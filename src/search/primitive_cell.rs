use std::collections::BTreeMap;

use nalgebra::{Dyn, Matrix3, OMatrix, Vector3, U3};

use super::solve::{
    pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
};
use crate::base::{
    orbits_from_permutations, Cell, Lattice, Linear, MoyoError, Permutation, Position, Rotation,
    Translation, UnimodularTransformation, EPS,
};
use crate::math::HNF;

#[derive(Debug)]
pub struct PrimitiveCell {
    /// Primitive cell
    pub cell: Cell,
    /// Transformation matrix from the primitive cell to the input cell
    pub linear: Linear,
    /// Mapping from sites of the input cell to those of the primitive cell (many-to-one)
    pub site_mapping: Vec<usize>,
    /// Translations in the **input** cell
    pub translations: Vec<Translation>,
    /// Permutations induced by translations in the input cell
    pub permutations: Vec<Permutation>,
}

impl PrimitiveCell {
    /// Return primitive cell and transformation matrix from the primitive cell to the input cell
    /// Possible replacements for spglib/src/primitive.h::prm_get_primitive
    pub fn new(cell: &Cell, symprec: f64) -> Result<Self, MoyoError> {
        // cell.lattice.basis * reduced_trans_mat = reduced_cell.lattice.basis
        let (reduced_lattice, reduced_trans_mat) = cell.lattice.minkowski_reduce()?;
        let reduced_cell =
            UnimodularTransformation::from_linear(reduced_trans_mat).transform_cell(cell);

        // Check if symprec is sufficiently small
        let minimum_basis_norm = reduced_lattice
            .basis
            .column_iter()
            .map(|v| v.norm())
            .fold(f64::INFINITY, f64::min);
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > minimum_basis_norm / 2.0 {
            return Err(MoyoError::TooSmallSymprecError);
        }

        // Try possible translations: overlap the `src`the site to the `dst`th site
        let pivot_site_indices = pivot_site_indices(&reduced_cell.numbers);
        let mut permutations_translations_tmp = vec![];
        let src = pivot_site_indices[0];
        for dst in pivot_site_indices.iter() {
            let translation = reduced_cell.positions[*dst] - reduced_cell.positions[src];
            let new_positions: Vec<Position> = reduced_cell
                .positions
                .iter()
                .map(|pos| pos + translation)
                .collect();

            // Because the translation may not be optimal to minimize distance between input and acted positions,
            // use a larger symprec (diameter of a Ball) for finding correspondence
            if let Some(permutation) =
                solve_correspondence(&reduced_cell, &new_positions, rough_symprec)
            {
                permutations_translations_tmp.push((permutation, translation));
            }
        }
        assert!(!permutations_translations_tmp.is_empty());

        // Purify translations by permutations
        let mut translations = vec![];
        let mut permutations = vec![];
        for (permutation, rough_translation) in permutations_translations_tmp.iter() {
            let (translation, distance) = symmetrize_translation_from_permutation(
                &reduced_cell,
                permutation,
                &Rotation::identity(),
                rough_translation,
            );
            if distance < symprec {
                translations.push(translation);
                permutations.push(permutation.clone());
            }
        }

        let size = translations.len() as i32;
        assert!(size > 0);
        if reduced_cell.num_atoms() % (size as usize) != 0 {
            return Err(MoyoError::PrimitiveCellError);
        }

        // Recover a transformation matrix from primitive to input cell
        let mut columns: Vec<Vector3<i32>> = vec![
            Vector3::new(size, 0, 0),
            Vector3::new(0, size, 0),
            Vector3::new(0, 0, size),
        ];
        for translation in translations.iter() {
            columns.push((translation * (size as f64)).map(|e| e.round() as i32));
        }
        let hnf = HNF::new(&OMatrix::<i32, U3, Dyn>::from_columns(&columns));
        let trans_mat_inv =
            Matrix3::<i32>::from_columns(&[hnf.h.column(0), hnf.h.column(1), hnf.h.column(2)])
                .map(|e| e as f64)
                / (size as f64);
        let trans_mat = trans_mat_inv
            .try_inverse()
            .ok_or(MoyoError::PrimitiveCellError)?;
        if relative_ne!(trans_mat.determinant(), size as f64, epsilon = EPS) {
            return Err(MoyoError::PrimitiveCellError);
        }

        // Primitive cell
        let (primitive_cell, site_mapping) = primitive_cell_from_transformation(
            &reduced_cell,
            &trans_mat,
            &translations,
            &permutations,
        )
        .ok_or(MoyoError::PrimitiveCellError)?;
        let (_, prim_trans_mat) = primitive_cell.lattice.minkowski_reduce()?;
        let reduced_prim_cell =
            UnimodularTransformation::from_linear(prim_trans_mat).transform_cell(&primitive_cell);

        // (input cell)
        //    -[reduced_trans_mat]-> (reduced cell)
        //    <-[trans_mat]- (primitive cell)
        //    -[prim_trans_mat]-> (reduced primitive cell)
        let inv_prim_trans_mat = prim_trans_mat
            .map(|e| e as f64)
            .try_inverse()
            .unwrap()
            .map(|e| e.round() as i32);
        let inv_reduced_trans_mat = reduced_trans_mat
            .map(|e| e as f64)
            .try_inverse()
            .ok_or(MoyoError::PrimitiveCellError)?;
        Ok(Self {
            cell: reduced_prim_cell,
            linear: inv_prim_trans_mat.map(|e| e as f64) * trans_mat * inv_reduced_trans_mat,
            site_mapping,
            translations: translations
                .iter()
                .map(|translation| reduced_trans_mat.map(|e| e as f64) * translation)
                .collect(),
            permutations,
        })
    }
}

fn primitive_cell_from_transformation(
    cell: &Cell,
    trans_mat: &Linear,
    translations: &Vec<Translation>,
    permutations: &[Permutation],
) -> Option<(Cell, Vec<usize>)> {
    let trans_mat_inv = trans_mat.try_inverse().unwrap();
    let new_lattice = Lattice::new(cell.lattice.basis * trans_mat_inv);

    let num_atoms = cell.num_atoms();
    let orbits = orbits_from_permutations(num_atoms, permutations);
    let representatives = (0..num_atoms)
        .filter(|&i| orbits[i] == i)
        .collect::<Vec<_>>();

    // Eq. (25) of https://arxiv.org/pdf/2211.15008.pdf
    let mut new_positions = vec![Vector3::zeros(); representatives.len()];
    let mut new_numbers = vec![0; representatives.len()];
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();
    for orbit in representatives.iter() {
        let mut acc = Vector3::zeros();
        for (inv_perm, translation) in inverse_permutations.iter().zip(translations.iter()) {
            let mut frac_displacements =
                cell.positions[inv_perm.apply(*orbit)] + translation - cell.positions[*orbit];
            frac_displacements -= frac_displacements.map(|e| e.round());
            acc += frac_displacements;
        }
        new_positions[*orbit] = cell.positions[*orbit] + acc / (translations.len() as f64);
        new_numbers[*orbit] = cell.numbers[*orbit];
    }

    let primitive_cell = Cell::new(new_lattice, new_positions, new_numbers);
    let site_mapping = site_mapping_from_orbits(&orbits);
    Some((primitive_cell, site_mapping))
}

fn site_mapping_from_orbits(orbits: &[usize]) -> Vec<usize> {
    let mut mapping = BTreeMap::new();
    for ri in orbits.iter() {
        if mapping.contains_key(&ri) {
            continue;
        }
        mapping.insert(ri, mapping.len());
    }

    orbits.iter().map(|ri| *mapping.get(&ri).unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Vector3};

    use crate::base::{Cell, Lattice, Translation};

    use super::{site_mapping_from_orbits, PrimitiveCell};

    #[test]
    fn test_site_mapping_from_orbits() {
        let orbits = vec![0, 0, 2, 2, 0, 6];
        assert_eq!(site_mapping_from_orbits(&orbits), vec![0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_search_primitive_cell() {
        let symprec = 1e-4;

        // Conventional fcc
        {
            let cell = Cell::new(
                Lattice::new(matrix![
                    1.0, 0.0, 0.0;
                    0.0, 1.0, 0.0;
                    0.0, 0.0, 1.0;
                ]),
                vec![
                    Vector3::new(0.5 * symprec, 0.0, 0.0),
                    Vector3::new(0.0, 0.5, 0.5 + 0.5 * symprec),
                    Vector3::new(0.5, 0.0, 0.5),
                    Vector3::new(0.5, 0.5, 0.0),
                ],
                vec![0, 0, 0, 0],
            );

            let result = PrimitiveCell::new(&cell, symprec).unwrap();
            assert_eq!(result.site_mapping, vec![0, 0, 0, 0]);
            assert_relative_eq!(
                result.cell.positions[0],
                Vector3::zeros(),
                epsilon = symprec
            );
            assert_eq!(result.cell.numbers[0], 0);
        }

        // bcc in non-minkowski-reduced cell
        {
            let cell = Cell::new(
                Lattice::new(matrix![
                    1.0, 0.0, 0.0;
                    1.0, 1.0, 0.0;
                    0.0, 0.0, 1.0;
                ]),
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.5, 0.0, 0.5)],
                vec![0, 0],
            );
            let result = PrimitiveCell::new(&cell, symprec).unwrap();
            assert_eq!(result.site_mapping, vec![0, 0]);
            assert_eq!(result.cell.numbers[0], 0);
            assert_relative_eq!(result.translations[0], Translation::new(0.0, 0.0, 0.0));
            assert_relative_eq!(result.translations[1], Translation::new(0.5, 0.0, 0.5));
        }
    }
}
