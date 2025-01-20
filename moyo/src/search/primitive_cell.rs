use std::collections::BTreeMap;

use log::debug;
use nalgebra::{Dyn, Matrix3, OMatrix, Vector3, U3};

use super::solve::{
    pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
    PeriodicKdTree,
};
use crate::base::{
    orbits_from_permutations, Cell, Linear, MagneticCell, MagneticMoment, MoyoError, Permutation,
    Position, Rotation, Transformation, Translation, UnimodularTransformation, EPS,
};
use crate::math::HNF;

#[derive(Debug)]
pub struct PrimitiveCell {
    /// Primitive cell
    pub cell: Cell,
    /// Transformation matrix from the **primitive** cell to the input cell
    pub linear: Linear,
    /// Mapping from sites of the input cell to those of the primitive cell (many-to-one).
    /// The `i`th atom in the input cell is equivalent to the `site_mapping[i]`th atom in the **primitive** cell.
    pub site_mapping: Vec<usize>,
    /// Translations in the **input** cell
    pub translations: Vec<Translation>,
    /// Permutations induced by translations in the input cell.
    /// `translations[i]` moves the `k`th site to the `permutations[i].apply(k)`th site.
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
            .reduce(f64::min)
            .unwrap();
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > minimum_basis_norm / 2.0 {
            debug!(
                "symprec is too large compared to the basis vectors. Consider reducing symprec."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        // Try possible translations: overlap the `src`the site to the `dst`th site
        let pkdtree = PeriodicKdTree::new(&reduced_cell, rough_symprec);
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
            if let Some(permutation) = solve_correspondence(&pkdtree, &reduced_cell, &new_positions)
            {
                permutations_translations_tmp.push((permutation, translation));
            }
        }

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

        // Check number of translations
        let size = translations.len() as i32;
        if (size == 0) || (reduced_cell.num_atoms() % (size as usize) != 0) {
            debug!("Failed to properly find translations: {} translations in {} atoms. Consider increasing symprec.", size, reduced_cell.num_atoms());
            return Err(MoyoError::TooSmallToleranceError);
        }
        debug!("Found {} pure translations", size);

        // Recover a transformation matrix from primitive to input cell
        let trans_mat = if let Some(trans_mat) =
            transformation_matrix_from_translations(&translations)
        {
            trans_mat
        } else {
            debug!("Failed to find a transformation matrix for a primitive cell. Consider increasing symprec.");
            return Err(MoyoError::TooSmallToleranceError);
        };

        // Primitive cell
        let (primitive_cell, site_mapping, _) = primitive_cell_from_transformation(
            &reduced_cell,
            &trans_mat,
            &translations,
            &permutations,
        );
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
        let inv_reduced_trans_mat = reduced_trans_mat.map(|e| e as f64).try_inverse().unwrap();
        Ok(Self {
            cell: reduced_prim_cell,
            linear: ((inv_prim_trans_mat * trans_mat).map(|e| e as f64) * inv_reduced_trans_mat)
                .map(|e| e.round() as i32),
            site_mapping,
            translations: translations
                .iter()
                .map(|translation| reduced_trans_mat.map(|e| e as f64) * translation)
                .collect(),
            permutations,
        })
    }
}

#[derive(Debug)]
pub struct PrimitiveMagneticCell<M: MagneticMoment> {
    pub magnetic_cell: MagneticCell<M>,
    /// Transformation matrix from the **primitive** magnetic cell to the input magnetic cell
    pub linear: Linear,
    /// Mapping from sites of the input magnetic cell to those of the primitive magnetic cell (many-to-one).
    /// The `i`th atom in the input cell is equivalent to the `site_mapping[i]`th atom in the **primitive** magnetic cell.
    pub site_mapping: Vec<usize>,
    /// Translations in the **input** magnetic cell
    pub translations: Vec<Translation>,
    /// Permutations induced by translations in the input magnetic cell
    #[allow(dead_code)]
    pub permutations: Vec<Permutation>,
}

impl<M: MagneticMoment> PrimitiveMagneticCell<M> {
    pub fn new(
        magnetic_cell: &MagneticCell<M>,
        symprec: f64,
        mag_symprec: f64,
    ) -> Result<Self, MoyoError> {
        // Prepare candidate translations from nonmagnetic cell
        let prim_nonmagnetic_cell = PrimitiveCell::new(&magnetic_cell.cell, symprec)?;

        // Filter translations that keep magnetic moments
        let mut translations = vec![];
        let mut permutations = vec![];
        for (translation, permutation) in prim_nonmagnetic_cell
            .translations
            .iter()
            .zip(prim_nonmagnetic_cell.permutations.iter())
        {
            let new_magnetic_moments = (0..magnetic_cell.cell.num_atoms())
                .map(|i| magnetic_cell.magnetic_moments[permutation.apply(i)].clone())
                .collect::<Vec<_>>();
            let take = magnetic_cell
                .magnetic_moments
                .iter()
                .zip(new_magnetic_moments.iter())
                .all(|(m1, m2)| m1.is_close(m2, mag_symprec));
            if take {
                translations.push(*translation);
                permutations.push(permutation.clone());
            }
        }

        // Check number of translations
        let size = translations.len() as i32;
        if (size == 0) || (magnetic_cell.cell.num_atoms() % (size as usize) != 0) {
            debug!("Failed to properly find translations: {} translations in {} atoms. Consider increasing symprec.", size, magnetic_cell.cell.num_atoms());
            return Err(MoyoError::TooSmallToleranceError);
        }
        debug!("Found {} pure translations", size);

        // Recover a transformation matrix from primitive to input cell
        // trans_mat: prim_mag_cell -> magnetic_cell (input)
        let trans_mat = if let Some(trans_mat) =
            transformation_matrix_from_translations(&translations)
        {
            trans_mat
        } else {
            debug!("Failed to find a transformation matrix for a primitive cell. Consider increasing symprec.");
            return Err(MoyoError::TooSmallToleranceError);
        };

        // Primitive magnetic cell
        let (prim_mag_cell, site_mapping) = primitive_magnetic_cell_from_transformation(
            magnetic_cell,
            &trans_mat,
            &translations,
            &permutations,
        );
        let (_, prim_trans_mat) = prim_mag_cell.cell.lattice.minkowski_reduce()?;
        // prim_trans_mat: prim_mag_cell -> reduced_prim_mag_cell
        let reduced_prim_mag_cell = UnimodularTransformation::from_linear(prim_trans_mat)
            .transform_magnetic_cell(&prim_mag_cell);

        let inv_prim_trans_mat = prim_trans_mat
            .map(|e| e as f64)
            .try_inverse()
            .unwrap()
            .map(|e| e.round() as i32);

        Ok(Self {
            magnetic_cell: reduced_prim_mag_cell,
            linear: inv_prim_trans_mat * trans_mat,
            site_mapping,
            translations,
            permutations,
        })
    }
}

fn transformation_matrix_from_translations(translations: &[Translation]) -> Option<Linear> {
    let size = translations.len() as i32;
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
        .unwrap()
        .map(|e| e.round() as i32);

    if relative_ne!(
        trans_mat.map(|e| e as f64).determinant(),
        size as f64,
        epsilon = EPS
    ) {
        return None;
    }
    Some(trans_mat)
}

/// Transform `cell` to a primitive cell by inverse of `trans_mat`
fn primitive_cell_from_transformation(
    cell: &Cell,
    trans_mat: &Linear,
    translations: &[Translation],
    permutations: &[Permutation],
) -> (Cell, Vec<usize>, Vec<usize>) {
    let new_lattice =
        Transformation::from_linear(*trans_mat).inverse_transform_lattice(&cell.lattice);

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
    for (i, &orbit_i) in representatives.iter().enumerate() {
        let mut acc = Vector3::zeros();
        for (inv_perm, translation) in inverse_permutations.iter().zip(translations.iter()) {
            let mut frac_displacements =
                cell.positions[inv_perm.apply(orbit_i)] + translation - cell.positions[orbit_i];
            frac_displacements -= frac_displacements.map(|e| e.round()); // in [-0.5, 0.5]
            acc += frac_displacements;
        }
        new_positions[i] = trans_mat.map(|e| e as f64)
            * (cell.positions[orbit_i] + acc / (translations.len() as f64));
        new_numbers[i] = cell.numbers[orbit_i];
    }

    let primitive_cell = Cell::new(new_lattice, new_positions, new_numbers);
    let site_mapping = site_mapping_from_orbits(&orbits);
    (primitive_cell, site_mapping, representatives)
}

fn primitive_magnetic_cell_from_transformation<M: MagneticMoment>(
    magnetic_cell: &MagneticCell<M>,
    trans_mat: &Linear,
    translations: &[Translation],
    permutations: &[Permutation],
) -> (MagneticCell<M>, Vec<usize>) {
    let (primitive_cell, site_mapping, representatives) = primitive_cell_from_transformation(
        &magnetic_cell.cell,
        trans_mat,
        translations,
        permutations,
    );
    let new_magnetic_moments = representatives
        .iter()
        .map(|&i| magnetic_cell.magnetic_moments[i].clone())
        .collect::<Vec<_>>();
    let primitive_magnetic_cell = MagneticCell::from_cell(primitive_cell, new_magnetic_moments);
    (primitive_magnetic_cell, site_mapping)
}

fn site_mapping_from_orbits(orbits: &[usize]) -> Vec<usize> {
    let mut mapping = BTreeMap::new();
    let mut count = 0;
    for ri in orbits.iter() {
        mapping.entry(ri).or_insert_with(|| {
            let value = count;
            count += 1;
            value
        });
    }

    orbits.iter().map(|ri| *mapping.get(&ri).unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Matrix3, Vector3};

    use crate::base::{
        Cell, Collinear, Lattice, MagneticCell, MagneticMoment, Transformation, Translation,
    };

    use super::{
        site_mapping_from_orbits, transformation_matrix_from_translations, PrimitiveCell,
        PrimitiveMagneticCell,
    };

    #[test]
    fn test_transformation_matrix_from_translations() {
        // bcc
        let translations = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.5, 0.5, 0.5)];
        transformation_matrix_from_translations(&translations).unwrap();
    }

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
                Lattice::new(Matrix3::identity()),
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
                    1.0, 1.0, 0.0;
                    0.0, 1.0, 0.0;
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

    #[test]
    fn test_rhombohedral_lattice() {
        let a = 4.0;
        let b = 7.0;
        let rhombohedral_lattice = Lattice::new(matrix![
            3.0_f64.sqrt() / 2.0 * a, 0.5 * a, b;
            -3.0_f64.sqrt() / 2.0 * a, 0.5 * a, b;
            0.0, -a, b;
        ]);
        let trans_mat = matrix![
            1, 0, 1;
            -1, 1, 1;
            0, -1, 1;
        ];

        let lattice =
            Lattice::new((rhombohedral_lattice.basis * trans_mat.map(|e| e as f64)).transpose());
        let cell = Cell::new(
            lattice,
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                Vector3::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
                Vector3::new(0.0, 0.0, 0.1),
                Vector3::new(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 + 0.1),
                Vector3::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0 + 0.1),
            ],
            vec![0, 0, 0, 0, 0, 0],
        );
        let symprec = 1e-4;
        let prim_cell = PrimitiveCell::new(&cell, symprec).unwrap();

        assert_relative_eq!(
            prim_cell.cell.lattice.basis * prim_cell.linear.map(|e| e as f64),
            cell.lattice.basis,
            epsilon = 1e-8
        );
    }

    #[test]
    fn test_magnetic_primitive_cell() {
        let symprec = 1e-4;
        let mag_symprec = 1e-4;

        let magnetic_cell = MagneticCell::new(
            Lattice::new(Matrix3::identity()),
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.5, 0.0, 0.5),
                Vector3::new(0.5, 0.5, 0.0),
            ],
            vec![0, 0, 0, 0],
            vec![
                Collinear(1.0),
                Collinear(1.0),
                Collinear(-1.0),
                Collinear(-1.0),
            ],
        );
        let result = PrimitiveMagneticCell::new(&magnetic_cell, symprec, mag_symprec).unwrap();
        let prim_mag_cell = &result.magnetic_cell;
        assert_eq!(prim_mag_cell.cell.num_atoms(), 2);
        assert!(prim_mag_cell.magnetic_moments[0].is_close(&Collinear(1.0), 1e-8));
        assert!(prim_mag_cell.magnetic_moments[1].is_close(&Collinear(-1.0), 1e-8));

        assert_eq!(result.site_mapping, vec![0, 0, 1, 1]);
        assert_eq!(result.permutations.len(), 2);

        assert_eq!(result.translations.len(), 2);
        assert_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(result.translations[1], Vector3::new(0.0, 0.5, 0.5));

        assert_eq!(result.linear.map(|e| e as f64).determinant() as i32, 2);
        let (recovered_mag_cell, _) =
            Transformation::from_linear(result.linear).transform_magnetic_cell(prim_mag_cell);
        assert_eq!(recovered_mag_cell.cell.num_atoms(), 4);
        assert_relative_eq!(
            recovered_mag_cell.cell.lattice.basis,
            magnetic_cell.cell.lattice.basis,
        );
    }
}
