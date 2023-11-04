use std::collections::HashMap;
use std::collections::HashSet;

use nalgebra::{Dyn, Matrix3, OMatrix, Vector3, U3};
use union_find::{QuickFindUf, UnionByRank, UnionFind};

use crate::base::cell::{Cell, Position, SiteMapping};
use crate::base::error::MoyoError;
use crate::base::lattice::Lattice;
use crate::base::operation::{Permutation, Rotation, Translation};
use crate::base::transformation::{OriginShift, Transformation, TransformationMatrix};
use crate::math::hnf::hnf;

const EPS: f64 = 1e-8;

#[derive(Debug)]
pub struct PrimitiveCellSearchResult {
    pub primitive_cell: Cell,
    /// Transformation matrix from the primitive cell to the input cell
    pub trans_mat: TransformationMatrix,
    /// Mapping from sites of the input cell to those of the primitive cell (many-to-one)
    pub site_mapping: SiteMapping,
    /// Translations in the **input** cell
    pub translations: Vec<Translation>,
    /// Permutations induced by translations in the input cell
    pub permutations: Vec<Permutation>,
}

impl PrimitiveCellSearchResult {
    pub fn new(
        primitive_cell: Cell,
        trans_mat: TransformationMatrix,
        site_mapping: SiteMapping,
        translations: Vec<Translation>,
        permutations: Vec<Permutation>,
    ) -> Self {
        Self {
            primitive_cell,
            trans_mat,
            site_mapping,
            translations,
            permutations,
        }
    }
}

/// Return primitive cell and transformation matrix from the primitive cell to the input cell
pub fn search_primitive_cell(
    cell: &Cell,
    symprec: f64,
) -> Result<PrimitiveCellSearchResult, MoyoError> {
    // cell.lattice.basis * reduced_trans_mat = reduced_cell.lattice.basis
    let (reduced_lattice, reduced_trans_mat) = cell.lattice.minkowski_reduce()?;
    let reduced_cell = cell.transform(&Transformation::new(
        reduced_trans_mat,
        OriginShift::zeros(),
    ));

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

    // Choose atomic specie with the smallest occurrence
    let mut counter = HashMap::new();
    for number in reduced_cell.numbers.iter() {
        let count = counter.entry(number).or_insert(0);
        *count += 1;
    }
    let pivot_atomic_specie = *counter.iter().min_by_key(|(_, count)| *count).unwrap().0;
    let pivot_site_indices = reduced_cell
        .numbers
        .iter()
        .enumerate()
        .filter(|(_, number)| *number == pivot_atomic_specie)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();

    // Try possible translations: overlap the `src`the site to the `dst`th site
    // TODO: this part takes O(num_atoms^3)
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
    assert!(permutations_translations_tmp.len() > 0);

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
        return Err(MoyoError::PrimitiveCellSearchError);
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
    let hnf = hnf(&OMatrix::<i32, U3, Dyn>::from_columns(&columns));
    let trans_mat_inv =
        Matrix3::<i32>::from_columns(&[hnf.h.column(0), hnf.h.column(1), hnf.h.column(2)])
            .map(|e| e as f64)
            / (size as f64);
    let trans_mat = trans_mat_inv
        .try_inverse()
        .ok_or(MoyoError::PrimitiveCellSearchError)?
        .map(|e| e.round() as i32);
    if relative_ne!(
        trans_mat.map(|e| e as f64).determinant(),
        size as f64,
        epsilon = EPS
    ) {
        return Err(MoyoError::PrimitiveCellSearchError);
    }

    // Primitive cell
    let (primitive_cell, site_mapping) =
        primitive_cell_from_transformation(&reduced_cell, &trans_mat, &translations, &permutations)
            .ok_or(MoyoError::PrimitiveCellSearchError)?;

    // (input cell) --(reduced_trans_mat)--> (Minkowski reduced cell) <--(trans_mat)-- (primitive cell)
    let inv_reduced_trans_mat = reduced_trans_mat
        .map(|e| e as f64)
        .try_inverse()
        .unwrap()
        .map(|e| e as i32);
    Ok(PrimitiveCellSearchResult::new(
        primitive_cell,
        trans_mat * inv_reduced_trans_mat,
        site_mapping,
        translations
            .iter()
            .map(|translation| reduced_trans_mat.map(|e| e as f64) * translation)
            .collect(),
        permutations,
    ))
}

/// Return correspondence between the input and acted positions.
/// Assume that reduced_cell is Minkowski reduced and symprec is sufficiently small for Babai's algorithm.
/// Search permutation such that new_positions[i] = reduced_cell.positions[permutation[i]].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation[i]-th site.
/// Be careful that the current implementation takes O(num_atoms^2) time!
fn solve_correspondence(
    reduced_cell: &Cell,
    new_positions: &Vec<Position>,
    symprec: f64,
) -> Option<Permutation> {
    let num_atoms = reduced_cell.num_atoms();

    let mut mapping = vec![0; num_atoms];
    let mut visited = vec![false; num_atoms];

    for i in 0..num_atoms {
        for j in 0..num_atoms {
            if visited[j] || reduced_cell.numbers[i] != reduced_cell.numbers[j] {
                continue;
            }

            let mut frac_displacement = new_positions[j] - reduced_cell.positions[i];
            frac_displacement -= frac_displacement.map(|e| e.round());
            let distance = reduced_cell
                .lattice
                .cartesian_coords(&frac_displacement)
                .norm();
            if distance < symprec {
                mapping[i] = j;
                visited[j] = true;
                break;
            }
        }
    }

    if visited.iter().all(|&v| v) {
        Some(Permutation::new(mapping))
    } else {
        None
    }
}

fn symmetrize_translation_from_permutation(
    reduced_cell: &Cell,
    permutation: &Permutation,
    rotation: &Rotation,
    rough_translation: &Translation,
) -> (Translation, f64) {
    // argmin_{t} sum_{i} | pbc(rotation * positions[i] + t - positions[permutation[i]]) |^2
    //   = 1/num_atoms * sum_{i} pbc(positions[permutation[i]] - rotation * positions[i])
    let num_atoms = reduced_cell.num_atoms();

    let translation = (0..num_atoms)
        .map(|i| {
            let mut frac_displacement = reduced_cell.positions[permutation.apply(i)]
                - rotation.map(|e| e as f64) * reduced_cell.positions[i];

            // To avoid rounding error, we first subtract rough translation. Then, the remainder should be almost zeros.
            frac_displacement -= rough_translation;
            frac_displacement -= frac_displacement.map(|e| e.round());
            frac_displacement += rough_translation;

            frac_displacement
        })
        .sum::<Vector3<_>>()
        / (num_atoms as f64);
    let new_positions = reduced_cell
        .positions
        .iter()
        .map(|pos| rotation.map(|e| e as f64) * pos + translation)
        .collect::<Vec<_>>();
    let distance = (0..num_atoms)
        .map(|i| {
            let mut frac_displacement =
                reduced_cell.positions[permutation.apply(i)] - new_positions[i];
            frac_displacement -= frac_displacement.map(|e| e.round());
            reduced_cell
                .lattice
                .cartesian_coords(&frac_displacement)
                .norm()
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    (translation, distance)
}

fn primitive_cell_from_transformation(
    cell: &Cell,
    trans_mat: &TransformationMatrix,
    translations: &Vec<Translation>,
    permutations: &Vec<Permutation>,
) -> Option<(Cell, SiteMapping)> {
    let trans_mat_inv = trans_mat.map(|e| e as f64).try_inverse().unwrap();
    let new_lattice = Lattice::new(cell.lattice.basis * trans_mat_inv);

    // Create site mapping from union-find tree
    let num_atoms = cell.num_atoms();
    let mut uf = QuickFindUf::<UnionByRank>::new(num_atoms);
    for permutation in permutations.iter() {
        for i in 0..num_atoms {
            uf.union(i, permutation.apply(i));
        }
    }
    let site_mapping: SiteMapping = (0..num_atoms).map(|i| uf.find(i)).collect();
    let mut orbits_set = HashSet::new();
    for i in 0..num_atoms {
        orbits_set.insert(uf.find(i));
    }
    let orbits = orbits_set.into_iter().collect::<Vec<_>>();

    // Eq. (25) of https://arxiv.org/pdf/2211.15008.pdf
    let mut new_positions = vec![Vector3::zeros(); orbits.len()];
    let mut new_numbers = vec![0; orbits.len()];
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();
    for orbit in orbits.iter() {
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
    Some((primitive_cell, site_mapping))
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Vector3};

    use crate::base::cell::Cell;
    use crate::base::lattice::Lattice;
    use crate::base::operation::{Permutation, Rotation, Translation};

    use super::{
        search_primitive_cell, solve_correspondence, symmetrize_translation_from_permutation,
    };

    #[test]
    fn test_solve_correspondence() {
        // Conventional fcc
        let reduced_cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 1.0;
            ]),
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.5, 0.0, 0.5),
                Vector3::new(0.5, 0.5, 0.0),
            ],
            vec![0, 0, 0, 0],
        );
        let symprec = 1e-4;

        {
            // Translation::new(0.0, 0.5, 0.5);
            let new_positions = vec![
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.0, 1.0, 1.0),
                Vector3::new(0.5, 0.5, 1.0),
                Vector3::new(0.5, 1.0, 0.5),
            ];
            let actual = solve_correspondence(&reduced_cell, &new_positions, symprec).unwrap();
            let expect = Permutation::new(vec![1, 0, 3, 2]);
            assert_eq!(actual, expect);
        }
        {
            // Translation::new(0.0, 0.5, 0.5 - 2 * symprec);
            let new_positions = vec![
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.0, 1.0, 1.0 - 2.0 * symprec),
                Vector3::new(0.5, 0.5, 1.0),
                Vector3::new(0.5, 1.0, 0.5),
            ];
            let actual = solve_correspondence(&reduced_cell, &new_positions, symprec);
            assert_eq!(actual, None);
        }
    }

    #[test]
    fn test_symmetrize_translation_from_permutation() {
        // Conventional fcc
        let symprec = 1e-2;
        let distorted_reduced_cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 1.0;
            ]),
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.5, 0.5 + 0.5 * symprec),
                Vector3::new(0.5, 0.0, 0.5),
                Vector3::new(0.5, 0.5, 0.0),
            ],
            vec![0, 0, 0, 0],
        );

        let permutation = Permutation::new(vec![1, 0, 3, 2]);
        let (actual, distance) = symmetrize_translation_from_permutation(
            &distorted_reduced_cell,
            &permutation,
            &Rotation::identity(),
            &Translation::new(0.0, 0.5, 0.5 + 0.5 * symprec),
        );
        let expect = Translation::new(0.0, 0.5, 0.5);
        assert_relative_eq!(actual, expect);
        assert_relative_eq!(distance, 0.5 * symprec);
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

            let result = search_primitive_cell(&cell, symprec).unwrap();
            assert_eq!(result.site_mapping, vec![0, 0, 0, 0]);
            assert_relative_eq!(
                result.primitive_cell.positions[0],
                Vector3::zeros(),
                epsilon = symprec
            );
            assert_eq!(result.primitive_cell.numbers[0], 0);
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
            let result = search_primitive_cell(&cell, symprec).unwrap();
            assert_eq!(result.site_mapping, vec![0, 0]);
            assert_eq!(result.primitive_cell.numbers[0], 0);
            assert_relative_eq!(result.translations[0], Translation::new(0.0, 0.0, 0.0));
            assert_relative_eq!(result.translations[1], Translation::new(0.5, 0.0, 0.5));
        }
    }
}
