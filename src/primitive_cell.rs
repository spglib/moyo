use nalgebra::{Dyn, OMatrix, Vector3, U3};
use std::collections::HashMap;

use crate::cell::{Cell, Position};
use crate::error::MoyoError;
use crate::math::hnf::hnf;
use crate::operation::{CellWithOperations, Permutation, Rotation, Translation};
use crate::transformation::Transformation;

/// Return primitive cell and transformation matrix from the primitive cell to the input cell
pub fn search_primitive_cell(
    cell: &Cell,
    symprec: f64,
) -> Result<(CellWithOperations, Transformation), MoyoError> {
    let (reduced_lattice, reduced_trans) = cell.lattice.minkowski_reduce()?;
    let reduced_cell = cell.transform(&reduced_trans);

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
    let mut permutations_translations = vec![];
    for (permutation, rough_translation) in permutations_translations_tmp.iter() {
        let (translation, distance) = symmetrize_translation_from_permutation(
            &reduced_cell,
            permutation,
            &Rotation::identity(),
            rough_translation,
        );
        if distance < symprec {
            permutations_translations.push((permutation, translation));
        }
    }

    let size = permutations_translations.len() as i32;
    assert!(size > 0);
    if reduced_cell.num_atoms % (size as usize) != 0 {
        return Err(MoyoError::PrimitiveCellSearchError);
    }

    // Recover a transformation matrix from translations
    let mut columns: Vec<Vector3<i32>> = vec![
        Vector3::new(size, 0, 0),
        Vector3::new(0, size, 0),
        Vector3::new(0, 0, size),
    ];
    for (_, translation) in permutations_translations {
        columns.push((translation * (size as f64)).map(|e| e.round() as i32));
    }
    let hnf = hnf(&OMatrix::<i32, U3, Dyn>::from_columns(&columns));

    unimplemented!();
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
    let num_atoms = reduced_cell.num_atoms;

    let mut permutation = vec![0; num_atoms];
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
                permutation[i] = j;
                visited[j] = true;
                break;
            }
        }
    }

    if visited.iter().all(|&v| v) {
        Some(permutation)
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
    let num_atoms = reduced_cell.num_atoms;

    let translation = (0..num_atoms)
        .map(|i| {
            let mut frac_displacement = reduced_cell.positions[permutation[i]]
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
            let mut frac_displacement = reduced_cell.positions[permutation[i]] - new_positions[i];
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

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Vector3};

    use crate::cell::Cell;
    use crate::lattice::Lattice;
    use crate::operation::{Rotation, Translation};

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
            let expect = vec![1, 0, 3, 2];
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

        let permutation = vec![1, 0, 3, 2];
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
        // Conventional fcc
        let cell = Cell::new(
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

        let _ = search_primitive_cell(&cell, symprec);
    }
}
