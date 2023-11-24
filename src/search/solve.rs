use std::collections::BTreeMap;

use nalgebra::Vector3;

use crate::base::cell::{AtomicSpecie, Cell, Position};
use crate::base::operation::{Permutation, Rotation, Translation};

/// Choose atomic specie with the smallest occurrence
pub fn pivot_site_indices(numbers: &Vec<AtomicSpecie>) -> Vec<usize> {
    let mut counter = BTreeMap::new();
    for number in numbers.iter() {
        let count = counter.entry(number).or_insert(0);
        *count += 1;
    }
    let pivot_atomic_specie = *counter.iter().min_by_key(|(_, count)| *count).unwrap().0;
    let pivot_site_indices = numbers
        .iter()
        .enumerate()
        .filter(|(_, number)| *number == pivot_atomic_specie)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    pivot_site_indices
}

/// Return correspondence between the input and acted positions.
/// Assume that reduced_cell is Minkowski reduced and symprec is sufficiently small for Babai's algorithm.
/// Search permutation such that new_positions[i] = reduced_cell.positions[permutation[i]].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation[i]-th site.
/// TODO: Be careful that the current implementation takes O(num_atoms^2) time!
pub fn solve_correspondence(
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

            let mut frac_displacement = reduced_cell.positions[j] - new_positions[i];
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

pub fn symmetrize_translation_from_permutation(
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

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Vector3};

    use crate::base::cell::Cell;
    use crate::base::lattice::Lattice;
    use crate::base::operation::{Permutation, Rotation, Translation};

    use super::{
        pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
    };

    #[test]
    fn test_pivot_site_indices() {
        let numbers = vec![0, 1, 1, 1, 2, 0, 2, 2];
        let actual = pivot_site_indices(&numbers);
        let expect = vec![0, 5];
        assert_eq!(actual, expect);
    }

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
}
