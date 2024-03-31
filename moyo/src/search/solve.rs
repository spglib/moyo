use std::collections::BTreeMap;

use kiddo::{KdTree, SquaredEuclidean};
use nalgebra::Vector3;

use crate::base::{AtomicSpecie, Cell, Lattice, Permutation, Position, Rotation, Translation};

pub struct PeriodicKdTree {
    num_sites: usize,
    lattice: Lattice,
    kdtree: KdTree<f64, 3>,
}

#[derive(Debug)]
pub struct PeriodicNeighbor {
    pub index: usize,
    pub image: [i32; 3],
    pub distance: f64,
}

impl PeriodicKdTree {
    /// Construct a periodic kd-tree from the given **Minkowski-reduced** cell.
    pub fn new(reduced_cell: &Cell) -> Self {
        let mut entries = vec![];
        // Append entries by (offset, position) order
        for nz in -1..=1 {
            for ny in -1..=1 {
                for nx in -1..=1 {
                    let offset = Vector3::new(nx as f64, ny as f64, nz as f64);
                    for position in reduced_cell.positions.iter() {
                        let cart_coords =
                            reduced_cell.lattice.cartesian_coords(&(position + offset));
                        entries.push([cart_coords.x, cart_coords.y, cart_coords.z]);
                    }
                }
            }
        }
        Self {
            num_sites: reduced_cell.num_atoms(),
            lattice: reduced_cell.lattice.clone(),
            kdtree: (&entries).into(),
        }
    }

    /// Return the nearest neighbor within the given radius if exists.
    pub fn nearest(&self, position: &Position, radius: f64) -> Option<PeriodicNeighbor> {
        let cart_coords = self.lattice.cartesian_coords(position);
        let within = self.kdtree.nearest_n_within::<SquaredEuclidean>(
            &[cart_coords.x, cart_coords.y, cart_coords.z],
            radius.powi(2), // squared distance for KdTree
            1,
            false,
        );
        if within.len() == 0 {
            return None;
        }

        let mut item = within[0].item as usize;
        let mut remainders = vec![];
        for d in [self.num_sites, 3, 3, 3] {
            remainders.push(item % d);
            item /= d;
        }

        Some(PeriodicNeighbor {
            index: remainders[0],
            image: [
                remainders[1] as i32 - 1,
                remainders[2] as i32 - 1,
                remainders[3] as i32 - 1,
            ],
            distance: within[0].distance.sqrt(),
        })
    }
}

/// Choose atomic specie with the smallest occurrence
pub fn pivot_site_indices(numbers: &[AtomicSpecie]) -> Vec<usize> {
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
/// Search permutation such that new_positions[i] = reduced_cell.positions[permutation[i]].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation[i]-th site.
/// This function takes O(num_atoms * log(num_atoms)) time.
pub fn solve_correspondence(
    pkdtree: &PeriodicKdTree,
    reduced_cell: &Cell,
    new_positions: &[Position],
    symprec: f64,
) -> Option<Permutation> {
    let num_atoms = pkdtree.num_sites;
    let mut mapping = vec![0; num_atoms];
    let mut visited = vec![false; num_atoms];

    for i in 0..num_atoms {
        let neighbor = pkdtree.nearest(&new_positions[i], symprec)?;
        let j = neighbor.index;
        if visited[j] || reduced_cell.numbers[i] != reduced_cell.numbers[j] {
            return None;
        }
        if neighbor.distance > symprec {
            return None;
        }

        mapping[i] = j;
        visited[j] = true;
    }

    if visited.iter().all(|&v| v) {
        Some(Permutation::new(mapping))
    } else {
        None
    }
}

/// Return correspondence between the input and acted positions.
/// Assume that reduced_cell is Minkowski reduced and symprec is sufficiently small for Babai's algorithm.
/// Search permutation such that new_positions[i] = reduced_cell.positions[permutation[i]].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation[i]-th site.
/// This function takes O(num_atoms^2) time.
#[allow(clippy::needless_range_loop)]
pub fn solve_correspondence_naive(
    reduced_cell: &Cell,
    new_positions: &[Position],
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
            frac_displacement -= frac_displacement.map(|e| e.round()); // in [-0.5, 0.5]
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
    let distance = (0..num_atoms)
        .map(|i| {
            let mut frac_displacement = rotation.map(|e| e as f64) * reduced_cell.positions[i]
                + translation
                - reduced_cell.positions[permutation.apply(i)];
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
    use nalgebra::{Matrix3, Vector3};

    use crate::base::{Cell, Lattice, Permutation, Rotation, Translation};

    use super::{
        pivot_site_indices, solve_correspondence, solve_correspondence_naive,
        symmetrize_translation_from_permutation, PeriodicKdTree,
    };

    #[test]
    fn test_periodic_kdtree() {
        let pkdtree = PeriodicKdTree::new(&Cell::new(
            Lattice::new(Matrix3::identity()),
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.5, 0.0, 0.5),
                Vector3::new(0.5, 0.5, 0.0),
            ],
            vec![0, 0, 0, 0],
        ));

        {
            let neighbor = pkdtree.nearest(&Vector3::new(0.0, 0.0, 0.0), 1e-4).unwrap();
            assert!(neighbor.index == 0);
            assert!(neighbor.image == [0, 0, 0]);
        }

        {
            let neighbor = pkdtree.nearest(&Vector3::new(1.0, 0.5, 0.5), 1e-4).unwrap();
            assert!(neighbor.index == 1);
            assert!(neighbor.image == [1, 0, 0]);
        }

        {
            let neighbor = pkdtree
                .nearest(&Vector3::new(1.5, -0.0, -0.5), 1e-4)
                .unwrap();
            assert!(neighbor.index == 2);
            assert!(neighbor.image == [1, 0, -1]);
        }
    }

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
            Lattice::new(Matrix3::identity()),
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.5, 0.0, 0.5),
                Vector3::new(0.5, 0.5, 0.0),
            ],
            vec![0, 0, 0, 0],
        );
        let symprec = 1e-4;
        let pkdtree = PeriodicKdTree::new(&reduced_cell);

        {
            // Translation::new(0.0, 0.5, 0.5);
            let new_positions = vec![
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.0, 1.0, 1.0),
                Vector3::new(0.5, 0.5, 1.0),
                Vector3::new(0.5, 1.0, 0.5),
            ];
            let expect = Permutation::new(vec![1, 0, 3, 2]);

            let actual_naive =
                solve_correspondence_naive(&reduced_cell, &new_positions, symprec).unwrap();
            assert_eq!(actual_naive, expect);

            let actual_kdtree =
                solve_correspondence(&pkdtree, &reduced_cell, &new_positions, symprec).unwrap();
            assert_eq!(actual_kdtree, expect);
        }
        {
            // Translation::new(0.0, 0.5, 0.5 - 2 * symprec);
            let new_positions = vec![
                Vector3::new(0.0, 0.5, 0.5),
                Vector3::new(0.0, 1.0, 1.0 - 2.0 * symprec),
                Vector3::new(0.5, 0.5, 1.0),
                Vector3::new(0.5, 1.0, 0.5),
            ];

            let actual_naive = solve_correspondence_naive(&reduced_cell, &new_positions, symprec);
            assert_eq!(actual_naive, None);

            let actual_kdtree =
                solve_correspondence(&pkdtree, &reduced_cell, &new_positions, symprec);
            assert_eq!(actual_kdtree, None);
        }
    }

    #[test]
    fn test_symmetrize_translation_from_permutation() {
        // Conventional fcc
        let symprec = 1e-2;
        let distorted_reduced_cell = Cell::new(
            Lattice::new(Matrix3::identity()),
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
