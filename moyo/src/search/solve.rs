use std::collections::BTreeMap;

use itertools::iproduct;
use nalgebra::Vector3;
use rustc_hash::FxHashMap;

use crate::base::{AtomicSpecie, Cell, Lattice, Permutation, Position, Rotation, Translation};

/// Nearest-neighbor search under periodic boundary conditions, backed by a
/// uniform bin grid over Cartesian coordinates. Each point is registered in
/// every bin its ball of radius `symprec` overlaps -- at most eight bins,
/// because the bin size is at least `2 * symprec`. A query then only
/// inspects the single bin containing the query position: any point within
/// `symprec` of the query has registered itself there.
#[doc(hidden)]
pub struct PeriodicNeighborSearch {
    num_sites: usize,
    lattice: Lattice,
    bin_size: f64,
    /// Bin index -> entries registered to the bin (indices into `cart_coords`
    /// and `indices`)
    bins: FxHashMap<(i64, i64, i64), Vec<usize>>,
    cart_coords: Vec<Vector3<f64>>,
    /// Entry -> site index in the original cell
    indices: Vec<usize>,
    symprec: f64,
}

#[doc(hidden)]
#[derive(Debug)]
pub struct PeriodicNeighbor {
    pub index: usize,
    pub distance: f64,
}

impl PeriodicNeighborSearch {
    /// Construct a periodic neighbor search from the given **Minkowski-reduced** cell.
    pub fn new(reduced_cell: &Cell, symprec: f64) -> Self {
        let lattice = reduced_cell.lattice.clone();

        // Twice the padding for safety
        let padding =
            2.0 * symprec / (3.0 * (lattice.basis * lattice.basis.transpose()).trace()).sqrt();

        // `f64::EPSILON` floor keeps the bin size positive for degenerate symprec;
        // queries then fall back to exact-match behavior instead of dividing by zero.
        let bin_size = (2.0 * symprec).max(f64::EPSILON);
        let mut bins: FxHashMap<(i64, i64, i64), Vec<usize>> = FxHashMap::default();
        let mut cart_coords = vec![];
        let mut indices = vec![];
        for offset in iproduct!(-1..=1, -1..=1, -1..=1) {
            for (index, position) in reduced_cell.positions.iter().enumerate() {
                let mut new_position = *position;
                new_position -= position.map(|e| e.floor()); // [0, 1)
                new_position += Vector3::new(offset.0 as f64, offset.1 as f64, offset.2 as f64);
                if new_position[0] < -padding
                    || new_position[0] > 1.0 + padding
                    || new_position[1] < -padding
                    || new_position[1] > 1.0 + padding
                    || new_position[2] < -padding
                    || new_position[2] > 1.0 + padding
                {
                    continue;
                }

                let cart = lattice.cartesian_coords(&new_position);
                // Register the entry in every bin overlapping its symprec-ball
                let lo = Self::bin_key(&cart.add_scalar(-symprec), bin_size);
                let hi = Self::bin_key(&cart.add_scalar(symprec), bin_size);
                for key in iproduct!(lo.0..=hi.0, lo.1..=hi.1, lo.2..=hi.2) {
                    bins.entry(key).or_default().push(cart_coords.len());
                }
                cart_coords.push(cart);
                indices.push(index);
            }
        }

        Self {
            num_sites: reduced_cell.num_atoms(),
            lattice,
            bin_size,
            bins,
            cart_coords,
            indices,
            symprec,
        }
    }

    fn bin_key(cart_coords: &Vector3<f64>, bin_size: f64) -> (i64, i64, i64) {
        // `as i64` saturates on overflow, which only degrades bin locality,
        // never memory safety.
        (
            (cart_coords.x / bin_size).floor() as i64,
            (cart_coords.y / bin_size).floor() as i64,
            (cart_coords.z / bin_size).floor() as i64,
        )
    }

    /// Return the nearest neighbor within symprec if exists.
    pub fn nearest(&self, position: &Position) -> Option<PeriodicNeighbor> {
        let mut wrapped_position = *position;
        wrapped_position -= wrapped_position.map(|e| e.floor()); // [0, 1)
        let cart = self.lattice.cartesian_coords(&wrapped_position);

        let mut nearest: Option<PeriodicNeighbor> = None;
        let entries = self.bins.get(&Self::bin_key(&cart, self.bin_size))?;
        for &entry in entries {
            let distance = (self.cart_coords[entry] - cart).norm();
            if distance <= self.symprec && nearest.as_ref().is_none_or(|n| distance < n.distance) {
                nearest = Some(PeriodicNeighbor {
                    index: self.indices[entry],
                    distance,
                });
            }
        }
        nearest
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
    numbers
        .iter()
        .enumerate()
        .filter(|(_, number)| *number == pivot_atomic_specie)
        .map(|(i, _)| i)
        .collect::<Vec<_>>()
}

/// Return correspondence between the input and acted positions.
/// Search permutation such that new_positions\[i\] = reduced_cell.positions\[permutation\[i\]\].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation\[i\]-th site.
/// This function takes O(num_atoms) expected time thanks to the bin-grid nearest-neighbor queries.
#[doc(hidden)]
pub fn solve_correspondence(
    neighbor_search: &PeriodicNeighborSearch,
    reduced_cell: &Cell,
    new_positions: &[Position],
) -> Option<Permutation> {
    let num_atoms = neighbor_search.num_sites;
    let mut mapping = vec![None; num_atoms];

    for i in 0..num_atoms {
        let neighbor = neighbor_search.nearest(&new_positions[i])?;
        let j = neighbor.index;
        if reduced_cell.numbers[i] != reduced_cell.numbers[j] {
            return None;
        }
        if mapping[i].is_some() {
            return None;
        }

        mapping[i] = Some(j);
    }

    let mapping = mapping.into_iter().map(|v| v.unwrap()).collect::<Vec<_>>();
    assert_eq!(mapping.len(), num_atoms);
    Some(Permutation::new(mapping))
}

/// Return correspondence between the input and acted positions.
/// Assume that reduced_cell is Minkowski reduced and symprec is sufficiently small for Babai's algorithm.
/// Search permutation such that new_positions\[i\] = reduced_cell.positions\[permutation\[i\]\].
/// Then, a corresponding symmetry operation moves the i-th site into the permutation\[i\]-th site.
/// This function takes O(num_atoms^2) time.
#[doc(hidden)]
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
        PeriodicNeighborSearch, pivot_site_indices, solve_correspondence,
        solve_correspondence_naive, symmetrize_translation_from_permutation,
    };

    #[test]
    fn test_periodic_neighbor_search() {
        let neighbor_search = PeriodicNeighborSearch::new(
            &Cell::new(
                Lattice::new(Matrix3::identity()),
                vec![
                    Vector3::new(0.0, 0.0, 0.0),
                    Vector3::new(0.0, 0.5, 0.5),
                    Vector3::new(0.5, 0.0, 0.5),
                    Vector3::new(0.5, 0.5, 0.0),
                ],
                vec![0, 0, 0, 0],
            ),
            1e-4,
        );

        {
            let neighbor = neighbor_search
                .nearest(&Vector3::new(0.0, 0.0, 0.0))
                .unwrap();
            assert!(neighbor.index == 0);
        }

        {
            let neighbor = neighbor_search
                .nearest(&Vector3::new(1.0, 0.5, 0.5))
                .unwrap();
            assert!(neighbor.index == 1);
        }

        {
            let neighbor = neighbor_search
                .nearest(&Vector3::new(1.5, -0.0, -0.5))
                .unwrap();
            assert!(neighbor.index == 2);
        }
    }

    #[test]
    fn test_periodic_neighbor_search_axis_aligned_grid() {
        // Many sites sharing coordinate values on every axis (an axis-aligned
        // n x n x n grid). Cross-check every translation against the naive
        // O(num_atoms^2) solver.
        let n = 4;
        let mut positions = vec![];
        let mut numbers = vec![];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    positions.push(Vector3::new(
                        i as f64 / n as f64,
                        j as f64 / n as f64,
                        k as f64 / n as f64,
                    ));
                    numbers.push(0);
                }
            }
        }
        let reduced_cell = Cell::new(
            Lattice::new(Matrix3::identity() * n as f64),
            positions,
            numbers,
        );
        let symprec = 1e-4;
        let neighbor_search = PeriodicNeighborSearch::new(&reduced_cell, symprec);

        for dst in 0..reduced_cell.num_atoms() {
            let translation = reduced_cell.positions[dst] - reduced_cell.positions[0];
            let new_positions = reduced_cell
                .positions
                .iter()
                .map(|pos| pos + translation)
                .collect::<Vec<_>>();

            let actual = solve_correspondence(&neighbor_search, &reduced_cell, &new_positions);
            let expect = solve_correspondence_naive(&reduced_cell, &new_positions, symprec);
            assert_eq!(actual, expect);
            assert!(actual.is_some());
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
        let neighbor_search = PeriodicNeighborSearch::new(&reduced_cell, symprec);

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

            let actual_neighbor_search =
                solve_correspondence(&neighbor_search, &reduced_cell, &new_positions).unwrap();
            assert_eq!(actual_neighbor_search, expect);
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

            let actual_neighbor_search =
                solve_correspondence(&neighbor_search, &reduced_cell, &new_positions);
            assert_eq!(actual_neighbor_search, None);
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
