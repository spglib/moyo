use itertools::iproduct;
use log::debug;
use nalgebra::Vector3;
use std::collections::HashMap;

use crate::base::{
    Lattice, MoyoError, Permutation, Position, UnimodularTransformation, orbits_from_permutations,
};
use crate::data::{HallNumber, WyckoffPosition, WyckoffPositionSpace, iter_wyckoff_positions};
use crate::math::SNF;

/// * `prim_num_atoms` - Number of atoms in the primitive cell
/// * `prim_permutations` - Permutations of the primitive cell
/// * `site_mapping` - Mapping from the site in cell to that in its primitive cell
pub fn orbits_in_cell(
    prim_num_atoms: usize,
    prim_permutations: &[Permutation],
    site_mapping: &[usize],
) -> Vec<usize> {
    // prim_orbits: [prim_num_atoms] -> [prim_num_atoms]
    let prim_orbits = orbits_from_permutations(prim_num_atoms, prim_permutations);

    let num_atoms = site_mapping.len();
    let mut map = HashMap::new();
    let mut orbits = vec![]; // [num_atoms] -> [num_atoms]
    for i in 0..num_atoms {
        // site_mapping: [num_atoms] -> [prim_num_atoms]
        let key = prim_orbits[site_mapping[i]]; // in [prim_num_atoms]
        map.entry(key).or_insert(i);
        orbits.push(*map.get(&key).unwrap());
    }
    orbits
}

/// Group sites of a conventional cell by crystallographic orbit. Shared
/// between the bulk and layer pipelines.
pub(super) struct OrbitGrouping {
    /// `[std_num_atoms] -> [num_orbits]`: orbit index for each site.
    pub mapping: Vec<usize>,
    /// `[num_orbits] -> [std_num_atoms]`: representative site for each orbit
    /// (the lowest-index site in the orbit). Used only for diagnostics.
    pub remapping: Vec<usize>,
    /// Number of sites per orbit, indexed by orbit. `multiplicities.len()` is
    /// the orbit count.
    pub multiplicities: Vec<usize>,
}

pub(super) fn group_sites_by_orbit(
    prim_num_atoms: usize,
    prim_permutations: &[Permutation],
    site_mapping: &[usize],
    std_num_atoms: usize,
) -> OrbitGrouping {
    let orbits = orbits_in_cell(prim_num_atoms, prim_permutations, site_mapping);

    let mut num_orbits = 0;
    let mut mapping = vec![0; std_num_atoms];
    let mut remapping = vec![];
    for i in 0..std_num_atoms {
        if orbits[i] == i {
            mapping[i] = num_orbits;
            remapping.push(i);
            num_orbits += 1;
        } else {
            mapping[i] = mapping[orbits[i]];
        }
    }

    let mut multiplicities = vec![0; num_orbits];
    for i in 0..std_num_atoms {
        multiplicities[mapping[i]] += 1;
    }

    OrbitGrouping {
        mapping,
        remapping,
        multiplicities,
    }
}

/// Run the per-orbit Wyckoff-assignment loop over a precomputed
/// `OrbitGrouping`. The closure produces the database-specific Wyckoff
/// representative for an `(orbit_position, multiplicity)` pair (returning
/// `None` if no orbit in the database matches that multiplicity at that
/// position). Shared between bulk and layer pipelines: the only thing that
/// differs is the database lookup, which the caller passes in as the
/// closure.
pub(super) fn assign_wyckoffs_by_orbit<W, F>(
    group: &OrbitGrouping,
    positions: &[Position],
    mut find_for_orbit: F,
) -> Result<Vec<W>, MoyoError>
where
    W: Clone,
    F: FnMut(&Position, usize) -> Option<W>,
{
    let mut representative_wyckoffs: Vec<Option<W>> = vec![None; group.multiplicities.len()];
    for (i, position) in positions.iter().enumerate() {
        let orbit = group.mapping[i];
        if representative_wyckoffs[orbit].is_some() {
            continue;
        }
        if let Some(w) = find_for_orbit(position, group.multiplicities[orbit]) {
            representative_wyckoffs[orbit] = Some(w);
        }
    }

    for (orbit, wyckoff) in representative_wyckoffs.iter().enumerate() {
        if wyckoff.is_none() {
            debug!(
                "Failed to assign Wyckoff position with multiplicity {} at representative site {}",
                group.multiplicities[orbit], group.remapping[orbit]
            );
        }
    }
    let representative_wyckoffs = representative_wyckoffs
        .into_iter()
        .map(|w| w.ok_or(MoyoError::WyckoffPositionAssignmentError))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(group
        .mapping
        .iter()
        .map(|&orbit| representative_wyckoffs[orbit].clone())
        .collect())
}

/// For each conventional-basis normalizer operation, apply it to `positions`
/// and recompute the Wyckoff positions of the structure.
///
/// Returns one inner `Vec` per operation (same order as `conventional_ops`),
/// each holding one [`WyckoffPosition`] per atom in `positions` order. The
/// result is RAW: duplicate assignments produced by different operations are
/// NOT removed (the caller deduplicates / canonicalizes).
///
/// Contract: `conventional_ops`, `lattice`, and `positions` must all be in the
/// conventional standardized setting (the Wyckoff database coordinates are);
/// `site_orbits[i]` is the orbit index of atom `i` and `orbit_multiplicities`
/// is indexed by orbit. The orbit structure is invariant under the normalizer,
/// so it is computed once by the caller and reused for every operation. The
/// lattice is normalizer-invariant (`P^T M P = M`), so it is used unchanged for
/// the Cartesian-distance tolerance inside [`match_wyckoff_coordinates`].
pub(crate) fn wyckoff_positions_under_normalizer(
    conventional_ops: &[UnimodularTransformation],
    lattice: &Lattice,
    positions: &[Position],
    hall_number: HallNumber,
    site_orbits: &[usize],
    orbit_multiplicities: &[usize],
    symprec: f64,
) -> Result<Vec<Vec<WyckoffPosition>>, MoyoError> {
    let num_orbits = orbit_multiplicities.len();

    let mut result = Vec::with_capacity(conventional_ops.len());
    for op in conventional_ops {
        let linear = op.linear.map(|e| e as f64);

        // Apply the active map `x -> P x + p` (reduced mod 1) to each atom and
        // assign a Wyckoff position per orbit. `match_wyckoff_coordinates` only
        // matches the database representative coordinate (up to lattice
        // translations and the free-variable span), so -- as in
        // `assign_wyckoffs_by_orbit` -- each orbit is tried against successive
        // atoms until one lands on that representative.
        let mut representative_wyckoffs: Vec<Option<WyckoffPosition>> = vec![None; num_orbits];
        for (i, &orbit) in site_orbits.iter().enumerate() {
            if representative_wyckoffs[orbit].is_some() {
                continue;
            }
            let transformed = (linear * positions[i] + op.origin_shift).map(|e| e.rem_euclid(1.0));
            if let Some(wyckoff) = iter_wyckoff_positions(hall_number, orbit_multiplicities[orbit])
                .find(|w| match_wyckoff_coordinates(&transformed, w.coordinates, lattice, symprec))
                .cloned()
            {
                representative_wyckoffs[orbit] = Some(wyckoff);
            }
        }

        let representative_wyckoffs = representative_wyckoffs
            .into_iter()
            .map(|w| w.ok_or(MoyoError::WyckoffPositionAssignmentError))
            .collect::<Result<Vec<_>, _>>()?;
        let assignment = site_orbits
            .iter()
            .map(|&orbit| representative_wyckoffs[orbit].clone())
            .collect();
        result.push(assignment);
    }

    Ok(result)
}

/// Test whether `position` matches a Wyckoff orbit described by `coordinates`
/// (the parsed-coordinates string from a Wyckoff entry, e.g. `"x,1/2,z"`)
/// modulo a bounded integer offset. Shared between the bulk and layer
/// pipelines: both Wyckoff databases store their representative coordinate as
/// the same `&str` shape, so the SNF-based offset search is identical.
///
/// The search probes offsets in `[-1, 1]^3` first, then the remaining shell
/// in `[-2, 2]^3`. The math (suppressing tolerances): we want `y, offset` so
/// that `lattice * (space.linear * y + space.origin - position - offset)` is
/// near zero. Decomposing `space.linear` as `D = L * space.linear * R` (SNF)
/// reduces the integer search to picking `D^{-1} L (offset + position -
/// space.origin)` and reading off whether the residual is small in Cartesian
/// distance.
pub(super) fn match_wyckoff_coordinates(
    position: &Position,
    coordinates: &str,
    lattice: &Lattice,
    symprec: f64,
) -> bool {
    let space = WyckoffPositionSpace::new(coordinates);
    let snf = SNF::new(&space.linear);

    let iter_multi_1 = iproduct!(-1..=1, -1..=1, -1..=1);
    let iter_multi_2 = iproduct!(-2_i32..=2_i32, -2_i32..=2_i32, -2_i32..=2_i32)
        .filter(|&(n1, n2, n3)| n1.abs() == 2 || n2.abs() == 2 || n3.abs() == 2);

    for offset in iter_multi_1.chain(iter_multi_2) {
        let offset = Vector3::new(offset.0 as f64, offset.1 as f64, offset.2 as f64);
        let b = snf.l.map(|e| e as f64) * (offset + position - space.origin);
        let mut rinvy = Vector3::zeros();
        for i in 0..3 {
            if snf.d[(i, i)] != 0 {
                rinvy[i] = b[i] / snf.d[(i, i)] as f64;
            }
        }

        let y = snf.r.map(|e| e as f64) * rinvy;
        let diff = space.linear.map(|e| e as f64) * y + space.origin - position - offset;
        if lattice.cartesian_coords(&diff).norm() < symprec {
            return true;
        }
    }
    false
}
