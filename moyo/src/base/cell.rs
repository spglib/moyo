use std::collections::BTreeMap;
use std::f64::consts::PI;

use log::debug;
use nalgebra::{Matrix3, Vector3};

use crate::utils::to_3_slice;
use serde::{Deserialize, Serialize};
use union_find::{QuickFindUf, UnionByRank, UnionFind};

use super::error::MoyoError;
use super::lattice::Lattice;
use super::permutation::Permutation;
use super::tolerance::{AngleTolerance, is_angle_within_tolerance};

/// Fractional coordinates
pub type Position = Vector3<f64>;
/// Atomic number
pub type AtomicSpecie = i32;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Representing a crystal structure
pub struct Cell {
    /// Lattice of the cell.
    pub lattice: Lattice,
    /// `positions[i]` is a fractional coordinates of the i-th site.
    pub positions: Vec<Position>,
    /// `numbers[i]` is an atomic number of the i-th site.
    pub numbers: Vec<AtomicSpecie>,
}

impl Cell {
    pub fn new(lattice: Lattice, positions: Vec<Position>, numbers: Vec<AtomicSpecie>) -> Self {
        if positions.len() != numbers.len() {
            panic!("positions and numbers should be the same length");
        }
        Self {
            lattice,
            positions,
            numbers,
        }
    }

    /// Return the number of atoms in the cell.
    pub fn num_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Returns fractional positions as a vector of `[f64; 3]` arrays.
    ///
    /// Each `[f64; 3]` contains the fractional coordinates `[x, y, z]`
    /// of a site in the cell.
    ///
    /// This is a convenience method for users who do not depend on `nalgebra`.
    /// It is equivalent to converting each `Vector3<f64>` in `self.positions`.
    pub fn positions_as_arrays(&self) -> Vec<[f64; 3]> {
        self.positions.iter().map(|p| to_3_slice(p)).collect()
    }

    /// Rotate the cell by the given rotation matrix.
    pub fn rotate(&self, rotation_matrix: &Matrix3<f64>) -> Self {
        Self::new(
            self.lattice.rotate(rotation_matrix),
            self.positions.clone(),
            self.numbers.clone(),
        )
    }
}

/// Validate the layer-group convention that the third basis vector `c` is
/// perpendicular to the in-plane axes `a` and `b` (paper Fu et al. 2024 eq. 5).
///
/// Defers to `is_angle_within_tolerance` so the `(symprec, angle_tolerance)`
/// pair is interpreted exactly as in the rest of the symmetry pipeline.
pub fn validate_aperiodic_axis(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<(), MoyoError> {
    validate_lattice_aperiodic_axis(&cell.lattice, symprec, angle_tolerance)
}

/// Lattice-level half of `validate_aperiodic_axis`. The cell-level wrapper
/// delegates here; `LayerLattice::new` calls this directly without needing a
/// `Cell` to wrap.
pub fn validate_lattice_aperiodic_axis(
    lattice: &Lattice,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<(), MoyoError> {
    let a = lattice.basis.column(0);
    let b = lattice.basis.column(1);
    let c = lattice.basis.column(2);

    let dev_ca = c.angle(&a) - PI / 2.0;
    let dev_cb = c.angle(&b) - PI / 2.0;

    let ok_ca = is_angle_within_tolerance(dev_ca, c.norm(), a.norm(), symprec, angle_tolerance);
    let ok_cb = is_angle_within_tolerance(dev_cb, c.norm(), b.norm(), symprec, angle_tolerance);

    if !ok_ca || !ok_cb {
        debug!(
            "Aperiodic axis is not orthogonal: dev(c,a)={:.6} rad, dev(c,b)={:.6} rad",
            dev_ca, dev_cb
        );
        return Err(MoyoError::AperiodicAxisNotOrthogonal {
            dev_ca: dev_ca.abs(),
            dev_cb: dev_cb.abs(),
        });
    }
    Ok(())
}

/// If and only if the `i`th and `j`th atoms are equivalent, `orbits[i] == orbits[j]`.
/// For each orbit, only one of them satisfies `orbits[i] == i`.
pub fn orbits_from_permutations(num_atoms: usize, permutations: &[Permutation]) -> Vec<usize> {
    let mut uf = QuickFindUf::<UnionByRank>::new(num_atoms);
    for permutation in permutations.iter() {
        for i in 0..num_atoms {
            uf.union(i, permutation.apply(i));
        }
    }
    let mut identifier_mapping = BTreeMap::new();
    for i in 0..num_atoms {
        identifier_mapping.entry(uf.find(i)).or_insert(i);
    }

    (0..num_atoms)
        .map(|i| *identifier_mapping.get(&uf.find(i)).unwrap())
        .collect()
}

#[cfg(test)]
mod tests {
    use std::panic;

    use nalgebra::{Matrix3, matrix, vector};

    use super::{Cell, orbits_from_permutations, validate_aperiodic_axis};
    use crate::base::error::MoyoError;
    use crate::base::lattice::Lattice;
    use crate::base::permutation::Permutation;
    use crate::base::tolerance::AngleTolerance;

    #[test]
    fn test_orbits_from_permutations() {
        {
            let num_atoms = 3;
            let permutations = vec![Permutation::new(vec![2, 1, 0])];
            assert_eq!(
                orbits_from_permutations(num_atoms, &permutations),
                vec![0, 1, 0]
            );
        }
        {
            let num_atoms = 3;
            let permutations = vec![Permutation::new(vec![1, 0, 2])];
            assert_eq!(
                orbits_from_permutations(num_atoms, &permutations),
                vec![0, 0, 2]
            );
        }
    }

    #[test]
    fn test_mismatched_length() {
        let lattice = Lattice::new(Matrix3::<f64>::identity());
        let positions = vec![vector![0.0, 0.0, 0.0], vector![0.5, 0.5, 0.5]];
        let numbers = vec![1];

        let result = panic::catch_unwind(|| Cell::new(lattice, positions, numbers));
        assert!(result.is_err());
    }

    const TEST_SYMPREC: f64 = 1e-4;

    #[test]
    fn test_validate_aperiodic_axis_orthogonal() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        assert!(validate_aperiodic_axis(&cell, TEST_SYMPREC, AngleTolerance::Default).is_ok());
    }

    #[test]
    fn test_validate_aperiodic_axis_oblique_in_plane_ok() {
        // Hexagonal-style oblique in-plane (a,b) with c perpendicular:
        // a, b in xy-plane; c along z. Validator must pass.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            -0.5, (3.0_f64).sqrt() / 2.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        assert!(validate_aperiodic_axis(&cell, TEST_SYMPREC, AngleTolerance::Default).is_ok());
    }

    #[test]
    fn test_validate_aperiodic_axis_tilted_c_rejected() {
        // c has a non-trivial in-plane component; should be rejected with default tol.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.5, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        assert!(matches!(
            validate_aperiodic_axis(&cell, TEST_SYMPREC, AngleTolerance::Default),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }

    #[test]
    fn test_validate_aperiodic_axis_explicit_radian() {
        // Small tilt (~1.7 degrees) accepted with a 5-degree tolerance, rejected with a 1-degree one.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.15, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        assert!(
            validate_aperiodic_axis(
                &cell,
                TEST_SYMPREC,
                AngleTolerance::Radian(5.0_f64.to_radians())
            )
            .is_ok()
        );
        assert!(matches!(
            validate_aperiodic_axis(
                &cell,
                TEST_SYMPREC,
                AngleTolerance::Radian(1.0_f64.to_radians())
            ),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }
}
