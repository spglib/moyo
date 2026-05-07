use std::f64::consts::PI;

use log::debug;
use nalgebra::Matrix3;
use serde::{Deserialize, Serialize};

use super::error::MoyoError;
use super::lattice::Lattice;
use super::tolerance::{AngleTolerance, is_angle_within_tolerance};

/// A `Lattice` whose third basis vector is the aperiodic stacking direction
/// of a layer system, with `c` perpendicular to `a, b` (paper Fu et al. 2024
/// eq. 5).
///
/// Construction via [`LayerLattice::new`] runs the perpendicularity check up
/// front, so any function taking `&LayerLattice` can rely on `c . a = 0` and
/// `c . b = 0` within tolerance, plus all three basis vectors being non-zero.
/// The newtype prevents bulk-pipeline `Lattice` values from being passed where
/// the layer pipeline expects them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerLattice {
    inner: Lattice,
}

impl LayerLattice {
    /// Validate `lattice` against the layer-group periodicity contract
    /// (`c` perpendicular to `a, b`, all basis vectors non-zero) and wrap it.
    pub fn new(
        lattice: Lattice,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        let a = lattice.basis.column(0);
        let b = lattice.basis.column(1);
        let c = lattice.basis.column(2);

        let na = a.norm();
        let nb = b.norm();
        let nc = c.norm();
        if na == 0.0 || nb == 0.0 || nc == 0.0 {
            return Err(MoyoError::DegenerateLattice);
        }

        let dev_ca = c.angle(&a) - PI / 2.0;
        let dev_cb = c.angle(&b) - PI / 2.0;

        let ok_ca = is_angle_within_tolerance(dev_ca, nc, na, symprec, angle_tolerance);
        let ok_cb = is_angle_within_tolerance(dev_cb, nc, nb, symprec, angle_tolerance);

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
        Ok(Self { inner: lattice })
    }

    /// Wrap `lattice` without running the perpendicularity check.
    /// Use only at internal construction sites that produce a lattice
    /// guaranteed to satisfy the layer contract by construction (e.g. the
    /// primitive layer lattice built from in-plane translations of an
    /// already-validated input).
    pub(crate) fn new_unchecked(lattice: Lattice) -> Self {
        Self { inner: lattice }
    }

    /// Lattice basis matrix `[a | b | c]` (column-vector convention) with `c`
    /// the aperiodic stacking direction.
    ///
    /// Deliberately exposes the bare matrix instead of a `&Lattice`: a
    /// `&Lattice` is exactly the input shape bulk-pipeline helpers consume
    /// (e.g. `search_bravais_group(lattice, ..)`), so handing one out from
    /// `LayerLattice` would re-introduce the conflation the newtype prevents.
    /// In-crate sites that need a `Lattice` reconstruct it explicitly via
    /// `Lattice { basis: *ll.basis() }`.
    pub fn basis(&self) -> &Matrix3<f64> {
        &self.inner.basis
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::*;

    const TEST_SYMPREC: f64 = 1e-4;

    #[test]
    fn test_layer_lattice_new_orthogonal() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let ll = LayerLattice::new(lattice, TEST_SYMPREC, AngleTolerance::Default).unwrap();
        assert_eq!(ll.basis().column(2)[2], 5.0);
    }

    #[test]
    fn test_layer_lattice_new_oblique_in_plane_ok() {
        // Hexagonal-style oblique in-plane (a,b) with c perpendicular.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            -0.5, (3.0_f64).sqrt() / 2.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        assert!(LayerLattice::new(lattice, TEST_SYMPREC, AngleTolerance::Default).is_ok());
    }

    #[test]
    fn test_layer_lattice_new_rejects_tilted_c() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.5, 0.0, 5.0;
        ]);
        assert!(matches!(
            LayerLattice::new(lattice, TEST_SYMPREC, AngleTolerance::Default),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }

    #[test]
    fn test_layer_lattice_new_explicit_radian() {
        // Small tilt (~1.7 deg) accepted with 5-deg tolerance, rejected with 1-deg.
        let lattice_factory = || {
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.15, 0.0, 5.0;
            ])
        };
        assert!(
            LayerLattice::new(
                lattice_factory(),
                TEST_SYMPREC,
                AngleTolerance::Radian(5.0_f64.to_radians()),
            )
            .is_ok()
        );
        assert!(matches!(
            LayerLattice::new(
                lattice_factory(),
                TEST_SYMPREC,
                AngleTolerance::Radian(1.0_f64.to_radians()),
            ),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }

    #[test]
    fn test_layer_lattice_new_rejects_zero_norm() {
        // c = 0 is a degenerate basis -- caught before the angle check (which
        // would otherwise NaN-out and silently pass under AngleTolerance::Default).
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 0.0;
        ]);
        assert!(matches!(
            LayerLattice::new(lattice, TEST_SYMPREC, AngleTolerance::Default),
            Err(MoyoError::DegenerateLattice)
        ));
    }
}
