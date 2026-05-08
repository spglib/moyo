use std::f64::consts::PI;

use log::debug;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use super::super::error::MoyoError;
use super::super::lattice::Lattice;
use super::super::tolerance::{AngleTolerance, is_angle_within_tolerance};

/// A `Lattice` whose third basis vector is the aperiodic stacking direction
/// of a layer system, with `c` along the z-axis and `a, b` in the xy-plane
/// (paper Fu et al. 2024 eq. 5).
///
/// Construction via [`LayerLattice::new`] enforces this orientation up front,
/// so any function taking `&LayerLattice` can rely on `a_z = b_z = 0` (and
/// hence `c . a = c . b = 0`) within tolerance, plus all three basis vectors
/// being non-zero. The newtype prevents bulk-pipeline `Lattice` values from
/// being passed where the layer pipeline expects them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerLattice {
    inner: Lattice,
}

impl LayerLattice {
    /// Validate `lattice` against the layer-group periodicity contract
    /// (`c` along z, `a, b` in the xy-plane, all basis vectors non-zero) and
    /// wrap it.
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

        // Reject globally-rotated layer cells: `c . a = c . b = 0` alone admits
        // any 3D rotation of the slab, but the in-plane block extracted by
        // `Lattice2D::from_inplane_of` only retains the upper-left 2x2 of the
        // basis. Require `a, b` in the xy-plane so that extraction is exact.
        let z_axis = Vector3::new(0.0_f64, 0.0, 1.0);
        let dev_az = a.angle(&z_axis) - PI / 2.0;
        let dev_bz = b.angle(&z_axis) - PI / 2.0;
        let ok_az = is_angle_within_tolerance(dev_az, na, 1.0, symprec, angle_tolerance);
        let ok_bz = is_angle_within_tolerance(dev_bz, nb, 1.0, symprec, angle_tolerance);
        if !ok_az || !ok_bz {
            debug!(
                "In-plane axes are not in the xy-plane: dev(a,xy)={:.6} rad, dev(b,xy)={:.6} rad",
                dev_az, dev_bz
            );
            return Err(MoyoError::InPlaneAxesNotInXY {
                dev_az: dev_az.abs(),
                dev_bz: dev_bz.abs(),
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
    /// Public callers go through this getter (rather than a `&Lattice`
    /// accessor) so a layer-pipeline value cannot silently flow into
    /// bulk-only entry points (e.g. `search_bravais_group(lattice, ..)`).
    /// In-crate sites that need a bulk `Lattice` use `Self::as_lattice`.
    pub fn basis(&self) -> &Matrix3<f64> {
        &self.inner.basis
    }

    /// Bulk `Lattice` view of the same basis, for in-crate helpers that
    /// consume a `Lattice` (the LG-restricted Bravais filter, the
    /// metric-tensor lattice symmetrizer). Returns a new `Lattice` rather
    /// than a borrow so a `LayerLattice` and a `Lattice` never share a
    /// backing object the caller could later mistake for the layer one.
    pub(crate) fn as_lattice(&self) -> Lattice {
        Lattice {
            basis: self.inner.basis,
        }
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
    fn test_layer_lattice_new_rejects_globally_rotated() {
        // A globally-rotated valid layer cell: c is still perpendicular to a,b
        // (so the orthogonality check passes) but a, b have non-zero z-components.
        // The 2x2 in-plane extraction would silently drop those z-components, so
        // construction must reject this case.
        let s = (0.5_f64).sqrt();
        let lattice = Lattice::new(matrix![
            s, 0.0, -s;
            0.0, 1.0, 0.0;
            s, 0.0, s;
        ]);
        assert!(matches!(
            LayerLattice::new(lattice, TEST_SYMPREC, AngleTolerance::Default),
            Err(MoyoError::InPlaneAxesNotInXY { .. })
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
