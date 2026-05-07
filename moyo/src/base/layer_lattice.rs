use nalgebra::Matrix3;
use serde::{Deserialize, Serialize};

use super::cell::validate_lattice_aperiodic_axis;
use super::error::MoyoError;
use super::lattice::Lattice;
use super::tolerance::AngleTolerance;

/// A `Lattice` whose third basis vector is the aperiodic stacking direction
/// of a layer system, with `c` perpendicular to `a, b` (paper Fu et al. 2024
/// eq. 5).
///
/// Construction via [`LayerLattice::new`] runs the perpendicularity check up
/// front, so any function taking `&LayerLattice` can rely on `c . a = 0` and
/// `c . b = 0` within tolerance. The newtype prevents bulk-pipeline `Lattice`
/// values from being passed where the layer pipeline expects them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerLattice {
    inner: Lattice,
}

impl LayerLattice {
    /// Validate `lattice` against the layer-group periodicity contract
    /// (`c` perpendicular to `a, b`) and wrap it. Returns
    /// `Err(MoyoError::AperiodicAxisNotOrthogonal { .. })` when `c` is not
    /// perpendicular to the in-plane axes within tolerance.
    pub fn new(
        lattice: Lattice,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        validate_lattice_aperiodic_axis(&lattice, symprec, angle_tolerance)?;
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

    #[test]
    fn test_layer_lattice_new_orthogonal() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let ll = LayerLattice::new(lattice, 1e-4, AngleTolerance::Default).unwrap();
        assert_eq!(ll.basis().column(2)[2], 5.0);
    }

    #[test]
    fn test_layer_lattice_new_rejects_tilted_c() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.5, 0.0, 5.0;
        ]);
        assert!(matches!(
            LayerLattice::new(lattice, 1e-4, AngleTolerance::Default),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }
}
