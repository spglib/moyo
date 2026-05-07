use serde::{Deserialize, Serialize};

use super::cell::{AtomicSpecie, Cell, Position, validate_aperiodic_axis};
use super::error::MoyoError;
use super::lattice::Lattice;
use super::tolerance::AngleTolerance;

/// A `Cell` whose third basis vector is the aperiodic stacking direction
/// of a layer system (paper Fu et al. 2024 §2).
///
/// Construction via [`LayerCell::new`] runs [`validate_aperiodic_axis`] so that
/// any function taking `&LayerCell` can rely on `c . a = 0` and `c . b = 0`
/// within tolerance. The newtype prevents bulk-pipeline `Cell` values from
/// being passed where the layer pipeline expects them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerCell {
    inner: Cell,
}

impl LayerCell {
    /// Validate `cell` against the layer-group periodicity contract (§3.2)
    /// and wrap it. Returns `Err(MoyoError::AperiodicAxisNotOrthogonal { .. })`
    /// when `c` is not perpendicular to the in-plane axes within tolerance.
    pub fn new(
        cell: Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        validate_aperiodic_axis(&cell, symprec, angle_tolerance)?;
        Ok(Self { inner: cell })
    }

    /// Wrap `cell` without running the perpendicularity check.
    /// Use only at internal construction sites that produce a cell guaranteed
    /// to satisfy the layer contract by construction (e.g. the primitive layer
    /// cell built from in-plane translations of an already-validated input).
    pub(crate) fn new_unchecked(cell: Cell) -> Self {
        Self { inner: cell }
    }

    /// Internal access to the wrapped `Cell`.
    ///
    /// Deliberately `pub(crate)`: a layer cell must not be downgradeable to a
    /// bulk `Cell` from outside the crate, otherwise the bulk pipeline could
    /// silently accept a layer-shaped input. Use the field-level getters
    /// (`lattice`, `positions`, `numbers`) for read-only inspection.
    pub(crate) fn cell(&self) -> &Cell {
        &self.inner
    }

    pub fn lattice(&self) -> &Lattice {
        &self.inner.lattice
    }

    pub fn positions(&self) -> &[Position] {
        &self.inner.positions
    }

    pub fn numbers(&self) -> &[AtomicSpecie] {
        &self.inner.numbers
    }

    pub fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};

    use super::*;
    use crate::base::Lattice;

    #[test]
    fn test_layer_cell_new_orthogonal() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        let layer = LayerCell::new(cell, 1e-4, AngleTolerance::Default).unwrap();
        assert_eq!(layer.num_atoms(), 1);
    }

    #[test]
    fn test_layer_cell_new_rejects_tilted_c() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.5, 0.0, 5.0;
        ]);
        let cell = Cell::new(lattice, vec![vector![0.0, 0.0, 0.0]], vec![1]);
        assert!(matches!(
            LayerCell::new(cell, 1e-4, AngleTolerance::Default),
            Err(MoyoError::AperiodicAxisNotOrthogonal { .. })
        ));
    }
}
