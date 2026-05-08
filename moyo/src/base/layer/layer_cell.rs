use serde::{Deserialize, Serialize};

use super::super::cell::{AtomicSpecie, Cell, Position};
use super::super::error::MoyoError;
use super::super::tolerance::AngleTolerance;
use super::layer_lattice::LayerLattice;

/// A cell whose third basis vector is the aperiodic stacking direction
/// of a layer system (paper Fu et al. 2024 §2).
///
/// Internally a tuple of `(LayerLattice, positions, numbers)` rather than a
/// wrapped `Cell`, so neither `LayerCell` nor `LayerLattice` hands out a borrow
/// of the bulk type at any visibility. In-crate sites that need a `Cell`
/// reconstruct one explicitly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerCell {
    lattice: LayerLattice,
    positions: Vec<Position>,
    numbers: Vec<AtomicSpecie>,
}

impl LayerCell {
    /// Validate `cell.lattice` against the layer-group periodicity contract
    /// (`c` perpendicular to `a, b`) and decompose `cell` into a `LayerCell`.
    pub fn new(
        cell: Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        let Cell {
            lattice,
            positions,
            numbers,
        } = cell;
        let layer_lattice = LayerLattice::new(lattice, symprec, angle_tolerance)?;
        Ok(Self {
            lattice: layer_lattice,
            positions,
            numbers,
        })
    }

    /// Wrap pre-validated parts without running the perpendicularity check.
    /// Use only at internal construction sites that produce parts guaranteed
    /// to satisfy the layer contract by construction.
    pub(crate) fn new_unchecked(
        lattice: LayerLattice,
        positions: Vec<Position>,
        numbers: Vec<AtomicSpecie>,
    ) -> Self {
        Self {
            lattice,
            positions,
            numbers,
        }
    }

    pub fn lattice(&self) -> &LayerLattice {
        &self.lattice
    }

    pub fn positions(&self) -> &[Position] {
        &self.positions
    }

    pub fn numbers(&self) -> &[AtomicSpecie] {
        &self.numbers
    }

    pub fn num_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Bulk `Cell` view of this layer cell, for in-crate helpers that
    /// consume a `Cell` (the bulk symmetry search, primitive-cell finder,
    /// transformation routines). Allocates a fresh `Cell` rather than
    /// handing out a borrow so the layer-to-bulk crossing remains a value
    /// boundary: bulk-only state can never be aliased back into the
    /// `LayerCell`.
    pub(crate) fn as_cell(&self) -> Cell {
        Cell::new(
            self.lattice.as_lattice(),
            self.positions.clone(),
            self.numbers.clone(),
        )
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
        assert_eq!(layer.lattice().basis().column(2)[2], 5.0);
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
