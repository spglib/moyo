use serde::{Deserialize, Serialize};

use super::cell::{AtomicSpecie, Cell, Position};
use super::error::MoyoError;
use super::layer_lattice::LayerLattice;
use super::tolerance::AngleTolerance;

/// A cell whose third basis vector is the aperiodic stacking direction
/// of a layer system (paper Fu et al. 2024 §2; plan §2.3).
///
/// Internally a tuple of `(LayerLattice, positions, numbers)` rather than a
/// wrapped `Cell`, so neither `LayerCell` nor `LayerLattice` hands out a borrow
/// of the bulk type at any visibility. In-crate sites that need a `Cell`
/// reconstruct one explicitly (see `to_cell`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerCell {
    lattice: LayerLattice,
    positions: Vec<Position>,
    numbers: Vec<AtomicSpecie>,
}

impl LayerCell {
    /// Validate `cell.lattice` against the layer-group periodicity contract
    /// (§3.2) and decompose `cell` into a `LayerCell`.
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
