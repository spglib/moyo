use moyo::base::{Cell, Lattice};
use moyo::utils::{to_3_slice, to_3x3_slice, to_vector3};

use crate::ffi::{free_slice, leak_slice};

/// Crystal structure.
/// All pointer fields are owned by the containing `MoyoDataset` and are freed
/// by `moyo_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoCell {
    /// Row-wise basis vectors
    pub basis: [[f64; 3]; 3],
    /// Fractional coordinates of the `num_atoms` sites
    pub positions: *const [f64; 3],
    /// Atomic numbers of the `num_atoms` sites
    pub numbers: *const i32,
    /// Number of atoms in the cell
    pub num_atoms: i32,
}

impl From<&Cell> for MoyoCell {
    fn from(cell: &Cell) -> Self {
        let num_atoms = cell.num_atoms() as i32;
        // Transpose to row-major basis
        let basis = to_3x3_slice(&cell.lattice.basis.transpose());
        let positions = cell.positions.iter().map(to_3_slice).collect::<Vec<_>>();
        let numbers = cell.numbers.clone();

        MoyoCell {
            basis,
            positions: leak_slice(positions),
            numbers: leak_slice(numbers),
            num_atoms,
        }
    }
}

impl From<&MoyoCell> for Cell {
    fn from(cell: &MoyoCell) -> Self {
        let lattice = Lattice::from_basis(cell.basis);
        let positions = unsafe {
            std::slice::from_raw_parts(cell.positions, cell.num_atoms as usize)
                .iter()
                .map(to_vector3)
                .collect()
        };
        let numbers =
            unsafe { std::slice::from_raw_parts(cell.numbers, cell.num_atoms as usize).to_vec() };
        Cell::new(lattice, positions, numbers)
    }
}

pub(crate) unsafe fn moyo_cell_free_members(cell: &MoyoCell) {
    unsafe {
        free_slice(cell.positions, cell.num_atoms as usize);
        free_slice(cell.numbers, cell.num_atoms as usize);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use moyo::base::{Cell, Lattice};
    use nalgebra::{matrix, vector};

    #[test]
    fn test_roundtrip_moyo_cell() {
        let original = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                1.0, 1.0, 0.0;
                1.0, 0.0, 1.0
            ]),
            vec![vector![0.0, 0.0, 0.0]],
            vec![1],
        );
        let moyoc = MoyoCell::from(&original);
        let reconstructed = Cell::from(&moyoc);
        assert_eq!(original.num_atoms(), reconstructed.num_atoms());
        assert_relative_eq!(original.lattice.basis, reconstructed.lattice.basis);
        unsafe { moyo_cell_free_members(&moyoc) }; // to avoid leaks
    }
}
