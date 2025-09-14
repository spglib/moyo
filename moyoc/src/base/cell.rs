use moyo::base::{Cell, Lattice};
use moyo::utils::{to_3_slice, to_3x3_slice, to_vector3};

#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoCell {
    /// Row-wise basis vectors
    pub basis: [[f64; 3]; 3],
    pub positions: *const [f64; 3],
    pub numbers: *const i32,
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
            positions: positions.leak().as_ptr(),
            numbers: numbers.leak().as_ptr(),
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
        Cell::new(lattice.into(), positions, numbers)
    }
}

#[no_mangle]
pub extern "C" fn free_moyo_cell(cell: MoyoCell) {
    unsafe {
        let _ = Vec::from_raw_parts(
            cell.positions as *mut f64,
            (cell.num_atoms as usize) * 3,
            (cell.num_atoms as usize) * 3,
        );
        let _ = Vec::from_raw_parts(
            cell.numbers as *mut i32,
            cell.num_atoms as usize,
            cell.num_atoms as usize,
        );
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
        free_moyo_cell(moyoc); // to avoid leaks
    }
}
