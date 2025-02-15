use nalgebra::vector;

use moyo::base::{Cell, Lattice};

#[repr(C)]
pub struct MoyoLattice {
    /// Row-wise basis vectors
    basis: [[f64; 3]; 3],
}

impl From<&MoyoLattice> for Lattice {
    fn from(lattice: &MoyoLattice) -> Self {
        Lattice::from_basis(lattice.basis)
    }
}

impl From<&Lattice> for MoyoLattice {
    fn from(lattice: &Lattice) -> Self {
        MoyoLattice {
            // Since nalgebra stores matrices in column-major order, we need to transpose them
            basis: *lattice.basis.transpose().as_ref(),
        }
    }
}

#[repr(C)]
pub struct MoyoCell {
    lattice: MoyoLattice,
    positions: *const [f64; 3],
    numbers: *const i32,
    num_atoms: i32,
}

impl From<&Cell> for MoyoCell {
    fn from(cell: &Cell) -> Self {
        let num_atoms = cell.num_atoms() as i32;
        let lattice = &cell.lattice;
        MoyoCell {
            lattice: lattice.into(),
            positions: {
                let positions = cell
                    .positions
                    .iter()
                    .map(|pos| [pos[0], pos[1], pos[2]])
                    .collect::<Vec<_>>();
                positions.as_ptr()
            },
            numbers: cell.numbers.as_ptr(),
            num_atoms,
        }
    }
}

impl From<&MoyoCell> for Cell {
    fn from(cell: &MoyoCell) -> Self {
        let lattice = &cell.lattice;
        let positions = unsafe {
            std::slice::from_raw_parts(cell.positions, cell.num_atoms as usize)
                .iter()
                .map(|&pos| vector![pos[0], pos[1], pos[2]])
                .collect()
        };
        let numbers = unsafe {
            std::slice::from_raw_parts(cell.numbers, cell.num_atoms as usize)
                .iter()
                .map(|&number| number)
                .collect()
        };
        Cell::new(lattice.into(), positions, numbers)
    }
}

// #[no_mangle]
// pub extern "C" fn moyo_cell_free(cell: *mut MoyoCell) {
//     if cell.is_null() {
//         return;
//     }
//
//     unsafe {
//         drop(Box::from_raw(cell));
//     }
// }
