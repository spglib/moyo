use moyo::base::{Collinear, MagneticCell, NonCollinear};
use moyo::utils::{to_3_slice, to_vector3};

use super::cell::{MoyoCell, moyo_cell_free_members};
use crate::ffi::{free_slice, leak_slice};

/// Collinear magnetic structure.
/// All pointer fields are owned by the containing dataset and are freed by
/// `moyo_collinear_magnetic_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoCollinearMagneticCell {
    /// Crystal structure
    pub cell: MoyoCell,
    /// Collinear magnetic moments of the `cell.num_atoms` sites
    pub magnetic_moments: *const f64,
}

impl From<&MagneticCell<Collinear>> for MoyoCollinearMagneticCell {
    fn from(magnetic_cell: &MagneticCell<Collinear>) -> Self {
        let magnetic_moments = magnetic_cell
            .magnetic_moments
            .iter()
            .map(|m| m.0)
            .collect::<Vec<_>>();
        Self {
            cell: (&magnetic_cell.cell).into(),
            magnetic_moments: leak_slice(magnetic_moments),
        }
    }
}

impl From<&MoyoCollinearMagneticCell> for MagneticCell<Collinear> {
    fn from(magnetic_cell: &MoyoCollinearMagneticCell) -> Self {
        let cell = (&magnetic_cell.cell).into();
        let magnetic_moments = unsafe {
            std::slice::from_raw_parts(
                magnetic_cell.magnetic_moments,
                magnetic_cell.cell.num_atoms as usize,
            )
            .iter()
            .map(|&m| Collinear(m))
            .collect()
        };
        MagneticCell::from_cell(cell, magnetic_moments)
    }
}

pub(crate) unsafe fn moyo_collinear_magnetic_cell_free_members(
    magnetic_cell: &MoyoCollinearMagneticCell,
) {
    unsafe {
        moyo_cell_free_members(&magnetic_cell.cell);
        free_slice(
            magnetic_cell.magnetic_moments,
            magnetic_cell.cell.num_atoms as usize,
        );
    }
}

/// Non-collinear magnetic structure.
/// All pointer fields are owned by the containing dataset and are freed by
/// `moyo_noncollinear_magnetic_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoNonCollinearMagneticCell {
    /// Crystal structure
    pub cell: MoyoCell,
    /// Non-collinear magnetic moments of the `cell.num_atoms` sites in Cartesian coordinates
    pub magnetic_moments: *const [f64; 3],
}

impl From<&MagneticCell<NonCollinear>> for MoyoNonCollinearMagneticCell {
    fn from(magnetic_cell: &MagneticCell<NonCollinear>) -> Self {
        let magnetic_moments = magnetic_cell
            .magnetic_moments
            .iter()
            .map(|m| to_3_slice(&m.0))
            .collect::<Vec<_>>();
        Self {
            cell: (&magnetic_cell.cell).into(),
            magnetic_moments: leak_slice(magnetic_moments),
        }
    }
}

impl From<&MoyoNonCollinearMagneticCell> for MagneticCell<NonCollinear> {
    fn from(magnetic_cell: &MoyoNonCollinearMagneticCell) -> Self {
        let cell = (&magnetic_cell.cell).into();
        let magnetic_moments = unsafe {
            std::slice::from_raw_parts(
                magnetic_cell.magnetic_moments,
                magnetic_cell.cell.num_atoms as usize,
            )
            .iter()
            .map(|m| NonCollinear(to_vector3(m)))
            .collect()
        };
        MagneticCell::from_cell(cell, magnetic_moments)
    }
}

pub(crate) unsafe fn moyo_noncollinear_magnetic_cell_free_members(
    magnetic_cell: &MoyoNonCollinearMagneticCell,
) {
    unsafe {
        moyo_cell_free_members(&magnetic_cell.cell);
        free_slice(
            magnetic_cell.magnetic_moments,
            magnetic_cell.cell.num_atoms as usize,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use moyo::base::Lattice;
    use nalgebra::{matrix, vector};

    #[test]
    fn test_roundtrip_moyo_collinear_magnetic_cell() {
        let original = MagneticCell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                1.0, 1.0, 0.0;
                1.0, 0.0, 1.0
            ]),
            vec![vector![0.0, 0.0, 0.0], vector![0.5, 0.5, 0.5]],
            vec![1, 1],
            vec![Collinear(0.7), Collinear(-0.7)],
        );
        let moyoc = MoyoCollinearMagneticCell::from(&original);
        let reconstructed = MagneticCell::<Collinear>::from(&moyoc);
        assert_eq!(original.num_atoms(), reconstructed.num_atoms());
        assert_relative_eq!(
            original.cell.lattice.basis,
            reconstructed.cell.lattice.basis
        );
        for (m1, m2) in original
            .magnetic_moments
            .iter()
            .zip(reconstructed.magnetic_moments.iter())
        {
            assert_relative_eq!(m1.0, m2.0);
        }
        unsafe { moyo_collinear_magnetic_cell_free_members(&moyoc) };
    }

    #[test]
    fn test_roundtrip_moyo_noncollinear_magnetic_cell() {
        let original = MagneticCell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                1.0, 1.0, 0.0;
                1.0, 0.0, 1.0
            ]),
            vec![vector![0.0, 0.0, 0.0], vector![0.5, 0.5, 0.5]],
            vec![1, 1],
            vec![
                NonCollinear(vector![0.0, 0.0, 0.7]),
                NonCollinear(vector![0.0, 0.0, -0.7]),
            ],
        );
        let moyoc = MoyoNonCollinearMagneticCell::from(&original);
        let reconstructed = MagneticCell::<NonCollinear>::from(&moyoc);
        assert_eq!(original.num_atoms(), reconstructed.num_atoms());
        assert_relative_eq!(
            original.cell.lattice.basis,
            reconstructed.cell.lattice.basis
        );
        for (m1, m2) in original
            .magnetic_moments
            .iter()
            .zip(reconstructed.magnetic_moments.iter())
        {
            assert_relative_eq!(m1.0, m2.0);
        }
        unsafe { moyo_noncollinear_magnetic_cell_free_members(&moyoc) };
    }
}
