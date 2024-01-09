use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{vector, Matrix3};

use crate::base::{
    Cell, Lattice, MoyoError, Rotation, Transformation, UnimodularTransformation, EPS,
};
use crate::data::{arithmetic_crystal_class_entry, hall_symbol_entry, HallSymbol, LatticeSystem};
use crate::identify::SpaceGroup;
use crate::search::PrimitiveCell;

pub struct StandardizedCell {
    pub prim_cell: Cell,
    /// Transformation from the input primitive cell to the primitive standardized cell.
    pub prim_transformation: UnimodularTransformation,
    pub cell: Cell,
    /// Transformation from the input primitive cell to the standardized cell.
    pub transformation: Transformation,
    /// Rotation matrix to map the lattice of the input primitive cell to that of the standardized cell.
    pub rotation_matrix: Matrix3<f64>,
}

impl StandardizedCell {
    /// Standardize the input **primitive** cell.
    /// For triclinic space groups, Niggli reduction is performed.
    /// Basis vectors are rotated to be a upper triangular matrix.
    pub fn new(prim_cell: &PrimitiveCell, space_group: &SpaceGroup) -> Result<Self, MoyoError> {
        let entry = hall_symbol_entry(space_group.hall_number);
        let arithmetic_number = entry.arithmetic_number;
        let bravais_class = arithmetic_crystal_class_entry(arithmetic_number).bravais_class;
        let lattice_system = LatticeSystem::from_bravais_class(bravais_class);

        // To standardized primitive cell
        let prim_transformation = match lattice_system {
            LatticeSystem::Triclinic => {
                let (_, linear) = prim_cell.cell.lattice.niggli_reduce()?;
                UnimodularTransformation::from_linear(linear)
            }
            _ => space_group.transformation.clone(),
        };
        let prim_std_cell = prim_transformation.transform_cell(&prim_cell.cell);

        // To (conventional) standardized cell
        let conv_trans = Transformation::from_linear(entry.centering.linear());
        let std_cell = conv_trans.transform_cell(&prim_std_cell);

        // Symmetrize lattice
        let std_rotations = HallSymbol::from_hall_number(entry.hall_number)
            .traverse()
            .rotations;
        let (_, rotation_matrix) = symmetrize_lattice(&std_cell.lattice, &std_rotations);

        // TODO: match Wyckoff positions

        Ok(StandardizedCell {
            prim_cell: prim_std_cell.rotate(&rotation_matrix),
            prim_transformation: prim_transformation.clone(),
            cell: std_cell.rotate(&rotation_matrix),
            // prim_transformation * (conv_trans.linear, 0)
            transformation: Transformation::new(
                prim_transformation.linear * conv_trans.linear,
                prim_transformation.origin_shift,
            ),
            rotation_matrix,
        })
    }
}

fn symmetrize_lattice(lattice: &Lattice, rotations: &Vec<Rotation>) -> (Lattice, Matrix3<f64>) {
    let metric_tensor = lattice.metric_tensor();
    let mut symmetrized_metric_tensor: Matrix3<f64> = rotations
        .iter()
        .map(|rotation| {
            rotation.transpose().map(|e| e as f64) * metric_tensor * rotation.map(|e| e as f64)
        })
        .sum();
    symmetrized_metric_tensor /= rotations.len() as f64;

    // upper-triangular basis
    let tri_basis = Cholesky::new_unchecked(symmetrized_metric_tensor)
        .l()
        .transpose();

    // tri_basis \approx rotation_matrix * lattice.basis
    // QR(tri_basis * lattice.basis^-1) = rotation_matrix * strain
    let qr = QR::new(tri_basis * lattice.basis.try_inverse().unwrap());
    let r = qr.r();
    let signs =
        Matrix3::<f64>::from_diagonal(&vector![sign(r[(0, 0)]), sign(r[(1, 1)]), sign(r[(2, 2)])]);
    let mut rotation_matrix = QR::new(tri_basis * lattice.basis.try_inverse().unwrap()).q();
    // Remove axis-direction freedom
    rotation_matrix *= signs;
    (Lattice::new(tri_basis), rotation_matrix)
}

fn sign(x: f64) -> f64 {
    if x > EPS {
        1.0
    } else if x < -EPS {
        -1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::symmetrize_lattice;
    use crate::base::{traverse, Lattice};
    use crate::data::{GeometricCrystalClass, PointGroupRepresentative};

    #[test]
    fn test_symmetrize_lattice_cubic() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0001;
            0.0, -0.999, 0.0;
            0.0, 0.0, -1.0001;
        ]);
        let rep = PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::Oh);
        let rotations = traverse(&rep.generators);

        let (new_lattice, rotation_matrix) = symmetrize_lattice(&lattice, &rotations);
        assert_relative_eq!(new_lattice.basis[(1, 1)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(2, 2)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(0, 1)], 0.0);
        assert_relative_eq!(new_lattice.basis[(0, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 1)], 0.0);

        assert_relative_eq!(
            rotation_matrix * lattice.basis,
            new_lattice.basis,
            epsilon = 1e-2
        );
    }
}
