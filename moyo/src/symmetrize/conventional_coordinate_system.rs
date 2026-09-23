use super::conventional_cell::{
    AXIS_PERMUTATIONS3, monoclinic_candidate_corrections, monoclinic_rank_key,
    orthorhombic_rank_key, select_conventional_correction,
};
use crate::base::{
    Lattice, Linear, MoyoError, Operations, Transformation, UnimodularTransformation,
};
use crate::data::{HallSymbol, LatticeSystem, arithmetic_crystal_class_entry, hall_symbol_entry};
use crate::identify::SpaceGroup;

/// A choice of conventional basis and origin for an identified space-group setting.
///
/// This is the first stage of standardization. It selects coordinate transformations
/// from the input primitive lattice and the identified space group. Position
/// refinement, Cartesian orientation, and construction of the output cells belong
/// to [`super::standardize::StandardizedCell`].
///
/// Both transformations use the column-vector convention `A_new = A_input * P`
/// and `x_new = P^-1 * (x_input - p)`, where `p` is expressed in the input basis.
/// The selected primitive and conventional systems share an origin and are related
/// by the fixed centering matrix `C = conv_trans_linear`: `P_conv = P_prim * C`.
/// The operation lists describe the target symmetry in these coordinate systems.
pub(super) struct ConventionalCoordinateSystem {
    /// From the input primitive system to the selected primitive system.
    pub prim_transformation: UnimodularTransformation,
    /// From the selected primitive basis to the conventional basis, with zero origin shift.
    pub conv_trans_linear: Linear,
    /// From the input primitive system to the selected conventional system.
    pub transformation: Transformation,
    /// Target operations in the selected conventional system, in Hall traversal order.
    pub conv_std_operations: Operations,
    /// Target operations in the selected primitive system, in Hall traversal order.
    pub prim_std_operations: Operations,
}

impl ConventionalCoordinateSystem {
    /// Select a basis and origin while keeping the identified Hall setting.
    ///
    /// Triclinic bases are Niggli-reduced. Monoclinic and orthorhombic bases are
    /// ranked among admissible affine-normalizer corrections; other lattice
    /// systems use the identified transformation. `epsilon` is the tolerance for
    /// checking centering translations and matching origin shifts.
    pub fn new(
        prim_lattice: &Lattice,
        space_group: &SpaceGroup,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        let entry =
            hall_symbol_entry(space_group.hall_number).ok_or(MoyoError::StandardizationError)?;

        // Prepare operations in primitive standard
        let hs = HallSymbol::from_hall_number(space_group.hall_number)
            .ok_or(MoyoError::StandardizationError)?;
        let (conv_std_operations, prim_std_operations) = hs.traverse_and_primitive_traverse();

        // To standardized primitive cell
        let lattice_system = arithmetic_crystal_class_entry(entry.arithmetic_number)
            .unwrap()
            .lattice_system();
        // For monoclinic and orthorhombic systems, the identified setting leaves some
        // freedom in the conventional basis (the affine normalizer). The chosen
        // correction is folded into `prim_transformation` so that the primitive
        // standardized cell is always related to the conventional one by the fixed
        // centering matrix `entry.centering.linear()`.
        let conv_lattice_tmp = Transformation::from_linear(
            space_group.transformation.linear * entry.centering.linear(),
        )
        .transform_lattice(prim_lattice);
        let (prim_transformation, conv_trans_linear) = match lattice_system {
            LatticeSystem::Triclinic => (
                standardize_triclinic_cell(prim_lattice, &space_group.transformation),
                Linear::identity(),
            ),
            LatticeSystem::Monoclinic => {
                let prim_correction = select_conventional_correction(
                    &conv_lattice_tmp,
                    entry.centering,
                    &prim_std_operations,
                    &hs.primitive_generators(),
                    &monoclinic_candidate_corrections(&conv_lattice_tmp),
                    monoclinic_rank_key,
                    epsilon,
                );
                (
                    space_group.transformation.clone() * prim_correction,
                    entry.centering.linear(),
                )
            }
            LatticeSystem::Orthorhombic => {
                let prim_correction = select_conventional_correction(
                    &conv_lattice_tmp,
                    entry.centering,
                    &prim_std_operations,
                    &hs.primitive_generators(),
                    &AXIS_PERMUTATIONS3,
                    orthorhombic_rank_key,
                    epsilon,
                );
                (
                    space_group.transformation.clone() * prim_correction,
                    entry.centering.linear(),
                )
            }
            _ => (space_group.transformation.clone(), entry.centering.linear()),
        };

        // prim_transformation * (conv_trans_linear, 0)
        let transformation = Transformation::new(
            prim_transformation.linear * conv_trans_linear,
            prim_transformation.origin_shift,
        );

        Ok(Self {
            prim_transformation,
            conv_trans_linear,
            transformation,
            conv_std_operations,
            prim_std_operations,
        })
    }
}

/// Niggli reduction for distorted triclinic lattice systems is numerically so challenging.
/// Thus, we skip checking reduction condition.
fn standardize_triclinic_cell(
    lattice: &Lattice,
    transformation_to_prim_std: &UnimodularTransformation,
) -> UnimodularTransformation {
    let lattice_prim_std_tmp = transformation_to_prim_std.transform_lattice(lattice);
    let (_, niggli_linear) = lattice_prim_std_tmp.unchecked_niggli_reduce();
    UnimodularTransformation::new(
        niggli_linear * transformation_to_prim_std.linear,
        transformation_to_prim_std.origin_shift,
    )
}

#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, matrix, vector};

    use super::ConventionalCoordinateSystem;
    use crate::base::{Lattice, UnimodularTransformation};
    use crate::data::{Centering, HallSymbol};
    use crate::identify::SpaceGroup;

    #[test]
    fn test_imma_basis_and_origin_selection() {
        // Imma permits swapping a and b only together with an origin shift.
        // Select the shorter a axis using only a lattice and an identified group.
        let conv_lattice = Lattice::new(matrix![
            7.0, 0.0, 0.0;
            0.0, 3.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let reference_prim_lattice =
            Lattice::new((conv_lattice.basis * Centering::I.inverse()).transpose());
        let identified_transformation = UnimodularTransformation::new(
            matrix![1, 1, 0; 0, 1, 0; 0, 0, 1],
            vector![0.17, 0.23, 0.31],
        );
        let space_group =
            SpaceGroup::from_hall_number_and_transformation(343, identified_transformation)
                .unwrap();
        let to_input = space_group.transformation.inverse();
        let prim_lattice = to_input.transform_lattice(&reference_prim_lattice);
        let prim_operations = to_input.transform_operations(
            &HallSymbol::from_hall_number(343)
                .unwrap()
                .primitive_traverse(),
        );

        let selected =
            ConventionalCoordinateSystem::new(&prim_lattice, &space_group, 1e-8).unwrap();
        let conv_basis = selected
            .transformation
            .transform_lattice(&prim_lattice)
            .basis;
        assert_relative_eq!(conv_basis.column(0).norm(), 3.0, epsilon = 1e-8);
        assert_relative_eq!(conv_basis.column(1).norm(), 7.0, epsilon = 1e-8);
        assert_relative_eq!(conv_basis.column(2).norm(), 5.0, epsilon = 1e-8);

        let correction = to_input * selected.prim_transformation.clone();
        assert!(correction.origin_shift.norm() > 1e-3);
        assert_eq!(selected.conv_trans_linear, Centering::I.linear());
        assert_eq!(
            selected.transformation.linear,
            selected.prim_transformation.linear * selected.conv_trans_linear
        );
        assert_relative_eq!(
            selected.transformation.origin_shift,
            selected.prim_transformation.origin_shift
        );

        let transformed = selected
            .prim_transformation
            .transform_operations(&prim_operations);
        assert_eq!(transformed.len(), selected.prim_std_operations.len());
        for operation in &transformed {
            let target = selected
                .prim_std_operations
                .iter()
                .find(|target| target.rotation == operation.rotation)
                .unwrap();
            let mut diff = operation.translation - target.translation;
            diff -= diff.map(|x| x.round());
            assert_relative_eq!(diff, Vector3::zeros(), epsilon = 1e-8);
        }
    }
}
