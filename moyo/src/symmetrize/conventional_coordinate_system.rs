use super::conventional_cell::{
    AXIS_PERMUTATIONS3, monoclinic_candidate_corrections, monoclinic_rank_key,
    orthorhombic_rank_key, select_conventional_correction,
};
use crate::base::{
    Cell, Lattice, Linear, MoyoError, Operations, Transformation, UnimodularTransformation,
};
use crate::data::{HallSymbol, LatticeSystem, arithmetic_crystal_class_entry, hall_symbol_entry};
use crate::identify::SpaceGroup;

pub(super) struct ConventionalCoordinateSystem {
    pub prim_transformation: UnimodularTransformation,
    pub conv_trans_linear: Linear,
    pub transformation: Transformation,
    pub conv_std_operations: Operations,
    pub prim_std_operations: Operations,
}

impl ConventionalCoordinateSystem {
    pub fn new(
        prim_cell: &Cell,
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
        .transform_lattice(&prim_cell.lattice);
        let (prim_transformation, conv_trans_linear) = match lattice_system {
            LatticeSystem::Triclinic => (
                standardize_triclinic_cell(&prim_cell.lattice, &space_group.transformation),
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
