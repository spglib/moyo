use itertools::izip;
use nalgebra::Matrix3;

use super::standardize::StandardizedCell;
use crate::base::{
    MagneticCell, MagneticMoment, MoyoError, Operations, Permutation, RotationMagneticMomentAction,
    Transformation, UnimodularTransformation,
};
use crate::data::{get_magnetic_space_group_type, ConstructType};
use crate::identify::{
    family_space_group_from_magnetic_space_group,
    primitive_maximal_space_subgroup_from_magnetic_space_group, MagneticSpaceGroup,
};
use crate::search::{PrimitiveMagneticCell, PrimitiveMagneticSymmetrySearch};

pub struct StandardizedMagneticCell<M: MagneticMoment> {
    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    pub prim_mag_cell: MagneticCell<M>,
    /// Transformation from the input primitive magnetic cell to the primitive standardized magnetic cell.
    pub prim_transformation: UnimodularTransformation,
    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    pub mag_cell: MagneticCell<M>,
    pub transformation: Transformation,
    // ------------------------------------------------------------------------
    // Miscellaneous
    // ------------------------------------------------------------------------
    pub rotation_matrix: Matrix3<f64>,
    /// Mapping from the site in the `mag_cell` to that in the `prim_mag_cell`
    #[allow(dead_code)]
    pub site_mapping: Vec<usize>,
}

impl<M: MagneticMoment> StandardizedMagneticCell<M> {
    /// Standardize the input **primitive** magnetic cell.
    /// For triclinic magnetic space groups, Niggli reduction is performed.
    /// Basis vectors are rotated to be a upper triangular matrix.
    pub fn new(
        prim_mag_cell: &PrimitiveMagneticCell<M>,
        magnetic_symmetry_search: &PrimitiveMagneticSymmetrySearch,
        magnetic_space_group: &MagneticSpaceGroup,
        symprec: f64,
        mag_symprec: f64,
        epsilon: f64,
        action: RotationMagneticMomentAction,
    ) -> Result<Self, MoyoError> {
        // Symmetrize positions and lattice by reference space group
        let ref_space_group = magnetic_space_group.reference_space_group();
        let msg_type = get_magnetic_space_group_type(magnetic_space_group.uni_number)
            .ok_or(MoyoError::MagneticStandardizationError)?;
        let (ref_prim_operations, ref_prim_permutations) =
            Self::reference_symmetry_operations_and_permutations(
                magnetic_symmetry_search,
                msg_type.construct_type,
                mag_symprec,
            );
        let ref_std_cell = StandardizedCell::new(
            &prim_mag_cell.magnetic_cell.cell,
            &ref_prim_operations,
            &ref_prim_permutations,
            &ref_space_group,
            symprec,
            epsilon,
        )?;

        // Need to rotate magnetic moments because the standardization rotates the ref cell.
        let prim_std_magnetic_moments_tmp = (0..prim_mag_cell.magnetic_cell.magnetic_moments.len())
            .map(|i| {
                prim_mag_cell.magnetic_cell.magnetic_moments[ref_std_cell.site_mapping[i]]
                    .act_rotation(&ref_std_cell.rotation_matrix, action)
            })
            .collect::<Vec<_>>();
        let cart_rotations = magnetic_symmetry_search
            .magnetic_operations
            .iter()
            .map(|mops| {
                mops.operation
                    .cartesian_rotation(&prim_mag_cell.magnetic_cell.cell.lattice)
            })
            .collect::<Vec<_>>();
        let time_reversals = magnetic_symmetry_search
            .magnetic_operations
            .iter()
            .map(|mops| mops.time_reversal)
            .collect::<Vec<_>>();

        Self::new_from_ref_cell(
            &ref_std_cell,
            &prim_std_magnetic_moments_tmp,
            &cart_rotations,
            &time_reversals,
            &magnetic_symmetry_search.permutations,
            action,
        )
    }

    fn new_from_ref_cell(
        ref_std_cell: &StandardizedCell,
        prim_std_magnetic_moments_tmp: &[M],
        cart_rotations: &[Matrix3<f64>],
        time_reversals: &[bool],
        permutations: &[Permutation],
        action: RotationMagneticMomentAction,
    ) -> Result<Self, MoyoError> {
        // Symmetrize magnetic moments by magnetic space group
        let prim_std_magnetic_moments = Self::symmetrize_magnetic_moments(
            &prim_std_magnetic_moments_tmp,
            &cart_rotations,
            &time_reversals,
            &permutations,
            action,
        );
        let prim_std_mag_cell =
            MagneticCell::from_cell(ref_std_cell.prim_cell.clone(), prim_std_magnetic_moments);

        // To (conventional) standardized magnetic cell
        let refined_prim_mag_cell = ref_std_cell
            .prim_transformation
            .inverse()
            .transform_magnetic_cell(&prim_std_mag_cell);
        let (conv_std_mag_cell, _) = ref_std_cell
            .transformation
            .transform_magnetic_cell(&refined_prim_mag_cell);

        Ok(Self {
            // Primitive standardized magnetic cell
            prim_mag_cell: prim_std_mag_cell,
            prim_transformation: ref_std_cell.prim_transformation.clone(),
            // Standardized magnetic cell
            mag_cell: conv_std_mag_cell,
            transformation: ref_std_cell.transformation.clone(),
            // Miscellaneous
            rotation_matrix: ref_std_cell.rotation_matrix,
            site_mapping: ref_std_cell.site_mapping.clone(),
        })
    }

    fn reference_symmetry_operations_and_permutations(
        magnetic_symmetry_search: &PrimitiveMagneticSymmetrySearch,
        construct_type: ConstructType,
        epsilon: f64,
    ) -> (Operations, Vec<Permutation>) {
        let (ref_operations, contained) = match construct_type {
            ConstructType::Type1 | ConstructType::Type2 | ConstructType::Type3 => {
                let (fsg, _, contained) = family_space_group_from_magnetic_space_group(
                    &magnetic_symmetry_search.magnetic_operations,
                    epsilon,
                );
                (fsg, contained)
            }
            ConstructType::Type4 => {
                let (xsg, contained) = primitive_maximal_space_subgroup_from_magnetic_space_group(
                    &magnetic_symmetry_search.magnetic_operations,
                );
                (xsg, contained)
            }
        };
        let ref_permutations = magnetic_symmetry_search
            .permutations
            .iter()
            .enumerate()
            .filter_map(|(i, perm)| {
                if contained[i] {
                    Some(perm.clone())
                } else {
                    None
                }
            })
            .collect();
        (ref_operations, ref_permutations)
    }

    fn symmetrize_magnetic_moments(
        magnetic_moments: &[M],
        cart_rotations: &[Matrix3<f64>],
        time_reversals: &[bool],
        permutations: &[Permutation],
        action: RotationMagneticMomentAction,
    ) -> Vec<M> {
        let inverse_permutations = permutations
            .iter()
            .map(|permutation| permutation.inverse())
            .collect::<Vec<_>>();
        (0..magnetic_moments.len())
            .map(|i| {
                let mut equiv_magmoms = vec![];
                for (inv_perm, cart_rotation, time_reversal) in izip!(
                    inverse_permutations.iter(),
                    cart_rotations.iter(),
                    time_reversals.iter()
                ) {
                    let magmom = magnetic_moments[inv_perm.apply(i)].act_magnetic_operation(
                        &cart_rotation,
                        *time_reversal,
                        action,
                    );
                    equiv_magmoms.push(magmom)
                }
                MagneticMoment::average(&equiv_magmoms)
            })
            .collect()
    }
}
