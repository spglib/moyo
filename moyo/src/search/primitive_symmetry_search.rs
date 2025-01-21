use itertools::iproduct;
use std::collections::{HashMap, HashSet, VecDeque};

use log::debug;
use nalgebra::{Matrix3, Vector3};

use super::{
    primitive_cell::PrimitiveMagneticCell,
    solve::{
        pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
        PeriodicKdTree,
    },
    PrimitiveCell,
};
use crate::base::{
    traverse, AngleTolerance, Cell, Lattice, MagneticCell, MagneticMoment, MagneticOperation,
    MagneticOperations, MoyoError, Operation, Operations, Permutation, Rotation,
    RotationMagneticMomentAction, Rotations, Transformation, EPS,
};

#[derive(Debug)]
pub struct PrimitiveSymmetrySearch {
    /// Operations in the given primitive cell
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
}

impl PrimitiveSymmetrySearch {
    /// Return coset representatives of the space group w.r.t. its translation subgroup.
    /// Assume `primitive_cell` is a primitive cell and its basis vectors are Minkowski reduced.
    /// Returned operations are guaranteed to form a group.
    /// If the group closure and tolerance (symprec and angle_tolerance) are incompatible, the former is prioritized.
    /// Possible replacements for spglib/src/spacegroup.h::spa_search_spacegroup
    pub fn new(
        primitive_cell: &Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        // Check if symprec is sufficiently small
        let minimum_basis_norm = primitive_cell.lattice.basis.column(0).norm();
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > minimum_basis_norm / 2.0 {
            debug!(
                "symprec is too large compared to the basis vectors. Consider reducing symprec."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        // Search symmetry operations
        let pkdtree = PeriodicKdTree::new(primitive_cell, rough_symprec);
        let bravais_group =
            search_bravais_group(&primitive_cell.lattice, symprec, angle_tolerance)?;
        let pivot_site_indices = pivot_site_indices(&primitive_cell.numbers);
        let mut symmetries_tmp = vec![];
        let src = pivot_site_indices[0];
        for rotation in bravais_group.iter() {
            let rotated_positions = primitive_cell
                .positions
                .iter()
                .map(|pos| rotation.map(|e| e as f64) * pos)
                .collect::<Vec<_>>();
            for dst in pivot_site_indices.iter() {
                // Try to overlap the `src`-th site to the `dst`-th site
                let translation = primitive_cell.positions[*dst] - rotated_positions[src];
                let new_positions = rotated_positions
                    .iter()
                    .map(|pos| pos + translation)
                    .collect::<Vec<_>>();

                if let Some(permutation) =
                    solve_correspondence(&pkdtree, primitive_cell, &new_positions)
                {
                    symmetries_tmp.push((*rotation, translation, permutation));
                    // Do not break here because there may be multiple translations with the rough tolerance
                }
            }
        }
        debug!(
            "Number of symmetry operation candidates: {}",
            symmetries_tmp.len()
        );

        // Purify symmetry operations by permutations
        let mut operations_and_permutations = vec![];
        for (rotation, rough_translation, permutation) in symmetries_tmp.iter() {
            let (translation, distance) = symmetrize_translation_from_permutation(
                primitive_cell,
                permutation,
                rotation,
                rough_translation,
            );
            if distance < symprec {
                operations_and_permutations
                    .push((Operation::new(*rotation, translation), permutation.clone()));
            }
        }
        if operations_and_permutations.is_empty() {
            debug!(
                "No symmetry operations are found. Consider increasing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooSmallToleranceError);
        }

        // Recover operations by group multiplication
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        let mut operations = vec![];
        let mut permutations = vec![];
        queue.push_back((
            Operation::identity(),
            Permutation::identity(primitive_cell.num_atoms()),
        ));
        while !queue.is_empty() {
            let (ops_lhs, permutation_lhs) = queue.pop_front().unwrap();
            if visited.contains(&ops_lhs.rotation) {
                continue;
            }
            visited.insert(ops_lhs.rotation);
            operations.push(ops_lhs.clone());
            permutations.push(permutation_lhs.clone());

            for (ops_rhs, permutation_rhs) in operations_and_permutations.iter() {
                let mut new_ops = ops_lhs.clone() * ops_rhs.clone();
                new_ops.translation -= new_ops.translation.map(|e| e.round()); // Consider up to the translation subgroup
                let new_permutation = permutation_lhs.clone() * permutation_rhs.clone();
                queue.push_back((new_ops, new_permutation));
            }
        }
        if operations.len() != operations_and_permutations.len() {
            debug!("Found operations do not form a group. Consider reducing symprec and angle_tolerance.");
            return Err(MoyoError::TooLargeToleranceError);
        }

        if !Self::check_closure(&operations, &primitive_cell.lattice, rough_symprec) {
            debug!("Some centering translations are missing. Consider reducing symprec and angle_tolerance.");
            return Err(MoyoError::TooLargeToleranceError);
        }

        debug!("Order of point group: {}", operations.len());
        Ok(Self {
            operations,
            permutations,
        })
    }

    fn check_closure(operations: &Operations, lattice: &Lattice, symprec: f64) -> bool {
        let mut translations_map = HashMap::new();
        for operation in operations.iter() {
            translations_map.insert(operation.rotation, operation.translation);
        }
        for ops1 in operations.iter() {
            for ops2 in operations.iter() {
                let ops12 = ops1.clone() * ops2.clone();
                let diff =
                    (translations_map[&ops12.rotation] - ops12.translation).map(|e| e - e.round());
                if lattice.cartesian_coords(&diff).norm() > symprec {
                    return false;
                }
            }
        }
        true
    }
}

#[derive(Debug)]
pub struct PrimitiveMagneticSymmetrySearch {
    /// Magnetic operations in the given primitive magnetic cell
    pub magnetic_operations: MagneticOperations,
    pub permutations: Vec<Permutation>,
}

impl PrimitiveMagneticSymmetrySearch {
    /// Return coset representatives of the magnetic space group w.r.t. its translation subgroup.
    /// Assume `primitive_magnetic_cell` is a primitive magnetic cell and its basis vectors are Minkowski reduced.
    /// Returned magnetic operations are guaranteed to form a group.
    /// If the group closure and tolerance (symprec, angle_tolerance, and mag_symprec) are incompatible, the former is prioritized.
    pub fn new<M: MagneticMoment>(
        primitive_magnetic_cell: &MagneticCell<M>,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        mag_symprec: f64,
        action: RotationMagneticMomentAction,
    ) -> Result<Self, MoyoError> {
        // Prepare candidate operations from primitive **nonmagnetic** cell
        let prim_nonmag_cell = PrimitiveCell::new(&primitive_magnetic_cell.cell, symprec)?;
        let prim_nonmag_symmetry =
            PrimitiveSymmetrySearch::new(&prim_nonmag_cell.cell, symprec, angle_tolerance)?;
        // Candidate operations in the primitive **magnetic** cell
        let candidate_operations =
            operations_in_cell(&prim_nonmag_cell, &prim_nonmag_symmetry.operations);

        // Find time reversal parts that keep magnetic moments
        let pkdtree = PeriodicKdTree::new(&primitive_magnetic_cell.cell, symprec);
        let mut magnetic_operations = vec![];
        let mut permutations = vec![];
        for operation in candidate_operations.iter() {
            let new_positions = primitive_magnetic_cell
                .cell
                .positions
                .iter()
                .map(|pos| operation.rotation.map(|e| e as f64) * pos + operation.translation)
                .collect::<Vec<_>>();
            if let Some(permutation) =
                solve_correspondence(&pkdtree, &primitive_magnetic_cell.cell, &new_positions)
            {
                let cartesian_rotation =
                    operation.cartesian_rotation(&primitive_magnetic_cell.cell.lattice);
                let rotated_magnetic_moments = primitive_magnetic_cell
                    .magnetic_moments
                    .iter()
                    .map(|m| m.act_rotation(&cartesian_rotation, action))
                    .collect::<Vec<_>>();

                let new_magnetic_moments = (0..primitive_magnetic_cell.num_atoms())
                    .map(|i| primitive_magnetic_cell.magnetic_moments[permutation.apply(i)].clone())
                    .collect::<Vec<_>>();

                // Find time_reversal s.t. time_reversal * rotated_magnetic_moments[:] = new_magnetic_moments[:]
                for time_reversal in [true, false] {
                    let acted_magnetic_moments = rotated_magnetic_moments
                        .iter()
                        .map(|m| m.act_time_reversal(time_reversal))
                        .collect::<Vec<_>>();
                    let take = acted_magnetic_moments
                        .iter()
                        .zip(new_magnetic_moments.iter())
                        .all(|(m1, m2)| m1.is_close(m2, mag_symprec));
                    if take {
                        magnetic_operations.push(MagneticOperation::from_operation(
                            operation.clone(),
                            time_reversal,
                        ));
                        permutations.push(permutation.clone());
                    }
                }
            }
        }

        // Check closure
        if !Self::check_closure(
            &magnetic_operations,
            &primitive_magnetic_cell.cell.lattice,
            symprec,
        ) {
            debug!("Some centering translations are missing. Consider reducing symprec and angle_tolerance.");
            return Err(MoyoError::TooLargeToleranceError);
        }

        Ok(Self {
            magnetic_operations,
            permutations,
        })
    }

    fn check_closure(
        magnetic_operations: &MagneticOperations,
        lattice: &Lattice,
        symprec: f64,
    ) -> bool {
        let mut translations_map = HashMap::new();
        for mops in magnetic_operations.iter() {
            translations_map.insert(
                (mops.operation.rotation, mops.time_reversal),
                mops.operation.translation,
            );
        }
        for mops1 in magnetic_operations.iter() {
            for mops2 in magnetic_operations.iter() {
                let mops12 = mops1.clone() * mops2.clone();
                let diff = (translations_map[&(mops12.operation.rotation, mops12.time_reversal)]
                    - mops12.operation.translation)
                    .map(|e| e - e.round());
                if lattice.cartesian_coords(&diff).norm() > symprec {
                    return false;
                }
            }
        }
        true
    }
}

pub fn operations_in_cell(prim_cell: &PrimitiveCell, prim_operations: &Operations) -> Operations {
    let input_operations =
        Transformation::from_linear(prim_cell.linear).transform_operations(prim_operations);
    let mut operations = vec![];
    for t1 in prim_cell.translations.iter() {
        for operation2 in input_operations.iter() {
            // (E, t1) (rotation, t2) = (rotation, t1 + t2)
            let t12 = (t1 + operation2.translation).map(|e| e % 1.);
            operations.push(Operation::new(operation2.rotation, t12));
        }
    }
    operations
}

pub fn magnetic_operations_in_magnetic_cell<M: MagneticMoment>(
    prim_mag_cell: &PrimitiveMagneticCell<M>,
    prim_mag_operations: &MagneticOperations,
) -> MagneticOperations {
    let input_mag_operations = Transformation::from_linear(prim_mag_cell.linear)
        .transform_magnetic_operations(prim_mag_operations);
    let mut mag_operations = vec![];
    for t1 in prim_mag_cell.translations.iter() {
        for ops2 in input_mag_operations.iter() {
            // (E, t1) (rotation, t2) theta = (rotation, t1 + t2) theta
            let t12 = (t1 + ops2.operation.translation).map(|e| e % 1.);
            mag_operations.push(MagneticOperation::new(
                ops2.operation.rotation,
                t12,
                ops2.time_reversal,
            ));
        }
    }
    mag_operations
}

/// Relevant to spglib.c/symmetry.c::get_lattice_symmetry
fn search_bravais_group(
    minkowski_lattice: &Lattice,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<Rotations, MoyoError> {
    // Candidate column vectors for rotation matrix
    let lengths = minkowski_lattice
        .basis
        .column_iter()
        .map(|v| v.norm())
        .collect::<Vec<_>>();
    let mut candidate_lattice_points = [vec![], vec![], vec![]];
    // It would be sufficient to search coeffs in [-1, 0, 1] because the third column of the Minkowski lattice is already in [-1, 0, 1]^3
    for coeffs in iproduct!(-1..=1, -1..=1, -1..=1) {
        let v = minkowski_lattice.basis
            * Vector3::new(coeffs.0 as f64, coeffs.1 as f64, coeffs.2 as f64);
        let v_length = v.norm();
        for (i, &length) in lengths.iter().enumerate() {
            if (v_length - length).abs() < symprec {
                candidate_lattice_points[i].push(coeffs);
            }
        }
    }

    let mut rotations = vec![];
    for (c0, c1) in iproduct!(
        candidate_lattice_points[0].iter(),
        candidate_lattice_points[1].iter()
    ) {
        let v0 = minkowski_lattice.basis * Vector3::new(c0.0 as f64, c0.1 as f64, c0.2 as f64);
        let v1 = minkowski_lattice.basis * Vector3::new(c1.0 as f64, c1.1 as f64, c1.2 as f64);

        // Check the (0, 1)-th element of metric tensor beforehand
        if !compare_nondiagonal_matrix_tensor_element(
            &minkowski_lattice.basis,
            &v0,
            &v1,
            0,
            1,
            symprec,
            angle_tolerance,
        ) {
            continue;
        }

        for c2 in candidate_lattice_points[2].iter() {
            let v2 = minkowski_lattice.basis * Vector3::new(c2.0 as f64, c2.1 as f64, c2.2 as f64);

            // Check determinant
            let rotation = Rotation::from_columns(&[
                Vector3::new(c0.0, c0.1, c0.2),
                Vector3::new(c1.0, c1.1, c1.2),
                Vector3::new(c2.0, c2.1, c2.2),
            ]);

            if relative_ne!(
                rotation.map(|e| e as f64).determinant().abs(),
                1.0,
                epsilon = EPS
            ) {
                continue;
            }

            // Check the (1, 2)th elements of metric tensor
            if !compare_nondiagonal_matrix_tensor_element(
                &minkowski_lattice.basis,
                &v1,
                &v2,
                1,
                2,
                symprec,
                angle_tolerance,
            ) {
                continue;
            }
            // Check the (2, 0)th elements of metric tensor
            if !compare_nondiagonal_matrix_tensor_element(
                &minkowski_lattice.basis,
                &v2,
                &v0,
                2,
                0,
                symprec,
                angle_tolerance,
            ) {
                continue;
            }

            rotations.push(rotation);
        }
    }

    // 48 for Oh
    if rotations.is_empty() || (48 % rotations.len() != 0) {
        debug!("Found automorphisms for the lattice do not form a group. Consider reducing symprec and angle_tolerance.");
        return Err(MoyoError::TooLargeToleranceError);
    }

    let complemented_rotations = traverse(&rotations);
    if complemented_rotations.len() != rotations.len() {
        debug!("Found automorphisms for the lattice do not form a group. Consider reducing symprec and angle_tolerance.");
        return Err(MoyoError::TooLargeToleranceError);
    }
    debug!("Order of Bravais group: {}", complemented_rotations.len());
    Ok(complemented_rotations)
}

/// Compare (basis.column(col1), basis.column(col2)) and (b1, b2)
/// Return true if these pairs of vectors are close enough
fn compare_nondiagonal_matrix_tensor_element(
    basis: &Matrix3<f64>,
    b1: &Vector3<f64>,
    b2: &Vector3<f64>,
    col1: usize,
    col2: usize,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> bool {
    let theta_org = basis.column(col1).angle(&basis.column(col2));
    let theta_new = b1.angle(b2);
    let cos_dtheta = theta_org.cos() * theta_new.cos() + theta_org.sin() * theta_new.sin();

    match angle_tolerance {
        AngleTolerance::Radian(angle_tolerance) => cos_dtheta.acos().abs() < angle_tolerance,
        AngleTolerance::Default => {
            // Eq.(7) of https://arxiv.org/pdf/1808.01590.pdf
            let sin_dtheta2 = 1.0 - cos_dtheta.powi(2);
            let length_ave2 = (basis.column(col1).norm() + b1.norm())
                * (basis.column(col2).norm() + b2.norm())
                / 4.0;
            sin_dtheta2 * length_ave2 < symprec * symprec
        }
    }
}

#[cfg(test)]
mod tests {
    use std::vec;

    use nalgebra::{matrix, Matrix3, Vector3};
    use test_log::test;

    use super::{search_bravais_group, PrimitiveMagneticSymmetrySearch};
    use crate::base::{
        AngleTolerance, Collinear, Lattice, MagneticCell, NonCollinear,
        RotationMagneticMomentAction,
    };

    #[test]
    fn test_search_bravais_group() {
        let symprec = 1e-4;

        {
            // m-3m
            let lattice = Lattice::new(matrix![
                0.0, 0.5, 0.5;
                0.5, 0.0, 0.5;
                0.5, 0.5, 0.0;
            ]);
            let rotations =
                search_bravais_group(&lattice, symprec, AngleTolerance::Radian(1e-2)).unwrap();
            assert_eq!(rotations.len(), 48);
        }

        {
            // 6/mmm
            let lattice = Lattice::new(matrix![
                1.0, 0.0, 0.0;
                -0.5, f64::sqrt(3.0) / 2.0, 0.0;
                0.0, 0.0, 1.0;
            ]);
            let rotations =
                search_bravais_group(&lattice, symprec, AngleTolerance::Default).unwrap();
            assert_eq!(rotations.len(), 24);
        }

        {
            let lattice = Lattice::new(matrix![
                0.5, 0.0, 0.5;
                0.5, 0.5, 0.0;
                0.0, 0.5, 0.5;
            ]);
            let rotations =
                search_bravais_group(&lattice, symprec, AngleTolerance::Default).unwrap();
            assert_eq!(rotations.len(), 48);
        }
    }

    #[test]
    fn test_primitive_magnetic_symmetry_search() {
        let symprec = 1e-4;
        let angle_tolerance = AngleTolerance::Default;
        let mag_symprec = 1e-4;

        {
            // AFM bcc (type 4)
            let prim_mag_cell = MagneticCell::new(
                Lattice::new(Matrix3::identity()),
                vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.5, 0.5, 0.5)],
                vec![0, 0],
                vec![Collinear(1.0), Collinear(-1.0)],
            );
            let result = PrimitiveMagneticSymmetrySearch::new(
                &prim_mag_cell,
                symprec,
                angle_tolerance,
                mag_symprec,
                RotationMagneticMomentAction::Axial,
            )
            .unwrap();
            assert_eq!(result.magnetic_operations.len(), 96);
            assert_eq!(result.permutations.len(), result.magnetic_operations.len());
            assert_eq!(
                result
                    .magnetic_operations
                    .iter()
                    .filter(|mops| mops.time_reversal)
                    .count(),
                48
            );
        }
        {
            // Primitive fcc (type 2)
            let prim_mag_cell = MagneticCell::new(
                Lattice::new(matrix![
                    0.0, 0.5, 0.5;
                    0.5, 0.0, 0.5;
                    0.5, 0.5, 0.0;
                ]),
                vec![Vector3::new(0.0, 0.0, 0.0)],
                vec![0],
                vec![NonCollinear(Vector3::new(0.0, 0.0, 0.0))],
            );
            let result = PrimitiveMagneticSymmetrySearch::new(
                &prim_mag_cell,
                symprec,
                angle_tolerance,
                mag_symprec,
                RotationMagneticMomentAction::Axial,
            )
            .unwrap();
            assert_eq!(result.magnetic_operations.len(), 96);
            assert_eq!(result.permutations.len(), result.magnetic_operations.len());
            assert_eq!(
                result
                    .magnetic_operations
                    .iter()
                    .filter(|mops| mops.time_reversal)
                    .count(),
                48
            );
        }
    }
}
