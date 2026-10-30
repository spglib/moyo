use itertools::iproduct;
use std::collections::{HashMap, VecDeque};

use log::debug;
use nalgebra::{Matrix3, Vector3};
use rustc_hash::{FxHashMap, FxHashSet};

use super::{
    PrimitiveCell,
    primitive_cell::PrimitiveMagneticCell,
    solve::{
        PeriodicNeighborSearch, pivot_site_indices, solve_correspondence,
        solve_correspondence_with, symmetrize_translation_from_permutation,
    },
};
use crate::base::{
    AngleTolerance, Cell, EPS, Lattice, Linear, MagneticCell, MagneticMoment, MagneticOperation,
    MagneticOperations, MoyoError, Operation, Operations, Permutation, Rotation,
    RotationMagneticMomentAction, Rotations, Transformation, Translation,
    is_angle_within_tolerance, traverse,
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
        let neighbor_search = PeriodicNeighborSearch::new(primitive_cell, rough_symprec);
        let bravais_group =
            search_bravais_group(&primitive_cell.lattice, symprec, angle_tolerance)?;
        let pivot_site_indices = pivot_site_indices(&primitive_cell.numbers);
        let mut symmetries_tmp = vec![];
        let src = pivot_site_indices[0];
        let mut mapping = Vec::with_capacity(primitive_cell.num_atoms());
        let mut visited = vec![false; primitive_cell.num_atoms()];
        for rotation in bravais_group.iter() {
            let rotated_positions = primitive_cell
                .positions
                .iter()
                .map(|pos| rotation.map(|e| e as f64) * pos)
                .collect::<Vec<_>>();
            for dst in pivot_site_indices.iter() {
                // Try to overlap the `src`-th site to the `dst`-th site
                let translation = primitive_cell.positions[*dst] - rotated_positions[src];
                if let Some(permutation) = solve_correspondence_with(
                    &neighbor_search,
                    primitive_cell,
                    |i| rotated_positions[i] + translation,
                    &mut mapping,
                    &mut visited,
                ) {
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
        let mut visited = FxHashSet::default();
        visited.insert(Rotation::identity());
        let mut operations = vec![];
        let mut permutations = vec![];
        queue.push_back((
            Operation::identity(),
            Permutation::identity(primitive_cell.num_atoms()),
        ));
        while !queue.is_empty() {
            let (ops_lhs, permutation_lhs) = queue.pop_front().unwrap();
            operations.push(ops_lhs.clone());
            permutations.push(permutation_lhs.clone());

            // Expand identity first to enqueue every candidate. Once all expected
            // rotations are queued, keep their BFS order and let check_closure
            // validate every product, including missing rotations.
            if operations.len() > 1 && visited.len() == operations_and_permutations.len() {
                continue;
            }

            for (ops_rhs, permutation_rhs) in operations_and_permutations.iter() {
                let mut new_ops = ops_lhs.clone() * ops_rhs.clone();
                if !visited.insert(new_ops.rotation) {
                    continue;
                }
                new_ops.translation -= new_ops.translation.map(|e| e.round()); // Consider up to the translation subgroup
                let new_permutation = permutation_lhs.clone() * permutation_rhs.clone();
                queue.push_back((new_ops, new_permutation));
            }
        }
        if operations.len() != operations_and_permutations.len() {
            debug!(
                "Found operations do not form a group. Consider reducing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        if !Self::check_closure(&operations, &primitive_cell.lattice, rough_symprec) {
            debug!(
                "Some centering translations are missing. Consider reducing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        debug!("Order of point group: {}", operations.len());
        Ok(Self {
            operations,
            permutations,
        })
    }

    fn check_closure(operations: &Operations, lattice: &Lattice, symprec: f64) -> bool {
        let mut translations_map = FxHashMap::default();
        for operation in operations.iter() {
            translations_map.insert(operation.rotation, operation.translation);
        }
        for ops1 in operations.iter() {
            for ops2 in operations.iter() {
                let ops12 = ops1.clone() * ops2.clone();
                let Some(translation) = translations_map.get(&ops12.rotation) else {
                    return false;
                };
                let diff = (translation - ops12.translation).map(|e| e - e.round());
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
        let neighbor_search = PeriodicNeighborSearch::new(&primitive_magnetic_cell.cell, symprec);
        let mut magnetic_operations = vec![];
        let mut permutations = vec![];
        for operation in candidate_operations.iter() {
            let new_positions = primitive_magnetic_cell
                .cell
                .positions
                .iter()
                .map(|pos| operation.rotation.map(|e| e as f64) * pos + operation.translation)
                .collect::<Vec<_>>();
            if let Some(permutation) = solve_correspondence(
                &neighbor_search,
                &primitive_magnetic_cell.cell,
                &new_positions,
            ) {
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
            debug!(
                "Some centering translations are missing. Consider reducing symprec and angle_tolerance."
            );
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
                let Some(expected_translation) =
                    translations_map.get(&(mops12.operation.rotation, mops12.time_reversal))
                else {
                    return false;
                };
                let diff =
                    (expected_translation - mops12.operation.translation).map(|e| e - e.round());
                if lattice.cartesian_coords(&diff).norm() > symprec {
                    return false;
                }
            }
        }
        true
    }
}

pub fn operations_in_cell(prim_cell: &PrimitiveCell, prim_operations: &Operations) -> Operations {
    operations_in_cell_from_linear_and_translations(
        &prim_cell.linear,
        &prim_cell.translations,
        prim_operations,
    )
}

pub(crate) fn operations_in_cell_from_linear_and_translations(
    linear: &Linear,
    translations: &[Translation],
    prim_operations: &Operations,
) -> Operations {
    let input_operations =
        Transformation::from_linear(*linear).transform_operations(prim_operations);
    let mut operations = vec![];
    for t1 in translations {
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
pub(crate) fn search_bravais_group(
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
    // For a Minkowski-reduced lattice, it is sufficient to search coeffs in
    // `[-1, 0, 1]` because the third column is already in `[-1, 0, 1]^3`.
    // This bounded integer search is exact under that precondition.
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
        debug!(
            "Found automorphisms for the lattice do not form a group. Consider reducing symprec and angle_tolerance."
        );
        return Err(MoyoError::TooLargeToleranceError);
    }

    let complemented_rotations = traverse(&rotations);
    if complemented_rotations.len() != rotations.len() {
        debug!(
            "Found automorphisms for the lattice do not form a group. Consider reducing symprec and angle_tolerance."
        );
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
    let len_u = (basis.column(col1).norm() + b1.norm()) / 2.0;
    let len_v = (basis.column(col2).norm() + b2.norm()) / 2.0;
    is_angle_within_tolerance(
        theta_new - theta_org,
        len_u,
        len_v,
        symprec,
        angle_tolerance,
    )
}

#[cfg(test)]
mod tests {
    use std::vec;

    use nalgebra::{Matrix3, Vector3, matrix};
    use test_log::test;

    use super::{PrimitiveMagneticSymmetrySearch, PrimitiveSymmetrySearch, search_bravais_group};
    use crate::base::{
        AngleTolerance, Cell, Collinear, Lattice, MagneticCell, NonCollinear, Operation,
        RotationMagneticMomentAction,
    };

    #[test]
    fn test_closure_checks_rotations_and_translations() {
        let lattice = Lattice::new(Matrix3::identity());
        let quarter_turn = matrix![0, -1, 0; 1, 0, 0; 0, 0, 1];
        let rotations = [
            Matrix3::identity(),
            quarter_turn,
            quarter_turn * quarter_turn,
            quarter_turn.transpose(),
        ];
        let operations = rotations
            .iter()
            .map(|&rotation| Operation::new(rotation, Vector3::zeros()))
            .collect::<Vec<_>>();
        assert!(PrimitiveSymmetrySearch::check_closure(
            &operations,
            &lattice,
            1e-5
        ));

        let missing_rotation = operations
            .iter()
            .filter(|op| op.rotation != quarter_turn.transpose())
            .cloned()
            .collect();
        assert!(!PrimitiveSymmetrySearch::check_closure(
            &missing_rotation,
            &lattice,
            1e-5
        ));

        let wrong_translation = operations
            .into_iter()
            .map(|mut op| {
                if op.rotation == quarter_turn {
                    op.translation = Vector3::new(0.1, 0.0, 0.0);
                }
                op
            })
            .collect();
        assert!(!PrimitiveSymmetrySearch::check_closure(
            &wrong_translation,
            &lattice,
            1e-5
        ));
    }

    #[test]
    fn test_closure_preserves_site_permutations() {
        let mut positions = vec![Vector3::zeros()];
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut position = Vector3::zeros();
                position[axis] = 0.2 * sign;
                positions.push(position);
            }
        }
        let cell = Cell::new(
            Lattice::new(Matrix3::identity()),
            positions,
            vec![0, 1, 1, 1, 1, 1, 1],
        );
        let search = PrimitiveSymmetrySearch::new(&cell, 1e-5, AngleTolerance::default()).unwrap();
        assert_eq!(search.operations.len(), 48);
        assert_eq!(search.permutations.len(), 48);
        assert_eq!(
            search
                .operations
                .iter()
                .map(|op| op.rotation)
                .collect::<std::collections::HashSet<_>>()
                .len(),
            48
        );
        for (operation, permutation) in search.operations.iter().zip(&search.permutations) {
            for (site, position) in cell.positions.iter().enumerate() {
                assert_relative_eq!(
                    operation.rotation.map(f64::from) * position + operation.translation,
                    cell.positions[permutation.apply(site)],
                    epsilon = 1e-12
                );
            }
            for (rhs, rhs_permutation) in search.operations.iter().zip(&search.permutations) {
                let rotation = operation.rotation * rhs.rotation;
                let (_, product) = search
                    .operations
                    .iter()
                    .zip(&search.permutations)
                    .find(|(candidate, _)| candidate.rotation == rotation)
                    .unwrap();
                for site in 0..cell.num_atoms() {
                    assert_eq!(
                        product.apply(site),
                        permutation.apply(rhs_permutation.apply(site))
                    );
                }
            }
        }
    }

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
                search_bravais_group(&lattice, symprec, AngleTolerance::default()).unwrap();
            assert_eq!(rotations.len(), 24);
        }

        {
            let lattice = Lattice::new(matrix![
                0.5, 0.0, 0.5;
                0.5, 0.5, 0.0;
                0.0, 0.5, 0.5;
            ]);
            let rotations =
                search_bravais_group(&lattice, symprec, AngleTolerance::default()).unwrap();
            assert_eq!(rotations.len(), 48);
        }
    }

    #[test]
    fn test_primitive_magnetic_symmetry_search() {
        let symprec = 1e-4;
        let angle_tolerance = AngleTolerance::default();
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
