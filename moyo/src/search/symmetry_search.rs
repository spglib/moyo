use itertools::iproduct;
use std::collections::{HashMap, HashSet, VecDeque};

use log::debug;
use nalgebra::{Matrix3, Vector3};

use super::solve::{
    pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
    PeriodicKdTree,
};
use crate::base::{
    traverse, AngleTolerance, Cell, Lattice, MoyoError, Operations, Permutation, Position,
    Rotation, Translation, EPS,
};

#[derive(Debug)]
pub struct PrimitiveSymmetrySearch {
    /// Operations in the given primitive cell
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
    pub bravais_group: Vec<Rotation>,
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
            for dst in pivot_site_indices.iter() {
                // Try to overlap the `src`-th site to the `dst`-th site
                let translation = primitive_cell.positions[*dst]
                    - rotation.map(|e| e as f64) * primitive_cell.positions[src];
                let new_positions: Vec<Position> = primitive_cell
                    .positions
                    .iter()
                    .map(|pos| rotation.map(|e| e as f64) * pos + translation)
                    .collect();

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
                operations_and_permutations.push((rotation, translation, permutation.clone()));
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
        let mut rotations = vec![];
        let mut translations = vec![];
        let mut permutations = vec![];
        queue.push_back((
            Rotation::identity(),
            Translation::zeros(),
            Permutation::identity(primitive_cell.num_atoms()),
        ));
        while !queue.is_empty() {
            let (rotation_lhs, translation_lhs, permutation_lhs) = queue.pop_front().unwrap();
            if visited.contains(&rotation_lhs) {
                continue;
            }
            visited.insert(rotation_lhs);
            rotations.push(rotation_lhs);
            translations.push(translation_lhs);
            permutations.push(permutation_lhs.clone());

            for (&rotation_rhs, translation_rhs, permutation_rhs) in
                operations_and_permutations.iter()
            {
                let new_rotation = rotation_lhs * rotation_rhs;
                let new_translation = (rotation_lhs.map(|e| e as f64) * translation_rhs
                    + translation_lhs)
                    .map(|e| e - e.round());
                let new_permutation = permutation_lhs.clone() * permutation_rhs.clone();
                queue.push_back((new_rotation, new_translation, new_permutation));
            }
        }
        if rotations.len() != operations_and_permutations.len() {
            debug!("Found operations do not form a group. Consider reducing symprec and angle_tolerance.");
            return Err(MoyoError::TooLargeToleranceError);
        }

        // Check closure
        let mut translations_map = HashMap::new();
        for (rotation, translation) in rotations.iter().zip(translations.iter()) {
            translations_map.insert(rotation.clone(), translation.clone());
        }
        let mut closed = true;
        for (r1, t1) in rotations.iter().zip(translations.iter()) {
            if !closed {
                break;
            }
            for (r2, t2) in rotations.iter().zip(translations.iter()) {
                // (r1, t1) * (r2, t2) = (r1 * r2, r1 * t2 + t1)
                let r = r1 * r2;
                let t = r1.map(|e| e as f64) * t2 + t1;
                let diff = (translations_map[&r] - t).map(|e| e - e.round());
                if primitive_cell.lattice.cartesian_coords(&diff).norm() > rough_symprec {
                    closed = false;
                    break;
                }
            }
        }
        if !closed {
            debug!("Some centering translations are missing. Consider reducing symprec and angle_tolerance.");
            return Err(MoyoError::TooLargeToleranceError);
        }

        debug!("Order of point group: {}", rotations.len());
        Ok(Self {
            operations: Operations::new(rotations, translations),
            permutations,
            bravais_group,
        })
    }
}

/// Relevant to spglib.c/symmetry.c::get_lattice_symmetry
fn search_bravais_group(
    minkowski_lattice: &Lattice,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<Vec<Rotation>, MoyoError> {
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

    // Recover rotations by group multiplication
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
    use nalgebra::matrix;

    use super::search_bravais_group;
    use crate::base::{AngleTolerance, Lattice};

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
}
