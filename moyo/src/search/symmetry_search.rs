use itertools::iproduct;
use std::collections::HashSet;

use nalgebra::{Matrix3, Vector3};

use super::solve::{
    pivot_site_indices, solve_correspondence, symmetrize_translation_from_permutation,
};
use crate::base::{
    AngleTolerance, Cell, Lattice, MoyoError, Operations, Permutation, Position, Rotation, EPS,
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
            return Err(MoyoError::TooSmallSymprecError);
        }

        // Search symmetry operations
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
                    solve_correspondence(primitive_cell, &new_positions, rough_symprec)
                {
                    symmetries_tmp.push((*rotation, translation, permutation));
                    // If a translation part is found, it should be unique (up to lattice translations)
                    break;
                }
            }
        }
        assert!(!symmetries_tmp.is_empty());

        // Purify symmetry operations by permutations
        let mut rotations = vec![];
        let mut translations = vec![];
        let mut permutations = vec![];
        for (rotation, rough_translation, permutation) in symmetries_tmp.iter() {
            let (translation, distance) = symmetrize_translation_from_permutation(
                primitive_cell,
                permutation,
                rotation,
                rough_translation,
            );
            if distance < symprec {
                rotations.push(*rotation);
                translations.push(translation);
                permutations.push(permutation.clone());
            }
        }
        if rotations.is_empty() {
            return Err(MoyoError::PrimitiveSymmetrySearchError);
        }

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

    if rotations.is_empty() {
        return Err(MoyoError::BravaisGroupSearchError);
    }

    // Check to reproduce rotation operations by group multiplication
    let mut rotation_set = HashSet::new();
    for rotation in rotations.iter() {
        rotation_set.insert(*rotation);
    }
    for (r1, r2) in iproduct!(rotation_set.iter(), rotation_set.iter()) {
        let r12 = r1 * r2;
        if !rotation_set.contains(&r12) {
            return Err(MoyoError::BravaisGroupSearchError);
        }
    }

    Ok(rotations)
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
