use itertools::iproduct;

use nalgebra::{Matrix3, Vector3};

use crate::base::cell::Cell;
use crate::base::error::MoyoError;
use crate::base::lattice::Lattice;
use crate::base::operation::Rotation;
use crate::base::tolerance::{AngleTolerance, EPS};

pub struct SymmetrySearchResult {}

/// Return coset representatives of the space group w.r.t. its translation subgroup.
/// Assume `primitive_cell` is a primitive cell and its basis vectors are Minkowski reduced.
/// Possible replacements for spglib/src/spacegroup.h::spa_search_spacegroup
pub fn search_symmetry_operations(
    primitive_cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<SymmetrySearchResult, MoyoError> {
    let bravais_group = search_bravais_group(&primitive_cell.lattice, symprec, angle_tolerance)?;

    // Search symmetry operations

    unimplemented!()
}

/// Relevant to spglib.c/symmetry.c::get_lattice_symmetry
fn search_bravais_group(
    minkowski_lattice: &Lattice,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<Vec<Rotation>, MoyoError> {
    assert!(minkowski_lattice.is_minkowski_reduced());

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
                candidate_lattice_points[i].push(coeffs.clone());
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

    if rotations.len() == 0 {
        return Err(MoyoError::BravaisGroupSearchError);
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
    let cos_org = basis.column(col1).dot(&basis.column(col2))
        / (basis.column(col1).norm() * basis.column(col2).norm());
    let cos_new = b1.dot(b2) / (b1.norm() * b2.norm());
    let cos_dtheta =
        cos_org * cos_new + (1.0 - cos_org.powi(2)).sqrt() * (1.0 - cos_new.powi(2)).sqrt();

    match angle_tolerance {
        AngleTolerance::Radian(angle_tolerance) => cos_dtheta.acos().abs() < angle_tolerance,
        AngleTolerance::Default => {
            // Eq.(7) of https://arxiv.org/pdf/1808.01590.pdf
            let sin_dtheta2 = 1.0 - cos_dtheta.powi(2);
            let length_ave2 = (basis.column(col1).norm() + b1.norm())
                * (basis.column(col2).norm() + b2.norm())
                / 4.0;
            (sin_dtheta2 * length_ave2).sqrt() < symprec
        }
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, Vector3};

    use super::{search_bravais_group, search_symmetry_operations};
    use crate::base::cell::Cell;
    use crate::base::lattice::Lattice;
    use crate::base::tolerance::AngleTolerance;

    #[test]
    fn test_search_symmetry_operations() {
        let symprec = 1e-4;

        // Primitive fcc
        let primitive_cell = Cell::new(
            Lattice::new(matrix![
                0.0, 0.5, 0.5;
                0.5, 0.0, 0.5;
                0.5, 0.5, 0.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.0)],
            vec![0],
        );
        let result =
            search_symmetry_operations(&primitive_cell, symprec, AngleTolerance::Default).unwrap();
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
                1.0, -0.5, 0.0;
                0.0, f64::sqrt(3.0) / 2.0, 0.0;
                0.0, 0.0, 1.0;
            ]);
            let rotations =
                search_bravais_group(&lattice, symprec, AngleTolerance::Default).unwrap();
            assert_eq!(rotations.len(), 24);
        }
    }
}
