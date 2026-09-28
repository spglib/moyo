use approx::assert_relative_eq;
use nalgebra::Matrix3;

pub fn assert_lattice_standardization(
    basis: &Matrix3<f64>,
    std_basis: &Matrix3<f64>,
    prim_std_basis: &Matrix3<f64>,
    std_linear: &Matrix3<f64>,
    prim_std_linear: &Matrix3<f64>,
    std_rotation_matrix: &Matrix3<f64>,
) {
    // prim_std_linear should be an inverse of an integer matrix
    let prim_std_linear_inv = prim_std_linear.map(|e| e as f64).try_inverse().unwrap();
    assert_relative_eq!(
        prim_std_linear_inv,
        prim_std_linear_inv.map(|e| e.round()),
        epsilon = 1e-8
    );

    // Refinement adds a symmetric positive stretch before the rigid rotation.
    let stretch =
        std_rotation_matrix.transpose() * std_basis * (basis * std_linear).try_inverse().unwrap();
    assert_relative_eq!(stretch, stretch.transpose(), epsilon = 1e-12);
    assert!(stretch.symmetric_eigen().eigenvalues.min() > 0.0);
    assert_relative_eq!(
        std_rotation_matrix.transpose() * std_rotation_matrix,
        Matrix3::identity(),
        epsilon = 1e-12
    );
    // The same stretch and rotation apply to the primitive cell.
    assert_relative_eq!(
        std_rotation_matrix * stretch * basis * prim_std_linear,
        *prim_std_basis,
        epsilon = 1e-10
    );
}
