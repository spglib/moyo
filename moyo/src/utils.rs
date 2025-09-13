use nalgebra::{Matrix3, Vector3};

pub fn to_3x3_slice<T: Copy>(mat: &Matrix3<T>) -> [[T; 3]; 3] {
    // Be careful that nalgebra stores matrices in column-major order
    [
        [mat[(0, 0)], mat[(0, 1)], mat[(0, 2)]],
        [mat[(1, 0)], mat[(1, 1)], mat[(1, 2)]],
        [mat[(2, 0)], mat[(2, 1)], mat[(2, 2)]],
    ]
}

pub fn to_matrix3<T: Copy + nalgebra::Scalar>(mat: &[[T; 3]; 3]) -> Matrix3<T> {
    // Be careful that nalgebra stores matrices in column-major order
    Matrix3::from_columns(&[
        nalgebra::Vector3::new(mat[0][0], mat[1][0], mat[2][0]),
        nalgebra::Vector3::new(mat[0][1], mat[1][1], mat[2][1]),
        nalgebra::Vector3::new(mat[0][2], mat[1][2], mat[2][2]),
    ])
}

pub fn to_3_slice<T: Copy>(vec: &Vector3<T>) -> [T; 3] {
    [vec[0], vec[1], vec[2]]
}

pub fn to_vector3<T: Copy + nalgebra::Scalar>(vec: &[T; 3]) -> Vector3<T> {
    Vector3::new(vec[0], vec[1], vec[2])
}
