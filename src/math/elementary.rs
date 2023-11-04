use nalgebra::SMatrix;

/// Return elementary matrix swapping the `col1`th and `col2`th columns
pub fn swapping_column_matrix<const DIM: usize>(
    col1: usize,
    col2: usize,
) -> SMatrix<i32, DIM, DIM> {
    let mut trans_mat = SMatrix::<i32, DIM, DIM>::zeros();
    for i in 0..DIM {
        if i == col1 {
            trans_mat[(col1, col2)] = 1;
        } else if i == col2 {
            trans_mat[(col2, col1)] = 1;
        } else {
            trans_mat[(i, i)] = 1
        }
    }
    trans_mat
}

/// Return elementary matrix adding the `k`-multiplied `col1`th column into the `col2`th column
#[allow(dead_code)]
pub fn adding_column_matrix<const DIM: usize>(
    col1: usize,
    col2: usize,
    k: i32,
) -> SMatrix<i32, DIM, DIM> {
    let mut trans_mat = SMatrix::<i32, DIM, DIM>::identity();
    for i in 0..DIM {
        if i == col1 {
            trans_mat[(col1, col2)] = k;
        }
    }
    trans_mat
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::{adding_column_matrix, swapping_column_matrix};

    #[test]
    fn test_swapping_column_matrix() {
        let actual = swapping_column_matrix::<3>(0, 1);
        let expect = matrix![
            0, 1, 0;
            1, 0, 0;
            0, 0, 1;
        ];
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_adding_column_matrix() {
        let actual = adding_column_matrix::<3>(0, 2, -1);
        let expect = matrix![
            1, 0, -1;
            0, 1, 0;
            0, 0, 1;
        ];
        assert_eq!(actual, expect);
    }
}
