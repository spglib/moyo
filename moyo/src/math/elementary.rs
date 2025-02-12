use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, OMatrix};

/// Return elementary matrix swapping the `col1`th and `col2`th columns
pub fn swapping_column_matrix<D: Dim>(dim: D, col1: usize, col2: usize) -> OMatrix<i32, D, D>
where
    DefaultAllocator: Allocator<D, D>,
{
    let mut trans_mat = OMatrix::zeros_generic(dim, dim);
    for i in 0..dim.value() {
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
pub fn adding_column_matrix<D: Dim>(dim: D, col1: usize, col2: usize, k: i32) -> OMatrix<i32, D, D>
where
    DefaultAllocator: Allocator<D, D>,
{
    let mut trans_mat = OMatrix::identity_generic(dim, dim);
    for i in 0..dim.value() {
        if i == col1 {
            trans_mat[(col1, col2)] = k;
        }
    }
    trans_mat
}

/// Return elementary matrix changing sign of the `col`th column
pub fn changing_column_sign_matrix<D: Dim>(dim: D, col: usize) -> OMatrix<i32, D, D>
where
    DefaultAllocator: Allocator<D, D>,
{
    let mut trans_mat = OMatrix::identity_generic(dim, dim);
    trans_mat[(col, col)] = -1;
    trans_mat
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, U3};

    use super::{adding_column_matrix, changing_column_sign_matrix, swapping_column_matrix};

    #[test]
    fn test_swapping_column_matrix() {
        let actual = swapping_column_matrix(U3, 0, 1);
        let expect = matrix![
            0, 1, 0;
            1, 0, 0;
            0, 0, 1;
        ];
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_adding_column_matrix() {
        let actual = adding_column_matrix(U3, 0, 2, -1);
        let expect = matrix![
            1, 0, -1;
            0, 1, 0;
            0, 0, 1;
        ];
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_changing_column_sign_matrix() {
        let actual = changing_column_sign_matrix(U3, 0);
        let expect = matrix![
            -1, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ];
        assert_eq!(actual, expect);
    }
}
