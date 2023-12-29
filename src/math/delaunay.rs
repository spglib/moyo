use nalgebra::{Matrix3, Vector3, U3};

use crate::base::transformation::TransformationMatrix;

use super::elementary::{adding_column_matrix, changing_column_sign_matrix};

const EPS: f64 = 1e-8;

/// basis is column-wise
pub fn delaunay_reduce(basis: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<i32>) {
    let mut reduced_basis = *basis;
    let mut trans_mat = Matrix3::<i32>::identity();

    loop {
        let mut superbase = vec![];
        let mut sum_vec = Vector3::<f64>::zeros();
        for base in basis.column_iter() {
            superbase.push(base.clone_owned());
            sum_vec += base;
        }
        superbase.push(-sum_vec);

        let mut update = false;
        for i in 0..3 {
            if update {
                break;
            }
            for j in i + 1..4 {
                if superbase[i].dot(&superbase[j]) < EPS {
                    let mut trans_mat_tmp = TransformationMatrix::identity();
                    for k in 0..3 {
                        if (k == i) || (k == j) {
                            continue;
                        }
                        trans_mat_tmp *= adding_column_matrix(U3, i, k, 1);
                    }
                    trans_mat_tmp *= changing_column_sign_matrix(U3, i);

                    reduced_basis *= trans_mat_tmp.map(|e| e as f64);
                    trans_mat *= trans_mat_tmp;
                    update = true;
                    break;
                }
            }
        }

        if update {
            break;
        }
    }

    (reduced_basis, trans_mat)
}

#[cfg(test)]
mod tests {
    use nalgebra::{vector, Matrix3};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::delaunay_reduce;

    #[test]
    fn test_delaunay_small() {
        {
            // https://github.com/spglib/spglib/blob/develop/test/functional/c/test_delaunay.cpp
            let basis = Matrix3::from_columns(&[
                vector![
                    -2.2204639179669590,
                    -4.4409278359339179,
                    179.8575773553236843
                ],
                vector![1.2819854407640749, 0.0, 103.8408207018900669],
                vector![10.5158083946732219, 0.0, 883.3279051525505565],
            ]);
            let (reduced_basis, trans_mat) = delaunay_reduce(&basis);
            assert_relative_eq!(basis * trans_mat.map(|e| e as f64), reduced_basis);
        }
    }

    #[test]
    fn test_delaunay_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let basis = Matrix3::<f64>::new(
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
            );
            let (reduced_basis, trans_mat) = delaunay_reduce(&basis);
            assert_relative_eq!(basis * trans_mat.map(|e| e as f64), reduced_basis);
        }

        for _ in 0..256 {
            let basis = Matrix3::<f64>::new(
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
            );
            let (reduced_basis, trans_mat) = delaunay_reduce(&basis);
            assert_relative_eq!(
                basis * trans_mat.map(|e| e as f64),
                reduced_basis,
                epsilon = 1e-4
            );
        }
    }
}
