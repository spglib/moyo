use nalgebra::{matrix, Matrix3};

use super::cycle_checker::CycleChecker;

const EPS: f64 = 1e-8;

pub fn niggli_reduce(basis: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<i32>) {
    let mut reduced_basis = *basis;
    let mut trans_mat = Matrix3::<i32>::identity();

    let mut cc = CycleChecker::new();
    let mut step = 1;

    while step <= 8 {
        let params = NiggliParameters::new(&reduced_basis);
        let branch = match step {
            1 => step1(&params, &mut trans_mat),
            2 => step2(&params, &mut trans_mat),
            3 => step3(&params, &mut trans_mat),
            4 => step4(&params, &mut trans_mat),
            5 => step5(&params, &mut trans_mat),
            6 => step6(&params, &mut trans_mat),
            7 => step7(&params, &mut trans_mat),
            8 => step8(&params, &mut trans_mat),
            _ => unreachable!(),
        };
        reduced_basis = basis * trans_mat.map(|e| e as f64);

        if branch && (step == 2 || step == 5 || step == 6 || step == 7 || step == 8) {
            step = 1;
        } else {
            step += 1;
        }

        // If the transformation matrix is already visited in Step 1, terminate
        if (step == 1) && !cc.insert(&trans_mat) {
            break;
        }
    }

    // Preserve parity
    if trans_mat.map(|e| e as f64).determinant() < 0. {
        reduced_basis *= -1.;
        trans_mat *= -1;
    }

    (reduced_basis, trans_mat)
}

#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[allow(clippy::collapsible_if)]
pub fn is_niggli_reduced(basis: &Matrix3<f64>) -> bool {
    let params = NiggliParameters::new(basis);

    // Common conditions:
    //   A <= B <= C, A >= |eta|, A >= |zeta|, B >= |xi|
    if !(params.b - params.a >= -EPS) {
        return false;
    }
    if !(params.c - params.b >= -EPS) {
        return false;
    }
    if !(params.a - params.eta.abs() >= -EPS) {
        return false;
    }
    if !(params.b - params.xi.abs() >= -EPS) {
        return false;
    }

    if params.sign_xi * params.sign_eta * params.sign_zeta > 0 {
        // Type-I cell
        // Main conditions:
        //    xi > 0, eta > 0, zeta > 0
        if !(params.xi > EPS) {
            return false;
        }
        if !(params.eta > EPS) {
            return false;
        }
        if !(params.zeta > EPS) {
            return false;
        }

        // Special conditions:
        //    If A = B, xi <= eta
        //    If B = C, eta <= zeta
        //    If B = |xi|, zeta <= 2 * eta
        //    If A = |eta|, zeta <= 2 * xi
        //    If A = |zeta|, eta <= 2 * xi
        if (params.a - params.b).abs() < EPS {
            if !(params.eta - params.xi >= -EPS) {
                return false;
            }
        }
        if (params.b - params.c).abs() < EPS {
            if !(params.zeta - params.eta >= -EPS) {
                return false;
            }
        }
        if (params.b - params.xi.abs()).abs() < EPS {
            if !(2.0 * params.eta - params.zeta >= -EPS) {
                return false;
            }
        }
        if (params.a - params.eta.abs()).abs() < EPS {
            if !(2.0 * params.xi - params.zeta >= -EPS) {
                return false;
            }
        }
        if (params.a - params.zeta.abs()).abs() < EPS {
            if !(2.0 * params.xi - params.eta >= -EPS) {
                return false;
            }
        }
    } else {
        // Type-II cell
        // Main conditions:
        //    xi <= 0, eta <= 0, zeta <= 0
        if !(params.xi <= EPS) {
            return false;
        }
        if !(params.eta <= EPS) {
            return false;
        }
        if !(params.zeta <= EPS) {
            return false;
        }

        // Special conditions:
        //    If A = B, |xi| <= |eta|
        //    If B = C, |eta| <= |zeta|
        //    If B = |xi|, zeta = 0
        //    If A = |eta|, zeta = 0
        //    If A = |zeta|, eta = 0
        //    If xi + eta + zeta = A + B, A <= |eta| + |zeta|
        if (params.a - params.b).abs() < EPS {
            if !(params.eta.abs() - params.xi.abs() >= -EPS) {
                return false;
            }
        }
        if (params.b - params.c).abs() < EPS {
            if !(params.zeta.abs() - params.eta.abs() >= -EPS) {
                return false;
            }
        }
        if (params.b - params.xi.abs()).abs() < EPS {
            if !(params.zeta.abs() < EPS) {
                return false;
            }
        }
        if (params.a - params.eta.abs()).abs() < EPS {
            if !(params.zeta.abs() < EPS) {
                return false;
            }
        }
        if (params.a - params.zeta.abs()).abs() < EPS {
            if !(params.eta.abs() < EPS) {
                return false;
            }
        }
        if (params.xi + params.eta + params.zeta - params.a - params.b).abs() < EPS {
            if !(params.eta.abs() + params.zeta.abs() - params.a >= -EPS) {
                return false;
            }
        }
    }
    true
}

/// If A > B or (A = B, |xi| > |eta|), swap a and b
fn step1(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.a - params.b > EPS)
        || ((params.a - params.b).abs() < EPS && params.xi.abs() > params.eta.abs())
    {
        let trans_mat_tmp = matrix![
            0, -1, 0;
            -1, 0, 0;
            0, 0, -1;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// If B > C or (B = C, |eta| > |zeta|), swap b and c
fn step2(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.b - params.c > EPS)
        || ((params.b - params.c).abs() < EPS && params.eta.abs() > params.zeta.abs())
    {
        let trans_mat_tmp = matrix![
            -1, 0, 0;
            0, 0, -1;
            0, -1, 0;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// Adjust axis directions for type I cell
fn step3(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if params.sign_eta * params.sign_xi * params.sign_zeta > 0 {
        // Here params.sign_* != 0
        let trans_mat_tmp = matrix![
            if params.sign_xi == -1 { -1 } else { 1 }, 0, 0;
            0, if params.sign_eta == -1 { -1 } else { 1 }, 0;
            0, 0, if params.sign_zeta == -1 { -1 } else { 1 };
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// Adjust axis directions for type II cell
fn step4(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if params.sign_xi == -1 && params.sign_eta == -1 && params.sign_zeta == -1 {
        return false;
    }

    if params.sign_xi * params.sign_eta * params.sign_zeta <= 0 {
        let mut i = 1;
        let mut j = 1;
        let mut k = 1;
        let mut p = -1;

        if params.sign_xi == 1 {
            i = -1;
        } else if params.sign_xi == 0 {
            p = 0;
        }
        if params.sign_eta == 1 {
            j = -1;
        } else if params.sign_eta == 0 {
            p = 1;
        }
        if params.sign_zeta == 1 {
            k = -1;
        } else if params.sign_zeta == 0 {
            p = 2;
        }
        if i * j * k == -1 {
            match p {
                0 => i = -1,
                1 => j = -1,
                2 => k = -1,
                _ => unreachable!(),
            }
        }

        let trans_mat_tmp = matrix![
            i, 0, 0;
            0, j, 0;
            0, 0, k;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// If |xi| > B or (xi = B, 2 eta < zeta) or (xi = -B, zeta < 0)
fn step5(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.xi.abs() - params.b > EPS)
        || ((params.xi - params.b).abs() < EPS && params.zeta - 2.0 * params.eta > EPS)
        || ((params.xi + params.b).abs() < EPS && -params.zeta > EPS)
    {
        let trans_mat_tmp = matrix![
            1, 0, 0;
            0, 1, -params.sign_xi;
            0, 0, 1;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// If |eta| > A or (eta = A, 2 * xi < zeta) or (eta = -A, zeta < 0)
fn step6(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.eta.abs() - params.a > EPS)
        || ((params.eta - params.a).abs() < EPS && params.zeta - 2.0 * params.xi > EPS)
        || ((params.eta + params.a).abs() < EPS && -params.zeta > EPS)
    {
        let trans_mat_tmp = matrix![
            1, 0, -params.sign_eta;
            0, 1, 0;
            0, 0, 1;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// If |zeta| > A or (zeta = A, 2 * xi < eta) or (zeta = -A, eta < 0)
fn step7(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.zeta.abs() - params.a > EPS)
        || ((params.zeta - params.a).abs() < EPS && params.eta - 2.0 * params.xi > EPS)
        || ((params.zeta + params.a).abs() < EPS && -params.eta > EPS)
    {
        let trans_mat_tmp = matrix![
            1, -params.sign_zeta, 0;
            0, 1, 0;
            0, 0, 1;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

/// If xi + eta + zeta + A + B < 0 or (xi + eta + zeta + A + B = 0, 2 * (A + eta) + zeta > 0)
fn step8(params: &NiggliParameters, trans_mat: &mut Matrix3<i32>) -> bool {
    if (params.xi + params.eta + params.zeta + params.a + params.b < -EPS)
        || ((params.xi + params.eta + params.zeta + params.a + params.b).abs() < EPS
            && 2.0 * (params.a + params.eta) + params.zeta > EPS)
    {
        let trans_mat_tmp = matrix![
            1, 0, 1;
            0, 1, 1;
            0, 0, 1;
        ];
        *trans_mat *= trans_mat_tmp;
        true
    } else {
        false
    }
}

#[derive(Debug)]
struct NiggliParameters {
    /// a * a
    a: f64,
    /// b * b
    b: f64,
    /// c * c
    c: f64,
    /// 2 * b * c * cos(alpha)
    xi: f64,
    /// 2 * c * a * cos(beta)
    eta: f64,
    /// 2 * a * b * cos(gamma)
    zeta: f64,
    sign_xi: i32,
    sign_eta: i32,
    sign_zeta: i32,
}

impl NiggliParameters {
    pub fn new(basis: &Matrix3<f64>) -> Self {
        let metric_tensor = basis.transpose() * basis;
        let xi = 2.0 * metric_tensor[(1, 2)];
        let eta = 2.0 * metric_tensor[(2, 0)];
        let zeta = 2.0 * metric_tensor[(0, 1)];

        Self {
            a: metric_tensor[(0, 0)],
            b: metric_tensor[(1, 1)],
            c: metric_tensor[(2, 2)],
            xi,
            eta,
            zeta,
            sign_xi: sign(xi),
            sign_eta: sign(eta),
            sign_zeta: sign(zeta),
        }
    }
}

fn sign(x: f64) -> i32 {
    if x > EPS {
        1
    } else if x < -EPS {
        -1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Matrix3;
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::{is_niggli_reduced, matrix, niggli_reduce, sign, NiggliParameters, EPS};

    #[test]
    fn test_sign() {
        assert_eq!(sign(2.0 * EPS), 1);
        assert_eq!(sign(-2.0 * EPS), -1);
        assert_eq!(sign(0.5 * EPS), 0);
        assert_eq!(sign(-0.5 * EPS), 0);
    }

    #[test]
    fn test_niggli_reduction_small() {
        // Example in Acta Cryst. (1976). A32, 297
        let g = matrix![
            9.0, -22.0 / 2.0, -4.0 / 2.0;
            -22.0 / 2.0, 27.0, -5.0 / 2.0;
            -4.0 / 2.0, -5.0 / 2.0, 4.0;
        ];
        let basis = g.cholesky().unwrap().l().transpose();

        let (reduced_basis, trans_mat) = niggli_reduce(&basis);
        assert_relative_eq!(basis * trans_mat.map(|e| e as f64), reduced_basis);
        assert_eq!(trans_mat.map(|e| e as f64).determinant().round() as i32, 1);

        let params = NiggliParameters::new(&reduced_basis);
        assert_relative_eq!(params.a, 4.0, epsilon = 1e-6);
        assert_relative_eq!(params.b, 9.0, epsilon = 1e-6);
        assert_relative_eq!(params.c, 9.0, epsilon = 1e-6);
        assert_relative_eq!(params.xi, 9.0, epsilon = 1e-6);
        assert_relative_eq!(params.eta, 3.0, epsilon = 1e-6);
        assert_relative_eq!(params.zeta, 4.0, epsilon = 1e-6);

        assert!(is_niggli_reduced(&reduced_basis));
    }

    #[test]
    fn test_niggli_reduction_small2() {
        let basis = matrix![
            -101.0,
            95.0,
            126.0;
            7.0,
            4.0,
            46.0;
            -127.0,
            73.0,
            5.0;
        ]
        .transpose();
        let (reduced_basis, _) = niggli_reduce(&basis);
        assert!(is_niggli_reduced(&reduced_basis));
    }

    #[test]
    fn test_niggli_reduction_small3() {
        let basis = matrix![
            17.0,
            -105.0,
            -117.0;
            105.0,
            -108.0,
            -113.0;
            85.0,
            2.0,
            2.0;
        ]
        .transpose();
        let (reduced_basis, _) = niggli_reduce(&basis);
        assert!(is_niggli_reduced(&reduced_basis));
    }

    #[test]
    fn test_niggli_reduction_random() {
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
            let (reduced_basis, trans_mat) = niggli_reduce(&basis);

            assert!(is_niggli_reduced(&reduced_basis));
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
            let (reduced_basis, trans_mat) = niggli_reduce(&basis);
            assert!(is_niggli_reduced(&reduced_basis));
            assert_relative_eq!(
                basis * trans_mat.map(|e| e as f64),
                reduced_basis,
                epsilon = 1e-4
            );
        }
    }
}
