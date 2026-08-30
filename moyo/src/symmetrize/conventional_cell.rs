use itertools::Itertools;
use log::debug;
use nalgebra::Matrix3;
use once_cell::sync::Lazy;

use crate::base::{EPS, Lattice, Operations, UnimodularLinear, UnimodularTransformation};
use crate::data::Centering;
use crate::identify::match_origin_shift;

/// The six axis permutations as proper (det = +1) signed-permutation matrices. The
/// three odd permutations (transpositions) are made right-handed by negating their
/// first axis, since [`Transformation`] requires a positive determinant.
///
/// One proper representative per permutation is sufficient: for every orthorhombic
/// space group, whether a permutation keeps the setting does not depend on the choice
/// of axis signs (the point-group rotations are diagonal, so conjugation by a sign
/// flip leaves them unchanged, and a sign flip only changes the sign of a translation
/// part, which is absorbed by the origin shift). The sign-flipped axes are re-oriented
/// away later by `symmetrize_lattice`.
pub(super) static AXIS_PERMUTATIONS3: Lazy<Vec<UnimodularLinear>> = Lazy::new(|| {
    vec![
        // Even permutations (already det = +1).
        Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1), // abc (identity)
        Matrix3::new(0, 1, 0, 0, 0, 1, 1, 0, 0), // cab
        Matrix3::new(0, 0, 1, 1, 0, 0, 0, 1, 0), // bca
        // Odd permutations, first axis negated to restore det = +1.
        Matrix3::new(0, 1, 0, -1, 0, 0, 0, 0, 1), // b(-a)c  (a <-> b)
        Matrix3::new(-1, 0, 0, 0, 0, 1, 0, 1, 0), // (-a)cb  (b <-> c)
        Matrix3::new(0, 0, 1, 0, 1, 0, -1, 0, 0), // c b(-a) (a <-> c)
    ]
});

/// Candidate unimodular corrections with entries in `[-1, 1]`.
///
/// This bounded search is used as a best-effort search space for monoclinic
/// conventional-cell standardization. Exhaustiveness is not claimed here.
pub(super) static UNIMODULAR3_RANGE1: Lazy<Vec<UnimodularLinear>> = Lazy::new(|| {
    (0..9)
        .map(|_| -1_i32..=1_i32)
        .multi_cartesian_product()
        .filter_map(|v| {
            let mat = Matrix3::new(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
            let det = mat.map(|e| e as f64).determinant().round() as i32;
            if det == 1 { Some(mat) } else { None }
        })
        .collect()
});

/// Ranking key for monoclinic conventional cells: closest to orthogonal first, then
/// (following the ITA convention) the non-acute monoclinic angle among the supplements,
/// which have the same `|cos|`.
pub(super) fn monoclinic_rank_key(lattice: &Lattice) -> Vec<f64> {
    let cos_angles = lattice.lattice_constant()[3..]
        .iter()
        .map(|angle_deg| angle_deg.to_radians().cos())
        .collect::<Vec<_>>();
    let skewness = cos_angles.iter().map(|cos| cos.abs()).sum::<f64>();
    let signed_cos_sum = cos_angles.iter().sum::<f64>();
    vec![skewness, signed_cos_sum]
}

/// Ranking key for orthorhombic conventional cells: lexicographically shortest
/// `(a, b, c)`, i.e. `a <= b <= c` as far as the setting allows.
pub(super) fn orthorhombic_rank_key(lattice: &Lattice) -> Vec<f64> {
    let lc = lattice.lattice_constant();
    vec![lc[0], lc[1], lc[2]]
}

/// Choose, among `candidates` (changes of basis of the conventional cell), the one that
/// keeps the Hall setting and minimizes `rank_key` (compared lexicographically with
/// `EPS` ties; the first candidate wins a tie, so candidate order makes the selection
/// stable).
///
/// A candidate keeps the setting when it preserves the centering translations as a set
/// and conjugates the space group onto itself up to an origin shift, i.e. when it
/// belongs to the affine normalizer. The origin shift is solved by `match_origin_shift`
/// exactly as in space-group identification, so settings such as `Imma`, whose
/// `a <-> b` swap needs an origin shift, are handled as well.
///
/// The correction is returned in the primitive standardized basis (`Q corr Q^-1` with
/// `Q = centering.linear()`, which is integral because `corr` maps the centering
/// lattice onto itself) so that the caller can fold it into the primitive
/// transformation and keep the fixed primitive-to-conventional matrix `Q`.
pub(super) fn select_conventional_correction<F>(
    conv_lattice: &Lattice,
    centering: Centering,
    prim_std_operations: &Operations,
    db_prim_generators: &Operations,
    candidates: &[UnimodularLinear],
    rank_key: F,
    epsilon: f64,
) -> UnimodularTransformation
where
    F: Fn(&Lattice) -> Vec<f64>,
{
    let q = centering.linear().map(|e| e as f64);
    let q_inv = q.try_inverse().unwrap();

    let mut best: Option<(Vec<f64>, UnimodularTransformation)> = None;
    for corr in candidates {
        let corr_f64 = corr.map(|e| e as f64);
        if corr_f64.determinant().round() as i32 != 1 {
            continue;
        }
        if !preserves_centering(corr, centering, epsilon) {
            continue;
        }

        let prim_corr_f64 = q * corr_f64 * q_inv;
        let prim_corr = prim_corr_f64.map(|e| e.round() as i32);
        if (prim_corr_f64 - prim_corr.map(|e| e as f64)).abs().max() > epsilon {
            continue;
        }
        let Some(origin_shift) =
            match_origin_shift(prim_std_operations, &prim_corr, db_prim_generators, epsilon)
        else {
            continue;
        };

        let key =
            rank_key(&UnimodularTransformation::from_linear(*corr).transform_lattice(conv_lattice));
        let is_better = match &best {
            Some((best_key, _)) => lexicographic_less(&key, best_key),
            None => true,
        };
        if is_better {
            best = Some((key, UnimodularTransformation::new(prim_corr, origin_shift)));
        }
    }

    match best {
        Some((_, correction)) => correction,
        None => {
            debug!("No admissible correction of the conventional cell; keep the identified one");
            UnimodularTransformation::from_linear(UnimodularLinear::identity())
        }
    }
}

/// Whether `linear` maps the set of centering translations onto itself modulo the
/// integer lattice. Checking the set (rather than each translation individually) is
/// required for `F`-centering, whose three face translations are permuted by an axis
/// permutation.
fn preserves_centering(linear: &UnimodularLinear, centering: Centering, epsilon: f64) -> bool {
    let lattice_points = centering.lattice_points();
    let linear_f64 = linear.map(|e| e as f64);
    lattice_points.iter().all(|translation| {
        let mapped = linear_f64 * translation;
        lattice_points.iter().any(|other| {
            let mut diff = mapped - other;
            diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
            diff.iter().all(|e| e.abs() <= epsilon)
        })
    })
}

/// Lexicographic `lhs < rhs` where entries within `EPS` are considered equal.
fn lexicographic_less(lhs: &[f64], rhs: &[f64]) -> bool {
    for (l, r) in lhs.iter().zip(rhs.iter()) {
        if (l - r).abs() < EPS {
            continue;
        }
        return l < r;
    }
    false
}
