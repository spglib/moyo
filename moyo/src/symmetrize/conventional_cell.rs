use itertools::Itertools;
use nalgebra::{Matrix3, Vector3};
use once_cell::sync::Lazy;
use std::cmp::Ordering;

use crate::base::{EPS, Lattice, Linear, Operations, Transformation, UnimodularTransformation};
use crate::data::Centering;

/// Candidate unimodular corrections with entries in `[-1, 1]`.
///
/// This bounded search is used as a best-effort search space for monoclinic
/// conventional-cell standardization. Exhaustiveness is not claimed here.
static UNIMODULAR3_RANGE1: Lazy<Vec<UnimodularTransformation>> = Lazy::new(|| {
    (0..9)
        .map(|_| -1_i32..=1_i32)
        .multi_cartesian_product()
        .filter_map(|v| {
            let mat = Matrix3::new(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
            let det = mat.map(|e| e as f64).determinant().round() as i32;
            if det == 1 {
                Some(UnimodularTransformation::from_linear(mat))
            } else {
                None
            }
        })
        .collect()
});

/// The six axis permutations as proper (det = +1) signed-permutation matrices. The
/// three odd permutations (transpositions) are made right-handed by negating their
/// first axis, since [`Transformation`] requires a positive determinant.
///
/// One proper representative per permutation is sufficient: for every orthorhombic
/// space group, whether a permutation preserves the operation set and centering does
/// not depend on the choice of axis signs (the point-group rotations are diagonal, so
/// conjugation by a sign flip leaves them unchanged, and centering half-translations
/// absorb the sign of any non-symmorphic translation). The sign-flipped axes are
/// re-oriented away later by `symmetrize_lattice`.
static AXIS_PERMUTATIONS3: Lazy<Vec<UnimodularTransformation>> = Lazy::new(|| {
    [
        // Even permutations (already det = +1).
        Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1), // abc (identity)
        Matrix3::new(0, 1, 0, 0, 0, 1, 1, 0, 0), // cab
        Matrix3::new(0, 0, 1, 1, 0, 0, 0, 1, 0), // bca
        // Odd permutations, first axis negated to restore det = +1.
        Matrix3::new(0, 1, 0, -1, 0, 0, 0, 0, 1), // b(-a)c  (a <-> b)
        Matrix3::new(-1, 0, 0, 0, 0, 1, 0, 1, 0), // (-a)cb  (b <-> c)
        Matrix3::new(0, 0, 1, 0, 1, 0, -1, 0, 0), // c b(-a) (a <-> c)
    ]
    .into_iter()
    .map(UnimodularTransformation::from_linear)
    .collect()
});

/// Standardize monoclinic cell by choosing reduced basis vectors for perpendicular plane to the unique axis while keeping matrix representations of `conv_std_generators`.
pub(super) fn standardize_monoclinic_conv_cell(
    prim_lattice: &Lattice,
    transformation_to_prim_std: &UnimodularTransformation,
    centering: Centering,
    conv_std_generators: &Operations,
    epsilon: f64,
) -> Linear {
    let mut candidate_conv_transformations = vec![];
    // Search candidate corrections in the bounded `[-1, 1]` window and choose
    // the least skewed valid one. This is a best-effort bounded search.
    for trans_corr in UNIMODULAR3_RANGE1.iter() {
        // trans_corr should keep centering translations
        if centering.lattice_points().iter().any(|translation| {
            let mut diff = trans_corr.linear.map(|e| e as f64) * translation - translation;
            diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
            diff.iter().any(|e| e.abs() > epsilon)
        }) {
            continue;
        }

        // trans_corr should keep the matrix representations of `conv_std_generators`.
        if conv_std_generators.iter().any(|ops| {
            let corr_ops = trans_corr.transform_operation(ops);
            if corr_ops.rotation != ops.rotation {
                true
            } else {
                let mut diff = corr_ops.translation - ops.translation;
                diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
                diff.iter().any(|e| e.abs() > epsilon)
            }
        }) {
            continue;
        }

        let refined_conv_trans = Transformation::new(
            transformation_to_prim_std.linear * centering.linear() * trans_corr.linear,
            transformation_to_prim_std.origin_shift,
        );
        let refined_conv_lattice = refined_conv_trans.transform_lattice(prim_lattice);
        let cos_angles = refined_conv_lattice.lattice_constant()[3..]
            .iter()
            .map(|angle_deg| angle_deg.to_radians().cos())
            .collect::<Vec<_>>();
        // `skewness` gets smaller as the basis vectors are closer to orthorhombic.
        let skewness = cos_angles.iter().map(|cos| cos.abs()).sum::<f64>();
        // The monoclinic angle and its supplement give the same `skewness`
        // (`|cos|` is symmetric about 90 deg), so the unimodular correction that
        // flips the angle to its supplement is an equally valid candidate. Among
        // such ties we follow the ITA convention and prefer a non-acute angle,
        // i.e. the smaller (more negative) signed cosine sum.
        let signed_cos_sum = cos_angles.iter().sum::<f64>();
        candidate_conv_transformations.push((
            skewness,
            signed_cos_sum,
            centering.linear() * trans_corr.linear,
        ));
    }

    candidate_conv_transformations
        .into_iter()
        .min_by(|(skewness_lhs, cos_lhs, _), (skewness_rhs, cos_rhs, _)| {
            if (skewness_lhs - skewness_rhs).abs() < EPS {
                cos_lhs.partial_cmp(cos_rhs).unwrap()
            } else {
                skewness_lhs.partial_cmp(skewness_rhs).unwrap()
            }
        })
        .unwrap()
        .2
}

/// Standardize an orthorhombic conventional cell by choosing the axis permutation
/// that orders the basis-vector lengths as `a <= b <= c` as much as possible, while
/// keeping the matrix representations of the symmetry operations (as a set) and the
/// centering.
///
/// Unlike the monoclinic case, an axis permutation relabels generators among
/// themselves, so the *set* of `conv_std_operations` -- not each individual generator
/// -- must be preserved. Permutations that would change the operation set or move a
/// centered face (e.g. for `C`-centering) are rejected, which reproduces seekpath's
/// "number of distinct projections" without a hardcoded per-space-group table.
pub(super) fn standardize_orthorhombic_conv_cell(
    prim_lattice: &Lattice,
    transformation_to_prim_std: &UnimodularTransformation,
    centering: Centering,
    conv_std_operations: &Operations,
    epsilon: f64,
) -> Linear {
    let mut candidate_conv_transformations = vec![];
    // Consider the six axis permutations and keep those that preserve the centering
    // and the operation set; pick the one giving the smallest lengths.
    for trans_perm in AXIS_PERMUTATIONS3.iter() {
        // trans_perm should keep centering translations.
        if centering.lattice_points().iter().any(|translation| {
            let mut diff = trans_perm.linear.map(|e| e as f64) * translation - translation;
            diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
            diff.iter().any(|e| e.abs() > epsilon)
        }) {
            continue;
        }

        // trans_perm should keep the *set* of conv_std_operations: each conjugated
        // operation must coincide with some operation in the set (rotation exactly,
        // translation modulo the centering lattice). `conv_std_operations` holds one
        // coset representative per rotation, so translations are compared modulo the
        // centering translations, not only modulo 1.
        if conv_std_operations.iter().any(|ops| {
            let perm_ops = trans_perm.transform_operation(ops);
            !conv_std_operations.iter().any(|other| {
                perm_ops.rotation == other.rotation
                    && translations_equal_mod_centering(
                        &perm_ops.translation,
                        &other.translation,
                        centering,
                        epsilon,
                    )
            })
        }) {
            continue;
        }

        let refined_conv_trans = Transformation::new(
            transformation_to_prim_std.linear * centering.linear() * trans_perm.linear,
            transformation_to_prim_std.origin_shift,
        );
        let lc = refined_conv_trans
            .transform_lattice(prim_lattice)
            .lattice_constant();
        let lengths = [lc[0], lc[1], lc[2]];
        candidate_conv_transformations.push((lengths, centering.linear() * trans_perm.linear));
    }

    // Choose the lexicographically smallest (a, b, c). `min_by` keeps the first
    // element among EPS-ties, so iteration order makes the selection stable.
    candidate_conv_transformations
        .into_iter()
        .min_by(|(lengths_lhs, _), (lengths_rhs, _)| {
            lengths_lhs
                .iter()
                .zip(lengths_rhs.iter())
                .find_map(|(lhs, rhs)| {
                    if (lhs - rhs).abs() < EPS {
                        None
                    } else {
                        Some(lhs.partial_cmp(rhs).unwrap())
                    }
                })
                .unwrap_or(Ordering::Equal)
        })
        .unwrap()
        .1
}

/// Whether two fractional translations are equal modulo the centering lattice (which
/// includes the integer lattice). Used to compare symmetry operations of a
/// conventional centered cell, where a given coset may be represented by translations
/// differing by a centering vector.
fn translations_equal_mod_centering(
    lhs: &Vector3<f64>,
    rhs: &Vector3<f64>,
    centering: Centering,
    epsilon: f64,
) -> bool {
    let mut diff = lhs - rhs;
    diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
    centering.lattice_points().iter().any(|lattice_point| {
        let mut residual = diff - lattice_point;
        residual -= residual.map(|e| e.round());
        residual.iter().all(|e| e.abs() <= epsilon)
    })
}
