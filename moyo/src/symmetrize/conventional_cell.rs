use itertools::iproduct;
use log::warn;
use nalgebra::{Matrix2, Matrix3, Vector2};
use once_cell::sync::Lazy;

use crate::base::{EPS, Lattice, Operations, UnimodularLinear, UnimodularTransformation};
use crate::data::Centering;
use crate::identify::match_origin_shift;
use crate::math::minkowski_reduce_2d;

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

/// Candidate changes of basis for a monoclinic conventional cell.
///
/// The two basis vectors perpendicular to the unique axis are Lagrange-reduced in their
/// plane and completed to the Delaunay triple `v1, v2, -(v1 + v2)`, whose members are
/// pairwise non-acute. The candidates are all ordered pairs of distinct triple vectors,
/// with every sign choice for the pair and for the unique axis.
///
/// Pairs from the triple are sufficient for the following reasons.
///
/// 1. Any two of the three vectors form a basis of the plane lattice, since
///    `v3 = -(v1 + v2)`.
/// 2. Whether a change of basis keeps the Hall setting depends only on the integer
///    coefficients modulo 2: the centering and glide translations are half-integer
///    vectors, so `preserves_centering` and `match_origin_shift` only see the class of
///    each new axis in `L / 2L`. The triple vectors realize all three non-zero classes,
///    hence every admissible basis has an admissible "parity twin" in the triple, and
///    the search never falls back to the identity.
/// 3. For an obtuse superbase in two dimensions, `+-v1, +-v2, +-v3` are the
///    Voronoi-relevant vectors, i.e. the shortest vectors of their classes modulo `2L`.
///    Among the equally admissible bases of a parity signature, the triple pair thus
///    has the shortest in-plane axes, which is the convention (ITA, Parthe-Gelato,
///    spglib): the shortest axes compatible with the setting first, and only then the
///    angle. Bases with longer axes in the same class, e.g. `(v1, v3 + 2 v1)`, may have
///    an angle closer to 90 deg but are deliberately excluded.
pub(super) fn monoclinic_candidate_corrections(conv_lattice: &Lattice) -> Vec<UnimodularLinear> {
    let unique_axis = monoclinic_unique_axis(conv_lattice);
    let i = (unique_axis + 1) % 3;
    let j = (unique_axis + 2) % 3;
    let basis = conv_lattice.basis; // column-wise
    let vi = basis.column(i).into_owned();
    let vj = basis.column(j).into_owned();

    // Coordinates of the in-plane basis in an orthonormal frame of the plane.
    let e1 = vi.normalize();
    let e2 = (vj - vj.dot(&e1) * e1).normalize();
    let basis_2d = Matrix2::new(vi.dot(&e1), vj.dot(&e1), vi.dot(&e2), vj.dot(&e2));
    let (reduced_2d, trans_2d) = minkowski_reduce_2d(&basis_2d);

    // Integer coefficients of the Delaunay triple with respect to `(vi, vj)`.
    let c1 = Vector2::new(trans_2d[(0, 0)], trans_2d[(1, 0)]);
    let mut c2 = Vector2::new(trans_2d[(0, 1)], trans_2d[(1, 1)]);
    if reduced_2d.column(0).dot(&reduced_2d.column(1)) > 0.0 {
        c2 = -c2;
    }
    let c3 = -(c1 + c2);
    let triple = [c1, c2, c3];

    let mut candidates = vec![];
    for (x, y) in iproduct!(0..3, 0..3).filter(|(x, y)| x != y) {
        for (sign_x, sign_y, sign_unique) in iproduct!([1, -1], [1, -1], [1, -1]) {
            let mut corr = UnimodularLinear::zeros();
            corr[(i, i)] = sign_x * triple[x][0];
            corr[(j, i)] = sign_x * triple[x][1];
            corr[(i, j)] = sign_y * triple[y][0];
            corr[(j, j)] = sign_y * triple[y][1];
            corr[(unique_axis, unique_axis)] = sign_unique;
            candidates.push(corr);
        }
    }
    candidates
}

/// Index of the unique axis of a monoclinic conventional cell, read off the metric: the
/// unique axis is perpendicular to the other two, so it is the axis opposite the angle
/// farthest from 90 deg (`alpha` <-> `a`, `beta` <-> `b`, `gamma` <-> `c`). This works for
/// every unique-axis setting of the Hall-symbol database. If all three angles are 90 deg
/// within noise, the in-plane basis is orthogonal, hence already reduced, and any choice
/// leaves the identity among the candidates.
fn monoclinic_unique_axis(conv_lattice: &Lattice) -> usize {
    let lc = conv_lattice.lattice_constant();
    (0..3)
        .max_by(|&i, &j| {
            (lc[3 + i] - 90.0)
                .abs()
                .partial_cmp(&(lc[3 + j] - 90.0).abs())
                .unwrap()
        })
        .unwrap()
}

/// Ranking key for monoclinic conventional cells: closest to orthogonal first, then
/// (following the ITA convention) the non-acute monoclinic angle among the supplements
/// -- they have the same `|cos|` -- and finally the lexicographically shortest
/// `(a, b, c)`. The last key orders `a <= c` when the setting allows the `a <-> c` swap
/// (P2, P2_1, Pm, P2/m, P2_1/m), the convention of E. Parthe and L. M. Gelato, "The
/// best unit cell for monoclinic structures compatible with b axis unique and cell
/// choice 1", Acta Cryst. A39, 169-173 (1983), which spglib follows as well.
pub(super) fn monoclinic_rank_key(lattice: &Lattice) -> Vec<f64> {
    let lc = lattice.lattice_constant();
    let cos_angles = lc[3..]
        .iter()
        .map(|angle_deg| angle_deg.to_radians().cos())
        .collect::<Vec<_>>();
    let skewness = cos_angles.iter().map(|cos| cos.abs()).sum::<f64>();
    let signed_cos_sum = cos_angles.iter().sum::<f64>();
    vec![skewness, signed_cos_sum, lc[0], lc[1], lc[2]]
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
            // Unreachable when the identity is among the candidates; a candidate set
            // that misses it (or a bug in the admissibility test) would otherwise pass
            // an uncorrected cell through unnoticed.
            warn!("No admissible correction of the conventional cell; keep the identified one");
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
