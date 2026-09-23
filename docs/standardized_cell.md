# `StandardizedCell` specification

Status: proposed contract for the three-dimensional `StandardizedCell` in
[`moyo/src/symmetrize/standardize.rs`](../moyo/src/symmetrize/standardize.rs).
This specifies the intended behavior for Point 2 of
[issue #451](https://github.com/spglib/moyo/issues/451); the current implementation
does not yet satisfy the lattice guarantee below.

## Core guarantee

On success, both `cell` and `prim_cell` represent a crystal invariant under the
supplied, identified symmetry, up to floating-point numerical error. This
requires exact lattice symmetry and exact site symmetry, independently of
`rotate_basis`.

Write lattice vectors as columns of a nonsingular matrix $A$, fractional
positions as $x_i$, and species as $Z_i$. For every target operation $(W,w)$ in
the basis of the returned cell, there must be a bijective site permutation
$\pi_{(W,w)}$ such that

$$
W^T A^T A W = A^T A,
\qquad
W x_i+w-x_{\pi_{(W,w)}(i)}\in\mathbb Z^3,
\qquad
Z_i=Z_{\pi_{(W,w)}(i)}.
$$

The equalities are evaluated within numerical error as defined below. The
permutations must realize the target group action. The output may have
additional symmetry; the guarantee is that every target operation is a
symmetry of the output.

## Inputs and the target symmetry

`prim_cell` is the input primitive cell. `prim_operations` and
`prim_permutations` describe its detected symmetry and site correspondence.
`space_group` fixes the identified group, Hall setting, and transformation to
that setting. Standardization keeps this identification fixed.

The exact target operations are the refined operations of `space_group`,
expressed in the respective standardized primitive and conventional bases.
The supplied permutations are aligned with these operations while retaining
input primitive-site indices. An operation's fractional matrix and translation
must always be interpreted in the basis and origin where they are defined.

Detected translations can contain search-level errors and need not close to
machine precision. Those translations guide the identification and site
correspondence; the invariant above is checked against the refined group.
An inconsistent group identification or permutation action is a
`StandardizationError`.

## Coordinate changes and field meanings

Let $A_0$ and $x_i^0$ denote the input primitive basis and positions. Let
$(P_p,p)$ be `prim_transformation`, and let $C$ be the fixed
primitive-to-conventional matrix for the selected centering. Then

$$
P_c=P_p C,
\qquad
\texttt{transformation}=(P_c,p),
\qquad
A_p^{(0)}=A_0P_p,
\qquad
A_c^{(0)}=A_0P_c.
$$

The primitive transformation is unimodular and orientation-preserving, and
$\det C>0$. Its fractional coordinate change before position refinement is

$$
y_i=P_p^{-1}(x_i^0-p).
$$

These transformations describe the crystallographic change of basis and origin.
Lattice and position refinement are additional operations. In particular,
`rotation_matrix * input_basis * transformation.linear` generally differs from
the returned basis, and $y_i$ generally differs from the refined position.

| Field | Required meaning |
| --- | --- |
| `prim_cell` | Refined primitive cell, with the input primitive-site count, order, and species preserved. |
| `prim_transformation` | Crystallographic basis and origin change $(P_p,p)$ before refinement. |
| `cell` | Conventional expansion of the same refined primitive crystal. |
| `transformation` | Crystallographic basis and origin change $(P_pC,p)$ before refinement. |
| `rotation_matrix` | Proper Cartesian rotation applied after lattice metric refinement; identity when `rotate_basis = false`. |
| `site_mapping` | For every conventional site, the corresponding index in `prim_cell`. |
| `wyckoffs` | One Wyckoff assignment per conventional site, consistent with the returned positions and target Hall setting. |

## Lattice refinement and Cartesian orientation

Form the conventional metric before refinement,
$G_c^{(0)}=(A_c^{(0)})^T A_c^{(0)}$. With $H_c$ the distinct rotation parts of the
target conventional operations, use the group average

$$
\bar G_c=\frac{1}{|H_c|}\sum_{W\in H_c}W^T G_c^{(0)}W.
$$

This positive-definite metric satisfies $W^T\bar G_cW=\bar G_c$ for every
$W\in H_c$. It also preserves a metric that already satisfies those equations.
The returned conventional basis must have metric $\bar G_c$ for either value of
`rotate_basis`.

To define the orientation unambiguously, choose the canonical upper-triangular
basis $B$ with $B^TB=\bar G_c$ and the same handedness as $A_c^{(0)}$. Concretely,
if $\bar G_c=LL^T$ is a Cholesky factorization and
$h=\operatorname{sign}(\det A_c^{(0)})$, take

$$
B=\operatorname{diag}(1,1,h)L^T.
$$

The handedness adjustment acts on Cartesian rows, so it preserves the metric.
Separate the deformation to this basis by its right polar decomposition:

$$
F=B(A_c^{(0)})^{-1}=QU,
\qquad
Q\in SO(3),
\qquad
U=U^T>0.
$$

| `rotate_basis` | Returned conventional basis $A_c$ | `rotation_matrix` |
| --- | --- | --- |
| `true` | $B=QUA_c^{(0)}$ | $Q$ |
| `false` | $UA_c^{(0)}=Q^TB$ | $I$ |

Thus `false` applies a pure symmetric strain correction without a superposed
rigid rotation. Both modes have the same metric, fractional positions, site
ordering, Wyckoff assignments, and crystallographic transformations. The
rotation must remain orthogonal; it must not absorb the strain correction.
For an already invariant input metric, $U=I$ within numerical error.

The full Cartesian deformation is recoverable as
$A_c(A_c^{(0)})^{-1}$; no additional stored matrix is required by this contract.
The polar decomposition specifies the intended meaning of rotation and replaces
the current QR-based separation.

## Position refinement and primitive/conventional consistency

Refine the transformed primitive positions $y_i$ using the aligned site
permutations and the target primitive operations. For a consistent choice of
nearby periodic images, project onto the affine site-symmetry constraints,
minimizing the sum of squared Cartesian displacements in the refined primitive
metric. Already invariant positions are fixed points of this projection.

Keep the primitive site indices and species. Periodic wrapping is a choice of
representation; all position comparisons are modulo lattice translations.

Construct both outputs from one refined crystal. Their bases and positions must
satisfy

$$
A_c=A_pC,
\qquad
C x_j^c-x_{\texttt{site\_mapping}[j]}^p\in\mathbb Z^3.
$$

Consequently $N_c=(\det C)N_p$, and every primitive site has exactly $\det C$
conventional representatives of the same species. The two outputs share an
origin and the same Cartesian strain correction and rotation. Derive the
primitive basis from $A_p=A_cC^{-1}$ so the two lattices cannot acquire
independent refinement errors.

Wyckoff assignment uses the refined crystal. It must agree with the site
orbits, multiplicities, and coordinates of that crystal in the chosen setting.

## Numerical error and success criteria

`symprec` and the existing `epsilon` argument control recognition and
correspondence of the approximate input. They do not define the allowed
symmetry residual of a successful output. Increasing a search tolerance must
not permit proportionately larger residuals in the refined cell.

Measure lattice error using

$$
r_G=\max_{W\in H}\frac{\|W^TGW-G\|_F}{\|G\|_F},
$$

and site error using the supplied, aligned permutations,

$$
r_x=\frac{1}{\ell}\max_{(W,w),i}\min_{n\in\mathbb Z^3}
\|A(Wx_i+w-x_{\pi_{(W,w)}(i)}-n)\|_2,
\qquad
\ell=|\det A|^{1/3}.
$$

Numerical error bounds must account for floating-point precision, lattice
conditioning, and operation magnitudes. For well-conditioned `f64` regression
fixtures, require both residuals to be at most $10^{-12}$. Any larger bound for
an ill-conditioned case needs an explicit numerical-error justification,
independent of the symmetry-search tolerance.

Successful output must satisfy the lattice, site, field, and centering
invariants together. If refinement cannot produce such a cell, return
`StandardizationError`. A positive-definite metric and finite coordinates are
required throughout.

## Acceptance cases and integration requirements

Validation must exercise the returned `StandardizedCell`, including both
primitive and conventional cells, for both rotation settings.

- Perturb lengths, angles, and atomic positions in fixtures spanning the
  lattice systems. Check target-operation residuals directly.
- Include the sheared NaCl example from issue #451 with the identified
  monoclinic group, and a cubic fixture with unequal input axis lengths.
- Include centered cells and nonsymmorphic operations to check the centering
  relation, site mapping, and fractional translations.
- Include globally rotated and left-handed inputs. Check proper rotations,
  preserved handedness, and the equivalence of the two orientation modes.
- Include already exact cells, general sites, and special sites. Refinement
  must be a fixed point when reapplied with the same target operations and
  correspondences, up to numerical error and periodic representatives.
- Vary the search tolerance while retaining the same target symmetry and
  correspondence; output residuals must remain at numerical precision.
- Reject inconsistent site-permutation actions. A successful test must check
  the supplied target operations, rather than relying on symmetry
  re-identification at the original search tolerance.

The shared dataset tests currently assert the pure-rotation lattice relation.
Replace that assertion for distorted inputs with the refined metric and
primitive/conventional invariants above. Keep the rigid relation as an exact
input control. Update the Rust and binding documentation of the public dataset
transformation fields when implementing this contract.

`StandardizedMagneticCell` consumes this struct. Its Cartesian symmetry
operations must be formed using the refined lattice in the same fractional
basis as those operations. Apply `rotation_matrix` to Cartesian magnetic
moments; the strain matrix is not a rotation of magnetic moments. The magnetic
standardization tests must continue to validate the resulting moment symmetry.

This document specifies the 3D cell contract. Layer-cell standardization has its
own lattice constraints and remains a separate specification.
