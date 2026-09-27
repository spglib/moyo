# Unimodular transformation types

This is the prerequisite for the handedness changes in [#473](https://github.com/spglib/moyo/issues/473).
Both types represent passive affine changes of crystallographic coordinates:
`B = A P`, `x_new = P^-1 (x_old - p)`.

| Rust type | Linear part | Intended use |
| --- | --- | --- |
| `UnimodularTransformation` | Integer `P`, determinant +1 or -1 | General changes of primitive basis, including handedness changes; Euclidean normalizer elements |
| `ProperUnimodularTransformation` | Integer `P`, determinant +1 | Orientation-preserving identification conjugators and setting corrections |

The existing general type keeps its name. The proper type wraps it privately and
shares its coordinate-transformation methods through immutable `Deref`. Neither
type exposes mutable access to its linear part or cached inverse. Rust callers
read the matrix through `linear()` instead of the former public `linear` field.
Constructors follow the existing panic-on-invalid-input convention; determinant
validation and inversion use integer arithmetic, and the inverse must fit `i32`.

Conversion from proper to general is infallible (`From`). Conversion from general
to proper is checked (`TryFrom`), rejecting orientation reversal. Inversion and
composition preserve the proper type when all operands are proper. Mixed
composition returns the general type, including when the particular product has
determinant +1. Origin shifts follow the existing affine composition convention.
Serialization retains the existing transformation fields for both types.

The bulk integral-normalizer search, type-IV magnetic conjugator search, and
conventional-cell correction selector already admit only determinant +1; their
return types become proper. Space-group result witnesses and lattice-aware
paths keep the general type because they must support a change of handedness.
Python continues to expose its existing frozen general transformation class;
the Rust boundary explicitly widens proper results for that binding.

This prerequisite leaves reduction outputs and layer-group behavior unchanged.
The subsequent reduction change will use general unimodular transformations to
make 3D output bases right-handed. General integer `Transformation` (possibly
changing the lattice index) is a separate contract, addressed in that change.

Acceptance checks cover constructor rejection, checked conversion, affine
composition order, inverse round trips, mixed-type products, immutable linear
access, and unchanged serialization. Existing identification, standardization,
normalizer, magnetic, and layer regressions must pass before this prerequisite
is committed.
