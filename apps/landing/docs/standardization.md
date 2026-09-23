# Moyo conventions of standardized cell

This document specifies the returned-cell contract for the three-dimensional standardization pipeline.
The core Rust result is `StandardizedCell`; its conventional and primitive cells become `MoyoDataset.std_cell` and `MoyoDataset.prim_std_cell`.
The same cell conventions apply to the language bindings.

!!! note "Implementation status"
    The exact-symmetry and Cartesian-orientation requirements below are the target contract.
    The current implementation refines positions but uses only the rotation returned by lattice symmetrization, retaining the input metric after the selected change of basis.
    Applying the refined lattice to both output cells, the polar-decomposition convention, and roundoff-level validation remain to be implemented.

**Standardization** comprises selecting a conventional coordinate system, symmetrizing the lattice and positions under fixed target operations, and constructing the returned cells and metadata.
`ConventionalCoordinateSystem` owns the basis and origin selection; `StandardizedCell` owns the complete result.

## `setting` option in `MoyoDataset::new`

We refer to the criteria to choose a representative of each space-group type as a setting.
After symmetry operations of an input crystal structure are determined, moyo transforms the space group into a representative of its space-group type defined by the chosen `setting`.
Here, the transformation comprises both a change of basis and an origin shift.

Moyo supports three settings;

- `Setting.Standard`: the so called ["ITA setting"](http://webbdcrista2.ehu.es/cgi-bin/cryst/programs/nph-def-ita_settings), which is one of the conventional descriptions for each space-group type used in the International Tables for Crystallography Volume A (2016). The ITA setting chooses unique axis b setting, cell choice 1 for monoclinic space groups, hexagonal axes for rhombohedral groups, and origin choice 2 for centrosymmetric groups.
- `Setting.Spglib`: the setting used in [spglib](https://spglib.readthedocs.io/en/stable/index.html), which chooses the smallest one in the serial numbers for Hall symbols described in [Prof. Seto's page](https://yseto.net/en/sg/sg1).
  Be sure that this serial number, so called "Hall number", would not be a standard crystallographic definition [^setting].
- `Setting.HallNumber`: allows users to specify a "Hall number" (an integer from 1 to 530) directly.

Moyo chooses `Setting.Standard` as the default setting, which is different from spglib's default `Setting.Spglib`.
This change of the default behavior affects in centrosymmetric groups: moyo chooses origin choice 2 by default, while spglib chooses origin choice 1 by default.

## Returned-cell contract

### Target symmetry and numerical accuracy

The target operations are those of the selected Hall symbol, expressed in the conventional and primitive coordinate systems chosen by `ConventionalCoordinateSystem`.
The detected operations and their atom permutations establish the correspondence with that target symmetry.
The returned cells must satisfy the target operations, including their refined translations.

For either returned cell, let $\mathbf{A}_*$ contain its basis vectors as columns and let $\mathbf{G}_* = \mathbf{A}_*^\mathsf{T}\mathbf{A}_*$.
Every target operation $(\mathbf{W}, \mathbf{w})$ must satisfy

$$
\mathbf{W}^\mathsf{T}\mathbf{G}_*\mathbf{W} = \mathbf{G}_*.
$$

It must also induce a species-preserving permutation $\pi$ of the returned sites such that, for every site $i$,

$$
\mathbf{W}\mathbf{x}_i + \mathbf{w} - \mathbf{x}_{\pi(i)} \in \mathbb{Z}^3.
$$

These are exact identities in real arithmetic.
In floating-point arithmetic, their residuals must be at roundoff scale, independently of the tolerances used to recognize symmetry in the input.
For example, dimensionless residuals can be measured as

$$
r_G =
\frac{
  \lVert \mathbf{W}^\mathsf{T}\mathbf{G}_*\mathbf{W} - \mathbf{G}_* \rVert_F
}{
  (1 + \lVert \mathbf{W} \rVert_F^2)\lVert \mathbf{G}_* \rVert_F
},
$$

$$
r_{x,i} =
\frac{
  \min_{\mathbf{n}\in\mathbb{Z}^3}
  \lVert \mathbf{A}_*(\mathbf{W}\mathbf{x}_i+\mathbf{w}-\mathbf{x}_{\pi(i)}-\mathbf{n}) \rVert_2
}{
  \lVert \mathbf{A}_* \rVert_F
  (1+\lVert\mathbf{W}\rVert_F\lVert\mathbf{x}_i\rVert_2+\lVert\mathbf{w}\rVert_2+\lVert\mathbf{x}_{\pi(i)}\rVert_2)
}.
$$

Here $\lVert\cdot\rVert_F$ is the Frobenius norm and $\lVert\cdot\rVert_2$ is the Euclidean vector norm.
The required scale is $O(\epsilon_{\mathrm{mach}})$ for these residuals, allowing for numerical conditioning and accumulated floating-point operations.
`symprec` and the fractional-coordinate tolerance `epsilon` are not the accuracy promised for the returned cells.
The concrete roundoff bounds must be established and tested with the numerical implementation.
If the target constraints cannot be met within those bounds, construction must return `MoyoError::StandardizationError`.

Both values of `rotate_basis` must satisfy this contract.
Changing Cartesian orientation must leave the refined metric, fractional positions, site correspondences, and Wyckoff assignments invariant.
If the input already satisfies the target symmetry, symmetrization preserves its metric and fractional positions up to roundoff; the selected coordinate change and Cartesian orientation still apply.

### Coordinate transformations and refinement

With column-wise input basis $\mathbf{A}$ and fractional coordinates $\mathbf{x}$, a selected transformation $(\mathbf{P},\mathbf{p})$ means

$$
\mathbf{A}' = \mathbf{A}\mathbf{P},\qquad
\mathbf{x}' = \mathbf{P}^{-1}(\mathbf{x}-\mathbf{p}).
$$

The origin shift $\mathbf{p}$ is expressed in the input basis.
These equations describe the coordinate change before lattice and position refinement.
The matrices and origin shifts returned by standardization record that coordinate change; they do not encode the subsequent corrections to the metric or fractional positions.

`StandardizedCell` takes a primitive input cell.
Its `transformation` and `prim_transformation` fields describe the selected conventional and primitive coordinate systems relative to that input.
`MoyoDataset` composes these transformations with the input-to-primitive transformation, so its corresponding fields refer to the original dataset input.
On the dataset, $(\mathbf{P}_{\mathrm{std}},\mathbf{p}_{\mathrm{std}})$ means `(std_linear, std_origin_shift)`, and $(\mathbf{P}_{\mathrm{prim}},\mathbf{p}_{\mathrm{prim}})$ means `(prim_std_linear, prim_std_origin_shift)`.
The returned Cartesian rotation is `std_rotation_matrix`.

| `StandardizedCell` field | Meaning |
| --- | --- |
| `cell` | Conventional output in the selected Hall setting; exposed as `MoyoDataset.std_cell`. |
| `prim_cell` | Primitive description of the same refined crystal; exposed as `MoyoDataset.prim_std_cell`. |
| `transformation` | $(\mathbf{P}_{\mathrm{std}},\mathbf{p}_{\mathrm{std}})$, selecting the conventional coordinates before refinement. |
| `prim_transformation` | $(\mathbf{P}_{\mathrm{prim}},\mathbf{p}_{\mathrm{prim}})$, selecting the primitive coordinates before refinement. |
| `rotation_matrix` | The common proper Cartesian rotation applied to both refined output lattices; identity when `rotate_basis=false`. |
| `site_mapping` | Maps each conventional output site to its corresponding primitive output site. |
| `wyckoffs` | One Wyckoff position for each conventional output site, in that site's order and in the selected setting. |

The coordinate transformations preserve the selected setting while refinement enforces its symmetry.
Recovering the original distorted lattice or positions requires retaining the input cell.

## `rotate_basis` option in `MoyoDataset::new`

This option controls the Cartesian orientation of the refined lattice.
Let $\mathbf{A}'_{\mathrm{std}}=\mathbf{A}\mathbf{P}_{\mathrm{std}}$ be the conventional basis before refinement and let $\mathbf{B}$ be the refined basis in the canonical orientation described below, with the same handedness.
Define the right polar decomposition

$$
\mathbf{F}=\mathbf{B}(\mathbf{A}'_{\mathrm{std}})^{-1}
=\mathbf{R}\mathbf{U},\qquad
\mathbf{R}^\mathsf{T}\mathbf{R}=\mathbf{I},\quad
\det\mathbf{R}=1,\quad
\mathbf{U}=\mathbf{U}^\mathsf{T}>0.
$$

The symmetric positive-definite stretch $\mathbf{U}$ describes lattice refinement in the input Cartesian frame.
The proper rotation $\mathbf{R}$ takes that refined lattice to the canonical orientation.

| Option | Returned conventional basis | Returned rotation matrix |
| --- | --- | --- |
| `rotate_basis=true` (default) | $\mathbf{A}_{\mathrm{std}}=\mathbf{B}=\mathbf{R}\mathbf{U}\mathbf{A}\mathbf{P}_{\mathrm{std}}$ | $\mathbf{R}$ |
| `rotate_basis=false` | $\mathbf{A}_{\mathrm{std}}=\mathbf{R}^\mathsf{T}\mathbf{B}=\mathbf{U}\mathbf{A}\mathbf{P}_{\mathrm{std}}$ | $\mathbf{I}$ |

The same stretch and rotation apply to the primitive output.
Fractional positions are refined in the selected coordinates and are identical for both orientation options.
When the input lattice already satisfies the target symmetry, $\mathbf{U}=\mathbf{I}$ up to roundoff, recovering the relation $\mathbf{A}_{\mathrm{std}}=\mathbf{R}\mathbf{A}\mathbf{P}_{\mathrm{std}}$ for `rotate_basis=true`.
The current implementation uses this last relation even for distorted lattices, as noted in the implementation status above.

## Standardized cell with `setting=Setting.Standard`, `rotate_basis=true`, and right-handed input basis vectors

When `setting=Setting.Standard`, the conventional output gives the space group in the ITA setting.
With `rotate_basis=true`, its refined basis $\mathbf{A}_{\mathrm{std}}=\mathbf{B}$ has the form below.
The parameters $a$, $b$, and $c$ denote positive lengths.

!!! caution
    The Rust field `Cell.lattice.basis` stores basis vectors as columns.
    The other language bindings expose basis vectors as rows; transpose those arrays when using the equations on this page.

| Crystal family | "Conventional" basis vectors $\mathbf{A}_{\mathrm{std}}$ <br> (`MoyoDataset.std_cell`) | Additional conditions                                             |
| -------------- | -------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Triclinic      | $\begin{pmatrix} a_x & b_x & c_x \\ 0 & b_y & c_y \\ 0 & 0 & c_z \end{pmatrix}$        | Niggli reduced [^std-cell-2]; $a_x, b_y, c_z \gt 0$ |
| Monoclinic     | $\begin{pmatrix} a & 0 & c \cos \beta \\ 0 & b & 0 \\ 0 & 0 & c \sin \beta \end{pmatrix}$    | $a, b, c \sin \beta \gt 0$; $\cos \beta \le 0$ [^std-cell-7] |
| Orthorhombic   | $\begin{pmatrix} a & 0 & 0 \\ 0 & b & 0 \\ 0 & 0 & c \end{pmatrix}$                          | $a, b, c \gt 0$; $a \le b \le c$ as far as possible [^std-cell-6] |
| Tetragonal     | $\begin{pmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & c \end{pmatrix}$                          | $a, c \gt 0$                                        |
| Hexagonal      | $\begin{pmatrix} a & -a / 2 & 0 \\0 & \sqrt{3} a / 2 & 0 \\ 0 & 0 & c \end{pmatrix}$         | $a, c > 0$                                          |
| Cubic          | $\begin{pmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & a \end{pmatrix}$                          | $a > 0$                                             |

For left-handed input, use the same metric and negate the final Cartesian row of the displayed upper-triangular basis.
Equivalently, left-multiply by $\operatorname{diag}(1,1,-1)$, making the final diagonal entry negative while preserving all lengths and angles.
This handedness convention also applies to other Hall settings; their canonical basis is upper triangular, with the first two diagonal entries positive.

## Primitive standardized cell with `setting=Setting.Standard`, `rotate_basis=true`, and right-handed input basis vectors

Let $\mathbf{P}_{\mathrm{prim}}$ be `MoyoDataset.prim_std_linear` and $\mathbf{p}_{\mathrm{prim}}$ be `MoyoDataset.prim_std_origin_shift`.
The transformation $(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}})$ selects the primitive coordinates before refinement.
Moyo chooses a transformation matrix $\mathbf{Q}$ from a primitive cell to the standardized cell as the following table.

| Crystal family | Bravais class | Transformation matrix from primitive to conventional, $\mathbf{Q}$                    | $\mathbf{Q}^{-1}$                                                                                            |
| -------------- | ------------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Triclinic      | aP            | $\mathbf{Q}_P = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$    | $\mathbf{Q}_P^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$                      |
| Monoclinic     | mP            | $\mathbf{Q}_P$                                                                        | $\mathbf{Q}_P^{-1}$                                                                                          |
|                | mC            | $\mathbf{Q}_C = \begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$   | $\mathbf{Q}_C^{-1} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ -1/2 & 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$             |
| Orthorhombic   | oP            | $\mathbf{Q}_P$                                                                        | $\mathbf{Q}_P^{-1}$                                                                                          |
|                | oS            | $\mathbf{Q}_C$                                                                        | $\mathbf{Q}_C^{-1}$                                                                                          |
|                | oF            | $\mathbf{Q}_F = \begin{pmatrix} -1 & 1 & 1 \\ 1 & -1 & 1 \\ 1 & 1 & -1 \end{pmatrix}$ | $\mathbf{Q}_F^{-1} = \begin{pmatrix} 0 & 1/2 & 1/2 \\ 1/2 & 0 & 1/2 \\ 1/2 & 1/2 & 0 \end{pmatrix}$          |
|                | oI            | $\mathbf{Q}_I = \begin{pmatrix} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix}$    | $\mathbf{Q}_I^{-1} = \begin{pmatrix} -1/2 & 1/2 & 1/2 \\ 1/2 & -1/2 & 1/2 \\ 1/2 & 1/2 & -1/2 \end{pmatrix}$ |
| Tetragonal     | tP            | $\mathbf{Q}_P$                                                                        | $\mathbf{Q}_P^{-1}$                                                                                          |
|                | tI            | $\mathbf{Q}_I$                                                                        | $\mathbf{Q}_I^{-1}$                                                                                          |
| Hexagonal      | hR            | $\mathbf{Q}_R = \begin{pmatrix} 1 & 0 & 1 \\ -1 & 1 & 1 \\ 0 & -1 & 1 \end{pmatrix}$  | $\mathbf{Q}_R^{-1} = \begin{pmatrix} 2/3 & -1/3 & -1/3 \\ 1/3 & 1/3 & -2/3 \\ 1/3 & 1/3 & 1/3 \end{pmatrix}$ |
|                | hP            | $\mathbf{Q}_P$                                                                        | $\mathbf{Q}_P^{-1}$                                                                                          |
| Cubic          | cP            | $\mathbf{Q}_P$                                                                        | $\mathbf{Q}_P^{-1}$                                                                                          |
|                | cF            | $\mathbf{Q}_F$                                                                        | $\mathbf{Q}_F^{-1}$                                                                                          |
|                | cI            | $\mathbf{Q}_I$                                                                        | $\mathbf{Q}_I^{-1}$                                                                                          |

The conventional and primitive outputs share an origin and satisfy the following relations for either orientation option:

$$
(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}}) = (\mathbf{P}_{\mathrm{std}}, \mathbf{p}_{\mathrm{std}}) (\mathbf{Q}, \mathbf{0})^{-1}
$$

$$
\mathbf{A}_{\mathrm{prim}} = \mathbf{A}_{\mathrm{std}} \mathbf{Q}^{-1}.
$$

Thus $\mathbf{P}_{\mathrm{std}}=\mathbf{P}_{\mathrm{prim}}\mathbf{Q}$ and $\mathbf{p}_{\mathrm{std}}=\mathbf{p}_{\mathrm{prim}}$.
The fixed centering matrix $\mathbf{Q}$ changes the description of the refined crystal; it introduces no further refinement.

For `StandardizedCell`, the primitive output preserves the input primitive cell's site order and species.
If $m(j)=$ `site_mapping[j]`, then every conventional output site $j$ has the same species as primitive output site $m(j)$ and satisfies

$$
\mathbf{Q}\mathbf{x}_{\mathrm{std},j}-\mathbf{x}_{\mathrm{prim},m(j)}\in\mathbb{Z}^3.
$$

Every primitive site has $|\det\mathbf{Q}|$ conventional copies, and `site_mapping` and `wyckoffs` each have one entry per conventional site.
`MoyoDataset.mapping_std_prim` maps original input sites to primitive output sites, and `MoyoDataset.wyckoffs` is indexed by original input sites.

[^setting]: That being said, the order of the Hall symbols are the same as Table A1.4.2.7 in International Tables for Crystallography Volume B (2010).

[^std-cell-2]: Applied regardless of `rotate_basis` value.

[^std-cell-7]: The basis vectors $\mathbf{a}$ and $\mathbf{c}$ are taken from the Delaunay-reduced triple $\mathbf{v}_1, \mathbf{v}_2, -(\mathbf{v}_1 + \mathbf{v}_2)$ of the lattice plane perpendicular to the unique axis, whose members are pairwise non-acute. Among the pairs that keep the Hall setting (the centering and the glide translations, up to an origin shift), moyo chooses the one with $\beta$ closest to $\pi / 2$, prefers the non-acute value ($\pi / 2 \le \beta \lt \pi$, i.e. $\cos \beta \le 0$) between the supplements following the ITA convention, and finally the lexicographically smallest $(a, b, c)$, which gives $a \le c$ for the settings that allow the $\mathbf{a} \leftrightarrow \mathbf{c}$ swap ($P2$, $P2_1$, $Pm$, $P2/m$, $P2_1/m$) following E. Parthe and L. M. Gelato, Acta Cryst. A**39**, 169-173 (1983), as spglib does. Consequently $\pi / 2 \le \beta \le 2\pi / 3$.

[^std-cell-6]: moyo orders the basis vectors as $a \le b \le c$ as far as the space-group setting allows. Among the six axis permutations, only those that preserve the centering and map the space group onto itself up to an origin shift (i.e. elements of the affine normalizer) are admissible, and moyo picks the admissible one with the lexicographically smallest $(a, b, c)$. Full ordering is not always attainable; for example, side-face-centered cells (oS) admit only the $\mathbf{a} \leftrightarrow \mathbf{b}$ swap, enforcing $a \le b$ alone.
