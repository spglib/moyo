# Standardization of crystal structures

Standardization describes a crystal in a conventional coordinate system and refines its lattice and atomic positions to satisfy the identified symmetry.
This page specifies the intended conventional and primitive standardized cells, including their basis, origin, and Cartesian orientation.

All basis matrices on this page store lattice vectors as columns.
The Python and C interfaces use row-wise lattice arrays; transpose them when applying these equations.

## Space-group setting

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

## Returned-cell specification

### Target symmetry and numerical accuracy

The conventional and primitive standardized cells must satisfy the space group of the selected Hall symbol, expressed in their respective coordinate systems.

For either returned cell, let $\mathbf{A}\sb{\ast}$ contain its basis vectors as columns and let $\mathbf{G}\sb{\ast} = \mathbf{A}\sb{\ast}^\mathsf{T}\mathbf{A}\sb{\ast}$.
Every target operation $(\mathbf{W}, \mathbf{w})$ must satisfy

$$
\mathbf{W}^\mathsf{T}\mathbf{G}_*\mathbf{W} = \mathbf{G}_*.
$$

It must also induce a species-preserving permutation $\pi$ of the returned sites such that, for every site $i$,

$$
\mathbf{W}\mathbf{x}_i + \mathbf{w} - \mathbf{x}_{\pi(i)} \in \mathbb{Z}^3.
$$

These are exact identities in real arithmetic.
In floating-point arithmetic, they must hold within numerical roundoff, independently of the tolerances used to recognize symmetry in the input.

Both values of `rotate_basis` must satisfy this specification.
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
They preserve the crystal's geometry, while refinement can change its metric and atomic positions to enforce symmetry.

## Cartesian orientation

The `rotate_basis` option controls the Cartesian orientation of the refined lattice.
Let $\mathbf{A}'\sb{\mathrm{std}}=\mathbf{A}\mathbf{P}\sb{\mathrm{std}}$ be the conventional basis before refinement and let $\mathbf{B}$ be the refined basis in the canonical orientation described below, with the same handedness.
Define the right polar decomposition

$$
\mathbf{F}=\mathbf{B}(\mathbf{A}'_{\mathrm{std}})^{-1}
=\mathbf{R}\mathbf{U},\qquad
\mathbf{R}^\mathsf{T}\mathbf{R}=\mathbf{I},\quad
\det\mathbf{R}=1,\quad
\mathbf{U}=\mathbf{U}^\mathsf{T}\gt 0.
$$

The symmetric positive-definite stretch $\mathbf{U}$ describes lattice refinement in the input Cartesian frame.
The proper rotation $\mathbf{R}$ takes that refined lattice to the canonical orientation.

With `rotate_basis=true` (default), the conventional basis is

$$
\mathbf{A}_{\mathrm{std}}=\mathbf{B}=\mathbf{R}\mathbf{U}\mathbf{A}\mathbf{P}_{\mathrm{std}}.
$$

With `rotate_basis=false`, the conventional basis is

$$
\mathbf{A}_{\mathrm{std}}=\mathbf{R}^\mathsf{T}\mathbf{B}=\mathbf{U}\mathbf{A}\mathbf{P}_{\mathrm{std}}.
$$

The same stretch and rotation apply to the primitive output.
Fractional positions are refined in the selected coordinates and are identical for both orientation options.
When the input lattice already satisfies the target symmetry, $\mathbf{U}=\mathbf{I}$ up to roundoff, recovering the relation $\mathbf{A}\sb{\mathrm{std}}=\mathbf{R}\mathbf{A}\mathbf{P}\sb{\mathrm{std}}$ for `rotate_basis=true`.

## Conventional standardized cell

With `Setting.Standard`, the conventional cell gives the space group in the ITA setting.
For `rotate_basis=true` and right-handed input, its refined basis $\mathbf{A}\sb{\mathrm{std}}=\mathbf{B}$ has the form below.
The parameters $a$, $b$, and $c$ denote positive lengths.

### Triclinic

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a_x & b_x & c_x \cr
0 & b_y & c_y \cr
0 & 0 & c_z
\end{pmatrix}.
$$

Niggli reduced [^std-cell-2]; $a\sb{x}, b\sb{y}, c\sb{z} \gt 0$.

### Monoclinic

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a & 0 & c \cos \beta \cr
0 & b & 0 \cr
0 & 0 & c \sin \beta
\end{pmatrix}.
$$

$a, b, c \sin \beta \gt 0$; $\cos \beta \le 0$ (see [basis selection](#monoclinic-basis-selection)).

### Orthorhombic

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a & 0 & 0 \cr
0 & b & 0 \cr
0 & 0 & c
\end{pmatrix}.
$$

$a, b, c \gt 0$; $a \le b \le c$ as far as possible (see [basis selection](#orthorhombic-basis-selection)).

### Tetragonal

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a & 0 & 0 \cr
0 & a & 0 \cr
0 & 0 & c
\end{pmatrix}.
$$

$a, c \gt 0$.

### Hexagonal

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a & -a / 2 & 0 \cr
0 & \sqrt{3} a / 2 & 0 \cr
0 & 0 & c
\end{pmatrix}.
$$

$a, c \gt 0$.

### Cubic

$$
\mathbf{A}_{\mathrm{std}} = \begin{pmatrix}
a & 0 & 0 \cr
0 & a & 0 \cr
0 & 0 & a
\end{pmatrix}.
$$

$a \gt 0$.

### Handedness

For left-handed input, use the same metric and negate the final Cartesian row of the displayed upper-triangular basis.
Equivalently, left-multiply by $\mathrm{diag}(1,1,-1)$, making the final diagonal entry negative while preserving all lengths and angles.
This handedness convention also applies to other Hall settings; their canonical basis is upper triangular, with the first two diagonal entries positive.

## Primitive standardized cell

The primitive and conventional cells describe the same refined crystal.
For `Setting.Standard`, the change of basis from primitive to conventional is given by the centering matrix $\mathbf{Q}$ below.

| Crystal family | Bravais class | Centering |
| --- | --- | --- |
| Triclinic | aP | [P](#p-centering) |
| Monoclinic | mP | [P](#p-centering) |
|  | mC | [C](#c-centering) |
| Orthorhombic | oP | [P](#p-centering) |
|  | oS | [C](#c-centering) |
|  | oF | [F](#f-centering) |
|  | oI | [I](#i-centering) |
| Tetragonal | tP | [P](#p-centering) |
|  | tI | [I](#i-centering) |
| Hexagonal | hR | [R](#r-centering) |
|  | hP | [P](#p-centering) |
| Cubic | cP | [P](#p-centering) |
|  | cF | [F](#f-centering) |
|  | cI | [I](#i-centering) |

The conventional and primitive outputs share an origin and satisfy the following relations for either orientation option:

$$
\mathbf{A}_{\mathrm{prim}} = \mathbf{A}_{\mathrm{std}} \mathbf{Q}^{-1},\qquad
\mathbf{P}_{\mathrm{prim}} = \mathbf{P}_{\mathrm{std}} \mathbf{Q}^{-1},\qquad
\mathbf{p}_{\mathrm{prim}} = \mathbf{p}_{\mathrm{std}}.
$$

Here $(\mathbf{P}\sb{\mathrm{std}},\mathbf{p}\sb{\mathrm{std}})$ and $(\mathbf{P}\sb{\mathrm{prim}},\mathbf{p}\sb{\mathrm{prim}})$ select the conventional and primitive coordinates before refinement.
The centering matrix changes the description of the refined crystal; it introduces no further refinement.

Each conventional site corresponds to a primitive site of the same species, with fractional positions satisfying

$$
\mathbf{Q}\mathbf{x}_{\mathrm{std}}-\mathbf{x}_{\mathrm{prim}}\in\mathbb{Z}^3.
$$

Each primitive site has $|\det\mathbf{Q}|$ conventional copies.

### P centering

$$
\mathbf{Q}_P = \begin{pmatrix}
1 & 0 & 0 \cr
0 & 1 & 0 \cr
0 & 0 & 1
\end{pmatrix},\qquad
\mathbf{Q}_P^{-1} = \begin{pmatrix}
1 & 0 & 0 \cr
0 & 1 & 0 \cr
0 & 0 & 1
\end{pmatrix}.
$$

### C centering

$$
\mathbf{Q}_C = \begin{pmatrix}
1 & -1 & 0 \cr
1 & 1 & 0 \cr
0 & 0 & 1
\end{pmatrix},\qquad
\mathbf{Q}_C^{-1} = \begin{pmatrix}
1/2 & 1/2 & 0 \cr
-1/2 & 1/2 & 0 \cr
0 & 0 & 1
\end{pmatrix}.
$$

### F centering

$$
\mathbf{Q}_F = \begin{pmatrix}
-1 & 1 & 1 \cr
1 & -1 & 1 \cr
1 & 1 & -1
\end{pmatrix},\qquad
\mathbf{Q}_F^{-1} = \begin{pmatrix}
0 & 1/2 & 1/2 \cr
1/2 & 0 & 1/2 \cr
1/2 & 1/2 & 0
\end{pmatrix}.
$$

### I centering

$$
\mathbf{Q}_I = \begin{pmatrix}
0 & 1 & 1 \cr
1 & 0 & 1 \cr
1 & 1 & 0
\end{pmatrix},\qquad
\mathbf{Q}_I^{-1} = \begin{pmatrix}
-1/2 & 1/2 & 1/2 \cr
1/2 & -1/2 & 1/2 \cr
1/2 & 1/2 & -1/2
\end{pmatrix}.
$$

### R centering

$$
\mathbf{Q}_R = \begin{pmatrix}
1 & 0 & 1 \cr
-1 & 1 & 1 \cr
0 & -1 & 1
\end{pmatrix},\qquad
\mathbf{Q}_R^{-1} = \begin{pmatrix}
2/3 & -1/3 & -1/3 \cr
1/3 & 1/3 & -2/3 \cr
1/3 & 1/3 & 1/3
\end{pmatrix}.
$$

[^setting]: That being said, the order of the Hall symbols are the same as Table A1.4.2.7 in International Tables for Crystallography Volume B (2010).

[^std-cell-2]: Applied regardless of `rotate_basis` value.

## Basis-selection notes

### Monoclinic basis selection

The basis vectors $\mathbf{a}$ and $\mathbf{c}$ are taken from the Delaunay-reduced triple $\mathbf{v}\sb{1}, \mathbf{v}\sb{2}, -(\mathbf{v}\sb{1} + \mathbf{v}\sb{2})$ of the lattice plane perpendicular to the unique axis, whose members are pairwise non-acute. Among the pairs that keep the Hall setting (the centering and the glide translations, up to an origin shift), moyo chooses the one with $\beta$ closest to $\pi / 2$, prefers the non-acute value ($\pi / 2 \le \beta \lt \pi$, i.e. $\cos \beta \le 0$) between the supplements following the ITA convention, and finally the lexicographically smallest $(a, b, c)$, which gives $a \le c$ for the settings that allow the $\mathbf{a} \leftrightarrow \mathbf{c}$ swap ($P2$, $P2\sb{1}$, $Pm$, $P2/m$, $P2\sb{1}/m$) following E. Parthe and L. M. Gelato, Acta Cryst. A**39**, 169-173 (1983), as spglib does. Consequently $\pi / 2 \le \beta \le 2\pi / 3$.

### Orthorhombic basis selection

moyo orders the basis vectors as $a \le b \le c$ as far as the space-group setting allows. Among the six axis permutations, only those that preserve the centering and map the space group onto itself up to an origin shift (i.e. elements of the affine normalizer) are admissible, and moyo picks the admissible one with the lexicographically smallest $(a, b, c)$. Full ordering is not always attainable; for example, side-face-centered cells (oS) admit only the $\mathbf{a} \leftrightarrow \mathbf{b}$ swap, enforcing $a \le b$ alone.
