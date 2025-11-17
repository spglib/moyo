# Moyo conventions of standardized cell

[!NOTE]
This document is based on moyo==0.7.3.

This document describes how the standardized cell `MoyoDataset.std_cell` and the primitive standardized cell `MoyoDataset.prim_std_cell` are specified in Moyo.
The conventions described here apply uniformly to the core Rust implementation as well as all other language bindings.

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

## `rotate_basis` option in `MoyoDataset::new`

Because the basis vectors of the input crystal structure are not assumed to align with Cartesian x, y, and z axes, moyo rotates the basis vectors at `rotate_basis=True` (default) as follows:

## Standardized cell with `setting=Setting.Standard`, `rotate_basis=true`, and right-handed input basis vectors

Let $\mathbf{P}_{\mathrm{std}}$ be `MoyoDataset.std_linear` and $\mathbf{p}_{\mathrm{std}}$ be `MoyoDataset.std_origin_shift`.
The transformation $(\mathbf{P}_{\mathrm{std}}, \mathbf{p}_{\mathrm{std}})$ brings the input crystal structure into the standardized cell (`MoyoDataset.std_cell`).
When `setting=Setting.Standard`, the standardized cell gives the space group in the ITA setting.

Let $\mathbf{A}$ be the column-wise basis vectors of the input crystal structure, `Cell.basis`, and $\mathbf{A}_{\mathbf{std}}$ be the column-wise basis vectors of the standardized cell, `MoyoDataset.std_cell.basis`.
When `rotate_basis=true`, the following relation holds:
$
\mathbf{A}_{\mathbf{std}} = \mathbf{R} \mathbf{A} \mathbf{P}_{\mathrm{std}},
$
where $\mathbf{R}$ is `MoyoDataset.std_rotation_matrix`.
Here, $\mathbf{R}$ is a proper rotation matrix that brings the input basis vectors into a certain orientation depending on the crystal family.

[!CAUTION]
The Rust-implementation attribute `Cell.basis` stores the basis vectors in column-wise manner, while the other language bindings store them in row-wise manner.

The standardized cell basis vectors $\mathbf{A}_{\mathrm{std}}$ are defined for each crystal family as the following table.
Note that `rotate_basis=true` is assumed, and that the input basis vectors are right-handed.

| Crystal family | "Conventional" basis vectors $\mathbf{A}_{\mathrm{std}}$ <br> (`MoyoDataset.std_cell.basis`) | Additional conditions                                             |
| -------------- | -------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Triclinic      | $\begin{pmatrix} a_x & b_x & c_x \\ a_y & b_y & c_y \\ a_z & b_z & c_z \end{pmatrix}$        | Niggli reduced [^std-cell-2]; $a_x, b_y, c_z \gt 0$ [^std-cell-3] |
| Monoclinic     | $\begin{pmatrix} a & 0 & c \cos \beta \\ 0 & b & 0 \\ 0 & 0 & c \sin \beta \end{pmatrix}$    | $a, b, c \sin \beta \gt 0$ [^std-cell-4] [^std-cell-7]            |
| Orthorhombic   | $\begin{pmatrix} a & 0 & 0 \\ 0 & b & 0 \\ 0 & 0 & c \end{pmatrix}$                          | $a, b, c \gt 0$ [^std-cell-4] [^std-cell-6]                       |
| Tetragonal     | $\begin{pmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & c \end{pmatrix}$                          | $a, c \gt 0$ [^std-cell-4]                                        |
| Hexagonal      | $\begin{pmatrix} a & -a / 2 & 0 \\0 & \sqrt{3} a / 2 & 0 \\ 0 & 0 & c \end{pmatrix}$         | $a, c > 0$ [^std-cell-4]                                          |
| Cubic          | $\begin{pmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & a \end{pmatrix}$                          | $a > 0$ [^std-cell-5]                                             |

## Primitive standardized cell with `setting=Setting.Standard`, `rotate_basis=true`, and right-handed input basis vectors

Let $\mathbf{P}_{\mathrm{prim}}$ be `MoyoDataset.prim_std_linear` and $\mathbf{p}_{\mathrm{prim}}$ be `MoyoDataset.prim_std_origin_shift`.
The transformation $(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}})$ brings the input crystal structure into the primitive standardized cell (`MoyoDataset.prim_std_cell`).
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

Depending on the Bravais class of the standardized cell, the following relation holds:

$
(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}}) = (\mathbf{P}_{\mathrm{std}}, \mathbf{p}_{\mathrm{std}}) (\mathbf{Q}, \mathbf{0})^{-1}
$

$
\mathbf{A}_{\mathrm{prim}} = \mathbf{A}_{\mathrm{std}} \mathbf{Q}^{-1}.
$

[^setting]: That being said, the order of the Hall symbols are the same as Table A1.4.2.7 in International Tables for Crystallography Volume B (2010).

[^std-cell-2]: Applied regardless of `rotate_basis` value.

[^std-cell-3]: $c_z < 0$ for left-handed input basis vectors.

[^std-cell-4]: $c < 0$ for left-handed input basis vectors.

[^std-cell-7]: moyo tries to bring $\beta$ closer to $\pi / 2$ as much as possible.

[^std-cell-6]: Currently, moyo does not enforce the order in $a$, $b$, and $c$.

[^std-cell-5]: Negative (z, z)-component for left-handed input basis vectors.
