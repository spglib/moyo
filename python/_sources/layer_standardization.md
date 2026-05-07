# Moyo conventions of standardized layer cell

[!NOTE]
This document describes how the standardized layer cell `MoyoLayerDataset.std_cell` and the primitive standardized layer cell `MoyoLayerDataset.prim_std_cell` are specified in Moyo. The conventions described here apply uniformly to the core Rust implementation as well as all other language bindings.

For the bulk space-group standardization, see [Moyo conventions of standardized cell](./standardization.md). The layer-group convention deviates from the bulk one in two places: the third basis vector is the aperiodic stacking direction (not a lattice translation), and the only centerings are `p` (primitive) and `c` (rectangular-centered).

## Basis vectors and the aperiodic axis

Three basis vectors $(\mathbf{a}_s, \mathbf{b}_s, \mathbf{c}_s)$ in the same column-vector convention as the input cell.

- $\mathbf{c}_s$ is the aperiodic stacking direction. By construction it is perpendicular to $\mathbf{a}_s$ and $\mathbf{b}_s$.
- The magnitude of $\mathbf{c}_s$ is preserved from the input: $|\mathbf{c}_s| = |\mathbf{c}|$. Moyo never rescales $\mathbf{c}$.

## `setting` option in `MoyoLayerDataset::new`

A *setting* selects the representative of each layer-group type. Moyo supports three settings, mirroring the bulk space-group case but using layer Hall numbers (1 - 116) instead of bulk Hall numbers (1 - 530):

- `LayerSetting.Standard`: the BCS / ITE convention. Uses cell choice 1 for monoclinic-oblique LGs with multiple cell choices in ITE (LG 5, 7), origin choice 2 for centrosymmetric LGs with two ITE origins (LG 52, 62, 64), monoclinic-rectangular `:a` axis labelling (LG 8-18), and orthorhombic `abc` axis labelling (LG 19-48).
- `LayerSetting.Spglib`: the smallest layer Hall number for each LG (spglib's first row per LG). Differs from `Standard` only at LG 52, 62, 64 (origin choice 1 vs 2).
- `LayerSetting.HallNumber`: allows users to specify a layer Hall number (an integer from 1 to 116) directly.

Moyo chooses `LayerSetting.Standard` as the default. This change of default behaviour from spglib's affects the centrosymmetric LGs 52, 62, 64 (origin choice 2 vs 1).

## `rotate_basis` option in `MoyoLayerDataset::new`

When `rotate_basis=true` (default), Moyo rotates the standardized basis so that:

- $\mathbf{c}_s$ lies along Cartesian $z$,
- $\mathbf{a}_s$ lies along Cartesian $x$,
- $\mathbf{b}_s$ lies in the $xy$-plane.

When `rotate_basis=false`, the basis directions inherit the input orientation, modulo the LG-canonical relabelling of $(\mathbf{a}, \mathbf{b})$.

## Standardized cell with `setting=LayerSetting.Standard`, `rotate_basis=true`, and right-handed input basis vectors

Let $\mathbf{P}_{\mathrm{std}}$ be `MoyoLayerDataset.std_linear` and $\mathbf{p}_{\mathrm{std}}$ be `MoyoLayerDataset.std_origin_shift`. The transformation $(\mathbf{P}_{\mathrm{std}}, \mathbf{p}_{\mathrm{std}})$ brings the input layer cell into the standardized layer cell (`MoyoLayerDataset.std_cell`).

Let $\mathbf{A}$ be the column-wise basis vectors of the input layer cell and $\mathbf{A}_{\mathrm{std}}$ be the column-wise basis vectors of the standardized layer cell. When `rotate_basis=true`, the following relation holds:

$$
\mathbf{A}_{\mathrm{std}} = \mathbf{R} \mathbf{A} \mathbf{P}_{\mathrm{std}},
$$

where $\mathbf{R}$ is `MoyoLayerDataset.std_rotation_matrix`, a proper rotation matrix that brings the input basis into the orientation described above.

The third column of $\mathbf{A}_{\mathrm{std}}$ — that is, $\mathbf{c}_s$ — is parallel to the input $\mathbf{c}$ and equal in length: $|\mathbf{c}_s| = |\mathbf{c}|$.

The standardized layer-cell basis vectors $\mathbf{A}_{\mathrm{std}}$ are defined for each LG crystal system as the following table.

| LG crystal system        | "Conventional" basis vectors $\mathbf{A}_{\mathrm{std}}$                            | Additional conditions                                         |
| ------------------------ | ----------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Triclinic / oblique      | $\begin{pmatrix} a_x & b_x & 0 \\ 0 & b_y & 0 \\ 0 & 0 & c \end{pmatrix}$           | 2D Minkowski reduced [^lstd-1]; $a_x, b_y, c \gt 0$ [^lstd-2] |
| Monoclinic / oblique     | $\begin{pmatrix} a_x & b_x & 0 \\ 0 & b_y & 0 \\ 0 & 0 & c \end{pmatrix}$           | 2D Minkowski reduced [^lstd-1]; $a_x, b_y, c \gt 0$ [^lstd-2] |
| Monoclinic / rectangular | $\begin{pmatrix} a & 0 & 0 \\ 0 & b & 0 \\ 0 & 0 & c \end{pmatrix}$                 | $a, b, c \gt 0$ [^lstd-2] [^lstd-3]                           |
| Orthorhombic             | $\begin{pmatrix} a & 0 & 0 \\ 0 & b & 0 \\ 0 & 0 & c \end{pmatrix}$                 | $a, b, c \gt 0$ [^lstd-2]                                     |
| Tetragonal               | $\begin{pmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & c \end{pmatrix}$                 | $a, c \gt 0$ [^lstd-2]                                        |
| Trigonal / hexagonal     | $\begin{pmatrix} a & -a/2 & 0 \\ 0 & \sqrt{3} a / 2 & 0 \\ 0 & 0 & c \end{pmatrix}$ | $a, c \gt 0$ [^lstd-2]                                        |

Which of $\mathbf{a}$ or $\mathbf{b}$ carries the 2-fold or mirror normal in the monoclinic-rectangular case is a labelling choice exposed via `LayerSetting` (paper Fu et al. 2024 Table 5 codes `:a`, `:b`; analogously `:abc` / `:b-ac` for orthorhombic). It is not part of the reduction step.

## Centering and the primitive standardized cell

Layer groups admit only two centerings: `p` (primitive) and `c` (rectangular-centered, in-plane only). There is no body, face, or rhombohedral centering. The conventional cell may double the primitive 2D cell area in the `c` case, but never extends along the aperiodic axis.

Let $\mathbf{Q}$ be the transformation matrix from the primitive layer cell to the conventional layer cell:

| LG Bravais class       | Centering | $\mathbf{Q}$                                                                        | $\mathbf{Q}^{-1}$                                                                                |
| ---------------------- | --------- | ----------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `mp`, `op`, `tp`, `hp` | `p`       | $\mathbf{Q}_p = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$  | $\mathbf{Q}_p^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$          |
| `oc`                   | `c`       | $\mathbf{Q}_c = \begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$ | $\mathbf{Q}_c^{-1} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ -1/2 & 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$ |

Note that $\mathbf{Q}_c$ and $\mathbf{Q}_c^{-1}$ leave the third column (the aperiodic axis) untouched.

Let $\mathbf{P}_{\mathrm{prim}}$ be `MoyoLayerDataset.prim_std_linear` and $\mathbf{p}_{\mathrm{prim}}$ be `MoyoLayerDataset.prim_std_origin_shift`. The transformation $(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}})$ brings the input layer cell into the primitive standardized layer cell (`MoyoLayerDataset.prim_std_cell`):

$$
(\mathbf{P}_{\mathrm{prim}}, \mathbf{p}_{\mathrm{prim}}) = (\mathbf{P}_{\mathrm{std}}, \mathbf{p}_{\mathrm{std}}) (\mathbf{Q}, \mathbf{0})^{-1}.
$$

## Pearson symbol

`MoyoLayerDataset.pearson_symbol` is built as `<2D Bravais type><N>`, where the Bravais type is one of `mp`, `op`, `oc`, `tp`, `hp` and `N` is the conventional-cell atom count.

## Atom positions

Atom positions in `MoyoLayerDataset.std_cell` and `MoyoLayerDataset.prim_std_cell` are symmetrised to the orbits prescribed by the identified layer group, mirroring the bulk space-group standardisation.

[^lstd-1]: Applied regardless of `rotate_basis` value.

[^lstd-2]: $c < 0$ for left-handed input basis vectors.

[^lstd-3]: For monoclinic-rectangular LGs (LG 8-18) the in-plane lattice is rectangular by Bravais constraint, so 2D Minkowski reduction is trivial (length sort / sign fixing only).
