use std::collections::HashMap;

use log::{debug, warn};
use nalgebra::linalg::QR;
use nalgebra::{Matrix3, Vector3};

use super::standardize::orbits_in_cell;
use crate::base::{
    Cell, EPS, Lattice, Lattice2D, LayerCell, LayerLattice, Linear, MoyoError, Operations,
    Permutation, Position, Rotations, Transformation, UnimodularTransformation, project_rotations,
};
use crate::data::{
    LayerArithmeticNumber, LayerHallNumber, LayerHallSymbol, LayerLatticeSystem,
    LayerWyckoffPosition, WyckoffPositionSpace, iter_layer_wyckoff_positions,
    layer_arithmetic_crystal_class_entry, layer_hall_symbol_entry,
};
use crate::identify::LayerGroup;
use crate::math::{SNF, lift_2d_to_3d};

/// Standardized cell for a layer-group setting.
///
/// Mirrors [`super::StandardizedCell`] for the bulk space-group case.
/// Output contract: `c_s` is the aperiodic axis along cartesian z
/// (when `rotate_basis = true`), `a_s` is along cartesian x, `b_s` lies
/// in the xy-plane, `|c_s| = |c|` is preserved from the input, and the
/// in-plane `(a_s, b_s)` block obeys the LG crystal system's metric
/// conditions (paper Fu et al. 2024 Figure 1 / Appendix C).
#[allow(dead_code)] // wired up by M5 (`MoyoLayerDataset`)
pub(crate) struct StandardizedLayerCell {
    /// Primitive standardized layer cell.
    pub prim_layer_cell: LayerCell,
    /// Transformation from the input primitive layer cell to the primitive
    /// standardized layer cell.
    pub prim_transformation: UnimodularTransformation,
    /// Conventional standardized layer cell.
    pub layer_cell: LayerCell,
    /// Wyckoff positions of sites in `layer_cell`.
    pub wyckoffs: Vec<LayerWyckoffPosition>,
    /// Transformation from the input primitive layer cell to the conventional
    /// standardized layer cell.
    pub transformation: Transformation,
    /// Rigid rotation taking the input primitive layer cell's basis to the
    /// standardized one. Identity when `rotate_basis = false`.
    pub rotation_matrix: Matrix3<f64>,
    /// Mapping from the site in `layer_cell` to that in `prim_layer_cell`.
    pub site_mapping: Vec<usize>,
}

#[allow(dead_code)] // wired up by M5 (`MoyoLayerDataset`)
impl StandardizedLayerCell {
    /// Standardize the input **primitive** layer cell.
    ///
    /// Produces the LG-canonical conventional layer cell and assigns a
    /// Wyckoff letter / site-symmetry symbol to every atom. Mirrors
    /// [`super::StandardizedCell::new`] but uses the layer-group databases
    /// (lowercase `p`/`c` Hall symbols, [`LayerCentering`], LG arithmetic
    /// classes, layer Wyckoff positions) and a layer-shaped lattice
    /// symmetrization (Cholesky on the layer-block metric tensor) so the
    /// output's third axis stays the aperiodic stacking direction.
    pub fn new(
        prim_layer_cell: &LayerCell,
        prim_layer_operations: &Operations,
        prim_layer_permutations: &[Permutation],
        layer_group: &LayerGroup,
        symprec: f64,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<Self, MoyoError> {
        let (
            prim_std_layer_cell,
            prim_std_permutations,
            prim_transformation,
            std_layer_cell,
            transformation,
            rotation_matrix,
            site_mapping,
        ) = Self::standardize_and_symmetrize_cell(
            prim_layer_cell,
            prim_layer_operations,
            prim_layer_permutations,
            layer_group,
            epsilon,
            rotate_basis,
        )?;

        let wyckoffs = Self::assign_wyckoffs(
            &prim_std_layer_cell,
            &prim_std_permutations,
            &std_layer_cell,
            &site_mapping,
            layer_group.hall_number,
            symprec,
        )?;

        Ok(StandardizedLayerCell {
            prim_layer_cell: prim_std_layer_cell,
            prim_transformation,
            layer_cell: std_layer_cell,
            wyckoffs,
            transformation,
            rotation_matrix,
            site_mapping,
        })
    }

    #[allow(clippy::type_complexity)]
    fn standardize_and_symmetrize_cell(
        prim_layer_cell: &LayerCell,
        prim_layer_operations: &Operations,
        prim_layer_permutations: &[Permutation],
        layer_group: &LayerGroup,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<
        (
            LayerCell,
            Vec<Permutation>,
            UnimodularTransformation,
            LayerCell,
            Transformation,
            Matrix3<f64>,
            Vec<usize>,
        ),
        MoyoError,
    > {
        let entry = layer_hall_symbol_entry(layer_group.hall_number)
            .ok_or(MoyoError::StandardizationError)?;

        // Layer-canonical operations from the Hall symbol (conv + prim).
        let lh_symbol = LayerHallSymbol::from_hall_number(layer_group.hall_number)
            .ok_or(MoyoError::StandardizationError)?;
        let (conv_std_operations, prim_std_operations) =
            lh_symbol.traverse_and_primitive_traverse();

        // Choose the prim_transformation. For oblique LGs (LG 1-7) the
        // in-plane lattice has full freedom and we 2D-Minkowski-reduce it
        // *after* applying `layer_group.transformation`; mirrors the bulk
        // triclinic case which composes Niggli with the SG transformation.
        // For the constrained crystal systems (rectangular, square,
        // hexagonal) the in-plane metric is already pinned by the LG, so
        // the bare LG transformation suffices.
        let lattice_system = layer_lattice_system_for(entry.arithmetic_number)
            .ok_or(MoyoError::StandardizationError)?;
        let prim_transformation = match lattice_system {
            LayerLatticeSystem::Oblique => {
                standardize_oblique_layer_cell(prim_layer_cell, &layer_group.transformation)
            }
            _ => layer_group.transformation.clone(),
        };

        let prim_std_cell_tmp =
            transform_layer_cell_unimodular(prim_layer_cell, &prim_transformation);

        // Symmetrize positions in the primitive standardized cell using the
        // refined LG operations from the Hall symbol. Permutations come
        // from the M2 search, which uses a different operation order than
        // the database traversal -- align them by rotation key.
        let mut permutation_mapping = HashMap::new();
        let prim_rotations =
            project_rotations(&prim_transformation.transform_operations(prim_layer_operations));
        for (prim_rotation, permutation) in
            prim_rotations.iter().zip(prim_layer_permutations.iter())
        {
            permutation_mapping.insert(*prim_rotation, permutation.clone());
        }
        let prim_std_permutations = prim_std_operations
            .iter()
            .map(|ops| {
                permutation_mapping
                    .get(&ops.rotation)
                    .cloned()
                    .ok_or(MoyoError::StandardizationError)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let new_prim_std_positions = symmetrize_positions(
            &prim_std_cell_tmp,
            &prim_std_operations,
            &prim_std_permutations,
            epsilon,
        );

        let prim_std_cell = LayerCell::new_unchecked(
            prim_std_cell_tmp.lattice().clone(),
            new_prim_std_positions,
            prim_std_cell_tmp.numbers().to_vec(),
        );

        // To (conventional) standardized cell. Layer centering only acts in
        // the in-plane block; the (1,1) entry on the aperiodic axis means
        // the conventional cell is not extended along c.
        let conv_trans_linear: Linear = entry.centering.linear();
        let (std_cell, site_mapping) = transform_layer_cell_centering(
            &prim_std_cell,
            &Transformation::from_linear(conv_trans_linear),
        );

        // prim_transformation * (conv_trans_linear, 0)
        let transformation = Transformation::new(
            prim_transformation.linear * conv_trans_linear,
            prim_transformation.origin_shift,
        );

        if rotate_basis {
            // Symmetrize the in-plane metric using the conventional
            // operations' rotations (all in layer-block form by
            // construction of `LayerHallSymbol`).
            let layer_rotations = project_rotations(&conv_std_operations);
            let (new_inner_lattice, rotation_matrix) =
                symmetrize_layer_lattice(&inner_lattice(&std_cell), &layer_rotations);

            // Apply the same rigid rotation to the primitive cell.
            let new_prim_inner_lattice = inner_lattice(&prim_std_cell).rotate(&rotation_matrix);

            let prim_rotated = LayerCell::new_unchecked(
                LayerLattice::new_unchecked(new_prim_inner_lattice),
                prim_std_cell.positions().to_vec(),
                prim_std_cell.numbers().to_vec(),
            );
            let std_rotated = LayerCell::new_unchecked(
                LayerLattice::new_unchecked(new_inner_lattice),
                std_cell.positions().to_vec(),
                std_cell.numbers().to_vec(),
            );
            Ok((
                prim_rotated,
                prim_std_permutations,
                prim_transformation,
                std_rotated,
                transformation,
                rotation_matrix,
                site_mapping,
            ))
        } else {
            Ok((
                prim_std_cell,
                prim_std_permutations,
                prim_transformation,
                std_cell,
                transformation,
                Matrix3::identity(),
                site_mapping,
            ))
        }
    }

    fn assign_wyckoffs(
        prim_std_cell: &LayerCell,
        prim_std_permutations: &[Permutation],
        std_cell: &LayerCell,
        site_mapping: &[usize],
        hall_number: LayerHallNumber,
        symprec: f64,
    ) -> Result<Vec<LayerWyckoffPosition>, MoyoError> {
        // Group sites in the conventional cell by crystallographic orbit.
        let orbits = orbits_in_cell(
            prim_std_cell.num_atoms(),
            prim_std_permutations,
            site_mapping,
        );
        let mut num_orbits = 0;
        let mut mapping = vec![0; std_cell.num_atoms()];
        let mut remapping = vec![];
        for i in 0..std_cell.num_atoms() {
            if orbits[i] == i {
                mapping[i] = num_orbits;
                remapping.push(i);
                num_orbits += 1;
            } else {
                mapping[i] = mapping[orbits[i]];
            }
        }

        let mut multiplicities = vec![0; num_orbits];
        for i in 0..std_cell.num_atoms() {
            multiplicities[mapping[i]] += 1;
        }

        let mut representative_wyckoffs: Vec<Option<LayerWyckoffPosition>> = vec![None; num_orbits];
        let lattice = inner_lattice(std_cell);
        for (i, position) in std_cell.positions().iter().enumerate() {
            let orbit = mapping[i];
            if representative_wyckoffs[orbit].is_some() {
                continue;
            }
            if let Ok(wyckoff) = assign_layer_wyckoff_position(
                position,
                multiplicities[orbit],
                hall_number,
                &lattice,
                symprec,
            ) {
                representative_wyckoffs[orbit] = Some(wyckoff);
            }
        }

        for (i, wyckoff) in representative_wyckoffs.iter().enumerate() {
            if wyckoff.is_none() {
                debug!(
                    "Failed to assign layer Wyckoff position with multiplicity {} (representative site {})",
                    multiplicities[i], remapping[i]
                );
            }
        }
        let representative_wyckoffs = representative_wyckoffs
            .into_iter()
            .map(|w| w.ok_or(MoyoError::WyckoffPositionAssignmentError))
            .collect::<Result<Vec<_>, _>>()?;

        let wyckoffs = (0..std_cell.num_atoms())
            .map(|i| representative_wyckoffs[mapping[orbits[i]]].clone())
            .collect::<Vec<_>>();
        Ok(wyckoffs)
    }
}

fn layer_lattice_system_for(
    arithmetic_number: LayerArithmeticNumber,
) -> Option<LayerLatticeSystem> {
    layer_arithmetic_crystal_class_entry(arithmetic_number).map(|e| e.layer_lattice_system())
}

/// 2D-Minkowski-reduce the in-plane block *after* applying
/// `transformation_to_prim_std`. Mirrors `standardize_triclinic_cell` for the
/// bulk case, except the third basis vector is left untouched.
fn standardize_oblique_layer_cell(
    prim_layer_cell: &LayerCell,
    transformation_to_prim_std: &UnimodularTransformation,
) -> UnimodularTransformation {
    let bulk_lattice = Lattice {
        basis: *prim_layer_cell.lattice().basis(),
    };
    let lattice_after = transformation_to_prim_std.transform_lattice(&bulk_lattice);
    let inplane = Lattice2D::from_inplane_of(&lattice_after.basis);
    let (_, trans_mat_2d) = match inplane.minkowski_reduce() {
        Ok(t) => t,
        // 2D Minkowski reduction is total over non-degenerate inputs; if it
        // fails, fall back to identity rather than aborting standardization.
        Err(_) => return transformation_to_prim_std.clone(),
    };
    let lifted: Linear = lift_2d_to_3d(&trans_mat_2d);
    UnimodularTransformation::new(
        lifted * transformation_to_prim_std.linear,
        transformation_to_prim_std.origin_shift,
    )
}

fn inner_lattice(cell: &LayerCell) -> Lattice {
    Lattice {
        basis: *cell.lattice().basis(),
    }
}

fn transform_layer_cell_unimodular(
    cell: &LayerCell,
    transformation: &UnimodularTransformation,
) -> LayerCell {
    let bulk = Cell::new(
        inner_lattice(cell),
        cell.positions().to_vec(),
        cell.numbers().to_vec(),
    );
    let new_bulk = transformation.transform_cell(&bulk);
    LayerCell::new_unchecked(
        LayerLattice::new_unchecked(new_bulk.lattice),
        new_bulk.positions,
        new_bulk.numbers,
    )
}

fn transform_layer_cell_centering(
    cell: &LayerCell,
    transformation: &Transformation,
) -> (LayerCell, Vec<usize>) {
    let bulk = Cell::new(
        inner_lattice(cell),
        cell.positions().to_vec(),
        cell.numbers().to_vec(),
    );
    let (new_bulk, site_mapping) = transformation.transform_cell(&bulk);
    let new_layer = LayerCell::new_unchecked(
        LayerLattice::new_unchecked(new_bulk.lattice),
        new_bulk.positions,
        new_bulk.numbers,
    );
    (new_layer, site_mapping)
}

/// Symmetrize the in-plane block of a layer lattice with `c` along z.
///
/// The averaged metric tensor of layer-block rotations decouples the
/// in-plane block from the aperiodic axis (`g_13 = g_23 = 0`). Cholesky on
/// the in-plane 2x2 block then produces the layer-canonical orientation
/// `a` along x, `b` in the xy-plane, with `|c|` and the c-axis direction
/// preserved up to the rigid rotation returned alongside.
fn symmetrize_layer_lattice(lattice: &Lattice, rotations: &Rotations) -> (Lattice, Matrix3<f64>) {
    let metric_tensor = lattice.metric_tensor();
    let mut sym_metric: Matrix3<f64> = rotations
        .iter()
        .map(|rotation| {
            rotation.transpose().map(|e| e as f64) * metric_tensor * rotation.map(|e| e as f64)
        })
        .sum();
    sym_metric /= rotations.len() as f64;

    // Force exact zeros on the aperiodic-axis off-diagonal block. The
    // averaging above already drives them to zero up to floating-point
    // noise; pinning them prevents Cholesky on the 3x3 from leaking
    // residual c-y / c-x components into the standardized basis.
    sym_metric[(0, 2)] = 0.0;
    sym_metric[(2, 0)] = 0.0;
    sym_metric[(1, 2)] = 0.0;
    sym_metric[(2, 1)] = 0.0;

    let g11 = sym_metric[(0, 0)];
    let g22 = sym_metric[(1, 1)];
    let g33 = sym_metric[(2, 2)];
    let g12 = sym_metric[(0, 1)];

    let ax = g11.sqrt();
    let bx = if ax.abs() > EPS { g12 / ax } else { 0.0 };
    let by = (g22 - bx * bx).max(0.0).sqrt();
    let cz = g33.sqrt();

    // Layer-canonical basis (column-vector convention):
    //   col 0 = a = (ax, 0, 0)   along x
    //   col 1 = b = (bx, by, 0)  in xy
    //   col 2 = c = (0,  0,  cz) along z
    let new_basis = Matrix3::new(
        ax, bx, 0.0, //
        0.0, by, 0.0, //
        0.0, 0.0, cz,
    );

    // Rotation = new_basis * old_basis^-1, projected to the orthogonal group.
    let mut rotation_matrix = QR::new(new_basis * lattice.basis.try_inverse().unwrap()).q();
    if rotation_matrix.determinant() < 0.0 {
        rotation_matrix *= -1.0;
    }

    (Lattice { basis: new_basis }, rotation_matrix)
}

/// Symmetrize fractional positions by their LG site-symmetry orbits. Same
/// formula as the bulk `symmetrize_positions` -- the layer pipeline just
/// passes layer operations.
fn symmetrize_positions(
    cell: &LayerCell,
    operations: &Operations,
    permutations: &[Permutation],
    epsilon: f64,
) -> Vec<Position> {
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();

    (0..cell.num_atoms())
        .map(|i| {
            let mut acc = Vector3::zeros();
            for (inv_perm, operation) in inverse_permutations.iter().zip(operations.iter()) {
                let mut frac_displacements = operation.rotation.map(|e| e as f64)
                    * cell.positions()[inv_perm.apply(i)]
                    + operation.translation
                    - cell.positions()[i];
                frac_displacements -= frac_displacements.map(|e| e.round());
                acc += frac_displacements;
            }
            acc /= permutations.len() as f64;
            if acc.abs().max() > epsilon {
                warn!(
                    "Large displacement during layer symmetrization: {:?} for site {}",
                    acc, i
                );
            }
            cell.positions()[i] + acc
        })
        .collect::<Vec<_>>()
}

/// Assign a layer Wyckoff position by solving for a compatible integer
/// offset. Mirrors the bulk `assign_wyckoff_position`; the only difference
/// is the database it consults (`iter_layer_wyckoff_positions`).
fn assign_layer_wyckoff_position(
    position: &Position,
    multiplicity: usize,
    hall_number: LayerHallNumber,
    lattice: &Lattice,
    symprec: f64,
) -> Result<LayerWyckoffPosition, MoyoError> {
    use itertools::iproduct;

    for wyckoff in iter_layer_wyckoff_positions(hall_number, multiplicity) {
        let space = WyckoffPositionSpace::new(wyckoff.coordinates);
        let snf = SNF::new(&space.linear);

        let iter_multi_1 = iproduct!(-1..=1, -1..=1, -1..=1);
        let iter_multi_2 = iproduct!(-2_i32..=2_i32, -2_i32..=2_i32, -2_i32..=2_i32)
            .filter(|&(n1, n2, n3)| n1.abs() == 2 || n2.abs() == 2 || n3.abs() == 2);

        for offset in iter_multi_1.chain(iter_multi_2) {
            let offset = Vector3::new(offset.0 as f64, offset.1 as f64, offset.2 as f64);
            let b = snf.l.map(|e| e as f64) * (offset + position - space.origin);
            let mut rinvy = Vector3::zeros();
            for i in 0..3 {
                if snf.d[(i, i)] != 0 {
                    rinvy[i] = b[i] / snf.d[(i, i)] as f64;
                }
            }

            let y = snf.r.map(|e| e as f64) * rinvy;
            let diff = space.linear.map(|e| e as f64) * y + space.origin - position - offset;
            if lattice.cartesian_coords(&diff).norm() < symprec {
                return Ok(wyckoff.clone());
            }
        }
    }
    Err(MoyoError::WyckoffPositionAssignmentError)
}

#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, matrix, vector};

    use super::*;
    use crate::base::{AngleTolerance, Cell, Lattice, traverse};
    use crate::data::{GeometricCrystalClass, LayerSetting, PointGroupRepresentative};
    use crate::search::{LayerPrimitiveCell, LayerPrimitiveSymmetrySearch};

    const SYMPREC: f64 = 1e-4;

    /// Build a layer cell, run the M2 + M3 + M4 pipeline, and return the
    /// resulting `StandardizedLayerCell` plus the identified `LayerGroup`.
    /// Tests below assert on the standardized output.
    fn run_layer_pipeline(
        cell: Cell,
        setting: LayerSetting,
    ) -> (StandardizedLayerCell, LayerGroup) {
        let layer = LayerCell::new(cell, SYMPREC, AngleTolerance::Default).unwrap();
        let primitive = LayerPrimitiveCell::new(&layer, SYMPREC).unwrap();
        let symmetry = LayerPrimitiveSymmetrySearch::new(
            &primitive.layer_cell,
            SYMPREC,
            AngleTolerance::Default,
        )
        .unwrap();
        let nc = primitive.layer_cell.lattice().basis().column(2).norm();
        let volume = primitive.layer_cell.lattice().basis().determinant().abs();
        let epsilon = SYMPREC / volume.powf(1.0 / 3.0);
        let _ = nc;
        let layer_group =
            LayerGroup::new(&symmetry.operations, setting, epsilon).expect("identification failed");
        let standardized = StandardizedLayerCell::new(
            &primitive.layer_cell,
            &symmetry.operations,
            &symmetry.permutations,
            &layer_group,
            SYMPREC,
            epsilon,
            true,
        )
        .expect("standardization failed");
        (standardized, layer_group)
    }

    /// LG 1 (p1): triclinic-oblique in-plane lattice with two generic atoms
    /// kills every higher-symmetry candidate (mirrors the LG 1 fixture in
    /// `tests/test_layer_symmetry.rs`). Standardized cell preserves the
    /// c-axis length and the two Wyckoff letters match (both 'a').
    #[test]
    fn test_layer_p1_standardize_round_trip() {
        let gamma = 80.0_f64.to_radians();
        let a = 1.0;
        let b = 1.5;
        let cell = Cell::new(
            Lattice::new(matrix![
                a, b * gamma.cos(), 0.0;
                0.0, b * gamma.sin(), 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.1, 0.2, 0.1), Vector3::new(0.3, 0.5, 0.2)],
            vec![1, 1],
        );
        let (std, lg) = run_layer_pipeline(cell, LayerSetting::Standard);
        assert_eq!(lg.number, 1);
        // |c| preserved.
        assert_relative_eq!(
            std.layer_cell.lattice().basis().column(2).norm(),
            5.0,
            epsilon = 1e-10
        );
        // c along z, a along x.
        assert_relative_eq!(
            std.layer_cell.lattice().basis()[(2, 2)],
            5.0,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            std.layer_cell.lattice().basis()[(0, 2)],
            0.0,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            std.layer_cell.lattice().basis()[(1, 2)],
            0.0,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            std.layer_cell.lattice().basis()[(1, 0)],
            0.0,
            epsilon = 1e-10
        );
        // LG 1 has only Wyckoff 'a' with multiplicity 1, site symmetry "1".
        assert_eq!(std.wyckoffs.len(), 2);
        for w in std.wyckoffs.iter() {
            assert_eq!(w.letter, 'a');
            assert_eq!(w.multiplicity, 1);
            assert_eq!(w.site_symmetry, "1");
        }
    }

    /// LG 3 (p112): primitive monoclinic-oblique with 2-fold along c. An
    /// atom on the c-axis with z != 0, ±1/2 sits at Wyckoff 'a'
    /// (multiplicity 1, site symmetry "2"). Choosing z != 0 avoids picking
    /// up an extra m_z that would promote LG 3 to LG 6 (p112/m).
    #[test]
    fn test_layer_p112_standardize_origin_atom() {
        let gamma = 80.0_f64.to_radians();
        let a = 1.0;
        let b = 1.5;
        let cell = Cell::new(
            Lattice::new(matrix![
                a, b * gamma.cos(), 0.0;
                0.0, b * gamma.sin(), 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.13)],
            vec![1],
        );
        let (std, lg) = run_layer_pipeline(cell, LayerSetting::Standard);
        assert_eq!(lg.number, 3);
        assert_eq!(std.wyckoffs.len(), 1);
        assert_eq!(std.wyckoffs[0].letter, 'a');
        assert_eq!(std.wyckoffs[0].site_symmetry, "2");
    }

    /// LG 55 (p4mm): single atom at the 4mm site (origin) with z != 0, ±1/2.
    /// The standardized cell must be square (a = b) with c along z preserved.
    /// The atom lands on Wyckoff 'a' (site symmetry "4mm").
    #[test]
    fn test_layer_p4mm_standardize_high_symmetry_site() {
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.1)],
            vec![1],
        );
        let (std, lg) = run_layer_pipeline(cell, LayerSetting::Standard);
        assert_eq!(lg.number, 55);
        // Square: a = b.
        let a = std.layer_cell.lattice().basis().column(0).norm();
        let b = std.layer_cell.lattice().basis().column(1).norm();
        assert_relative_eq!(a, b, epsilon = 1e-10);
        assert_relative_eq!(
            std.layer_cell.lattice().basis().column(2).norm(),
            5.0,
            epsilon = 1e-10
        );
        // High-symmetry site: site_symmetry contains "4mm".
        assert_eq!(std.wyckoffs.len(), 1);
        assert_eq!(std.wyckoffs[0].letter, 'a');
        assert!(std.wyckoffs[0].site_symmetry.contains("4mm"));
    }

    /// In-plane skew should be removed by the standardize pass: feeding a
    /// skewed-but-equivalent oblique cell yields a basis whose
    /// `c`-component on `a` and `b` is zero (layer-block form).
    #[test]
    fn test_standardize_orthogonalizes_inplane_block_against_c() {
        // Same physical structure as `test_layer_p1_standardize_round_trip`
        // but with the in-plane basis pre-rotated within the xy-plane: the
        // standardized output must still satisfy a_z = b_z = 0.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.3, 0.0;
                0.0, 1.05, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.2, 0.3, 0.1)],
            vec![1],
        );
        let (std, _) = run_layer_pipeline(cell, LayerSetting::Standard);
        let basis = std.layer_cell.lattice().basis();
        assert_relative_eq!(basis[(2, 0)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(basis[(2, 1)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(basis[(0, 2)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(basis[(1, 2)], 0.0, epsilon = 1e-10);
    }

    /// Two equivalent atoms (related by 2-fold along c) collapse to one
    /// orbit and share a Wyckoff letter.
    #[test]
    fn test_layer_p112_two_atoms_share_wyckoff() {
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.2, 0.0;
                0.0, 1.3, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![
                Vector3::new(0.31, 0.42, 0.07),
                Vector3::new(-0.31, -0.42, 0.07),
            ],
            vec![1, 1],
        );
        let (std, lg) = run_layer_pipeline(cell, LayerSetting::Standard);
        assert_eq!(lg.number, 3);
        assert_eq!(std.layer_cell.num_atoms(), 2);
        assert_eq!(std.wyckoffs.len(), 2);
        // Both atoms in the same orbit (one Wyckoff letter, multiplicity 2).
        assert_eq!(std.wyckoffs[0].letter, std.wyckoffs[1].letter);
        assert_eq!(std.wyckoffs[0].multiplicity, 2);
        let _ = vector![0.0_f64, 0.0, 0.0];
    }

    #[test]
    fn test_symmetrize_layer_lattice_square() {
        // Square lattice with small distortion should snap to a perfect square
        // (a = b, gamma = 90 deg) with a along x and c along z.
        let lattice = Lattice::new(matrix![
            1.001, 0.0, 0.0;
            0.0, 0.999, 0.0;
            0.0, 0.0, 7.0;
        ]);
        let rep =
            PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::C4v);
        let rotations = traverse(&rep.generators);

        let (new_lattice, rotation_matrix) = symmetrize_layer_lattice(&lattice, &rotations);
        // a = b, in-plane rectangular, c along z, |c| preserved.
        assert_relative_eq!(new_lattice.basis[(0, 0)], new_lattice.basis[(1, 1)]);
        assert_relative_eq!(new_lattice.basis[(2, 2)], 7.0, epsilon = 1e-12);
        assert_relative_eq!(new_lattice.basis[(0, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 1)], 0.0);
        // a along x.
        assert_relative_eq!(new_lattice.basis[(1, 0)], 0.0, epsilon = 1e-12);
        // Rotation matrix is essentially identity (diagonal-only distortion).
        assert_relative_eq!(rotation_matrix.determinant(), 1.0, epsilon = 1e-8);
    }

    #[test]
    fn test_symmetrize_layer_lattice_hexagonal_preserves_c() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            -0.5, (3.0_f64).sqrt() / 2.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        let rep =
            PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::C6v);
        let rotations = traverse(&rep.generators);

        let (new_lattice, _) = symmetrize_layer_lattice(&lattice, &rotations);
        // |c| preserved.
        assert_relative_eq!(new_lattice.basis[(2, 2)], 5.0, epsilon = 1e-12);
        // |a| = |b|, gamma = 120 deg => g12 = -|a|^2/2.
        let a = new_lattice.basis[(0, 0)];
        assert!(a > 0.0);
        let bx = new_lattice.basis[(0, 1)];
        let by = new_lattice.basis[(1, 1)];
        assert_relative_eq!((bx * bx + by * by).sqrt(), a, epsilon = 1e-10);
        assert_relative_eq!(bx / a, -0.5, epsilon = 1e-10);
    }
}
