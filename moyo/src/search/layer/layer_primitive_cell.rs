use log::debug;

use super::super::primitive_cell::{
    compose_input_to_reduced_prim_linear, primitive_cell_from_pure_translations,
    search_pure_translations,
};
use crate::base::{
    Cell, Lattice, Lattice2D, LayerCell, LayerLattice, Linear, MoyoError, Permutation, Translation,
    UnimodularTransformation,
};

/// 2D-Minkowski-reduce `cell`'s in-plane block, leaving `c` untouched.
/// Returns the reduced cell and the lifted 3D unimodular `T` such that
/// `cell.lattice.basis * T == reduced.lattice.basis`. The bulk's 3D Minkowski
/// reduction would mix `c` into a/b and break the aperiodic-axis convention.
fn minkowski_reduce_inplane(cell: &Cell) -> Result<(Cell, Linear), MoyoError> {
    let trans_mat = Lattice2D::lift_inplane_minkowski_reduce(&cell.lattice.basis)?;
    let reduced = UnimodularTransformation::from_linear(trans_mat).transform_cell(cell);
    Ok((reduced, trans_mat))
}

/// Result of a 2D primitive cell search for a layer system.
/// Mirrors `PrimitiveCell` but the discovered translations are constrained to
/// the in-plane sublattice (no `c`-direction translations).
#[derive(Debug)]
pub(crate) struct LayerPrimitiveCell {
    /// Primitive layer cell whose third basis vector equals the input `c`.
    pub layer_cell: LayerCell,
    /// Transformation matrix from the **primitive** cell to the input cell.
    pub linear: Linear,
    /// Mapping from sites of the input cell to those of the primitive cell (many-to-one).
    pub site_mapping: Vec<usize>,
    /// Pure translations in the **input** cell (all have zero `c`-component).
    pub translations: Vec<Translation>,
    /// Permutations induced by the translations. Mirrors `PrimitiveCell` for
    /// API parity; the layer pipeline uses `LayerPrimitiveSymmetrySearch`'s
    /// own permutations downstream rather than this field.
    #[allow(dead_code)]
    pub permutations: Vec<Permutation>,
}

impl LayerPrimitiveCell {
    /// Find the primitive cell of a 2D-periodic (layer) system.
    ///
    /// In-plane axes are 2D-Minkowski-reduced before the kd-tree-based
    /// translation search (the bulk `PeriodicKdTree` requires a Minkowski-
    /// reduced basis for its `[-1, 1]^3` periodic-image enumeration to be
    /// exact). The aperiodic `c` axis is left untouched, so the layer
    /// contract is preserved. Candidate translations whose fractional
    /// `c`-component is non-zero (within `symprec / |c|`) are discarded
    /// after the search -- the layer pipeline only consumes in-plane
    /// translations. The `LayerLattice::new` constructor enforces the
    /// non-degenerate / orthogonal-`c` invariants up front, so basis norms
    /// are not re-checked here.
    pub fn new(layer_cell: &LayerCell, symprec: f64) -> Result<Self, MoyoError> {
        // Reconstruct a bulk `Cell` once: explicit at the call site so the
        // layer-to-bulk crossing is visible.
        let owned_cell = Cell::new(
            Lattice {
                basis: *layer_cell.lattice().basis(),
            },
            layer_cell.positions().to_vec(),
            layer_cell.numbers().to_vec(),
        );

        let (reduced_cell, reduced_trans_mat) = minkowski_reduce_inplane(&owned_cell)?;

        // Sanity-check tolerance against the *reduced* in-plane vectors.
        let reduced_basis = &reduced_cell.lattice.basis;
        let nc = reduced_basis.column(2).norm();
        let min_inplane_norm = reduced_basis
            .column(0)
            .norm()
            .min(reduced_basis.column(1).norm());
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > min_inplane_norm / 2.0 {
            debug!(
                "symprec is too large compared to the in-plane basis vectors. Consider reducing symprec."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        // Reuse the bulk pure-translation search (kd-tree + correspondence +
        // symmetrize), which is correct on a Minkowski-reduced basis.
        let (raw_translations, raw_permutations) = search_pure_translations(&reduced_cell, symprec);

        // Drop translations with non-zero `c`-component: they are not
        // in-plane lattice translations, so they do not belong to the layer
        // group. Snap the surviving ones to `tz = 0`.
        let z_tol = symprec / nc;
        let mut translations = vec![];
        let mut permutations = vec![];
        for (mut translation, permutation) in raw_translations
            .into_iter()
            .zip(raw_permutations.into_iter())
        {
            let mut tz = translation[2];
            tz -= tz.round();
            if tz.abs() > z_tol {
                debug!(
                    "Skipping translation with non-zero c-component: tz={:.6} (tol={:.6})",
                    tz, z_tol
                );
                continue;
            }
            translation[2] = 0.0;
            translations.push(translation);
            permutations.push(permutation);
        }

        // Layer translations have `tz = 0`, so the shared post-search core
        // (size check + 3D HNF + primitive-cell build) yields a `trans_mat`
        // in layer-block form (`W_33 = 1`, `W_3i = W_i3 = 0`) automatically,
        // and the primitive-cell builder preserves the `c` column of the
        // basis and the `z` of fractional positions by construction.
        let (primitive_cell, trans_mat, site_mapping) =
            primitive_cell_from_pure_translations(&reduced_cell, &translations, &permutations)?;

        // Re-reduce so the returned layer cell satisfies the Minkowski-reduced
        // precondition of downstream bulk routines (`PrimitiveSymmetrySearch`,
        // `PeriodicKdTree`, `search_bravais_group`). The HNF-derived primitive
        // basis is upper-triangular, not Minkowski-reduced, so this step is
        // required even though the input was already reduced.
        let (reduced_prim_cell, prim_trans_mat) = minkowski_reduce_inplane(&primitive_cell)?;

        // (input cell)
        //    -[reduced_trans_mat]-> (reduced cell)
        //    <-[trans_mat]- (primitive cell)
        //    -[prim_trans_mat]-> (reduced primitive cell)
        let linear: Linear =
            compose_input_to_reduced_prim_linear(reduced_trans_mat, trans_mat, prim_trans_mat);

        let translations_in_input = translations
            .iter()
            .map(|t| reduced_trans_mat.map(|e| e as f64) * t)
            .collect::<Vec<_>>();

        let Cell {
            lattice: prim_lattice,
            positions: prim_positions,
            numbers: prim_numbers,
        } = reduced_prim_cell;
        Ok(Self {
            layer_cell: LayerCell::new_unchecked(
                LayerLattice::new_unchecked(prim_lattice),
                prim_positions,
                prim_numbers,
            ),
            linear,
            site_mapping,
            translations: translations_in_input,
            permutations,
        })
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, matrix};

    use super::LayerPrimitiveCell;
    use crate::base::{AngleTolerance, Cell, Lattice, LayerCell};

    fn make_layer(cell: Cell) -> LayerCell {
        LayerCell::new(cell, 1e-4, AngleTolerance::Default).unwrap()
    }

    #[test]
    fn test_layer_p1_single_atom() {
        // p1 square lattice, one atom: no contraction expected.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.3, 0.4, 0.2)],
            vec![1],
        );
        let input_basis = cell.lattice.basis;
        let layer = make_layer(cell);
        let result = LayerPrimitiveCell::new(&layer, 1e-4).unwrap();
        assert_eq!(result.layer_cell.num_atoms(), 1);
        assert_eq!(result.translations.len(), 1);
        assert_relative_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        // Lattice is unchanged (input was already 2D-Minkowski-reduced).
        assert_relative_eq!(
            *result.layer_cell.lattice().basis(),
            input_basis,
            epsilon = 1e-8
        );
        // The atom z-coordinate carries through unchanged.
        assert_relative_eq!(result.layer_cell.positions()[0][2], 0.2, epsilon = 1e-12);
    }

    #[test]
    fn test_layer_doubled_in_plane() {
        // Two atoms related by a (1/2, 0, 0) in-plane translation. The
        // primitive cell halves along `a`, with `c` left intact.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.1), Vector3::new(0.5, 0.0, 0.1)],
            vec![1, 1],
        );
        let input_c = cell.lattice.basis.column(2).into_owned();
        let layer = make_layer(cell);
        let result = LayerPrimitiveCell::new(&layer, 1e-4).unwrap();
        assert_eq!(result.translations.len(), 2);
        // One of the translations must be (1/2, 0, 0); both must have zero z.
        let has_half = result
            .translations
            .iter()
            .any(|t| (t[0] - 0.5).abs() < 1e-8 && t[1].abs() < 1e-8 && t[2].abs() < 1e-12);
        assert!(has_half);
        for t in result.translations.iter() {
            assert_eq!(t[2], 0.0);
        }
        // Primitive cell halves along `a` and preserves `c`.
        assert_relative_eq!(
            result.layer_cell.lattice().basis().column(0).into_owned(),
            Vector3::new(0.5, 0.0, 0.0),
            epsilon = 1e-8
        );
        assert_relative_eq!(
            result.layer_cell.lattice().basis().column(2).into_owned(),
            input_c,
            epsilon = 1e-12
        );
        assert_eq!(result.layer_cell.num_atoms(), 1);
        assert_relative_eq!(result.layer_cell.positions()[0][2], 0.1, epsilon = 1e-12);
    }

    #[test]
    fn test_layer_along_c_translation_discarded() {
        // Atoms at (0,0,0) and (0,0,1/2): the (0,0,1/2) candidate aligns atoms
        // but is not an in-plane lattice translation, so it is silently discarded.
        // The cell is treated as already primitive with two distinct atoms.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.5)],
            vec![1, 1],
        );
        let input_basis = cell.lattice.basis;
        let layer = make_layer(cell);
        let result = LayerPrimitiveCell::new(&layer, 1e-4).unwrap();
        assert_eq!(result.translations.len(), 1);
        assert_relative_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(result.layer_cell.num_atoms(), 2);
        assert_relative_eq!(
            *result.layer_cell.lattice().basis(),
            input_basis,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_layer_honeycomb_like() {
        // Hexagonal in-plane cell with two atoms at (1/3, 2/3, 1/2) and
        // (2/3, 1/3, 1/2). They are not related by an in-plane lattice
        // translation, so the cell is already primitive.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                -0.5, (3.0_f64).sqrt() / 2.0, 0.0;
                0.0, 0.0, 6.0;
            ]),
            vec![
                Vector3::new(1.0 / 3.0, 2.0 / 3.0, 0.5),
                Vector3::new(2.0 / 3.0, 1.0 / 3.0, 0.5),
            ],
            vec![1, 1],
        );
        let input_basis = cell.lattice.basis;
        let layer = make_layer(cell);
        let result = LayerPrimitiveCell::new(&layer, 1e-4).unwrap();
        assert_eq!(result.translations.len(), 1);
        assert_relative_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(result.layer_cell.num_atoms(), 2);
        assert_relative_eq!(
            *result.layer_cell.lattice().basis(),
            input_basis,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_layer_skewed_in_plane_basis_is_reduced() {
        // Skewed in-plane basis: a = (1, 0, 0), b = (4, 1, 0), c along z.
        // Without 2D Minkowski reduction the kd-tree's [-1, 1]^3 image search
        // can miss the true nearest periodic image. After reduction the
        // primitive cell finder still returns a single trivial translation
        // for a one-atom cell, and the linear field correctly maps the input
        // basis back to the primitive (= input, here) basis.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 4.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.1)],
            vec![1],
        );
        let input_basis = cell.lattice.basis;
        let layer = make_layer(cell);
        let result = LayerPrimitiveCell::new(&layer, 1e-4).unwrap();
        assert_eq!(result.layer_cell.num_atoms(), 1);
        assert_eq!(result.translations.len(), 1);
        // primitive_basis * linear == input_basis (reconstruct the input basis
        // through the reported transform; no claim on whether basis equals input).
        let reconstructed = result.layer_cell.lattice().basis() * result.linear.map(|e| e as f64);
        assert_relative_eq!(reconstructed, input_basis, epsilon = 1e-8);
    }
}
