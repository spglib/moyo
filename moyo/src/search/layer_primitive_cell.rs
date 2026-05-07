use log::debug;

use super::primitive_cell::{
    primitive_cell_from_transformation, transformation_matrix_from_translations,
};
use super::solve::{
    PeriodicKdTree, pivot_site_indices, solve_correspondence,
    symmetrize_translation_from_permutation,
};
use crate::base::{
    Cell, Lattice, LayerCell, LayerLattice, Linear, MoyoError, Permutation, Position, Rotation,
    Translation,
};

/// Result of a 2D primitive cell search for a layer system.
/// Mirrors `PrimitiveCell` but the discovered translations are constrained to
/// the in-plane sublattice (no `c`-direction translations).
#[derive(Debug)]
#[allow(dead_code)] // consumed by later layer-group milestones
pub(crate) struct LayerPrimitiveCell {
    /// Primitive layer cell whose third basis vector equals the input `c`.
    pub layer_cell: LayerCell,
    /// Transformation matrix from the **primitive** cell to the input cell.
    pub linear: Linear,
    /// Mapping from sites of the input cell to those of the primitive cell (many-to-one).
    pub site_mapping: Vec<usize>,
    /// Pure translations in the **input** cell (all have zero `c`-component).
    pub translations: Vec<Translation>,
    /// Permutations induced by the translations.
    pub permutations: Vec<Permutation>,
}

#[allow(dead_code)] // consumed by later layer-group milestones
impl LayerPrimitiveCell {
    /// Find the primitive cell of a 2D-periodic (layer) system.
    ///
    /// Candidate translations whose fractional `c`-component is non-zero
    /// (within `symprec / |c|`) are silently discarded -- the layer-group
    /// pipeline only cares about in-plane translations. The lattice contract
    /// (perpendicular, non-degenerate basis) is already enforced by
    /// `LayerLattice::new`, so basis norms are not re-checked here.
    pub fn new(layer_cell: &LayerCell, symprec: f64) -> Result<Self, MoyoError> {
        // Reconstruct a bulk `Cell` once: explicit at the call site so the
        // layer-to-bulk crossing is visible. Downstream helpers
        // (`PeriodicKdTree`, `solve_correspondence`) take `&Cell`.
        let owned_cell = Cell::new(
            Lattice {
                basis: *layer_cell.lattice().basis(),
            },
            layer_cell.positions().to_vec(),
            layer_cell.numbers().to_vec(),
        );
        let cell = &owned_cell;
        // We deliberately skip the 3D Minkowski reduction used by `PrimitiveCell::new`:
        // mixing `c` into the in-plane basis would break the convention that `c` is
        // the aperiodic axis. 2D Minkowski reduction belongs to standardization,
        // not the search step.
        let basis = &cell.lattice.basis;
        let na = basis.column(0).norm();
        let nb = basis.column(1).norm();
        let nc = basis.column(2).norm();

        // Sanity-check tolerance against in-plane vectors only -- `c` is not a
        // candidate lattice translation direction for layer systems.
        let min_inplane_norm = na.min(nb);
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > min_inplane_norm / 2.0 {
            debug!(
                "symprec is too large compared to the in-plane basis vectors. Consider reducing symprec."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        let pkdtree = PeriodicKdTree::new(cell, rough_symprec);
        let pivots = pivot_site_indices(&cell.numbers);
        let src = pivots[0];

        // Tolerance on the fractional `c`-component, derived from the cartesian
        // tolerance `symprec` divided by `|c|`.
        let z_tol = symprec / nc;

        let mut permutations_translations_tmp = vec![];
        for dst in pivots.iter() {
            let translation = cell.positions[*dst] - cell.positions[src];

            // Wrap the c-component into [-0.5, 0.5] and skip candidates with a
            // non-zero c-component before doing the (expensive) correspondence
            // solve -- they cannot be in-plane lattice translations.
            let mut tz = translation[2];
            tz -= tz.round();
            if tz.abs() > z_tol {
                debug!(
                    "Skipping translation with non-zero c-component: tz={:.6} (tol={:.6})",
                    tz, z_tol
                );
                continue;
            }

            let new_positions: Vec<Position> =
                cell.positions.iter().map(|pos| pos + translation).collect();

            if let Some(permutation) = solve_correspondence(&pkdtree, cell, &new_positions) {
                permutations_translations_tmp.push((permutation, translation));
            }
        }

        // Purify translations by permutations and snap the c-component to zero.
        let mut translations = vec![];
        let mut permutations = vec![];
        for (permutation, rough_translation) in permutations_translations_tmp.iter() {
            let (mut translation, distance) = symmetrize_translation_from_permutation(
                cell,
                permutation,
                &Rotation::identity(),
                rough_translation,
            );
            if distance < symprec {
                translation[2] = 0.0;
                translations.push(translation);
                permutations.push(permutation.clone());
            }
        }

        let size = translations.len() as i32;
        if (size == 0) || (cell.num_atoms() % (size as usize) != 0) {
            debug!(
                "Failed to properly find layer translations: {} translations in {} atoms.",
                size,
                cell.num_atoms()
            );
            return Err(MoyoError::TooSmallToleranceError);
        }
        debug!("Found {} pure layer translations", size);

        // Build the transformation matrix from primitive (a', b', c) to input
        // (a, b, c). Reuse the bulk 3D HNF helper: layer translations have
        // `tz = 0` (snapped above), so the third row of the resulting basis is
        // locked to `(0, 0, size)`, giving a `trans_mat` with the layer-group
        // block form (`W_33 = 1`, `W_3i = W_i3 = 0`) automatically.
        let trans_mat =
            transformation_matrix_from_translations(&translations).ok_or_else(|| {
                debug!("Failed to find a transformation matrix for the primitive layer cell.");
                MoyoError::TooSmallToleranceError
            })?;

        // Reuse the bulk primitive-cell builder. `trans_mat` preserves the
        // `c` column of the basis and the `z` component of fractional
        // positions by construction, so the layer contract is preserved.
        let (primitive_cell, site_mapping, _) =
            primitive_cell_from_transformation(cell, &trans_mat, &translations, &permutations);

        let Cell {
            lattice: prim_lattice,
            positions: prim_positions,
            numbers: prim_numbers,
        } = primitive_cell;
        Ok(Self {
            layer_cell: LayerCell::new_unchecked(
                LayerLattice::new_unchecked(prim_lattice),
                prim_positions,
                prim_numbers,
            ),
            linear: trans_mat,
            site_mapping,
            translations,
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
        // Lattice is unchanged.
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
}
