use std::collections::BTreeMap;

use log::debug;
use nalgebra::{Dyn, Matrix2, Matrix3, OMatrix, U2, Vector2, Vector3};

use super::solve::{
    PeriodicKdTree, pivot_site_indices, solve_correspondence,
    symmetrize_translation_from_permutation,
};
use crate::base::{
    Cell, EPS, Lattice, Linear, MoyoError, Permutation, Position, Rotation, Translation,
    orbits_from_permutations,
};
use crate::math::HNF;

/// Result of a 2D primitive cell search for a layer system.
/// Mirrors `PrimitiveCell` but the discovered translations are constrained to
/// the in-plane sublattice (no `c`-direction translations).
#[derive(Debug)]
#[allow(dead_code)] // consumed by later layer-group milestones
pub(crate) struct LayerPrimitiveCell {
    /// Primitive layer cell whose third basis vector equals the input `c`.
    pub cell: Cell,
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
    /// Per the layer-group plan §4.2, candidate translations are filtered to
    /// those whose fractional `c`-component is zero within `symprec / |c|`.
    /// A non-trivial candidate (one that aligns atoms) with a non-zero
    /// `c`-component falsifies the layer-group hypothesis and yields
    /// `MoyoError::SpuriousAperiodicTranslation`.
    pub fn new(cell: &Cell, symprec: f64) -> Result<Self, MoyoError> {
        // We deliberately skip the 3D Minkowski reduction used by `PrimitiveCell::new`:
        // mixing `c` into the in-plane basis would invalidate the §3.1 axis convention,
        // and the in-plane 2D Minkowski reduction is part of the M4 standardization
        // pipeline, not the search step.
        let basis = &cell.lattice.basis;
        let na = basis.column(0).norm();
        let nb = basis.column(1).norm();
        let nc = basis.column(2).norm();
        if na == 0.0 || nb == 0.0 || nc == 0.0 {
            return Err(MoyoError::PrimitiveCellError);
        }

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

            // Wrap the c-component into [-0.5, 0.5] before testing the layer constraint.
            let mut tz = translation[2];
            tz -= tz.round();

            let new_positions: Vec<Position> =
                cell.positions.iter().map(|pos| pos + translation).collect();

            if let Some(permutation) = solve_correspondence(&pkdtree, cell, &new_positions) {
                if tz.abs() > z_tol {
                    // The candidate aligns atoms but has a non-zero c-component:
                    // this means the input has a true lattice translation along the
                    // aperiodic axis, contradicting the layer-group hypothesis.
                    debug!(
                        "Spurious aperiodic translation found: tz={:.6} (tol={:.6})",
                        tz, z_tol
                    );
                    return Err(MoyoError::SpuriousAperiodicTranslation);
                }
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

        // Build the 2D transformation matrix from primitive (a', b') to input (a, b),
        // then lift to 3D with the c-axis fixed.
        let trans_mat_2d =
            if let Some(m) = transformation_matrix_from_inplane_translations(&translations) {
                m
            } else {
                debug!("Failed to find a 2D transformation matrix for the primitive layer cell.");
                return Err(MoyoError::TooSmallToleranceError);
            };
        let trans_mat: Linear = lift_2d_to_3d(&trans_mat_2d);

        let (primitive_cell, site_mapping, _) = primitive_layer_cell_from_transformation(
            cell,
            &trans_mat,
            &translations,
            &permutations,
        );

        // The relation `cell.lattice.basis * trans_mat^{-1}.cols * size = primitive.basis`
        // (modulo the in-plane block) is satisfied by construction; no further reduction
        // is performed here -- standardization (M4) handles 2D Minkowski reduction.
        Ok(Self {
            cell: primitive_cell,
            linear: trans_mat,
            site_mapping,
            translations,
            permutations,
        })
    }
}

/// Build the 2D HNF transformation matrix from a set of in-plane translations.
/// The c-component of each translation is assumed to already be zero (or
/// negligible).
fn transformation_matrix_from_inplane_translations(
    translations: &[Translation],
) -> Option<Matrix2<i32>> {
    let size = translations.len() as i32;
    let mut columns: Vec<Vector2<i32>> = vec![Vector2::new(size, 0), Vector2::new(0, size)];
    for translation in translations.iter() {
        let scaled = Vector2::new(
            (translation[0] * size as f64).round() as i32,
            (translation[1] * size as f64).round() as i32,
        );
        columns.push(scaled);
    }
    let basis = OMatrix::<i32, U2, Dyn>::from_columns(&columns);
    let hnf = HNF::new(&basis);
    let trans_mat_inv = Matrix2::<i32>::from_columns(&[hnf.h.column(0), hnf.h.column(1)])
        .map(|e| e as f64)
        / (size as f64);
    let trans_mat = trans_mat_inv.try_inverse()?.map(|e| e.round() as i32);

    if relative_ne!(
        trans_mat.map(|e| e as f64).determinant(),
        size as f64,
        epsilon = EPS
    ) {
        return None;
    }
    Some(trans_mat)
}

/// Lift a 2x2 in-plane integer transformation into the 3x3 layer-group block
/// form `[[A, 0], [0, 1]]` (paper eq. 4 with `W_33 = 1`).
fn lift_2d_to_3d(m: &Matrix2<i32>) -> Linear {
    Matrix3::<i32>::new(
        m[(0, 0)],
        m[(0, 1)],
        0, //
        m[(1, 0)],
        m[(1, 1)],
        0, //
        0,
        0,
        1,
    )
}

/// Transform `cell` to a primitive layer cell by inverse of `trans_mat`,
/// preserving the input `c` axis and the fractional `z` coordinates.
fn primitive_layer_cell_from_transformation(
    cell: &Cell,
    trans_mat: &Linear,
    translations: &[Translation],
    permutations: &[Permutation],
) -> (Cell, Vec<usize>, Vec<usize>) {
    // Build the new lattice with the in-plane block transformed and `c` carried through.
    let new_basis_inplane = cell.lattice.basis * trans_mat.map(|e| e as f64).try_inverse().unwrap();
    let new_basis = Matrix3::from_columns(&[
        new_basis_inplane.column(0),
        new_basis_inplane.column(1),
        cell.lattice.basis.column(2),
    ]);
    // The c-column inferred from the transform agrees with the input c by construction
    // because `trans_mat` has the layer-block form (paper eq. 4 with `W_33 = 1`).
    let new_lattice = Lattice { basis: new_basis };

    let num_atoms = cell.num_atoms();
    let orbits = orbits_from_permutations(num_atoms, permutations);
    let representatives = (0..num_atoms)
        .filter(|&i| orbits[i] == i)
        .collect::<Vec<_>>();

    let mut new_positions = vec![Vector3::zeros(); representatives.len()];
    let mut new_numbers = vec![0; representatives.len()];
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();
    // Eq. (25) of https://arxiv.org/pdf/2211.15008.pdf, with `c`-axis pinned: the
    // averaged in-plane shift is transformed but `z` stays equal to the input.
    for (i, &orbit_i) in representatives.iter().enumerate() {
        let mut acc = Vector3::zeros();
        for (inv_perm, translation) in inverse_permutations.iter().zip(translations.iter()) {
            let mut frac_displacements =
                cell.positions[inv_perm.apply(orbit_i)] + translation - cell.positions[orbit_i];
            frac_displacements -= frac_displacements.map(|e| e.round());
            acc += frac_displacements;
        }
        let averaged = cell.positions[orbit_i] + acc / (translations.len() as f64);
        let mut new_pos = trans_mat.map(|e| e as f64) * averaged;
        // Force the z-component to equal the input z exactly (no rescaling along c).
        new_pos[2] = cell.positions[orbit_i][2];
        new_positions[i] = new_pos;
        new_numbers[i] = cell.numbers[orbit_i];
    }

    let primitive_cell = Cell::new(new_lattice, new_positions, new_numbers);
    let site_mapping = site_mapping_from_orbits(&orbits);
    (primitive_cell, site_mapping, representatives)
}

fn site_mapping_from_orbits(orbits: &[usize]) -> Vec<usize> {
    let mut mapping = BTreeMap::new();
    let mut count = 0;
    for ri in orbits.iter() {
        mapping.entry(ri).or_insert_with(|| {
            let value = count;
            count += 1;
            value
        });
    }
    orbits.iter().map(|ri| *mapping.get(&ri).unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, matrix};

    use super::LayerPrimitiveCell;
    use crate::base::{Cell, Lattice, MoyoError};

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
        let result = LayerPrimitiveCell::new(&cell, 1e-4).unwrap();
        assert_eq!(result.cell.num_atoms(), 1);
        assert_eq!(result.translations.len(), 1);
        assert_relative_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        // Lattice is unchanged.
        assert_relative_eq!(
            result.cell.lattice.basis,
            cell.lattice.basis,
            epsilon = 1e-8
        );
        // The atom z-coordinate carries through unchanged.
        assert_relative_eq!(result.cell.positions[0][2], 0.2, epsilon = 1e-12);
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
        let result = LayerPrimitiveCell::new(&cell, 1e-4).unwrap();
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
            result.cell.lattice.basis.column(0).into_owned(),
            Vector3::new(0.5, 0.0, 0.0),
            epsilon = 1e-8
        );
        assert_relative_eq!(
            result.cell.lattice.basis.column(2).into_owned(),
            cell.lattice.basis.column(2).into_owned(),
            epsilon = 1e-12
        );
        assert_eq!(result.cell.num_atoms(), 1);
        assert_relative_eq!(result.cell.positions[0][2], 0.1, epsilon = 1e-12);
    }

    #[test]
    fn test_layer_spurious_stacking_rejected() {
        // Atoms at (0,0,0) and (0,0,1/2): a (0,0,1/2) translation aligns atoms,
        // so the layer-group hypothesis is falsified.
        let cell = Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.5)],
            vec![1, 1],
        );
        let err = LayerPrimitiveCell::new(&cell, 1e-4).unwrap_err();
        assert_eq!(err, MoyoError::SpuriousAperiodicTranslation);
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
        let result = LayerPrimitiveCell::new(&cell, 1e-4).unwrap();
        assert_eq!(result.translations.len(), 1);
        assert_relative_eq!(result.translations[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(result.cell.num_atoms(), 2);
        assert_relative_eq!(
            result.cell.lattice.basis,
            cell.lattice.basis,
            epsilon = 1e-12
        );
    }
}
