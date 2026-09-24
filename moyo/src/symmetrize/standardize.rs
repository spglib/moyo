use nalgebra::Matrix3;
use std::collections::HashMap;

use super::conventional_coordinate_system::ConventionalCoordinateSystem;
use super::symmetrization::{symmetrize_lattice, symmetrize_positions};
use super::wyckoff::{assign_wyckoffs_by_orbit, group_sites_by_orbit, match_wyckoff_coordinates};
use crate::base::{
    Cell, MoyoError, Operations, Permutation, Transformation, UnimodularTransformation,
    project_rotations,
};
use crate::data::{HallNumber, WyckoffPosition, iter_wyckoff_positions};
use crate::identify::SpaceGroup;

/// Result of the full standardization pipeline for a primitive cell.
///
/// [`ConventionalCoordinateSystem`] selects the basis and origin. This type then
/// refines the lattice and atomic positions, constructs the primitive and conventional cells,
/// applies the requested Cartesian orientation, and assigns Wyckoff positions.
///
/// The [returned-cell specification](https://spglib.github.io/moyo/standardization/#returned-cell-specification)
/// describes the target symmetry, coordinate transformations, and Cartesian orientation.
///
/// Both output cells use the refined metric.
/// Coordinate transformations describe the selected basis and origin before
/// refinement; the rotation matrix describes only the final Cartesian rotation.
pub struct StandardizedCell {
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    /// Primitive output, preserving the input primitive cell's site order and species.
    pub prim_cell: Cell,
    /// Coordinate transformation from the input primitive cell to the selected
    /// primitive system, before lattice and position refinement.
    pub prim_transformation: UnimodularTransformation,
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    /// Conventional output in the selected Hall setting.
    pub cell: Cell,
    /// One Wyckoff position per conventional site, in `cell` order and the selected setting.
    pub wyckoffs: Vec<WyckoffPosition>,
    /// Coordinate transformation from the input primitive cell to the selected
    /// conventional system, before lattice and position refinement.
    pub transformation: Transformation,
    // ------------------------------------------------------------------------
    // Miscellaneous
    // ------------------------------------------------------------------------
    /// Proper Cartesian rotation applied to both output lattices.
    /// Identity when `rotate_basis` is false; it does not encode lattice refinement.
    pub rotation_matrix: Matrix3<f64>,
    /// One primitive site index per conventional site: `site_mapping[j]` gives
    /// the site in `prim_cell` corresponding to site `j` in `cell`.
    pub site_mapping: Vec<usize>,
}

impl StandardizedCell {
    /// Standardize the input **primitive** cell.
    /// See [`Self`] for the returned-cell specification.
    ///
    /// Returns [`MoyoError::StandardizationError`] if lattice refinement fails.
    pub fn new(
        prim_cell: &Cell,
        prim_operations: &Operations,
        prim_permutations: &[Permutation],
        space_group: &SpaceGroup,
        symprec: f64,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<Self, MoyoError> {
        let coordinate_system =
            ConventionalCoordinateSystem::new(&prim_cell.lattice, space_group, epsilon)?;
        let (
            prim_std_cell,
            prim_std_permutations,
            prim_transformation,
            std_cell,
            transformation,
            rotation_matrix,
            site_mapping,
        ) = Self::standardize_and_symmetrize_cell(
            prim_cell,
            prim_operations,
            prim_permutations,
            coordinate_system,
            epsilon,
            rotate_basis,
        )?;

        let wyckoffs = Self::assign_wyckoffs(
            &prim_std_cell,
            &prim_std_permutations,
            &std_cell,
            &site_mapping,
            space_group.hall_number,
            symprec,
        )?;

        Ok(StandardizedCell {
            // Primitive standardized cell
            prim_cell: prim_std_cell,
            prim_transformation,
            // Standardized cell
            cell: std_cell,
            wyckoffs,
            transformation,
            // Miscellaneous
            rotation_matrix,
            site_mapping,
        })
    }

    #[allow(clippy::type_complexity)]
    fn standardize_and_symmetrize_cell(
        prim_cell: &Cell,
        prim_operations: &Operations,
        prim_permutations: &[Permutation],
        coordinate_system: ConventionalCoordinateSystem,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<
        (
            Cell,
            Vec<Permutation>,
            UnimodularTransformation,
            Cell,
            Transformation,
            Matrix3<f64>,
            Vec<usize>,
        ),
        MoyoError,
    > {
        let ConventionalCoordinateSystem {
            prim_transformation,
            conv_trans_linear,
            transformation,
            conv_std_operations,
            prim_std_operations,
        } = coordinate_system;

        let prim_std_cell_tmp = prim_transformation.transform_cell(prim_cell);

        // Symmetrize positions of prim_std_cell by refined symmetry operations.
        // Reorder permutations because prim_std_operations (from the Hall-symbol
        // traversal) is in a different order than `prim_operations` (from the
        // symmetry search).
        let prim_std_permutations = align_primitive_permutations(
            &prim_transformation,
            prim_operations,
            prim_permutations,
            &prim_std_operations,
        )?;
        let new_prim_std_positions = symmetrize_positions(
            &prim_std_cell_tmp.positions,
            &prim_std_operations,
            &prim_std_permutations,
            epsilon,
        );

        let centering = Transformation::from_linear(conv_trans_linear);
        let conv_lattice = centering.transform_lattice(&prim_std_cell_tmp.lattice);
        let (refined_lattice, rotation_matrix) =
            symmetrize_lattice(&conv_lattice, &project_rotations(&conv_std_operations))?;
        let (std_lattice, rotation_matrix) = if rotate_basis {
            (refined_lattice, rotation_matrix)
        } else {
            // Undo only the polar rotation, retaining the symmetric stretch.
            (
                refined_lattice.rotate(&rotation_matrix.transpose()),
                Matrix3::identity(),
            )
        };

        // Note: prim_transformation.transform_cell does not change the order of sites
        let prim_std_cell = Cell::new(
            centering.inverse_transform_lattice(&std_lattice),
            new_prim_std_positions,
            prim_std_cell_tmp.numbers.clone(),
        );

        // To (conventional) standardized cell
        let (std_cell, site_mapping) = centering.transform_cell(&prim_std_cell);

        Ok((
            prim_std_cell,
            prim_std_permutations,
            prim_transformation,
            std_cell,
            transformation,
            rotation_matrix,
            site_mapping,
        ))
    }

    fn assign_wyckoffs(
        prim_std_cell: &Cell,
        prim_std_permutations: &[Permutation],
        std_cell: &Cell,
        site_mapping: &[usize],
        hall_number: HallNumber,
        symprec: f64,
    ) -> Result<Vec<WyckoffPosition>, MoyoError> {
        let group = group_sites_by_orbit(
            prim_std_cell.num_atoms(),
            prim_std_permutations,
            site_mapping,
            std_cell.num_atoms(),
        );
        assign_wyckoffs_by_orbit(&group, &std_cell.positions, |position, multiplicity| {
            iter_wyckoff_positions(hall_number, multiplicity)
                .find(|w| {
                    match_wyckoff_coordinates(position, w.coordinates, &std_cell.lattice, symprec)
                })
                .cloned()
        })
    }
}

/// Align a set of primitive-cell permutations (as produced by the symmetry
/// search) with a target operation list (as produced by traversing a Hall
/// symbol) by matching rotation matrices. Both pipelines need this because
/// the search and the database traversal generate operations in different
/// orders. Errors with `StandardizationError` if any target rotation is
/// missing from the input — i.e. the input set is not closed under the
/// transformation, which would indicate a primitive-cell or
/// hall-symbol-database bug rather than a user error.
pub(super) fn align_primitive_permutations(
    prim_transformation: &UnimodularTransformation,
    prim_operations: &Operations,
    prim_permutations: &[Permutation],
    target_operations: &Operations,
) -> Result<Vec<Permutation>, MoyoError> {
    let mut permutation_mapping = HashMap::new();
    let prim_rotations =
        project_rotations(&prim_transformation.transform_operations(prim_operations));
    for (rot, perm) in prim_rotations.iter().zip(prim_permutations.iter()) {
        permutation_mapping.insert(*rot, perm.clone());
    }
    target_operations
        .iter()
        .map(|ops| {
            permutation_mapping
                .get(&ops.rotation)
                .cloned()
                .ok_or(MoyoError::StandardizationError)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Rotation3, Vector3, matrix, vector};

    use super::StandardizedCell;
    use crate::base::{Cell, Lattice, Operations, Permutation, UnimodularTransformation};
    use crate::data::{HallNumber, HallSymbol, hall_symbol_entry};
    use crate::identify::SpaceGroup;

    fn fixture(hall_number: HallNumber) -> (Cell, Operations, Vec<Permutation>, SpaceGroup) {
        let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        let conventional_operations = hall_symbol.traverse();
        // Analytic metrics, selected by invariance in this Hall setting. Include
        // all unique axes for monoclinic groups and primitive rhombohedral axes.
        let conventional_lattice = [
            matrix![4.0, 0.0, 0.0; 0.7, 5.0, 0.0; 0.2, 0.4, 6.0],
            matrix![4.0, 0.0, 0.0; 0.0, 5.0, 0.0; -1.0, 0.0, 6.0],
            matrix![4.0, 0.0, 0.0; 0.0, 5.0, 0.0; 0.0, -1.0, 6.0],
            matrix![4.0, 0.0, 0.0; -1.0, 5.0, 0.0; 0.0, 0.0, 6.0],
            Matrix3::from_diagonal(&vector![4.0, 5.0, 6.0]),
            Matrix3::from_diagonal(&vector![4.0, 4.0, 6.0]),
            matrix![4.0, 0.0, 0.0; -2.0, 12.0_f64.sqrt(), 0.0; 0.0, 0.0, 6.0],
            matrix![0.0, 2.0, 2.0; 2.0, 0.0, 2.0; 2.0, 2.0, 0.0],
            4.0 * Matrix3::identity(),
        ]
        .into_iter()
        .map(Lattice::new)
        .find(|lattice| {
            let metric = lattice.metric_tensor();
            conventional_operations.iter().all(|operation| {
                let w = operation.rotation.map(f64::from);
                (w.transpose() * metric * w - metric).norm() < 1e-12 * metric.norm()
            })
        })
        .expect("an analytic metric must cover each Hall setting");
        let operations = hall_symbol.primitive_traverse();
        let mut positions = Vec::<Vector3<f64>>::new();
        let mut numbers = vec![];
        for (species, seed) in [(1, vector![0.137, 0.239, 0.371]), (2, Vector3::zeros())] {
            for operation in &operations {
                let position = (operation.rotation.map(f64::from) * seed + operation.translation)
                    .map(|x| x.rem_euclid(1.0));
                if positions.iter().zip(&numbers).any(|(other, &number)| {
                    let delta = position - other;
                    number == species && (delta - delta.map(f64::round)).norm() < 1e-10
                }) {
                    continue;
                }
                positions.push(position);
                numbers.push(species);
            }
        }
        let cell = Cell::new(
            Lattice {
                basis: conventional_lattice.basis
                    * hall_symbol_entry(hall_number).unwrap().centering.inverse(),
            },
            positions,
            numbers,
        );
        let permutations = operations
            .iter()
            .map(|operation| {
                Permutation::new(
                    cell.positions
                        .iter()
                        .zip(&cell.numbers)
                        .map(|(position, species)| {
                            let image = operation.rotation.map(f64::from) * position
                                + operation.translation;
                            cell.positions
                                .iter()
                                .zip(&cell.numbers)
                                .position(|(other, other_species)| {
                                    let delta = image - other;
                                    species == other_species
                                        && (delta - delta.map(f64::round)).norm() < 1e-10
                                })
                                .unwrap()
                        })
                        .collect(),
                )
            })
            .collect();
        // Also exercise a nontrivial input basis and origin. Site order is unchanged.
        let to_input = UnimodularTransformation::new(
            matrix![1, 1, 0; 0, 1, 0; 0, 0, 1],
            vector![0.13, 0.17, 0.19],
        );
        (
            to_input.transform_cell(&cell),
            to_input.transform_operations(&operations),
            permutations,
            SpaceGroup::from_hall_number_and_transformation(hall_number, to_input.inverse())
                .unwrap(),
        )
    }

    fn assert_symmetry(cell: &Cell, operations: &Operations) {
        let metric = cell.lattice.metric_tensor();
        for operation in operations {
            let w = operation.rotation.map(f64::from);
            assert!((w.transpose() * metric * w - metric).norm() < 1e-12 * metric.norm());
            let mut visited = vec![false; cell.num_atoms()];
            for (position, species) in cell.positions.iter().zip(&cell.numbers) {
                let image = w * position + operation.translation;
                let mapped = cell
                    .positions
                    .iter()
                    .zip(&cell.numbers)
                    .position(|(other, other_species)| {
                        let delta = image - other;
                        species == other_species && (delta - delta.map(f64::round)).norm() < 1e-12
                    })
                    .expect("target symmetry must map to a site of the same species");
                assert!(!visited[mapped], "site action must be a permutation");
                visited[mapped] = true;
            }
        }
    }

    fn check_hall_setting(hall_number: HallNumber) {
        let (exact, operations, permutations, space_group) = fixture(hall_number);
        let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        let centering = hall_symbol_entry(hall_number).unwrap().centering;
        let rotation = Rotation3::from_euler_angles(0.37, -0.21, 0.13).into_inner();
        for handedness in [1.0, -1.0] {
            for perturbed in [false, true] {
                let mut input = exact.clone();
                let distortion = if perturbed {
                    matrix![1.0007, 0.0002, -0.0003; 0.0004, 0.9995, 0.0001; 0.0002, -0.0003, 1.0003]
                } else {
                    Matrix3::identity()
                };
                input.lattice.basis = rotation
                    * Matrix3::from_diagonal(&vector![1.0, 1.0, handedness])
                    * distortion
                    * input.lattice.basis;
                if perturbed {
                    for (i, position) in input.positions.iter_mut().enumerate() {
                        *position += 1e-5
                            * vector![
                                (i as f64).sin(),
                                (2.0 * i as f64).cos(),
                                (3.0 * i as f64).sin()
                            ];
                    }
                }
                let standardize = |rotate_basis| {
                    StandardizedCell::new(&input, &operations, &permutations, &space_group, 1e-2, 1e-2, rotate_basis)
                        .unwrap_or_else(|error| panic!("Hall {hall_number}, handedness {handedness}, perturbed {perturbed}: {error}"))
                };
                let unrotated = standardize(false);
                let rotated = standardize(true);
                for standardized in [&unrotated, &rotated] {
                    assert_symmetry(&standardized.prim_cell, &hall_symbol.primitive_traverse());
                    assert_symmetry(&standardized.cell, &hall_symbol.traverse());
                    assert_eq!(standardized.prim_cell.numbers, input.numbers);
                    assert_eq!(
                        standardized.cell.num_atoms(),
                        input.num_atoms() * centering.order()
                    );
                    assert_relative_eq!(
                        standardized.prim_cell.lattice.basis * centering.linear().map(f64::from),
                        standardized.cell.lattice.basis,
                        epsilon = 1e-12
                    );
                    assert_eq!(
                        standardized.cell.lattice.basis.determinant().signum(),
                        handedness
                    );
                    let mut copies = vec![0; input.num_atoms()];
                    for (j, &i) in standardized.site_mapping.iter().enumerate() {
                        copies[i] += 1;
                        assert_eq!(standardized.cell.numbers[j], input.numbers[i]);
                        let delta = centering.linear().map(f64::from)
                            * standardized.cell.positions[j]
                            - standardized.prim_cell.positions[i];
                        assert!((delta - delta.map(f64::round)).norm() < 1e-12);
                    }
                    assert!(copies.into_iter().all(|count| count == centering.order()));
                    if !perturbed {
                        let selected = standardized.prim_transformation.transform_cell(&input);
                        assert_relative_eq!(
                            standardized.rotation_matrix * selected.lattice.basis,
                            standardized.prim_cell.lattice.basis,
                            epsilon = 1e-12
                        );
                        for (position, expected) in standardized
                            .prim_cell
                            .positions
                            .iter()
                            .zip(&selected.positions)
                        {
                            let delta = position - expected;
                            assert!((delta - delta.map(f64::round)).norm() < 1e-12);
                        }
                    }
                }
                assert_relative_eq!(
                    unrotated.rotation_matrix,
                    Matrix3::identity(),
                    epsilon = 1e-12
                );
                assert_relative_eq!(
                    rotated.rotation_matrix.transpose() * rotated.rotation_matrix,
                    Matrix3::identity(),
                    epsilon = 1e-12
                );
                assert_relative_eq!(rotated.rotation_matrix.determinant(), 1.0, epsilon = 1e-12);
                assert_relative_eq!(
                    rotated.rotation_matrix * unrotated.cell.lattice.basis,
                    rotated.cell.lattice.basis,
                    epsilon = 1e-12
                );
                assert_eq!(unrotated.prim_cell.positions, rotated.prim_cell.positions);
                assert_eq!(unrotated.cell.positions, rotated.cell.positions);
                assert_eq!(unrotated.site_mapping, rotated.site_mapping);
                assert_eq!(
                    unrotated
                        .wyckoffs
                        .iter()
                        .map(|w| w.letter)
                        .collect::<Vec<_>>(),
                    rotated
                        .wyckoffs
                        .iter()
                        .map(|w| w.letter)
                        .collect::<Vec<_>>(),
                );
                let selected_basis = unrotated
                    .transformation
                    .transform_lattice(&input.lattice)
                    .basis;
                let stretch = unrotated.cell.lattice.basis * selected_basis.try_inverse().unwrap();
                assert_relative_eq!(stretch, stretch.transpose(), epsilon = 1e-12);
                assert!(stretch.symmetric_eigen().eigenvalues.min() > 0.0);
                for (row, col) in [(1, 0), (2, 0), (2, 1)] {
                    assert!(rotated.cell.lattice.basis[(row, col)].abs() < 1e-12);
                }
                assert!(rotated.cell.lattice.basis[(0, 0)] > 0.0);
                assert!(rotated.cell.lattice.basis[(1, 1)] > 0.0);
            }
        }
    }

    #[test]
    fn test_standardized_cell_cubic_specification() {
        check_hall_setting(523); // Fm-3m, including a general orbit and a special orbit.
    }

    #[test]
    fn test_standardized_cell_all_hall_settings() {
        for hall_number in 1..=530 {
            check_hall_setting(hall_number);
        }
    }
}
