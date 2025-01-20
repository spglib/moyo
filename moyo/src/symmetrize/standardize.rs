use itertools::iproduct;
use log::debug;
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{vector, Matrix3, Vector3};
use std::collections::HashMap;

use crate::base::{
    orbits_from_permutations, project_rotations, Cell, Lattice, MoyoError, Operations, Permutation,
    Position, Rotations, Transformation, UnimodularTransformation, EPS,
};
use crate::data::{
    arithmetic_crystal_class_entry, hall_symbol_entry, iter_wyckoff_positions, HallNumber,
    HallSymbol, LatticeSystem, WyckoffPosition, WyckoffPositionSpace,
};
use crate::identify::SpaceGroup;
use crate::math::SNF;

pub struct StandardizedCell {
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    pub prim_cell: Cell,
    /// Transformation from the input primitive cell to the primitive standardized cell.
    pub prim_transformation: UnimodularTransformation,
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    pub cell: Cell,
    /// Wyckoff positions of sites in the `cell`
    pub wyckoffs: Vec<WyckoffPosition>,
    /// Transformation from the input primitive cell to the standardized cell.
    pub transformation: Transformation,
    /// Rotation matrix to map the lattice of the input primitive cell to that of the standardized cell.
    // ------------------------------------------------------------------------
    // Miscellaneous
    // ------------------------------------------------------------------------
    pub rotation_matrix: Matrix3<f64>,
    /// Mapping from the site in the `cell` to that in the `prim_cell`
    pub site_mapping: Vec<usize>,
}

impl StandardizedCell {
    /// Standardize the input **primitive** cell.
    /// For triclinic space groups, Niggli reduction is performed.
    /// Basis vectors are rotated to be a upper triangular matrix.
    /// TODO: option not to rotate basis vectors
    pub fn new(
        prim_cell: &Cell,
        prim_operations: &Operations,
        prim_permutations: &[Permutation],
        space_group: &SpaceGroup,
        symprec: f64,
    ) -> Result<Self, MoyoError> {
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
            space_group,
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
        space_group: &SpaceGroup,
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
        let entry =
            hall_symbol_entry(space_group.hall_number).ok_or(MoyoError::StandardizationError)?;

        // To standardized primitive cell
        let arithmetic_number = entry.arithmetic_number;
        let lattice_system = arithmetic_crystal_class_entry(arithmetic_number).lattice_system();
        let prim_transformation = match lattice_system {
            LatticeSystem::Triclinic => standardize_triclinic_cell(&prim_cell.lattice),
            _ => space_group.transformation.clone(),
        };
        let prim_std_cell_tmp = prim_transformation.transform_cell(&prim_cell);

        // Symmetrize positions of prim_std_cell by refined symmetry operations
        // 1. Prepare operations in primitive standard
        let hs = HallSymbol::from_hall_number(space_group.hall_number)
            .ok_or(MoyoError::StandardizationError)?;
        let conv_std_operations = hs.traverse();
        let prim_std_operations = Transformation::from_linear(entry.centering.linear())
            .inverse_transform_operations(&conv_std_operations);

        // 2. Reorder permutations because prim_std_operations may have different order from symmetry_search.operations!
        let mut permutation_mapping = HashMap::new();
        let prim_operations = prim_transformation.transform_operations(&prim_operations);
        let prim_rotations = project_rotations(&prim_operations);
        for (prim_rotation, permutation) in prim_rotations.iter().zip(prim_permutations.iter()) {
            permutation_mapping.insert(*prim_rotation, permutation.clone());
        }
        let prim_std_permutations = prim_std_operations
            .iter()
            .map(|ops| permutation_mapping.get(&ops.rotation).unwrap().clone())
            .collect::<Vec<_>>();
        let new_prim_std_positions = symmetrize_positions(
            &prim_std_cell_tmp,
            &prim_std_operations,
            &prim_std_permutations,
        );

        // Note: prim_transformation.transform_cell does not change the order of sites
        let prim_std_cell = Cell::new(
            prim_std_cell_tmp.lattice.clone(),
            new_prim_std_positions,
            prim_std_cell_tmp.numbers.clone(),
        );

        // To (conventional) standardized cell
        let conv_trans = Transformation::from_linear(entry.centering.linear());
        let (std_cell, site_mapping) = conv_trans.transform_cell(&prim_std_cell);

        // Symmetrize lattice
        let (_, rotation_matrix) =
            symmetrize_lattice(&std_cell.lattice, &project_rotations(&conv_std_operations));

        Ok((
            prim_std_cell.rotate(&rotation_matrix),
            prim_std_permutations,
            prim_transformation.clone(),
            std_cell.rotate(&rotation_matrix),
            // prim_transformation * (conv_trans.linear, 0)
            Transformation::new(
                prim_transformation.linear * conv_trans.linear,
                prim_transformation.origin_shift,
            ),
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
        // Group sites in std_cell by crystallographic orbits
        let orbits = orbits_in_cell(
            prim_std_cell.num_atoms(),
            prim_std_permutations,
            site_mapping,
        );
        let mut num_orbits = 0;
        let mut mapping = vec![0; std_cell.num_atoms()]; // [std_cell.num_atoms()] -> [num_orbits]
        let mut remapping = vec![]; // [num_orbits] -> [std_cell.num_atoms()]
        for i in 0..std_cell.num_atoms() {
            if orbits[i] == i {
                mapping[i] = num_orbits;
                remapping.push(i);
                num_orbits += 1;
            } else {
                mapping[i] = mapping[orbits[i]];
            }
        }

        // multiplicities: orbit -> multiplicity
        let mut multiplicities = vec![0; num_orbits];
        for i in 0..std_cell.num_atoms() {
            multiplicities[mapping[i]] += 1;
        }
        // Assign Wyckoff positions to representative sites: orbit -> WyckoffPosition
        let mut representative_wyckoffs = vec![None; num_orbits];
        for (i, position) in std_cell.positions.iter().enumerate() {
            let orbit = mapping[i];
            if representative_wyckoffs[orbit].is_some() {
                continue;
            }

            if let Ok(wyckoff) = assign_wyckoff_position(
                position,
                multiplicities[orbit],
                hall_number,
                &std_cell.lattice,
                symprec,
            ) {
                representative_wyckoffs[orbit] = Some(wyckoff);
            }
        }

        for (i, wyckoff) in representative_wyckoffs.iter().enumerate() {
            if wyckoff.is_none() {
                debug!(
                    "Failed to assign Wyckoff positions with multiplicity {}: {:?}",
                    multiplicities[i], std_cell.positions[i]
                );
            }
        }
        let representative_wyckoffs = representative_wyckoffs
            .into_iter()
            .map(|wyckoff| wyckoff.ok_or(MoyoError::WyckoffPositionAssignmentError))
            .collect::<Result<Vec<_>, _>>()?;

        let wyckoffs = (0..std_cell.num_atoms())
            .map(|i| representative_wyckoffs[mapping[orbits[i]]].clone())
            .collect::<Vec<_>>();
        Ok(wyckoffs)
    }
}

/// * `prim_num_atoms` - Number of atoms in the primitive cell
/// * `prim_permutations` - Permutations of the primitive cell
/// * `site_mapping` - Mapping from the site in cell to that in its primitive cell
pub fn orbits_in_cell(
    prim_num_atoms: usize,
    prim_permutations: &[Permutation],
    site_mapping: &[usize],
) -> Vec<usize> {
    // prim_orbits: [prim_num_atoms] -> [prim_num_atoms]
    let prim_orbits = orbits_from_permutations(prim_num_atoms, prim_permutations);

    let num_atoms = site_mapping.len();
    let mut map = HashMap::new();
    let mut orbits = vec![]; // [num_atoms] -> [num_atoms]
    for i in 0..num_atoms {
        // site_mapping: [num_atoms] -> [prim_num_atoms]
        let key = prim_orbits[site_mapping[i]]; // in [prim_num_atoms]
        map.entry(key).or_insert(i);
        orbits.push(*map.get(&key).unwrap());
    }
    orbits
}

/// Niggli reduction for distorted triclinic lattice systems is numerically so challenging.
/// Thus, we skip checking reduction condition.
fn standardize_triclinic_cell(lattice: &Lattice) -> UnimodularTransformation {
    let (_, linear) = lattice.unchecked_niggli_reduce();
    UnimodularTransformation::from_linear(linear)
}

fn assign_wyckoff_position(
    position: &Position,
    multiplicity: usize,
    hall_number: HallNumber,
    lattice: &Lattice,
    symprec: f64,
) -> Result<WyckoffPosition, MoyoError> {
    for wyckoff in iter_wyckoff_positions(hall_number, multiplicity) {
        // Find variable `y` and integers offset `offset` such that
        //    | lattice * (space.linear * y + space.origin - position - offset) | < symprec.
        // Let SNF decomposition of space.linear be
        //    D = L * space.linear * R.
        // lattice * (space.linear * y + space.origin - position - offset)
        //    = lattice * (L^-1 * D * R^-1 * y + space.origin - position - offset)
        //    = lattice * L^-1 * (D * R^-1 * y + L * (space.origin - position - offset))

        //    = lattice * (D * R^-1 * y - L * (offset + position - space.origin)
        let space = WyckoffPositionSpace::new(wyckoff.coordinates);
        let snf = SNF::new(&space.linear);
        for offset in iproduct!(-1..=1, -1..=1, -1..=1) {
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

/// Symmetrize positions by site symmetry groups
fn symmetrize_positions(
    cell: &Cell,
    operations: &Operations,
    permutations: &[Permutation],
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
                    * cell.positions[inv_perm.apply(i)]
                    + operation.translation
                    - cell.positions[i];
                frac_displacements -= frac_displacements.map(|e| e.round()); // in [-0.5, 0.5]
                acc += frac_displacements;
            }
            cell.positions[i] + acc / (permutations.len() as f64)
        })
        .collect::<Vec<_>>()
}

fn symmetrize_lattice(lattice: &Lattice, rotations: &Rotations) -> (Lattice, Matrix3<f64>) {
    let metric_tensor = lattice.metric_tensor();
    let mut symmetrized_metric_tensor: Matrix3<f64> = rotations
        .iter()
        .map(|rotation| {
            rotation.transpose().map(|e| e as f64) * metric_tensor * rotation.map(|e| e as f64)
        })
        .sum();
    symmetrized_metric_tensor /= rotations.len() as f64;

    // upper-triangular basis
    let tri_basis = Cholesky::new_unchecked(symmetrized_metric_tensor)
        .l()
        .transpose();

    // tri_basis \approx rotation_matrix * lattice.basis
    // QR(tri_basis * lattice.basis^-1) = rotation_matrix * strain
    let qr = QR::new(tri_basis * lattice.basis.try_inverse().unwrap());
    let r = qr.r();
    let signs =
        Matrix3::<f64>::from_diagonal(&vector![sign(r[(0, 0)]), sign(r[(1, 1)]), sign(r[(2, 2)])]);
    let mut rotation_matrix = QR::new(tri_basis * lattice.basis.try_inverse().unwrap()).q();
    // Remove axis-direction freedom
    rotation_matrix *= signs;
    (Lattice::new(tri_basis.transpose()), rotation_matrix)
}

fn sign(x: f64) -> f64 {
    if x > EPS {
        1.0
    } else if x < -EPS {
        -1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::symmetrize_lattice;
    use crate::base::{traverse, Lattice};
    use crate::data::{GeometricCrystalClass, PointGroupRepresentative};

    #[test]
    fn test_symmetrize_lattice_cubic() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0001;
            0.0, -0.999, 0.0;
            0.0, 0.0, -1.0001;
        ]);
        let rep = PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::Oh);
        let rotations = traverse(&rep.generators);

        let (new_lattice, rotation_matrix) = symmetrize_lattice(&lattice, &rotations);
        assert_relative_eq!(new_lattice.basis[(1, 1)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(2, 2)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(0, 1)], 0.0);
        assert_relative_eq!(new_lattice.basis[(0, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 1)], 0.0);

        assert_relative_eq!(
            rotation_matrix * lattice.basis,
            new_lattice.basis,
            epsilon = 1e-2
        );
    }
}
