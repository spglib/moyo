use itertools::{iproduct, izip};
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{vector, Matrix3, Vector3};
use std::collections::HashMap;

use crate::base::{
    orbits_from_permutations, Cell, Lattice, MoyoError, Operations, Permutation, Position,
    Rotation, Transformation, UnimodularTransformation, EPS,
};
use crate::data::{
    arithmetic_crystal_class_entry, hall_symbol_entry, iter_wyckoff_positions, HallNumber,
    HallSymbol, LatticeSystem, WyckoffPosition, WyckoffPositionSpace,
};
use crate::identify::SpaceGroup;
use crate::math::SNF;
use crate::search::{PrimitiveCell, PrimitiveSymmetrySearch};

pub struct StandardizedCell {
    pub prim_cell: Cell,
    /// Transformation from the input primitive cell to the primitive standardized cell.
    pub prim_transformation: UnimodularTransformation,
    pub cell: Cell,
    /// Wyckoff positions of sites in the `cell`
    pub wyckoffs: Vec<WyckoffPosition>,
    /// Transformation from the input primitive cell to the standardized cell.
    pub transformation: Transformation,
    /// Rotation matrix to map the lattice of the input primitive cell to that of the standardized cell.
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
        prim_cell: &PrimitiveCell,
        symmetry_search: &PrimitiveSymmetrySearch,
        space_group: &SpaceGroup,
        symprec: f64,
    ) -> Result<Self, MoyoError> {
        let entry = hall_symbol_entry(space_group.hall_number);
        let arithmetic_number = entry.arithmetic_number;
        let bravais_class = arithmetic_crystal_class_entry(arithmetic_number).bravais_class;
        let lattice_system = LatticeSystem::from_bravais_class(bravais_class);

        // To standardized primitive cell
        let prim_transformation = match lattice_system {
            LatticeSystem::Triclinic => {
                let (_, linear) = prim_cell.cell.lattice.niggli_reduce()?;
                UnimodularTransformation::from_linear(linear)
            }
            _ => space_group.transformation.clone(),
        };
        let new_prim_positions = symmetrize_positions(
            &prim_cell.cell,
            &symmetry_search.operations,
            &symmetry_search.permutations,
        );
        // Note: prim_transformation.transform_cell does not change the order of sites
        let prim_std_cell = prim_transformation.transform_cell(&Cell::new(
            prim_cell.cell.lattice.clone(),
            new_prim_positions,
            prim_cell.cell.numbers.clone(),
        ));

        // To (conventional) standardized cell
        let conv_trans = Transformation::from_linear(entry.centering.linear());
        let (std_cell, site_mapping) = conv_trans.transform_cell(&prim_std_cell);

        // Group sites in std_cell by crystallographic orbits
        let orbits = orbits_in_cell(
            prim_cell.cell.num_atoms(),
            &symmetry_search.permutations,
            &site_mapping,
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
        let representative_wyckoffs: Result<Vec<_>, MoyoError> = multiplicities
            .iter()
            .enumerate()
            .map(|(orbit, multiplicity)| {
                let i = remapping[orbit];
                assign_wyckoff_position(
                    &std_cell.positions[i],
                    *multiplicity,
                    space_group.hall_number,
                    &std_cell.lattice,
                    symprec,
                )
            })
            .collect();
        let representative_wyckoffs = representative_wyckoffs?;

        let wyckoffs = (0..std_cell.num_atoms())
            .map(|i| representative_wyckoffs[mapping[orbits[i]]].clone())
            .collect::<Vec<_>>();

        // Symmetrize lattice
        let std_rotations = HallSymbol::from_hall_number(entry.hall_number)
            .traverse()
            .rotations;
        let (_, rotation_matrix) = symmetrize_lattice(&std_cell.lattice, &std_rotations);

        Ok(StandardizedCell {
            prim_cell: prim_std_cell.rotate(&rotation_matrix),
            prim_transformation: prim_transformation.clone(),
            cell: std_cell.rotate(&rotation_matrix),
            wyckoffs,
            // prim_transformation * (conv_trans.linear, 0)
            transformation: Transformation::new(
                prim_transformation.linear * conv_trans.linear,
                prim_transformation.origin_shift,
            ),
            rotation_matrix,
            site_mapping,
        })
    }
}

/// * `prim_num_atoms` - Number of atoms in the primitive cell
/// * `prim_permutations` - Permutations of the primitive cell
/// * `site_mapping` - Mapping from the site in cell to that in its primitive cell
pub fn orbits_in_cell(
    prim_num_atoms: usize,
    prim_permutations: &Vec<Permutation>,
    site_mapping: &Vec<usize>,
) -> Vec<usize> {
    // prim_orbits: [prim_num_atoms] -> [prim_num_atoms]
    let prim_orbits = orbits_from_permutations(prim_num_atoms, &prim_permutations);

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

fn assign_wyckoff_position(
    position: &Position,
    multiplicity: usize,
    hall_number: HallNumber,
    lattice: &Lattice,
    symprec: f64,
) -> Result<WyckoffPosition, MoyoError> {
    for wyckoff in iter_wyckoff_positions(hall_number, multiplicity) {
        // Find variable `y` and integers offset `offset` such that
        //    | space.linear * y + space.origin - position - offset | < symprec.
        // Let SNF decomposition of space.linear be
        //    D = L * space.linear * R.
        // argmin_{y in Q^3, n in Z^3} | space.linear * y + space.origin - position - offset |
        //    = argmin_{y in Q^3, n in Z^3} | L^-1 * D * R^-1 * y + space.origin - position - offset |
        //    = argmin_{y in Q^3, n in Z^3} | D * R^-1 * y - L * (offset + position - space.origin) |
        let space = WyckoffPositionSpace::new(&wyckoff.coordinates);
        let snf = SNF::new(&space.linear);
        for offset in iproduct!(-1..=1, -1..=1, -1..=1) {
            let offset = Vector3::new(offset.0 as f64, offset.1 as f64, offset.2 as f64);
            let b = snf.l.map(|e| e as f64) * (offset + position - space.origin);
            let mut rinvy = Vector3::zeros();
            let mut valid = true;
            for i in 0..3 {
                if snf.d[(i, i)] == 0 {
                    let mut bi = Vector3::zeros();
                    bi[i] = b[i];
                    if lattice.cartesian_coords(&bi).norm() > symprec {
                        valid = false;
                        break;
                    }
                } else {
                    rinvy[i] = b[i] / snf.d[(i, i)] as f64;
                }
            }

            if !valid {
                continue;
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
    permutations: &Vec<Permutation>,
) -> Vec<Position> {
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();
    (0..cell.num_atoms())
        .map(|i| {
            let mut acc = Vector3::zeros();
            for (inv_perm, rotation, translation) in izip!(
                inverse_permutations.iter(),
                operations.rotations.iter(),
                operations.translations.iter(),
            ) {
                if inv_perm.apply(i) == i {
                    continue;
                }
                let mut frac_displacements =
                    rotation.map(|e| e as f64) * cell.positions[inv_perm.apply(i)] + translation
                        - cell.positions[i];
                frac_displacements -= frac_displacements.map(|e| e.round()); // in [-0.5, 0.5]
                acc += frac_displacements;
            }
            cell.positions[i] + acc / (permutations.len() as f64)
        })
        .collect::<Vec<_>>()
}

fn symmetrize_lattice(lattice: &Lattice, rotations: &Vec<Rotation>) -> (Lattice, Matrix3<f64>) {
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
    (Lattice::new(tri_basis), rotation_matrix)
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
