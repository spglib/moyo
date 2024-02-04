use itertools::izip;
use nalgebra::Vector3;

use crate::base::{orbits_from_permutations, Cell, MoyoError, Transformation};
use crate::search::PrimitiveSymmetrySearch;

pub struct WyckoffPositionAssignment {}

impl WyckoffPositionAssignment {
    // 0. We already know crystallographic orbits for the std conventional cell
    // 1. For each site in std conventional cell, find the site symmetry group
    //      { (Ri, vi + v) | Ri * x + vi + v = x, v in Z^3, (Ri, vi) in coset representatives }
    //      Here, rotation parts should be unique
    /// Assign Wyckoff positions to **primitive** cell
    pub fn new(
        prim_cell: &Cell,
        symmetry_search: &PrimitiveSymmetrySearch,
        to_conv: &Transformation,
        symprec: f64,
    ) -> Result<Self, MoyoError> {
        // Determine a site symmetry group for each site and (crystallographic) orbits
        let orbits = orbits_from_permutations(prim_cell.num_atoms(), &symmetry_search.permutations);

        unimplemented!()
    }
}

#[cfg(test)]
mod tests {}
