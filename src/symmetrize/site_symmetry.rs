// 0. We already know crystallographic orbits for the std conventional cell
// 1. For each site in std conventional cell, find the site symmetry group
//      { (Ri, vi + v) | Ri * x + vi + v = x, v in Z^3, (Ri, vi) in coset representatives }
//      Here, rotation parts should be unique
// How to store Wyckoff positions?
// "multiplicity" is defined on conventional cell!

#[cfg(test)]
mod tests {}
