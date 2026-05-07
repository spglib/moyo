use super::primitive_symmetry_search::search_bravais_group;
use crate::base::{AngleTolerance, Lattice, MoyoError, Rotation, Rotations};

/// Search the Bravais group restricted to layer-group form.
///
/// Reuses [`search_bravais_group`] and post-filters to rotations `W` satisfying
/// the layer-group block form (paper eq. 4):
/// `W[0,2] = W[1,2] = W[2,0] = W[2,1] = 0` and `W[2,2] = ±1`.
/// This excludes cubic point groups and any in-plane / out-of-plane mixing.
#[allow(dead_code)] // wired up by later layer-group milestones
pub(crate) fn search_layer_bravais_group(
    minkowski_lattice: &Lattice,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<Rotations, MoyoError> {
    let rotations = search_bravais_group(minkowski_lattice, symprec, angle_tolerance)?;
    Ok(rotations.into_iter().filter(is_layer_block_form).collect())
}

#[allow(dead_code)] // wired up by later layer-group milestones
fn is_layer_block_form(rotation: &Rotation) -> bool {
    rotation[(0, 2)] == 0
        && rotation[(1, 2)] == 0
        && rotation[(2, 0)] == 0
        && rotation[(2, 1)] == 0
        && rotation[(2, 2)].abs() == 1
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;
    use test_log::test;

    use super::{is_layer_block_form, search_layer_bravais_group};
    use crate::base::{AngleTolerance, Lattice};

    #[test]
    fn test_square_lattice() {
        // Square lattice (a=b along x,y; c along z, c != a). Holohedry 4/mmm.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 2.0;
        ]);
        let rotations =
            search_layer_bravais_group(&lattice, 1e-4, AngleTolerance::Radian(1e-2)).unwrap();
        // WHY: 4/mmm has order 16 and every element already fixes the c-axis (W_33 = ±1) with no in/out mixing.
        assert_eq!(rotations.len(), 16);
        for rotation in &rotations {
            assert!(is_layer_block_form(rotation));
        }
    }

    #[test]
    fn test_cubic_lattice() {
        // Cubic lattice (a=b=c, orthogonal). Full Bravais group is m-3m (order 48).
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 1.0;
        ]);
        let rotations =
            search_layer_bravais_group(&lattice, 1e-4, AngleTolerance::Radian(1e-2)).unwrap();
        // WHY: m-3m restricted to elements fixing the z-axis is the 4/mmm site stabilizer (order 16).
        assert_eq!(rotations.len(), 16);
        for rotation in &rotations {
            assert!(is_layer_block_form(rotation));
        }
    }

    #[test]
    fn test_hexagonal_lattice() {
        // Hexagonal lattice (gamma = 120°, c along z). Holohedry 6/mmm.
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            -0.5, f64::sqrt(3.0) / 2.0, 0.0;
            0.0, 0.0, 1.0;
        ]);
        let rotations =
            search_layer_bravais_group(&lattice, 1e-4, AngleTolerance::default()).unwrap();
        // WHY: 6/mmm has order 24 and every element fixes the c-axis (the 6-fold is along c), so the LG filter is a no-op.
        assert_eq!(rotations.len(), 24);
        for rotation in &rotations {
            assert!(is_layer_block_form(rotation));
        }
    }
}
