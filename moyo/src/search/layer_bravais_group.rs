use crate::base::Rotation;

/// True iff `rotation` has the layer-group block form (paper Fu et al. 2024
/// eq. 4): `W[0,2] = W[1,2] = W[2,0] = W[2,1] = 0` and `W[2,2] = +/-1`.
/// Used to filter bulk space-group rotations down to layer-group ones.
pub(crate) fn is_layer_block_form(rotation: &Rotation) -> bool {
    rotation[(0, 2)] == 0
        && rotation[(1, 2)] == 0
        && rotation[(2, 0)] == 0
        && rotation[(2, 1)] == 0
        && rotation[(2, 2)].abs() == 1
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::is_layer_block_form;

    #[test]
    fn test_layer_block_form_predicate() {
        // Identity, 4-fold along c, mirror through xy plane: layer block form.
        assert!(is_layer_block_form(&matrix![1, 0, 0; 0, 1, 0; 0, 0, 1]));
        assert!(is_layer_block_form(&matrix![0, -1, 0; 1, 0, 0; 0, 0, 1]));
        assert!(is_layer_block_form(&matrix![1, 0, 0; 0, 1, 0; 0, 0, -1]));

        // 3-fold along [111]: mixes in/out, not layer block form.
        assert!(!is_layer_block_form(&matrix![0, 0, 1; 1, 0, 0; 0, 1, 0]));
        // Off-block element on row 0: not layer block form.
        assert!(!is_layer_block_form(&matrix![1, 0, 1; 0, 1, 0; 0, 0, 1]));
        // W[2,2] = 0: not layer block form.
        assert!(!is_layer_block_form(&matrix![1, 0, 0; 0, 1, 0; 0, 0, 0]));
    }
}
