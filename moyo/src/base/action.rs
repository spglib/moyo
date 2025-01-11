/// How rotation acts on magnetic moments.
#[derive(Debug, Clone, Copy)]
pub enum RotationMagneticMomentAction {
    /// m -> R @ m
    Polar,
    /// m -> (det R) R @ m
    Axial,
}
