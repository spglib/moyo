#[derive(Debug, Copy, Clone)]
pub enum AngleTolerance {
    Radian(f64),
    Default,
}

pub const EPS: f64 = 1e-8;
