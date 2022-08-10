// A list of commonly used constants
use super::*;

// F32 Constants
pub const SQRT_2: f32 = 1.4142135623730951_f32;
pub const SQRT_3: f32 = 1.7320508075688772_f32;
pub const PI: f32 = 3.1415926535897_f32;
pub const ONE_OVER_SQRT_2: f32 = 1. / SQRT_2;
pub const TWO_PI_OVER_THREE: f32 = 2. * PI / 3.;
pub const SQRT_3_OVER_3: f32 = SQRT_3 / 3.;
pub const ONE_THIRD: f32 = 1. / 3.;
pub const TWO_THIRDS: f32 = 2. / 3.;
pub const ONE_HALF: f32 = 0.5;
pub const SQRT_3_OVER_2: f32 = SQRT_3 / 2.;

// Fixed-point Constants
pub const SQRT_2_FXD: Fxd = Fxd::unwrapped_from_str("1.4142135623730951");
pub const SQRT_3_FXD: Fxd = Fxd::unwrapped_from_str("1.7320508075688772");
pub const PI_FXD: Fxd = Fxd::unwrapped_from_str("3.1415926535897");
pub const ONE_OVER_SQRT_2_FXD: Fxd = Fxd::unwrapped_from_str("0.7071067811865475");
pub const TWO_PI_OVER_THREE_FXD: Fxd = Fxd::unwrapped_from_str("2.0943951023931953");
pub const SQRT_3_OVER_3_FXD: Fxd = Fxd::unwrapped_from_str("0.5773502691896257");
pub const ONE_THIRD_FXD: Fxd = Fxd::unwrapped_from_str("0.3333333333333333");
pub const TWO_THIRDS_FXD: Fxd = Fxd::unwrapped_from_str("0.6666666666666666");
pub const ONE_HALF_FXD: Fxd = Fxd::unwrapped_from_str("0.5");
pub const SQRT_3_OVER_2_FXD: Fxd = Fxd::unwrapped_from_str("0.8660254037844386");
pub const ZERO: Fxd = Fxd::unwrapped_from_str("0");
pub const ONE: Fxd = Fxd::unwrapped_from_str("1");
pub const TWO: Fxd = Fxd::unwrapped_from_str("2");
