#![cfg_attr(not(test), no_std)]

// External crate imports
extern crate cordic as cd;
extern crate fixed as fx;
extern crate nalgebra as na;
pub use nalgebra::SVector as Vec;  //TODO: Replace new nalgebra vector creation with generic new
pub use fixed::types::I2F30;
pub use fixed::FixedI32;
pub use fixed::types::extra::U29;

type Fxd = fixed::types::I3F29;
//type Fxd = fixed::FixedI32::<U29>;
type Flt = f32;
//TODO: Switch to using fixed crate for fixed-integer math

// unifi-gfm crate modules
pub mod calc;
pub mod constants;
pub mod dvoc;
pub mod droop;
pub mod refs;
pub mod sims;
mod tests;
