#![cfg_attr(not(test), no_std)]

// External crate imports
extern crate nalgebra as na;
pub use nalgebra::SVector as Vec;  //TODO: Replace new nalgebra vector creation with generic new
//TODO: Switch to using fixed crate for fixed-integer math

// unifi-gfm crate modules
pub mod calc;
pub mod constants;
pub mod dvoc;
pub mod droop;
pub mod refs;
pub mod sims;
mod tests;
