#![cfg_attr(not(test), no_std)]

/* External crate imports */
extern crate nalgebra as na;
pub use nalgebra::SVector as Vec;
use core::ops::{Add, Sub, Mul, Div, Rem};

/* unifi-gfm crate modules */
pub mod calculations;
pub mod constants;
pub mod dc;
pub mod dynamics;
pub mod gfm;
pub mod gfl;
pub mod inverter;
pub mod osg;
pub mod pll;
pub mod reference_frames;
pub mod simulations;

#[cfg(test)]
mod tests;  // Meaningful tests are not implemented yet...

/* Generic type with a trait bound for acceptable number types for use with reference frame structures */
pub trait Num<Rhs = Self, Output = Self>: Add<Rhs, Output = Output>
                                        + Sub<Rhs, Output = Output>
                                        + Mul<Rhs, Output = Output>
                                        + Div<Rhs, Output = Output>
                                        + Rem<Rhs, Output = Output> 
                                        + core::ops::Neg
                                        + Copy
{}  // (https://stackoverflow.com/questions/40776020/is-there-any-way-to-restrict-a-generic-type-to-one-of-several-types)
impl Num for f32 {}
impl Num for f64 {}  // TODO: Determine if implementations can be generic such that f32 or f64 can be used? The issue is caused by adding or multiplying by floats within dynamics steps (https://stackoverflow.com/questions/36013519/how-can-i-make-a-rust-function-accept-any-floating-type-as-an-argument, https://docs.rs/num/latest/num/trait.Float.html)
