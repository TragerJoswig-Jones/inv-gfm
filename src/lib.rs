#![cfg_attr(not(test), no_std)]

/* External crate imports */
extern crate nalgebra as na;
pub use nalgebra::SVector as Vec;  //TODO: Replace new nalgebra vector creation with generic new
use core::ops::{Add, Sub, Mul, Div, Rem};


/* unifi-gfm crate modules */
pub mod calc;
pub mod constants;
pub mod gfm;
pub mod refs;
pub mod sims;
#[cfg(test)]
mod tests;

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
impl Num for f64 {}
