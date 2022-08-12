#![cfg_attr(not(test), no_std)]

/* External crate imports */
extern crate cordic as cd;
extern crate fixed as fx;
extern crate nalgebra as na;
pub use fixed::types::I2F30;
pub use fixed::{FixedI32, FixedI64};
pub use fixed::traits::Fixed;
pub use nalgebra::SVector as Vec;  //TODO: Replace new nalgebra vector creation with generic new
use cd::CordicNumber;
use fixed::types::extra::{IsLessOrEqual, True, Unsigned, U29, U30, U32, U61, U62, U64};

/* unifi-gfm crate modules */
pub mod calc;
pub mod constants;
pub mod dvoc;
pub mod droop;
pub mod refs;
pub mod sims;
#[cfg(test)]
mod tests;

/* Generic type with a trait bound for acceptable number types for use with reference frame structures */
pub trait Num<Rhs = Self, Output = Self>: Copy
                                        + fx::traits::Fixed
                                        + CordicNumber
{}  // Creates a trait to restrict the types for structs and functions (https://stackoverflow.com/questions/40776020/is-there-any-way-to-restrict-a-generic-type-to-one-of-several-types)
impl<Fract> Num for FixedI32<Fract>   // Restricts the useable fixed-point numbers to 32-bit with fractional less than 29 which aligns with CordicNumber type (https://github.com/sebcrozet/cordic/blob/0cb0773e879721ad8c72cd36dcb7eb27bd2f83a4/cordic/src/cordic_number.rs#L97-L103)
where 
    Fract: 'static 
          + Unsigned 
          + IsLessOrEqual<U32, Output = True>
          + IsLessOrEqual<U30, Output = True>
          + IsLessOrEqual<U29, Output = True>
{}
impl<Fract> Num for FixedI64<Fract>   // Restricts the useable fixed-point numbers to 64-bit with fractional less than 61 which aligns with CordicNumber type (https://github.com/sebcrozet/cordic/blob/0cb0773e879721ad8c72cd36dcb7eb27bd2f83a4/cordic/src/cordic_number.rs#L205-L211)
where 
    Fract: 'static 
          + Unsigned 
          + IsLessOrEqual<U64, Output = True>
          + IsLessOrEqual<U62, Output = True>
          + IsLessOrEqual<U61, Output = True>
{}
// TODO: Add support for FixedI64 or other fixed-point bit sizes?
