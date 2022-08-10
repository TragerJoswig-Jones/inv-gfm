/* Structures for reference frame transformations */
use super::constants::*;
use super::*;
use core::ops::{Add, Sub, Mul, Div, Rem};
// use idsp::cossin;  // TODO: Determine if we want to use fixed point math?

// Enumeration of Acceptable Numbers for Reference Frame Structures
// pub enum Num {
//     F32(f32),
//     FXD(Fxd),
// }  // TODO: Should all <T> and T, be replaced with <Num> and Num,  (https://stackoverflow.com/questions/40776020/is-there-any-way-to-restrict-a-generic-type-to-one-of-several-types) 

// Generic type with a trait bound for acceptable number types for use with reference frame structures.
// Want a similar trait to Num here (https://github.com/rust-num/num-traits/blob/master/src/lib.rs)
pub trait Num<Rhs = Self, Output = Self>: Add<Rhs, Output = Output>
                                        + Sub<Rhs, Output = Output>
                                        + Mul<Rhs, Output = Output>
                                        + Div<Rhs, Output = Output>
                                        + Rem<Rhs, Output = Output> 
                                        + core::ops::Neg
                                        + Copy
{}  // (https://stackoverflow.com/questions/40776020/is-there-any-way-to-restrict-a-generic-type-to-one-of-several-types)
impl Num for f32 {}
impl Num for Fxd {}

// TODO: If we want the crate to work for f32 and Fxd numbers we may need to implement certain functions seperately for each enumeration of Num, but ultimately it should just be two seperate crates
// TODO: Look into using macros to streamline / reduce code repetition if including f32 and fxd in the same crate

/*
ABC Three-Phase Values 
*/
pub struct ABC<T> {
    pub a: T,
    pub b: T,
    pub c: T,
}
pub trait ToFromABC<T> {
    fn to_abc(&self) -> ABC<T>;
    fn from_abc(a: T, b: T, c: T) -> Self;
}
impl ToFromAlphaBeta<Flt> for ABC<Flt> {
    fn to_ab(&self) -> AlphaBeta<Flt> {
        let alpha = TWO_THIRDS * self.a - ONE_THIRD * self.b - ONE_THIRD * self.c;
        let beta = (self.b - self.c) * SQRT_3_OVER_3;
        let gamma = ONE_THIRD * (self.a + self.b + self.c);
        AlphaBeta{ alpha, beta, gamma }
    }
    fn from_ab(alpha: Flt, beta: Flt, gamma: Flt) -> Self {
        let a = alpha + gamma;
        let b = -ONE_HALF * alpha + SQRT_3_OVER_2 * beta + gamma;
        let c = -ONE_HALF * alpha - SQRT_3_OVER_2 * beta + gamma;
        ABC{ a, b, c }
    }
}
impl ToFromDQZ<Flt> for ABC<Flt> {
    fn to_dqz(&self, sin_cos: SinCos<Flt>) -> DQZ<Flt> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: TWO_THIRDS * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: ONE_THIRD * (self.a + self.b + self.c),
        }
    }
    fn from_dqz(d: Flt, q: Flt, z: Flt, sin_cos: SinCos<Flt>) -> Self {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * d + sin_cos.cos_val * q + z,
              b: left.sin_val * d + left.cos_val * q + z, 
              c: right.sin_val * d + right.cos_val * q + z,
            }
    }
}
impl ToFromPolar<Flt> for ABC<Flt> {
    fn to_polar(&self) -> Polar<Flt> {
        // Implemented as ABC to AB then AB to Polar, but throws an error if gamme != 0 (This indicates the ABC signal is not balanced)
        let ab = self.to_ab();
        if ab.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: Flt, theta: Flt) -> Self {  // TODO: Make a macro for this subfunction as it is the same as for Fxd except the SinCos call?
        let sin_cos = SinCos::<Flt>::from_theta(theta);
        ABC{ a: r * sin_cos.sin_val, 
             b: r * sin_cos.rotate_left_120().sin_val, 
             c: r * sin_cos.rotate_right_120().sin_val} // TODO: double check the left/right diretions here
    }
}
impl ToFromAlphaBeta<Fxd> for ABC<Fxd> {
    fn to_ab(&self) -> AlphaBeta<Fxd> {
        let alpha = TWO_THIRDS_FXD * self.a - ONE_THIRD_FXD * self.b - ONE_THIRD_FXD * self.c;
        let beta = (self.b - self.c) * SQRT_3_OVER_3_FXD;
        let gamma = ONE_THIRD_FXD * (self.a + self.b + self.c);
        AlphaBeta{ alpha, beta, gamma }
    }
    fn from_ab(alpha: Fxd, beta: Fxd, gamma: Fxd) -> Self {
        let a = alpha + gamma;
        let b = -ONE_HALF_FXD * alpha + SQRT_3_OVER_2_FXD * beta + gamma;
        let c = -ONE_HALF_FXD * alpha - SQRT_3_OVER_2_FXD* beta + gamma;
        ABC{ a, b, c }
    }
}
impl ToFromDQZ<Fxd> for ABC<Fxd> {
    fn to_dqz(&self, sin_cos: SinCos<Fxd>) -> DQZ<Fxd> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS_FXD * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: TWO_THIRDS_FXD * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: ONE_THIRD_FXD * (self.a + self.b + self.c),
        }
    }
    fn from_dqz(d: Fxd, q: Fxd, z: Fxd, sin_cos: SinCos<Fxd>) -> Self {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * d + sin_cos.cos_val * q + z,
              b: left.sin_val * d + left.cos_val * q + z, 
              c: right.sin_val * d + right.cos_val * q + z,
            }
    }
}
impl ToFromPolar<Fxd> for ABC<Fxd> {
    fn to_polar(&self) -> Polar<Fxd> {
        // Implemented as ABC to AB then AB to Polar, but throws an error if gamme != 0 (This indicates the ABC signal is not balanced)
        let ab = self.to_ab();
        if ab.gamma != ZERO {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: Fxd, theta: Fxd) -> Self {
        let sin_cos = SinCos::<Fxd>::from_theta(theta);
        ABC{ a: r * sin_cos.sin_val, 
             b: r * sin_cos.rotate_left_120().sin_val, 
             c: r * sin_cos.rotate_right_120().sin_val} // TODO: double check the left/right diretions here
    }
}

/* 
Polar Coordinate Values 
*/
pub struct Polar<T> {
    pub r: T,
    pub theta: T,
}
pub trait ToFromPolar<T> {
    fn to_polar(&self) -> Polar<T>;
    fn from_polar(r: T, theta: T) -> Self;
}
impl ToFromABC<Flt> for Polar<Flt> {
    fn to_abc(&self) -> ABC<Flt> {
        ABC::from_polar(self.r, self.theta)
    }
    fn from_abc(a: Flt, b: Flt, c: Flt) -> Self {
        // Creates ABC object and calls its to_polar function
        // TODO: Determine if this should avoid the intermediate ABC object creation
        let abc: ABC<Flt> = ABC { a, b, c };
        abc.to_polar()
    }
}
impl ToFromAlphaBeta<Flt> for Polar<Flt> {
    fn to_ab(&self) -> AlphaBeta<Flt> {
        AlphaBeta::from_polar(self.r, self.theta)
    }
    fn from_ab(alpha: Flt, beta: Flt, gamma: Flt) -> Self {
        if gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", gamma);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(alpha, 2.) + libm::powf(beta, 2.)) / SQRT_2, 
                   theta: libm::atan2f(beta, alpha)}
        }
    }
}
impl ToFromDQZ<Flt> for Polar<Flt> {
    fn to_dqz(&self, sin_cos: SinCos<Flt>) -> DQZ<Flt> {
        let sin_cos_thetas = SinCos::<Flt>::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: 0. 
            }
    }
    fn from_dqz(d: Flt, q: Flt, z: Flt, sin_cos: SinCos<Flt>) -> Self {
        if z != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the DQZ zero value.", z);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(d, 2.) + libm::powf(q, 2.)) / SQRT_2, 
                   theta: libm::atan2f(d * sin_cos.sin_val + q * sin_cos.cos_val, d * sin_cos.cos_val - q * sin_cos.sin_val)
            }
        }
    }
}
impl ToFromABC<Fxd> for Polar<Fxd> {
    fn to_abc(&self) -> ABC<Fxd> {
        ABC::from_polar(self.r, self.theta)
    }
    fn from_abc(a: Fxd, b: Fxd, c: Fxd) -> Self {
        // Creates ABC object and calls its to_polar function
        // TODO: Determine if this should avoid the intermediate ABC object creation
        let abc: ABC<Fxd> = ABC { a, b, c };
        abc.to_polar()
    }
}
impl ToFromAlphaBeta<Fxd> for Polar<Fxd> {
    fn to_ab(&self) -> AlphaBeta<Fxd> {
        AlphaBeta::from_polar(self.r, self.theta)
    }
    fn from_ab(alpha: Fxd, beta: Fxd, gamma: Fxd) -> Self {
        if gamma != ZERO {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", gamma);
        } else {
            Polar{ r: cd::sqrt(alpha * alpha + beta * beta) / SQRT_2_FXD, 
                   theta: cd::atan2(beta, alpha)}
        }
    }
}
impl ToFromDQZ<Fxd> for Polar<Fxd> {
    fn to_dqz(&self, sin_cos: SinCos<Fxd>) -> DQZ<Fxd> {
        let sin_cos_thetas = SinCos::<Fxd>::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: ZERO,
            }
    }
    fn from_dqz(d: Fxd, q: Fxd, z: Fxd, sin_cos: SinCos<Fxd>) -> Self {
        if z != ZERO {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the DQZ zero value.", z);
        } else {
            Polar{ r: cd::sqrt(d * d + q * q) / SQRT_2_FXD, 
                   theta: cd::atan2(d * sin_cos.sin_val + q * sin_cos.cos_val, d * sin_cos.cos_val - q * sin_cos.sin_val)
            }
        }
    }
}

/* 
Alpha-Beta Reference Frame Values 
*/
pub struct AlphaBeta<T> {
    pub alpha: T,
    pub beta: T,
    pub gamma: T,
}
pub trait ToFromAlphaBeta<T> {
    fn to_ab(&self) -> AlphaBeta<T>;
    fn from_ab(alpha: T, beta: T, gamma: T) -> Self;
}
impl AlphaBeta<Flt> {
    pub fn from_ab_(alpha: Flt, beta: Flt) -> AlphaBeta<Flt> {
        AlphaBeta{
            alpha,
            beta,
            gamma: 0.,
        }
    }
}
impl AlphaBeta<Fxd> {
    pub fn from_ab_(alpha: Fxd, beta: Fxd) -> AlphaBeta<Fxd> {
        AlphaBeta{
            alpha,
            beta,
            gamma: ZERO,
        }
    }
}
impl ToFromABC<Flt> for AlphaBeta<Flt> {
    fn to_abc(&self) -> ABC<Flt> {
        ABC{ a: self.alpha + self.gamma, 
             b: -ONE_HALF * self.alpha + SQRT_3_OVER_2 * self.beta + self.gamma, 
             c: -ONE_HALF * self.alpha - SQRT_3_OVER_2 * self.beta + self.gamma 
            }
    }
    fn from_abc(a: Flt, b: Flt, c: Flt) -> Self {
        AlphaBeta{ alpha: TWO_THIRDS * a - ONE_THIRD * b - ONE_THIRD * c,
                   beta: SQRT_3_OVER_3 * (b - c),
                   gamma:  ONE_THIRD * (a + b + c), 
                }
    }
}
impl ToFromABC<Fxd> for AlphaBeta<Fxd> {
    fn to_abc(&self) -> ABC<Fxd> {
        ABC{ a: self.alpha + self.gamma, 
             b: -ONE_HALF_FXD * self.alpha + SQRT_3_OVER_2_FXD * self.beta + self.gamma, 
             c: -ONE_HALF_FXD * self.alpha - SQRT_3_OVER_2_FXD * self.beta + self.gamma 
            }
    }
    fn from_abc(a: Fxd, b: Fxd, c: Fxd) -> Self {
        AlphaBeta{ alpha: TWO_THIRDS_FXD * a - ONE_THIRD_FXD * b - ONE_THIRD_FXD * c,
                   beta: SQRT_3_OVER_3_FXD * (b - c),
                   gamma:  ONE_THIRD_FXD * (a + b + c), 
                }
    }
}

impl ToFromPolar<Flt> for AlphaBeta<Flt> {
    fn to_polar(&self) -> Polar<Flt> {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {  // TODO: Implement for both Fxd and f32
            Polar{ r: libm::sqrtf(libm::powf(self.alpha, 2.) + libm::powf(self.beta, 2.)) / SQRT_2, 
                theta: libm::atan2f(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: Flt, theta: Flt) -> Self {  // TODO: Implement using SinCos
        AlphaBeta {
            alpha: SQRT_2 * r * libm::cosf(theta),
            beta: SQRT_2 * r * libm::sinf(theta),
            gamma: 0.,
        }
    }
}
impl ToFromPolar<Fxd> for AlphaBeta<Fxd> {
    fn to_polar(&self) -> Polar<Fxd> {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {  // TODO: Implement for both Fxd and f32
            Polar{ r: cd::sqrt(self.alpha * self.alpha + self.beta * self.beta) / SQRT_2_FXD, 
                theta: cd::atan2(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: Fxd, theta: Fxd) -> Self {  // TODO: Implement using SinCos
        AlphaBeta {
            alpha: SQRT_2_FXD * r * cd::cos(theta),
            beta: SQRT_2_FXD * r * cd::sin(theta),
            gamma: ZERO,
        }
    }
}
impl<T: Num> ToFromDQZ<T> for AlphaBeta<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let d = sin_cos.cos_val * self.alpha + sin_cos.sin_val * self.beta;
        let q = sin_cos.cos_val * self.beta - (sin_cos.sin_val * self.alpha);
        let z = self.gamma;
        DQZ{ d, q, z }
    }
    fn from_dqz(d: T, q: T, z: T, sin_cos: SinCos<T>) -> Self {
        return AlphaBeta { alpha: sin_cos.cos_val * d - sin_cos.sin_val * q, 
                           beta: sin_cos.sin_val * d + sin_cos.cos_val * q, 
                           gamma: z }
    }
}

/* 
Direct-Quadrature-Zero (DQZ) Reference Frame Values 
*/
pub struct DQZ<T: Num> {
    pub d: T,
    pub q: T,
    pub z: T,
}
pub trait ToFromDQZ<T: Num> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T>;
    fn from_dqz(d: T, q: T, z: T, sin_cos: SinCos<T>) -> Self;
}
// TODO: Should we implement ToFrom functions for the DQZ structure? 
//       The issue is that DQZ transforms require a reference angle.
//       If not constructing a DQZ struct will have to be done from 
//       an existing reference frame struct.

/* 
Stores sine and cosine pair values 
*/
pub struct SinCos<T: Num> {
    sin_val: T,
    cos_val: T,
    theta: T,
}

pub trait Trig<T: Num> {
    fn from_theta(theta: T) -> SinCos<T>;
    fn rotate_right_120(&self) -> SinCos<T>;
    fn rotate_left_120(&self) -> SinCos<T>;
}

impl Trig<Flt> for SinCos<Flt> {
    fn from_theta(theta: Flt) -> Self {
        let sin_val = libm::sinf(theta);
        let cos_val = libm::cosf(theta);
        Self{ sin_val, cos_val, theta }
    }
    // Rotate the reference angle, theta, by 120 degrees counter-clockwise
    fn rotate_right_120(&self) -> SinCos<Flt> {  // TODO: Should we have a function that modifies the values of the SinCos object instead of returning a new one?
        return SinCos{ sin_val: -ONE_HALF * self.sin_val + SQRT_3_OVER_2 * self.cos_val, 
                        cos_val: -ONE_HALF * self.cos_val - SQRT_3_OVER_2 * self.sin_val,
                        theta: self.theta + TWO_PI_OVER_THREE,
            }
    }

    // Rotate the reference angle, theta, by 120 degrees clockwise
    fn rotate_left_120(&self) -> SinCos<Flt> {
        return SinCos{ sin_val: -ONE_HALF * self.sin_val - SQRT_3_OVER_2 * self.cos_val, 
                        cos_val: -ONE_HALF * self.cos_val + SQRT_3_OVER_2 * self.sin_val,
                        theta: self.theta - TWO_PI_OVER_THREE,
                    }
    }
}

impl Trig<Fxd> for SinCos<Fxd> {
    fn from_theta(theta: Fxd) -> SinCos<Fxd> {
        let sin_cos = cd::sin_cos(theta);
        Self{ sin_val: sin_cos.0, cos_val: sin_cos.1, theta }
    }
    // Rotate the reference angle, theta, by 120 degrees counter-clockwise
    fn rotate_right_120(&self) -> SinCos<Fxd> {  // TODO: Should we have a function that modifies the values of the SinCos object instead of returning a new one?
        return SinCos{ sin_val: -ONE_HALF_FXD * self.sin_val + SQRT_3_OVER_2_FXD * self.cos_val, 
                        cos_val: -ONE_HALF_FXD * self.cos_val - SQRT_3_OVER_2_FXD * self.sin_val,
                        theta: self.theta + TWO_PI_OVER_THREE_FXD,
            }
    }

    // Rotate the reference angle, theta, by 120 degrees clockwise
    fn rotate_left_120(&self) -> SinCos<Fxd> {
        return SinCos{ sin_val: -ONE_HALF_FXD * self.sin_val - SQRT_3_OVER_2_FXD * self.cos_val, 
                        cos_val: -ONE_HALF_FXD * self.cos_val + SQRT_3_OVER_2_FXD * self.sin_val,
                        theta: self.theta - TWO_PI_OVER_THREE_FXD,
                    }
    }
}

pub trait SinCosValues<T> {
    fn sin_value(&self) -> T;
    fn cos_value(&self) -> T;
    fn theta_value(&self) -> T;
}

impl<T: Num> SinCosValues<T> for SinCos<T> {
    // Get the sin_val assosiated with the SinCos object
    fn sin_value(&self) -> T {
        self.sin_val
    }

    // Get the cos_val assosiated with the SinCos object
    fn cos_value(&self) -> T {
        self.cos_val
    }

    // Get the theta assosiated with the SinCos object
    fn theta_value(&self) -> T {
        self.theta
    }
}
