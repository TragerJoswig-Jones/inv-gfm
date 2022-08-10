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
                                        //+ Mul<i32, Output = Output>
{}  // pub trait Num: core::ops::Add<Output = Self> + core::ops::Mul<Output = Self>  {}; (https://stackoverflow.com/questions/40776020/is-there-any-way-to-restrict-a-generic-type-to-one-of-several-types)
// TODO: Add multiply by integer scalar to above list???
// impl Mul<i32> for f32 {
//     type Output = f32;
    
//     fn mul(self, rhs: i32) -> f32 {
//         return self * rhs as f32
//     }
// } 
impl Num for f32 {}
impl Num for Fxd {}

// TODO: If we want the crate to work for f32 and Fxd numbers we may need to implement certain functions seperately for each enumeration of Num

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
impl<T: Num> ToFromAlphaBeta<T> for ABC<T> {
    fn to_ab(&self) -> AlphaBeta<T> {
        let alpha = TWO_THIRDS * self.a - ONE_THIRD * self.b - ONE_THIRD * self.c;
        let beta = (self.b - self.c) * SQRT_3_OVER_3;
        let gamma = ONE_THIRD * (self.a + self.b + self.c);
        AlphaBeta{ alpha, beta, gamma }
    }
    fn from_ab(alpha: T, beta: T, gamma: T) -> Self {
        let a = alpha + gamma;
        let b = -ONE_HALF * alpha + SQRT_3_OVER_2 * beta + gamma;
        let c = -ONE_HALF * alpha - SQRT_3_OVER_2 * beta + gamma;
        ABC{ a, b, c }
    }
}
impl<T: Num> ToFromDQZ<T> for ABC<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: TWO_THIRDS * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: ONE_THIRD * (self.a + self.b + self.c),
        }
    }
    fn from_dqz(d: T, q: T, z: T, sin_cos: SinCos<T>) -> Self {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * d + sin_cos.cos_val * q + z,
              b: left.sin_val * d + left.cos_val * q + z, 
              c: right.sin_val * d + right.cos_val * q + z,
            }
    }
}
impl<T: Num> ToFromPolar<T> for ABC<T> {
    fn to_polar(&self) -> Polar<T> {
        // Implemented as ABC to AB then AB to Polar, but throws an error if gamme != 0 (This indicates the ABC signal is not balanced)
        let ab = self.to_ab();
        if ab.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: T, theta: T) -> Self {
        let sin_cos = SinCos::from_theta(theta);
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
impl<T: Num> ToFromABC<T> for Polar<T> {
    fn to_abc(&self) -> ABC<T> {
        ABC::from_polar(self.r, self.theta)
    }
    fn from_abc(a: T, b: T, c: T) -> Self {
        // Creates ABC object and calls its to_polar function
        // TODO: Determine if this should avoid the intermediate ABC object creation
        let abc: ABC<T> = ABC { a, b, c };
        abc.to_polar()
    }
}
impl<T: Num> ToFromAlphaBeta<T> for Polar<T> {
    fn to_ab(&self) -> AlphaBeta<T> {
        AlphaBeta::from_polar(self.r, self.theta)
    }
    fn from_ab(alpha: T, beta: T, gamma: T) -> Self {
        if gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", gamma);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(alpha, 2.) + libm::powf(beta, 2.)) / SQRT_2, 
                theta: libm::atan2f(beta, alpha)}
        }
    }
}
impl<T: Num> ToFromDQZ<T> for Polar<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let sin_cos_thetas = SinCos::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: 0. 
            }
    }
    fn from_dqz(d: T, q: T, z: T, sin_cos: SinCos<T>) -> Self {
        if z != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the DQZ zero value.", z);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(d, 2.) + libm::powf(q, 2.)) / SQRT_2, 
                   theta: libm::atan2f(d * sin_cos.sin_val + q * sin_cos.cos_val, d * sin_cos.cos_val - q * sin_cos.sin_val)
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
impl<T: Num> AlphaBeta<T> {
    pub fn from_ab_(alpha: T, beta: T) -> AlphaBeta<T> {
        AlphaBeta{
            alpha,
            beta,
            gamma: 0.,
        }
    }
}
impl<T: Num> ToFromABC<T> for AlphaBeta<T> {
    fn to_abc(&self) -> ABC<T> {
        ABC{ a: self.alpha + self.gamma, 
             b: -ONE_HALF * self.alpha + SQRT_3_OVER_2 * self.beta + self.gamma, 
             c: -ONE_HALF * self.alpha - SQRT_3_OVER_2 * self.beta + self.gamma 
            }
    }
    fn from_abc(a: T, b: T, c: T) -> Self {
        AlphaBeta{ alpha: TWO_THIRDS * a - ONE_THIRD * b - ONE_THIRD * c,
                   beta: SQRT_3_OVER_3 * (b - c),
                   gamma:  ONE_THIRD * (a + b + c), 
                }
    }
}
impl<T: Num> ToFromPolar<T> for AlphaBeta<T> {
    fn to_polar(&self) -> Polar<T> {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {  // TODO: Implement for both Fxd and f32
            Polar{ r: libm::sqrtf(libm::powf(self.alpha, 2.) + libm::powf(self.beta, 2.)) / SQRT_2, 
                theta: libm::atan2f(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: T, theta: T) -> Self {  // TODO: Implement using SinCos
        AlphaBeta {
            alpha: SQRT_2 * r * libm::cosf(theta),
            beta: SQRT_2 * r * libm::sinf(theta),
            gamma: 0.,
        }
    }
}
impl<T: Num> ToFromDQZ<T> for AlphaBeta<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let d = sin_cos.cos_val * self.alpha + sin_cos.sin_val * self.beta;
        let q = -sin_cos.sin_val * self.alpha + sin_cos.cos_val * self.beta;
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
pub struct DQZ<T> {
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

impl Trig<f32> for SinCos<f32> {
    fn from_theta(theta: f32) -> Self {
        let sin_val = libm::sinf(theta);
        let cos_val = libm::cosf(theta);
        Self{ sin_val, cos_val, theta }
    }
    // Rotate the reference angle, theta, by 120 degrees counter-clockwise
    fn rotate_right_120(&self) -> SinCos<f32> {  // TODO: Should we have a function that modifies the values of the SinCos object instead of returning a new one?
        return SinCos{ sin_val: -ONE_HALF * self.sin_val + SQRT_3_OVER_2 * self.cos_val, 
                        cos_val: -ONE_HALF * self.cos_val - SQRT_3_OVER_2 * self.sin_val,
                        theta: self.theta + TWO_PI_OVER_THREE,
            }
    }

    // Rotate the reference angle, theta, by 120 degrees clockwise
    fn rotate_left_120(&self) -> SinCos<f32> {
        return SinCos{ sin_val: -ONE_HALF * self.sin_val - SQRT_3_OVER_2 * self.cos_val, 
                        cos_val: -ONE_HALF * self.cos_val + SQRT_3_OVER_2 * self.sin_val,
                        theta: self.theta - TWO_PI_OVER_THREE,
                    }
    }
}

impl Trig<Fxd> for SinCos<Fxd> {
    fn from_theta(theta: Fxd) -> Self {
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
