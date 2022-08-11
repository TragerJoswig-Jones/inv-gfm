/* Structures for reference frame transformations */
use super::constants::*;
use super::*;

// TODO: If we want the crate to work for f32 and Fxd numbers we may need to implement certain functions seperately for each enumeration of Num, but ultimately it should just be two seperate crates
// TODO: Look into using macros to streamline / reduce code repetition if including f32 and fxd in the same crate

/*
ABC Three-Phase Values 
*/
pub struct ABC<T: Num> {
    pub a: T,
    pub b: T,
    pub c: T,
}
pub trait ToFromABC<T: Num> {
    fn to_abc(&self) -> ABC<T>;
    fn from_abc(a: T, b: T, c: T) -> Self;
}
impl<T: Num> ToFromAlphaBeta<T> for ABC<T> {
    fn to_ab(&self) -> AlphaBeta<T> {
        let alpha = T::from_fixed(TWO_THIRDS) * self.a - T::from_fixed(ONE_THIRD) * self.b - T::from_fixed(ONE_THIRD) * self.c;
        let beta = (self.b - self.c) * T::from_fixed(SQRT_3_OVER_3);
        let gamma = T::from_fixed(ONE_THIRD) * (self.a + self.b + self.c);
        AlphaBeta{ alpha, beta, gamma }
    }
    fn from_ab(alpha: T, beta: T, gamma: T) -> Self {
        let a = alpha + gamma;
        let b = T::from_fixed(-ONE_HALF) * alpha + T::from_fixed(SQRT_3_OVER_2) * beta + gamma;
        let c = T::from_fixed(-ONE_HALF) * alpha - T::from_fixed(SQRT_3_OVER_2)* beta + gamma;
        ABC{ a, b, c }
    }
}
impl<T: Num> ToFromDQZ<T> for ABC<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: T::from_fixed(TWO_THIRDS) * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: T::from_fixed(TWO_THIRDS) * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: T::from_fixed(ONE_THIRD) * (self.a + self.b + self.c),
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
        if ab.gamma != T::from_fixed(ZERO) {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: T, theta: T) -> Self {
        let sin_cos = SinCos::<T>::from_theta(theta);
        ABC{ a: r * sin_cos.sin_val, 
             b: r * sin_cos.rotate_left_120().sin_val, 
             c: r * sin_cos.rotate_right_120().sin_val} // TODO: double check the left/right diretions here
    }
}

/* 
Polar Coordinate Values 
*/
pub struct Polar<T: Num> {
    pub r: T,
    pub theta: T,
}
pub trait ToFromPolar<T: Num> {
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
        if gamma != T::from_fixed(ZERO) {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", gamma);
        } else {
            Polar{ r: cd::sqrt(alpha * alpha + beta * beta) / T::from_fixed(SQRT_2), 
                   theta: cd::atan2(beta, alpha)}
        }
    }
}
impl<T: Num> ToFromDQZ<T> for Polar<T> {
    fn to_dqz(&self, sin_cos: SinCos<T>) -> DQZ<T> {
        let sin_cos_thetas = SinCos::<T>::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: T::from_fixed(ZERO),
            }
    }
    fn from_dqz(d: T, q: T, z: T, sin_cos: SinCos<T>) -> Self {
        if z != T::from_fixed(ZERO) {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the DQZ zero value.", z);
        } else {
            Polar{ r: cd::sqrt(d * d + q * q) / T::from_fixed(SQRT_2), 
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
impl<T: Num> AlphaBeta<T> {
    pub fn from_ab_(alpha: T, beta: T) -> AlphaBeta<T> {
        AlphaBeta{
            alpha,
            beta,
            gamma: T::from_fixed(ZERO),
        }
    }
}
impl<T: Num> ToFromABC<T> for AlphaBeta<T> {
    fn to_abc(&self) -> ABC<T> {
        ABC{ a: self.alpha + self.gamma, 
             b: -T::from_fixed(ONE_HALF) * self.alpha + T::from_fixed(SQRT_3_OVER_2) * self.beta + self.gamma, 
             c: -T::from_fixed(ONE_HALF) * self.alpha - T::from_fixed(SQRT_3_OVER_2) * self.beta + self.gamma 
            }
    }
    fn from_abc(a: T, b: T, c: T) -> Self {
        AlphaBeta{ alpha: T::from_fixed(TWO_THIRDS) * a - T::from_fixed(ONE_THIRD) * b - T::from_fixed(ONE_THIRD) * c,
                   beta: T::from_fixed(SQRT_3_OVER_3) * (b - c),
                   gamma:  T::from_fixed(ONE_THIRD) * (a + b + c), 
                }
    }
}
impl<T: Num> ToFromPolar<T> for AlphaBeta<T> {
    fn to_polar(&self) -> Polar<T> {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {  // TODO: Implement for both T and f32
            Polar{ r: cd::sqrt(self.alpha * self.alpha + self.beta * self.beta) / T::from_fixed(SQRT_2), 
                theta: cd::atan2(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: T, theta: T) -> Self {  // TODO: Implement using SinCos
        AlphaBeta {
            alpha: T::from_fixed(SQRT_2) * r * cd::cos(theta),
            beta: T::from_fixed(SQRT_2) * r * cd::sin(theta),
            gamma: T::from_fixed(ZERO),
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

impl<T: Num> Trig<T> for SinCos<T> {
    fn from_theta(theta: T) -> SinCos<T> {
        let sin_cos = cd::sin_cos(theta);
        Self{ sin_val: sin_cos.0, cos_val: sin_cos.1, theta }
    }
    // Rotate the reference angle, theta, by 120 degrees counter-clockwise
    fn rotate_right_120(&self) -> SinCos<T> {  // TODO: Should we have a function that modifies the values of the SinCos object instead of returning a new one?
        return SinCos{ sin_val: -T::from_fixed(ONE_HALF) * self.sin_val + T::from_fixed(SQRT_3_OVER_2) * self.cos_val, 
                        cos_val: -T::from_fixed(ONE_HALF) * self.cos_val - T::from_fixed(SQRT_3_OVER_2) * self.sin_val,
                        theta: self.theta +  T::from_fixed(TWO_PI_OVER_THREE),
            }
    }

    // Rotate the reference angle, theta, by 120 degrees clockwise
    fn rotate_left_120(&self) -> SinCos<T> {
        return SinCos{ sin_val: -T::from_fixed(ONE_HALF) * self.sin_val - T::from_fixed(SQRT_3_OVER_2) * self.cos_val, 
                        cos_val: -T::from_fixed(ONE_HALF) * self.cos_val + T::from_fixed(SQRT_3_OVER_2) * self.sin_val,
                        theta: self.theta -  T::from_fixed(TWO_PI_OVER_THREE),
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
