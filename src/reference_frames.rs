/* Structures for reference frame transformations */
use super::constants::*;
use super::*;
// use idsp::cossin;  // TODO: Determine if we want to use floating-point math?

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
impl ToFromAlphaBeta<f32> for ABC<f32> {
    fn to_ab(&self) -> AlphaBeta<f32> {
        let alpha = TWO_THIRDS * self.a - ONE_THIRD * self.b - ONE_THIRD * self.c;
        let beta = (self.b - self.c) * SQRT_3_OVER_3;
        let gamma = ONE_THIRD * (self.a + self.b + self.c);
        AlphaBeta{ alpha, beta, gamma }
    }
    fn from_ab(alpha: f32, beta: f32, gamma: f32) -> Self {
        let a = alpha + gamma;
        let b = -ONE_HALF * alpha + SQRT_3_OVER_2 * beta + gamma;
        let c = -ONE_HALF * alpha - SQRT_3_OVER_2 * beta + gamma;
        ABC{ a, b, c }
    }
}
impl ToFromDQZ<f32> for ABC<f32> {
    fn to_dqz(&self, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: TWO_THIRDS * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: ONE_THIRD * (self.a + self.b + self.c),
        }
    }
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: &SinCos<f32>) -> Self {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * d + sin_cos.cos_val * q + z,
              b: left.sin_val * d + left.cos_val * q + z, 
              c: right.sin_val * d + right.cos_val * q + z,
            }
    }
}
impl ToFromPolar<f32> for ABC<f32> {
    fn to_polar(&self) -> Polar<f32> {
        // Implemented as ABC to AB then AB to Polar, but throws an error if gamme != 0 (This indicates the ABC signal is not balanced)
        let ab = self.to_ab();
        if ab.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: f32, theta: f32) -> Self {  // TODO: Make a macro for this subfunction as it is the same as for Fxd except the SinCos call?
        let sin_cos = SinCos::<f32>::from_theta(theta);
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
impl ToFromABC<f32> for Polar<f32> {
    fn to_abc(&self) -> ABC<f32> {
        ABC::from_polar(self.r, self.theta)
    }
    fn from_abc(a: f32, b: f32, c: f32) -> Self {
        // Creates ABC object and calls its to_polar function
        // TODO: Determine if this should avoid the intermediate ABC object creation
        let abc: ABC<f32> = ABC { a, b, c };
        abc.to_polar()
    }
}
impl ToFromAlphaBeta<f32> for Polar<f32> {
    fn to_ab(&self) -> AlphaBeta<f32> {
        AlphaBeta::from_polar(self.r, self.theta)
    }
    fn from_ab(alpha: f32, beta: f32, gamma: f32) -> Self {
        if gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", gamma);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(alpha, 2.) + libm::powf(beta, 2.)) / SQRT_2, 
                   theta: libm::atan2f(beta, alpha)}
        }
    }
}
impl ToFromDQZ<f32> for Polar<f32> {
    fn to_dqz(&self, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        let sin_cos_thetas = SinCos::<f32>::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: 0. 
            }
    }
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: &SinCos<f32>) -> Self {
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
impl AlphaBeta<f32> {
    pub fn from_ab_(alpha: f32, beta: f32) -> AlphaBeta<f32> {
        AlphaBeta{
            alpha,
            beta,
            gamma: 0.,
        }
    }
}
impl ToFromABC<f32> for AlphaBeta<f32> {
    fn to_abc(&self) -> ABC<f32> {
        ABC{ a: self.alpha + self.gamma, 
             b: -ONE_HALF * self.alpha + SQRT_3_OVER_2 * self.beta + self.gamma, 
             c: -ONE_HALF * self.alpha - SQRT_3_OVER_2 * self.beta + self.gamma 
            }
    }
    fn from_abc(a: f32, b: f32, c: f32) -> Self {
        AlphaBeta{ alpha: TWO_THIRDS * a - ONE_THIRD * b - ONE_THIRD * c,
                   beta: SQRT_3_OVER_3 * (b - c),
                   gamma:  ONE_THIRD * (a + b + c), 
                }
    }
}
impl ToFromPolar<f32> for AlphaBeta<f32> {
    fn to_polar(&self) -> Polar<f32> {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {  // TODO: Implement for both Fxd and f32
            Polar{ r: libm::sqrtf(libm::powf(self.alpha, 2.) + libm::powf(self.beta, 2.)) / SQRT_2, 
                theta: libm::atan2f(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: f32, theta: f32) -> Self {  // TODO: Implement using SinCos
        AlphaBeta {
            alpha: SQRT_2 * r * libm::cosf(theta),
            beta: SQRT_2 * r * libm::sinf(theta),
            gamma: 0.,
        }
    }
}
impl<T: Num> ToFromDQZ<T> for AlphaBeta<T> {
    fn to_dqz(&self, sin_cos: &SinCos<T>) -> DQZ<T> {
        let d = sin_cos.cos_val * self.alpha + sin_cos.sin_val * self.beta;
        let q = sin_cos.cos_val * self.beta - (sin_cos.sin_val * self.alpha);
        let z = self.gamma;
        DQZ{ d, q, z }
    }
    fn from_dqz(d: T, q: T, z: T, sin_cos: &SinCos<T>) -> Self {
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
    fn to_dqz(&self, sin_cos: &SinCos<T>) -> DQZ<T>;
    fn from_dqz(d: T, q: T, z: T, sin_cos: &SinCos<T>) -> Self;
}
// Implement ToFrom functions with all other ref frames for the DQZ structure
impl DQZ<f32> {
    /* ABC */
    pub fn from_abc(a: f32, b: f32, c: f32, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS * (sin_cos.sin_val * a + left.sin_val * b + right.sin_val * c),
            q: TWO_THIRDS * (sin_cos.cos_val * a + left.cos_val * b + right.cos_val * c),
            z: ONE_THIRD * (a + b + c),
        }
    }
    pub fn to_abc(&self, sin_cos: &SinCos<f32>) -> ABC<f32> {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * self.d + sin_cos.cos_val * self.q + self.z,
              b: left.sin_val * self.d + left.cos_val * self.q + self.z, 
              c: right.sin_val * self.d + right.cos_val * self.q + self.z,
            }
    }
    /* Polar */
    pub fn to_polar(&self, sin_cos: &SinCos<f32>) -> Polar<f32> {
        if self.z != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the DQZ zero value.", self.z);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(self.d, 2.) + libm::powf(self.q, 2.)) / SQRT_2, 
                   theta: libm::atan2f(self.d * sin_cos.sin_val + self.q * sin_cos.cos_val, self.d * sin_cos.cos_val - self.q * sin_cos.sin_val)
            }
        }
    }
    pub fn from_polar(r: f32, theta: f32, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        let sin_cos_thetas = SinCos::<f32>::from_theta(theta + sin_cos.theta);
        DQZ { d: r * sin_cos_thetas.cos_val,
              q: r * sin_cos_thetas.sin_val, 
              z: 0. 
            }
    }
    /* AlphaBeta */
    pub fn to_ab(&self, sin_cos: &SinCos<f32>) -> AlphaBeta<f32> {
        return AlphaBeta { alpha: sin_cos.cos_val * self.d - sin_cos.sin_val * self.q, 
                           beta: sin_cos.sin_val * self.d + sin_cos.cos_val * self.q, 
                           gamma: self.z }
    }
    pub fn from_ab(alpha: f32, beta: f32, gamma: f32, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        let d = sin_cos.cos_val * alpha + sin_cos.sin_val * beta;
        let q = sin_cos.cos_val * beta - (sin_cos.sin_val * alpha);
        let z = gamma;
        DQZ{ d, q, z }
    }
    pub fn from_ab_(alpha: f32, beta: f32, sin_cos: &SinCos<f32>) -> DQZ<f32> {
        return DQZ::from_ab(alpha, beta, 0., sin_cos)
    }

}

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
    fn flip_theta(&self) -> SinCos<T>; 
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

    // Flip the sign of the reference angle, theta
    fn flip_theta(&self) -> SinCos<f32> {
        return SinCos{ sin_val: -self.sin_val,
                       cos_val: self.cos_val,
                       theta: -self.theta,        
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