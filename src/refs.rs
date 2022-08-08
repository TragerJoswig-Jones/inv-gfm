/* Structures for reference frame transformations */
use super::constants::*;
// use idsp::cossin;  // TODO: Determine if we want to use fixed point math?

/*
ABC Three-Phase Values 
*/
pub struct ABC {
    pub a: f32,
    pub b: f32,
    pub c: f32,
}
pub trait ToFromABC {
    fn to_abc(&self) -> ABC;
    fn from_abc(a: f32, b: f32, c: f32) -> Self;
}
impl ToFromAlphaBeta for ABC {
    fn to_ab(&self) -> AlphaBeta {
        let alpha = TWO_THIRDS * self.a - ONE_THIRD * self.b - ONE_THIRD * self.c;
        let beta = SQRT_3_OVER_3 * (self.b - self.c);
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
impl ToFromDQZ for ABC {
    fn to_dqz(&self, sin_cos: SinCos) -> DQZ {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_right_120();
        DQZ {
            d: TWO_THIRDS * (sin_cos.sin_val * self.a + left.sin_val * self.b + right.sin_val * self.c),
            q: TWO_THIRDS * (sin_cos.cos_val * self.a + left.cos_val * self.b + right.cos_val * self.c),
            z: ONE_THIRD * (self.a + self.b + self.c),
        }
    }
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: SinCos) -> Self {
        let left = sin_cos.rotate_left_120();
        let right = sin_cos.rotate_left_120();
        ABC { a: sin_cos.sin_val * d + sin_cos.cos_val * q + z,
              b: left.sin_val * d + left.cos_val * q + z, 
              c: right.sin_val * d + right.cos_val * q + z,
            }
    }
}
impl ToFromPolar for ABC {
    fn to_polar(&self) -> Polar {
        // Implemented as ABC to AB then AB to Polar, but throws an error if gamme != 0 (This indicates the ABC signal is not balanced)
        let ab = self.to_ab();
        if ab.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", ab.gamma);
        } else {
            ab.to_polar()
        }
    }
    fn from_polar(r: f32, theta: f32) -> Self {
        let sin_cos = SinCos::from_theta(theta);
        ABC{ a: r * sin_cos.sin_val, 
             b: r * sin_cos.rotate_left_120().sin_val, 
             c: r * sin_cos.rotate_right_120().sin_val} // TODO: double check the left/right diretions here
    }
}

/* 
Polar Coordinate Values 
*/
pub struct Polar {
    pub r: f32,
    pub theta: f32,
}
pub trait ToFromPolar {
    fn to_polar(&self) -> Polar;
    fn from_polar(r: f32, theta: f32) -> Self;
}
impl ToFromABC for Polar {
    fn to_abc(&self) -> ABC {
        ABC::from_polar(self.r, self.theta)
    }
    fn from_abc(a: f32, b: f32, c: f32) -> Self {
        // Creates ABC object and calls its to_polar function
        // TODO: Determine if this should avoid the intermediate ABC object creation
        let abc: ABC = ABC { a, b, c };
        abc.to_polar()
    }
}
impl ToFromAlphaBeta for Polar {
    fn to_ab(&self) -> AlphaBeta {
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
impl ToFromDQZ for Polar {
    fn to_dqz(&self, sin_cos: SinCos) -> DQZ {
        let sin_cos_thetas = SinCos::from_theta(self.theta + sin_cos.theta);
        DQZ { d: self.r * sin_cos_thetas.cos_val,
              q: self.r * sin_cos_thetas.sin_val, 
              z: 0. 
            }
    }
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: SinCos) -> Self {
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
pub struct AlphaBeta {
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
}
pub trait ToFromAlphaBeta {
    fn to_ab(&self) -> AlphaBeta;
    fn from_ab(alpha: f32, beta: f32, gamma: f32) -> Self;
}
impl AlphaBeta {
    pub fn from_ab_(alpha: f32, beta: f32) -> AlphaBeta {
        AlphaBeta{
            alpha,
            beta,
            gamma: 0.,
        }
    }
}
impl ToFromABC for AlphaBeta {
    fn to_abc(&self) -> ABC {
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
impl ToFromPolar for AlphaBeta {
    fn to_polar(&self) -> Polar {
        if self.gamma != 0. {
            panic!("Cannot create a Polar value from unbalanced three-phase ABC values! Got {} as the AlphaBeta gamma value.", self.gamma);
        } else {
            Polar{ r: libm::sqrtf(libm::powf(self.alpha, 2.) + libm::powf(self.beta, 2.)) / SQRT_2, 
                theta: libm::atan2f(self.beta, self.alpha)}
        }
    }
    fn from_polar(r: f32, theta: f32) -> Self {
        AlphaBeta {
            alpha: SQRT_2 * r * libm::cosf(theta),
            beta: SQRT_2 * r * libm::sinf(theta),
            gamma: 0.,
        }
    }
}
impl ToFromDQZ for AlphaBeta {
    fn to_dqz(&self, sin_cos: SinCos) -> DQZ {
        let d = sin_cos.cos_val * self.alpha + sin_cos.sin_val * self.beta;
        let q = -sin_cos.sin_val * self.alpha + sin_cos.cos_val * self.beta;
        let z = self.gamma;
        DQZ{ d, q, z }
    }
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: SinCos) -> Self {
        return AlphaBeta { alpha: sin_cos.cos_val * d - sin_cos.sin_val * q, 
                           beta: sin_cos.sin_val * d + sin_cos.cos_val * q, 
                           gamma: z }
    }
}

/* 
Direct-Quadrature-Zero (DQZ) Reference Frame Values 
*/
pub struct DQZ {
    pub d: f32,
    pub q: f32,
    pub z: f32,
}
pub trait ToFromDQZ {
    fn to_dqz(&self, sin_cos: SinCos) -> DQZ;
    fn from_dqz(d: f32, q: f32, z: f32, sin_cos: SinCos) -> Self;
}
// TODO: Should we implement ToFrom functions for the DQZ structure? 
//       The issue is that DQZ transforms require a reference angle.
//       If not constructing a DQZ struct will have to be done from 
//       an existing reference frame struct.

/* 
Stores sine and cosine pair values 
*/
pub struct SinCos {
    sin_val: f32,
    cos_val: f32,
    theta: f32,
}
impl SinCos {
    pub fn from_theta(theta: f32) -> Self {
        let sin_val = libm::sinf(theta);
        let cos_val = libm::cosf(theta);
        Self{ sin_val, cos_val, theta }
    }

    // Rotate the reference angle, theta, by 120 degrees counter-clockwise
    pub fn rotate_right_120(&self) -> Self {
        return SinCos { sin_val: -ONE_HALF * self.sin_val + SQRT_3_OVER_2 * self.cos_val,
                        cos_val: -ONE_HALF * self.cos_val - SQRT_3_OVER_2 * self.sin_val,
                        theta: self.theta + TWO_PI_OVER_THREE,
                    }
    }

    // Rotate the reference angle, theta, by 120 degrees clockwise
    pub fn rotate_left_120(&self) -> Self {
        return SinCos{ sin_val: -ONE_HALF * self.sin_val - SQRT_3_OVER_2 * self.cos_val, 
                       cos_val: -ONE_HALF * self.cos_val + SQRT_3_OVER_2 * self.sin_val,
                       theta: self.theta - TWO_PI_OVER_THREE,
                    }
    }

    // Get the sin_val assosiated with the SinCos object
    pub fn sin_val(&self) -> f32 {
        self.sin_val
    }

    // Get the cos_val assosiated with the SinCos object
    pub fn cos_val(&self) -> f32 {
        self.cos_val
    }

    // Get the theta assosiated with the SinCos object
    pub fn theta(&self) -> f32 {
        self.theta
    }
}
