// Structures for reference frame transformations
use super::constants::SQRT_2;
// use idsp::cossin;  // TODO: Use fixed point math?

pub struct AlphaBeta {
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
}

pub fn alpha_beta_fr_polar(v: f32, theta: f32) -> AlphaBeta {
    // let (cos_val, sin_val) = cossin(theta);
    let alpha = SQRT_2 * v * libm::cosf(theta);
    let beta = SQRT_2 * v * libm::sinf(theta);
    AlphaBeta {
        alpha,
        beta,
        gamma: 0.,
    }
}

pub fn alpha_beta_fr_ab(alpha: f32, beta: f32) -> AlphaBeta {
    AlphaBeta {
        alpha,
        beta,
        gamma: 0.,
    }
}

// TODO: Implement DQ reference frame & DQ reference frame power calculation? 