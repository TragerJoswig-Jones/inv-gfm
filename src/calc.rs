// Functions for standard calculations
use super::refs::*;

/* FLOATING-POINT IMPLEMENTATIONS */
/// Calculates power from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_power(v: AlphaBeta<f32>, i: AlphaBeta<f32>) -> (f32, f32) {
    let p = 3. / 2. * (v.alpha * i.alpha + v.beta * i.beta);
    let q = 3. / 2. * (v.beta * i.alpha - v.alpha * i.beta);
    return (p, q)
}

/// Calculates active power, p, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_p(v: AlphaBeta<f32>, i: AlphaBeta<f32>) -> f32 {
    return 3. / 2. * (v.alpha * i.alpha + v.beta * i.beta)
}

/// Calculates the reactive power, q, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_q(v: AlphaBeta<f32>, i: AlphaBeta<f32>) -> f32 {
    return 3. / 2. * (v.beta * i.alpha - v.alpha * i.beta)
}

/// Calculates power from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_power(v: DQZ<f32>, i: DQZ<f32>) -> (f32, f32) {
    let p = 3. / 2. * (v.d * i.d + v.q * i.q);
    let q = 3. / 2. * (v.q * i.d - v.d * i.q);
    return (p, q)
}

/// Calculates active power, p, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_p(v: DQZ<f32>, i: DQZ<f32>) -> f32{
    return 3. / 2. * (v.d * i.d + v.q * i.q)
}

/// Calculates reactove power, q, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_q(v: DQZ<f32>, i: DQZ<f32>) -> f32 {
    return 3. / 2. * (v.q * i.d - v.d * i.q)
}