// Functions for standard calculations
use super::reference_frames::*;

/* FLOATING-POINT IMPLEMENTATIONS */
/// Calculates power from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_power(v: &AlphaBeta<f32>, i: &AlphaBeta<f32>, n_phase: f32) -> (f32, f32) {
    let p = n_phase * (v.alpha * i.alpha + v.beta * i.beta) / 2.;
    let q = n_phase * (v.beta * i.alpha - v.alpha * i.beta) / 2.;
    return (p, q)
}

/// Calculates active power, p, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_p(v: &AlphaBeta<f32>, i: &AlphaBeta<f32>, n_phase: f32) -> f32 {
    return n_phase * (v.alpha * i.alpha + v.beta * i.beta) / 2.
}

/// Calculates the reactive power, q, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_q(v: &AlphaBeta<f32>, i: &AlphaBeta<f32>, n_phase: f32) -> f32 {
    return n_phase * (v.beta * i.alpha - v.alpha * i.beta) / 2.
}

/// Calculates power from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_power(v: &DQZ<f32>, i: &DQZ<f32>, n_phase: f32) -> (f32, f32) {
    let p = n_phase * (v.d * i.d + v.q * i.q) / 2.;
    let q = n_phase * (v.q * i.d - v.d * i.q) / 2.;
    return (p, q)
}

/// Calculates active power, p, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_p(v: &DQZ<f32>, i: &DQZ<f32>, n_phase: f32) -> f32{
    return n_phase * (v.d * i.d + v.q * i.q) / 2.
}

/// Calculates reactove power, q, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_q(v: &DQZ<f32>, i: &DQZ<f32>, n_phase: f32) -> f32 {
    return n_phase * (v.q * i.d - v.d * i.q) / 2.
}