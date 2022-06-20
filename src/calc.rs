// Functions for standard calculations
use super::refs::AlphaBeta;

/// Calculates power from the given alpha-beta voltage, x, and alpha-beta current, u
pub fn calc_power(v: AlphaBeta, i: AlphaBeta) -> (f32, f32) {
    let p = 1.5 * (v.alpha * i.alpha + v.beta * i.beta);
    let q = 1.5 * (v.beta * i.alpha - v.alpha * i.beta);
    return (p, q)
}