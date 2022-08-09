// Functions for standard calculations
use super::refs::*;

/// Calculates power from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_power(v: AlphaBeta, i: AlphaBeta) -> (f32, f32) {
    let p = 1.5 * (v.alpha * i.alpha + v.beta * i.beta);
    let q = 1.5 * (v.beta * i.alpha - v.alpha * i.beta);
    return (p, q)
}

/// Calculates active power, p, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_p(v: AlphaBeta, i: AlphaBeta) -> f32 {
    return 1.5 * (v.alpha * i.alpha + v.beta * i.beta)
}

/// Calculates the reactive power, q, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_q(v: AlphaBeta, i: AlphaBeta) -> f32 {
    return 1.5 * (v.beta * i.alpha - v.alpha * i.beta)
}

/// Calculates power from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_power(v: DQZ, i: DQZ) -> (f32, f32) {
    let p = 1.5 * (v.d * i.d + v.q * i.q);
    let q = 1.5 * (v.q * i.d - v.d * i.q);
    return (p, q)
}

/// Calculates active power, p, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_p(v: DQZ, i: DQZ) -> f32{
    return 1.5 * (v.d * i.d + v.q * i.q)
}

/// Calculates reactove power, q, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_q(v: DQZ, i: DQZ) -> f32 {
    return 1.5 * (v.q * i.d - v.d * i.q)
}

// Trying to implement add and multiply for arrays
// type Coord = [i32; 2];
// fn new_coord_from<F: Iterator<Item=i32>>(src: F) -> Coord {
//     let mut result = [0; 2];
//     for (rref, val) in result.iter_mut().zip(src) {
//         *rref = val;
//     }
//     result
// }

// fn add(a: Coord, b: Coord) -> Coord {
//     new_coord_from(a.iter().zip(&b).map(|(a, b)| a + b))
// }