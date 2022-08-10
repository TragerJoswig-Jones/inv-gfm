// Functions for standard calculations
use super::refs::*;

/// Calculates power from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_power<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> (T, T) {
    let p = 3 / 2 * (v.alpha * i.alpha + v.beta * i.beta);
    let q = 3 / 2 * (v.beta * i.alpha - v.alpha * i.beta);
    return (p, q)
}

/// Calculates active power, p, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_p<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> T {
    return 3 / 2 * (v.alpha * i.alpha + v.beta * i.beta)
}

/// Calculates the reactive power, q, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_q<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> T {
    return 3 / 2 * (v.beta * i.alpha - v.alpha * i.beta)
}

/// Calculates power from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_power<T: Num>(v: DQZ<T>, i: DQZ<T>) -> (T, T) {
    let p = 3 / 2 * (v.d * i.d + v.q * i.q);
    let q = 3 / 2 * (v.q * i.d - v.d * i.q);
    return (p, q)
}

/// Calculates active power, p, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_p<T: Num>(v: DQZ<T>, i: DQZ<T>) -> T{
    return 3 / 2 * (v.d * i.d + v.q * i.q)
}

/// Calculates reactove power, q, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_q<T: Num>(v: DQZ<T>, i: DQZ<T>) -> T {
    return 3 / 2 * (v.q * i.d - v.d * i.q)
}
