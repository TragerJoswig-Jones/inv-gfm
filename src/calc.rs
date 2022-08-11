// Functions for standard calculations
use super::refs::*;
use super::constants::*;
use super::*;

/* FIXED-POINT IMPLEMENTATIONS */
/// Calculates power from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_power<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> (T, T) {
    let p = (v.alpha * i.alpha + v.beta * i.beta) * T::from_num(1.5);
    let q = (v.beta * i.alpha - v.alpha * i.beta) * T::from_num(1.5);
    return (p, q)
}

/// Calculates active power, p, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_p<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> T {
    return (v.alpha * i.alpha + v.beta * i.beta) * T::from_num(1.5)
}

/// Calculates the reactive power, q, from the given alpha-beta voltage, v, and alpha-beta current, i
pub fn calc_ab_q<T: Num>(v: AlphaBeta<T>, i: AlphaBeta<T>) -> T {
    return (v.beta * i.alpha - v.alpha * i.beta) * T::from_num(1.5)
}

/// Calculates power from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_power<T: Num>(v: DQZ<T>, i: DQZ<T>) -> (T, T) {
    let p = (v.d * i.d + v.q * i.q) * T::from_num(1.5);
    let q = (v.q * i.d - v.d * i.q) * T::from_num(1.5);
    return (p, q)
}

/// Calculates active power, p, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_p<T: Num>(v: DQZ<T>, i: DQZ<T>) -> T{
    return (v.d * i.d + v.q * i.q) * T::from_num(1.5)
}

/// Calculates reactove power, q, from the given DQZ voltage, v, and DQZ current, i
pub fn calc_dq_q<T: Num>(v: DQZ<T>, i: DQZ<T>) -> T {
    return (v.q * i.d - v.d * i.q) * T::from_num(1.5)
}
