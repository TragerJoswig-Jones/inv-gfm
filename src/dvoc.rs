/* Dispatchable Virtual Oscillator Controller (dVOC) implementation */
use super::calc::*;
use super::constants::*;
use super::refs::*;
use super::sims::*;
use super::*;

/* Define a dVOC controller */
const DVOC_STATES: usize = 2;
const DVOC_INPUTS: usize = 2;
type DvocStates<T> =  Vec<T, DVOC_STATES>;
pub struct DvocController<T: Num> {
    // Internal States
    pub v: T,  // voltage state (p.u.)
    pub theta: T, // angle state (p.u.)
    pub x: DvocStates<T>, // array of states; [v, theta]
    theta_idx: ThetaIdx, // index of theta value; 1 for dVOC states

    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    x_nom: T, // nominal voltage (p.u.)
    pub w_nom: T, // nominal frequency (rad)
    pub s_rated: T,  // maximum expected power output (VA)
    pub kv: T, // Base voltage (V)
    xi: T,
    c: T,  // Oscillator capacitance (F)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)
}

impl<T: Num> Dynamics<T, DVOC_STATES, DVOC_INPUTS> for DvocController<T> {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as an array of T fixed-point values: [v, theta]
    // * 'u' - alpha-beta current (p.u.) as an array of T fixed-point values: [ialpha, ibeta]
    fn dynamics(&self, x:  &DvocStates<T>, u: [T; DVOC_INPUTS]) -> DvocStates<T> {
        let (v, theta) = (x[0], x[1] * self.w_nom);
        let v_dq = DQZ{ d: v * T::from_fixed(SQRT_2), q: T::from_fixed(ZERO), z: T::from_fixed(ZERO)};  // TODO: Determine the best way to handle multiplying by a constant
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<T>::from_theta(theta));
        let (p, q) = calc_dq_power(v_dq, i_dq);

        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = T::from_fixed(ONE) / (T::from_fixed(SQRT_2) * self.c * x[0]);
        let dv_dt = T::from_fixed(TWO) * self.xi * x[0] * ((self.x_nom * self.x_nom ) - (x[0] * x[0])) - _sqrt2cv * (q - self.q_ref);
        let dtheta_dt = T::from_fixed(ONE) - _sqrt2cv / (x[0] * self.w_nom) * (p - self.p_ref); 
        return na::Vector2::new(dv_dt, dtheta_dt)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, DVOC_STATES, DVOC_INPUTS> for DvocController<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, DVOC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, DVOC_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

impl<T: Num> DvocController<T> {
    // Sets the reference active power within the dVOC controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in Watts
    pub fn set_p_ref(&mut self, p_ref: T) {
        self.p_ref = p_ref / self.s_rated;
    }
}

pub fn build_dvoc_controller<T: Num>(v_nom: T, w_nom: T, s_rated: T, xi: T, c: T) -> DvocController<T> {
    DvocController {
        v_nom,
        x_nom: T::from_fixed(ONE),
        w_nom,
        s_rated,
        v: T::from_fixed(ONE),  
        theta: T::from_fixed(ZERO),
        x: na::Vector2::new(T::from_fixed(ONE), T::from_fixed(ZERO)),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        kv: v_nom,
        xi,
        c,
        p_ref: T::from_fixed(ZERO),
        q_ref: T::from_fixed(ZERO),
    }
}

pub fn build_default_dvoc_controller<T: Num>(v_nom: T, f_nom: T) -> DvocController<T> {
    DvocController {
        v_nom,
        x_nom: T::from_fixed(ONE),
        w_nom: T::from_fixed(TWO*PI) * f_nom,
        s_rated: T::from_num(1000),
        v: T::from_fixed(ONE),  
        theta: T::from_fixed(ZERO),
        x: na::Vector2::new(T::from_fixed(ONE), T::from_fixed(ZERO)),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        kv: v_nom,
        xi: T::from_num(15),
        c: T::from_num(0.2679),
        p_ref: T::from_fixed(ZERO),
        q_ref: T::from_fixed(ZERO),
    }
}

pub fn build_dvoc_controller_from_flt<T: Num>(v_nom: f32, f_nom: f32, s_rated: f32, xi: f32, c: f32) -> DvocController<T> {
    let w_nom = T::from_num(2. * f_nom) * T::from_fixed(PI);
    DvocController {
        v_nom: T::from_num(v_nom),
        x_nom: T::from_fixed(ONE),
        w_nom,
        s_rated: T::from_num(s_rated),
        v: T::from_fixed(ONE),  
        theta: T::from_fixed(ZERO),
        x: na::Vector2::new(T::from_fixed(ONE), T::from_fixed(ZERO)),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        kv: T::from_num(v_nom),
        xi: T::from_num(xi),
        c: T::from_num(c),
        p_ref: T::from_fixed(ZERO),
        q_ref: T::from_fixed(ZERO),
    }
}