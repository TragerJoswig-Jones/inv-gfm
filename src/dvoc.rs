use super::calc::*;
use super::constants::*;
use super::refs::*;
use super::sims::*;
use super::*;

/* Dispatchable Virtual Oscillator Controller (dVOC) implementation */
const DVOC_STATES: usize = 2;
const DVOC_INPUTS: usize = 2;
type DvocStates<T> =  Vec<T, DVOC_STATES>;
/* Define a dVOC controller */
pub struct DvocController<T: Num> {
    // Internal States
    pub v: T,  // voltage state (p.u.)
    pub theta: T, // angle state (p.u.)
    pub x: DvocStates<T>, // array of states; [v, theta]
    theta_idx: ThetaIdx, // index of theta value; 1 for dVOC states

    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    x_nom: T, // nominal voltage (p.u.)
    pub w_nom: T, // nominal frequency (rad/s)
    pub kv: T, // Base voltage (V)
    xi: T,
    c: T,  // Oscillator capacitance (F)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)
}

impl Dynamics<f32, DVOC_STATES, DVOC_INPUTS> for DvocController<f32> {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as an array of T fixed-point values: [v, theta]
    // * 'u' - alpha-beta current (p.u.) as an array of T fixed-point values: [ialpha, ibeta]
    fn dynamics(&self, x:  &DvocStates<f32>, u: [f32; DVOC_INPUTS]) -> DvocStates<f32> {
        let (v, theta) = (x[0], x[1] * self.w_nom);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0. };  // TODO: Determine the best way to handle multiplying by a constant
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(v_dq, i_dq);

        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = 1. / (SQRT_2 * self.c * x[0]);
        let dv_dt = 2. * self.xi * x[0] * ((self.x_nom * self.x_nom ) - (x[0] * x[0])) - _sqrt2cv * (q - self.q_ref);
        let dtheta_dt = 1. - _sqrt2cv / x[0] / self.w_nom * (p - self.p_ref); 
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
    // * 'p_ref' - The desired active power reference in p.u.
    pub fn set_p_ref(&mut self, p_ref: T) {
        self.p_ref = p_ref;
    }
}

pub fn build_dvoc_controller(v_nom: f32, w_nom: f32, xi: f32, c: f32) -> DvocController<f32> {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom,
        v: 1.,  
        theta: 0.,
        x: na::Vector2::new(1., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        kv: v_nom,
        xi,
        c,
        p_ref: 0.,
        q_ref: 0.,
    }
}

pub fn build_default_dvoc_controller(v_nom: f32, f_nom: f32) -> DvocController<f32> {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom: 2.*PI * f_nom,
        v: 1.,  
        theta: 0.,
        x: na::Vector2::new(1., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        kv: v_nom,
        xi: 15.,
        c: 0.2679,
        p_ref: 0.,
        q_ref: 0.,
    }
}