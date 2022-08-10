/* Droop controller implementation */
use super::calc::*;
use super::constants::{SQRT_2, SQRT_3, PI};
use super::refs::*;
use super::sims::*;
use super::*;

/* Define a droop controller */
const DROOP_STATES: usize = 4;
const DROOP_INPUTS: usize = 2;
type DroopStates<T: Num> =  Vec<T, DROOP_STATES>;
pub struct DroopController<T: Num> {
    // Internal States
    pub v: T,  // voltage state (p.u.)
    pub theta: T, // angle state (p.u.)
    pub p_filt: T, // low-pass filter active power (p.u.)
    pub q_filt: T, // low-pass filter reactive power (p.u.)
    pub x: DroopStates<T>,
    theta_idx: ThetaIdx,
    
    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    x_nom: T, // nominal voltage (p.u.)
    pub w_nom: T, // nominal frequency (rad)
    pub w_c: T, // power low-pass filter cutoff frequency (rad)
    pub s_rated: T,  // maximum expected power output (VA)
    pub mp: T, // frequency droop slope (rad/s)
    pub mq: T, // voltage droop slope (V)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)
}

impl<T: Num> Dynamics<T, DROOP_STATES, DROOP_INPUTS> for DroopController<T> {
    // Calculates the voltage dynamics of the droop controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) and filtered powers as a tuple of f32 values: (v, theta, p_filt, q_filt)
    // * 'u' - alpha-beta current (A) as a tuple of f32 values: (ialpha, ibeta)
    fn dynamics(&self, x: &DroopStates<T>, u: [T; DROOP_INPUTS]) -> DroopStates<T> {
        let (v, theta, p_filt, q_filt) = (x[0], x[1], x[2], x[3]);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<T>::from_theta(theta));
        let (p, q) = calc_dq_power(v_dq, i_dq);  // TODO: Determine is this calculation can be done in p.u.

        // Unit dynamics (eq.13 & eq.17 from 'Control of Parallel Connected Inverters in Standalone ac Supply Systems' by Chandorkar M., Et al.)
        let dp_filt_dt = self.w_c * (p - p_filt);
        let dq_filt_dt = self.w_c * (q - q_filt);
        let dv_dt = - self.mq * dq_filt_dt;
        let dtheta_dt = self.w_nom - self.mp * (p_filt - self.p_ref);
        return na::Vector4::new(dv_dt, dtheta_dt, dp_filt_dt, dq_filt_dt)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, DROOP_STATES, DROOP_INPUTS> for DroopController<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, DROOP_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, DROOP_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
}

impl<T: Num> DroopController<T> {
    // Sets the reference active power within the droop controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in Watts
    pub fn set_p_ref(&mut self, p_ref: T) {
        self.p_ref = p_ref;
    }
}

pub fn build_droop_controller(v_nom: f32, w_nom: f32, s_rated: f32, mp: f32, mq: f32, w_c: f32) -> DroopController<f32> {
    DroopController {
        v_nom,
        x_nom: 1.,
        w_nom,
        s_rated,
        v: v_nom,  
        theta: 0.,
        p_filt: 0.,
        q_filt: 0.,
        x: na::Vector4::new(v_nom, 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        mp,
        mq,
        w_c,
        p_ref: 0.,
        q_ref: 0.,
    }
}

pub fn build_default_droop_controller(v_nom: f32, f_nom: f32) -> DroopController<f32> {
    DroopController {
        v_nom,
        x_nom: 1.,
        w_nom: 2. * PI * f_nom,
        s_rated: 1000.,
        v: v_nom,  
        theta: 0.,
        p_filt: 0.,
        q_filt: 0.,
        x: na::Vector4::new(v_nom, 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        mp: 0.0026,
        mq: 0.005,
        w_c: 2.*PI*30.,
        p_ref: 0.,
        q_ref: 0.,
    }
}
