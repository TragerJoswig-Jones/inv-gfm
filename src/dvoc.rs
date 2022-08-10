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
    ki: T, // Base current (A)
    pub kv: T, // Base voltage (V)
    xi: T,
    c: T,  // Oscillator capacitance (F)
    l: T,  // Oscillator inductance (H)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)
}

impl Dynamics<Flt, DVOC_STATES, DVOC_INPUTS> for DvocController<Flt> {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of Flt values: (v, theta)
    // * 'u' - alpha-beta current (A) as a tuple of Flt values: (ialpha, ibeta)
    fn dynamics(&self, x:  &DvocStates<Flt>, u: [Flt; DVOC_INPUTS]) -> DvocStates<Flt> {
        let (v, theta) = (x[0], x[1]);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::<Flt>::from_ab_(u[0], u[1]).to_dqz(SinCos::<Flt>::from_theta(theta));
        let (p, q) = calc_dq_power(v_dq, i_dq);  // TODO: Determine is this calculation can be done in p.u.

        // Unit dynamics (eq.11-12 from 'A Grid-compatible Virtual Oscillator Controller' by Lu M., Et al.)
        let kvki_3cv = self.kv * self.ki / (3. * self.c * v);
        let dv_dt = self.xi / (self.kv * self.kv) * v * (2. * (self.v_nom * self.v_nom ) - 2. * (v * v)) - kvki_3cv * (q - self.q_ref);
        let dtheta_dt = self.w_nom - kvki_3cv / v * (p - self.p_ref);
        return na::Vector2::new(dv_dt, dtheta_dt)
    }
}

// impl Dynamics<Fxd, DVOC_STATES, DVOC_INPUTS> for DvocController<Fxd> {
//     // Calculates the voltage dynamics of the dVOC controller using the given input, u.
//     // # Arguments    
//     // * 'x' - polar voltage (p.u.) as an array of fixed-point values: [v, theta]
//     // * 'u' - alpha-beta current (A) as an array of fixed-point values: [ialpha, ibeta]
//     fn dynamics(&self, x:  &DvocStates<Fxd>, u: [Fxd; DVOC_INPUTS]) -> DvocStates<Fxd> {
//         let (v, theta) = (x[0], x[1]);
//         let v_dq = DQZ{ d: v * SQRT_2_FXD, q: ZERO, z: ZERO};
//         let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<Fxd>::from_theta(theta));
//         let (p, q) = calc_dq_power(v_dq, i_dq);  // TODO: Determine is this calculation can be done in p.u.

//         // Unit dynamics (eq.11-12 from 'A Grid-compatible Virtual Oscillator Controller' by Lu M., Et al.)
//         let kvki_3cv = self.kv * self.ki / (3 * self.c * v);
//         let dv_dt = self.xi / (self.kv * self.kv) * v * (2 * (self.v_nom * self.v_nom ) - 2 * (v * v)) - kvki_3cv * (q - self.q_ref);
//         let dtheta_dt = self.w_nom - kvki_3cv / v * (p - self.p_ref);
//         return na::Vector2::new(dv_dt, dtheta_dt)
//     }
// }

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
}

impl<T: Num> DvocController<T> {
    // Sets the reference active power within the dVOC controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in Watts
    pub fn set_p_ref(&mut self, p_ref: T) {
        self.p_ref = p_ref;
    }
}

pub fn build_dvoc_controller(v_nom: Flt, w_nom: Flt, s_rated: Flt, xi: Flt, c: Flt) -> DvocController<Flt> {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom,
        s_rated,
        v: v_nom,  
        theta: 0.,
        x: na::Vector2::new(v_nom, 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        ki: 3. * v_nom / s_rated,
        kv: v_nom,
        xi,
        c,
        l: 1. / (w_nom * w_nom * c),
        p_ref: 0.,
        q_ref: 0.,
    }
}

pub fn build_default_dvoc_controller(v_nom: Flt, f_nom: Flt) -> DvocController<Flt> {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom: 2. * PI * f_nom,
        s_rated: 1000.,
        v: v_nom,  
        theta: 0.,
        x: na::Vector2::new(v_nom, 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
        ki: 3. * v_nom / 1000.,
        kv: v_nom,
        xi: 15.,
        c: 0.2679,
        l: 1. / (4. * PI * PI * 3600. * 0.2679),
        p_ref: 0.,
        q_ref: 0.,
    }
}
