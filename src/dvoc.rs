/* Dispatchable Virtual Oscillator Controller (dVOC) implementation */
use super::calc::*;
use super::constants::*;
use super::refs::*;
use super::sims::*;

/* Define a dVOC controller */
const DVOC_STATES: usize = 2;
const DVOC_INPUTS: usize = 2;
pub struct DvocController {
    // Internal States
    pub v: f32,  // voltage state (p.u.)
    pub theta: f32, // angle state (p.u.)
    
    // Other Parameters
    pub v_nom: f32, // nominal voltage (V)
    x_nom: f32, // nominal voltage (p.u.)
    pub w_nom: f32, // nominal frequency (rad)
    pub s_rated: f32,  // maximum expected power output (VA)
    ki: f32, // Base current (A)
    pub kv: f32, // Base voltage (V)
    xi: f32,
    c: f32,  // Oscillator capacitance (F)
    l: f32,  // Oscillator inductance (H)
    pub p_ref: f32,  // Active power reference (p.u.)
    pub q_ref: f32,  // Reactive power reference (p.u.)
}

impl Dynamics<DVOC_STATES, DVOC_INPUTS> for DvocController {
    // Steps the dVOC dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'u' - alpha-beta current as a tuple of f32 values: (ialpha, ibeta)
    fn step(&mut self, dt: f32, u: [f32; DVOC_INPUTS]) {
        let dx_dt1 = self.dynamics([self.v, self.theta], u);
        let dx_dt2 = self.dynamics([self.v + dt * dx_dt1[0], self.theta + dt * dx_dt1[1]], u);
        let dv_dt = (dx_dt1[0] + dx_dt2[0]) * 0.5;
        let dtheta_dt = (dx_dt1[1] + dx_dt2[1]) * 0.5;
        self.v = self.v + dt * dv_dt;
        self.theta = (self.theta + dt * dtheta_dt) % (2.*PI);
    }

    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of f32 values: (v, theta)
    // * 'u' - alpha-beta current (A) as a tuple of f32 values: (ialpha, ibeta)
    fn dynamics(&self, x: [f32; DVOC_STATES], u: [f32; DVOC_INPUTS]) -> [f32; DVOC_STATES] {
        let (v, theta) = (x[0], x[1]);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::from_theta(theta));
        let (p, q) = calc_dq_power(v_dq, i_dq);  // TODO: Determine is this calculation can be done in p.u.

        // Unit dynamics (eq.11-12 from 'A Grid-compatible Virtual Oscillator Controller' by Lu M., Et al.)
        let kvki_3cv = self.kv * self.ki / (3. * self.c * v);
        let dv_dt = self.xi / (self.kv * self.kv) * v * (2. * (self.v_nom * self.v_nom ) - 2. * (v * v)) - kvki_3cv * (q - self.q_ref);
        let dtheta_dt = self.w_nom - kvki_3cv / v * (p - self.p_ref);
        return [dv_dt, dtheta_dt]
    }
}

impl DvocController {
    // Sets the reference active power within the dVOC controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in Watts
    pub fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
}

pub fn build_dvoc_controller(v_nom: f32, w_nom: f32, s_rated: f32, xi: f32, c: f32) -> DvocController {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom,
        s_rated,
        v: v_nom,  
        theta: 0.,
        ki: 3. * v_nom / s_rated,
        kv: v_nom,
        xi,
        c,
        l: 1. / (w_nom * w_nom * c),
        p_ref: 0.,
        q_ref: 0.,
    }
}

pub fn build_default_dvoc_controller(v_nom: f32, f_nom: f32) -> DvocController {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom: 2. * PI * f_nom,
        s_rated: 1000.,
        v: v_nom,  
        theta: 0.,
        ki: 3. * v_nom / 1000.,
        kv: v_nom,
        xi: 15.,
        c: 0.2679,
        l: 1. / (4. * PI * PI * 3600. * 0.2679),
        p_ref: 0.,
        q_ref: 0.,
    }
}
