use super::constants::{SQRT_2, PI};
use super::refs::{alpha_beta_fr_polar, alpha_beta_fr_ab};

// TODO: Determine what should be in per unit and how to handle unit conversions

/* Implement dynamics for all available elements */
pub trait Dynamics<const X: usize, const U: usize>{
    fn dynamics(&self, x: [f32; X], u: [f32; U]) -> (f32, f32);
    fn step(&mut self, dt: f32, u: [f32; U]) -> ();
}

/* Implement a dynamic step function with no inputs */
pub trait NoInputStep<const X: usize, const U: usize>{
    fn step_(&mut self, dt: f32) -> ();
}

impl<T, const X: usize, const U: usize> NoInputStep<X, U> for T where T: Dynamics<X, U>{
    // Steps the dVOC dynamics
    // # Arguments
    // * 'dt' - The step period in seconds (1 / fs)
    // * 'u' - An empty array as there are no inputs
    fn step_(&mut self, dt: f32) {
        self.step(dt, [0.; U]);
    }
}

/* Define an RL Filter object */
const LINE_STATES: usize = 2;
const LINE_INPUTS: usize = 4;
pub struct RLFilter {
    // RL Filter Parameters
    pub w_nom: f32, // nominal frequency (rad)
    pub s_rated: f32,  // maximum expected power output (VA)
    rf: f32,  // Line resistance (Ohms)
    lf: f32,  // Line inductance (H)

    // Internal States
    pub i_alpha: f32,  // alpha current state (p.u.)
    pub i_beta: f32, // beta current state (p.u.)
}

impl Dynamics<LINE_STATES, LINE_INPUTS> for RLFilter {
    // Steps the dVOC dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'dt' - The step period in seconds (1 / fs)
    // * 'u' - input voltages as an array of f32 values: (v1, theta1, v2, theta2)
    fn step(&mut self, dt: f32, u: [f32; LINE_INPUTS]) {
        let (di_alpha_dt1, di_beta_dt1) = self.dynamics([self.i_alpha, self.i_beta], u);
        let (di_alpha_dt2, di_beta_dt2) = self.dynamics([self.i_alpha + dt * di_alpha_dt1, self.i_beta + dt * di_beta_dt1], u);
        let di_alpha_dt = (di_alpha_dt1 + di_alpha_dt2) * 0.5;
        let di_beta_dt = (di_beta_dt1 + di_beta_dt2) * 0.5;
        self.i_alpha = self.i_alpha + dt * di_alpha_dt;
        self.i_alpha = (self.i_beta + dt * di_beta_dt) % 1.;
    }

    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of f32 values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of f32 values: (v1, theta1, v2, theta2)
    fn dynamics(&self, x: [f32; LINE_STATES], u: [f32; LINE_INPUTS]) -> (f32, f32) {
        let (v1, theta1, v2, theta2) = (u[0], u[1], u[2], u[3]);
        let v1_ab = alpha_beta_fr_polar(v1, theta1);
        let v2_ab = alpha_beta_fr_polar(v2, theta2);
        let i_ab = alpha_beta_fr_ab(x[0], x[1]);
        
        // Unitc line dynamics
        let di_alpha_dt = 1./self.lf*(v1_ab.alpha - v2_ab.alpha - self.rf*i_ab.alpha);
        let di_beta_dt = 1./self.lf*(v1_ab.beta - v2_ab.beta - self.rf*i_ab.beta);

        (di_alpha_dt, di_beta_dt)
    }
}

pub fn build_rl_line(w_nom: f32, s_rated: f32, rf: f32, lf: f32) -> RLFilter {
    RLFilter {
        // RL Filter Parameters
        w_nom,
        s_rated,
        rf,
        lf,

        // Internal States
        i_alpha: 0.,
        i_beta: 0.,
    }
}

/* Define a stiff AC voltage source object */
const ACVS_STATES: usize = 2;
const ACVS_INPUTS: usize = 0;
pub struct ACVoltSrc {
    // Parameters
    pub v_nom: f32, // nominal RMS LN voltage (V)
    pub w_nom: f32, // nominal frequency (rad)
    pub s_rated: f32,  // maximum expected power output (VA)
    
    
    // Internal States
    pub v: f32,  // alpha current state (p.u.)
    pub theta: f32, // beta current state (p.u.)
}

impl Dynamics<ACVS_STATES, ACVS_INPUTS> for ACVoltSrc {
    // Steps the dVOC dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'dt' - The step period in seconds (1 / fs)
    // * 'u' - An empty array as there are no inputs
    fn step(&mut self, dt: f32, u: [f32; ACVS_INPUTS]) {
        let (dv_dt1, dtheta_dt1) = self.dynamics([self.v, self.theta], u);
        let (dv_dt2, dtheta_dt2) = self.dynamics([self.v + dt * dv_dt1, self.theta + dt * dtheta_dt1], u);
        let dv_dt = (dv_dt1 + dv_dt2) * 0.5;
        let dtheta_dt = (dtheta_dt1 + dtheta_dt2) * 0.5;
        self.v = self.v + dt * dv_dt;
        self.theta = (self.theta + dt * dtheta_dt) % (2.*PI);
    }

    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of f32 values: (v, theta)
    // * 'u' - An empty array as there are no inputs
    fn dynamics(&self, _x: [f32; ACVS_STATES], _u: [f32; ACVS_INPUTS]) -> (f32, f32) {
        return (0., self.w_nom)
    }
}

pub fn build_ac_volt_src(v_nom: f32, w_nom: f32, s_rated: f32) -> ACVoltSrc {
    ACVoltSrc {
        // Parameters
        v_nom,
        w_nom,
        s_rated,
        
        // Internal States
        v: v_nom,
        theta: 0.,
    }
}

