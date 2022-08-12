/* Dynamical objects for simulating power systems */
use crate::constants::*;
use crate::refs::*;
use crate::*;

// TODO: Determine what should be in per unit and how to handle unit conversions

/* Dynamics trait that should be implemented for all available elements */
// * 'X' - Number of states, 
// * 'U' - Number of inputs
pub trait Dynamics<T: Num, const X: usize, const U: usize>{
    fn dynamics(&self, x:  &Vec<T, X>, u: [T; U]) ->  Vec<T, X>;
}

pub trait XState<T: Num, const X: usize, const U: usize>{
    fn get_x(&self) -> &Vec<T, X>;
    fn set_x(&mut self, x: Vec<T, X>);
    fn get_theta_idx(&self) -> &ThetaIdx;
    fn get_w_nom(&self) -> T;
}

pub trait RK2Step<T: Num, const X: usize, const U: usize>{
    fn step(&mut self, dt: T, u: [T; U]) -> ();
}

pub struct ThetaIdx {
    pub has_theta: bool,
    pub theta_idx: usize,
}

impl<D, const X: usize, const U: usize> RK2Step<f32, X, U> for D where D: XState<f32, X, U> + Dynamics<f32, X, U>{
    // Steps the dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    fn step(&mut self, dt: f32, u: [f32; U]) {
        let x = self.get_x();
        let dx_dt1 = self.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
        let dx_dt2 = self.dynamics(&(x + dx_dt1 * dt), u);  
        let dx_dt = (dx_dt1 + dx_dt2) * 0.5; // TODO: Use a smaller fractional fixed-point number for dx_dt values?
        let mut x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
        let theta_idx = self.get_theta_idx();
        if  theta_idx.has_theta {
            x1[(theta_idx.theta_idx)] = x1[(theta_idx.theta_idx)] % (2. * PI / self.get_w_nom());
        }
        self.set_x(x1);
    }
}

/* Implement a dynamic step function with no inputs */
pub trait NoInputStep<T, const X: usize, const U: usize>{
    fn step_(&mut self, dt: T) -> ();
}

impl<D, const X: usize, const U: usize> NoInputStep<f32, X, U> for D where D: RK2Step<f32, X, U> + Dynamics<f32, X, U>{
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
type LineStates<T> =  Vec<T, LINE_STATES>;
pub struct RLFilter<T: Num> {
    // RL Filter Parameters
    pub w_nom: T, // nominal frequency (rad)
    rf: T,  // Line resistance (p.u.)
    lf: T,  // Line inductance (p.u.)

    // Internal States
    pub i_alpha: T,  // alpha current state (p.u.)
    pub i_beta: T, // beta current state (p.u.)
    pub x: LineStates<T>,
    theta_idx: ThetaIdx, // -1 for RLFilter
}

impl Dynamics<f32, LINE_STATES, LINE_INPUTS> for RLFilter<f32> {
    // Calculates the p.u. current dynamics for the RL line using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1, theta1, v2, theta2)
    fn dynamics(&self, x: &LineStates<f32>, u: [f32; LINE_INPUTS]) -> LineStates<f32> {
        let (v1, theta1, v2, theta2) = (u[0], u[1], u[2], u[3]);
        let v1_ab = AlphaBeta::from_polar(v1, theta1);
        let v2_ab = AlphaBeta::from_polar(v2, theta2);
        let i_ab = AlphaBeta::from_ab_(x[0], x[1]);
        
        // Unitc line dynamics
        let di_alpha_dt = (v1_ab.alpha - v2_ab.alpha - self.rf * i_ab.alpha) / self.lf;
        let di_beta_dt  = (v1_ab.beta  - v2_ab.beta  - self.rf * i_ab.beta) / self.lf;

        na::Vector2::new(di_alpha_dt, di_beta_dt)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, LINE_STATES, LINE_INPUTS> for RLFilter<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, LINE_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, LINE_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

pub fn build_rl_line(w_nom: f32, rf: f32, lf: f32) -> RLFilter<f32> {
    RLFilter {
        // RL Filter Parameters
        w_nom,
        rf,
        lf,

        // Internal States
        i_alpha: 0.,
        i_beta: 0.,
        x: na::Vector2::new(0., 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 1},
    }
}

/* Define a stiff AC voltage source object */
const ACVS_STATES: usize = 2;
const ACVS_INPUTS: usize = 0;
type ACVSStates<T> =  Vec<T, ACVS_STATES>;
pub struct ACVoltSrc<T: Num> {
    // Parameters
    pub v_nom: T, // nominal RMS LN voltage (V)
    pub w_nom: T, // nominal frequency (rad)
    
    // Internal States
    pub v: T,  // alpha current state (p.u.)
    pub theta: T, // beta current state (p.u.)
    pub x: ACVSStates<T>,
    theta_idx: ThetaIdx,
}

impl Dynamics<f32, ACVS_STATES, ACVS_INPUTS> for ACVoltSrc<f32> {
    // Calculates the p.u. voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of T values: (v, theta)
    // * 'u' - An empty array as there are no inputs
    fn dynamics(&self, _x: &ACVSStates<f32>, _u: [f32; ACVS_INPUTS]) -> ACVSStates<f32> {
        return na::Vector2::new(0., 1.)
    }
}

// Implement functions for getting and setting the states of an ACVoltSrc object
impl<T: Num> XState<T, ACVS_STATES, ACVS_INPUTS> for ACVoltSrc<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, ACVS_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, ACVS_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

pub fn build_ac_volt_src(v_nom: f32, w_nom: f32) -> ACVoltSrc<f32> {
    ACVoltSrc {
        // Parameters
        v_nom,
        w_nom,
        
        // Internal States
        v: 1.,
        theta: 0.,
        x: na::Vector2::new(1., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
    }
}
