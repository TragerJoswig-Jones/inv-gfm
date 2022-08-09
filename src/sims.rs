/* Dynamical objects for simulating power systems */
use super::constants::*;
use super::refs::*;
use super::*;

// TODO: Determine what should be in per unit and how to handle unit conversions

pub type State = f32;  // Type alias to indicate a state variable
pub type Input = f32;  // Type alias to indicate an input variable
pub type Param = f32;  // Type alias to indicate a parameter

/* Implement dynamics for all available elements */
// * 'X' - Number of states, 
// * 'U' - Number of inputs
pub trait Dynamics<const X: usize, const U: usize>{
    fn dynamics(&self, x:  &Vec<State, X>, u: [Input; U]) ->  Vec<State, X>;
}

pub trait XState<const X: usize, const U: usize>{
    fn get_x(&self) -> &Vec<State, X>;
    fn set_x(&mut self, x: Vec<State, X>);
    fn get_theta_idx(&self) -> &ThetaIdx;
}

pub trait RK2Step<const X: usize, const U: usize>{
    fn step(&mut self, dt: Param, u: [Input; U]) -> ();
}

pub struct ThetaIdx {
    pub has_theta: bool,
    pub theta_idx: usize,
}

impl<T, const X: usize, const U: usize> RK2Step<X, U> for T where T: XState<X, U> + Dynamics<X, U>{
    // Steps the dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'u' - alpha-beta current as a tuple of f32 values: (ialpha, ibeta)
    fn step(&mut self, dt: Param, u: [Input; U]) {
        let x = self.get_x();
        let dx_dt1 = self.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
        let dx_dt2 = self.dynamics(&(x + dt * dx_dt1), u);
        let dx_dt = (dx_dt1 + dx_dt2) * 0.5;
        let mut x1 = x + dt * dx_dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
        let theta_idx = self.get_theta_idx();
        if  theta_idx.has_theta {
            x1[(theta_idx.theta_idx)] = x1[(theta_idx.theta_idx)] % (2.*PI);
        }
        self.set_x(x1);
    }
}

/* Implement a dynamic step function with no inputs */
pub trait NoInputStep<const X: usize, const U: usize>{
    fn step_(&mut self, dt: Param) -> ();
}

impl<T, const X: usize, const U: usize> NoInputStep<X, U> for T where T: RK2Step<X, U> + Dynamics<X, U>{
    // Steps the dVOC dynamics
    // # Arguments
    // * 'dt' - The step period in seconds (1 / fs)
    // * 'u' - An empty array as there are no inputs
    fn step_(&mut self, dt: Param) {
        self.step(dt, [0.; U]);
    }
}

/* Define an RL Filter object */
const LINE_STATES: usize = 2;
const LINE_INPUTS: usize = 4;
type LineStates =  Vec<State, LINE_STATES>;
pub struct RLFilter {
    // RL Filter Parameters
    pub w_nom: Param, // nominal frequency (rad)
    pub s_rated: Param,  // maximum expected power output (VA)
    rf: Param,  // Line resistance (Ohms)
    lf: Param,  // Line inductance (H)

    // Internal States
    pub i_alpha: State,  // alpha current state (p.u.)
    pub i_beta: State, // beta current state (p.u.)
    pub x: LineStates,
    theta_idx: ThetaIdx, // -1 for RLFilter
}

impl Dynamics<LINE_STATES, LINE_INPUTS> for RLFilter {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of f32 values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of f32 values: (v1, theta1, v2, theta2)
    fn dynamics(&self, x: &LineStates, u: [Input; LINE_INPUTS]) -> LineStates {
        let (v1, theta1, v2, theta2) = (u[0], u[1], u[2], u[3]);
        let v1_ab = AlphaBeta::from_polar(v1, theta1);
        let v2_ab = AlphaBeta::from_polar(v2, theta2);
        let i_ab = AlphaBeta::from_ab_(x[0], x[1]);
        
        // Unitc line dynamics
        let di_alpha_dt = 1./self.lf*(v1_ab.alpha - v2_ab.alpha - self.rf*i_ab.alpha);
        let di_beta_dt  = 1./self.lf*(v1_ab.beta - v2_ab.beta - self.rf*i_ab.beta);

        na::Vector2::new(di_alpha_dt, di_beta_dt)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl XState<LINE_STATES, LINE_INPUTS> for RLFilter {
    fn get_x(&self) -> &nalgebra::SVector<State, LINE_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<State, LINE_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
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
        x: na::Vector2::new(0., 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 1},
    }
}

/* Define a stiff AC voltage source object */
const ACVS_STATES: usize = 2;
const ACVS_INPUTS: usize = 0;
type ACVSStates=  Vec<State, ACVS_STATES>;
pub struct ACVoltSrc {
    // Parameters
    pub v_nom: Param, // nominal RMS LN voltage (V)
    pub w_nom: Param, // nominal frequency (rad)
    pub s_rated: Param,  // maximum expected power output (VA)
    
    // Internal States
    pub v: State,  // alpha current state (p.u.)
    pub theta: State, // beta current state (p.u.)
    pub x: ACVSStates,
    theta_idx: ThetaIdx,
}

impl Dynamics<ACVS_STATES, ACVS_INPUTS> for ACVoltSrc {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of f32 values: (v, theta)
    // * 'u' - An empty array as there are no inputs
    fn dynamics(&self, _x: &ACVSStates, _u: [Input; ACVS_INPUTS]) -> ACVSStates {
        return na::Vector2::new(0., self.w_nom)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl XState<ACVS_STATES, ACVS_INPUTS> for ACVoltSrc {
    fn get_x(&self) -> &nalgebra::SVector<State, ACVS_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<State, ACVS_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
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
        x: na::Vector2::new(v_nom, 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
    }
}

