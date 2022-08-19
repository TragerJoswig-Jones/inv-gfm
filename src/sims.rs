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
    fn step(&mut self, dt: T, u: [T; U]) -> Vec<f32, X>;
}

pub struct ThetaIdx {
    pub has_theta: bool,
    pub theta_idx: usize,
}

impl<D, const X: usize, const U: usize> RK2Step<f32, X, U> for D where D: XState<f32, X, U> + Dynamics<f32, X, U>{
    // Steps the dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; U]) -> Vec<f32, X>{
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
        return dx_dt  // TODO: Determine if this should be returned. So far this is only used for the double-loop voltage controller as it requires the dtheta_dt value from a gfm controller
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

/* 
Voltage Node Interface
*/
/// The 'NodeInterface' trait is used to indicate an object that has a voltage state with dynamics.
pub trait NodeInterface<T: Num, const X: usize>: RK2Step<f32, X, 2> + 
                                                 NoInputStep<f32, X, 2> + 
                                                 XState<f32, X, 2> 
{
    /// Returns the voltage state values of the voltage node object
    fn get_voltage(&self) -> [T; 2];
    fn get_pu_voltage(&self) -> [T; 2];
    fn set_voltage(&mut self, v: [T; 2]) -> ();
    // TODO: Should we add get_w_nom and get_v_nom here as well or should this scaling be build into the functions (possibly add get_voltage_pu and set_voltage_pu) 
}

/* 
Current Edge Interface
*/
/// The 'EdgeInterface' trait is used to indicate an object that has a current state with dynamics.
pub trait EdgeInterface<T: Num, const X: usize>: RK2Step<f32, X, 2> + 
                                                 NoInputStep<f32, X, 2> + 
                                                 XState<f32, X, 2> 
{
    /// Returns the current state values of the current edge object
    fn get_current(&self) -> [T; 2];
    fn get_pu_current(&self) -> [T; 2];
    fn set_current(&mut self, v: [T; 2]) -> ();
}


/* 
Define an RL Branch object 
*/
const RL_STATES: usize = 2;
const RL_INPUTS: usize = 4;
type RlStates<T> =  Vec<T, RL_STATES>;
pub struct RlBranch<T: Num> {
    // RL Filter Parameters
    pub w_nom: T, // nominal frequency (rad)
    rf: T,  // Line resistance (p.u.)
    lf: T,  // Line inductance (p.u.)

    // Internal States
    pub i_alpha: T,  // alpha current state (p.u.)
    pub i_beta: T, // beta current state (p.u.)
    pub x: RlStates<T>,
    theta_idx: ThetaIdx, // -1 for RlBranch
}

impl Dynamics<f32, RL_STATES, RL_INPUTS> for RlBranch<f32> {
    // Calculates the p.u. current dynamics for the RL branch using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1, theta1, v2, theta2)
    fn dynamics(&self, x: &RlStates<f32>, u: [f32; RL_INPUTS]) -> RlStates<f32> {
        let (v1_alpha, v1_beta, v2_alpha, v2_beta) = (u[0], u[1], u[2], u[3]);
        let v1_ab = AlphaBeta::from_ab_(v1_alpha, v1_beta);
        let v2_ab = AlphaBeta::from_ab_(v2_alpha, v2_beta);
        let i_ab = AlphaBeta::from_ab_(x[0], x[1]);
        
        // Unitc line dynamics
        let di_alpha_dt = (v1_ab.alpha - v2_ab.alpha - self.rf * i_ab.alpha) / self.lf;
        let di_beta_dt  = (v1_ab.beta  - v2_ab.beta  - self.rf * i_ab.beta) / self.lf;

        na::Vector2::new(di_alpha_dt, di_beta_dt)
    }
}

// Implement functions for getting and setting the states of the RlBranch object
impl<T: Num> XState<T, RL_STATES, RL_INPUTS> for RlBranch<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, RL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, RL_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

pub fn build_rl_branch(w_nom: f32, rf: f32, lf: f32) -> RlBranch<f32> {
    RlBranch {
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

/* 
Define an RC Branch object 
*/
const RC_STATES: usize = 2;
const RC_INPUTS: usize = 2;
type RcStates<T> =  Vec<T, RC_STATES>;
pub struct RcBranch<T: Num> {
    // Rc Branch Parameters
    pub v_nom: T, // nominal frequency (rad)
    pub rc: T,  // Capacitor resistance (p.u.)
    pub cf: T,  // Capacitor capacitance (p.u.)
    
    // Internal States
    pub x: RcStates<T>,  // [v_alpha, alpha voltage state (p.u.); v_beta, beta voltage state (p.u.)]
    theta_idx: ThetaIdx, // No angular terms for RcBranch, but required for generic RK2Step definition
}

impl Dynamics<f32, RC_STATES, RC_INPUTS> for RcBranch<f32> {
    // Calculates the p.u. current dynamics for the RC branch using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1, theta1, v2, theta2)
    fn dynamics(&self, _x: &RlStates<f32>, u: [f32; RC_INPUTS]) -> RlStates<f32> {
        // Unit line dynamics
        let dv_alpha_dt = u[0] / self.cf;
        let dv_beta_dt  = u[1] / self.cf;

        na::Vector2::new(dv_alpha_dt, dv_beta_dt)
    }
}

// Implement functions for getting and setting the states of the RcBranch object
impl<T: Num> XState<T, RC_STATES, RC_INPUTS> for RcBranch<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, RC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, RC_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.v_nom  // TODO: Replace this function somehow?
    }
}

pub fn build_rc_branch(v_nom: f32, rc: f32, cf: f32) -> RcBranch<f32> {
    RcBranch {
        // RC Filter Parameters
        v_nom,
        rc,
        cf,

        // Internal States
        x: na::Vector2::new(SQRT_2, 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 1},
    }
}

/* 
Define an LCL Filter object 
*/
const LCL_STATES: usize = 6;
const LCL_INPUTS: usize = 4;
type LclStates<T> =  Vec<T, LCL_STATES>;
pub struct LclFilter<T: Num> {
    // LCL Filter Components
    pub from_rl_branch: RlBranch<T>,  // TODO: Determine if there is a better way to handle the states of these components. Currently, they just sit idle as the LclFitler states are stepped.
    pub rc_branch: RcBranch<T>, 
    pub to_rl_branch: RlBranch<T>,

    // Internal States
    pub x: LclStates<T>,
    theta_idx: ThetaIdx, // -1 for LCLFilter
}

impl Dynamics<f32, LCL_STATES, LCL_INPUTS> for LclFilter<f32> {
    // Calculates the p.u. current dynamics for the LCL line using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1_alpha, v1_beta, v2_alpha, v2_beta)
    fn dynamics(&self, x: &LclStates<f32>, u: [f32; LCL_INPUTS]) -> LclStates<f32> {    // TODO: Clean up this function to reduce the number of new vectors being created
        let x1 = na::Vector2::new(x[(0)], x[(1)]);
        let x2 = na::Vector2::new(x[(2)], x[(3)]);
        let x3 = na::Vector2::new(x[(4)], x[(5)]);
        let i_cap_alpha = x[(0)] - x[(4)];  // Find cap net currents by taking the difference of the inductor currents
        let i_cap_beta = x[(1)] - x[(5)];
        let v_node_alpha = x[(2)] + i_cap_alpha * self.rc_branch.rc;  // Find voltage at node by taking the sum of the cap voltage and rc voltage
        let v_node_beta = x[(3)] + i_cap_beta * self.rc_branch.rc;
        let u1 = [u[0], u[1], v_node_alpha, v_node_beta];
        let u2 = [i_cap_alpha, i_cap_beta];
        let u3 = [v_node_alpha, v_node_beta, u[2], u[3]];
        let from_di_dt = self.from_rl_branch.dynamics(&x1, u1);  // TODO: make the dynamics trait take a vector or slice for x, such that we can pass a slice of x here to the rl dynamics function (&x.fixed_rows::<2>(0)). Unsure how to make the S term of Matrix generic though (https://stackoverflow.com/questions/60885237/nalgebra-implementing-a-function-for-a-generic-matrixmn)
        let dv_dt = self.rc_branch.dynamics(&x2, u2);
        let to_di_dt = self.to_rl_branch.dynamics(&x3, u3);

        na::Vector6::new(from_di_dt[(0)], from_di_dt[(1)], dv_dt[(0)], dv_dt[(1)], to_di_dt[(0)], to_di_dt[(1)])
    }
}

// Implement functions for getting and setting the states of the LCLFilter object
impl<T: Num> XState<T, LCL_STATES, LCL_INPUTS> for LclFilter<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, LCL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, LCL_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.from_rl_branch.w_nom  // TODO: Get rid of this function from the XState trait?
    }
}

impl LclFilter<f32> {
    pub fn get_from_current(&self) -> [f32; 2] {  // TODO: Determine if this should be a slice of a vector or a new vector or remain an array
        return [self.x[(0)], self.x[(1)]]  // TODO: Determine if this should return the LCL filter states or make calls to the LCL filter components to get their states
    }
    pub fn get_to_current(&self) -> [f32; 2] {
        return [self.x[(4)], self.x[(5)]]
    }
    pub fn get_voltage(&self) -> [f32; 2] {
        return [self.x[(2)], self.x[(3)]]
    }
}

pub fn build_lcl_filter(w_nom: f32, v_nom: f32, rf: f32, lf: f32, rc: f32, cf: f32, rg: f32, lg: f32) -> LclFilter<f32> {
    let from_rl_branch = build_rl_branch(w_nom, rf, lf);
    let rc_branch = build_rc_branch(v_nom, rc, cf);
    let to_rl_branch = build_rl_branch(w_nom, rg, lg);
    LclFilter {
        // RL Filter Parameters
        from_rl_branch,
        rc_branch, 
        to_rl_branch,

        // Internal States
        x: na::Vector6::new(0., 0., SQRT_2, 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 0},
    }
}

/* 
Define a stiff AC voltage source object 
*/
const ACVS_STATES: usize = 2;
const ACVS_INPUTS: usize = 0;
type AcvsStates<T> =  Vec<T, ACVS_STATES>;
pub struct AcVoltSrc<T: Num> {
    // Parameters
    pub v_nom: T, // nominal RMS LN voltage (V)
    pub w_nom: T, // nominal frequency (rad)
    
    // Internal States
    pub v: T,  // alpha current state (p.u.)
    pub theta: T, // beta current state (p.u.)
    pub x: AcvsStates<T>,
    theta_idx: ThetaIdx,
}

impl Dynamics<f32, ACVS_STATES, ACVS_INPUTS> for AcVoltSrc<f32> {
    // Calculates the p.u. voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of T values: (v, theta)
    // * 'u' - An empty array as there are no inputs
    fn dynamics(&self, _x: &AcvsStates<f32>, _u: [f32; ACVS_INPUTS]) -> AcvsStates<f32> {
        return na::Vector2::new(0., 1.)
    }
}

// Implement functions for getting and setting the states of an AcVoltSrc object
impl<T: Num> XState<T, ACVS_STATES, ACVS_INPUTS> for AcVoltSrc<T> {
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

pub fn build_ac_volt_src(v_nom: f32, w_nom: f32) -> AcVoltSrc<f32> {
    AcVoltSrc {
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

