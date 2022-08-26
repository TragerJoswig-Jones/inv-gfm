
/* Dynamical objects for simulating power systems */
use crate::constants::*;
use crate::dynamics::*;
use crate::reference_frames::*;
use crate::*;

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
const LINE_INPUTS: usize = 4;
/// The 'Line' trait is used to indicate an object that has a current state with dynamics.
pub trait Line<T: Num, const X: usize>: RK2Step<f32, X, 4> + 
                                        NoInputStep<f32, X, 4> + 
                                        XState<f32, X, 4> +
                                        Dynamics<f32, X, 4>
{
    /// Returns the current state values of the current edge object
    fn get_fr_current(&self) -> [T; 2];
    fn get_to_current(&self) -> [T; 2];
    fn get_fr_pu_current(&self) -> [T; 2];
    fn get_to_pu_current(&self) -> [T; 2];
    fn set_fr_current(&mut self, i: [T; 2]) -> ();
    fn set_to_current(&mut self, i: [T; 2]) -> ();
    // Opens and closes line switch
    fn open_switch(&mut self);
    fn close_switch(&mut self);
    fn switch_is_closed(&self) -> bool;
}

/* 
Define an RL Branch object 
*/
const RL_STATES: usize = 2;
const RL_INPUTS: usize = 4;
type RlStates<T> =  Vec<T, RL_STATES>;
pub struct RlBranch<T: Num> {
    // RL Filter Parameters
    pub i_base: T, // base current (A)
    pub w_nom: T, // nominal frequency (rad)
    rf: T,  // Line resistance (p.u.)
    lf: T,  // Line inductance (p.u.)

    switch_closed: bool, // to-side switch state

    // Internal States
    pub i_alpha: T,  // alpha current state (p.u.)
    pub i_beta: T, // beta current state (p.u.)
    pub x: RlStates<T>,
    theta_idx: ThetaIdx, // -1 for RlBranch
}
impl<T: Num> RlBranch<T> {
    pub fn get_current(&self) -> [T; 2] {
        [self.i_base * self.x[(0)], self.i_base * self.x[(1)]]
    }
    pub fn get_pu_current(&self) -> [T; 2] {
        [self.x[(0)], self.x[(1)]]
    }
    pub fn set_current(&mut self, i: [T; 2]) -> () {
        self.x[(0)] = i[0] / self.i_base;
        self.x[(1)] = i[1] / self.i_base;
    }
    pub fn set_pu_current(&mut self, i: [T; 2]) -> () {
        self.x[(0)] = i[0];
        self.x[(1)] = i[1];
    }
}

impl Dynamics<f32, RL_STATES, RL_INPUTS> for RlBranch<f32> {
    // Calculates the p.u. current dynamics for the RL branch using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1, theta1, v2, theta2)
    fn dynamics(&self, x: &RlStates<f32>, u: [f32; RL_INPUTS]) -> RlStates<f32> {
        let (v1_alpha, v1_beta, v2_alpha, v2_beta) = (u[0], u[1], u[2], u[3]);
        let mut v1_ab = AlphaBeta::from_ab_(v1_alpha, v1_beta);
        if !self.switch_closed {
            v1_ab.alpha = v2_alpha;
            v1_ab.beta = v2_beta;
        }
        let v2_ab = AlphaBeta::from_ab_(v2_alpha, v2_beta);
        let i_ab = AlphaBeta::from_ab_(x[0], x[1]);
        
        // calculate dynamics
        let di_alpha_dt = (v1_ab.alpha - v2_ab.alpha - self.rf * i_ab.alpha) / self.lf;
        let di_beta_dt  = (v1_ab.beta  - v2_ab.beta  - self.rf * i_ab.beta) / self.lf;

        na::Vector2::new(di_alpha_dt, di_beta_dt)
    }
}

// Implement functions for getting and setting the states of the RlBranch object
impl<T: Num> XState<T, RL_STATES, RL_INPUTS> for RlBranch<T> {
    fn get_x(&self) -> &Vec<T, RL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, RL_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

// Implement Line trait for RLBranch
impl Line<f32, RL_STATES> for RlBranch<f32> {
    fn get_fr_current(&self) -> [f32; 2] {
        self.get_current()
    }
    fn get_fr_pu_current(&self) -> [f32; 2] {
        self.get_pu_current()
    }
    fn get_to_current(&self) -> [f32; 2] {
        self.get_current()
    }
    fn get_to_pu_current(&self) -> [f32; 2] {
        self.get_pu_current()
    }
    fn set_fr_current(&mut self, i: [f32; 2]) -> () {
        self.set_current(i)
    }
    fn set_to_current(&mut self, i: [f32; 2]) -> () {
        self.set_current(i)
    }
    fn open_switch(&mut self) {
        self.switch_closed = false;
    }
    fn close_switch(&mut self) {
        self.switch_closed = true;
    }
    fn switch_is_closed(&self) -> bool {
        self.switch_closed
    }
}

pub fn build_rl_branch(i_base: f32, w_nom: f32, rf: f32, lf: f32) -> RlBranch<f32> {
    RlBranch {
        // RL Filter Parameters
        i_base,
        w_nom,
        rf,
        lf,

        switch_closed: false,

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
    // Calculates the voltage dynamics for the RC branch using the given input, u.
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
    fn get_x(&self) -> &Vec<T, RC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, RC_STATES>) {
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
    // Parameters
    v_nom: f32,
    i_base: f32,
    
    // LCL Filter Components
    pub from_rl_branch: RlBranch<T>,  // TODO: Determine if there is a better way to handle the states of these components. Currently, they just sit idle as the LclFitler states are stepped.
    pub rc_branch: RcBranch<T>, 
    pub to_rl_branch: RlBranch<T>,

    switch_closed: bool, // to-side switch state

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
        let mut u3 = [v_node_alpha, v_node_beta, u[2], u[3]];
        if !self.switch_closed {
            u3[2] = v_node_alpha;
            u3[3] = v_node_beta;
        }
        let from_di_dt = self.from_rl_branch.dynamics(&x1, u1);  // TODO: make the dynamics trait take a vector or slice for x, such that we can pass a slice of x here to the rl dynamics function (&x.fixed_rows::<2>(0)). Unsure how to make the S term of Matrix generic though (https://stackoverflow.com/questions/60885237/nalgebra-implementing-a-function-for-a-generic-matrixmn)
        let dv_dt = self.rc_branch.dynamics(&x2, u2);
        let to_di_dt = self.to_rl_branch.dynamics(&x3, u3);

        na::Vector6::new(from_di_dt[(0)], from_di_dt[(1)], dv_dt[(0)], dv_dt[(1)], to_di_dt[(0)], to_di_dt[(1)])
    }
}

// Implement functions for getting and setting the states of the LCLFilter object
impl<T: Num> XState<T, LCL_STATES, LCL_INPUTS> for LclFilter<T> {
    fn get_x(&self) -> &Vec<T, LCL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, LCL_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.from_rl_branch.w_nom  // TODO: Get rid of this function from the XState trait?
    }
}

// Implement Line trait for LclFilter
impl Line<f32, LCL_STATES> for LclFilter<f32> {
    fn get_fr_current(&self) -> [f32; 2] {
        [self.x[(0)] * self.i_base, self.x[(1)] * self.i_base]
    }
    fn get_fr_pu_current(&self) -> [f32; 2] {
        [self.x[(0)], self.x[(1)]]
    }
    fn get_to_current(&self) -> [f32; 2] {
        [self.x[(4)] * self.i_base, self.x[(5)] * self.i_base]
    }
    fn get_to_pu_current(&self) -> [f32; 2] {
        [self.x[(4)], self.x[(5)]]
    }
    fn set_fr_current(&mut self, i: [f32; 2]) -> () {
        self.x[(0)] = i[0];
        self.x[(1)] = i[1];
    }
    fn set_to_current(&mut self, i: [f32; 2]) -> () {
        self.x[(4)] = i[0];
        self.x[(5)] = i[1];
    }
    fn open_switch(&mut self) {
        self.switch_closed = false;
    }
    fn close_switch(&mut self) {
        self.switch_closed = true;
    }
    fn switch_is_closed(&self) -> bool {
        self.switch_closed
    }
}

impl LclFilter<f32> {
    pub fn get_fr_current(&self) -> [f32; 2] {  // TODO: Determine if this should be a slice of a vector or a new vector or remain an array
        return [self.x[(0)], self.x[(1)]]  // TODO: Determine if this should return the LCL filter states or make calls to the LCL filter components to get their states
    }
    pub fn get_to_current(&self) -> [f32; 2] {
        return [self.x[(4)], self.x[(5)]]
    }
    pub fn get_voltage(&self) -> [f32; 2] {
        // Calculate the node voltage as sum of cap and rc voltages
        return [self.x[(2)] + self.rc_branch.rc * (self.x[(0)] - self.x[(4)]), self.x[(3)] + self.rc_branch.rc * (self.x[(1)] - self.x[(5)])]  
    }
}

pub fn build_lcl_filter(w_nom: f32, i_base: f32, v_nom: f32, rf: f32, lf: f32, rc: f32, cf: f32, rg: f32, lg: f32) -> LclFilter<f32> {
    let from_rl_branch = build_rl_branch(i_base, w_nom, rf, lf);
    let rc_branch = build_rc_branch(v_nom, rc, cf);
    let to_rl_branch = build_rl_branch(i_base, w_nom, rg, lg);
    LclFilter {
        // Parameters
        i_base,
        v_nom,

        // Components
        from_rl_branch,
        rc_branch, 
        to_rl_branch,

        switch_closed: false,

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
    fn get_x(&self) -> &Vec<T, ACVS_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, ACVS_STATES>) {
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


/* 
Define a LineToBus object 
*/
const LTB_INPUTS: usize = 2;
type LtbStates<T, const X: usize> =  Vec<T, X>;
pub struct LineToBus<'a, T: Num, const X: usize, const N: usize> {
    // Components
    pub line: &'a mut dyn Line<T, N>,  // TODO: Should the LineToBus directly own these components? Could make access to internal states and functions easier... Look at using Boxed for size error https://docs.rust-embedded.org/book/collections/, https://stackoverflow.com/questions/25818082/vector-of-objects-belonging-to-a-trait
    pub bus: AcVoltSrc<T>, 

    // Internal States
    pub x: LtbStates<T, X>,  // TODO: Determine if there is a better way to handle the states of these components. Currently, they just sit idle as the LineToBus states are stepped. Possible to do, but would need to change Xstate trait or implement step seperately for this struct
}

impl<'a, T: Num, const X: usize, const N: usize> LineToBus<'a, T, X, N> {
    pub fn open_switch(&mut self) {
        self.line.open_switch();
    }
    pub fn close_switch(&mut self) {
        self.line.close_switch();
    }
    pub fn switch_is_closed(&self) -> bool {
        self.line.switch_is_closed()
    }
}

impl<'a, const X: usize, const L: usize> Dynamics<f32, X, LTB_INPUTS> for LineToBus<'a, f32, X, L> {
    // Calculates the p.u. current dynamics for the LineToBus using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an vector of T values: (i_alpha, i_beta)
    // * 'u' - input voltages as an array of T values: (v1_alpha, v1_beta)
    fn dynamics(&self, x: &LtbStates<f32, X>, u: [f32; LTB_INPUTS]) -> LtbStates<f32, X> {    // TODO: Clean up this function to reduce the number of new vectors being created
        let x_bus = x.fixed_slice::<2, 1>(0, 0);
        let x_line = x.fixed_slice::<L, 1>(X-L, 0);
        let x_bus_alpha_beta = AlphaBeta::from_polar(x_bus[(0)], x_bus[(1)] * self.bus.w_nom);
        let u_line = [u[0], u[1], x_bus_alpha_beta.alpha, x_bus_alpha_beta.beta];
        let dx_dt_bus = self.bus.dynamics(&x_bus.into(), []);
        let dx_dt_line = self.line.dynamics(&x_line.into(), u_line);  // TODO: make the dynamics trait take a vector or slice for x, such that we can pass a slice of x here to the rl dynamics function (&x.fixed_rows::<2>(0)). Unsure how to make the S term of Matrix generic though (https://stackoverflow.com/questions/60885237/nalgebra-implementing-a-function-for-a-generic-matrixmn)

        let mut dx_dt: Vec<f32, X> = na::zero();
        let mut dx_dt_bus_slice = dx_dt.fixed_slice_mut::<2, 1>(0, 0); 
        dx_dt_bus_slice.copy_from(&dx_dt_bus);
        let mut dx_dt_line_slice = dx_dt.fixed_slice_mut::<L, 1>(X-L, 0);
        dx_dt_line_slice.copy_from(&dx_dt_line);
        return dx_dt
    }
}

impl<'a, const X: usize, const N: usize> XState<f32, X, LTB_INPUTS> for LineToBus<'a, f32, X, N> {
    fn get_x(&self) -> &Vec<f32, X> {
        //let mut x: Vec<f32, X> = zero();
        //let x_bus = self.bus.get_x();
        //let x_line = self.line.get_x();
        //let mut x_l = x.fixed_slice_mut::<2, 1>(0, 0); 
        //x_l.copy_from(x_bus);
        //let mut x_r = x.fixed_slice_mut::<N, 1>(N-1, 0);
        //x_r.copy_from(x_line);
        //self.x = x;  // TODO: Determine how to do with without storing the state in the LineToBus object here... The point of all these calls to elements is to avoid storing the state twice. If I change this back, make sure to make the &mut self param just &self again
        return &self.x
    }
    fn set_x(&mut self, x: Vec<f32, X>) {
        //let x_l = x.fixed_slice::<2, 1>(0, 0);
        //let x_r = x.fixed_slice::<N, 1>(N-1, 0);
        //self.line.set_x(x_r.into());
        //self.bus.set_x(x_l.into());
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.bus.theta_idx
    }
    fn get_w_nom(&self) -> f32 {
        return self.bus.w_nom
    }
}

pub fn build_line_to_bus<'a, const X: usize, const N: usize>(line: &'a mut dyn Line<f32, N>, bus: AcVoltSrc<f32>) -> LineToBus<'a, f32, X, N> {
    // Initialize x to the states of the line and bus
    let mut x: Vec<f32, X> = na::zero();
    let x_bus = bus.get_x();
    let x_line = line.get_x();
    let mut x_l = x.fixed_slice_mut::<2, 1>(0, 0); 
    x_l.copy_from(x_bus);
    let mut x_r = x.fixed_slice_mut::<N, 1>(N-1, 0);
    x_r.copy_from(x_line);

    // Construct and return the LineToBus struct 
    LineToBus { 
        line, 
        bus, 
        
        x,
    }
}