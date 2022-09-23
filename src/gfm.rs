use crate::calculations::*;
use crate::constants::*;
use crate::dynamics::*;
use crate::inverter::*;
use crate::pll::*;
use crate::reference_frames::*;
use crate::*;

//// Grid Forming Controllers //// 

/* 
Dispatchable Virtual Oscillator Controller (dVOC)
*/
const DVOC_STATES: usize = 2;
const DVOC_INPUTS: usize = 2;
type DvocStates<T> =  Vec<T, DVOC_STATES>;
/// A dispatchable Virtual Oscillator Controller (dVOC) based on ['A Grid-compatible Virtual Oscillator Controller' by Lu M., Et al](https://doi.org/10.1109/ECCE.2019.8913128)
pub struct DvocController<T: Num> {
    // Internal States
    pub x: DvocStates<T>, // array of states; [v, theta] in p.u.
    theta_idx: StateLimits<T, 1>,  // Theta wraps index is 1

    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    pub w_nom: T, // nominal frequency (rad/s)
    xi: T,  // speed constant
    c: T,  // oscillator capacitance (F)
    pub p_ref: T,  // active power reference (p.u.)
    pub q_ref: T,  // reactive power reference (p.u.)

    n_phase: T, // # of phases for power calculation (e.g. Single-Phase, 1., or Three-Phase, 3.)
    step_method: fn(&mut dyn StepDynamics<T, DVOC_STATES, DVOC_INPUTS>, T, [T; DVOC_INPUTS])-> Vec<T, DVOC_STATES>,
}

impl DvocController<f32> {
    /// Constructs a dispatchable virtual oscillator controller from the given controller parameters
    /// # Arguments
    /// * 'v_nom' - nominal voltage (V)
    /// * 'w_nom' - nominal frequency (rad/s)
    /// * 'xi' - speed gain
    /// * 'c' - virtual oscillator capacitance (C)
    /// * 'n_phase' - number of electrical phases of the inverter
    /// * 'step_method' - the method to be used to step the controller (e.g. forward_euler_step, rk2_step)
    pub fn new(v_nom: f32, w_nom: f32, xi: f32, c: f32, n_phase: f32, 
                                 step_method: fn(&mut dyn StepDynamics<f32, DVOC_STATES, DVOC_INPUTS>, f32, [f32; DVOC_INPUTS])-> Vec<f32, DVOC_STATES>
                                ) -> DvocController<f32> {
        DvocController {
            v_nom,
            w_nom,
            x: na::Vector2::new(1., 0.),
            theta_idx: StateLimits::new_theta_wrap([1], w_nom),
            xi,
            c,
            p_ref: 0.,
            q_ref: 0.,
            n_phase,
            step_method,
        }
    }
}

impl Dynamics<f32, DVOC_STATES, DVOC_INPUTS> for DvocController<f32> {
    /// Calculates the voltage dynamics of the dVOC controller using the given input, u.
    /// # Arguments    
    /// * 'x' - polar voltage (p.u.) as an array: [v, theta]
    /// * 'u' - alpha-beta current (p.u.) as an array: [i_alpha, i_beta]
    fn dynamics(&self, x:  &DvocStates<f32>, u: [f32; DVOC_INPUTS]) -> DvocStates<f32> {
        let (v, theta) = (x[0], x[1] * self.w_nom);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0. };  // TODO: Determine the best way to handle multiplying by a constant
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq, self.n_phase);

        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = 1. / (SQRT_2 * self.c * x[0]);
        let dv_dt = 2. * self.xi * x[0] * (1. - (x[0] * x[0])) - _sqrt2cv * (q - self.q_ref);
        let dtheta_dt = 1. - _sqrt2cv / x[0] / self.w_nom * (p - self.p_ref); 
        return na::Vector2::new(dv_dt, dtheta_dt)
    }
}

impl StepDynamics<f32, DVOC_STATES, DVOC_INPUTS> for DvocController<f32> {
    /// Steps the dynamics
    /// # Arguments
    /// * 'dt' - The time period to step the system in seconds
    /// * 'u' - inputs (p.u.) as an array of T values: (i_alpha, i_beta)
    /// 
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; DVOC_INPUTS]) -> Vec<f32, DVOC_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        wrap_angle(&mut self.x, &self.theta_idx);
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, DVOC_STATES, DVOC_INPUTS> for DvocController<T> {
    fn get_x(&self) -> &Vec<T, DVOC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, DVOC_STATES>) {
        self.x = x;
    }
}

impl InvInterface<f32, DVOC_STATES> for DvocController<f32> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of DvocController must be implemented on f32 to use scalars and constants)   
    fn get_voltage(&self) -> [f32; 2] {
        return [self.x[(0)] * self.v_nom, self.x[(1)] * self.w_nom]
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        return [self.x[(0)], self.x[(1)]] 
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.x[(0)] = v[0];  // Sets the voltage magnitude
        self.x[(1)] = v[1];  // Sets the voltage angle
    }
    fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
    fn set_q_ref(&mut self, q_ref: f32) {
        self.q_ref = q_ref;
    }
    fn output(&self) -> [f32; 2] {
        self.get_voltage()
    }
    fn get_w_nom(&self) -> f32 {
        return self.w_nom
    }
}
impl InvController<f32, DVOC_STATES, DVOC_INPUTS> for DvocController<f32> {}

/// Constructs a dispatchable virtual oscillator controller with a default set of controller parameters
/// # Arguments
/// * 'v_nom' - nominal voltage (V)
/// * 'w_nom' - nominal oscillator (rad/s)
pub fn build_default_dvoc_controller(v_nom: f32, f_nom: f32) -> DvocController<f32> {
    DvocController {
        v_nom,
        w_nom: 2.*PI * f_nom,
        x: na::Vector2::new(1., 0.),
        theta_idx: StateLimits::new_theta_wrap([1], 2.*PI * f_nom),
        xi: 15.,
        c: 0.2679,
        p_ref: 0.,
        q_ref: 0.,
        n_phase: 3.,
        step_method: rk2_step,
    }
}

/* 
Droop controller
*/
const DROOP_STATES: usize = 3;
const DROOP_INPUTS: usize = 2;  // 
type DroopStates<T> =  Vec<T, DROOP_STATES>;
/// A droop controller based on ['Control of Parallel Connected Inverters in Standalone ac Supply Systems' by Chandorkar M., Et al.](https://doi.org/10.1109/28.195899)
pub struct DroopController<T: Num> {
    // Internal States
    pub x: DroopStates<T>,  // [theta: angle state (p.u.), p_filt: low-pass filter active power (p.u.), q_filt: low-pass filter reactive power (p.u.)]
    theta_idx: StateLimits<T, 1>,  // Theta wraps index is 0
    
    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    pub w_nom: T, // nominal frequency (rad)
    pub w_c: T, // power low-pass filter cutoff frequency (rad)
    pub mp: T, // frequency droop slope (rad/s)
    pub mq: T, // voltage droop slope (V)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)

    n_phase: T, // # of phases for power calculation (e.g. Single-Phase, 1., or Three-Phase, 3.)
    step_method: fn(&mut dyn StepDynamics<T, DROOP_STATES, DROOP_INPUTS>, T, [T; DROOP_INPUTS])-> Vec<T, DROOP_STATES>,
}

impl DroopController<f32> {
    /// Constructs a droop controller from the given controller parameters
    /// # Arguments
    /// * 'v_nom' - nominal voltage (V)
    /// * 'w_nom' - nominal oscillator (rad/s)
    /// * 'mp' - active power/frequency droop slope
    /// * 'mq' - reactive power/voltage droop slope
    /// * 'n_phase' - number of electrical phases of the inverter
    /// * 'step_method' - the method to be used to step the controller (e.g. forward_euler_step, rk2_step)
    pub fn new(v_nom: f32, w_nom: f32, mp: f32, mq: f32, w_c: f32, n_phase: f32, 
                                  step_method: fn(&mut dyn StepDynamics<f32, DROOP_STATES, DROOP_INPUTS>, f32, [f32; DROOP_INPUTS])-> Vec<f32, DROOP_STATES>
                                 ) -> DroopController<f32> {
        DroopController {
            v_nom,
            w_nom,
            x: na::Vector3::new(0., 0., 0.),
            theta_idx: StateLimits::new_theta_wrap([0], w_nom),
            mp,
            mq,
            w_c,
            p_ref: 0.,
            q_ref: 0.,
            n_phase,
            step_method,
        }
    }
}

impl Dynamics<f32, DROOP_STATES, DROOP_INPUTS> for DroopController<f32> {
    // Calculates the voltage dynamics of the droop controller using the given input, u.
    // # Arguments
    // * 'x' - polar voltage (p.u.) and filtered powers as a tuple of f32 values: (v, theta, p_filt, q_filt)
    // * 'u' - alpha-beta current (A) as a tuple of f32 values: (ialpha, ibeta)
    fn dynamics(&self, x: &DroopStates<f32>, u: [f32; DROOP_INPUTS]) -> DroopStates<f32> {
        let (theta, p_filt, q_filt) = (x[0] * self.w_nom, x[1], x[2]);
        let v = 1. - self.mq * (q_filt - self.q_ref);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq, self.n_phase);

        // Per-unit dynamics (based on eq.13 & eq.17 from 'Control of Parallel Connected Inverters in Standalone ac Supply Systems' by Chandorkar M., Et al.)
        let dp_filt_dt = self.w_c * (p - p_filt);
        let dq_filt_dt = self.w_c * (q - q_filt);
        let dtheta_dt = 1. - self.mp * (p_filt - self.p_ref);
        return na::Vector3::new(dtheta_dt, dp_filt_dt, dq_filt_dt)
    }
}

impl StepDynamics<f32, DROOP_STATES, DROOP_INPUTS> for DroopController<f32> {
    /// Steps the dynamics
    /// # Arguments
    /// * 'dt' - The time period to step the system in seconds
    /// * 'u' - inputs (p.u.) as an array of T values: (i_alpha, i_beta)
    /// 
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; DROOP_INPUTS]) -> Vec<f32, DROOP_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        wrap_angle(&mut self.x, &self.theta_idx);
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, DROOP_STATES, DROOP_INPUTS> for DroopController<T> {  // TODO: Implement XState trait with a macro as it is the same for each object
    fn get_x(&self) -> &Vec<T, DROOP_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, DROOP_STATES>) {
        self.x = x;
    }
}

// TODO: Determine if we want to make this implementation a macro as it will be the same for each GFM controller (Note that we cannot implement it on a generic type that includes all GFM controllers unless we want to access the internal parameters through functions which may be slower...)
impl InvInterface<f32, DROOP_STATES> for DroopController<f32> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of DroopController must be implemented on f32 to use scalars and constants)
    fn get_voltage(&self) -> [f32; 2] {
        [self.compute_voltage(), self.x[(0)] * self.w_nom]
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        [self.compute_pu_voltage(), self.x[(0)]] 
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.x[(0)] = v[1];  // Sets the voltage angle
        // Note the voltage magnitude is not set as it is directly calculated from Q,filt
    }
    // Sets the active power reference within the droop controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in p.u.
    fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
    // Sets the reactive power reference within the droop controller
    // # Arguments
    // * 'q_ref' - The desired reactive power reference in p.u.
    fn set_q_ref(&mut self, q_ref: f32) {
        self.q_ref = q_ref;
    }
    // Returns the reference voltage for the droop controller
    fn output(&self) -> [f32; 2] {
        return [self.compute_voltage(), self.x[(0)]]  // [v, theta]
    }
    fn get_w_nom(&self) -> f32 {
        return self.w_nom
    }
}
impl InvController<f32, DROOP_STATES, DROOP_INPUTS> for DroopController<f32> {}

impl DroopController<f32> {
    /// Returns the p.u. voltage calculated from the Q,filt state
    fn compute_pu_voltage(&self) -> f32 {
        1. - self.mq * (self.x[(2)] - self.q_ref)
    }
    /// Returns the unit voltage calculated from the Q,filt state
    fn compute_voltage(&self) -> f32 {
        self.v_nom * self.compute_pu_voltage()
    }
}

/// Constructs a droop controller from a default set of controller parameters
/// # Arguments
/// * 'v_nom' - nominal voltage (V)
/// * 'f_nom' - nominal frequency (Hz)
pub fn build_default_droop_controller(v_nom: f32, f_nom: f32) -> DroopController<f32> {
    let w_nom = 2. * PI * f_nom;
    DroopController {
        v_nom,
        w_nom,
        x: na::Vector3::new(0., 0., 0.),
        theta_idx: StateLimits::new_theta_wrap([0], w_nom),
        mp: 0.0026,
        mq: 0.005,
        w_c: 2.*PI*30.,
        p_ref: 0.,
        q_ref: 0.,
        n_phase: 3.,
        step_method: rk2_step,
    }
}

/* 
Virtual Synchronous Machine (VSM) controller
*/
const VSM_STATES: usize = 3;
const VSM_INPUTS: usize = 4;  // [ia, ib, vg_alpha , vg_beta]
type VsmStates<T> =  Vec<T, VSM_STATES>;
/// A Virtual Synchronous Machine controller based on ['A Virtual Synchronous Machine implementation for distributed control of power converters in SmartGrid' by D'Arco S., Suul J.A., and Fosso O.B.](https://doi.org/10.1016/j.epsr.2015.01.001)
pub struct VsmController<'a, T: Num, const P: usize> {
    // Internal States
    pub x: VsmStates<T>,  // [theta: angle state (p.u.), omega: angular frequency (p.u.), q_filt: low-pass filter reactive power (p.u.)]
    theta_idx: StateLimits<T, 1>, // Theta, [0], wraps within -pi to pi
    pub pll: &'a mut dyn PhaseLockLoop<T, P>,  // Phase-lock loop object used to track the grid voltage angle
    
    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    pub w_nom: T, // nominal frequency (rad)
    pub w_c: T, // power low-pass filter cutoff frequency (rad)
    pub mp: T, // frequency droop slope (rad/s)
    pub mq: T, // voltage droop slope (V)
    pub j: T, // inertia constant
    pub d: T, // damping coefficient
    pub p_ref: T,  // active power reference (p.u.)
    pub q_ref: T,  // reactive power reference (p.u.)

    n_phase: T, // # of phases for power calculation (e.g. Single-Phase, 1., or Three-Phase, 3.)
    step_method: fn(&mut dyn StepDynamics<T, VSM_STATES, VSM_INPUTS>, T, [T; VSM_INPUTS])-> Vec<T, VSM_STATES>,
}

impl<'a, const P: usize> Dynamics<f32, VSM_STATES, VSM_INPUTS> for VsmController<'a, f32, P> {
    /// Calculates the dynamics of the virtual synchronous machine controller using the given input, u.
    /// The dynamics are based on equations 2, 4, 5, 8, and 9 from 'A Virtual Synchronous Machine implementation for distributed control of power converters in SmartGrid' by D'Arco S., Suul J.A., and Fosso O.B.
    /// # Arguments
    /// * 'x' - voltage angle, omega, and filtered reactive power as a vector of f32 values: &(theta, omega, q_filt)
    /// * 'u' - alpha-beta current and alpha-beta grid-voltage as an array of f32 values: (i_alpha, i_beta, vg_alpha, vg_beta)
    fn dynamics(&self, x: &VsmStates<f32>, u: [f32; VSM_INPUTS]) -> VsmStates<f32> {
        let (theta, omega, q_filt) = (x[0], x[1], x[2]);
        let v = 1. - self.mq * (q_filt - self.q_ref);

        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta * self.w_nom));
        let (p, q) = calc_dq_power(&v_dq, &i_dq, self.n_phase);

        let dq_filt_dt = self.w_c * (q - q_filt);
        let freq_droop = self.mp * (self.p_ref - p);
        let pll_omega = self.pll.get_pu_omega();
        let domega_dt = 1. / (self.j * self.w_nom) * (1. + freq_droop - x[(1)] - self.d * (x[(1)] - pll_omega));  // TODO: Make self.j = self.j * self.w_nom?
        return na::Vector3::new(omega, domega_dt, dq_filt_dt)
    }
}

impl<'a, const P: usize> StepDynamics<f32, VSM_STATES, VSM_INPUTS> for VsmController<'a, f32, P> {
    /// Steps the dynamics
    /// # Arguments
    /// * 'dt' - The time period to step the system in seconds
    /// * 'u' - inputs (p.u.) as an array of T values: (i_alpha, i_beta, vg_alpha, vg_beta)
    /// 
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; VSM_INPUTS]) -> Vec<f32, VSM_STATES> {
        let u_pll = [u[2], u[3]];  // [vg_alpha, vg_beta, vg_gamme]
        self.pll.step(dt, u_pll);  // Steps the pll and stores the omega value  // TODO: Determine if we should have the user step the pll themselves outside of the vsm step function
        let dx_dt = (self.step_method)(self, dt, u);  // Note that the vsm dynamics method takes the grid voltage inputs here and does not use them; This is to satify the input size constraints of the dynamics trait
        wrap_angle(&mut self.x, &self.theta_idx);
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<'a, T: Num, const P: usize> XState<T, VSM_STATES, VSM_INPUTS> for VsmController<'a, T, P> {  // TODO: Implement XState trait with a macro as it is the same for each object
    fn get_x(&self) -> &Vec<T, VSM_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, VSM_STATES>) {
        self.x = x;
    }
}

impl<'a, const P: usize> InvInterface<f32, VSM_STATES> for VsmController<'a, f32, P> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of VsmController must be implemented on f32 to use scalars and constants)
    fn get_voltage(&self) -> [f32; 2] {
        [self.compute_voltage(), self.x[(0)] * self.w_nom]
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        [self.compute_pu_voltage(), self.x[(0)]] 
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.x[(0)] = v[1];  // Sets the voltage angle
        // Note the voltage magnitude is not set as it is directly calculated from Q,filt
    }
    fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
    fn set_q_ref(&mut self, q_ref: f32) {
        self.q_ref = q_ref;
    }
    fn output(&self) -> [f32; 2] {
        return [self.compute_voltage(), self.x[(0)]]
    }
    fn get_w_nom(&self) -> f32 {
        return self.w_nom
    }
}
impl<'a, const P: usize> InvController<f32, VSM_STATES, VSM_INPUTS> for VsmController<'a, f32, P> {} 

impl<'a, const P: usize> VsmController<'a, f32, P> {
    /// Constructs a virtual synchronous machine controller from the given controller parameters
    /// # Arguments
    /// * 'v_nom' - nominal voltage (V)
    /// * 'w_nom' - nominal frequency (rad/s)
    /// * 'mp' - active power/frequency droop slope
    /// * 'mq' - reactive power/voltage droop slope
    /// * 'pll' - reference to a phase-lock loop object
    /// * 'n_phase' - number of electrical phases of the inverter
    /// * 'step_method' - the method to be used to step the controller (e.g. forward_euler_step, rk2_step)
    pub fn new(v_nom: f32, w_nom: f32, mp: f32, mq: f32, j: f32, d: f32, w_c: f32, pll: &'a mut dyn PhaseLockLoop<f32, P>, n_phase: f32,
                                                    step_method: fn(&mut dyn StepDynamics<f32, VSM_STATES, VSM_INPUTS>, f32, [f32; VSM_INPUTS])-> Vec<f32, VSM_STATES>) -> VsmController<'a, f32, P> {
        VsmController {
            v_nom,
            w_nom,
            x: na::Vector3::new(0., 1., 0.),
            theta_idx: StateLimits::new_theta_wrap([0], w_nom),
            pll,
            mp,
            mq,
            j,
            d,
            w_c,
            p_ref: 0.,
            q_ref: 0.,
            n_phase,
            step_method,
        }
    }

    /// Returns the p.u. voltage calculated from the Q,filt state
    fn compute_pu_voltage(&self) -> f32 {
        1. - self.mq * (self.x[(2)] - self.q_ref)
    }
    /// Returns the unit voltage calculated from the Q,filt state
    fn compute_voltage(&self) -> f32 {
        self.v_nom * self.compute_pu_voltage()
    }
}

/// Constructs a virtual synchronous machine controller from a default set of controller parameters
/// # Arguments
/// * 'v_nom' - nominal voltage (V)
/// * 'f_nom' - nominal frequency (Hz)
/// * 'pll' - reference to a phase-lock loop object
pub fn build_default_vsm_controller<'a, const P: usize>(v_nom: f32, f_nom: f32, pll: &'a mut dyn PhaseLockLoop<f32, P>) -> VsmController<'a, f32, P> {
    let w_nom = 2. * PI * f_nom;
    let h = 0.002;
    VsmController {
        v_nom,
        w_nom,
        x: na::Vector3::new(0., 1., 0.),
        theta_idx: StateLimits::new_theta_wrap([0], w_nom),
        pll,
        mp: 0.0026,
        mq: 0.005,
        j: 2. * h / (w_nom*w_nom),
        d: 100. / w_nom,
        w_c: 2.*PI*30.,
        p_ref: 0.,
        q_ref: 0.,
        n_phase: 3.,
        step_method: rk2_step,
    }
}
