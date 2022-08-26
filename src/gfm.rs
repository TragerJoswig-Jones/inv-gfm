use crate::calculations::*;
use crate::constants::*;
use crate::dynamics::*;
use crate::inverter::*;
use crate::reference_frames::*;
use crate::*;

//// Grid Forming Controllers //// 

/* 
Dispatchable Virtual Oscillator Controller (dVOC) implementation 
*/
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
    pub kv: T, // base voltage scalar (V)
    xi: T,
    c: T,  // oscillator capacitance (F)
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
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq);

        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = 1. / (SQRT_2 * self.c * x[0]);
        let dv_dt = 2. * self.xi * x[0] * ((self.x_nom * self.x_nom) - (x[0] * x[0])) - _sqrt2cv * (q - self.q_ref);
        let dtheta_dt = 1. - _sqrt2cv / x[0] / self.w_nom * (p - self.p_ref); 
        return na::Vector2::new(dv_dt, dtheta_dt)
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
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}

// TODO: Determine if we want to make this implementation a macro as it will be the same for each GFM controller (Note that we cannot implement it on a generic type that includes all GFM controllers unless we want to access the internal parameters through functions which may be slower...)
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
    // Sets the active power reference within the dVOC controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in p.u.
    fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
    // Sets the reactive power reference within the dVOC controller
    // # Arguments
    // * 'q_ref' - The desired reactive power reference in p.u.
    fn set_q_ref(&mut self, q_ref: f32) {
        self.q_ref = q_ref;
    }
    // Returns the reference voltage for the dVOC controller
    fn output(&self) -> [f32; 2] {
        self.get_voltage()
    }
}
impl InvController<f32, DVOC_STATES> for DvocController<f32> {}

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

/* 
Droop controller implementation 
*/
const DROOP_STATES: usize = 3;
const DROOP_INPUTS: usize = 2;
type DroopStates<T> =  Vec<T, DROOP_STATES>;
/* Define a droop controller */
pub struct DroopController<T: Num> {
    // Internal States
    pub x: DroopStates<T>,  // [theta: angle state (p.u.), p_filt: low-pass filter active power (p.u.), q_filt: low-pass filter reactive power (p.u.)]
    theta_idx: ThetaIdx,  // Theta index is 0
    
    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    pub w_nom: T, // nominal frequency (rad)
    pub w_c: T, // power low-pass filter cutoff frequency (rad)
    pub mp: T, // frequency droop slope (rad/s)
    pub mq: T, // voltage droop slope (V)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)
}

impl Dynamics<f32, DROOP_STATES, DROOP_INPUTS> for DroopController<f32> {
    // Calculates the voltage dynamics of the droop controller using the given input, u.
    // # Arguments
    // * 'x' - polar voltage (p.u.) and filtered powers as a tuple of f32 values: (v, theta, p_filt, q_filt)
    // * 'u' - alpha-beta current (A) as a tuple of f32 values: (ialpha, ibeta)
    fn dynamics(&self, x: &DroopStates<f32>, u: [f32; DROOP_INPUTS]) -> DroopStates<f32> {
        let (theta, p_filt, q_filt) = (x[0] * self.w_nom, x[1], x[2]);
        let v = 1. - self.mq * (p_filt - self.p_ref);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq);

        // Per-unit dynamics (based on eq.13 & eq.17 from 'Control of Parallel Connected Inverters in Standalone ac Supply Systems' by Chandorkar M., Et al.)
        let dp_filt_dt = self.w_c * (p - p_filt);
        let dq_filt_dt = self.w_c * (q - q_filt);
        let dtheta_dt = 1. - self.mp * (p_filt - self.p_ref);
        return na::Vector3::new(dtheta_dt, dp_filt_dt, dq_filt_dt)
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
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
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
}
impl InvController<f32, DROOP_STATES> for DroopController<f32> {}

impl DroopController<f32> {
    /// Returns the p.u. voltage calculated from the Q,filt state
    fn compute_pu_voltage(&self) -> f32 {
        1. - self.mq * (self.x[(3)] - self.q_ref)
    }
    /// Returns the unit voltage calculated from the Q,filt state
    fn compute_voltage(&self) -> f32 {
        self.v_nom * self.compute_pu_voltage()
    }
}

pub fn build_droop_controller(v_nom: f32, w_nom: f32, mp: f32, mq: f32, w_c: f32) -> DroopController<f32> {
    DroopController {
        v_nom,
        w_nom,
        x: na::Vector3::new(0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 0},
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
        w_nom: 2. * PI * f_nom,
        x: na::Vector3::new(0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 0},
        mp: 0.0026,
        mq: 0.005,
        w_c: 2.*PI*30.,
        p_ref: 0.,
        q_ref: 0.,
    }
}

/* 
Define a double-loop voltage control object 
*/
const DLVC_STATES: usize = 4;
// [vd int, vq int, id int, iq int]
const DLVC_INPUTS: usize = 8;
// [E; ref voltage magnitude (p.u.), omega; ref voltage frequency (rad), 
//  Vd_c; cap direct-axis voltage; , Vq_c; cap quad-axis voltage, 
//  id_f; filter inductor direct-axis current, iq_f; filter inductor quad-axis current,
//  id_g; grid-side direct-axis current, iq_g; grid-side quad-axis current]
const DLVC_OUTPUTS: usize = 2;
// [ud, uq]
type DlvcStates<T> =  Vec<T, DLVC_STATES>;

/// Defines a double-loop voltage controller based on Fig.4b from ('Control of Power Converters in AC Microgrids' by Rocabert J., Et. al)
pub struct DlvController<T: Num> {
    // Parameters
    pub x_nom: T, // nominal unit value
    pub kp_v: T, // voltage proportional gain
    pub ki_v: T, // voltage integral gain
    pub kp_i: T, // current proportional gain
    pub ki_i: T, // current integral gain
    pub lf: T, // filter-side inductance value
    pub cf: T, // filter capacitance value

    pub i_max: T,  // Maximum current to saturate i_ref at
    pub i_min: T,  // Minimum current to saturate i_ref at

    // Internal States
    pub x: DlvcStates<T>,  // [vd int, vq int, id int, iq int]
    theta_idx: ThetaIdx,  // No angular states
}

impl DlvController<f32> {
    // Calculates the dynamics of the double-loop voltage controller's integrators using the given input, u.
    // # Arguments    
    // * 'x' - An array of state values;    [vd int, vq int, id int, iq int].
    // * 'u' - An array of input values;    [E; ref voltage magnitude (p.u.), omega; ref voltage frequency (rad), 
    //                                       Vd_c; cap direct-axis voltage; , Vq_c; cap quad-axis voltage, 
    //                                       id_f; filter inductor direct-axis current, iq_f; filter inductor quad-axis current,
    //                                       id_g; grid-side direct-axis current, iq_g; grid-side quad-axis current].
    pub fn output(&self, u: [f32; DLVC_INPUTS]) -> [f32; DLVC_OUTPUTS] {
        let x = &self.x;
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = u[0] - u[2] - v_virtual_impedance[0];
        let vq_err = - u[3] - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - u[4];
        let iq_err = iq_ref - u[5];

        let vd_ref = self.kp_i * id_err + self.ki_i * x[(2)];
        let vq_ref = self.kp_i * iq_err + self.ki_i * x[(3)];

        let ud = u[2] + vd_ref + u[1] * self.lf  * iq_ref;
        let uq = u[3] + vq_ref - u[1] * self.lf  * id_ref;

        return [ud, uq];
    }

    pub fn step_output(&mut self, dt: f32, u: [f32; DLVC_INPUTS]) -> [f32; DLVC_OUTPUTS] {
        let x = self.get_x();
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = u[0] - u[2] - v_virtual_impedance[0];
        let vq_err = - u[3] - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - u[4];
        let iq_err = iq_ref - u[5];

        let vd_ref = self.kp_i * id_err + self.ki_i * x[(2)];
        let vq_ref = self.kp_i * iq_err + self.ki_i * x[(3)];

        let ud = u[2] + vd_ref + u[1] * self.lf  * iq_ref;
        let uq = u[3] + vq_ref - u[1] * self.lf  * id_ref;

        let dx_dt = na::Vector4::new(vd_err, vq_err, id_err, iq_err);
        self.set_x(x + dx_dt * dt);

        return [ud, uq];
    }
}

// TODO: DETERMINE IF THESE INTEGRATOR DYNAMICS FUNCTIONS SHOULD REMAIN OR IF WE SHOULD SWITCH TO ANOTHER FORMAT FOR THE CONTROLLER TO AVOID DOUBLE CALCULATIONS???
impl Dynamics<f32, DLVC_STATES, DLVC_INPUTS> for DlvController<f32> {  // TODO: DETERMINE IF THIS SHOULD BE REPRESENTED IN A DIFFERENT WAY. THE DYNAMICS FUNCTION MAY NOT BE ABLE TO PROPERLY STORE THE OUTPUT Udq THIS WAY
    // Calculates the dynamics of the double-loop voltage controller's integrators using the given input, u.
    // # Arguments    
    // * 'x' - An array of state values;    [vd int, vq int, id int, iq int].
    // * 'u' - An array of input values;    [E; ref voltage magnitude (p.u.), omega; ref voltage frequency (rad), 
    //                                       Vd_c; cap direct-axis voltage; , Vq_c; cap quad-axis voltage, 
    //                                       id_f; filter inductor direct-axis current, iq_f; filter inductor quad-axis current,
    //                                       id_g; grid-side direct-axis current, iq_g; grid-side quad-axis current].
    fn dynamics(&self, x: &DlvcStates<f32>, u: [f32; DLVC_INPUTS]) -> DlvcStates<f32> {
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = u[0] - u[2] - v_virtual_impedance[0];
        let vq_err = - u[3] - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - u[4]; 
        let iq_err = iq_ref - u[5];
        return na::Vector4::new(vd_err, vq_err, id_err, iq_err)  // TODO: Should the integrator track x_err or ki * x_err?
    }
}

// Implement functions for getting and setting the states of an DLVController object
impl<T: Num> XState<T, DLVC_STATES, DLVC_INPUTS> for DlvController<T> {
    fn get_x(&self) -> &Vec<T, DLVC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, DLVC_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {  // TODO: Possibly change this to a limiter so it can be used for PI integrator saturation as well as angle wrapping???
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.x_nom;
    }
}

pub fn build_double_loop_voltage_controller(x_nom: f32, kp_v: f32, ki_v: f32, kp_i: f32, ki_i: f32, lf: f32, cf: f32, i_max: f32, i_min: f32) -> DlvController<f32> {
    DlvController {
        // Parameters
        x_nom,
        kp_v,
        ki_v,
        kp_i,
        ki_i,
        lf,
        cf,

        i_max,
        i_min,
        
        // Internal States
        x: na::Vector4::new(0., 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 0},
    }
}
