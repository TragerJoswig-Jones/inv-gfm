use super::calc::*;
use super::constants::*;
use super::refs::*;
use super::sims::*;
use super::*;

/* 
Grid-Forming Controller Interface
*/
/// The 'GFMController' trait is used to indicate a control object that can be used within a GFM object.
/// -   This trait requires the control object to have already implemented functions to set and get the 
///     state and step the dynamics of the control object with and without input.
/// -   Note this trait is seperate from the GFMInterface so that the GFMInterface can also be implemented 
///     on the GFM struct.
pub trait GFMController<T: Num, const X: usize>: RK2Step<f32, X, 2> + // TODO: Base/Require the GFMController to have implemented the NodeInterface trait
                                                 NoInputStep<f32, X, 2> + 
                                                 XState<f32, X, 2> + 
                                                 GFMInterface<f32, X> 
{}

pub trait GFMInterface<T: Num, const X: usize> {
    /// Returns the voltage magnitude (p.u.) and angle (rad) state values of the GFM controller
    fn get_voltage(&self) -> [T; 2];
    fn get_pu_voltage(&self) -> [T; 2];
    fn set_voltage(&mut self, v: [T; 2]) -> ();
    fn set_p_ref(&mut self, p_ref: T) -> ();
    fn set_q_ref(&mut self, q_ref: T) -> ();
    // TODO: Should we add get_w_nom and get_v_nom here as well or should this scaling be build into the functions (possibly add get_voltage_pu and set_voltage_pu) 
}

pub struct GFM<'a, T: Num, const X: usize> {
    pub ctrl: &'a mut dyn GFMController<T, X>,  // A reference to a grid forming control object
    synced: bool,  // If false, then presynchronization dynamics will be used
    gamma: T,  // A scalar used for presynchronization dynamics
    pub sync_tol: T,  // The tolerance for the difference between the grid and GFM voltages to be considered synchronized during presynchronization
}

/// The 'Presync' trait allows a GFM controller to step its dynamics such that it can synchronize to the input voltage
pub trait Presync<T: Num, const X: usize> {
    /// Steps the GFM controller by timestep 'dt' using dynamics according to 'synced' given the inputs 
    /// 'u' (An array of alpha beta currents in p.u.: [i.alpha, i.beta]) and 'vg' (An array of alpha beta voltages in p.u.: [vg.alpha, vg.beta]).
    fn gfm_step(&mut self, dt: T, u: [T; 2], vg: [T; 2]) -> Vec<f32, X>;
    /// Steps the GFM controller by timestep 'dt' using dynamics such that it synchronizes to the input voltage,
    /// 'vg' (An array of alpha beta voltages in p.u.: [vg.alpha, vg.beta]).
    fn presync_step(&mut self, dt: T, vg: [T; 2]) -> Vec<f32, X>;
    /// Sets the synced parameter of the GFM to true
    fn disable_presync(&mut self) -> ();
    /// Sets the synced parameter of the GFM to false
    fn enable_presync(&mut self) -> ();
    /// Checks if the GFM voltage and the given grid voltage, 'vg', are within the 'sync_tol' of the GFM
    fn check_sync(&self, vg: [T; 2]) -> bool;
}

impl<'a, const X: usize> Presync<f32, X> for GFM<'a, f32, X>{
    fn gfm_step(&mut self, dt: f32, u: [f32; 2], vg: [f32; 2]) -> Vec<f32, X> {
        if self.synced {
            return self.ctrl.step(dt, u);
        }
        else {
            return self.presync_step(dt, vg);
        }
    }
    fn presync_step(&mut self, dt: f32, vg: [f32; 2]) -> Vec<f32, X> {
        let v_inv = self.get_pu_voltage();
        let sin_cos = SinCos::<f32>::from_theta(v_inv[1] * self.ctrl.get_w_nom());
        // Presynchronization dynamics based on ('A Pre-synchronization Strategy for Grid-forming Virtual Oscillator Controlled Inverters' by Lu, M., Et al.)
        // dv_sync = -self.gamma * (v_inv[0] - v_grid[0] * cos(v_inv[1] - v_grid[1])) // Where v_inv and v_grid are polar voltages
        // Use trig identities cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y), and cos(theta) = v_alpha / (SQRT_2 * V) to obtain:
        let dv_sync = -self.gamma * (v_inv[0] - ((vg[0] * sin_cos.cos_value() + vg[1] * sin_cos.sin_value()))  / SQRT_2) * dt;  // TODO: Check the conversion to alpha-beta/polar mixed
        let w_sync = self.gamma / v_inv[0] * (vg[1] * sin_cos.cos_value() - vg[0] * sin_cos.sin_value()) * dt;  // TODO: Determine if this needs to be scaled by 1 / w_nom
        self.ctrl.set_voltage([v_inv[0] + dv_sync, v_inv[1] + w_sync]);
        let dx_dt = self.ctrl.step(dt, [0.; 2]);   // Zero-input inverter dynamics
        // Combine the default dynamics and presync dynamics to return
        let mut dx_dt_sync: Vec<f32, X> = na::zero();
        dx_dt_sync[(0)] = dv_sync; dx_dt_sync[(1)] = w_sync;
        return dx_dt + dx_dt_sync; 
    }
    fn disable_presync(&mut self) -> () {
        let v = self.get_pu_voltage();
        if v[0] < 0. {  // If the voltage magnitude is negative flip the magnitude to be positive and rotate the increment the angle by 180 deg
           self.set_voltage([-v[0], v[1] + PI / self.ctrl.get_w_nom()]) 
        }
        self.synced = true;
    }
    fn enable_presync(&mut self) -> () {
        self.synced = false;
    }
    fn check_sync(&self, vg: [f32; 2]) -> bool {
        // TODO: Determine the most computational efficient way to check the distance between the two voltages. Note that vg should be an alpha-beta voltage, while v_inv is in polar coords
        let v_inv = self.get_pu_voltage();
        let sin_cos = SinCos::<f32>::from_theta(v_inv[1] * self.ctrl.get_w_nom());
        let dv_sync = -self.gamma * (v_inv[0] - ((vg[0] * sin_cos.cos_value() + vg[1] * sin_cos.sin_value()))  / SQRT_2);  
        let w_sync = self.gamma / v_inv[0] * (vg[1] * sin_cos.cos_value() - vg[0] * sin_cos.sin_value()); 
        if (dv_sync*dv_sync + w_sync*w_sync) < self.sync_tol {
            return true
        }
        return false
    }
}

// Implement GFMInterface for GFM object such that users can make calls to functions directly from the GFM object
// TODO: With this should we make the ctrl parameter private???
// TODO: Crate a seperate trait 'NodeInterface' with voltage calls and make it a requirement of GFMInterface (Implement NodeInterface on ACVoltSrc)
impl<'a, const X: usize> GFMInterface<f32, X> for GFM<'a, f32, X> {
    fn get_voltage(&self) -> [f32; 2] {
        self.ctrl.get_voltage()
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        self.ctrl.get_pu_voltage()
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.ctrl.set_voltage(v)
    }
    fn set_p_ref(&mut self, p_ref: f32) -> () {
        self.ctrl.set_p_ref(p_ref)
    }
    fn set_q_ref(&mut self, q_ref: f32) -> () {
        self.ctrl.set_q_ref(q_ref)
    }
}

pub fn build_gfm<'a, const X: usize>(ctrl: &'a mut dyn GFMController<f32, X>, gamma: f32) -> GFM<'a, f32, X> {
    GFM { ctrl, synced: false, gamma, sync_tol: 1e-3 }
}

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
    pub kv: T, // Base voltage (V)
    xi: T,
    c: T,  // Oscillator capacitance (F)
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
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<f32>::from_theta(theta));
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
    fn get_x(&self) -> &nalgebra::SVector<T, DVOC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, DVOC_STATES>) {
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
impl GFMInterface<f32, DVOC_STATES> for DvocController<f32> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of DvocController must be implemented on f32 to use scalars and constants)
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
}
impl GFMController<f32, DVOC_STATES> for DvocController<f32> {}

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
const DROOP_STATES: usize = 4;
const DROOP_INPUTS: usize = 2;
type DroopStates<T> =  Vec<T, DROOP_STATES>;
/* Define a droop controller */
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
        let (v, theta, p_filt, q_filt) = (x[0], x[1] * self.w_nom, x[2], x[3]);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0.};
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq);

        // Per-unit dynamics (based on eq.13 & eq.17 from 'Control of Parallel Connected Inverters in Standalone ac Supply Systems' by Chandorkar M., Et al.)
        let dp_filt_dt = self.w_c * (p - p_filt);
        let dq_filt_dt = self.w_c * (q - q_filt);
        let dv_dt = - self.mq * dq_filt_dt;
        let dtheta_dt = 1. - self.mp * (p_filt - self.p_ref);
        return na::Vector4::new(dv_dt, dtheta_dt, dp_filt_dt, dq_filt_dt)
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, DROOP_STATES, DROOP_INPUTS> for DroopController<T> {  // TODO: Implement XState trait with a macro as it is the same for each object
    fn get_x(&self) -> &nalgebra::SVector<T, DROOP_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, DROOP_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
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

pub fn build_droop_controller(v_nom: f32, w_nom: f32, mp: f32, mq: f32, w_c: f32) -> DroopController<f32> {
    DroopController {
        v_nom,
        w_nom,
        v: 1.,  
        theta: 0.,
        p_filt: 0.,
        q_filt: 0.,
        x: na::Vector4::new(1., 0., 0., 0.),
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
        w_nom: 2. * PI * f_nom,
        v: 1.,  
        theta: 0.,
        p_filt: 0.,
        q_filt: 0.,
        x: na::Vector4::new(1., 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: true, theta_idx: 1},
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

    // Internal States
    pub x: DlvcStates<T>,  // [vd int, vq int, id int, iq int]
    theta_idx: ThetaIdx,  // No angular states
}

impl DlvController<f32> {
    pub fn output(&self, u: [f32; DLVC_INPUTS]) -> [f32; DLVC_OUTPUTS] {
        let x = &self.x;
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = u[0] - u[2] - v_virtual_impedance[0];
        let vq_err = - u[3] - v_virtual_impedance[1];

        let id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];

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

        let id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];

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
        let vq_err = - u[2] - v_virtual_impedance[1];

        let id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];
        let iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];

        let id_err = id_ref - u[4]; 
        let iq_err = iq_ref - u[5];
        return na::Vector4::new(vd_err, vq_err, id_err, iq_err)
    }
}

// Implement functions for getting and setting the states of an DLVController object
impl<T: Num> XState<T, DLVC_STATES, DLVC_INPUTS> for DlvController<T> {
    fn get_x(&self) -> &nalgebra::SVector<T, DLVC_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: nalgebra::SVector<T, DLVC_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {  // TODO: Possibly change this to a limiter so it can be used for PI integrator saturation as well as angle wrapping???
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.x_nom;
    }
}

pub fn build_double_loop_voltage_controller(x_nom: f32, kp_v: f32, ki_v: f32, kp_i: f32, ki_i: f32, lf: f32, cf: f32) -> DlvController<f32> {
    DlvController {
        // Parameters
        x_nom,
        kp_v,
        ki_v,
        kp_i,
        ki_i,
        lf,
        cf,
        
        // Internal States
        x: na::Vector4::new(0., 0., 0., 0.),
        theta_idx: ThetaIdx {has_theta: false, theta_idx: 0},
    }
}
