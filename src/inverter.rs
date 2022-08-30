use super::constants::*;
use crate::dynamics::*;
use super::reference_frames::*;
use super::*;

/* 
Inverter Controller Interface
*/
/// The 'InvController' trait is used to indicate a control object that can be used within a Inverter object.
/// -   This trait requires the control object to have already implemented functions to set and get the 
///     state and step the dynamics of the control object with and without input.
/// -   Note this trait is seperate from the InvInterface so that the Inverter can also be implemented 
///     on the Inv struct.
pub trait InvController<T: Num, const X: usize>: StepDynamics<f32, X, 2> + // TODO: Base/Require the InvController to have implemented the NodeInterface trait
                                                 NoInputStep<f32, X, 2> + 
                                                 XState<f32, X, 2> + 
                                                 InvInterface<f32, X> 
{}

pub trait InvInterface<T: Num, const X: usize> {
    /// Returns the voltage magnitude (p.u.) and angle (rad) state values of the Inverter controller
    fn get_voltage(&self) -> [T; 2];
    fn get_pu_voltage(&self) -> [T; 2];
    fn set_voltage(&mut self, v: [T; 2]) -> ();
    fn set_p_ref(&mut self, p_ref: T) -> ();
    fn set_q_ref(&mut self, q_ref: T) -> ();
    // TODO: Should we add get_w_nom and get_v_nom here as well or should this scaling be build into the functions (possibly add get_voltage_pu and set_voltage_pu) 
    /// Get the voltage reference from the Inverter controller
    fn output(&self) -> [T; 2];
    fn get_w_nom(&self) -> T;
}

/* 
Inverter Controller Presynchronization
*/
/// Defines a 'PresyncInvController' structure that wraps an 'InvController' object and adds presynchronization capabilities
pub struct PresyncInvController<'a, T: Num, const X: usize> {
    pub ctrl: &'a mut dyn InvController<T, X>,  // A reference to a grid forming control object
    synced: bool,  // If false, then presynchronization dynamics will be used
    gamma: T,  // A scalar used for presynchronization dynamics
    pub sync_tol: T,  // The tolerance for the difference between the grid and inverter voltages to be considered synchronized during presynchronization
}

/// The 'Presync' trait allows a Inverter controller to step its dynamics such that it can synchronize to an AC voltage
pub trait Presync<T: Num, const X: usize> {
    /// Steps the Inverter controller by timestep 'dt' using dynamics according to 'synced' given the inputs 
    /// 'u' (An array of alpha beta currents in p.u.: [i.alpha, i.beta]) and 'vg' (An array of alpha beta voltages in p.u.: [vg.alpha, vg.beta]).
    fn inv_step(&mut self, dt: T, u: [T; 2], vg: [T; 2]) -> Vec<f32, X>;
    /// Steps the Inverter controller by timestep 'dt' using dynamics such that it synchronizes to the input voltage,
    /// 'vg' (An array of alpha beta voltages in p.u.: [vg.alpha, vg.beta]).
    fn presync_step(&mut self, dt: T, vg: [T; 2]) -> Vec<f32, X>;
    /// Sets the synced parameter of the inverter to true
    fn disable_presync(&mut self) -> ();
    /// Sets the synced parameter of the inverter to false
    fn enable_presync(&mut self) -> ();
    /// Checks if the inverter voltage and the given grid voltage, 'vg', are within the 'sync_tol' of the inverter controller voltage
    fn check_sync(&self, vg: [T; 2]) -> bool;
}

impl<'a, const X: usize> Presync<f32, X> for PresyncInvController<'a, f32, X>{
    fn inv_step(&mut self, dt: f32, u: [f32; 2], vg: [f32; 2]) -> Vec<f32, X> {
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

// Implement InvInterface for Inverter object such that users can make calls to functions directly from the Inverter object
// TODO: With this should we make the ctrl parameter private???
// TODO: Crate a seperate trait 'NodeInterface' with voltage calls and make it a requirement of GfmInterface (Implement NodeInterface on ACVoltSrc)
impl<'a, const X: usize> InvInterface<f32, X> for PresyncInvController<'a, f32, X> {
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
    fn output(&self) -> [f32; 2] {
        self.ctrl.output()
    }
    fn get_w_nom(&self) -> f32 {
        self.ctrl.get_w_nom()
    }
}

pub fn add_presynch<'a, const X: usize>(ctrl: &'a mut dyn InvController<f32, X>, gamma: f32) -> PresyncInvController<'a, f32, X> {
    PresyncInvController { ctrl, synced: false, gamma, sync_tol: 1e-3 }
}


/* 
Double-loop Voltage Controller (DVLC)
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

    step_method: fn(&mut dyn StepDynamics<T, DLVC_STATES, DLVC_INPUTS>, T, [T; DLVC_INPUTS])-> Vec<T, DLVC_STATES>,
    //sat_idx: StateLimits<T, 1>,  //  Saturate current PI controller integral, x[1]
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

impl StepDynamics<f32, DLVC_STATES, DLVC_INPUTS> for DlvController<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; DLVC_INPUTS]) -> Vec<f32, DLVC_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        return dx_dt 
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
        step_method: rk2_step,
    }
}
