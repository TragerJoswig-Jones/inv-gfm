use super::constants::*;
use crate::dynamics::*;
use super::reference_frames::*;
use super::*;

/*
Inverter Controller Interface
*/
/// The 'InvController' trait is used to indicate a control object that can be used within a Presyncrhonization object.
/// -   This trait requires the control object to have already implemented functions to set and get the 
///     state and step the dynamics of the control object with and without input.
/// -   Note this trait is seperate from the InvInterface so that the Inverter can also be implemented 
///     on the Inv struct.
pub trait InvController<T: Num, const X: usize, const U: usize>: StepDynamics<f32, X, U> + // TODO: Base/Require the InvController to have implemented the NodeInterface trait
                                                                 NoInputStep<f32, X, U> + 
                                                                 XState<f32, X, U> + 
                                                                 InvInterface<f32, X> 
{}

pub trait InvInterface<T: Num, const X: usize> {
    /// Returns the voltage magnitude (V) and angle (rad) state values of the Inverter controller
    fn get_voltage(&self) -> [T; 2];

    /// Returns the voltage magnitude (p.u.) and angle (p.u.) state values of the Inverter controller
    fn get_pu_voltage(&self) -> [T; 2];

    /// Sets the voltage magnitude and angle of the controller from the given vector, v
    /// # Arguments
    /// * 'v' - A vector containing the voltage magnitude and angle in p.u.: [v, theta]
    fn set_voltage(&mut self, v: [T; 2]) -> ();  // TODO: Remove this as it is not universal to the controller to have the voltage magnitude as a state

    /// Sets the active power reference within the controller
    /// # Arguments
    /// * 'p_ref' - The desired active power reference in p.u.
    fn set_p_ref(&mut self, p_ref: T) -> ();

    /// Sets the reactive power reference within the controller
    /// # Arguments
    /// * 'q_ref' - The desired reactive power reference in p.u.
    fn set_q_ref(&mut self, q_ref: T) -> ();

    /// Returns the reference voltage for the inverter controller
    fn output(&self) -> [T; 2];

    /// Returns the nominal frequency of the controller
    fn get_w_nom(&self) -> T;
}

/* 
Inverter Controller Presynchronization
*/
/// Defines a 'PresyncInvController' structure that wraps an 'InvController' object and adds presynchronization capabilities.
/// This can be added to any 'InvController' objects, but was designed to be used with grid-forming controllers (e.g. dVOC, droop, VSM)
pub struct PresyncInvController<'a, T: Num, const X: usize, const U: usize> {
    synced: bool,  // If false, then presynchronization dynamics will be used
    gamma: T,  // A scalar used for presynchronization dynamics
    pub sync_tol: T,  // The tolerance for the difference between the grid and inverter voltages to be considered synchronized during presynchronization
    pub ctrl: &'a mut dyn InvController<T, X, U>,  // A reference to a grid forming control object  // TODO: CHECK IF THIS CAN/SHOULD BE IMPLEMENTED WITHOUT THE REFERENCE
}

/// The 'Presync' trait allows a Inverter controller to step its dynamics such that it can synchronize to an AC voltage
pub trait Presync<T: Num, const X: usize, const U: usize> {
    /// Steps the Inverter controller by timestep 'dt' using dynamics according to 'synced' given the inputs 
    /// * 'u' - an array of alpha beta currents in p.u.: (i_alpha, i_beta) 
    /// * 'vg' - an array of alpha beta voltages in p.u.: (vg_alpha, vg_beta).
    fn inv_step(&mut self, dt: T, u: [T; U], vg: [T; 2]) -> Vec<f32, X>;
    /// Steps the Inverter controller by timestep 'dt' using dynamics such that it synchronizes to the input voltage,
    /// * 'vg' - An array of alpha beta voltages in p.u.: (vg_alpha, vg_beta).
    fn presync_step(&mut self, dt: T, u: [T; U], vg: [T; 2]) -> Vec<f32, X>;
    /// Sets the synced parameter of the inverter to true and flips the sign of the voltage vector's magnitude if it is negative
    fn disable_presync(&mut self) -> ();
    /// Sets the synced parameter of the inverter to false
    fn enable_presync(&mut self) -> ();
    /// Checks if the inverter voltage and the given grid voltage, 'vg', are within the 'sync_tol' of the inverter controller voltage
    fn check_sync(&self, vg: [T; 2]) -> bool;
}

impl<'a, const X: usize, const U: usize> Presync<f32, X, U> for PresyncInvController<'a, f32, X, U>{
    fn inv_step(&mut self, dt: f32, u: [f32; U], vg: [f32; 2]) -> Vec<f32, X> {
        if self.synced {
            return self.ctrl.step(dt, u);
        }
        else {
            let mut u_sync = u;  // TODO: Determine a better way to step during presynch for the controllers with PLLs. Here I am zeroing out the measured currents. This is a somewhat ugly approach though...
            u_sync[0] = 0.; u_sync[1] = 0.;
            return self.presync_step(dt, u_sync, vg);
        }
    }
    fn presync_step(&mut self, dt: f32, u: [f32; U], vg: [f32; 2]) -> Vec<f32, X> {
        let v_inv = self.get_pu_voltage();
        let sin_cos = SinCos::<f32>::from_theta(v_inv[1] * self.ctrl.get_w_nom());
        // Presynchronization dynamics based on ('A Pre-synchronization Strategy for Grid-forming Virtual Oscillator Controlled Inverters' by Lu, M., Et al.)
        // dv_sync = -self.gamma * (v_inv[0] - v_grid[0] * cos(v_inv[1] - v_grid[1])) // Where v_inv and v_grid are polar voltages
        // Use trig identities cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y), and cos(theta) = v_alpha / (SQRT_2 * V) to obtain:
        let dv_sync = -self.gamma * (v_inv[0] - ((vg[0] * sin_cos.cos_value() + vg[1] * sin_cos.sin_value()))  / SQRT_2) * dt;  // TODO: Check the conversion to alpha-beta/polar mixed
        let w_sync = self.gamma / v_inv[0] * (vg[1] * sin_cos.cos_value() - vg[0] * sin_cos.sin_value()) * dt;  // TODO: Determine if this needs to be scaled by 1 / w_nom
        self.ctrl.set_voltage([v_inv[0] + dv_sync, v_inv[1] + w_sync]);
        let dx_dt: Vec<f32, X>;
        dx_dt = self.ctrl.step(dt, u);   // Zero current input inverter dynamics; Note that u is required only for PLL based controllers
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
impl<'a, const X: usize, const U: usize> InvInterface<f32, X> for PresyncInvController<'a, f32, X, U> {
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

/// Wraps an inverter controller object in a PresyncInvController object adding presynchronization capabilities to the controller
pub fn add_presynch<'a, const X: usize, const U: usize>(ctrl: &'a mut dyn InvController<f32, X, U>, gamma: f32) -> PresyncInvController<'a, f32, X, U> {
    PresyncInvController { ctrl, synced: false, gamma, sync_tol: 1e-3 }
}

/* 
Double-loop Voltage Controller (DVLC)
*/
const DLVC_STATES: usize = 4;  // [vd int, vq int, id int, iq int]
const DLVC_INPUTS: usize = 8;  
// [E; ref voltage magnitude (p.u.), omega; ref voltage frequency (rad), 
//  Vd_c; cap direct-axis voltage; , Vq_c; cap quad-axis voltage, 
//  id_f; filter inductor direct-axis current, iq_f; filter inductor quad-axis current,
//  id_g; grid-side direct-axis current, iq_g; grid-side quad-axis current]
const DLVC_OUTPUTS: usize = 2;  // [ud, uq]
type DlvcStates<T> =  Vec<T, DLVC_STATES>;

/// Defines a double-loop voltage controller based on Fig.4b from ('Control of Power Converters in AC Microgrids' by Rocabert J., Et. al)
pub struct DlvController<T: Num> {
    // Parameters
    pub kp_v: T, // voltage proportional gain
    pub ki_v: T, // voltage integral gain
    pub kp_i: T, // current proportional gain
    pub ki_i: T, // current integral gain
    pub lf: T, // filter-side inductance value (p.u.)
    pub cf: T, // filter capacitance value (p.u.)

    pub i_max: T,  // Maximum current to saturate i_ref
    pub i_min: T,  // Minimum current to saturate i_ref

    // Internal States
    pub x: DlvcStates<T>,  // [vd int, vq int, id int, iq int]

    step_method: fn(&mut dyn StepDynamics<T, DLVC_STATES, DLVC_INPUTS>, T, [T; DLVC_INPUTS])-> Vec<T, DLVC_STATES>,
    //sat_idx: StateLimits<T, 1>,  //  Saturate current PI controller integral, x[1]
}

impl DlvController<f32> {
    /// Constructs a double-loop voltage controller from the given controller parameters
    /// # Arguments
    /// * 'kp_v' - voltage-loop proportional gain
    /// * 'ki_v' - voltage-loop integral gain
    /// * 'kp_i' - current-loop proportional gain
    /// * 'ki_i' - current-loop integral gain
    /// * 'lf' - filter inductance value (p.u.)
    /// * 'cf' - filter capacitance value (p.u.)
    /// * 'i_max' - maximum reference current value (p.u.)
    /// * 'i_min' - minimum reference current value (p.u.)
    pub fn new(kp_v: f32, ki_v: f32, kp_i: f32, ki_i: f32, lf: f32, cf: f32, i_max: f32, i_min: f32) -> DlvController<f32> {
        DlvController {
            // Parameters
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

    /// Calculates the outputs of the double-loop voltage controller using the given input, u.
    /// # Arguments    
    /// * 'x' - An array of state values: 
    ///     * vd int; voltage direct-axis integrator, 
    ///     * vq int; voltage quadrature-axis integrator, 
    ///     * id int; current direct-axis integrator, , 
    ///     * iq int; current quadrature-axis integrator, .
    /// * 'u' - An array of input values:   
    ///     * E; ref voltage magnitude (p.u.), 
    ///     * omega; ref voltage frequency (rad),
    ///     * vd_c; cap direct-axis voltage; , 
    ///     * vq_c; cap quad-axis voltage,
    ///     * id_f; filter inductor direct-axis current, 
    ///     * iq_f; filter inductor quad-axis current,
    ///     * id_g; grid-side direct-axis current, 
    ///     * iq_g; grid-side quad-axis current.
    pub fn output(&self, u: [f32; DLVC_INPUTS]) -> [f32; DLVC_OUTPUTS] {
        let x = &self.x;
        let (v_ref, omega, vc_d, vc_q, if_d, if_q, ig_d, ig_q) = (u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7]);
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = v_ref - vc_d - v_virtual_impedance[0];
        let vq_err = - vc_q - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - if_d; 
        let iq_err = iq_ref - if_q;

        let vd_ref = self.kp_i * id_err + self.ki_i * x[(2)];
        let vq_ref = self.kp_i * iq_err + self.ki_i * x[(3)];

        let ud = vc_d + vd_ref + omega * self.lf  * iq_ref;
        let uq = vc_q + vq_ref - omega * self.lf  * id_ref;

        return [ud, uq];
    }

    /// Calculates the dynamics of the double-loop voltage controller's integrators and steps their states using the given input, u.
    /// # Arguments    
    /// * 'x' - An array of state values: 
    ///     * vd int; voltage direct-axis integrator, 
    ///     * vq int; voltage quadrature-axis integrator, 
    ///     * id int; current direct-axis integrator, , 
    ///     * iq int; current quadrature-axis integrator, .
    /// * 'u' - An array of input values:   
    ///     * E; ref voltage magnitude (p.u.), 
    ///     * omega; ref voltage frequency (rad),
    ///     * vd_c; cap direct-axis voltage; , 
    ///     * vq_c; cap quad-axis voltage,
    ///     * id_f; filter inductor direct-axis current, 
    ///     * iq_f; filter inductor quad-axis current,
    ///     * id_g; grid-side direct-axis current, 
    ///     * iq_g; grid-side quad-axis current.
    pub fn step_output(&mut self, dt: f32, u: [f32; DLVC_INPUTS]) -> [f32; DLVC_OUTPUTS] {
        let x = self.get_x();
        let (v_ref, omega, vc_d, vc_q, if_d, if_q, ig_d, ig_q) = (u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7]);
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = v_ref - vc_d - v_virtual_impedance[0];
        let vq_err = - vc_q - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - if_d; 
        let iq_err = iq_ref - if_q;

        let vd_ref = self.kp_i * id_err + self.ki_i * x[(2)];
        let vq_ref = self.kp_i * iq_err + self.ki_i * x[(3)];

        let ud = vc_d + vd_ref + omega * self.lf  * iq_ref;
        let uq = vc_q + vq_ref - omega * self.lf  * id_ref;

        let dx_dt = na::Vector4::new(vd_err, vq_err, id_err, iq_err);
        self.set_x(x + dx_dt * dt);

        return [ud, uq];
    }
}

impl Dynamics<f32, DLVC_STATES, DLVC_INPUTS> for DlvController<f32> { 
    /// Calculates the dynamics of the double-loop voltage controller's integrators using the given input, u.
    /// # Arguments    
    /// * 'x' - An array of state values: 
    ///     * vd int; voltage direct-axis integrator, 
    ///     * vq int; voltage quadrature-axis integrator, 
    ///     * id int; current direct-axis integrator, , 
    ///     * iq int; current quadrature-axis integrator, .
    /// * 'u' - An array of input values:   
    ///     * E; ref voltage magnitude (p.u.), 
    ///     * omega; ref voltage frequency (rad),
    ///     * vd_c; cap direct-axis voltage; , 
    ///     * vq_c; cap quad-axis voltage,
    ///     * id_f; filter inductor direct-axis current, 
    ///     * iq_f; filter inductor quad-axis current,
    ///     * id_g; grid-side direct-axis current, 
    ///     * iq_g; grid-side quad-axis current.
    fn dynamics(&self, x: &DlvcStates<f32>, u: [f32; DLVC_INPUTS]) -> DlvcStates<f32> {
        let (v_ref, _omega, vc_d, vc_q, if_d, if_q, ig_d, ig_q) = (u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7]);
        let v_virtual_impedance = [0., 0.];  // TODO: Determine how to implement the virtual impedance / grid-side compensation and get this value here
        let vd_err = v_ref - vc_d - v_virtual_impedance[0];
        let vq_err = - vc_q - v_virtual_impedance[1];

        let mut id_ref = self.kp_v * vd_err + self.ki_v * x[(0)];  // TODO: Missing FF componenet here using LCL cap and capacitor dq voltage?
        let mut iq_ref = self.kp_v * vq_err + self.ki_v * x[(1)];
        
        id_ref = id_ref.clamp(self.i_min, self.i_max);
        iq_ref = iq_ref.clamp(self.i_min, self.i_max);

        let id_err = id_ref - if_d; 
        let iq_err = iq_ref - if_q;
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
