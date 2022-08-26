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
pub trait InvController<T: Num, const X: usize>: RK2Step<f32, X, 2> + // TODO: Base/Require the InvController to have implemented the NodeInterface trait
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
    fn output(&self) -> [f32; 2];
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
}

pub fn build_gfm<'a, const X: usize>(ctrl: &'a mut dyn InvController<f32, X>, gamma: f32) -> PresyncInvController<'a, f32, X> {
    PresyncInvController { ctrl, synced: false, gamma, sync_tol: 1e-3 }
}