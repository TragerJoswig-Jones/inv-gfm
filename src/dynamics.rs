/* Traits for objects with dynamics */
use crate::constants::*;
use crate::*;

/// Generic trait for defining the dynamics of an object.
/// * 'X' - Number of states, 
/// * 'U' - Number of inputs
pub trait Dynamics<T: Num, const X: usize, const U: usize> {
    /// Returns a vector containing the state dynamics of the object for the given input, u, and states, x.
    fn dynamics(&self, x:  &Vec<T, X>, u: [T; U]) ->  Vec<T, X>;
}

/// Generic trait for accessing the states of an object.
pub trait XState<T: Num, const X: usize, const U: usize> {
    fn get_x(&self) -> &Vec<T, X>;
    fn set_x(&mut self, x: Vec<T, X>);
}

/// Generic trait for adding a function to step the states of an object.
///
/// This trait requires the XState and Dynamics traits to be implemented.
pub trait StepDynamics<T: Num, const X: usize, const U: usize>: XState<T, X, U> + Dynamics<T, X, U> {
    /// Steps the dynamics
    /// # Arguments
    /// * 'dt' - The time period to step the system in seconds
    /// * 'u' - inputs (p.u.) as an array of T values: (input_1, input_2, ...)
    /// 
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: T, u: [T; U]) -> Vec<T, X>;
}

/// Steps a component's states using the given step size (s), dt, and input, u, using the 1st-order Forward Euler method and returns the calculated dynamics
pub fn forward_euler_step<const X: usize, const U: usize>(component: &mut dyn StepDynamics<f32, X, U>, dt: f32, u: [f32; U]) -> Vec<f32, X> {
    let x = component.get_x();
    let dx_dt = component.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
    let x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
    component.set_x(x1);
    return dx_dt
} 

/// Steps a component's states using the given step size (s), dt, and input, u, using a 2nd-order Runge-Kutta method and returns the calculated dynamics
pub fn rk2_step<const X: usize, const U: usize>(component: &mut dyn StepDynamics<f32, X, U>, dt: f32, u: [f32; U]) -> Vec<f32, X> {
    let x = component.get_x();
    let dx_dt1 = component.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
    let dx_dt2 = component.dynamics(&(x + dx_dt1 * dt), u);  
    let dx_dt = (dx_dt1 + dx_dt2) * 0.5; // TODO: Use a smaller fractional fixed-point number for dx_dt values?
    let x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
    component.set_x(x1);
    return dx_dt
}

/// Steps a component's states using the given step size (s), dt, and input, u, using the 4th-order Runge-Kutta method and returns the calculated dynamics
pub fn rk4_step<const X: usize, const U: usize>(component: &mut dyn StepDynamics<f32, X, U>, dt: f32, u: [f32; U]) -> Vec<f32, X> {
    let x = component.get_x();
    let dx_dt1 = component.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
    let dx_dt2 = component.dynamics(&(x + dx_dt1 * dt * 0.5), u);  
    let dx_dt3 = component.dynamics(&(x + dx_dt2 * dt * 0.5), u);  
    let dx_dt4 = component.dynamics(&(x + dx_dt3 * dt), u);  
    let dx_dt = (dx_dt1 + 2.*dx_dt2 + 2.*dx_dt3 + dx_dt4) * ONE_SIXTH; // TODO: Use a smaller fractional fixed-point number for dx_dt values?
    let x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
    component.set_x(x1);
    return dx_dt
}

/* Implement a dynamic step function with no inputs */
pub trait NoInputStep<T, const X: usize, const U: usize>{
    /// Steps the control object with all zero-valued inputs
    fn step_(&mut self, dt: T) -> ();
}

impl<D, const X: usize, const U: usize> NoInputStep<f32, X, U> for D where D: StepDynamics<f32, X, U> + Dynamics<f32, X, U>{
    /// Steps the dVOC dynamics with all zero-valued inputs
    /// # Arguments
    /// * 'dt' - The step period in seconds
    fn step_(&mut self, dt: f32) {
        self.step(dt, [0.; U]);
    }
}

/* 
State Wrapping and Saturation
*/
// A structure for storing information about states that require wrapping or saturation
pub struct StateLimits<T: Num, const N: usize> {
    pub idxs: [usize; N],
    pub mins: [T; N],
    pub maxs: [T; N],
}

impl<const N: usize> StateLimits<f32, N> {
    pub fn new(idxs: [usize; N], mins: [f32; N], maxs: [f32; N]) -> StateLimits<f32, N> {
        StateLimits { idxs, mins, maxs }
    }
    
    pub fn new_theta_wrap(idxs: [usize; N], w_nom: f32) -> StateLimits<f32, N> {
        StateLimits { idxs, mins: [0.; N], maxs: [2.*PI/w_nom; N] }
    }
}

pub fn wrap_states<const X: usize, const N: usize>(x: &mut Vec<f32, X>, wrap_info: &StateLimits<f32, N>) -> () {
    // Wraps the specified states around to be between their corresponding max and min values in 'wrap_info'
    for (i, indx) in wrap_info.idxs.iter().enumerate() {
        if x[(*indx)] > wrap_info.maxs[i] {
            x[(*indx)] = (x[(*indx)] - wrap_info.mins[i]) % (wrap_info.maxs[i] - wrap_info.mins[i]) + wrap_info.mins[i];
        } else if x[(*indx)] < wrap_info.maxs[i] {
            x[(*indx)] = (x[(*indx)] - wrap_info.mins[i]) % (wrap_info.maxs[i] - wrap_info.mins[i]) + wrap_info.maxs[i];
        }
    }
}

pub fn wrap_angle<const X: usize, const N: usize>(x: &mut Vec<f32, X>, wrap_info: &StateLimits<f32, N>) -> () {
    // Computationally cheap implementation of 'wrap_states' for wrapping an angular states that only increments
    for (i, indx) in wrap_info.idxs.iter().enumerate() {
        if x[(*indx)] > wrap_info.maxs[i] {
            x[(*indx)] = x[(*indx)] % wrap_info.maxs[i];
        }
    }
}

pub fn saturate_states<const X: usize, const N: usize>(x: &mut Vec<f32, X>, sat_info: &StateLimits<f32, N>) -> () {
    // Saturates the specified states at the specified max and min values
    for (i, indx) in sat_info.idxs.iter().enumerate() {
        if x[(*indx)] > sat_info.maxs[i] {
            x[(*indx)] = sat_info.maxs[i];
        } else if x[(*indx)] < sat_info.mins[i] {
            x[(*indx)] = sat_info.mins[i]       
        }
    }
}