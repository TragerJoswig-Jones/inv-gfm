/* Traits for objects with dynamics */
use crate::constants::*;
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
}

pub trait StepDynamics<T: Num, const X: usize, const U: usize>: XState<T, X, U> + Dynamics<T, X, U> {
    fn step(&mut self, dt: T, u: [T; U]) -> Vec<T, X>;
}

// TODO: Determine the cost of the dynamic dispatch here for StepDynamics...
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

pub trait RK2Step<T: Num, const X: usize, const U: usize>{
    fn step(&mut self, dt: T, u: [T; U]) -> Vec<f32, X>;
}

// impl<D, const X: usize, const U: usize> RK2Step<f32, X, U> for D where D: XState<f32, X, U> + Dynamics<f32, X, U>{
//     // Steps the dynamics using a 2nd-order Runge-Kutta method
//     // # Arguments
//     // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
//     // # Returns the dynamics, 'dx_dt' used to step the states
//     fn step(&mut self, dt: f32, u: [f32; U]) -> Vec<f32, X> {
//         let x = self.get_x();
//         let dx_dt1 = self.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
//         let dx_dt2 = self.dynamics(&(x + dx_dt1 * dt), u);  
//         let dx_dt = (dx_dt1 + dx_dt2) * 0.5; // TODO: Use a smaller fractional fixed-point number for dx_dt values?
//         let mut x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
//         let theta_idx = self.get_theta_idx();
//         if  theta_idx.has_theta {
//             x1[(theta_idx.theta_idx)] = x1[(theta_idx.theta_idx)] % (2. * PI / self.get_w_nom());  // TODO: Consider wrapping from PI to -PI instead of from 2PI to 0. This may reduce the largest number encountered here for fixed-point nums
//         }
//         self.set_x(x1);
//         return dx_dt  // TODO: Determine if this should be returned. So far this is only used for the double-loop voltage controller as it requires the dtheta_dt value from a gfm controller
//     }
// }

/* Implement a dynamic step function with no inputs */
pub trait NoInputStep<T, const X: usize, const U: usize>{
    fn step_(&mut self, dt: T) -> ();
}

impl<D, const X: usize, const U: usize> NoInputStep<f32, X, U> for D where D: StepDynamics<f32, X, U> + Dynamics<f32, X, U>{
    /// Steps the dVOC dynamics
    /// # Arguments
    /// * 'dt' - The step period in seconds (1 / fs)
    /// * 'u' - An empty array as there are no inputs
    fn step_(&mut self, dt: f32) {
        self.step(dt, [0.; U]);
    }
}

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