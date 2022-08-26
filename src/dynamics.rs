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
    fn step(&mut self, dt: f32, u: [f32; U]) -> Vec<f32, X> {
        let x = self.get_x();
        let dx_dt1 = self.dynamics(x, u);  //TODO: Replace dx_dt with reference to vector within the object?
        let dx_dt2 = self.dynamics(&(x + dx_dt1 * dt), u);  
        let dx_dt = (dx_dt1 + dx_dt2) * 0.5; // TODO: Use a smaller fractional fixed-point number for dx_dt values?
        let mut x1 = x + dx_dt * dt;  //TODO: Test if &mut x would allow for direct modification of elements of the state vector allowing us to avoid creating x1 here
        let theta_idx = self.get_theta_idx();
        if  theta_idx.has_theta {
            x1[(theta_idx.theta_idx)] = x1[(theta_idx.theta_idx)] % (2. * PI / self.get_w_nom());  // TODO: Consider wrapping from PI to -PI instead of from 2PI to 0. This may reduce the largest number encountered here for fixed-point nums
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