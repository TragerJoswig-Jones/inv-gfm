use crate::dynamics::*;
use super::reference_frames::*;
use super::*;

const OSG_INPUTS: usize = 1;
pub trait OrthogonalSignalGenerator<T: Num, const X: usize>: StepDynamics<f32, X, OSG_INPUTS> {
    // /// Function description
    /// Returns the orthogonal signals, [x_alpha, x_beta], of the associated orthogonal signal generator
    fn get_signals(&self) -> [T; 2];
}

/* 
Orthogonal Signal Generator with Second-order Generalized Integrator (SOGI)
*/
// Based on "A New Single-Phase PLL Structure Based on Second Order Generalized Integrator" by Ciobotaru M., Et. al
const SOGI_STATES: usize = 2;
const SOGI_INPUTS: usize = 1;
type SogiStates<T> =  Vec<T, SOGI_STATES>;
pub struct OrthSigGenSogi<T: Num> {
    // parameters
    pub w_nom: T, // nominal frequency
    pub w_res: T, // resonant frequency  // TODO: Create function to set this dynamically based on PLL omega
    pub k: T, // bandwidth scalar

    // states
    pub x: SogiStates<T>, // array of states; [v, theta] //// TODO: Are these right? 
    step_method: fn(&mut dyn StepDynamics<T, SOGI_STATES, SOGI_INPUTS>, T, [T; SOGI_INPUTS])-> Vec<T, SOGI_STATES>,
}

impl OrthSigGenSogi<f32> {
    pub fn new(w_nom: f32, k: f32, step_method: fn(&mut dyn StepDynamics<f32, SOGI_STATES, SOGI_INPUTS>, f32, [f32; SOGI_INPUTS])-> Vec<f32, SOGI_STATES>) -> Self {  // TODO: Determine if we want to implement a constructor like this for each struct?
        OrthSigGenSogi {
            w_nom,
            w_res: w_nom,
            k,
            step_method,
            
            x: na::Vector2::new(0., 0.),
        }
    }
}

impl OrthogonalSignalGenerator<f32, 2> for OrthSigGenSogi<f32> {
    fn get_signals(&self) -> [f32; 2] {
        [self.x[(0)], self.x[(1)]]
    }
}


impl Dynamics<f32, SOGI_STATES, SOGI_INPUTS> for OrthSigGenSogi<f32> {
    fn dynamics(&self, x:  &SogiStates<f32>, u: [f32; SOGI_INPUTS]) -> SogiStates<f32> {
        let v_err = u[0] - x[(0)];
        let dv_dt = self.w_res * ( self.k * v_err - x[(1)] );
        let dqv_dt = self.w_res * x[(0)];
        na::Vector2::new(dv_dt, dqv_dt)
    } 
}

impl StepDynamics<f32, SOGI_STATES, SOGI_INPUTS> for OrthSigGenSogi<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; SOGI_INPUTS]) -> Vec<f32, SOGI_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        return dx_dt
     }
}

// Implement functions for getting and setting the states of the OrthogonalSysGenSogi object
impl<T: Num> XState<T, SOGI_STATES, SOGI_INPUTS> for OrthSigGenSogi<T> {
    fn get_x(&self) -> &Vec<T, SOGI_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, SOGI_STATES>) {
        self.x = x;
    }
}