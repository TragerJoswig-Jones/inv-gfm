use crate::dynamics::*;
use super::*;

const OSG_INPUTS: usize = 1;
pub trait OrthogonalSignalGenerator<T: Num, const X: usize>: StepDynamics<f32, X, OSG_INPUTS> {
    /// Returns the orthogonal signals, [x_alpha, x_beta], of the associated orthogonal signal generator
    fn get_signals(&self) -> [T; 2];
    /// Sets the resonant frequency of the orthogonal signal generator. Used to set the OSG frequency based on a PLL.
    fn set_resonant_frequency(&mut self, omega: T) -> ();
}

/* 
Orthogonal Signal Generator with Second-order Generalized Integrator (SOGI)
*/
const SOGI_STATES: usize = 2;
const SOGI_INPUTS: usize = 1;
type SogiStates<T> =  Vec<T, SOGI_STATES>;
/// SOGI Orthogonal Signal Generator based on ["A New Single-Phase PLL Structure Based on Second Order Generalized Integrator" by Ciobotaru M., Et. al](https://doi.org/10.1109/pesc.2006.1711988)
pub struct OrthSigGenSogi<T: Num> {
    // parameters
    pub w_nom: T, // nominal frequency
    pub w_res: T, // resonant frequency  
    pub k: T, // bandwidth scalar

    // states
    pub x: SogiStates<T>, // array of states; [v, theta] //// TODO: Are these right? 
    step_method: fn(&mut dyn StepDynamics<T, SOGI_STATES, SOGI_INPUTS>, T, [T; SOGI_INPUTS])-> Vec<T, SOGI_STATES>,
}

impl OrthSigGenSogi<f32> {
    /// Constructs a SOGI Orthogonal Signal Generator from the given controller parameters
    /// # Arguments
    /// * 'w_nom' - nominal frequency (rad/s)
    /// * 'k' - bandwidth scalar
    /// * 'step_method' - the method to be used to step the controller (e.g. forward_euler_step, rk2_step)
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
    fn set_resonant_frequency(&mut self, omega: f32) -> () {
        self.w_res = omega;
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
    /// Steps the dynamics of the SOGI-OSG
    /// # Arguments
    /// * 'dt' - time step period (s)
    /// * 'u' - input (p.u.) wrapped in an array: (u)
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; SOGI_INPUTS]) -> Vec<f32, SOGI_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        return dx_dt
     }
}

impl<T: Num> XState<T, SOGI_STATES, SOGI_INPUTS> for OrthSigGenSogi<T> {
    fn get_x(&self) -> &Vec<T, SOGI_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, SOGI_STATES>) {
        self.x = x;
    }
}

/* 
Hilbert Transform Orthogonal Signal Generator
*/
// TODO: Implement this OSG flavor based on ['Single-phase synchronisation with Hilbert transformers: a linear and frequency independent orthogonal system generator' by Føyen S., Et al.](http://doi.org/10.13140/RG.2.2.28306.07361)