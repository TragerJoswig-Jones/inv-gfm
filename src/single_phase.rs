use crate::dynamics::*;
use super::reference_frames::*;
use super::*;

/* 
Orthogonal System Generator with Second-order Generalized Integrator (SOGI)
*/
// Based on "A New Single-Phase PLL Structure Based on Second Order Generalized Integrator" by Ciobotaru M., Et. al
const SOGI_STATES: usize = 2;
const SOGI_INPUTS: usize = 1;
type SogiStates<T> =  Vec<T, SOGI_STATES>;
pub struct OrthogonalSysGenSogi<T: Num> {
    // parameters
    pub w_res: T, // resonant frequency
    pub k: T, // bandwidth scalar

    // states
    pub x: SogiStates<T>, // array of states; [v, theta]
}

impl Dynamics<f32, SOGI_STATES, SOGI_INPUTS> for OrthogonalSysGenSogi<f32> {
    fn dynamics(&self, x:  &SogiStates<f32>, u: [f32; SOGI_INPUTS]) -> SogiStates<f32> {
        let v_err = u[0] - x[(0)];
        let dv_dt = self.w_res * ( self.k * v_err - x[(1)] );
        let dqv_dt = self.w_res * x[(0)];
    } 
}

// Implement functions for getting and setting the states of the OrthogonalSysGenSogi object
impl<T: Num> XState<T, SOGI_STATES, SOGI_INPUTS> for OrthogonalSysGenSogi<T> {
    fn get_x(&self) -> &Vec<T, GFL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, GFL_STATES>) {
        self.x = x;
    }
    fn get_theta_idx(&self) -> &ThetaIdx {
        return &self.theta_idx
    }
    fn get_w_nom(&self) -> T {
        return self.w_nom
    }
}