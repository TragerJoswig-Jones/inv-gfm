/* DC-side controllers */
use crate::dynamics::*;
use crate::*;

/* 
PI Controller
*/
const PI_STATES: usize = 1;
const PI_INPUTS: usize = 1;
const PI_OUTPUTS: usize = 1;
type PiControllerStates<T> =  Vec<T, PI_STATES>;
pub struct PiController<T: Num> {
    // parameters
    pub kp: T, // proptional scalar
    pub ki: T, // integral scalar

    u_ref: T, // reference value

    // states
    pub x: PiControllerStates<T>, // array of states; [integral state]
    pub y: T, // output value
}

impl PiController<f32> {
    pub fn new(kp: f32, ki: f32, ref_val: f32) -> Self {
        PiController {
            kp,
            ki,
            u_ref: ref_val,
            
            x: na::Vector1::new(0.),
            y: 0.,
        }
    }

    pub fn output(&self, u: [f32; PI_INPUTS])  -> [f32; PI_OUTPUTS] {
        let u_err = u[0] - self.u_ref;
        [self.kp * u_err + self.ki * self.x[(0)]]
    }

    pub fn output_step(&mut self, dt: f32, u: [f32; PI_INPUTS])  -> ([f32; PI_OUTPUTS], [f32; PI_STATES]) {
        let u_err = u[0] - self.u_ref;

        // Update the integrator state
        self.x[(0)] = self.x[(0)] + dt * u_err;

        ([self.kp * u_err + self.ki * self.x[(0)]], [u_err])
    }

    pub fn set_ref(&mut self, ref_val: f32) {
        self.u_ref = ref_val;
    }

    pub fn get_output(&self) -> f32 {
        return self.y
    }
}

impl Dynamics<f32, PI_STATES, PI_INPUTS> for PiController<f32> {  // NOTE: Dynamics are not used for this objects step function, but are required by the StepDynamics trait.
    fn dynamics(&self, x:  &PiControllerStates<f32>, u: [f32; PI_INPUTS]) -> PiControllerStates<f32> {
        let u_err = u[0] - self.u_ref;
        na::Vector1::new(u_err)
    } 
}

impl StepDynamics<f32, PI_STATES, PI_INPUTS> for PiController<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; PI_INPUTS]) -> Vec<f32, PI_STATES> {
        let (output, dx_dt) = self.output_step(dt, u);
        self.y = output[0];
        return na::Vector1::new(dx_dt[0])
     }
}

// Implement functions for getting and setting the states of the OrthogonalSysGenSogi object
impl<T: Num> XState<T, PI_STATES, PI_INPUTS> for PiController<T> {
    fn get_x(&self) -> &Vec<T, PI_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, PI_STATES>) {
        self.x = x;
    }
}

/*
Define a Capacitor object 
*/
const CAP_STATES: usize = 1;
const CAP_INPUTS: usize = 1;
const CAP_OUTPUTS: usize = 1;
type CapStates<T> =  Vec<T, CAP_STATES>;
pub struct Capacitor<T: Num> {
    // Capacitor Parameters
    pub v_nom: T, // nominal frequency (rad)
    pub rc: T,  // leakage resistance (p.u.)
    pub c: T,  // capacitance (p.u.)
    
    // Internal States
    pub x: CapStates<T>,  // [v_alpha, alpha voltage state (p.u.); v_beta, beta voltage state (p.u.)]
    step_method: fn(&mut dyn StepDynamics<T, CAP_STATES, CAP_INPUTS>, T, [T; CAP_INPUTS])-> Vec<T, CAP_STATES>,
}

impl Dynamics<f32, CAP_STATES, CAP_INPUTS> for Capacitor<f32> {
    // Calculates the voltage dynamics for the capacitor using the given input, u.
    // # Arguments    
    // * 'x' - internal states as an array of T values: (i_alpha, i_beta)
    // * 'u' - input current as an array of T values: (v1, theta1, v2, theta2)
    fn dynamics(&self, _x: &CapStates<f32>, u: [f32; CAP_INPUTS]) -> CapStates<f32> {
        let dv_dt = u[0] / self.c;
        na::Vector1::new(dv_dt)
    }
}

impl StepDynamics<f32, CAP_STATES, CAP_INPUTS> for Capacitor<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; CAP_INPUTS]) -> Vec<f32, CAP_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the Capacitor object
impl<T: Num> XState<T, CAP_STATES, CAP_INPUTS> for Capacitor<T> {
    fn get_x(&self) -> &Vec<T, CAP_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, CAP_STATES>) {
        self.x = x;
    }
}

impl Capacitor<f32> {
    pub fn new(v_nom: f32, rc: f32, c: f32) -> Self {
        Capacitor {
            // RC Filter Parameters
            v_nom,
            rc,
            c,
    
            // Internal States
            x: na::Vector1::new(1.),
    
            step_method: forward_euler_step,
        }
    }

    pub fn get_voltage(&self, u: [f32; CAP_INPUTS]) -> [f32; CAP_OUTPUTS]{
        [self.x[(0)] + self.rc * u[0]]
    }
}