use crate::calculations::*;
use crate::constants::*;
use crate::dynamics::*;
use crate::inverter::*;
use crate::reference_frames::*;
use crate::*;

// TODO: Make the input to the PLL in the alpha-beta reference frame so it can be generalized to single-phase or three-phase inverters?

/* 
Synchronous Reference Frame (SRF) Phase-locked Loop (PLL)
*/
const SRF_PLL_STATES: usize = 2;
const SRF_PLL_INPUTS: usize = 3;
const SRF_PLL_OUTPUTS: usize = 3;
type SrfPllStates<T> =  Vec<T, SRF_PLL_STATES>;
pub struct SrfPhaseLockedLoop<T: Num> {
    // parameters
    pub w_nom: T, // nominal frequency
    pub kp: T, // proportional gain
    pub ki: T, // integral gain

    // states
    pub x: SrfPllStates<T>, // array of states; [theta, PI controller integral]
    pub dx_dt: SrfPllStates<T>, // TODO: DETERMINE IF I WANT TO STORE THIS VALUE FOR COMPUTING THE OUTPUT?

    wrap_idx: StateLimits<T, 1>,  // Wrap Theta, x[0]
}

impl SrfPhaseLockedLoop<f32> {  
    pub fn new(w_nom: f32, kp: f32, ki: f32) -> Self {  // TODO: Determine if we want to implement a constructor like this for each struct?
        SrfPhaseLockedLoop { 
            w_nom, 
            kp, 
            ki, 
            x: na::Vector2::new(0., 0.),
            dx_dt: na::Vector2::new(w_nom, 0.),
            wrap_idx: StateLimits::new_theta_wrap([0], w_nom),
        }
    }

    // // Steps the PLL dynamics and stores dx_dt
    // pub fn step_dx_dt(&mut self, dt: f32, u: [f32; SRF_PLL_INPUTS]) -> () {  // TODO: DETERMINE IF I WANT TO STORE THIS VALUE FOR COMPUTING THE OUTPUT?
    //     self.dx_dt = self.step(dt, u);
    // }

    // pub fn output_from_dx_dt(&self) -> [f32; SRF_PLL_OUTPUTS] {
    //     let v_mag = libm::fabsf(vd);
    //     [v_mag, self.x[(0)], self.dx_dt[(0)]]  // TODO: DETERMINE IF I WANT TO STORE THESE VALUES FOR COMPUTING THE OUTPUT? Need vd here for output
    // }

    pub fn output(&self, u: [f32; SRF_PLL_INPUTS]) -> [f32; SRF_PLL_OUTPUTS] {
        let x = self.x;
        let sin_cos = SinCos::from_theta(x[(0)]);
        let v_in_dq = DQZ::from_abc(u[0], u[1], u[2], &sin_cos);
        let omega = self.w_nom + self.kp * v_in_dq.q + self.ki * x[(1)];
        let v_mag = libm::fabsf(v_in_dq.d);
        [v_mag, self.x[(0)], omega]
    }
}

impl Dynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<f32> {
    fn dynamics(&self, x:  &SrfPllStates<f32>, u: [f32; SRF_PLL_INPUTS]) ->  SrfPllStates<f32> {
        let sin_cos = SinCos::from_theta(x[(0)]);
        let v_in_dq = DQZ::from_abc(u[0], u[1], u[2], &sin_cos);  // TODO: Possibly make, v_in_dq, the input u to reduce calculations for the output?
        let omega = self.w_nom + self.kp * v_in_dq.q + self.ki * x[(1)];
        return na::Vector2::new(omega, v_in_dq.q);
    }
}

impl<T: Num> XState<T, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<T> {
    fn get_x(&self) -> &Vec<T, SRF_PLL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, SRF_PLL_STATES>) {
        self.x = x;
    }
}

/* 
Double Synchronous Reference Frame (DSRF) Phase-locked Loop (PLL)
*/
// Based on "Double Synchronous Reference Frame PLL for Power Converters Control" by Rodríguez P., Et. al
const DSRF_PLL_STATES: usize = 6;
const DSRF_PLL_INPUTS: usize = 3;
const DSRF_PLL_OUTPUTS: usize = 3;
type DsrfPllStates<T> =  Vec<T, DSRF_PLL_STATES>;
pub struct DsrfPhaseLockedLoop<T: Num> {
    // parameters
    pub w_nom: T, // nominal frequency
    pub w_filt: T, // low-pass filter frequency
    pub kp: T, // proportional gain
    pub ki: T, // integral gain

    // states
    pub x: DsrfPllStates<T>, // array of states; [theta, d^(+1) filt, q^(+1) filt, d^(-1) filt, q^(-1) filt, PI controller integral]
    pub dx_dt: DsrfPllStates<T>, // TODO: DETERMINE IF I WANT TO STORE THIS VALUE FOR COMPUTING THE OUTPUT?
}

impl DsrfPhaseLockedLoop<f32> {
    // // Steps the PLL dynamics and stores dx_dt
    // pub fn step_dx_dt(&mut self, dt: f32, u: [f32; DSRF_PLL_INPUTS]) -> () {  // TODO: DETERMINE IF I WANT TO STORE THIS VALUE FOR COMPUTING THE OUTPUT?
    //     self.dx_dt = self.step(dt, u);
    // }

    pub fn output(&self, u: [f32; DSRF_PLL_INPUTS]) -> [f32; DSRF_PLL_OUTPUTS] {  // TODO: Remove inputs here?
        let x = self.x;
        let omega = self.w_nom + self.kp * x[(2)] + self.ki * x[(5)];
        [x[(1)], x[(0)], omega]  // [v_mag, theta, omega]
    }

    fn decoupling_cell(&self, theta: f32, d_n: f32, q_n: f32, d_m: f32, q_m: f32, d_n_filt: f32, q_n_filt: f32, d_m_filt: f32, q_m_filt: f32) -> [f32; 4] {
        // Based on eq. 11 from "Double Synchronous Reference Frame PLL for Power Converters Control" by Rodríguez P., Et. al
        // n is dq^(+1) reference frame
        let sin_cos = SinCos::from_theta(2. * theta);
        let d_n_star = d_n - d_m_filt * sin_cos.cos_value() - q_m_filt * sin_cos.sin_value();
        let q_n_star = q_n - q_m_filt * sin_cos.cos_value() + d_m_filt * sin_cos.sin_value(); 
        let d_m_star = d_m - d_n_filt * sin_cos.cos_value() + q_n_filt * sin_cos.sin_value();
        let q_m_star = q_m - q_n_filt * sin_cos.cos_value() - d_n_filt * sin_cos.sin_value(); 
        [d_n_star, q_n_star, d_m_star, q_m_star]
    }
}

impl Dynamics<f32, DSRF_PLL_STATES, DSRF_PLL_INPUTS> for DsrfPhaseLockedLoop<f32> {
    // Based on Fig.5 from "Double Synchronous Reference Frame PLL for Power Converters Control" by Rodríguez P., Et. al
    fn dynamics(&self, x:  &DsrfPllStates<f32>, u: [f32; DSRF_PLL_INPUTS]) ->  DsrfPllStates<f32> {
        let sin_cos = SinCos::from_theta(x[(0)]);
        let v_alpha_beta = AlphaBeta::from_abc(u[0], u[1], u[2]);
        let v_dq_pos = v_alpha_beta.to_dqz(&sin_cos); 
        let v_dq_neg = v_alpha_beta.to_dqz(&sin_cos.flip_theta()); 
        
        let dc_values = self.decoupling_cell(x[(0)], v_dq_pos.d, v_dq_pos.q, v_dq_neg.d, v_dq_neg.q, 
                                                       x[(1)], x[(2)], x[(3)], x[(4)]);
        let omega = self.w_nom + self.kp * dc_values[1] + self.ki * x[(5)];  // TODO: Should the addition of w_nom be here? not shown in fig.5
        return na::Vector6::new(omega, //  anglular velocity
                                self.w_filt * (dc_values[0] - x[(1)]), // v_s_d_pos lpf
                                self.w_filt * (dc_values[1] - x[(2)]), // v_s_q_pos lpf
                                self.w_filt * (dc_values[2] - x[(3)]), // v_s_d_neg lpf
                                self.w_filt * (dc_values[3] - x[(4)]), // v_s_q_neg lpf
                                dc_values[1]); // PI controller integral 
    }
}

impl<T: Num> XState<T, DSRF_PLL_STATES, DSRF_PLL_INPUTS> for DsrfPhaseLockedLoop<T> {
    fn get_x(&self) -> &Vec<T, DSRF_PLL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, DSRF_PLL_STATES>) {
        self.x = x;
    }
}

/* 
Frequency-locked Loop (FLL)
*/



/* 
Grid-Following Controller (GFL)
*/
const GFL_STATES: usize = 2;
const GFL_INPUTS: usize = 2;
type GflStates<T> =  Vec<T, GFL_STATES>;
/* Define a GFL controller */
pub struct GflController<T: Num> {
    // Internal States
    pub v: T,  // voltage state (p.u.)
    pub theta: T, // angle state (p.u.)
    pub x: GflStates<T>, // array of states; [v, theta]
    theta_idx: StateLimits<T, 1>, // index of theta value; 1 for GFL states

    // Other Parameters
    pub v_nom: T, // nominal voltage (V)
    x_nom: T, // nominal voltage (p.u.)
    pub w_nom: T, // nominal frequency (rad/s)
    pub kv: T, // Base voltage (V)
    xi: T,
    c: T,  // Oscillator capacitance (F)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)

    n_phase: T, // # of phases for power calculation (e.g. Single-Phase, 1., or Three-Phase, 3.)
    step_method: fn(&mut dyn StepDynamics<T, GFL_STATES, GFL_INPUTS>, T, [T; GFL_INPUTS])-> Vec<T, GFL_STATES>,
}

impl Dynamics<f32, GFL_STATES, GFL_INPUTS> for GflController<f32> {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as an array of T fixed-point values: [v, theta]
    // * 'u' - alpha-beta current (p.u.) as an array of T fixed-point values: [ialpha, ibeta]
    fn dynamics(&self, x:  &GflStates<f32>, u: [f32; GFL_INPUTS]) -> GflStates<f32> {
        let (v, theta) = (x[0], x[1] * self.w_nom);
        let v_dq = DQZ{ d: v * SQRT_2, q: 0., z: 0. };  // TODO: Determine the best way to handle multiplying by a constant
        let i_dq = AlphaBeta::from_ab_(u[0], u[1]).to_dqz(&SinCos::<f32>::from_theta(theta));
        let (p, q) = calc_dq_power(&v_dq, &i_dq, self.n_phase);

        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = 1. / (SQRT_2 * self.c * x[0]);
        let dv_dt = 2. * self.xi * x[0] * ((self.x_nom * self.x_nom) - (x[0] * x[0])) - _sqrt2cv * (q - self.q_ref);
        let dtheta_dt = 1. - _sqrt2cv / x[0] / self.w_nom * (p - self.p_ref); 
        return na::Vector2::new(dv_dt, dtheta_dt)
    }
}

impl StepDynamics<f32, GFL_STATES, GFL_INPUTS> for GflController<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; GFL_INPUTS]) -> Vec<f32, GFL_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        wrap_angle(&mut self.x, &self.theta_idx);
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<T: Num> XState<T, GFL_STATES, GFL_INPUTS> for GflController<T> {
    fn get_x(&self) -> &Vec<T, GFL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, GFL_STATES>) {
        self.x = x;
    }
}

// TODO: Determine if we want to make this implementation a macro as it will be the same for each GFM controller (Note that we cannot implement it on a generic type that includes all GFM controllers unless we want to access the internal parameters through functions which may be slower...)
impl InvInterface<f32, GFL_STATES> for GflController<f32> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of GflController must be implemented on f32 to use scalars and constants)
    fn get_voltage(&self) -> [f32; 2] {
        return [self.x[(0)] * self.v_nom, self.x[(1)] * self.w_nom]
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        return [self.x[(0)], self.x[(1)]] 
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.x[(0)] = v[0];  // Sets the voltage magnitude
        self.x[(1)] = v[1];  // Sets the voltage angle
    }
    // Sets the active power reference within the dVOC controller
    // # Arguments
    // * 'p_ref' - The desired active power reference in p.u.
    fn set_p_ref(&mut self, p_ref: f32) {
        self.p_ref = p_ref;
    }
    // Sets the reactive power reference within the dVOC controller
    // # Arguments
    // * 'q_ref' - The desired reactive power reference in p.u.
    fn set_q_ref(&mut self, q_ref: f32) {
        self.q_ref = q_ref;
    }
    // Returns the reference voltage for the dVOC controller
    fn output(&self) -> [f32; 2] {
        self.get_voltage()
    }
    fn get_w_nom(&self) -> f32 {
        return self.w_nom
    }
}
impl InvController<f32, GFL_STATES> for GflController<f32> {}

pub fn build_gfl_controller(v_nom: f32, w_nom: f32, xi: f32, c: f32, n_phase: f32) -> GflController<f32> {
    GflController {
        v_nom,
        x_nom: 1.,
        w_nom,
        v: 1.,  
        theta: 0.,
        x: na::Vector2::new(1., 0.),
        theta_idx: StateLimits::new_theta_wrap([0], w_nom),
        kv: v_nom,
        xi,
        c,
        p_ref: 0.,
        q_ref: 0.,
        n_phase,
        step_method: rk2_step,
    }
}