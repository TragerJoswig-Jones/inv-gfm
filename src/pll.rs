use crate::constants::*;
use crate::dynamics::*;
use crate::reference_frames::*;
use crate::*;

const PLL_INPUTS: usize = 2;
pub trait PhaseLockLoop<T: Num, const X: usize>: StepDynamics<f32, X, PLL_INPUTS> {
    /// Returns the voltage magnitude, V , of the associated PLL
    fn get_voltage_magnitude(&self, u: [T; PLL_INPUTS]) -> T;
    /// Returns the voltage angle, theta (p.u.), of the associated PLL
    fn get_pu_voltage_angle(&self) -> T;
    /// Returns the angle frequency, omega (p.u.), of the associated PLL
    fn get_pu_omega(&self) -> T;
    /// Returns the voltage angle, theta (rad), of the associated PLL
    fn get_voltage_angle(&self) -> T;
    /// Returns the angle frequency, omega (rad/s), of the associated PLL
    fn get_omega(&self) -> T;
    /// Returns the nominal angle frequency, omega_nom (rad/s), of the associated PLL
    fn get_omega_nom(&self) -> T;
}

/* 
Synchronous Reference Frame (SRF) Phase-locked Loop (PLL)
*/
pub const SRF_PLL_STATES: usize = 2;
const SRF_PLL_INPUTS: usize = 2;
const SRF_PLL_OUTPUTS: usize = 3;
type SrfPllStates<T> =  Vec<T, SRF_PLL_STATES>;
/// A Synchronous Refrence Frame Phase-Lock Loop (SRF-PLL) 
pub struct SrfPhaseLockedLoop<T: Num> {
    // parameters
    pub w_nom: T, // nominal frequency
    pub kp: T, // proportional gain
    pub ki: T, // integral gain

    // states
    pub x: SrfPllStates<T>, // array of states; [theta, PI controller integral]
    pub dx_dt: SrfPllStates<T>, // TODO: DETERMINE IF I WANT TO STORE THIS VALUE FOR COMPUTING THE OUTPUT?
    pub omega: T, // stores the previously calculated omega value

    theta_idx: StateLimits<T, 1>,  // Wrap Theta, x[0]
    step_method: fn(&mut dyn StepDynamics<T, SRF_PLL_STATES, SRF_PLL_INPUTS>, T, [T; SRF_PLL_INPUTS])-> Vec<T, SRF_PLL_STATES>,
}

impl SrfPhaseLockedLoop<f32> {  
    /// Constructs a synchronous reference frame phase-lock loop from the given controller parameters
    /// # Arguments
    /// * 'w_nom' - nominal frequency (rad/s)
    /// * 'kp' - speed gain
    /// * 'ki' - virtual oscillator capacitance (C)
    /// * 'step_method' - the method to be used to step the controller (e.g. forward_euler_step, rk2_step)
    pub fn new(w_nom: f32, kp: f32, ki: f32, 
               step_method: fn(&mut dyn StepDynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS>, f32, [f32; SRF_PLL_INPUTS])-> Vec<f32, SRF_PLL_STATES>
              ) -> Self {
        SrfPhaseLockedLoop { 
            w_nom, 
            kp, 
            ki, 
            x: na::Vector2::new(0., 0.),
            dx_dt: na::Vector2::new(1., 0.),
            omega: 1.,
            theta_idx: StateLimits::new_theta_wrap([0], w_nom),
            step_method,
        }
    }

    /// Calculates the output of the PLL using the given input, u.
    /// # Arguments    
    /// * 'u' - alpha-beta reference values (p.u.) as an array: [u_alpha, u_beta]
    pub fn output(&self, u: [f32; SRF_PLL_INPUTS]) -> [f32; SRF_PLL_OUTPUTS] { 
        let x = self.x;
        let sin_cos = SinCos::from_theta(x[(0)]);
        let v_in_dq = DQZ::from_ab(u[0], u[1], u[2], &sin_cos);
        let omega = 1. + self.kp * v_in_dq.q + self.ki * x[(1)];
        let v_mag = libm::sqrtf(v_in_dq.d*v_in_dq.d + v_in_dq.q*v_in_dq.q) / SQRT_2;
        [v_mag, self.x[(0)], omega]
    }
}

impl Dynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<f32> {
    /// Calculates the dynamics of the PLL using the given input, u.
    /// # Arguments    
    /// * 'u' - alpha-beta reference values (p.u.) as an array: [u_alpha, u_beta]
    fn dynamics(&self, x:  &SrfPllStates<f32>, u: [f32; SRF_PLL_INPUTS]) ->  SrfPllStates<f32> {
        let sin_cos = SinCos::from_theta(x[(0)] * self.w_nom);
        let v_in_dq = DQZ::from_ab_(u[0], u[1], &sin_cos);  // TODO: Possibly make, v_in_dq, the input u to reduce calculations for the output?
        let omega = 1. + self.kp * v_in_dq.q + self.ki * x[(1)];
        return na::Vector2::new(omega, v_in_dq.q);
    }
}

impl StepDynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<f32> {
    /// Steps the dynamics of the SRF-PLL
    /// # Arguments
    /// * 'dt' - step time period (s)
    /// * 'u' - inputs (p.u.) as an array of T values: [x_alpha, x_beta]
    /// Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; SRF_PLL_INPUTS]) -> Vec<f32, SRF_PLL_STATES> {
        let dx_dt = (self.step_method)(self, dt, u);
        self.omega = dx_dt[0];
        wrap_angle(&mut self.x, &self.theta_idx);
        return dx_dt
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

impl PhaseLockLoop<f32, SRF_PLL_STATES> for SrfPhaseLockedLoop<f32> {
    fn get_voltage_magnitude(&self, u: [f32; PLL_INPUTS]) -> f32 {
        libm::sqrtf(u[0]*u[0] + u[1]*u[1])  // Caculates voltage magnitude from alpha-beta voltages
    } 
    fn get_pu_voltage_angle(&self) -> f32 {
        self.x[(0)]
    }
    fn get_pu_omega(&self) -> f32 {
        self.omega
    }
    fn get_voltage_angle(&self) -> f32 {
        self.x[(0)] * self.w_nom
    }
    fn get_omega(&self) -> f32 {
        self.omega * self.w_nom
    }
    fn get_omega_nom(&self) -> f32 {
        self.w_nom
    }
}

/* 
Double Synchronous Reference Frame (DSRF) Phase-locked Loop (PLL)
*/
// Based on "Double Synchronous Reference Frame PLL for Power Converters Control" by Rodríguez P., Et. al
pub const DSRF_PLL_STATES: usize = 6;
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

    pub fn output(&self, _u: [f32; DSRF_PLL_INPUTS]) -> [f32; DSRF_PLL_OUTPUTS] {  // TODO: Remove inputs here?
        let x = self.x;
        let omega = self.w_nom + self.kp * x[(2)] + self.ki * x[(5)];
        [x[(1)], x[(0)], omega]  // [v_mag, theta, omega]
    }

    fn decoupling_cell(&self, theta: f32, d_n: f32, q_n: f32, d_m: f32, q_m: f32, d_n_filt: f32, q_n_filt: f32, d_m_filt: f32, q_m_filt: f32) -> [f32; 4] {
        // Based on eq. 11 from "Double Synchronous Reference Frame PLL for Power Converters Control" by Rodríguez P., Et. al
        // n is dq^(+1) reference frame
        let sin_cos = SinCos::from_theta(2. * theta * self.w_nom);
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
        let v_alpha_beta = AlphaBeta::from_ab(u[0], u[1], u[2]);
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
// TODO: Implement a FLL