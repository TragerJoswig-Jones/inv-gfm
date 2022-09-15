use crate::calculations::*;
use crate::constants::*;
use crate::dynamics::*;
use crate::inverter::*;
use crate::reference_frames::*;
use crate::*;

// TODO: Make the input to the PLL in the alpha-beta reference frame so it can be generalized to single-phase or three-phase inverters?
const PLL_INPUTS: usize = 3;
pub trait PhaseLockLoop<T: Num, const X: usize>: StepDynamics<f32, X, PLL_INPUTS> {
    // /// Function description
    // Returns the voltage magnitude, V (p.u.), of the associated PLL
    fn get_voltage_magnitude(&self, u: [T; PLL_INPUTS]) -> T;
    // Returns the voltage angle, theta (rad), of the associated PLL
    fn get_voltage_angle(&self) -> T;
    // Returns the angle frequency, omega (rad/s), of the associated PLL
    fn get_omega(&self) -> T;
}

/* 
Synchronous Reference Frame (SRF) Phase-locked Loop (PLL)
*/
pub const SRF_PLL_STATES: usize = 2;
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
    pub omega: T, // stores the previously calculated omega value

    theta_idx: StateLimits<T, 1>,  // Wrap Theta, x[0]
    step_method: fn(&mut dyn StepDynamics<T, SRF_PLL_STATES, SRF_PLL_INPUTS>, T, [T; SRF_PLL_INPUTS])-> Vec<T, SRF_PLL_STATES>,
}

impl SrfPhaseLockedLoop<f32> {  
    pub fn new(w_nom: f32, kp: f32, ki: f32) -> Self {  // TODO: Determine if we want to implement a constructor like this for each struct?
        SrfPhaseLockedLoop { 
            w_nom, 
            kp, 
            ki, 
            x: na::Vector2::new(0., 0.),
            dx_dt: na::Vector2::new(1., 0.),
            omega: 1.,
            theta_idx: StateLimits::new_theta_wrap([0], w_nom),
            step_method: rk2_step,
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

    pub fn output(&self, u: [f32; SRF_PLL_INPUTS]) -> [f32; SRF_PLL_OUTPUTS] {  // TODO: Make this function standard for all controllers? At least dynamical objects with outputs that are based on the input as well
        let x = self.x;
        let sin_cos = SinCos::from_theta(x[(0)]);
        let v_in_dq = DQZ::from_ab(u[0], u[1], u[2], &sin_cos);
        let omega = 1. + self.kp * v_in_dq.q + self.ki * x[(1)];
        let v_mag = libm::sqrtf(v_in_dq.d*v_in_dq.d + v_in_dq.q*v_in_dq.q) / SQRT_2;
        [v_mag, self.x[(0)], omega]
    }
}

impl Dynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<f32> {
    fn dynamics(&self, x:  &SrfPllStates<f32>, u: [f32; SRF_PLL_INPUTS]) ->  SrfPllStates<f32> {
        let sin_cos = SinCos::from_theta(x[(0)] * self.w_nom);
        let v_in_dq = DQZ::from_ab(u[0], u[1], u[2], &sin_cos);  // TODO: Possibly make, v_in_dq, the input u to reduce calculations for the output?
        let omega = 1. + self.kp * v_in_dq.q + self.ki * x[(1)];
        return na::Vector2::new(omega, v_in_dq.q);
    }
}

impl StepDynamics<f32, SRF_PLL_STATES, SRF_PLL_INPUTS> for SrfPhaseLockedLoop<f32> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
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
    fn get_voltage_angle(&self) -> f32 {
        self.x[(0)] * self.w_nom
    }
    fn get_omega(&self) -> f32 {
        self.omega * self.w_nom
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

    pub fn output(&self, u: [f32; DSRF_PLL_INPUTS]) -> [f32; DSRF_PLL_OUTPUTS] {  // TODO: Remove inputs here?
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



/* 
Grid-Following Controller (GFL)
*/
const GFL_STATES: usize = 2;
const GFL_INPUTS: usize = 4;
const GFL_OUTPUTS: usize = 2;
type GflStates<T> =  Vec<T, GFL_STATES>;
/* Define a GFL controller */
pub struct GflController<'a, T: Num, const P: usize> {
    // Internal States
    pub x: GflStates<T>, // array of states; [vd PI integrator, vq PI integrator]
    sat_idx: StateLimits<T, 2>, // indeces and limit values for integrators
    pub pll: &'a mut dyn PhaseLockLoop<T, P>,
    pub v: T,  // previously output voltage magnitude (p.u.)
    pub theta: T, // previouly output voltage angle (p.u.)

    // Parameters
    pub v_nom: T, // nominal voltage (V)
    pub w_nom: T, // nominal frequency (rad/s)
    kp_d: T, // vd PI controller proportional gain
    ki_d: T, // vd PI controller integral gain
    kp_q: T, // vq PI controller proportional gain
    ki_q: T, // vq PI controller integral gain
    l: T,  // filter inductance (p.u.)
    pub p_ref: T,  // Active power reference (p.u.)
    pub q_ref: T,  // Reactive power reference (p.u.)

    n_phase: T, // # of phases for power calculation (e.g. Single-Phase, 1., or Three-Phase, 3.)
    step_method: fn(&mut dyn StepDynamics<T, GFL_STATES, GFL_INPUTS>, T, [T; GFL_INPUTS])-> Vec<T, GFL_STATES>,
}

impl<'a, const P: usize> GflController<'a, f32, P> {
    pub fn output_step(&mut self, dt: f32, u: [f32; GFL_INPUTS]) -> ([f32; GFL_OUTPUTS], [f32; GFL_STATES]) {
        // TODO: REPEATED CODE FROM DYNAMICS: This function steps the integrators and returns the output, but does not fit within the dynamics framework used for other controllers 
        // Step the PLL
        let u_pll = [u[2], u[3], 0.];  // [vg_alpha, vg_beta, vg_gamme]
        self.pll.step(dt, u_pll);  // Steps the pll and stores the omega value
        // Calculate GFL output
        let theta = self.pll.get_voltage_angle();  // Get the angle (rad) of the PLL
        let sin_cos = SinCos::from_theta(theta);
        let i_dq = DQZ::from_ab_(u[0], u[1], &sin_cos);
        let vg_dq = DQZ::from_ab_(u[2], u[3], &sin_cos);
        
        let id_ref = self.p_ref / self.n_phase / vg_dq.d;
        let iq_ref = self.q_ref / self.n_phase / vg_dq.d;
        let id_err = id_ref - i_dq.d;
        let iq_err = iq_ref - i_dq.q;

		let pi_d = self.kp_d * id_err + self.ki_d * self.x[(0)];
        let pi_q = self.kp_q * iq_err + self.ki_q * self.x[(1)];
        
        let u_d = pi_d + vg_dq.d + self.pll.get_omega() * self.l * iq_ref;
        let u_q = pi_q + vg_dq.q - self.pll.get_omega() * self.l * id_ref;
        let v = Polar::from_dqz(u_d, u_q, 0., &sin_cos);

        // Step the integrator states
        self.x[(0)] = self.x[(0)] + id_err;
        self.x[(1)] = self.x[(1)] + iq_err;
        
        ([v.r, v.theta], [id_err, iq_err])
    }

    pub fn output(&self, u: [f32; GFL_INPUTS]) -> [f32; GFL_OUTPUTS] {
        // TODO: REPEATED CODE FROM DYNAMICS: Determine a good way to not repeat this calculation from stepping the dynamics of the integrators. Note the integrators are linear and fully dependant on sampled inputs so only a first-order method is needed.
            let theta = self.pll.get_voltage_angle();  // Get the angle (rad) of the PLL
            let sin_cos = SinCos::from_theta(theta);
            let i_dq = DQZ::from_ab_(u[0], u[1], &sin_cos);
            let vg_dq = DQZ::from_ab_(u[2], u[3], &sin_cos);  // USED FOR OUTPUT
            
            let id_ref = self.p_ref / self.n_phase / vg_dq.d;  // USED FOR OUTPUT
            let iq_ref = self.q_ref / self.n_phase / vg_dq.d;  // USED FOR OUTPUT
            let id_err = id_ref - i_dq.d;  // USED FOR OUTPUT
            let iq_err = iq_ref - i_dq.q;  // USED FOR OUTPUT
        // END REPEATED CODE

		let pi_d = self.kp_d * id_err + self.ki_d * self.x[(0)];
        let pi_q = self.kp_q * iq_err + self.ki_q * self.x[(1)];
        
        let u_d = pi_d + vg_dq.d + self.pll.get_omega() * self.l * iq_ref;
        let u_q = pi_q + vg_dq.q - self.pll.get_omega() * self.l * id_ref;
        let v = Polar::from_dqz(u_d, u_q, 0., &sin_cos);
        [v.r, v.theta]
    }
}

impl<'a, const P: usize> Dynamics<f32, GFL_STATES, GFL_INPUTS> for GflController<'a, f32, P> {
    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - PI integrator states as an array of T fixed-point values: [PI_i_d, PI_i_q]
    // * 'u' - alpha-beta current (p.u.) and alpha-beta grid voltage (p.u.) as an array of T fixed-point values: [i_alpha, i_beta, vg_alpha, vg_beta]
    fn dynamics(&self, _x:  &GflStates<f32>, u: [f32; GFL_INPUTS]) -> GflStates<f32> {
        let theta = self.pll.get_voltage_angle();  // Get the angle (rad) of the PLL
        let sin_cos = SinCos::from_theta(theta);
        let i_dq = DQZ::from_ab_(u[0], u[1], &sin_cos);
        let vg_dq = DQZ::from_ab_(u[2], u[3], &sin_cos);
        
        let id_ref = self.p_ref / self.n_phase / vg_dq.d;
        let iq_ref = self.q_ref / self.n_phase / vg_dq.d;
        let id_err = id_ref - i_dq.d;
        let iq_err = iq_ref - i_dq.q;

        return na::Vector2::new(id_err, iq_err)
    }
}

impl<'a, const P: usize> StepDynamics<f32, GFL_STATES, GFL_INPUTS> for GflController<'a, f32, P> {
    // Steps the dynamics
    // # Arguments
    // * 'u' - inputs (p.u.) as an array of T values: [input1, input2, ...]
    // # Returns the dynamics, 'dx_dt' used to step the states
    fn step(&mut self, dt: f32, u: [f32; GFL_INPUTS]) -> Vec<f32, GFL_STATES> {
        let u_pll = [u[2], u[3], 0.];  // [vg_alpha, vg_beta, vg_gamme]
        self.pll.step(dt, u_pll);  // Steps the pll and stores the omega value
        let dx_dt = (self.step_method)(self, dt, u);
        saturate_states(&mut self.x, &self.sat_idx);
        let output = self.output(u);  // Update the output parameters of the controller
        self.v = output[0];
        self.theta = output[1];
        return dx_dt 
    }
}

// Implement functions for getting and setting the states of the dVOC object
impl<'a, T: Num, const P: usize> XState<T, GFL_STATES, GFL_INPUTS> for GflController<'a, T, P> {
    fn get_x(&self) -> &Vec<T, GFL_STATES> {
        return &self.x
    }
    fn set_x(&mut self, x: Vec<T, GFL_STATES>) {
        self.x = x;
    }
}

// TODO: Determine if we want to make this implementation a macro as it will be the same for each GFM controller (Note that we cannot implement it on a generic type that includes all GFM controllers unless we want to access the internal parameters through functions which may be slower...)
impl<'a, const P: usize> InvInterface<f32, GFL_STATES> for GflController<'a, f32, P> {  // TODO: Determine if this can remain based on generic num type T (Issue arises as dynamics of GflController must be implemented on f32 to use scalars and constants)
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
impl<'a, const P: usize> InvController<f32, GFL_STATES, GFL_INPUTS> for GflController<'a, f32, P> {}

pub fn build_gfl_controller<'a, const P: usize>(v_nom: f32, w_nom: f32, kp_d: f32, ki_d: f32, kp_q: f32, ki_q: f32, l: f32, pll: &'a mut dyn PhaseLockLoop<f32, P>, n_phase: f32) -> GflController<'a, f32, P> {
    GflController {
        x: na::Vector2::new(0., 0.),
        sat_idx: StateLimits{ idxs: [0, 1], mins: [-1., -1.], maxs: [1., 1.] },
        pll,
        v: 1.,
        theta: 0.,
        v_nom,
        w_nom,
        kp_d,
        ki_d,
        kp_q,
        ki_q,
        l,
        p_ref: 0.,
        q_ref: 0.,
        n_phase,
        step_method: rk2_step,
    }
}
