use crate::dynamics::*;
use crate::inverter::*;
use crate::pll::*;
use crate::reference_frames::*;
use crate::*;

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
}

impl<'a, const P: usize> GflController<'a, f32, P> {
    pub fn output_step(&mut self, dt: f32, u: [f32; GFL_INPUTS]) -> ([f32; GFL_OUTPUTS], [f32; GFL_STATES]) {
        // TODO: REPEATED CODE FROM DYNAMICS: This function steps the integrators and returns the output, but does not fit within the dynamics framework used for other controllers 
        // Step the PLL  // TODO: Determine if the PLL dynamics should be built in here somehow to reduce repeated calculations (vg_dq calculated twice atm)
        let u_pll = [u[2], u[3], 0.];  // [vg_alpha, vg_beta, vg_gamme]
        self.pll.step(dt, u_pll);  // Steps the pll and stores the omega value
        // Calculate GFL output
        let theta = self.pll.get_voltage_angle();  // Get the angle (rad) of the PLL
        let sin_cos = SinCos::from_theta(theta);
        let i_dq = DQZ::from_ab_(u[0], u[1], &sin_cos);
        let vg_dq = DQZ::from_ab_(u[2], u[3], &sin_cos);
        
        let id_ref =  self.p_ref * 2. / self.n_phase / vg_dq.d;
        let iq_ref = -self.q_ref * 2. / self.n_phase / vg_dq.d;
        let id_err = id_ref - i_dq.d;
        let iq_err = iq_ref - i_dq.q;

		let pi_d = self.kp_d * id_err + self.ki_d * self.x[(0)];
        let pi_q = self.kp_q * iq_err + self.ki_q * self.x[(1)];
        
        let u_d = pi_d + vg_dq.d + self.pll.get_pu_omega() * self.l * iq_ref;
        let u_q = pi_q + vg_dq.q - self.pll.get_pu_omega() * self.l * id_ref;
        let v = Polar::from_dqz(u_d, u_q, 0., &sin_cos);

        // Step the integrator states
        self.x[(0)] = self.x[(0)] + id_err * dt;
        self.x[(1)] = self.x[(1)] + iq_err * dt;
        
        ([v.r, v.theta / self.w_nom], [id_err, iq_err])
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
        
        let u_d = pi_d + vg_dq.d + self.pll.get_pu_omega() * self.l * iq_ref;
        let u_q = pi_q + vg_dq.q - self.pll.get_pu_omega() * self.l * id_ref;
        let v = Polar::from_dqz(u_d, u_q, 0., &sin_cos);
        [v.r, v.theta / self.w_nom]
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
        //let u_pll = [u[2], u[3], 0.];  // [vg_alpha, vg_beta, vg_gamme]
        //self.pll.step(dt, u_pll);  // Steps the pll and stores the omega value
        let (output, dx_dt) = self.output_step(dt, u);
        saturate_states(&mut self.x, &self.sat_idx);
        self.v = output[0];
        self.theta = output[1];
        return na::Vector2::new(dx_dt[0], dx_dt[1]) 
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
        return [self.v * self.v_nom, self.theta * self.w_nom]
    }
    fn get_pu_voltage(&self) -> [f32; 2] {
        return [self.v, self.theta] 
    }
    fn set_voltage(&mut self, v: [f32; 2]) -> () {
        self.v = v[0];  // Sets the voltage magnitude
        self.theta = v[1];  // Sets the voltage angle
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
    fn output(&self) -> [f32; 2] {  // TODO: Determine if this conflicts with the gfl specific output function. Possibly remove fn output from InvInterface trait
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
        sat_idx: StateLimits{ idxs: [0, 1], mins: [-100., -100.], maxs: [100., 100.] },
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
    }
}
