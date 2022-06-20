use super::calc::calc_power;
use super::constants::{SQRT_2, SQRT_3};
use super::refs::{alpha_beta_fr_polar, alpha_beta_fr_ab};

#[test]
fn test_dvoc_params() {
    let v_nom: f32 = 120.;
    let w_nom: f32 = 60.;
    let s_rated: f32 = 500.;
    let dt: f32 = 1.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let dvoc = build_dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);
    assert_eq!(dvoc.v_nom, v_nom);
    assert_eq!(dvoc.w_nom, w_nom);
    assert_eq!(dvoc.s_rated, s_rated);
}

pub struct DvocController {
    v_nom: f32, // nominal voltage (V)
    x_nom: f32, // nominal voltage (p.u.)
    w_nom: f32, // nominal frequency (rad)
    s_rated: f32,  // maximum expected power output (VA)
    dt: f32,  // step size (s)
    pub v: f32,  // voltage state (p.u.)
    pub theta: f32, // angle state (p.u.)
    ki: f32, // Base current (A)
    kv: f32, // Base voltage (V)
    xi: f32,
    c: f32,  // Oscillator capacitance (F)
    l: f32,  // Oscillator inductance (H)
    pub p_ref: f32,  // Active power reference (p.u.)
    pub q_ref: f32,  // Reactive power reference (p.u.)
}

impl DvocController {
    // Steps the dVOC dynamics using a 2nd-order Runge-Kutta method
    // # Arguments
    // * 'u' - alpha-beta current as a tuple of f32 values: (ialpha, ibeta)
    pub fn step(&mut self, u: (f32, f32)) {
        let (dv_dt1, dtheta_dt1) = self.dynamics((self.v, self.theta), u);
        let (dv_dt2, dtheta_dt2) = self.dynamics((self.v + self.dt * dv_dt1, self.theta + self.dt * dtheta_dt1), u);
        let dv_dt = (dv_dt1 + dv_dt2) * 0.5;
        let dtheta_dt = (dtheta_dt1 + dtheta_dt2) * 0.5;
        self.v = self.v + self.dt * dv_dt;
        self.theta = self.theta + self.dt * dtheta_dt;
    }

    // Calculates the voltage dynamics of the dVOC controller using the given input, u.
    // # Arguments    
    // * 'x' - polar voltage (p.u.) as a tuple of f32 values: (v, theta)
    // * 'u' - alpha-beta current (A) as a tuple of f32 values: (ialpha, ibeta)
    fn dynamics(&mut self, x: (f32, f32), u: (f32, f32)) -> (f32, f32) {
        let (v, theta) = x;
        let x_ab = alpha_beta_fr_polar(v * self.kv, theta);
        let u_ab = alpha_beta_fr_ab(u.0, u.1);
        let (p, q) = calc_power(x_ab, u_ab);  // TODO: Determine is this calculation can be done in p.u.
        
        // Per unit dynamics (eq.26 from 'A Grid-compatible Virtual Oscillator Controller')
        let _sqrt2cv = 1. / (SQRT_2 * self.c * v);
        let dv_dt = 2. * self.xi * v * ((self.x_nom * self.x_nom ) - (v * v)) - _sqrt2cv * (q / self.s_rated - self.q_ref);
        let dtheta_dt = 1. - _sqrt2cv / (v * self.w_nom) * (p / self.s_rated - self.p_ref);

        // Unit dynamics (eq.11-12 from 'A Grid-compatible Virtual Oscillator Controller')
        //let kvki_3cv = self.kv * self.ki / (3. * self.c * v);
        //let dv_dt = self.xi / (self.kv * self.kv) * v * (2. * (self.x_nom * self.x_nom ) - 2. * (v * v)) - kvki_3cv * (q - self.q_ref);
        //let dtheta_dt = self.w_nom - kvki_3cv / v * (p - self.p_ref);
        return (dv_dt, dtheta_dt)
    }
}

pub fn build_dvoc_controller(v_nom: f32, w_nom: f32, s_rated: f32, dt: f32, xi: f32, c: f32) -> DvocController {
    DvocController {
        v_nom,
        x_nom: 1.,
        w_nom,
        s_rated,
        dt,
        v: 1.,  
        theta: 0.,
        ki: 3. * v_nom / s_rated,
        kv: v_nom,
        xi,
        c,
        l: 1. / (w_nom * w_nom * c),
        p_ref: 0.,
        q_ref: 0.,
    }
}
