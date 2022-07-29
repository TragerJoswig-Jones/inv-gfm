use super::dvoc::{build_dvoc_controller, build_default_dvoc_controller};

#[test]
fn it_works() {
    let result = 2 + 2;
    assert_eq!(result, 4);
}

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

#[test]
fn test_default_dvoc() {
    let v_nom = 80.;
    let f_nom = 60.;
    let dvoc = build_default_dvoc_controller(v_nom, f_nom);
    assert_eq!(dvoc.v_nom, v_nom);
}