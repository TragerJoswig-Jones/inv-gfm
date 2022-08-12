extern crate fixed as fx;
use crate::dvoc::*;

#[test]
fn test_dvoc_params() {
    let v_nom: f32 = 80.;
    let f_nom: f32 = 60.;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let dvoc = build_dvoc_controller_from_flt::<fx::types::I16F16>(v_nom, f_nom, xi, c);
    assert_eq!(dvoc.v_nom, v_nom);
    assert_eq!(dvoc.f_nom, f_nom);
}

#[test]
fn test_default_dvoc() {
    let v_nom = 80.;
    let f_nom = 60.;
    let dvoc = build_default_dvoc_controller::<fx::types::I16F16>(fx::types::I16F16::from_num(v_nom), fx::types::I16F16::from_num(f_nom));
    assert_eq!(dvoc.v_nom, v_nom);
}