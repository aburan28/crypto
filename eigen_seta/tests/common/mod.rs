#![allow(dead_code)]
#[allow(unused_imports)]
use eigen_seta::field::Field;
use eigen_seta::fp::Fp;
use eigen_seta::fp2::Fp2;
use eigen_seta::params::{ParamSet, Toy16};
use eigen_seta::weier::WCurve;
use num_bigint::BigUint;

pub type Cfg = <Toy16 as ParamSet<4>>::Cfg;
pub type F = Fp2<Cfg, 4>;

pub fn fp_hex(h: &str) -> Fp<Cfg, 4> {
    let n = BigUint::parse_bytes(h.as_bytes(), 16).unwrap();
    let mut l = [0u64; 4];
    for (i, d) in n.to_u64_digits().iter().enumerate() {
        l[i] = *d;
    }
    Fp::from_canonical_limbs(l)
}

pub fn ed_curve() -> WCurve<F> {
    WCurve::new(F::zero(), Fp2::from_fp(fp_hex(Toy16::ED_A4_HEX)), Fp2::from_fp(fp_hex(Toy16::ED_A6_HEX)))
}
