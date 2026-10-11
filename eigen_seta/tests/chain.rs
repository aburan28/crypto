//! The fixed-schedule Kani chain against a direct evaluation of the embedded
//! endomorphism, on the CM(-7) curve E_D with alpha = 2 omega - t + c.
use eigen_seta::chain::decaps_chain;
use eigen_seta::cm::{chain3_vartime, eigenvector_2k, mul_vartime, point_of_order, torsion_basis, Omega};
use eigen_seta::ntheory::{modpos, sqrt_mod_2power, two_squares, v2};
use eigen_seta::params::{ParamSet, Toy16};
use eigen_seta::weier::{WIso, WPoint};
use num_bigint::{BigInt, BigUint};
use rand_core::{OsRng, RngCore};

mod common;
use common::*;

#[test]
fn chain_matches_direct_endomorphism() {
    let mut rng = OsRng;
    let ed = ed_curve();
    let om = Omega::find(&ed, &mut rng);
    let big_e = Toy16::E;
    let t = om.trace;
    // eigenvalues mu, nu of omega mod 2^E: roots of x^2 - t x + 2
    let mod_e = BigUint::from(1u32) << big_e;
    let s_d = sqrt_mod_2power(&((BigUint::from(1u32) << (big_e + 2)) - 7u32), big_e as u64 + 2);
    let s_d = BigInt::from(s_d);
    let mu = modpos(&((BigInt::from(t) + &s_d) / 2), &mod_e);
    let nu = modpos(&((BigInt::from(t) - &s_d) / 2), &mod_e);
    let w = |x: &WPoint<F>| om.eval(x);
    let p = eigenvector_2k(&ed, &w, &nu, big_e, &mut rng);
    let q = eigenvector_2k(&ed, &w, &mu, big_e, &mut rng);
    assert!(bool::from(om.eval(&p).ct_eq(&mul_vartime(&ed, &p, &mu))), "P not an eigenvector");
    assert!(bool::from(om.eval(&q).ct_eq(&mul_vartime(&ed, &q, &nu))), "Q not an eigenvector");
    // phi : E_D -> E_C of degree d = 3^k; on E_C, psi' = phi theta' phi^ has eigenvalues +-d s_D
    let kdeg = 5usize;
    let d = BigUint::from(3u32).pow(kdeg as u32);
    let k3 = point_of_order(&ed, 3, kdeg, &mut rng);
    let (y1, y2) = torsion_basis(&ed, 3, Toy16::F3, &mut rng);
    let mut pts = [p, q, y1, y2];
    let ec = chain3_vartime(&ed, &k3, kdeg, &mut pts);
    let [pc, qc, x1, x2] = pts;
    let sd_mod = modpos(&(BigInt::from(d.clone()) * &s_d), &mod_e); // d * s_D
    // chain length e' with 2^e' - 7 = c^2 + c'^2 and the canonical schedule
    let mut found = None;
    'search: for ep in (40..=big_e).rev() {
        let m7 = &d * &d * 7u32;
        if (BigUint::from(1u32) << ep) <= m7 { continue; }
        let m = (BigUint::from(1u32) << ep) - m7;
        if let Some((x, y)) = two_squares(&m, &mut rng, 1000) {
            let (c, cp) = if &x % 2u32 == BigUint::from(0u32) { (x, y) } else { (y, x) };
            let n2 = BigInt::from(1u8) << ep;
            for sc in [1i32, -1] {
                for scp in [1i32, -1] {
                    let c = BigInt::from(c.clone()) * sc;
                    let cp = BigInt::from(cp.clone()) * scp;
                    let lam_s = (BigInt::from(sd_mod.clone()) + &c) % &n2;
                    if v2(&(&cp - &lam_s)).max(v2(&(&cp + &lam_s))) == 2 {
                        found = Some((ep, c, cp));
                        break 'search;
                    }
                }
            }
        }
    }
    let (ep, c, cp) = found.expect("no admissible chain length");
    let n2 = BigUint::from(1u32) << ep;
    let lam_s = modpos(&(BigInt::from(sd_mod.clone()) + &c), &n2);
    let lam_t = modpos(&(&c - BigInt::from(sd_mod.clone())), &n2);
    let cpm = modpos(&cp, &n2);
    let shift = BigUint::from(1u32) << (big_e - ep);
    let s = mul_vartime(&ec, &pc, &shift);
    let tt = mul_vartime(&ec, &qc, &shift);
    let k1 = (mul_vartime(&ec, &s, &cpm), mul_vartime(&ec, &s, &lam_s).neg());
    let k2 = (mul_vartime(&ec, &tt, &cpm), mul_vartime(&ec, &tt, &lam_t).neg());
    let o = WPoint::infinity();
    let t0 = std::time::Instant::now();
    let (fc, imgs) = decaps_chain(&ec, &k1, &k2, ep, &[(x1, o), (o, x1), (x2, o)]);
    eprintln!("chain length {ep}: {:?}", t0.elapsed());
    assert_eq!(fc.j_invariant(), ec.j_invariant(), "codomain is not E_C");
    // direct: alpha(phi Y) = phi(d theta'(Y) + c Y),  theta' = 2 omega - t
    let n3 = BigUint::from(3u32).pow(Toy16::F3 as u32);
    let dt = modpos(&(&c - BigInt::from(d.clone()) * t), &n3);
    let dd = d.clone();
    let alpha_y = |y: &WPoint<F>| {
        let w2 = ed.add(&om.eval(y), &om.eval(y));
        ed.add(&mul_vartime(&ed, &w2, &dd), &mul_vartime(&ed, y, &dt))
    };
    let mut ims = [alpha_y(&y1), alpha_y(&y2)];
    let _ = chain3_vartime(&ed, &k3, kdeg, &mut ims);
    let alpha = |x: &WPoint<F>| if bool::from(x.ct_eq(&x1)) { ims[0] } else { ims[1] };
    let cpx1 = mul_vartime(&ec, &x1, &modpos(&cp, &n3));
    let cpx2 = mul_vartime(&ec, &x2, &modpos(&cp, &n3));
    let sig = WIso::find(&fc, &ec);
    let mut ok = false;
    for swap in [false, true] {
        let pick = |pr: &(WPoint<F>, WPoint<F>), first: bool| if first ^ swap { pr.0 } else { pr.1 };
        for s_c in [sig, sig.negate()] {
            if !(bool::from(s_c.eval(&pick(&imgs[0], false)).ct_eq(&cpx1.neg()))
                && bool::from(s_c.eval(&pick(&imgs[2], false)).ct_eq(&cpx2.neg())))
            {
                continue;
            }
            for s_a in [sig, sig.negate()] {
                if bool::from(s_a.eval(&pick(&imgs[1], true)).ct_eq(&cpx1)) {
                    let a1 = s_a.eval(&pick(&imgs[0], true));
                    let a2 = s_a.eval(&pick(&imgs[2], true));
                    assert!(bool::from(a1.ct_eq(&alpha(&x1))), "alpha(X1) mismatch");
                    assert!(bool::from(a2.ct_eq(&alpha(&x2))), "alpha(X2) mismatch");
                    ok = true;
                }
            }
        }
    }
    assert!(ok, "no consistent identification of the codomain factors");
    let _ = rng.next_u32();
}
