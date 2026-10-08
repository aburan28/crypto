use isogeny_algos::curve::{jinv, Curve, Pt};
use isogeny_algos::field::{factor_u64, Field, Rng, Zp};
use isogeny_algos::gf3n::GF3n;
use isogeny_algos::curve::Isogeny;
use isogeny_algos::weier::*;

/// General Weierstrass Vélu agrees with the short-form Vélu in large characteristic.
#[test]
fn gw_velu_matches_short_form_large_char() {
    let fp = Zp::new(1_000_003);
    let mut rng = Rng::new(9);
    let mut tested = 0;
    for _ in 0..200 {
        let a = fp.random(&mut rng);
        let b = fp.random(&mut rng);
        let sc = Curve::new(a, b);
        if !isogeny_algos::curve::is_smooth(&fp, &sc) {
            continue;
        }
        let n = isogeny_algos::curve::order(&fp, &sc, &mut rng);
        let gw = GWCurve::new([0, 0, 0, a, b]);
        for ell in [3u64, 5, 7] {
            if n % ell != 0 {
                continue;
            }
            // point of order ell
            let mut lv = 1;
            while n % (lv * ell) == 0 {
                lv *= ell;
            }
            let k = loop {
                let r = sc.random_point(&fp, &mut rng);
                let mut t = isogeny_algos::curve::pmul(&fp, &sc, &r, (n / lv) as u128);
                let mut ok = None;
                loop {
                    let nt = isogeny_algos::curve::pmul(&fp, &sc, &t, ell as u128);
                    if nt == Pt::Inf {
                        ok = Some(t);
                        break;
                    }
                    if nt == t {
                        break;
                    }
                    t = nt;
                }
                if let Some(k) = ok {
                    if k != Pt::Inf {
                        break k;
                    }
                }
            };
            let short = isogeny_algos::kernel::velu::velu_cyclic(&fp, &sc, &k, ell);
            let gwv = gw_velu_cyclic(&fp, &gw, &k, ell);
            // same codomain j, and the maps agree on a random point
            assert_eq!(jinv(&fp, &short.cod), gwv.cod.j(&fp));
            let r = sc.random_point(&fp, &mut rng);
            if let (Pt::Aff(sx, _), Pt::Aff(gx, _)) = (short.eval(&fp, &r), gwv.eval(&fp, &r)) {
                assert_eq!(sx, gx);
                tested += 1;
            }
        }
    }
    assert!(tested >= 10, "only {tested}");
}

/// Ordinary curves y^2 = x^3 + a2 x^2 + a6 over GF(3^n): Vélu codomain is a homomorphic image
/// with equal point count; the kernel maps to the identity.
#[test]
fn gw_velu_characteristic_three() {
    for n in [3u32, 4, 5] {
        let f = GF3n::new(n);
        let mut rng = Rng::new(300 + n as u64);
        let q = 3u64.pow(n);
        let points = |e: &GWCurve<u64>| -> Vec<Pt<u64>> {
            let mut v = vec![Pt::Inf];
            for x in 0..q {
                for y in 0..q {
                    let pt = Pt::Aff(x, y);
                    if e.on_curve(&f, &pt) {
                        v.push(pt);
                    }
                }
            }
            v
        };
        let mut done = 0;
        for _ in 0..200 {
            if done >= 3 {
                break;
            }
            let a2 = f.random(&mut rng);
            let a6 = f.random(&mut rng);
            if a2 == 0 || a6 == 0 {
                continue;
            }
            let e = GWCurve::new([0, a2, 0, 0, a6]);
            if f.is_zero(e.disc(&f)) {
                continue;
            }
            let pts = points(&e);
            let order = pts.len() as u64;
            for ell in [5u64, 7, 11, 13] {
                if order % ell != 0 || done >= 3 {
                    continue;
                }
                let k = loop {
                    let idx = rng.below(order);
                    let kk = e.mul(&f, &pts[idx as usize], order / ell);
                    if kk != Pt::Inf {
                        break kk;
                    }
                };
                let v = gw_velu_cyclic(&f, &e, &k, ell);
                assert!(v.cod.on_curve(&f, &Pt::Inf));
                assert_eq!(v.eval(&f, &k), Pt::Inf, "kernel -> O");
                let cpts = points(&v.cod);
                assert_eq!(cpts.len() as u64, order, "isogenous: equal point count");
                let mut checks = 0;
                for _ in 0..120 {
                    let p = pts[rng.below(order) as usize];
                    let qq = pts[rng.below(order) as usize];
                    let (ip, iq) = (v.eval(&f, &p), v.eval(&f, &qq));
                    if ip == Pt::Inf || iq == Pt::Inf {
                        continue;
                    }
                    assert!(v.cod.on_curve(&f, &ip), "image on codomain");
                    let s = v.eval(&f, &e.add(&f, &p, &qq));
                    assert_eq!(s, v.cod.add(&f, &ip, &iq), "homomorphism");
                    checks += 1;
                }
                assert!(checks >= 15, "only {checks} non-kernel samples");
                let _ = factor_u64;
                done += 1;
            }
        }
        assert!(done >= 1, "n = {n}: no ell-divisible order found");
    }
}
