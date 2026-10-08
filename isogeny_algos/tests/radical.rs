use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::montgomery;
use isogeny_algos::kernel::radical::*;
use isogeny_algos::path::csidh::Csidh;
use isogeny_algos::weier::GWCurve;

/// Radical chains on CSIDH-512 equal the CSIDH action of l_3^k and l_5^k.
#[test]
fn radical_chains_match_csidh_action_on_csidh512() {
    let cs = Csidh::<FpM<8>>::csidh512();
    let f = &cs.fp;
    let mut rng = Rng::new(1200);
    let p1 = f.modulus().add_small(1);
    let e3 = root_exponent(f, 3).expect("p = 2 mod 3");
    let e5 = root_exponent(f, 5).expect("p = 4 mod 5");
    let e7 = root_exponent(f, 7).expect("p != 1 mod 7");
    for (ell, idx) in [(3u64, 0usize), (5, 1), (7, 2)] {
        // start curve: [l_3] or [l_5] applied once from a random curve in the class
        let mut ex = vec![0i32; 74];
        ex[10] = 1;
        ex[20] = -1;
        let a0 = cs.action_fast(f.zero(), &ex, &mut rng);
        let ew = montgomery::to_weierstrass(f, a0);
        let p = loop {
            let r = random_point_f(f, &ew, &mut rng);
            let p = pmul_big(f, &ew, &r, &p1.divrem_small(ell).0);
            if p != Pt::Inf {
                break p;
            }
        };
        let (a1, a2, a3) = tangent_form(f, &ew, &p).unwrap();
        let k = 6;
        let js: Vec<_> = if ell == 3 {
            assert!(f.is_zero(a2), "flex point: a2 = 0");
            let mut cur = (a1, a3);
            let mut js = vec![weierstrass_j(f, model3(f, cur.0, cur.1))];
            for _ in 0..k {
                cur = step3(f, cur.0, cur.1, &e3);
                js.push(weierstrass_j(f, model3(f, cur.0, cur.1)));
            }
            js
        } else if ell == 5 {
            let (b, c) = tate_bc(f, a1, a2, a3);
            assert_eq!(b, c, "order 5: c = b");
            let mut b = b;
            let mut js = vec![weierstrass_j(f, tate_model(f, b, b))];
            for _ in 0..k {
                b = step5(f, b, &e5);
                js.push(weierstrass_j(f, tate_model(f, b, b)));
            }
            js
        } else {
            // order 7: Tate (b,c) lies on b^2 - bc - c^3 = 0, so the X_1(7) parameter is t = b/c
            let (b, c) = tate_bc(f, a1, a2, a3);
            let mut t = f.div(b, c);
            let mut js = vec![weierstrass_j(f, model7(f, t))];
            for _ in 0..k {
                t = step7(f, t, &e7);
                js.push(weierstrass_j(f, model7(f, t)));
            }
            js
        };
        assert_eq!(js[0], montgomery::j_invariant(f, a0));
        // CSIDH: l_ell^s for s = 1..k. The radical chain walks one of l_ell, l_ell^{-1};
        // fix the sign from step 1, then the whole chain must match that direction.
        let dir = {
            let mut e = ex.clone();
            e[idx] = 1;
            let plus = cs.action_fast(f.zero(), &e, &mut rng);
            if js[1] == montgomery::j_invariant(f, plus) {
                1i32
            } else {
                -1i32
            }
        };
        for s in 1..=k {
            let mut e = ex.clone();
            e[idx] = dir * s as i32;
            let a = cs.action_fast(f.zero(), &e, &mut rng);
            assert_eq!(js[s], montgomery::j_invariant(f, a), "l = {ell}, step {s}, dir {dir}");
        }
    }
}

/// Over an ordinary-curve field: consecutive radical steps are l-isogenous and non-backtracking.
#[test]
fn radical_steps_are_isogenies() {
    // a 40-bit prime p = 2 mod 3 and p != 1 mod 5, so cube and fifth roots are unique
    let mut p = (1u64 << 40) + 1;
    while !(is_prime(p) && p % 3 == 2 && p % 5 != 1) {
        p += 2;
    }
    let fp = Zp::new(p);
    let mut rng = Rng::new(1201);
    let mut ran = 0;
    for ell in [3u64, 5] {
        let e_root = root_exponent(&fp, ell).expect("unique roots");
        ran += 1;
        let phi = Phi::compute(&fp, ell as usize);
        let (e, p) = isogeny_algos::testdata::curve_with_point(&fp, ell, &mut rng);
        let (a1, a2, a3) = tangent_form(&fp, &e, &p).unwrap();
        let mut js = vec![];
        if ell == 3 {
            let mut cur = (a1, a3);
            js.push(weierstrass_j(&fp, model3(&fp, cur.0, cur.1)));
            for _ in 0..8 {
                cur = step3(&fp, cur.0, cur.1, &e_root);
                js.push(weierstrass_j(&fp, model3(&fp, cur.0, cur.1)));
            }
        } else {
            let (b, _) = tate_bc(&fp, a1, a2, a3);
            let mut b = b;
            js.push(weierstrass_j(&fp, tate_model(&fp, b, b)));
            for _ in 0..8 {
                b = step5(&fp, b, &e_root);
                js.push(weierstrass_j(&fp, tate_model(&fp, b, b)));
            }
        }
        assert_eq!(js[0], jinv(&fp, &e));
        for i in 0..js.len() - 1 {
            assert_eq!(phi.eval(&fp, js[i], js[i + 1]), 0, "l = {ell}, step {i}");
            if i + 2 < js.len() {
                assert_ne!(js[i], js[i + 2], "non-backtracking");
            }
        }
    }
    assert_eq!(ran, 2);
}

/// Radical 7-isogenies on the X_1(7) family: consecutive curves in the chain are 7-isogenous
/// (independently checked by Phi_7) and the chain is non-backtracking. The 7-th root is unique
/// because the prime is != 1 mod 7.
#[test]
fn radical7_steps_are_isogenies() {
    // a 40-bit prime with p != 1 mod 7 (so gcd(7, p-1) = 1 and the 7th root is unique)
    let mut p = (1u64 << 40) + 1;
    while !(is_prime(p) && p % 7 != 1) {
        p += 2;
    }
    let fp = Zp::new(p);
    let e7 = root_exponent(&fp, 7).expect("unique 7th root");
    let phi = Phi::compute(&fp, 7);
    let mut rng = Rng::new(7007);
    let origin = Pt::Aff(0u64, 0u64);
    let mut ran = 0;
    'outer: for _ in 0..400 {
        let t0 = fp.random(&mut rng);
        // skip degenerate parameters (t in {0,1} and singular members)
        if t0 == 0 || t0 == 1 {
            continue;
        }
        let gw0 = GWCurve::new(model7(&fp, t0));
        if fp.is_zero(gw0.disc(&fp)) || !gw0.on_curve(&fp, &origin) {
            continue;
        }
        // (0,0) must have order exactly 7 on the family member
        if gw0.mul(&fp, &origin, 7) != Pt::Inf {
            continue;
        }
        // build the j-invariant chain
        let mut cur = t0;
        let mut js = vec![jinv_tate7(&fp, cur)];
        for _ in 0..7 {
            cur = step7(&fp, cur, &e7);
            let gw = GWCurve::new(model7(&fp, cur));
            if fp.is_zero(gw.disc(&fp)) {
                continue 'outer; // hit a cusp/degenerate step; try another start
            }
            js.push(jinv_tate7(&fp, cur));
        }
        for i in 0..js.len() - 1 {
            assert_eq!(phi.eval(&fp, js[i], js[i + 1]), 0, "p={p} step {i}: not 7-isogenous");
            if i + 2 < js.len() {
                assert_ne!(js[i], js[i + 2], "p={p} step {i}: backtracking");
            }
        }
        ran += 1;
        if ran >= 4 {
            break;
        }
    }
    assert!(ran >= 4, "expected >= 4 usable radical-7 chains, got {ran}");
}

fn jinv_tate7(fp: &Zp, t: u64) -> u64 {
    weierstrass_j(fp, model7(fp, t))
}
