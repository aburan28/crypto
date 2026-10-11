mod common;
use common::*;
use eigen_seta::chain::{graph_identity, map_pair, prodker_step, two_torsion_xs};
use eigen_seta::cm::{eigenvector_2k, mul_vartime, torsion_basis, Omega};
use eigen_seta::kani::{div_dbl, div_dbl_iter, div_on_curve, Glue, Richelot, Split};
use eigen_seta::ntheory::{modpos, sqrt_mod_2power, two_squares, v2};
use eigen_seta::params::{ParamSet, Toy16};
use eigen_seta::weier::{WCurve, WPoint};
use num_bigint::{BigInt, BigUint};
use rand_core::OsRng;

#[test]
fn chain_stepwise() {
    let mut rng = OsRng;
    let ed = ed_curve();
    let om = Omega::find(&ed, &mut rng);
    let big_e = Toy16::E;
    let t = om.trace;
    let mod_e = BigUint::from(1u32) << big_e;
    let s_d = BigInt::from(sqrt_mod_2power(
        &((BigUint::from(1u32) << (big_e + 2)) - 7u32),
        big_e as u64 + 2,
    ));
    let mu = modpos(&((BigInt::from(t) + &s_d) / 2), &mod_e);
    let nu = modpos(&((BigInt::from(t) - &s_d) / 2), &mod_e);
    let w = |x: &WPoint<F>| om.eval(x);
    let p = eigenvector_2k(&ed, &w, &nu, big_e, &mut rng);
    let q = eigenvector_2k(&ed, &w, &mu, big_e, &mut rng);
    let sd_mod = modpos(&s_d, &mod_e);
    let mut found = None;
    'search: for ep in (40..=big_e).rev() {
        let m = (BigUint::from(1u32) << ep) - 7u32;
        if let Some((x, y)) = two_squares(&m, &mut rng, 1000) {
            let (c, cp) = if &x % 2u32 == BigUint::from(0u32) {
                (x, y)
            } else {
                (y, x)
            };
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
    let (e, c, cp) = found.unwrap();
    let n2 = BigUint::from(1u32) << e;
    let lam_s = modpos(&(BigInt::from(sd_mod.clone()) + &c), &n2);
    let lam_t = modpos(&(&c - BigInt::from(sd_mod.clone())), &n2);
    let cpm = modpos(&cp, &n2);
    let shift = BigUint::from(1u32) << (big_e - e);
    let s = mul_vartime(&ed, &p, &shift);
    let tt = mul_vartime(&ed, &q, &shift);
    let k1 = (
        mul_vartime(&ed, &s, &cpm),
        mul_vartime(&ed, &s, &lam_s).neg(),
    );
    let k2 = (
        mul_vartime(&ed, &tt, &cpm),
        mul_vartime(&ed, &tt, &lam_t).neg(),
    );
    let (x1, _x2) = torsion_basis(&ed, 3, Toy16::F3, &mut rng);
    let o = WPoint::infinity();
    let onp = |c1: &WCurve<F>, c2: &WCurve<F>, pr: &(WPoint<F>, WPoint<F>)| {
        bool::from(c1.on_curve(&pr.0)) && bool::from(c2.on_curve(&pr.1))
    };
    // step 1
    let k1a = graph_identity(&ed, &k1);
    let k2a = graph_identity(&ed, &k2);
    let pa = graph_identity(&ed, &(x1, o));
    eprintln!(
        "after graph: on curve {} {} {}",
        onp(&ed, &ed, &k1a),
        onp(&ed, &ed, &k2a),
        onp(&ed, &ed, &pa)
    );
    // check 2-torsion structure of the kernel at step 2
    let a1 = ed.dbl_iter(&k1a.0, e - 2);
    let a2 = ed.dbl_iter(&k2a.0, e - 2);
    let b1 = ed.dbl_iter(&k1a.1, e - 2);
    let b2 = ed.dbl_iter(&k2a.1, e - 2);
    eprintln!(
        "step2 2-torsion inf flags: K1=({},{}) K2=({},{})",
        bool::from(a1.inf),
        bool::from(b1.inf),
        bool::from(a2.inf),
        bool::from(b2.inf)
    );
    let f = prodker_step(&ed, &ed, &k1a, &k2a, e - 2);
    let k1b = map_pair(&f, &k1a);
    let k2b = map_pair(&f, &k2a);
    let pb = map_pair(&f, &pa);
    let (c1, c2) = (f.0.codomain, f.1.codomain);
    eprintln!(
        "after prodker: on curve {} {} {}",
        onp(&c1, &c2, &k1b),
        onp(&c1, &c2, &k2b),
        onp(&c1, &c2, &pb)
    );
    // orders: kernel gens should now have order 2^(e-1)
    let z1 = c1.dbl_iter(&k1b.0, e - 1);
    let z2 = c2.dbl_iter(&k1b.1, e - 1);
    eprintln!(
        "after prodker: [2^(e-1)]K1 inf = ({}, {})",
        bool::from(z1.inf),
        bool::from(z2.inf)
    );
    // glue
    let a1 = c1.dbl_iter(&k1b.0, e - 3);
    let a2 = c1.dbl_iter(&k2b.0, e - 3);
    let b1 = c2.dbl_iter(&k1b.1, e - 3);
    let b2 = c2.dbl_iter(&k2b.1, e - 3);
    eprintln!(
        "glue 2-torsion: y zero {} {} {} {}; inf {} {} {} {}",
        bool::from(a1.y.is_zero()),
        bool::from(a2.y.is_zero()),
        bool::from(b1.y.is_zero()),
        bool::from(b2.y.is_zero()),
        bool::from(a1.inf),
        bool::from(a2.inf),
        bool::from(b1.inf),
        bool::from(b2.inf)
    );
    let ta = two_torsion_xs(&c1, &a1, &a2);
    let tb = two_torsion_xs(&c2, &b1, &b2);
    eprintln!(
        "a-xs distinct: {} {} {}; b-xs distinct: {} {} {}",
        ta[0] != ta[1],
        ta[1] != ta[2],
        ta[0] != ta[2],
        tb[0] != tb[1],
        tb[1] != tb[2],
        tb[0] != tb[2]
    );
    eprintln!(
        "a3 is a root: {}  b3 is a root: {}",
        bool::from(c1.rhs(&ta[2]).is_zero()),
        bool::from(c2.rhs(&tb[2]).is_zero())
    );
    let g = Glue::new(ta, tb);
    eprintln!(
        "glue h zero: {}",
        g.h.iter().all(|c| bool::from(c.is_zero()))
    );
    let d1 = g.map(&k1b.0, &k1b.1);
    let d2 = g.map(&k2b.0, &k2b.1);
    let dp = g.map(&pb.0, &pb.1);
    eprintln!(
        "after glue: divisors on curve {} {} {}",
        div_on_curve(&g.h, &d1),
        div_on_curve(&g.h, &d2),
        div_on_curve(&g.h, &dp)
    );
    // order check of glued kernel: [2^(e-3)] D should be zero divisor: [2^(e-4)]D has order 2
    let d1h = div_dbl_iter(&g.h, &d1, e - 4);
    let d2h = div_dbl_iter(&g.h, &d2, e - 4);
    eprintln!(
        "glued order-2 elements on curve {} {}",
        div_on_curve(&g.h, &d1h),
        div_on_curve(&g.h, &d2h)
    );
    {
        let mut dk = d1;
        for k in 0..(e - 2) {
            let ok = div_on_curve(&g.h, &dk);
            let uz = bool::from(dk.u1.is_zero() & dk.u0.is_zero());
            let vz = bool::from(dk.v1.is_zero() & dk.v0.is_zero());
            if !ok || uz || vz || k < 3 {
                eprintln!(
                    "  dbl {k}: on curve {ok} u zero {uz} v zero {vz} u1 zero {}",
                    bool::from(dk.u1.is_zero())
                );
            }
            if !ok || uz {
                break;
            }
            dk = div_dbl(&g.h, &dk);
        }
    }
    for k in [e - 6, e - 5, e - 4, e - 3] {
        let dk = div_dbl_iter(&g.h, &d1, k);
        eprintln!(
            "  [2^{k}]D1: u zero {} v zero {}",
            bool::from(dk.u1.is_zero() & dk.u0.is_zero()),
            bool::from(dk.v1.is_zero() & dk.v0.is_zero())
        );
    }
    // one Richelot step from those
    let r = Richelot::new(&g.h, &d1h, &d2h);
    let d1r = r.map(&d1);
    let dpr = r.map(&dp);
    eprintln!(
        "after 1st Richelot: on curve {} {}",
        div_on_curve(&r.hnew, &d1r),
        div_on_curve(&r.hnew, &dpr)
    );
    // naive Jacobian phase: r = e - 6 Richelot steps, leaf by direct doubling
    let r_steps = e - 6;
    let mut h = g.h;
    let mut gens = [d1, d2];
    let mut pay = dp;
    for i in 0..r_steps {
        // kernel order now 2^(e-3-i): leaf = [2^(e-4-i)] gens
        let l1 = div_dbl_iter(&h, &gens[0], e - 4 - i);
        let l2 = div_dbl_iter(&h, &gens[1], e - 4 - i);
        let rr = Richelot::new(&h, &l1, &l2);
        gens = [rr.map(&gens[0]), rr.map(&gens[1])];
        pay = rr.map(&pay);
        h = rr.hnew;
        let hz = h.iter().all(|c| bool::from(c.is_zero()));
        if hz || !div_on_curve(&h, &gens[0]) || !div_on_curve(&h, &pay) {
            eprintln!(
                "naive Richelot step {i} (of {r_steps}): h zero {hz}; leaf u1 zero {}",
                bool::from(l1.u1.is_zero() & l1.u0.is_zero())
            );
            break;
        }
    }
    eprintln!(
        "naive Jacobian phase done; gens on curve {}",
        div_on_curve(&h, &gens[0])
    );
    let s1 = div_dbl_iter(&h, &gens[0], 2);
    let s2 = div_dbl_iter(&h, &gens[1], 2);
    eprintln!(
        "split inputs on curve {} {}",
        div_on_curve(&h, &s1),
        div_on_curve(&h, &s2)
    );
    eprintln!("G1 = x^2 + {:?} x + {:?}", s1.u1, s1.u0);
    eprintln!("G2 = x^2 + {:?} x + {:?}", s2.u1, s2.u0);
    eprintln!("h = {:?}", h);
    // the doubled kernel elements should be 2-torsion: [2] s1 = 0  <=>  v = 0 mod u (u | h)
    eprintln!("v of s1: {:?} {:?}", s1.v1, s1.v0);
    let sp = Split::new(&h, &s1, &s2);
    let kk1 = sp.map(&gens[0]);
    let kk2 = sp.map(&gens[1]);
    let pp = sp.map(&pay);
    eprintln!(
        "after split: on curve {} {} {}",
        onp(&sp.e1, &sp.e2, &kk1),
        onp(&sp.e1, &sp.e2, &kk2),
        onp(&sp.e1, &sp.e2, &pp)
    );
    eprintln!(
        "split curves j: {:?} {:?}",
        sp.e1.j_invariant(),
        sp.e2.j_invariant()
    );
    // strategy version
    let mut gens2 = [d1, d2];
    let mut extra = vec![dp];
    let h2 = eigen_seta::chain::jacobian_phase(g.h, &mut gens2, e - 3, r_steps, &mut extra);
    eprintln!(
        "strategy phase: same h as naive {}, gens on curve {}",
        h2 == h,
        div_on_curve(&h2, &gens2[0])
    );
    let _ = div_dbl::<F>;
}
