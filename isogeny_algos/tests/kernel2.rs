use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{divpoly, modpoly::Phi};
use isogeny_algos::kernel::kohel::kohel;
use isogeny_algos::kernel::velu::*;
use isogeny_algos::kernel::xonly::*;
use isogeny_algos::poly;
use isogeny_algos::testdata::*;

/// A prime p with a full rational 2-torsion curve is easy: pick roots of a cubic.
fn curve_with_full_2torsion(fp: &Zp, rng: &mut Rng) -> (Curve<u64>, [u64; 3]) {
    // y^2 = (x-r1)(x-r2)(x-r3) with r1+r2+r3 = 0
    loop {
        let (r1, r2) = (fp.random(rng), fp.random(rng));
        let r3 = fp.neg(fp.add(r1, r2));
        if r1 == r2 || r1 == r3 || r2 == r3 {
            continue;
        }
        // x^3 + a x + b with a = r1r2+r1r3+r2r3, b = -r1r2r3
        let a = fp.add(fp.add(fp.mul(r1, r2), fp.mul(r1, r3)), fp.mul(r2, r3));
        let b = fp.neg(fp.mul(r1, fp.mul(r2, r3)));
        let c = Curve::new(a, b);
        if is_smooth(fp, &c) {
            return (c, [r1, r2, r3]);
        }
    }
}

#[test]
fn general_velu_even_degree() {
    let mut rng = Rng::new(21);
    let fp = Zp::new(next_prime(1u64 << 30));
    let phi2 = Phi::compute(&fp, 2);
    for _ in 0..3 {
        let (e, roots) = curve_with_full_2torsion(&fp, &mut rng);
        let j = jinv(&fp, &e);
        // 2-isogenies, one per 2-torsion point
        for &r in &roots {
            let p = Pt::Aff(r, 0);
            let v = velu_general(&fp, &e, &cyclic_reps(&fp, &e, &p, 2));
            assert_eq!(v.deg, 2);
            assert_eq!(phi2.eval(&fp, j, jinv(&fp, &v.cod)), 0, "Phi_2(j, j') != 0");
            assert!(check_homomorphism(&fp, &v, &mut rng, 5));
            // Kohel with kernel polynomial (x - r)
            let h = vec![fp.neg(r), 1];
            let k = kohel(&fp, &e, &h, 2);
            assert_eq!(k.cod, v.cod);
            assert!(check_x_identity(&fp, &k, &mut rng, 5));
            for _ in 0..4 {
                let x = fp.random(&mut rng);
                assert_eq!(k.eval_x(&fp, x), v.eval_x(&fp, x));
            }
        }
        // full 2-torsion kernel: E/E[2] is isomorphic to E (same j), degree 4, non-cyclic
        let gens: Vec<_> = roots.iter().map(|&r| Pt::Aff(r, 0)).collect();
        let reps = subgroup_reps(&fp, &e, &gens, 16).unwrap();
        assert_eq!(reps.len(), 3);
        let v = velu_general(&fp, &e, &reps);
        assert_eq!(v.deg, 4);
        assert_eq!(
            jinv(&fp, &v.cod),
            j,
            "E/E[2] must have the same j-invariant as E"
        );
        assert!(check_homomorphism(&fp, &v, &mut rng, 5));
        let fx: Vec<u64> = vec![e.b, e.a, 0, 1];
        let k = kohel(&fp, &e, &fx, 4);
        assert_eq!(k.cod, v.cod);
        assert!(check_x_identity(&fp, &k, &mut rng, 5));
        for _ in 0..4 {
            let x = fp.random(&mut rng);
            assert_eq!(k.eval_x(&fp, x), v.eval_x(&fp, x));
        }
    }
}

#[test]
fn general_velu_cyclic_order_4_and_6() {
    let mut rng = Rng::new(22);
    let fp = Zp::new(next_prime(1u64 << 28));
    for &n in &[4u64, 6, 8, 9, 10, 12, 15] {
        // random curve with a rational point of exact order n (cyclic)
        let (e, p) = loop {
            let c = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
            if !is_smooth(&fp, &c) {
                continue;
            }
            let ord = order(&fp, &c, &mut rng);
            if !ord.is_multiple_of(n) {
                continue;
            }
            let r = c.random_point(&fp, &mut rng);
            let q = pmul(&fp, &c, &r, (ord / n) as u128);
            // exact order n: n*q = O and (n/p)*q != O for every prime p | n
            let mut ok = pmul(&fp, &c, &q, n as u128) == Pt::Inf;
            for pr in [2u64, 3, 5] {
                if n % pr == 0 && pmul(&fp, &c, &q, (n / pr) as u128) == Pt::Inf {
                    ok = false;
                }
            }
            if ok {
                break (c, q);
            }
        };
        let v = velu_cyclic(&fp, &e, &p, n);
        assert_eq!(v.deg, n);
        assert!(check_homomorphism(&fp, &v, &mut rng, 6), "n = {n}");
        // kernel polynomial = prod over reps
        let xs: Vec<u64> = cyclic_reps(&fp, &e, &p, n)
            .iter()
            .map(|q| if let Pt::Aff(x, _) = q { *x } else { 0 })
            .collect();
        let h = poly::from_roots(&fp, &xs);
        let k = kohel(&fp, &e, &h, n);
        assert_eq!(k.cod, v.cod, "kohel vs velu codomain, n = {n}");
        for _ in 0..4 {
            let x = fp.random(&mut rng);
            assert_eq!(k.eval_x(&fp, x), v.eval_x(&fp, x), "n = {n}");
        }
    }
}

#[test]
fn xonly_velu_matches_and_handles_extension_kernels() {
    let mut rng = Rng::new(23);
    let fp = Zp::new(next_prime(1u64 << 30));
    // (1) rational generator: x-only agrees with point-based Velu
    for &ell in &[3u64, 5, 7, 11, 13, 31] {
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let (_, pts) = kernel_poly_from_point(&fp, &e, &p, ell);
        let v = velu(&fp, &e, &pts, ell);
        let x0 = if let Pt::Aff(x, _) = p {
            x
        } else {
            unreachable!()
        };
        let xo = velu_xonly(&fp, &e, x0, ell).unwrap();
        assert_eq!(xo.cod, v.cod, "l = {ell}");
        for _ in 0..3 {
            let q = e.random_point(&fp, &mut rng);
            assert_eq!(xo.eval(&fp, &q), v.eval(&fp, &q));
        }
    }
    // (2) Galois-stable kernel with non-rational points: kernel polynomial irreducible of degree 2
    //     (l = 5, Frobenius acts on the kernel by a scalar of order 2 modulo +-1)
    let f2 = Zp2::new(fp.p);
    let mut found = 0;
    for _ in 0..200 {
        let c = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
        if !is_smooth(&fp, &c) {
            continue;
        }
        let ks = divpoly::kernel_polys(&fp, &c, 5, &mut rng);
        for h in ks {
            if poly::deg(&fp, &poly::factor_squarefree(&fp, &h, &mut rng)[0]) != 2 {
                continue; // want an irreducible quadratic kernel polynomial
            }
            if poly::factor_squarefree(&fp, &h, &mut rng).len() != 1 {
                continue;
            }
            // roots in F_{p^2}
            let h2: Vec<_> = h.iter().map(|&a| f2.from_u64(a)).collect();
            let roots = poly::roots(&f2, &h2, &mut rng);
            assert_eq!(roots.len(), 2);
            let c2 = Curve::new(f2.from_u64(c.a), f2.from_u64(c.b));
            let xo = velu_xonly(&f2, &c2, roots[0], 5).unwrap();
            let k = kohel(&fp, &c, &h, 5);
            // codomain over F_p^2 must lie in F_p and equal Kohel's
            assert_eq!(
                (xo.cod.a, xo.cod.b),
                ((k.cod.a, 0), (k.cod.b, 0)).into_pair()
            );
            found += 1;
        }
        if found >= 2 {
            break;
        }
    }
    assert!(found >= 1, "no irreducible-quadratic kernel found");
}

trait IntoPair {
    fn into_pair(self) -> ((u64, u64), (u64, u64));
}
impl IntoPair for ((u64, u64), (u64, u64)) {
    fn into_pair(self) -> ((u64, u64), (u64, u64)) {
        self
    }
}

#[test]
fn fast_xonly_velu_matches() {
    use isogeny_algos::fpm::FpM;
    let mut rng = Rng::new(24);
    let fp = Zp::new(next_prime(1u64 << 40));
    for &ell in &[3u64, 5, 7, 11, 31, 101, 401] {
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let x0 = if let Pt::Aff(x, _) = p {
            x
        } else {
            unreachable!()
        };
        let slow = velu_xonly(&fp, &e, x0, ell).unwrap();
        let fast = velu_xonly_fast(&fp, &e, x0, ell).unwrap();
        assert_eq!(slow.cod, fast.cod, "l={ell}");
        let xs1: Vec<u64> = slow.reps.iter().map(|r| r.0).collect();
        let xs2: Vec<u64> = fast.reps.iter().map(|r| r.0).collect();
        assert_eq!(xs1, xs2);
        // FpM<1> agrees too
        let fm = FpM::<1>::from_u64_modulus(fp.p);
        let em = Curve::new(fm.from_u64(e.a), fm.from_u64(e.b));
        let fm_iso = velu_xonly_fast(&fm, &em, fm.from_u64(x0), ell).unwrap();
        assert_eq!(fm.to_canonical(&fm_iso.cod.a)[0], fast.cod.a);
        assert!(check_homomorphism(&fp, &fast, &mut rng, 3));
    }
}

#[test]
fn fast_sqrt_velu_matches_velu() {
    use isogeny_algos::kernel::sqrt_velu::sqrt_velu_fast;
    let mut rng = Rng::new(25);
    let fp = Zp::new(next_prime(1u64 << 32));
    for &ell in &[3u64, 5, 7, 11, 13, 31, 61, 101, 211, 401, 1009, 2003, 10007] {
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let x0 = if let Pt::Aff(x, _) = p {
            x
        } else {
            unreachable!()
        };
        let v = velu_xonly_fast(&fp, &e, x0, ell).unwrap();
        let s = sqrt_velu_fast(&fp, &e, x0, ell).expect("sqrt-velu");
        assert_eq!(s.cod, v.cod, "codomain l={ell}");
        for _ in 0..3 {
            let q = e.random_point(&fp, &mut rng);
            assert_eq!(s.eval(&fp, &q), v.eval(&fp, &q), "eval l={ell}");
        }
    }
}
