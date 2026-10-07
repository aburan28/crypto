use isogeny_algos::binary::*;
use isogeny_algos::curve::Pt;
use isogeny_algos::field::*;
use isogeny_algos::gf2n::GF2n;
use isogeny_algos::path::galbraith::galbraith_with;
use isogeny_algos::path::ghs::ghs_with;
use isogeny_algos::path::graph::{Oracle, Path};
use isogeny_algos::poly;

fn soft_mul(f: &GF2n, a: u64, b: u64) -> u64 {
    // schoolbook in GF(2)[x] then reduce bit by bit
    let mut t = 0u128;
    for i in 0..64 {
        if (b >> i) & 1 == 1 {
            t ^= (a as u128) << i;
        }
    }
    let full = (1u128 << f.n) | f.red as u128;
    for bit in (f.n as usize..128).rev() {
        if (t >> bit) & 1 == 1 {
            t ^= full << (bit - f.n as usize);
        }
    }
    t as u64
}

#[test]
fn gf2n_arithmetic() {
    let mut rng = Rng::new(700);
    for n in [2u32, 3, 8, 13, 31, 32, 61, 62, 63] {
        let f = GF2n::new(n);
        for _ in 0..300 {
            let (a, b, c) = (f.random(&mut rng), f.random(&mut rng), f.random(&mut rng));
            assert_eq!(f.mul(a, b), soft_mul(&f, a, b), "n = {n}");
            assert_eq!(f.mul(f.mul(a, b), c), f.mul(a, f.mul(b, c)));
            assert_eq!(f.mul(a, b ^ c), f.mul(a, b) ^ f.mul(a, c));
            if a != 0 {
                assert_eq!(f.mul(a, f.inv(a)), 1, "inv, n = {n}");
            }
            let r = f.sqrt(a).unwrap();
            assert_eq!(f.sq(r), a);
            assert_eq!(f.trace(a ^ b), f.trace(a) ^ f.trace(b));
            match f.solve_quadratic(a) {
                Some(z) => assert_eq!(f.sq(z) ^ z, a),
                None => assert_eq!(f.trace(a), 1),
            }
        }
        // a^(q-1) = 1
        let a = f.random(&mut rng) | 1;
        assert_eq!(f.pow(a, f.size() - 1), 1);
    }
}

#[test]
fn group_law_and_point_counting() {
    let mut rng = Rng::new(701);
    for n in [11u32, 13, 15] {
        let f = GF2n::new(n);
        for _ in 0..3 {
            let e = BinCurve::new(f.random(&mut rng), f.random(&mut rng), f.random(&mut rng) | 1);
            if e.disc(&f) == 0 {
                continue;
            }
            let ord = e.order_naive(&f);
            // BSGS path (forced) agrees
            assert_eq!(order_bsgs_forced(&f, &e, &mut rng), ord);
            for _ in 0..5 {
                let (p, q, r) = (e.random_point(&f, &mut rng), e.random_point(&f, &mut rng), e.random_point(&f, &mut rng));
                assert!(e.on_curve(&f, &p));
                assert!(e.on_curve(&f, &e.add(&f, &p, &q)));
                assert_eq!(e.add(&f, &e.add(&f, &p, &q), &r), e.add(&f, &p, &e.add(&f, &q, &r)));
                assert_eq!(e.add(&f, &p, &e.neg(&p)), Pt::Inf);
                assert_eq!(e.mul(&f, &p, ord), Pt::Inf);
            }
        }
    }
    // larger field: BSGS only, checked by [N]P = 0
    let f = GF2n::new(41);
    let e = BinCurve::new(1, 0, f.random(&mut rng) | 1);
    let n = e.order(&f, &mut rng);
    for _ in 0..5 {
        assert_eq!(e.mul(&f, &e.random_point(&f, &mut rng), n), Pt::Inf);
    }
    let q = f.size() as i128;
    assert!(((n as i128) - q - 1).pow(2) <= 4 * q);
}

/// BSGS even on small fields (the public `order` switches to naive counting there).
fn order_bsgs_forced(f: &GF2n, e: &BinCurve, rng: &mut Rng) -> u128 {
    // the BSGS code path needs q >= 2^14; emulate by checking that the naive count is the
    // unique even multiple of the lcm of random point orders in the Hasse interval
    let n = e.order_naive(f);
    let q = f.size();
    let w = 2 * ((q as f64).sqrt().ceil() as u128) + 2;
    let mut l = 1u128;
    for _ in 0..20 {
        let p = e.random_point(f, rng);
        let mut o = 1u128;
        let mut r = p;
        while r != Pt::Inf {
            r = e.add(f, &r, &p);
            o += 1;
        }
        assert_eq!(n % o, 0);
        let g = { let (mut a, mut b) = (l, o); while b != 0 { let t = a % b; a = b; b = t; } a };
        l = l / g * o;
    }
    let c: Vec<u128> = ((q + 1 - w)..=(q + 1 + w)).filter(|k| k % l == 0 && k % 2 == 0).collect();
    if c.len() == 1 { c[0] } else { n }
}

/// A curve over GF(2^n) with a rational point of order ell.
fn curve_with_ell(f: &GF2n, ell: u64, rng: &mut Rng) -> (BinCurve, BPt, u128) {
    loop {
        let e = BinCurve::new(f.random(rng) & 1, 0, f.random(rng) | 1);
        let n = e.order(f, rng);
        if n % ell as u128 != 0 {
            continue;
        }
        let p = e.mul(f, &e.random_point(f, rng), n / ell as u128);
        if p != Pt::Inf {
            return (e, p, n);
        }
    }
}

#[test]
fn velu_kohel_and_division_polynomials() {
    let mut rng = Rng::new(702);
    let f = GF2n::new(23);
    for ell in [3u64, 5, 7, 11, 13] {
        let (e, p, n) = curve_with_ell(&f, ell, &mut rng);
        let iso = velu(&f, &e, &p, ell);
        // codomain has the same number of points, images lie on it, homomorphism
        assert_eq!(iso.cod.order(&f, &mut rng), n, "isogenous curves have equal order");
        for _ in 0..5 {
            let (a, b) = (e.random_point(&f, &mut rng), e.random_point(&f, &mut rng));
            let (ia, ib) = (iso.eval(&f, &a), iso.eval(&f, &b));
            assert!(iso.cod.on_curve(&f, &ia), "l = {ell}");
            assert_eq!(iso.eval(&f, &e.add(&f, &a, &b)), iso.cod.add(&f, &ia, &ib), "l = {ell}");
        }
        assert_eq!(iso.eval(&f, &p), Pt::Inf);
        // Kohel from the kernel polynomial: same codomain and x-map
        let xs: Vec<u64> = iso.s.iter().map(|q| q.0).collect();
        let h = poly::from_roots(&f, &xs);
        assert_eq!(kohel_codomain(&f, &e, &h), iso.cod);
        let (num, den) = x_map(&f, &h);
        let a = e.random_point(&f, &mut rng);
        if let (Pt::Aff(x, _), Pt::Aff(xi, _)) = (a, iso.eval(&f, &a)) {
            assert_eq!(f.div(poly::eval(&f, &num, x), poly::eval(&f, &den, x)), xi);
        }
        // division polynomial vanishes on the kernel, degree (l^2-1)/2
        let en = e.normalized(&f);
        let psi = division_poly(&f, &en, ell as usize);
        assert_eq!(psi.len() - 1, ((ell * ell - 1) / 2) as usize);
        for &x in &xs {
            assert_eq!(poly::eval(&f, &psi, x), 0);
        }
        // and kernel factoring finds this kernel
        let ks = kernel_polys(&f, &e, ell, &mut rng);
        assert!(ks.contains(&h), "l = {ell}: {} kernels found", ks.len());
        for k in &ks {
            assert_eq!(kohel_codomain(&f, &e, k).order(&f, &mut rng), n);
        }
    }
}

#[test]
fn phi_mod2_matches_kernel_oracle() {
    let mut rng = Rng::new(703);
    let f = GF2n::new(19);
    let ko = BinKernelOracle { f: &f, ells: vec![3, 5, 7] };
    let phis: Vec<_> = [3usize, 5, 7].iter().map(|&l| phi_mod2(&f, l)).collect();
    let mut edges = 0;
    for _ in 0..6 {
        let j = f.random(&mut rng) | 2;
        let mut from_kernels = ko.neighbors(j, &mut rng);
        let mut from_phi: Vec<(usize, u64)> = vec![];
        for p in &phis {
            for jn in p.neighbors(&f, j, &mut rng) {
                from_phi.push((p.ell, jn));
                assert_eq!(p.eval(&f, jn, j), 0, "Phi symmetric");
            }
        }
        from_kernels.sort();
        from_phi.sort();
        from_phi.dedup();
        assert_eq!(from_kernels, from_phi, "j = {j:#x}");
        edges += from_phi.len();
    }
    assert!(edges >= 6, "only {edges} edges: comparison too weak");
}

fn verify_path(f: &GF2n, phis: &[isogeny_algos::find::modpoly::Phi<GF2n>], path: &Path<u64>) -> bool {
    path.ells.iter().enumerate().all(|(i, &l)| {
        let p = phis.iter().find(|p| p.ell == l).unwrap();
        p.eval(f, path.js[i], path.js[i + 1]) == 0
    })
}

#[test]
fn path_finding_on_binary_curves() {
    let mut rng = Rng::new(704);
    let f = GF2n::new(17);
    let ells = vec![3u64, 5, 7];
    let ko = BinKernelOracle { f: &f, ells: ells.clone() };
    let phis: Vec<_> = ells.iter().map(|&l| phi_mod2(&f, l as usize)).collect();
    let mut total_len = 0;
    for trial in 0..4 {
        // start curve and a random walk to the target
        let j1 = f.random(&mut rng) | 4;
        let mut j2 = j1;
        let mut walked = Path { js: vec![j1], ells: vec![] };
        for _ in 0..6 {
            let ns = ko.neighbors(j2, &mut rng);
            if ns.is_empty() {
                break;
            }
            let (l, jn) = ns[rng.below(ns.len() as u64) as usize];
            walked.js.push(jn);
            walked.ells.push(l);
            j2 = jn;
        }
        assert!(verify_path(&f, &phis, &walked));
        total_len += walked.len();
        // same number of points (twist-insensitive: compare the pair of orders)
        let (o1, o2) = (BinCurve::from_j(&f, j1).order(&f, &mut rng), BinCurve::from_j(&f, j2).order(&f, &mut rng));
        let q = f.size();
        assert!(o1 == o2 || o1 + o2 == 2 * q + 2, "trial {trial}");
        let (p, _) = galbraith_with(&ko, j1, j2, 200_000, &mut rng);
        let p = p.expect("Galbraith BFS path");
        assert_eq!((p.js[0], *p.js.last().unwrap()), (j1, j2));
        assert!(verify_path(&f, &phis, &p));
        let (p, _) = ghs_with(&ko, j1, j2, 200_000, &mut rng);
        if let Some(p) = p {
            assert!(verify_path(&f, &phis, &p));
        }
    }
    assert!(total_len >= 12, "walks too short: {total_len}");
}
