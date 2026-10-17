//! Tests for the parts of Bernstein–De Feo–Leroux–Smith 2020 added in `strassen`,
//! `kernel::sqrt_mult`, `kernel::sqrt_velu_irr` and `path::csurf`. Every check is against an
//! independent computation: naive products, trial division / Pollard rho (`field::factor_u64`),
//! integer cyclotomic polynomials, Vélu over the extension field, Weierstrass Vélu, and the
//! modular polynomials Φ_l. Timings are printed (run with `--nocapture`).
use isogeny_algos::bigint::Big;
use isogeny_algos::curve::{jinv, Curve, Pt};
use isogeny_algos::ext::Fp4;
use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::fpm::FpM;
use isogeny_algos::fpr::FpR;
use isogeny_algos::kernel::montgomery::*;
use isogeny_algos::kernel::sqrt_mult::*;
use isogeny_algos::kernel::sqrt_velu_irr::*;
use isogeny_algos::kernel::sqrt_velu_mont::SqrtVeluMont;
use isogeny_algos::kernel::velu::velu_general;
use isogeny_algos::kernel::xonly::velu_xonly;
use isogeny_algos::path::csidh::Csidh;
use isogeny_algos::path::csurf::*;
use isogeny_algos::poly;
use isogeny_algos::strassen::*;
use std::time::Instant;

// ------------------------------------------------------------------------------------------
// §2 Strassen

#[test]
fn strassen_factorial_matches_naive() {
    let mut rng = Rng::new(1);
    let mut moduli: Vec<u64> = vec![
        2,
        3,
        4,
        6,
        8,
        12,
        1 << 20,
        1000003,
        1000003 * 999983,
        (1u64 << 61) - 1,
        3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43,
        (1u64 << 62) - 57,
    ];
    for _ in 0..12 {
        moduli.push(2 + rng.below((1u64 << 62) - 2));
    }
    let mut ls: Vec<u64> = vec![
        0, 1, 2, 3, 4, 15, 16, 17, 63, 64, 65, 99, 100, 101, 1000, 4096,
    ];
    for _ in 0..6 {
        ls.push(rng.below(20000));
    }
    let mut checked = 0;
    for &n in &moduli {
        for &l in &ls {
            let want = factorial_naive(l, n);
            // min_b = 2 forces the trees even for small l
            assert_eq!(factorial_mod_with(l, n, 2), want, "l = {l}, n = {n}");
            assert_eq!(factorial_mod(l, n), want, "l = {l}, n = {n}");
            checked += 1;
        }
    }
    // Wilson: (p-1)! = -1 mod p, an independent closed form (p ~ 2^20)
    for &p in &[1048573u64, 1048571] {
        assert!(is_prime(p));
        assert_eq!(factorial_mod(p - 1, p), p - 1, "Wilson p = {p}");
    }
    println!("strassen factorial: {checked} (l, n) pairs = naive");
}

#[test]
fn strassen_factorial_timing() {
    let n = (1u64 << 61) - 1;
    println!("l          naive (ms)   strassen (ms)");
    for &l in &[1_000u64, 10_000, 100_000, 1_000_000, 4_000_000] {
        let t = Instant::now();
        let a = factorial_naive(l, n);
        let tn = t.elapsed().as_secs_f64() * 1e3;
        let t = Instant::now();
        let b = factorial_mod_with(l, n, 2);
        let ts = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!(a, b);
        println!("{l:<10} {tn:>10.3}   {ts:>10.3}");
    }
}

#[test]
fn rem_tree_matches_horner_and_evaltree() {
    let fp = Zp::new((1u64 << 61) - 1);
    let mut rng = Rng::new(2);
    for &m in &[1usize, 2, 3, 7, 64, 300, 2048] {
        let pts: Vec<u64> = (0..m).map(|_| fp.random(&mut rng)).collect();
        let a: Vec<u64> = (0..2 * m + 5).map(|_| fp.random(&mut rng)).collect();
        let t0 = Instant::now();
        let rt = RemTree::from_roots(&fp, &pts);
        let v1 = rt.eval(&fp, &a);
        let t_rem = t0.elapsed().as_secs_f64() * 1e3;
        let t0 = Instant::now();
        let et = poly::EvalTree::new(&fp, &pts);
        let a_red = poly::rem(&fp, &a, &et.tree.last().unwrap()[0]);
        let v2 = et.eval(&fp, &a_red);
        let t_ev = t0.elapsed().as_secs_f64() * 1e3;
        for k in 0..m {
            assert_eq!(v1[k], poly::eval(&fp, &a, pts[k]), "m = {m}, k = {k}");
            assert_eq!(v1[k], v2[k]);
        }
        if m >= 300 {
            println!(
                "m = {m}: RemTree {t_rem:.2} ms, poly::EvalTree {t_ev:.2} ms (build + one eval)"
            );
        }
    }
}

#[test]
fn strassen_factoring_matches_trial_division() {
    let mut rng = Rng::new(3);
    let mut ns: Vec<u64> = vec![
        2,
        3,
        4,
        9,
        49,
        1 << 39,
        3u64.pow(25),
        1000003 * 999983,
        1048573 * 1048571,
        2 * 3 * 5 * 7 * 11 * 13 * 1000003,
        104729 * 104729 * 13,
    ];
    // balanced semiprimes (smallest factor ~ 2^19-2^20: the slowest case) and random numbers
    for _ in 0..3 {
        let p = next_prime((1 << 19) + rng.below(1 << 19));
        let q = next_prime((1 << 19) + rng.below(1 << 19));
        ns.push(p * q);
    }
    for _ in 0..15 {
        ns.push(2 + rng.below((1u64 << 40) - 2));
    }
    let t = Instant::now();
    for &n in &ns {
        assert_eq!(factor(n), factor_u64(n), "n = {n}");
        // the first prime by trial division, independently of factor_u64
        let mut d = 2u64;
        while d * d <= n && n % d != 0 {
            d += 1;
        }
        let spf = if d * d <= n { d } else { n };
        assert_eq!(smallest_factor(n), spf, "n = {n}");
    }
    println!(
        "strassen factoring: {} numbers <= 2^40 = factor_u64 in {:.2} s",
        ns.len(),
        t.elapsed().as_secs_f64()
    );
}

// ------------------------------------------------------------------------------------------
// §3 Algorithm 1

fn check_alg1<F: Field>(f: &F, rng: &mut Rng, label: &str) {
    for &(m, r, n) in &[
        (0i64, 1i64, 1u64),
        (1, 1, 2),
        (1, 1, 10),
        (0, 1, 100),
        (5, 3, 257),
        (-40, 7, 1000),
        (1, 2, 4097),
    ] {
        let zeta = f.random(rng);
        let alg = Alg1::progression(f, zeta, m, r, n);
        let s = alg.set();
        let want_s: Vec<i64> = (0..n as i64).map(|t| m + t * r).collect();
        let mut ws = want_s.clone();
        ws.sort_unstable();
        assert_eq!(s, ws, "{label}: S = I + J ∪ K");
        assert!(no_common_differences(&alg.i, &alg.j));
        for _ in 0..3 {
            let alpha = f.random(rng);
            let d = hs_direct(f, zeta, &want_s, alpha);
            assert_eq!(
                alg.eval(f, alpha),
                d,
                "{label}: Alg 1, (m, r, n) = ({m}, {r}, {n})"
            );
            assert_eq!(alg.eval_pollard(f, alpha), d, "{label}: §3.2, n = {n}");
        }
    }
    // an arbitrary set: random gaps
    let mut s = vec![];
    let mut x = 0i64;
    for _ in 0..600 {
        x += 1 + (rng.below(5) == 0) as i64 * rng.below(9) as i64;
        s.push(x);
    }
    let zeta = f.random(rng);
    let alg = Alg1::for_set(f, zeta, &s, &[1, 2, 3]);
    assert_eq!(alg.set(), s);
    let alpha = f.random(rng);
    assert_eq!(
        alg.eval(f, alpha),
        hs_direct(f, zeta, &s, alpha),
        "{label}: general S"
    );
}

#[test]
fn alg1_matches_direct_product() {
    let mut rng = Rng::new(4);
    check_alg1(&Zp::new((1u64 << 61) - 1), &mut rng, "Zp 61");
    check_alg1(&Zp::new(1009), &mut rng, "Zp 1009");
    let m127 = Big::from_u64(1).shl(127).sub_small(1);
    check_alg1(&FpM::<2>::new(&m127), &mut rng, "FpM<2> 2^127-1");
    let cs = Csidh::<FpM<8>>::csidh512();
    check_alg1(&cs.fp, &mut rng, "FpM<8> CSIDH-512");
}

#[test]
fn alg1_cyclotomic() {
    // Phi_n(X) = h_S(X) for S = (Z/nZ)^*, zeta of exact order n in F_p (p = 1 mod n)
    let mut rng = Rng::new(5);
    for &n in &[1u64, 2, 7, 12, 64, 101, 105, 210, 360, 1024, 1155] {
        let mut p = n * ((1u64 << 40) / n) + 1;
        while !is_prime(p) {
            p += n;
        }
        let fp = Zp::new(p);
        let fac = factor_u64(n);
        let zeta = loop {
            let g = 2 + rng.below(p - 3);
            let z = fp.pow(g, ((p - 1) / n) as u128);
            if fac.iter().all(|&(q, _)| fp.pow(z, (n / q) as u128) != 1) {
                break z;
            }
        };
        let units: Vec<i64> = (0..n as i64)
            .filter(|&i| gcd(i as u64, n) == 1 || n == 1)
            .collect();
        let alg = Alg1::for_set(&fp, zeta, &units, &[1, 2, 3, 6]);
        assert_eq!(alg.set(), units);
        let phi = cyclotomic_int(n);
        assert_eq!(phi.len() - 1, units.len(), "deg Phi_n = phi(n)");
        let phi_p: Vec<u64> = phi
            .iter()
            .map(|&c| c.rem_euclid(p as i128) as u64)
            .collect();
        for _ in 0..4 {
            let alpha = fp.random(&mut rng);
            assert_eq!(
                alg.eval(&fp, alpha),
                poly::eval(&fp, &phi_p, alpha),
                "Phi_{n}"
            );
        }
        // roots: Phi_n(zeta^u) = 0 for units u
        assert_eq!(
            alg.eval(&fp, fp.pow(zeta, units[units.len() / 2] as u128)),
            0
        );
        println!(
            "Phi_{n}: #S = {}, #I = {}, #J = {}, #K = {}",
            units.len(),
            alg.i.len(),
            alg.j.len(),
            alg.k.len()
        );
    }
}

#[test]
fn alg1_timing() {
    let fp = Zp::new((1u64 << 61) - 1);
    let mut rng = Rng::new(6);
    for &n in &[10_000u64, 100_000, 1_000_000] {
        let zeta = fp.random(&mut rng);
        let alpha = fp.random(&mut rng);
        let t = Instant::now();
        let alg = Alg1::progression(&fp, zeta, 1, 1, n);
        let t_setup = t.elapsed().as_secs_f64() * 1e3;
        let t = Instant::now();
        let v = alg.eval(&fp, alpha);
        let t_eval = t.elapsed().as_secs_f64() * 1e3;
        // direct product by successive powers (n multiplications for powers, n for the product)
        let t = Instant::now();
        let (mut acc, mut z) = (1u64, zeta);
        for _ in 0..n {
            acc = fp.mul(acc, fp.sub(alpha, z));
            z = fp.mul(z, zeta);
        }
        let t_dir = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!(v, acc);
        println!(
            "n = {n}: Alg 1 setup {t_setup:.2} ms + eval {t_eval:.2} ms; direct {t_dir:.2} ms"
        );
    }
}

// ------------------------------------------------------------------------------------------
// §4.14 irrational generators

/// Order class of lambda in (Z/lZ)^* / {±1}.
fn class_order(lambda: u64, ell: u64) -> usize {
    let mut x = lambda % ell;
    let mut k = 1;
    while x != 1 && x != ell - 1 {
        x = x * lambda % ell;
        k += 1;
    }
    k
}

fn powmod(b: u64, mut e: u64, m: u64) -> u64 {
    let (mut r, mut b) = (1u128, (b % m) as u128);
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m as u128;
        }
        b = b * b % m as u128;
        e >>= 1;
    }
    r as u64
}

/// A Montgomery curve E_A/F_p whose Frobenius has an eigenvalue lambda mod l of class order e
/// (the other eigenvalue mu with mu^e != lambda^e up to sign, so the lambda-line is the whole
/// l-torsion of E or its twist over F_{p^e}), with its trace t.
fn curve_with_eigen(fp: &Zp, ell: u64, e: usize, rng: &mut Rng) -> (u64, u64, u64, i64) {
    let p = fp.p;
    loop {
        let a = fp.random(rng);
        if a == 2 || a == p - 2 {
            continue;
        }
        let w = to_weierstrass(fp, a);
        let n = isogeny_algos::curve::order(fp, &w, rng);
        let t = p as i64 + 1 - n as i64;
        let tm = t.rem_euclid(ell as i64) as u64;
        for lam in 2..ell - 1 {
            let v = (lam * lam % ell + ell * ell - tm * lam % ell + p % ell) % ell;
            if v != 0 || class_order(lam, ell) != e {
                continue;
            }
            let mu = p % ell * powmod(lam, ell - 2, ell) % ell;
            let (le, me) = (powmod(lam, e as u64, ell), powmod(mu, e as u64, ell));
            if le == me {
                continue;
            }
            return (a, lam, mu, t);
        }
    }
}

/// x(P) in G for P of order l in the lambda-eigenline (on E or the twist over F_{p^e}).
fn eigen_point<G: ExtField<B = Zp>>(
    g: &G,
    p: u64,
    a: u64,
    t: i64,
    ell: u64,
    lam: u64,
    e: usize,
    rng: &mut Rng,
) -> G::E {
    // t_k: trace of pi^k
    let (mut t0, mut t1) = (2i128, t as i128);
    for _ in 1..e {
        (t0, t1) = (t1, t as i128 * t1 - p as i128 * t0);
    }
    let pe = (p as i128).pow(e as u32);
    let on_e = powmod(lam, e as u64, ell) == 1;
    let n = if on_e { pe + 1 - t1 } else { pe + 1 + t1 } as u128;
    let mut cof = n;
    while cof.is_multiple_of(ell as u128) {
        cof /= ell as u128;
    }
    let al = g.lift(a);
    let e24 = a24(g, al);
    let half = g.q().sub_small(1).shr(1);
    loop {
        let x = g.random(rng);
        let r = g.add(g.mul(g.add(g.mul(x, x), g.mul(al, x)), x), x);
        if g.is_zero(r) || (g.pow_big(r, &half) == g.one()) != on_e {
            continue;
        }
        let mut k = ladder_xz(g, e24, (x, g.one()), &Big::from_u128(cof));
        if g.is_zero(k.1) {
            continue;
        }
        loop {
            let k2 = ladder_xz(g, e24, k, &Big::from_u64(ell));
            if g.is_zero(k2.1) {
                break;
            }
            k = k2;
        }
        return g.div(k.0, k.1);
    }
}

/// Checks the §4.14 √élu on <P> against Montgomery Vélu over G and Weierstrass x-only Vélu over G.
fn check_irr<G: ExtField<B = Zp>>(
    g: &G,
    a: u64,
    xp: G::E,
    ell: u64,
    lam: u64,
    rng: &mut Rng,
    label: &str,
) {
    let fp = g.base_field();
    let found = find_eigenvalue(g, a, xp, ell).expect("eigenvalue");
    assert!(
        found == lam || found == ell - lam,
        "{label}: eigenvalue class"
    );
    assert!(g.down(xp).is_none(), "{label}: x(P) must be irrational");
    let t = Instant::now();
    let sv = SqrtVeluIrr::new(g, a, xp, ell, canon(lam as i64, ell)).expect("setup");
    let a_new = sv.codomain();
    let t_sv = t.elapsed().as_secs_f64() * 1e3;
    // reference: Vélu over G from all (l-1)/2 multiples
    let t = Instant::now();
    let al = g.lift(a);
    let ms = multiples(g, a24(g, al), (xp, g.one()), ((ell - 1) / 2) as usize);
    let a_ref = g
        .down(velu_codomain(g, al, &ms, ell))
        .expect("codomain in F_p");
    let t_ref = t.elapsed().as_secs_f64() * 1e3;
    assert_eq!(a_new, a_ref, "{label}: codomain");
    assert_eq!(sv.codomain_orbits(), a_ref, "{label}: orbit-norm codomain");
    for _ in 0..3 {
        let u = fp.random(rng);
        let want = g.down(isog_x(g, &ms, g.lift(u))).expect("image in F_p");
        assert_eq!(sv.eval(u), want, "{label}: x-map");
    }
    // Weierstrass x-only Vélu over G (division-polynomial multiples), compared by j
    let w = to_weierstrass(g, al);
    let xw = g.add(xp, g.div(al, g.from_u64(3)));
    let vw = velu_xonly(g, &w, xw, ell).expect("velu_xonly");
    assert_eq!(
        jinv(g, &vw.cod),
        g.lift(j_invariant(&fp, a_new)),
        "{label}: j vs Weierstrass Vélu"
    );
    let (ni, nj, nk) = sv.sizes();
    println!(
        "{label}: l = {ell}, e = {}, #I = {ni}, #J = {nj}, #K = {nk} (#S = {}); √élu-irr {t_sv:.2} ms, Vélu over F_p^e {t_ref:.2} ms",
        sv.sys.e,
        (ell - 1) / 2
    );
}

#[test]
fn sqrt_velu_irrational_generators() {
    let mut rng = Rng::new(7);
    let p = next_prime(1 << 19);
    let fp = Zp::new(p);
    // e = 2 over Zp2 (lambda^2 = -1 needs l = 1 mod 4)
    let g2 = Zp2::new(p);
    for &ell in &[101u64, 197] {
        let (a, lam, _mu, t) = curve_with_eigen(&fp, ell, 2, &mut rng);
        let xp = eigen_point(&g2, p, a, t, ell, lam, 2, &mut rng);
        check_irr(&g2, a, xp, ell, lam, &mut rng, "Zp2");
    }
    // e = 3 over FpR (l = 1 mod 3)
    let g3 = FpR::new(p, 3);
    for &ell in &[103u64, 151] {
        let (a, lam, _mu, t) = curve_with_eigen(&fp, ell, 3, &mut rng);
        let xp = eigen_point(&g3, p, a, t, ell, lam, 3, &mut rng);
        check_irr(&g3, a, xp, ell, lam, &mut rng, "FpR r=3");
    }
    // e = 4 over Fp4 (lambda of order 8: l = 1 mod 8)
    let g4 = Fp4::new(p);
    let ell = 113u64;
    let (a, lam, _mu, t) = curve_with_eigen(&fp, ell, 4, &mut rng);
    let xp = eigen_point(&g4, p, a, t, ell, lam, 4, &mut rng);
    check_irr(&g4, a, xp, ell, lam, &mut rng, "Fp4");
}

#[test]
fn irrational_index_systems_partition_s() {
    // I ± J and K partition (Z/lZ)^*/{±1} exactly, for every class order e of lambda
    for &ell in &[101u64, 103, 211, 1009, 2003] {
        let mut seen_e = std::collections::BTreeSet::new();
        for lam in 2..ell - 1 {
            let e = class_order(lam, ell);
            if !seen_e.insert(e) {
                continue;
            }
            let sys = index_system(ell, lam);
            let orbs = orbits(ell, lam);
            let orbit_of = |c: u64| orbs.iter().find(|o| o.contains(&c)).unwrap().clone();
            let full =
                |reps: &[u64]| -> Vec<u64> { reps.iter().flat_map(|&c| orbit_of(c)).collect() };
            let (fi, fj, fk) = (full(&sys.i), full(&sys.j), full(&sys.k));
            let mut cover: Vec<u64> = fk.clone();
            for &i in &fi {
                for &j in &fj {
                    cover.push(canon(i as i64 + j as i64, ell));
                    cover.push(canon(i as i64 - j as i64, ell));
                }
            }
            cover.sort_unstable();
            let want: Vec<u64> = (1..=(ell - 1) / 2).collect();
            assert_eq!(cover, want, "l = {ell}, lambda = {lam}, e = {e}");
        }
    }
}

// ------------------------------------------------------------------------------------------
// Appendix A.5

#[test]
fn a5_fused_codomain_and_image() {
    let cs = Csidh::<FpM<8>>::csidh512();
    let f = &cs.fp;
    let mut rng = Rng::new(8);
    let a = cs.action_fast(
        f.zero(),
        &{
            let mut e = vec![0i32; 74];
            e[1] = 1;
            e
        },
        &mut rng,
    );
    for &ell in &[3u64, 7, 101, 373, 587] {
        let cof = cs.p1.divrem_small(ell).0;
        let k = loop {
            let x = f.random(&mut rng);
            let k = ladder_xz(f, a24(f, a), (x, f.one()), &cof);
            if !f.is_zero(k.1) {
                break k;
            }
        };
        let ms = multiples(f, a24(f, a), k, ((ell - 1) / 2) as usize);
        let sv = SqrtVeluMont::new(f, a, k, ell).unwrap();
        let u = f.random(&mut rng);
        let (a_new, x_new) = fused_codomain_eval(f, &sv, u);
        assert_eq!(
            a_new,
            velu_codomain(f, a, &ms, ell),
            "fused codomain, l = {ell}"
        );
        assert_eq!(x_new, isog_x(f, &ms, u), "fused image, l = {ell}");
        if ell == 587 {
            let reps = 20;
            let t = Instant::now();
            for _ in 0..reps {
                std::hint::black_box(fused_codomain_eval(f, &sv, u));
            }
            let tf = t.elapsed().as_secs_f64() * 1e3 / reps as f64;
            let t = Instant::now();
            for _ in 0..reps {
                std::hint::black_box((sv.codomain(f), sv.eval(f, u)));
            }
            let ts = t.elapsed().as_secs_f64() * 1e3 / reps as f64;
            println!("A.5 merge at l = 587 (CSIDH-512): fused {tf:.3} ms, separate {ts:.3} ms");
        }
    }
}

// ------------------------------------------------------------------------------------------
// §6 CSURF

/// Whether two short Weierstrass curves are isomorphic over F_p (not just twists).
fn fp_isomorphic<F: Field>(f: &F, c1: &Curve<F::E>, c2: &Curve<F::E>) -> bool {
    if jinv(f, c1) != jinv(f, c2) {
        return false;
    }
    let qm1 = f.q().sub_small(1);
    // r is a k-th power in F_q^* iff r^((q-1)/gcd(k, q-1)) = 1
    let is_kth_power = |r: F::E, k: u64| {
        let g = gcd(k, qm1.rem_small(k));
        f.pow_big(r, &qm1.divrem_small(g).0) == f.one()
    };
    if f.is_zero(c1.b) {
        return is_kth_power(f.div(c2.a, c1.a), 4);
    }
    if f.is_zero(c1.a) {
        return is_kth_power(f.div(c2.b, c1.b), 6);
    }
    // u^2 = (b2 a1) / (b1 a2), u^4 = a2 / a1
    let u2 = f.div(f.mul(c2.b, c1.a), f.mul(c1.b, c2.a));
    f.sq(u2) == f.div(c2.a, c1.a) && is_kth_power(u2, 2)
}

fn small_csurf() -> Csurf<Zp> {
    (10..16)
        .rev()
        .find_map(Csurf::<Zp>::with_n_primes)
        .expect("a CSURF prime below 2^62")
}

/// A 2-isogeny from E_A^- by Vélu from the explicit kernel point (x0, 0), as a Weierstrass curve.
fn velu2<F: Field>(f: &F, a: F::E, x0: F::E) -> Curve<F::E> {
    let w = minus_to_weierstrass(f, a);
    let xw = f.add(x0, f.div(a, f.from_u64(3)));
    velu_general(f, &w, &[Pt::Aff(xw, f.zero())]).cod
}

fn check_two_steps<F: Field>(
    cs: &Csurf<F>,
    phi2: Option<&Phi<F>>,
    start: F::E,
    n: usize,
    label: &str,
) {
    let f = &cs.fp;
    let mut a = start;
    let mut alpha = cs.alpha_of(a);
    let mut prev_j = None;
    for step in 0..n {
        // the other non-zero 2-torsion point is not halvable
        let other = f.neg(f.inv(alpha));
        assert!(
            !(cs.is_square(other) && cs.is_square(f.add(other, f.inv(other)))),
            "{label}: both (alpha,0), (-1/alpha,0) halvable"
        );
        let (a2, alpha2) = cs.two_step(alpha);
        assert_eq!(
            alpha2,
            cs.alpha_of(a2),
            "{label}: next kernel from the radical, step {step}"
        );
        // generic 2-isogeny by Vélu with kernel (alpha, 0): F_p-isomorphic to E_{A'}^-
        let cod = velu2(f, a, alpha);
        assert!(
            fp_isomorphic(f, &cod, &minus_to_weierstrass(f, a2)),
            "{label}: radical 2-step = Vélu, step {step}"
        );
        let (j, j2) = (j_minus(f, a), j_minus(f, a2));
        if let Some(phi) = phi2 {
            assert!(f.is_zero(phi.eval(f, j, j2)), "{label}: Phi_2(j, j') = 0");
        }
        // non-backtracking (the l_2 direction is consistent along the chain)
        if let Some(pj) = prev_j {
            assert_ne!(j2, pj, "{label}: l_2 chain backtracks at step {step}");
        }
        // the vertical isogeny (kernel (0,0)) leaves the surface: codomain has one rational
        // 2-torsion point
        let down = velu2(f, a, f.zero());
        let cubic = vec![down.b, down.a, f.zero(), f.one()];
        let mut r2 = Rng::new(step as u64);
        let roots = if f.char() < (1u64 << 62) && f.q().bits() <= 62 {
            Some(poly::roots(f, &cubic, &mut r2).len())
        } else {
            None
        };
        if let Some(nr) = roots {
            assert_eq!(nr, 1, "{label}: kernel (0,0) is vertical");
        }
        prev_j = Some(j);
        (a, alpha) = (a2, alpha2);
    }
    // l_2 l_2-bar = (2): back to the start
    let fwd = cs.two_power(start, n as i32);
    assert_eq!(
        cs.two_power(fwd, -(n as i32)),
        start,
        "{label}: l_2^n l_2bar^n = 1"
    );
}

#[test]
fn csurf_two_steps_radical_vs_velu() {
    let cs = small_csurf();
    let f = &cs.fp;
    assert_eq!(cs.p() % 8, 7);
    println!(
        "CSURF test prime p = {} ({} bits), odd primes {:?}",
        cs.p(),
        64 - cs.p().leading_zeros(),
        cs.primes
    );
    let phi2 = Phi::compute(f, 2);
    check_two_steps(&cs, Some(&phi2), f.zero(), 30, "small p from A = 0");
    // from a random surface curve
    let mut rng = Rng::new(9);
    let mut e = vec![0i32; cs.primes.len() + 1];
    e[1] = 1;
    e[3] = -2;
    let a = cs.action(f.zero(), &e, &mut rng);
    check_two_steps(&cs, Some(&phi2), a, 30, "small p random");
}

/// p = 8 * (consecutive odd primes) * l_last - 1 with the base just below `min_bits` - 12 bits
/// and l_last < 2^12 the first larger prime that makes p prime (at most 64 N bits).
fn find_big_csurf<const N: usize>(min_bits: usize) -> Csurf<FpM<N>> {
    let odd: Vec<u64> = (3u64..4096).filter(|&q| is_prime(q)).collect();
    let mut base = vec![];
    let mut bits = 3.0f64;
    let mut idx = 0;
    while bits + (odd[idx] as f64).log2() < (min_bits - 12) as f64 {
        base.push(odd[idx]);
        bits += (odd[idx] as f64).log2();
        idx += 1;
    }
    for &last in &odd[idx..] {
        let mut v = base.clone();
        v.push(last);
        if let Some(c) = Csurf::<FpM<N>>::new_big(&v) {
            return c;
        }
    }
    panic!("no CSURF prime");
}

#[test]
fn csurf_odd_steps_match_velu_and_phi() {
    let cs = small_csurf();
    let f = &cs.fp;
    let mut rng = Rng::new(10);
    let phis: Vec<Phi<Zp>> = [3usize, 5, 7, 11, 13]
        .iter()
        .map(|&l| Phi::compute(f, l))
        .collect();
    let mut a = f.zero();
    for round in 0..3 {
        // to_plus / from_plus round trip and isomorphism
        let m = cs.to_plus(a);
        assert_eq!(cs.from_plus(m), a, "from_plus(to_plus(A)) = A");
        assert!(fp_isomorphic(
            f,
            &to_weierstrass(f, m),
            &minus_to_weierstrass(f, a)
        ));
        for (idx, &ell) in cs.primes.iter().enumerate() {
            for positive in [true, false] {
                let (m2, k) = cs.odd_step_plus(m, ell, positive, &mut rng);
                // Weierstrass x-only Vélu on the same kernel
                let w = to_weierstrass(f, m);
                let xw = f.add(f.div(k.0, k.1), f.div(m, f.from_u64(3)));
                let vw = velu_xonly(f, &w, xw, ell).expect("velu_xonly");
                assert!(
                    fp_isomorphic(f, &vw.cod, &to_weierstrass(f, m2)),
                    "odd step = Weierstrass Vélu, l = {ell}"
                );
                // √élu on the same kernel
                let sv = SqrtVeluMont::new(f, m, k, ell).unwrap();
                assert_eq!(sv.codomain(f), m2, "√élu = Vélu, l = {ell}");
                // modular polynomial
                if idx < phis.len() && phis[idx].ell as u64 == ell {
                    let (j, j2) = (j_invariant(f, m), j_invariant(f, m2));
                    assert!(f.is_zero(phis[idx].eval(f, j, j2)), "Phi_{ell}");
                }
                // the codomain is on the surface with a unique Montgomery^- coefficient
                let a2 = cs.from_plus(m2);
                assert!(fp_isomorphic(
                    f,
                    &to_weierstrass(f, m2),
                    &minus_to_weierstrass(f, a2)
                ));
            }
        }
        let mut e = vec![0i32; cs.primes.len() + 1];
        e[0] = 1;
        e[1 + round] = 1;
        a = cs.action(a, &e, &mut rng);
    }
}

fn random_key(n: usize, bound: u64, rng: &mut Rng) -> Vec<i32> {
    (0..n)
        .map(|_| rng.below(2 * bound + 1) as i32 - bound as i32)
        .collect()
}

fn check_action<F: Field>(cs: &mut Csurf<F>, bound: u64, seed: u64, label: &str) {
    let f = cs.fp.clone();
    let mut rng = Rng::new(seed);
    let n = cs.primes.len() + 1;
    let ka = random_key(n, bound, &mut rng);
    let kb = random_key(n, bound, &mut rng);
    let e0 = f.zero();
    let t = Instant::now();
    let a = cs.action(e0, &ka, &mut rng);
    let t_step = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let a_b = cs.action_batched(e0, &ka, &mut rng);
    let t_batch = t.elapsed().as_secs_f64();
    assert_eq!(a, a_b, "{label}: stepwise = batched");
    let b = cs.action_batched(e0, &kb, &mut rng);
    let ab = cs.action_batched(a, &kb, &mut rng);
    let ba = cs.action_batched(b, &ka, &mut rng);
    assert_eq!(ab, ba, "{label}: actions commute");
    let neg: Vec<i32> = ka.iter().map(|&x| -x).collect();
    assert_eq!(
        cs.action_batched(a, &neg, &mut rng),
        e0,
        "{label}: inverse key returns"
    );
    // twisting inverts the action: [a] E^t = ([a]^{-1} E)^t, E_A^t = E_{-A}
    let mut rng2 = Rng::new(seed ^ 1);
    assert_eq!(
        cs.action_batched(f.neg(b), &ka, &mut rng2),
        f.neg(cs.action_batched(b, &neg, &mut rng2)),
        "{label}: twist"
    );
    // the same with √élu for every odd degree
    let old = cs.sqrt_velu_min;
    cs.sqrt_velu_min = 3;
    let t = Instant::now();
    let a_sv = cs.action_batched(e0, &ka, &mut rng);
    let t_sv = t.elapsed().as_secs_f64();
    cs.sqrt_velu_min = old;
    assert_eq!(a_sv, a, "{label}: √élu action = Vélu action");
    println!(
        "{label}: key bound {bound}, {n} exponents: stepwise {:.3} s, batched {:.3} s, batched+√élu {:.3} s",
        t_step, t_batch, t_sv
    );
}

#[test]
fn csurf_actions_commute_and_invert() {
    let mut cs = small_csurf();
    check_action(&mut cs, 3, 11, "CSURF 62-bit");
    let mut cs2 = find_big_csurf::<2>(118);
    println!(
        "CSURF ~128-bit prime: {} bits, primes {:?}",
        cs2.p1.bits(),
        cs2.primes
    );
    let phi2 = Phi::compute(&cs2.fp, 2);
    check_two_steps(&cs2, Some(&phi2), cs2.fp.zero(), 10, "128-bit");
    check_action(&mut cs2, 2, 12, "CSURF 128-bit");
}

#[test]
fn csurf512_prime_and_action() {
    let primes = csurf512_primes();
    assert_eq!(
        primes.iter().filter(|&&q| (5..=337).contains(&q)).count(),
        66
    );
    // the only reading of the paper's formula that is prime; neighbours are not
    for extra in [359u64, 397] {
        let mut v = primes.clone();
        v.push(extra);
        v.sort_unstable();
        assert!(Csurf::<FpM<9>>::new_big(&v).is_none(), "with {extra}");
    }
    let mut cs = Csurf::<FpM<9>>::csurf512();
    assert_eq!(cs.p1.bits(), 513);
    assert_eq!(cs.p1.sub_small(1).rem_small(8), 7);
    let f = cs.fp.clone();
    // l_2 steps radical vs Vélu, and a small key: commutation and inverse
    check_two_steps(&cs, None, f.zero(), 4, "CSURF-512");
    let mut rng = Rng::new(13);
    let n = cs.primes.len() + 1;
    let mut ka = vec![0i32; n];
    let mut kb = vec![0i32; n];
    ka[0] = 2; // l_2^2
    ka[2] = 1; // l_5
    ka[n - 1] = -1; // l_389^-1
    kb[1] = -1; // l_3^-1
    kb[n - 2] = 1; // l_383
    cs.sqrt_velu_min = 300;
    let t = Instant::now();
    let a = cs.action_batched(f.zero(), &ka, &mut rng);
    let ta = t.elapsed().as_secs_f64();
    let b = cs.action_batched(f.zero(), &kb, &mut rng);
    let ab = cs.action_batched(a, &kb, &mut rng);
    let ba = cs.action_batched(b, &ka, &mut rng);
    assert_eq!(ab, ba, "CSURF-512 commute");
    cs.sqrt_velu_min = u64::MAX;
    let neg: Vec<i32> = ka.iter().map(|&x| -x).collect();
    assert_eq!(cs.action(a, &neg, &mut rng), f.zero(), "CSURF-512 inverse");
    println!("CSURF-512: key with 4 isogenies (two of degree >= 383 by √élu) in {ta:.3} s");
}

#[test]
fn a5_remtree_sqrt_velu_large_degree() {
    // p = 4 * 3 * 1009 * 100003 * c - 1 (256 bits): kernels of degree 1009 and 100003
    let mut rng = Rng::new(14);
    let w = isogeny_algos::testdata::supersingular_workload::<4>(&[3, 1009, 100003], 256, &mut rng);
    let f = &w.f;
    let a = f.zero();
    let p1 = f.modulus().add_small(1);
    for &ell in &[1009u64, 100003] {
        let cof = p1.divrem_small(ell).0;
        let k = loop {
            let x = f.random(&mut rng);
            let k = ladder_xz(f, a24(f, a), (x, f.one()), &cof);
            if !f.is_zero(k.1) {
                break k;
            }
        };
        let t = Instant::now();
        let ms = multiples(f, a24(f, a), k, ((ell - 1) / 2) as usize);
        let a_ref = velu_codomain(f, a, &ms, ell);
        let t_velu = t.elapsed().as_secs_f64() * 1e3;
        let u = f.random(&mut rng);
        let x_ref = isog_x(f, &ms, u);
        let sv = SqrtVeluMont::new(f, a, k, ell).unwrap();
        let t = Instant::now();
        let rt = RemTreeSqrtVelu::new(f, &sv);
        let (a_rt, x_rt) = (rt.codomain(f, &sv), rt.eval(f, &sv, u));
        let t_rt = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!(a_rt, a_ref, "RemTree codomain, l = {ell}");
        assert_eq!(x_rt, x_ref, "RemTree x-map, l = {ell}");
        let mut sv_t = SqrtVeluMont::new(f, a, k, ell).unwrap();
        sv_t.force_tree = Some(true);
        let t = Instant::now();
        let (a_et, x_et) = (sv_t.codomain(f), sv_t.eval(f, u));
        let t_et = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!((a_et, x_et), (a_ref, x_ref));
        let t = Instant::now();
        let (a_pp, x_pp) = (sv.codomain(f), sv.eval(f, u));
        let t_pp = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!((a_pp, x_pp), (a_ref, x_ref));
        let t = Instant::now();
        let fused = fused_codomain_eval(f, &sv, u);
        let t_fu = t.elapsed().as_secs_f64() * 1e3;
        assert_eq!(fused, (a_ref, x_ref));
        println!(
            "l = {ell} (b = {}, b' = {}), codomain + 1 image: RemTree {t_rt:.1} ms, EvalTree {t_et:.1} ms, pair-by-pair {t_pp:.1} ms, fused pair-by-pair {t_fu:.1} ms; Vélu {t_velu:.1} ms",
            sv.b, sv.bp
        );
    }
}
