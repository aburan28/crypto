//! Tests for `fixed_degree::{factor, cornacchia, baselines}`.
use isogeny_algos::curve::*;
use isogeny_algos::field::{next_prime, Field, Rng, Zp, Zp2};
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::fixed_degree::baselines::*;
use isogeny_algos::fixed_degree::cornacchia::*;
use isogeny_algos::fixed_degree::factor::*;
use isogeny_algos::int::Int;
use isogeny_algos::path::graph::Path;
use isogeny_algos::quat::{klpt::ideal_from, lll_gram, Alg};
use std::collections::BTreeSet;

fn int(v: i64) -> Int {
    Int::from(v)
}

// ------------------------------------------------------------------ factor

#[test]
fn factor_random_products() {
    let mut rng = Rng::new(1);
    for t in 0..12 {
        let mut n = Int::one();
        let mut want: Vec<(Int, u32)> = vec![];
        let k = 2 + t % 3;
        for _ in 0..k {
            let bits = 12 + rng.below(20);
            let q = next_prime((1u64 << bits) + rng.below(1 << bits));
            let e = 1 + rng.below(2) as u32;
            n = &n * &Int::from(q).pow(e);
            want.push((Int::from(q), e));
        }
        n = &n * &int(1 << (t % 4));
        if t % 4 > 0 {
            want.push((int(2), t % 4));
        }
        want.sort();
        // merge equal primes (unlikely)
        let mut merged: Vec<(Int, u32)> = vec![];
        for (q, e) in want {
            match merged.last_mut() {
                Some((r, f)) if *r == q => *f += e,
                _ => merged.push((q, e)),
            }
        }
        let got = factor(&n);
        assert_eq!(got, merged, "n = {}", n);
        assert_eq!(product(&got), n);
        assert_eq!(omega(&n), merged.len());
    }
    // a bounded run that must give up: product of two ~40-bit primes with a tiny budget
    let a = Int::from(next_prime(1 << 40));
    let b = Int::from(next_prime((1 << 41) + 12345));
    let c = Int::from(next_prime((1 << 42) + 999));
    let n = &(&a * &b) * &c;
    assert!(factor_bounded(&n, 100).is_none());
    assert_eq!(factor(&n).len(), 3);
}

#[test]
fn sqrt_mod_matches_brute_force() {
    for n in 1u64..=200 {
        let fac: Vec<(Int, u32)> = if n == 1 {
            vec![]
        } else {
            factor(&Int::from(n))
        };
        for a in 0..n {
            let brute: Vec<Int> = (0..n).filter(|x| x * x % n == a).map(Int::from).collect();
            let got = sqrt_mod_composite(&Int::from(a), &Int::from(n), &fac);
            assert_eq!(got, brute, "a = {a}, n = {n}");
        }
    }
    // prime powers incl. 2^k and p | a
    for (p, k) in [(2u64, 9u32), (3, 6), (5, 4), (7, 3), (11, 2)] {
        let pk = p.pow(k);
        for a in 0..pk {
            let brute: Vec<Int> = (0..pk).filter(|x| x * x % pk == a).map(Int::from).collect();
            assert_eq!(sqrt_mod_prime_power(&Int::from(a), &Int::from(p), k), brute);
        }
    }
}

// ------------------------------------------------------------------ Cornacchia

#[test]
fn cornacchia_all_matches_brute_force() {
    let mut total = 0usize;
    for d in [1i64, 2, 3, 4, 5, 7, 11, 12, 15, 20, 27, 163] {
        for m in 1i64..=2500 {
            let mut brute = vec![];
            let mut y = 0i64;
            while d * y * y <= m {
                let r = m - d * y * y;
                let x = (r as f64).sqrt() as i64;
                for xx in [x - 1, x, x + 1] {
                    if xx >= 0 && xx * xx == r {
                        brute.push((int(xx), int(y)));
                    }
                }
                y += 1;
            }
            brute.sort();
            brute.dedup();
            let got = cornacchia_all(&int(d), &int(m), &factor(&int(m)));
            assert_eq!(got, brute, "d = {d}, m = {m}");
            total += got.len();
        }
    }
    eprintln!("cornacchia_all: {total} solutions checked against brute force");
}

#[test]
fn binary_quadratic_matches_brute_force() {
    let mut rng = Rng::new(7);
    let mut nsol = 0;
    let mut cases = 0;
    while cases < 300 {
        let r = |rng: &mut Rng, k: i64| rng.below(2 * k as u64 + 1) as i64 - k;
        let a = 1 + rng.below(6) as i64;
        let c = 1 + rng.below(6) as i64;
        let b = r(&mut rng, 6);
        if b * b - 4 * a * c >= 0 {
            continue;
        }
        let (d, e) = (r(&mut rng, 25), r(&mut rng, 25));
        let f = -(rng.below(600) as i64);
        cases += 1;
        let lam = ((a + c) as f64 - (((a - c) * (a - c) + b * b) as f64).sqrt()) / 2.0;
        let s = (d.abs() + e.abs()) as f64;
        let rr = ((s + (s * s + 4.0 * lam * (-f) as f64).sqrt()) / (2.0 * lam)).ceil() as i64 + 2;
        let mut brute = BTreeSet::new();
        for x in -rr..=rr {
            for y in -rr..=rr {
                if a * x * x + b * x * y + c * y * y + d * x + e * y + f == 0 {
                    brute.insert((int(x), int(y)));
                }
            }
        }
        match solve_binary_quadratic(&int(a), &int(b), &int(c), &int(d), &int(e), &int(f)) {
            BqOutcome::Solved { sols, .. } => {
                let got: BTreeSet<_> = sols.into_iter().collect();
                assert_eq!(got, brute, "{a} {b} {c} {d} {e} {f}");
                nsol += got.len();
            }
            BqOutcome::Abandoned(r) => panic!("abandoned small case: {r:?}"),
        }
    }
    eprintln!("binary quadratic: {cases} random equations, {nsol} solutions, all = brute force");
    // the omega heuristic: N = 2 * 3 * 5 * ... * 47 (15 primes) vs B ln ln N with B = 11 -> fine,
    // B = 1 -> abandoned. x^2 + y^2 = N: a=c=1, f=-N; D = -4, E = 0, N' = 16 N.
    let n: Int = [2i64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        .iter()
        .fold(Int::one(), |acc, &q| &acc * &int(q));
    let coef = [int(1), int(0), int(1), int(0), int(0), -n.clone()];
    let strict = BqParams {
        b: 1.0,
        ..Default::default()
    };
    assert!(matches!(
        solve_binary_quadratic_with(&coef, &strict),
        BqOutcome::Abandoned(AbandonReason::TooManyPrimeFactors { .. })
    ));
}

fn chol_inv_frob2(g: &[[Int; 4]; 4]) -> f64 {
    // G = L L^T; |x|^2 <= ||L^{-1}||_F^2 Q(x)
    let a: Vec<Vec<f64>> = g
        .iter()
        .map(|r| r.iter().map(|v| v.to_f64()).collect())
        .collect();
    let mut l = vec![vec![0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            let mut s = a[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            l[i][j] = if i == j { s.sqrt() } else { s / l[j][j] };
        }
    }
    let mut inv = vec![vec![0f64; 4]; 4];
    for c in 0..4 {
        for i in 0..4 {
            let mut s = if i == c { 1.0 } else { 0.0 };
            for k in 0..i {
                s -= l[i][k] * inv[k][c];
            }
            inv[i][c] = s / l[i][i];
        }
    }
    inv.iter().flatten().map(|v| v * v).sum()
}

fn transform(h: &[[Int; 4]; 4], g: &[[Int; 4]; 4]) -> [[Int; 4]; 4] {
    std::array::from_fn(|r| {
        std::array::from_fn(|s| {
            let mut acc = Int::zero();
            for i in 0..4 {
                for j in 0..4 {
                    acc = &acc + &(&h[r][i] * &(&g[i][j] * &h[s][j]));
                }
            }
            acc
        })
    })
}

/// Random positive-definite Gram: U^T diag U for a random unimodular U, then LLL.
fn synthetic_gram(rng: &mut Rng, diag: [i64; 4], mix_steps: usize) -> [[Int; 4]; 4] {
    let mut u: [[Int; 4]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| int((i == j) as i64)));
    for _ in 0..mix_steps {
        let i = rng.below(4) as usize;
        let mut j = rng.below(4) as usize;
        if j == i {
            j = (i + 1) % 4;
        }
        let c = int(rng.below(7) as i64 - 3);
        for t in 0..4 {
            u[i][t] = &u[i][t] + &(&c * &u[j][t]);
        }
    }
    let dg: [[Int; 4]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| int(if i == j { diag[i] } else { 0 })));
    let g = transform(&u, &dg);
    transform(&lll_gram(&g), &g)
}

#[test]
fn norm_equation_brute_force_per_guess() {
    let mut rng = Rng::new(11);
    let mut checked = 0;
    for inst in 0..6 {
        let g = synthetic_gram(&mut rng, [1 + inst % 3, 2, 3 + inst, 5], 12);
        let x0: [Int; 4] = std::array::from_fn(|_| int(rng.below(9) as i64 - 4));
        let d = eval_form(&g, &x0);
        if d.is_zero() {
            continue;
        }
        let r = ((chol_inv_frob2(&g) * d.to_f64()).sqrt()).ceil() as i64 + 1;
        let bound = 2i64;
        for (x3, x4) in guess_order(bound) {
            let mut brute = BTreeSet::new();
            for x1 in -r..=r {
                for x2 in -r..=r {
                    let x = [int(x1), int(x2), int(x3), int(x4)];
                    if eval_form(&g, &x) == d {
                        brute.insert((int(x1), int(x2)));
                    }
                }
            }
            let coef = guess_quadratic(&g, &d, &int(x3), &int(x4));
            let BqOutcome::Solved { sols, .. } =
                solve_binary_quadratic_with(&coef, &BqParams::default())
            else {
                panic!("abandoned");
            };
            let got: BTreeSet<_> = sols.into_iter().collect();
            assert_eq!(got, brute, "instance {inst} guess ({x3}, {x4})");
            checked += 1;
        }
        // and the driver finds a solution whenever the planted one is inside the box
        let planted_in_box = x0[2].abs() <= int(bound) && x0[3].abs() <= int(bound);
        match solve_norm_eq_cornacchia(&g, &d, &int(bound)) {
            CornacchiaOutcome::Found(x, _) => assert_eq!(eval_form(&g, &x), d),
            CornacchiaOutcome::NotFound(_) => assert!(!planted_in_box),
        }
    }
    eprintln!("norm equation: {checked} (x3, x4) guesses agree with brute force over (x1, x2)");
}

#[test]
fn norm_equation_planted_synthetic_and_quaternion() {
    let mut rng = Rng::new(21);
    // synthetic, larger entries
    let mut tot = CornacchiaStats::default();
    for _ in 0..10 {
        let diag = [
            1 + rng.below(1000) as i64,
            1 + rng.below(100_000) as i64,
            1 + rng.below(10_000_000) as i64,
            1 + rng.below(1_000_000_000) as i64,
        ];
        let g = synthetic_gram(&mut rng, diag, 20);
        let mut x0: [Int; 4] = std::array::from_fn(|_| int(rng.below(2001) as i64 - 1000));
        x0[2] = int(rng.below(7) as i64 - 3);
        x0[3] = int(rng.below(7) as i64 - 3);
        let d = eval_form(&g, &x0);
        match solve_norm_eq_cornacchia(&g, &d, &int(3)) {
            CornacchiaOutcome::Found(x, st) => {
                assert_eq!(eval_form(&g, &x), d);
                tot.guesses += st.guesses;
                tot.factorisations += st.factorisations;
                tot.abandoned += st.abandoned;
            }
            CornacchiaOutcome::NotFound(st) => panic!("planted solution missed: {st:?}"),
        }
    }
    eprintln!("synthetic (planted, 10 instances): {tot:?}");
    // Gram matrices of random left O_0-ideals of prime norm, LLL-reduced
    for p in [10007i64, 1_000_003] {
        let alg = Alg::new(&int(p));
        let o0 = alg.o0();
        let mut tot = CornacchiaStats::default();
        for _ in 0..5 {
            let n = Int::from(next_prime(100 + rng.below(900)));
            let a = isogeny_algos::quat::Quat::from_ints([
                rng.below(1000) as i64,
                rng.below(1000) as i64,
                rng.below(1000) as i64,
                rng.below(1000) as i64,
            ]);
            let lat = ideal_from(&alg, &o0, &n, &a);
            let g0 = lat.gram(&alg);
            let g = transform(&lll_gram(&g0), &g0);
            let x0 = [
                int(rng.below(41) as i64 - 20),
                int(rng.below(41) as i64 - 20),
                int(rng.below(5) as i64 - 2),
                int(rng.below(5) as i64 - 2),
            ];
            let d = eval_form(&g, &x0);
            if d.is_zero() {
                continue;
            }
            match solve_norm_eq_cornacchia(&g, &d, &int(2)) {
                CornacchiaOutcome::Found(x, st) => {
                    assert_eq!(eval_form(&g, &x), d);
                    tot.guesses += st.guesses;
                    tot.factorisations += st.factorisations;
                    tot.abandoned += st.abandoned;
                }
                CornacchiaOutcome::NotFound(st) => panic!("planted solution missed: {st:?}"),
            }
        }
        eprintln!("quaternion ideals p = {p}: {tot:?}");
    }
}

// ------------------------------------------------------------------ baselines

fn prime_3mod4_near(bits: u32, extra: u64) -> u64 {
    let mut p = next_prime((1u64 << bits) + extra);
    while p % 4 != 3 {
        p = next_prime(p + 1);
    }
    p
}

fn phi2(p: u64, ell: usize) -> (Zp2, Phi<Zp2>) {
    let fp = Zp::new(p);
    let f2 = Zp2::new(p);
    let phi = Phi::compute(&fp, ell).lift(|x| (x, 0u64));
    (f2, phi)
}

/// Non-backtracking random j-walk of `steps` l-steps.
fn walk(f2: &Zp2, phi: &Phi<Zp2>, j: (u64, u64), steps: usize, rng: &mut Rng) -> Path<(u64, u64)> {
    let mut js = vec![j];
    let mut prev = None;
    for _ in 0..steps {
        let cur = *js.last().unwrap();
        let mut ns = phi.neighbors(f2, cur, rng);
        ns.retain(|&n| Some(n) != prev);
        let nx = ns[rng.below(ns.len() as u64) as usize];
        prev = Some(cur);
        js.push(nx);
    }
    Path {
        ells: vec![phi.ell; steps],
        js,
    }
}

#[test]
fn exhaustive_prime_degree_finds_planted() {
    let mut rng = Rng::new(31);
    // p + 1 divisible by 4 * 3 * 5 * 7 * 11 * 13
    let base = 4u64 * 3 * 5 * 7 * 11 * 13;
    let mut c = 3u64;
    let p = loop {
        let p = base * c - 1;
        if isogeny_algos::field::is_prime(p) {
            break p;
        }
        c += 1;
    };
    let f2 = Zp2::new(p);
    let mut e = Curve::new(f2.one(), f2.zero());
    for _ in 0..4 {
        let k = isogeny_algos::testdata::ss_point_of_order(&f2, &e, p, 3, &mut rng);
        e = isogeny_algos::kernel::velu::velu_from_point(&f2, &e, &k, 3).cod;
    }
    let fp = Zp::new(p);
    for d in [5u64, 7, 11, 13] {
        let k = isogeny_algos::testdata::ss_point_of_order(&f2, &e, p, d, &mut rng);
        let e2 = isogeny_algos::kernel::velu::velu_from_point(&f2, &e, &k, d).cod;
        let (j1, j2) = (jinv(&f2, &e), jinv(&f2, &e2));
        let phi = Phi::compute(&fp, d as usize).lift(|x| (x, 0u64));
        let r = exhaustive_prime_degree_phi(&f2, &phi, j1, j2);
        assert!(r.match_index.is_some(), "Phi_{d} search missed");
        let m = exhaustive_prime_degree_kernel(&f2, &e, p, d, j2, &mut rng)
            .expect("kernel search missed");
        assert_eq!(jinv(&f2, &m.codomain), j2);
        assert_eq!(pmul(&f2, &e, &m.kernel_gen, d as u128), Pt::Inf);
        eprintln!(
            "exhaustive d = {d} (p = {p}): {} distinct Phi roots, kernel match after {} / {} subgroups",
            r.candidates,
            m.tried,
            d + 1
        );
    }
}

#[test]
fn mitm_and_vow_find_planted_paths() {
    let mut rng = Rng::new(41);
    // ~2^22 / 12 ~ 350k supersingular j's, 3 * 2^13 ~ 24.6k cyclic 2^14-isogenies from j1:
    // the planted claw is unique with probability ~ 0.93.
    let p = prime_3mod4_near(22, 12345);
    let (f2, phi) = phi2(p, 2);
    let e = 14usize;
    let mut summary = vec![];
    for inst in 0..3 {
        let j0 = (1728u64, 0u64);
        let j1 = *walk(&f2, &phi, j0, 40 + inst, &mut rng).js.last().unwrap();
        let planted = walk(&f2, &phi, j1, e, &mut rng);
        let j2 = *planted.js.last().unwrap();
        assert!(verify_cyclic_path(&f2, &phi, &planted));

        let (m, ms) = mitm_smooth(&f2, &phi, e, j1, j2);
        let m = m.expect("MITM failed");
        assert_eq!(m.js[0], j1);
        assert_eq!(*m.js.last().unwrap(), j2);
        assert_eq!(m.len(), e);
        assert!(verify_cyclic_path(&f2, &phi, &m));

        let prm = VowParams {
            seed: 100 + inst as u64,
            ..VowParams::new(16)
        };
        let (v, vs) = vow_claw(&f2, &phi, e, j1, j2, &prm);
        let v = v.expect("vOW failed");
        assert_eq!(v.js[0], j1);
        assert_eq!(*v.js.last().unwrap(), j2);
        assert_eq!(v.len(), e);
        assert!(verify_cyclic_path(&f2, &phi, &v));
        let theory = vow_expected_evals(vs.space, prm.memory);
        eprintln!(
            "inst {inst}: MITM expansions {} leaves {}+{} (theory 2*(l+1)l^(e/2-1) = {}); vOW f-evals {} \
             (theory 2.5 sqrt(|S|^3/w) = {:.0}, |S| = {}, w = {}), dist {}, collisions {} (useless {}), versions {}",
            ms.expansions,
            ms.leaves1,
            ms.leaves2,
            2 * 3 * (1 << (e / 2 - 1)),
            vs.f_evals,
            theory,
            vs.space,
            prm.memory,
            vs.distinguished,
            vs.collisions,
            vs.useless,
            vs.versions
        );
        summary.push((vs.f_evals as f64) / theory);
    }
    eprintln!("vOW measured/theory ratios: {summary:?}");
    let q = grover_tani_cost_model(2, e as u32);
    eprintln!(
        "quantum cost model l = 2, e = {e}: Grover {:.1}, Tani {:.1} queries",
        q.grover_queries, q.tani_queries
    );
    assert!((q.tani_queries - (192f64 * 192.0).cbrt()).abs() < 1e-6);
}

#[test]
fn mitm_odd_exponent_ell3() {
    let mut rng = Rng::new(51);
    let p = prime_3mod4_near(20, 777);
    let (f2, phi) = phi2(p, 3);
    let e = 7usize;
    let j1 = *walk(&f2, &phi, (1728, 0), 30, &mut rng).js.last().unwrap();
    let planted = walk(&f2, &phi, j1, e, &mut rng);
    let j2 = *planted.js.last().unwrap();
    let (m, ms) = mitm_smooth(&f2, &phi, e, j1, j2);
    let m = m.expect("MITM failed");
    assert!(verify_cyclic_path(&f2, &phi, &m) && m.len() == e && *m.js.last().unwrap() == j2);
    let (v, vs) = vow_claw(&f2, &phi, e, j1, j2, &VowParams::new(8));
    let v = v.expect("vOW failed");
    assert!(verify_cyclic_path(&f2, &phi, &v) && v.len() == e && *v.js.last().unwrap() == j2);
    eprintln!(
        "l = 3, e = 7: MITM leaves {}+{} (theory 4*3^3 = 108, 4*3^2 = 36); vOW f-evals {} (theory {:.0})",
        ms.leaves1,
        ms.leaves2,
        vs.f_evals,
        vow_expected_evals(vs.space, 8)
    );
}

/// vOW scaling: e = 18 (|S| = 2 * 3 * 2^8 = 1536), memory w in {16, 64}, mean over instances.
/// p ~ 2^28 (~2.2e7 supersingular j's vs 3 * 2^17 cyclic 2^18-isogenies): the planted claw is
/// unique with probability ~ 0.98, so vOW is looking for a single golden collision.
#[test]
fn vow_scaling_vs_theory() {
    let mut rng = Rng::new(61);
    let p = prime_3mod4_near(28, 4321);
    let (f2, phi) = phi2(p, 2);
    let e = 18usize;
    let inst = 4;
    let mut targets = vec![];
    for _ in 0..inst {
        let j1 = *walk(&f2, &phi, (1728, 0), 50, &mut rng).js.last().unwrap();
        let j2 = *walk(&f2, &phi, j1, e, &mut rng).js.last().unwrap();
        targets.push((j1, j2));
    }
    let mut mitm_leaves = 0;
    for &(j1, j2) in &targets {
        let (m, ms) = mitm_smooth(&f2, &phi, e, j1, j2);
        assert!(verify_cyclic_path(&f2, &phi, &m.expect("MITM failed")));
        mitm_leaves = ms.leaves1 + ms.leaves2;
    }
    eprintln!("MITM e = {e}: {mitm_leaves} leaves per instance (theory 2 * 3 * 2^8 = 1536)");
    for w in [16usize, 64] {
        let mut sum = 0u64;
        let mut space = 0;
        for (i, &(j1, j2)) in targets.iter().enumerate() {
            let prm = VowParams {
                seed: 1000 + i as u64,
                ..VowParams::new(w)
            };
            let (v, vs) = vow_claw(&f2, &phi, e, j1, j2, &prm);
            let v = v.expect("vOW failed");
            assert!(verify_cyclic_path(&f2, &phi, &v) && v.len() == e);
            assert_eq!((v.js[0], *v.js.last().unwrap()), (j1, j2));
            sum += vs.f_evals;
            space = vs.space;
        }
        let mean = sum as f64 / inst as f64;
        let th = vow_expected_evals(space, w);
        eprintln!(
            "vOW e = {e}, w = {w}: mean f-evals {mean:.0} over {inst} instances, theory {th:.0}, ratio {:.2}",
            mean / th
        );
    }
}
