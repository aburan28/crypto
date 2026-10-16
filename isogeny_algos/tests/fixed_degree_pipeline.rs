//! End-to-end tests of the fixed-degree isogeny algorithm (Benčina–Kutas–Merz–Petit–Stopar–
//! Weitkämper): End(E) from a path, Kirschmer–Voight connecting ideals, norm forms, the
//! norm-equation solvers inside the pipeline, ideal -> isogeny (kernel and torsion/Kani
//! representations) and the hybrid method. Every check compares against an independent
//! computation: Phi_l neighbours on the j-line, Deuring images, Vélu codomains, Kani splits with
//! a twisted negative control, and Weil pairings.
use isogeny_algos::curve::{from_j, jinv, Curve};
use isogeny_algos::ext::Fp4;
use isogeny_algos::field::Rng;
use isogeny_algos::fixed_degree::connecting::{
    combine, connecting_ideal, connecting_via_o0, connects, norm_form,
};
use isogeny_algos::fixed_degree::endring::{end_ring, random_ss_j, sub_ideals, JGraph, J};
use isogeny_algos::fixed_degree::hybrid::hybrid_search;
use isogeny_algos::fixed_degree::pipeline::{
    find_isogeny_of_degree, plant_ideal, Instance, Solver,
};
use isogeny_algos::fixed_degree::represent::{
    fast_ideal_to_j, ideal_to_isogeny, is_deuring_smooth, two_adic, Repr, ReprOptions, E4,
};
use isogeny_algos::int::Int;
use isogeny_algos::quat::deuring::Deuring;

/// p = 2^15 * 49 - 1: v_2(p + 1) = 15 (torsion level a = 12, 2^a > sqrt(p) ~ 1267), and
/// p^2 - 1 has a large 256-smooth odd part (searched for; see the report).
const P: u64 = 1605631;

fn setup(seed: u64) -> (Deuring, JGraph, Rng) {
    setup_p(P, seed)
}

fn setup_p(p: u64, seed: u64) -> (Deuring, JGraph, Rng) {
    let mut rng = Rng::new(seed);
    let dg = Deuring::new(p, &mut rng);
    let g = JGraph::new(p, 2);
    (dg, g, rng)
}

fn model(f: &Fp4, j: J) -> Curve<E4> {
    from_j(f, f.embed(j))
}

fn is_prime(n: u64) -> bool {
    isogeny_algos::field::is_prime(n)
}

/// The l + 1 sub-ideals of an ideal are exactly the l-isogenous curves (Phi_l roots).
#[test]
fn sub_ideals_are_the_phi_neighbours() {
    let (dg, _, mut rng) = setup(1);
    let (alg, o0) = (&dg.alg, &dg.o0);
    for ell in [2u64, 3] {
        let mut g = JGraph::new(P, ell as usize);
        let mut ideal = o0.clone();
        for step in 0..3 {
            let j = fast_ideal_to_j(&dg, &ideal).unwrap();
            let subs = sub_ideals(alg, o0, &ideal, ell);
            assert_eq!(subs.len(), ell as usize + 1, "ell = {ell}, step {step}");
            let mut js: Vec<J> = subs
                .iter()
                .map(|s| fast_ideal_to_j(&dg, s).unwrap())
                .collect();
            js.sort();
            js.dedup();
            assert_eq!(js, g.neighbours(j), "ell = {ell}, step {step}");
            ideal = subs[rng.below(subs.len() as u64) as usize].clone();
        }
    }
}

/// Step 1: End(E) of random supersingular curves (independent random walks on the j-line).
#[test]
fn end_ring_of_random_curves() {
    let (dg, mut g, mut rng) = setup(2);
    let o0 = &dg.o0;
    for t in 0..4 {
        let j = random_ss_j(&mut g, P, 25 + 5 * t, &mut rng);
        let er = end_ring(&dg, &mut g, j).expect("end ring");
        assert_eq!(fast_ideal_to_j(&dg, &er.ideal), Some(j));
        // the path is a Phi_2 path
        for w in er.path.windows(2) {
            assert!(g.neighbours(w[0]).contains(&w[1]));
        }
        // a maximal order: same discriminant (covolume) as O_0
        let or = er.ideal.right_order(&dg.alg, o0);
        assert_eq!(or.det(), o0.det());
        eprintln!(
            "End(E): path length {}, {} graph vertices, {} Deuring evaluations",
            er.path.len() - 1,
            er.graph_vertices,
            er.deuring_evals
        );
    }
}

/// Steps 2–3: Kirschmer–Voight connecting ideal and its norm form.
#[test]
fn connecting_ideal_and_norm_form() {
    let (dg, mut g, mut rng) = setup(3);
    let (alg, o0) = (&dg.alg, &dg.o0);
    for _ in 0..3 {
        let j1 = random_ss_j(&mut g, P, 20, &mut rng);
        let j2 = random_ss_j(&mut g, P, 20, &mut rng);
        let i1 = end_ring(&dg, &mut g, j1).unwrap().ideal;
        let i2 = end_ring(&dg, &mut g, j2).unwrap().ideal;
        let (o1, o2) = (i1.right_order(alg, o0), i2.right_order(alg, o0));
        let kv = connecting_ideal(alg, &o1, &o2);
        assert!(connects(alg, &kv, &o1, &o2));
        let viao0 = connecting_via_o0(alg, &i1, &i2);
        assert!(connects(alg, &viao0, &o1, &o2));
        // the Deuring image of I_1 * (connecting ideal) is E_2
        assert_eq!(dg.ideal_to_j(&i1.mul(alg, &kv), &mut rng), Some(j2));
        let nf = norm_form(alg, &o1, &kv);
        for _ in 0..20 {
            let x: [Int; 4] = std::array::from_fn(|_| Int::from(rng.below(41) as i64 - 20));
            let a = combine(&nf.basis, &x);
            let mut q = Int::zero();
            for r in 0..4 {
                for s in 0..4 {
                    q = &q + &(&(&nf.gram[r][s] * &x[r]) * &x[s]);
                }
            }
            let (n, d) = alg.nrd(&a);
            assert_eq!(d, Int::one());
            assert_eq!(&q * &nf.norm, &n * &Int::from(2i64));
        }
    }
}

fn planted(dg: &Deuring, g: &mut JGraph, rng: &mut Rng, primes: &[u64]) -> (J, J, Int) {
    let j1 = random_ss_j(g, dg.p, 30, rng);
    let i1 = end_ring(dg, g, j1).unwrap().ideal;
    let qs: Vec<Int> = primes.iter().map(|&q| Int::from(q)).collect();
    let jd = plant_ideal(&dg.alg, &dg.o0, &i1, &qs, rng);
    let j2 = fast_ideal_to_j(dg, &i1.mul(&dg.alg, &jd)).unwrap();
    let d = qs.iter().fold(Int::one(), |a, q| &a * q);
    (j1, j2, d)
}

/// Steps 1–6 with Cornacchia, non-smooth prime d with sqrt(p) < d < 2^a: torsion representation
/// certified by Kani (and the twisted control), Weil pairing, codomain j.
#[test]
fn pipeline_cornacchia_prime_degree_kani() {
    let (dg, mut g, mut rng) = setup(4);
    let f = dg.f2.clone();
    let a = two_adic(P) - 3;
    let sqrt_p = (P as f64).sqrt() as u64;
    let primes: Vec<u64> = (sqrt_p..(1 << a)).filter(|&q| is_prime(q)).collect();
    let mut done = 0;
    for t in 0..3 {
        let q = primes[rng.below(primes.len() as u64) as usize];
        assert!(!is_deuring_smooth(&dg, &Int::from(q)));
        let (j1, j2, d) = planted(&dg, &mut g, &mut rng, &[q]);
        let inst = Instance {
            j1,
            j2,
            e1: Some(model(&f, j1)),
            e2: Some(model(&f, j2)),
            d: d.clone(),
        };
        let sol = find_isogeny_of_degree(
            &dg,
            &mut g,
            &inst,
            Solver::Cornacchia,
            &ReprOptions::default(),
            &mut rng,
        )
        .expect("pipeline");
        let rep = &sol.isogeny;
        assert!(rep.j_ok);
        assert_eq!(jinv(&f, &rep.codomain), f.embed(j2));
        match &rep.repr {
            Repr::Torsion { certificate, .. } => {
                assert!(certificate.certified(), "Kani certificate: {certificate:?}")
            }
            _ => panic!("expected the torsion representation"),
        }
        assert_eq!(rep.pairing_ok, Some(true));
        eprintln!(
            "trial {t}: d = {d} (log_p d = {:.3}), N(conn) = {}, x = {:?}, times {:?}",
            (q as f64).ln() / (P as f64).ln(),
            sol.connecting_norm,
            sol.degree.x,
            sol.times
        );
        done += 1;
    }
    assert_eq!(done, 3);
}

/// Smooth d (every prime power has rational torsion over F_{p^4}, so `Deuring::kernel_of_ideal`
/// applies): kernel representation; the Vélu chain from the caller's E_1 with the transported
/// kernel lands on a curve with j(E_2).
#[test]
fn pipeline_smooth_degree_kernel() {
    // p = 2^10 * 1443 - 1: p + 1 = 2^10 3 13 37, p - 1 = 2 5 7 11 19 101 (t_odd ~ p^1.47)
    let (dg, mut g, mut rng) = setup_p(1477631, 5);
    let f = dg.f2.clone();
    // planted smooth degrees above sqrt(p) ~ 1215 that leave t_odd room for N(I_1')
    for primes in [vec![3u64, 13, 37], vec![5, 7, 37]] {
        let (j1, j2, d) = planted(&dg, &mut g, &mut rng, &primes);
        assert!(is_deuring_smooth(&dg, &d));
        let inst = Instance {
            j1,
            j2,
            e1: Some(model(&f, j1)),
            e2: Some(model(&f, j2)),
            d: d.clone(),
        };
        let sol = find_isogeny_of_degree(
            &dg,
            &mut g,
            &inst,
            Solver::Cornacchia,
            &ReprOptions::default(),
            &mut rng,
        )
        .expect("pipeline");
        match &sol.isogeny.repr {
            Repr::Kernel { gens } => {
                let order: u64 = gens.iter().map(|g| g.0.pow(g.1)).product();
                assert_eq!(Int::from(order), d, "kernel order");
            }
            _ => panic!("expected the kernel representation"),
        }
        assert!(sol.isogeny.j_ok, "Vélu codomain j");
        eprintln!(
            "d = {d} = {primes:?}: x = {:?}, times {:?}",
            sol.degree.x, sol.times
        );
    }
}

/// Random E_1, E_2 (independent random walks, nothing planted) and a prime d ~ 8p > p, so that
/// degree-d isogenies E_1 -> E_2 exist generically (about 12 (d + 1) / p of them): the pipeline
/// finds one. d > 2^a, so the torsion representation carries no Kani certificate; it is checked
/// by the Deuring codomain j and the Weil pairing.
#[test]
fn pipeline_random_target_large_degree() {
    let (dg, mut g, mut rng) = setup(8);
    let f = dg.f2.clone();
    let mut ok = 0;
    for t in 0..2 {
        let j1 = random_ss_j(&mut g, P, 30, &mut rng);
        let j2 = random_ss_j(&mut g, P, 31, &mut rng);
        let mut d = 8 * P + rng.below(P);
        while !is_prime(d) {
            d += 1;
        }
        let inst = Instance {
            j1,
            j2,
            e1: Some(model(&f, j1)),
            e2: Some(model(&f, j2)),
            d: Int::from(d),
        };
        match find_isogeny_of_degree(
            &dg,
            &mut g,
            &inst,
            Solver::Cornacchia,
            &ReprOptions::default(),
            &mut rng,
        ) {
            Ok(sol) => {
                assert!(sol.isogeny.j_ok, "Deuring codomain j");
                assert_eq!(sol.isogeny.pairing_ok, Some(true));
                assert_eq!(jinv(&f, &sol.isogeny.codomain), f.embed(j2));
                eprintln!(
                    "trial {t}: d = {d} (log_p d = {:.3}), x = {:?}, guesses {:?}, times {:?}",
                    (d as f64).ln() / (P as f64).ln(),
                    sol.degree.x,
                    sol.degree.cornacchia_stats.as_ref().map(|s| s.guesses),
                    sol.times
                );
                ok += 1;
            }
            Err(e) => eprintln!("trial {t}: d = {d}: {e}"),
        }
    }
    assert!(ok >= 1);
}

/// Coron's bivariate Coppersmith method inside the full pipeline (planted d).
#[test]
fn pipeline_coppersmith_bivariate() {
    let (dg, mut g, mut rng) = setup(6);
    let f = dg.f2.clone();
    let sqrt_p = (P as f64).sqrt() as u64;
    let mut q = sqrt_p + rng.below(sqrt_p);
    while !is_prime(q) {
        q += 1;
    }
    let (j1, j2, d) = planted(&dg, &mut g, &mut rng, &[q]);
    let inst = Instance {
        j1,
        j2,
        e1: Some(model(&f, j1)),
        e2: Some(model(&f, j2)),
        d,
    };
    let sol = find_isogeny_of_degree(
        &dg,
        &mut g,
        &inst,
        Solver::CoppersmithBivariate,
        &ReprOptions::default(),
        &mut rng,
    )
    .expect("pipeline");
    assert!(sol.isogeny.j_ok);
    assert_eq!(sol.isogeny.pairing_ok, Some(true));
    match &sol.isogeny.repr {
        Repr::Torsion { certificate, .. } => assert!(certificate.certified(), "{certificate:?}"),
        _ => panic!("expected the torsion representation"),
    }
    eprintln!(
        "Coron: d = {q}, x = {:?}, times {:?}",
        sol.degree.x, sol.times
    );
}

/// Hybrid: d = 3^7 = 2187 (sqrt(p) < d < 2^a), guess a 3^2 prefix, solve the 3^5 rest,
/// compose; the composite's Deuring image is E_2 and its isogeny is Kani-certified.
#[test]
fn hybrid_three_power_degree() {
    let (dg, mut g, mut rng) = setup(7);
    let (alg, o0) = (&dg.alg, &dg.o0);
    let e = 7u32;
    let (j1, j2, _) = planted(&dg, &mut g, &mut rng, &vec![3; e as usize]);
    let i1 = end_ring(&dg, &mut g, j1).unwrap().ideal;
    let i2 = end_ring(&dg, &mut g, j2).unwrap().ideal;
    let res = hybrid_search(alg, o0, &i1, &i2, 3, e, 2, Solver::Cornacchia).expect("hybrid");
    assert_eq!(
        fast_ideal_to_j(&dg, &i1.mul(alg, &res.ideal)),
        Some(j2),
        "hybrid ideal lands on E_2"
    );
    eprintln!("hybrid: {} prefixes tried", res.guesses);
    let rep = ideal_to_isogeny(
        &dg,
        &i1,
        &res.ideal,
        None,
        &ReprOptions::default(),
        &mut rng,
    )
    .expect("isogeny");
    assert!(rep.j_ok);
    assert_eq!(rep.pairing_ok, Some(true));
    match &rep.repr {
        Repr::Torsion { certificate, .. } => assert!(certificate.certified(), "{certificate:?}"),
        _ => panic!("expected the torsion representation"),
    }
}

/// Steps 2–5 at cryptographic size (no Deuring step, so p is not limited): O_1 = O_R(I_1) for a
/// random ideal I_1 of O_0, O_2 = O_R(I_2) with I_2 a reduced ideal equivalent to I_1 J_0 for a
/// random J_0 of prime norm d (so the Kirschmer–Voight ideal does not reveal d); every solver
/// must return alpha with Nrd(alpha) = d N(I), and J = I conj(alpha) / N(I) must have norm d,
/// left order O_1 and right order beta^{-1} O_2 beta (beta = conj(alpha) / N(I)).
#[test]
fn steps_two_to_five_cryptographic_size() {
    use isogeny_algos::fixed_degree::endring::reduce_ideal;
    use isogeny_algos::fixed_degree::pipeline::{random_prime_norm_ideal, solve_degree};
    use isogeny_algos::quat::Alg;
    let mut rng = Rng::new(9);
    // (p, solver, d / sqrt(p) range): the trivariate solver only reaches |x_i| of a few units,
    // so it gets the smallest epsilon
    for (p, cases) in [
        (
            Int::from((1i128 << 61) - 1),
            vec![
                (Solver::Cornacchia, 100u64),
                (Solver::CoppersmithBivariate, 100),
                (Solver::CoppersmithTrivariate, 4),
            ],
        ),
        (
            &Int::from(1i64).shl(127) - &Int::one(),
            // (Cornacchia here needs N ~ 2^200 factored: beyond Pollard rho, so it is not run)
            vec![(Solver::CoppersmithBivariate, 30)],
        ),
    ] {
        let alg = Alg::new(&p);
        let o0 = alg.o0();
        let next_prime = |mut v: Int| {
            while !v.is_probable_prime() {
                v = &v + &Int::one();
            }
            v
        };
        for (solver, mult) in cases {
            let q1 = next_prime(&p.isqrt() * &Int::from(3i64));
            let i1 = reduce_ideal(
                &alg,
                &o0,
                &random_prime_norm_ideal(&alg, &o0, &q1, &mut rng),
            )
            .0;
            let o1 = i1.right_order(&alg, &o0);
            // d ~ p^{1/2 + eps} with a small eps (planted)
            let d = next_prime(&p.isqrt() * &Int::from((mult + rng.below(mult)) as i64));
            let j0 = random_prime_norm_ideal(&alg, &o1, &d, &mut rng);
            let i2 = reduce_ideal(&alg, &o0, &i1.mul(&alg, &j0)).0;
            let o2 = i2.right_order(&alg, &o0);
            let conn = connecting_ideal(&alg, &o1, &o2);
            assert!(connects(&alg, &conn, &o1, &o2));
            let nconn = conn.norm_int(&o1);
            let t = std::time::Instant::now();
            let sol = solve_degree(&alg, &o1, &conn, &d, solver, None)
                .unwrap_or_else(|| panic!("{solver:?} at {} bits", p.bits()));
            let el = t.elapsed();
            assert_eq!(sol.ideal.norm(&o1), (d.clone(), Int::one()));
            assert_eq!(sol.ideal.left_order(&alg, &o1), o1);
            let beta = sol.alpha.conj().scale(&Int::one(), &nconn);
            let (bn, bd) = alg.nrd(&beta);
            let beta_inv = beta.conj().scale(&bd, &bn);
            let conj_o2 = o2.lmul(&alg, &beta_inv).rmul(&alg, &beta);
            assert_eq!(sol.ideal.right_order(&alg, &o1), conj_o2);
            eprintln!(
                "{} bits, {solver:?}: d = {d}, N(conn) = {nconn}, x = {:?}, {el:?}",
                p.bits(),
                sol.x
            );
        }
    }
}
