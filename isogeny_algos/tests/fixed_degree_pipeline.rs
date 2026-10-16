//! End-to-end tests of the fixed-degree isogeny algorithm (Benčina–Kutas–Merz–Petit–Stopar–
//! Weitkämper): End(E) from a path, Kirschmer–Voight connecting ideals, norm forms, the
//! norm-equation solvers inside the pipeline, ideal -> isogeny (kernel and torsion/Kani
//! representations) and the hybrid method. Every check compares against an independent
//! computation: Phi_l neighbours on the j-line, Deuring images, Vélu codomains, Kani splits with
//! a twisted negative control, and Weil pairings.
use isogeny_algos::curve::{from_j, jinv, Curve};
use isogeny_algos::ext::Fp4;
use isogeny_algos::field::{factor_u64, Field, Rng};
use isogeny_algos::fixed_degree::connecting::{
    combine, connecting_ideal, connecting_via_o0, connects, norm_form,
};
use isogeny_algos::fixed_degree::endring::{end_ring, random_ss_j, sub_ideals, JGraph, J};
use isogeny_algos::fixed_degree::hybrid::hybrid_search;
use isogeny_algos::fixed_degree::pipeline::{
    find_isogeny_of_degree, plant_ideal, Instance, Solver,
};
use isogeny_algos::fixed_degree::represent::{
    fast_ideal_to_j, ideal_to_isogeny, is_deuring_smooth, two_adic, Repr,
    ReprOptions, E4,
};
use isogeny_algos::int::Int;
use isogeny_algos::quat::deuring::Deuring;

/// p = 2^15 * 49 - 1: v_2(p + 1) = 15 (torsion level a = 12, 2^a > sqrt(p) ~ 1267), and
/// p^2 - 1 has a large 256-smooth odd part (searched for; see the report).
const P: u64 = 1605631;

fn setup(seed: u64) -> (Deuring, JGraph, Rng) {
    let mut rng = Rng::new(seed);
    let dg = Deuring::new(P, &mut rng);
    let g = JGraph::new(P, 2);
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
        let er = end_ring(&dg, &mut g, j, &mut rng).expect("end ring");
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
        let i1 = end_ring(&dg, &mut g, j1, &mut rng).unwrap().ideal;
        let i2 = end_ring(&dg, &mut g, j2, &mut rng).unwrap().ideal;
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

fn planted(
    dg: &Deuring,
    g: &mut JGraph,
    rng: &mut Rng,
    primes: &[u64],
) -> (J, J, Int) {
    let j1 = random_ss_j(g, P, 30, rng);
    let i1 = end_ring(dg, g, j1, rng).unwrap().ideal;
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

/// Smooth d (every prime power has rational torsion): kernel representation, Vélu chain from
/// the caller's E_1 lands on E_2; random E_2 (no planting) with d > p.
#[test]
fn pipeline_smooth_degree_kernel_random_target() {
    let (dg, mut g, mut rng) = setup(5);
    let f = dg.f2.clone();
    // a divisor of t_odd a little above 4p (the isogeny count between two curves is ~ 12 d / p)
    let mut ds: Vec<u64> = vec![1];
    for t in &dg.tors {
        let mut nd = vec![];
        for &x in &ds {
            let mut y = x;
            for _ in 0..=t.e {
                nd.push(y);
                y = y.saturating_mul(t.ell);
            }
        }
        ds = nd;
    }
    ds.sort();
    ds.dedup();
    let dsel: Vec<u64> = ds
        .into_iter()
        .filter(|&x| x > 4 * P && x < 64 * P)
        .collect();
    assert!(!dsel.is_empty());
    let mut ok = 0;
    for t in 0..2 {
        let j1 = random_ss_j(&mut g, P, 30, &mut rng);
        let j2 = random_ss_j(&mut g, P, 31, &mut rng);
        let d = dsel[rng.below(dsel.len() as u64) as usize];
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
                assert!(matches!(sol.isogeny.repr, Repr::Kernel { .. }));
                assert!(sol.isogeny.j_ok, "Vélu codomain j");
                eprintln!(
                    "trial {t}: d = {d} = {:?}, x = {:?}, times {:?}",
                    factor_u64(d),
                    sol.degree.x,
                    sol.times
                );
                ok += 1;
            }
            Err(e) => eprintln!("trial {t}: d = {d}: {e}"),
        }
    }
    assert!(ok >= 1);
}

/// The Coppersmith solvers inside the pipeline (bivariate: Coron; trivariate), planted d.
#[test]
fn pipeline_coppersmith_solvers() {
    let (dg, mut g, mut rng) = setup(6);
    let f = dg.f2.clone();
    let sqrt_p = (P as f64).sqrt() as u64;
    for (solver, lo, hi) in [
        (Solver::CoppersmithBivariate, sqrt_p, 2 * sqrt_p),
        (Solver::CoppersmithTrivariate, sqrt_p / 2, sqrt_p),
    ] {
        let mut q = lo + rng.below(hi - lo);
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
        let sol = find_isogeny_of_degree(&dg, &mut g, &inst, solver, &ReprOptions::default(), &mut rng)
            .unwrap_or_else(|e| panic!("{solver:?}: {e}"));
        assert!(sol.isogeny.j_ok);
        assert_eq!(sol.isogeny.pairing_ok, Some(true));
        if let Repr::Torsion { certificate, .. } = &sol.isogeny.repr {
            assert!(certificate.certified(), "{certificate:?}");
        }
        eprintln!("{solver:?}: d = {q}, x = {:?}, times {:?}", sol.degree.x, sol.times);
    }
}

/// Hybrid: d = 3^e, guess a 3^{e_1} prefix, solve the rest, compose; Deuring image and
/// isogeny checked.
#[test]
fn hybrid_three_power_degree() {
    let (dg, mut g, mut rng) = setup(7);
    let (alg, o0) = (&dg.alg, &dg.o0);
    // 3^9 = 19683 > sqrt(p) ~ 863; 3^9 exceeds both the rational 3-torsion and 2^a, so the
    // representation is torsion images checked by the pairing (no Kani certificate)
    let e = 9u32;
    let (j1, j2, _) = planted(&dg, &mut g, &mut rng, &vec![3; e as usize]);
    let i1 = end_ring(&dg, &mut g, j1, &mut rng).unwrap().ideal;
    let i2 = end_ring(&dg, &mut g, j2, &mut rng).unwrap().ideal;
    let res = hybrid_search(alg, o0, &i1, &i2, 3, e, 3, Solver::Cornacchia).expect("hybrid");
    assert_eq!(
        dg.ideal_to_j(&i1.mul(alg, &res.ideal), &mut rng),
        Some(j2),
        "hybrid ideal lands on E_2"
    );
    eprintln!("hybrid: {} prefixes tried", res.guesses);
    let rep = ideal_to_isogeny(&dg, &i1, &res.ideal, None, &ReprOptions::default(), &mut rng)
        .expect("isogeny");
    assert!(rep.j_ok);
}
