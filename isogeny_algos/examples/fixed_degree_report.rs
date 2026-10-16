//! Measurement harness for the fixed-degree isogeny study (research/fixed_degree_isogenies_20261016).
//!
//! Prints one JSON object per measurement (stdout) for:
//!  * `pipeline`: the end-to-end algorithm (steps 1–6) at p = 1605631 on planted prime degrees
//!    sqrt(p) < d < 2^a (Cornacchia, and Coron's bivariate method) and on random targets with
//!    d ~ 8p (Cornacchia), with per-step wall times and every check's outcome;
//!  * `steps25`: steps 2–5 at p = 2^61 - 1 and 2^127 - 1 for each norm-equation solver against
//!    d / sqrt(p);
//!  * `baseline`: MITM and van Oorschot–Wiener claw search for d = 2^e (counts);
//!  * `strassen`: l! mod (2^61 - 1), naive vs Strassen, wall time;
//!  * `csurf`: one CSURF-512 action.
//!
//! ```bash
//! cargo run --release --example fixed_degree_report > records.jsonl
//! ```
//! Every value is computed natively here; the plots only lay out these records.
use isogeny_algos::curve::from_j;
use isogeny_algos::field::{is_prime, Rng, Zp};
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::fixed_degree::baselines::{
    mitm_smooth, vow_claw, vow_expected_evals, VowParams,
};
use isogeny_algos::fixed_degree::connecting::{connecting_ideal, connects};
use isogeny_algos::fixed_degree::endring::{end_ring, random_ss_j, reduce_ideal, JGraph, J};
use isogeny_algos::fixed_degree::pipeline::{
    find_isogeny_of_degree, plant_ideal, random_prime_norm_ideal, solve_degree, Instance, Solver,
};
use isogeny_algos::fixed_degree::represent::{fast_ideal_to_j, two_adic, Repr, ReprOptions};
use isogeny_algos::fpm::FpM;
use isogeny_algos::int::Int;
use isogeny_algos::path::csurf::Csurf;
use isogeny_algos::quat::deuring::Deuring;
use isogeny_algos::quat::Alg;
use std::time::Instant;

const P: u64 = 1605631;

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

fn pipeline_records() {
    let mut rng = Rng::new(20261016);
    let dg = Deuring::new(P, &mut rng);
    let mut g = JGraph::new(P, 2);
    let f = dg.f2.clone();
    let a = two_adic(P) - 3;
    let sqrt_p = (P as f64).sqrt() as u64;
    let lp = (P as f64).ln();
    // planted degrees spread over (sqrt p, 2^a), and random targets with d ~ 8p
    let mut cases: Vec<(Solver, u64, bool)> = vec![];
    for i in 0..6u64 {
        let lo = sqrt_p + i * ((1u64 << a) - sqrt_p) / 6;
        let mut q = lo + rng.below(((1u64 << a) - sqrt_p) / 6);
        while !is_prime(q) {
            q += 1;
        }
        cases.push((Solver::Cornacchia, q, true));
    }
    for _ in 0..2 {
        let mut q = sqrt_p + rng.below(sqrt_p);
        while !is_prime(q) {
            q += 1;
        }
        cases.push((Solver::CoppersmithBivariate, q, true));
    }
    for _ in 0..3 {
        let mut q = 8 * P + rng.below(4 * P);
        while !is_prime(q) {
            q += 1;
        }
        cases.push((Solver::Cornacchia, q, false));
    }
    for (idx, (solver, d, planted)) in cases.into_iter().enumerate() {
        let j1 = random_ss_j(&mut g, P, 30, &mut rng);
        let j2: J = if planted {
            let i1 = end_ring(&dg, &mut g, j1).unwrap().ideal;
            let jd = plant_ideal(&dg.alg, &dg.o0, &i1, &[Int::from(d)], &mut rng);
            fast_ideal_to_j(&dg, &i1.mul(&dg.alg, &jd)).unwrap()
        } else {
            random_ss_j(&mut g, P, 31, &mut rng)
        };
        let inst = Instance {
            j1,
            j2,
            e1: Some(from_j(&f, f.embed(j1))),
            e2: Some(from_j(&f, f.embed(j2))),
            d: Int::from(d),
        };
        let t = Instant::now();
        let res = find_isogeny_of_degree(
            &dg,
            &mut g,
            &inst,
            solver,
            &ReprOptions::default(),
            &mut rng,
        );
        let total = ms(t);
        match res {
            Ok(sol) => {
                let (repr, cert) = match &sol.isogeny.repr {
                    Repr::Kernel { .. } => ("kernel".to_string(), "n/a".to_string()),
                    Repr::Torsion { certificate, .. } => (
                        "torsion".to_string(),
                        if certificate.certified() {
                            "kani2_certified".to_string()
                        } else {
                            format!("{certificate:?}").replace('"', "'")
                        },
                    ),
                };
                println!(
                    "{{\"kind\":\"pipeline\",\"case\":{idx},\"p\":{P},\"solver\":\"{solver:?}\",\"planted\":{planted},\"d\":{d},\"log_p_d\":{:.4},\"ok\":true,\"end_ring_ms\":{:.1},\"connecting_ms\":{:.2},\"solve_ms\":{:.2},\"isogeny_ms\":{:.1},\"total_ms\":{:.1},\"path1\":{},\"path2\":{},\"graph_vertices\":{},\"deuring_evals\":{},\"n_conn\":{},\"x\":\"{:?}\",\"repr\":\"{repr}\",\"certificate\":\"{cert}\",\"pairing_ok\":{},\"j_ok\":{}}}",
                    (d as f64).ln() / lp,
                    sol.times.end_ring_ms,
                    sol.times.connecting_ms,
                    sol.times.solve_ms,
                    sol.times.isogeny_ms,
                    total,
                    sol.end1.path.len() - 1,
                    sol.end2.path.len() - 1,
                    sol.end1.graph_vertices + sol.end2.graph_vertices,
                    sol.end1.deuring_evals + sol.end2.deuring_evals,
                    sol.connecting_norm,
                    sol.degree.x,
                    sol.isogeny.pairing_ok.map_or("null".to_string(), |b| b.to_string()),
                    sol.isogeny.j_ok
                );
            }
            Err(e) => {
                // diagnostic: does Cornacchia find a solution on the same instance?
                let alt = if solver != Solver::Cornacchia {
                    find_isogeny_of_degree(
                        &dg,
                        &mut g,
                        &inst,
                        Solver::Cornacchia,
                        &ReprOptions::default(),
                        &mut rng,
                    )
                    .map(|s| format!("{:?}", s.degree.x))
                    .unwrap_or_else(|e| e.to_string())
                } else {
                    "n/a".to_string()
                };
                println!(
                    "{{\"kind\":\"pipeline\",\"case\":{idx},\"p\":{P},\"solver\":\"{solver:?}\",\"planted\":{planted},\"d\":{d},\"log_p_d\":{:.4},\"ok\":false,\"error\":\"{e}\",\"total_ms\":{total:.1},\"cornacchia_same_instance\":\"{alt}\"}}",
                    (d as f64).ln() / lp
                )
            }
        }
    }
}

fn steps25_records() {
    let mut rng = Rng::new(61);
    let next_prime = |mut v: Int| {
        while !v.is_probable_prime() {
            v = &v + &Int::one();
        }
        v
    };
    let p61 = Int::from((1i128 << 61) - 1);
    let p127 = &Int::from(1i64).shl(127) - &Int::one();
    let mut cases: Vec<(Int, Solver, u64)> = vec![];
    for m in [4u64, 16, 64, 256] {
        cases.push((p61.clone(), Solver::Cornacchia, m));
        cases.push((p61.clone(), Solver::CoppersmithBivariate, m));
    }
    for m in [2u64, 4, 8] {
        cases.push((p61.clone(), Solver::CoppersmithTrivariate, m));
    }
    for m in [4u64, 16, 64] {
        cases.push((p127.clone(), Solver::CoppersmithBivariate, m));
    }
    for (p, solver, mult) in cases {
        let alg = Alg::new(&p);
        let o0 = alg.o0();
        for rep in 0..2 {
            let q1 = next_prime(&p.isqrt() * &Int::from(3i64));
            let i1 = reduce_ideal(
                &alg,
                &o0,
                &random_prime_norm_ideal(&alg, &o0, &q1, &mut rng),
            )
            .0;
            let o1 = i1.right_order(&alg, &o0);
            let d = next_prime(&p.isqrt() * &Int::from((mult + rng.below(mult)) as i64));
            let j0 = random_prime_norm_ideal(&alg, &o1, &d, &mut rng);
            let i2 = reduce_ideal(&alg, &o0, &i1.mul(&alg, &j0)).0;
            let o2 = i2.right_order(&alg, &o0);
            let conn = connecting_ideal(&alg, &o1, &o2);
            assert!(connects(&alg, &conn, &o1, &o2));
            let t = Instant::now();
            let sol = solve_degree(&alg, &o1, &conn, &d, solver, None);
            let el = ms(t);
            let xmax = sol
                .as_ref()
                .map(|s| s.x.iter().map(|v| v.abs()).max().unwrap().to_string())
                .unwrap_or("null".into());
            println!(
                "{{\"kind\":\"steps25\",\"p_bits\":{},\"solver\":\"{solver:?}\",\"rep\":{rep},\"d_over_sqrtp\":{:.2},\"log_p_d\":{:.4},\"ok\":{},\"ms\":{el:.1},\"x_max\":{xmax}}}",
                p.bits(),
                d.to_f64() / p.isqrt().to_f64(),
                d.to_f64().ln() / p.to_f64().ln(),
                sol.is_some()
            );
        }
    }
}

fn baseline_records() {
    // p = 2^22-ish, p = 3 mod 4: planted 2^e walks
    let p = (4_000_000u64..)
        .find(|&q| q % 4 == 3 && is_prime(q))
        .unwrap();
    let fp = Zp::new(p);
    let phi = Phi::compute(&fp, 2).lift(|x| (x, 0u64));
    let mut g = JGraph::new(p, 2);
    let mut rng = Rng::new(22);
    for e in [12usize, 14, 16, 18] {
        for rep in 0..2 {
            let j1 = random_ss_j(&mut g, p, 25, &mut rng);
            // planted non-backtracking walk of length e
            let mut prev = None;
            let mut cur = j1;
            for _ in 0..e {
                let nb: Vec<J> = g
                    .neighbours(cur)
                    .into_iter()
                    .filter(|&v| Some(v) != prev)
                    .collect();
                let nx = nb[rng.below(nb.len() as u64) as usize];
                prev = Some(cur);
                cur = nx;
            }
            let j2 = cur;
            let t = Instant::now();
            let (path, st) = mitm_smooth(&g.f2, &phi, e, j1, j2);
            let mitm_ms = ms(t);
            let w = 16usize;
            let t = Instant::now();
            let (vpath, vst) = vow_claw(&g.f2, &phi, e, j1, j2, &VowParams::new(w));
            let vow_ms = ms(t);
            let space = 2 * 3 * (1u64 << (e / 2 - 1));
            println!(
                "{{\"kind\":\"baseline\",\"p\":{p},\"e\":{e},\"rep\":{rep},\"mitm_found\":{},\"mitm_leaves\":{},\"mitm_ms\":{mitm_ms:.1},\"vow_found\":{},\"vow_f_evals\":{},\"vow_expected\":{:.0},\"vow_w\":{w},\"vow_ms\":{vow_ms:.1}}}",
                path.is_some(),
                st.leaves1 + st.leaves2,
                vpath.is_some(),
                vst.f_evals,
                vow_expected_evals(space, w)
            );
        }
    }
}

fn strassen_records() {
    use isogeny_algos::strassen::{factorial_mod, factorial_naive};
    let n = (1u64 << 61) - 1;
    for l in [1_000u64, 10_000, 100_000, 1_000_000, 4_000_000] {
        let reps = (4_000_000 / l).clamp(1, 50);
        let t = Instant::now();
        let mut a = 0;
        for _ in 0..reps {
            a = factorial_naive(l, n);
        }
        let naive = ms(t) / reps as f64;
        let t = Instant::now();
        let mut b = 0;
        for _ in 0..reps {
            b = factorial_mod(l, n);
        }
        let fast = ms(t) / reps as f64;
        assert_eq!(a, b);
        println!(
            "{{\"kind\":\"strassen\",\"l\":{l},\"n_bits\":61,\"reps\":{reps},\"naive_ms\":{naive:.4},\"strassen_ms\":{fast:.4}}}"
        );
    }
}

fn csurf_record() {
    let mut cs = Csurf::<FpM<9>>::csurf512();
    let mut rng = Rng::new(512);
    let n = cs.primes.len() + 1;
    let key: Vec<i32> = (0..n).map(|_| rng.below(3) as i32 - 1).collect();
    use isogeny_algos::field::Field;
    let e0 = cs.fp.zero();
    let t = Instant::now();
    let a = cs.action_batched(e0, &key, &mut rng);
    let el = ms(t);
    let neg: Vec<i32> = key.iter().map(|&x| -x).collect();
    let back = cs.action_batched(a, &neg, &mut rng);
    let _ = &mut cs;
    println!(
        "{{\"kind\":\"csurf\",\"prime\":\"CSURF-512\",\"key_bound\":1,\"ms\":{el:.1},\"inverse_returns\":{}}}",
        back == e0
    );
}

fn main() {
    let which: Vec<String> = std::env::args().skip(1).collect();
    let all = which.is_empty();
    let on = |k: &str| all || which.iter().any(|w| w == k);
    if on("strassen") {
        strassen_records();
    }
    if on("csurf") {
        csurf_record();
    }
    if on("steps25") {
        steps25_records();
    }
    if on("baseline") {
        baseline_records();
    }
    if on("pipeline") {
        pipeline_records();
    }
}
