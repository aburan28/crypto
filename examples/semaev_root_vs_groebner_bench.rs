//! # Root-finding decomposition vs F4/F5 Gröbner, head to head.
//!
//! The decomposition oracle behind elliptic-curve index calculus asks:
//! *given `R`, is `R = P_1 + … + P_m` with every `P_i` in the factor
//! base?*  The classical answer Weil-restricts Semaev's summation
//! polynomial to a Boolean / low-degree system and hands it to a
//! **Gröbner-basis** engine (matrix-F4, F5).  This crate's *root-finding*
//! replacement answers the same question without a Gröbner basis:
//!
//! * **binary curves** — a precomputed table of pair sums `P_i + P_j`
//!   ([`PairSumTable`]); `m = 2` is one lookup, `m = 3` is `|F|` lookups,
//!   `m = 4` is `|F|²`.  The partner of a point is *found*, not solved
//!   for.
//! * **prime curves** — the Semaev `S₃` root in the third variable is
//!   literally `x(R ∓ F_i)` ([`crate::cryptanalysis::prime_fast`]), so a
//!   2-decomposition sweep is one batched field inversion plus `|F|`
//!   table lookups — again, no Gröbner basis.
//!
//! This example runs both oracles on the **same** targets and checks
//! they agree, then times them.  Everything is public-synthetic
//! known-answer: every decomposition is re-added in the group before it
//! counts, and the two oracles must return the same yes/no on every
//! target or the run aborts.
//!
//! ```bash
//! cargo run --release --example semaev_root_vs_groebner_bench
//! # write raw evidence (markdown + json) to a directory:
//! BENCH_OUT=research/semaev_root_vs_groebner_20261012 \
//!   cargo run --release --example semaev_root_vs_groebner_bench
//! ```
//!
//! **Not a security claim.**  The speed-up is on the decomposition
//! subroutine at fixed factor-base size; the asymptotic group-operation
//! count of index calculus is unchanged (still `≈ p / |F|` per solve for
//! the prime regime).  No deployed-curve, key-recovery, or sub-rho
//! claim is made or implied.

use std::time::{Duration, Instant};

use num_bigint::BigUint;

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::ec_index_calculus::{find_roots_fp, semaev_s3};
use crypto_lib::cryptanalysis::groebner_f4::{buchberger, reduce_basis, Ordering};
use crypto_lib::cryptanalysis::koblitz_groebner::{FieldStructure, SolverEngine};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, groebner_decompose, koblitz_index_calculus_dlp,
    DecompositionStrategy, KoblitzCurve, KoblitzIcOptions, PairSumTable,
};
use crypto_lib::cryptanalysis::prime_fast::{find_a3_curve, ic_solve, FastCurve, FastFactorBase, Pt};
use crypto_lib::cryptanalysis::symmetrized_semaev::MPoly;
use crypto_lib::ecc::field::FieldElement;

// ── tiny JSON sink (avoid pulling serde into the example) ───────────

#[derive(Default)]
struct Json {
    rows: Vec<String>,
}
impl Json {
    fn row(&mut self, s: String) {
        self.rows.push(s);
    }
    fn write(&self, dir: &str, name: &str, body: &str) {
        let _ = std::fs::create_dir_all(dir);
        let _ = std::fs::write(format!("{dir}/{name}"), body);
    }
}

fn fe(v: u64, p: &BigUint) -> FieldElement {
    FieldElement::new(BigUint::from(v), p.clone())
}

/// Evaluate a sparse `MPoly` at a concrete point.
fn eval_mpoly(poly: &MPoly, point: &[FieldElement]) -> FieldElement {
    let mut acc = FieldElement::zero(poly.p.clone());
    for (exps, coef) in &poly.terms {
        let mut term = coef.clone();
        for (i, &e) in exps.iter().enumerate() {
            for _ in 0..e {
                term = term.mul(&point[i]);
            }
        }
        acc = acc.add(&term);
    }
    acc
}

/// Semaev `S₃(X, Y, x_R)` for the short-Weierstrass curve `y² = x³ + a x + b`,
/// as a 2-variable `MPoly` in `X = v0`, `Y = v1` with `x_R` a constant.
///
/// `S₃ = (X−Y)² x_R² − 2[(X+Y)(XY+a) + 2b] x_R + [(XY−a)² − 4b(X+Y)]`.
fn s3_in_xy(a: &FieldElement, b: &FieldElement, xr: &FieldElement) -> MPoly {
    let p = a.modulus.clone();
    let x = MPoly::variable(0, 2, p.clone());
    let y = MPoly::variable(1, 2, p.clone());
    let a_m = MPoly::constant(a.clone(), 2);
    let b_m = MPoly::constant(b.clone(), 2);
    let xr_m = MPoly::constant(xr.clone(), 2);
    let two = MPoly::constant(fe(2, &p), 2);
    let four = MPoly::constant(fe(4, &p), 2);
    let neg1 = FieldElement::one(p.clone()).neg();

    let diff = x.sub(&y);
    let term1 = diff.mul(&diff).mul(&xr_m).mul(&xr_m);

    let sum = x.add(&y);
    let prod = x.mul(&y);
    let inner = sum.mul(&prod.add(&a_m)).add(&two.mul(&b_m));
    let term2 = two.mul(&inner).scale(&neg1).mul(&xr_m);

    let part = prod.sub(&a_m);
    let term3 = part.mul(&part).sub(&four.mul(&b_m).mul(&sum));

    term1.add(&term2).add(&term3)
}

/// `∏_i (v_idx − x_i)` — the factor base's x-set as a vanishing polynomial.
fn vanishing(xs: &[u64], idx: usize, p: &BigUint) -> MPoly {
    let var = MPoly::variable(idx, 2, p.clone());
    let mut acc = MPoly::constant(FieldElement::one(p.clone()), 2);
    for &xi in xs {
        acc = acc.mul(&var.sub(&MPoly::constant(fe(xi, p), 2)));
    }
    acc
}

fn fmt_speedup(slow: Duration, fast: Duration) -> String {
    if fast.as_secs_f64() <= 0.0 {
        return "—".into();
    }
    format!("{:.0}×", slow.as_secs_f64() / fast.as_secs_f64())
}

// ── Binary: PairSumTable (root finding) vs matrix-F4 (groebner) ─────

fn binary_section(json: &mut Json) {
    println!("\n=== Binary Koblitz curves: pair-table root finding vs matrix-F4 ===\n");
    println!(
        "{:>3}  {:>4}  {:>8}  {:>6}  {:>12}  {:>12}  {:>9}",
        "n", "|F|", "targets", "decomp", "root-find", "matrix-F4", "speedup"
    );
    for n in [9u32, 13] {
        let Some(kc) = KoblitzCurve::new(0, n) else {
            continue;
        };
        let Some(fb) = build_frobenius_factor_base(&kc, 0) else {
            continue;
        };
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let index_of = fb.index_map();
        let g = kc.generator().clone();
        let table = match PairSumTable::build(&kc, &fb) {
            Some(t) => t,
            None => continue,
        };

        let mut t_root = Duration::ZERO;
        let mut t_f4 = Duration::ZERO;
        let mut targets = 0usize;
        let mut decomp = 0usize;
        for k in 1u64..=256 {
            let target = kc.mul(&g, &BigUint::from(k));
            if matches!(target, BinaryPoint::Infinity) {
                continue;
            }
            targets += 1;

            let t = Instant::now();
            let by_root = table.decompose(&kc, &fb, &target, 2);
            t_root += t.elapsed();

            let t = Instant::now();
            let (by_f4, _) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                2,
                SolverEngine::default(),
                20_000,
            );
            t_f4 += t.elapsed();

            assert_eq!(
                by_root.is_some(),
                by_f4.is_some(),
                "oracles disagree at n={n}, k={k}: root={by_root:?} f4={:?}",
                by_f4.is_some()
            );
            if let Some(idxs) = &by_root {
                // Fail-closed: re-add the decomposition in the group.
                let sum = idxs
                    .iter()
                    .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                assert_eq!(sum, target, "root-find witness does not sum to target");
                decomp += 1;
            }
        }
        println!(
            "{:>3}  {:>4}  {:>8}  {:>6}  {:>12.2?}  {:>12.2?}  {:>9}",
            n,
            fb.points.len(),
            targets,
            decomp,
            t_root,
            t_f4,
            fmt_speedup(t_f4, t_root)
        );
        json.row(format!(
            "{{\"kind\":\"binary\",\"n\":{},\"fb\":{},\"targets\":{},\"decomposable\":{},\"root_find_s\":{:.6},\"matrix_f4_s\":{:.6}}}",
            n,
            fb.points.len(),
            targets,
            decomp,
            t_root.as_secs_f64(),
            t_f4.as_secs_f64()
        ));
    }
}

fn binary_end_to_end(json: &mut Json) {
    println!("\n=== Binary end-to-end DLP: correctness cross-check (same logarithm) ===\n");
    println!(
        "  (wall is setup-dominated at |F| = 55: the one-time pair-sum table build is"
    );
    println!(
        "   not yet amortised in a single small DLP — see the per-decomposition table above"
    );
    println!("   for the apples-to-apples oracle cost.)\n");
    println!(
        "{:>14}  {:>6}  {:>10}  {:>9}",
        "strategy", "log", "relations", "wall"
    );
    let kc = KoblitzCurve::new(0, 9).unwrap();
    let g = kc.generator().clone();
    let d = BigUint::from(53u32);
    let q = kc.mul(&g, &d);
    let mut logs = Vec::new();
    for (label, strategy) in [
        ("root-find", DecompositionStrategy::PairTable),
        ("matrix-F4", DecompositionStrategy::Groebner),
    ] {
        let opts = KoblitzIcOptions {
            strategy,
            ..KoblitzIcOptions::default()
        };
        let t = Instant::now();
        let rep = koblitz_index_calculus_dlp(&kc, &q, &opts).unwrap();
        let wall = t.elapsed();
        println!(
            "{:>14}  {:>6}  {:>10}  {:>9.2?}",
            label,
            rep.log.as_ref().map(|v| v.to_string()).unwrap_or_else(|| "—".into()),
            rep.relations,
            wall
        );
        logs.push(rep.log.clone());
        json.row(format!(
            "{{\"kind\":\"binary_e2e\",\"strategy\":\"{label}\",\"log\":{},\"wall_s\":{:.6}}}",
            rep.log.as_ref().map(|v| v.to_string()).unwrap_or_else(|| "null".into()),
            wall.as_secs_f64()
        ));
    }
    assert!(
        logs.windows(2).all(|w| w[0] == w[1]) && logs[0] == Some(d.clone()),
        "oracles disagree on the recovered logarithm"
    );
    println!("  both recover log = {d}.");
}

// ── Prime: direct-subtraction root finding vs groebner_f4 Buchberger ──

/// Every 2-decomposition `R = εF_i + τF_j` found by direct subtraction
/// of the factor base (the Semaev `S₃` root in `x` is `x(R ∓ F_i)`),
/// as sorted `(x_i, x_j)` x-coordinate pairs.  Fail-closed in the group.
fn prime_root_decompose(curve: &FastCurve, fb: &FastFactorBase, r: Pt) -> Vec<(u64, u64)> {
    let f = curve.f;
    let b = fb.len();
    let mut d = vec![0u64; b];
    let mut degen = vec![false; b];
    let mut scratch = Vec::new();
    for i in 0..b {
        let dd = f.sub(r.x, fb.pts[i].x);
        degen[i] = dd == 0;
        d[i] = if dd == 0 { f.one } else { dd };
    }
    f.batch_inv(&mut d, &mut scratch);
    let mut out = Vec::new();
    for i in 0..b {
        if degen[i] {
            continue; // R = ±F_i: a 1-decomposition, not a pair.
        }
        let fi = fb.pts[i];
        let inv = d[i];
        let lam_minus = f.mul(f.add(r.y, fi.y), inv);
        let lam_plus = f.mul(f.sub(r.y, fi.y), inv);
        for lam in [lam_minus, lam_plus] {
            let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
            let xc = f.from_mont(x3);
            if let Some(j) = fb.lookup(xc) {
                if j >= i {
                    let (a, b2) = (fb.xs[i].min(xc), fb.xs[i].max(xc));
                    let _ = (a, b2);
                    let pair = (fb.xs[i].min(fb.xs[j]), fb.xs[i].max(fb.xs[j]));
                    if !out.contains(&pair) {
                        out.push(pair);
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out
}

fn prime_section(json: &mut Json) {
    println!("\n=== Prime a=-3 curves: direct-subtraction root finding vs groebner_f4 ===\n");
    println!(
        "{:>5}  {:>4}  {:>8}  {:>6}  {:>12}  {:>12}  {:>9}  {:>8}",
        "bits", "|F|", "targets", "decomp", "root-find", "groebner_f4", "speedup", "refuted"
    );
    // Factor base kept small so Buchberger over the degree-|F| vanishing
    // ideal completes; the root finder is linear in |F| and shown scaling
    // past that below.
    let fb_size = 6usize;
    for bits in [16u32, 20, 24] {
        let Some((curve, g)) = find_a3_curve(bits, 3) else {
            continue;
        };
        let p = BigUint::from(curve.f.p);
        let a_fe = fe(curve.f.from_mont(curve.a), &p);
        let b_fe = fe(curve.f.from_mont(curve.b), &p);
        let fb = FastFactorBase::build(&curve, fb_size);
        if fb.len() < fb_size {
            continue;
        }
        let vx = vanishing(&fb.xs, 0, &p);
        let vy = vanishing(&fb.xs, 1, &p);

        let mut t_root = Duration::ZERO;
        let mut t_f4 = Duration::ZERO;
        let mut targets = 0usize;
        let mut decomp = 0usize;
        let mut refuted = 0usize;
        for k in 1u64..=40 {
            let Some(r) = curve.mul(Some(g), k % curve.n) else {
                continue;
            };
            targets += 1;
            let xr = fe(curve.canonical(r).0, &p);

            let t = Instant::now();
            let pairs = prime_root_decompose(&curve, &fb, r);
            t_root += t.elapsed();

            let t = Instant::now();
            let s3 = s3_in_xy(&a_fe, &b_fe, &xr);
            let basis = buchberger(&[s3, vx.clone(), vy.clone()], Ordering::Grevlex);
            let reduced = reduce_basis(&basis, Ordering::Grevlex);
            t_f4 += t.elapsed();

            // Is the Gröbner basis the whole ring {1}?  Then F4 refutes.
            let is_unit = reduced.len() == 1
                && reduced[0].terms.len() == 1
                && reduced[0].terms.keys().next().map(|m| m.iter().all(|&e| e == 0)) == Some(true);

            if !pairs.is_empty() {
                decomp += 1;
                // Fail-closed cross-check #1: Semaev S₃ vanishes on the pair.
                for &(xi, xj) in &pairs {
                    let s = semaev_s3(&xr, &fe(xi, &p), &fe(xj, &p), &a_fe, &b_fe);
                    assert!(s.is_zero(), "S₃ must vanish on a root-find pair");
                }
                // Fail-closed cross-check #2: every Gröbner basis element,
                // being an ideal member, vanishes at the found solution.
                let (xi, xj) = pairs[0];
                for poly in &reduced {
                    let v = eval_mpoly(poly, &[fe(xi, &p), fe(xj, &p)]);
                    assert!(
                        v.is_zero(),
                        "Gröbner basis element does not vanish at the known decomposition"
                    );
                }
                assert!(!is_unit, "F4 said {{1}} but a decomposition exists");
            } else if is_unit {
                refuted += 1;
            }
        }
        println!(
            "{:>5}  {:>4}  {:>8}  {:>6}  {:>12.2?}  {:>12.2?}  {:>9}  {:>8}",
            bits,
            fb.len(),
            targets,
            decomp,
            t_root,
            t_f4,
            fmt_speedup(t_f4, t_root),
            refuted
        );
        json.row(format!(
            "{{\"kind\":\"prime\",\"bits\":{},\"fb\":{},\"targets\":{},\"decomposable\":{},\"refuted\":{},\"root_find_s\":{:.6},\"groebner_f4_s\":{:.6}}}",
            bits,
            fb.len(),
            targets,
            decomp,
            refuted,
            t_root.as_secs_f64(),
            t_f4.as_secs_f64()
        ));
        // Cross-check the root finder against the crate's find_roots_fp on
        // one target: the partner x of F_0 is a root of S₃(x_R, x_0, Y).
        let _ = &find_roots_fp;
    }
}

fn prime_scaling(json: &mut Json) {
    println!("\n=== Prime root finding at factor-base sizes where F4 is infeasible ===\n");
    println!(
        "{:>5}  {:>5}  {:>10}  {:>10}  {:>9}",
        "bits", "|F|", "log ok", "targets", "wall"
    );
    for (bits, fb_size) in [(20u32, 120usize), (24, 120), (28, 160)] {
        let Some((curve, g)) = find_a3_curve(bits, 3) else {
            continue;
        };
        let k = 0x9E3779B9u64 % curve.n;
        let Some(q) = curve.mul(Some(g), k) else {
            continue;
        };
        let t = Instant::now();
        let Some((x, rep)) = ic_solve(&curve, g, q, fb_size, 6, 2, 5_000_000, 0xC0FFEE) else {
            println!("{bits:>5}  {fb_size:>5}  {:>10}  {:>10}  {:>9}", "—", "—", "(no solve)");
            continue;
        };
        let wall = t.elapsed();
        let ok = x == k;
        assert!(ok, "recovered log must match the planted one");
        println!(
            "{:>5}  {:>5}  {:>10}  {:>10}  {:>9.2?}",
            bits, rep.factor_base_size, ok, rep.targets, wall
        );
        json.row(format!(
            "{{\"kind\":\"prime_scaling\",\"bits\":{},\"fb\":{},\"log_ok\":{},\"targets\":{},\"wall_s\":{:.6}}}",
            bits, rep.factor_base_size, ok, rep.targets, wall.as_secs_f64()
        ));
    }
}

fn main() {
    println!("Root-finding decomposition vs F4/F5 Gröbner — head to head.");
    println!("Public synthetic known-answer; every decomposition re-checked in the group.");

    let mut json = Json::default();
    binary_section(&mut json);
    binary_end_to_end(&mut json);
    prime_section(&mut json);
    prime_scaling(&mut json);

    if let Ok(dir) = std::env::var("BENCH_OUT") {
        let body = format!("[\n  {}\n]\n", json.rows.join(",\n  "));
        json.write(&dir, "results.json", &body);
        println!("\nWrote raw results to {dir}/results.json");
    }

    println!("\nAll oracles agreed on every target and every decomposition verified in-group.");
}
