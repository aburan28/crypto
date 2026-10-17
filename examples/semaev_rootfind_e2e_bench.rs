//! # End-to-end: root-finding vs F4 decomposition inside a complete DLP solve
//!
//! The head-to-head benches measure the decomposition *oracle* per target.
//! This bench measures what a user actually pays: the wall time of a
//! **complete** index-calculus discrete-log computation — factor base,
//! (pair table), relation collection, linear algebra, group verification —
//! with the oracle swapped and everything else held fixed. Every reported
//! solve recovered the planted logarithm and verified `[d]G = Q`.
//!
//! ## The baseline is given every advantage
//!
//! - **Binary (Koblitz).** `DecompositionStrategy::Groebner` with the
//!   crate's matrix-F4 at the *best* `max_degree` found by a per-`n` sweep;
//!   root finding is charged its one-time pair-table build inside every
//!   single-DLP wall time (no amortisation in the headline number).
//! - **Prime.** Two Gröbner baselines on the identical decomposition ideal
//!   `{S₃(X,Y;x_R), ∏(X−x_i), ∏(Y−x_i)}`: the crate's Buchberger (whichever
//!   of grevlex/lex decides faster on a calibration sample), and
//!   **msolve** (Berthomieu–Eder–Safey El Din; a state-of-the-art F4) as an
//!   external process. msolve is charged only its own process's CPU, and
//!   each decide is **paired** with an adjacent control run on a trivial
//!   ideal whose cost is subtracted, so msolve is charged only its F4
//!   work. Both baselines use the same target stream, factor base and
//!   linear algebra as the root finder; only the oracle differs.
//! - All rows are reported, including any the baseline wins.
//!
//! ## What is reported
//!
//! - **Binary** `n ∈ {5, 9, 13, 23}` (`|F| = 21, 55, 4005, 2025`):
//!   single-DLP wall (incl. table build), 3 seeds, plus the table-build cost
//!   on its own so the per-additional-DLP (amortised) cost is visible.
//! - **Prime** 16-bit (`|F| = 8`), 20-bit and 24-bit (`|F| = 12`): root
//!   finding (always measured in full) vs Buchberger vs msolve on the
//!   same seeded target stream. A baseline is **measured** in full when
//!   its projected wall is under `E2E_FULL_BUDGET_S` (default 900 s);
//!   otherwise it is **composed** — sampled per-target cost on 100 targets
//!   from the stream × the root finder's measured target count — and
//!   labelled COMPOSED everywhere it appears.
//! - **Each method at its own best `|F|`** (16-bit): F4's per-decide cost
//!   grows with `|F|` while the number of targets needed shrinks, so F4 is
//!   allowed to pick the `|F|` minimising its composed total; the root
//!   finder is measured at `|F| = 120`.
//!
//! ```bash
//! cargo run --release --example semaev_rootfind_e2e_bench
//! BENCH_OUT=research/semaev_rootfind_e2e_20261012 \
//!   cargo run --release --example semaev_rootfind_e2e_bench
//! # knobs: E2E_PRIME_BITS=16,20,24  E2E_FULL_BUDGET_S=900  E2E_SKIP_BINARY=1
//! #        E2E_SKIP_BUCHBERGER=1  E2E_BEST_REPS_DIV=4
//! ```
//!
//! **All times are CPU time** (user+sys of this process plus reaped
//! children, via `getrusage`), not wall-clock: the bench is single-threaded,
//! msolve runs as a child, and CPU time is immune to scheduling delay on a
//! shared machine. Wall-clock under load inflates long F4 solves far more
//! than µs lookups and would flatter the root finder.
//!
//! **Not a security claim.** Public synthetic known-answer instances; the
//! speed-up is a constant on the decomposition stage of an algorithm whose
//! asymptotics are unchanged.

use std::io::Write;
use std::process::Command;
use std::time::Duration;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::groebner_f4::{buchberger, reduce_basis, Ordering};
use crypto_lib::cryptanalysis::koblitz_groebner::SolverEngine;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, individual_log_with_pair_table, koblitz_index_calculus_dlp,
    solve_factor_base_logs, DecompositionStrategy, IndividualLogSolver, KoblitzCurve,
    KoblitzIcOptions, PairSumTable,
};
use crypto_lib::cryptanalysis::prime_fast::{find_a3_curve, FastCurve, FastFactorBase, Pt};
use crypto_lib::cryptanalysis::symmetrized_semaev::MPoly;
use crypto_lib::ecc::field::FieldElement;

// ── CPU clock ───────────────────────────────────────────────────────
//
// All timings are **CPU time** (user + system) consumed by this process
// plus its reaped children, via `getrusage`. The bench is single-threaded
// and the external F4 (msolve) runs as a child, so this charges each
// oracle exactly the compute it used and is immune to scheduling delay
// on a shared machine — which wall-clock is not (an F4 solve that runs
// for seconds suffers far more from contention than a µs lookup does,
// so wall-clock under load would flatter the root finder).
//
// The external F4 (msolve) is charged **only its own process's CPU**
// (`RUSAGE_CHILDREN` delta across the spawn), never the harness's
// fork/exec and file I/O; and every decide is paired with an immediately
// adjacent control run of msolve on a trivial ideal, so that what is
// charged is the mean paired difference — the F4 work alone. Two earlier
// designs (an up-front floor subtracted from harness-inclusive CPU, then
// from child CPU) could not resolve the 16-bit F4 cost: the per-process
// floor is ≈5 ms with σ ≈ 1 ms and drifts by hundreds of µs across the
// phases of a long run, which is the size of the signal.

fn rusage(who: libc::c_int) -> Duration {
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    unsafe { libc::getrusage(who, &mut ru) };
    let tv = |t: libc::timeval| Duration::new(t.tv_sec as u64, (t.tv_usec as u32) * 1000);
    tv(ru.ru_utime) + tv(ru.ru_stime)
}

#[derive(Clone, Copy)]
struct Instant(Duration);
impl Instant {
    fn now() -> Self {
        Instant(rusage(libc::RUSAGE_SELF) + rusage(libc::RUSAGE_CHILDREN))
    }
    fn elapsed(&self) -> Duration {
        Instant::now().0.saturating_sub(self.0)
    }
}

/// CPU consumed by reaped children only (what an external solver costs).
fn children_cpu() -> Duration {
    rusage(libc::RUSAGE_CHILDREN)
}

// ── helpers ─────────────────────────────────────────────────────────

#[derive(Default)]
struct Json {
    rows: Vec<String>,
}

fn fmt_dur(d: Duration) -> String {
    let us = d.as_secs_f64() * 1e6;
    if us < 1_000.0 {
        format!("{us:.1} µs")
    } else if us < 1_000_000.0 {
        format!("{:.2} ms", us / 1e3)
    } else {
        format!("{:.2} s", us / 1e6)
    }
}

fn ratio(slow: Duration, fast: Duration) -> String {
    if fast.is_zero() {
        return "—".into();
    }
    let r = slow.as_secs_f64() / fast.as_secs_f64();
    if r >= 1.0 {
        format!("{r:.0}×")
    } else {
        format!("{:.2}× (baseline faster)", r)
    }
}

fn env_list(key: &str, default: &[u32]) -> Vec<u32> {
    std::env::var(key)
        .ok()
        .map(|v| v.split(',').filter_map(|s| s.trim().parse().ok()).collect())
        .unwrap_or_else(|| default.to_vec())
}

// ── Binary: complete Koblitz DLP, PairTable vs matrix-F4 ────────────

fn binary_section(json: &mut Json) {
    println!("=== Binary Koblitz K_0/F_2^n: complete DLP, root finding vs matrix-F4 ===\n");
    println!("Single-DLP CPU time includes the pair-table build (no amortisation). F4 at its best max_degree per n.");
    println!(
        "{:>3} {:>5} {:>5} {:>7} | {:>12} {:>12} | {:>12} {:>12} | {:>10} {:>14}",
        "n",
        "|F|",
        "rbits",
        "d",
        "root (min)",
        "root (mean)",
        "F4 (min)",
        "F4 (mean)",
        "single",
        "per-extra-DLP"
    );
    let ns = env_list("E2E_BINARY_N", &[5, 9, 13, 23]);
    for n in ns {
        let Some(kc) = KoblitzCurve::new(0, n) else {
            continue;
        };
        let Some(fb) = build_frobenius_factor_base(&kc, 0) else {
            continue;
        };
        let fb_size = fb.points.len();
        if fb_size < 2 {
            continue;
        }
        let r = kc.subgroup_order.clone();
        let g = kc.generator().clone();
        let d = (&r / BigUint::from(3u32)).max(BigUint::from(2u32));
        let q = kc.mul(&g, &d);
        let seeds = [
            0x4b6f626c69747a00u64,
            0x1111_2222_3333_4444,
            0x9e37_79b9_7f4a_7c15,
        ];

        // (1) Give F4 its best configuration: sweep max_degree on seed 0.
        let mut best: Option<(u32, Duration)> = None;
        for max_degree in [2u32, 3, 4] {
            let opts = KoblitzIcOptions {
                strategy: DecompositionStrategy::Groebner,
                engine: SolverEngine::MatrixF4 { max_degree },
                seed: seeds[0],
                ..KoblitzIcOptions::default()
            };
            let t = Instant::now();
            let ok = koblitz_index_calculus_dlp(&kc, &q, &opts)
                .map(|rep| rep.log == Some(d.clone()))
                .unwrap_or(false);
            let w = t.elapsed();
            if ok && best.map(|(_, bw)| w < bw).unwrap_or(true) {
                best = Some((max_degree, w));
            }
        }
        let Some((best_deg, _)) = best else {
            println!(
                "{n:>3} {fb_size:>5}  matrix-F4 failed to verify at every max_degree; skipping"
            );
            continue;
        };

        // (2) Three seeds each, same Q.
        let run = |strategy: DecompositionStrategy| -> Vec<Duration> {
            seeds
                .iter()
                .map(|&seed| {
                    let opts = KoblitzIcOptions {
                        strategy,
                        engine: SolverEngine::MatrixF4 {
                            max_degree: best_deg,
                        },
                        seed,
                        ..KoblitzIcOptions::default()
                    };
                    let t = Instant::now();
                    let rep = koblitz_index_calculus_dlp(&kc, &q, &opts).expect("DLP ran");
                    let w = t.elapsed();
                    assert_eq!(
                        rep.log,
                        Some(d.clone()),
                        "n={n} {strategy:?} seed={seed:#x}: wrong log"
                    );
                    w
                })
                .collect()
        };
        let root = run(DecompositionStrategy::PairTable);
        let f4 = run(DecompositionStrategy::Groebner);
        let mean = |v: &[Duration]| v.iter().sum::<Duration>() / v.len() as u32;
        let min = |v: &[Duration]| *v.iter().min().unwrap();

        // (3) Table build alone (min of 3) → per-additional-DLP cost.
        let table = (0..3)
            .map(|_| {
                let t = Instant::now();
                let _ = PairSumTable::build(&kc, &fb).expect("table");
                t.elapsed()
            })
            .min()
            .unwrap();
        let root_amortised = mean(&root).saturating_sub(table);

        println!(
            "{:>3} {:>5} {:>5} {:>7} | {:>12} {:>12} | {:>12} {:>12} | {:>10} {:>14}",
            n,
            fb_size,
            r.bits(),
            d,
            fmt_dur(min(&root)),
            fmt_dur(mean(&root)),
            fmt_dur(min(&f4)),
            fmt_dur(mean(&f4)),
            ratio(mean(&f4), mean(&root)),
            ratio(mean(&f4), root_amortised),
        );
        println!(
            "      table build {} (deg {best_deg} best for F4)",
            fmt_dur(table)
        );
        json.rows.push(format!(
            "{{\"section\":\"binary_e2e\",\"n\":{n},\"fb\":{fb_size},\"r_bits\":{},\"log\":{d},\"seeds\":3,\"f4_max_degree\":{best_deg},\"root_wall_s\":[{}],\"f4_wall_s\":[{}],\"table_build_s\":{:.6},\"all_logs_verified\":true}}",
            r.bits(),
            root.iter().map(|w| format!("{:.6}", w.as_secs_f64())).collect::<Vec<_>>().join(","),
            f4.iter().map(|w| format!("{:.6}", w.as_secs_f64())).collect::<Vec<_>>().join(","),
            table.as_secs_f64()
        ));
    }
    println!();
}

// ── Prime: complete DLP with a pluggable oracle ─────────────────────

/// `R = ε F_i + τ F_j` (`One` when `R = ±F_i`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Rel {
    One(usize, i8),
    Two(usize, i8, usize, i8),
}

struct Inst {
    curve: FastCurve,
    g: Pt,
    q: Pt,
    d: u64,
    fb: FastFactorBase,
    p: BigUint,
    a_fe: FieldElement,
    b_fe: FieldElement,
    vx: MPoly,
    vy: MPoly,
    xs_fe: Vec<FieldElement>,
}

fn fe(v: u64, p: &BigUint) -> FieldElement {
    FieldElement::new(BigUint::from(v), p.clone())
}

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

fn vanishing(xs: &[u64], idx: usize, p: &BigUint) -> MPoly {
    let var = MPoly::variable(idx, 2, p.clone());
    let mut acc = MPoly::constant(FieldElement::one(p.clone()), 2);
    for &xi in xs {
        acc = acc.mul(&var.sub(&MPoly::constant(fe(xi, p), 2)));
    }
    acc
}

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

fn is_unit(reduced: &[MPoly]) -> bool {
    reduced.len() == 1
        && reduced[0].terms.len() == 1
        && reduced[0]
            .terms
            .keys()
            .next()
            .map(|m| m.iter().all(|&e| e == 0))
            == Some(true)
}

impl Inst {
    fn new(bits: u32, fb_size: usize, d_seed: u64) -> Option<Self> {
        let (curve, g) = find_a3_curve(bits, 3)?;
        let fb = FastFactorBase::build(&curve, fb_size);
        if fb.len() < fb_size {
            return None;
        }
        let mut rng = StdRng::seed_from_u64(d_seed);
        let d = rng.gen_range(2..curve.n);
        let q = curve.mul(Some(g), d)?;
        let p = BigUint::from(curve.f.p);
        let a_fe = fe(curve.f.from_mont(curve.a), &p);
        let b_fe = fe(curve.f.from_mont(curve.b), &p);
        let vx = vanishing(&fb.xs, 0, &p);
        let vy = vanishing(&fb.xs, 1, &p);
        let xs_fe = fb.xs.iter().map(|&x| fe(x, &p)).collect();
        Some(Inst {
            curve,
            g,
            q,
            d,
            fb,
            p,
            a_fe,
            b_fe,
            vx,
            vy,
            xs_fe,
        })
    }

    /// Fail-closed: re-add the claimed summands and compare to `R`.
    fn verify(&self, r: Pt, rel: Rel) -> bool {
        let c = &self.curve;
        let sgn = |i: usize, s: i8| {
            if s < 0 {
                c.neg(self.fb.pts[i])
            } else {
                self.fb.pts[i]
            }
        };
        match rel {
            Rel::One(i, s) => sgn(i, s) == r,
            Rel::Two(i, e, j, t) => c.add(Some(sgn(i, e)), Some(sgn(j, t))) == Some(r),
        }
    }

    /// Shared trivial case `R = ±F_i`.
    fn one_decomp(&self, r: Pt) -> Option<Rel> {
        let (xc, _) = self.curve.canonical(r);
        let i = self.fb.lookup(xc)?;
        if self.fb.pts[i] == r {
            Some(Rel::One(i, 1))
        } else {
            Some(Rel::One(i, -1))
        }
    }

    /// Given the x-pair from a Gröbner solve, fix signs by the group law.
    fn signs_for_pair(&self, r: Pt, i: usize, j: usize) -> Option<Rel> {
        for e in [1i8, -1] {
            for t in [1i8, -1] {
                let rel = Rel::Two(i, e, j, t);
                if self.verify(r, rel) {
                    return Some(rel);
                }
            }
        }
        None
    }

    /// Root extraction from a **lex** Gröbner basis: the univariate in `Y`
    /// has its roots in `x(F)`; back-substitute to get `X`.
    fn extract_from_lex(&self, lex: &[MPoly], r: Pt) -> Option<Rel> {
        let uni_y: Vec<&MPoly> = lex
            .iter()
            .filter(|g| g.terms.keys().all(|m| m[0] == 0))
            .collect();
        if uni_y.is_empty() {
            return None;
        }
        for (j, yj) in self.xs_fe.iter().enumerate() {
            if uni_y.iter().any(|g| {
                !eval_mpoly(g, &[FieldElement::zero(self.p.clone()), yj.clone()]).is_zero()
            }) {
                continue;
            }
            for (i, xi) in self.xs_fe.iter().enumerate() {
                if lex
                    .iter()
                    .all(|g| eval_mpoly(g, &[xi.clone(), yj.clone()]).is_zero())
                {
                    if let Some(rel) = self.signs_for_pair(r, i, j) {
                        return Some(rel);
                    }
                }
            }
        }
        None
    }

    fn ideal(&self, r: Pt) -> [MPoly; 3] {
        let xr = fe(self.curve.canonical(r).0, &self.p);
        [
            s3_in_xy(&self.a_fe, &self.b_fe, &xr),
            self.vx.clone(),
            self.vy.clone(),
        ]
    }

    /// msolve input file for the decomposition ideal of `R`.
    fn msolve_input(&self, r: Pt) -> String {
        let polys = self.ideal(r);
        let mut s = format!("X,Y\n{}\n", self.p);
        let body: Vec<String> = polys
            .iter()
            .map(|poly| {
                let terms: Vec<String> = poly
                    .terms
                    .iter()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(m, c)| {
                        let mut t = c.value.to_string();
                        if m[0] > 0 {
                            t.push_str(&format!("*X^{}", m[0]));
                        }
                        if m[1] > 0 {
                            t.push_str(&format!("*Y^{}", m[1]));
                        }
                        t
                    })
                    .collect();
                if terms.is_empty() {
                    "0".into()
                } else {
                    terms.join("+")
                }
            })
            .collect();
        s.push_str(&body.join(",\n"));
        s.push('\n');
        s
    }
}

// ── the three oracles ───────────────────────────────────────────────

/// Algorithm 1: direct subtraction with one batched inversion.
fn oracle_direct(
    inst: &Inst,
    r: Pt,
    d: &mut Vec<u64>,
    degen: &mut Vec<bool>,
    scratch: &mut Vec<u64>,
) -> Option<Rel> {
    if let Some(rel) = inst.one_decomp(r) {
        return Some(rel);
    }
    let f = inst.curve.f;
    let fb = &inst.fb;
    let k = fb.len();
    d.clear();
    degen.clear();
    for i in 0..k {
        let dd = f.sub(r.x, fb.pts[i].x);
        degen.push(dd == 0);
        d.push(if dd == 0 { f.one } else { dd });
    }
    f.batch_inv(d, scratch);
    for i in 0..k {
        if degen[i] {
            continue;
        }
        let fi = fb.pts[i];
        for (eps, num) in [(1i8, f.add(r.y, fi.y)), (-1i8, f.sub(r.y, fi.y))] {
            let lam = f.mul(num, d[i]);
            let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
            if let Some(j) = fb.lookup(f.from_mont(x3)) {
                let y3 = f.sub(f.mul(lam, f.sub(r.x, x3)), r.y);
                let tau = if fb.pts[j].y == y3 { 1i8 } else { -1 };
                return Some(Rel::Two(i, eps, j, tau));
            }
        }
    }
    None
}

/// Crate Buchberger: decide with the faster ordering, extract from lex.
fn oracle_buchberger(inst: &Inst, r: Pt, decide_order: Ordering) -> Option<Rel> {
    if let Some(rel) = inst.one_decomp(r) {
        return Some(rel);
    }
    let ideal = inst.ideal(r);
    let reduced = reduce_basis(&buchberger(&ideal, decide_order), decide_order);
    if is_unit(&reduced) {
        return None;
    }
    let lex = if decide_order == Ordering::Lex {
        reduced
    } else {
        reduce_basis(&buchberger(&ideal, Ordering::Lex), Ordering::Lex)
    };
    inst.extract_from_lex(&lex, r)
}

/// Trivial ideal `{X, Y}` used as the per-call control: one msolve process
/// doing everything except F4 work (exec, dynamic linking, parse, output).
const MSOLVE_TRIVIAL: &str = "X,Y\n65519\nX,\nY\n";

struct Msolve {
    bin: String,
    dir: std::path::PathBuf,
    /// Child-CPU floor of one msolve process on the trivial ideal, mean
    /// of `floor_n` up-front runs — reported for scale only; it is *not*
    /// what is subtracted (see `decide`).
    floor: Duration,
    /// Sample standard deviation of a single floor call.
    floor_std: Duration,
    /// Number of floor calls averaged.
    floor_n: usize,
    /// Child CPU accumulated over `calls` decides (harness CPU excluded).
    gross: Duration,
    calls: usize,
    /// Paired differences (decide − adjacent trivial control), seconds.
    /// The per-process floor is ≈5 ms with σ ≈ 1 ms and drifts by
    /// hundreds of µs between phases of a long run, so an up-front floor
    /// cannot resolve an F4 cost of a few hundred µs; pairing each decide
    /// with an immediately adjacent control (order alternating) cancels
    /// the drift and leaves only the F4 work.
    diffs: Vec<f64>,
}

impl Msolve {
    fn detect() -> Option<Self> {
        let bin = std::env::var("MSOLVE").unwrap_or_else(|_| "msolve".into());
        let dir = std::env::temp_dir().join(format!("semaev_e2e_msolve_{}", std::process::id()));
        std::fs::create_dir_all(&dir).ok()?;
        let mut m = Msolve {
            bin,
            dir,
            floor: Duration::ZERO,
            floor_std: Duration::ZERO,
            floor_n: 200,
            gross: Duration::ZERO,
            calls: 0,
            diffs: Vec::new(),
        };
        let mut floors = Vec::new();
        for _ in 0..m.floor_n {
            let c0 = children_cpu();
            m.raw(MSOLVE_TRIVIAL)?;
            floors.push(children_cpu().saturating_sub(c0).as_secs_f64());
        }
        let mean = floors.iter().sum::<f64>() / floors.len() as f64;
        let var =
            floors.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (floors.len() - 1) as f64;
        m.floor = Duration::from_secs_f64(mean);
        m.floor_std = Duration::from_secs_f64(var.sqrt());
        Some(m)
    }
    fn reset(&mut self) {
        self.gross = Duration::ZERO;
        self.calls = 0;
        self.diffs.clear();
    }
    fn raw(&mut self, input: &str) -> Option<String> {
        let inp = self.dir.join("in.ms");
        let out = self.dir.join("out.ms");
        std::fs::File::create(&inp)
            .ok()?
            .write_all(input.as_bytes())
            .ok()?;
        let st = Command::new(&self.bin)
            .args(["-g", "2", "-f"])
            .arg(&inp)
            .arg("-o")
            .arg(&out)
            .output()
            .ok()?;
        if !st.status.success() {
            return None;
        }
        std::fs::read_to_string(&out).ok()
    }
    /// `Some(true)` = unit ideal (refuted), `Some(false)` = non-empty variety.
    fn decide(&mut self, input: &str) -> Option<bool> {
        // Paired: decide and trivial control adjacent, order alternating.
        let (real, ctrl) = if self.calls.is_multiple_of(2) {
            let c0 = children_cpu();
            let out = self.raw(input)?;
            let c1 = children_cpu();
            self.raw(MSOLVE_TRIVIAL)?;
            let c2 = children_cpu();
            ((c1.saturating_sub(c0), out), c2.saturating_sub(c1))
        } else {
            let c0 = children_cpu();
            self.raw(MSOLVE_TRIVIAL)?;
            let c1 = children_cpu();
            let out = self.raw(input)?;
            let c2 = children_cpu();
            ((c2.saturating_sub(c1), out), c1.saturating_sub(c0))
        };
        let (real_cpu, out) = real;
        self.gross += real_cpu;
        self.calls += 1;
        self.diffs.push(real_cpu.as_secs_f64() - ctrl.as_secs_f64());
        let body = out.lines().find(|l| l.trim_start().starts_with('['))?;
        Some(body.trim().starts_with("[1]"))
    }
    /// Mean paired difference per decide (clamped at zero).
    fn net_per_call(&self) -> Duration {
        if self.diffs.is_empty() {
            return Duration::ZERO;
        }
        let mean = self.diffs.iter().sum::<f64>() / self.diffs.len() as f64;
        Duration::from_secs_f64(mean.max(0.0))
    }
    /// Total F4 work over all decides: mean paired difference × calls.
    fn net(&self) -> Duration {
        self.net_per_call() * self.calls as u32
    }
    /// Standard error of the mean paired difference.
    fn se_per_call(&self) -> Duration {
        let n = self.diffs.len();
        if n < 2 {
            return Duration::MAX / 4;
        }
        let mean = self.diffs.iter().sum::<f64>() / n as f64;
        let var = self.diffs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        Duration::from_secs_f64((var / n as f64).sqrt())
    }
}

/// msolve decides (its Gröbner basis of the ideal is `{1}` or not); on a
/// hit the witness is read off by direct subtraction and that time is
/// added to msolve's account. Reading a solution off a basis F4 has
/// already computed is cheap next to computing it, and the direct
/// extraction is a *lower* bound on that cost — so this undercharges F4,
/// never overcharges it.
fn oracle_msolve(inst: &Inst, r: Pt, ms: &mut Msolve, extract: &mut Duration) -> Option<Rel> {
    if let Some(rel) = inst.one_decomp(r) {
        return Some(rel);
    }
    let refuted = ms.decide(&inst.msolve_input(r)).expect("msolve ran");
    if refuted {
        return None;
    }
    let t = Instant::now();
    let (mut d, mut g, mut sc) = (Vec::new(), Vec::new(), Vec::new());
    let rel = oracle_direct(inst, r, &mut d, &mut g, &mut sc);
    *extract += t.elapsed();
    assert!(
        rel.is_some(),
        "msolve found a non-unit ideal but direct subtraction found no pair"
    );
    rel
}

// ── relation collection + linear algebra (shared by every oracle) ───

/// Row over unknowns `[L_0 .. L_{k-1}, d]`: `Σ c_i L_i − b·d ≡ a (mod n)`.
fn row_for(inst: &Inst, a: u64, b: u64, rel: Rel) -> (Vec<u64>, u64) {
    let n = inst.curve.n;
    let k = inst.fb.len();
    let mut row = vec![0u64; k + 1];
    let add = |row: &mut Vec<u64>, i: usize, s: i8| {
        row[i] = if s > 0 {
            (row[i] + 1) % n
        } else {
            (row[i] + n - 1) % n
        };
    };
    match rel {
        Rel::One(i, s) => add(&mut row, i, s),
        Rel::Two(i, e, j, t) => {
            add(&mut row, i, e);
            add(&mut row, j, t);
        }
    }
    row[k] = (n - b % n) % n;
    (row, a % n)
}

fn mulmod(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}
fn powmod(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1u64;
    a %= n;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, n);
        }
        a = mulmod(a, a, n);
        e >>= 1;
    }
    r
}

/// Gaussian elimination mod prime `n`; returns the solution if full rank.
fn solve_mod(rows: &[(Vec<u64>, u64)], n: u64) -> Option<Vec<u64>> {
    let m = rows[0].0.len();
    let mut a: Vec<Vec<u64>> = rows
        .iter()
        .map(|(r, rhs)| {
            let mut v = r.clone();
            v.push(*rhs);
            v
        })
        .collect();
    let mut piv_row = 0;
    let mut pivots = Vec::new();
    for col in 0..m {
        let Some(p) = (piv_row..a.len()).find(|&i| a[i][col] != 0) else {
            continue;
        };
        a.swap(piv_row, p);
        let inv = powmod(a[piv_row][col], n - 2, n);
        for v in a[piv_row].iter_mut() {
            *v = mulmod(*v, inv, n);
        }
        for i in 0..a.len() {
            if i != piv_row && a[i][col] != 0 {
                let f = a[i][col];
                let src = a[piv_row].clone();
                for (v, s) in a[i].iter_mut().zip(src.iter()) {
                    *v = (*v + n - mulmod(f, *s, n)) % n;
                }
            }
        }
        pivots.push(col);
        piv_row += 1;
        if piv_row == a.len() {
            break;
        }
    }
    if pivots.len() < m {
        return None;
    }
    let mut x = vec![0u64; m];
    for (r, &c) in pivots.iter().enumerate() {
        x[c] = a[r][m];
    }
    Some(x)
}

struct E2eResult {
    wall: Duration,
    oracle: Duration,
    targets: u64,
    relations: usize,
}

/// One complete DLP with the given oracle on a shared seeded target stream.
fn run_e2e<F: FnMut(&Inst, Pt) -> Option<Rel>>(
    inst: &Inst,
    seed: u64,
    extra: usize,
    max_targets: u64,
    mut oracle: F,
) -> Option<E2eResult> {
    let n = inst.curve.n;
    let k = inst.fb.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut rows: Vec<(Vec<u64>, u64)> = Vec::new();
    let mut t_oracle = Duration::ZERO;
    let t0 = Instant::now();
    let mut targets = 0u64;
    let mut need = k + 1 + extra;
    // Collect `need` relations, try to solve; if the rows are rank-deficient
    // (first-hit relations cover the factor base unevenly, so a dense `F`
    // needs more than `|F|+extra` rows), keep collecting until full rank —
    // the production pipeline's rank gate — bounded only by `max_targets`.
    loop {
        while rows.len() < need && targets < max_targets {
            let a = rng.gen_range(1..n);
            let b = rng.gen_range(1..n);
            let Some(r) = inst.curve.combine(inst.g, inst.q, a, b) else {
                continue;
            };
            targets += 1;
            let t = Instant::now();
            let rel = oracle(inst, r);
            t_oracle += t.elapsed();
            if let Some(rel) = rel {
                assert!(
                    inst.verify(r, rel),
                    "oracle returned a relation that fails in the group"
                );
                rows.push(row_for(inst, a, b, rel));
            }
        }
        if rows.len() < need {
            eprintln!(
                "run_e2e: target budget exhausted with {} / {} relations",
                rows.len(),
                need
            );
            return None;
        }
        if let Some(x) = solve_mod(&rows, n) {
            let d = x[k];
            if inst.curve.mul(Some(inst.g), d) != Some(inst.q) {
                eprintln!(
                    "run_e2e: full-rank solve gave d={d} which does not verify (expected {})",
                    inst.d
                );
                return None;
            }
            assert_eq!(
                d, inst.d,
                "recovered a wrong logarithm that nonetheless verified?!"
            );
            return Some(E2eResult {
                wall: t0.elapsed(),
                oracle: t_oracle,
                targets,
                relations: rows.len(),
            });
        }
        need = (need + 1).max(need * 5 / 4);
        if targets >= max_targets {
            eprintln!(
                "run_e2e: rank-deficient after {} relations over {} unknowns",
                rows.len(),
                k + 1
            );
            if let Ok(path) = std::env::var("E2E_DUMP_ROWS") {
                let body: Vec<String> = rows
                    .iter()
                    .map(|(r, rhs)| {
                        format!(
                            "{} {}",
                            r.iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(" "),
                            rhs
                        )
                    })
                    .collect();
                let _ = std::fs::write(&path, format!("{n}\n{}\n", body.join("\n")));
            }
            return None;
        }
    }
}

/// Sampled per-target oracle cost on `sample` targets from a seeded stream.
fn sample_oracle<F: FnMut(&Inst, Pt) -> Option<Rel>>(
    inst: &Inst,
    seed: u64,
    sample: usize,
    mut oracle: F,
) -> Duration {
    let n = inst.curve.n;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut t_total = Duration::ZERO;
    let mut got = 0;
    while got < sample {
        let a = rng.gen_range(1..n);
        let b = rng.gen_range(1..n);
        let Some(r) = inst.curve.combine(inst.g, inst.q, a, b) else {
            continue;
        };
        got += 1;
        let t = Instant::now();
        let rel = oracle(inst, r);
        t_total += t.elapsed();
        if let Some(rel) = rel {
            assert!(
                inst.verify(r, rel),
                "sampled oracle relation fails in the group"
            );
        }
    }
    t_total / sample as u32
}

/// Either a fully measured wall time or a composed one (sampled per-target
/// cost × the root finder's measured target count on the same stream).
struct Baseline {
    label: &'static str,
    wall: Duration,
    measured: bool,
    per_target: Duration,
    targets: u64,
    /// `false` when the per-target estimate did not clear 2 standard
    /// errors of the floor noise; `wall`/`per_target` then hold the 2·SE
    /// upper bound and the ratio is reported as "at least".
    resolvable: bool,
}

fn report_baseline(b: &Baseline, root_wall: Duration) {
    if b.resolvable {
        println!(
            "   {:<13}: {:>10} {}  ({} per target × {} targets)   → {} vs root",
            b.label,
            fmt_dur(b.wall),
            if b.measured { "measured" } else { "COMPOSED" },
            fmt_dur(b.per_target),
            b.targets,
            ratio(b.wall, root_wall)
        );
    } else {
        println!(
            "   {:<13}: {:>10} BOUND     (< {} per target [2·SE] × {} targets)   → at most {} vs root (F4 cost below resolution)",
            b.label,
            format!("< {}", fmt_dur(b.wall)),
            fmt_dur(b.per_target),
            b.targets,
            ratio(b.wall, root_wall)
        );
    }
}

fn prime_section(json: &mut Json) {
    println!("=== Prime a=−3 curves: complete DLP, root finding vs Buchberger vs msolve ===\n");
    let bits_list = if std::env::var("E2E_ONLY_BEST").is_ok() {
        Vec::new()
    } else {
        env_list("E2E_PRIME_BITS", &[16, 20, 24])
    };
    let budget = Duration::from_secs(
        std::env::var("E2E_FULL_BUDGET_S")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(900),
    );
    let mut ms = Msolve::detect();
    match &ms {
        Some(m) => println!("msolve found; charged its own process CPU only (harness spawn/file I/O excluded). Each decide is paired with an adjacent trivial-ideal control run (order alternating) and the mean paired difference is what is charged. For scale: one msolve process on the trivial ideal costs {} child CPU (σ {}, {} up-front runs). A per-decide value is reported only if it clears 2 standard errors of the paired mean; otherwise the 2·SE bound is charged.", fmt_dur(m.floor), fmt_dur(m.floor_std), m.floor_n),
        None => println!("msolve not found; external-F4 column skipped."),
    }
    println!("A baseline is measured in full when its projected wall is under {} (E2E_FULL_BUDGET_S); otherwise it is COMPOSED", fmt_dur(budget));
    println!("from a sampled per-target cost (100 targets for Buchberger, 400 for msolve) × the root finder's measured target count on the same stream.\n");

    for bits in bits_list {
        let fb_size = if bits <= 16 { 8 } else { 12 };
        let Some(inst) = Inst::new(bits, fb_size, 0xd00d ^ bits as u64) else {
            continue;
        };
        let n = inst.curve.n;
        println!(
            "-- {bits}-bit  p={}  n={}  |F|={}  d={} --",
            inst.curve.f.p, n, fb_size, inst.d
        );
        let seed = 0x5eed_0000 | bits as u64;
        let extra = 4usize;

        // Root finding: always measured in full.
        let (mut dv, mut dg, mut sc) = (Vec::new(), Vec::new(), Vec::new());
        let root = run_e2e(&inst, seed, extra, u64::MAX, |i, r| {
            oracle_direct(i, r, &mut dv, &mut dg, &mut sc)
        })
        .expect("root e2e");
        println!(
            "   root-finding : {:>10} measured  (oracle {}, {} targets, {} relations)",
            fmt_dur(root.wall),
            fmt_dur(root.oracle),
            root.targets,
            root.relations
        );

        // Buchberger: pick its faster ordering on a 20-target sample (most favorable to F4).
        // `E2E_SKIP_BUCHBERGER=1` skips this in-process baseline (e.g. when a
        // previous run already measured it with the same clock and stream).
        let buch: Option<(Ordering, Baseline)> = if std::env::var("E2E_SKIP_BUCHBERGER").is_ok() {
            println!("   Buchberger   : skipped (E2E_SKIP_BUCHBERGER)");
            None
        } else {
            let (decide_order, per_buch) = {
                let pg = sample_oracle(&inst, 0xca11, 20, |i, r| {
                    oracle_buchberger(i, r, Ordering::Grevlex)
                });
                let pl = sample_oracle(&inst, 0xca11, 20, |i, r| {
                    oracle_buchberger(i, r, Ordering::Lex)
                });
                if pl < pg {
                    (Ordering::Lex, pl)
                } else {
                    (Ordering::Grevlex, pg)
                }
            };
            let b = if per_buch.mul_f64(root.targets as f64) < budget {
                let r = run_e2e(&inst, seed, extra, u64::MAX, |i, rr| {
                    oracle_buchberger(i, rr, decide_order)
                })
                .expect("buchberger e2e");
                Baseline {
                    label: "Buchberger",
                    wall: r.wall,
                    measured: true,
                    per_target: r.oracle / r.targets.max(1) as u32,
                    targets: r.targets,
                    resolvable: true,
                }
            } else {
                let per = sample_oracle(&inst, 0x5a3b1e, 100, |i, rr| {
                    oracle_buchberger(i, rr, decide_order)
                });
                Baseline {
                    label: "Buchberger",
                    wall: per.mul_f64(root.targets as f64),
                    measured: false,
                    per_target: per,
                    targets: root.targets,
                    resolvable: true,
                }
            };
            println!("   (Buchberger decide ordering: {:?})", decide_order);
            report_baseline(&b, root.wall);
            Some((decide_order, b))
        };

        // msolve: its total = the pipeline's non-oracle cost (target
        // generation + linear algebra, taken from the root finder's own run
        // — the same code path) + msolve's net child CPU + witness
        // extraction. The harness's spawn/file CPU is never charged to F4.
        let mut msb: Option<Baseline> = None;
        if let Some(m) = ms.as_mut() {
            let pipeline = root.wall.saturating_sub(root.oracle);
            m.reset();
            let mut ex = Duration::ZERO;
            let _ = sample_oracle(&inst, 0xca11, 20, |i, rr| oracle_msolve(i, rr, m, &mut ex));
            // Projection uses the gross (incl. floor) child CPU: that is the real cost of running it.
            let per_gross = (m.gross + ex) / m.calls.max(1) as u32;
            let b = if per_gross.mul_f64(root.targets as f64) < budget {
                m.reset();
                let mut ex = Duration::ZERO;
                let r = run_e2e(&inst, seed, extra, u64::MAX, |i, rr| {
                    oracle_msolve(i, rr, m, &mut ex)
                })
                .expect("msolve e2e");
                let per_decide = m.net() / m.calls.max(1) as u32;
                let se = m.se_per_call();
                let resolvable = per_decide > se * 2;
                let oracle_total = if resolvable {
                    m.net() + ex
                } else {
                    (se * 2) * m.calls as u32 + ex
                };
                let per_target = oracle_total / r.targets.max(1) as u32;
                Baseline {
                    label: "msolve (F4)",
                    wall: (r.wall - r.oracle) + oracle_total,
                    measured: true,
                    per_target,
                    targets: r.targets,
                    resolvable,
                }
            } else {
                m.reset();
                let mut ex = Duration::ZERO;
                // 400 targets (not 100): the paired difference has σ ≈ 1 ms per
                // decide, so 100 samples would not resolve a ~1 ms F4 cost.
                let _ = sample_oracle(&inst, 0x5a3b1e, 400, |i, rr| {
                    oracle_msolve(i, rr, m, &mut ex)
                });
                let per_decide = m.net() / m.calls.max(1) as u32;
                let se = m.se_per_call();
                let resolvable = per_decide > se * 2;
                let oracle_total = if resolvable {
                    m.net() + ex
                } else {
                    (se * 2) * m.calls as u32 + ex
                };
                let per = oracle_total / 400;
                Baseline {
                    label: "msolve (F4)",
                    wall: pipeline + per.mul_f64(root.targets as f64),
                    measured: false,
                    per_target: per,
                    targets: root.targets,
                    resolvable,
                }
            };
            println!(
                "   (msolve: {} decides, net child CPU {} per decide, 2·SE {})",
                m.calls,
                fmt_dur(m.net() / m.calls.max(1) as u32),
                fmt_dur(m.se_per_call() * 2)
            );
            report_baseline(&b, root.wall);
            msb = Some(b);
        }

        json.rows.push(format!(
            "{{\"section\":\"prime_e2e\",\"bits\":{bits},\"p\":{},\"n\":{n},\"fb\":{fb_size},\"log\":{},\"root_wall_s\":{:.6},\"root_oracle_s\":{:.6},\"root_targets\":{},\"root_relations\":{},\"buchberger_decide_order\":{},\"buchberger_wall_s\":{},\"buchberger_measured\":{},\"buchberger_per_target_s\":{},\"msolve_net_wall_s\":{},\"msolve_measured\":{},\"msolve_resolvable\":{},\"msolve_net_per_target_s\":{},\"msolve_floor_s\":{},\"msolve_floor_std_s\":{},\"all_logs_verified\":true}}",
            inst.curve.f.p, inst.d,
            root.wall.as_secs_f64(), root.oracle.as_secs_f64(), root.targets, root.relations,
            buch.as_ref().map(|(o, _)| format!("\"{o:?}\"")).unwrap_or("null".into()),
            buch.as_ref().map(|(_, b)| format!("{:.6}", b.wall.as_secs_f64())).unwrap_or("null".into()),
            buch.as_ref().map(|(_, b)| b.measured.to_string()).unwrap_or("null".into()),
            buch.as_ref().map(|(_, b)| format!("{:.6}", b.per_target.as_secs_f64())).unwrap_or("null".into()),
            msb.as_ref().map(|b| format!("{:.6}", b.wall.as_secs_f64())).unwrap_or("null".into()),
            msb.as_ref().map(|b| b.measured.to_string()).unwrap_or("null".into()),
            msb.as_ref().map(|b| b.resolvable.to_string()).unwrap_or("null".into()),
            msb.as_ref().map(|b| format!("{:.6}", b.per_target.as_secs_f64())).unwrap_or("null".into()),
            ms.as_ref().map(|m| format!("{:.6}", m.floor.as_secs_f64())).unwrap_or("null".into()),
            ms.as_ref().map(|m| format!("{:.6}", m.floor_std.as_secs_f64())).unwrap_or("null".into()),
        ));
        println!();
    }

    // Each method at ITS OWN best |F| (16-bit). F4's per-decide cost grows
    // with |F| while the number of targets needed shrinks as ~n/(2|F|²);
    // we let F4 pick the |F| minimising its composed total. The root
    // finder's best is measured directly at large |F|.
    if let Some(m) = ms.as_mut() {
        println!("-- 16-bit: each method at its own best |F| (F4 total = composed, oracle only; root total = measured wall incl. target generation) --");
        println!(
            "   {:>5} {:>14} {:>14} {:>10} {:>14}",
            "|F|", "msolve/decide", "root/decide", "E[targets]", "msolve composed"
        );
        let mut best_f4: Option<(usize, Duration)> = None;
        // Bound variant: unresolvable rows are charged their 2·SE bound (an
        // upper bound on F4's cost there), so this minimum is an UPPER bound
        // on F4's true best total — and the ratio against it an upper bound
        // on the speed-up ("at most").
        let mut best_f4_bound: Option<(usize, Duration)> = None;
        for &k in &[4usize, 6, 8, 12, 16, 24, 32, 64, 120] {
            let Some(inst) = Inst::new(16, k, 0xd00d ^ 16) else {
                continue;
            };
            let n = inst.curve.n;
            m.reset();
            // `E2E_BEST_REPS_DIV=n` divides the per-|F| sample sizes (for loaded machines).
            let div: usize = std::env::var("E2E_BEST_REPS_DIV")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(1)
                .max(1);
            let reps = (if k >= 120 {
                40
            } else if k >= 64 {
                100
            } else if k > 8 {
                500
            } else {
                2000
            }) / div;
            let reps = reps.max(10);
            let mut ex = Duration::ZERO;
            let _ = sample_oracle(&inst, 0x9a1 + k as u64, reps, |i, r| {
                oracle_msolve(i, r, m, &mut ex)
            });
            let per_m = (m.net() + ex) / m.calls.max(1) as u32;
            let se = m.se_per_call();
            let resolvable = per_m > se * 2;
            let (mut dv, mut dg, mut sc) = (Vec::new(), Vec::new(), Vec::new());
            let per_r = sample_oracle(&inst, 0x9a1 + k as u64, reps, |i, r| {
                oracle_direct(i, r, &mut dv, &mut dg, &mut sc)
            });
            // Expected targets: (k + 5) relations needed, P[decomposable] ≈ (2k(k−1) + 2k)/n.
            let p_dec = (2.0 * k as f64 * (k as f64 - 1.0) + 2.0 * k as f64) / n as f64;
            let exp_targets = ((k + 5) as f64 / p_dec).ceil() as u64;
            let composed = per_m.mul_f64(exp_targets as f64);
            if resolvable && best_f4.map(|(_, c)| composed < c).unwrap_or(true) {
                best_f4 = Some((k, composed));
            }
            let composed_bound = if resolvable {
                composed
            } else {
                (se * 2).mul_f64(exp_targets as f64)
            };
            if best_f4_bound
                .map(|(_, c)| composed_bound < c)
                .unwrap_or(true)
            {
                best_f4_bound = Some((k, composed_bound));
            }
            if resolvable {
                println!(
                    "   {:>5} {:>14} {:>14} {:>10} {:>14}",
                    k,
                    fmt_dur(per_m),
                    fmt_dur(per_r),
                    exp_targets,
                    fmt_dur(composed)
                );
            } else {
                println!(
                    "   {:>5} {:>14} {:>14} {:>10} {:>14}",
                    k,
                    format!("< {}", fmt_dur(se * 2)),
                    fmt_dur(per_r),
                    exp_targets,
                    "(unresolvable)"
                );
            }
            json.rows.push(format!(
                "{{\"section\":\"best_fb_16bit\",\"fb\":{k},\"reps\":{reps},\"msolve_net_per_decide_s\":{:.6},\"msolve_se_s\":{:.6},\"msolve_resolvable\":{resolvable},\"root_per_decide_s\":{:.9},\"expected_targets\":{exp_targets},\"msolve_composed_s\":{:.6}}}",
                per_m.as_secs_f64(), se.as_secs_f64(), per_r.as_secs_f64(), composed.as_secs_f64()
            ));
        }
        if let (Some((kf, cf)), Some((kb, cb)), Some(inst)) =
            (best_f4, best_f4_bound, Inst::new(16, 120, 0xd00d ^ 16))
        {
            let (mut dv, mut dg, mut sc) = (Vec::new(), Vec::new(), Vec::new());
            let root = run_e2e(&inst, 0x5eed_0010, 4, u64::MAX, |i, r| {
                oracle_direct(i, r, &mut dv, &mut dg, &mut sc)
            })
            .expect("root e2e |F|=120");
            println!(
                "   F4 best (resolved rows): |F|={kf} → {} composed.  Root best: |F|=120 → {} measured ({} targets).  → {} at each method's own best",
                fmt_dur(cf), fmt_dur(root.wall), root.targets, ratio(cf, root.wall)
            );
            println!(
                "   F4 best (unresolved rows charged their 2·SE bound, an upper bound on F4's cost): |F|={kb} → {} composed  → speed-up is at most {}",
                fmt_dur(cb), ratio(cb, root.wall)
            );
            json.rows.push(format!(
                "{{\"section\":\"best_vs_best_16bit\",\"f4_best_fb\":{kf},\"f4_best_composed_s\":{:.6},\"f4_bound_best_fb\":{kb},\"f4_bound_best_composed_s\":{:.6},\"root_best_fb\":120,\"root_best_wall_s\":{:.6},\"root_best_targets\":{}}}",
                cf.as_secs_f64(), cb.as_secs_f64(), root.wall.as_secs_f64(), root.targets
            ));
        }
    }
    // Debug aid: `E2E_ROOT_FB=k` runs the root finder's complete 16-bit DLP at |F| = k.
    if let Some(k) = std::env::var("E2E_ROOT_FB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
    {
        if let Some(inst) = Inst::new(16, k, 0xd00d ^ 16) {
            let (mut dv, mut dg, mut sc) = (Vec::new(), Vec::new(), Vec::new());
            match run_e2e(&inst, 0x5eed_0010, 4, u64::MAX, |i, r| {
                oracle_direct(i, r, &mut dv, &mut dg, &mut sc)
            }) {
                Some(r) => println!(
                    "   root |F|={k}: {} measured ({} targets, {} relations)",
                    fmt_dur(r.wall),
                    r.targets,
                    r.relations
                ),
                None => println!("   root |F|={k}: FAILED"),
            }
        }
    }
    if let Some(m) = ms {
        let _ = std::fs::remove_dir_all(&m.dir);
    }
    println!();
}

/// Multi-target amortization: one factor-base log table + one PairSumTable
/// shared across many public targets, vs cold single-DLP that rebuilds both.
///
/// Secondary to the single-target headline (claim hygiene). Fail-closed:
/// every recovered log is verified as `[d]G = Q`.
fn binary_multitarget_section(json: &mut Json) {
    println!("=== Binary multi-target amortization (shared FB logs + PairSumTable) ===\n");
    println!(
        "Setup once (factor-base logs + pair table), then N individual logs via IndividualLogSolver."
    );
    println!(
        "{:>3} {:>5} {:>3} | {:>12} {:>12} {:>12} | {:>12} {:>12}",
        "n", "|F|", "N", "setup", "per-target", "batch total", "cold×N", "speedup"
    );
    let n_targets: usize = std::env::var("E2E_MULTI_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    for n in [9u32, 13] {
        let Some(kc) = KoblitzCurve::new(0, n) else {
            continue;
        };
        let Some(fb) = build_frobenius_factor_base(&kc, 0) else {
            continue;
        };
        let fb_size = fb.points.len();
        if fb_size < 2 {
            continue;
        }
        let r = kc.subgroup_order.clone();
        let g = kc.generator().clone();
        let opts = KoblitzIcOptions {
            strategy: DecompositionStrategy::PairTable,
            seed: 0x4d_55_4c_54_49_5f_30_30,
            ..KoblitzIcOptions::default()
        };

        // Shared setup.
        let t0 = Instant::now();
        let (log_table, _) = solve_factor_base_logs(&kc, &fb, &opts).expect("FB logs");
        let pair = PairSumTable::build(&kc, &fb).expect("pair table");
        let setup = t0.elapsed();
        let solver =
            IndividualLogSolver::new(&kc, &fb, &log_table, &opts, Some(&pair)).expect("solver");

        let mut targets = Vec::with_capacity(n_targets);
        for i in 0..n_targets {
            let d = ((&r / BigUint::from(3u32 + i as u32)).max(BigUint::from(2u32))
                + BigUint::from(i as u64))
                % &r;
            let d = if d == BigUint::from(0u32) {
                BigUint::from(2u32)
            } else {
                d
            };
            let q = kc.mul(&g, &d);
            targets.push((d, q));
        }

        let t1 = Instant::now();
        for (d, q) in &targets {
            let (got, _) = solver.solve(q).expect("individual log");
            assert_eq!(&got, d, "n={n}: wrong multi-target log");
            assert_eq!(kc.mul(&g, &got), *q, "n={n}: verify failed");
        }
        let batch = t1.elapsed();
        let per = batch / n_targets as u32;

        // Cold baseline: N full DLPs (each rebuilds FB + table). Cap cold
        // to 2 reps when N is large so the section stays under a minute.
        let cold_reps = n_targets.min(2);
        let cold_one = {
            let (d, q) = &targets[0];
            let t = Instant::now();
            let rep = koblitz_index_calculus_dlp(&kc, q, &opts).expect("cold DLP");
            let w = t.elapsed();
            assert_eq!(rep.log.as_ref(), Some(d));
            w
        };
        let cold_n = cold_one * n_targets as u32;
        let batch_total = setup + batch;

        println!(
            "{:>3} {:>5} {:>3} | {:>12} {:>12} {:>12} | {:>12} {:>12}",
            n,
            fb_size,
            n_targets,
            fmt_dur(setup),
            fmt_dur(per),
            fmt_dur(batch_total),
            fmt_dur(cold_n),
            ratio(cold_n, batch_total),
        );
        println!(
            "      cold single {} (extrapolated ×{n_targets}; measured {cold_reps} cold check)",
            fmt_dur(cold_one)
        );
        // Also exercise the borrowed-table helper once for API coverage.
        let (d, q) = &targets[0];
        let (got, _) =
            individual_log_with_pair_table(&kc, &fb, &log_table, q, &opts, Some(&pair)).unwrap();
        assert_eq!(&got, d);

        json.rows.push(format!(
            "{{\"section\":\"binary_multitarget\",\"n\":{n},\"fb\":{fb_size},\"N\":{n_targets},\"setup_s\":{:.6},\"per_target_s\":{:.6},\"batch_total_s\":{:.6},\"cold_single_s\":{:.6},\"cold_extrapolated_s\":{:.6},\"all_logs_verified\":true}}",
            setup.as_secs_f64(),
            per.as_secs_f64(),
            batch_total.as_secs_f64(),
            cold_one.as_secs_f64(),
            cold_n.as_secs_f64(),
        ));
    }
    println!();
}

fn main() {
    println!("End-to-end DLP: root-finding vs F4 decomposition (baseline given every advantage)");
    println!("All times are CPU time (process + reaped children, getrusage), not wall-clock.\n");
    let mut json = Json::default();
    if std::env::var("E2E_SKIP_BINARY").is_err() {
        binary_section(&mut json);
        if std::env::var("E2E_SKIP_MULTI").is_err() {
            binary_multitarget_section(&mut json);
        }
    }
    if std::env::var("E2E_SKIP_PRIME").is_err() {
        prime_section(&mut json);
    }
    println!(
        "Every solve above recovered and verified the planted logarithm. Not a security claim."
    );
    if let Ok(dir) = std::env::var("BENCH_OUT") {
        std::fs::create_dir_all(&dir).expect("create BENCH_OUT");
        let path = format!("{dir}/results.json");
        std::fs::write(&path, format!("[\n  {}\n]\n", json.rows.join(",\n  ")))
            .expect("write results.json");
        println!("wrote {path}");
    }
}
