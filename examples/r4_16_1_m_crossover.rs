//! **Round 4, lead 16.1: where kernel saturation pays — the feasibility
//! crossover in `m`.**
//!
//! `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §16.1.  The affine saturation law
//! (Prop. 15.2) charges a base closed under the rational 2-torsion kernel
//! `M = ⟨T⟩` a yield factor `|M|^{m−1}` and merges its columns by `|M|`, so
//! the mitm cost factor is `|M|^{m−2}`.  The `u`-frame (`u = 1/(x + 1)`,
//! `F_u = {u ∈ V}`, closed under `+T` because `1 ∈ V`) buys symmetrised
//! summation polynomials in `w = u² + u`, `s = Σu` of lower degree.  This
//! driver asks, on one shared `V ∋ 1` per `(n, m)` and one shared target
//! list per `n`, which arm the repository's own Boolean engine can solve
//! under one counted cap, and what a relation costs.
//!
//! Arms (all through the library's oracle functions, engine
//! `InheritedF4 { max_degree: 3 }` with the default split rule, the same
//! `node_budget` for every solve):
//!
//! | m | plain `x`-frame over `F_x = {x ∈ V}` | kernel-saturated `u`-frame over `F_u = {u ∈ V}` |
//! |---|---|---|
//! | 2 | `x-direct`: Semaev `S₃` | `sym`: symmetrised `S₃` |
//! | 3 | `x-direct`: `S₄` (degree 6), `x-chained`: two `S₃` links | `sym`: symmetrised `S₄` (degree 4) |
//! | 4 | `x-chained`: three `S₃` links, `4ℓ + 2n` unknowns | `sym-chained`: three symmetrised `S₃` links, `4(ℓ−1) + 1 + 2n` unknowns |
//!
//! Every verdict is gated against exact enumeration on the arm's own base
//! (a `found` relation must re-sum to `R`, or to `R` with `T` appended; a
//! `refuted` target must have no decomposition), so the hit rate and the
//! completeness are measured, not assumed.
//!
//! Unit: the engine's counted **word XORs** (`f4_word_ops_thread` deltas:
//! elimination, specialisation and linear elimination; the Macaulay build
//! and the symbolic substitutions are not counted, so every oracle cost is a
//! lower bound).  A per-cell wall guard exists only as a safety stop and is
//! reported if it fires.  The conversion to group-addition equivalents (GAE)
//! uses the raw word-XOR price of `ic_boundary::calibrate_word_xor` over one
//! single-word curve addition, both measured in this process on a contended
//! host.
//!
//! ```text
//! cargo run --release --example r4_16_1_m_crossover -- \
//!     --curves 17:1,19:0,23:0 --ms 2,3,4 --targets 24 --node-budget 4096 \
//!     --json research/isogeny_fibers_20261010/round4/16_1_m_crossover/panel_b4096.json
//! ```

use std::collections::HashSet;
use std::env;
use std::fs;
use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::ic_boundary::calibrate_word_xor;
use crypto_lib::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, FieldStructure, SolverEngine, MAX_VARS,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{groebner_decompose, KoblitzCurve};
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    build_symmetrised_factor_base, build_symmetrised_factor_base_from_basis,
    chained_symmetrised_groebner_decompose, direct_x_groebner_decompose,
    plain_subspace_factor_base, random_subspace_containing_one, symmetrised_groebner_decompose,
    OracleOutcome,
};

/// Engine knobs that change what the solver builds or counts; the run
/// refuses to start with any of them set, so every arm runs the default.
const ENGINE_ENV: &[&str] = &[
    "IC_REDUCTION_CACHE",
    "KIC_F4_INHERIT",
    "SOLVER_SPLIT_RULE",
    "KIC_F4_MAX_DEGREE_ONLY",
    "KIC_F4_CLOSURE_ROUNDS",
    "KIC_F4_SOLVER_FULL_READBACK",
    "F4_F2_MAX_ROWS",
    "F4_F2_MAX_COLS",
    "KIC_LINEAR_ELIM",
];

struct Opts {
    curves: Vec<(u32, u8)>,
    ms: Vec<usize>,
    targets: usize,
    node_budget: usize,
    seed: u64,
    cell_wall_s: f64,
    json: Option<String>,
    /// Extra cells `(n, m, spec)`: `spec` is a dimension `ℓ` (random
    /// non-invariant `V ∋ 1`) or `div=i.j` (the Frobenius-invariant `V` of
    /// that divisor of `xⁿ − 1`, as the `ic` bases use).
    extra: Vec<(u32, usize, VSpec)>,
    /// Arms not run, as `(arm, m)`: recorded in the JSON as skipped.
    skip: Vec<(String, usize)>,
}

#[derive(Clone, Debug)]
enum VSpec {
    /// `ℓ = ⌈(n + log₂ m!)/m⌉`, the repository's E1 convention.
    Rule,
    Dim(usize),
    Div(Vec<usize>),
}

fn rule_ell(n: u32, m: usize) -> usize {
    let fact: f64 = (1..=m).map(|k| k as f64).product();
    ((n as f64 + fact.log2()) / m as f64).ceil() as usize
}

fn parse_args() -> Opts {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut o = Opts {
        curves: vec![(17, 1), (19, 0), (23, 0)],
        ms: vec![2, 3, 4],
        targets: 24,
        node_budget: 4096,
        seed: 0x16_1C_0550,
        cell_wall_s: 600.0,
        json: None,
        extra: Vec::new(),
        skip: Vec::new(),
    };
    let mut i = 0;
    while i < args.len() {
        let v = args.get(i + 1).cloned().unwrap_or_default();
        match args[i].as_str() {
            "--curves" => {
                o.curves = v
                    .split(',')
                    .map(|c| {
                        let (n, a) = c.split_once(':').expect("--curves n:a,…");
                        (n.parse().expect("n"), a.parse().expect("a"))
                    })
                    .collect()
            }
            "--ms" => o.ms = v.split(',').map(|x| x.parse().expect("--ms")).collect(),
            "--targets" => o.targets = v.parse().expect("--targets"),
            "--node-budget" => o.node_budget = v.parse().expect("--node-budget"),
            "--seed" => o.seed = v.parse().expect("--seed"),
            "--cell-wall-s" => o.cell_wall_s = v.parse().expect("--cell-wall-s"),
            "--json" => o.json = Some(v),
            "--skip" => {
                for c in v.split(',') {
                    let (arm, m) = c.split_once('@').expect("--skip arm@m");
                    o.skip.push((arm.to_string(), m.parse().expect("m")));
                }
            }
            "--v-spec" => {
                for c in v.split(',') {
                    let parts: Vec<&str> = c.split(':').collect();
                    assert_eq!(parts.len(), 3, "--v-spec n:m:ell or n:m:div=i.j");
                    let spec = if let Some(d) = parts[2].strip_prefix("div=") {
                        VSpec::Div(d.split('.').map(|x| x.parse().expect("div")).collect())
                    } else {
                        VSpec::Dim(parts[2].parse().expect("ell"))
                    };
                    o.extra.push((
                        parts[0].parse().expect("n"),
                        parts[1].parse().expect("m"),
                        spec,
                    ));
                }
            }
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    o
}

fn bits(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

/// Exact decomposition check over a signed base (closed under negation):
/// is `target` a sum of exactly `m` base points (with repetition)?
struct Exact {
    fc: FastCurve,
    points: Vec<FastPoint>,
    set: HashSet<u64>,
    pairs: HashSet<u64>,
}

impl Exact {
    fn new(fc: &FastCurve, pts: &[BinaryPoint], m: usize) -> Self {
        let points: Vec<FastPoint> = pts.iter().map(|p| fc.lift(p)).collect();
        let set: HashSet<u64> = points.iter().map(|p| p.pack()).collect();
        let mut pairs = HashSet::new();
        if m >= 3 {
            for i in 0..points.len() {
                for j in i..points.len() {
                    pairs.insert(fc.add(points[i], points[j]).pack());
                }
            }
        }
        Self {
            fc: fc.clone(),
            points,
            set,
            pairs,
        }
    }

    fn decomposes(&self, target: FastPoint, m: usize) -> bool {
        let fc = &self.fc;
        match m {
            2 => self
                .points
                .iter()
                .any(|&p| self.set.contains(&fc.add(target, fc.neg(p)).pack())),
            3 => self
                .points
                .iter()
                .any(|&p| self.pairs.contains(&fc.add(target, fc.neg(p)).pack())),
            4 => {
                for i in 0..self.points.len() {
                    let a = fc.add(target, fc.neg(self.points[i]));
                    for j in i..self.points.len() {
                        let b = fc.add(a, fc.neg(self.points[j]));
                        if self.pairs.contains(&b.pack()) {
                            return true;
                        }
                    }
                }
                false
            }
            _ => unreachable!(),
        }
    }

    /// First-hit enumeration steps (one curve addition and one lookup each)
    /// to decide one target: the combinatorial reference on the same base.
    fn first_hit_steps(&self, target: FastPoint, m: usize) -> (u64, bool) {
        let fc = &self.fc;
        let mut steps = 0u64;
        match m {
            2 => {
                for &p in &self.points {
                    steps += 1;
                    if self.set.contains(&fc.add(target, fc.neg(p)).pack()) {
                        return (steps, true);
                    }
                }
            }
            3 => {
                for i in 0..self.points.len() {
                    let a = fc.add(target, fc.neg(self.points[i]));
                    for j in i..self.points.len() {
                        steps += 1;
                        let b = fc.add(a, fc.neg(self.points[j]));
                        if self.set.contains(&b.pack()) {
                            return (steps, true);
                        }
                    }
                }
            }
            _ => return (0, false),
        }
        (steps, false)
    }
}

#[derive(Default)]
struct CellAcc {
    attempts: u64,
    found: u64,
    refuted: u64,
    budget: u64,
    unbuildable: u64,
    truth: u64,
    found_on_truth: u64,
    gate_failures: u64,
    via_t: u64,
    word_xors: u64,
    word_xors_found: u64,
    splits: u64,
    wall_ns: u128,
    n_vars: usize,
    n_equations: usize,
    degree: u32,
    built_degree: u32,
    oversize_targets: u64,
    wall_truncated: bool,
    per_attempt: Vec<Value>,
}

fn measure_ns_per_add(fc: &FastCurve, pts: &[FastPoint], reps: usize) -> f64 {
    let mut acc = pts[0];
    let t0 = Instant::now();
    for k in 0..reps {
        acc = fc.add(acc, pts[k % pts.len()]);
    }
    let ns = t0.elapsed().as_nanos() as f64 / reps as f64;
    std::hint::black_box(acc);
    ns
}

/// Number of unknowns an arm's system would have, for the layout check
/// that the library applies (`MAX_VARS`).
fn layout_vars(arm: &str, m: usize, ell: usize, n: usize) -> usize {
    match arm {
        "x-direct" => m * ell,
        "x-chained" => m * ell + (m - 2) * n,
        "sym" => m * (ell - 1) + 1,
        "sym-chained" => 4 * (ell - 1) + 1 + 2 * n,
        _ => 0,
    }
}

fn arms_for(m: usize) -> Vec<&'static str> {
    match m {
        2 => vec!["x-direct", "sym"],
        3 => vec!["x-direct", "x-chained", "sym"],
        4 => vec!["x-chained", "sym-chained"],
        _ => vec![],
    }
}

fn main() {
    for var in ENGINE_ENV {
        assert!(
            env::var(var).is_err(),
            "the run needs {var} unset (default engine)"
        );
    }
    let o = parse_args();
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    let ns_per_word_xor = calibrate_word_xor();
    let mut out_curves = Vec::new();

    for &(n, a) in &o.curves {
        let kc = KoblitzCurve::new(a, n).expect("Koblitz curve with a usable subgroup");
        let fc = FastCurve::new(&kc.curve).expect("single-word field");
        let st = FieldStructure::new(n, &kc.curve.irreducible);
        let r: u64 = kc.subgroup_order.to_u64_digits()[0];
        let h: u64 = kc.cofactor.to_u64_digits().first().copied().unwrap_or(1);
        let irr_bits: Vec<u32> = kc.curve.irreducible.low_terms.clone();
        let mut rng = StdRng::seed_from_u64(o.seed ^ ((n as u64) << 32) ^ a as u64);
        // One target list per curve, shared by every m and every arm:
        // uniform in the prime-order subgroup, u(R) finite and non-zero.
        let g = kc.generator().clone();
        let mut targets: Vec<BinaryPoint> = Vec::new();
        while targets.len() < o.targets {
            let p = kc.mul(&g, &BigUint::from(rng.gen_range(1..r)));
            if let BinaryPoint::Affine { x, .. } = &p {
                if *x != F2mElement::one(n) {
                    targets.push(p);
                }
            }
        }
        let fast_targets: Vec<FastPoint> = targets.iter().map(|t| fc.lift(t)).collect();
        let mut sample: Vec<FastPoint> = fast_targets.clone();
        sample.push(fc.lift(&g));
        let ns_per_add = measure_ns_per_add(&fc, &sample, 2_000_000);
        let gae_per_xor = ns_per_word_xor / ns_per_add;
        eprintln!(
            "n={n} a={a} r={r} h={h}  ns/XOR {ns_per_word_xor:.3} ns/add {ns_per_add:.1}  GAE/XOR {gae_per_xor:.5}"
        );

        let mut cells = Vec::new();
        let mut plan: Vec<(usize, VSpec)> = o.ms.iter().map(|&m| (m, VSpec::Rule)).collect();
        plan.extend(
            o.extra
                .iter()
                .filter(|e| e.0 == n)
                .map(|e| (e.1, e.2.clone())),
        );
        for (m, spec) in plan {
            // One V ∋ 1 per (n, m), shared by both frames.
            let mut vrng = StdRng::seed_from_u64(o.seed ^ ((n as u64) << 32) ^ ((m as u64) << 8));
            let (v, invariant) = match &spec {
                VSpec::Rule => (
                    random_subspace_containing_one(&kc, rule_ell(n, m), &mut vrng),
                    false,
                ),
                VSpec::Dim(l) => (random_subspace_containing_one(&kc, *l, &mut vrng), false),
                VSpec::Div(idx) => (
                    build_symmetrised_factor_base(&kc, idx)
                        .expect("invariant V ∋ 1 for that divisor")
                        .v_basis,
                    true,
                ),
            };
            let ell = v.len();
            let fb_x = plain_subspace_factor_base(&kc, &v).expect("F_x");
            let index_x = fb_x.index_map();
            let fb_u = build_symmetrised_factor_base_from_basis(&kc, v.clone()).expect("F_u");
            let ex_x = Exact::new(&fc, &fb_x.points, m);
            let ex_u = Exact::new(&fc, &fb_u.points, m);
            let t_fast = fc.lift(&fb_u.two_torsion);
            // Columns of the relation matrix: signed classes of [h]P.  On F_x
            // (no Frobenius orbits on a non-invariant V) that is one per ±P;
            // on F_u, P, −P, P + T, −P + T share one (Prop. 4.1), T itself none.
            let cols_x = column_count(&fc, &ex_x.points, h);
            let cols_u = column_count(&fc, &ex_u.points, h);
            // Truth per target on each base (R + T counts on F_u: the oracle
            // may return a relation through R + T with T appended).
            let truth_x: Vec<bool> = fast_targets
                .iter()
                .map(|&t| ex_x.decomposes(t, m))
                .collect();
            let truth_u: Vec<bool> = fast_targets
                .iter()
                .map(|&t| ex_u.decomposes(t, m) || ex_u.decomposes(fc.add(t, t_fast), m))
                .collect();
            // Combinatorial reference (first-hit enumeration) on each base.
            let enum_ref = |ex: &Exact| -> Value {
                if m > 3 {
                    return Value::Null;
                }
                let (mut steps, mut hits) = (0u64, 0u64);
                for &t in &fast_targets {
                    let (s, hit) = ex.first_hit_steps(t, m);
                    steps += s;
                    hits += hit as u64;
                }
                json!({
                    "steps": steps, "hits": hits,
                    "gae_per_relation": if hits > 0 { steps as f64 / hits as f64 } else { f64::INFINITY },
                })
            };
            let enum_x = enum_ref(&ex_x);
            let enum_u = enum_ref(&ex_u);

            for arm in arms_for(m) {
                if o.skip.iter().any(|(a, sm)| a == arm && *sm == m) {
                    eprintln!("  n={n} m={m} ell={ell} {arm}: skipped (--skip)");
                    cells.push(json!({
                        "m": m, "ell": ell, "arm": arm, "v_invariant": invariant, "skipped": true,
                    }));
                    continue;
                }
                let on_u = arm.starts_with("sym");
                let points = if on_u { &fb_u.points } else { &fb_x.points };
                let truth = if on_u { &truth_u } else { &truth_x };
                let cols = if on_u { cols_u } else { cols_x };
                let vars = layout_vars(arm, m, ell, n as usize);
                let mut acc = CellAcc::default();
                let cell_start = Instant::now();
                if vars > MAX_VARS {
                    eprintln!(
                        "  n={n} m={m} {arm}: {vars} unknowns > {MAX_VARS}: the engine's layout refuses it"
                    );
                    acc.unbuildable = targets.len() as u64;
                    acc.attempts = targets.len() as u64;
                    acc.n_vars = vars;
                } else {
                    for (ti, target) in targets.iter().enumerate() {
                        if cell_start.elapsed().as_secs_f64() > o.cell_wall_s {
                            acc.wall_truncated = true;
                            break;
                        }
                        let before = f4_word_ops_thread();
                        let t0 = Instant::now();
                        let outcome: Option<OracleOutcome> = match arm {
                            "x-direct" => direct_x_groebner_decompose(
                                &kc,
                                &fb_x.subspace_basis,
                                &fb_x.points,
                                &index_x,
                                &st,
                                target,
                                m,
                                engine,
                                o.node_budget,
                            ),
                            "x-chained" => {
                                let (rel, stats) = groebner_decompose(
                                    &kc,
                                    &fb_x,
                                    &index_x,
                                    &st,
                                    target,
                                    m,
                                    engine,
                                    o.node_budget,
                                );
                                if stats.unsupported {
                                    None
                                } else {
                                    Some(OracleOutcome {
                                        relation: rel,
                                        used_t: false,
                                        complete: !stats.exhausted,
                                        n_vars: vars,
                                        n_equations: (m - 1) * n as usize,
                                        degree: 0,
                                        effort: stats.splits as u64,
                                        built_degree: stats.max_degree_built,
                                        oversize: stats.oversize,
                                    })
                                }
                            }
                            "sym" => symmetrised_groebner_decompose(
                                &kc,
                                &fb_u,
                                &st,
                                target,
                                m,
                                engine,
                                o.node_budget,
                            ),
                            "sym-chained" => chained_symmetrised_groebner_decompose(
                                &kc,
                                &fb_u,
                                &st,
                                target,
                                engine,
                                o.node_budget,
                            ),
                            _ => unreachable!(),
                        };
                        let wall = t0.elapsed().as_nanos();
                        let xors = f4_word_ops_thread().wrapping_sub(before);
                        acc.attempts += 1;
                        acc.word_xors += xors;
                        acc.wall_ns += wall;
                        if truth[ti] {
                            acc.truth += 1;
                        }
                        let verdict = match &outcome {
                            None => {
                                acc.unbuildable += 1;
                                "unbuildable"
                            }
                            Some(oc) => {
                                acc.n_vars = oc.n_vars;
                                acc.n_equations = oc.n_equations;
                                acc.degree = acc.degree.max(oc.degree);
                                acc.built_degree = acc.built_degree.max(oc.built_degree);
                                acc.splits += oc.effort;
                                if oc.oversize > 0 {
                                    acc.oversize_targets += 1;
                                }
                                match (&oc.relation, oc.complete) {
                                    (Some(idx), _) => {
                                        // Re-sum in the group: R, or R with T appended.
                                        let sum = idx.iter().fold(FastPoint::INFINITY, |s, &i| {
                                            fc.add(s, fc.lift(&points[i]))
                                        });
                                        let tgt = fast_targets[ti];
                                        let ok = sum == tgt || (on_u && sum == fc.add(tgt, t_fast));
                                        if !ok || !truth[ti] {
                                            acc.gate_failures += 1;
                                        }
                                        if oc.used_t {
                                            acc.via_t += 1;
                                        }
                                        acc.found += 1;
                                        acc.word_xors_found += xors;
                                        if truth[ti] {
                                            acc.found_on_truth += 1;
                                        }
                                        "found"
                                    }
                                    (None, true) => {
                                        if truth[ti] {
                                            acc.gate_failures += 1;
                                        }
                                        acc.refuted += 1;
                                        "refuted"
                                    }
                                    (None, false) => {
                                        acc.budget += 1;
                                        "budget"
                                    }
                                }
                            }
                        };
                        acc.per_attempt.push(json!({
                            "target": ti, "verdict": verdict, "truth": truth[ti],
                            "word_xors": xors, "wall_ms_contended": wall as f64 / 1e6,
                        }));
                    }
                }
                let solved = acc.found + acc.refuted;
                let solved_frac = if acc.attempts > 0 {
                    solved as f64 / acc.attempts as f64
                } else {
                    0.0
                };
                let feasible = acc.found > 0 && solved_frac >= 0.5 && acc.gate_failures == 0;
                let xor_per_attempt = if acc.attempts > 0 {
                    acc.word_xors as f64 / acc.attempts as f64
                } else {
                    f64::NAN
                };
                let xor_per_rel = if acc.found > 0 {
                    acc.word_xors as f64 / acc.found as f64
                } else {
                    f64::INFINITY
                };
                let gae_per_rel = xor_per_rel * gae_per_xor;
                // Relation phase to a full-rank matrix: cols + 1 relations.
                let rel_phase_gae = gae_per_rel * (cols as f64 + 1.0);
                let s_rel = rel_phase_gae / (r as f64).sqrt();
                eprintln!(
                    "  n={n} m={m} ell={ell}{} {arm:11} |F|={:5} cols={:4} vars={:2} found {:2} refuted {:2} budget {:2} unbuildable {:2} truth {:2} gate {} | XOR/att {:.3e} XOR/rel {:.3e} S_rel {:.3e} wall {:.1}s{}",
                    if invariant { "inv" } else { "" },
                    points.len(),
                    cols,
                    acc.n_vars,
                    acc.found,
                    acc.refuted,
                    acc.budget,
                    acc.unbuildable,
                    acc.truth,
                    acc.gate_failures,
                    xor_per_attempt,
                    xor_per_rel,
                    s_rel,
                    acc.wall_ns as f64 / 1e9,
                    if acc.wall_truncated { " WALL-TRUNCATED" } else { "" }
                );
                cells.push(json!({
                    "m": m, "ell": ell, "arm": arm, "v_invariant": invariant,
                    "frame": if on_u { "u (kernel-saturated, +T-closed)" } else { "x (plain)" },
                    "base_points": points.len(), "columns": cols,
                    "n_vars": acc.n_vars, "n_equations": acc.n_equations,
                    "boolean_degree": acc.degree, "built_degree": acc.built_degree,
                    "oversize_targets": acc.oversize_targets,
                    "attempts": acc.attempts, "found": acc.found, "refuted": acc.refuted,
                    "budget_exhausted": acc.budget, "unbuildable": acc.unbuildable,
                    "truth_decomposable": acc.truth, "found_on_truth": acc.found_on_truth,
                    "gate_failures": acc.gate_failures, "via_t": acc.via_t,
                    "solved_fraction": solved_frac,
                    "hit_rate": if acc.attempts > 0 { acc.found as f64 / acc.attempts as f64 } else { 0.0 },
                    "exact_hit_rate": if acc.attempts > 0 { acc.truth as f64 / acc.attempts as f64 } else { 0.0 },
                    "completeness": if acc.truth > 0 { acc.found_on_truth as f64 / acc.truth as f64 } else { f64::NAN },
                    "feasible": feasible,
                    "word_xors": acc.word_xors, "word_xors_on_found": acc.word_xors_found,
                    "xor_per_attempt": xor_per_attempt, "xor_per_relation": xor_per_rel,
                    "gae_per_relation": gae_per_rel,
                    "relation_phase_gae": rel_phase_gae,
                    "S_relation_phase": s_rel,
                    "splits": acc.splits,
                    "wall_s_contended": acc.wall_ns as f64 / 1e9,
                    "wall_truncated": acc.wall_truncated,
                    "enumeration_reference": if on_u { enum_u.clone() } else { enum_x.clone() },
                    "per_attempt": acc.per_attempt,
                }));
            }
            cells.push(json!({
                "m": m, "ell": ell, "arm": "_base", "v_invariant": invariant,
                "v_basis": v.iter().map(bits).collect::<Vec<u64>>(),
                "fx_points": fb_x.points.len(), "fu_points": fb_u.points.len(),
                "cols_x": cols_x, "cols_u": cols_u,
                "exact_hits_x": truth_x.iter().filter(|&&b| b).count(),
                "exact_hits_u": truth_u.iter().filter(|&&b| b).count(),
                "enumeration_x": enum_x, "enumeration_u": enum_u,
            }));
        }
        out_curves.push(json!({
            "n": n, "a": a, "r": r, "h": h, "irreducible_low_terms": irr_bits,
            "targets": o.targets,
            "ns_per_word_xor_contended": ns_per_word_xor,
            "ns_per_add_contended": ns_per_add,
            "gae_per_word_xor": gae_per_xor,
            "rho_gae_derived": (std::f64::consts::PI * r as f64 / 4.0).sqrt(),
            "cells": cells,
        }));
        if let Some(path) = &o.json {
            let doc = json!({
                "experiment": "round4 lead 16.1: m crossover, plain x-frame vs kernel-saturated u-frame",
                "engine": "InheritedF4 { max_degree: 3 }, default split rule",
                "node_budget": o.node_budget, "seed": o.seed, "targets": o.targets,
                "feasible_rule": "found >= 1, (found + refuted) / attempts >= 0.5, zero gate failures",
                "unit": "word XORs counted by the engine (lower bound: Macaulay build and substitutions uncounted); GAE via measured ns/XOR over ns/add (contended host)",
                "curves": out_curves,
            });
            fs::write(path, serde_json::to_string_pretty(&doc).unwrap()).expect("write json");
        }
    }
}

/// Columns of the relation matrix: distinct `±[h]P` classes over the base
/// (excluding points with `[h]P = O`).
fn column_count(fc: &FastCurve, pts: &[FastPoint], h: u64) -> usize {
    let mut seen = HashSet::new();
    for &p in pts {
        let q = fc.mul_u64(p, h);
        if q.infinity {
            continue;
        }
        let nq = fc.neg(q);
        seen.insert(q.pack().min(nq.pack()));
    }
    seen.len()
}
