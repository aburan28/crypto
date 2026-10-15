//! IC improvement round 3, item 02: a better gate for the in-engine
//! linear-row harvest of the Koblitz Gröbner oracle.
//!
//! Two modes, one process per cell, every arm on the same system:
//!
//! * `--mode bench` (default): the full root-set oracle solve item 05
//!   benched (`solve_boolean_system_filtered`, `SolverEngine::default()`,
//!   `split_rule_default()`, no early stop).  Per target it also records
//!   the gate features: `n_vars`, generators, occurring variables, the free
//!   variables and the basis rank at the root's first reduction, the root
//!   reduction's word operations, and the rows the root harvest takes.
//! * `--mode prod`: the production decomposition call,
//!   `groebner_decompose` (the path `DecompositionStrategy::Groebner` takes
//!   with the default options: node budget 4096, stop at the first root
//!   that lifts), on the same targets.  Every decomposition is re-checked in
//!   the group (`Σ P_i = R`).
//!
//! `--arms` lists `name:spec` separated by `/`; spec is `off`, `env` (no
//! override: whatever `KIC_LINEAR_HARVEST` says) or a
//! `parse_linear_harvest` string.  The arm order rotates with the target
//! index.  The first arm is the reference whose roots (bench) or outcome
//! (prod) every other arm must equal.
//!
//! See `research/ic_improve3_20261015/02-harvest-gate/PROTOCOL.md`.
//!
//! ```bash
//! cargo run --release --example ic_improve3_02_harvest_gate -- --n 23 --ell 12 --m 2 \
//!     --targets 4 --arms i:off/g:root,gpv=0.9607/r:root/i2:off
//! ```
#![allow(clippy::too_many_arguments, clippy::type_complexity)]

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_profile, f4_profile_reset, f4_word_ops_thread, linear_harvest_counts_reset,
    linear_harvest_counts_thread, parse_linear_harvest, set_linear_harvest, set_node_trace,
    solve_boolean_system_filtered, split_rule_default, take_node_trace, FieldStructure,
    SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, decomposition_solver_equations, groebner_decompose,
    FrobeniusFactorBase, KoblitzCurve,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;

fn tcpu() -> f64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    unsafe {
        libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
}

fn set_arm(spec: &str) {
    match spec {
        "env" => set_linear_harvest(None),
        "off" => set_linear_harvest(Some(None)),
        s => set_linear_harvest(Some(Some(
            parse_linear_harvest(s).unwrap_or_else(|| panic!("harvest spec {s}")),
        ))),
    }
}

#[derive(Default, Clone)]
struct Arm {
    roots: Vec<u64>,
    cpu: f64,
    splits: usize,
    reductions: usize,
    eliminated: usize,
    word_ops: u64,
    calls: u64,
    rank_sum: u64,
    spec_ops: u64,
    harvests: u64,
    harvested_rows: u64,
    exhausted: bool,
    /// prod mode: the decomposition (sorted factor-base indices), if any.
    found: Option<Vec<usize>>,
    verified: bool,
}

fn opts(budget: usize) -> SolveOptions {
    SolveOptions {
        engine: SolverEngine::default(),
        max_solutions: usize::MAX,
        node_budget: budget,
        split_rule: split_rule_default(),
    }
}

fn bench_arm(eqs: &[F2BoolPoly], nv: usize, spec: &str) -> Arm {
    set_arm(spec);
    f4_profile_reset();
    linear_harvest_counts_reset();
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (mut roots, st) = solve_boolean_system_filtered(eqs, nv, &opts(usize::MAX), |_| false);
    let cpu = tcpu() - c0;
    let p = f4_profile();
    let (harvests, harvested_rows) = linear_harvest_counts_thread();
    set_linear_harvest(None);
    roots.sort_unstable();
    roots.dedup();
    Arm {
        roots,
        cpu,
        splits: st.splits,
        reductions: st.reductions,
        eliminated: st.eliminated,
        word_ops: f4_word_ops_thread() - w0,
        calls: p.calls,
        rank_sum: p.rows,
        spec_ops: p.specialise_word_ops,
        harvests,
        harvested_rows,
        exhausted: st.exhausted,
        verified: true,
        ..Default::default()
    }
}

fn prod_arm(cell: &PipeCell, index_of: &HashMap<(BigUint, BigUint), usize>, target: &BinaryPoint, spec: &str) -> Arm {
    set_arm(spec);
    f4_profile_reset();
    linear_harvest_counts_reset();
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (found, st) = groebner_decompose(
        &cell.kc,
        &cell.fb,
        index_of,
        &cell.st,
        target,
        cell.m,
        SolverEngine::default(),
        4096,
    );
    let cpu = tcpu() - c0;
    let p = f4_profile();
    let (harvests, harvested_rows) = linear_harvest_counts_thread();
    set_linear_harvest(None);
    let verified = match &found {
        Some(idxs) => {
            let sum = idxs
                .iter()
                .fold(BinaryPoint::Infinity, |acc, &i| cell.kc.add(&acc, &cell.fb.points[i]));
            &sum == target
        }
        None => true,
    };
    Arm {
        cpu,
        splits: st.splits,
        reductions: st.reductions,
        eliminated: st.eliminated,
        word_ops: f4_word_ops_thread() - w0,
        calls: p.calls,
        rank_sum: p.rows,
        spec_ops: p.specialise_word_ops,
        harvests,
        harvested_rows,
        exhausted: st.exhausted,
        found,
        verified,
        ..Default::default()
    }
}

fn hash64(xs: impl Iterator<Item = u64>) -> u64 {
    xs.fold(0xcbf29ce484222325u64, |h, r| (h ^ r).wrapping_mul(0x100000001b3))
}

fn arm_json(name: &str, a: &Arm) -> String {
    let per = |x: u64| {
        if a.calls > 0 {
            x as f64 / a.calls as f64
        } else {
            0.0
        }
    };
    let found = match &a.found {
        Some(v) => format!("{:016x}", hash64(v.iter().map(|&i| i as u64))),
        None => "none".into(),
    };
    format!(
        r#""{name}_cpu":{:.6},"{name}_splits":{},"{name}_reductions":{},"{name}_elim":{},"{name}_word_ops":{},"{name}_calls":{},"{name}_rank_per_call":{:.1},"{name}_spec_ops":{},"{name}_harvests":{},"{name}_harvested_rows":{},"{name}_roots":{},"{name}_root_hash":"{:016x}","{name}_exhausted":{},"{name}_found":"{found}","{name}_verified":{}"#,
        a.cpu,
        a.splits,
        a.reductions,
        a.eliminated,
        a.word_ops,
        a.calls,
        per(a.rank_sum),
        a.spec_ops,
        a.harvests,
        a.harvested_rows,
        a.roots.len(),
        hash64(a.roots.iter().copied()),
        a.exhausted,
        a.verified,
    )
}

fn support(p: &F2BoolPoly) -> u64 {
    p.terms.iter().fold(0u64, |a, t| a | t.mask)
}

/// Root features with the harvest off: free variables at the root's first
/// reduction, that reduction's word operations and basis rank (one F4 call),
/// from a run stopped after the first reduction.
fn root_features(eqs: &[F2BoolPoly], nv: usize) -> (u32, u64, f64, u64) {
    set_linear_harvest(Some(None));
    set_node_trace(true);
    f4_profile_reset();
    let w0 = f4_word_ops_thread();
    let _ = solve_boolean_system_filtered(eqs, nv, &opts(1), |_| false);
    let w1 = f4_word_ops_thread();
    let p = f4_profile();
    let trace = take_node_trace();
    set_node_trace(false);
    set_linear_harvest(None);
    let free = trace
        .first()
        .map(|e| nv as u32 - e.assigned.count_ones())
        .unwrap_or(0);
    let rank = if p.calls > 0 {
        p.rows as f64 / p.calls as f64
    } else {
        0.0
    };
    (free, w1 - w0, rank, p.calls)
}

struct PipeCell {
    kc: KoblitzCurve,
    fb: FrobeniusFactorBase,
    st: FieldStructure,
    m: usize,
}

impl PipeCell {
    fn new(a: u8, n: u32, ell: u32, m: usize) -> Self {
        let kc = KoblitzCurve::new(a, n)
            .or_else(|| KoblitzCurve::new(1 - a, n))
            .expect("curve");
        let fb = build_standard_subspace_factor_base(&kc, ell).expect("factor base");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        PipeCell { kc, fb, st, m }
    }
    /// Target `t`: a uniform multiple of the generator, seeded per
    /// (n, ℓ, m, t) exactly as item 05 seeded it.
    fn system(&self, t: u64, ell: u32) -> (Vec<F2BoolPoly>, usize, BinaryPoint) {
        let ord = self.kc.subgroup_order.to_u64_digits()[0];
        let mut rng = StdRng::seed_from_u64(
            20261015 ^ (self.kc.n as u64 * 1_000_000 + ell as u64 * 1000 + self.m as u64 * 100 + t),
        );
        loop {
            let k = rng.gen_range(1..ord);
            let target = self.kc.mul(self.kc.generator(), &BigUint::from(k));
            if let Some((eqs, nv)) = decomposition_solver_equations(
                &self.kc,
                &self.fb,
                &self.st,
                &target,
                self.m,
                SolverEngine::default(),
            ) {
                return (eqs, nv, target);
            }
        }
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let flag = |k: &str| {
        argv.iter()
            .position(|x| x == k)
            .map(|i| argv[i + 1].clone())
    };
    let n: u32 = flag("--n").map(|v| v.parse().unwrap()).unwrap_or(13);
    let ell: u32 = flag("--ell").map(|v| v.parse().unwrap()).unwrap_or(4);
    let m: usize = flag("--m").map(|v| v.parse().unwrap()).unwrap_or(3);
    let a: u8 = flag("--a").map(|v| v.parse().unwrap()).unwrap_or(0);
    let targets: u64 = flag("--targets").map(|v| v.parse().unwrap()).unwrap_or(4);
    let first: u64 = flag("--first").map(|v| v.parse().unwrap()).unwrap_or(0);
    let phase = flag("--phase").unwrap_or_else(|| "pilot".into());
    let mode = flag("--mode").unwrap_or_else(|| "bench".into());
    let arms: Vec<(String, String)> = flag("--arms")
        .unwrap_or_else(|| "i:off/g:root,gpv=0.9607/r:root/i2:off".into())
        .split('/')
        .map(|s| {
            let (a, b) = s.split_once(':').expect("name:spec");
            (a.to_string(), b.to_string())
        })
        .collect();
    let cell = PipeCell::new(a, n, ell, m);
    let index_of = cell.fb.index_map();
    let label = format!(
        "pipe-n{n}-l{ell}-m{m}-a{}",
        cell.kc.curve.a.raw_bits().first().copied().unwrap_or(0)
    );
    for t in first..first + targets {
        let (eqs, nv, target) = cell.system(t, ell);
        let ngens = eqs.iter().filter(|p| !p.is_zero()).count();
        let occ = eqs.iter().fold(0u64, |a, p| a | support(p)).count_ones();
        let (root_free, root_word_ops, root_rank, root_calls) = if mode == "bench" {
            root_features(&eqs, nv)
        } else {
            (0, 0, 0.0, 0)
        };
        let k = arms.len();
        let mut results: Vec<Option<Arm>> = vec![None; k];
        for j in 0..k {
            let idx = (j + t as usize) % k;
            results[idx] = Some(if mode == "bench" {
                bench_arm(&eqs, nv, &arms[idx].1)
            } else {
                prod_arm(&cell, &index_of, &target, &arms[idx].1)
            });
        }
        let results: Vec<Arm> = results.into_iter().map(Option::unwrap).collect();
        let r0 = &results[0];
        let ok = results.iter().all(|r| {
            r.verified
                && if mode == "bench" {
                    !r.exhausted && r.roots == r0.roots
                } else {
                    r.found.is_some() == r0.found.is_some() && r.exhausted == r0.exhausted
                }
        });
        let same_decomp = results.iter().all(|r| r.found == r0.found);
        let mut out = format!(
            r#"{{"cell":"{label}","mode":"{mode}","phase":"{phase}","n":{n},"ell":{ell},"m":{m},"target":{t},"n_vars":{nv},"gens":{ngens},"occ_vars":{occ},"root_free":{root_free},"root_word_ops":{root_word_ops},"root_rank":{root_rank:.1},"root_calls":{root_calls},"order_offset":{}"#,
            t as usize % k
        );
        for ((name, _), r) in arms.iter().zip(&results) {
            out += ",";
            out += &arm_json(name, r);
        }
        out += &format!(r#","gate_ok":{ok},"same_decomposition":{same_decomp}}}"#);
        println!("{out}");
    }
}
