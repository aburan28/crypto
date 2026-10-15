//! IC improvement round 2, item 05: in-engine harvest of the multi-variable
//! linear rows the inherited-F4 solver computes at every node and discards.
//!
//! The incumbent is the pipeline's oracle solve: the equations
//! `decomposition_solver_equations` hands the solver, solved by
//! `solve_boolean_system_filtered` with `SolverEngine::default()` and
//! `split_rule_default()`, harvest off.  Each candidate arm is the same call
//! with a `LinearHarvest` configuration set on the thread
//! (`koblitz_groebner::set_linear_harvest`).  All arms run in one process on
//! the same system; the arm order rotates with the target index so no arm
//! always runs first.  `--arms` lists `name:spec` separated by `/`, where spec is `off`,
//! `env` (no override: whatever `KIC_LINEAR_HARVEST` says) or a
//! `parse_linear_harvest` string such as `1` (the default), `every`, `root`,
//! `k=2`, `root,gpv=0.9607`.
//! The first arm is the reference whose roots every other arm must equal.
//!
//! See `research/ic_improve_20261015/05-linear-harvest/PROTOCOL.md`.
//!
//! ```bash
//! cargo run --release --example ic_improve2_05_linear_harvest -- --n 17 --ell 6 --m 3 --targets 4 \
//!     --arms i:off/h:every/r:root/i2:off
//! ```
#![allow(clippy::too_many_arguments, clippy::type_complexity, dead_code)]

use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, f4_profile, f4_profile_reset, f4_word_ops_thread,
    linear_harvest_counts_reset, linear_harvest_counts_thread, parse_linear_harvest,
    set_linear_harvest, solve_boolean_system_filtered, split_rule_default, FieldStructure,
    SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, decomposition_solver_equations, KoblitzCurve,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

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
    cols_sum: u64,
    spec_ops: u64,
    harvests: u64,
    harvested_rows: u64,
    exhausted: bool,
}

fn solve_arm(eqs: &[F2BoolPoly], nv: usize, spec: &str) -> Arm {
    match spec {
        "env" => set_linear_harvest(None),
        "off" => set_linear_harvest(Some(None)),
        s => set_linear_harvest(Some(Some(
            parse_linear_harvest(s).expect("harvest spec"),
        ))),
    }
    let opts = SolveOptions {
        engine: SolverEngine::default(),
        max_solutions: usize::MAX,
        node_budget: usize::MAX,
        split_rule: split_rule_default(),
    };
    f4_profile_reset();
    linear_harvest_counts_reset();
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (mut roots, st) = solve_boolean_system_filtered(eqs, nv, &opts, |_| false);
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
        cols_sum: p.cols,
        spec_ops: p.specialise_word_ops,
        harvests,
        harvested_rows,
        exhausted: st.exhausted,
    }
}

fn arm_json(name: &str, a: &Arm) -> String {
    let per = |x: u64| {
        if a.calls > 0 {
            x as f64 / a.calls as f64
        } else {
            0.0
        }
    };
    format!(
        r#""{name}_cpu":{:.6},"{name}_splits":{},"{name}_reductions":{},"{name}_elim":{},"{name}_word_ops":{},"{name}_calls":{},"{name}_rank_per_call":{:.1},"{name}_cols_per_call":{:.1},"{name}_spec_ops":{},"{name}_harvests":{},"{name}_harvested_rows":{},"{name}_roots":{},"{name}_root_hash":"{:016x}""#,
        a.cpu,
        a.splits,
        a.reductions,
        a.eliminated,
        a.word_ops,
        a.calls,
        per(a.rank_sum),
        per(a.cols_sum),
        a.spec_ops,
        a.harvests,
        a.harvested_rows,
        a.roots.len(),
        a.roots
            .iter()
            .fold(0xcbf29ce484222325u64, |h, &r| (h ^ r).wrapping_mul(0x100000001b3)),
    )
}

fn support(p: &F2BoolPoly) -> u64 {
    p.terms.iter().fold(0u64, |a, t| a | t.mask)
}

struct PipeCell {
    kc: KoblitzCurve,
    fb: crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    st: FieldStructure,
    m: usize,
}

impl PipeCell {
    fn new(a: u8, n: u32, ell: u32, m: usize) -> Self {
        // a = 0 first; a = 1 where a = 0 has no usable subgroup at this n.
        let kc = KoblitzCurve::new(a, n)
            .or_else(|| KoblitzCurve::new(1 - a, n))
            .expect("curve");
        let fb = build_standard_subspace_factor_base(&kc, ell).expect("factor base");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        PipeCell { kc, fb, st, m }
    }
    /// Target `t`: a uniform multiple of the generator, seeded per (n, ℓ, m, t).
    fn system(&self, t: u64, ell: u32, plain: bool) -> (Vec<F2BoolPoly>, usize, u64) {
        let ord = self.kc.subgroup_order.to_u64_digits()[0];
        let mut rng = StdRng::seed_from_u64(
            20261015 ^ (self.kc.n as u64 * 1_000_000 + ell as u64 * 1000 + self.m as u64 * 100 + t),
        );
        loop {
            let k = rng.gen_range(1..ord);
            let target = self.kc.mul(self.kc.generator(), &BigUint::from(k));
            if plain {
                // Item 02's lead-probe layout: the system in its own
                // (non-interleaved) variable order.
                let x_r = match &target {
                    crypto_lib::binary_ecc::BinaryPoint::Affine { x, .. } => x.clone(),
                    _ => continue,
                };
                if let Some(sys) = build_decomposition_system(
                    &self.fb.subspace_basis,
                    &x_r,
                    &self.kc.curve.b,
                    self.m,
                    &self.st,
                ) {
                    return (sys.equations, sys.n_vars, k);
                }
                continue;
            }
            if let Some((eqs, nv)) = decomposition_solver_equations(
                &self.kc,
                &self.fb,
                &self.st,
                &target,
                self.m,
                SolverEngine::default(),
            ) {
                return (eqs, nv, k);
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
    let arms: Vec<(String, String)> = flag("--arms")
        .unwrap_or_else(|| "i:off/h:every/i2:off".into())
        .split('/')
        .map(|s| {
            let (a, b) = s.split_once(':').expect("name:spec");
            (a.to_string(), b.to_string())
        })
        .collect();
    let cell = PipeCell::new(a, n, ell, m);
    let label = format!(
        "pipe-n{n}-l{ell}-m{m}-a{}",
        cell.kc.curve.a.raw_bits().first().copied().unwrap_or(0)
    );
    for t in first..first + targets {
        let (eqs, nv, _k) = cell.system(t, ell, false);
        let ngens = eqs.iter().filter(|p| !p.is_zero()).count();
        let occ = eqs.iter().fold(0u64, |a, p| a | support(p)).count_ones();
        // Rotate the arm order with the target so no arm always runs first.
        let k = arms.len();
        let mut results: Vec<Option<Arm>> = vec![None; k];
        for j in 0..k {
            let idx = (j + t as usize) % k;
            results[idx] = Some(solve_arm(&eqs, nv, &arms[idx].1));
        }
        let results: Vec<Arm> = results.into_iter().map(Option::unwrap).collect();
        let reference = &results[0].roots;
        let ok = results
            .iter()
            .all(|r| !r.exhausted && &r.roots == reference);
        let mut out = format!(
            r#"{{"cell":"{label}","phase":"{phase}","target":{t},"n_vars":{nv},"gens":{ngens},"occ_vars":{occ},"order_offset":{}"#,
            t as usize % k
        );
        for ((name, _), r) in arms.iter().zip(&results) {
            out += ",";
            out += &arm_json(name, r);
        }
        out += &format!(r#","gate_ok":{ok}}}"#);
        println!("{out}");
    }
}
