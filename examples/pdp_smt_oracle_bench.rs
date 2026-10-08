//! Point-decomposition oracle comparison: native CDCL+XOR against external
//! SMT solvers (cvc5, Z3) on identical Semaev instances.
//!
//! Instrument of `research/pdp_smt_oracle_20261008/PROTOCOL.md`.  Stage
//! diagnostic: it prices the decomposition oracle only; `S`, end-to-end
//! cost and speedup stay unset.
//!
//! ```bash
//! SMT_CVC5_BIN=/path/to/cvc5 SMT_Z3_BIN=/path/to/z3 \
//! cargo run --release --example pdp_smt_oracle_bench -- \
//!     --n 13 --a 0 --m 3 --fb standard --ell 5 --targets 40 --seed 20261008 \
//!     --solvers cvc5,z3 --encodings bool,bv1 --timeout-s 120 \
//!     --out research/pdp_smt_oracle_20261008/results/n13-m3.jsonl
//! ```
//!
//! For every target (a random multiple of the subgroup generator), the
//! verdict of each arm is checked against exhaustive enumeration and every
//! returned decomposition is re-added in the group.  A disagreement is
//! recorded and counted; it is a bug report, not a result.  One JSON line
//! per (target, arm), then one summary line per arm.
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::koblitz_index_calculus::*;
use crypto_lib::cryptanalysis::smt_oracle::{
    find_solver, smt_decompose_detailed, SmtEncoding, SmtSolveOptions, SmtSolverKind,
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::io::Write;
use std::time::{Duration, Instant};

struct Args {
    n: u32,
    a: u8,
    m: usize,
    /// `standard` (polynomial-basis subspace of dimension `ell`) or `frobenius` (orbit base, index 0).
    fb: String,
    ell: Option<u32>,
    targets: usize,
    seed: u64,
    solvers: Vec<SmtSolverKind>,
    encodings: Vec<SmtEncoding>,
    timeout: Duration,
    max_models: usize,
    macaulay_degree: Option<u32>,
    out: Option<String>,
}

fn parse_args() -> Args {
    let mut args = Args {
        n: 13,
        a: 1,
        m: 3,
        fb: "standard".to_string(),
        ell: None,
        targets: 20,
        seed: 20261008,
        solvers: vec![SmtSolverKind::Cvc5, SmtSolverKind::Z3],
        encodings: vec![SmtEncoding::Bool, SmtEncoding::BitVec1],
        timeout: Duration::from_secs(120),
        max_models: 64,
        macaulay_degree: Some(2),
        out: None,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let value = argv.get(i + 1).cloned().unwrap_or_default();
        match argv[i].as_str() {
            "--n" => args.n = value.parse().expect("--n"),
            "--a" => args.a = value.parse().expect("--a"),
            "--m" => args.m = value.parse().expect("--m"),
            "--fb" => args.fb = value,
            "--ell" => args.ell = Some(value.parse().expect("--ell")),
            "--targets" => args.targets = value.parse().expect("--targets"),
            "--seed" => args.seed = value.parse().expect("--seed"),
            "--timeout-s" => args.timeout = Duration::from_secs(value.parse().expect("--timeout-s")),
            "--max-models" => args.max_models = value.parse().expect("--max-models"),
            "--macaulay" => {
                args.macaulay_degree = match value.as_str() {
                    "none" | "0" => None,
                    d => Some(d.parse().expect("--macaulay")),
                }
            }
            "--solvers" => {
                args.solvers = value
                    .split(',')
                    .map(|s| SmtSolverKind::parse(s).unwrap_or_else(|| panic!("solver {s}")))
                    .collect()
            }
            "--encodings" => {
                args.encodings = value
                    .split(',')
                    .map(|s| SmtEncoding::parse(s).unwrap_or_else(|| panic!("encoding {s}")))
                    .collect()
            }
            "--out" => args.out = Some(value),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    args
}

#[derive(Default)]
struct ArmTotals {
    cells: usize,
    found: usize,
    refuted: usize,
    exhausted: usize,
    disagreements: usize,
    spurious: usize,
    native_conflicts: u64,
    solver_conflicts: Option<u64>,
    solver_decisions: Option<u64>,
    solver_calls: usize,
    wall_ms: Vec<f64>,
}

impl ArmTotals {
    fn add_counter(slot: &mut Option<u64>, value: Option<u64>, first: bool) {
        *slot = match (if first { Some(0) } else { *slot }, value) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        };
    }

    fn median_ms(&self) -> f64 {
        let mut v = self.wall_ms.clone();
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    }
}

fn verify_sum(kc: &KoblitzCurve, fb: &FrobeniusFactorBase, ids: &[usize], target: &BinaryPoint) -> bool {
    let sum = ids
        .iter()
        .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
    &sum == target
}

fn main() {
    let args = parse_args();
    let kc = KoblitzCurve::new(args.a, args.n).expect("usable curve");
    let fb = match args.fb.as_str() {
        "standard" => {
            let ell = args.ell.unwrap_or((args.n as usize).div_ceil(args.m) as u32);
            build_standard_subspace_factor_base(&kc, ell).expect("standard-subspace factor base")
        }
        "frobenius" => build_frobenius_factor_base(&kc, 0).expect("frobenius factor base"),
        other => panic!("unknown factor base {other}"),
    };
    let st = FieldStructure::new(args.n, &kc.curve.irreducible);
    let index = fb.index_map();
    let r = kc.subgroup_order.to_u64().expect("subgroup order fits u64");
    let ell = fb.subspace_basis.len();

    let mut rng = rand::rngs::StdRng::seed_from_u64(args.seed);
    let targets: Vec<(u64, BinaryPoint)> = (0..args.targets)
        .map(|_| {
            let k = rng.gen_range(1..r);
            (k, kc.mul(kc.generator(), &BigUint::from(k)))
        })
        .collect();

    let mut arms: Vec<(String, Option<(SmtSolverKind, SmtEncoding, std::path::PathBuf)>)> =
        vec![("native-cdcl-xor".to_string(), None)];
    for &kind in &args.solvers {
        let Some(binary) = find_solver(kind) else {
            eprintln!("{} not found (PATH or SMT_{}_BIN); arm skipped", kind.name(), kind.name().to_ascii_uppercase());
            continue;
        };
        for &encoding in &args.encodings {
            arms.push((
                format!("{}-{}", kind.name(), encoding.name()),
                Some((kind, encoding, binary.clone())),
            ));
        }
    }

    let mut out: Box<dyn Write> = match &args.out {
        Some(path) => {
            if let Some(parent) = std::path::Path::new(path).parent() {
                std::fs::create_dir_all(parent).expect("create output directory");
            }
            Box::new(std::fs::File::create(path).expect("create output file"))
        }
        None => Box::new(std::io::stdout()),
    };
    let mut emit = |v: Value| {
        writeln!(out, "{v}").expect("write");
    };

    let header = json!({
        "kind": "pdp_smt_oracle_header",
        "protocol": "research/pdp_smt_oracle_20261008/PROTOCOL.md",
        "n": args.n, "a": args.a, "m": args.m, "factor_base": args.fb, "ell": ell,
        "factor_base_points": fb.points.len(),
        "subgroup_order": r,
        "targets": args.targets, "seed": args.seed,
        "max_models": args.max_models, "macaulay_degree": args.macaulay_degree,
        "trace_constraint": true,
        "timeout_s": args.timeout.as_secs(),
        "arms": arms.iter().map(|(name, _)| name.clone()).collect::<Vec<_>>(),
        "wall_time_note": "not isolated; practicality note only, counts are the metric",
    });
    emit(header);

    let mut totals: Vec<ArmTotals> = arms.iter().map(|_| ArmTotals::default()).collect();
    for (k, target) in &targets {
        let reference = enumerate_decompose(&kc, &fb, &index, target, args.m);
        for (arm_idx, (arm, smt)) in arms.iter().enumerate() {
            let start = Instant::now();
            let (found, stats, solver_calls, conflicts, decisions, rows, error) = match smt {
                None => {
                    let (out, stats) = sat_decompose_with(
                        &kc,
                        &fb,
                        &index,
                        &st,
                        target,
                        args.m,
                        args.max_models,
                        args.macaulay_degree,
                        SatDecompositionOptions::default(),
                    );
                    let (calls, conflicts) = (stats.solver_calls, stats.conflicts);
                    (out, stats, calls, Some(conflicts), None, None, None)
                }
                Some((kind, encoding, binary)) => {
                    let mut options = SmtSolveOptions::new(*kind, binary.clone());
                    options.encoding = *encoding;
                    options.timeout = args.timeout;
                    options.max_models = args.max_models;
                    options.macaulay_degree = args.macaulay_degree;
                    let (out, stats, report) =
                        smt_decompose_detailed(&kc, &fb, &index, &st, target, args.m, &options);
                    (
                        out,
                        stats,
                        report.calls,
                        report.conflicts,
                        report.decisions,
                        Some(report.rows),
                        report.error,
                    )
                }
            };
            let ms = start.elapsed().as_secs_f64() * 1000.0;
            let decided = found.is_some() || stats.refuted;
            let agrees = if decided {
                found.is_some() == reference.is_some()
            } else {
                true // an exhausted attempt makes no claim
            };
            let sum_ok = found
                .as_ref()
                .map(|ids| verify_sum(&kc, &fb, ids, target))
                .unwrap_or(true);
            let t = &mut totals[arm_idx];
            let first = t.cells == 0;
            t.cells += 1;
            t.found += usize::from(found.is_some());
            t.refuted += usize::from(stats.refuted);
            t.exhausted += usize::from(stats.exhausted);
            t.disagreements += usize::from(!agrees || !sum_ok);
            t.spurious += stats.spurious;
            t.solver_calls += solver_calls;
            if smt.is_none() {
                t.native_conflicts += stats.conflicts;
            } else {
                ArmTotals::add_counter(&mut t.solver_conflicts, conflicts, first);
                ArmTotals::add_counter(&mut t.solver_decisions, decisions, first);
            }
            t.wall_ms.push(ms);
            emit(json!({
                "kind": "pdp_smt_oracle_cell",
                "n": args.n, "a": args.a, "m": args.m, "target_scalar": k, "arm": arm,
                "reference_found": reference.is_some(),
                "found": found.is_some(), "refuted": stats.refuted, "exhausted": stats.exhausted,
                "agrees": agrees, "sum_verified": sum_ok, "spurious": stats.spurious,
                "models": stats.models, "solver_calls": solver_calls,
                "implied_rows": stats.implied_rows, "rows": rows,
                "conflicts": conflicts, "decisions": decisions,
                "wall_ms": ms, "error": error,
            }));
        }
    }
    for (arm_idx, (arm, _)) in arms.iter().enumerate() {
        let t = &totals[arm_idx];
        emit(json!({
            "kind": "pdp_smt_oracle_summary",
            "n": args.n, "a": args.a, "m": args.m, "arm": arm,
            "cells": t.cells, "found": t.found, "refuted": t.refuted,
            "decided": t.found + t.refuted, "exhausted": t.exhausted,
            "disagreements": t.disagreements, "spurious": t.spurious,
            "solver_calls": t.solver_calls,
            "native_conflicts": if arm_idx == 0 { Some(t.native_conflicts) } else { None },
            "solver_conflicts": t.solver_conflicts, "solver_decisions": t.solver_decisions,
            "total_wall_ms": t.wall_ms.iter().sum::<f64>(), "median_wall_ms": t.median_ms(),
        }));
    }
}
