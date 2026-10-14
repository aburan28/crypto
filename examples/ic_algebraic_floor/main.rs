//! `ic_algebraic_floor`: cost per relation of algebraic and SAT
//! decompositions for `m ≥ 3` summands against the meet-in-the-middle
//! (pair-table) floor, on the same binary Koblitz cells and targets
//! (`research/ic_m3_algebraic_vs_mitm_20261014/`).
//!
//! ```bash
//! cargo run --release --example ic_algebraic_floor -- run --m 3 \
//!     --cells 1:13:6,1:17:7 --arms f4-chain,f6-ic,sat,s4-hybrid-g4 \
//!     --pool 2048 --kyes 6 --kno 6 --budget-secs 120 --seed 20261014 --out DIR
//! cargo run --release --example ic_algebraic_floor -- analyze --dir DIR --out ANALYSIS.json
//! ```
//!
//! A cell is `a:n:d`: the Koblitz curve `K_a / GF(2^n)` with the standard
//! subspace factor base `V = span(1, z, …, z^{d−1})` (every point with
//! abscissa in `V`, both signs).  A relation is a target `R = [k]G` that
//! decomposes as `R = P_1 + … + P_m` over the base; every arm answers the
//! decision "does `R` decompose, and with which summands", stopping at the
//! first decomposition.
//!
//! * The floor (`truth`) is `PairSumTable` probing (`m = 3`: `|F|` probes a
//!   target, `m = 4`: `|F|²`), run over a pool of seeded uniform targets.
//!   It is exact and complete, so it is also the ground truth; on small
//!   cells it is itself checked against exhaustive enumeration
//!   (`pdp_split::exhaustive_decompositions`).
//! * The algebraic arms run on `kyes` decomposable and `kno` non-decomposable
//!   targets drawn in pool order, each solve in its own child process under
//!   a wall budget and an address-space cap, so a solve that does not finish
//!   is recorded as censored at the budget, never dropped.
//!
//! A stage diagnostic (AGENTS.md §8): decomposition cost per relation,
//! never an end-to-end ECDLP.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, BinaryInstance};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, FieldStructure, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, groebner_decompose_f6_ic,
    sat_decompose_with, FrobeniusFactorBase, KoblitzCurve, PairSumTable,
    SatDecompositionOptions,
};
use crypto_lib::cryptanalysis::pdp_split::{
    build_for_instance, exhaustive_decompositions, solve_branch, summands_of,
    verify_decomposition,
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn num<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    flag(args, name)
        .map(|v| {
            v.parse()
                .unwrap_or_else(|_| panic!("{name} takes a number, got {v}"))
        })
        .unwrap_or(default)
}

/// One cell's curve, base and fast views, built the same way in the
/// orchestrator and in every child.
struct Cell {
    a: u8,
    n: u32,
    d: u32,
    kc: KoblitzCurve,
    fb: FrobeniusFactorBase,
    inst: BinaryInstance,
    r: u64,
}

impl Cell {
    fn new(a: u8, n: u32, d: u32) -> Self {
        let kc = KoblitzCurve::new(a, n).expect("Koblitz curve with a usable subgroup");
        let fb = build_standard_subspace_factor_base(&kc, d).expect("standard subspace base");
        let inst = koblitz_instance(a, n).expect("binary instance");
        assert!(
            inst.irreducible.degree == kc.curve.irreducible.degree
                && inst.irreducible.low_terms == kc.curve.irreducible.low_terms,
            "instance and curve must share the field representation"
        );
        let r = kc.subgroup_order.to_u64().expect("r fits u64");
        Self {
            a,
            n,
            d,
            kc,
            fb,
            inst,
            r,
        }
    }

    fn target(&self, k: u64) -> BinaryPoint {
        self.kc.mul(self.kc.generator(), &BigUint::from(k))
    }

    fn x_of(&self, p: &BinaryPoint) -> u64 {
        self.inst.fast.lift(p).x
    }

    fn basis_words(&self) -> Vec<u64> {
        (0..self.d).map(|i| 1u64 << i).collect()
    }

    /// Columns of the relation matrix: points up to sign.
    fn columns(&self) -> usize {
        self.fb.signed_orbits.len()
    }

    fn check_witness(&self, idx: &[usize], target: &BinaryPoint) -> bool {
        let sum = idx
            .iter()
            .fold(BinaryPoint::Infinity, |s, &i| self.kc.add(&s, &self.fb.points[i]));
        &sum == target
    }

    /// Sorted abscissae of a decomposition, the form `pdp_split` compares.
    fn abscissae(&self, idx: &[usize]) -> Vec<u64> {
        let mut xs: Vec<u64> = idx.iter().map(|&i| self.x_of(&self.fb.points[i])).collect();
        xs.sort_unstable();
        xs
    }
}

fn parse_cell(s: &str) -> (u8, u32, u32) {
    let p: Vec<&str> = s.split(':').collect();
    assert_eq!(p.len(), 3, "a cell is a:n:d, got {s}");
    (
        p[0].parse().unwrap(),
        p[1].parse().unwrap(),
        p[2].parse().unwrap(),
    )
}

fn cell_seed(seed: u64, m: usize, n: u32, d: u32) -> u64 {
    seed ^ (m as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
        ^ (n as u64).wrapping_mul(0xc2b2_ae3d_27d4_eb4f)
        ^ (d as u64).wrapping_mul(0x1656_67b1_9e37_79f9)
}

/// Group additions the sorted pair-table scan performs up to and
/// including outer row `row` (`m = 3`: one `R − P_k` a row; `m = 4`: the
/// row's `l + 1` pair sums and as many rests).
fn floor_adds(m: usize, rows_done: u64) -> u64 {
    match m {
        3 => rows_done,
        4 => rows_done * (rows_done + 1),
        _ => unreachable!(),
    }
}

// ── truth: the floor and the ground truth ─────────────────────────────

fn truth(
    cell: &Cell,
    m: usize,
    pool: usize,
    kyes: usize,
    kno: usize,
    seed: u64,
    exhaustive_cap: f64,
) -> Value {
    let t0 = Instant::now();
    let table = PairSumTable::build_full(&cell.kc, &cell.fb).expect("pair table");
    let build_ns = t0.elapsed().as_nanos() as u64;
    let pts = cell.fb.points.len();
    let build_adds = (pts as u64) * (pts as u64 + 1) / 2;
    let mut rng = StdRng::seed_from_u64(cell_seed(seed, m, cell.n, cell.d));
    let mut rows = Vec::with_capacity(pool);
    let mut yes_sel = Vec::new();
    let mut no_sel = Vec::new();
    let (mut sum_ns, mut sum_adds, mut decomposable, mut total_decomps) = (0u64, 0u64, 0u64, 0u64);
    for _ in 0..pool {
        let k = rng.gen_range(1..cell.r);
        let target = cell.target(k);
        let fast = table.curve().lift(&target);
        // Decision mode, timed: the sorted scan, stopping at the first
        // witness, exactly what a relation search pays.
        let mut first: Option<Vec<usize>> = None;
        let t = Instant::now();
        table.witnesses_fast(fast, m, &mut |w| {
            first = Some(w.to_vec());
            false
        });
        let ns = t.elapsed().as_nanos() as u64;
        let rows_done = match &first {
            Some(w) => *w.iter().max().unwrap() as u64 + 1,
            None => pts as u64,
        };
        let adds = floor_adds(m, rows_done);
        // Counting mode, untimed: every sorted decomposition.
        let mut all: Vec<Vec<usize>> = Vec::new();
        table.witnesses_fast(fast, m, &mut |w| {
            let mut v = w.to_vec();
            v.sort_unstable();
            all.push(v);
            true
        });
        all.sort();
        all.dedup();
        for w in &all {
            assert!(cell.check_witness(w, &target), "pair-table witness does not add up");
        }
        assert_eq!(first.is_some(), !all.is_empty(), "decision and count disagree");
        sum_ns += ns;
        sum_adds += adds;
        total_decomps += all.len() as u64;
        let yes = !all.is_empty();
        if yes {
            decomposable += 1;
            if yes_sel.len() < kyes {
                yes_sel.push(k);
            }
        } else if no_sel.len() < kno {
            no_sel.push(k);
        }
        rows.push(json!({"k": k, "yes": yes, "decomps": all.len(), "ns": ns, "adds": adds}));
    }
    // Exhaustive cross-check of the floor on the selected targets, where
    // enumeration is affordable.
    let liftable = (0..1u64 << cell.d)
        .filter(|&x| crypto_lib::cryptanalysis::pdp_split::lift_x(&cell.inst, x).is_some())
        .count();
    let tuples = (0..m).fold(1f64, |acc, i| acc * (liftable + i) as f64 / (i + 1) as f64);
    let mut exhaustive = Vec::new();
    if tuples <= exhaustive_cap {
        for &k in yes_sel.iter().chain(no_sel.iter()) {
            let target = cell.target(k);
            let xr = cell.x_of(&target);
            let ex = exhaustive_decompositions(&cell.inst, &cell.basis_words(), m, xr);
            let fast = table.curve().lift(&target);
            let mut mitm: BTreeSet<Vec<u64>> = BTreeSet::new();
            table.witnesses_fast(fast, m, &mut |w| {
                mitm.insert(cell.abscissae(w));
                true
            });
            exhaustive.push(json!({"k": k, "exhaustive": ex.len(), "mitm": mitm.len(), "equal": ex == mitm}));
        }
    }
    // Interleave yes/no so the early-stop rule sees both strata first.
    let mut selected = Vec::new();
    for i in 0..yes_sel.len().max(no_sel.len()) {
        if let Some(&k) = yes_sel.get(i) {
            selected.push(json!({"k": k, "yes": true}));
        }
        if let Some(&k) = no_sel.get(i) {
            selected.push(json!({"k": k, "yes": false}));
        }
    }
    json!({
        "cmd": "truth", "m": m, "a": cell.a, "n": cell.n, "d": cell.d,
        "curve": cell.kc.label(), "r": cell.r, "seed": seed,
        "fb_points": pts, "fb_columns": cell.columns(), "liftable_abscissae": liftable,
        "table_build_ns": build_ns, "table_build_adds": build_adds,
        "pool": pool, "decomposable": decomposable, "total_decomps": total_decomps,
        "floor_sum_ns": sum_ns, "floor_sum_adds": sum_adds,
        "exhaustive_tuples": tuples, "exhaustive_check": exhaustive,
        "selected": selected, "rows": rows,
    })
}

// ── solve: one algebraic decision in a child process ──────────────────

fn limit_address_space(gib: u64) {
    let bytes = gib << 30;
    let lim = libc::rlimit {
        rlim_cur: bytes,
        rlim_max: bytes,
    };
    // SAFETY: setrlimit reads the struct it is handed.
    unsafe {
        libc::setrlimit(libc::RLIMIT_AS, &lim);
    }
}

/// `…-gG` guesses `G` unknowns; `…-fK` guesses all but `K` of the `m·d`
/// summand unknowns (`g = max(0, m·d − K)`), so the branch systems keep
/// the same size as `d` grows.
fn hybrid_g(arm: &str, m: usize, d: u32) -> Option<usize> {
    if let Some((_, g)) = arm.rsplit_once("-g") {
        return g.parse().ok();
    }
    let (_, k) = arm.rsplit_once("-f")?;
    let k: usize = k.parse().ok()?;
    Some((m * d as usize).saturating_sub(k))
}

fn solve(cell: &Cell, m: usize, arm: &str, k: u64, budget: Duration) -> Value {
    let target = cell.target(k);
    let index_of = cell.fb.index_map();
    let st = FieldStructure::new(cell.kc.n, &cell.kc.curve.irreducible);
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    let node_budget = usize::MAX / 2;
    let x0 = f4_word_ops_thread();
    let t0 = Instant::now();
    let mut extra = json!({});
    let (hit, exhausted, unsupported): (Option<Vec<usize>>, bool, bool) = match arm {
        "f4-chain" | "f6-ic" => {
            let (h, s) = if arm == "f4-chain" {
                groebner_decompose(&cell.kc, &cell.fb, &index_of, &st, &target, m, engine, node_budget)
            } else {
                groebner_decompose_f6_ic(
                    &cell.kc, &cell.fb, &index_of, &st, &target, m, engine, node_budget,
                )
            };
            extra = json!({"reductions": s.reductions, "splits": s.splits,
                "max_degree_built": s.max_degree_built});
            (h, s.exhausted, s.unsupported)
        }
        "sat" => {
            let options = SatDecompositionOptions {
                restrict_to_factor_base: true,
                wall_budget: Some(budget),
                ..SatDecompositionOptions::default()
            };
            let (h, s) = sat_decompose_with(
                &cell.kc, &cell.fb, &index_of, &st, &target, m, 1 << 24, None, options,
            );
            extra = json!({"conflicts": s.conflicts, "models": s.models,
                "solver_calls": s.solver_calls, "spurious": s.spurious, "refuted": s.refuted});
            (h, s.exhausted, s.unsupported)
        }
        a if a.starts_with("s4-hybrid") || a.starts_with("s22-hybrid") => {
            let g = hybrid_g(a, m, cell.d).expect("hybrid arm names its g");
            let plan: Vec<usize> = if a.starts_with("s4") { vec![3] } else { vec![2, 2] };
            assert_eq!(plan.iter().sum::<usize>(), m, "{a} does not match m = {m}");
            let xr = cell.x_of(&target);
            match build_for_instance(&cell.inst, &plan, cell.d as usize, xr) {
                Err(e) => {
                    extra = json!({"error": e});
                    (None, true, true)
                }
                Ok((sys, _)) => {
                    let words = cell.basis_words();
                    let (mut xors, mut branches, mut unresolved, mut timed_out) = (0u64, 0u64, 0u64, false);
                    let mut found: Option<Vec<u64>> = None;
                    let mut max_sd = 0u32;
                    'branches: for guess in 0..(1u64 << g) {
                        let left = budget.saturating_sub(t0.elapsed());
                        if left.is_zero() {
                            timed_out = true;
                            break;
                        }
                        let bo = solve_branch(&sys, g, guess, Some(left), 24);
                        branches += 1;
                        xors += bo.word_xors;
                        max_sd = max_sd.max(bo.solving_degree);
                        if bo.timed_out {
                            timed_out = true;
                            break;
                        }
                        match bo.solutions {
                            None => unresolved += 1,
                            Some(sols) => {
                                for a in sols {
                                    let xs = summands_of(&sys, &words, a);
                                    if verify_decomposition(&cell.inst, &xs, xr) == Some(true) {
                                        found = Some(xs);
                                        break 'branches;
                                    }
                                }
                            }
                        }
                    }
                    extra = json!({"g": g, "n_vars": sys.n_vars, "branches": branches,
                        "unresolved": unresolved, "split_word_xors": xors,
                        "solving_degree_max": max_sd, "witness_abscissae": found});
                    // The abscissae verified on the curve; map back to base
                    // indices only to report a hit.
                    let hit = found.map(|_| Vec::new());
                    (hit, timed_out || (unresolved > 0), false)
                }
            }
        }
        other => panic!("unknown arm {other}"),
    };
    let wall_ns = t0.elapsed().as_nanos() as u64;
    let word_xors = f4_word_ops_thread() - x0;
    let witness_ok = match &hit {
        Some(idx) if !idx.is_empty() => cell.check_witness(idx, &target),
        Some(_) => true, // hybrid: verified on the curve by verify_decomposition
        None => true,
    };
    // A completed "no" is a refutation; an exhausted "no" is no verdict.
    let verdict = if hit.is_some() {
        "yes"
    } else if exhausted || unsupported {
        "none"
    } else {
        "no"
    };
    json!({"cmd": "solve", "arm": arm, "m": m, "a": cell.a, "n": cell.n, "d": cell.d, "k": k,
        "verdict": verdict, "unsupported": unsupported, "wall_ns": wall_ns, "word_xors": word_xors,
        "witness_ok": witness_ok, "extra": extra})
}

// ── run: orchestrate cells, arms and per-target timeouts ──────────────

fn run_child(
    exe: &Path,
    m: usize,
    arm: &str,
    cell: (u8, u32, u32),
    k: u64,
    budget_secs: u64,
    mem_gib: u64,
) -> Value {
    let (a, n, d) = cell;
    let grace = Duration::from_secs(budget_secs + 15);
    let t0 = Instant::now();
    let mut child = Command::new(exe)
        .args([
            "solve", "--m", &m.to_string(), "--arm", arm, "--cell", &format!("{a}:{n}:{d}"),
            "--k", &k.to_string(), "--budget-secs", &budget_secs.to_string(),
            "--mem-gib", &mem_gib.to_string(),
        ])
        // One core per solve, so wall time compares like with like.
        .env("RAYON_NUM_THREADS", "1")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn solve child");
    loop {
        match child.try_wait().expect("wait") {
            Some(status) => {
                let out = child.wait_with_output().expect("collect");
                let text = String::from_utf8_lossy(&out.stdout);
                if let Some(line) = text.lines().rev().find(|l| l.starts_with('{')) {
                    if let Ok(v) = serde_json::from_str::<Value>(line) {
                        return v;
                    }
                }
                let err = String::from_utf8_lossy(&out.stderr);
                let oom = err.contains("memory allocation") || err.contains("alloc");
                return json!({"cmd": "solve", "arm": arm, "m": m, "a": a, "n": n, "d": d, "k": k,
                    "verdict": "crash", "oom": oom, "status": status.to_string(),
                    "wall_ns": t0.elapsed().as_nanos() as u64,
                    "stderr_tail": err.chars().rev().take(400).collect::<String>().chars().rev().collect::<String>()});
            }
            None => {
                if t0.elapsed() > grace {
                    let _ = child.kill();
                    let _ = child.wait();
                    return json!({"cmd": "solve", "arm": arm, "m": m, "a": a, "n": n, "d": d, "k": k,
                        "verdict": "timeout", "wall_ns": t0.elapsed().as_nanos() as u64});
                }
                std::thread::sleep(Duration::from_millis(20));
            }
        }
    }
}

fn censored(v: &Value) -> bool {
    if v["unsupported"].as_bool() == Some(true) {
        return false;
    }
    matches!(v["verdict"].as_str(), Some("timeout") | Some("crash") | Some("none"))
        || v["wall_ns"].as_u64().unwrap_or(0) as f64 / 1e9 >= v["budget_secs"].as_f64().unwrap_or(f64::MAX)
}

fn cmd_run(args: &[String]) {
    let m: usize = num(args, "--m", 3);
    let cells: Vec<(u8, u32, u32)> = flag(args, "--cells")
        .expect("--cells a:n:d,...")
        .split(',')
        .map(parse_cell)
        .collect();
    let arms: Vec<String> = flag(args, "--arms")
        .unwrap_or_default()
        .split(',')
        .filter(|s| !s.is_empty())
        .map(String::from)
        .collect();
    let pool: usize = num(args, "--pool", 2048);
    let kyes: usize = num(args, "--kyes", 6);
    let kno: usize = num(args, "--kno", 6);
    let budget: u64 = num(args, "--budget-secs", 120);
    let mem: u64 = num(args, "--mem-gib", 4);
    let seed: u64 = num(args, "--seed", 20261014);
    let stop_after: usize = num(args, "--stop-after", 4);
    let exhaustive_cap: f64 = num(args, "--exhaustive-cap", 5e6);
    let out = PathBuf::from(flag(args, "--out").expect("--out DIR"));
    std::fs::create_dir_all(&out).unwrap();
    let exe = std::env::current_exe().unwrap();
    let mut stopped: BTreeSet<String> = BTreeSet::new();
    for &(a, n, d) in &cells {
        let cell = Cell::new(a, n, d);
        let tr = truth(&cell, m, pool, kyes, kno, seed, exhaustive_cap);
        let tpath = out.join(format!("truth-m{m}-n{n}-d{d}.json"));
        std::fs::write(&tpath, serde_json::to_string(&tr).unwrap()).unwrap();
        eprintln!(
            "m={m} n={n} d={d}: |F|={} decomposable {}/{} floor {:.1} µs/target, exhaustive checks {}",
            tr["fb_points"], tr["decomposable"], pool,
            tr["floor_sum_ns"].as_u64().unwrap() as f64 / pool as f64 / 1e3,
            tr["exhaustive_check"].as_array().unwrap().len()
        );
        drop(cell);
        let selected: Vec<(u64, bool)> = tr["selected"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| (s["k"].as_u64().unwrap(), s["yes"].as_bool().unwrap()))
            .collect();
        for arm in &arms {
            let path = out.join(format!("arm-{arm}-m{m}-n{n}-d{d}.jsonl"));
            let mut f = std::fs::File::create(&path).unwrap();
            if stopped.contains(arm) {
                writeln!(f, "{}", json!({"cmd": "skipped", "arm": arm, "m": m, "n": n, "d": d,
                    "reason": "an earlier cell of this arm met the early-stop rule"}))
                .unwrap();
                eprintln!("  {arm}: skipped (stopped at a smaller cell)");
                continue;
            }
            let mut lead_censored = 0usize;
            for (i, &(k, yes)) in selected.iter().enumerate() {
                let mut v = run_child(&exe, m, arm, (a, n, d), k, budget, mem);
                v["truth_yes"] = json!(yes);
                v["budget_secs"] = json!(budget);
                v["censored"] = json!(censored(&v));
                let agrees = match v["verdict"].as_str() {
                    Some("yes") => Some(yes),
                    Some("no") => Some(!yes),
                    _ => None,
                };
                v["agrees"] = json!(agrees);
                writeln!(f, "{v}").unwrap();
                f.flush().unwrap();
                eprintln!(
                    "  {arm} k={k} truth={yes} verdict={} {:.3}s",
                    v["verdict"], v["wall_ns"].as_u64().unwrap_or(0) as f64 / 1e9
                );
                if v["unsupported"].as_bool() == Some(true) {
                    eprintln!("  {arm}: unsupported at n={n} (encoder or variable limit)");
                    break;
                }
                if i < stop_after && v["censored"].as_bool().unwrap() {
                    lead_censored += 1;
                }
                if i + 1 == stop_after && lead_censored == stop_after {
                    writeln!(f, "{}", json!({"cmd": "stopped", "arm": arm, "m": m, "n": n, "d": d,
                        "reason": format!("the first {stop_after} targets were all censored")}))
                    .unwrap();
                    stopped.insert(arm.clone());
                    eprintln!("  {arm}: early stop at n={n}");
                    break;
                }
            }
        }
    }
}

fn cmd_solve(args: &[String]) {
    let m: usize = num(args, "--m", 3);
    let arm = flag(args, "--arm").expect("--arm");
    let (a, n, d) = parse_cell(&flag(args, "--cell").expect("--cell"));
    let k: u64 = num(args, "--k", 1);
    let budget = Duration::from_secs(num(args, "--budget-secs", 120));
    let mem: u64 = num(args, "--mem-gib", 4);
    let cell = Cell::new(a, n, d);
    limit_address_space(mem);
    let v = solve(&cell, m, &arm, k, budget);
    println!("{v}");
}

fn cmd_truth(args: &[String]) {
    let m: usize = num(args, "--m", 3);
    let (a, n, d) = parse_cell(&flag(args, "--cell").expect("--cell"));
    let cell = Cell::new(a, n, d);
    let mut v = truth(
        &cell,
        m,
        num(args, "--pool", 512),
        num(args, "--kyes", 6),
        num(args, "--kno", 6),
        num(args, "--seed", 1),
        num(args, "--exhaustive-cap", 5e6),
    );
    v.as_object_mut().unwrap().remove("rows");
    println!("{v}");
}

// ── analyze ───────────────────────────────────────────────────────────

mod analysis;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("run") => cmd_run(&args),
        Some("solve") => cmd_solve(&args),
        Some("truth") => cmd_truth(&args),
        Some("analyze") => analysis::cmd_analyze(&args),
        _ => eprintln!("usage: ic_algebraic_floor run|solve|truth|analyze …"),
    }
}
