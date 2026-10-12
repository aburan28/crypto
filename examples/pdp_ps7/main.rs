//! PS-7 harness: CDCL with an algebraic leaf, and dynamic branching
//! (`research/pdp_oracle_speedup_20261011/protocols/PS-7.md`).
//!
//! ```bash
//! cargo run --release --example pdp_ps7 -- run  --arm a2 --l 8 --corpus corpus_l8.txt --gates gates.jsonl
//! cargo run --release --example pdp_ps7 -- reg  --arm a2 [--leaf quartic|group]
//! cargo run --release --example pdp_ps7 -- fidelity --l 8 --corpus corpus_l8.txt
//! cargo run --release --example pdp_ps7 -- graph-stats --l 5 --n 15
//! cargo run --release --example pdp_ps7 -- summarize FILE.jsonl ...
//! ```
//!
//! Arms (all decide the same question, whether `x_R` has a decomposition
//! `±P₁ ± P₂ ± P₃ = ±R` over the base):
//!
//! | arm | what it is |
//! |:--|:--|
//! | `a1`  | PS-2's plain loop (`kernel.rs`, row-batched inversions) |
//! | `a1u` | the same pair loop over the shell's one-pair-at-a-time leaf (control) |
//! | `a2`  | the library CDCL solver with the leaf as lazy-clause theory, native branching order |
//! | `a2p` | `a2` with the bare guard clause (the library's level-0 path; control) |
//! | `a2l` | the same solver code, lexicographic static order (order-matched control) |
//! | `a2h` | the example-local solver copy with its hooks inactive (fidelity control) |
//! | `a3`  | the same solver code, re-evaluated minimum-vertex-cover order |
//!
//! This is an oracle-stage measurement (`AGENTS.md` section 8): nothing here is an end-to-end ECDLP cost.

#[allow(dead_code)]
#[path = "../pdp_speedup/common.rs"]
mod common;
#[allow(dead_code)]
#[path = "../pdp_speedup/corpus.rs"]
mod corpus;
#[allow(dead_code)]
#[path = "../pdp_speedup/group.rs"]
mod group;
#[allow(dead_code)]
#[path = "../pdp_speedup/kernel.rs"]
mod kernel;

mod acct;
mod dynorder;
mod leaf;
#[allow(dead_code, clippy::all)]
mod sat_hooked;
mod shell;

use common::Rung;
use corpus::{make_targets, read_corpus, Record};
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::binary_semaev_s4::weil_descend_s4;
use crypto_lib::cryptanalysis::semaev_corpus::CORPUS;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use dynorder::XGraph;
use group::{base_points, degenerate_target, lift, ArtinSchreier};
use kernel::ScalarKernel;
use leaf::{plain_loop_counts, rows_entered, BaseTable, Gates, GroupLeaf, Leaf, LeafCounts, QuarticLeaf};
use serde_json::{json, Value};
use shell::{plain_loop, run_shell, Encoding, Engine, Order, Outcome, ShellCfg, ShellRun, ShellState, Stats};
use std::collections::BTreeMap;
use std::io::Write;
use std::time::Instant;

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}

fn arg_s(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn flag(args: &[String], name: &str) -> bool {
    args.iter().any(|a| a == name)
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Arm {
    A1,
    A1u,
    A2,
    A2p,
    A2l,
    A2h,
    A3,
}

impl Arm {
    fn parse(s: &str) -> Arm {
        match s {
            "a1" => Arm::A1,
            "a1u" => Arm::A1u,
            "a2" => Arm::A2,
            "a2p" => Arm::A2p,
            "a2l" => Arm::A2l,
            "a2h" => Arm::A2h,
            "a3" => Arm::A3,
            other => panic!("unknown arm {other}"),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Arm::A1 => "a1",
            Arm::A1u => "a1u",
            Arm::A2 => "a2",
            Arm::A2p => "a2p",
            Arm::A2l => "a2l",
            Arm::A2h => "a2h",
            Arm::A3 => "a3",
        }
    }

    /// (engine, order, encoding) of a shell arm.
    fn shell(self) -> Option<(Engine, Order, Encoding)> {
        match self {
            Arm::A1 | Arm::A1u => None,
            Arm::A2 => Some((Engine::Library, Order::Native, Encoding::Aux)),
            Arm::A2p => Some((Engine::Library, Order::Native, Encoding::Plain)),
            Arm::A2h => Some((Engine::Hooked, Order::Native, Encoding::Aux)),
            Arm::A2l => Some((Engine::Hooked, Order::Lex, Encoding::Aux)),
            Arm::A3 => Some((Engine::Hooked, Order::Mvc, Encoding::Aux)),
        }
    }
}

fn stats_json(s: &Stats) -> Value {
    json!({
        "decisions": s.decisions, "conflicts": s.conflicts, "restarts": s.restarts,
        "propagations": s.propagations, "learnt_clauses": s.learnt_clauses,
        "clause_visits": s.clause_visits, "clause_lit_visits": s.clause_lit_visits,
        "analyze_lit_visits": s.analyze_lit_visits,
        "learnt_lits_raw": s.learnt_lits_raw, "learnt_lits_kept": s.learnt_lits_kept,
        "max_level": s.max_level, "conflict_level_sum": s.conflict_level_sum,
    })
}

fn leaf_json(c: &LeafCounts) -> Value {
    json!({
        "mul": c.mul, "sqr": c.sqr, "inv": c.inv, "add_n": c.add_n, "flip1": c.flip1,
        "zero_tests": c.zero_tests, "member_tests": c.member_tests,
        "exact_fallbacks": c.exact_fallbacks, "evals": c.evals,
    })
}

/// Shell work in counted steps (one step = one literal, clause, assignment or monomial examined).
fn shell_words<L: Leaf>(run: &ShellRun, st: &ShellState<L>, graph: Option<&XGraph>) -> u64 {
    let s = &run.stats;
    let mut w = s.clause_visits
        + s.clause_lit_visits
        + s.analyze_lit_visits
        + s.propagations
        + s.decisions
        + st.literals_returned
        + st.decide_scans;
    if let (Some(d), Some(g)) = (&run.dyn_stats, graph) {
        w += d.decisions * g.n_monomials as u64 + d.mis_steps;
    }
    w
}

struct Ref<'a> {
    found: bool,
    roots: &'a [[u64; 3]],
}

#[allow(clippy::too_many_arguments)]
fn shell_record<L: Leaf>(
    arm: Arm,
    run: &ShellRun,
    st: &ShellState<L>,
    graph: Option<&XGraph>,
    gates: Option<&Gates>,
    refr: &Ref,
    all: bool,
    base: Value,
) -> Value {
    let found = run.outcome == Outcome::Sat;
    let known = run.outcome != Outcome::Unknown;
    let total = st.acct.total_pairs;
    let words = shell_words(run, st, graph);
    let lc = st.leaf.counts();
    let leaf_bits = gates.map(|g| g.price(&lc));
    let ok = known
        && found == refr.found
        && (if all {
            run.roots == refr.roots
        } else {
            run.witness.map(|w| refr.roots.binary_search(&w).is_ok()).unwrap_or(!found)
        });
    let a = &st.acct;
    let mut v = base;
    let o = v.as_object_mut().unwrap();
    o.insert("arm".into(), json!(arm.name()));
    o.insert("outcome".into(), json!(format!("{:?}", run.outcome)));
    o.insert("found".into(), json!(found));
    o.insert("ref_found".into(), json!(refr.found));
    o.insert("ok".into(), json!(ok));
    o.insert("all_roots".into(), json!(all));
    o.insert("roots".into(), json!(run.roots.len()));
    o.insert("pairs_total".into(), json!(total));
    o.insert("leaf_visits".into(), json!(st.leaf_visits));
    o.insert("pairs_saved".into(), json!(total as i64 - st.leaf_visits as i64));
    o.insert("solver".into(), stats_json(&run.stats));
    o.insert(
        "theory".into(),
        json!({
            "frontier_calls": st.frontier_calls, "refuted": st.refuted, "feasible": st.feasible,
            "clauses_returned": st.clauses_returned, "literals_returned": st.literals_returned,
            "first_feasible_at": st.first_feasible_at, "decide_scans": st.decide_scans,
        }),
    );
    o.insert(
        "acct".into(),
        json!({
            "learned": a.learned, "pair_only": a.learned_pair_only, "other": a.learned_other,
            "units": a.learned_units, "cube_sum": a.cube_sum, "marginal_sum": a.marginal_sum,
            "saved_sum": a.saved_sum, "max_cube": a.max_cube, "max_saved": a.max_saved,
            "shorter_than_pair": a.shorter_than_pair, "theory_excluded": a.theory_excluded,
            "evaluated_pairs": a.evaluated_pairs, "excluded_pairs": a.excluded_pairs,
            "len_hist": a.len_hist,
        }),
    );
    if let Some(d) = &run.dyn_stats {
        o.insert(
            "dyn".into(),
            json!({
                "decisions": d.decisions, "differs_from_lex": d.differs_from_lex,
                "no_edge_fallbacks": d.no_edge_fallbacks, "cover_size_sum": d.cover_size_sum,
                "mis_steps": d.mis_steps, "max_mis_steps_one_node": d.max_mis_steps_one_node,
            }),
        );
    }
    o.insert("leaf".into(), leaf_json(&lc));
    o.insert("leaf_bitops".into(), json!(leaf_bits));
    o.insert("shell_words".into(), json!(words));
    o.insert("shell_bitops64".into(), json!(words * 64));
    o.insert(
        "total_bitops".into(),
        json!(leaf_bits.map(|b| b + words * 64)),
    );
    o.insert("wall_ns".into(), json!(run.wall_ns as u64));
    v
}

fn plain_record(
    arm: Arm,
    found: bool,
    pairs: u64,
    total: u64,
    counts: LeafCounts,
    roots: &[[u64; 3]],
    gates: Option<&Gates>,
    refr: &Ref,
    all: bool,
    wall_ns: u128,
    base: Value,
) -> Value {
    let ok = found == refr.found
        && (if all {
            roots == refr.roots
        } else {
            roots.first().map(|w| refr.roots.binary_search(w).is_ok()).unwrap_or(!found)
        });
    let leaf_bits = gates.map(|g| g.price(&counts));
    let mut v = base;
    let o = v.as_object_mut().unwrap();
    o.insert("arm".into(), json!(arm.name()));
    o.insert("outcome".into(), json!(if found { "Sat" } else { "Unsat" }));
    o.insert("found".into(), json!(found));
    o.insert("ref_found".into(), json!(refr.found));
    o.insert("ok".into(), json!(ok));
    o.insert("all_roots".into(), json!(all));
    o.insert("roots".into(), json!(roots.len()));
    o.insert("pairs_total".into(), json!(total));
    o.insert("leaf_visits".into(), json!(pairs));
    o.insert("pairs_saved".into(), json!(total as i64 - pairs as i64));
    o.insert("leaf".into(), leaf_json(&counts));
    o.insert("leaf_bitops".into(), json!(leaf_bits));
    o.insert("shell_words".into(), json!(0));
    o.insert("shell_bitops64".into(), json!(0));
    o.insert("total_bitops".into(), json!(leaf_bits));
    o.insert("wall_ns".into(), json!(wall_ns as u64));
    v
}

fn emit(out: &mut Option<std::fs::File>, v: &Value) {
    let line = v.to_string();
    println!("{line}");
    if let Some(f) = out {
        writeln!(f, "{line}").unwrap();
    }
}

fn open_out(args: &[String]) -> Option<std::fs::File> {
    arg_s(args, "--out").map(|p| {
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(p)
            .unwrap()
    })
}

fn load_gates(args: &[String], rung: &Rung) -> Option<Gates> {
    let path = arg_s(args, "--gates")?;
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    Some(
        Gates::from_ps1(&text, rung.n, rung.l)
            .unwrap_or_else(|| panic!("{path} has no row for n={} l={}", rung.n, rung.l)),
    )
}

/// One target of a rung through one arm with the group-law leaf.
#[allow(clippy::too_many_arguments)]
fn run_group_target(
    arm: Arm,
    rung: &Rung,
    tab: &BaseTable,
    kernel: &ScalarKernel,
    gates: Option<&Gates>,
    graph: Option<&XGraph>,
    r: (u64, u64),
    refr: &Ref,
    all: bool,
    log: bool,
    max_conflicts: u64,
    base: Value,
) -> Value {
    let m = tab.xs.len() as u64;
    let total = m * (m + 1) / 2;
    match arm {
        Arm::A1 => {
            let mut out = Vec::new();
            let t0 = Instant::now();
            let (found, pairs) = kernel.run(r, all, &mut out);
            let wall = t0.elapsed().as_nanos();
            out.sort_unstable();
            out.dedup();
            let counts = plain_loop_counts(m, rows_entered(m, pairs), pairs);
            plain_record(arm, found, pairs, total, counts, &out, gates, refr, all, wall, base)
        }
        Arm::A1u => {
            let mut leaf = GroupLeaf::new(&rung.gf, rung.l, tab, r);
            let t0 = Instant::now();
            let (found, pairs, roots) = plain_loop(&mut leaf, all);
            let wall = t0.elapsed().as_nanos();
            plain_record(arm, found, pairs, total, leaf.counts(), &roots, gates, refr, all, wall, base)
        }
        _ => {
            let (engine, order, encoding) = arm.shell().unwrap();
            let leaf = GroupLeaf::new(&rung.gf, rung.l, tab, r);
            let cfg = ShellCfg {
                l: rung.l,
                engine,
                order,
                encoding,
                all_roots: all,
                graph,
                log_learnt: log && engine == Engine::Hooked,
                max_conflicts,
            };
            let (run, st) = run_shell(&cfg, leaf);
            shell_record(arm, &run, &st, graph, gates, refr, all, base)
        }
    }
}

fn cmd_run(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let path = arg_s(args, "--corpus").expect("--corpus");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let loaded = read_corpus(&text);
    let rung = Rung::new(l);
    assert_eq!(loaded.l, l, "corpus is for a different rung");
    assert_eq!(loaded.modulus, rung.modulus(), "corpus modulus differs");
    let arm = Arm::parse(&arg_s(args, "--arm").expect("--arm"));
    let k: usize = arg(args, "--targets").unwrap_or(16);
    let all = flag(args, "--all");
    let log = !flag(args, "--no-log");
    let summary_only = flag(args, "--summary-only");
    let max_conflicts: u64 = arg(args, "--max-conflicts").unwrap_or(u64::MAX);
    let label = arg_s(args, "--label").unwrap_or_default();
    let gates = load_gates(args, &rung);
    let mut out = open_out(args);

    // The first k uniform targets of the PS-0 corpus, re-derived from the frozen generator.
    let recs: Vec<&Record> = loaded.recs.iter().filter(|r| !r.t.planted).take(k).collect();
    let gen = make_targets(&rung, k, 0);
    for (r, g) in recs.iter().zip(&gen) {
        assert_eq!((r.t.xr, r.t.yr), (g.xr, g.yr), "corpus target differs from the frozen generator");
    }

    let sol = ArtinSchreier::new(&rung.gf);
    let kernel = ScalarKernel::new(&rung.gf, &sol, l);
    let tab = BaseTable::new(base_points(&rung.gf, &sol, l), l);
    let graph = if arm == Arm::A3 { Some(XGraph::new(l)) } else { None };
    let t_all = Instant::now();
    let (mut n_ok, mut n) = (0u64, 0u64);
    for rec in recs {
        let base = json!({
            "kind": "target", "label": label, "l": l, "n": rung.n, "idx": rec.t.idx,
            "xr": format!("{:x}", rec.t.xr), "space": tab.xs.len(),
            "leaf": "group", "wall_includes_accounting": log,
        });
        let refr = Ref { found: !rec.roots.is_empty(), roots: &rec.roots };
        let v = run_group_target(
            arm,
            &rung,
            &tab,
            &kernel,
            gates.as_ref(),
            graph.as_ref(),
            (rec.t.xr, rec.t.yr),
            &refr,
            all,
            log,
            max_conflicts,
            base,
        );
        n += 1;
        n_ok += v["ok"].as_bool().unwrap() as u64;
        if !summary_only {
            emit(&mut out, &v);
        }
    }
    let wall = t_all.elapsed().as_nanos() as u64;
    emit(
        &mut out,
        &json!({"kind": "summary", "label": label, "arm": arm.name(), "l": l, "targets": n, "ok": n_ok,
                "all_roots": all, "wall_ns": wall, "log": log, "rustc_profile": "release"}),
    );
    if n_ok != n {
        eprintln!("VERDICT MISMATCH: {} of {} targets disagree with the corpus", n - n_ok, n);
        std::process::exit(2);
    }
}

/// The rung of a regression instance: its own modulus (not the bench's).
fn rung_of(n: u32, l: u32, low_terms: &[u32]) -> Rung {
    use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, IrreduciblePoly};
    use num_bigint::BigUint;
    let irr = IrreduciblePoly {
        degree: n,
        low_terms: low_terms.to_vec(),
    };
    let gf = Gf2::new(&irr);
    let one = F2mElement::one(n);
    let curve = BinaryCurve {
        m: n,
        irreducible: irr.clone(),
        a: one.clone(),
        b: one,
        generator: BinaryPoint::Infinity,
        order: BigUint::from(1u32),
        cofactor: BigUint::from(1u32),
    };
    Rung { l, n, irr, gf, curve }
}

fn cmd_reg(args: &[String]) {
    let arm = Arm::parse(&arg_s(args, "--arm").expect("--arm"));
    let leaf_kind = arg_s(args, "--leaf").unwrap_or_else(|| "quartic".into());
    let k: usize = arg(args, "--targets").unwrap_or(CORPUS.len());
    let all = flag(args, "--all");
    let log = !flag(args, "--no-log");
    let max_conflicts: u64 = arg(args, "--max-conflicts").unwrap_or(u64::MAX);
    let label = arg_s(args, "--label").unwrap_or_default();
    let mut out = open_out(args);
    let (mut n, mut n_ok, mut skipped) = (0u64, 0u64, 0u64);
    for inst in CORPUS.iter().take(k) {
        let rung = rung_of(inst.n, inst.l, inst.irr_low_terms);
        let xr = rung.gf.from_element(&inst.x_r());
        // The reference: every sorted triple of V^3 through the crate's general S4 evaluator,
        // sharing nothing with the leaf or the shell.
        let reference: Vec<[u64; 3]> = corpus::brute_roots(&rung, xr, 1);
        let refr = Ref {
            found: inst.truly_sat,
            roots: &reference,
        };
        assert_eq!(!reference.is_empty(), inst.truly_sat, "{}: exhaustive roots disagree with truly_sat", inst.name);
        let base = json!({
            "kind": "target", "label": label, "l": inst.l, "n": inst.n, "name": inst.name,
            "xr": format!("{:x}", xr), "labelled_sat": inst.labelled_sat, "truly_sat": inst.truly_sat,
            "wall_includes_accounting": log,
        });
        let graph = if arm == Arm::A3 { Some(XGraph::new(inst.l)) } else { None };
        let v = match leaf_kind.as_str() {
            "group" => {
                let sol = ArtinSchreier::new(&rung.gf);
                let Some(y) = lift(&rung.gf, &sol, xr) else {
                    skipped += 1;
                    continue;
                };
                if degenerate_target(&rung.gf, inst.l, xr) {
                    skipped += 1;
                    continue;
                }
                let kernel = ScalarKernel::new(&rung.gf, &sol, inst.l);
                let tab = BaseTable::new(base_points(&rung.gf, &sol, inst.l), inst.l);
                let mut b = base.clone();
                b["leaf"] = json!("group");
                run_group_target(
                    arm, &rung, &tab, &kernel, None, graph.as_ref(), (xr, y), &refr, all, log, max_conflicts, b,
                )
            }
            _ => {
                let mut b = base.clone();
                b["leaf"] = json!("quartic");
                let mut leaf = QuarticLeaf::new(&rung, xr);
                match arm {
                    Arm::A1 | Arm::A1u => {
                        let m = leaf.space().len() as u64;
                        let t0 = Instant::now();
                        let (found, pairs, roots) = plain_loop(&mut leaf, all);
                        let wall = t0.elapsed().as_nanos();
                        plain_record(arm, found, pairs, m * (m + 1) / 2, leaf.counts(), &roots, None, &refr, all, wall, b)
                    }
                    _ => {
                        let (engine, order, encoding) = arm.shell().unwrap();
                        let cfg = ShellCfg {
                            l: inst.l,
                            engine,
                            order,
                            encoding,
                            all_roots: all,
                            graph: graph.as_ref(),
                            log_learnt: log && engine == Engine::Hooked,
                            max_conflicts,
                        };
                        let (run, st) = run_shell(&cfg, leaf);
                        shell_record(arm, &run, &st, graph.as_ref(), None, &refr, all, b)
                    }
                }
            }
        };
        n += 1;
        n_ok += v["ok"].as_bool().unwrap() as u64;
        emit(&mut out, &v);
    }
    emit(
        &mut out,
        &json!({"kind": "summary", "label": label, "arm": arm.name(), "set": "regression-60", "leaf": leaf_kind,
                "targets": n, "ok": n_ok, "skipped_not_applicable": skipped, "all_roots": all}),
    );
    if n_ok != n {
        eprintln!("VERDICT MISMATCH: {} of {} inputs disagree with the corpus", n - n_ok, n);
        std::process::exit(2);
    }
}

/// Library solver against the example-local copy with its hooks inactive: every counter must agree.
fn cmd_fidelity(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let path = arg_s(args, "--corpus").expect("--corpus");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let loaded = read_corpus(&text);
    let rung = Rung::new(l);
    let k: usize = arg(args, "--targets").unwrap_or(16);
    let sol = ArtinSchreier::new(&rung.gf);
    let tab = BaseTable::new(base_points(&rung.gf, &sol, l), l);
    let mut out = open_out(args);
    let mut bad = 0u64;
    for rec in loaded.recs.iter().filter(|r| !r.t.planted).take(k) {
        let mut results = Vec::new();
        for engine in [Engine::Library, Engine::Hooked] {
            let leaf = GroupLeaf::new(&rung.gf, l, &tab, (rec.t.xr, rec.t.yr));
            let cfg = ShellCfg {
                l,
                engine,
                order: Order::Native,
                encoding: Encoding::Aux,
                all_roots: true,
                graph: None,
                log_learnt: false,
                max_conflicts: u64::MAX,
            };
            let (run, st) = run_shell(&cfg, leaf);
            results.push((run, st.leaf_visits, st.frontier_calls, st.clauses_returned));
        }
        let (a, b) = (&results[0], &results[1]);
        let same = a.0.stats == b.0.stats
            && a.0.outcome == b.0.outcome
            && a.0.roots == b.0.roots
            && a.1 == b.1
            && a.2 == b.2
            && a.3 == b.3;
        bad += (!same) as u64;
        emit(
            &mut out,
            &json!({"kind": "fidelity", "l": l, "idx": rec.t.idx, "identical": same,
                    "library": stats_json(&a.0.stats), "hooked": stats_json(&b.0.stats),
                    "leaf_visits": [a.1, b.1], "frontier_calls": [a.2, b.2], "roots": a.0.roots.len()}),
        );
    }
    emit(&mut out, &json!({"kind": "fidelity_summary", "l": l, "targets": k, "differing": bad}));
    if bad != 0 {
        std::process::exit(2);
    }
}

fn cmd_graph_stats(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let n: u32 = arg(args, "--n").unwrap_or(3 * l);
    let low = arg_s(args, "--low-terms").map(|s| s.split(',').map(|t| t.parse().unwrap()).collect::<Vec<u32>>());
    let rung = match low {
        Some(t) => rung_of(n, l, &t),
        None => Rung::with(n, l),
    };
    let xr = rung.gf.to_element(0x1234_5678 & rung.mask());
    let sys = weil_descend_s4(n, l, &rung.irr, &F2mElement::one(n), &xr);
    let shape = dynorder::shape(&sys);
    let g_sys = XGraph::from_system(&sys);
    let g_new = XGraph::new(l);
    let (mut a1, mut a2) = (vec![0u64; g_sys.nx], vec![0u64; g_new.nx]);
    let mut buf = Vec::new();
    g_sys.adjacency(0, 0, &mut a1, &mut buf);
    g_new.adjacency(0, 0, &mut a2, &mut buf);
    println!(
        "{}",
        json!({
            "kind": "graph-stats", "n": n, "l": l,
            "correspondence_equations": g_sys.n_equations, "correspondence_monomials_deg_ge2": g_sys.n_monomials,
            "semaev_equations": sys.semaev.len(),
            "semaev_monomials": sys.semaev.iter().map(|p| p.len()).sum::<usize>(),
            "x_vertices": shape.x_vertices, "x_edges": shape.x_edges,
            "x_edges_within_summand": shape.x_edges_within_summand,
            "x_edges_across_summands": shape.x_edges_across_summands,
            "x_complete_graph_edges": shape.x_vertices * (shape.x_vertices - 1) / 2,
            "x_min_vertex_cover": shape.x_min_vertex_cover,
            "e_vertices": shape.e_vertices, "e_edges": shape.e_edges,
            "xe_edges": 0,
            "xgraph_new_equals_from_system": a1 == a2 && g_sys.n_monomials == g_new.n_monomials,
        })
    );
}

// ── summaries ────────────────────────────────────────────────────────

#[derive(Default)]
struct Agg {
    n: u64,
    ok: u64,
    sums: BTreeMap<&'static str, f64>,
}

fn getf(v: &Value, path: &[&str]) -> f64 {
    let mut cur = v;
    for p in path {
        cur = &cur[*p];
    }
    cur.as_f64().unwrap_or(0.0)
}

fn cmd_summarize(args: &[String]) {
    let mut groups: BTreeMap<(String, String, String, String), Agg> = BTreeMap::new();
    for path in args.iter().filter(|a| !a.starts_with("--")) {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        for line in text.lines() {
            let Ok(v) = serde_json::from_str::<Value>(line) else { continue };
            if v["kind"] != "target" || v["all_roots"] == true {
                continue;
            }
            let set = v["name"].as_str().map(|_| "reg60".to_string()).unwrap_or_else(|| format!("l={}", v["l"]));
            let class = if v["ref_found"] == true { "found" } else { "refuted" };
            let key = (set, v["leaf"].as_str().unwrap_or("").to_string(), v["arm"].as_str().unwrap().to_string(), class.to_string());
            let g = groups.entry(key).or_default();
            g.n += 1;
            g.ok += (v["ok"] == true) as u64;
            let mut add = |k: &'static str, x: f64| *g.sums.entry(k).or_insert(0.0) += x;
            add("pairs_total", getf(&v, &["pairs_total"]));
            add("leaf_visits", getf(&v, &["leaf_visits"]));
            add("pairs_saved", getf(&v, &["pairs_saved"]));
            add("conflicts", getf(&v, &["solver", "conflicts"]));
            add("decisions", getf(&v, &["solver", "decisions"]));
            add("restarts", getf(&v, &["solver", "restarts"]));
            add("theory_refuted", getf(&v, &["theory", "refuted"]));
            add("acct_learned", getf(&v, &["acct", "learned"]));
            add("acct_pair_only", getf(&v, &["acct", "pair_only"]));
            add("acct_cube", getf(&v, &["acct", "cube_sum"]));
            add("acct_marginal", getf(&v, &["acct", "marginal_sum"]));
            add("acct_saved", getf(&v, &["acct", "saved_sum"]));
            add("shorter", getf(&v, &["acct", "shorter_than_pair"]));
            add("leaf_bitops", getf(&v, &["leaf_bitops"]));
            add("shell_bitops64", getf(&v, &["shell_bitops64"]));
            add("total_bitops", getf(&v, &["total_bitops"]));
            add("wall_ns", getf(&v, &["wall_ns"]));
            add("dyn_decisions", getf(&v, &["dyn", "decisions"]));
            add("dyn_differs", getf(&v, &["dyn", "differs_from_lex"]));
        }
    }
    println!("| set | leaf | arm | class | targets | ok | pairs/target | leaves/target | saved/target | conflicts/target | learned/target | saved per learned | cube per learned | marginal per learned | shorter than pair | leaf bitops/target | shell bitops/target | total bitops/target | ms/target |");
    println!("|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((set, leaf, arm, class), g) in &groups {
        let n = g.n as f64;
        let s = |k: &str| g.sums.get(k).copied().unwrap_or(0.0);
        let learned = s("conflicts");
        let per = |a: f64, b: f64| if b > 0.0 { format!("{:.3}", a / b) } else { "-".to_string() };
        println!(
            "| {set} | {leaf} | {arm} | {class} | {} | {} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {} | {} | {} | {:.0} | {:.0} | {:.0} | {:.0} | {:.2} |",
            g.n,
            g.ok,
            s("pairs_total") / n,
            s("leaf_visits") / n,
            s("pairs_saved") / n,
            s("conflicts") / n,
            learned / n,
            per(s("pairs_saved"), learned),
            per(s("acct_cube"), s("acct_pair_only")),
            per(s("acct_marginal"), s("acct_pair_only")),
            s("shorter"),
            s("leaf_bitops") / n,
            s("shell_bitops64") / n,
            s("total_bitops") / n,
            s("wall_ns") / n / 1e6,
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(|s| s.as_str()) {
        Some("run") => cmd_run(&args),
        Some("reg") => cmd_reg(&args),
        Some("fidelity") => cmd_fidelity(&args),
        Some("graph-stats") => cmd_graph_stats(&args),
        Some("summarize") => cmd_summarize(&args),
        _ => {
            eprintln!("usage: pdp_ps7 run|reg|fidelity|graph-stats|summarize ... (see the module docs)");
            std::process::exit(1);
        }
    }
}
