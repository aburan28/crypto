//! Refutation-degree trace of the PDP3 Gröbner decomposers (lit2026 W4,
//! arXiv 2603.24890): for every decomposition attempt, where the inherited
//! F4 tree decides its nodes — refuted at degree 0 (a constant), 1 (linear
//! elimination), 2 or the top degree, forced, or stalled — and the first
//! degree at which the attempt's root system is inconsistent.
//!
//! ```bash
//! cargo run --release --example f4_refutation_degree_probe -- [--a 0] [--n 17] [--dim 6] \
//!     [--m 3] [--targets 8] [--budget 20000] [--engines f4,f6p]
//! ```
//!
//! Workload and arms as in `f6_ic_probe` (standard-subspace base, targets
//! `[k]G` for `k = 1 …`, `InheritedF4 { max_degree: 3 }`).  One JSON line per
//! arm and target: the verdict, the solver's counters and the
//! `RefutationDegreeTrace` of that attempt.  The trace reduces every node's
//! system again from scratch, so this probe's wall times mean nothing; time
//! with `f6_ic_probe`.  Set `KIC_F4_EARLY_ABORT=1` to trace the early-abort
//! arm (the verdicts and the node counters must not move).

use crypto_lib::cryptanalysis::koblitz_groebner::{
    early_abort_mode, refutation_degree_trace_take, set_refutation_degree_trace, FieldStructure,
    SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, groebner_decompose_f6_ic,
    groebner_decompose_f6_ic_pair_index, KoblitzCurve,
};
use num_bigint::BigUint;

fn flag(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn read_jsonl(path: &std::path::Path) -> Vec<serde_json::Value> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

/// `--summarize <raw dir>`: tabulate a run's recorded probe outputs.
///
/// - `trace_<cell>_mode<k>.jsonl` (this probe): per cell, engine and mode,
///   the sums over targets of the solver counters and of every trace field.
/// - `time_<cell>_<arm>.round<r>.jsonl` (`f6_ic_probe`): per cell, engine
///   and arm, `total_ms` of every round, its median, min and max; the
///   counters (`reductions`, `splits`, `hits`, `word_xors`) and whether
///   they agree across rounds; and per-round ratios `orig / arm` and
///   `prev / arm` with their median, min and max (rounds paired by index).
///
/// Prints one JSON document.  Every value is copied or derived (sum,
/// median, ratio) from the files named; nothing is measured here.
fn summarize(dir: &str) {
    use serde_json::{json, Map, Value};
    let mut entries: Vec<std::path::PathBuf> = std::fs::read_dir(dir)
        .expect("raw dir")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    entries.sort();
    // Traces.
    let mut traces = Map::new();
    for path in entries.iter().filter(|p| {
        let n = p.file_name().unwrap().to_string_lossy();
        n.starts_with("trace_") && n.ends_with(".jsonl")
    }) {
        let name = path.file_stem().unwrap().to_string_lossy().to_string();
        let mut per_engine: Map<String, Value> = Map::new();
        for line in read_jsonl(path) {
            let arm = line["arm"].as_str().unwrap().to_string();
            let slot = per_engine
                .entry(arm)
                .or_insert_with(|| json!({"attempts": 0, "hits": 0, "witnesses": []}));
            slot["attempts"] = json!(slot["attempts"].as_u64().unwrap() + 1);
            slot["hits"] =
                json!(slot["hits"].as_u64().unwrap() + line["hit"].as_bool().unwrap() as u64);
            slot["witnesses"]
                .as_array_mut()
                .unwrap()
                .push(line["witness"].clone());
            for key in [
                "reductions",
                "splits",
                "propagations",
                "infeasible",
                "eliminated",
                "geo_refute",
            ] {
                let add = line[key].as_u64().unwrap_or(0);
                slot[key] = json!(slot[key].as_u64().unwrap_or(0) + add);
            }
            for (key, value) in line["trace"].as_object().unwrap() {
                match value {
                    Value::Number(n) => {
                        slot[key] = json!(slot[key].as_u64().unwrap_or(0) + n.as_u64().unwrap());
                    }
                    Value::Array(xs) => {
                        let mut acc: Vec<u64> = slot[key]
                            .as_array()
                            .map(|a| a.iter().map(|v| v.as_u64().unwrap()).collect())
                            .unwrap_or_else(|| vec![0; xs.len()]);
                        for (i, x) in xs.iter().enumerate() {
                            acc[i] += x.as_u64().unwrap();
                        }
                        slot[key] = json!(acc);
                    }
                    _ => {}
                }
            }
        }
        // Derived: shares of the refuted nodes and of degree-2 refutations.
        for slot in per_engine.values_mut() {
            let r = slot["word_ops_refuted_nodes"].as_f64().unwrap_or(0.0);
            let s = slot["word_ops_split_nodes"].as_f64().unwrap_or(0.0);
            if r + s > 0.0 {
                slot["refuted_share_of_node_word_ops"] = json!(r / (r + s));
            }
            let deg = |k: &str, d: usize| slot[k][d].as_u64().unwrap_or(0) as f64;
            let at2 =
                deg("refuted_first_degree_fresh", 2) + deg("refuted_first_degree_inherited", 2);
            let at3 =
                deg("refuted_first_degree_fresh", 3) + deg("refuted_first_degree_inherited", 3);
            if at2 + at3 > 0.0 {
                slot["refuted_at_degree2_share"] = json!(at2 / (at2 + at3));
            }
        }
        traces.insert(name, Value::Object(per_engine));
    }
    // Timing rounds: (cell, arm, round) -> lines.
    let mut timing: std::collections::BTreeMap<(String, String), Vec<(u32, Vec<Value>)>> =
        Default::default();
    for path in entries.iter() {
        let n = path.file_name().unwrap().to_string_lossy().to_string();
        let Some(rest) = n
            .strip_prefix("time_")
            .and_then(|r| r.strip_suffix(".jsonl"))
        else {
            continue;
        };
        let (stem, round) = rest
            .split_once(".round")
            .expect("time_<cell>_<arm>.round<r>");
        let (cell, arm) = stem.rsplit_once('_').expect("<cell>_<arm>");
        timing
            .entry((cell.to_string(), arm.to_string()))
            .or_default()
            .push((round.parse().unwrap(), read_jsonl(path)));
    }
    let mut cells = Map::new();
    let totals = |rounds: &Vec<(u32, Vec<Value>)>, engine: &str| -> Vec<(u32, Value)> {
        let mut out: Vec<(u32, Value)> = rounds
            .iter()
            .filter_map(|(r, lines)| {
                lines
                    .iter()
                    .find(|l| l["arm"] == engine && l.get("total_ms").is_some())
                    .map(|l| (*r, l.clone()))
            })
            .collect();
        out.sort_by_key(|(r, _)| *r);
        out
    };
    let cell_names: std::collections::BTreeSet<String> =
        timing.keys().map(|(c, _)| c.clone()).collect();
    for cell in cell_names {
        let mut per_engine = Map::new();
        for engine in ["f4", "f6p"] {
            let mut arms = Map::new();
            let ms_of = |arm: &str| -> Vec<(u32, f64)> {
                timing
                    .get(&(cell.clone(), arm.to_string()))
                    .map(|rounds| {
                        totals(rounds, engine)
                            .into_iter()
                            .map(|(r, l)| (r, l["total_ms"].as_f64().unwrap()))
                            .collect()
                    })
                    .unwrap_or_default()
            };
            for arm in ["orig", "prev", "cand1", "cand2"] {
                let Some(rounds) = timing.get(&(cell.clone(), arm.to_string())) else {
                    continue;
                };
                let lines = totals(rounds, engine);
                if lines.is_empty() {
                    continue;
                }
                let ms: Vec<f64> = lines
                    .iter()
                    .map(|(_, l)| l["total_ms"].as_f64().unwrap())
                    .collect();
                let mut sorted = ms.clone();
                let counters: Vec<Value> = lines
                    .iter()
                    .map(|(_, l)| {
                        json!({"hits": l["hits"], "reductions": l["reductions"], "splits": l["splits"],
                               "geo_refute": l["geo_refute"], "word_xors": l["word_xors"]})
                    })
                    .collect();
                let consistent = counters.windows(2).all(|w| w[0] == w[1]);
                // Per-target verdicts of every round, for the identity check.
                let verdicts: Vec<Vec<Value>> = rounds
                    .iter()
                    .map(|(_, ls)| {
                        ls.iter()
                            .filter(|l| l["arm"] == engine && l.get("target").is_some())
                            .map(|l| {
                                json!([
                                    l["target"],
                                    l["hit"],
                                    l["exhausted"],
                                    l["reductions"],
                                    l["splits"],
                                    l["infeasible"]
                                ])
                            })
                            .collect()
                    })
                    .collect();
                let verdicts_consistent = verdicts.windows(2).all(|w| w[0] == w[1]);
                let mut entry = json!({
                    "rounds": lines.iter().map(|(r, _)| *r).collect::<Vec<_>>(),
                    "total_ms": ms,
                    "median_ms": median(&mut sorted),
                    "min_ms": sorted.first(),
                    "max_ms": sorted.last(),
                    "counters": counters[0],
                    "counters_identical_across_rounds": consistent,
                    "per_target": verdicts[0],
                    "per_target_identical_across_rounds": verdicts_consistent,
                });
                for reference in ["orig", "prev"] {
                    if reference == arm {
                        continue;
                    }
                    // The stage unit (`word_xors`) is deterministic: one
                    // ratio, from round 1 of each arm.
                    if let Some(base) = timing
                        .get(&(cell.clone(), reference.to_string()))
                        .and_then(|rounds| totals(rounds, engine).first().cloned())
                    {
                        let b = base.1["word_xors"].as_f64().unwrap();
                        let c = counters[0]["word_xors"].as_f64().unwrap();
                        entry[format!("word_xors_ratio_{reference}_over_this")] = json!(b / c);
                    }
                    let base = ms_of(reference);
                    let mine = ms_of(arm);
                    let mut ratios: Vec<f64> = base
                        .iter()
                        .filter_map(|(r, b)| mine.iter().find(|(q, _)| q == r).map(|(_, c)| b / c))
                        .collect();
                    if ratios.is_empty() {
                        continue;
                    }
                    let per_round = ratios.clone();
                    let med = median(&mut ratios);
                    entry[format!("ratio_{reference}_over_this")] = json!({
                        "per_round": per_round, "median": med,
                        "min": ratios.first(), "max": ratios.last()
                    });
                }
                arms.insert(arm.to_string(), entry);
            }
            per_engine.insert(engine.to_string(), Value::Object(arms));
        }
        cells.insert(cell, Value::Object(per_engine));
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({"traces": traces, "timing": cells})).unwrap()
    );
}

fn main() {
    if let Some(dir) = flag("--summarize") {
        summarize(&dir);
        return;
    }
    let a: u8 = flag("--a").and_then(|v| v.parse().ok()).unwrap_or(0);
    let n: u32 = flag("--n").and_then(|v| v.parse().ok()).unwrap_or(17);
    let dim: u32 = flag("--dim").and_then(|v| v.parse().ok()).unwrap_or(6);
    let m: usize = flag("--m").and_then(|v| v.parse().ok()).unwrap_or(3);
    let targets: u64 = flag("--targets").and_then(|v| v.parse().ok()).unwrap_or(8);
    let budget: usize = flag("--budget")
        .and_then(|v| v.parse().ok())
        .unwrap_or(20_000);
    let engines = flag("--engines").unwrap_or_else(|| "f4,f6p".to_string());
    let engines: Vec<&str> = engines.split(',').collect();

    let kc = KoblitzCurve::new(a, n).expect("curve has a usable subgroup");
    let fb = build_standard_subspace_factor_base(&kc, dim).expect("factor base");
    let index_of = fb.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    eprintln!(
        "K({a}, {n}), standard subspace dim {dim}: {} base points, m = {m}, {targets} targets, budget {budget}, early abort {}",
        fb.points.len(),
        early_abort_mode()
    );
    for arm in &engines {
        for k in 1..=targets {
            let target = kc.mul(kc.generator(), &BigUint::from(k));
            set_refutation_degree_trace(true);
            let (hit, stats) = match *arm {
                "f4" => groebner_decompose(&kc, &fb, &index_of, &st, &target, m, engine, budget),
                "f6" => {
                    groebner_decompose_f6_ic(&kc, &fb, &index_of, &st, &target, m, engine, budget)
                }
                "f6p" => groebner_decompose_f6_ic_pair_index(
                    &kc, &fb, &index_of, &st, &target, m, engine, budget,
                ),
                other => panic!("unknown engine {other}; use f4, f6 or f6p"),
            };
            let trace = refutation_degree_trace_take().expect("trace is on");
            set_refutation_degree_trace(false);
            if let Some(indices) = &hit {
                let sum = indices
                    .iter()
                    .fold(crypto_lib::binary_ecc::BinaryPoint::Infinity, |acc, &i| {
                        kc.add(&acc, &fb.points[i])
                    });
                assert_eq!(sum, target, "{arm}: witness for target {k} does not add up");
            }
            println!(
                "{{\"arm\":\"{arm}\",\"a\":{a},\"n\":{n},\"dim\":{dim},\"m\":{m},\"target\":{k},\"early_abort\":{},\"hit\":{},\"witness\":{},\"exhausted\":{},\"reductions\":{},\"splits\":{},\"propagations\":{},\"infeasible\":{},\"eliminated\":{},\"max_degree_built\":{},\"geo_refute\":{},\"trace\":{}}}",
                early_abort_mode(),
                hit.is_some(),
                serde_json::to_string(&hit).unwrap(),
                stats.exhausted,
                stats.reductions,
                stats.splits,
                stats.propagations,
                stats.infeasible_branches,
                stats.eliminated,
                stats.max_degree_built,
                stats.geometric_refutations,
                serde_json::to_string(&trace).unwrap()
            );
        }
    }
}
