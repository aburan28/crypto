//! # Tables for research/isogeny_conductor_gap_ic_20261010
//!
//! Reads the committed `ecbench` sessions of the round and the certified
//! pairs, and prints the Markdown tables of the README and a JSON twin:
//!
//! * **part A** (fiber-aware relation generation): per curve and arm, the
//!   mean `S`, the relation yield, the relation-phase group-operation
//!   equivalents per relation, the oracle set-up cost, the solver's
//!   operations per relation and solving degree, the online wall, and
//!   each fiber-aware arm's ratio to the blind arm of its oracle family;
//! * **part B** (isogenous pairs): per pair and arm, the crater's and the
//!   floor's figures side by side with their ratio.
//!
//! Every figure is a mean over the session's measured (non-warm-up) runs
//! of one arm on one curve; nothing is computed that the records do not
//! carry.  Warm-up runs and unverified runs are counted and excluded.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_analyze -- \
//!   --sessions research/isogeny_conductor_gap_ic_20261010/sessions \
//!   --pairs research/isogeny_conductor_gap_ic_20261010/pairs \
//!   --out research/isogeny_conductor_gap_ic_20261010/TABLES.md \
//!   --json research/isogeny_conductor_gap_ic_20261010/tables.json
//! ```

use std::collections::BTreeMap;
use std::fmt::Write as _;

use serde_json::{json, Value};

#[derive(Default, Clone)]
struct Agg {
    runs: u64,
    verified: u64,
    warmups: u64,
    s_sum: f64,
    s_min: f64,
    s_max: f64,
    lower_bound: bool,
    online_ns: f64,
    trials: f64,
    relations: f64,
    lookups: f64,
    rel_gae: f64,
    setup_gae: f64,
    fb_gae: f64,
    total_gae: f64,
    fiber_size: Option<u64>,
    fiber_complete: Option<bool>,
    solver_ops: f64,
    solver_calls: f64,
    degree_mean_sum: f64,
    degree_mean_n: u64,
    degree_max: u32,
    n_vars: Option<u64>,
    levels: BTreeMap<String, u64>,
    pinned_by_repeated_row: u64,
    exhausted: u64,
    family: String,
    cofactor: u64,
    r: u64,
    slug: String,
    method: String,
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or(0.0)
}

fn phase<'a>(r: &'a Value, name: &str) -> Option<&'a Value> {
    r["phases"].as_array()?.iter().find(|p| p["name"] == name)
}

fn absorb(a: &mut Agg, r: &Value) {
    if r["warmup"].as_bool().unwrap_or(false) {
        a.warmups += 1;
        return;
    }
    a.runs += 1;
    if r["outcome"]["status"] == "exhausted" {
        a.exhausted += 1;
    }
    if r["outcome"]["status"] != "verified" {
        return;
    }
    if let Some(p) = phase(r, "linear_algebra") {
        a.pinned_by_repeated_row += p["native"]["pinned_by_repeated_row"].as_u64().unwrap_or(0);
    }
    a.verified += 1;
    let s = f(&r["cost"]["s"]);
    a.s_sum += s;
    a.s_min = if a.verified == 1 { s } else { a.s_min.min(s) };
    a.s_max = a.s_max.max(s);
    a.lower_bound |= r["cost"]["lower_bound"].as_bool().unwrap_or(false);
    a.total_gae += f(&r["cost"]["total_gae"]);
    a.online_ns += f(&r["online"]["wall_ns"]);
    if let Some(p) = phase(r, "relations") {
        a.trials += f(&p["native"]["trials"]);
        a.relations += f(&p["native"]["relations"]);
        a.lookups += f(&p["native"]["lookups"]);
        a.rel_gae += f(&p["gae"]);
    }
    if let Some(p) = phase(r, "oracle_setup") {
        a.setup_gae += f(&p["gae"]);
        if let Some(fs) = p["native"]["fiber_size"].as_u64() {
            a.fiber_size = Some(fs);
            a.fiber_complete = Some(p["native"]["fiber_complete"].as_u64() == Some(1));
        }
    }
    if let Some(p) = phase(r, "factor_base") {
        a.fb_gae += f(&p["gae"]);
    }
    if r["solver"].is_object() {
        a.solver_ops += f(&r["solver"]["ops"]);
        a.solver_calls += f(&r["solver"]["calls"]);
        if let Some(d) = r["solver"]["solving_degree_mean"].as_f64() {
            a.degree_mean_sum += d;
            a.degree_mean_n += 1;
        }
        a.degree_max = a
            .degree_max
            .max(r["solver"]["solving_degree_max"].as_u64().unwrap_or(0) as u32);
        a.n_vars = r["solver"]["n_vars"].as_u64();
    }
    *a.levels
        .entry(r["isolation"]["level"].as_str().unwrap_or("?").to_string())
        .or_default() += 1;
}

impl Agg {
    fn s_mean(&self) -> f64 {
        if self.verified == 0 {
            0.0
        } else {
            self.s_sum / self.verified as f64
        }
    }
    fn yield_(&self) -> f64 {
        if self.trials == 0.0 {
            0.0
        } else {
            self.relations / self.trials
        }
    }
    fn rel_gae_per_relation(&self) -> f64 {
        if self.relations == 0.0 {
            0.0
        } else {
            self.rel_gae / self.relations
        }
    }
    fn solver_ops_per_relation(&self) -> Option<f64> {
        (self.solver_calls > 0.0 && self.relations > 0.0).then(|| self.solver_ops / self.relations)
    }
    fn degree_mean(&self) -> Option<f64> {
        (self.degree_mean_n > 0).then(|| self.degree_mean_sum / self.degree_mean_n as f64)
    }
    fn online_ms(&self) -> f64 {
        if self.verified == 0 {
            0.0
        } else {
            self.online_ns / self.verified as f64 / 1e6
        }
    }
    fn level_string(&self) -> String {
        self.levels
            .iter()
            .map(|(k, v)| format!("{k}×{v}"))
            .collect::<Vec<_>>()
            .join(" ")
    }
    fn json(&self) -> Value {
        json!({
            "runs": self.runs, "verified": self.verified, "warmups": self.warmups,
            "s_mean": self.s_mean(), "s_min": self.s_min, "s_max": self.s_max, "s_lower_bound": self.lower_bound,
            "total_gae_mean": if self.verified == 0 { 0.0 } else { self.total_gae / self.verified as f64 },
            "online_ms_mean": self.online_ms(),
            "trials": self.trials, "relations": self.relations, "yield": self.yield_(),
            "lookups": self.lookups,
            "relation_phase_gae_per_relation": self.rel_gae_per_relation(),
            "oracle_setup_gae_mean": if self.verified == 0 { 0.0 } else { self.setup_gae / self.verified as f64 },
            "factor_base_gae_mean": if self.verified == 0 { 0.0 } else { self.fb_gae / self.verified as f64 },
            "fiber_size": self.fiber_size, "fiber_complete": self.fiber_complete,
            "solver_ops_per_relation": self.solver_ops_per_relation(),
            "solving_degree_mean": self.degree_mean(), "solving_degree_max": (self.degree_max > 0).then_some(self.degree_max),
            "solver_n_vars": self.n_vars,
            "levels": self.levels, "pinned_by_repeated_row": self.pinned_by_repeated_row, "exhausted": self.exhausted,
            "family": self.family, "cofactor": self.cofactor, "r": self.r,
            "slug": self.slug, "method": self.method,
        })
    }
}

fn spec_key(spec: &Value) -> String {
    match spec["kind"].as_str().unwrap_or("") {
        "prime_explicit" => format!("p:{}:{}:{}", spec["p"], spec["a"], spec["b"]),
        "binary_explicit" => format!("b:{}:{}:{}", spec["n"], spec["a"], spec["b"]),
        "koblitz" => format!("k:{}:{}", spec["a"], spec["n"]),
        "binary_random" => format!("r:{}:{}", spec["n"], spec["seed"]),
        "prime_search" => format!("s:{}:{}", spec["bits"], spec["seed"]),
        other => format!("?:{other}"),
    }
}

#[derive(Clone)]
struct PairRole {
    pair: String,
    role: String,
    conductor: u64,
    gap: u64,
    j: Option<u64>,
    aut: Option<u64>,
}

fn load_pairs(dir: &str) -> BTreeMap<String, PairRole> {
    let mut out = BTreeMap::new();
    for file in [
        "prime_pairs.json",
        "koblitz_k0_pairs.json",
        "koblitz_k1_pairs.json",
        "binary_pairs.json",
    ] {
        let Ok(text) = std::fs::read_to_string(format!("{dir}/{file}")) else {
            continue;
        };
        let v: Value = serde_json::from_str(&text).expect("pairs JSON");
        for pr in v["pairs"].as_array().unwrap_or(&Vec::new()) {
            let family = pr["family"].as_str().unwrap_or("?");
            let pair = match family {
                "prime" => format!(
                    "prime D={} {}^{} {}b",
                    pr["d_k"], pr["ell"], pr["h"], pr["bits"]
                ),
                "koblitz" => format!("K_{} n={}", pr["a"], pr["n"]),
                _ => format!("binary n={} ell={}", pr["n"], pr["ell"]),
            };
            let gap = pr["gap"].as_u64().unwrap_or(0);
            for c in pr["curves"].as_array().unwrap_or(&Vec::new()) {
                out.insert(
                    spec_key(&c["spec"]),
                    PairRole {
                        pair: pair.clone(),
                        role: c["role"].as_str().unwrap_or("?").to_string(),
                        conductor: c["conductor"].as_u64().unwrap_or(0),
                        gap,
                        j: c["j"].as_u64(),
                        aut: c["aut"].as_u64(),
                    },
                );
            }
            if pr["floor_isomorphic_model"].is_object() {
                out.insert(
                    spec_key(&pr["floor_isomorphic_model"]["spec"]),
                    PairRole {
                        pair: pair.clone(),
                        role: "floor-iso".into(),
                        conductor: gap,
                        gap,
                        j: None,
                        aut: None,
                    },
                );
            }
        }
    }
    out
}

fn fmt_opt(v: Option<f64>, prec: usize) -> String {
    v.map(|x| format!("{x:.prec$}")).unwrap_or("–".into())
}

fn ratio(a: f64, b: f64) -> String {
    if b == 0.0 {
        "–".into()
    } else {
        format!("{:.3}", a / b)
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let sessions = get("--sessions").expect("--sessions DIR");
    let pairs_dir = get("--pairs").expect("--pairs DIR");
    let out = get("--out").expect("--out FILE");
    let json_out = get("--json");
    let roles = load_pairs(&pairs_dir);

    // session name → (curve key, arm) → Agg
    let mut by_session: BTreeMap<String, BTreeMap<(String, String), Agg>> = BTreeMap::new();
    let mut session_ids: BTreeMap<String, String> = BTreeMap::new();
    let mut entries: Vec<_> = std::fs::read_dir(&sessions)
        .expect("sessions dir")
        .flatten()
        .collect();
    entries.sort_by_key(|e| e.file_name());
    for e in entries {
        let path = e.path().join("records.jsonl");
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        let name = e.file_name().to_string_lossy().to_string();
        let groups = by_session.entry(name.clone()).or_default();
        for line in text.lines() {
            let Ok(r) = serde_json::from_str::<Value>(line) else {
                continue;
            };
            if let Some(id) = r["session_id"].as_str() {
                session_ids.insert(name.clone(), id.to_string());
            }
            let key = spec_key(&r["workload"]["curve_spec"]);
            let arm = r["arm"].as_str().unwrap_or("?").to_string();
            let a = groups.entry((key, arm)).or_default();
            a.family = r["workload"]["curve"]["family"]
                .as_str()
                .unwrap_or("?")
                .to_string();
            a.cofactor = r["workload"]["curve"]["cofactor"].as_u64().unwrap_or(0);
            a.r = r["workload"]["curve"]["r"].as_u64().unwrap_or(0);
            a.slug = r["workload"]["curve"]["slug"]
                .as_str()
                .unwrap_or("?")
                .to_string();
            a.method = r["method"]["id"].as_str().unwrap_or("?").to_string();
            absorb(a, &r);
        }
    }

    let mut md = String::new();
    let mut js = json!({"schema": "isogeny_gap_tables/v1", "sessions": {}, "pairs": []});
    writeln!(
        md,
        "# Tables of research/isogeny_conductor_gap_ic_20261010\n"
    )
    .unwrap();
    // The certified pairs, from the constructor's own records.
    writeln!(md, "## Certified pairs\n").unwrap();
    writeln!(md, "| family | pair | field | #E | r | cofactor | f_π | gap (crater→floor) | crater j / Aut | floor j | edges | certificate |").unwrap();
    writeln!(
        md,
        "|---|---|---|---:|---:|---:|---:|---:|---|---|---:|---|"
    )
    .unwrap();
    for file in [
        "prime_pairs.json",
        "koblitz_k0_pairs.json",
        "koblitz_k1_pairs.json",
        "binary_pairs.json",
    ] {
        let Ok(text) = std::fs::read_to_string(format!("{pairs_dir}/{file}")) else {
            continue;
        };
        let v: Value = serde_json::from_str(&text).expect("pairs JSON");
        for pr in v["pairs"].as_array().unwrap_or(&Vec::new()) {
            let family = pr["family"].as_str().unwrap_or("?");
            let (pair, field, order, r, cof, fpi, cj, fj, edges, cert) = match family {
                "prime" => (
                    format!("D={} ℓ^h={}^{}", pr["d_k"], pr["ell"], pr["h"]),
                    format!("p = {} ({} b)", pr["p"], pr["bits"]),
                    pr["group_order"].to_string(), pr["r"].to_string(), pr["cofactor"].to_string(),
                    pr["frobenius_conductor"].to_string(),
                    format!("{} / {}", pr["crater"]["j"], pr["crater"]["aut"]),
                    pr["floor"]["j"].to_string(),
                    pr["edges"].as_array().map(|e| e.len()).unwrap_or(0),
                    "j(O_K) crater; rank-2 E[ℓ] above, rank-1 at the floor; #E′, order-r image, planted scalar transported".to_string(),
                ),
                "koblitz" => (
                    format!("K_{} descent", pr["a"]),
                    format!("F_2^{}", pr["n"]),
                    pr["group_order"].to_string(), pr["r"].to_string(), pr["cofactor"].to_string(),
                    pr["frobenius_conductor"].to_string(),
                    "1 / 2 (τ)".to_string(), "–".to_string(),
                    pr["edges"].as_array().map(|e| e.len()).unwrap_or(0),
                    format!("kernels in F_(q^m); rank-1 E[ℓ] at the floor: {}", pr["floor_certificates"]),
                ),
                _ => (
                    format!("random b, ℓ={}", pr["ell"]),
                    format!("F_2^{} (seed {})", pr["n"], pr["seed_found"]),
                    pr["group_order"].to_string(), pr["r"].to_string(), pr["cofactor"].to_string(),
                    pr["frobenius_conductor"].to_string(),
                    "– / 2".to_string(), "–".to_string(), 1,
                    format!("{}; D_K = {}, m = {}", pr["direction"], pr["fundamental_discriminant"], pr["ext_degree"]),
                ),
            };
            writeln!(md, "| {family} | {pair} | {field} | {order} | {r} | {cof} | {fpi} | {} | {cj} | {fj} | {edges} | {cert} |", pr["gap"]).unwrap();
            js["pairs"].as_array_mut().unwrap().push(json!({"family": family, "pair": pair, "gap": pr["gap"], "r": pr["r"], "cofactor": pr["cofactor"], "frobenius_conductor": pr["frobenius_conductor"]}));
        }
    }
    writeln!(md).unwrap();
    writeln!(md, "Generated by `examples/isogeny_gap_analyze.rs` from the committed sessions; every figure is a mean over the measured runs of one arm on one curve (warm-ups and unverified runs excluded and counted).  `S` is the cold whole-pipeline cost over √r in ecbench's unit; a `≥` marks a lower bound (an algebraic arm's solver work is reported in its own unit beside it, never folded in).  Yield is relations per trial.  `rel/rel` is relation-phase group-operation equivalents per relation.  Wall is the online window, a practicality note.\n").unwrap();

    for (name, groups) in &by_session {
        let sid = session_ids.get(name).cloned().unwrap_or_default();
        let mut sess_json = json!({"session_id": sid, "curves": {}});
        let is_fiber = name.starts_with("fiber_");
        writeln!(md, "## `{name}` (`{sid}`)\n").unwrap();
        if is_fiber {
            writeln!(md, "| curve | family | h | fiber | arm | runs ok (exhausted, collision-pinned) | S mean | S/blind | yield | yield/blind | rel/rel | rel/rel ÷ blind | setup GAE | solver ops/rel | ops ÷ blind | deg mean | deg max | vars | wall ms | levels |").unwrap();
            writeln!(md, "|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|").unwrap();
        } else {
            writeln!(md, "| pair | role | conductor | j | Aut | curve | arm | runs ok | S mean | S min–max | yield | rel/rel | deg mean | deg max | wall ms | levels |").unwrap();
            writeln!(
                md,
                "|---|---|---:|---:|---:|---|---|---:|---:|---|---:|---:|---:|---:|---:|---|"
            )
            .unwrap();
        }
        // curve order: by key; arm order: as in the spec (first seen)
        let mut curves: Vec<String> = groups.keys().map(|(c, _)| c.clone()).collect();
        curves.dedup();
        for ckey in curves {
            let arms: Vec<(&(String, String), &Agg)> =
                groups.iter().filter(|((c, _), _)| c == &ckey).collect();
            let blind_of = |arm: &str| -> Option<&Agg> {
                let fam = arm.split('-').next().unwrap_or("");
                arms.iter()
                    .find(|((_, a), _)| a == &format!("{fam}-blind"))
                    .map(|(_, g)| *g)
            };
            let role = roles.get(&ckey);
            for ((_, arm), a) in &arms {
                let blind = blind_of(arm);
                let s_mark = if a.lower_bound { "≥" } else { "" };
                let s_ratio = blind
                    .map(|b| ratio(a.s_mean(), b.s_mean()))
                    .unwrap_or("–".into());
                let y_ratio = blind
                    .map(|b| ratio(a.yield_(), b.yield_()))
                    .unwrap_or("–".into());
                let rr_ratio = blind
                    .map(|b| ratio(a.rel_gae_per_relation(), b.rel_gae_per_relation()))
                    .unwrap_or("–".into());
                let ops_ratio = match (
                    a.solver_ops_per_relation(),
                    blind.and_then(|b| b.solver_ops_per_relation()),
                ) {
                    (Some(x), Some(y)) => ratio(x, y),
                    _ => "–".into(),
                };
                if is_fiber {
                    writeln!(
                        md,
                        "| `{}` | {} | {} | {} | {} | {}/{} | {}{:.3} | {} | {:.4} | {} | {:.1} | {} | {:.0} | {} | {} | {} | {} | {} | {:.2} | {} |",
                        a.slug, a.family, a.cofactor,
                        a.fiber_size.map(|v| v.to_string()).unwrap_or("–".into()),
                        arm, format!("{}{}{}", a.verified, if a.exhausted > 0 { format!(" ({} exhausted)", a.exhausted) } else { String::new() }, if a.pinned_by_repeated_row > 0 { format!(" ({} collision-pinned)", a.pinned_by_repeated_row) } else { String::new() }), a.runs, s_mark, a.s_mean(), s_ratio, a.yield_(), y_ratio,
                        a.rel_gae_per_relation(), rr_ratio,
                        if a.verified == 0 { 0.0 } else { a.setup_gae / a.verified as f64 },
                        fmt_opt(a.solver_ops_per_relation(), 0), ops_ratio,
                        fmt_opt(a.degree_mean(), 3),
                        if a.degree_max > 0 { a.degree_max.to_string() } else { "–".into() },
                        a.n_vars.map(|v| v.to_string()).unwrap_or("–".into()),
                        a.online_ms(), a.level_string()
                    )
                    .unwrap();
                } else {
                    writeln!(
                        md,
                        "| {} | {} | {} | {} | {} | `{}` | {} | {}/{} | {}{:.3} | {:.3}–{:.3} | {:.4} | {:.1} | {} | {} | {:.2} | {} |",
                        role.map(|r| r.pair.clone()).unwrap_or("–".into()),
                        role.map(|r| r.role.clone()).unwrap_or("–".into()),
                        role.map(|r| r.conductor.to_string()).unwrap_or("–".into()),
                        role.and_then(|r| r.j).map(|j| j.to_string()).unwrap_or("–".into()),
                        role.and_then(|r| r.aut).map(|j| j.to_string()).unwrap_or("–".into()),
                        a.slug, arm, format!("{}{}{}", a.verified, if a.exhausted > 0 { format!(" ({} exhausted)", a.exhausted) } else { String::new() }, if a.pinned_by_repeated_row > 0 { format!(" ({} collision-pinned)", a.pinned_by_repeated_row) } else { String::new() }), a.runs, s_mark, a.s_mean(), a.s_min, a.s_max, a.yield_(),
                        a.rel_gae_per_relation(), fmt_opt(a.degree_mean(), 3),
                        if a.degree_max > 0 { a.degree_max.to_string() } else { "–".into() },
                        a.online_ms(), a.level_string()
                    )
                    .unwrap();
                }
                let mut aj = a.json();
                if let Some(r) = role {
                    aj["pair"] = json!(r.pair);
                    aj["role"] = json!(r.role);
                    aj["conductor"] = json!(r.conductor);
                    aj["gap"] = json!(r.gap);
                }
                sess_json["curves"][&ckey][arm] = aj;
            }
        }
        writeln!(md).unwrap();
        // Part B: pair ratios floor / crater per arm.
        if !is_fiber {
            let mut pairs: BTreeMap<String, BTreeMap<String, BTreeMap<String, &Agg>>> =
                BTreeMap::new();
            for ((ckey, arm), a) in groups {
                if let Some(r) = roles.get(ckey) {
                    pairs
                        .entry(r.pair.clone())
                        .or_default()
                        .entry(arm.clone())
                        .or_default()
                        .insert(r.role.clone(), a);
                }
            }
            if !pairs.is_empty() {
                writeln!(md, "Pair ratios (floor ÷ crater, and the isomorphic floor model ÷ floor where present):\n").unwrap();
                writeln!(md, "| pair | gap | arm | S crater | S floor | S floor/crater | S iso/floor | yield crater | yield floor | yield floor/crater | deg crater | deg floor | wall crater ms | wall floor ms |").unwrap();
                writeln!(
                    md,
                    "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
                )
                .unwrap();
                for (pair, arms) in &pairs {
                    for (arm, by_role) in arms {
                        let (Some(c), Some(fl)) = (by_role.get("crater"), by_role.get("floor"))
                        else {
                            continue;
                        };
                        let gap = roles
                            .values()
                            .find(|r| r.pair == *pair)
                            .map(|r| r.gap)
                            .unwrap_or(0);
                        let iso = by_role
                            .get("floor-iso")
                            .map(|i| ratio(i.s_mean(), fl.s_mean()))
                            .unwrap_or("–".into());
                        writeln!(
                            md,
                            "| {} | {} | {} | {:.3} | {:.3} | {} | {} | {:.4} | {:.4} | {} | {} | {} | {:.2} | {:.2} |",
                            pair, gap, arm, c.s_mean(), fl.s_mean(), ratio(fl.s_mean(), c.s_mean()), iso,
                            c.yield_(), fl.yield_(), ratio(fl.yield_(), c.yield_()),
                            fmt_opt(c.degree_mean(), 3), fmt_opt(fl.degree_mean(), 3), c.online_ms(), fl.online_ms()
                        )
                        .unwrap();
                        sess_json["pairs"][pair][arm] = json!({
                            "gap": gap,
                            "s_floor_over_crater": if c.s_mean() == 0.0 { Value::Null } else { json!(fl.s_mean() / c.s_mean()) },
                            "yield_floor_over_crater": if c.yield_() == 0.0 { Value::Null } else { json!(fl.yield_() / c.yield_()) },
                            "s_iso_over_floor": by_role.get("floor-iso").map(|i| json!(i.s_mean() / fl.s_mean())).unwrap_or(Value::Null),
                            "degree_crater": c.degree_mean(), "degree_floor": fl.degree_mean(),
                            "degree_max_crater": c.degree_max, "degree_max_floor": fl.degree_max,
                        });
                    }
                }
                writeln!(md).unwrap();
            }
        }
        js["sessions"][name] = sess_json;
    }
    std::fs::write(&out, md).expect("write tables");
    if let Some(p) = json_out {
        std::fs::write(&p, serde_json::to_string_pretty(&js).unwrap()).expect("write json");
    }
    eprintln!("tables written to {out}");
}
