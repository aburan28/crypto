//! # Tables for research/isogeny_conductor_gap_round2_20261011
//!
//! Reads the round's committed `ecbench` sessions and round 1's certified
//! pairs, and prints the README's tables and a JSON twin:
//!
//! * **A0**: per pair and role (crater, floor, isomorphic floor model),
//!   the oracle set-up, `fiber_samples` and mean `S` of the closed table
//!   under the legacy and the coset closure;
//! * **E1**: per pair and role, for each `fiber_scale=d`, the subgroup size
//!   `|[d]E(F_p)[h]|`, set-up, relation-phase cost and mean `S`, with the
//!   best `d` marked;
//! * **E2**: per pair and role, mean `S` of the automorphism-folded base
//!   (`group=auto`) and the same-column negation base, beside rho.
//!
//! Every figure is a mean over measured (non-warm-up), verified runs of one
//! arm on one curve.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_round2_analyze -- \
//!   --sessions research/isogeny_conductor_gap_round2_20261011/sessions \
//!   --pairs research/isogeny_conductor_gap_ic_20261010/pairs/prime_pairs.json \
//!   --out research/isogeny_conductor_gap_round2_20261011/TABLES.md \
//!   --json research/isogeny_conductor_gap_round2_20261011/tables.json
//! ```

use std::collections::BTreeMap;
use std::fmt::Write as _;

use serde_json::{json, Value};

#[derive(Default, Clone)]
struct Agg {
    runs: u64,
    verified: u64,
    s_sum: f64,
    s_min: f64,
    s_max: f64,
    setup_gae: f64,
    rel_gae: f64,
    relations: f64,
    online_ns: f64,
    fiber_size: Option<u64>,
    fiber_samples: Option<u64>,
    table_entries: Option<u64>,
}

impl Agg {
    fn absorb(&mut self, r: &Value) {
        if r["warmup"].as_bool().unwrap_or(false) {
            return;
        }
        self.runs += 1;
        if r["outcome"]["status"] != "verified" {
            return;
        }
        self.verified += 1;
        let s = r["cost"]["s"].as_f64().unwrap_or(0.0);
        self.s_min = if self.verified == 1 {
            s
        } else {
            self.s_min.min(s)
        };
        self.s_max = self.s_max.max(s);
        self.s_sum += s;
        self.online_ns += r["online"]["wall_ns"].as_f64().unwrap_or(0.0);
        for p in r["phases"].as_array().into_iter().flatten() {
            match p["name"].as_str() {
                Some("oracle_setup") => {
                    self.setup_gae += p["gae"].as_f64().unwrap_or(0.0);
                    let n = &p["native"];
                    self.fiber_size = n["fiber_size"].as_u64().or(self.fiber_size);
                    self.fiber_samples = n["fiber_samples"].as_u64().or(self.fiber_samples);
                    self.table_entries = n["fiber_table_entries"].as_u64().or(self.table_entries);
                }
                Some("relations") => {
                    self.rel_gae += p["gae"].as_f64().unwrap_or(0.0);
                    self.relations += p["native"]["relations"].as_f64().unwrap_or(0.0);
                }
                _ => {}
            }
        }
    }
    fn mean(&self, x: f64) -> f64 {
        if self.verified == 0 {
            0.0
        } else {
            x / self.verified as f64
        }
    }
    fn s(&self) -> f64 {
        self.mean(self.s_sum)
    }
    fn setup(&self) -> f64 {
        self.mean(self.setup_gae)
    }
    fn rel(&self) -> f64 {
        self.mean(self.rel_gae)
    }
    fn json(&self) -> Value {
        json!({
            "runs": self.runs, "verified": self.verified,
            "s_mean": self.s(), "s_min": self.s_min, "s_max": self.s_max,
            "oracle_setup_gae_mean": self.setup(), "relation_phase_gae_mean": self.rel(),
            "online_ms_mean": self.mean(self.online_ns) / 1e6,
            "fiber_size": self.fiber_size, "fiber_samples": self.fiber_samples,
            "fiber_table_entries": self.table_entries,
        })
    }
}

fn spec_key(spec: &Value) -> String {
    format!("p:{}:{}:{}", spec["p"], spec["a"], spec["b"])
}

#[derive(Clone)]
struct Role {
    pair: String,
    order: (i64, u64, u64, u64),
    role: &'static str,
}

fn load_pairs(path: &str) -> BTreeMap<String, Role> {
    let v: Value =
        serde_json::from_str(&std::fs::read_to_string(path).expect("pairs file")).expect("JSON");
    let mut out = BTreeMap::new();
    for pr in v["pairs"].as_array().expect("pairs") {
        let pair = format!(
            "D={} {}^{} {}b h={}",
            pr["d_k"], pr["ell"], pr["h"], pr["bits"], pr["cofactor"]
        );
        let order = (
            pr["d_k"].as_i64().unwrap_or(0),
            pr["bits"].as_u64().unwrap_or(0),
            pr["cofactor"].as_u64().unwrap_or(0),
            pr["ell"].as_u64().unwrap_or(0),
        );
        for c in pr["curves"].as_array().expect("curves") {
            let role = match c["role"].as_str() {
                Some("crater") => "crater",
                Some("floor") => "floor",
                _ => continue,
            };
            out.insert(
                spec_key(&c["spec"]),
                Role {
                    pair: pair.clone(),
                    order,
                    role,
                },
            );
        }
        out.insert(
            spec_key(&pr["floor_isomorphic_model"]["spec"]),
            Role {
                pair,
                order,
                role: "floor-iso",
            },
        );
    }
    out
}

type Table = BTreeMap<((i64, u64, u64, u64), String, &'static str), BTreeMap<String, Agg>>;

fn load_session(dir: &str, roles: &BTreeMap<String, Role>) -> Option<Table> {
    let text = std::fs::read_to_string(format!("{dir}/records.jsonl")).ok()?;
    let mut t: Table = BTreeMap::new();
    for line in text.lines() {
        let r: Value = serde_json::from_str(line).expect("record");
        let key = spec_key(&r["workload"]["curve_spec"]);
        let role = roles
            .get(&key)
            .expect("every curve is a certified pair curve");
        let arm = r["arm"].as_str().unwrap_or("?").to_string();
        t.entry((role.order, role.pair.clone(), role.role))
            .or_default()
            .entry(arm)
            .or_default()
            .absorb(&r);
    }
    Some(t)
}

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|a| a == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    args[i + 1].clone()
}

fn table_json(t: &Table) -> Value {
    let rows: Vec<Value> = t
        .iter()
        .map(|((_, pair, role), arms)| {
            json!({
                "pair": pair, "role": role,
                "arms": arms.iter().map(|(k, a)| (k.clone(), a.json())).collect::<BTreeMap<_, _>>(),
            })
        })
        .collect();
    Value::Array(rows)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let sessions = arg(&args, "--sessions");
    let roles = load_pairs(&arg(&args, "--pairs"));
    let mut md = String::new();
    let mut js = serde_json::Map::new();

    // A0
    let mut a0: Table = BTreeMap::new();
    for name in ["a0_coset_prime26", "a0_coset_h5408"] {
        if let Some(t) = load_session(&format!("{sessions}/{name}"), &roles) {
            a0.extend(t);
        }
    }
    writeln!(md, "## A0. Legacy against coset fiber closure\n").unwrap();
    writeln!(
        md,
        "| pair | role | samples | set-up legacy | set-up coset | legacy / coset | S legacy | S coset | S rho |"
    )
    .unwrap();
    writeln!(md, "|:--|:--|--:|--:|--:|--:|--:|--:|--:|").unwrap();
    for ((_, pair, role), arms) in &a0 {
        let (Some(l), Some(c)) = (arms.get("ic-closed-legacy"), arms.get("ic-closed-coset")) else {
            continue;
        };
        let rho = arms.get("rho-neg").map(|a| a.s()).unwrap_or(0.0);
        writeln!(
            md,
            "| {pair} | {role} | {} | {:.0} | {:.0} | {:.2} | {:.1} | {:.1} | {:.2} |",
            l.fiber_samples.unwrap_or(0),
            l.setup(),
            c.setup(),
            l.setup() / c.setup().max(1.0),
            l.s(),
            c.s(),
            rho
        )
        .unwrap();
    }
    js.insert("a0".into(), table_json(&a0));

    // E1
    let mut e1: Table = BTreeMap::new();
    for name in ["e1_scale_h1568", "e1_scale_h726", "e1_scale_h200"] {
        if let Some(t) = load_session(&format!("{sessions}/{name}"), &roles) {
            e1.extend(t);
        }
    }
    writeln!(md, "\n## E1. Partial fiber `[d]E(F_p)[h]`\n").unwrap();
    writeln!(
        md,
        "Each cell is `|H_d|` : mean `S`; **bold** is the cheapest `d` on that curve.\n"
    )
    .unwrap();
    let mut last_pair = String::new();
    for ((_, pair, role), arms) in &e1 {
        let mut ds: Vec<(u64, &Agg)> = arms
            .iter()
            .filter_map(|(k, a)| {
                k.strip_prefix("ic-closed-d")
                    .and_then(|d| d.parse().ok())
                    .map(|d| (d, a))
            })
            .collect();
        ds.sort_by_key(|(d, _)| *d);
        if *pair != last_pair {
            write!(md, "\n| {pair} | rho |").unwrap();
            for (d, _) in &ds {
                write!(md, " d={d} |").unwrap();
            }
            write!(md, "\n|:--|--:|").unwrap();
            for _ in &ds {
                write!(md, "--:|").unwrap();
            }
            writeln!(md).unwrap();
            last_pair = pair.clone();
        }
        let best = ds.iter().map(|(_, a)| a.s()).fold(f64::INFINITY, f64::min);
        let rho = arms.get("rho-neg").map(|a| a.s()).unwrap_or(0.0);
        write!(md, "| {role} | {rho:.2} |").unwrap();
        for (_, a) in &ds {
            let cell = format!("{} : {:.1}", a.fiber_size.unwrap_or(0), a.s());
            if a.s() == best {
                write!(md, " **{cell}** |").unwrap();
            } else {
                write!(md, " {cell} |").unwrap();
            }
        }
        writeln!(md).unwrap();
    }
    // Matched-size comparison across roles: every (pair, |H|) with more
    // than one role measured.
    writeln!(
        md,
        "\n### E1 at matched subgroup size\n\n| pair | \\|H\\| | role : d : S |\n|:--|--:|:--|"
    )
    .unwrap();
    let mut by_size: BTreeMap<(String, u64), Vec<String>> = BTreeMap::new();
    for ((_, pair, role), arms) in &e1 {
        for (k, a) in arms {
            if let (Some(d), Some(s)) = (k.strip_prefix("ic-closed-d"), a.fiber_size) {
                by_size
                    .entry((pair.clone(), s))
                    .or_default()
                    .push(format!("{role} : {d} : {:.1}", a.s()));
            }
        }
    }
    for ((pair, s), cells) in &by_size {
        if cells.iter().any(|c| c.starts_with("crater"))
            && cells.iter().any(|c| c.starts_with("floor :"))
        {
            writeln!(md, "| {pair} | {s} | {} |", cells.join("; ")).unwrap();
        }
    }
    js.insert("e1".into(), table_json(&e1));

    // E2
    let e2 = load_session(&format!("{sessions}/e2_glv_crater_floor"), &roles).unwrap_or_default();
    writeln!(
        md,
        "\n## E2. Automorphism fold on the crater against the floor\n"
    )
    .unwrap();
    writeln!(
        md,
        "| pair | role | S rho | S glv auto | S glv negation | auto / negation | set-up auto | set-up negation |"
    )
    .unwrap();
    writeln!(md, "|:--|:--|--:|--:|--:|--:|--:|--:|").unwrap();
    for ((_, pair, role), arms) in &e2 {
        let (Some(a), Some(n)) = (arms.get("ic-glv-auto"), arms.get("ic-glv-neg")) else {
            continue;
        };
        let rho = arms.get("rho-neg").map(|a| a.s()).unwrap_or(0.0);
        writeln!(
            md,
            "| {pair} | {role} | {rho:.2} | {:.1} | {:.1} | {:.2} | {:.0} | {:.0} |",
            a.s(),
            n.s(),
            a.s() / n.s().max(1e-9),
            a.setup(),
            n.setup()
        )
        .unwrap();
    }
    js.insert("e2".into(), table_json(&e2));

    std::fs::write(arg(&args, "--out"), &md).expect("write tables");
    std::fs::write(
        arg(&args, "--json"),
        serde_json::to_string_pretty(&Value::Object(js)).unwrap() + "\n",
    )
    .expect("write json");
    print!("{md}");
}
