//! # Tables for research/isogeny_conductor_gap_round3_20261011
//!
//! Reads the E4 records (`e4/records.jsonl`, written by
//! `examples/isogeny_gap_two_torsion_f4.rs`) and the E5 `ecbench` sessions,
//! and prints the README's tables and a JSON twin.
//!
//! * **E4 shape**: for each factor-base size `B`, subgroup order `|H|` and
//!   target kind, the range over every curve and target of the solving
//!   degree, widest matrix and `F_p` row-reduction multiplications.
//! * **E4 gain**: the direct system at `B` against the quotient by `|H|`,
//!   as ratios of the random-target means.
//! * **E4 gap**: on every pair whose crater has rational `E[2]`, the
//!   crater's `E[2]` beside the floor's and the floor model's `C4`, at each
//!   `B`.
//! * **E5**: per session and arm, the verified runs, mean `S`, solver
//!   operations and maximum solving degree, with the combined arm's solver
//!   work over the blind arm's.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_round3_tables -- \
//!   --dir research/isogeny_conductor_gap_round3_20261011 \
//!   --out research/isogeny_conductor_gap_round3_20261011/TABLES.md \
//!   --json research/isogeny_conductor_gap_round3_20261011/tables.json
//! ```

use std::collections::BTreeMap;
use std::fmt::Write as _;

use serde_json::{json, Value};

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|a| a == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    args[i + 1].clone()
}

#[derive(Default, Clone)]
struct Range {
    n: u64,
    min: f64,
    max: f64,
    sum: f64,
}

impl Range {
    fn add(&mut self, v: f64) {
        if self.n == 0 {
            self.min = v;
            self.max = v;
        } else {
            self.min = self.min.min(v);
            self.max = self.max.max(v);
        }
        self.n += 1;
        self.sum += v;
    }
    fn mean(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.sum / self.n as f64
        }
    }
    fn cell(&self) -> String {
        if self.min == self.max {
            format!("{:.0}", self.min)
        } else {
            format!("{:.0}–{:.0}", self.min, self.max)
        }
    }
    fn json(&self) -> Value {
        json!({"n": self.n, "min": self.min, "max": self.max, "mean": self.mean()})
    }
}

#[derive(Default, Clone)]
struct Shape {
    systems: u64,
    base_degree: Range,
    degree: Range,
    cols: Range,
    ops: Range,
    ms: Range,
}

impl Shape {
    fn add(&mut self, r: &Value) {
        self.systems += 1;
        self.base_degree.add(r["base_degree"].as_f64().unwrap());
        self.degree.add(r["solving_degree_max"].as_f64().unwrap());
        self.cols.add(r["max_cols_to_solution"].as_f64().unwrap());
        self.ops.add(r["field_ops"].as_f64().unwrap());
        self.ms.add(r["ms"].as_f64().unwrap());
    }
    fn json(&self) -> Value {
        json!({"systems": self.systems, "base_degree": self.base_degree.json(),
               "solving_degree_max": self.degree.json(), "max_cols_to_solution": self.cols.json(),
               "field_ops": self.ops.json(), "ms": self.ms.json()})
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = arg(&args, "--dir");
    let mut md = String::new();
    let mut js = serde_json::Map::new();

    // ── E4 ──
    let text = std::fs::read_to_string(format!("{dir}/e4/records.jsonl")).expect("E4 records");
    let recs: Vec<Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).expect("record"))
        .filter(|r: &Value| r["kind"] == "system")
        .collect();
    let failures = recs
        .iter()
        .filter(|r| {
            r["f4_agrees_with_enumeration"] != true || r["planted_annihilates_basis"] == false
        })
        .count();
    let timeouts = recs.iter().filter(|r| r["timed_out"] == true).count();
    let curves: std::collections::BTreeSet<String> = recs
        .iter()
        .map(|r| r["curve_key"].as_str().unwrap().to_string())
        .collect();
    let planted = recs
        .iter()
        .filter(|r| r["target_kind"] == "planted")
        .count();
    writeln!(
        md,
        "## E4. Two-torsion quotients of the S3 system\n\n{} systems on {} curves; {} planted targets. F4 and enumeration disagree, or a planted solution fails the basis, on {} systems; {} timed out.\n",
        recs.len(),
        curves.len(),
        planted,
        failures,
        timeouts
    )
    .unwrap();

    // Shape by (B, |H|, target kind).
    let mut shape: BTreeMap<(u64, u64, String), Shape> = BTreeMap::new();
    for r in &recs {
        shape
            .entry((
                r["base_size"].as_u64().unwrap(),
                r["subgroup_order"].as_u64().unwrap(),
                r["target_kind"].as_str().unwrap().to_string(),
            ))
            .or_default()
            .add(r);
    }
    writeln!(
        md,
        "### Shape: every curve and role, by base size B and subgroup order |H|\n"
    )
    .unwrap();
    writeln!(
        md,
        "| B | \\|H\\| | targets | systems | base polynomial degree | solving degree | widest matrix (cols) | F_p row-reduction mults |"
    )
    .unwrap();
    writeln!(md, "|--:|--:|:--|--:|--:|--:|--:|--:|").unwrap();
    for ((b, h, kind), s) in &shape {
        writeln!(
            md,
            "| {b} | {h} | {kind} | {} | {} | {} | {} | {} |",
            s.systems,
            s.base_degree.cell(),
            s.degree.cell(),
            s.cols.cell(),
            s.ops.cell()
        )
        .unwrap();
    }
    js.insert(
        "e4_shape".into(),
        Value::Array(
            shape
                .iter()
                .map(|((b, h, k), s)| json!({"base_size": b, "subgroup_order": h, "target_kind": k, "shape": s.json()}))
                .collect(),
        ),
    );

    // Gain: direct at B over quotient by |H| at B (random targets).
    writeln!(
        md,
        "\n### Gain at matched B: direct system over the quotient by |H| (random targets, means)\n\n| B | \\|H\\| | solving degree direct ÷ quotient | widest matrix ÷ | row-reduction mults ÷ |\n|--:|--:|--:|--:|--:|"
    )
    .unwrap();
    let mut gain = Vec::new();
    for ((b, h, kind), s) in &shape {
        if kind != "random" || *h == 1 {
            continue;
        }
        let Some(d) = shape.get(&(*b, 1, "random".to_string())) else {
            continue;
        };
        let (rd, rc, ro) = (
            d.degree.mean() / s.degree.mean(),
            d.cols.mean() / s.cols.mean(),
            d.ops.mean() / s.ops.mean(),
        );
        writeln!(md, "| {b} | {h} | {rd:.2} | {rc:.2} | {ro:.2} |").unwrap();
        gain.push(json!({"base_size": b, "subgroup_order": h, "degree_ratio": rd, "cols_ratio": rc, "ops_ratio": ro}));
    }
    js.insert("e4_gain".into(), Value::Array(gain));

    // Gap: crater E[2] against floor / floor-iso C4.
    let mut by: BTreeMap<(String, u64, String, String), Shape> = BTreeMap::new();
    let mut rank: BTreeMap<String, u64> = BTreeMap::new();
    for r in &recs {
        if r["target_kind"] != "random" {
            continue;
        }
        let pair = r["pair"].as_str().unwrap().to_string();
        let role = r["role"].as_str().unwrap().to_string();
        rank.insert(format!("{pair}|{role}"), r["two_rank"].as_u64().unwrap());
        by.entry((
            pair,
            r["base_size"].as_u64().unwrap(),
            role,
            r["subgroup"].as_str().unwrap().to_string(),
        ))
        .or_default()
        .add(r);
    }
    writeln!(
        md,
        "\n### Gap: the crater's E[2] against the floor's C4, matched B (random targets, means)\n\n| pair | B | crater E[2]: degree / cols / mults | floor C4 | floor model C4 | crater C4 | floor ÷ crater mults |\n|:--|--:|:--|:--|:--|:--|--:|"
    )
    .unwrap();
    let mut gap = Vec::new();
    let pairs: std::collections::BTreeSet<(String, u64)> =
        by.keys().map(|(p, b, _, _)| (p.clone(), *b)).collect();
    for (pair, b) in pairs {
        if rank.get(&format!("{pair}|crater")) != Some(&2) {
            continue;
        }
        let get = |role: &str, h: &str| by.get(&(pair.clone(), b, role.to_string(), h.to_string()));
        let cell = |s: Option<&Shape>| {
            s.map(|s| {
                format!(
                    "{:.0} / {:.0} / {:.0}",
                    s.degree.mean(),
                    s.cols.mean(),
                    s.ops.mean()
                )
            })
            .unwrap_or("–".into())
        };
        let (c, f, fi, cc) = (
            get("crater", "E[2]"),
            get("floor", "C4"),
            get("floor-iso", "C4"),
            get("crater", "C4"),
        );
        let ratio = match (c, f) {
            (Some(c), Some(f)) => format!("{:.3}", f.ops.mean() / c.ops.mean()),
            _ => "–".into(),
        };
        writeln!(
            md,
            "| {pair} | {b} | {} | {} | {} | {} | {ratio} |",
            cell(c),
            cell(f),
            cell(fi),
            cell(cc)
        )
        .unwrap();
        gap.push(json!({"pair": pair, "base_size": b,
            "crater_E2": c.map(|s| s.json()), "floor_C4": f.map(|s| s.json()),
            "floor_iso_C4": fi.map(|s| s.json()), "crater_C4": cc.map(|s| s.json())}));
    }
    js.insert("e4_gap".into(), Value::Array(gap));
    let mut ranks = BTreeMap::new();
    for (k, v) in &rank {
        ranks.insert(k.clone(), *v);
    }
    js.insert("e4_two_rank".into(), json!(ranks));

    // ── E5 ──
    writeln!(md, "\n## E5. Combined Weil descent, audited sessions\n").unwrap();
    let mut e5 = serde_json::Map::new();
    let mut names: Vec<String> = std::fs::read_dir(format!("{dir}/sessions"))
        .map(|d| {
            d.filter_map(|e| e.ok())
                .filter(|e| e.path().join("records.jsonl").exists())
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    if names.is_empty() {
        writeln!(md, "No completed E5 session.").unwrap();
    }
    for name in names {
        let text = std::fs::read_to_string(format!("{dir}/sessions/{name}/records.jsonl")).unwrap();
        #[derive(Default)]
        struct Arm {
            runs: u64,
            verified: u64,
            s: Range,
            ops: Range,
            deg_max: u64,
            wall_s: Range,
            statuses: BTreeMap<String, u64>,
        }
        let mut arms: BTreeMap<String, Arm> = BTreeMap::new();
        for line in text.lines() {
            let r: Value = serde_json::from_str(line).unwrap();
            if r["warmup"] == true {
                continue;
            }
            let a = arms
                .entry(r["arm"].as_str().unwrap().to_string())
                .or_default();
            a.runs += 1;
            *a.statuses
                .entry(r["outcome"]["status"].as_str().unwrap_or("?").to_string())
                .or_default() += 1;
            if r["outcome"]["status"] != "verified" {
                continue;
            }
            a.verified += 1;
            a.s.add(r["cost"]["s"].as_f64().unwrap_or(0.0));
            a.wall_s
                .add(r["online"]["wall_ns"].as_f64().unwrap_or(0.0) / 1e9);
            if r["solver"].is_object() {
                a.ops.add(r["solver"]["ops"].as_f64().unwrap_or(0.0));
                a.deg_max = a
                    .deg_max
                    .max(r["solver"]["solving_degree_max"].as_u64().unwrap_or(0));
            }
        }
        writeln!(
            md,
            "### `{name}`\n\n| arm | runs | verified | S | solver ops | max solving degree | online wall (s) |\n|:--|--:|--:|--:|--:|--:|--:|"
        )
        .unwrap();
        let blind_ops = arms.get("f4-blind").map(|a| a.ops.mean());
        let mut sj = serde_json::Map::new();
        for (arm, a) in &arms {
            writeln!(
                md,
                "| {arm} | {} | {} | {} | {} | {} | {} |",
                a.runs,
                a.verified,
                if a.verified > 0 {
                    format!("{:.2}", a.s.mean())
                } else {
                    "–".into()
                },
                if a.ops.n > 0 {
                    a.ops.cell()
                } else {
                    "–".into()
                },
                if a.deg_max > 0 {
                    a.deg_max.to_string()
                } else {
                    "–".into()
                },
                if a.verified > 0 {
                    format!("{:.1}", a.wall_s.mean())
                } else {
                    "–".into()
                },
            )
            .unwrap();
            sj.insert(arm.clone(), json!({"runs": a.runs, "verified": a.verified, "statuses": a.statuses,
                "s": a.s.json(), "solver_ops": a.ops.json(), "solving_degree_max": a.deg_max,
                "online_wall_s": a.wall_s.json(),
                "solver_ops_over_blind": blind_ops.filter(|b| *b > 0.0 && a.ops.n > 0).map(|b| a.ops.mean() / b)}));
        }
        if let (Some(b), Some(c)) = (blind_ops, arms.get("f4-combined")) {
            if c.ops.n > 0 && b > 0.0 {
                writeln!(
                    md,
                    "\nCombined over blind solver operations: {:.1}×.",
                    c.ops.mean() / b
                )
                .unwrap();
            }
        }
        writeln!(md).unwrap();
        e5.insert(name, Value::Object(sj));
    }
    js.insert("e5".into(), Value::Object(e5));

    std::fs::write(arg(&args, "--out"), &md).expect("write tables");
    std::fs::write(
        arg(&args, "--json"),
        serde_json::to_string_pretty(&Value::Object(js)).unwrap() + "\n",
    )
    .expect("write json");
    print!("{md}");
}
