//! Score `research/pdp_smt_oracle_20261008/PROTOCOL.md` from the frozen
//! JSONL written by `pdp_smt_oracle_bench`.
//!
//! ```bash
//! cargo run --release --example pdp_smt_oracle_score -- \
//!     research/pdp_smt_oracle_20261008/results/n13-m3-ell5.jsonl \
//!     research/pdp_smt_oracle_20261008/results/n19-m3-ell7.jsonl \
//!     research/pdp_smt_oracle_20261008/results/n23-m3-ell8.jsonl
//! ```
//!
//! Prints one Markdown table per `n` (every arm: decided, found, refuted,
//! exhausted, disagreements, spurious; on the cells every arm refuted:
//! native conflicts, the arm's own counters, median wall) and the P1–P4
//! verdicts, computed only from the cell lines.  Counts are compared
//! within an arm; cross-arm columns are decided cells and wall time.
use serde_json::Value;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct Cell {
    arm: String,
    target: u64,
    found: bool,
    refuted: bool,
    exhausted: bool,
    agrees: bool,
    sum_verified: bool,
    spurious: u64,
    conflicts: Option<u64>,
    decisions: Option<u64>,
    wall_ms: f64,
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(v[v.len() / 2])
}

fn fmt_opt(v: Option<u64>) -> String {
    v.map_or_else(|| "n/a".to_string(), |x| x.to_string())
}

fn fmt_ms(v: Option<f64>) -> String {
    v.map_or_else(|| "n/a".to_string(), |x| format!("{x:.0}"))
}

struct Run {
    n: u64,
    arms: Vec<String>,
    cells: Vec<Cell>,
}

fn load(path: &str) -> Run {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut n = 0;
    let mut arms = Vec::new();
    let mut cells = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let v: Value = serde_json::from_str(line).unwrap_or_else(|e| panic!("{path}: {e}"));
        match v["kind"].as_str().unwrap_or("") {
            "pdp_smt_oracle_header" => {
                n = v["n"].as_u64().expect("n");
                arms = v["arms"]
                    .as_array()
                    .expect("arms")
                    .iter()
                    .map(|a| a.as_str().expect("arm").to_string())
                    .collect();
            }
            "pdp_smt_oracle_cell" => cells.push(Cell {
                arm: v["arm"].as_str().expect("arm").to_string(),
                target: v["target_scalar"].as_u64().expect("target"),
                found: v["found"].as_bool().unwrap_or(false),
                refuted: v["refuted"].as_bool().unwrap_or(false),
                exhausted: v["exhausted"].as_bool().unwrap_or(false),
                agrees: v["agrees"].as_bool().unwrap_or(false),
                sum_verified: v["sum_verified"].as_bool().unwrap_or(false),
                spurious: v["spurious"].as_u64().unwrap_or(0),
                conflicts: v["conflicts"].as_u64(),
                decisions: v["decisions"].as_u64(),
                wall_ms: v["wall_ms"].as_f64().unwrap_or(0.0),
            }),
            _ => {}
        }
    }
    Run { n, arms, cells }
}

/// Per-arm figures for one run.
struct ArmRow {
    cells: usize,
    found: usize,
    refuted: usize,
    exhausted: usize,
    disagreements: usize,
    spurious: u64,
    /// On cells refuted by both this arm and the native arm.
    common_refuted: usize,
    native_conflicts_on_common: u64,
    own_conflicts_on_common: Option<u64>,
    own_decisions_on_common: Option<u64>,
    median_wall_refuting: Option<f64>,
    native_median_wall_on_common: Option<f64>,
}

fn arm_rows(run: &Run) -> BTreeMap<String, ArmRow> {
    let native = &run.arms[0];
    let by_target_native: BTreeMap<u64, &Cell> = run
        .cells
        .iter()
        .filter(|c| &c.arm == native)
        .map(|c| (c.target, c))
        .collect();
    let mut rows = BTreeMap::new();
    for arm in &run.arms {
        let mine: Vec<&Cell> = run.cells.iter().filter(|c| &c.arm == arm).collect();
        let common: Vec<(&Cell, &Cell)> = mine
            .iter()
            .filter(|c| c.refuted)
            .filter_map(|c| {
                by_target_native
                    .get(&c.target)
                    .filter(|nc| nc.refuted)
                    .map(|nc| (*c, *nc))
            })
            .collect();
        let sum_opt = |f: &dyn Fn(&Cell) -> Option<u64>| -> Option<u64> {
            common.iter().map(|(c, _)| f(c)).sum::<Option<u64>>()
        };
        rows.insert(
            arm.clone(),
            ArmRow {
                cells: mine.len(),
                found: mine.iter().filter(|c| c.found).count(),
                refuted: mine.iter().filter(|c| c.refuted).count(),
                exhausted: mine.iter().filter(|c| c.exhausted).count(),
                disagreements: mine.iter().filter(|c| !c.agrees || !c.sum_verified).count(),
                spurious: mine.iter().map(|c| c.spurious).sum(),
                common_refuted: common.len(),
                native_conflicts_on_common: common
                    .iter()
                    .map(|(_, nc)| nc.conflicts.unwrap_or(0))
                    .sum(),
                own_conflicts_on_common: if common.is_empty() {
                    None
                } else {
                    sum_opt(&|c| c.conflicts)
                },
                own_decisions_on_common: if common.is_empty() {
                    None
                } else {
                    sum_opt(&|c| c.decisions)
                },
                median_wall_refuting: median(common.iter().map(|(c, _)| c.wall_ms).collect()),
                native_median_wall_on_common: median(
                    common.iter().map(|(_, nc)| nc.wall_ms).collect(),
                ),
            },
        );
    }
    rows
}

fn main() {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    assert!(!paths.is_empty(), "pass one or more result JSONL files");
    let runs: Vec<Run> = paths.iter().map(|p| load(p)).collect();
    let mut p1 = true;
    let mut p2 = true;
    let mut p2_scored = false;
    let mut p3 = true;
    // arm -> (n, wall ratio to native on common refuted cells)
    let mut ratios: BTreeMap<String, Vec<(u64, f64)>> = BTreeMap::new();

    for run in &runs {
        let rows = arm_rows(run);
        let native = &rows[&run.arms[0]];
        println!("### n = {}\n", run.n);
        println!("| arm | cells | found | refuted | exhausted | disagree | spurious | common refuted | native conflicts (common) | own conflicts (common) | own decisions (common) | median wall refuting, ms | native median wall (common), ms | wall ratio |");
        println!("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
        for arm in &run.arms {
            let r = &rows[arm];
            let ratio = match (r.median_wall_refuting, r.native_median_wall_on_common) {
                (Some(a), Some(b)) if b > 0.0 => Some(a / b),
                _ => None,
            };
            println!(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                arm,
                r.cells,
                r.found,
                r.refuted,
                r.exhausted,
                r.disagreements,
                r.spurious,
                r.common_refuted,
                r.native_conflicts_on_common,
                fmt_opt(r.own_conflicts_on_common),
                fmt_opt(r.own_decisions_on_common),
                fmt_ms(r.median_wall_refuting),
                fmt_ms(r.native_median_wall_on_common),
                ratio.map_or_else(|| "n/a".to_string(), |x| format!("{x:.2}")),
            );
            if r.disagreements > 0 || r.spurious > 0 {
                p1 = false;
            }
            if arm != &run.arms[0] {
                if r.found + r.refuted > native.found + native.refuted {
                    p3 = false;
                }
                if let Some(x) = ratio {
                    ratios.entry(arm.clone()).or_default().push((run.n, x));
                }
                if arm.starts_with("z3") && r.common_refuted > 0 {
                    p2_scored = true;
                    match r.own_conflicts_on_common {
                        Some(own) if (own as f64) >= 5.0 * r.native_conflicts_on_common as f64 => {}
                        _ => p2 = false,
                    }
                }
            }
        }
        println!();
    }

    // P4: for arms with a ratio at both the smallest and largest n, the ratio
    // at the largest n is at least the ratio at the smallest.
    let n_min = runs.iter().map(|r| r.n).min().unwrap();
    let n_max = runs.iter().map(|r| r.n).max().unwrap();
    let mut p4 = true;
    let mut p4_scored = false;
    println!("### wall ratio to native on common refuted cells, by n\n");
    for (arm, v) in &ratios {
        let at = |n: u64| v.iter().find(|(m, _)| *m == n).map(|(_, x)| *x);
        let line: Vec<String> = v.iter().map(|(n, x)| format!("n={n}: {x:.2}")).collect();
        println!("- {arm}: {}", line.join(", "));
        if let (Some(a), Some(b)) = (at(n_min), at(n_max)) {
            p4_scored = true;
            if b < a {
                p4 = false;
            }
        }
    }
    println!();
    let verdict = |ok: bool, scored: bool| match (scored, ok) {
        (false, _) => "not scored (no qualifying cells)",
        (true, true) => "PASS",
        (true, false) => "FAIL",
    };
    println!("P1 (correctness): {}", verdict(p1, true));
    println!(
        "P2 (Z3 parity cost ≥ 5× native conflicts on common refuted cells): {}",
        verdict(p2, p2_scored)
    );
    println!(
        "P3 (no SMT arm decides more cells than native): {}",
        verdict(p3, true)
    );
    println!(
        "P4 (wall ratio to native at n={n_max} ≥ ratio at n={n_min}): {}",
        verdict(p4, p4_scored && n_max > n_min)
    );
}
