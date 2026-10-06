//! Readout of `research/ic_gb_ladder_20261003`: the external-engine refutation-degree
//! ladder.  Reads `runs/<tag>/dump/*.catalogue.jsonl` (the draws) and
//! `runs/<tag>/<cell>.jsonl` (one line per Singular process) and prints, per cell and
//! arm, every draw's reading, the cell median, per-arm slopes over `ℓ`, the paired
//! `rr − x4` difference, the excess of `x4` over the sharp first-fall bound 7, and the
//! calibration against the in-tree ladder's readings where the draws are shared.
//!
//!     gb_ladder_analyze research/ic_gb_ladder_20261003/runs/registered \
//!         [research/ic_rr_norm_ladder_20260930/runs/registered]
//!
//! A draw's reading is the least `D` whose process reported `refuted_at = D`; otherwise
//! `≥ D + 1` for the largest `D` that finished unrefuted, or `≥ D` when the process at
//! `D` was killed by its limit (censored, never negative).  `Dp` marks a satisfiable
//! system resolved by pinning every variable at `D`; `triv` marks a system refuted at
//! degree 1, i.e. one that contains a constant equation (the trace condition of the
//! `rr` form on a draw with `V ⊂ ker Tr`), which is excluded from medians and slopes.

use serde_json::Value;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Reading {
    Exact(u32),
    Pinned(u32),
    Trivial,
    AtLeast(u32),
    None,
}

impl Reading {
    fn text(self) -> String {
        match self {
            Reading::Exact(d) => d.to_string(),
            Reading::Pinned(d) => format!("{d}p"),
            Reading::Trivial => "triv".to_string(),
            Reading::AtLeast(d) => format!("≥{d}"),
            Reading::None => "—".to_string(),
        }
    }
}

type Key = (String, String, u32); // cell, arm, draw

fn read_lines(path: &std::path::Path) -> Vec<Value> {
    std::fs::read_to_string(path)
        .map(|s| {
            s.lines()
                .filter(|l| !l.trim().is_empty())
                .map(|l| serde_json::from_str(l).expect("json line"))
                .collect()
        })
        .unwrap_or_default()
}

fn cell_ell(cell: &str) -> u32 {
    cell.split('l').nth(1).unwrap().parse().unwrap()
}

fn median_rule(readings: &[Reading]) -> Reading {
    // majority exact -> median of exact readings; else the least lower bound among bounds
    let mut exact: Vec<u32> = readings
        .iter()
        .filter_map(|r| match r {
            Reading::Exact(d) | Reading::Pinned(d) => Some(*d),
            _ => None,
        })
        .collect();
    let bounds: Vec<u32> = readings
        .iter()
        .filter_map(|r| {
            if let Reading::AtLeast(d) = r {
                Some(*d)
            } else {
                None
            }
        })
        .collect();
    if exact.is_empty() && bounds.is_empty() {
        return Reading::None;
    }
    if exact.len() * 2 > exact.len() + bounds.len() {
        exact.sort_unstable();
        let m = exact.len();
        let med = if m % 2 == 1 {
            exact[m / 2] as f64
        } else {
            (exact[m / 2 - 1] + exact[m / 2]) as f64 / 2.0
        };
        Reading::Exact(med.round() as u32)
    } else {
        // every reading is at least this
        let all_min = readings
            .iter()
            .filter_map(|r| match r {
                Reading::Exact(d) | Reading::Pinned(d) | Reading::AtLeast(d) => Some(*d),
                Reading::None | Reading::Trivial => None,
            })
            .min()
            .unwrap();
        Reading::AtLeast(all_min)
    }
}

fn ols(points: &[(f64, f64)]) -> Option<f64> {
    if points.len() < 2 {
        return None;
    }
    let n = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / n;
    let my = points.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = points.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx == 0.0 {
        return None;
    }
    Some(points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>() / sxx)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = std::path::Path::new(args.get(1).expect("runs dir"));
    let intree = args.get(2).map(std::path::PathBuf::from);
    // draws
    let mut draws: BTreeMap<Key, (String, u64)> = BTreeMap::new(); // -> (dump, roots)
    let mut cells: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(dir.join("dump")).expect("dump dir") {
        let p = entry.unwrap().path();
        if !p.to_string_lossy().ends_with(".catalogue.jsonl") {
            continue;
        }
        for v in read_lines(&p) {
            let cell = v["cell"].as_str().unwrap().to_string();
            if !cells.contains(&cell) {
                cells.push(cell.clone());
            }
            let roots = v["roots"].as_u64().unwrap_or(0);
            draws.insert(
                (
                    cell,
                    v["arm"].as_str().unwrap().to_string(),
                    v["draw"].as_u64().unwrap() as u32,
                ),
                (v["dump"].as_str().unwrap().to_string(), roots),
            );
        }
    }
    cells.sort_by_key(|c| (cell_ell(c), c.clone()));
    // process lines -> readings and cost
    let mut lines: BTreeMap<Key, Vec<Value>> = BTreeMap::new();
    for cell in &cells {
        for v in read_lines(&dir.join(format!("{cell}.jsonl"))) {
            lines
                .entry((
                    cell.clone(),
                    v["arm"].as_str().unwrap().to_string(),
                    v["draw"].as_u64().unwrap() as u32,
                ))
                .or_default()
                .push(v);
        }
    }
    let mut readings: BTreeMap<Key, (Reading, f64)> = BTreeMap::new(); // reading, wall of the deciding/last run
    for (key, (_, roots)) in &draws {
        if *roots > 0 {
            continue;
        }
        // the last line per degree wins (a phase-2 retry supersedes a phase-1 memory kill)
        let mut by_d: BTreeMap<u32, Value> = BTreeMap::new();
        for l in lines.get(key).cloned().unwrap_or_default() {
            by_d.insert(l["dmax"].as_u64().unwrap() as u32, l);
        }
        let ls: Vec<Value> = by_d.into_values().collect();
        let mut r = Reading::None;
        let mut wall = 0.0;
        let mut max_unref = 0u32;
        for l in &ls {
            let d = l["dmax"].as_u64().unwrap() as u32;
            wall += l["wall"].as_f64().unwrap_or(0.0);
            if l["kind"] == "killed" {
                r = Reading::AtLeast(d);
                break;
            }
            let ra = l["refuted_at"].as_u64().unwrap() as u32;
            if ra == 1 {
                r = Reading::Trivial;
                break;
            }
            if ra > 0 {
                r = Reading::Exact(ra);
                break;
            }
            if l["pinned"].as_u64().unwrap_or(0) > 0 && l["pinned"] == l["n_vars"] {
                r = Reading::Pinned(d);
                break;
            }
            max_unref = max_unref.max(d);
            r = Reading::AtLeast(d + 1);
        }
        let _ = max_unref;
        readings.insert(key.clone(), (r, wall));
    }
    // in-tree readings (last final line per draw, arm)
    let mut intree_r: BTreeMap<Key, Reading> = BTreeMap::new();
    if let Some(it) = &intree {
        for cell in &cells {
            for v in read_lines(&it.join(format!("{cell}.jsonl"))) {
                let o = &v["outcome"];
                if o["partial"].as_bool().unwrap_or(false) {
                    continue;
                }
                let r = match o["kind"].as_str().unwrap_or("") {
                    "resolved" => Reading::Exact(o["degree"].as_u64().unwrap() as u32),
                    "at_least" => Reading::AtLeast(o["degree"].as_u64().unwrap() as u32),
                    _ => Reading::None,
                };
                intree_r.insert(
                    (
                        cell.clone(),
                        v["arm"].as_str().unwrap().to_string(),
                        v["draw"].as_u64().unwrap() as u32,
                    ),
                    r,
                );
            }
        }
    }
    println!("# External-engine refutation-degree ladder: readout\n");
    println!("## Per draw (reading; Singular wall seconds summed over the draw's degree runs; in-tree reading where shared)\n");
    println!("| cell | arm | draw | n_vars | reading | wall s | in-tree |");
    println!("|:--|:--|--:|--:|:--|--:|:--|");
    let mut per_cell: BTreeMap<(String, String), Vec<Reading>> = BTreeMap::new();
    for (key, (r, wall)) in &readings {
        let nv = {
            let (dump, _) = &draws[key];
            // n_vars from the dump's ring line
            std::fs::read_to_string(dir.join("dump").join(dump))
                .ok()
                .and_then(|s| {
                    s.lines()
                        .next()
                        .and_then(|l| l.split("x(1..").nth(1))
                        .and_then(|t| t.split(')').next())
                        .and_then(|n| n.parse::<u32>().ok())
                })
                .unwrap_or(0)
        };
        println!(
            "| {} | {} | {} | {} | {} | {:.0} | {} |",
            key.0,
            key.1,
            key.2,
            nv,
            r.text(),
            wall,
            intree_r.get(key).map_or(String::new(), |x| x.text())
        );
        per_cell
            .entry((key.0.clone(), key.1.clone()))
            .or_default()
            .push(*r);
    }
    println!("\n## Per cell\n");
    println!("| cell | arm | readings | median | x4 − 7 |");
    println!("|:--|:--|:--|:--|:--|");
    let mut medians: BTreeMap<(String, String), Reading> = BTreeMap::new();
    for ((cell, arm), rs) in &per_cell {
        let m = median_rule(rs);
        medians.insert((cell.clone(), arm.clone()), m);
        let excess = match (arm.as_str(), m) {
            ("x4", Reading::Exact(d)) => format!("{}", d as i64 - 7),
            ("x4", Reading::AtLeast(d)) => format!("≥{}", d as i64 - 7),
            _ => String::new(),
        };
        println!(
            "| {cell} | {arm} | {} | {} | {excess} |",
            rs.iter().map(|r| r.text()).collect::<Vec<_>>().join(" "),
            m.text()
        );
    }
    println!("\n## Slopes (OLS of exact cell medians against ℓ, per curve and arm)\n");
    let mut curves: Vec<String> = cells
        .iter()
        .map(|c| c.split('l').next().unwrap().to_string())
        .collect();
    curves.sort();
    curves.dedup();
    for arm in ["rr", "x4", "ctrl"] {
        // pool K1 cells across n: the curve label is K<a>n<n>; slope over all K1 cells with exact medians
        for a in ["K0", "K1"] {
            let pts: Vec<(f64, f64)> = medians
                .iter()
                .filter(|((c, ar), m)| {
                    ar == arm && c.starts_with(a) && matches!(m, Reading::Exact(_))
                })
                .map(|((c, _), m)| {
                    (
                        cell_ell(c) as f64,
                        if let Reading::Exact(d) = m {
                            *d as f64
                        } else {
                            0.0
                        },
                    )
                })
                .collect();
            if let Some(s) = ols(&pts) {
                println!(
                    "- {arm} on {a}: slope {s:.3} over {} cells {:?}",
                    pts.len(),
                    pts
                );
            }
        }
    }
    println!("\n## Paired rr − x4 on draws exact on both arms\n");
    let mut diffs: Vec<i64> = Vec::new();
    for ((cell, arm, draw), (r, _)) in &readings {
        if arm != "rr" {
            continue;
        }
        if let (Reading::Exact(a), Some((Reading::Exact(b), _))) =
            (r, readings.get(&(cell.clone(), "x4".to_string(), *draw)))
        {
            diffs.push(*a as i64 - *b as i64);
        }
    }
    if diffs.is_empty() {
        println!("none");
    } else {
        let mean = diffs.iter().sum::<i64>() as f64 / diffs.len() as f64;
        println!(
            "{} draws: mean {mean:.2}, min {}, max {}",
            diffs.len(),
            diffs.iter().min().unwrap(),
            diffs.iter().max().unwrap()
        );
    }
    println!("\n## Calibration against the in-tree ladder (shared draws, both exact)\n");
    let (mut agree, mut differ) = (0, Vec::new());
    for (key, (r, _)) in &readings {
        if let (Reading::Exact(a), Some(Reading::Exact(b))) = (r, intree_r.get(key)) {
            if a == b {
                agree += 1;
            } else {
                differ.push(format!(
                    "{} {} d{}: external {a}, in-tree {b}",
                    key.0, key.1, key.2
                ));
            }
        }
    }
    println!("agree on {agree} draws; differ on {}:", differ.len());
    for d in differ {
        println!("- {d}");
    }
    let _ = curves;
}
