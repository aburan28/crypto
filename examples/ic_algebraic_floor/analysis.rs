//! `analyze`: cost per relation per cell, growth slopes in `n` with
//! bootstrap intervals, and the frozen decision rule
//! (`research/ic_m3_algebraic_vs_mitm_20261014/PROTOCOL.md` §5).
//!
//! Cost per relation of an arm on a cell is the stratified estimator
//! `(p·c_yes + (1 − p)·c_no) / p`: `p` is the decomposable fraction of the
//! cell's uniform target pool (from the exact floor), `c_yes` and `c_no`
//! the arm's mean wall time on the decomposable and non-decomposable
//! targets it ran.  A censored target enters at the wall it was stopped at
//! (at least the budget), so a cell with any censored target, or one the
//! early-stop rule cut short, is a **lower bound**, never a measurement.
//! The floor's cost per relation is its mean decision-mode wall over the
//! pool divided by `p`, plus the pair-table build spread over the `|F|`
//! relations a solve needs (`amortized`, primary) or without it (`online`).

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Clone, Debug)]
struct Row {
    yes: bool,
    wall: f64,
    censored: bool,
    xors: f64,
    conflicts: f64,
    nodes: f64,
    unsupported: bool,
}

#[derive(Clone, Debug)]
struct ArmCell {
    n: u32,
    d: u32,
    rows: Vec<Row>,
    stopped: bool,
    skipped: bool,
    disagreements: usize,
    bad_witness: usize,
}

impl ArmCell {
    fn status(&self) -> &'static str {
        if self.skipped {
            "not_run"
        } else if self.rows.iter().any(|r| r.unsupported) {
            "unsupported"
        } else if self.disagreements > 0 || self.bad_witness > 0 {
            "inexact"
        } else if self.rows.is_empty() {
            "not_run"
        } else if self.stopped || self.rows.iter().any(|r| r.censored) {
            "lower_bound"
        } else {
            "measured"
        }
    }
}

#[derive(Clone, Debug)]
struct Truth {
    n: u32,
    d: u32,
    pool_yes: Vec<bool>,
    pool_ns: Vec<f64>,
    pool_adds: Vec<f64>,
    build_ns: f64,
    build_adds: f64,
    columns: f64,
    fb_points: u64,
    total_decomps: u64,
    exhaustive_checked: usize,
    exhaustive_equal: usize,
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn resample<T: Clone>(rng: &mut StdRng, v: &[T]) -> Vec<T> {
    (0..v.len())
        .map(|_| v[rng.gen_range(0..v.len())].clone())
        .collect()
}

/// Floor cost per relation, wall ns (amortized, online) and adds.
/// Sensitivity switch (`--floor online`, post-hoc, not the frozen primary):
/// read the floor without its table build.
static FLOOR_ONLINE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn floor_cost(t: &Truth, yes: &[bool], ns: &[f64], adds: &[f64]) -> Option<(f64, f64, f64)> {
    let p = yes.iter().filter(|&&y| y).count() as f64 / yes.len() as f64;
    if p == 0.0 {
        return None;
    }
    let online = mean(ns) / p;
    let build = !FLOOR_ONLINE.load(std::sync::atomic::Ordering::Relaxed);
    let amort = online + if build { t.build_ns / t.columns } else { 0.0 };
    let adds_amort = mean(adds) / p + if build { t.build_adds / t.columns } else { 0.0 };
    Some((amort, online, adds_amort))
}

fn arm_cost(p: f64, rows: &[Row], key: fn(&Row) -> f64) -> Option<f64> {
    if p == 0.0 {
        return None;
    }
    let y: Vec<f64> = rows.iter().filter(|r| r.yes).map(key).collect();
    let n: Vec<f64> = rows.iter().filter(|r| !r.yes).map(key).collect();
    Some((p * mean(&y) + (1.0 - p) * mean(&n)) / p)
}

fn ols_slope(xs: &[f64], ys: &[f64]) -> f64 {
    let mx = mean(xs);
    let my = mean(ys);
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    sxy / sxx
}

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() - 1) as f64 * q).round() as usize;
    v[i]
}

fn load(
    dirs: &[PathBuf],
) -> (
    BTreeMap<(usize, u32), Truth>,
    BTreeMap<(usize, String), Vec<ArmCell>>,
) {
    let mut truths = BTreeMap::new();
    let mut arms: BTreeMap<(usize, String), Vec<ArmCell>> = BTreeMap::new();
    for dir in dirs {
        let mut entries: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .flatten()
            .map(|e| e.path())
            .collect();
        entries.sort();
        for path in entries {
            let name = path.file_name().unwrap().to_string_lossy().to_string();
            if name.starts_with("truth-") && name.ends_with(".json") {
                let v: Value =
                    serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
                let rows = v["rows"].as_array().unwrap();
                let ex = v["exhaustive_check"].as_array().unwrap();
                let t = Truth {
                    n: v["n"].as_u64().unwrap() as u32,
                    d: v["d"].as_u64().unwrap() as u32,
                    pool_yes: rows.iter().map(|r| r["yes"].as_bool().unwrap()).collect(),
                    pool_ns: rows.iter().map(|r| r["ns"].as_f64().unwrap()).collect(),
                    pool_adds: rows.iter().map(|r| r["adds"].as_f64().unwrap()).collect(),
                    build_ns: v["table_build_ns"].as_f64().unwrap(),
                    build_adds: v["table_build_adds"].as_f64().unwrap(),
                    columns: v["fb_columns"].as_f64().unwrap(),
                    fb_points: v["fb_points"].as_u64().unwrap(),
                    total_decomps: v["total_decomps"].as_u64().unwrap(),
                    exhaustive_checked: ex.len(),
                    exhaustive_equal: ex
                        .iter()
                        .filter(|e| e["equal"].as_bool() == Some(true))
                        .count(),
                };
                truths.insert((v["m"].as_u64().unwrap() as usize, t.n), t);
            } else if name.starts_with("arm-") && name.ends_with(".jsonl") {
                let text = std::fs::read_to_string(&path).unwrap();
                let mut cell: Option<ArmCell> = None;
                let mut key: Option<(usize, String)> = None;
                for line in text.lines().filter(|l| !l.trim().is_empty()) {
                    let v: Value = serde_json::from_str(line).unwrap();
                    let m = v["m"].as_u64().unwrap() as usize;
                    let arm = v["arm"].as_str().unwrap().to_string();
                    key = Some((m, arm));
                    let c = cell.get_or_insert(ArmCell {
                        n: v["n"].as_u64().unwrap() as u32,
                        d: v["d"].as_u64().unwrap() as u32,
                        rows: Vec::new(),
                        stopped: false,
                        skipped: false,
                        disagreements: 0,
                        bad_witness: 0,
                    });
                    match v["cmd"].as_str().unwrap() {
                        "skipped" => c.skipped = true,
                        "stopped" => c.stopped = true,
                        "solve" => {
                            let budget = v["budget_secs"].as_f64().unwrap() * 1e9;
                            let censored = v["censored"].as_bool().unwrap()
                                && !v["unsupported"].as_bool().unwrap_or(false);
                            let mut wall = v["wall_ns"].as_f64().unwrap();
                            if censored {
                                wall = wall.max(budget);
                            }
                            if v["agrees"].as_bool() == Some(false) {
                                c.disagreements += 1;
                            }
                            if v["witness_ok"].as_bool() == Some(false) {
                                c.bad_witness += 1;
                            }
                            let ex = &v["extra"];
                            c.rows.push(Row {
                                yes: v["truth_yes"].as_bool().unwrap(),
                                wall,
                                censored,
                                xors: v["word_xors"].as_f64().unwrap_or(0.0)
                                    + ex["split_word_xors"].as_f64().unwrap_or(0.0),
                                conflicts: ex["conflicts"].as_f64().unwrap_or(0.0),
                                nodes: ex["nodes"].as_f64().unwrap_or(0.0),
                                unsupported: v["unsupported"].as_bool().unwrap_or(false),
                            });
                        }
                        other => panic!("unknown record {other}"),
                    }
                }
                if let (Some(c), Some(k)) = (cell, key) {
                    arms.entry(k).or_default().push(c);
                }
            }
        }
    }
    for v in arms.values_mut() {
        v.sort_by_key(|c| c.n);
    }
    (truths, arms)
}

pub fn cmd_analyze(args: &[String]) {
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let dirs: Vec<PathBuf> = flag("--dir")
        .expect("--dir DIR[,DIR]")
        .split(',')
        .map(PathBuf::from)
        .collect();
    let boot: usize = flag("--boot").and_then(|v| v.parse().ok()).unwrap_or(2000);
    let seed: u64 = flag("--seed")
        .and_then(|v| v.parse().ok())
        .unwrap_or(20261014);
    FLOOR_ONLINE.store(
        flag("--floor").as_deref() == Some("online"),
        std::sync::atomic::Ordering::Relaxed,
    );
    // Sensitivity switch (`--exclude-n 11,…`, post-hoc): drop whole cells.
    let exclude: Vec<u32> = flag("--exclude-n")
        .map(|s| s.split(',').map(|v| v.parse().unwrap()).collect())
        .unwrap_or_default();
    let out = PathBuf::from(flag("--out").expect("--out ANALYSIS.json"));
    let (truths, arms) = load(&dirs);
    let truths: BTreeMap<(usize, u32), Truth> = truths
        .into_iter()
        .filter(|(k, _)| !exclude.contains(&k.1))
        .collect();
    let arms: BTreeMap<(usize, String), Vec<ArmCell>> = arms
        .into_iter()
        .map(|(k, v)| {
            (
                k,
                v.into_iter().filter(|c| !exclude.contains(&c.n)).collect(),
            )
        })
        .collect();
    let mut rng = StdRng::seed_from_u64(seed);

    let ms: Vec<usize> = truths
        .keys()
        .map(|k| k.0)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut report = json!({"protocol": "research/ic_m3_algebraic_vs_mitm_20261014/PROTOCOL.md",
        "boot": boot, "seed": seed, "by_m": {}});
    let mut md = String::new();
    md.push_str("# ANALYSIS (written by `ic_algebraic_floor analyze`)\n\n");
    md.push_str(
        "Cost per relation in wall-clock ns (log2), one core, same binary and host. \
        `measured` = every target finished; `lower_bound` = some target censored at its budget or \
        the cell stopped early; intervals are 95% bootstrap percentiles.\n\n",
    );
    let mut csv = String::from("m,arm,n,d,status,log2_cost_ns,lo,hi,log2_counted,counted_unit\n");

    for &m in &ms {
        let cells: Vec<&Truth> = truths
            .iter()
            .filter(|(k, _)| k.0 == m)
            .map(|(_, t)| t)
            .collect();
        // Point estimates.
        let mut cell_rows = Vec::new();
        for t in &cells {
            let p = t.pool_yes.iter().filter(|&&y| y).count() as f64 / t.pool_yes.len() as f64;
            let fc = floor_cost(t, &t.pool_yes, &t.pool_ns, &t.pool_adds);
            cell_rows.push(json!({"n": t.n, "d": t.d, "fb_points": t.fb_points, "columns": t.columns,
                "pool": t.pool_yes.len(), "p_decomposable": p,
                "decomps_per_target": t.total_decomps as f64 / t.pool_yes.len() as f64,
                "floor_amortized_ns": fc.map(|c| c.0), "floor_online_ns": fc.map(|c| c.1),
                "floor_amortized_adds": fc.map(|c| c.2),
                "table_build_ns": t.build_ns,
                "exhaustive_checked": t.exhaustive_checked, "exhaustive_equal": t.exhaustive_equal}));
        }
        // Bootstrap replicates: per cell, floor (amortized, online) and per arm.
        let arm_names: Vec<String> = arms
            .keys()
            .filter(|k| k.0 == m)
            .map(|k| k.1.clone())
            .collect();
        let mut floor_reps: Vec<Vec<Option<(f64, f64)>>> =
            vec![Vec::with_capacity(boot); cells.len()];
        let mut arm_reps: BTreeMap<(String, u32), Vec<Option<f64>>> = BTreeMap::new();
        for _ in 0..boot {
            for (ci, t) in cells.iter().enumerate() {
                let idx: Vec<usize> = (0..t.pool_yes.len())
                    .map(|_| rng.gen_range(0..t.pool_yes.len()))
                    .collect();
                let yes: Vec<bool> = idx.iter().map(|&i| t.pool_yes[i]).collect();
                let ns: Vec<f64> = idx.iter().map(|&i| t.pool_ns[i]).collect();
                let adds: Vec<f64> = idx.iter().map(|&i| t.pool_adds[i]).collect();
                let fc = floor_cost(t, &yes, &ns, &adds);
                floor_reps[ci].push(fc.map(|c| (c.0, c.1)));
                let p = yes.iter().filter(|&&y| y).count() as f64 / yes.len() as f64;
                for arm in &arm_names {
                    if let Some(ac) = arms[&(m, arm.clone())].iter().find(|c| c.n == t.n) {
                        if ac.rows.is_empty() {
                            continue;
                        }
                        let ys: Vec<Row> = ac.rows.iter().filter(|r| r.yes).cloned().collect();
                        let no: Vec<Row> = ac.rows.iter().filter(|r| !r.yes).cloned().collect();
                        let mut rs = resample(&mut rng, &ys);
                        rs.extend(resample(&mut rng, &no));
                        arm_reps
                            .entry((arm.clone(), t.n))
                            .or_default()
                            .push(arm_cost(p, &rs, |r| r.wall));
                    }
                }
            }
        }
        let ns_of: Vec<f64> = cells.iter().map(|t| t.n as f64).collect();
        let log2 = |x: f64| x.log2();
        // Floor slopes (all cells with p > 0).
        let floor_ok: Vec<usize> = (0..cells.len())
            .filter(|&i| cell_rows[i]["floor_amortized_ns"].is_f64())
            .collect();
        let fit = |sel: &[usize], ys: &dyn Fn(usize) -> f64| -> f64 {
            let xs: Vec<f64> = sel.iter().map(|&i| ns_of[i]).collect();
            let y: Vec<f64> = sel.iter().map(|&i| ys(i)).collect();
            ols_slope(&xs, &y)
        };
        let floor_point = fit(&floor_ok, &|i| {
            log2(cell_rows[i]["floor_amortized_ns"].as_f64().unwrap())
        });
        let floor_online_point = fit(&floor_ok, &|i| {
            log2(cell_rows[i]["floor_online_ns"].as_f64().unwrap())
        });
        let floor_adds_point = fit(&floor_ok, &|i| {
            log2(cell_rows[i]["floor_amortized_adds"].as_f64().unwrap())
        });
        let floor_slope_rep = |b: usize, sel: &[usize]| -> Option<f64> {
            let ys: Option<Vec<f64>> = sel
                .iter()
                .map(|&i| floor_reps[i][b].map(|c| log2(c.0)))
                .collect();
            ys.map(|y| ols_slope(&sel.iter().map(|&i| ns_of[i]).collect::<Vec<_>>(), &y))
        };
        let mut fs: Vec<f64> = (0..boot)
            .filter_map(|b| floor_slope_rep(b, &floor_ok))
            .collect();
        let floor_ci = (pct(&mut fs, 0.025), pct(&mut fs, 0.975));
        for (ci, t) in cells.iter().enumerate() {
            let mut v: Vec<f64> = floor_reps[ci]
                .iter()
                .filter_map(|c| c.map(|c| log2(c.0)))
                .collect();
            if !v.is_empty() {
                let (lo, hi) = (pct(&mut v, 0.025), pct(&mut v, 0.975));
                cell_rows[ci]["floor_log2_ci"] = json!([lo, hi]);
                csv.push_str(&format!(
                    "{m},floor-mitm,{},{},measured,{:.4},{lo:.4},{hi:.4},{:.4},group_adds\n",
                    t.n,
                    t.d,
                    log2(cell_rows[ci]["floor_amortized_ns"].as_f64().unwrap()),
                    log2(cell_rows[ci]["floor_amortized_adds"].as_f64().unwrap())
                ));
            }
        }

        md.push_str(&format!("## m = {m}\n\n### Cells and the floor\n\n"));
        md.push_str("| n | d | \\|F\\| points | p (decomposable) | decomps/target | exhaustive = MITM | floor amortized log2 ns [95%] | floor online log2 ns |\n|--:|--:|--:|--:|--:|--:|--:|--:|\n");
        for r in &cell_rows {
            let ci = r["floor_log2_ci"].as_array();
            md.push_str(&format!(
                "| {} | {} | {} | {:.3} | {:.3} | {}/{} | {} | {} |\n",
                r["n"],
                r["d"],
                r["fb_points"],
                r["p_decomposable"].as_f64().unwrap(),
                r["decomps_per_target"].as_f64().unwrap(),
                r["exhaustive_equal"],
                r["exhaustive_checked"],
                match (r["floor_amortized_ns"].as_f64(), ci) {
                    (Some(c), Some(ci)) => format!(
                        "{:.2} [{:.2}, {:.2}]",
                        log2(c),
                        ci[0].as_f64().unwrap(),
                        ci[1].as_f64().unwrap()
                    ),
                    _ => "p = 0".into(),
                },
                r["floor_online_ns"]
                    .as_f64()
                    .map(|c| format!("{:.2}", log2(c)))
                    .unwrap_or_default()
            ));
        }
        md.push_str(&format!(
            "\nFloor slope (log2 ns per unit n, amortized): **{floor_point:.3}** [{:.3}, {:.3}]; online {floor_online_point:.3}; counted (group additions, amortized) {floor_adds_point:.3}.\n\n",
            floor_ci.0, floor_ci.1
        ));

        // Arms.
        let mut arm_out = serde_json::Map::new();
        md.push_str("### Arms\n\n| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |\n|---|--:|--:|---|--:|--:|--:|--:|--:|\n");
        for arm in &arm_names {
            let acs = &arms[&(m, arm.clone())];
            let mut cells_json = Vec::new();
            let mut measured: Vec<usize> = Vec::new();
            let mut bounded: Vec<usize> = Vec::new();
            let mut point: BTreeMap<usize, f64> = BTreeMap::new();
            for ac in acs {
                let Some(ci) = cells.iter().position(|t| t.n == ac.n) else {
                    continue;
                };
                let p = cell_rows[ci]["p_decomposable"].as_f64().unwrap();
                let status = ac.status();
                let c = arm_cost(p, &ac.rows, |r| r.wall);
                let xors = arm_cost(p, &ac.rows, |r| r.xors);
                let confl = arm_cost(p, &ac.rows, |r| r.conflicts);
                let (counted, unit) = if arm == "sat" {
                    (confl, "conflicts")
                } else if arm.starts_with("wide") {
                    (arm_cost(p, &ac.rows, |r| r.nodes), "search_nodes")
                } else {
                    (xors, "word_xors")
                };
                let reps = arm_reps.get(&(arm.clone(), ac.n));
                let ci95 = reps.map(|r| {
                    let mut v: Vec<f64> = r.iter().filter_map(|c| c.map(log2)).collect();
                    (pct(&mut v, 0.025), pct(&mut v, 0.975))
                });
                if let Some(c) = c {
                    if !ac.rows.is_empty() {
                        point.insert(ci, log2(c));
                        match status {
                            "measured" => measured.push(ci),
                            "lower_bound" => bounded.push(ci),
                            _ => {}
                        }
                    }
                }
                let floor = cell_rows[ci]["floor_amortized_ns"].as_f64();
                let yes = ac.rows.iter().filter(|r| r.yes).count();
                let cens = ac.rows.iter().filter(|r| r.censored).count();
                md.push_str(&format!(
                    "| {arm} | {} | {} | {status} | {}/{} | {cens} | {} | {} | {} |\n",
                    ac.n,
                    ac.d,
                    yes,
                    ac.rows.len() - yes,
                    match (c, ci95, ac.rows.is_empty() || status == "unsupported") {
                        (Some(c), Some(ci), false) if status == "measured" =>
                            format!("{:.2} [{:.2}, {:.2}]", log2(c), ci.0, ci.1),
                        (Some(c), _, false) => format!("≥ {:.2}", log2(c)),
                        _ => "—".into(),
                    },
                    match (c, floor, ac.rows.is_empty() || status == "unsupported") {
                        (Some(c), Some(f), false) => format!(
                            "{}{:.3e}",
                            if status == "measured" { "" } else { "≥ " },
                            c / f
                        ),
                        _ => "—".into(),
                    },
                    match (counted, status) {
                        (Some(x), "measured") if x > 0.0 => format!("{x:.3e} {unit}"),
                        _ => "—".into(),
                    }
                ));
                if !ac.rows.is_empty() && status != "unsupported" {
                    if let Some(c) = c {
                        let (lo, hi) = ci95.unwrap_or((f64::NAN, f64::NAN));
                        csv.push_str(&format!(
                            "{m},{arm},{},{},{status},{:.4},{lo:.4},{hi:.4},{},{unit}\n",
                            ac.n,
                            ac.d,
                            log2(c),
                            counted
                                .filter(|x| *x > 0.0)
                                .map(|x| format!("{:.4}", log2(x)))
                                .unwrap_or_default()
                        ));
                    }
                }
                cells_json.push(json!({"n": ac.n, "d": ac.d, "status": status,
                    "targets": ac.rows.len(), "yes_targets": yes, "censored": cens,
                    "stopped": ac.stopped, "disagreements": ac.disagreements, "bad_witness": ac.bad_witness,
                    "cost_per_relation_ns": c, "log2_ci": ci95.map(|c| [c.0, c.1]),
                    "ratio_to_floor": c.zip(floor).map(|(c, f)| c / f),
                    "counted_per_relation": counted, "counted_unit": unit}));
            }
            // Slopes and the decision rule.
            let slope_of = |sel: &[usize]| -> Option<f64> {
                (sel.len() >= 2).then(|| {
                    ols_slope(
                        &sel.iter().map(|&i| ns_of[i]).collect::<Vec<_>>(),
                        &sel.iter().map(|&i| point[&i]).collect::<Vec<_>>(),
                    )
                })
            };
            let boot_delta = |sel: &[usize]| -> Option<(f64, f64, f64, f64)> {
                if sel.len() < 2 {
                    return None;
                }
                let mut ds = Vec::new();
                let mut ss = Vec::new();
                for b in 0..boot {
                    let ys: Option<Vec<f64>> = sel
                        .iter()
                        .map(|&i| {
                            arm_reps
                                .get(&(arm.clone(), cells[i].n))
                                .and_then(|r| r[b])
                                .map(log2)
                        })
                        .collect();
                    let Some(ys) = ys else { continue };
                    let s = ols_slope(&sel.iter().map(|&i| ns_of[i]).collect::<Vec<_>>(), &ys);
                    // The floor fitted on the same cells, same replicate.
                    if let Some(f) = floor_slope_rep(b, sel) {
                        ds.push(s - f);
                        ss.push(s);
                    }
                }
                if ds.is_empty() {
                    return None;
                }
                Some((
                    pct(&mut ss, 0.025),
                    pct(&mut ss, 0.975),
                    pct(&mut ds, 0.025),
                    pct(&mut ds, 0.975),
                ))
            };
            let floor_on = |sel: &[usize]| -> Option<f64> {
                (sel.len() >= 2).then(|| {
                    fit(sel, &|i| {
                        log2(cell_rows[i]["floor_amortized_ns"].as_f64().unwrap())
                    })
                })
            };
            let mut all_sel: Vec<usize> = measured.iter().chain(bounded.iter()).copied().collect();
            all_sel.sort_unstable();
            let mean_measured_n = mean(&measured.iter().map(|&i| ns_of[i]).collect::<Vec<_>>());
            // Substituting lower bounds gives a lower bound on the slope
            // only when every bounded cell lies above the fit's mean n and
            // the measured cells sit below it (OLS weights).
            let lb_valid = !bounded.is_empty() && {
                let mean_all = mean(&all_sel.iter().map(|&i| ns_of[i]).collect::<Vec<_>>());
                bounded.iter().all(|&i| ns_of[i] > mean_all)
            };
            let m_slope = slope_of(&measured);
            let m_boot = boot_delta(&measured);
            let lb_slope = slope_of(&all_sel);
            let lb_boot = boot_delta(&all_sel);
            let verdict = if measured.len() >= 3 {
                match m_boot {
                    Some((_, _, _, hi)) if hi < 0.0 => "grows more slowly than the floor",
                    Some((_, _, lo, _)) if lo > 0.0 => "grows faster than the floor",
                    Some(_) => "indistinguishable from the floor on this range",
                    None => "undetermined",
                }
            } else if all_sel.len() >= 3 && lb_valid {
                match lb_boot {
                    Some((_, _, lo, _)) if lo > 0.0 => {
                        "grows faster than the floor (lower-bound fit)"
                    }
                    _ => "undetermined (lower-bound fit does not separate)",
                }
            } else {
                "undetermined (fewer than 3 usable cells)"
            };
            md.push_str(&format!(
                "\n**{arm}**: measured cells {}, lower-bound cells {}; slope on measured cells {} (floor on the same cells {}), Δ = arm − floor 95% {}; lower-bound fit slope {} ({}), Δ 95% {}. **Verdict: {verdict}.**\n\n",
                measured.len(), bounded.len(),
                m_slope.map(|s| format!("{s:.3}")).unwrap_or("—".into()),
                floor_on(&measured).map(|s| format!("{s:.3}")).unwrap_or("—".into()),
                m_boot.map(|b| format!("[{:.3}, {:.3}]", b.2, b.3)).unwrap_or("—".into()),
                lb_slope.map(|s| format!("{s:.3}")).unwrap_or("—".into()),
                if lb_valid { "valid as a lower bound" } else { "not a valid bound" },
                lb_boot.map(|b| format!("[{:.3}, {:.3}]", b.2, b.3)).unwrap_or("—".into()),
            ));
            if arm_names.last() != Some(arm) {
                md.push_str("| arm | n | d | status | targets run (yes/no) | censored | log2 ns per relation [95%] | ÷ floor | counted per relation |\n|---|--:|--:|---|--:|--:|--:|--:|--:|\n");
            }
            let _ = mean_measured_n;
            arm_out.insert(
                arm.clone(),
                json!({"cells": cells_json,
                "measured_cells": measured.iter().map(|&i| cells[i].n).collect::<Vec<_>>(),
                "lower_bound_cells": bounded.iter().map(|&i| cells[i].n).collect::<Vec<_>>(),
                "slope_measured": m_slope, "slope_measured_ci": m_boot.map(|b| [b.0, b.1]),
                "floor_slope_same_cells": floor_on(&measured),
                "delta_ci": m_boot.map(|b| [b.2, b.3]),
                "slope_lower_bound_fit": lb_slope, "lower_bound_fit_valid": lb_valid,
                "delta_lower_bound_ci": lb_boot.map(|b| [b.2, b.3]),
                "verdict": verdict}),
            );
        }
        report["by_m"][m.to_string()] = json!({"cells": cell_rows,
            "floor_slope_amortized": floor_point, "floor_slope_amortized_ci": [floor_ci.0, floor_ci.1],
            "floor_slope_online": floor_online_point, "floor_slope_counted_adds": floor_adds_point,
            "arms": arm_out});
    }
    std::fs::write(&out, serde_json::to_string_pretty(&report).unwrap()).unwrap();
    std::fs::write(out.with_extension("md"), md).unwrap();
    std::fs::write(out.with_file_name("cells.csv"), csv).unwrap();
    println!("wrote {}", out.display());
}
