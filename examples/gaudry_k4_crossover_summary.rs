//! Summarise `research/gaudry_k4_char2_crossover_20261015` (PROTOCOL.md
//! Parts D, E and F) from its ecbench sessions' records:
//!
//! - every run: `S`, total GAE, the relation and linear-algebra phases,
//!   trials, relations, the matrix's `row_ops`;
//! - Part D: per workload, `ic-dense` against `ic-wied` and against the
//!   committed #1883 value of the same workload;
//! - Part E: per target, IC `S` ÷ rho's theoretical `S` (0.886) and ÷ the
//!   measured `rho-neg` run on the same workload (matched by workload id,
//!   across sessions where rho and IC ran apart); per size the mean; the
//!   crossover `m*` and `log₂ n*` by linear interpolation of
//!   `log₂(IC ÷ rho theory)` in `log₂ n` between adjacent sizes; the
//!   matrix's share of IC's cost;
//! - Part F: per workload, each large-prime arm's `S` ÷ `ic-wied`'s, and
//!   its large-prime counters.
//!
//! Writes `raw/summary.json`.
//!
//! ```sh
//! cargo run --release --example gaudry_k4_crossover_summary -- research/gaudry_k4_char2_crossover_20261015
//! ```

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::Path;

use serde_json::{json, Value};

const RHO_THEORY_S: f64 = 0.886;

fn num(v: &Value) -> f64 {
    match v {
        Value::Number(n) => n.as_f64().unwrap(),
        Value::String(s) => s.parse::<f64>().unwrap(),
        _ => f64::NAN,
    }
}

fn records(dir: &Path) -> Vec<(String, Value)> {
    let mut out = Vec::new();
    let Ok(entries) = fs::read_dir(dir) else {
        return out;
    };
    let mut names: Vec<_> = entries
        .filter_map(|e| e.ok())
        .filter(|e| e.path().join("records.jsonl").exists())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    names.sort();
    for name in names {
        let text = fs::read_to_string(dir.join(&name).join("records.jsonl")).unwrap();
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            out.push((name.clone(), serde_json::from_str(line).unwrap()));
        }
    }
    out
}

fn phase_gae(r: &Value, name: &str) -> f64 {
    r["phases"]
        .as_array()
        .and_then(|ps| ps.iter().find(|p| p["name"] == name))
        .map_or(0.0, |p| num(&p["gae"]))
}

fn row(session: &str, r: &Value) -> Value {
    let curve = &r["workload"]["curve"];
    let n = num(&curve["group_order"]);
    let rr = num(&curve["r"]);
    let c = &r["counters"];
    json!({
        "session": session,
        "arm": r["arm"],
        "workload": r["workload"]["workload_id"],
        "slug": curve["slug"],
        "m": r["workload"]["curve_spec"]["m"],
        "h": num(&curve["cofactor"]),
        "log2_n": n.log2(),
        "log2_r": rr.log2(),
        "status": r["outcome"]["status"],
        "s": num(&r["cost"]["s"]),
        "total_gae": num(&r["cost"]["total_gae"]),
        "relations_gae": phase_gae(r, "relations"),
        "linear_algebra_gae": phase_gae(r, "linear_algebra"),
        "trials": c.get("targets_tried").cloned().unwrap_or(Value::Null),
        "relations": c.get("relations_found").cloned().unwrap_or(Value::Null),
        "matrix_work": c.get("matrix_work").cloned().unwrap_or(Value::Null),
        "solve_wall_s": num(&r["time"]["solve_wall_ns"]) / 1e9,
        "lp": r["phases"].as_array().and_then(|ps| ps.iter().find(|p| p["name"] == "relations"))
            .map(|p| p["native"].as_object().map(|o| o.iter().filter(|(k, _)| k.starts_with("lp_")).map(|(k, v)| (k.clone(), v.clone())).collect::<serde_json::Map<_, _>>()))
            .unwrap_or_default(),
    })
}

fn main() {
    let dir = env::args().nth(1).expect("the study directory");
    let dir = Path::new(&dir);
    let rows: Vec<Value> = records(&dir.join("sessions"))
        .iter()
        .map(|(s, r)| row(s, r))
        .collect();
    let by = |sess_prefix: &str, arm: &str| -> BTreeMap<String, &Value> {
        rows.iter()
            .filter(|r| r["session"].as_str().unwrap().starts_with(sess_prefix) && r["arm"] == arm)
            .map(|r| (r["workload"].as_str().unwrap().to_string(), r))
            .collect()
    };

    // Part D, with #1883's committed ic-frobenius rows.
    let committed: BTreeMap<String, f64> =
        records(&dir.join("../gaudry_quartic_d4_char2_20261015/sessions"))
            .iter()
            .filter(|(_, r)| r["arm"] == "ic-frobenius")
            .map(|(_, r)| {
                (
                    r["workload"]["workload_id"].as_str().unwrap().to_string(),
                    num(&r["cost"]["s"]),
                )
            })
            .collect();
    let dense = by("d-", "ic-dense");
    let wied = by("d-", "ic-wied");
    let part_d: Vec<Value> = dense
        .iter()
        .filter_map(|(w, d)| {
            let x = wied.get(w)?;
            let f = |v: &Value, k: &str| num(&v[k]);
            Some(json!({
                "workload": w,
                "m": d["m"],
                "dense_s": d["s"],
                "committed_s": committed.get(w),
                "reproduces_committed": committed.get(w).map(|c| *c == f(d, "s")),
                "wied_s": x["s"],
                "s_ratio_wied_over_dense": f(x, "s") / f(d, "s"),
                "trials_ratio_wied_over_dense": f(x, "trials") / f(d, "trials"),
                "row_ops_dense": d["matrix_work"],
                "row_ops_wied": x["matrix_work"],
                "row_ops_ratio_dense_over_wied": f(d, "matrix_work") / f(x, "matrix_work"),
                "la_share_dense": f(d, "linear_algebra_gae") / f(d, "total_gae"),
                "la_share_wied": f(x, "linear_algebra_gae") / f(x, "total_gae"),
            }))
        })
        .collect();

    // Part E.
    let ic = by("e-", "ic-wied");
    let rho = by("e-", "rho-neg");
    let mut part_e_targets = Vec::new();
    let mut per_m: BTreeMap<i64, Vec<(f64, f64, f64)>> = BTreeMap::new();
    for (w, i) in &ic {
        let s_ic = num(&i["s"]);
        let ratio_theory = s_ic / RHO_THEORY_S;
        let rho_s = rho.get(w).map(|r| num(&r["s"]));
        let m = i["m"].as_i64().unwrap();
        per_m.entry(m).or_default().push((
            num(&i["log2_n"]),
            ratio_theory,
            rho_s.map_or(f64::NAN, |r| s_ic / r),
        ));
        part_e_targets.push(json!({
            "workload": w,
            "m": m,
            "h": i["h"],
            "log2_n": i["log2_n"],
            "log2_r": i["log2_r"],
            "status": i["status"],
            "rho_status": rho.get(w).map(|r| r["status"].clone()),
            "ic_s": s_ic,
            "rho_s": rho_s,
            "ic_over_rho_theory": ratio_theory,
            "ic_over_rho_measured": rho_s.map(|r| s_ic / r),
            "la_share": num(&i["linear_algebra_gae"]) / num(&i["total_gae"]),
            "trials": i["trials"],
            "relations": i["relations"],
            "trials_per_relation": num(&i["trials"]) / num(&i["relations"]),
            "ic_wall_s": i["solve_wall_s"],
            "rho_wall_s": rho.get(w).map(|r| r["solve_wall_s"].clone()),
        }));
    }
    let means: Vec<(i64, f64, f64, f64)> = per_m
        .iter()
        .map(|(m, v)| {
            let n = v.len() as f64;
            // geometric means over the targets
            let ln = v.iter().map(|x| x.0).sum::<f64>() / n;
            let rt = (v.iter().map(|x| x.1.ln()).sum::<f64>() / n).exp();
            let rm = (v.iter().map(|x| x.2.ln()).sum::<f64>() / n).exp();
            (*m, ln, rt, rm)
        })
        .collect();
    let crossover = means.windows(2).find_map(|w| {
        let (a, b) = (w[0], w[1]);
        (a.2 >= 1.0 && b.2 < 1.0).then(|| {
            let (ya, yb) = (a.2.log2(), b.2.log2());
            let t = ya / (ya - yb);
            json!({
                "between_m": [a.0, b.0],
                "m_star": a.0 as f64 + t * (b.0 - a.0) as f64,
                "log2_n_star": a.1 + t * (b.1 - a.1),
            })
        })
    });
    let slope = {
        let pts: Vec<(f64, f64)> = means.iter().map(|x| (x.1, x.2.log2())).collect();
        let n = pts.len() as f64;
        if n >= 2.0 {
            let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
            let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
            let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
            let sxx: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
            Some(sxy / sxx)
        } else {
            None
        }
    };

    // Part F.
    let base = by("f-", "ic-wied");
    let mut part_f = Vec::new();
    for arm in ["lp1-700", "lp2-700", "lp2-500"] {
        for (w, x) in by("f-", arm) {
            let Some(b) = base.get(&w) else { continue };
            part_f.push(json!({
                "arm": arm,
                "workload": w,
                "m": x["m"],
                "s": x["s"],
                "base_s": b["s"],
                "s_ratio_over_base": num(&x["s"]) / num(&b["s"]),
                "trials_ratio_over_base": num(&x["trials"]) / num(&b["trials"]),
                "la_gae_ratio_over_base": num(&x["linear_algebra_gae"]) / num(&b["linear_algebra_gae"]),
                "lp": x["lp"],
            }));
        }
    }

    let out = json!({
        "schema": "gaudry_k4_crossover.summary/v1",
        "rho_theory_s": RHO_THEORY_S,
        "runs": rows,
        "part_d": part_d,
        "part_e": {
            "targets": part_e_targets,
            "per_m": means.iter().map(|x| json!({"m": x.0, "log2_n": x.1, "ic_over_rho_theory": x.2, "ic_over_rho_measured": x.3})).collect::<Vec<_>>(),
            "slope_log2_ratio_per_log2_n": slope,
            "crossover": crossover,
        },
        "part_f": part_f,
    });
    let path = dir.join("raw/summary.json");
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!("wrote {}", path.display());
}
