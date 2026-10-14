//! `analyze`: the one table, the ratios to each boundary and the fitted
//! growth of every arm, from the raw records of the study directory.
//! Reads what the other subcommands wrote (and the committed F4 rows the
//! `prime-rows` pass copied), writes `ANALYSIS.json` and `ANALYSIS.md`.
//! Medians, least-squares slopes in bits per rung with their standard
//! error, and ratios; nothing is re-measured.

use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::path::Path;

fn read_jsonl(path: &Path) -> Vec<Value> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    text.lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

fn f(v: &Value, k: &str) -> Option<f64> {
    v.get(k).and_then(|x| x.as_f64())
}

fn s<'a>(v: &'a Value, k: &str) -> &'a str {
    v.get(k).and_then(|x| x.as_str()).unwrap_or("")
}

fn median(xs: &mut [f64]) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = xs.len();
    Some(if n % 2 == 1 {
        xs[n / 2]
    } else {
        0.5 * (xs[n / 2 - 1] + xs[n / 2])
    })
}

/// Least squares of `y = a + b x`: `(b, standard error of b, r², n)`.
fn slope(points: &[(f64, f64)]) -> Option<(f64, f64, f64, usize)> {
    let n = points.len();
    if n < 2 {
        return None;
    }
    let nf = n as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / nf;
    let my = points.iter().map(|p| p.1).sum::<f64>() / nf;
    let sxx: f64 = points.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let sxy: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let syy: f64 = points.iter().map(|p| (p.1 - my).powi(2)).sum();
    if sxx == 0.0 {
        return None;
    }
    let b = sxy / sxx;
    let a = my - b * mx;
    let ss_res: f64 = points.iter().map(|p| (p.1 - (a + b * p.0)).powi(2)).sum();
    let r2 = if syy > 0.0 { 1.0 - ss_res / syy } else { 1.0 };
    let se = if n > 2 {
        (ss_res / (nf - 2.0) / sxx).sqrt()
    } else {
        f64::NAN
    };
    Some((b, se, r2, n))
}

#[derive(Default, Clone)]
struct Cell {
    ns: Vec<f64>,
    ns_contended: Vec<f64>,
    mul: Vec<f64>,
    sqr: Vec<f64>,
    inv: Vec<f64>,
    lookups: Vec<f64>,
    candidates: Vec<f64>,
    table_bytes: f64,
    setup_secs: Vec<f64>,
    base_points: f64,
    n: f64,
}

pub fn run(args: &[String]) {
    let dir = args
        .iter()
        .position(|a| a == "--dir")
        .and_then(|i| args.get(i + 1))
        .cloned()
        .unwrap_or_else(|| "research/pdp_root_finder_20261012".into());
    let raw = Path::new(&dir).join("raw");

    // ── isolation records: command → contended flag ───────────────────
    let mut contended: BTreeMap<String, bool> = BTreeMap::new();
    for rec in read_jsonl(&raw.join("isolated.jsonl")) {
        let c = rec
            .get("contended")
            .and_then(|x| x.as_bool())
            .unwrap_or(false);
        if let Some(lab) = rec.get("label").and_then(|x| x.as_str()) {
            if !lab.is_empty() {
                contended.insert(lab.to_string(), c);
                continue;
            }
        }
        let label = rec
            .get("argv")
            .or_else(|| rec.get("command"))
            .map(|c| c.to_string())
            .unwrap_or_default();
        // Key on the `--label` value inside the recorded command.
        if let Some(i) = label.find("--label") {
            let rest = &label[i + 7..];
            let tok: String = rest
                .trim_start_matches(['"', ',', ' '])
                .chars()
                .take_while(|c| c.is_alphanumeric() || *c == '-' || *c == '_')
                .collect();
            contended.insert(tok, c);
        }
    }

    // ── binary timing rows ────────────────────────────────────────────
    let mut bin: BTreeMap<(String, u64, u64), Cell> = BTreeMap::new();
    for row in read_jsonl(&raw.join("bin_time.jsonl")) {
        let arm = s(&row, "arm").to_string();
        let l = f(&row, "l").unwrap_or(0.0) as u64;
        let n = f(&row, "n").unwrap_or(0.0) as u64;
        let c = bin.entry((arm, n, l)).or_default();
        c.n = n as f64;
        let label = s(&row, "label").to_string();
        let is_contended = contended.get(&label).copied().unwrap_or(false);
        if let Some(ns) = f(&row, "ns_per_target") {
            if is_contended {
                c.ns_contended.push(ns);
            } else {
                c.ns.push(ns);
            }
        }
        for (k, v) in [
            ("mul_per_target", &mut c.mul),
            ("sqr_per_target", &mut c.sqr),
            ("inv_per_target", &mut c.inv),
            ("lookups_per_target", &mut c.lookups),
            ("candidates_per_target", &mut c.candidates),
            ("setup_secs", &mut c.setup_secs),
        ] {
            if let Some(x) = f(&row, k) {
                v.push(x);
            }
        }
        c.table_bytes = f(&row, "table_bytes").unwrap_or(0.0);
        c.base_points = f(&row, "base_points").unwrap_or(0.0);
    }
    // Counts rows (the PS arms carry no counters in the timing rows).
    for row in read_jsonl(&raw.join("bin_count.jsonl")) {
        let arm = s(&row, "arm").to_string();
        let l = f(&row, "l").unwrap_or(0.0) as u64;
        let n = f(&row, "n").unwrap_or(0.0) as u64;
        let c = bin.entry((arm, n, l)).or_default();
        if let Some(pt) = row.get("per_target") {
            c.mul = vec![f(pt, "mul").unwrap_or(0.0)];
            c.sqr = vec![f(pt, "sqr").unwrap_or(0.0)];
            c.inv = vec![f(pt, "inv").unwrap_or(0.0)];
            c.lookups = vec![f(pt, "lookups").unwrap_or(0.0)];
            c.candidates = vec![f(pt, "candidates").unwrap_or(0.0)];
        }
        c.base_points = f(&row, "base_points").unwrap_or(c.base_points);
        c.n = n as f64;
        if let Some(pr) = row.get("prepare") {
            c.table_bytes = f(pr, "table_bytes").unwrap_or(c.table_bytes);
        }
    }

    // ── prime ladder rows ─────────────────────────────────────────────
    let mut prime: BTreeMap<(u64, String, u64), Cell> = BTreeMap::new();
    for name in [
        "prime_ladder_m2.jsonl",
        "prime_ladder_m3.jsonl",
        "prime_ladder_m4.jsonl",
        "prime_time.jsonl",
    ] {
        for row in read_jsonl(&raw.join(name)) {
            let m = f(&row, "m").unwrap_or(0.0) as u64;
            let t = f(&row, "t").unwrap_or(0.0) as u64;
            let arm = s(&row, "arm").to_string();
            let label = s(&row, "label").to_string();
            let c = prime.entry((m, arm, t)).or_default();
            // The counts pass is contended by construction (other jobs ran); only
            // rows whose label is a timing round carry wall evidence.
            let timing_row = label.starts_with("round");
            let is_contended = contended.get(&label).copied().unwrap_or(!timing_row);
            if let Some(ns) = f(&row, "ns_per_target") {
                if is_contended {
                    c.ns_contended.push(ns);
                } else {
                    c.ns.push(ns);
                }
            }
            for (k, v) in [
                ("mul_per_target", &mut c.mul),
                ("sqr_per_target", &mut c.sqr),
                ("inv_per_target", &mut c.inv),
                ("lookups_per_target", &mut c.lookups),
                ("candidates_per_target", &mut c.candidates),
            ] {
                if let Some(x) = f(&row, k) {
                    v.push(x);
                }
            }
            if let Some(tb) = row.get("table") {
                c.table_bytes = f(tb, "bytes").unwrap_or(0.0);
                if let Some(sec) = f(tb, "secs") {
                    c.setup_secs.push(sec);
                }
            }
            c.base_points = f(&row, "base_points").unwrap_or(0.0);
            c.n = f(&row, "N").unwrap_or(0.0);
        }
    }

    // ── exactness census of the committed rows ────────────────────────
    let mut census = json!({"rows": 0, "brute_checked": 0, "brute_agree": 0, "f4_finished_rows": 0, "f4_verdict_agree": 0, "f4_verdict_disagree": 0, "f4_unfinished_rows": 0, "algebraic_checked": 0, "algebraic_but_no_rational": 0});
    for row in read_jsonl(&raw.join("prime_rows.jsonl")) {
        census["rows"] = json!(census["rows"].as_u64().unwrap() + 1);
        if let Some(arms) = row.get("arms").and_then(|a| a.as_array()) {
            for a in arms {
                if let Some(b) = a.get("brute_agree").and_then(|x| x.as_bool()) {
                    census["brute_checked"] = json!(census["brute_checked"].as_u64().unwrap() + 1);
                    if b {
                        census["brute_agree"] = json!(census["brute_agree"].as_u64().unwrap() + 1);
                    }
                }
            }
        }
        let alg = row
            .get("algebraic_ordered_solutions")
            .and_then(|x| x.as_u64());
        let br = row.get("brute_rational").and_then(|x| x.as_u64());
        if let (Some(alg), Some(br)) = (alg, br) {
            census["algebraic_checked"] = json!(census["algebraic_checked"].as_u64().unwrap() + 1);
            if alg > 0 && br == 0 {
                census["algebraic_but_no_rational"] =
                    json!(census["algebraic_but_no_rational"].as_u64().unwrap() + 1);
            }
        }
        if let Some(f4r) = row.get("f4") {
            let finished = f4r
                .get("finished")
                .and_then(|x| x.as_bool())
                .unwrap_or(false);
            if !finished {
                census["f4_unfinished_rows"] =
                    json!(census["f4_unfinished_rows"].as_u64().unwrap() + 1);
            } else if let (Some(inc), Some(alg)) =
                (f4r.get("inconsistent").and_then(|x| x.as_bool()), alg)
            {
                census["f4_finished_rows"] =
                    json!(census["f4_finished_rows"].as_u64().unwrap() + 1);
                if inc == (alg == 0) {
                    census["f4_verdict_agree"] =
                        json!(census["f4_verdict_agree"].as_u64().unwrap() + 1);
                } else {
                    census["f4_verdict_disagree"] =
                        json!(census["f4_verdict_disagree"].as_u64().unwrap() + 1);
                }
            }
        }
    }
    // ── committed F4 rows (prime) ─────────────────────────────────────
    let mut f4: BTreeMap<(String, u64, u64, u64), Vec<(f64, f64, f64)>> = BTreeMap::new(); // (kind,p,m,t) → (muladds, ms, degree)
    let mut f4_unfinished: BTreeMap<(String, u64, u64, u64), u64> = BTreeMap::new();
    let mut rf_at_rows: BTreeMap<(String, u64, u64, u64, String), Vec<(f64, f64)>> =
        BTreeMap::new(); // → (mul, ns)
    for row in read_jsonl(&raw.join("prime_rows.jsonl")) {
        let key = (
            s(&row, "kind").to_string(),
            f(&row, "p").unwrap_or(0.0) as u64,
            f(&row, "m").unwrap_or(0.0) as u64,
            f(&row, "t").unwrap_or(0.0) as u64,
        );
        if let Some(f4r) = row.get("f4") {
            let finished = f4r
                .get("finished")
                .and_then(|x| x.as_bool())
                .unwrap_or(false);
            if !finished {
                *f4_unfinished.entry(key.clone()).or_default() += 1;
            } else if let (Some(mm), Some(ms)) = (f(f4r, "muladds"), f(f4r, "ms")) {
                f4.entry(key.clone()).or_default().push((
                    mm,
                    ms,
                    f(f4r, "degree").unwrap_or(f64::NAN),
                ));
            }
        }
        if let Some(arms) = row.get("arms").and_then(|a| a.as_array()) {
            for a in arms {
                let k2 = (key.0.clone(), key.1, key.2, key.3, s(a, "arm").to_string());
                rf_at_rows
                    .entry(k2)
                    .or_default()
                    .push((f(a, "mul").unwrap_or(0.0), f(a, "ns_median").unwrap_or(0.0)));
            }
        }
    }

    // ── binary Gröbner references on (23, 8) ──────────────────────────
    let mut f6_xors = Vec::new();
    let mut f6_ms = Vec::new();
    for row in read_jsonl(&raw.join("f6_ic_probe_n23_dim8.jsonl")) {
        if row.get("target").is_some() {
            if let (Some(x), Some(ms)) = (f(&row, "word_xors"), f(&row, "ms")) {
                f6_xors.push(x);
                f6_ms.push(ms);
            }
        }
    }
    let mut split_rows = Vec::new();
    for row in read_jsonl(&raw.join("f4_split_3-1-23-8.jsonl")) {
        if row
            .get("kind")
            .and_then(|k| k.as_str())
            .map(|k| k == "cell" || k == "config")
            .unwrap_or(false)
        {
            continue;
        }
        split_rows.push(row);
    }

    // ── assemble ──────────────────────────────────────────────────────
    let mut out = Map::new();
    let mut md = String::new();
    md.push_str("# Analysis (written by `pdp_rootfind analyze`; medians, least-squares slopes, ratios of recorded values)\n\n");

    // Binary table.
    let mut bin_rows = Vec::new();
    md.push_str("## Binary, per target (uniform refuted targets)\n\n");
    md.push_str("| arm | n | ℓ | base points | mul | sqr | inv | lookups | candidates | ns median (uncontended) | ns min | runs | contended runs | table MB | ÷ rf-pairs mul | ÷ rf-mitm mul | ÷ baseline ns |\n|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let get_mul = |arm: &str, n: u64, l: u64| -> Option<f64> {
        bin.get(&(arm.to_string(), n, l))
            .and_then(|c| median(&mut c.mul.clone()))
    };
    let get_ns = |arm: &str, n: u64, l: u64| -> Option<f64> {
        bin.get(&(arm.to_string(), n, l))
            .and_then(|c| median(&mut c.ns.clone()))
    };
    for ((arm, n, l), c) in &bin {
        let ns_med = median(&mut c.ns.clone());
        let ns_min = c.ns.iter().cloned().fold(f64::INFINITY, f64::min);
        let mul = median(&mut c.mul.clone());
        let r_pairs = mul.and_then(|m| get_mul("rf-pairs", *n, *l).map(|p| p / m));
        let r_mitm = mul.and_then(|m| get_mul("rf-mitm", *n, *l).map(|p| m / p));
        let r_base = ns_med.and_then(|t| get_ns("baseline", *n, *l).map(|b| b / t));
        let fmt = |x: Option<f64>| x.map(|v| format!("{v:.1}")).unwrap_or("—".into());
        let fmt3 = |x: Option<f64>| x.map(|v| format!("{v:.3}")).unwrap_or("—".into());
        md.push_str(&format!(
            "| {arm} | {n} | {l} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.1} | {} | {} | {} |\n",
            c.base_points, fmt(mul), fmt(median(&mut c.sqr.clone())), fmt3(median(&mut c.inv.clone())), fmt(median(&mut c.lookups.clone())), fmt(median(&mut c.candidates.clone())),
            fmt(ns_med), if ns_min.is_finite() { format!("{ns_min:.1}") } else { "—".into() }, c.ns.len(), c.ns_contended.len(), c.table_bytes / 1e6,
            r_pairs.map(|v| format!("{v:.2}")).unwrap_or("—".into()), r_mitm.map(|v| format!("{v:.2}")).unwrap_or("—".into()), r_base.map(|v| format!("{v:.1}")).unwrap_or("—".into())
        ));
        let t_enum = 2f64.powi(2 * *l as i32) / 2.0;
        let pair_floor = 2f64.powi(*l as i32);
        let cand = median(&mut c.candidates.clone());
        let look = median(&mut c.lookups.clone());
        let r_enum = cand.map(|v| v / t_enum);
        let r_pair = look.map(|v| v / pair_floor);
        bin_rows.push(json!({"arm": arm, "n": n, "l": l, "base_points": c.base_points, "T_enum": t_enum, "pair_floor": pair_floor, "ratio_candidates_over_T_enum": r_enum, "ratio_lookups_over_pair_floor": r_pair, "mul": mul, "sqr": median(&mut c.sqr.clone()), "inv": median(&mut c.inv.clone()), "lookups": median(&mut c.lookups.clone()), "candidates": median(&mut c.candidates.clone()), "ns_median": ns_med, "ns_min": if ns_min.is_finite() { Some(ns_min) } else { None }, "ns_runs": c.ns.len(), "ns_contended_runs": c.ns_contended.len(), "ns_all": c.ns, "table_bytes": c.table_bytes, "setup_secs": median(&mut c.setup_secs.clone()), "ratio_mul_vs_rf_pairs": r_pairs, "ratio_mul_vs_rf_mitm": r_mitm, "ratio_ns_vs_baseline": r_base}));
    }
    out.insert("binary".into(), Value::Array(bin_rows));
    // S·T² accounting of the table arm (G5): S = stored keys, T = per-target
    // work in group-operation equivalents (mul + sqr over 6, the cost of one
    // batched affine addition in field operations), r ≈ 2^(n−1) (cofactor 2),
    // |Aut| = 2 (negation only; the standard subspace is not Frobenius-stable).
    let mut st2 = Vec::new();
    md.push_str("\n## S·T² of the table arm against the generic line (model: T = (mul + sqr)/6 group operations, S = stored keys = |F|², r = 2^(n−1), |Aut| = 2)\n\n| arm | n | ℓ | S (keys) | T (group-op eq.) | log₂(S·T²) | log₂(r/|Aut|) | gap (bits) |\n|:--|--:|--:|--:|--:|--:|--:|--:|\n");
    for ((arm, n, l), c) in &bin {
        if arm != "rf-mitm" || c.table_bytes == 0.0 {
            continue;
        }
        let (Some(mul), Some(sqr)) = (median(&mut c.mul.clone()), median(&mut c.sqr.clone()))
        else {
            continue;
        };
        // Stored keys: one per doubling and two per unordered pair, m² in all.
        let s_keys = c.base_points * c.base_points;
        let t_gae = (mul + sqr) / 6.0;
        let lhs = (s_keys * t_gae * t_gae).log2();
        let rhs = (*n as f64 - 1.0) - 1.0;
        md.push_str(&format!(
            "| {arm} | {n} | {l} | {s_keys:.3e} | {t_gae:.1} | {lhs:.1} | {rhs:.1} | {:+.1} |\n",
            lhs - rhs
        ));
        st2.push(json!({"arm": arm, "n": n, "l": l, "S_keys": s_keys, "T_gae": t_gae, "log2_ST2": lhs, "log2_r_over_aut": rhs, "gap_bits": lhs - rhs}));
    }
    out.insert("st2_binary".into(), Value::Array(st2));

    // Binary slopes per arm (bits per rung of ℓ), over the n = 3ℓ rungs only.
    let mut slopes = Vec::new();
    md.push_str("\n## Growth, bits per rung (least squares over the n = 3ℓ rungs; standard error; r²)\n\n| track | arm | quantity | slope | s.e. | r² | rungs |\n|:--|:--|:--|--:|--:|--:|--:|\n");
    let arms: Vec<String> = bin
        .keys()
        .map(|k| k.0.clone())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    for arm in &arms {
        for (qname, pick) in [("mul", 0), ("ns", 1)] {
            let pts: Vec<(f64, f64)> = bin
                .iter()
                .filter(|((a, n, l), _)| a == arm && *n == 3 * l)
                .filter_map(|((_, _, l), c)| {
                    let v = if pick == 0 {
                        median(&mut c.mul.clone())
                    } else {
                        median(&mut c.ns.clone())
                    };
                    v.filter(|v| *v > 0.0).map(|v| (*l as f64, v.log2()))
                })
                .collect();
            if let Some((b, se, r2, n)) = slope(&pts) {
                md.push_str(&format!(
                    "| binary | {arm} | {qname} | {b:.3} | {se:.3} | {r2:.4} | {n} |\n"
                ));
                slopes.push(json!({"track": "binary", "arm": arm, "quantity": qname, "slope_bits_per_rung": b, "se": se, "ci95_halfwidth": 1.96 * se, "r2": r2, "rungs": n, "points": pts}));
            }
        }
    }

    // Prime table.
    let mut prime_rows = Vec::new();
    md.push_str("\n## Prime, per target (uniform refuted targets, Kummer towers at p₁)\n\n| m | arm | t | N | base points | mul | sqr | inv | lookups | candidates | ns median (uncontended) | ns min | runs | contended runs | table MB | ÷ pairs3 mul | ÷ mitm3 mul |\n|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let pget = |m: u64, arm: &str, t: u64| -> Option<f64> {
        prime
            .get(&(m, arm.to_string(), t))
            .and_then(|c| median(&mut c.mul.clone()))
    };
    for ((m, arm, t), c) in &prime {
        let ns_med = median(&mut c.ns.clone());
        let ns_c = median(&mut c.ns_contended.clone());
        let ns_min =
            c.ns.iter()
                .chain(c.ns_contended.iter())
                .cloned()
                .fold(f64::INFINITY, f64::min);
        let mul = median(&mut c.mul.clone());
        let r_pairs = mul.and_then(|x| pget(*m, "pairs3", *t).map(|p| p / x));
        let r_mitm = mul.and_then(|x| pget(*m, "mitm3", *t).map(|p| x / p));
        let fmt = |x: Option<f64>| x.map(|v| format!("{v:.1}")).unwrap_or("—".into());
        md.push_str(&format!(
            "| {m} | {arm} | {t} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.1} | {} | {} |\n",
            c.n, c.base_points, fmt(mul), fmt(median(&mut c.sqr.clone())), fmt(median(&mut c.inv.clone())), fmt(median(&mut c.lookups.clone())), fmt(median(&mut c.candidates.clone())),
            fmt(ns_med), if ns_min.is_finite() { format!("{ns_min:.1}") } else { "—".into() }, c.ns.len(), c.ns_contended.len(), c.table_bytes / 1e6,
            r_pairs.map(|v| format!("{v:.2}")).unwrap_or("—".into()), r_mitm.map(|v| format!("{v:.2}")).unwrap_or("—".into())
        ));
        // Floors: T_enum(m, t) = 2^((m−1)t)/(m−1)! candidate sub-tuples, the
        // pair table 2^(⌊m/2⌋ t) lookups.
        let fact = |k: u64| (1..=k).map(|x| x as f64).product::<f64>();
        let t_enum = 2f64.powi(((*m - 1) * *t) as i32) / fact(*m - 1);
        let pair_floor = 2f64.powi(((*m / 2) * *t) as i32);
        let cand = median(&mut c.candidates.clone());
        let look = median(&mut c.lookups.clone());
        let r_enum = cand.map(|v| v / t_enum);
        let r_pair = look.and_then(|v| if v > 0.0 { Some(v / pair_floor) } else { None });
        prime_rows.push(json!({"m": m, "arm": arm, "t": t, "N": c.n, "base_points": c.base_points, "T_enum": t_enum, "pair_floor": pair_floor, "ratio_candidates_over_T_enum": r_enum, "ratio_lookups_over_pair_floor": r_pair, "mul": mul, "sqr": median(&mut c.sqr.clone()), "inv": median(&mut c.inv.clone()), "lookups": median(&mut c.lookups.clone()), "candidates": median(&mut c.candidates.clone()), "ns_median": ns_med, "ns_median_contended": ns_c, "ns_min": if ns_min.is_finite() { Some(ns_min) } else { None }, "ns_runs": c.ns.len(), "ns_contended_runs": c.ns_contended.len(), "ns_all": c.ns, "table_bytes": c.table_bytes, "setup_secs": median(&mut c.setup_secs.clone()), "ratio_mul_vs_pairs3": r_pairs, "ratio_mul_vs_mitm3": r_mitm}));
    }
    out.insert("prime".into(), Value::Array(prime_rows));
    let parms: Vec<(u64, String)> = prime
        .keys()
        .map(|k| (k.0, k.1.clone()))
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    for (m, arm) in &parms {
        for (qname, pick) in [("mul", 0), ("ns", 1)] {
            let pts: Vec<(f64, f64)> = prime
                .iter()
                .filter(|((mm, a, _), _)| mm == m && a == arm)
                .filter_map(|((_, _, t), c)| {
                    let v = if pick == 0 {
                        median(&mut c.mul.clone())
                    } else {
                        median(&mut c.ns.clone())
                    };
                    v.filter(|v| *v > 0.0).map(|v| (*t as f64, v.log2()))
                })
                .collect();
            if let Some((b, se, r2, n)) = slope(&pts) {
                md.push_str(&format!(
                    "| prime m={m} | {arm} | {qname} | {b:.3} | {se:.3} | {r2:.4} | {n} |\n"
                ));
                slopes.push(json!({"track": format!("prime m={m}"), "arm": arm, "quantity": qname, "slope_bits_per_rung": b, "se": se, "ci95_halfwidth": 1.96 * se, "r2": r2, "rungs": n, "points": pts}));
            }
        }
    }
    // F4 growth on the committed cells, finished rows only: log₂ of the
    // median multiply-adds (and ms) of a cell against t, per (kind, p, m).
    {
        let mut groups: BTreeMap<(String, u64, u64), Vec<(f64, f64, f64)>> = BTreeMap::new();
        for ((kind, p, m, t), rows) in &f4 {
            let mut mm: Vec<f64> = rows.iter().map(|r| r.0).collect();
            let mut ms: Vec<f64> = rows.iter().map(|r| r.1).collect();
            if let (Some(med), Some(msm)) = (median(&mut mm), median(&mut ms)) {
                if med > 0.0 && msm > 0.0 {
                    groups.entry((kind.clone(), *p, *m)).or_default().push((
                        *t as f64,
                        med.log2(),
                        msm.log2(),
                    ));
                }
            }
        }
        for ((kind, p, m), pts) in &groups {
            if pts.len() < 4 {
                continue;
            }
            for (qname, sel) in [("muladds", 1usize), ("ms", 2usize)] {
                let xy: Vec<(f64, f64)> = pts
                    .iter()
                    .map(|q| (q.0, if sel == 1 { q.1 } else { q.2 }))
                    .collect();
                if let Some((b, se, r2, n)) = slope(&xy) {
                    md.push_str(&format!(
                        "| prime m={m} (committed {kind} cells, p={p}) | f4_fp_tower | {qname} | {b:.3} | {se:.3} | {r2:.4} | {n} |\n"
                    ));
                    slopes.push(json!({"track": format!("prime m={m}"), "arm": format!("f4_fp_tower {kind} p={p}"), "quantity": qname, "slope_bits_per_rung": b, "se": se, "ci95_halfwidth": 1.96 * se, "r2": r2, "rungs": n, "points": xy}));
                }
            }
        }
    }
    out.insert("slopes".into(), Value::Array(slopes));

    md.push_str(&format!(
        "\n## Exactness census of the committed PKM rows\n\n`{}`\n",
        census
    ));
    out.insert("prime_rows_census".into(), census);

    // F4 ratios on the committed prime rows (finished F4 runs only; the
    // unfinished ones are censored lower bounds and are counted apart).
    let mut f4_rows = Vec::new();
    md.push_str("\n## Committed F4 rows against the root-finders on the same cells (medians over the rows of a cell)\n\n| kind | p | m | t | N | finished rows | unfinished rows | F4 multiply-adds | F4 ms | F4 degree | arm | mul per target | ns (row pass, contended) | F4 multiply-adds ÷ mul |\n|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|\n");
    for ((kind, p, m, t), rows) in &f4 {
        let mut mm: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let mut ms: Vec<f64> = rows.iter().map(|r| r.1).collect();
        let mut dd: Vec<f64> = rows.iter().map(|r| r.2).collect();
        let (mmed, msmed, dmed) = (
            median(&mut mm).unwrap(),
            median(&mut ms).unwrap(),
            median(&mut dd).unwrap_or(f64::NAN),
        );
        for ((k2, p2, m2, t2, arm), vals) in &rf_at_rows {
            if k2 == kind && p2 == p && m2 == m && t2 == t {
                let mut muls: Vec<f64> = vals.iter().map(|v| v.0).collect();
                let mut nss: Vec<f64> = vals.iter().map(|v| v.1).collect();
                let (mu, nsm) = (median(&mut muls).unwrap(), median(&mut nss).unwrap());
                let unf = f4_unfinished
                    .get(&(kind.clone(), *p, *m, *t))
                    .cloned()
                    .unwrap_or(0);
                md.push_str(&format!("| {kind} | {p} | {m} | {t} | {} | {} | {unf} | {mmed:.3e} | {msmed:.1} | {dmed} | {arm} | {mu:.0} | {nsm:.0} | {:.3e} |\n", m * t, rows.len(), mmed / mu.max(1.0)));
                f4_rows.push(json!({"kind": kind, "p": p, "m": m, "t": t, "N": m * t, "rows": rows.len(), "unfinished_rows": unf, "f4_muladds": mmed, "f4_ms": msmed, "f4_degree": dmed, "arm": arm, "mul": mu, "ns_contended": nsm, "ratio_muladds_over_mul": mmed / mu.max(1.0)}));
            }
        }
    }
    out.insert("f4_prime".into(), Value::Array(f4_rows));

    // Binary Gröbner references on (23, 8).
    let mut bref = Map::new();
    if !f6_xors.is_empty() {
        let mut x = f6_xors.clone();
        let mut m = f6_ms.clone();
        bref.insert(
            "inherited_f4_word_xors_median".into(),
            json!(median(&mut x)),
        );
        bref.insert("inherited_f4_ms_median".into(), json!(median(&mut m)));
        bref.insert("inherited_f4_targets".into(), json!(f6_xors.len()));
        bref.insert("inherited_f4_word_xors_all".into(), json!(f6_xors));
        bref.insert("inherited_f4_ms_all".into(), json!(f6_ms));
        md.push_str(&format!("\n## Binary Gröbner references on icv1-f2m23-tm5197-1f85e9e1, ℓ = 8, m = 3\n\n- inherited F4 (`f6_ic_probe`, default decomposition oracle): {} targets, median {:.3e} word XORs and {:.1} ms per target.\n", f6_xors.len(), median(&mut f6_xors.clone()).unwrap(), median(&mut f6_ms.clone()).unwrap()));
    }
    bref.insert("split_dreg_rows".into(), Value::Array(split_rows.clone()));
    if split_rows.is_empty() {
        md.push_str("- split-dreg F4 (unsplit S₄, g = 0): no solve row (see the run's stderr: killed or budget).\n");
    } else {
        md.push_str(&format!("- split-dreg F4 rows: {}\n", split_rows.len()));
    }
    out.insert("binary_groebner_refs".into(), Value::Object(bref));

    let json_out = Value::Object(out);
    std::fs::write(
        Path::new(&dir).join("ANALYSIS.json"),
        serde_json::to_string_pretty(&json_out).unwrap(),
    )
    .unwrap();
    std::fs::write(Path::new(&dir).join("ANALYSIS.md"), md).unwrap();
    println!("wrote {}/ANALYSIS.json and ANALYSIS.md", dir);
}
