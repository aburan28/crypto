//! Apply the frozen n41/n53 shared-rank decision to sealed ecbench records.
//!
//! An operation-count diagnostic (`research/ecbench_n41_n53_shared_rank_20261005/PROTOCOL.md`).
//! Counted quotients are lower bounds while native work is unpriced; every
//! wall figure from the macOS session is L0 and stays exploratory. Nothing
//! here can promote an IC/rho speedup.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::record::Record;
use crypto_lib::cryptanalysis::ecbench::runner::{read_records, read_session};
use crypto_lib::cryptanalysis::ecbench::stats::cluster_bootstrap_ci;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

const SPEC_ID: &str = "ECS1hb5e99ec6937e";
const ARMS: [&str; 4] = ["rho-strong", "ic-k8", "ic-k16", "ic-k16-control"];
const IC_ARMS: [&str; 3] = ["ic-k8", "ic-k16", "ic-k16-control"];
const CURVES: [(&str, u32, u32, &str); 2] = [
    ("icv1-f2m41-tm2308219-7f48b14a", 41, 82, "W15dcf1a74cb8"),
    ("icv1-f2m53-tm56619371-dac20a85", 53, 106, "W32013b060ceb"),
];
const TARGETS_PER_CURVE: usize = 16;
const ROUNDS: u32 = 5;
const BOOT_SEED: u64 = 202610055103;
const RESAMPLES: usize = 20000;
const COMPARISONS: [(&str, &str); 4] = [
    ("rho-strong", "ic-k8"),
    ("rho-strong", "ic-k16"),
    ("ic-k8", "ic-k16"),
    ("ic-k16", "ic-k16-control"),
];
const ONLINE_PHASES: [&str; 5] = [
    "target_query",
    "target_PDP",
    "target_relation_check",
    "target_descent",
    "target_recovery_check",
];

fn bytes(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn value(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&bytes(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn ratio_of_sums(pairs: &[(f64, f64)]) -> Option<f64> {
    let (sa, sb) = pairs
        .iter()
        .fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
    (!pairs.is_empty() && sa > 0.0).then(|| sb / sa)
}

/// Ratio of sums `Σ b / Σ a` over workload strata with the frozen
/// 20,000-resample two-stage percentile interval.
fn strata_ratio(strata: &BTreeMap<String, Vec<(f64, f64)>>) -> Value {
    let s: Vec<Vec<(f64, f64)>> = strata.values().filter(|v| !v.is_empty()).cloned().collect();
    let flat: Vec<(f64, f64)> = s.iter().flatten().copied().collect();
    let ratio = ratio_of_sums(&flat);
    let ci = if s.len() >= 2 {
        cluster_bootstrap_ci(&s, RESAMPLES, BOOT_SEED, |x| {
            ratio_of_sums(&x.iter().flatten().copied().collect::<Vec<_>>())
        })
    } else {
        None
    };
    json!({
        "ratio_of_sums": ratio,
        "ci95": ci.map(|(lo, hi)| json!([lo, hi])),
        "ci_method": if s.len() >= 2 { "cluster" } else { "none" },
        "pairs": flat.len(),
        "workloads": s.len(),
        "resamples": RESAMPLES,
        "seed": BOOT_SEED,
    })
}

fn verdict(ci: Option<(f64, f64)>, keep_below_one: bool) -> &'static str {
    match ci {
        None => "not_evaluable",
        Some((lo, hi)) => {
            if lo > 1.0 {
                if keep_below_one { "falsified" } else { "retained" }
            } else if hi < 1.0 {
                if keep_below_one { "retained" } else { "falsified" }
            } else {
                "undecided"
            }
        }
    }
}

fn ci_of(v: &Value) -> Option<(f64, f64)> {
    let a = v.as_array()?;
    Some((a.first()?.as_f64()?, a.get(1)?.as_f64()?))
}

fn checked_receipt(path: &Path, session_dir: &Path, session: &str, expected_replays: usize) -> Result<Value, String> {
    let receipt = value(path)?;
    if receipt["ok"] != true
        || receipt["session_status"] != "complete"
        || receipt["session_id"] != session
        || receipt["records"] != 768
    {
        return Err(format!("{} does not audit the frozen session", path.display()));
    }
    let replays = receipt["replays"].as_array().ok_or("audit replays missing")?;
    if replays.len() != expected_replays || replays.iter().any(|r| r["reproduced"] != true) {
        return Err(format!(
            "{} replayed {} runs ({} expected) or failed to reproduce one",
            path.display(),
            replays.len(),
            expected_replays
        ));
    }
    for name in ["host.json", "plan.json", "records.jsonl", "session.json", "spec.json"] {
        let actual = sha256_hex(&bytes(&session_dir.join(name))?);
        if receipt["files"][name] != actual {
            return Err(format!("{} has a wrong {name} hash", path.display()));
        }
    }
    Ok(receipt)
}

fn verified(r: &Record) -> bool {
    r.outcome.status == "verified" && r.outcome.matches_target == Some(true)
}

fn online_gae(r: &Record) -> f64 {
    r.phases.iter().filter(|p| p.name.starts_with("target_")).map(|p| p.gae).sum()
}

fn checked_record(r: &Record, session: &str) -> Result<(), String> {
    if r.session_id != session || r.warmup != (r.round == 0) || !ARMS.contains(&r.arm.as_str()) {
        return Err(format!("{} is not a record of the frozen session", r.run_id));
    }
    if !verified(r) {
        if r.cost.s.is_some() && r.outcome.status != "verified" {
            return Err(format!("{} has a cost without a verified answer", r.run_id));
        }
        return Ok(());
    }
    if !r.cost.lower_bound || r.cost.unpriced.is_empty() || !r.cost.deterministic {
        return Err(format!("{} fails the accounting gate", r.run_id));
    }
    let charged = r.cost.total_gae.ok_or("missing cold GAE")?;
    let phase_sum: f64 = r.phases.iter().map(|p| p.gae).sum();
    if (charged - phase_sum).abs() > 1e-7_f64.max(charged * 1e-11) {
        return Err(format!("{} has nonexclusive charged phases", r.run_id));
    }
    let online = r.online.as_ref().ok_or_else(|| format!("{} lacks an online interval", r.run_id))?;
    if online.wall_ns == 0 || online.phases_ns.values().sum::<u64>() != online.wall_ns {
        return Err(format!("{} has incomplete exclusive online phases", r.run_id));
    }
    if r.arm.starts_with("ic-") {
        let names: BTreeSet<&str> = online.phases_ns.keys().map(String::as_str).collect();
        let expected: BTreeSet<&str> = ONLINE_PHASES.into_iter().collect();
        if names != expected || r.detail["rank"]["verified"] != true {
            return Err(format!("{} lacks full-rank IC phases", r.run_id));
        }
    }
    Ok(())
}

fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

fn arm_summary(arm: &str, measured: &[&Record]) -> Value {
    let ok: Vec<&Record> = measured.iter().copied().filter(|r| verified(r)).collect();
    let mut statuses = BTreeMap::<String, u64>::new();
    let mut errors = BTreeMap::<String, u64>::new();
    let mut levels = BTreeMap::<String, u64>::new();
    for r in measured {
        *statuses.entry(r.outcome.status.clone()).or_default() += 1;
        *levels.entry(r.isolation.level.name().to_string()).or_default() += 1;
        if !verified(r) {
            *errors
                .entry(r.outcome.error.clone().unwrap_or_else(|| r.outcome.status.clone()))
                .or_default() += 1;
        }
    }
    let gae: Vec<f64> = ok.iter().filter_map(|r| r.cost.total_gae).collect();
    let s: Vec<f64> = ok.iter().filter_map(|r| r.cost.s).collect();
    let online_ns: Vec<f64> = ok.iter().filter_map(|r| r.online.as_ref().map(|o| o.wall_ns as f64)).collect();
    let solve_ns: Vec<f64> = ok.iter().filter_map(|r| r.time.solve_wall_ns.map(|x| x as f64)).collect();
    let failed_process_ns: Vec<f64> = measured
        .iter()
        .filter(|r| !verified(r))
        .map(|r| r.time.process_wall_ns as f64)
        .collect();
    let max_rss = measured.iter().map(|r| r.time.max_rss_kib).max().unwrap_or(0);
    let mut phase_mean = Map::new();
    if !ok.is_empty() {
        let mut sums = BTreeMap::<String, f64>::new();
        for r in &ok {
            for p in &r.phases {
                *sums.entry(p.name.clone()).or_default() += p.gae;
            }
        }
        for (k, v) in sums {
            phase_mean.insert(k, json!(v / ok.len() as f64));
        }
    }
    let mut row = json!({
        "arm": arm,
        "method_id": measured.first().map(|r| r.method.method_id.clone()),
        "measured": measured.len(),
        "verified": ok.len(),
        "statuses": statuses,
        "errors": errors,
        "isolation_levels": levels,
        "timeouts": statuses.get("timeout").copied().unwrap_or(0),
        "total_cold_gae_lower_bound": (!gae.is_empty()).then(|| gae.iter().sum::<f64>()),
        "mean_cold_gae_lower_bound": mean(&gae),
        "mean_cold_s_lower_bound": mean(&s),
        "mean_online_wall_ns_l0_exploratory": mean(&online_ns),
        "mean_solve_wall_ns_l0_exploratory": mean(&solve_ns),
        "mean_failed_process_wall_ns_l0_exploratory": mean(&failed_process_ns),
        "max_rss_kib": max_rss,
        "phase_mean_gae_lower_bound": phase_mean,
    });
    if arm.starts_with("ic-") {
        let first = ok.first();
        let fb = first.and_then(|r| r.factor_base.as_ref());
        let attempts: Vec<f64> = ok.iter().map(|r| *r.counters.get("target_attempts").unwrap_or(&0) as f64).collect();
        let pair_entries = first.and_then(|r| {
            r.phases
                .iter()
                .find(|p| p.name == "oracle_setup")
                .and_then(|p| p.native.get("pair_table_entries").copied())
        });
        let target: Vec<f64> = ok.iter().map(|r| online_gae(r)).collect();
        let reusable: Vec<f64> = ok
            .iter()
            .map(|r| r.cost.total_gae.unwrap_or(0.0) - online_gae(r))
            .collect();
        let extra = json!({
            "usable_points": fb.map(|f| f.signed_points),
            "folded_columns": fb.map(|f| f.columns),
            "factor_base_id": fb.map(|f| f.fb_id.clone()),
            "factor_base_points_sha256": fb.map(|f| f.points_sha256.clone()),
            "pair_table_entries": pair_entries,
            "rank_trials": first.and_then(|r| r.counters.get("rank_trials").copied()),
            "rank_hits": first.and_then(|r| r.counters.get("rank_hits").copied()),
            "matrix_rank": first.and_then(|r| r.counters.get("matrix_rank").copied()),
            "mean_target_attempts_including_failed": mean(&attempts),
            "target_attempts_total": attempts.iter().sum::<f64>(),
            "mean_target_online_gae_lower_bound": mean(&target),
            "mean_reusable_setup_gae_lower_bound": mean(&reusable),
        });
        for (k, v) in extra.as_object().unwrap() {
            row[k] = v.clone();
        }
    }
    row
}

fn comparison_row(session_dir: &Path, a: &str, b: &str) -> Result<Value, String> {
    let path = session_dir.join(format!("comparisons/{a}__{b}.json"));
    let row = value(&path)?;
    if row["ops"]["resamples"] != RESAMPLES
        || row["ops"]["same_session"] != true
        || row["ops"]["bounded"] != true
        || row["wall"]["status"] != "descriptive"
    {
        return Err(format!("{} is not the frozen 20,000-resample bounded L0 comparison", path.display()));
    }
    Ok(row)
}

fn curve_block(
    session_dir: &Path,
    slug: &str,
    n: u32,
    automorphisms: u32,
    primary: &str,
    by_pair: &BTreeMap<(String, u32), BTreeMap<String, &Record>>,
) -> Result<Value, String> {
    let mut by_arm = BTreeMap::<String, Vec<&Record>>::new();
    for arms in by_pair.values() {
        for (arm, r) in arms {
            by_arm.entry(arm.clone()).or_default().push(r);
        }
    }
    let floor_s = by_arm["rho-strong"][0].boundaries.floor_s;
    let r_order = by_arm["rho-strong"][0].workload.curve.r;
    for r in by_arm.values().flatten() {
        if r.boundaries.floor_s != floor_s
            || r.boundaries.automorphisms_available != automorphisms
            || r.workload.curve.field_degree != Some(n)
            || r.workload.curve.slug != slug
        {
            return Err(format!("{slug}: a record disagrees on floor, automorphisms or curve"));
        }
    }
    let arms: Vec<Value> = ARMS.iter().map(|a| arm_summary(a, &by_arm[*a])).collect();

    // Primary target-zero row: the mean over five paired rounds.
    let mut rounds = Vec::new();
    let mut prim = BTreeMap::<&str, Vec<(f64, f64)>>::new();
    let mut prim_wall = BTreeMap::<&str, Vec<(f64, f64)>>::new();
    for round in 1..=ROUNDS {
        let arms_r = by_pair
            .get(&(primary.to_string(), round))
            .ok_or_else(|| format!("{slug}: primary workload {primary} lacks round {round}"))?;
        let g = |arm: &str| arms_r[arm].cost.total_gae.filter(|_| verified(arms_r[arm]));
        let w = |arm: &str| arms_r[arm].online.as_ref().filter(|_| verified(arms_r[arm])).map(|o| o.wall_ns as f64);
        for (key, a, b) in [("k8_over_rho", "rho-strong", "ic-k8"), ("k16_over_rho", "rho-strong", "ic-k16"), ("k16_over_k8", "ic-k8", "ic-k16"), ("k16_over_control", "ic-k16-control", "ic-k16")] {
            if let (Some(x), Some(y)) = (g(a), g(b)) {
                prim.entry(key).or_default().push((x, y));
            }
            if let (Some(x), Some(y)) = (w(a), w(b)) {
                prim_wall.entry(key).or_default().push((x, y));
            }
        }
        rounds.push(json!({
            "round": round,
            "algorithm_seed": arms_r["rho-strong"].algorithm_seed,
            "status": ARMS.iter().map(|a| (a.to_string(), json!(arms_r[*a].outcome.status))).collect::<Map<_, _>>(),
            "cold_gae_lower_bound": ARMS.iter().map(|a| (a.to_string(), json!(g(a)))).collect::<Map<_, _>>(),
            "online_wall_ns_l0_exploratory": ARMS.iter().map(|a| (a.to_string(), json!(w(a)))).collect::<Map<_, _>>(),
        }));
    }
    let prim_json = |m: &BTreeMap<&str, Vec<(f64, f64)>>| -> Value {
        ["k8_over_rho", "k16_over_rho", "k16_over_k8", "k16_over_control"]
            .iter()
            .map(|k| {
                let pairs = m.get(k).cloned().unwrap_or_default();
                (
                    k.to_string(),
                    json!({
                        "paired_rounds": pairs.len(),
                        "complete": pairs.len() as u32 == ROUNDS,
                        "ratio_of_sums": ratio_of_sums(&pairs),
                    }),
                )
            })
            .collect::<Map<_, _>>()
            .into()
    };

    // Secondary: the harness's own per-curve comparison rows, plus the
    // target-only (H2) and L0 online-wall ratios with the same bootstrap.
    let mut secondary = Map::new();
    for (a, b) in COMPARISONS {
        let row = comparison_row(session_dir, a, b)?;
        let curve = row["curves"]
            .as_array()
            .and_then(|cs| cs.iter().find(|c| c["slug"] == slug))
            .cloned()
            .unwrap_or(Value::Null);
        secondary.insert(
            format!("{b}_over_{a}"),
            json!({
                "numerator": b,
                "denominator": a,
                "unit": "cold counted S lower bound (ecbench.gae / sqrt r)",
                "overall_status": row["ops"]["status"],
                "ratio_of_sums": curve["ratio_b_over_a"],
                "ci95": curve["ci95"],
                "ci_method": curve["ci_method"],
                "pairs": curve["pairs"],
                "workloads": curve["workloads"],
                "resamples": RESAMPLES,
                "seed": BOOT_SEED,
                "comparison_id": row["comparison_id"],
            }),
        );
    }
    let mut h2_strata = BTreeMap::<String, Vec<(f64, f64)>>::new();
    let mut wall_strata = BTreeMap::<&str, BTreeMap<String, Vec<(f64, f64)>>>::new();
    for ((wid, _), arms_r) in by_pair {
        let k8 = arms_r["ic-k8"];
        let k16 = arms_r["ic-k16"];
        if verified(k8) && verified(k16) {
            h2_strata.entry(wid.clone()).or_default().push((online_gae(k8), online_gae(k16)));
        }
        for (key, a, b) in [("k8_over_k16", "ic-k16", "ic-k8"), ("rho_over_k8", "ic-k8", "rho-strong"), ("rho_over_k16", "ic-k16", "rho-strong"), ("k16_over_control", "ic-k16-control", "ic-k16")] {
            let (ra, rb) = (arms_r[a], arms_r[b]);
            if verified(ra) && verified(rb) {
                let (oa, ob) = (ra.online.as_ref().unwrap().wall_ns as f64, rb.online.as_ref().unwrap().wall_ns as f64);
                wall_strata.entry(key).or_default().entry(wid.clone()).or_default().push((oa, ob));
            }
        }
    }
    let h2 = strata_ratio(&h2_strata);
    let h2_verdict = verdict(ci_of(&h2["ci95"]), true);
    let wall_secondary: Map<String, Value> = ["k8_over_k16", "rho_over_k8", "rho_over_k16", "k16_over_control"]
        .iter()
        .map(|k| {
            let mut v = strata_ratio(wall_strata.get(k).unwrap_or(&BTreeMap::new()));
            v["unit"] = json!("online_wall_ns; L0 exploratory, never admitted");
            (k.to_string(), v)
        })
        .collect();
    let aa_max_dev = wall_strata
        .get("k16_over_control")
        .map(|m| {
            m.values()
                .flatten()
                .map(|(a, b)| ((b - a) / a).abs())
                .fold(0.0_f64, f64::max)
        });

    let h1: Map<String, Value> = [("ic-k8", "ic-k8_over_rho-strong"), ("ic-k16", "ic-k16_over_rho-strong")]
        .iter()
        .map(|(arm, key)| {
            let ci = ci_of(&secondary[*key]["ci95"]);
            (arm.to_string(), json!({"verdict": verdict(ci, false), "ci95": secondary[*key]["ci95"], "ratio_of_sums": secondary[*key]["ratio_of_sums"]}))
        })
        .collect();

    Ok(json!({
        "curve": slug,
        "field_degree": n,
        "subgroup_order": r_order,
        "log2_r": (r_order as f64).log2(),
        "automorphisms_available": automorphisms,
        "generic_floor_s": floor_s,
        "primary_workload_id": primary,
        "arms": arms,
        "primary_one_target": {
            "workload_id": primary,
            "target_index": 0,
            "rounds": rounds,
            "counted": prim_json(&prim),
            "online_wall_l0_exploratory": prim_json(&prim_wall),
        },
        "secondary_counted": secondary,
        "h1_counted_ic_over_rho_above_one": h1,
        "h2_target_only_counted_k16_over_k8": {"verdict": h2_verdict, "statistic": h2},
        "secondary_online_wall_l0_exploratory": wall_secondary,
        "k16_aa_max_relative_online_wall_deviation_l0": aa_max_dev,
    }))
}

fn run(round_dir: &Path, session_dir: &Path, audit_path: &Path, out_path: &Path) -> Result<(), String> {
    let freeze = value(&round_dir.join("FREEZE.json"))?;
    let spec_sha = sha256_hex(&bytes(&round_dir.join("SPEC.json"))?);
    let plan_sha = sha256_hex(&bytes(&round_dir.join("PLAN.json"))?);
    if freeze["spec_file_sha256"] != spec_sha || freeze["plan_file_sha256"] != plan_sha || freeze["spec_id"] != SPEC_ID {
        return Err("SPEC.json or PLAN.json differs from FREEZE.json".into());
    }
    let session = read_session(session_dir)?;
    if session.spec_id != SPEC_ID || session.status != "complete" {
        return Err("session is not a complete run of the frozen spec".into());
    }
    let session_spec_sha = sha256_hex(&bytes(&session_dir.join("spec.json"))?);
    if session_spec_sha != spec_sha {
        return Err("the session's spec.json is not the frozen SPEC.json byte for byte".into());
    }
    let records = read_records(session_dir)?;
    if records.len() != 768 || records.iter().filter(|r| r.warmup).count() != 128 {
        return Err("expected 128 warm-ups and 640 measured executions".into());
    }
    for r in &records {
        checked_record(r, &session.session_id)?;
    }
    let verified_measured = records.iter().filter(|r| !r.warmup && verified(r)).count();
    let audit = checked_receipt(audit_path, session_dir, &session.session_id, verified_measured)?;

    let mut by_curve = BTreeMap::<String, BTreeMap<(String, u32), BTreeMap<String, &Record>>>::new();
    for r in records.iter().filter(|r| !r.warmup) {
        by_curve
            .entry(r.workload.curve.slug.clone())
            .or_default()
            .entry((r.workload.workload_id.clone(), r.round))
            .or_default()
            .insert(r.arm.clone(), r);
    }
    let expected_arms: BTreeSet<&str> = ARMS.into_iter().collect();
    for (slug, pairs) in &by_curve {
        if pairs.len() != TARGETS_PER_CURVE * ROUNDS as usize {
            return Err(format!("{slug}: expected {} (workload, round) pairs", TARGETS_PER_CURVE * ROUNDS as usize));
        }
        for ((wid, round), arms) in pairs {
            if arms.keys().map(String::as_str).collect::<BTreeSet<_>>() != expected_arms {
                return Err(format!("{wid} round {round} lacks an arm"));
            }
            let rho = arms["rho-strong"];
            for r in arms.values() {
                if r.workload.target != rho.workload.target || r.algorithm_seed != rho.algorithm_seed {
                    return Err(format!("{wid} round {round} is not a paired solve"));
                }
                if verified(r) && verified(rho) && r.outcome.recovered != rho.outcome.recovered {
                    return Err(format!("{wid} round {round}: two verified arms disagree on the logarithm"));
                }
            }
        }
    }
    let mut curves = Vec::new();
    for (slug, n, autos, primary) in CURVES {
        let pairs = by_curve.get(slug).ok_or_else(|| format!("{slug} missing from the session"))?;
        curves.push(curve_block(session_dir, slug, n, autos, primary, pairs)?);
    }

    // Decision, exactly as preregistered.
    let mut any_cell_zero = false;
    let mut all_arms_all_verified = true;
    let mut intervals = Vec::new();
    let mut aa_counted_ok = true;
    let mut cells = Vec::new();
    for c in &curves {
        for arm in c["arms"].as_array().unwrap() {
            let name = arm["arm"].as_str().unwrap();
            if IC_ARMS.contains(&name) {
                let v = arm["verified"].as_u64().unwrap();
                let m = arm["measured"].as_u64().unwrap();
                if v == 0 { any_cell_zero = true; }
                if v != m { all_arms_all_verified = false; }
            }
            cells.push(json!({"curve": c["curve"], "arm": name, "measured": arm["measured"], "verified": arm["verified"], "timeouts": arm["timeouts"], "errors": arm["errors"]}));
        }
        for key in ["ic-k8_over_rho-strong", "ic-k16_over_rho-strong"] {
            if let Some(ci) = ci_of(&c["secondary_counted"][key]["ci95"]) {
                intervals.push(ci);
            }
        }
        let aa = &c["secondary_counted"]["ic-k16-control_over_ic-k16"];
        if aa["pairs"].as_u64().unwrap_or(0) > 0 && aa["ratio_of_sums"] != 1.0 {
            aa_counted_ok = false;
        }
    }
    if !aa_counted_ok {
        return Err("the identical K16 control changed counted work".into());
    }
    let any_below = intervals.iter().any(|(_, hi)| *hi < 1.0);
    let all_above = !intervals.is_empty() && intervals.iter().all(|(lo, _)| *lo > 1.0);
    let decision = if any_below {
        "counted_lead_candidate"
    } else if all_arms_all_verified && all_above {
        "counted_quotients_above_one"
    } else if any_cell_zero && all_above {
        "frozen_parameters_do_not_transfer"
    } else {
        "undecided"
    };
    let mut levels = BTreeMap::<String, u64>::new();
    for r in records.iter().filter(|r| !r.warmup) {
        *levels.entry(r.isolation.level.name().to_string()).or_default() += 1;
    }
    let output = json!({
        "schema": "ecbench.n41_n53_shared_rank_decision/v1",
        "status": "complete_l0_counted_diagnostic",
        "class": "accounting",
        "decision": decision,
        "protocol_source_commit": freeze["source_commit_at_freeze"],
        "spec_id": SPEC_ID,
        "spec_sha256": spec_sha,
        "plan_sha256": plan_sha,
        "session_id": session.session_id,
        "session_binary_sha256": session.binary_sha256,
        "session_git_commit": session.git_commit,
        "session_env_class_id": session.env_class_id,
        "records_sha256": sha256_hex(&bytes(&session_dir.join("records.jsonl"))?),
        "audit_receipt_sha256": sha256_hex(&bytes(audit_path)?),
        "auditor_env_class_id": audit["auditor_env_class_id"],
        "replays_reproduced": audit["replays"].as_array().map(Vec::len),
        "measured_runs": 640,
        "measured_verified_runs": verified_measured,
        "isolation_levels": levels,
        "all_measured_L2": false,
        "cells": cells,
        "curves": curves,
        "online_speedup": Value::Null,
        "online_speedup_reason": "every run earned L0 on macOS (no affinity control); wall figures are exploratory and the L2 host run in L2_RUNBOOK.md is pending",
        "fully_priced_cold_speedup": Value::Null,
        "limits": [
            "counted IC and rho totals are lower bounds: field arithmetic, hashing, allocation, table probes and modular combination are counted but unpriced",
            "a ratio of two lower bounds does not bound the true IC/rho cost",
            "all wall intervals are L0 exploratory on a shared, loaded host",
            "the n37 method parameters were carried unchanged; a cell whose arm never verifies is unknown, not a cost",
            "n41 and n53 are not n83 or ECC2K-130 results",
        ],
    });
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out_path)
        .map_err(|e| format!("{}: {e}", out_path.display()))?;
    serde_json::to_writer_pretty(&mut out, &output).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let usage = "usage: ecbench_n41_n53_shared_rank_analyze ROUND_DIR SESSION_DIR AUDIT_JSON OUT_JSON";
    let round = args.next().ok_or(usage)?;
    let session = args.next().ok_or(usage)?;
    let audit = args.next().ok_or(usage)?;
    let out = args.next().ok_or(usage)?;
    if args.next().is_some() {
        return Err(usage.into());
    }
    run(Path::new(&round), Path::new(&session), Path::new(&audit), Path::new(&out))
}
