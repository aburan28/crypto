//! Apply the frozen n37 rank-column decision to sealed ecbench records.
//! This is an operation-count diagnostic. L0 wall ratios and unpriced native
//! work can never be promoted to an IC speedup by this analyzer.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::claim::id_sha256;
use crypto_lib::cryptanalysis::ecbench::record::Record;
use crypto_lib::cryptanalysis::ecbench::runner::{read_records, read_session};
use crypto_lib::cryptanalysis::ecbench::spec::Level;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

const KS: [u64; 6] = [4, 8, 12, 16, 24, 42];
const ARMS: [&str; 8] = [
    "rho-strong",
    "ic-k4",
    "ic-k8",
    "ic-k12",
    "ic-k16",
    "ic-k24",
    "ic-k42",
    "ic-k42-control",
];
const SESSION: &str = "ECBS1h9a1d373bbe54";
const SPEC: &str = "ECS1h76742ce139b5";
const SOURCE: &str = "694ac8706d6b5d7ee4083b9fe195f9a25efeecd6";

fn bytes(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn value(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&bytes(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn required_f64(v: &Value, key: &str) -> Result<f64, String> {
    v[key].as_f64().ok_or_else(|| format!("missing {key}"))
}

fn checked_receipt(path: &Path, session_dir: &Path, session: &str) -> Result<Value, String> {
    let receipt = value(path)?;
    if receipt["ok"] != true
        || receipt["session_status"] != "complete"
        || receipt["session_id"] != session
        || receipt["records"] != 384
        || receipt["verified_records"] != 384
    {
        return Err(format!(
            "{} does not audit the frozen session",
            path.display()
        ));
    }
    let replays = receipt["replays"]
        .as_array()
        .ok_or("audit replays missing")?;
    if replays.len() != 320 || replays.iter().any(|r| r["reproduced"] != true) {
        return Err(format!(
            "{} did not reproduce all 320 measured runs",
            path.display()
        ));
    }
    for name in [
        "host.json",
        "plan.json",
        "records.jsonl",
        "session.json",
        "spec.json",
    ] {
        let actual = sha256_hex(&bytes(&session_dir.join(name))?);
        if receipt["files"][name] != actual {
            return Err(format!("{} has a wrong {name} hash", path.display()));
        }
    }
    Ok(receipt)
}

fn comparison(session_dir: &Path, a: &str, b: &str) -> Result<Value, String> {
    let path = session_dir.join(format!("comparisons/{a}__{b}.json"));
    let row = value(&path)?;
    if row["ops"]["status"] != "ok"
        || row["ops"]["bounded"] != true
        || row["ops"]["pairs"] != 40
        || row["ops"]["same_seeds"] != true
        || row["ops"]["same_session"] != true
        || row["wall"]["status"] != "descriptive"
    {
        return Err(format!(
            "{} is not a paired bounded L0 comparison",
            path.display()
        ));
    }
    Ok(row)
}

fn checked_record(record: &Record, session: &str) -> Result<(), String> {
    if record.session_id != session
        || record.outcome.status != "verified"
        || record.outcome.matches_target != Some(true)
        || !record.cost.lower_bound
        || record.cost.unpriced.is_empty()
        || !record.cost.deterministic
        || record.isolation.level != Level::L0
        || record.warmup != (record.round == 0)
    {
        return Err(format!("{} fails status/accounting gate", record.run_id));
    }
    let online = record
        .online
        .as_ref()
        .ok_or_else(|| format!("{} lacks an online interval", record.run_id))?;
    if online.wall_ns == 0 || online.phases_ns.values().sum::<u64>() != online.wall_ns {
        return Err(format!(
            "{} has incomplete exclusive online phases",
            record.run_id
        ));
    }
    if record.arm.starts_with("ic-") {
        let names: BTreeSet<&str> = online.phases_ns.keys().map(String::as_str).collect();
        let expected: BTreeSet<&str> = [
            "target_query",
            "target_PDP",
            "target_relation_check",
            "target_descent",
            "target_recovery_check",
        ]
        .into_iter()
        .collect();
        if names != expected || record.detail["rank"]["verified"] != true {
            return Err(format!("{} lacks full-rank IC phases", record.run_id));
        }
    }
    let charged = record.cost.total_gae.ok_or("missing cold GAE")?;
    let phase_sum: f64 = record.phases.iter().map(|p| p.gae).sum();
    if (charged - phase_sum).abs() > 1e-7_f64.max(charged * 1e-11) {
        return Err(format!("{} has nonexclusive charged phases", record.run_id));
    }
    Ok(())
}

fn arm_row(
    session_dir: &Path,
    claims_dir: &Path,
    records: &[&Record],
    k: u64,
    floor_s: f64,
) -> Result<Value, String> {
    let name = format!("ic-k{k}");
    if records.len() != 40 || records.iter().any(|r| r.arm != name) {
        return Err(format!("{name} does not have 40 measured runs"));
    }
    let first = records[0];
    let fb = first
        .factor_base
        .as_ref()
        .ok_or("factor-base facts missing")?;
    if fb.columns != k || fb.signed_points != 74 * k {
        return Err(format!("{name} did not build the predicted usable support"));
    }
    let base_digest = first.detail["rank"]["base_sha256"]
        .as_str()
        .ok_or("base digest missing")?;
    for r in records {
        if r.factor_base.as_ref() != Some(fb)
            || r.detail["rank"]["base_sha256"] != base_digest
            || r.detail["rank"]["rank"] != k
            || r.detail["rank"]["column_logs"].as_array().map(Vec::len) != Some(k as usize)
            || r.detail["targets"][0]["verified"] != true
        {
            return Err(format!("{name} changed base, rank, or target status"));
        }
    }
    let claim_path = claims_dir.join(format!("{name}-first.json"));
    let claim_bytes = bytes(&claim_path)?;
    let claim: Value = serde_json::from_slice(&claim_bytes).map_err(|e| e.to_string())?;
    let manifest_hash = id_sha256(&claim["candidate_manifest"])?;
    let candidate_id = claim["candidate_id"]
        .as_str()
        .ok_or("candidate ID missing")?;
    if claim["candidate_manifest_sha256"] != manifest_hash
        || !candidate_id.ends_with(&format!("h{}", &manifest_hash[..12]))
        || claim["candidate_manifest"]["factor_base"]["inventory"]["usable_point_count"]
            != fb.signed_points
        || claim["candidate_manifest"]["factor_base"]["inventory"]["effective_columns"] != k
        || claim["ecbench_records"]["session"] != SESSION
        || claim["verdict"]
            .as_str()
            .is_none_or(|s| !s.starts_with("descriptive only"))
    {
        return Err(format!("{name} candidate manifest does not match session"));
    }

    let mut phase_sum = BTreeMap::<String, f64>::new();
    let mut total_gae = 0.0;
    let mut total_s = 0.0;
    let mut attempts = 0usize;
    let mut peak_rss_kib = 0u64;
    for r in records {
        total_gae += r.cost.total_gae.ok_or("missing cold GAE")?;
        total_s += r.cost.s.ok_or("missing S")?;
        peak_rss_kib = peak_rss_kib.max(r.time.max_rss_kib);
        attempts += r.detail["targets"][0]["attempts"]
            .as_array()
            .ok_or("target attempts missing")?
            .len();
        for phase in &r.phases {
            *phase_sum.entry(phase.name.clone()).or_default() += phase.gae;
        }
    }
    let mean_s = total_s / 40.0;
    let phase_mean: BTreeMap<String, f64> = phase_sum
        .into_iter()
        .map(|(name, sum)| (name, sum / 40.0))
        .collect();
    let reusable_setup_gae: f64 = phase_mean
        .iter()
        .filter(|(name, _)| !name.starts_with("target_"))
        .map(|(_, cost)| cost)
        .sum();
    let online_gae: f64 = phase_mean
        .iter()
        .filter(|(name, _)| name.starts_with("target_"))
        .map(|(_, cost)| cost)
        .sum();
    let rho = comparison(session_dir, "rho-strong", &name)?;
    let baseline = if k == 42 {
        json!({"ratio_b_over_a":1.0,"ci95":[1.0,1.0]})
    } else {
        comparison(session_dir, "ic-k42", &name)?["ops"].clone()
    };
    Ok(json!({
        "arm": name,
        "candidate_id": candidate_id,
        "candidate_manifest_sha256": manifest_hash,
        "diagnostic_claim_sha256": sha256_hex(&claim_bytes),
        "factor_base_id": fb.fb_id,
        "factor_base_points_sha256": fb.points_sha256,
        "rank_base_sha256": base_digest,
        "usable_points": fb.signed_points,
        "folded_columns": fb.columns,
        "rank_trials": first.detail["rank"]["trials"],
        "rank_hits": first.detail["rank"]["hits"],
        "measured_runs": 40,
        "verified_runs": 40,
        "total_cold_gae_lower_bound": total_gae,
        "mean_cold_s_lower_bound": mean_s,
        "s_over_generic_floor_lower_bound": mean_s / floor_s,
        "mean_reusable_setup_gae_lower_bound": reusable_setup_gae,
        "mean_target_online_gae_lower_bound": online_gae,
        "mean_target_attempts": attempts as f64 / 40.0,
        "target_attempts_total": attempts,
        "peak_rss_kib_max": peak_rss_kib,
        "phase_mean_gae_lower_bound": phase_mean,
        "cold_counted_over_rho_diagnostic": rho["ops"]["ratio_b_over_a"],
        "cold_counted_over_rho_ci95_diagnostic": rho["ops"]["ci95"],
        "cold_counted_over_k42_diagnostic": baseline["ratio_b_over_a"],
        "cold_counted_over_k42_ci95_diagnostic": baseline["ci95"],
    }))
}

fn run(
    session_dir: &Path,
    claims_dir: &Path,
    local_audit_path: &Path,
    independent_audit_path: &Path,
    out_path: &Path,
) -> Result<(), String> {
    let session = read_session(session_dir)?;
    if session.session_id != SESSION || session.spec_id != SPEC || session.status != "complete" {
        return Err("session identity/status differs from frozen protocol".into());
    }
    let record_bytes = bytes(&session_dir.join("records.jsonl"))?;
    let record_hash = sha256_hex(&record_bytes);
    let local = checked_receipt(local_audit_path, session_dir, SESSION)?;
    let independent = checked_receipt(independent_audit_path, session_dir, SESSION)?;
    if local["auditor_env_class_id"] != session.env_class_id
        || independent["auditor_env_class_id"] == session.env_class_id
        || independent["auditor_env_class_id"].is_null()
    {
        return Err("independent receipt is not from another host class".into());
    }
    let records = read_records(session_dir)?;
    if records.len() != 384 || records.iter().filter(|r| r.warmup).count() != 64 {
        return Err("expected 64 warm-ups and 320 measured executions".into());
    }
    let expected_arms: BTreeSet<&str> = ARMS.into_iter().collect();
    let mut by_pair = BTreeMap::<(String, u32), BTreeMap<String, &Record>>::new();
    let mut by_arm = BTreeMap::<String, Vec<&Record>>::new();
    for record in &records {
        checked_record(record, SESSION)?;
        if !expected_arms.contains(record.arm.as_str()) {
            return Err(format!("unexpected arm {}", record.arm));
        }
        if !record.warmup {
            by_pair
                .entry((record.workload.workload_id.clone(), record.round))
                .or_default()
                .insert(record.arm.clone(), record);
            by_arm.entry(record.arm.clone()).or_default().push(record);
        }
    }
    if by_pair.len() != 40 || by_arm.len() != 8 {
        return Err("expected eight targets in five complete measured rounds".into());
    }
    for ((workload, round), arms) in &by_pair {
        if arms.len() != 8
            || arms.keys().map(String::as_str).collect::<BTreeSet<_>>() != expected_arms
        {
            return Err(format!("{workload} round {round} lacks an arm"));
        }
        let rho = arms["rho-strong"];
        for record in arms.values() {
            if record.workload.target != rho.workload.target
                || record.algorithm_seed != rho.algorithm_seed
                || record.outcome.recovered != rho.outcome.recovered
            {
                return Err(format!("{workload} round {round} is not a paired solve"));
            }
        }
    }
    let rho_records = &by_arm["rho-strong"];
    let rho_total_gae: f64 = rho_records
        .iter()
        .map(|r| r.cost.total_gae.expect("checked above"))
        .sum();
    let rho_mean_s: f64 = rho_records
        .iter()
        .map(|r| r.cost.s.expect("checked above"))
        .sum::<f64>()
        / 40.0;
    let floor_s = rho_records[0].boundaries.floor_s;
    if rho_records
        .iter()
        .any(|r| r.boundaries.floor_s != floor_s || r.boundaries.automorphisms_available != 74)
    {
        return Err("generic floor or automorphism count changed".into());
    }
    let mut candidates = Vec::new();
    for k in KS {
        candidates.push(arm_row(
            session_dir,
            claims_dir,
            &by_arm[&format!("ic-k{k}")],
            k,
            floor_s,
        )?);
    }
    let selected = candidates
        .iter()
        .min_by(|a, b| {
            required_f64(a, "total_cold_gae_lower_bound")
                .unwrap()
                .total_cmp(&required_f64(b, "total_cold_gae_lower_bound").unwrap())
                .then_with(|| {
                    a["folded_columns"]
                        .as_u64()
                        .cmp(&b["folded_columns"].as_u64())
                })
        })
        .ok_or("no candidates")?;
    let selected_k = selected["folded_columns"].as_u64().ok_or("selected K")?;
    let selected_ratio = required_f64(selected, "cold_counted_over_k42_diagnostic")?;
    let upper_ci = selected["cold_counted_over_k42_ci95_diagnostic"][1]
        .as_f64()
        .ok_or("missing paired confidence limit")?;
    let aa = comparison(session_dir, "ic-k42", "ic-k42-control")?;
    if aa["ops"]["ratio_b_over_a"] != 1.0 {
        return Err("K42 A/A control changed counted work".into());
    }
    let decision = if selected_ratio <= 0.8 && upper_ci < 1.0 {
        "COUNTED_ENGINEERING_LEAD"
    } else {
        "NO_COUNTED_COLUMN_LEAD"
    };
    let output = json!({
        "schema": "n37-rank-columns-analysis/v1",
        "status": "independently_replayed_l0_bounded_diagnostic",
        "protocol_source_commit": SOURCE,
        "session_id": SESSION,
        "spec_id": SPEC,
        "measuring_binary_sha256": session.binary_sha256,
        "source_env_class_id": session.env_class_id,
        "independent_env_class_id": independent["auditor_env_class_id"],
        "local_audit_sha256": sha256_hex(&bytes(local_audit_path)?),
        "independent_audit_sha256": sha256_hex(&bytes(independent_audit_path)?),
        "records_sha256": record_hash,
        "target_count": 8,
        "measured_rounds_per_target": 5,
        "all_measured_verified": true,
        "generic_floor_s": floor_s,
        "reference": {
            "arm": "rho-strong",
            "measured_runs": 40,
            "verified_runs": 40,
            "total_cold_gae_lower_bound": rho_total_gae,
            "mean_cold_s_lower_bound": rho_mean_s,
            "s_over_generic_floor_lower_bound": rho_mean_s / floor_s,
        },
        "candidates": candidates,
        "aa_counted_ratio": aa["ops"]["ratio_b_over_a"],
        "selection": {
            "selected_k": selected_k,
            "selected_candidate_id": selected["candidate_id"],
            "selected_cold_counted_over_k42": selected_ratio,
            "selected_cold_counted_over_k42_ci95": selected["cold_counted_over_k42_ci95_diagnostic"],
            "selected_cold_counted_over_rho": selected["cold_counted_over_rho_diagnostic"],
            "decision": decision,
        },
        "online_speedup": Value::Null,
        "fully_priced_cold_speedup": Value::Null,
        "limits": [
            "all CPU wall intervals are L0 exploratory",
            "field arithmetic, hashing, allocation, and modular combination remain unpriced",
            "a ratio of two operation lower bounds does not bound the true speedup",
            "n37 is not an n41/n53/n83/n131 scaling result",
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
    let usage = "usage: n37_rank_columns_analyze SESSION CLAIMS LOCAL_AUDIT INDEPENDENT_AUDIT OUT";
    let session = args.next().ok_or(usage)?;
    let claims = args.next().ok_or(usage)?;
    let local = args.next().ok_or(usage)?;
    let independent = args.next().ok_or(usage)?;
    let out = args.next().ok_or(usage)?;
    if args.next().is_some() {
        return Err(usage.into());
    }
    run(
        Path::new(&session),
        Path::new(&claims),
        Path::new(&local),
        Path::new(&independent),
        Path::new(&out),
    )
}
