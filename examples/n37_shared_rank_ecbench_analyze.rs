//! Audit-ready one-target online summary from sealed native ecbench records.
//! Cold S remains the harness's comparison; this file never promotes L0 time.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::claim::id_sha256;
use crypto_lib::cryptanalysis::ecbench::record::Record;
use crypto_lib::cryptanalysis::ecbench::runner::{read_records, read_session};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

fn median(mut values: Vec<u64>) -> u64 {
    values.sort_unstable();
    values[values.len() / 2]
}

fn range(values: &[u64]) -> Value {
    json!({
        "min_ns": values.iter().min().copied(),
        "median_ns": median(values.to_vec()),
        "max_ns": values.iter().max().copied(),
    })
}

fn checked_online(record: &Record) -> Result<u64, String> {
    let online = record
        .online
        .as_ref()
        .ok_or_else(|| format!("{} has no online interval", record.run_id))?;
    if online.phases_ns.values().sum::<u64>() != online.wall_ns {
        return Err(format!("{} has nonexclusive online phases", record.run_id));
    }
    if record.outcome.status != "verified" || record.outcome.matches_target != Some(true) {
        return Err(format!("{} did not verify", record.run_id));
    }
    Ok(online.wall_ns)
}

fn run(session_dir: &Path, claim_path: &Path, out_path: &Path) -> Result<(), String> {
    let session = read_session(session_dir)?;
    if session.status != "complete" {
        return Err("session is incomplete".into());
    }
    let claim_bytes = fs::read(claim_path).map_err(|e| e.to_string())?;
    let claim: Value = serde_json::from_slice(&claim_bytes).map_err(|e| e.to_string())?;
    if claim["ecbench_records"]["session"] != session.session_id {
        return Err("claim names a different session".into());
    }
    for name in ["candidate", "workload"] {
        let record_key = format!("{name}_manifest");
        let digest_key = format!("{name}_manifest_sha256");
        let record = &claim[record_key.as_str()];
        let actual = id_sha256(record)?;
        if claim[digest_key.as_str()].as_str() != Some(actual.as_str()) {
            return Err(format!(
                "{name} manifest hash does not match its canonical record"
            ));
        }
    }
    let records = read_records(session_dir)?;
    let mut by_round: BTreeMap<u32, BTreeMap<String, &Record>> = BTreeMap::new();
    for record in records.iter().filter(|r| !r.warmup) {
        by_round
            .entry(record.round)
            .or_default()
            .insert(record.arm.clone(), record);
    }
    if by_round.len() != 5 || records.len() != 18 {
        return Err("expected one warm-up and five complete three-arm rounds".into());
    }
    let mut ic_online = Vec::new();
    let mut rho_online = Vec::new();
    let mut control_online = Vec::new();
    let mut paired_ratios = Vec::new();
    let mut aa_ratios = Vec::new();
    let mut rows = Vec::new();
    let mut phase_samples: BTreeMap<String, Vec<u64>> = BTreeMap::new();
    let mut ic_s = Vec::new();
    let mut rho_s = Vec::new();
    let mut online_gae = Vec::new();
    let mut setup_gae = Vec::new();
    let mut rank_trials = Vec::new();
    let mut attempts = Vec::new();
    for (round, arms) in by_round {
        let get = |arm: &str| {
            arms.get(arm)
                .copied()
                .ok_or_else(|| format!("round {round} lacks {arm}"))
        };
        let ic = get("ic-shared")?;
        let rho = get("rho-strong")?;
        let control = get("ic-control")?;
        if ic.workload.workload_id != rho.workload.workload_id
            || ic.workload.workload_id != control.workload.workload_id
            || ic.workload.target != rho.workload.target
            || ic.workload.target != control.workload.target
            || ic.algorithm_seed != rho.algorithm_seed
            || ic.algorithm_seed != control.algorithm_seed
            || ic.outcome.recovered != rho.outcome.recovered
            || ic.outcome.recovered != control.outcome.recovered
        {
            return Err(format!(
                "round {round} is not same-point, same-seed recovery"
            ));
        }
        let (i, r, c) = (
            checked_online(ic)?,
            checked_online(rho)?,
            checked_online(control)?,
        );
        if i == 0 || r == 0 || c == 0 {
            return Err(format!("round {round} has a zero online interval"));
        }
        let rank = &ic.detail["rank"];
        if rank["base_sha256"] != "8460ac4c28515db701c3897a03b4ce0f28abf7cd56fd759ad095f98436a76dcf"
            || rank["rank"] != 42
            || rank["verified"] != true
        {
            return Err(format!("round {round} rank/base differs from frozen setup"));
        }
        let n_attempts = ic.detail["targets"][0]["attempts"]
            .as_array()
            .ok_or("target attempts missing")?
            .len();
        attempts.push(n_attempts);
        rank_trials.push(rank["trials"].as_u64().ok_or("rank trials missing")?);
        let online_cost: f64 = ic
            .phases
            .iter()
            .filter(|p| p.name.starts_with("target_"))
            .map(|p| p.gae)
            .sum();
        let cold = ic.cost.total_gae.ok_or("IC cold GAE missing")?;
        online_gae.push(online_cost);
        setup_gae.push(cold - online_cost);
        ic_s.push(ic.cost.s.ok_or("IC S missing")?);
        rho_s.push(rho.cost.s.ok_or("rho S missing")?);
        for (name, ns) in &ic.online.as_ref().unwrap().phases_ns {
            phase_samples.entry(name.clone()).or_default().push(*ns);
        }
        ic_online.push(i);
        rho_online.push(r);
        control_online.push(c);
        paired_ratios.push(r as f64 / i as f64);
        aa_ratios.push(c as f64 / i as f64);
        rows.push(json!({
            "round": round,
            "ic_record": ic.record_id,
            "rho_record": rho.record_id,
            "control_record": control.record_id,
            "ic_online_ns": i,
            "rho_online_ns": r,
            "control_online_ns": c,
            "rho_over_ic_online_descriptive": r as f64 / i as f64,
            "control_over_ic_online_descriptive": c as f64 / i as f64,
            "ic_target_attempts": n_attempts,
            "rank_trials": rank["trials"],
        }));
    }
    paired_ratios.sort_by(f64::total_cmp);
    aa_ratios.sort_by(f64::total_cmp);
    let phases: BTreeMap<String, Value> = phase_samples
        .iter()
        .map(|(name, values)| (name.clone(), range(values)))
        .collect();
    let comparison_path = session_dir.join("comparisons/rho-strong__ic-shared.json");
    let comparison_bytes = fs::read(&comparison_path).map_err(|e| e.to_string())?;
    let comparison: Value = serde_json::from_slice(&comparison_bytes).map_err(|e| e.to_string())?;
    if comparison["ops"]["bounded"] != true || comparison["ops"]["pairs"] != 5 {
        return Err("cold comparison is not the expected bounded five-pair row".into());
    }
    let output = json!({
        "schema": "n37-shared-rank-ecbench-analysis/v1",
        "status": "verified_bounded_l0_diagnostic",
        "session_id": session.session_id,
        "candidate_id": claim["candidate_id"],
        "candidate_manifest_sha256": claim["candidate_manifest_sha256"],
        "workload_id": claim["workload_id"],
        "ecbench_workload_id": records[0].workload.workload_id,
        "public_target": records[0].workload.target,
        "target_count": 1,
        "measured_rounds": 5,
        "all_measured_verified": true,
        "ic_rank_trials": rank_trials,
        "ic_target_attempts": attempts,
        "ic_online_ns": range(&ic_online),
        "rho_online_ns": range(&rho_online),
        "control_online_ns": range(&control_online),
        "ic_online_phase_ns": phases,
        "rho_over_ic_online_ratio_descriptive_median": paired_ratios[2],
        "rho_over_ic_online_ratio_descriptive_range": [paired_ratios[0], paired_ratios[4]],
        "control_over_ic_online_ratio_descriptive_median": aa_ratios[2],
        "control_over_ic_online_ratio_descriptive_range": [aa_ratios[0], aa_ratios[4]],
        "ic_online_gae_lower_bound": online_gae,
        "ic_setup_gae_lower_bound": setup_gae,
        "ic_cold_s_lower_bound": ic_s,
        "rho_cold_s_lower_bound": rho_s,
        "cold_counted_ic_over_rho_bounded": comparison["ops"]["ratio_b_over_a"],
        "cold_counted_ic_over_rho_ci95_bounded": comparison["ops"]["ci95"],
        "online_speedup_admissible": false,
        "online_speedup": null,
        "limits": ["L0 macOS: no CPU isolation", "one target: no target-population interval", "native field, hash, allocation and modular work incompletely priced", "other-host replay receipt not yet attached"],
        "session_records_sha256": sha256_hex(&fs::read(session_dir.join("records.jsonl")).map_err(|e| e.to_string())?),
        "comparison_sha256": sha256_hex(&comparison_bytes),
        "diagnostic_claim_sha256": sha256_hex(&claim_bytes),
        "rounds": rows,
    });
    serde_json::to_writer_pretty(
        fs::File::create(out_path).map_err(|e| e.to_string())?,
        &output,
    )
    .map_err(|e| e.to_string())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let session = args.next().ok_or("usage: analyze SESSION CLAIM OUT")?;
    let claim = args.next().ok_or("missing claim")?;
    let out = args.next().ok_or("missing output")?;
    if args.next().is_some() {
        return Err("too many arguments".into());
    }
    run(Path::new(&session), Path::new(&claim), Path::new(&out))
}
