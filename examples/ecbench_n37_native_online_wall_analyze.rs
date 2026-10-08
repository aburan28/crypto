//! Recompute the fresh n37 native online-wall screen from sealed raw records.
//! Usage: ecbench_n37_native_online_wall_analyze RESEARCH_DIR SESSION_DIR RECEIPT.json OUT.json

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crypto_lib::cryptanalysis::ecbench::{audit, canonical, runner};
use serde_json::{json, Value};

const ARMS: [&str; 4] = ["rho-strong", "ic-k8", "ic-k16", "ic-k16-control"];
const BOOTSTRAP_SAMPLES: usize = 20_000;
const BOOTSTRAP_SEED: u64 = 202_610_046_103;
const SOURCE_COMMIT: &str = "9e15353bf83b15e4cb8332519d8515613ac4dfea";
const K8_CANDIDATE: &str =
    "IC1N37Ckb0fb592PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hd7ebd078ae0c";
const K16_CANDIDATE: &str =
    "IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hbbfdf029e5a2";
const K8_BASE_SHA256: &str = "8a0c8a2892c2579002215083812ae57ba536bf8eaf1a243b18a1152ecd317dbf";
const K16_BASE_SHA256: &str = "4dbbb8c436f34319fc5971b31245063bfbbca619a29254e64afc75e924e74b4f";

fn read(path: &Path) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn json_file(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn ratio(a: u128, b: u128) -> Result<f64, String> {
    if b == 0 {
        return Err("zero online denominator".into());
    }
    Ok(a as f64 / b as f64)
}

fn sums(group: &BTreeMap<String, u128>, numerator: &str, denominator: &str) -> Result<f64, String> {
    ratio(group[numerator], group[denominator])
}

/// Resample whole targets, retaining all five paired rounds in each block.
fn paired_interval(
    targets: &[BTreeMap<String, u128>],
    numerator: &str,
    denominator: &str,
) -> Result<Value, String> {
    if targets.len() != 16 {
        return Err("target-block bootstrap needs the frozen 16 targets".into());
    }
    let mut numerator_total = 0u128;
    let mut denominator_total = 0u128;
    for target in targets {
        numerator_total += target[numerator];
        denominator_total += target[denominator];
    }
    let observed = ratio(numerator_total, denominator_total)?;
    let mut rng = crypto_lib::cryptanalysis::ecbench::stats::Resampler::new(BOOTSTRAP_SEED);
    let mut draws = Vec::with_capacity(BOOTSTRAP_SAMPLES);
    for _ in 0..BOOTSTRAP_SAMPLES {
        let mut a = 0u128;
        let mut b = 0u128;
        for _ in 0..targets.len() {
            let target = &targets[rng.below(targets.len())];
            a += target[numerator];
            b += target[denominator];
        }
        draws.push(ratio(a, b)?);
    }
    draws.sort_by(f64::total_cmp);
    Ok(json!({
        "numerator": numerator,
        "denominator": denominator,
        "ratio_of_sums": observed,
        "target_block_bootstrap_95": [draws[499], draws[19_499]],
        "resamples": BOOTSTRAP_SAMPLES,
        "seed": BOOTSTRAP_SEED,
        "unit": "online_wall_ns; descriptive until L2 admission",
    }))
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: ecbench_n37_native_online_wall_analyze RESEARCH_DIR SESSION_DIR RECEIPT.json OUT.json".into());
    }
    let (research, session_dir, receipt_path, output) = (
        PathBuf::from(&args[1]),
        PathBuf::from(&args[2]),
        PathBuf::from(&args[3]),
        PathBuf::from(&args[4]),
    );
    let freeze = json_file(&research.join("FREEZE.json"))?;
    if freeze["source_commit"] != SOURCE_COMMIT
        || freeze["public_target_count"] != 16
        || freeze["primary_target_index"] != 0
        || freeze["planned_executions_including_warmups"] != 384
        || freeze["target_overlap_with_two_prior_n37_k8_k16_panels"] != 0
    {
        return Err("preregistered freeze fields changed".into());
    }
    for (name, key) in [
        ("SPEC.json", "spec_file_sha256"),
        ("PLAN.json", "plan_file_sha256"),
    ] {
        let got = canonical::sha256_hex(&read(&research.join(name))?);
        if freeze[key].as_str() != Some(&got) {
            return Err(format!("{name} differs from the preregistered SHA-256"));
        }
    }
    let plan = json_file(&research.join("PLAN.json"))?;
    let workloads = plan["workloads"]
        .as_array()
        .ok_or("frozen plan has no workloads")?;
    if workloads.len() != 16 || plan["executions"] != 384 {
        return Err("frozen plan has the wrong target or execution count".into());
    }
    let by_index: BTreeMap<u64, (String, Value)> = workloads
        .iter()
        .map(|w| {
            Ok((
                w["target_index"]
                    .as_u64()
                    .ok_or("plan target index missing")?,
                (
                    w["workload_id"]
                        .as_str()
                        .ok_or("plan workload ID missing")?
                        .to_string(),
                    w["target"].clone(),
                ),
            ))
        })
        .collect::<Result<_, String>>()?;
    if by_index.len() != 16
        || by_index.keys().copied().collect::<Vec<_>>() != (0..16).collect::<Vec<_>>()
    {
        return Err("plan target indices are not exactly 0..15".into());
    }
    let session = runner::read_session(&session_dir)?;
    let local = audit::audit(&session_dir, 0)?;
    if !local.ok || session.status != "complete" || session.records_written != 384 {
        return Err(format!(
            "session audit or completeness failed: {}",
            local.problems.join("; ")
        ));
    }
    if session.spec_id != freeze["spec_id"].as_str().unwrap_or("")
        || session.spec_sha256 != freeze["canonical_spec_sha256"].as_str().unwrap_or("")
        || read(&session_dir.join("spec.json"))? != read(&research.join("SPEC.json"))?
    {
        return Err("session used a different frozen spec".into());
    }
    let records_bytes = read(&session_dir.join("records.jsonl"))?;
    let records_sha256 = canonical::sha256_hex(&records_bytes);
    let receipt_bytes = read(&receipt_path)?;
    let receipt_sha256 = canonical::sha256_hex(&receipt_bytes);
    let receipt: audit::AuditReport =
        serde_json::from_slice(&receipt_bytes).map_err(|e| e.to_string())?;
    if !receipt.ok
        || receipt.session_id != session.session_id
        || receipt.files.get("records.jsonl") != Some(&records_sha256)
        || receipt.replays.len() != 320
        || receipt.replays.iter().any(|r| !r.reproduced)
        || receipt.auditor_env_class_id.as_deref() == Some(&session.env_class_id)
        || receipt.auditor_env_class_id.is_none()
    {
        return Err("independent receipt did not replay the exact measured session".into());
    }
    let records = runner::read_records(&session_dir)?;
    if records.len() != 384 {
        return Err("frozen 384-record session is incomplete".into());
    }
    let mut groups: BTreeMap<
        (u64, u32),
        BTreeMap<String, &crypto_lib::cryptanalysis::ecbench::record::Record>,
    > = BTreeMap::new();
    let mut levels = BTreeMap::<String, u64>::new();
    let mut arm_stats = BTreeMap::<String, Value>::new();
    let mut measured_counts = BTreeMap::<String, u64>::new();
    let mut failures = Vec::new();
    let mut all_l2 = true;
    for record in &records {
        if record.session_id != session.session_id
            || record.workload.planted.is_some()
            || record.workload.target_seed != 202610046101
            || record.binary_sha256 != session.binary_sha256
        {
            return Err(format!(
                "record {} differs from the frozen session",
                record.seq
            ));
        }
        let (wid, point) = by_index
            .get(&record.workload.target_index)
            .ok_or("record has an unplanned target index")?;
        if &record.workload.workload_id != wid
            || serde_json::to_value(&record.workload.target).map_err(|e| e.to_string())? != *point
        {
            return Err(format!(
                "record {} has a different public point",
                record.seq
            ));
        }
        if !ARMS.contains(&record.arm.as_str()) {
            return Err(format!("unplanned arm in record {}", record.seq));
        }
        if record.warmup {
            continue;
        }
        *measured_counts.entry(record.arm.clone()).or_default() += 1;
        *levels
            .entry(record.isolation.level.name().into())
            .or_default() += 1;
        all_l2 &= record.isolation.level.name() == "L2";
        if record.outcome.status != "verified" || record.outcome.matches_target != Some(true) {
            failures.push(json!({"seq":record.seq,"arm":record.arm,"workload_id":wid,
                "status":record.outcome.status,"error":record.outcome.error}));
            continue;
        }
        let online = record
            .online
            .as_ref()
            .ok_or("verified record has no online interval")?;
        if online.wall_ns == 0 || online.phases_ns.values().sum::<u64>() != online.wall_ns {
            return Err(format!(
                "record {} has nonexclusive online timing",
                record.seq
            ));
        }
        if record
            .time
            .solve_wall_ns
            .is_none_or(|ns| ns < online.wall_ns)
        {
            return Err(format!(
                "record {} has invalid cold solve timing",
                record.seq
            ));
        }
        if record.arm.starts_with("ic-") {
            let fb = record
                .factor_base
                .as_ref()
                .ok_or("IC record has no factor base")?;
            let (count, columns, sha) = if record.arm == "ic-k8" {
                (592, 8, K8_BASE_SHA256)
            } else {
                (1184, 16, K16_BASE_SHA256)
            };
            if fb.signed_points != count || fb.columns != columns || fb.fb_sha256 != sha {
                return Err(format!(
                    "record {} changed the frozen factor base",
                    record.seq
                ));
            }
        }
        let group = groups
            .entry((record.workload.target_index, record.round))
            .or_default();
        if group.insert(record.arm.clone(), record).is_some() {
            return Err(format!("duplicate target/round/arm at seq {}", record.seq));
        }
    }
    if !failures.is_empty() {
        let partial = json!({"schema":"ecbench.native_online_wall_decision/v1",
            "status":"incomplete", "failures":failures,"online_speedup":null,
            "records_sha256":records_sha256,"independent_receipt_sha256":receipt_sha256});
        std::fs::write(
            &output,
            serde_json::to_vec_pretty(&partial).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        return Ok(());
    }
    if measured_counts.values().copied().collect::<BTreeSet<_>>() != BTreeSet::from([80])
        || measured_counts.len() != 4
        || groups.len() != 80
    {
        return Err("measured arms or paired target/round groups are incomplete".into());
    }
    let mut target_sums = Vec::with_capacity(16);
    let mut per_target = Vec::with_capacity(16);
    let mut aa_max = 0.0f64;
    for index in 0..16u64 {
        let mut totals = BTreeMap::<String, u128>::new();
        let mut rounds = Vec::with_capacity(5);
        let mut workload_id = String::new();
        for round in 1..=5 {
            let group = groups
                .get(&(index, round))
                .ok_or("missing paired target round")?;
            if group.len() != 4 {
                return Err(format!("target {index} round {round} lacks an arm"));
            }
            let first = group[ARMS[0]];
            workload_id = first.workload.workload_id.clone();
            for arm in ARMS {
                let record = group[arm];
                if record.algorithm_seed != first.algorithm_seed {
                    return Err(format!("target {index} round {round} has mismatched seeds"));
                }
                *totals.entry(arm.into()).or_default() +=
                    record.online.as_ref().unwrap().wall_ns as u128;
                let row = arm_stats.entry(arm.into()).or_insert_with(|| {
                    json!({
                        "online_ns_total":0u64,"cold_solve_ns_total":0u64,
                        "max_rss_kib":0u64,"target_attempts_total":0u64,
                        "rank_trials_total":0u64,"rank_hits_total":0u64,
                    })
                });
                row["online_ns_total"] = json!(
                    row["online_ns_total"].as_u64().unwrap()
                        + record.online.as_ref().unwrap().wall_ns
                );
                row["cold_solve_ns_total"] = json!(
                    row["cold_solve_ns_total"].as_u64().unwrap()
                        + record.time.solve_wall_ns.unwrap()
                );
                row["max_rss_kib"] = json!(row["max_rss_kib"]
                    .as_u64()
                    .unwrap()
                    .max(record.time.max_rss_kib));
                for (name, key) in [
                    ("target_attempts_total", "target_attempts"),
                    ("rank_trials_total", "rank_trials"),
                    ("rank_hits_total", "rank_hits"),
                ] {
                    row[name] = json!(
                        row[name].as_u64().unwrap()
                            + record.counters.get(key).copied().unwrap_or(0)
                    );
                }
            }
            let k16 = group["ic-k16"].online.as_ref().unwrap().wall_ns as f64;
            let control = group["ic-k16-control"].online.as_ref().unwrap().wall_ns as f64;
            aa_max = aa_max.max((k16 / control - 1.0).abs());
            rounds.push(json!({"round":round,"algorithm_seed":first.algorithm_seed,
                "online_ns":ARMS.iter().map(|arm| ((*arm).to_string(), group[*arm].online.as_ref().unwrap().wall_ns)).collect::<BTreeMap<_,_>>(),
                "rho_over_k8":ratio(group["rho-strong"].online.as_ref().unwrap().wall_ns as u128,
                    group["ic-k8"].online.as_ref().unwrap().wall_ns as u128)?,
                "rho_over_k16":ratio(group["rho-strong"].online.as_ref().unwrap().wall_ns as u128,
                    group["ic-k16"].online.as_ref().unwrap().wall_ns as u128)?}));
        }
        per_target.push(json!({"target_index":index,"workload_id":workload_id,
            "rounds":rounds,"online_ns_totals":totals,
            "descriptive_rho_over_k8":sums(&totals,"rho-strong","ic-k8")?,
            "descriptive_rho_over_k16":sums(&totals,"rho-strong","ic-k16")?,
            "descriptive_k8_over_k16":sums(&totals,"ic-k8","ic-k16")?}));
        target_sums.push(totals);
    }
    let k8_over_k16 = paired_interval(&target_sums, "ic-k8", "ic-k16")?;
    let k8_ci = k8_over_k16["target_block_bootstrap_95"]
        .as_array()
        .ok_or("K8/K16 interval absent")?;
    let decision = if aa_max >= 0.05 {
        "carry_both_host_noise_exceeds_gate"
    } else if k8_ci[0].as_f64().unwrap() > 1.10 {
        "prioritize_k16_for_L2_wall_gate"
    } else if k8_ci[1].as_f64().unwrap() < 0.90 {
        "prioritize_k8_for_L2_wall_gate"
    } else {
        "carry_both_native_wall_inconclusive"
    };
    let result = json!({
        "schema":"ecbench.native_online_wall_decision/v1",
        "status":"complete_exploratory_hosted",
        "source_commit":SOURCE_COMMIT,
        "session_id":session.session_id,
        "session_binary_sha256":session.binary_sha256,
        "session_env_class_id":session.env_class_id,
        "auditor_env_class_id":receipt.auditor_env_class_id,
        "records_sha256":records_sha256,
        "independent_receipt_sha256":receipt_sha256,
        "measured_verified_runs":320,
        "same_target_paired_rounds":80,
        "isolation_levels":levels,
        "all_measured_L2":all_l2,
        "online_speedup":Value::Null,
        "online_speedup_reason":"hosted VM has no auditable host-wide L2 receipt; ratios are exploratory diagnostics",
        "primary_one_target":per_target[0],
        "secondary_target_rows":&per_target[1..],
        "panel_descriptive_ratios":{
            "rho_over_k8":paired_interval(&target_sums,"rho-strong","ic-k8")?,
            "rho_over_k16":paired_interval(&target_sums,"rho-strong","ic-k16")?,
            "k8_over_k16":k8_over_k16,
            "k16_over_k16_control":paired_interval(&target_sums,"ic-k16","ic-k16-control")?,
        },
        "k16_aa_max_relative_deviation":aa_max,
        "arm_stage_totals":arm_stats,
        "candidate_ids":{"ic-k8":K8_CANDIDATE,"ic-k16":K16_CANDIDATE,
            "ic-k16-control":K16_CANDIDATE},
        "decision":decision,
        "failures":[],
    });
    let mut bytes = serde_json::to_vec_pretty(&result).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    std::fs::write(&output, bytes).map_err(|e| format!("{}: {e}", output.display()))
}

fn main() {
    if let Err(error) = run() {
        eprintln!("ecbench_n37_native_online_wall_analyze: {error}");
        std::process::exit(1);
    }
}
