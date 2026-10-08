//! Apply the preregistered n37 complete-table setup-floor decision to sealed
//! ecbench records. This deliberately does not compute an IC/rho speedup.

use crypto_lib::cryptanalysis::ecbench::workload::PUBLIC_TARGET_LAW;
use crypto_lib::cryptanalysis::ecbench::{canonical::sha256_hex, record, runner};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;

const R: u64 = 230_603_167;
const PAIRS: u64 = 4_831_386;
const SPEC_ID: &str = "ECS1h999d25a394f1";
const SPEC_SHA: &str = "999d25a394f18734933455821ec41d1e699c8d1b69e890a83473bd33bcc3820b";
const CURVE: &str = "icv1-f2m37-tm534059-32aad96b";
const METHOD: &str = "rho.signed_frobenius_strong";
const SUPPORT_SHA: &str = "8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5";

fn read_json(path: &Path) -> Result<Value, String> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|e| format!("{} JSON: {e}", path.display()))
}

fn analyze(
    count_a: &Path,
    count_b: &Path,
    session_dir: &Path,
    audit_path: &Path,
) -> Result<Value, String> {
    let a_bytes = fs::read(count_a).map_err(|e| format!("{}: {e}", count_a.display()))?;
    let b_bytes = fs::read(count_b).map_err(|e| format!("{}: {e}", count_b.display()))?;
    if a_bytes != b_bytes {
        return Err("the two independent count receipts differ".into());
    }
    let count: Value = serde_json::from_slice(&a_bytes).map_err(|e| e.to_string())?;
    if count["schema"] != "n37-complete-table-floor/v1"
        || count["status"] != "PASS"
        || count["support_gzip_sha256"] != SUPPORT_SHA
        || count["r"].as_u64() != Some(R)
        || count["pairs_each"].as_u64() != Some(PAIRS)
        || count["cold_setup_gae_floor_each"].as_u64() != Some(PAIRS)
        || count["policies"].as_array().map(Vec::len) != Some(4)
        || count["policies"].as_array().is_none_or(|p| {
            p.iter()
                .any(|x| x["pair_gae"].as_f64() != Some(PAIRS as f64))
        })
    {
        return Err("count receipt does not certify four complete tables".into());
    }

    let session = runner::read_session(session_dir)?;
    if session.status != "complete"
        || session.spec_id != SPEC_ID
        || session.spec_sha256 != SPEC_SHA
        || session.records_written != 24
        || session.status_counts.get("verified") != Some(&24)
    {
        return Err("rho session does not match the frozen spec or completion gate".into());
    }
    let raw_records =
        fs::read(session_dir.join("records.jsonl")).map_err(|e| format!("records.jsonl: {e}"))?;
    let records_sha = sha256_hex(&raw_records);
    if session.records_sha256.as_deref() != Some(records_sha.as_str()) {
        return Err("rho records hash differs from sealed session".into());
    }
    // This receipt was migrated after serde_json's exact float parser was
    // enabled: preserve the value already written in the sealed source line.
    let lines = runner::read_record_lines_exact(session_dir)?;
    if lines.len() != 24 {
        return Err("rho session does not contain 24 warmup plus measured records".into());
    }

    let audit_bytes = fs::read(audit_path).map_err(|e| format!("audit: {e}"))?;
    let audit: Value = read_json(audit_path)?;
    if audit["schema"] != "ecbench.audit/v1"
        || audit["ok"] != true
        || audit["session_id"] != session.session_id
        || audit["verified_records"].as_u64() != Some(24)
        || audit["replays"].as_array().map(Vec::len) != Some(16)
        || audit["replays"].as_array().is_none_or(|replays| {
            replays
                .iter()
                .any(|r| r["reproduced"] != true || r["detail"] != "identical")
        })
        || audit["files"]["records.jsonl"].as_str() != Some(records_sha.as_str())
    {
        return Err("full ecbench audit and replay gate failed".into());
    }

    let mut rows = Vec::with_capacity(16);
    let mut workloads = BTreeSet::new();
    let mut targets = BTreeSet::new();
    let mut per_workload = BTreeMap::<String, u32>::new();
    let mut measured_failures = 0u64;
    let mut max_rho_gae = 0.0f64;
    let mut min_floor_ratio = f64::INFINITY;
    let mut rho_gae_sum = 0.0f64;
    for (line, rec) in &lines {
        if !record::check_line_seal(line, &rec.record_id)
            || rec.session_id != session.session_id
            || rec.arm != "rho-strong"
            || rec.method.id != METHOD
            || rec.workload.curve.slug != CURVE
            || rec.workload.curve.r != u128::from(R)
            || rec.workload.curve.group_order != 137_439_487_532
            || rec.workload.curve.generator[0] != "0x17262ad4f3"
            || rec.workload.curve.generator[1] != "0x82fe10673"
            || !rec.workload.curve.registered
            || rec.workload.target_law != PUBLIC_TARGET_LAW
            || rec.workload.planted.is_some()
            || rec.workload.target_seed != 202610030037
            || rec.workload.target_index >= 8
        {
            return Err(format!(
                "record {} differs from frozen target/method",
                rec.seq
            ));
        }
        if rec.warmup {
            continue;
        }
        workloads.insert(rec.workload.workload_id.clone());
        targets.insert(rec.workload.target.clone());
        *per_workload
            .entry(rec.workload.workload_id.clone())
            .or_default() += 1;
        let verified = rec.outcome.status == "verified"
            && rec.outcome.matches_target == Some(true)
            && rec.cost.unit == record::UNIT;
        let gae = if verified { rec.cost.total_gae } else { None };
        let ratio = gae.filter(|x| x.is_finite() && *x > 0.0).map(|x| {
            max_rho_gae = max_rho_gae.max(x);
            rho_gae_sum += x;
            let ratio = PAIRS as f64 / x;
            min_floor_ratio = min_floor_ratio.min(ratio);
            ratio
        });
        if ratio.is_none() {
            measured_failures += 1;
        }
        rows.push(json!({
            "run_id": rec.run_id,
            "record_id": rec.record_id,
            "workload_id": rec.workload.workload_id,
            "target_index": rec.workload.target_index,
            "target": rec.workload.target,
            "round": rec.round,
            "status": rec.outcome.status,
            "matches_target": rec.outcome.matches_target,
            "rho_total_gae": gae,
            "rho_s": rec.cost.s,
            "table_setup_floor_over_rho_gae": ratio,
            "rho_unpriced": rec.cost.unpriced,
            "isolation": rec.isolation.level,
        }));
    }
    if rows.len() != 16
        || workloads.len() != 8
        || targets.len() != 8
        || per_workload.values().any(|&count| count != 2)
    {
        return Err("expected two measured rounds on eight distinct workloads".into());
    }
    let decision = if measured_failures == 0 && min_floor_ratio >= 100.0 {
        "REJECT_COMPLETE_TABLE_COLD_N37"
    } else {
        "INCONCLUSIVE"
    };
    Ok(json!({
        "schema": "n37-complete-table-floor-decision/v1",
        "decision": decision,
        "scope": "charged-GAE necessary-cost gate for the fixed complete m<=3 pair-table oracle at n37; no runtime speedup, full IC cost, scaling claim, or general IC no-go",
        "spec_id": SPEC_ID,
        "session_id": session.session_id,
        "session_records_sha256": records_sha,
        "audit_sha256": sha256_hex(&audit_bytes),
        "count_receipt_sha256": sha256_hex(&a_bytes),
        "public_workloads": workloads.len(),
        "measured_runs": rows.len(),
        "measured_failures": measured_failures,
        "table_setup_floor_gae_each": PAIRS,
        "table_setup_floor_s_each": PAIRS as f64 / (R as f64).sqrt(),
        "rho_gae_sum": rho_gae_sum,
        "rho_gae_max": max_rho_gae,
        "minimum_table_floor_over_rho_gae": if measured_failures == 0 {Some(min_floor_ratio)} else {None},
        "preregistered_threshold": 100.0,
        "rows": rows,
    }))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert!(
        args.len() == 5,
        "usage: n37_complete_table_floor_analyze COUNT_A COUNT_B SESSION_DIR AUDIT OUTPUT"
    );
    let result = analyze(
        Path::new(&args[0]),
        Path::new(&args[1]),
        Path::new(&args[2]),
        Path::new(&args[3]),
    )
    .expect("analyze preregistered n37 floor gate");
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&args[4])
        .expect("refuse to overwrite decision");
    serde_json::to_writer_pretty(&mut file, &result).expect("write decision");
    file.write_all(b"\n").expect("newline");
    println!("{}", args[4]);
}
