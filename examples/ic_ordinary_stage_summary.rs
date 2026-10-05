//! Read-only paired analysis of two independently audited, target-free n17 panels.
//! This report is stage evidence only; it never computes an IC/rho speedup.
use crypto_lib::cryptanalysis::{
    ecbench::stats::{bootstrap_ci, mean},
    prepared_sat_control::sha256,
};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap, env, fs, fs::OpenOptions, io::Write, path::Path, process::ExitCode,
};

const QUERIES: usize = 512;
const BOOTSTRAPS: usize = 10_000;
const BOOTSTRAP_SEED: u64 = 20261005;
const STATUS: &str = "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT";

#[derive(Clone)]
struct Pair {
    f5_witness: bool,
    cms_witness: bool,
    cost_delta_ms: f64,
}

fn require(ok: bool, why: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(why.into())
    }
}

fn read(path: &Path, limit: u64) -> Result<(Vec<u8>, Value), String> {
    require(path.is_absolute(), "input path is not absolute")?;
    let meta = fs::symlink_metadata(path).map_err(|e| e.to_string())?;
    require(
        meta.is_file() && !meta.file_type().is_symlink() && meta.len() <= limit,
        "input is not a bounded regular file",
    )?;
    let bytes = fs::read(path).map_err(|e| e.to_string())?;
    require(bytes.len() as u64 <= limit, "input grew beyond limit")?;
    let value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    Ok((bytes, value))
}

fn number(value: &Value, what: &str) -> Result<u64, String> {
    value.as_u64().ok_or_else(|| format!("missing {what}"))
}

fn wilson95(successes: usize, count: usize) -> [f64; 2] {
    let z = 1.959963984540054_f64;
    let n = count as f64;
    let p = successes as f64 / n;
    let denom = 1.0 + z * z / n;
    let center = (p + z * z / (2.0 * n)) / denom;
    let radius = z * ((p * (1.0 - p) + z * z / (4.0 * n)) / n).sqrt() / denom;
    [center - radius, center + radius]
}

fn phase_ms(audit: &Value) -> Result<Value, String> {
    let phases = audit["exclusive_preparation_phases_ns"]
        .as_object()
        .ok_or("missing exclusive preparation phases")?;
    let mut output = serde_json::Map::new();
    let mut sum = 0_u64;
    for (name, value) in phases {
        let ns = number(value, name)?;
        sum = sum.checked_add(ns).ok_or("phase sum overflow")?;
        output.insert(name.clone(), json!(ns as f64 / 1e6));
    }
    require(
        sum == number(&audit["preparation_wall_ns"], "preparation wall")?,
        "preparation phases do not close",
    )?;
    Ok(Value::Object(output))
}

fn analyze(
    f5_audit_path: &Path,
    cms_audit_path: &Path,
    f5_producer_path: &Path,
    cms_producer_path: &Path,
) -> Result<Value, String> {
    let (f5_audit_bytes, f5_audit) = read(f5_audit_path, 2 * 1024 * 1024)?;
    let (cms_audit_bytes, cms_audit) = read(cms_audit_path, 2 * 1024 * 1024)?;
    let (f5_producer_bytes, f5_producer) = read(f5_producer_path, 16 * 1024 * 1024)?;
    let (cms_producer_bytes, cms_producer) = read(cms_producer_path, 16 * 1024 * 1024)?;
    for (name, audit, producer, bytes) in [
        ("matrix_f5", &f5_audit, &f5_producer, &f5_producer_bytes),
        (
            "cryptominisat",
            &cms_audit,
            &cms_producer,
            &cms_producer_bytes,
        ),
    ] {
        require(
            audit["status"] == STATUS
                && audit["source_bound_execution_admitted"] == true
                && audit["audited_ordinary_queries"] == QUERIES
                && audit["mathematics"]["declared_family"] == name
                && audit["mathematics"]["query_law"] == "independent-probe-scalar-rand08-v1"
                && audit["mathematics"]["panel_complete"] == true
                && audit["producer_sha256"] == sha256(bytes)
                && producer["panel_complete"] == true
                && producer["plan"]["family"] == name
                && producer["plan"]["algorithm_seed"] == 2026100311_u64
                && producer["plan"]["planned_queries"] == QUERIES,
            "original audit or paired producer identity differs",
        )?;
    }
    let f5_rows = f5_producer["query_records"]
        .as_array()
        .filter(|rows| rows.len() == QUERIES)
        .ok_or("F5 attempt count differs")?;
    let cms_rows = cms_producer["query_records"]
        .as_array()
        .filter(|rows| rows.len() == QUERIES)
        .ok_or("CMS attempt count differs")?;
    let mut pairs = Vec::with_capacity(QUERIES);
    let mut attempts = Vec::with_capacity(QUERIES);
    let mut f5_mix = BTreeMap::<String, usize>::new();
    let mut cms_mix = BTreeMap::<String, usize>::new();
    let mut cms_native_mix = BTreeMap::<String, usize>::new();
    let mut both = 0;
    let mut f5_only = 0;
    let mut cms_only = 0;
    let mut neither = 0;
    for trial in 0..QUERIES {
        let f5 = &f5_rows[trial];
        let cms = &cms_rows[trial];
        let scalar = number(&f5["attempt"]["scalar"], "F5 scalar")?;
        require(
            f5["attempt"]["trial"] == trial
                && cms["attempt"]["trial"] == trial
                && cms["attempt"]["scalar"] == scalar,
            "paired trial or public-point fixture differs",
        )?;
        let f5_outcome = f5["attempt"]["outcome"]
            .as_str()
            .ok_or("missing F5 outcome")?;
        let cms_outcome = cms["attempt"]["outcome"]
            .as_str()
            .ok_or("missing CMS outcome")?;
        let cms_native_status = cms["evidence"]["frontend"]["native_status"]
            .as_str()
            .ok_or("missing CMS native status")?;
        *f5_mix.entry(f5_outcome.to_string()).or_default() += 1;
        *cms_mix.entry(cms_outcome.to_string()).or_default() += 1;
        *cms_native_mix
            .entry(cms_native_status.to_string())
            .or_default() += 1;
        let f5_witness = f5_outcome == "witness";
        let cms_witness = cms_outcome == "witness";
        match (f5_witness, cms_witness) {
            (true, true) => both += 1,
            (true, false) => f5_only += 1,
            (false, true) => cms_only += 1,
            (false, false) => neither += 1,
        }
        let f5_ns = number(&f5["costs"]["observed_wall_ns"], "F5 query wall")?;
        let cms_ns = number(&cms["costs"]["observed_wall_ns"], "CMS query wall")?;
        require(f5_ns > 0 && cms_ns > 0, "nonpositive attempt wall")?;
        pairs.push(Pair {
            f5_witness,
            cms_witness,
            cost_delta_ms: cms_ns as f64 / 1e6 - f5_ns as f64 / 1e6,
        });
        attempts.push(json!({
            "trial":trial,"fixture_scalar":scalar,
            "f5_outcome":f5_outcome,"cms_outcome":cms_outcome,
            "cms_native_status":cms_native_status,
            "f5_observed_wall_ns":f5_ns,"cms_observed_wall_ns":cms_ns
        }));
    }
    require(
        both + f5_only + cms_only + neither == QUERIES,
        "paired contingency does not close",
    )?;
    let f5_witnesses = both + f5_only;
    let cms_witnesses = both + cms_only;
    require(
        f5_audit["mathematics"]["verified_witnesses"] == f5_witnesses
            && cms_audit["mathematics"]["verified_witnesses"] == cms_witnesses,
        "paired witness counts differ from independent audits",
    )?;
    let strata = vec![pairs];
    let witness_delta = (f5_witnesses as f64 - cms_witnesses as f64) / QUERIES as f64;
    let witness_delta_ci = bootstrap_ci(&strata, BOOTSTRAPS, BOOTSTRAP_SEED, |draw| {
        mean(
            &draw[0]
                .iter()
                .map(|p| (p.f5_witness as u8 as f64) - (p.cms_witness as u8 as f64))
                .collect::<Vec<_>>(),
        )
    })
    .ok_or("paired witness bootstrap failed")?;
    let cost_delta_ms = mean(
        &strata[0]
            .iter()
            .map(|pair| pair.cost_delta_ms)
            .collect::<Vec<_>>(),
    )
    .ok_or("paired cost mean failed")?;
    let cost_delta_ci = bootstrap_ci(&strata, BOOTSTRAPS, BOOTSTRAP_SEED ^ 0xC057, |draw| {
        mean(&draw[0].iter().map(|p| p.cost_delta_ms).collect::<Vec<_>>())
    })
    .ok_or("paired cost bootstrap failed")?;
    let f5_accepted = number(&f5_audit["mathematics"]["accepted_rows"], "F5 rows")?;
    let cms_accepted = number(&cms_audit["mathematics"]["accepted_rows"], "CMS rows")?;
    let f5_wall = number(&f5_audit["preparation_wall_ns"], "F5 wall")?;
    let cms_wall = number(&cms_audit["preparation_wall_ns"], "CMS wall")?;
    Ok(json!({
        "schema_version":1,
        "question":"paired-source-bound-n17-ordinary-stage-only-v1",
        "stage_only":true,"candidate_id":null,"workload_id":null,"run_id":null,
        "curve_id":"EC1N17Ce1hdfbf24105ef5",
        "ordinary_queries":QUERIES,"target_count":0,
        "f5_audit_sha256":sha256(&f5_audit_bytes),
        "cms_audit_sha256":sha256(&cms_audit_bytes),
        "f5_producer_sha256":sha256(&f5_producer_bytes),
        "cms_producer_sha256":sha256(&cms_producer_bytes),
        "f5_registration_sha256":f5_audit["registration_sha256"],
        "cms_registration_sha256":cms_audit["registration_sha256"],
        "f5":{
            "family":"matrix_f5","actual_usable_base_points":62,"folded_columns":29,
            "verified_witnesses":f5_witnesses,"witness_rate":f5_witnesses as f64 / QUERIES as f64,
            "witness_rate_wilson95":wilson95(f5_witnesses, QUERIES),
            "accepted_rows":f5_accepted,"final_rank":f5_audit["mathematics"]["rank"],
            "independently_recovered_logs":f5_audit["mathematics"]["independently_recovered_column_logs"].as_array().map(Vec::len),
            "outcome_mix":f5_mix,"preparation_wall_ns":f5_wall,
            "preparation_wall_seconds_per_accepted_row":(f5_accepted > 0).then(|| f5_wall as f64 / f5_accepted as f64 / 1e9),
            "exclusive_preparation_phases_ms":phase_ms(&f5_audit)?,
            "peak_memory":f5_audit["peak_memory"]
        },
        "cms":{
            "family":"cryptominisat","actual_usable_base_points":62,"folded_columns":29,
            "verified_witnesses":cms_witnesses,"witness_rate":cms_witnesses as f64 / QUERIES as f64,
            "witness_rate_wilson95":wilson95(cms_witnesses, QUERIES),
            "accepted_rows":cms_accepted,"final_rank":cms_audit["mathematics"]["rank"],
            "independently_recovered_logs":cms_audit["mathematics"]["independently_recovered_column_logs"].as_array().map(Vec::len),
            "feasible_inconclusive_queries":cms_audit["mathematics"]["feasible_inconclusive_queries"],
            "outcome_mix":cms_mix,"native_status_mix":cms_native_mix,
            "preparation_wall_ns":cms_wall,
            "preparation_wall_seconds_per_accepted_row":(cms_accepted > 0).then(|| cms_wall as f64 / cms_accepted as f64 / 1e9),
            "exclusive_preparation_phases_ms":phase_ms(&cms_audit)?,
            "peak_memory":cms_audit["peak_memory"]
        },
        "paired_witness_contingency":{"both":both,"f5_only":f5_only,"cms_only":cms_only,"neither":neither},
        "paired_f5_minus_cms_witness_rate":witness_delta,
        "paired_f5_minus_cms_witness_rate_bootstrap95":witness_delta_ci,
        "paired_cms_minus_f5_attempt_wall_mean_ms":cost_delta_ms,
        "paired_cms_minus_f5_attempt_wall_mean_bootstrap95_ms":cost_delta_ci,
        "bootstrap_resamples":BOOTSTRAPS,"bootstrap_seed":BOOTSTRAP_SEED,
        "wall_timing_class":"ordinary-host-exploratory",
        "cpu_speedup_claim":null,"online_speedup":null,
        "complete_single_target_result":false,
        "attempts":attempts
    }))
}

fn main() -> ExitCode {
    let args = env::args_os().collect::<Vec<_>>();
    if args.len() != 6 {
        eprintln!("usage: ic_ordinary_stage_summary F5_AUDIT CMS_AUDIT F5_PRODUCER CMS_PRODUCER NEW_OUTPUT_JSON");
        return ExitCode::FAILURE;
    }
    let out = Path::new(&args[5]);
    if !out.is_absolute() || out.exists() {
        eprintln!("new absolute output JSON required");
        return ExitCode::FAILURE;
    }
    let result = analyze(
        Path::new(&args[1]),
        Path::new(&args[2]),
        Path::new(&args[3]),
        Path::new(&args[4]),
    );
    let value = match result {
        Ok(value) => value,
        Err(error) => {
            eprintln!("{error}");
            return ExitCode::FAILURE;
        }
    };
    let mut file = match OpenOptions::new().write(true).create_new(true).open(out) {
        Ok(file) => file,
        Err(error) => {
            eprintln!("{error}");
            return ExitCode::FAILURE;
        }
    };
    let bytes = serde_json::to_vec_pretty(&value).expect("stage summary JSON");
    if let Err(error) = file.write_all(&bytes).and_then(|_| file.sync_all()) {
        eprintln!("{error}");
        return ExitCode::FAILURE;
    }
    println!("{}", out.display());
    ExitCode::SUCCESS
}
