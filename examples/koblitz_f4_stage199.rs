//! Native composer/verifier for the Stage 199 selected-default replay.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

const REPLAY_SCHEMA: &str = "koblitz_stage199_default_replay.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage199_default_replay_verification.v1";
const FINAL_SCHEMA: &str = "koblitz_stage199_final.v1";
const FINAL_VERIFICATION_SCHEMA: &str = "koblitz_stage199_final_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";
const STAGE198_SHA256: &str = "38dd39b13cbcf26ac0fd6e5f43f27bf6e733488d349cc24c197c4842047c5198";

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Artifact {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct ProcessMetrics {
    schema: String,
    meter: String,
    command: Vec<String>,
    environment: BTreeMap<String, String>,
    pid: u32,
    returncode: i32,
    timed_out: bool,
    watchdog_seconds: u64,
    wall_seconds: f64,
    user_seconds: Option<f64>,
    system_seconds: Option<f64>,
    total_core_seconds: Option<f64>,
    single_core_seconds: Option<f64>,
    peak_rss_bytes: Option<u64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct MeterReceipt {
    schema: String,
    process: ProcessMetrics,
    stdout: Artifact,
    stderr: Artifact,
    metrics: Artifact,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct DefaultReplay {
    schema: String,
    inherited_stage198_result_sha256: String,
    correct: bool,
    errors: Vec<String>,
    receipt: MeterReceipt,
    report: Value,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct Verification {
    schema: String,
    result_sha256: String,
    checks: usize,
    passed: usize,
    failed: Vec<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Charge {
    components: usize,
    wall_seconds_sum: f64,
    total_core_seconds_sum: f64,
    peak_rss_bytes_max: u64,
    single_core_seconds: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct FinalResult {
    schema: String,
    selected_commit: String,
    finalizer_commit: String,
    decision_status: String,
    runtime_default: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    single_core_seconds: Option<f64>,
    default_replay: DefaultReplay,
    measured_stage_charge: Charge,
    cumulative_measured_lower_bound: Charge,
    complete_campaign_cost: Option<Charge>,
    artifacts: BTreeMap<String, Artifact>,
    all_seven_gates_passed: bool,
    koblitz_index_calculus_sota: bool,
}

fn main() {
    if let Err(error) = real_main() {
        eprintln!("{error}");
        std::process::exit(2);
    }
}

fn real_main() -> AnyResult<()> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    match args.first().map(String::as_str) {
        Some("compose") if args.len() == 4 => {
            let result = compose_replay(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &result)?;
            let verification = verify_replay(Path::new(&args[2]))?;
            write_json(Path::new(&args[3]), &verification)?;
            println!(
                "{}",
                serde_json::to_string(&json!({
                    "schema":"koblitz_stage199_terminal.v1",
                    "correct":result.correct,
                    "verification_checks":verification.checks,
                    "verification_passed":verification.passed,
                    "verification_failed":verification.failed,
                }))
                .map_err(|e| e.to_string())?
            );
            if result.correct && verification.failed.is_empty() {
                Ok(())
            } else {
                Err("default replay failed".into())
            }
        }
        Some("verify") if args.len() == 3 => {
            let verification = verify_replay(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("default replay verification failed".into())
            }
        }
        Some("finalize") if args.len() == 5 => finalize(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &args[3],
            &args[4],
        ),
        Some("verify-final") if args.len() == 3 => {
            let verification = verify_final(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("final verification failed".into())
            }
        }
        _ => Err("usage: koblitz_f4_stage199 compose RUN_DIR RESULT VERIFICATION\n       koblitz_f4_stage199 verify RESULT VERIFICATION\n       koblitz_f4_stage199 finalize STAGE RESULT SELECTED_COMMIT FINALIZER_COMMIT\n       koblitz_f4_stage199 verify-final RESULT VERIFICATION".into()),
    }
}

fn compose_replay(directory: &Path) -> AnyResult<DefaultReplay> {
    let receipt: MeterReceipt = read_json(&directory.join("receipt.json"))?;
    let mut errors = validate_receipt(directory, &receipt);
    errors.extend(validate_command(&receipt.process));
    let report = match read_json::<Value>(&directory.join(&receipt.stdout.path)) {
        Ok(report) => {
            errors.extend(validate_report(&report));
            report
        }
        Err(error) => {
            errors.push(error);
            Value::Null
        }
    };
    let stage = find_stage_root(directory)?;
    if let Err(error) = verify_stage198_parent(stage) {
        errors.push(error);
    }
    Ok(DefaultReplay {
        schema: REPLAY_SCHEMA.into(),
        inherited_stage198_result_sha256: STAGE198_SHA256.into(),
        correct: errors.is_empty(),
        errors,
        receipt,
        report,
    })
}

fn validate_receipt(directory: &Path, receipt: &MeterReceipt) -> Vec<String> {
    let mut failures = Vec::new();
    if receipt.schema != RECEIPT_SCHEMA || receipt.process.schema != METER_SCHEMA {
        failures.push("meter schema".into());
    }
    if receipt.process.meter != "fresh child process wait4/getrusage; all child threads included"
        || !receipt.process.environment.is_empty()
        || receipt.process.returncode != 0
        || receipt.process.timed_out
        || receipt.process.watchdog_seconds != 360
        || receipt.process.single_core_seconds.is_some()
        || receipt.process.total_core_seconds.is_none()
        || receipt.process.peak_rss_bytes.is_none()
    {
        failures.push("process accounting contract".into());
    }
    if receipt.stdout.path != "stdout.txt"
        || receipt.stderr.path != "stderr.txt"
        || receipt.metrics.path != "metrics.json"
    {
        failures.push("receipt paths".into());
    }
    for (name, artifact) in [
        ("stdout", &receipt.stdout),
        ("stderr", &receipt.stderr),
        ("metrics", &receipt.metrics),
    ] {
        match artifact_at(directory, artifact) {
            Ok(actual) if actual == *artifact => {}
            Ok(_) => failures.push(format!("{name} artifact mismatch")),
            Err(error) => failures.push(format!("{name}: {error}")),
        }
    }
    match read_json::<ProcessMetrics>(&directory.join("metrics.json")) {
        Ok(metrics) if metrics == receipt.process => {}
        Ok(_) => failures.push("embedded process differs from metrics".into()),
        Err(error) => failures.push(format!("metrics replay: {error}")),
    }
    failures
}

fn validate_command(process: &ProcessMetrics) -> Vec<String> {
    let mut failures = Vec::new();
    let required = [
        "RAYON_NUM_THREADS=12",
        "PQ_F4_X1_BATCH=512",
        "F4_F2_DENSE_PAIR_SELECT=1",
        "F4_F2_DISABLE_INNER_BUILD_PARALLEL=1",
        "F4_F2_MATRIX_PROFILE=0",
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ];
    if process.command.first().map(String::as_str) != Some("/usr/bin/env")
        || required
            .iter()
            .any(|required| !process.command.iter().any(|arg| arg == required))
    {
        failures.push("frozen environment".into());
    }
    if process.command.iter().any(|arg| {
        arg.starts_with("F4_F2_FULL_M4RI=") || arg.starts_with("F4_F2_FULL_M4RI_MIN_ROWS=")
    }) {
        failures.push("default replay contains M4RI override".into());
    }
    if process.command.len() < 4
        || !process.command[process.command.len() - 4].ends_with("koblitz_pdp_backend")
        || process.command[process.command.len() - 3] != "native-f4"
        || !process.command[process.command.len() - 2]
            .ends_with("stage-175-current-f4-single-target-20261001/input/manifest.json")
        || process.command.last().map(String::as_str) != Some("300")
    {
        failures.push("backend invocation".into());
    }
    failures
}

fn validate_report(report: &Value) -> Vec<String> {
    let mut failures = Vec::new();
    for (pointer, expected) in [
        ("/status", json!("unsat")),
        ("/exhaustive", json!(true)),
        ("/source_instance_verified", json!(true)),
        ("/regenerated_source_exact", json!(true)),
        ("/source_instance_id", json!(SOURCE_INSTANCE_ID)),
        ("/solver_equations_blake3", json!(EQUATION_FINGERPRINT)),
        ("/fixed_x1_values", json!(512)),
        ("/fixed_x1_masks_visited", json!(512)),
        ("/fixed_x1_nonrational_skipped", json!(270)),
        ("/fixed_x1_systems_constructed", json!(242)),
        ("/fixed_x1_systems_completed", json!(242)),
        ("/algebraic_roots", json!(0)),
        ("/conflicts", Value::Null),
        ("/cost/timed_out", json!(false)),
        ("/cost/extra/full_m4ri_matrices", json!(481)),
        ("/cost/extra/full_m4ri_blocks", json!(351164)),
        ("/cost/extra/word_xors_performed", json!(102707985015u64)),
        ("/solver_full_m4ri_min_rows", json!(4096)),
        (
            "/solver_echelon_policy",
            json!("selected hybrid: full-matrix block-8 M4RI on eligible shapes from 4096 rows; current BlockTables elsewhere"),
        ),
        (
            "/factor_base_contract/target_subgroup_enumerated",
            json!(false),
        ),
        (
            "/factor_base_contract/discrete_log_labels_used",
            json!(false),
        ),
    ] {
        if report.pointer(pointer) != Some(&expected) {
            failures.push(format!("{pointer}: mismatch"));
        }
    }
    let profiled = [
        "lt256",
        "r256_511",
        "r512_1023",
        "r1024_2047",
        "r2048_4095",
        "ge4096",
    ]
    .iter()
    .any(|label| {
        pointer_u64(report, &format!("/cost/extra/matrix_profile_{label}_count")).unwrap_or(0) != 0
    });
    if profiled
        || report
            .pointer("/solver_matrix_profile")
            .and_then(Value::as_bool)
            != Some(false)
    {
        failures.push("profile must be disabled".into());
    }
    failures
}

fn verify_replay(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: DefaultReplay =
        serde_json::from_slice(&bytes).map_err(|e| format!("replay JSON: {e}"))?;
    let directory = path
        .parent()
        .ok_or_else(|| "replay has no parent".to_string())?;
    let replay = compose_replay(directory)?;
    let mut checks = 0usize;
    let mut failed = Vec::new();
    check(
        &mut checks,
        &mut failed,
        result.schema == REPLAY_SCHEMA,
        "schema",
    );
    check(
        &mut checks,
        &mut failed,
        result.inherited_stage198_result_sha256 == STAGE198_SHA256,
        "Stage 198 binding",
    );
    check(
        &mut checks,
        &mut failed,
        result.correct && replay.correct,
        "correctness",
    );
    check(
        &mut checks,
        &mut failed,
        result.errors == replay.errors,
        "error replay",
    );
    check(
        &mut checks,
        &mut failed,
        result.receipt == replay.receipt,
        "receipt replay",
    );
    check(
        &mut checks,
        &mut failed,
        result.report == replay.report,
        "report replay",
    );
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failed.len()),
        failed,
    })
}

fn finalize(
    stage: &Path,
    out: &Path,
    selected_commit: &str,
    finalizer_commit: &str,
) -> AnyResult<()> {
    validate_commit("selected", selected_commit)?;
    validate_commit("finalizer", finalizer_commit)?;
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    verify_stage198_parent(&stage)?;
    let replay_path = stage.join("development/default-replay/result.json");
    let replay: DefaultReplay = read_json(&replay_path)?;
    let verification = verify_replay(&replay_path)?;
    if !replay.correct || !verification.failed.is_empty() {
        return Err("default replay is not verified".into());
    }
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 682 + measured_stage_charge.components,
        wall_seconds_sum: 25_322.254_375_247 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 67_684.185_543 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    let mut artifacts = BTreeMap::new();
    for (name, relative) in [
        ("protocol", "PROTOCOL.md"),
        ("replay_result", "development/default-replay/result.json"),
        (
            "replay_verification",
            "development/default-replay/verification.json",
        ),
    ] {
        artifacts.insert(
            name.into(),
            artifact_relative(&stage, &stage.join(relative))?,
        );
    }
    write_json(
        out,
        &FinalResult {
            schema: FINAL_SCHEMA.into(),
            selected_commit: selected_commit.into(),
            finalizer_commit: finalizer_commit.into(),
            decision_status: "SELECTED_FOR_REPOSITORY_DEFAULT".into(),
            runtime_default: "unset selects full M4RI on eligible matrices from 4096 rows and five-column BlockTables elsewhere; F4_F2_FULL_M4RI=0 preserves the prior control; =1 preserves broad full M4RI from 128 rows".into(),
            claim_boundary: "One opened n=59 F4 default selection and unset replay only; no relation-yield, full-DLP, rho-crossover, external-reproduction, novelty, or SOTA claim.".into(),
            conflicts: None,
            single_core_seconds: None,
            default_replay: replay,
            measured_stage_charge,
            cumulative_measured_lower_bound,
            complete_campaign_cost: None,
            artifacts,
            all_seven_gates_passed: false,
            koblitz_index_calculus_sota: false,
        },
    )
}

fn verify_final(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: FinalResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("final JSON: {e}"))?;
    let stage_root = find_stage_root(path)?;
    let stage = stage_root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage_root.display()))?;
    let mut checks = 0usize;
    let mut failed = Vec::new();
    check(
        &mut checks,
        &mut failed,
        result.schema == FINAL_SCHEMA,
        "schema",
    );
    check(
        &mut checks,
        &mut failed,
        result.decision_status == "SELECTED_FOR_REPOSITORY_DEFAULT",
        "decision",
    );
    check(
        &mut checks,
        &mut failed,
        result.conflicts.is_none()
            && result.single_core_seconds.is_none()
            && result.complete_campaign_cost.is_none(),
        "null fields",
    );
    check(
        &mut checks,
        &mut failed,
        !result.all_seven_gates_passed && !result.koblitz_index_calculus_sota,
        "claim boundary",
    );
    check(
        &mut checks,
        &mut failed,
        result.default_replay.correct
            && result.default_replay.inherited_stage198_result_sha256 == STAGE198_SHA256,
        "default replay projection",
    );
    check(
        &mut checks,
        &mut failed,
        verify_stage198_parent(&stage).is_ok(),
        "Stage 198 parent",
    );
    let replay_path = stage.join("development/default-replay/result.json");
    let replay_verification = verify_replay(&replay_path)?;
    check(
        &mut checks,
        &mut failed,
        replay_verification.failed.is_empty(),
        "replay verification",
    );
    let replay: DefaultReplay = read_json(&replay_path)?;
    check(
        &mut checks,
        &mut failed,
        replay == result.default_replay,
        "replay projection",
    );
    let charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failed,
        charges_match(&charge, &result.measured_stage_charge),
        "stage charge",
    );
    let cumulative = Charge {
        components: 682 + charge.components,
        wall_seconds_sum: 25_322.254_375_247 + charge.wall_seconds_sum,
        total_core_seconds_sum: 67_684.185_543 + charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    check(
        &mut checks,
        &mut failed,
        charges_match(&cumulative, &result.cumulative_measured_lower_bound),
        "cumulative charge",
    );
    for (name, artifact) in &result.artifacts {
        checks += 2;
        match artifact_relative(&stage, &stage.join(&artifact.path)) {
            Ok(actual) => {
                if actual.bytes != artifact.bytes {
                    failed.push(format!("{name} bytes"));
                }
                if actual.sha256 != artifact.sha256 {
                    failed.push(format!("{name} SHA-256"));
                }
            }
            Err(error) => failed.push(format!("{name}: {error}")),
        }
    }
    let source = fs::read_to_string("src/cryptanalysis/pq_f4_f2.rs")
        .map_err(|e| format!("current F4 source: {e}"))?;
    let backend = fs::read_to_string("examples/koblitz_pdp_backend.rs")
        .map_err(|e| format!("backend source: {e}"))?;
    check(
        &mut checks,
        &mut failed,
        source.contains("full_m4ri_enabled_for(None)")
            && source.contains("full_m4ri_min_rows_for(None, None), 4096")
            && backend.contains("selected hybrid: full-matrix block-8 M4RI"),
        "selected source policy",
    );
    Ok(Verification {
        schema: FINAL_VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failed.len()),
        failed,
    })
}

fn measured_charge(development: &Path) -> AnyResult<Charge> {
    let mut metrics_files = Vec::new();
    collect_named_files(development, "metrics.json", &mut metrics_files)?;
    metrics_files.sort();
    let mut receipt_files = Vec::new();
    collect_named_files(development, "receipt.json", &mut receipt_files)?;
    receipt_files.sort();
    let mut covered = BTreeSet::new();
    for receipt_file in receipt_files {
        let directory = receipt_file
            .parent()
            .ok_or_else(|| format!("{} has no parent", receipt_file.display()))?;
        let receipt: MeterReceipt = read_json(&receipt_file)?;
        let errors = validate_receipt(directory, &receipt);
        if !errors.is_empty() {
            return Err(format!(
                "{} failed receipt replay: {}",
                receipt_file.display(),
                errors.join("; ")
            ));
        }
        covered.insert(directory.join(&receipt.metrics.path));
    }
    if covered != metrics_files.iter().cloned().collect::<BTreeSet<_>>() {
        return Err("development receipts and metrics differ".into());
    }
    let mut charge = Charge {
        components: metrics_files.len(),
        wall_seconds_sum: 0.0,
        total_core_seconds_sum: 0.0,
        peak_rss_bytes_max: 0,
        single_core_seconds: None,
    };
    for file in metrics_files {
        let metrics: ProcessMetrics = read_json(&file)?;
        charge.wall_seconds_sum += metrics.wall_seconds;
        charge.total_core_seconds_sum += metrics
            .total_core_seconds
            .ok_or_else(|| format!("{} lacks core time", file.display()))?;
        charge.peak_rss_bytes_max = charge.peak_rss_bytes_max.max(
            metrics
                .peak_rss_bytes
                .ok_or_else(|| format!("{} lacks RSS", file.display()))?,
        );
    }
    Ok(charge)
}

fn collect_named_files(directory: &Path, name: &str, out: &mut Vec<PathBuf>) -> AnyResult<()> {
    for entry in fs::read_dir(directory).map_err(|e| format!("{}: {e}", directory.display()))? {
        let entry = entry.map_err(|e| e.to_string())?;
        let kind = entry.file_type().map_err(|e| e.to_string())?;
        if kind.is_symlink() {
            return Err(format!("refusing symlink {}", entry.path().display()));
        }
        if kind.is_dir() {
            collect_named_files(&entry.path(), name, out)?;
        } else if kind.is_file() && entry.file_name() == name {
            out.push(entry.path());
        }
    }
    Ok(())
}

fn verify_stage198_parent(stage: &Path) -> AnyResult<()> {
    let path = stage
        .parent()
        .ok_or_else(|| "Stage 199 has no continuation parent".to_string())?
        .join("stage-198-hybrid-m4ri-confirm-20261002/result.json");
    let bytes = fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let actual = hex::encode(sha256(&bytes));
    if actual == STAGE198_SHA256 {
        Ok(())
    } else {
        Err(format!("Stage 198 result hash {actual}"))
    }
}

fn find_stage_root(path: &Path) -> AnyResult<&Path> {
    path.ancestors()
        .skip(1)
        .find(|candidate| {
            candidate.join("PROTOCOL.md").is_file() && candidate.join("development").is_dir()
        })
        .ok_or_else(|| format!("{} has no Stage 199 root", path.display()))
}

fn artifact_relative(root: &Path, path: &Path) -> AnyResult<Artifact> {
    let relative = path
        .strip_prefix(root)
        .map_err(|_| format!("{} escapes {}", path.display(), root.display()))?;
    artifact_from_file(path, relative.display().to_string())
}

fn artifact_at(directory: &Path, artifact: &Artifact) -> AnyResult<Artifact> {
    let path = Path::new(&artifact.path);
    if path.is_absolute()
        || path
            .components()
            .any(|component| matches!(component, std::path::Component::ParentDir))
    {
        return Err(format!("artifact path escapes: {}", artifact.path));
    }
    artifact_from_file(&directory.join(path), artifact.path.clone())
}

fn artifact_from_file(path: &Path, shown: String) -> AnyResult<Artifact> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: shown,
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
}

fn pointer_u64(value: &Value, pointer: &str) -> Option<u64> {
    value.pointer(pointer).and_then(Value::as_u64)
}

fn charges_match(left: &Charge, right: &Charge) -> bool {
    left.components == right.components
        && close(left.wall_seconds_sum, right.wall_seconds_sum)
        && close(left.total_core_seconds_sum, right.total_core_seconds_sum)
        && left.peak_rss_bytes_max == right.peak_rss_bytes_max
        && left.single_core_seconds == right.single_core_seconds
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn validate_commit(name: &str, commit: &str) -> AnyResult<()> {
    if commit.len() == 40 && commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        Ok(())
    } else {
        Err(format!("{name} commit must be a full SHA"))
    }
}

fn check(checks: &mut usize, failures: &mut Vec<String>, condition: bool, name: &str) {
    *checks += 1;
    if !condition {
        failures.push(name.into());
    }
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> AnyResult<T> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

fn write_json(path: &Path, value: &impl Serialize) -> AnyResult<()> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let mut file = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::to_writer_pretty(&mut file, value).map_err(|e| e.to_string())?;
    file.write_all(b"\n").map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strict_artifact_paths_reject_escape() {
        let artifact = Artifact {
            path: "../outside".into(),
            bytes: 0,
            sha256: String::new(),
        };
        assert!(artifact_at(Path::new("."), &artifact).is_err());
    }

    #[test]
    fn selected_policy_markers_are_specific() {
        let source = include_str!("../src/cryptanalysis/pq_f4_f2.rs");
        assert!(source.contains("full_m4ri_enabled_for(None)"));
        assert!(source.contains("full_m4ri_min_rows_for(None, None), 4096"));
    }
}
