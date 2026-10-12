//! Native composer/verifier for the Stage 200 parallel full-M4RI screen.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

const SCREEN_SCHEMA: &str = "koblitz_stage200_screen.v1";
const SCREEN_VERIFICATION_SCHEMA: &str = "koblitz_stage200_screen_verification.v1";
const FINAL_SCHEMA: &str = "koblitz_stage200_final.v1";
const FINAL_VERIFICATION_SCHEMA: &str = "koblitz_stage200_final_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";
const STAGE199_SHA256: &str = "2f4551a1644b5fc49d15a8be07a13ff4a5b56219595d496ef9e12b2077902edd";
const INITIAL_CANDIDATE_COMMIT: &str = "8fe333238";
const CORRECTED_TARGET: &str =
    "/Volumes/SSD990/kic-stage200-target-c2397b083/release/examples/koblitz_pdp_backend";

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Arm {
    SerialControl,
    ParallelCandidate,
}

impl Arm {
    fn directory(self) -> &'static str {
        match self {
            Self::SerialControl => "serial-control",
            Self::ParallelCandidate => "parallel-candidate",
        }
    }

    fn mode(self) -> &'static str {
        match self {
            Self::SerialControl => "F4_F2_FULL_M4RI_PARALLEL=0",
            Self::ParallelCandidate => "F4_F2_FULL_M4RI_PARALLEL=1",
        }
    }
}

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

#[derive(Clone, Debug, Deserialize, Serialize)]
struct RunRecord {
    order: usize,
    arm: Arm,
    correct: bool,
    errors: Vec<String>,
    receipt: MeterReceipt,
    report: Value,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Ratios {
    wall: f64,
    total_core: f64,
    peak_rss: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Assessment {
    all_correct: bool,
    candidate_over_control: Ratios,
    threshold: f64,
    passed_performance_gate: bool,
    status: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ScreenResult {
    schema: String,
    inherited_stage199_result_sha256: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    single_core_seconds: Option<f64>,
    runs: Vec<RunRecord>,
    assessment: Assessment,
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
    initial_candidate_commit: String,
    corrected_candidate_commit: String,
    revert_commit: String,
    finalizer_commit: String,
    decision_status: String,
    candidate_selected: bool,
    runtime_default: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    single_core_seconds: Option<f64>,
    screen: Assessment,
    failed_candidate_preserved_and_charged: bool,
    measured_stage_charge: Charge,
    cumulative_measured_lower_bound: Charge,
    complete_campaign_cost: Option<Charge>,
    unmetered_components: Vec<String>,
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
            let result = compose_screen(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &result)?;
            let verification = verify_screen(Path::new(&args[2]))?;
            write_json(Path::new(&args[3]), &verification)?;
            println!(
                "{}",
                serde_json::to_string(&json!({
                    "schema":"koblitz_stage200_terminal.v1",
                    "status":result.assessment.status,
                    "all_correct":result.assessment.all_correct,
                    "ratios":result.assessment.candidate_over_control,
                    "passed_performance_gate":result.assessment.passed_performance_gate,
                    "verification_checks":verification.checks,
                    "verification_passed":verification.passed,
                    "verification_failed":verification.failed,
                }))
                .map_err(|e| e.to_string())?
            );
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("screen verification failed".into())
            }
        }
        Some("verify") if args.len() == 3 => {
            let verification = verify_screen(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("screen verification failed".into())
            }
        }
        Some("finalize") if args.len() == 6 => finalize(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &args[3],
            &args[4],
            &args[5],
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
        _ => Err("usage: koblitz_f4_stage200 compose SCREEN_DIR RESULT VERIFICATION\n       koblitz_f4_stage200 verify RESULT VERIFICATION\n       koblitz_f4_stage200 finalize STAGE RESULT CORRECTED_COMMIT REVERT_COMMIT FINALIZER_COMMIT\n       koblitz_f4_stage200 verify-final RESULT VERIFICATION".into()),
    }
}

fn compose_screen(directory: &Path) -> AnyResult<ScreenResult> {
    let schedule = [Arm::SerialControl, Arm::ParallelCandidate];
    let mut runs = Vec::new();
    for (offset, arm) in schedule.into_iter().enumerate() {
        let run_dir = directory.join(arm.directory());
        let receipt: MeterReceipt = read_json(&run_dir.join("receipt.json"))?;
        let mut errors = validate_receipt(&run_dir, &receipt);
        errors.extend(validate_solver_process(&receipt.process, arm));
        let report = match read_json::<Value>(&run_dir.join(&receipt.stdout.path)) {
            Ok(report) => {
                errors.extend(validate_report(&report, arm));
                report
            }
            Err(error) => {
                errors.push(error);
                Value::Null
            }
        };
        runs.push(RunRecord {
            order: offset + 1,
            arm,
            correct: errors.is_empty(),
            errors,
            receipt,
            report,
        });
    }
    if runs.len() == 2 && normalized_report(&runs[0].report) != normalized_report(&runs[1].report) {
        runs[1]
            .errors
            .push("normalized solver report mismatch".into());
        runs[1].correct = false;
    }
    let assessment = assess(&runs)?;
    Ok(ScreenResult {
        schema: SCREEN_SCHEMA.into(),
        inherited_stage199_result_sha256: STAGE199_SHA256.into(),
        claim_boundary: "One opened n=59 F4 row-scheduling screen only; no relation-yield, full-DLP, rho-crossover, external-reproduction, novelty, or SOTA claim.".into(),
        conflicts: None,
        single_core_seconds: None,
        runs,
        assessment,
    })
}

fn validate_receipt(directory: &Path, receipt: &MeterReceipt) -> Vec<String> {
    let mut failures = Vec::new();
    if receipt.schema != RECEIPT_SCHEMA || receipt.process.schema != METER_SCHEMA {
        failures.push("meter schema".into());
    }
    if receipt.process.meter != "fresh child process wait4/getrusage; all child threads included"
        || !receipt.process.environment.is_empty()
        || receipt.process.watchdog_seconds == 0
        || receipt.process.single_core_seconds.is_some()
        || receipt.process.total_core_seconds.is_none()
        || receipt.process.peak_rss_bytes.is_none()
        || !receipt.process.wall_seconds.is_finite()
        || receipt.process.wall_seconds < 0.0
        || receipt
            .process
            .total_core_seconds
            .is_some_and(|value| !value.is_finite() || value < 0.0)
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

fn validate_solver_process(process: &ProcessMetrics, arm: Arm) -> Vec<String> {
    let mut failures = Vec::new();
    if process.returncode != 0 || process.timed_out || process.watchdog_seconds != 360 {
        failures.push("solver process terminal".into());
    }
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
    let modes = process
        .command
        .iter()
        .filter(|arg| arg.starts_with("F4_F2_FULL_M4RI_PARALLEL="))
        .collect::<Vec<_>>();
    if modes.len() != 1 || modes[0].as_str() != arm.mode() {
        failures.push("parallel-clear mode".into());
    }
    if process.command.iter().any(|arg| {
        arg.starts_with("F4_F2_FULL_M4RI=") || arg.starts_with("F4_F2_FULL_M4RI_MIN_ROWS=")
    }) {
        failures.push("selected-default M4RI override".into());
    }
    if process.command.len() < 4
        || process.command[process.command.len() - 4] != CORRECTED_TARGET
        || process.command[process.command.len() - 3] != "native-f4"
        || !process.command[process.command.len() - 2]
            .ends_with("stage-175-current-f4-single-target-20261001/input/manifest.json")
        || process.command.last().map(String::as_str) != Some("300")
    {
        failures.push("backend invocation".into());
    }
    failures
}

fn validate_report(report: &Value, arm: Arm) -> Vec<String> {
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
        ("/cost/ops", json!(318703372596u64)),
        ("/cost/extra/full_m4ri_matrices", json!(481)),
        ("/cost/extra/full_m4ri_blocks", json!(351164)),
        (
            "/cost/extra/full_m4ri_table_word_xors",
            json!(9222665017u64),
        ),
        ("/cost/extra/word_xors_performed", json!(102707985015u64)),
        ("/solver_full_m4ri_min_rows", json!(4096)),
        ("/solver_full_m4ri_parallel_min_words", json!(1048576)),
        ("/solver_matrix_profile", json!(false)),
        ("/single_thread_requested", json!(false)),
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
            failures.push(format!("{pointer}: terminal mismatch"));
        }
    }
    let expected_policy = match arm {
        Arm::SerialControl => "serial target-row clearing",
        Arm::ParallelCandidate => {
            "existing Rayon pool for target-row regions of at least 1048576 words"
        }
    };
    if report
        .pointer("/solver_full_m4ri_parallel_clear")
        .and_then(Value::as_str)
        != Some(expected_policy)
    {
        failures.push("parallel-clear policy report".into());
    }
    let expected_parallel = match arm {
        Arm::SerialControl => (0, 0, 0),
        Arm::ParallelCandidate => (75803, 611986047, 106712818981),
    };
    for (pointer, expected) in [
        ("/cost/extra/full_m4ri_parallel_blocks", expected_parallel.0),
        (
            "/cost/extra/full_m4ri_parallel_target_rows",
            expected_parallel.1,
        ),
        (
            "/cost/extra/full_m4ri_parallel_target_row_words",
            expected_parallel.2,
        ),
    ] {
        if pointer_u64(report, pointer) != Some(expected) {
            failures.push(format!("{pointer}: scheduler counter"));
        }
    }
    if report
        .pointer("/cost/extra")
        .and_then(Value::as_object)
        .is_none_or(|extra| {
            extra
                .iter()
                .any(|(key, value)| key.starts_with("matrix_profile_") && value.as_u64() != Some(0))
        })
    {
        failures.push("non-zero profiling counter".into());
    }
    failures
}

fn normalized_report(report: &Value) -> Value {
    let mut normalized = report.clone();
    if let Some(root) = normalized.as_object_mut() {
        root.remove("timing_ns");
        root.remove("solver_full_m4ri_parallel_clear");
        if let Some(cost) = root.get_mut("cost").and_then(Value::as_object_mut) {
            cost.remove("wall_ns");
            if let Some(extra) = cost.get_mut("extra").and_then(Value::as_object_mut) {
                for key in [
                    "build_ns",
                    "eliminate_ns",
                    "f4_wall_ns",
                    "full_m4ri_parallel_blocks",
                    "full_m4ri_parallel_target_rows",
                    "full_m4ri_parallel_target_row_words",
                ] {
                    extra.remove(key);
                }
            }
        }
    }
    normalized
}

fn assess(runs: &[RunRecord]) -> AnyResult<Assessment> {
    let control = runs
        .iter()
        .find(|run| run.arm == Arm::SerialControl)
        .ok_or_else(|| "serial control missing".to_string())?;
    let candidate = runs
        .iter()
        .find(|run| run.arm == Arm::ParallelCandidate)
        .ok_or_else(|| "parallel candidate missing".to_string())?;
    let ratios = Ratios {
        wall: candidate.receipt.process.wall_seconds / control.receipt.process.wall_seconds,
        total_core: candidate
            .receipt
            .process
            .total_core_seconds
            .ok_or_else(|| "candidate core unavailable".to_string())?
            / control
                .receipt
                .process
                .total_core_seconds
                .ok_or_else(|| "control core unavailable".to_string())?,
        peak_rss: candidate
            .receipt
            .process
            .peak_rss_bytes
            .ok_or_else(|| "candidate RSS unavailable".to_string())? as f64
            / control
                .receipt
                .process
                .peak_rss_bytes
                .ok_or_else(|| "control RSS unavailable".to_string())? as f64,
    };
    let all_correct = runs.len() == 2 && runs.iter().all(|run| run.correct);
    let threshold = 0.98;
    let passed = all_correct && ratios.wall < threshold && ratios.total_core < threshold;
    Ok(Assessment {
        all_correct,
        candidate_over_control: ratios,
        threshold,
        passed_performance_gate: passed,
        status: if all_correct && passed {
            "CONTINUE_TO_CONFIRMATION"
        } else if all_correct {
            "REJECTED_SCREEN"
        } else {
            "REJECTED_CORRECTNESS"
        }
        .into(),
    })
}

fn verify_screen(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: ScreenResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("screen JSON: {e}"))?;
    let directory = path
        .parent()
        .ok_or_else(|| "screen result has no parent".to_string())?;
    let mut checks = 0usize;
    let mut failed = Vec::new();
    check(
        &mut checks,
        &mut failed,
        result.schema == SCREEN_SCHEMA,
        "schema",
    );
    check(
        &mut checks,
        &mut failed,
        result.inherited_stage199_result_sha256 == STAGE199_SHA256,
        "Stage 199 projection",
    );
    check(
        &mut checks,
        &mut failed,
        result.conflicts.is_none() && result.single_core_seconds.is_none(),
        "null fields",
    );
    for run in &result.runs {
        let run_dir = directory.join(run.arm.directory());
        let mut errors = validate_receipt(&run_dir, &run.receipt);
        errors.extend(validate_solver_process(&run.receipt.process, run.arm));
        errors.extend(validate_report(&run.report, run.arm));
        check(
            &mut checks,
            &mut failed,
            errors == run.errors && run.correct == errors.is_empty(),
            &format!("{:?} replay", run.arm),
        );
    }
    check(
        &mut checks,
        &mut failed,
        result.runs.len() == 2
            && normalized_report(&result.runs[0].report)
                == normalized_report(&result.runs[1].report),
        "normalized report equality",
    );
    let assessment = assess(&result.runs)?;
    check(
        &mut checks,
        &mut failed,
        assessments_match(&assessment, &result.assessment),
        "assessment replay",
    );
    check(
        &mut checks,
        &mut failed,
        result.assessment.status == "REJECTED_SCREEN"
            && result.assessment.all_correct
            && !result.assessment.passed_performance_gate,
        "frozen screen decision",
    );
    Ok(Verification {
        schema: SCREEN_VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failed.len()),
        failed,
    })
}

fn finalize(
    stage: &Path,
    out: &Path,
    corrected_commit: &str,
    revert_commit: &str,
    finalizer_commit: &str,
) -> AnyResult<()> {
    for (name, commit) in [
        ("corrected", corrected_commit),
        ("revert", revert_commit),
        ("finalizer", finalizer_commit),
    ] {
        validate_commit(name, commit)?;
    }
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    verify_stage199_parent(&stage)?;
    let screen_path = stage.join("development/corrected-screen/result.json");
    let screen: ScreenResult = read_json(&screen_path)?;
    let replay = verify_screen(&screen_path)?;
    if !replay.failed.is_empty() || screen.assessment.status != "REJECTED_SCREEN" {
        return Err("screen is not a verified rejection".into());
    }
    let failed_receipt: MeterReceipt =
        read_json(&stage.join("development/screen/parallel-candidate/receipt.json"))?;
    if failed_receipt.process.returncode != 101
        || failed_receipt.process.timed_out
        || failed_receipt.process.wall_seconds <= 0.0
    {
        return Err("first candidate failure is not preserved".into());
    }
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 696 + measured_stage_charge.components,
        wall_seconds_sum: 25_527.893_758_747_003 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 68_523.607_063 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    let artifact_paths = [
        ("protocol", "PROTOCOL.md"),
        ("candidate_failure", "CANDIDATE_FAILURE.md"),
        ("candidate_patch", "candidate-parallel-m4ri.patch"),
        ("screen_result", "development/corrected-screen/result.json"),
        (
            "screen_verification",
            "development/corrected-screen/verification.json",
        ),
        (
            "failed_candidate_receipt",
            "development/screen/parallel-candidate/receipt.json",
        ),
        (
            "failed_candidate_stderr",
            "development/screen/parallel-candidate/stderr.txt",
        ),
    ];
    let mut artifacts = BTreeMap::new();
    for (name, relative) in artifact_paths {
        artifacts.insert(
            name.into(),
            artifact_relative(&stage, &stage.join(relative))?,
        );
    }
    write_json(
        out,
        &FinalResult {
            schema: FINAL_SCHEMA.into(),
            initial_candidate_commit: INITIAL_CANDIDATE_COMMIT.into(),
            corrected_candidate_commit: corrected_commit.into(),
            revert_commit: revert_commit.into(),
            finalizer_commit: finalizer_commit.into(),
            decision_status: screen.assessment.status.clone(),
            candidate_selected: false,
            runtime_default: "Stage 199 selected hybrid: full-matrix block-8 M4RI on eligible matrices from 4096 rows with serial target-row clearing; five-column BlockTables elsewhere".into(),
            claim_boundary: "One opened n=59 F4 row-scheduling screen. Direct MITM and full automorphism-aware rho remain unchanged; no relation-yield, full-DLP, external-reproduction, novelty, or SOTA claim.".into(),
            conflicts: None,
            single_core_seconds: None,
            screen: screen.assessment,
            failed_candidate_preserved_and_charged: true,
            measured_stage_charge,
            cumulative_measured_lower_bound,
            complete_campaign_cost: None,
            unmetered_components: vec![
                "final result write and final verification after the charged composition control"
                    .into(),
            ],
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
        result.decision_status == "REJECTED_SCREEN"
            && !result.candidate_selected
            && result.screen.status == result.decision_status,
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
        result.failed_candidate_preserved_and_charged,
        "failed candidate accounting",
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
        verify_stage199_parent(&stage).is_ok(),
        "Stage 199 parent",
    );
    for (name, commit) in [
        ("corrected", result.corrected_candidate_commit.as_str()),
        ("revert", result.revert_commit.as_str()),
        ("finalizer", result.finalizer_commit.as_str()),
    ] {
        check(
            &mut checks,
            &mut failed,
            validate_commit(name, commit).is_ok(),
            &format!("{name} commit"),
        );
    }
    let screen_path = stage.join("development/corrected-screen/result.json");
    let screen: ScreenResult = read_json(&screen_path)?;
    let screen_verification = verify_screen(&screen_path)?;
    check(
        &mut checks,
        &mut failed,
        screen_verification.failed.is_empty()
            && assessments_match(&screen.assessment, &result.screen),
        "screen replay",
    );
    let charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failed,
        charges_match(&charge, &result.measured_stage_charge),
        "stage charge",
    );
    let cumulative = Charge {
        components: 696 + charge.components,
        wall_seconds_sum: 25_527.893_758_747_003 + charge.wall_seconds_sum,
        total_core_seconds_sum: 68_523.607_063 + charge.total_core_seconds_sum,
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
        !source.contains("F4_F2_FULL_M4RI_PARALLEL")
            && !source.contains("full_m4ri_parallel_blocks")
            && !backend.contains("solver_full_m4ri_parallel_clear")
            && source.contains("full_m4ri_min_rows_for(None, None), 4096")
            && backend.contains("selected hybrid: full-matrix block-8 M4RI"),
        "runtime candidate reverted",
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

fn verify_stage199_parent(stage: &Path) -> AnyResult<()> {
    let path = stage
        .parent()
        .ok_or_else(|| "Stage 200 has no continuation parent".to_string())?
        .join("stage-199-hybrid-m4ri-default-20261002/result.json");
    let bytes = fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let actual = hex::encode(sha256(&bytes));
    if actual == STAGE199_SHA256 {
        Ok(())
    } else {
        Err(format!("Stage 199 result hash {actual}"))
    }
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

fn find_stage_root(path: &Path) -> AnyResult<&Path> {
    path.ancestors()
        .skip(1)
        .find(|candidate| {
            candidate.join("PROTOCOL.md").is_file() && candidate.join("development").is_dir()
        })
        .ok_or_else(|| format!("{} has no Stage 200 root", path.display()))
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

fn assessments_match(left: &Assessment, right: &Assessment) -> bool {
    left.all_correct == right.all_correct
        && close(
            left.candidate_over_control.wall,
            right.candidate_over_control.wall,
        )
        && close(
            left.candidate_over_control.total_core,
            right.candidate_over_control.total_core,
        )
        && close(
            left.candidate_over_control.peak_rss,
            right.candidate_over_control.peak_rss,
        )
        && close(left.threshold, right.threshold)
        && left.passed_performance_gate == right.passed_performance_gate
        && left.status == right.status
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
    fn rejected_candidate_is_absent_from_runtime_source() {
        let source = include_str!("../src/cryptanalysis/pq_f4_f2.rs");
        let backend = include_str!("koblitz_pdp_backend.rs");
        assert!(!source.contains("F4_F2_FULL_M4RI_PARALLEL"));
        assert!(!source.contains("full_m4ri_parallel_blocks"));
        assert!(!backend.contains("solver_full_m4ri_parallel_clear"));
        assert!(source.contains("full_m4ri_min_rows_for(None, None), 4096"));
    }
}
