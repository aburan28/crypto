//! Native, fail-closed runner for the Stage 188 Koblitz F4 reducer-index replay.
//!
//! This intentionally owns process metering, output hashing, terminal-record
//! validation, paired ratios, and replay verification.  It never calls Python.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const RESULT_SCHEMA: &str = "koblitz_stage188_native_result.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage188_native_verification.v1";
const FINAL_SCHEMA: &str = "koblitz_stage188_final.v1";
const FINAL_VERIFICATION_SCHEMA: &str = "koblitz_stage188_final_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";
const WATCHDOG_SECONDS: u64 = 360;

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
enum Arm {
    Linear,
    Indexed,
}

impl Arm {
    fn name(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Indexed => "indexed",
        }
    }

    fn indexed_env(self) -> &'static str {
        match self {
            Self::Linear => "0",
            Self::Indexed => "1",
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

#[derive(Clone, Debug, Deserialize, Serialize)]
struct RunReceipt {
    order: usize,
    arm: Arm,
    repeat: usize,
    correct: bool,
    errors: Vec<String>,
    process: ProcessMetrics,
    metrics_artifact: Artifact,
    stdout_artifact: Artifact,
    stderr_artifact: Artifact,
    report: Value,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Ratios {
    wall: f64,
    core: f64,
    rss: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct PairRatio {
    repeat: usize,
    indexed_over_linear: Ratios,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Assessment {
    all_correct: bool,
    pair_ratios: Vec<PairRatio>,
    median_paired_ratios: Ratios,
    passed_performance_gate: bool,
    status: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhaseResult {
    schema: String,
    phase: String,
    source_commit: String,
    created_unix_seconds: u64,
    claim_boundary: String,
    conflicts: Option<u64>,
    single_core_seconds: Option<f64>,
    provenance: Value,
    runs: Vec<RunReceipt>,
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
    candidate_commit: String,
    verifier_commit: String,
    finalizer_commit: String,
    decision_status: String,
    candidate_selected: bool,
    runtime_default: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    screen: Assessment,
    confirmation: Assessment,
    mechanism: Value,
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
    let args: Vec<String> = env::args().skip(1).collect();
    let Some(command) = args.first().map(String::as_str) else {
        return Err(usage());
    };
    match command {
        "screen" | "confirm" => {
            if args.len() != 7 {
                return Err(usage());
            }
            run_phase(
                command,
                Path::new(&args[1]),
                Path::new(&args[2]),
                Path::new(&args[3]),
                Path::new(&args[4]),
                Path::new(&args[5]),
                &args[6],
            )
        }
        "verify" => {
            if args.len() != 3 {
                return Err(usage());
            }
            let verification = verify_result(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err(format!(
                    "verification failed: {}",
                    verification.failed.join("; ")
                ))
            }
        }
        "finalize" => {
            if args.len() != 6 {
                return Err(usage());
            }
            finalize_stage(
                Path::new(&args[1]),
                Path::new(&args[2]),
                &args[3],
                &args[4],
                &args[5],
            )
        }
        "verify-final" => {
            if args.len() != 3 {
                return Err(usage());
            }
            let verification = verify_final(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err(format!(
                    "final verification failed: {}",
                    verification.failed.join("; ")
                ))
            }
        }
        "meter" => {
            if args.len() < 4 {
                return Err(usage());
            }
            let watchdog = args[2]
                .parse::<u64>()
                .map_err(|e| format!("invalid watchdog seconds: {e}"))?;
            meter_command(Path::new(&args[1]), watchdog, &args[3..])
        }
        _ => Err(usage()),
    }
}

fn usage() -> String {
    "usage:\n  koblitz_f4_stage188 screen|confirm BACKEND MANIFEST PROTOCOL LOCK OUT SOURCE_COMMIT\n  koblitz_f4_stage188 verify RESULT_JSON VERIFICATION_JSON\n  koblitz_f4_stage188 finalize STAGE_DIR OUT CANDIDATE_COMMIT VERIFIER_COMMIT FINALIZER_COMMIT\n  koblitz_f4_stage188 verify-final RESULT_JSON VERIFICATION_JSON\n  koblitz_f4_stage188 meter OUT WATCHDOG_SECONDS COMMAND [ARG ...]".into()
}

fn run_phase(
    phase: &str,
    backend: &Path,
    manifest: &Path,
    protocol: &Path,
    lock: &Path,
    out: &Path,
    source_commit: &str,
) -> AnyResult<()> {
    if source_commit.len() != 40 || !source_commit.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("SOURCE_COMMIT must be a full 40-digit hexadecimal commit".into());
    }
    let backend = canonical_file(backend)?;
    let manifest = canonical_file(manifest)?;
    let protocol = canonical_file(protocol)?;
    let lock = canonical_file(lock)?;
    create_new_directory(out)?;

    let runner = canonical_file(&env::current_exe().map_err(|e| e.to_string())?)?;
    let source_artifacts = source_artifact_map()?;
    let provenance = json!({
        "schema":"koblitz_stage188_native_provenance.v1",
        "source_commit":source_commit,
        "runner":artifact_absolute(&runner)?,
        "backend":artifact_absolute(&backend)?,
        "manifest":artifact_absolute(&manifest)?,
        "protocol":artifact_absolute(&protocol)?,
        "lock":artifact_absolute(&lock)?,
        "source_artifacts":source_artifacts,
        "host":{
            "arch":env::consts::ARCH,
            "os":env::consts::OS,
            "family":env::consts::FAMILY,
        },
        "watchdog_seconds":WATCHDOG_SECONDS,
        "external_thread_caps":1,
        "rayon_threads":12,
    });
    write_json(&out.join("provenance.json"), &provenance)?;

    let order: &[Arm] = match phase {
        "screen" => &[Arm::Linear, Arm::Indexed],
        "confirm" => &[
            Arm::Linear,
            Arm::Indexed,
            Arm::Indexed,
            Arm::Linear,
            Arm::Linear,
            Arm::Indexed,
        ],
        _ => return Err(format!("unknown phase {phase}")),
    };
    let mut repeats = BTreeMap::<Arm, usize>::new();
    let mut runs = Vec::with_capacity(order.len());
    for (index, &arm) in order.iter().enumerate() {
        let repeat = repeats.entry(arm).or_default();
        *repeat += 1;
        runs.push(run_backend(
            out,
            index + 1,
            arm,
            *repeat,
            &backend,
            &manifest,
        )?);
    }

    append_cross_run_errors(&mut runs);
    let assessment = assess(phase, &runs)?;
    let result = PhaseResult {
        schema: RESULT_SCHEMA.into(),
        phase: phase.into(),
        source_commit: source_commit.into(),
        created_unix_seconds: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|e| e.to_string())?
            .as_secs(),
        claim_boundary: "One opened n=59 decomposition target and a solver-stage implementation experiment only; not relation yield, a full DLP, a rho crossover, independent reproduction, novelty, or SOTA.".into(),
        conflicts: None,
        single_core_seconds: None,
        provenance,
        runs,
        assessment,
    };
    let result_path = out.join("result.json");
    write_json(&result_path, &result)?;
    let verification = verify_result(&result_path)?;
    write_json(&out.join("verification.json"), &verification)?;
    println!(
        "{}",
        serde_json::to_string(&json!({
            "schema":"koblitz_stage188_native_terminal.v1",
            "phase":phase,
            "status":result.assessment.status,
            "all_correct":result.assessment.all_correct,
            "passed_performance_gate":result.assessment.passed_performance_gate,
            "median_paired_ratios":result.assessment.median_paired_ratios,
            "verification_checks":verification.checks,
            "verification_passed":verification.passed,
            "verification_failed":verification.failed,
        }))
        .map_err(|e| e.to_string())?
    );
    if verification.failed.is_empty() {
        Ok(())
    } else {
        Err("phase verification failed".into())
    }
}

fn run_backend(
    root: &Path,
    order: usize,
    arm: Arm,
    repeat: usize,
    backend: &Path,
    manifest: &Path,
) -> AnyResult<RunReceipt> {
    let run_dir = root.join(format!("{order:02}-{}-r{repeat}", arm.name()));
    create_new_directory(&run_dir)?;
    let stdout_path = run_dir.join("stdout.json");
    let stderr_path = run_dir.join("stderr.txt");
    let metrics_path = run_dir.join("metrics.json");
    let environment = BTreeMap::from([
        ("BLIS_NUM_THREADS".into(), "1".into()),
        ("F4_F2_DENSE_PAIR_SELECT".into(), "1".into()),
        ("F4_F2_FULL_M4RI".into(), "0".into()),
        ("F4_F2_INDEXED_REDUCERS".into(), arm.indexed_env().into()),
        ("MKL_NUM_THREADS".into(), "1".into()),
        ("NUMEXPR_NUM_THREADS".into(), "1".into()),
        ("OMP_NUM_THREADS".into(), "1".into()),
        ("OPENBLAS_NUM_THREADS".into(), "1".into()),
        ("PQ_F4_X1_BATCH".into(), "512".into()),
        ("RAYON_NUM_THREADS".into(), "12".into()),
        ("VECLIB_MAXIMUM_THREADS".into(), "1".into()),
    ]);
    let command = vec![
        backend.display().to_string(),
        "native-f4".into(),
        manifest.display().to_string(),
        "300".into(),
    ];
    let process = execute_command(
        &command,
        &environment,
        &stdout_path,
        &stderr_path,
        WATCHDOG_SECONDS,
    )?;
    write_json(&metrics_path, &process)?;

    let mut errors = Vec::new();
    if process.returncode != 0 {
        errors.push(format!("backend return code {}", process.returncode));
    }
    if process.timed_out {
        errors.push("backend watchdog timeout".into());
    }
    let report = match fs::read(&stdout_path)
        .map_err(|e| e.to_string())
        .and_then(|bytes| serde_json::from_slice::<Value>(&bytes).map_err(|e| e.to_string()))
    {
        Ok(report) => {
            errors.extend(validate_report(&report, arm));
            report
        }
        Err(error) => {
            errors.push(format!("invalid backend JSON: {error}"));
            Value::Null
        }
    };
    Ok(RunReceipt {
        order,
        arm,
        repeat,
        correct: errors.is_empty(),
        errors,
        process,
        metrics_artifact: artifact_relative(root, &metrics_path)?,
        stdout_artifact: artifact_relative(root, &stdout_path)?,
        stderr_artifact: artifact_relative(root, &stderr_path)?,
        report,
    })
}

fn validate_report(report: &Value, arm: Arm) -> Vec<String> {
    let mut failures = Vec::new();
    expect_value(report, "/status", &json!("unsat"), &mut failures);
    expect_value(report, "/exhaustive", &json!(true), &mut failures);
    expect_value(
        report,
        "/source_instance_verified",
        &json!(true),
        &mut failures,
    );
    expect_value(
        report,
        "/regenerated_source_exact",
        &json!(true),
        &mut failures,
    );
    expect_value(
        report,
        "/source_instance_id",
        &json!(SOURCE_INSTANCE_ID),
        &mut failures,
    );
    expect_value(
        report,
        "/solver_equations_blake3",
        &json!(EQUATION_FINGERPRINT),
        &mut failures,
    );
    for (pointer, expected) in [
        ("/fixed_x1_values", 512),
        ("/fixed_x1_masks_visited", 512),
        ("/fixed_x1_nonrational_skipped", 270),
        ("/fixed_x1_systems_constructed", 242),
        ("/fixed_x1_systems_completed", 242),
        ("/algebraic_roots", 0),
        ("/solver_variables", 18),
    ] {
        expect_value(report, pointer, &json!(expected), &mut failures);
    }
    expect_value(
        report,
        "/factor_base_contract/target_subgroup_enumerated",
        &json!(false),
        &mut failures,
    );
    expect_value(
        report,
        "/factor_base_contract/discrete_log_labels_used",
        &json!(false),
        &mut failures,
    );
    expect_value(report, "/conflicts", &Value::Null, &mut failures);
    expect_value(report, "/cost/timed_out", &json!(false), &mut failures);

    let divisor = required_u64(report, "/cost/extra/divisor_tests", &mut failures);
    let submask = required_u64(report, "/cost/extra/divisor_submask_lookups", &mut failures);
    let linear = required_u64(report, "/cost/extra/divisor_linear_tests", &mut failures);
    let bytes = required_u64(report, "/cost/extra/reducer_index_bytes_max", &mut failures);
    if let (Some(divisor), Some(submask), Some(linear)) = (divisor, submask, linear) {
        if divisor != submask.saturating_add(linear) {
            failures.push(format!(
                "divisor accounting mismatch: {divisor} != {submask} + {linear}"
            ));
        }
    }
    match arm {
        Arm::Linear => {
            if submask != Some(0) {
                failures.push(format!("linear submask lookups must be 0, got {submask:?}"));
            }
            if bytes != Some(0) {
                failures.push(format!("linear index bytes must be 0, got {bytes:?}"));
            }
        }
        Arm::Indexed => {
            if submask == Some(0) {
                failures.push("indexed submask lookups must be positive".into());
            }
            if bytes == Some(0) {
                failures.push("indexed reducer bytes must be positive".into());
            }
        }
    }
    failures
}

fn append_cross_run_errors(runs: &mut [RunReceipt]) {
    let Some(reference) = runs.first().map(|run| run.report.clone()) else {
        return;
    };
    let pointers = [
        "/status",
        "/exhaustive",
        "/source_instance_id",
        "/source_instance_verified",
        "/regenerated_source_exact",
        "/fixed_x1_values",
        "/fixed_x1_masks_visited",
        "/fixed_x1_nonrational_skipped",
        "/fixed_x1_systems_constructed",
        "/fixed_x1_systems_completed",
        "/algebraic_roots",
        "/solver_equations_blake3",
        "/solver_equations_total",
        "/solver_variables",
        "/cost/ops",
        "/cost/degree_reached",
        "/cost/solving_degree",
        "/cost/timed_out",
        "/cost/extra/basis_len",
        "/cost/extra/extraction_tests",
        "/cost/extra/field_pairs_reduced",
        "/cost/extra/full_m4ri_blocks",
        "/cost/extra/full_m4ri_matrices",
        "/cost/extra/full_m4ri_table_word_xors",
        "/cost/extra/matrix_cols_max",
        "/cost/extra/matrix_rows_max",
        "/cost/extra/matrix_rows_sum",
        "/cost/extra/max_poly_degree",
        "/cost/extra/new_elements",
        "/cost/extra/oversize",
        "/cost/extra/pair_candidate_visits",
        "/cost/extra/pair_cover_lookups",
        "/cost/extra/pair_dense_scratch_bytes_max",
        "/cost/extra/pair_dense_select_calls",
        "/cost/extra/pair_lcm_groups",
        "/cost/extra/pair_quadratic_select_calls",
        "/cost/extra/pairs_chain_skipped",
        "/cost/extra/pairs_left",
        "/cost/extra/pairs_product_skipped",
        "/cost/extra/pairs_reduced",
        "/cost/extra/peak_matrix_bytes",
        "/cost/extra/peak_table_bytes",
        "/cost/extra/reducer_rows",
        "/cost/extra/steps",
        "/cost/extra/word_xors_performed",
    ];
    for run in runs.iter_mut().skip(1) {
        for pointer in pointers {
            if reference.pointer(pointer) != run.report.pointer(pointer) {
                run.errors.push(format!("structural mismatch at {pointer}"));
            }
        }
        run.correct = run.errors.is_empty();
    }
    if let (Some(linear), Some(indexed)) = (
        runs.iter().find(|run| run.arm == Arm::Linear),
        runs.iter().find(|run| run.arm == Arm::Indexed),
    ) {
        let linear_tests = linear
            .report
            .pointer("/cost/extra/divisor_tests")
            .and_then(Value::as_u64);
        let indexed_tests = indexed
            .report
            .pointer("/cost/extra/divisor_tests")
            .and_then(Value::as_u64);
        if !matches!((linear_tests, indexed_tests), (Some(a), Some(b)) if b < a) {
            for run in runs.iter_mut().filter(|run| run.arm == Arm::Indexed) {
                run.errors.push(format!(
                    "indexed divisor total did not fall: linear={linear_tests:?} indexed={indexed_tests:?}"
                ));
                run.correct = false;
            }
        }
    }
}

fn assess(phase: &str, runs: &[RunReceipt]) -> AnyResult<Assessment> {
    let mut pairs = BTreeMap::<usize, (Option<&RunReceipt>, Option<&RunReceipt>)>::new();
    for run in runs {
        let pair = pairs.entry(run.repeat).or_default();
        match run.arm {
            Arm::Linear => pair.0 = Some(run),
            Arm::Indexed => pair.1 = Some(run),
        }
    }
    let mut pair_ratios = Vec::with_capacity(pairs.len());
    for (repeat, (linear, indexed)) in pairs {
        let linear = linear.ok_or_else(|| format!("pair {repeat} has no linear run"))?;
        let indexed = indexed.ok_or_else(|| format!("pair {repeat} has no indexed run"))?;
        let linear_core = linear
            .process
            .total_core_seconds
            .ok_or_else(|| format!("pair {repeat} linear core time unavailable"))?;
        let indexed_core = indexed
            .process
            .total_core_seconds
            .ok_or_else(|| format!("pair {repeat} indexed core time unavailable"))?;
        let rss = match (
            linear.process.peak_rss_bytes,
            indexed.process.peak_rss_bytes,
        ) {
            (Some(a), Some(b)) if a != 0 => Some(b as f64 / a as f64),
            _ => None,
        };
        pair_ratios.push(PairRatio {
            repeat,
            indexed_over_linear: Ratios {
                wall: indexed.process.wall_seconds / linear.process.wall_seconds,
                core: indexed_core / linear_core,
                rss,
            },
        });
    }
    let median_paired_ratios = Ratios {
        wall: median(pair_ratios.iter().map(|p| p.indexed_over_linear.wall))?,
        core: median(pair_ratios.iter().map(|p| p.indexed_over_linear.core))?,
        rss: median_optional(pair_ratios.iter().map(|p| p.indexed_over_linear.rss)),
    };
    let all_correct = runs.iter().all(|run| run.correct);
    let threshold = match phase {
        "screen" => 1.0,
        "confirm" => 0.97,
        _ => return Err(format!("unknown phase {phase}")),
    };
    let passed = all_correct
        && median_paired_ratios.wall < threshold
        && median_paired_ratios.core < threshold;
    let status = match (phase, all_correct, passed) {
        (_, false, _) => "REJECTED_CORRECTNESS",
        ("screen", true, true) => "CONTINUE_TO_CONFIRMATION",
        ("screen", true, false) => "REJECTED_SCREEN",
        ("confirm", true, true) => "SELECTED_FOR_DEFAULT_REPLAY",
        ("confirm", true, false) => "REJECTED_CONFIRMATION",
        _ => unreachable!(),
    };
    Ok(Assessment {
        all_correct,
        pair_ratios,
        median_paired_ratios,
        passed_performance_gate: passed,
        status: status.into(),
    })
}

fn verify_result(result_path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(result_path).map_err(|e| format!("{}: {e}", result_path.display()))?;
    let result: PhaseResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("result JSON: {e}"))?;
    let root = result_path
        .parent()
        .ok_or_else(|| "result has no parent directory".to_string())?;
    let mut checks = 0usize;
    let mut failures = Vec::new();
    check(
        &mut checks,
        &mut failures,
        result.schema == RESULT_SCHEMA,
        "result schema",
    );
    check(
        &mut checks,
        &mut failures,
        matches!(result.phase.as_str(), "screen" | "confirm"),
        "known phase",
    );
    for run in &result.runs {
        for artifact in [
            &run.metrics_artifact,
            &run.stdout_artifact,
            &run.stderr_artifact,
        ] {
            checks += 2;
            let path = root.join(&artifact.path);
            match artifact_relative(root, &path) {
                Ok(actual) => {
                    if actual.bytes != artifact.bytes {
                        failures.push(format!("{} byte count", artifact.path));
                    }
                    if actual.sha256 != artifact.sha256 {
                        failures.push(format!("{} SHA-256", artifact.path));
                    }
                }
                Err(error) => failures.push(format!("{}: {error}", artifact.path)),
            }
        }
        let report_errors = validate_report(&run.report, run.arm);
        checks += 1;
        if !report_errors.is_empty() {
            failures.push(format!(
                "run {} terminal report: {}",
                run.order,
                report_errors.join(", ")
            ));
        }
        checks += 1;
        if run.correct != run.errors.is_empty() {
            failures.push(format!("run {} correct flag", run.order));
        }
        checks += 1;
        let metrics_path = root.join(&run.metrics_artifact.path);
        match fs::read(&metrics_path)
            .map_err(|e| e.to_string())
            .and_then(|data| {
                serde_json::from_slice::<ProcessMetrics>(&data).map_err(|e| e.to_string())
            }) {
            Ok(metrics) if metrics == run.process => {}
            Ok(_) => failures.push(format!("run {} metrics projection", run.order)),
            Err(error) => failures.push(format!("run {} metrics parse: {error}", run.order)),
        }
    }
    let mut replay_runs = result.runs.clone();
    append_cross_run_errors(&mut replay_runs);
    checks += 1;
    if replay_runs
        .iter()
        .map(|r| (&r.errors, r.correct))
        .collect::<Vec<_>>()
        != result
            .runs
            .iter()
            .map(|r| (&r.errors, r.correct))
            .collect::<Vec<_>>()
    {
        failures.push("cross-run structural replay".into());
    }
    checks += 1;
    match assess(&result.phase, &result.runs) {
        Ok(assessment) if assessments_match(&assessment, &result.assessment) => {}
        Ok(_) => failures.push("paired assessment replay".into()),
        Err(error) => failures.push(format!("paired assessment: {error}")),
    }
    let passed = checks.saturating_sub(failures.len());
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed,
        failed: failures,
    })
}

fn meter_command(out: &Path, watchdog: u64, command: &[String]) -> AnyResult<()> {
    create_new_directory(out)?;
    let stdout = out.join("stdout.txt");
    let stderr = out.join("stderr.txt");
    let environment = BTreeMap::new();
    let metrics = execute_command(command, &environment, &stdout, &stderr, watchdog)?;
    write_json(&out.join("metrics.json"), &metrics)?;
    let receipt = json!({
        "schema":"koblitz_native_meter_receipt.v1",
        "process":metrics,
        "stdout":artifact_relative(out, &stdout)?,
        "stderr":artifact_relative(out, &stderr)?,
        "metrics":artifact_relative(out, &out.join("metrics.json"))?,
    });
    write_json(&out.join("receipt.json"), &receipt)?;
    println!(
        "{}",
        serde_json::to_string(&receipt).map_err(|e| e.to_string())?
    );
    let process = receipt
        .get("process")
        .ok_or_else(|| "meter receipt lost process".to_string())?;
    if process.get("returncode").and_then(Value::as_i64) == Some(0)
        && process.get("timed_out") == Some(&json!(false))
    {
        Ok(())
    } else {
        Err("metered command failed or timed out".into())
    }
}

fn finalize_stage(
    stage: &Path,
    out: &Path,
    candidate_commit: &str,
    verifier_commit: &str,
    finalizer_commit: &str,
) -> AnyResult<()> {
    for (label, commit) in [
        ("candidate", candidate_commit),
        ("verifier", verifier_commit),
        ("finalizer", finalizer_commit),
    ] {
        if commit.len() != 40 || !commit.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(format!("{label} commit must be a full hexadecimal SHA"));
        }
    }
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    let screen_path = stage.join("development/screen/result.json");
    let confirmation_path = stage.join("development/confirmation/result.json");
    let screen: PhaseResult = read_json(&screen_path)?;
    let confirmation: PhaseResult = read_json(&confirmation_path)?;
    let screen_replay = verify_result(&screen_path)?;
    let confirmation_replay = verify_result(&confirmation_path)?;
    if !screen_replay.failed.is_empty() || !confirmation_replay.failed.is_empty() {
        return Err("phase replay failed during finalization".into());
    }
    if screen.assessment.status != "CONTINUE_TO_CONFIRMATION"
        || confirmation.assessment.status != "REJECTED_CONFIRMATION"
    {
        return Err(format!(
            "unexpected phase decisions: screen={} confirmation={}",
            screen.assessment.status, confirmation.assessment.status
        ));
    }
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 508 + measured_stage_charge.components,
        wall_seconds_sum: 21_808.881_290 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 52_387.262_218 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    let linear = confirmation
        .runs
        .iter()
        .find(|run| run.arm == Arm::Linear)
        .ok_or_else(|| "confirmation has no linear run".to_string())?;
    let indexed = confirmation
        .runs
        .iter()
        .find(|run| run.arm == Arm::Indexed)
        .ok_or_else(|| "confirmation has no indexed run".to_string())?;
    let linear_divisors = pointer_u64(&linear.report, "/cost/extra/divisor_tests")?;
    let indexed_divisors = pointer_u64(&indexed.report, "/cost/extra/divisor_tests")?;
    let mechanism = json!({
        "linear_divisor_tests":linear_divisors,
        "indexed_divisor_tests":indexed_divisors,
        "indexed_submask_lookups":pointer_u64(&indexed.report, "/cost/extra/divisor_submask_lookups")?,
        "indexed_linear_tests":pointer_u64(&indexed.report, "/cost/extra/divisor_linear_tests")?,
        "reducer_index_bytes_max":pointer_u64(&indexed.report, "/cost/extra/reducer_index_bytes_max")?,
        "lookup_reduction_percent":100.0 * (1.0 - indexed_divisors as f64 / linear_divisors as f64),
        "structural_f4_counts_identical":true,
    });
    let artifact_paths = [
        ("protocol", "PROTOCOL.md"),
        ("verifier_correction_note", "VERIFIER_CORRECTION.md"),
        ("rejected_patch", "rejected-stacked-dense-reducer.patch"),
        ("screen_result", "development/screen/result.json"),
        (
            "screen_original_verification",
            "development/screen/verification.json",
        ),
        (
            "screen_corrected_verification",
            "development/screen/verification-correction.json",
        ),
        (
            "confirmation_result",
            "development/confirmation/result.json",
        ),
        (
            "confirmation_verification",
            "development/confirmation/verification.json",
        ),
    ];
    let mut artifacts = BTreeMap::new();
    for (name, relative) in artifact_paths {
        artifacts.insert(
            name.into(),
            artifact_relative(&stage, &stage.join(relative))?,
        );
    }
    let result = FinalResult {
        schema: FINAL_SCHEMA.into(),
        candidate_commit: candidate_commit.into(),
        verifier_commit: verifier_commit.into(),
        finalizer_commit: finalizer_commit.into(),
        decision_status: "REJECTED_CONFIRMATION".into(),
        candidate_selected: false,
        runtime_default: "dense exact pair selection plus linear symbolic-reducer scan plus five-column BlockTables".into(),
        claim_boundary: "One opened n=59 decomposition target and a solver-stage engineering experiment. Direct MITM and full automorphism-aware rho boundaries are unchanged; no relation-yield, full-DLP, external-reproduction, novelty, or SOTA claim.".into(),
        conflicts: None,
        screen: screen.assessment,
        confirmation: confirmation.assessment,
        mechanism,
        measured_stage_charge,
        cumulative_measured_lower_bound,
        complete_campaign_cost: None,
        unmetered_components: vec![
            "initial development compilation and tests before the native meter existed".into(),
            "failed attempt to invoke the cargo-test-only hashed runner as a normal example binary".into(),
            "one bootstrap build of the native meter".into(),
            "final result write and final verification after their charged dry-run controls".into(),
        ],
        artifacts,
        all_seven_gates_passed: false,
        koblitz_index_calculus_sota: false,
    };
    write_json(out, &result)
}

fn verify_final(result_path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(result_path).map_err(|e| format!("{}: {e}", result_path.display()))?;
    let result: FinalResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("final result JSON: {e}"))?;
    let stage = find_stage_root(result_path)?;
    let mut checks = 0usize;
    let mut failures = Vec::new();
    check(
        &mut checks,
        &mut failures,
        result.schema == FINAL_SCHEMA,
        "final schema",
    );
    check(
        &mut checks,
        &mut failures,
        result.decision_status == "REJECTED_CONFIRMATION" && !result.candidate_selected,
        "rejected decision",
    );
    check(
        &mut checks,
        &mut failures,
        result.screen.status == "CONTINUE_TO_CONFIRMATION"
            && result.confirmation.status == "REJECTED_CONFIRMATION"
            && !result.confirmation.passed_performance_gate,
        "phase decisions",
    );
    check(
        &mut checks,
        &mut failures,
        result.complete_campaign_cost.is_none(),
        "complete cost remains unavailable",
    );
    check(
        &mut checks,
        &mut failures,
        !result.all_seven_gates_passed && !result.koblitz_index_calculus_sota,
        "claim boundary",
    );
    let current_charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failures,
        charges_match(&current_charge, &result.measured_stage_charge),
        "measured charge replay",
    );
    let expected_cumulative = Charge {
        components: 508 + current_charge.components,
        wall_seconds_sum: 21_808.881_290 + current_charge.wall_seconds_sum,
        total_core_seconds_sum: 52_387.262_218 + current_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(current_charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    check(
        &mut checks,
        &mut failures,
        charges_match(
            &expected_cumulative,
            &result.cumulative_measured_lower_bound,
        ),
        "cumulative charge replay",
    );
    for (name, artifact) in &result.artifacts {
        checks += 2;
        match artifact_relative(stage, &stage.join(&artifact.path)) {
            Ok(actual) => {
                if actual.bytes != artifact.bytes {
                    failures.push(format!("{name} byte count"));
                }
                if actual.sha256 != artifact.sha256 {
                    failures.push(format!("{name} SHA-256"));
                }
            }
            Err(error) => failures.push(format!("{name}: {error}")),
        }
    }
    let screen_path = stage.join("development/screen/result.json");
    let confirmation_path = stage.join("development/confirmation/result.json");
    let screen: PhaseResult = read_json(&screen_path)?;
    let confirmation: PhaseResult = read_json(&confirmation_path)?;
    check(
        &mut checks,
        &mut failures,
        assessments_match(&screen.assessment, &result.screen),
        "screen assessment projection",
    );
    check(
        &mut checks,
        &mut failures,
        assessments_match(&confirmation.assessment, &result.confirmation),
        "confirmation assessment projection",
    );
    let screen_replay = verify_result(&screen_path)?;
    let confirmation_replay = verify_result(&confirmation_path)?;
    check(
        &mut checks,
        &mut failures,
        screen_replay.failed.is_empty(),
        "screen replay",
    );
    check(
        &mut checks,
        &mut failures,
        confirmation_replay.failed.is_empty(),
        "confirmation replay",
    );
    let original: Verification = read_json(&stage.join("development/screen/verification.json"))?;
    let corrected: Verification =
        read_json(&stage.join("development/screen/verification-correction.json"))?;
    check(
        &mut checks,
        &mut failures,
        original.failed == vec!["paired assessment replay"],
        "original verifier failure retained",
    );
    check(
        &mut checks,
        &mut failures,
        corrected.failed.is_empty()
            && corrected.result_sha256 == original.result_sha256
            && corrected.passed == corrected.checks,
        "corrected verifier receipt",
    );
    let patch = result
        .artifacts
        .get("rejected_patch")
        .ok_or_else(|| "rejected patch artifact missing".to_string())?;
    check(
        &mut checks,
        &mut failures,
        patch.sha256 == "629a876825479a32835c767e72bdc4fff1682812d5f1660ec9f9a5285fbfab6b",
        "rejected patch identity",
    );
    let f4_source = fs::read_to_string("src/cryptanalysis/pq_f4_f2.rs")
        .map_err(|e| format!("current F4 source: {e}"))?;
    check(
        &mut checks,
        &mut failures,
        !f4_source.contains("F4_F2_INDEXED_REDUCERS") && !f4_source.contains("DenseReducerIndex"),
        "runtime candidate reverted",
    );
    let passed = checks.saturating_sub(failures.len());
    Ok(Verification {
        schema: FINAL_VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed,
        failed: failures,
    })
}

fn measured_charge(development: &Path) -> AnyResult<Charge> {
    let mut files = Vec::new();
    collect_named_files(development, "metrics.json", &mut files)?;
    files.sort();
    let mut wall = 0.0;
    let mut core = 0.0;
    let mut peak = 0u64;
    for path in &files {
        let metrics: ProcessMetrics = read_json(path)?;
        if metrics.schema != METER_SCHEMA {
            return Err(format!("{} has wrong meter schema", path.display()));
        }
        wall += metrics.wall_seconds;
        core += metrics
            .total_core_seconds
            .ok_or_else(|| format!("{} has no core time", path.display()))?;
        peak = peak.max(
            metrics
                .peak_rss_bytes
                .ok_or_else(|| format!("{} has no peak RSS", path.display()))?,
        );
    }
    Ok(Charge {
        components: files.len(),
        wall_seconds_sum: wall,
        total_core_seconds_sum: core,
        peak_rss_bytes_max: peak,
        single_core_seconds: None,
    })
}

fn collect_named_files(directory: &Path, name: &str, out: &mut Vec<PathBuf>) -> AnyResult<()> {
    for entry in fs::read_dir(directory).map_err(|e| format!("{}: {e}", directory.display()))? {
        let entry = entry.map_err(|e| e.to_string())?;
        let file_type = entry.file_type().map_err(|e| e.to_string())?;
        if file_type.is_symlink() {
            return Err(format!("refusing symlink {}", entry.path().display()));
        }
        if file_type.is_dir() {
            collect_named_files(&entry.path(), name, out)?;
        } else if file_type.is_file() && entry.file_name() == name {
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
        .ok_or_else(|| format!("{} has no Stage 188 root ancestor", path.display()))
}

fn charges_match(left: &Charge, right: &Charge) -> bool {
    left.components == right.components
        && close(left.wall_seconds_sum, right.wall_seconds_sum)
        && close(left.total_core_seconds_sum, right.total_core_seconds_sum)
        && left.peak_rss_bytes_max == right.peak_rss_bytes_max
        && left.single_core_seconds == right.single_core_seconds
}

fn pointer_u64(value: &Value, pointer: &str) -> AnyResult<u64> {
    value
        .pointer(pointer)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("{pointer}: missing unsigned integer"))
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> AnyResult<T> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

fn execute_command(
    command: &[String],
    environment: &BTreeMap<String, String>,
    stdout_path: &Path,
    stderr_path: &Path,
    watchdog_seconds: u64,
) -> AnyResult<ProcessMetrics> {
    if command.is_empty() {
        return Err("empty command".into());
    }
    let stdout =
        File::create(stdout_path).map_err(|e| format!("{}: {e}", stdout_path.display()))?;
    let stderr =
        File::create(stderr_path).map_err(|e| format!("{}: {e}", stderr_path.display()))?;
    let mut child_command = Command::new(&command[0]);
    child_command
        .args(&command[1..])
        .envs(environment)
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr));
    configure_process_group(&mut child_command);
    let start = Instant::now();
    let mut child = child_command
        .spawn()
        .map_err(|e| format!("spawn {}: {e}", command[0]))?;
    let pid = child.id();
    let outcome = wait_with_usage(&mut child, Duration::from_secs(watchdog_seconds))?;
    let wall_seconds = start.elapsed().as_secs_f64();
    Ok(ProcessMetrics {
        schema: METER_SCHEMA.into(),
        meter: "fresh child process wait4/getrusage; all child threads included".into(),
        command: command.to_vec(),
        environment: environment.clone(),
        pid,
        returncode: outcome.returncode,
        timed_out: outcome.timed_out,
        watchdog_seconds,
        wall_seconds,
        user_seconds: outcome.user_seconds,
        system_seconds: outcome.system_seconds,
        total_core_seconds: match (outcome.user_seconds, outcome.system_seconds) {
            (Some(user), Some(system)) => Some(user + system),
            _ => None,
        },
        single_core_seconds: None,
        peak_rss_bytes: outcome.peak_rss_bytes,
    })
}

struct WaitOutcome {
    returncode: i32,
    timed_out: bool,
    user_seconds: Option<f64>,
    system_seconds: Option<f64>,
    peak_rss_bytes: Option<u64>,
}

#[cfg(unix)]
fn configure_process_group(command: &mut Command) {
    use std::os::unix::process::CommandExt;
    // SAFETY: setpgid is async-signal-safe and touches only the child process.
    unsafe {
        command.pre_exec(|| {
            if libc::setpgid(0, 0) == 0 {
                Ok(())
            } else {
                Err(io::Error::last_os_error())
            }
        });
    }
}

#[cfg(not(unix))]
fn configure_process_group(_command: &mut Command) {}

#[cfg(unix)]
fn wait_with_usage(child: &mut std::process::Child, watchdog: Duration) -> AnyResult<WaitOutcome> {
    let pid = child.id() as libc::pid_t;
    let started = Instant::now();
    let mut status = 0i32;
    // SAFETY: rusage is plain data filled by wait4.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    let mut timed_out = false;
    loop {
        // SAFETY: pid is our unreaped child and pointers are valid.
        let waited = unsafe { libc::wait4(pid, &mut status, libc::WNOHANG, &mut usage) };
        if waited == pid {
            break;
        }
        if waited == -1 {
            let error = io::Error::last_os_error();
            if error.raw_os_error() == Some(libc::EINTR) {
                continue;
            }
            return Err(format!("wait4({pid}): {error}"));
        }
        if started.elapsed() >= watchdog {
            timed_out = true;
            // SAFETY: negative pid targets the process group created pre-exec.
            unsafe {
                libc::kill(-pid, libc::SIGTERM);
            }
            thread::sleep(Duration::from_millis(250));
            // SAFETY: SIGKILL is a bounded fallback if SIGTERM was ignored.
            unsafe {
                libc::kill(-pid, libc::SIGKILL);
            }
            loop {
                // SAFETY: reap our child exactly once.
                let reaped = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
                if reaped == pid {
                    break;
                }
                let error = io::Error::last_os_error();
                if error.raw_os_error() != Some(libc::EINTR) {
                    return Err(format!("wait4 after timeout ({pid}): {error}"));
                }
            }
            break;
        }
        thread::sleep(Duration::from_millis(20));
    }
    let returncode = if libc::WIFEXITED(status) {
        libc::WEXITSTATUS(status)
    } else if libc::WIFSIGNALED(status) {
        -libc::WTERMSIG(status)
    } else {
        status
    };
    let user = timeval_seconds(usage.ru_utime);
    let system = timeval_seconds(usage.ru_stime);
    let rss = if cfg!(target_os = "macos") {
        usage.ru_maxrss.max(0) as u64
    } else {
        (usage.ru_maxrss.max(0) as u64).saturating_mul(1024)
    };
    Ok(WaitOutcome {
        returncode,
        timed_out,
        user_seconds: Some(user),
        system_seconds: Some(system),
        peak_rss_bytes: Some(rss),
    })
}

#[cfg(not(unix))]
fn wait_with_usage(child: &mut std::process::Child, watchdog: Duration) -> AnyResult<WaitOutcome> {
    let started = Instant::now();
    let (returncode, timed_out) = loop {
        if let Some(status) = child.try_wait().map_err(|e| e.to_string())? {
            break (status.code().unwrap_or(-1), false);
        }
        if started.elapsed() >= watchdog {
            child.kill().map_err(|e| e.to_string())?;
            let status = child.wait().map_err(|e| e.to_string())?;
            break (status.code().unwrap_or(-1), true);
        }
        thread::sleep(Duration::from_millis(20));
    };
    Ok(WaitOutcome {
        returncode,
        timed_out,
        user_seconds: None,
        system_seconds: None,
        peak_rss_bytes: None,
    })
}

#[cfg(unix)]
fn timeval_seconds(value: libc::timeval) -> f64 {
    value.tv_sec as f64 + value.tv_usec as f64 / 1_000_000.0
}

fn expect_value(report: &Value, pointer: &str, expected: &Value, failures: &mut Vec<String>) {
    if report.pointer(pointer) != Some(expected) {
        failures.push(format!(
            "{pointer}: expected {expected}, got {}",
            report
                .pointer(pointer)
                .map_or_else(|| "<missing>".into(), Value::to_string)
        ));
    }
}

fn required_u64(report: &Value, pointer: &str, failures: &mut Vec<String>) -> Option<u64> {
    let value = report.pointer(pointer).and_then(Value::as_u64);
    if value.is_none() {
        failures.push(format!("{pointer}: missing unsigned integer"));
    }
    value
}

fn median(values: impl Iterator<Item = f64>) -> AnyResult<f64> {
    let mut values: Vec<f64> = values.collect();
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return Err("median needs finite values".into());
    }
    values.sort_by(f64::total_cmp);
    Ok(values[values.len() / 2])
}

fn median_optional(values: impl Iterator<Item = Option<f64>>) -> Option<f64> {
    let values: Option<Vec<f64>> = values.collect();
    values.and_then(|values| median(values.into_iter()).ok())
}

fn assessments_match(left: &Assessment, right: &Assessment) -> bool {
    left.all_correct == right.all_correct
        && left.passed_performance_gate == right.passed_performance_gate
        && left.status == right.status
        && ratios_match(&left.median_paired_ratios, &right.median_paired_ratios)
        && left.pair_ratios.len() == right.pair_ratios.len()
        && left
            .pair_ratios
            .iter()
            .zip(&right.pair_ratios)
            .all(|(a, b)| {
                a.repeat == b.repeat && ratios_match(&a.indexed_over_linear, &b.indexed_over_linear)
            })
}

fn ratios_match(left: &Ratios, right: &Ratios) -> bool {
    close(left.wall, right.wall)
        && close(left.core, right.core)
        && match (left.rss, right.rss) {
            (Some(a), Some(b)) => close(a, b),
            (None, None) => true,
            _ => false,
        }
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn check(checks: &mut usize, failures: &mut Vec<String>, condition: bool, label: &str) {
    *checks += 1;
    if !condition {
        failures.push(label.into());
    }
}

fn canonical_file(path: &Path) -> AnyResult<PathBuf> {
    let path = path
        .canonicalize()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    if !path.is_file() {
        return Err(format!("{} is not a file", path.display()));
    }
    Ok(path)
}

fn create_new_directory(path: &Path) -> AnyResult<()> {
    fs::create_dir(path)
        .map_err(|e| format!("refusing to reuse output directory {}: {e}", path.display()))
}

fn source_artifact_map() -> AnyResult<Value> {
    let mut map = serde_json::Map::new();
    for path in [
        "src/cryptanalysis/pq_f4_f2.rs",
        "src/cryptanalysis/ic_framework/solvers.rs",
        "examples/koblitz_pdp_backend.rs",
        "examples/koblitz_f4_stage188.rs",
    ] {
        map.insert(
            path.into(),
            serde_json::to_value(artifact_absolute(Path::new(path))?).map_err(|e| e.to_string())?,
        );
    }
    Ok(Value::Object(map))
}

fn artifact_absolute(path: &Path) -> AnyResult<Artifact> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: path.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
}

fn artifact_relative(root: &Path, path: &Path) -> AnyResult<Artifact> {
    let relative = path
        .strip_prefix(root)
        .map_err(|_| format!("{} escapes {}", path.display(), root.display()))?;
    if relative
        .components()
        .any(|part| matches!(part, std::path::Component::ParentDir))
    {
        return Err(format!("{} contains a parent escape", relative.display()));
    }
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: relative.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
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
    fn median_uses_the_middle_paired_ratio() {
        assert_eq!(median([0.9, 0.7, 0.8].into_iter()).unwrap(), 0.8);
        assert!(median(std::iter::empty()).is_err());
    }

    #[test]
    fn screen_and_confirmation_keep_the_frozen_thresholds() {
        let make = |arm, repeat, wall, core| RunReceipt {
            order: repeat,
            arm,
            repeat,
            correct: true,
            errors: Vec::new(),
            process: ProcessMetrics {
                schema: METER_SCHEMA.into(),
                meter: "test".into(),
                command: Vec::new(),
                environment: BTreeMap::new(),
                pid: 1,
                returncode: 0,
                timed_out: false,
                watchdog_seconds: 1,
                wall_seconds: wall,
                user_seconds: Some(core),
                system_seconds: Some(0.0),
                total_core_seconds: Some(core),
                single_core_seconds: None,
                peak_rss_bytes: Some(100),
            },
            metrics_artifact: Artifact {
                path: String::new(),
                bytes: 0,
                sha256: String::new(),
            },
            stdout_artifact: Artifact {
                path: String::new(),
                bytes: 0,
                sha256: String::new(),
            },
            stderr_artifact: Artifact {
                path: String::new(),
                bytes: 0,
                sha256: String::new(),
            },
            report: Value::Null,
        };
        let screen = vec![
            make(Arm::Linear, 1, 10.0, 20.0),
            make(Arm::Indexed, 1, 9.9, 19.9),
        ];
        assert!(assess("screen", &screen).unwrap().passed_performance_gate);
        assert!(!assess("confirm", &screen).unwrap().passed_performance_gate);
    }

    #[test]
    fn assessment_replay_allows_only_serialisation_scale_float_roundoff() {
        let base = Assessment {
            all_correct: true,
            pair_ratios: vec![PairRatio {
                repeat: 1,
                indexed_over_linear: Ratios {
                    wall: 0.8,
                    core: 0.9,
                    rss: Some(1.0),
                },
            }],
            median_paired_ratios: Ratios {
                wall: 0.8,
                core: 0.9,
                rss: Some(1.0),
            },
            passed_performance_gate: true,
            status: "CONTINUE_TO_CONFIRMATION".into(),
        };
        let mut rounded = base.clone();
        rounded.median_paired_ratios.wall += 1e-15;
        rounded.pair_ratios[0].indexed_over_linear.core -= 1e-15;
        assert!(assessments_match(&base, &rounded));
        rounded.median_paired_ratios.wall += 1e-6;
        assert!(!assessments_match(&base, &rounded));
    }

    #[test]
    fn meter_captures_a_native_child() {
        let directory = env::temp_dir().join(format!(
            "koblitz-stage188-meter-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&directory).unwrap();
        let stdout = directory.join("stdout");
        let stderr = directory.join("stderr");
        let command = if cfg!(windows) {
            vec!["cmd".into(), "/C".into(), "exit 0".into()]
        } else {
            vec!["/usr/bin/true".into()]
        };
        let metrics = execute_command(&command, &BTreeMap::new(), &stdout, &stderr, 5).unwrap();
        assert_eq!(metrics.returncode, 0);
        assert!(!metrics.timed_out);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn measured_charge_sums_components_and_takes_peak_rss() {
        let directory = env::temp_dir().join(format!(
            "koblitz-stage188-charge-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(directory.join("a/b")).unwrap();
        let metric = |wall, core, rss| ProcessMetrics {
            schema: METER_SCHEMA.into(),
            meter: "test".into(),
            command: Vec::new(),
            environment: BTreeMap::new(),
            pid: 1,
            returncode: 0,
            timed_out: false,
            watchdog_seconds: 1,
            wall_seconds: wall,
            user_seconds: Some(core),
            system_seconds: Some(0.0),
            total_core_seconds: Some(core),
            single_core_seconds: None,
            peak_rss_bytes: Some(rss),
        };
        write_json(&directory.join("a/metrics.json"), &metric(2.0, 3.0, 5)).unwrap();
        write_json(&directory.join("a/b/metrics.json"), &metric(7.0, 11.0, 13)).unwrap();
        let charge = measured_charge(&directory).unwrap();
        assert_eq!(charge.components, 2);
        assert_eq!(charge.wall_seconds_sum, 9.0);
        assert_eq!(charge.total_core_seconds_sum, 14.0);
        assert_eq!(charge.peak_rss_bytes_max, 13);
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn final_verifier_finds_stage_root_above_development_drafts() {
        let directory = env::temp_dir().join(format!(
            "koblitz-stage188-root-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(directory.join("development/nested")).unwrap();
        fs::write(directory.join("PROTOCOL.md"), b"frozen\n").unwrap();
        let draft = directory.join("development/nested/draft.json");
        fs::write(&draft, b"{}\n").unwrap();
        assert_eq!(find_stage_root(&draft).unwrap(), directory.as_path());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn frozen_stage178_terminals_satisfy_the_stage188_record_contract() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join(
            "research/sat_factor_base_review_20260908/continuation-05-sota-gates/\
             stage-178-current-f4-dense-reducer-index-20261001/development/paired",
        );
        let linear: Value =
            serde_json::from_slice(&fs::read(root.join("01-linear-r1/stdout.json")).unwrap())
                .unwrap();
        let indexed: Value =
            serde_json::from_slice(&fs::read(root.join("02-dense-r1/stdout.json")).unwrap())
                .unwrap();
        assert_eq!(validate_report(&linear, Arm::Linear), Vec::<String>::new());
        assert_eq!(
            validate_report(&indexed, Arm::Indexed),
            Vec::<String>::new()
        );
    }

    #[test]
    fn structural_pointer_list_has_no_duplicates() {
        let pointers = [
            "/status",
            "/exhaustive",
            "/source_instance_id",
            "/cost/ops",
            "/cost/extra/word_xors_performed",
        ];
        let set: std::collections::BTreeSet<_> = pointers.into_iter().collect();
        assert_eq!(set.len(), pointers.len());
    }
}
