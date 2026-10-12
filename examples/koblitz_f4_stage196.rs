//! Native verifier/composer for the Stage 196 full-m4ri-profile F4 experiment.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

const RESULT_SCHEMA: &str = "koblitz_stage196_phase_result.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage196_phase_verification.v1";
const FINAL_SCHEMA: &str = "koblitz_stage196_final.v1";
const FINAL_VERIFICATION_SCHEMA: &str = "koblitz_stage196_final_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";
const BIN_LABELS: [&str; 6] = [
    "lt256",
    "r256_511",
    "r512_1023",
    "r1024_2047",
    "r2048_4095",
    "ge4096",
];
const SUFFIX_THRESHOLDS: [(u64, usize); 5] = [(256, 1), (512, 2), (1024, 3), (2048, 4), (4096, 5)];

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Arm {
    Current,
    FullM4ri,
}

impl Arm {
    fn name(self) -> &'static str {
        match self {
            Self::Current => "current-profile",
            Self::FullM4ri => "full-m4ri-profile",
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
    repeat: usize,
    correct: bool,
    errors: Vec<String>,
    receipt: MeterReceipt,
    report: Value,
    bins: Vec<BinProfile>,
    residual: ResidualProfile,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct BinProfile {
    label: String,
    count: u64,
    rows_sum: u64,
    cols_sum: u64,
    eliminate_ns: u64,
    logical_xors: u64,
    performed_xors: u64,
    full_m4ri_count: u64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct ResidualProfile {
    eliminate_ns: u64,
    logical_xors: u64,
    performed_xors: u64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct SuffixAssessment {
    min_rows: u64,
    current_matrices: u64,
    full_matrices: u64,
    current_eliminate_ns: u64,
    full_eliminate_ns: u64,
    current_performed_xors: u64,
    full_performed_xors: u64,
    eliminate_ns_ratio: f64,
    performed_xors_ratio: f64,
    current_time_share: f64,
    qualifies: bool,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct ProfileDecision {
    all_correct: bool,
    recommended_min_rows: Option<u64>,
    suffixes: Vec<SuffixAssessment>,
    status: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhaseResult {
    schema: String,
    phase: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    single_core_seconds: Option<f64>,
    runs: Vec<RunRecord>,
    decision: ProfileDecision,
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
    finalizer_commit: String,
    decision_status: String,
    runtime_default: String,
    claim_boundary: String,
    conflicts: Option<u64>,
    profile: ProfileDecision,
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
    match args.first().map(String::as_str) {
        Some("compose") if args.len() == 4 => {
            let result = compose_profile(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &result)?;
            let verification = verify_result(Path::new(&args[2]))?;
            write_json(Path::new(&args[3]), &verification)?;
            println!(
                "{}",
                serde_json::to_string(&json!({
                    "schema":"koblitz_stage196_terminal.v1",
                    "phase":result.phase,
                    "status":result.decision.status,
                    "all_correct":result.decision.all_correct,
                    "recommended_min_rows":result.decision.recommended_min_rows,
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
        Some("verify") if args.len() == 3 => {
            let verification = verify_result(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("phase verification failed".into())
            }
        }
        Some("finalize") if args.len() == 5 => finalize_stage(
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
        _ => Err("usage: koblitz_f4_stage196 compose RUN_DIR RESULT VERIFICATION\n       koblitz_f4_stage196 verify RESULT VERIFICATION\n       koblitz_f4_stage196 finalize STAGE RESULT SOURCE_COMMIT FINALIZER_COMMIT\n       koblitz_f4_stage196 verify-final RESULT VERIFICATION".into()),
    }
}

fn compose_profile(directory: &Path) -> AnyResult<PhaseResult> {
    let schedule = [(Arm::Current, 1usize), (Arm::FullM4ri, 1usize)];
    let mut runs = Vec::with_capacity(schedule.len());
    for (offset, &(arm, repeat)) in schedule.iter().enumerate() {
        let order = offset + 1;
        let run_dir = directory.join(format!("{order:02}-{}-r{repeat}", arm.name()));
        let receipt_path = run_dir.join("receipt.json");
        let receipt: MeterReceipt = read_json(&receipt_path)?;
        let mut errors = validate_receipt(&run_dir, &receipt);
        errors.extend(validate_command(&receipt.process, arm));
        let report_path = run_dir.join(&receipt.stdout.path);
        let report = match read_json::<Value>(&report_path) {
            Ok(report) => {
                errors.extend(validate_report(&report, arm));
                report
            }
            Err(error) => {
                errors.push(error);
                Value::Null
            }
        };
        let (bins, residual) = match extract_profile(&report) {
            Ok(profile) => profile,
            Err(error) => {
                errors.push(error);
                empty_profile()
            }
        };
        runs.push(RunRecord {
            order,
            arm,
            repeat,
            correct: errors.is_empty(),
            errors,
            receipt,
            report,
            bins,
            residual,
        });
    }
    append_structural_errors(&mut runs);
    let decision = assess(&runs)?;
    Ok(PhaseResult {
        schema: RESULT_SCHEMA.into(),
        phase: "matrix-shape-profile".into(),
        claim_boundary: "One opened n=59 decomposition target and an opt-in elimination-shape profile only; not a selected speedup, relation yield, a full DLP, rho crossover, independent reproduction, novelty, or SOTA.".into(),
        conflicts: None,
        single_core_seconds: None,
        runs,
        decision,
    })
}

fn validate_receipt(directory: &Path, receipt: &MeterReceipt) -> Vec<String> {
    let mut failures = Vec::new();
    if receipt.schema != RECEIPT_SCHEMA {
        failures.push("meter receipt schema".into());
    }
    if receipt.process.schema != METER_SCHEMA {
        failures.push("process meter schema".into());
    }
    if receipt.process.returncode != 0 {
        failures.push(format!("return code {}", receipt.process.returncode));
    }
    if receipt.process.timed_out {
        failures.push("watchdog timeout".into());
    }
    if receipt.process.total_core_seconds.is_none() || receipt.process.peak_rss_bytes.is_none() {
        failures.push("CPU or RSS unavailable".into());
    }
    if !receipt.process.wall_seconds.is_finite()
        || receipt.process.wall_seconds < 0.0
        || receipt
            .process
            .total_core_seconds
            .is_some_and(|value| !value.is_finite() || value < 0.0)
    {
        failures.push("invalid resource accounting".into());
    }
    if receipt.stdout.path != "stdout.txt"
        || receipt.stderr.path != "stderr.txt"
        || receipt.metrics.path != "metrics.json"
    {
        failures.push("receipt artifact paths".into());
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

fn validate_command(process: &ProcessMetrics, arm: Arm) -> Vec<String> {
    let mut failures = Vec::new();
    if process.meter != "fresh child process wait4/getrusage; all child threads included"
        || !process.environment.is_empty()
        || process.watchdog_seconds != 360
        || process.single_core_seconds.is_some()
    {
        failures.push("process accounting contract".into());
    }
    let required = [
        "RAYON_NUM_THREADS=12",
        "PQ_F4_X1_BATCH=512",
        "F4_F2_DENSE_PAIR_SELECT=1",
        "F4_F2_DISABLE_INNER_BUILD_PARALLEL=1",
        "F4_F2_MATRIX_PROFILE=1",
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
    let expected_mode = match arm {
        Arm::Current => "F4_F2_FULL_M4RI=0",
        Arm::FullM4ri => "F4_F2_FULL_M4RI=1",
    };
    let modes = process
        .command
        .iter()
        .filter(|arg| arg.starts_with("F4_F2_FULL_M4RI="))
        .collect::<Vec<_>>();
    if modes.len() != 1 || modes[0].as_str() != expected_mode {
        failures.push("full-M4RI mode".into());
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
    if report
        .pointer("/cost/extra/inner_build_parallel_disabled_calls")
        .and_then(Value::as_u64)
        != Some(242)
    {
        failures.push("selected build-serial policy".into());
    }
    if report
        .pointer("/solver_matrix_profile")
        .and_then(Value::as_bool)
        != Some(true)
    {
        failures.push("matrix profile disabled".into());
    }
    let expected_bins = Value::Array(
        BIN_LABELS
            .iter()
            .map(|label| Value::String((*label).into()))
            .collect(),
    );
    if report.pointer("/solver_matrix_profile_bins") != Some(&expected_bins) {
        failures.push("matrix profile bins".into());
    }
    let full_matrices = pointer_u64(report, "/cost/extra/full_m4ri_matrices");
    match (arm, full_matrices) {
        (Arm::Current, Some(0)) => {}
        (Arm::FullM4ri, Some(value)) if value > 0 => {}
        _ => failures.push("full-M4RI routing".into()),
    }
    failures
}

fn extract_profile(report: &Value) -> AnyResult<(Vec<BinProfile>, ResidualProfile)> {
    let mut bins = Vec::with_capacity(BIN_LABELS.len());
    for label in BIN_LABELS {
        let value = |metric: &str| -> AnyResult<u64> {
            let pointer = format!("/cost/extra/matrix_profile_{label}_{metric}");
            pointer_u64(report, &pointer).ok_or_else(|| format!("missing {pointer}"))
        };
        bins.push(BinProfile {
            label: label.into(),
            count: value("count")?,
            rows_sum: value("rows_sum")?,
            cols_sum: value("cols_sum")?,
            eliminate_ns: value("eliminate_ns")?,
            logical_xors: value("logical_xors")?,
            performed_xors: value("performed_xors")?,
            full_m4ri_count: value("full_m4ri_count")?,
        });
    }
    let residual_value = |metric: &str| -> AnyResult<u64> {
        let pointer = format!("/cost/extra/matrix_profile_unbinned_{metric}");
        pointer_u64(report, &pointer).ok_or_else(|| format!("missing {pointer}"))
    };
    let residual = ResidualProfile {
        eliminate_ns: residual_value("eliminate_ns")?,
        logical_xors: residual_value("logical_xors")?,
        performed_xors: residual_value("performed_xors")?,
    };
    let count = bins.iter().map(|bin| bin.count).sum::<u64>();
    let expected_count = pointer_u64(report, "/cost/extra/f4_calls")
        .and_then(|calls| pointer_u64(report, "/cost/extra/steps").map(|steps| calls + steps))
        .ok_or_else(|| "missing F4 call/step counters".to_string())?;
    if count != expected_count {
        return Err(format!("profile matrix count {count} != {expected_count}"));
    }
    let full = bins.iter().map(|bin| bin.full_m4ri_count).sum::<u64>();
    let expected_full = pointer_u64(report, "/cost/extra/full_m4ri_matrices")
        .ok_or_else(|| "missing full-M4RI count".to_string())?;
    if full != expected_full {
        return Err(format!("profile full count {full} != {expected_full}"));
    }
    let profile_eliminate = bins.iter().map(|bin| bin.eliminate_ns).sum::<u64>();
    let profile_logical = bins.iter().map(|bin| bin.logical_xors).sum::<u64>();
    let profile_performed = bins.iter().map(|bin| bin.performed_xors).sum::<u64>();
    let eliminate = pointer_u64(report, "/cost/extra/eliminate_ns")
        .ok_or_else(|| "missing eliminate_ns".to_string())?;
    let logical = pointer_u64(report, "/cost/ops").ok_or_else(|| "missing ops".to_string())?;
    let performed = pointer_u64(report, "/cost/extra/word_xors_performed")
        .ok_or_else(|| "missing performed XORs".to_string())?;
    if profile_eliminate + residual.eliminate_ns != eliminate
        || profile_logical + residual.logical_xors != logical
        || profile_performed + residual.performed_xors != performed
    {
        return Err("profile bins and residual do not reproduce totals".into());
    }
    Ok((bins, residual))
}

fn empty_profile() -> (Vec<BinProfile>, ResidualProfile) {
    (
        BIN_LABELS
            .iter()
            .map(|label| BinProfile {
                label: (*label).into(),
                count: 0,
                rows_sum: 0,
                cols_sum: 0,
                eliminate_ns: 0,
                logical_xors: 0,
                performed_xors: 0,
                full_m4ri_count: 0,
            })
            .collect(),
        ResidualProfile {
            eliminate_ns: 0,
            logical_xors: 0,
            performed_xors: 0,
        },
    )
}

fn append_structural_errors(runs: &mut [RunRecord]) {
    let Some(reference) = runs.first().map(|run| run.report.clone()) else {
        return;
    };
    let pointers = [
        "/status",
        "/exhaustive",
        "/source_instance_id",
        "/solver_equations_blake3",
        "/solver_equations_total",
        "/solver_terms_total",
        "/solver_terms_max_per_system",
        "/fixed_x1_values",
        "/fixed_x1_masks_visited",
        "/fixed_x1_nonrational_skipped",
        "/fixed_x1_systems_constructed",
        "/fixed_x1_systems_completed",
        "/algebraic_roots",
        "/cost/degree_reached",
        "/cost/solving_degree",
        "/cost/timed_out",
        "/cost/extra/basis_len",
        "/cost/extra/extraction_tests",
        "/cost/extra/field_pairs_reduced",
        "/cost/extra/inner_build_parallel_disabled_calls",
        "/cost/extra/matrix_cols_max",
        "/cost/extra/matrix_rows_max",
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
        "/cost/extra/steps",
    ];
    for run in runs.iter_mut().skip(1) {
        for pointer in pointers {
            if reference.pointer(pointer) != run.report.pointer(pointer) {
                run.errors.push(format!("structural mismatch at {pointer}"));
            }
        }
        run.correct = run.errors.is_empty();
    }
}

fn assess(runs: &[RunRecord]) -> AnyResult<ProfileDecision> {
    let current = runs
        .iter()
        .find(|run| run.arm == Arm::Current)
        .ok_or_else(|| "profile lacks current arm".to_string())?;
    let full = runs
        .iter()
        .find(|run| run.arm == Arm::FullM4ri)
        .ok_or_else(|| "profile lacks full-M4RI arm".to_string())?;
    let all_correct = runs.iter().all(|run| run.correct);
    let total_current_ns = current.bins.iter().map(|bin| bin.eliminate_ns).sum::<u64>();
    if total_current_ns == 0 {
        return Ok(ProfileDecision {
            all_correct: false,
            recommended_min_rows: None,
            suffixes: Vec::new(),
            status: "INVALID_PROFILE".into(),
        });
    }
    let mut suffixes = Vec::with_capacity(SUFFIX_THRESHOLDS.len());
    let mut recommended = None;
    for (min_rows, start) in SUFFIX_THRESHOLDS {
        let current_matrices = current.bins[start..]
            .iter()
            .map(|bin| bin.count)
            .sum::<u64>();
        let full_matrices = full.bins[start..].iter().map(|bin| bin.count).sum::<u64>();
        let current_eliminate_ns = current.bins[start..]
            .iter()
            .map(|bin| bin.eliminate_ns)
            .sum::<u64>();
        let full_eliminate_ns = full.bins[start..]
            .iter()
            .map(|bin| bin.eliminate_ns)
            .sum::<u64>();
        let current_performed_xors = current.bins[start..]
            .iter()
            .map(|bin| bin.performed_xors)
            .sum::<u64>();
        let full_performed_xors = full.bins[start..]
            .iter()
            .map(|bin| bin.performed_xors)
            .sum::<u64>();
        let ratio = |numerator: u64, denominator: u64| {
            if denominator == 0 {
                f64::MAX
            } else {
                numerator as f64 / denominator as f64
            }
        };
        let eliminate_ns_ratio = ratio(full_eliminate_ns, current_eliminate_ns);
        let performed_xors_ratio = ratio(full_performed_xors, current_performed_xors);
        let current_time_share = current_eliminate_ns as f64 / total_current_ns as f64;
        let qualifies = all_correct
            && current_matrices > 0
            && full_matrices > 0
            && eliminate_ns_ratio < 0.90
            && performed_xors_ratio < 0.90
            && current_time_share >= 0.05;
        if qualifies {
            recommended = Some(min_rows);
        }
        suffixes.push(SuffixAssessment {
            min_rows,
            current_matrices,
            full_matrices,
            current_eliminate_ns,
            full_eliminate_ns,
            current_performed_xors,
            full_performed_xors,
            eliminate_ns_ratio,
            performed_xors_ratio,
            current_time_share,
            qualifies,
        });
    }
    let status = if !all_correct {
        "INVALID_PROFILE"
    } else if recommended.is_some() {
        "RECOMMEND_HYBRID_THRESHOLD"
    } else {
        "NO_HYBRID_THRESHOLD"
    };
    Ok(ProfileDecision {
        all_correct,
        recommended_min_rows: recommended,
        suffixes,
        status: status.into(),
    })
}

fn verify_result(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: PhaseResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("result JSON: {e}"))?;
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
        result.conflicts.is_none() && result.single_core_seconds.is_none(),
        "null conflict and single-core fields",
    );
    for run in &result.runs {
        let mut receipt_errors = validate_receipt(
            path.parent()
                .ok_or_else(|| "result has no parent".to_string())?
                .join(format!(
                    "{:02}-{}-r{}",
                    run.order,
                    run.arm.name(),
                    run.repeat
                ))
                .as_path(),
            &run.receipt,
        );
        receipt_errors.extend(validate_command(&run.receipt.process, run.arm));
        checks += 1;
        if !receipt_errors.is_empty() {
            failures.push(format!("run {} receipt replay", run.order));
        }
        let report_errors = validate_report(&run.report, run.arm);
        checks += 1;
        if !report_errors.is_empty() {
            failures.push(format!("run {} terminal replay", run.order));
        }
        checks += 1;
        match extract_profile(&run.report) {
            Ok((bins, residual)) if bins == run.bins && residual == run.residual => {}
            _ => failures.push(format!("run {} profile replay", run.order)),
        }
    }
    let mut replay = result.runs.clone();
    append_structural_errors(&mut replay);
    checks += 1;
    if replay
        .iter()
        .map(|run| (&run.errors, run.correct))
        .collect::<Vec<_>>()
        != result
            .runs
            .iter()
            .map(|run| (&run.errors, run.correct))
            .collect::<Vec<_>>()
    {
        failures.push("structural replay".into());
    }
    checks += 1;
    match assess(&result.runs) {
        Ok(decision) if decisions_match(&decision, &result.decision) => {}
        _ => failures.push("decision replay".into()),
    }
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failures.len()),
        failed: failures,
    })
}

fn finalize_stage(
    stage: &Path,
    out: &Path,
    candidate_commit: &str,
    finalizer_commit: &str,
) -> AnyResult<()> {
    for (name, commit) in [
        ("candidate", candidate_commit),
        ("finalizer", finalizer_commit),
    ] {
        if commit.len() != 40 || !commit.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(format!("{name} commit must be a full SHA"));
        }
    }
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    let profile_path = stage.join("development/profile/result.json");
    let profile: PhaseResult = read_json(&profile_path)?;
    let replay = verify_result(&profile_path)?;
    if !replay.failed.is_empty() || !valid_profile_terminal(&profile.decision) {
        return Err("profile is not a verified terminal decision".into());
    }
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 646 + measured_stage_charge.components,
        wall_seconds_sum: 24_599.183_897_831 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 64_219.173_358 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    let artifact_paths = [
        ("protocol", "PROTOCOL.md"),
        ("profile_result", "development/profile/result.json"),
        (
            "profile_verification",
            "development/profile/verification.json",
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
        finalizer_commit: finalizer_commit.into(),
        decision_status: profile.decision.status.clone(),
        runtime_default: "selected build-serial inner phases, unconditional Rayon BlockTables construction with fresh per-run word buffers, parallel row reduction, and the twelve-worker fixed-X1 outer batch".into(),
        claim_boundary: "One opened n=59 decomposition target and an opt-in elimination-shape profile. No runtime is selected here; direct MITM and full automorphism-aware rho boundaries are unchanged, with no relation-yield, full-DLP, external-reproduction, novelty, or SOTA claim.".into(),
        conflicts: None,
        profile: profile.decision,
        measured_stage_charge,
        cumulative_measured_lower_bound,
        complete_campaign_cost: None,
        unmetered_components: vec![
            "development compilation and tests before the exact native build".into(),
            "final result write and final verification after a charged composition control".into(),
        ],
        artifacts,
        all_seven_gates_passed: false,
        koblitz_index_calculus_sota: false,
    };
    write_json(out, &result)
}

fn verify_final(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: FinalResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("final result: {e}"))?;
    let stage_root = find_stage_root(path)?;
    let stage = stage_root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage_root.display()))?;
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
        result.decision_status == result.profile.status,
        "profile decision",
    );
    check(
        &mut checks,
        &mut failures,
        valid_profile_terminal(&result.profile),
        "profile projection",
    );
    check(
        &mut checks,
        &mut failures,
        result.complete_campaign_cost.is_none(),
        "complete cost unavailable",
    );
    check(
        &mut checks,
        &mut failures,
        !result.all_seven_gates_passed && !result.koblitz_index_calculus_sota,
        "claim boundary",
    );
    let charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failures,
        charges_match(&charge, &result.measured_stage_charge),
        "stage charge replay",
    );
    let cumulative = Charge {
        components: 646 + charge.components,
        wall_seconds_sum: 24_599.183_897_831 + charge.wall_seconds_sum,
        total_core_seconds_sum: 64_219.173_358 + charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(charge.peak_rss_bytes_max),
        single_core_seconds: None,
    };
    check(
        &mut checks,
        &mut failures,
        charges_match(&cumulative, &result.cumulative_measured_lower_bound),
        "cumulative charge replay",
    );
    for (name, artifact) in &result.artifacts {
        checks += 2;
        match artifact_relative(&stage, &stage.join(&artifact.path)) {
            Ok(actual) => {
                if actual.bytes != artifact.bytes {
                    failures.push(format!("{name} bytes"));
                }
                if actual.sha256 != artifact.sha256 {
                    failures.push(format!("{name} SHA-256"));
                }
            }
            Err(error) => failures.push(format!("{name}: {error}")),
        }
    }
    let profile_path = stage.join("development/profile/result.json");
    let profile: PhaseResult = read_json(&profile_path)?;
    let profile_replay = verify_result(&profile_path)?;
    check(
        &mut checks,
        &mut failures,
        profile_replay.failed.is_empty(),
        "profile replay",
    );
    check(
        &mut checks,
        &mut failures,
        decisions_match(&profile.decision, &result.profile),
        "profile decision replay",
    );
    let source = fs::read_to_string("src/cryptanalysis/pq_f4_f2.rs")
        .map_err(|e| format!("current F4 source: {e}"))?;
    check(
        &mut checks,
        &mut failures,
        source.contains("F4_F2_MATRIX_PROFILE")
            && source.contains("matrix_profile_bin")
            && source.contains("if !matrix_profile_enabled()"),
        "opt-in profile retained",
    );
    Ok(Verification {
        schema: FINAL_VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failures.len()),
        failed: failures,
    })
}

fn measured_charge(development: &Path) -> AnyResult<Charge> {
    let mut files = Vec::new();
    collect_named_files(development, "metrics.json", &mut files)?;
    files.sort();
    let mut receipt_files = Vec::new();
    collect_named_files(development, "receipt.json", &mut receipt_files)?;
    receipt_files.sort();
    let mut covered_metrics = BTreeSet::new();
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
        covered_metrics.insert(directory.join(&receipt.metrics.path));
    }
    let charged_metrics = files.iter().cloned().collect::<BTreeSet<_>>();
    if covered_metrics != charged_metrics {
        return Err("development receipts and charged metrics differ".into());
    }
    let mut wall = 0.0;
    let mut core = 0.0;
    let mut peak = 0u64;
    for file in &files {
        let metrics: ProcessMetrics = read_json(file)?;
        if metrics.schema != METER_SCHEMA {
            return Err(format!("{} has wrong schema", file.display()));
        }
        wall += metrics.wall_seconds;
        core += metrics
            .total_core_seconds
            .ok_or_else(|| format!("{} lacks core time", file.display()))?;
        peak = peak.max(
            metrics
                .peak_rss_bytes
                .ok_or_else(|| format!("{} lacks peak RSS", file.display()))?,
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

fn collect_named_files(
    directory: &Path,
    name: &str,
    out: &mut Vec<std::path::PathBuf>,
) -> AnyResult<()> {
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

fn charges_match(left: &Charge, right: &Charge) -> bool {
    left.components == right.components
        && close(left.wall_seconds_sum, right.wall_seconds_sum)
        && close(left.total_core_seconds_sum, right.total_core_seconds_sum)
        && left.peak_rss_bytes_max == right.peak_rss_bytes_max
        && left.single_core_seconds == right.single_core_seconds
}

fn find_stage_root(path: &Path) -> AnyResult<&Path> {
    path.ancestors()
        .skip(1)
        .find(|candidate| {
            candidate.join("PROTOCOL.md").is_file() && candidate.join("development").is_dir()
        })
        .ok_or_else(|| format!("{} has no Stage 196 root", path.display()))
}

fn artifact_relative(root: &Path, path: &Path) -> AnyResult<Artifact> {
    let relative = path
        .strip_prefix(root)
        .map_err(|_| format!("{} escapes {}", path.display(), root.display()))?;
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: relative.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
}

fn decisions_match(left: &ProfileDecision, right: &ProfileDecision) -> bool {
    left.all_correct == right.all_correct
        && left.recommended_min_rows == right.recommended_min_rows
        && left.status == right.status
        && left.suffixes.len() == right.suffixes.len()
        && left.suffixes.iter().zip(&right.suffixes).all(|(a, b)| {
            a.min_rows == b.min_rows
                && a.current_matrices == b.current_matrices
                && a.full_matrices == b.full_matrices
                && a.current_eliminate_ns == b.current_eliminate_ns
                && a.full_eliminate_ns == b.full_eliminate_ns
                && a.current_performed_xors == b.current_performed_xors
                && a.full_performed_xors == b.full_performed_xors
                && close(a.eliminate_ns_ratio, b.eliminate_ns_ratio)
                && close(a.performed_xors_ratio, b.performed_xors_ratio)
                && close(a.current_time_share, b.current_time_share)
                && a.qualifies == b.qualifies
        })
}

fn valid_profile_terminal(decision: &ProfileDecision) -> bool {
    match decision.status.as_str() {
        "RECOMMEND_HYBRID_THRESHOLD" => {
            decision.all_correct
                && decision.recommended_min_rows.is_some()
                && decision.suffixes.iter().any(|suffix| suffix.qualifies)
        }
        "NO_HYBRID_THRESHOLD" => {
            decision.all_correct
                && decision.recommended_min_rows.is_none()
                && decision.suffixes.iter().all(|suffix| !suffix.qualifies)
        }
        "INVALID_PROFILE" => !decision.all_correct,
        _ => false,
    }
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn check(checks: &mut usize, failures: &mut Vec<String>, condition: bool, name: &str) {
    *checks += 1;
    if !condition {
        failures.push(name.into());
    }
}

fn pointer_u64(value: &Value, pointer: &str) -> Option<u64> {
    value.pointer(pointer).and_then(Value::as_u64)
}

fn artifact_at(directory: &Path, artifact: &Artifact) -> AnyResult<Artifact> {
    let relative = Path::new(&artifact.path);
    if relative.is_absolute()
        || relative
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
    {
        return Err(format!("artifact path escapes: {}", artifact.path));
    }
    let path = directory.join(relative);
    let bytes = fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: artifact.path.clone(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
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
    fn decision_float_replay_is_tight() {
        assert!(close(1.0, 1.0 + 1e-13));
        assert!(!close(1.0, 1.0 + 1e-8));
    }

    #[test]
    fn report_validation_distinguishes_full_m4ri() {
        let report = |full: bool| {
            json!({
                "status":"unsat",
                "exhaustive":true,
                "source_instance_verified":true,
                "regenerated_source_exact":true,
                "source_instance_id":SOURCE_INSTANCE_ID,
                "solver_equations_blake3":EQUATION_FINGERPRINT,
                "fixed_x1_values":512,
                "fixed_x1_masks_visited":512,
                "fixed_x1_nonrational_skipped":270,
                "fixed_x1_systems_constructed":242,
                "fixed_x1_systems_completed":242,
                "algebraic_roots":0,
                "conflicts":null,
                "solver_matrix_profile":true,
                "solver_matrix_profile_bins":BIN_LABELS,
                "factor_base_contract":{
                    "target_subgroup_enumerated":false,
                    "discrete_log_labels_used":false
                },
                "cost":{
                    "timed_out":false,
                    "extra":{
                        "inner_build_parallel_disabled_calls":242,
                        "full_m4ri_matrices":if full { 10 } else { 0 }
                    }
                }
            })
        };
        assert!(validate_report(&report(false), Arm::Current).is_empty());
        assert!(validate_report(&report(true), Arm::FullM4ri).is_empty());
        assert!(!validate_report(&report(false), Arm::FullM4ri).is_empty());
    }

    #[test]
    fn profile_terminal_requires_consistent_recommendation() {
        let suffix = SuffixAssessment {
            min_rows: 4096,
            current_matrices: 2,
            full_matrices: 2,
            current_eliminate_ns: 100,
            full_eliminate_ns: 80,
            current_performed_xors: 100,
            full_performed_xors: 80,
            eliminate_ns_ratio: 0.8,
            performed_xors_ratio: 0.8,
            current_time_share: 0.1,
            qualifies: true,
        };
        let decision = ProfileDecision {
            all_correct: true,
            recommended_min_rows: Some(4096),
            suffixes: vec![suffix],
            status: "RECOMMEND_HYBRID_THRESHOLD".into(),
        };
        assert!(valid_profile_terminal(&decision));
    }
}
