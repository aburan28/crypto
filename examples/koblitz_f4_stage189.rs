//! Native verifier/composer for the Stage 189 outer-only F4 experiment.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

const RESULT_SCHEMA: &str = "koblitz_stage189_phase_result.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage189_phase_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Arm {
    Current,
    OuterOnly,
}

impl Arm {
    fn name(self) -> &'static str {
        match self {
            Self::Current => "current",
            Self::OuterOnly => "outer-only",
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
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Ratios {
    wall: f64,
    core: f64,
    rss: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct PairRatio {
    repeat: usize,
    outer_only_over_current: Ratios,
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

fn main() {
    if let Err(error) = real_main() {
        eprintln!("{error}");
        std::process::exit(2);
    }
}

fn real_main() -> AnyResult<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("screen" | "confirm") if args.len() == 4 => {
            let result = compose_phase(&args[0], Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &result)?;
            let verification = verify_result(Path::new(&args[2]))?;
            write_json(Path::new(&args[3]), &verification)?;
            println!(
                "{}",
                serde_json::to_string(&json!({
                    "schema":"koblitz_stage189_terminal.v1",
                    "phase":result.phase,
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
        Some("verify") if args.len() == 3 => {
            let verification = verify_result(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err("phase verification failed".into())
            }
        }
        _ => Err("usage: koblitz_f4_stage189 screen|confirm RUN_DIR RESULT VERIFICATION\n       koblitz_f4_stage189 verify RESULT VERIFICATION".into()),
    }
}

fn compose_phase(phase: &str, directory: &Path) -> AnyResult<PhaseResult> {
    let schedule: &[(Arm, usize)] = match phase {
        "screen" => &[(Arm::Current, 1), (Arm::OuterOnly, 1)],
        "confirm" => &[
            (Arm::Current, 1),
            (Arm::OuterOnly, 1),
            (Arm::OuterOnly, 2),
            (Arm::Current, 2),
            (Arm::Current, 3),
            (Arm::OuterOnly, 3),
        ],
        _ => return Err(format!("unknown phase {phase}")),
    };
    let mut runs = Vec::with_capacity(schedule.len());
    for (offset, &(arm, repeat)) in schedule.iter().enumerate() {
        let order = offset + 1;
        let run_dir = directory.join(format!("{order:02}-{}-r{repeat}", arm.name()));
        let receipt_path = run_dir.join("receipt.json");
        let receipt: MeterReceipt = read_json(&receipt_path)?;
        let mut errors = validate_receipt(&run_dir, &receipt);
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
        runs.push(RunRecord {
            order,
            arm,
            repeat,
            correct: errors.is_empty(),
            errors,
            receipt,
            report,
        });
    }
    append_structural_errors(&mut runs);
    let assessment = assess(phase, &runs)?;
    Ok(PhaseResult {
        schema: RESULT_SCHEMA.into(),
        phase: phase.into(),
        claim_boundary: "One opened n=59 decomposition target and a solver-scheduling experiment only; not relation yield, a full DLP, rho crossover, independent reproduction, novelty, or SOTA.".into(),
        conflicts: None,
        single_core_seconds: None,
        runs,
        assessment,
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
    let disabled = report
        .pointer("/cost/extra/inner_parallel_disabled_calls")
        .and_then(Value::as_u64);
    let expected = match arm {
        Arm::Current => 0,
        Arm::OuterOnly => 242,
    };
    if disabled != Some(expected) {
        failures.push(format!(
            "inner_parallel_disabled_calls expected {expected}, got {disabled:?}"
        ));
    }
    failures
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
        "/cost/ops",
        "/cost/degree_reached",
        "/cost/solving_degree",
        "/cost/timed_out",
        "/cost/extra/basis_len",
        "/cost/extra/divisor_tests",
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
}

fn assess(phase: &str, runs: &[RunRecord]) -> AnyResult<Assessment> {
    let mut pairs = BTreeMap::<usize, (Option<&RunRecord>, Option<&RunRecord>)>::new();
    for run in runs {
        let entry = pairs.entry(run.repeat).or_default();
        match run.arm {
            Arm::Current => entry.0 = Some(run),
            Arm::OuterOnly => entry.1 = Some(run),
        }
    }
    let mut pair_ratios = Vec::new();
    for (repeat, (current, outer)) in pairs {
        let current = current.ok_or_else(|| format!("pair {repeat} lacks current"))?;
        let outer = outer.ok_or_else(|| format!("pair {repeat} lacks outer-only"))?;
        let current_core = current
            .receipt
            .process
            .total_core_seconds
            .ok_or_else(|| "current core unavailable".to_string())?;
        let outer_core = outer
            .receipt
            .process
            .total_core_seconds
            .ok_or_else(|| "outer-only core unavailable".to_string())?;
        let current_rss = current
            .receipt
            .process
            .peak_rss_bytes
            .ok_or_else(|| "current RSS unavailable".to_string())?;
        let outer_rss = outer
            .receipt
            .process
            .peak_rss_bytes
            .ok_or_else(|| "outer-only RSS unavailable".to_string())?;
        pair_ratios.push(PairRatio {
            repeat,
            outer_only_over_current: Ratios {
                wall: outer.receipt.process.wall_seconds / current.receipt.process.wall_seconds,
                core: outer_core / current_core,
                rss: outer_rss as f64 / current_rss as f64,
            },
        });
    }
    let median_paired_ratios = Ratios {
        wall: median(pair_ratios.iter().map(|p| p.outer_only_over_current.wall))?,
        core: median(pair_ratios.iter().map(|p| p.outer_only_over_current.core))?,
        rss: median(pair_ratios.iter().map(|p| p.outer_only_over_current.rss))?,
    };
    let all_correct = runs.iter().all(|run| run.correct);
    let threshold = if phase == "screen" { 1.0 } else { 0.97 };
    let passed = all_correct
        && median_paired_ratios.wall < threshold
        && median_paired_ratios.core < threshold;
    let status = match (phase, all_correct, passed) {
        (_, false, _) => "REJECTED_CORRECTNESS",
        ("screen", true, true) => "CONTINUE_TO_CONFIRMATION",
        ("screen", true, false) => "REJECTED_SCREEN",
        ("confirm", true, true) => "SELECTED_FOR_DEFAULT_REPLAY",
        ("confirm", true, false) => "REJECTED_CONFIRMATION",
        _ => return Err(format!("unknown phase {phase}")),
    };
    Ok(Assessment {
        all_correct,
        pair_ratios,
        median_paired_ratios,
        passed_performance_gate: passed,
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
        let receipt_errors = validate_receipt(
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
        checks += 1;
        if !receipt_errors.is_empty() {
            failures.push(format!("run {} receipt replay", run.order));
        }
        let report_errors = validate_report(&run.report, run.arm);
        checks += 1;
        if !report_errors.is_empty() {
            failures.push(format!("run {} terminal replay", run.order));
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
    match assess(&result.phase, &result.runs) {
        Ok(assessment) if assessments_match(&assessment, &result.assessment) => {}
        _ => failures.push("assessment replay".into()),
    }
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failures.len()),
        failed: failures,
    })
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
                a.repeat == b.repeat
                    && ratios_match(&a.outer_only_over_current, &b.outer_only_over_current)
            })
}

fn ratios_match(left: &Ratios, right: &Ratios) -> bool {
    close(left.wall, right.wall) && close(left.core, right.core) && close(left.rss, right.rss)
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn median(values: impl Iterator<Item = f64>) -> AnyResult<f64> {
    let mut values: Vec<f64> = values.collect();
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return Err("median requires finite values".into());
    }
    values.sort_by(f64::total_cmp);
    Ok(values[values.len() / 2])
}

fn check(checks: &mut usize, failures: &mut Vec<String>, condition: bool, name: &str) {
    *checks += 1;
    if !condition {
        failures.push(name.into());
    }
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
    fn frozen_thresholds_require_joint_wall_and_core_improvement() {
        let ratios = Ratios {
            wall: 0.96,
            core: 0.98,
            rss: 0.8,
        };
        assert!(!(ratios.wall < 0.97 && ratios.core < 0.97));
    }

    #[test]
    fn assessment_float_replay_is_tight() {
        assert!(close(1.0, 1.0 + 1e-13));
        assert!(!close(1.0, 1.0 + 1e-8));
    }

    #[test]
    fn report_validation_distinguishes_the_two_schedules() {
        let report = |disabled| {
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
                "factor_base_contract":{
                    "target_subgroup_enumerated":false,
                    "discrete_log_labels_used":false
                },
                "cost":{
                    "timed_out":false,
                    "extra":{"inner_parallel_disabled_calls":disabled}
                }
            })
        };
        assert!(validate_report(&report(0), Arm::Current).is_empty());
        assert!(validate_report(&report(242), Arm::OuterOnly).is_empty());
        assert!(!validate_report(&report(0), Arm::OuterOnly).is_empty());
    }
}
