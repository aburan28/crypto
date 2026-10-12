//! Native composer/verifier for the Stage 191 single-core boundary replay.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

const RESULT_SCHEMA: &str = "koblitz_stage191_boundary_result.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage191_boundary_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const EQUATION_FINGERPRINT: &str =
    "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb";

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

#[derive(Clone, Debug, Deserialize, Serialize)]
struct MeterReceipt {
    schema: String,
    process: ProcessMetrics,
    stdout: Artifact,
    stderr: Artifact,
    metrics: Artifact,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct ResourceSummary {
    wall_seconds: f64,
    total_core_seconds: f64,
    single_core_seconds: f64,
    peak_rss_bytes: f64,
    conflicts: Option<u64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Ratios {
    wall: f64,
    total_core: f64,
    single_core: f64,
    peak_rss: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Charge {
    components: usize,
    wall_seconds_sum: f64,
    total_core_seconds_sum: f64,
    peak_rss_bytes_max: u64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct BoundaryResult {
    schema: String,
    source_commit: String,
    finalizer_commit: String,
    claim_boundary: String,
    direct_before: MeterReceipt,
    direct_before_report: Value,
    f4_single_core: MeterReceipt,
    f4_report: Value,
    direct_after: MeterReceipt,
    direct_after_report: Value,
    direct_reference_mean: ResourceSummary,
    f4_resources: ResourceSummary,
    f4_over_direct: Ratios,
    structural_checks_passed: bool,
    measured_stage_charge: Charge,
    cumulative_measured_lower_bound: Charge,
    complete_campaign_cost: Option<Charge>,
    artifacts: BTreeMap<String, Artifact>,
    all_seven_gates_passed: bool,
    koblitz_index_calculus_sota: bool,
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
        Some("compose") if args.len() == 5 => compose(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &args[3],
            &args[4],
        ),
        Some("verify") if args.len() == 3 => {
            let verification = verify(Path::new(&args[1]))?;
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
        _ => Err("usage: koblitz_f4_stage191 compose STAGE RESULT SOURCE_COMMIT FINALIZER_COMMIT\n       koblitz_f4_stage191 verify RESULT VERIFICATION".into()),
    }
}

fn compose(stage: &Path, out: &Path, source_commit: &str, finalizer_commit: &str) -> AnyResult<()> {
    validate_commit("source", source_commit)?;
    validate_commit("finalizer", finalizer_commit)?;
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    let before_dir = stage.join("development/01-direct-before");
    let f4_dir = stage.join("development/02-f4-single-core");
    let after_dir = stage.join("development/03-direct-after");
    let before = load_run(&before_dir)?;
    let f4 = load_run(&f4_dir)?;
    let after = load_run(&after_dir)?;
    let mut errors = Vec::new();
    errors.extend(validate_direct(&before.1, &before.0, "before"));
    errors.extend(validate_f4(&f4.1, &f4.0));
    errors.extend(validate_direct(&after.1, &after.0, "after"));
    for pointer in direct_identity_pointers() {
        if before.1.pointer(pointer) != after.1.pointer(pointer) {
            errors.push(format!("direct controls differ at {pointer}"));
        }
    }
    if !errors.is_empty() {
        return Err(format!("boundary records failed: {}", errors.join("; ")));
    }
    let before_resources = resources(&before.0)?;
    let after_resources = resources(&after.0)?;
    let f4_resources = resources(&f4.0)?;
    let direct_reference_mean = ResourceSummary {
        wall_seconds: mean(before_resources.wall_seconds, after_resources.wall_seconds),
        total_core_seconds: mean(
            before_resources.total_core_seconds,
            after_resources.total_core_seconds,
        ),
        single_core_seconds: mean(
            before_resources.single_core_seconds,
            after_resources.single_core_seconds,
        ),
        peak_rss_bytes: mean(
            before_resources.peak_rss_bytes,
            after_resources.peak_rss_bytes,
        ),
        conflicts: None,
    };
    let f4_over_direct = Ratios {
        wall: f4_resources.wall_seconds / direct_reference_mean.wall_seconds,
        total_core: f4_resources.total_core_seconds / direct_reference_mean.total_core_seconds,
        single_core: f4_resources.single_core_seconds / direct_reference_mean.single_core_seconds,
        peak_rss: f4_resources.peak_rss_bytes / direct_reference_mean.peak_rss_bytes,
    };
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 564 + measured_stage_charge.components,
        wall_seconds_sum: 23_476.882_993_374 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 59_943.278_38 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
    };
    let mut artifacts = BTreeMap::new();
    for (name, relative) in [
        ("protocol", "PROTOCOL.md"),
        (
            "direct_before_receipt",
            "development/01-direct-before/receipt.json",
        ),
        (
            "direct_before_stdout",
            "development/01-direct-before/stdout.txt",
        ),
        ("f4_receipt", "development/02-f4-single-core/receipt.json"),
        ("f4_stdout", "development/02-f4-single-core/stdout.txt"),
        (
            "direct_after_receipt",
            "development/03-direct-after/receipt.json",
        ),
        (
            "direct_after_stdout",
            "development/03-direct-after/stdout.txt",
        ),
    ] {
        artifacts.insert(
            name.into(),
            artifact_relative(&stage, &stage.join(relative))?,
        );
    }
    write_json(
        out,
        &BoundaryResult {
            schema: RESULT_SCHEMA.into(),
            source_commit: source_commit.into(),
            finalizer_commit: finalizer_commit.into(),
            claim_boundary: "One opened n=59 decomposition target. Single-core F4 versus same-binary direct MITM only; relation collection, linear algebra, recovery, full rho, external reproduction, novelty and SOTA remain outside this ratio.".into(),
            direct_before: before.0,
            direct_before_report: before.1,
            f4_single_core: f4.0,
            f4_report: f4.1,
            direct_after: after.0,
            direct_after_report: after.1,
            direct_reference_mean,
            f4_resources,
            f4_over_direct,
            structural_checks_passed: true,
            measured_stage_charge,
            cumulative_measured_lower_bound,
            complete_campaign_cost: None,
            artifacts,
            all_seven_gates_passed: false,
            koblitz_index_calculus_sota: false,
        },
    )
}

fn verify(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: BoundaryResult =
        serde_json::from_slice(&bytes).map_err(|e| format!("result JSON: {e}"))?;
    let stage = find_stage_root(path)?;
    let mut checks = 0usize;
    let mut failures = Vec::new();
    check(
        &mut checks,
        &mut failures,
        result.schema == RESULT_SCHEMA,
        "result schema",
    );
    let before_dir = stage.join("development/01-direct-before");
    let f4_dir = stage.join("development/02-f4-single-core");
    let after_dir = stage.join("development/03-direct-after");
    let before = load_run(&before_dir)?;
    let f4 = load_run(&f4_dir)?;
    let after = load_run(&after_dir)?;
    check(
        &mut checks,
        &mut failures,
        before.0.process == result.direct_before.process && before.1 == result.direct_before_report,
        "direct-before projection",
    );
    check(
        &mut checks,
        &mut failures,
        f4.0.process == result.f4_single_core.process && f4.1 == result.f4_report,
        "F4 projection",
    );
    check(
        &mut checks,
        &mut failures,
        after.0.process == result.direct_after.process && after.1 == result.direct_after_report,
        "direct-after projection",
    );
    check(
        &mut checks,
        &mut failures,
        validate_direct(&before.1, &before.0, "before").is_empty()
            && validate_direct(&after.1, &after.0, "after").is_empty()
            && validate_f4(&f4.1, &f4.0).is_empty(),
        "terminal and command contracts",
    );
    check(
        &mut checks,
        &mut failures,
        direct_identity_pointers()
            .iter()
            .all(|pointer| before.1.pointer(pointer) == after.1.pointer(pointer)),
        "direct-control identity",
    );
    let before_resources = resources(&before.0)?;
    let after_resources = resources(&after.0)?;
    let f4_resources = resources(&f4.0)?;
    let reference = ResourceSummary {
        wall_seconds: mean(before_resources.wall_seconds, after_resources.wall_seconds),
        total_core_seconds: mean(
            before_resources.total_core_seconds,
            after_resources.total_core_seconds,
        ),
        single_core_seconds: mean(
            before_resources.single_core_seconds,
            after_resources.single_core_seconds,
        ),
        peak_rss_bytes: mean(
            before_resources.peak_rss_bytes,
            after_resources.peak_rss_bytes,
        ),
        conflicts: None,
    };
    let ratios = Ratios {
        wall: f4_resources.wall_seconds / reference.wall_seconds,
        total_core: f4_resources.total_core_seconds / reference.total_core_seconds,
        single_core: f4_resources.single_core_seconds / reference.single_core_seconds,
        peak_rss: f4_resources.peak_rss_bytes / reference.peak_rss_bytes,
    };
    check(
        &mut checks,
        &mut failures,
        summary_matches(&reference, &result.direct_reference_mean)
            && summary_matches(&f4_resources, &result.f4_resources)
            && ratios_match(&ratios, &result.f4_over_direct),
        "resource and ratio replay",
    );
    let charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failures,
        charges_match(&charge, &result.measured_stage_charge),
        "stage charge replay",
    );
    let cumulative = Charge {
        components: 564 + charge.components,
        wall_seconds_sum: 23_476.882_993_374 + charge.wall_seconds_sum,
        total_core_seconds_sum: 59_943.278_38 + charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(charge.peak_rss_bytes_max),
    };
    check(
        &mut checks,
        &mut failures,
        charges_match(&cumulative, &result.cumulative_measured_lower_bound),
        "cumulative charge replay",
    );
    for (name, artifact) in &result.artifacts {
        checks += 2;
        match artifact_relative(stage, &stage.join(&artifact.path)) {
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
    check(
        &mut checks,
        &mut failures,
        result.complete_campaign_cost.is_none()
            && result.structural_checks_passed
            && !result.all_seven_gates_passed
            && !result.koblitz_index_calculus_sota,
        "claim boundary",
    );
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failures.len()),
        failed: failures,
    })
}

fn load_run(directory: &Path) -> AnyResult<(MeterReceipt, Value)> {
    let receipt: MeterReceipt = read_json(&directory.join("receipt.json"))?;
    let errors = validate_receipt(directory, &receipt);
    if !errors.is_empty() {
        return Err(format!("{}: {}", directory.display(), errors.join("; ")));
    }
    let report = read_json(&directory.join(&receipt.stdout.path))?;
    Ok((receipt, report))
}

fn validate_receipt(directory: &Path, receipt: &MeterReceipt) -> Vec<String> {
    let mut errors = Vec::new();
    if receipt.schema != RECEIPT_SCHEMA || receipt.process.schema != METER_SCHEMA {
        errors.push("meter schema".into());
    }
    if receipt.process.returncode != 0 || receipt.process.timed_out {
        errors.push("process failed or timed out".into());
    }
    if receipt.process.total_core_seconds.is_none() || receipt.process.peak_rss_bytes.is_none() {
        errors.push("CPU or RSS unavailable".into());
    }
    for artifact in [&receipt.stdout, &receipt.stderr, &receipt.metrics] {
        match artifact_at(directory, artifact) {
            Ok(actual) if actual == *artifact => {}
            Ok(_) => errors.push(format!("artifact mismatch: {}", artifact.path)),
            Err(error) => errors.push(error),
        }
    }
    errors
}

fn validate_direct(report: &Value, receipt: &MeterReceipt, label: &str) -> Vec<String> {
    let mut errors = Vec::new();
    for (pointer, expected) in [
        ("/schema", json!("koblitz_pdp_isolated_backend.v1")),
        ("/backend", json!("direct-mitm")),
        ("/status", json!("unsat")),
        ("/exhaustive", json!(true)),
        ("/source_instance_verified", json!(true)),
        ("/regenerated_source_exact", json!(true)),
        ("/source_instance_id", json!(SOURCE_INSTANCE_ID)),
    ] {
        if report.pointer(pointer) != Some(&expected) {
            errors.push(format!("{label} {pointer}"));
        }
    }
    if receipt.process.command.iter().any(|arg| {
        arg.starts_with("F4_F2_DISABLE_INNER_BUILD_PARALLEL=")
            || arg.starts_with("RAYON_NUM_THREADS=")
    }) {
        errors.push(format!("{label} direct command has an F4/thread override"));
    }
    if !receipt
        .process
        .command
        .iter()
        .any(|arg| arg == "direct-mitm")
    {
        errors.push(format!("{label} is not direct-mitm"));
    }
    for expected in [
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
    ] {
        if !receipt.process.command.iter().any(|arg| arg == expected) {
            errors.push(format!("{label} direct command lacks {expected}"));
        }
    }
    if report
        .pointer("/conflicts")
        .is_some_and(|value| !value.is_null())
    {
        errors.push(format!(
            "{label} direct conflicts must be null or unavailable"
        ));
    }
    errors
}

fn validate_f4(report: &Value, receipt: &MeterReceipt) -> Vec<String> {
    let mut errors = Vec::new();
    for (pointer, expected) in [
        ("/backend", json!("native-f4")),
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
        ("/single_thread_requested", json!(true)),
        ("/conflicts", Value::Null),
        (
            "/factor_base_contract/target_subgroup_enumerated",
            json!(false),
        ),
        (
            "/factor_base_contract/discrete_log_labels_used",
            json!(false),
        ),
        (
            "/cost/extra/inner_build_parallel_disabled_calls",
            json!(242),
        ),
    ] {
        if report.pointer(pointer) != Some(&expected) {
            errors.push(format!("F4 {pointer}"));
        }
    }
    for expected in [
        "RAYON_NUM_THREADS=1",
        "PQ_F4_X1_BATCH=1",
        "VECLIB_MAXIMUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "BLIS_NUM_THREADS=1",
        "NUMEXPR_NUM_THREADS=1",
        "F4_F2_DENSE_PAIR_SELECT=1",
        "F4_F2_FULL_M4RI=0",
    ] {
        if !receipt.process.command.iter().any(|arg| arg == expected) {
            errors.push(format!("F4 command lacks {expected}"));
        }
    }
    if receipt
        .process
        .command
        .iter()
        .any(|arg| arg.starts_with("F4_F2_DISABLE_INNER_BUILD_PARALLEL="))
    {
        errors.push("F4 command overrides selected build scheduling".into());
    }
    if !receipt.process.command.iter().any(|arg| arg == "native-f4") {
        errors.push("F4 command does not select native-f4".into());
    }
    errors
}

fn direct_identity_pointers() -> &'static [&'static str] {
    &[
        "/schema",
        "/backend",
        "/status",
        "/exhaustive",
        "/source_instance_id",
        "/source_instance_verified",
        "/regenerated_source_exact",
        "/factor_points",
        "/pair_entries",
        "/group_additions",
        "/witness_indices",
        "/source_witness_valid",
    ]
}

fn resources(receipt: &MeterReceipt) -> AnyResult<ResourceSummary> {
    let core = receipt
        .process
        .total_core_seconds
        .ok_or_else(|| "core time unavailable".to_string())?;
    Ok(ResourceSummary {
        wall_seconds: receipt.process.wall_seconds,
        total_core_seconds: core,
        single_core_seconds: core,
        peak_rss_bytes: receipt
            .process
            .peak_rss_bytes
            .ok_or_else(|| "peak RSS unavailable".to_string())? as f64,
        conflicts: None,
    })
}

fn measured_charge(development: &Path) -> AnyResult<Charge> {
    let mut files = Vec::new();
    collect_named_files(development, "metrics.json", &mut files)?;
    files.sort();
    let mut charge = Charge {
        components: files.len(),
        wall_seconds_sum: 0.0,
        total_core_seconds_sum: 0.0,
        peak_rss_bytes_max: 0,
    };
    for file in files {
        let metrics: ProcessMetrics = read_json(&file)?;
        if metrics.schema != METER_SCHEMA {
            return Err(format!("{} has wrong meter schema", file.display()));
        }
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

fn artifact_at(directory: &Path, artifact: &Artifact) -> AnyResult<Artifact> {
    let relative = Path::new(&artifact.path);
    if relative.is_absolute()
        || relative
            .components()
            .any(|part| matches!(part, std::path::Component::ParentDir))
    {
        return Err(format!("artifact path escapes: {}", artifact.path));
    }
    let bytes = fs::read(directory.join(relative)).map_err(|e| e.to_string())?;
    Ok(Artifact {
        path: artifact.path.clone(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
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

fn find_stage_root(path: &Path) -> AnyResult<&Path> {
    path.ancestors()
        .skip(1)
        .find(|candidate| {
            candidate.join("PROTOCOL.md").is_file() && candidate.join("development").is_dir()
        })
        .ok_or_else(|| format!("{} has no Stage 191 root", path.display()))
}

fn validate_commit(name: &str, commit: &str) -> AnyResult<()> {
    if commit.len() == 40 && commit.bytes().all(|b| b.is_ascii_hexdigit()) {
        Ok(())
    } else {
        Err(format!("{name} commit must be a full SHA"))
    }
}

fn mean(left: f64, right: f64) -> f64 {
    (left + right) / 2.0
}

fn summary_matches(left: &ResourceSummary, right: &ResourceSummary) -> bool {
    close(left.wall_seconds, right.wall_seconds)
        && close(left.total_core_seconds, right.total_core_seconds)
        && close(left.single_core_seconds, right.single_core_seconds)
        && close(left.peak_rss_bytes, right.peak_rss_bytes)
        && left.conflicts == right.conflicts
}

fn ratios_match(left: &Ratios, right: &Ratios) -> bool {
    close(left.wall, right.wall)
        && close(left.total_core, right.total_core)
        && close(left.single_core, right.single_core)
        && close(left.peak_rss, right.peak_rss)
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn charges_match(left: &Charge, right: &Charge) -> bool {
    left.components == right.components
        && close(left.wall_seconds_sum, right.wall_seconds_sum)
        && close(left.total_core_seconds_sum, right.total_core_seconds_sum)
        && left.peak_rss_bytes_max == right.peak_rss_bytes_max
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
    fn bracketing_reference_is_the_arithmetic_mean() {
        assert_eq!(mean(2.0, 4.0), 3.0);
    }

    #[test]
    fn ratios_are_compared_with_serialisation_scale_tolerance() {
        assert!(close(1.0, 1.0 + 1e-13));
        assert!(!close(1.0, 1.0 + 1e-8));
    }

    #[test]
    fn direct_identity_covers_counted_search_structure() {
        let pointers = direct_identity_pointers();
        assert!(pointers.contains(&"/factor_points"));
        assert!(pointers.contains(&"/pair_entries"));
        assert!(pointers.contains(&"/group_additions"));
    }
}
