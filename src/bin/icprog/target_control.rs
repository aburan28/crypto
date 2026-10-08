//! Original-capsule controller and data-only runtime audit for the bounded n17 worker.
use super::ordinary_control::capsule as preparation;
use super::sat_control::native;
use super::{
    f5_source_custody, ordinary_control, ordinary_preparation, target_build, target_custody,
    target_math,
};
#[path = "../prepared_target_worker/contract.rs"]
pub(super) mod capsule;
#[path = "../prepared_target_worker/journal.rs"]
#[allow(dead_code)] // The writer is shared with the worker; this controller inspects it.
pub(super) mod journal;
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path, time::Instant};

fn load(path: &Path) -> Result<Value, String> {
    target_math::parse(&read(path, 16 * 1024 * 1024)?)
}
/// Verify original preparation bytes/build/math without launching any executable.
pub(super) fn original_preparation(binding: &capsule::PreparationBinding) -> Result<Value, String> {
    let (root, execution, _) = capsule::preparation_paths(binding)?;
    let (record, seal) = preparation::check_capsule(&root)?;
    require(
        seal == binding.registration_sha256,
        "original preparation seal differs",
    )?;
    ordinary_control::verify_build_from(&root, &record, &root)?;
    let bytes = read(&execution.join("producer.json"), 16 * 1024 * 1024)?;
    let producer = target_math::parse_ordinary_producer(&bytes)?;
    let math = producer["mathematical_input"].clone();
    let checked = ordinary_preparation::audit(&math)?;
    require(
        sha256(&bytes) == binding.producer_sha256
            && canonical_sha(&math)? == binding.mathematics_sha256
            && checked["declared_family"] == "matrix_f5"
            && checked["rank"] == 29
            && checked["folded_columns"] == 29
            && checked["usable_points"] == 62
            && checked["geometric_points"] == 63
            && checked["panel_complete"] == true
            && checked["mathematical_preparation_complete"] == true,
        "original preparation lacks exact complete n17 F5 mathematics",
    )?;
    capsule::preparation_paths(binding)?;
    require(
        read(&execution.join("producer.json"), 16 * 1024 * 1024)? == bytes,
        "original preparation producer changed during verification",
    )?;
    Ok(math)
}
fn admissible(record: &capsule::Registration, expected: &str) -> Result<(), String> {
    require(
        journal::digest(expected) && !record.validation_only,
        "target registration is validation-only or external seal is invalid",
    )
}
fn frozen_self(root: &Path, record: &capsule::Registration) -> Result<String, String> {
    capsule::require_identity(&capsule::identity(), &record.worker_build_identity)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let hash = sha256(&read(&own, 128 * 1024 * 1024)?);
    require(
        own.canonicalize().map_err(|e| e.to_string())?
            == root
                .join("immutable/bin/icprog")
                .canonicalize()
                .map_err(|e| e.to_string())?
            && hash == record.auditor_sha256,
        "target action requires original frozen controller",
    )?;
    Ok(hash)
}
fn worker_args(root: &Path, execution: &Path, seal: &str) -> Vec<String> {
    vec![
        "run".into(),
        "--capsule".into(),
        root.to_string_lossy().into_owned(),
        "--execution".into(),
        execution.to_string_lossy().into_owned(),
        "--registration-sha256".into(),
        seal.into(),
    ]
}
fn preparation_args(root: &Path, execution: &Path, seal: &str, out: &Path) -> Vec<String> {
    vec![
        "ordinary-control-audit".into(),
        "--capsule".into(),
        root.to_string_lossy().into_owned(),
        "--execution".into(),
        execution.to_string_lossy().into_owned(),
        "--registration-sha256".into(),
        seal.into(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
    ]
}
fn exposure(cfg: &capsule::Config, record: &capsule::Registration, seal: &str) -> Value {
    json!({"schema_version":1,"scope":capsule::SCOPE,"registration_sha256":seal,
        "config_sha256":record.config_sha256,"target":cfg.target,"target_count":1,
        "input_kind":"supplied-public-point","fixture_scalar_generated":false,
        "exposure_stage":"before-sole-worker-launch; point already present in registration",
        "fresh_target_certified":false})
}
fn consume(
    root: &Path,
    execution: &Path,
    record: &capsule::Registration,
    cfg: &capsule::Config,
    seal: &str,
) -> Result<(), String> {
    admissible(record, seal)?;
    save(
        &root.join("consumed.json"),
        &json!(capsule::Claim {
            schema_version: 1,
            scope: capsule::SCOPE.into(),
            registration_sha256: seal.into(),
            execution: execution.into(),
            worker_sha256: record.worker_sha256.clone()
        }),
    )?;
    // A failure after consumption is terminal for this registration, even before spawn.
    save(
        &execution.join("target-exposure.json"),
        &exposure(cfg, record, seal),
    )
}
fn output_outside(execution: &Path, root: &Path, out: &Path) -> Result<(), String> {
    let parent = out
        .parent()
        .ok_or("target output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !parent.starts_with(execution) && !parent.starts_with(root),
        "target audit/inspection output must be outside original evidence",
    )
}

/// Consume before the sole launch; preserve every transport/source/drain outcome.
pub(super) fn execute(
    root: &Path,
    publication: &Path,
    execution: &Path,
    expected: &str,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let record = capsule::check_capsule(&root, expected)?;
    admissible(&record, expected)?;
    let own = frozen_self(&root, &record)?;
    target_build::verify(&root, &record)?;
    original_preparation(&record.preparation)?;
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "target registration already consumed; never retry",
    )?;
    let publication_preflight = target_custody::scientific_preflight(publication, &root, expected)?;
    let adopted = root
        .join("adoption-start.json")
        .try_exists()
        .map_err(|e| e.to_string())?;
    let adoption_preflight = if adopted {
        Some(f5_source_custody::verify_adoption(&root, expected)?)
    } else {
        None
    };
    let cfg = capsule::config(&root)?;
    fs::create_dir(execution).map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let preflight_path = execution.join("publication-preflight.json");
    save(&preflight_path, &publication_preflight)?;
    let preflight_sha256 = sha256(&read(&preflight_path, 16 * 1024 * 1024)?);
    if let Some(adoption) = &adoption_preflight {
        save(&execution.join("adoption-preflight.json"), adoption)?;
    }
    let adoption_preflight_sha256 = adoption_preflight
        .as_ref()
        .map(|_| read(&execution.join("adoption-preflight.json"), 65536).map(|b| sha256(&b)))
        .transpose()?;
    consume(&root, &execution, &record, &cfg, expected)?;
    let env = capsule::environment().into_iter().collect::<Vec<_>>();
    let args = worker_args(&root, &execution, expected);
    let child = native::measured_child_request(native::ChildRequest {
        program: &root.join("immutable/bin").join(capsule::WORKER),
        args: &args,
        cwd: &execution,
        stem: &execution.join("worker"),
        deadline_ms: cfg.worker_timeout_ms,
        ledger: &execution.join("worker-pids"),
        helper: false,
        input: None,
        environment: &env,
    });
    let worker_drain = native::drain_ledger(&execution.join("worker-pids"));
    // The nested auditor has its own process group; draining the outer group is insufficient.
    let prep_drain = native::drain_ledger(&execution.join("preparation-auditor-pids"));
    let gate = capsule::check_capsule(&root, expected);
    let adoption_gate = if adopted {
        Some(f5_source_custody::verify_adoption(&root, expected))
    } else {
        None
    };
    let evidence = native::inventory(&execution);
    let own_after = std::env::current_exe()
        .map_err(|e| e.to_string())
        .and_then(|p| read(&p, 128 * 1024 * 1024))
        .map(|b| sha256(&b));
    let unchanged = own_after.as_ref().is_ok_and(|v| v == &own);
    let mut terminal = json!({"schema_version":1,"scope":capsule::SCOPE,"registration_sha256":expected,
        "publication_preflight_sha256":preflight_sha256,
        "adoption_preflight_sha256":adoption_preflight_sha256,
        "adoption_gate_passed":adoption_gate.as_ref().map(Result::is_ok),
        "adoption_gate_error":adoption_gate.as_ref().and_then(|result| result.as_ref().err()),
        "controller_sha256":own,"controller_unchanged":unchanged,"controller_read_error":own_after.as_ref().err(),
        "source_gate_passed":gate.is_ok(),"source_gate_error":gate.as_ref().err(),
        "worker_drain_passed":worker_drain.is_ok(),"worker_drain_error":worker_drain.as_ref().err(),
        "preparation_auditor_drain_passed":prep_drain.is_ok(),"preparation_auditor_drain_error":prep_drain.as_ref().err(),
        "execution_files":evidence.as_ref().ok(),"execution_inventory_error":evidence.as_ref().err(),
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"full_goal_complete":false,"online_speedup":null});
    let mut completed = false;
    match child {
        Ok(child) => {
            completed = child.exit_code == Some(0) && !child.timed_out;
            terminal["worker"] = json!({"exit_code":child.exit_code,"timed_out":child.timed_out,"receipt":child.receipt});
        }
        Err(error) => terminal["worker_transport_error"] = json!(error),
    }
    save(&execution.join("terminal.json"), &terminal)?;
    require(completed && unchanged && gate.is_ok()
        && adoption_gate.as_ref().is_none_or(Result::is_ok)
        && worker_drain.is_ok() && prep_drain.is_ok() && evidence.is_ok(),
        "target worker failed/interrupted or source/drain gate failed; terminal retained; never retry")?;
    serde_json::to_string_pretty(&terminal).map_err(|e| e.to_string())
}

fn prefix(
    execution: &Path,
    cfg: &capsule::Config,
    record: &capsule::Registration,
    seal: &str,
) -> Result<Value, String> {
    if execution
        .join("attempts")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        let retained = journal::inspect(execution, seal)?;
        require(
            retained["binding"] == json!(cfg.binding(seal, &record.worker_sha256)),
            "target prefix does not match original target/plan/worker",
        )?;
        Ok(retained)
    } else {
        Ok(
            json!({"status":"NO_DURABLE_TARGET_ATTEMPT_RECORD","completed_attempts":[],"pending_start":null,
            "searches_executed_by_inspection":0,"source_bound_execution_admitted":false,
            "full_goal_complete":false,"online_speedup":null}),
        )
    }
}
pub(super) fn inspect(
    root: &Path,
    execution: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    output_outside(&execution, &root, out)?;
    let record = capsule::check_capsule(&root, expected)?;
    admissible(&record, expected)?;
    let own = frozen_self(&root, &record)?;
    capsule::claim(&root, &execution, &record, expected)?;
    let before = native::inventory(&execution)?;
    let retained = prefix(&execution, &capsule::config(&root)?, &record, expected)?;
    let terminal = if execution
        .join("terminal.json")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        Some(load(&execution.join("terminal.json"))?)
    } else {
        None
    };
    require(
        native::inventory(&execution)? == before,
        "target evidence changed during inspection",
    )?;
    capsule::check_capsule(&root, expected)?;
    require(
        frozen_self(&root, &record)? == own,
        "target inspector changed",
    )?;
    let result = json!({"schema_version":1,"status":"RETAINED_ORIGINAL_TARGET_PREFIX","registration_sha256":expected,
        "checker_sha256":own,"prefix":retained,"terminal":terminal,"execution_files":before,
        "native_children_executed_by_inspector":0,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}
fn ledger(path: &Path, pid: u32) -> Result<(), String> {
    require(
        read(path, 65536)? == format!("start {pid}\ndone {pid}\n").as_bytes(),
        "original target role PID ledger differs",
    )
}
fn status(receipt: &Value, recovered: bool) -> Result<(), String> {
    require(
        receipt["timed_out"] == false && receipt["exit_code"] == if recovered { 0 } else { 1 },
        "worker exit differs from verified recovery/budget-exhausted result",
    )
}
fn interval(
    producer: &Value,
    preparation_receipt: &Value,
    worker_receipt: &Value,
) -> Result<u64, String> {
    let mut total = 0u64;
    for n in [
        &producer["costs"]["observed_wall_ns"],
        &producer["reusable_validation_ns"],
        &producer["solver_setup_ns"],
        &preparation_receipt["child_wall_ns"],
    ] {
        total = total
            .checked_add(
                n.as_u64()
                    .ok_or("missing original target/setup/preparation interval")?,
            )
            .ok_or("target whole-interval overflow")?;
    }
    let outer = worker_receipt["child_wall_ns"]
        .as_u64()
        .ok_or("missing original worker interval")?;
    require(
        total > 0 && total <= outer,
        "target/setup/preparation intervals exceed original worker interval",
    )?;
    Ok(outer - total)
}
fn terminal_binding(value: &Value, seal: &str, own: &str) -> Result<(), String> {
    require(
        value["schema_version"] == 1
            && value["scope"] == capsule::SCOPE
            && value["registration_sha256"] == seal
            && value["controller_sha256"] == own
            && value["controller_unchanged"] == true
            && value["controller_read_error"].is_null()
            && value["source_gate_passed"] == true
            && value["source_gate_error"].is_null()
            && (value["adoption_gate_passed"].is_null() || value["adoption_gate_passed"] == true)
            && value["adoption_gate_error"].is_null()
            && value["worker_drain_passed"] == true
            && value["worker_drain_error"].is_null()
            && value["preparation_auditor_drain_passed"] == true
            && value["preparation_auditor_drain_error"].is_null()
            && value["worker_transport_error"].is_null(),
        "target terminal source/controller/drain gate failed",
    )?;
    provisional(value)
}
fn provisional(value: &Value) -> Result<(), String> {
    require(
        value["source_bound_execution_admitted"] == false
            && value["fresh_paired_qualification"] == false
            && value["headline_eligible"] == false
            && value["promotion_eligible"] == false
            && value["full_goal_complete"] == false
            && value["online_speedup"].is_null(),
        "target producer/controller must not self-admit scientific claims",
    )
}
fn verified_runtime(root: &Path, execution: &Path, seal: &str, own: &str) -> Result<Value, String> {
    let record = capsule::check_capsule(root, seal)?;
    admissible(&record, seal)?;
    require(
        frozen_self(root, &record)? == own,
        "target original auditor identity differs",
    )?;
    capsule::claim(root, execution, &record, seal)?;
    target_build::verify(root, &record)?;
    let cfg = capsule::config(root)?;
    let terminal = load(&execution.join("terminal.json"))?;
    terminal_binding(&terminal, seal, own)?;
    ordinary_control::verify_execution_inventory(execution, &terminal)?;
    let preflight_path = execution.join("publication-preflight.json");
    let preflight_bytes = read(&preflight_path, 16 * 1024 * 1024)?;
    let preflight = target_math::parse(&preflight_bytes)?;
    let publication = Path::new(
        preflight["publication_root"]
            .as_str()
            .ok_or("target publication preflight lacks root")?,
    );
    require(
        sha256(&preflight_bytes) == terminal["publication_preflight_sha256"]
            && preflight == target_custody::scientific_preflight(publication, root, seal)?,
        "original target publication preflight differs from registered bytes",
    )?;
    let adoption = if root
        .join("adoption-start.json")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        let bytes = read(&execution.join("adoption-preflight.json"), 65536)?;
        let checked = f5_source_custody::verify_adoption(root, seal)?;
        require(
            sha256(&bytes) == terminal["adoption_preflight_sha256"]
                && target_math::parse(&bytes)? == checked
                && terminal["adoption_gate_passed"] == true
                && terminal["adoption_gate_error"].is_null(),
            "original F5 source/card adoption preflight differs",
        )?;
        Some(checked)
    } else {
        require(
            terminal["adoption_preflight_sha256"].is_null()
                && terminal["adoption_gate_passed"].is_null()
                && !execution
                    .join("adoption-preflight.json")
                    .try_exists()
                    .map_err(|e| e.to_string())?,
            "unadopted F5 target has an adoption preflight",
        )?;
        None
    };
    require(
        load(&execution.join("target-exposure.json"))? == exposure(&cfg, &record, seal),
        "original target exposure differs",
    )?;
    let receipt = load(&execution.join("worker.receipt.json"))?;
    let argv = std::iter::once(
        root.join("immutable/bin")
            .join(capsule::WORKER)
            .to_string_lossy()
            .into_owned(),
    )
    .chain(worker_args(root, execution, seal))
    .collect::<Vec<_>>();
    let pid = ordinary_control::verify_receipt(
        execution,
        &execution.join("worker"),
        &receipt,
        &argv,
        &record.worker_sha256,
        cfg.worker_timeout_ms,
        &json!(capsule::environment()),
    )?;
    require(
        terminal["worker"]
            == json!({"exit_code":receipt["exit_code"],"timed_out":receipt["timed_out"],"receipt":receipt}),
        "original target controller/worker receipt differs",
    )?;
    ledger(&execution.join("worker-pids"), pid)?;
    require(
        load(&execution.join("worker-started.json"))?
            == json!({"schema_version":1,"scope":capsule::SCOPE,
        "registration_sha256":seal,"worker_sha256":record.worker_sha256,"worker_build_identity":record.worker_build_identity,
        "config_sha256":record.config_sha256,"pid":pid,"source_bound_execution_admitted":false}),
        "target worker start identity differs",
    )?;
    let worker_terminal = load(&execution.join("worker-terminal.json"))?;
    require(
        worker_terminal["schema_version"] == 1
            && worker_terminal["scope"] == capsule::SCOPE
            && worker_terminal["registration_sha256"] == seal
            && worker_terminal["preparation_auditor_drain_confirmed"] == true
            && worker_terminal["preparation_auditor_drain_error"].is_null()
            && worker_terminal["immutable_binding_unchanged"] == true
            && worker_terminal["immutable_binding_error"].is_null(),
        "target worker terminal source/drain differs",
    )?;
    provisional(&worker_terminal)?;

    let (prep_root, prep_execution, checker) = capsule::preparation_paths(&record.preparation)?;
    let prep_receipt = load(&execution.join("preparation-auditor.receipt.json"))?;
    let prep_argv = std::iter::once(checker.to_string_lossy().into_owned())
        .chain(preparation_args(
            &prep_root,
            &prep_execution,
            &record.preparation.registration_sha256,
            &execution.join("preparation-recheck.json"),
        ))
        .collect::<Vec<_>>();
    let prep_pid = ordinary_control::verify_receipt(
        &prep_root,
        &execution.join("preparation-auditor"),
        &prep_receipt,
        &prep_argv,
        &record.preparation.auditor_sha256,
        cfg.worker_timeout_ms,
        &json!(capsule::environment()),
    )?;
    status(&prep_receipt, true)?;
    ledger(&execution.join("preparation-auditor-pids"), prep_pid)?;
    let prep_audit = load(&execution.join("preparation-recheck.json"))?;
    require(
        target_math::parse(&read(
            &execution.join("preparation-auditor.stdout"),
            8 * 1024 * 1024,
        )?)? == prep_audit,
        "original preparation audit stdout/result differs",
    )?;
    let prep_bytes = read(&prep_execution.join("producer.json"), 16 * 1024 * 1024)?;
    let math = capsule::preparation_receipt(
        &record.preparation,
        &prep_audit,
        &target_math::parse_ordinary_producer(&prep_bytes)?,
        &prep_bytes,
    )?;
    require(
        math == original_preparation(&record.preparation)?,
        "original preparation changed",
    )?;
    let producer_bytes = read(&execution.join("producer.json"), 16 * 1024 * 1024)?;
    let producer = target_math::parse(&producer_bytes)?;
    let math_cfg: target_math::Config =
        serde_json::from_value(json!(cfg)).map_err(|e| e.to_string())?;
    let checked = target_math::verify(&math, &math_cfg, &producer)?;
    let attempt_files =
        target_math::verify_records(execution, &math_cfg, &producer, seal, &record.worker_sha256)?;
    let recovered = checked["verified_recovery"]
        .as_bool()
        .ok_or("missing target recovery verdict")?;
    status(&receipt, recovered)?;
    require(
        worker_terminal["producer_verified_recovery"] == recovered,
        "target worker recovery verdict differs from independent replay",
    )?;
    let tail = interval(&producer, &prep_receipt, &receipt)?;
    require(
        capsule::check_capsule(root, seal)?.auditor_sha256 == own
            && frozen_self(root, &record)? == own,
        "target original source/checker changed during audit",
    )?;
    require(
        original_preparation(&record.preparation)? == math
            && read(&prep_execution.join("producer.json"), 16 * 1024 * 1024)? == prep_bytes,
        "original preparation changed during target audit",
    )?;
    ordinary_control::verify_execution_inventory(execution, &terminal)?;
    Ok(
        json!({"schema_version":1,"status":"PASS_NATIVE_SOURCE_BOUND_TARGET_EXECUTION_AUDIT",
        "registration_sha256":seal,"source_bound_execution_admitted":true,"verified_recovery":recovered,
        "target":cfg.target,"target_count":1,"worker_pid":pid,"preparation_auditor_pid":prep_pid,
        "audited_worker_calls":1,"audited_preparation_auditor_calls":1,"producer_sha256":sha256(&producer_bytes),
        "preparation_registration_sha256":record.preparation.registration_sha256,
        "source_adoption":adoption,"mathematics":checked,
        "attempt_files":attempt_files,"recorded_online_interval_ns":producer["costs"]["online_wall_ns"],
        "five_phase_recorded_total_ns":checked["timing"]["five_phase_recorded_total_ns"],
        "exclusive_ic_online_phases_ns":checked["timing"]["exclusive_ic_online_phases_ns"],
        "outer_worker_wall_ns":receipt["child_wall_ns"],"unassigned_worker_setup_and_tail_ns":tail,
        "source_timing_attested":true,"wall_timing_class":"ordinary-host-exploratory; isolation unverified",
        "peak_memory":null,"operation_counts":null,"candidate_id":null,"workload_id":null,"run_id":null,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null}),
    )
}
/// Audit original files only. Never invoke, restore, retry, or drain any process here.
pub(super) fn audit(
    root: &Path,
    execution: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    output_outside(&execution, &root, out)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let own_sha = sha256(&read(&own, 128 * 1024 * 1024)?);
    let result = verified_runtime(&root, &execution, expected, &own_sha);
    let error = result.as_ref().err().cloned();
    let mut receipt = result.unwrap_or_else(|error|json!({"schema_version":1,"status":"TERMINAL_TARGET_EXECUTION_NOT_ADMITTED",
        "registration_sha256":expected,"error":error,"source_bound_execution_admitted":false,"verified_recovery":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null}));
    require(
        sha256(&read(&own, 128 * 1024 * 1024)?) == own_sha,
        "target audit checker changed",
    )?;
    receipt["checker_sha256"] = json!(own_sha);
    receipt["independent_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["native_children_executed_by_auditor"] = json!(0);
    save(out, &receipt)?;
    if let Some(error) = error {
        return Err(format!(
            "target original execution audit failed; receipt retained: {error}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
#[path = "target_control_tests.rs"]
mod tests;
