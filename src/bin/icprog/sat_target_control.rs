//! Frozen one-use controller and original data-only audit for the n17 SAT arm.
//! A publication/card preflight precedes consumption; the consumed capsule
//! cannot launch another worker even when the first launch fails.
use super::{
    ordinary_control, sat_control::native, sat_target_build, sat_target_build::capsule,
    sat_target_build::journal, sat_target_custody, target_math, target_sat_transport,
};
use crypto_lib::cryptanalysis::prepared_sat_control::sha256;
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path, time::Instant};

const SOURCE_QUESTION: &str = "fresh-paired-n17-source-publication-v1";
const CURVE_ID: &str = "EC1N17Ckb1hbbe2b5b6b1e6";

fn admissible(record: &capsule::Registration) -> Result<(), String> {
    require(
        !record.validation_only
            && record.exporter_control_sha256.is_some()
            && record.validation_archive_sha256.is_some(),
        "SAT target registration is validation-only or lacks exporter parity",
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
        "SAT target action requires the original frozen controller",
    )?;
    Ok(hash)
}
fn outside(root: &Path, execution: &Path, out: &Path) -> Result<(), String> {
    let parent = out
        .parent()
        .ok_or("SAT audit output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !parent.starts_with(root) && !parent.starts_with(execution),
        "SAT audit/inspection output must be outside original evidence",
    )
}
fn source_descriptor(
    path: &Path,
    record: &capsule::Registration,
    seal: &str,
    archive: &Value,
) -> Result<(Value, Vec<u8>), String> {
    let bytes = read(path, 65536)?;
    let value = target_math::parse(&bytes)?;
    let fields = [
        "schema_version",
        "question",
        "role",
        "curve_id",
        "registration_sha256",
        "source_manifest_sha256",
        "archive_sha256",
        "target_free",
        "source_bound_execution_admitted",
    ];
    require(
        value
            .as_object()
            .is_some_and(|o| o.len() == fields.len() && fields.iter().all(|&k| o.contains_key(k)))
            && value["schema_version"] == 1
            && value["question"] == SOURCE_QUESTION
            && value["role"] == "cms"
            && value["curve_id"] == CURVE_ID
            && value["registration_sha256"] == seal
            && value["source_manifest_sha256"] == record.source_manifest_sha256
            && value["archive_sha256"] == archive["archive_sha256"]
            && value["target_free"] == true
            && value["source_bound_execution_admitted"] == false
            && read(path, 65536)? == bytes,
        "SAT source descriptor differs from independent scientific publication",
    )?;
    Ok((value, bytes))
}
fn worker_args(root: &Path, execution: &Path, card: &Path, seal: &str) -> Vec<String> {
    vec![
        "run".into(),
        "--capsule".into(),
        root.to_string_lossy().into_owned(),
        "--execution".into(),
        execution.to_string_lossy().into_owned(),
        "--card".into(),
        card.to_string_lossy().into_owned(),
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
fn exposure(card: &capsule::Card, card_sha: &str, descriptor_sha: &str, seal: &str) -> Value {
    json!({"schema_version":1,"scope":capsule::SCOPE,"registration_sha256":seal,
        "card_sha256":card_sha,"source_descriptor_sha256":descriptor_sha,
        "target":card.target,"target_count":1,"input_kind":"supplied-public-point",
        "fixture_scalar_generated":false,
        "exposure_stage":"after four source descriptors; before sole SAT worker launch",
        "fresh_target_certified":false})
}
fn ledger(path: &Path, pid: u32) -> Result<(), String> {
    require(
        read(path, 65536)? == format!("start {pid}\ndone {pid}\n").as_bytes(),
        "SAT worker/preparation PID ledger differs",
    )
}

pub(super) fn execute(
    root: &Path,
    publication: &Path,
    source: &Path,
    card_path: &Path,
    execution: &Path,
    seal: &str,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    let source = source.canonicalize().map_err(|e| e.to_string())?;
    let card_path = card_path.canonicalize().map_err(|e| e.to_string())?;
    let execution_parent = execution
        .parent()
        .ok_or("SAT execution lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !execution_parent.starts_with(&root)
            && !execution_parent.starts_with(&publication)
            && !card_path.starts_with(&root)
            && !card_path.starts_with(&publication),
        "SAT target execution/card must stay outside sealed source and publication",
    )?;
    let record = capsule::check_capsule(&root, seal)?;
    admissible(&record)?;
    let own = frozen_self(&root, &record)?;
    sat_target_build::verify(&root, &record)?;
    sat_target_build::original_preparation(&record.preparation)?;
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "SAT target registration already consumed; never retry",
    )?;
    let publication_preflight =
        sat_target_custody::scientific_preflight(&publication, &root, seal)?;
    let (_, descriptor) = source_descriptor(&source, &record, seal, &publication_preflight)?;
    let descriptor_sha = sha256(&descriptor);
    let (card, card_sha) = capsule::card(&card_path)?;
    require(
        card.source_publications["cms"] == descriptor_sha,
        "SAT public point card does not bind the scientific source descriptor",
    )?;
    let cfg = capsule::config(&root)?;
    fs::create_dir(execution).map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    save(
        &execution.join("publication-preflight.json"),
        &publication_preflight,
    )?;
    native::create(&execution.join("source-descriptor.json"), &descriptor)?;
    let card_bytes = read(&card_path, 65536)?;
    native::create(&execution.join("card.json"), &card_bytes)?;
    require(
        sha256(&card_bytes) == card_sha,
        "SAT public point card changed before one-use claim",
    )?;
    save(
        &root.join("consumed.json"),
        &json!(capsule::Claim {
            schema_version: 1,
            scope: capsule::SCOPE.into(),
            registration_sha256: seal.into(),
            execution: execution.clone(),
            worker_sha256: record.worker_sha256.clone(),
            card: card_path.clone(),
            card_sha256: card_sha.clone(),
            publication_sha256: descriptor_sha.clone(),
        }),
    )?;
    // Any failure after this point consumes the registration permanently.
    save(
        &execution.join("target-exposure.json"),
        &exposure(&card, &card_sha, &descriptor_sha, seal),
    )?;
    let env = capsule::environment().into_iter().collect::<Vec<_>>();
    let args = worker_args(&root, &execution, &card_path, seal);
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
    let prep_drain = native::drain_ledger(&execution.join("preparation-auditor-pids"));
    let source_gate = capsule::check_capsule(&root, seal);
    let card_gate = capsule::card(&card_path).map(|(_, sha)| sha == card_sha);
    let descriptor_gate = read(&source, 65536).map(|bytes| bytes == descriptor);
    let evidence = native::inventory(&execution);
    let own_after = std::env::current_exe()
        .map_err(|e| e.to_string())
        .and_then(|p| read(&p, 128 * 1024 * 1024))
        .map(|bytes| sha256(&bytes));
    let unchanged = own_after.as_ref().is_ok_and(|hash| hash == &own);
    let mut terminal = json!({"schema_version":1,"scope":capsule::SCOPE,
        "registration_sha256":seal,"publication_preflight_sha256":sha256(&read(&execution.join("publication-preflight.json"),65536)?),
        "card_sha256":card_sha,"source_descriptor_sha256":descriptor_sha,
        "controller_sha256":own,"controller_unchanged":unchanged,
        "controller_read_error":own_after.as_ref().err(),
        "source_gate_passed":source_gate.is_ok(),"source_gate_error":source_gate.as_ref().err(),
        "card_gate_passed":card_gate.as_ref().is_ok_and(|v| *v),"card_gate_error":card_gate.as_ref().err(),
        "descriptor_gate_passed":descriptor_gate.as_ref().is_ok_and(|v| *v),"descriptor_gate_error":descriptor_gate.as_ref().err(),
        "worker_drain_passed":worker_drain.is_ok(),"worker_drain_error":worker_drain.as_ref().err(),
        "preparation_auditor_drain_passed":prep_drain.is_ok(),"preparation_auditor_drain_error":prep_drain.as_ref().err(),
        "execution_files":evidence.as_ref().ok(),"execution_inventory_error":evidence.as_ref().err(),
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "headline_eligible":false,"promotion_eligible":false,"online_speedup":null});
    let mut completed = false;
    match child {
        Ok(child) => {
            completed = child.exit_code == Some(0) && !child.timed_out;
            terminal["worker"] = json!({"exit_code":child.exit_code,"timed_out":child.timed_out,"receipt":child.receipt});
        }
        Err(error) => terminal["worker_transport_error"] = json!(error),
    }
    save(&execution.join("terminal.json"), &terminal)?;
    require(
        completed && unchanged && source_gate.is_ok() && card_gate == Ok(true)
            && descriptor_gate == Ok(true) && worker_drain.is_ok() && prep_drain.is_ok()
            && evidence.is_ok(),
        "SAT target worker failed/interrupted or source/drain gate failed; terminal retained; never retry",
    )?;
    serde_json::to_string_pretty(&terminal).map_err(|e| e.to_string())
}

pub(super) fn inspect(
    root: &Path,
    execution: &Path,
    card_path: &Path,
    seal: &str,
    out: &Path,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    outside(&root, &execution, out)?;
    let record = capsule::check_capsule(&root, seal)?;
    admissible(&record)?;
    let own = frozen_self(&root, &record)?;
    let (card, _) = capsule::claim(&root, &execution, card_path, &record, seal)?;
    let before = native::inventory(&execution)?;
    let prefix = if execution
        .join("attempts")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        let prefix = journal::inspect(&execution, seal)?;
        let cfg = capsule::config(&root)?;
        require(
            prefix["binding"]
                == json!(cfg.journal_binding(seal, &record.worker_sha256, card.target)),
            "SAT target prefix differs from frozen card/plan/worker",
        )?;
        prefix
    } else {
        json!({"status":"NO_DURABLE_TARGET_ATTEMPT_RECORD","completed_attempts":[],
            "pending_start":null,"searches_executed_by_inspection":0,
            "source_bound_execution_admitted":false,"online_speedup":null})
    };
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
        "SAT prefix changed during inspection",
    )?;
    capsule::check_capsule(&root, seal)?;
    require(frozen_self(&root, &record)? == own, "SAT inspector changed")?;
    let result = json!({"schema_version":1,"status":"RETAINED_ORIGINAL_SAT_TARGET_PREFIX",
        "registration_sha256":seal,"checker_sha256":own,"prefix":prefix,
        "terminal":terminal,"execution_files":before,
        "native_children_executed_by_inspector":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "online_wall_ns":null,"online_speedup":null});
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

fn verified_runtime(
    root: &Path,
    execution: &Path,
    card_path: &Path,
    seal: &str,
    own: &str,
    transport_out: &Path,
) -> Result<Value, String> {
    let record = capsule::check_capsule(root, seal)?;
    admissible(&record)?;
    require(
        frozen_self(root, &record)? == own,
        "SAT original auditor changed",
    )?;
    let (card, card_sha) = capsule::claim(root, execution, card_path, &record, seal)?;
    sat_target_build::verify(root, &record)?;
    let terminal = load(&execution.join("terminal.json"))?;
    require(
        terminal["schema_version"] == 1
            && terminal["scope"] == capsule::SCOPE
            && terminal["registration_sha256"] == seal
            && terminal["controller_sha256"] == own
            && terminal["controller_unchanged"] == true
            && terminal["controller_read_error"].is_null()
            && terminal["source_gate_passed"] == true
            && terminal["source_gate_error"].is_null()
            && terminal["card_gate_passed"] == true
            && terminal["card_gate_error"].is_null()
            && terminal["descriptor_gate_passed"] == true
            && terminal["descriptor_gate_error"].is_null()
            && terminal["worker_drain_passed"] == true
            && terminal["worker_drain_error"].is_null()
            && terminal["preparation_auditor_drain_passed"] == true
            && terminal["preparation_auditor_drain_error"].is_null()
            && terminal["execution_inventory_error"].is_null()
            && terminal["worker_transport_error"].is_null()
            && terminal["source_bound_execution_admitted"] == false
            && terminal["online_speedup"].is_null(),
        "SAT target terminal source/controller/drain gate failed",
    )?;
    ordinary_control::verify_execution_inventory(execution, &terminal)?;
    let preflight_bytes = read(&execution.join("publication-preflight.json"), 65536)?;
    let preflight = target_math::parse(&preflight_bytes)?;
    let publication = Path::new(
        preflight["publication_root"]
            .as_str()
            .ok_or("SAT publication root missing")?,
    );
    require(
        sha256(&preflight_bytes) == terminal["publication_preflight_sha256"]
            && preflight == sat_target_custody::scientific_preflight(publication, root, seal)?,
        "SAT original scientific publication preflight differs",
    )?;
    let (_, descriptor_bytes) = source_descriptor(
        &execution.join("source-descriptor.json"),
        &record,
        seal,
        &preflight,
    )?;
    let descriptor_sha = sha256(&descriptor_bytes);
    require(
        descriptor_sha == terminal["source_descriptor_sha256"]
            && descriptor_sha == card.source_publications["cms"],
        "SAT original source descriptor/card differs",
    )?;
    let claim = load(&root.join("consumed.json"))?;
    require(
        claim["publication_sha256"] == descriptor_sha,
        "SAT one-use claim differs",
    )?;
    require(
        read(&execution.join("card.json"), 65536)? == read(card_path, 65536)?
            && load(&execution.join("target-exposure.json"))?
                == exposure(&card, &card_sha, &descriptor_sha, seal),
        "SAT original public point exposure differs",
    )?;
    let receipt = load(&execution.join("worker.receipt.json"))?;
    let argv = std::iter::once(
        root.join("immutable/bin")
            .join(capsule::WORKER)
            .to_string_lossy()
            .into_owned(),
    )
    .chain(worker_args(root, execution, card_path, seal))
    .collect::<Vec<_>>();
    let cfg = capsule::config(root)?;
    let pid = ordinary_control::verify_receipt(
        execution,
        &execution.join("worker"),
        &receipt,
        &argv,
        &record.worker_sha256,
        cfg.worker_timeout_ms,
        &json!(capsule::environment()),
    )?;
    ledger(&execution.join("worker-pids"), pid)?;
    require(
        terminal["worker"]
            == json!({"exit_code":receipt["exit_code"],
            "timed_out":receipt["timed_out"],"receipt":receipt})
            && load(&execution.join("worker-started.json"))?
                == json!({"schema_version":1,"scope":capsule::SCOPE,
                    "registration_sha256":seal,"worker_sha256":record.worker_sha256,
                    "worker_build_identity":record.worker_build_identity,
                    "config_sha256":record.config_sha256,"card_sha256":card_sha,
                    "target":card.target,"pid":pid,"source_bound_execution_admitted":false}),
        "SAT original worker identity/launch differs",
    )?;
    let worker_terminal = load(&execution.join("worker-terminal.json"))?;
    require(
        worker_terminal["schema_version"] == 1
            && worker_terminal["scope"] == capsule::SCOPE
            && worker_terminal["registration_sha256"] == seal
            && worker_terminal["card_sha256"] == card_sha
            && worker_terminal["child_drain_confirmed"] == true
            && worker_terminal["child_drain_error"].is_null()
            && worker_terminal["preparation_auditor_drain_confirmed"] == true
            && worker_terminal["preparation_auditor_drain_error"].is_null()
            && worker_terminal["immutable_binding_unchanged"] == true
            && worker_terminal["immutable_binding_error"].is_null()
            && worker_terminal["card_unchanged"] == true
            && worker_terminal["card_error"].is_null()
            && worker_terminal["source_bound_execution_admitted"] == false
            && worker_terminal["online_speedup"].is_null(),
        "SAT worker terminal child/source/card gate failed",
    )?;
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
    ledger(&execution.join("preparation-auditor-pids"), prep_pid)?;
    let prep_audit = load(&execution.join("preparation-recheck.json"))?;
    let prep_bytes = read(&prep_execution.join("producer.json"), 16 * 1024 * 1024)?;
    let prep_producer = target_math::parse_ordinary_producer(&prep_bytes)?;
    let replayed_math = capsule::preparation_receipt(
        &record.preparation,
        &prep_audit,
        &prep_producer,
        &prep_bytes,
    )?;
    require(
        prep_receipt["exit_code"] == 0
            && prep_receipt["timed_out"] == false
            && target_math::parse(&read(
                &execution.join("preparation-auditor.stdout"),
                8 * 1024 * 1024,
            )?)? == prep_audit
            && sat_target_build::original_preparation(&record.preparation)? == replayed_math,
        "SAT original preparation recheck differs",
    )?;
    let transport = target_math::parse(
        target_sat_transport::run(root, execution, card_path, seal, transport_out)?.as_bytes(),
    )?;
    let recovered = transport["mathematics"]["verified_recovery"] == true;
    require(
        transport["status"] == "PASS_SAT_TRANSPORT_DATA_ONLY_RUNTIME_NOT_ADMITTED"
            && receipt["exit_code"] == if recovered { 0 } else { 1 }
            && receipt["timed_out"] == false
            && worker_terminal["producer_verified_recovery"] == recovered,
        "SAT worker exit/recovery differs from independent replay",
    )?;
    ordinary_control::verify_execution_inventory(execution, &terminal)?;
    Ok(
        json!({"schema_version":1,"status":"PASS_NATIVE_SOURCE_BOUND_SAT_TARGET_EXECUTION_AUDIT",
        "registration_sha256":seal,"card_sha256":card_sha,"target":card.target,
        "target_count":1,"source_bound_execution_admitted":true,
        "verified_recovery":recovered,"verified_scalar":transport["mathematics"]["verified_scalar"],
        "mathematics":transport["mathematics"],"native_source":transport["native_source"],
        "worker_pid":pid,"preparation_auditor_pid":prep_pid,
        "worker_sha256":record.worker_sha256,"controller_sha256":own,
        "source_manifest_sha256":record.source_manifest_sha256,
        "publication_archive_sha256":preflight["archive_sha256"],
        "recorded_online_interval_ns":transport["mathematics"]["timing"]["producer_recorded_online_interval_ns"],
        "wall_timing_class":"ordinary-host-exploratory; isolation unverified",
        "candidate_id":null,"workload_id":null,"run_id":null,
        "fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"online_speedup":null}),
    )
}

pub(super) fn audit(
    root: &Path,
    execution: &Path,
    card_path: &Path,
    seal: &str,
    out: &Path,
    transport_out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    outside(&root, &execution, out)?;
    outside(&root, &execution, transport_out)?;
    require(out != transport_out, "SAT audit outputs must be distinct")?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let own_sha = sha256(&read(&own, 128 * 1024 * 1024)?);
    let result = verified_runtime(&root, &execution, card_path, seal, &own_sha, transport_out);
    let error = result.as_ref().err().cloned();
    let mut receipt = result.unwrap_or_else(|reason| {
        json!({"schema_version":1,
        "status":"TERMINAL_SAT_TARGET_EXECUTION_NOT_ADMITTED","registration_sha256":seal,
        "error":reason,"source_bound_execution_admitted":false,"verified_recovery":false,
        "fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"online_speedup":null})
    });
    require(
        sha256(&read(&own, 128 * 1024 * 1024)?) == own_sha,
        "SAT target auditor changed",
    )?;
    receipt["checker_sha256"] = json!(own_sha);
    receipt["independent_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["native_children_executed_by_auditor"] = json!(0);
    save(out, &receipt)?;
    if let Some(error) = error {
        return Err(format!(
            "SAT original execution audit failed; receipt retained: {error}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn old_validation_archive_cannot_dispatch_and_descriptor_cannot_carry_target() {
        let publication = Path::new("research/ic_candidate_tournament_20260915/goal_20260924/native-sat-target-transport-v1/source-freeze-v1/result-v1");
        let seal = "df824d634e829032cb209374d1f395ced43d0604496414cb5d67f464c3e5c865";
        let record = capsule::registration(publication, seal).unwrap();
        assert!(admissible(&record).is_err());
        let archive = json!({"archive_sha256":"d".repeat(64)});
        let mut descriptor = json!({"schema_version":1,"question":SOURCE_QUESTION,
            "role":"cms","curve_id":CURVE_ID,"registration_sha256":seal,
            "source_manifest_sha256":record.source_manifest_sha256,
            "archive_sha256":"d".repeat(64),"target_free":true,
            "source_bound_execution_admitted":false});
        let path = std::env::temp_dir().join(format!(
            "sat-descriptor-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, serde_json::to_vec(&descriptor).unwrap()).unwrap();
        assert!(source_descriptor(&path, &record, seal, &archive).is_ok());
        descriptor["target"] = json!([40991, 73355]);
        fs::write(&path, serde_json::to_vec(&descriptor).unwrap()).unwrap();
        assert!(source_descriptor(&path, &record, seal, &archive).is_err());
        fs::remove_file(path).unwrap();
    }
}
