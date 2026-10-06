//! Data-only replay of the prepared SAT target's original native role files.
//! This is a prerequisite for, never a substitute for, a frozen one-use
//! controller audit and postpublication fresh paired comparison.
use super::{
    json as strict_json, oracle::Curve, ordinary_control, ordinary_preparation,
    sat_control::native, sat_source, target_math, target_sat_math,
};
use crypto_lib::cryptanalysis::prepared_sat_control::sha256;
use serde_json::{json, Value};
use std::{collections::BTreeSet, path::Path, time::Instant};

#[path = "../prepared_sat_target_worker/contract.rs"]
#[allow(dead_code)] // The worker owns some methods; this replay reads their contract.
mod capsule;
#[path = "../prepared_target_worker/journal.rs"]
#[allow(dead_code)]
mod journal;
use super::ordinary_control::capsule as preparation;

fn load(path: &Path) -> Result<Value, String> {
    target_math::parse(&native::read(path, 16 * 1024 * 1024)?)
}
fn require(condition: bool, message: &str) -> Result<(), String> {
    native::require(condition, message)
}
fn marker_prefix(bytes: &[u8], marker: &str) -> Result<String, String> {
    let line = format!("{marker}\n").into_bytes();
    let mut end = 0;
    let mut found = None;
    for part in bytes.split_inclusive(|&byte| byte == b'\n') {
        end += part.len();
        if part == line {
            require(
                found.replace(end).is_none(),
                "prepared readiness marker repeated",
            )?;
        }
    }
    Ok(sha256(
        &bytes[..found.ok_or("missing prepared readiness marker")?],
    ))
}
fn classify_native(receipt: &Value, stdout: &str) -> Result<&'static str, String> {
    let lines = stdout
        .lines()
        .filter(|line| line.starts_with("s "))
        .collect::<Vec<_>>();
    if receipt["timed_out"] == true {
        return Ok("TIMEOUT");
    }
    match (receipt["exit_code"].as_i64(), lines.as_slice()) {
        (Some(10), ["s SATISFIABLE"]) => Ok("SAT_MODEL"),
        (Some(20), ["s UNSATISFIABLE"]) => Ok("SOURCE_UNSAT"),
        (Some(15), ["s INDETERMINATE"]) => Ok("CONFLICT_BUDGET_INCONCLUSIVE"),
        (Some(0), ["s UNKNOWN"]) => Ok("UNKNOWN_INCONCLUSIVE"),
        _ => Err("SAT native exit/status lines are unrecognized or contradictory".into()),
    }
}

struct RoleEvidence {
    ready: Value,
    receipt: Value,
    stdout: Vec<u8>,
}
fn role_evidence(
    root: &Path,
    dir: &Path,
    cfg: &capsule::Config,
    record: &capsule::Registration,
    trial: usize,
    role: &str,
    input: Option<&[u8]>,
    deadline_ms: u64,
) -> Result<RoleEvidence, String> {
    let spec = capsule::prepared_role(
        root,
        cfg,
        trial,
        role,
        &record.prepared_exporter_sha256,
        &record.prepared_cms_sha256,
    )?;
    let ready_path = dir.join(format!("{role}.ready.json"));
    let ready_bytes = native::read(&ready_path, 65536)?;
    let ready = target_math::parse(&ready_bytes)?;
    let receipt = load(&dir.join(format!("{role}.receipt.json")))?;
    let stdout = native::read(&dir.join(format!("{role}.stdout")), 8 * 1024 * 1024)?;
    let stderr = native::read(&dir.join(format!("{role}.stderr")), 8 * 1024 * 1024)?;
    let argv = std::iter::once(spec.program.to_string_lossy().into_owned())
        .chain(spec.args)
        .collect::<Vec<_>>();
    let pid = ready["pid"]
        .as_u64()
        .filter(|&pid| pid > 1)
        .ok_or("prepared role lacks valid PID")?;
    require(
        ready["schema_version"] == 1
            && ready["state"] == "ready-without-input"
            && ready["argv"] == json!(argv)
            && ready["cwd"] == json!(dir)
            && ready["environment"] == json!({"LC_ALL":"C"})
            && ready["executable_sha256"] == spec.pin
            && ready["marker"] == spec.marker
            && ready["stdin_written_bytes"] == 0
            && ready["ready_after_launch_ns"].as_u64().is_some()
            && ready["stdout_at_ready_sha256"] == marker_prefix(&stdout, spec.marker)?
            && receipt["schema_version"] == 1
            && receipt["pid"] == pid
            && receipt["launch_to_ready_ns"] == ready["ready_after_launch_ns"]
            && receipt["executable_sha256_before"] == spec.pin
            && receipt["executable_sha256_after"] == spec.pin
            && receipt["readiness_marker_count"] == 1
            && receipt["process_group_drain_confirmed"] == true
            && receipt["stdout_bytes"] == stdout.len()
            && receipt["stdout_sha256"] == sha256(&stdout)
            && receipt["stderr_bytes"] == stderr.len()
            && receipt["stderr_sha256"] == sha256(&stderr),
        "prepared role identity/readiness/original output differs",
    )?;
    match input {
        Some(bytes) => require(
            receipt["state"] == "used"
                && receipt["argv"] == json!(argv)
                && receipt["cwd"] == json!(dir)
                && receipt["environment"] == json!({"LC_ALL":"C"})
                && receipt["marker"] == spec.marker
                && receipt["stdin_bytes"] == bytes.len()
                && receipt["stdin_sha256"] == sha256(bytes)
                && receipt["stdin_written_bytes"] == bytes.len()
                && receipt["stdin_write_error"].is_null()
                && receipt["deadline_ms"] == deadline_ms
                && receipt["timed_out"].is_boolean()
                && receipt["output_limit"] == false,
            "prepared role original stdin/deadline differs",
        )?,
        None => require(
            receipt["state"] == "cancelled-unused"
                && receipt["stdin_written_bytes"] == 0
                && receipt["ready_to_cancel_ns"].as_u64().is_some(),
            "unused prepared role received input or lacks cancellation",
        )?,
    }
    require(
        native::read(&ready_path, 65536)? == ready_bytes,
        "prepared role readiness changed during replay",
    )?;
    Ok(RoleEvidence {
        ready,
        receipt,
        stdout,
    })
}

fn geometric_base(math: &Value) -> Result<(Curve, Vec<super::oracle::Point>), String> {
    let curve = Curve::new(&strict_json::parse(&math["fixture"].to_string())?)?;
    let base = math["geometric_base"]
        .as_array()
        .ok_or("SAT preparation lacks geometric base")?
        .iter()
        .map(|point| curve.decode(&strict_json::parse(&point.to_string())?))
        .collect::<Result<Vec<_>, _>>()?;
    require(base.len() == 63, "SAT geometric base count differs")?;
    Ok((curve, base))
}
fn model_xs(model: &[bool]) -> Result<[u64; 3], String> {
    require(model.len() >= 51, "SAT model lacks source variables")?;
    Ok(std::array::from_fn(|i| {
        (0..6).fold(0, |value, bit| value | ((model[i * 6 + bit] as u64) << bit))
    }))
}
fn has_lift_for_xs(
    curve: &Curve,
    base: &[super::oracle::Point],
    point: super::oracle::Point,
    xs: [u64; 3],
) -> bool {
    base.iter()
        .filter(|p| p.is_some_and(|p| p.0 == xs[0]))
        .any(|a| {
            base.iter()
                .filter(|p| p.is_some_and(|p| p.0 == xs[1]))
                .any(|b| {
                    base.iter()
                        .filter(|p| p.is_some_and(|p| p.0 == xs[2]))
                        .any(|c| curve.add(curve.add(*a, *b), *c) == point)
                })
        })
}

fn verify_sources(
    root: &Path,
    execution: &Path,
    cfg: &capsule::Config,
    record: &capsule::Registration,
    math: &Value,
    producer: &Value,
) -> Result<Value, String> {
    let (curve, base) = geometric_base(math)?;
    let attempts = producer["attempts"]
        .as_array()
        .ok_or("missing SAT target attempts")?;
    let pool = load(&execution.join("pool-ready.json"))?;
    let mut expected_ready = Vec::new();
    let mut cancelled_exporters = Vec::new();
    let mut cancelled_cms = Vec::new();
    let mut pids = Vec::new();
    let mut checked = Vec::new();
    for trial in 0..cfg.plan.max_queries {
        let dir = execution.join(format!("query-{trial:03}"));
        let row = attempts.get(trial);
        let used = row.is_some_and(|row| row["backend_called"] == true);
        if used {
            let row = row.ok_or("missing used SAT target attempt")?;
            let point: [u64; 2] =
                serde_json::from_value(row["public_query"].clone()).map_err(|e| e.to_string())?;
            require(
                load(&dir.join("query.json"))? == json!({"trial":trial,"point":point}),
                "SAT original query file differs from mathematical attempt",
            )?;
            let mut request = serde_json::to_vec(&json!({
                "target_x":point[0].to_string(),"target_y":point[1].to_string(),
                "blind_instance_id":format!("native-sat-target-{trial:03}")
            }))
            .map_err(|e| e.to_string())?;
            request.push(b'\n');
            let exporter = role_evidence(
                root,
                &dir,
                cfg,
                record,
                trial,
                "exporter",
                Some(&request),
                cfg.plan.exporter_timeout_ms,
            )?;
            require(
                exporter.receipt["exit_code"] == 0 && exporter.receipt["timed_out"] == false,
                "SAT exporter failure cannot yield a complete source attempt",
            )?;
            let (manifest, anf, cnf, files) =
                native::validate_exports(&dir.join("instance"), point)?;
            let cms = role_evidence(
                root,
                &dir,
                cfg,
                record,
                trial,
                "cms",
                Some(cnf.as_bytes()),
                cfg.plan.solver_timeout_ms,
            )?;
            let stdout = std::str::from_utf8(&cms.stdout).map_err(|e| e.to_string())?;
            let status = classify_native(&cms.receipt, stdout)?;
            let source = &row["source_evidence"];
            let expected_source = json!({
                "manifest":manifest,"source_receipt":{
                    "files":files.clone(),"files_unchanged_before_after":true,
                    "exporter":exporter.receipt.clone(),"exporter_ready":exporter.ready.clone(),
                    "cms_ready":cms.ready.clone()},
                "anf_bytes":anf.len(),"anf_sha256":sha256(anf.as_bytes()),
                "cnf_bytes":cnf.len(),"cnf_sha256":sha256(cnf.as_bytes()),
                "stdout_bytes":cms.stdout.len(),"stdout_sha256":sha256(&cms.stdout),
                "native":{"exit_code":cms.receipt["exit_code"],
                    "timed_out":cms.receipt["timed_out"],"receipt":cms.receipt.clone()},
                "payload_policy":"outer executor retains original files"
            });
            require(
                source == &expected_source && row["native_status"] == status,
                "SAT target original source/model/role receipt differs",
            )?;
            if status == "SAT_MODEL" {
                let model = sat_source::verify_native_model(&anf, &cnf, stdout)?;
                require(
                    row["source_model_valid"] == true
                        && row["expanded_model_sha256"]
                            == sha256(&model.iter().map(|&bit| u8::from(bit)).collect::<Vec<_>>()),
                    "SAT source-valid model digest differs",
                )?;
                let xs = model_xs(&model)?;
                let query = curve.decode(&strict_json::parse(&json!(point).to_string())?)?;
                match row["outcome"].as_str() {
                    Some("verified_point_witness") => {
                        let indices: [usize; 3] =
                            serde_json::from_value(row["witness_indices"].clone())
                                .map_err(|e| e.to_string())?;
                        require(
                            indices.iter().enumerate().all(|(i, &index)| {
                                base.get(index)
                                    .is_some_and(|p| p.is_some_and(|p| p.0 == xs[i]))
                            }) && has_lift_for_xs(&curve, &base, query, xs),
                            "SAT model does not lift to the reported full-point witness",
                        )?;
                    }
                    Some("nonlifting_model") => require(
                        !has_lift_for_xs(&curve, &base, query, xs),
                        "SAT model labeled nonlifting despite a full-point lift",
                    )?,
                    _ => return Err("SAT model has unadmitted target disposition".into()),
                }
            }
            pids.push(
                exporter.ready["pid"]
                    .as_u64()
                    .ok_or("missing exporter PID")?,
            );
            pids.push(cms.ready["pid"].as_u64().ok_or("missing CMS PID")?);
            expected_ready.push(
                json!({"trial":trial,"role":"exporter","receipt":exporter.ready,
                "receipt_sha256":sha256(&native::read(&dir.join("exporter.ready.json"),65536)?)}),
            );
            expected_ready.push(json!({"trial":trial,"role":"cms","receipt":cms.ready,
                "receipt_sha256":sha256(&native::read(&dir.join("cms.ready.json"),65536)?)}));
            checked.push(json!({"trial":trial,"status":status,"outcome":row["outcome"],
                "source_files":files,"exporter_pid":pids[pids.len()-2],"cms_pid":pids[pids.len()-1]}));
        } else {
            require(
                !dir.join("query.json")
                    .try_exists()
                    .map_err(|e| e.to_string())?
                    && ![
                        "manifest.json",
                        "instance.anf",
                        "instance.xor.cnf",
                        "instance.magma",
                    ]
                    .iter()
                    .any(|name| dir.join("instance").join(name).exists()),
                "unused SAT target slot received a query or exported source",
            )?;
            for role in ["exporter", "cms"] {
                let evidence = role_evidence(root, &dir, cfg, record, trial, role, None, 0)?;
                require(
                    !evidence
                        .stdout
                        .split(|&byte| byte == b'\n')
                        .any(|line| line.starts_with(b"s ")),
                    "unused SAT role emitted a solve status",
                )?;
                pids.push(evidence.ready["pid"].as_u64().ok_or("missing unused PID")?);
                expected_ready.push(json!({"trial":trial,"role":role,"receipt":evidence.ready,
                    "receipt_sha256":sha256(&native::read(&dir.join(format!("{role}.ready.json")),65536)?)}));
                let cancelled = json!({"role":role,"trial":trial,"receipt":evidence.receipt});
                if role == "exporter" {
                    cancelled_exporters.push(cancelled);
                } else {
                    cancelled_cms.push(cancelled);
                }
            }
        }
    }
    require(
        pool == json!({"schema_version":1,"scope":capsule::SCOPE,
        "max_queries":cfg.plan.max_queries,"prepared_roles":expected_ready,
        "all_roles_ready_without_input":true,"source_bound_execution_admitted":false}),
        "SAT pool readiness boundary differs",
    )?;
    let mut cancelled = cancelled_exporters;
    cancelled.extend(cancelled_cms);
    require(
        load(&execution.join("prepared-pool-cancellation.json"))?
            == json!({"schema_version":1,"cancelled":cancelled,"errors":[],
            "all_unused_children_drained":true}),
        "SAT unused prepared roles were not all cancelled",
    )?;
    let ledger_bytes = native::read(&execution.join("child-pids"), 65536)?;
    let ledger = std::str::from_utf8(&ledger_bytes)
        .map_err(|e| e.to_string())?
        .lines()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    require(
        ledger.len() == 2 * pids.len()
            && ledger[..pids.len()]
                == pids
                    .iter()
                    .map(|pid| format!("start {pid}"))
                    .collect::<Vec<_>>(),
        "SAT native PID ledger launch sequence differs",
    )?;
    let expected = pids.iter().copied().collect::<BTreeSet<_>>();
    let done = ledger[pids.len()..]
        .iter()
        .map(|line| {
            line.strip_prefix("done ")
                .ok_or("SAT PID ledger has non-done tail")?
                .parse::<u64>()
                .map_err(|e| e.to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    require(
        expected.len() == pids.len()
            && done.len() == pids.len()
            && done.iter().copied().collect::<BTreeSet<_>>() == expected,
        "SAT native PID ledger has missing, repeated or extra drain",
    )?;
    Ok(
        json!({"schema_version":1,"status":"PASS_PREPARED_SAT_SOURCE_BYTES_DATA_ONLY",
        "audited_attempts":checked,"prepared_roles":pids.len(),
        "cancelled_unused_roles":cancelled.len(),"native_solver_calls_by_auditor":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "online_speedup":null}),
    )
}

pub(super) fn run(
    capsule_path: &Path,
    execution_path: &Path,
    card_path: &Path,
    seal: &str,
    out: &Path,
) -> Result<String, String> {
    let start = Instant::now();
    let root = capsule_path.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution_path.canonicalize().map_err(|e| e.to_string())?;
    let card_path = card_path.canonicalize().map_err(|e| e.to_string())?;
    let parent = out
        .parent()
        .ok_or("SAT transport audit output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !parent.starts_with(&root) && !parent.starts_with(&execution),
        "SAT transport audit output would alter original evidence",
    )?;
    let before = native::inventory(&execution)?;
    let result = (|| -> Result<Value, String> {
        let record = capsule::check_capsule(&root, seal)?;
        let (card, card_sha) = capsule::claim(&root, &execution, &card_path, &record, seal)?;
        let cfg = capsule::config(&root)?;
        let (prep_root, prep_execution, _) = capsule::preparation_paths(&record.preparation)?;
        let (prep_record, _) = preparation::check_capsule(&prep_root)?;
        ordinary_control::verify_build_from(&prep_root, &prep_record, &prep_root)?;
        let prep_bytes = native::read(&prep_execution.join("producer.json"), 16 * 1024 * 1024)?;
        let prep_producer = target_math::parse_ordinary_producer(&prep_bytes)?;
        let math = capsule::preparation_receipt(
            &record.preparation,
            &load(&execution.join("preparation-recheck.json"))?,
            &prep_producer,
            &prep_bytes,
        )?;
        let checked_preparation = ordinary_preparation::audit(&math)?;
        require(
            checked_preparation["declared_family"] == "cryptominisat"
                && checked_preparation["rank"] == 29
                && checked_preparation["panel_complete"] == true,
            "SAT original preparation failed independent mathematical replay",
        )?;
        let producer = load(&execution.join("producer.json"))?;
        let math_cfg: target_sat_math::Config = serde_json::from_value(json!({
            "schema_version":cfg.schema_version,"question":cfg.question,
            "plan":&cfg.plan,"worker_timeout_ms":cfg.worker_timeout_ms,
            "target":card.target}))
        .map_err(|e| e.to_string())?;
        let checked_math = target_sat_math::verify(&math, &math_cfg, &producer)?;
        let checked_sources = verify_sources(&root, &execution, &cfg, &record, &math, &producer)?;
        require(
            capsule::card(&card_path)?.1 == card_sha,
            "SAT public card changed during transport replay",
        )?;
        Ok(
            json!({"schema_version":1,"status":"PASS_SAT_TRANSPORT_DATA_ONLY_RUNTIME_NOT_ADMITTED",
            "registration_sha256":seal,"card_sha256":card_sha,
            "mathematics":checked_math,"native_source":checked_sources,
            "producer_sha256":sha256(&native::read(&execution.join("producer.json"),16*1024*1024)?),
            "native_solver_calls_by_auditor":0,"source_bound_execution_admitted":false,
            "fresh_paired_qualification":false,"headline_eligible":false,
            "promotion_eligible":false,"full_goal_complete":false,"online_speedup":null}),
        )
    })();
    let result = result.and_then(|value| {
        capsule::check_capsule(&root, seal)?;
        Ok(value)
    });
    require(
        native::inventory(&execution)? == before,
        "SAT execution changed during data-only transport audit",
    )?;
    let error = result.as_ref().err().cloned();
    let mut receipt = result.unwrap_or_else(|reason| {
        json!({"schema_version":1,
        "status":"SAT_TRANSPORT_NOT_VERIFIED","error":reason,
        "native_solver_calls_by_auditor":0,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"full_goal_complete":false,"online_speedup":null})
    });
    receipt["independent_audit_wall_ns"] = json!(start.elapsed().as_nanos().to_string());
    native::save(out, &receipt)?;
    if let Some(error) = error {
        return Err(format!(
            "SAT transport audit failed; receipt retained: {error}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn marker_prefix_requires_one_exact_complete_line() {
        let bytes = b"c banner\nc PREPARED_STDIN_READY_v1\ns SATISFIABLE\n";
        assert_eq!(
            marker_prefix(bytes, "c PREPARED_STDIN_READY_v1").unwrap(),
            sha256(b"c banner\nc PREPARED_STDIN_READY_v1\n")
        );
        assert!(marker_prefix(
            b"c PREPARED_STDIN_READY_v1\nc PREPARED_STDIN_READY_v1\n",
            "c PREPARED_STDIN_READY_v1"
        )
        .is_err());
        assert!(marker_prefix(
            b"c PREPARED_STDIN_READY_v1-extra\n",
            "c PREPARED_STDIN_READY_v1"
        )
        .is_err());
    }
    #[test]
    fn native_status_requires_exact_exit_and_one_status_line() {
        assert_eq!(
            classify_native(
                &json!({"timed_out":false,"exit_code":10}),
                "c PREPARED_STDIN_READY_v1\ns SATISFIABLE\nv 1 0\n"
            )
            .unwrap(),
            "SAT_MODEL"
        );
        assert_eq!(
            classify_native(
                &json!({"timed_out":true,"exit_code":null}),
                "c PREPARED_STDIN_READY_v1\n"
            )
            .unwrap(),
            "TIMEOUT"
        );
        assert!(classify_native(
            &json!({"timed_out":false,"exit_code":10}),
            "s UNSATISFIABLE\n"
        )
        .is_err());
        assert!(classify_native(
            &json!({"timed_out":false,"exit_code":20}),
            "s UNSATISFIABLE\ns SATISFIABLE\n"
        )
        .is_err());
    }
}
