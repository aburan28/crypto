//! Retain and replay the sole consumed native SAT control, without solver execution.
//! This is postexecution publication context, never a replacement frozen admission.
use super::sat_control::native;
use super::{json as records, prepared_sat, sat_control};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256, ControlConfig};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path};

pub const REGISTRATION: &str = "b098ebbe135abd4706e62d6286f513fe1c4f07fad61f7250487bd2d86ad317e4";
const ORIGINAL_AUDIT: &str = "3570722341b8ca2cd8794d51829cac1a1f3b8c1ff02a07a6af9377a36d2cbab7";
const CHECKER: &str = "c5338662f688c15f52e2b02debf5e7f82ee9c5450a8de7b1ad0961225a27bf41";
const CLAIM: &str = "2d5f2089fb40c228df368640f34b62890653d65a9645251e0e668247535726bb";
const TERMINAL: &str = "fac80e95c7e98e3758dacf5edef3e92290fd89cbae6e8172b445bd00c8b00f7e";

/// Exact original files are copied with create-only writes and retained names.
pub fn publish(
    capsule: &Path,
    execution: &Path,
    audit: &Path,
    context: &Path,
    out: &Path,
) -> Result<String, String> {
    let reg = native::check_capsule(capsule)?;
    require(
        canonical_sha(&reg)? == REGISTRATION,
        "not the accepted sole registration",
    )?;
    let claim = load(&capsule.join("consumed.json"))?;
    require(
        claim["execution"] == json!(execution.canonicalize().map_err(|e| e.to_string())?),
        "execution does not match consumed claim",
    )?;
    require(
        sha256(&read(audit, 1024 * 1024)?) == ORIGINAL_AUDIT,
        "original frozen audit differs",
    )?;
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let data = out.join("data");
    fs::create_dir(&data).map_err(|e| e.to_string())?;
    native::copy_tree(execution, &data.join("execution"))?;
    for name in [
        "registration.json",
        "seal.json",
        "config.json",
        "preparation.json",
        "consumed.json",
    ] {
        native::create(
            &data.join(name),
            &read(&capsule.join(name), 2 * 1024 * 1024)?,
        )?;
    }
    native::create(
        &data.join("independent-audit.json"),
        &read(audit, 1024 * 1024)?,
    )?;
    native::create(&data.join("host-context.txt"), &read(context, 64 * 1024)?)?;
    // Context includes the original operator streams, not a reconstruction of them.
    for suffix in [
        "dispatch.stdout",
        "dispatch.stderr",
        "audit.stdout",
        "audit.stderr",
    ] {
        let source = context
            .parent()
            .ok_or("missing context parent")?
            .join(format!("ic-native-sat-control-20261003-{suffix}"));
        native::create(&data.join(suffix), &read(&source, 1024 * 1024)?)?;
    }
    let manifest = json!({"schema_version":1,"scope":"closed-disclosed-n17-native-sat-control",
        "registration_sha256":REGISTRATION,"original_frozen_audit_sha256":ORIGINAL_AUDIT,
        "files":native::inventory(&data)?,"postexecution_publication_context":true,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,"online_speedup":null});
    save(&out.join("manifest.json"), &manifest)?;
    save(
        &out.join("seal.json"),
        &json!({"manifest_sha256":canonical_sha(&manifest)?}),
    )?;
    native::check_capsule(capsule)?;
    let result = replay(out)?;
    save(&out.join("publication-replay.json"), &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

/// Portable data/math replay. Do not execute archived Mac binaries or redo search.
pub fn replay(root: &Path) -> Result<Value, String> {
    let manifest = load(&root.join("manifest.json"))?;
    require(
        load(&root.join("seal.json"))?["manifest_sha256"] == canonical_sha(&manifest)?
            && manifest["schema_version"] == 1
            && manifest["scope"] == "closed-disclosed-n17-native-sat-control"
            && manifest["registration_sha256"] == REGISTRATION
            && manifest["original_frozen_audit_sha256"] == ORIGINAL_AUDIT
            && manifest["postexecution_publication_context"] == true
            && manifest["fresh_paired_qualification"] == false
            && manifest["headline_eligible"] == false
            && manifest["promotion_eligible"] == false
            && manifest["online_speedup"].is_null(),
        "publication seal/scope differs",
    )?;
    let data = root.join("data");
    native::check_tree(&data, &manifest["files"])?;
    let reg = load(&data.join("registration.json"))?;
    require(
        canonical_sha(&reg)? == REGISTRATION
            && load(&data.join("seal.json"))?["registration_sha256"] == REGISTRATION,
        "original registration differs",
    )?;
    for name in ["config.json", "preparation.json"] {
        require(
            reg[name] == sha256(&read(&data.join(name), 2 * 1024 * 1024)?),
            "frozen input bytes differ",
        )?;
    }
    let cfg: ControlConfig =
        serde_json::from_value(load(&data.join("config.json"))?).map_err(|e| e.to_string())?;
    let prep = load(&data.join("preparation.json"))?;
    let preparation = prepared_sat::verify(&records::parse(&prep.to_string())?)?;
    let audit_bytes = read(&data.join("independent-audit.json"), 1024 * 1024)?;
    require(
        sha256(&audit_bytes) == ORIGINAL_AUDIT,
        "original frozen audit replaced",
    )?;
    let frozen: Value = serde_json::from_slice(&audit_bytes).map_err(|e| e.to_string())?;
    require(
        frozen["status"] == "PASS_NATIVE_SOURCE_BOUND_CONTROL_AUDIT"
            && frozen["auditor_binary_sha256"] == CHECKER
            && reg["auditor_sha256"] == CHECKER
            && frozen["registration_sha256"] == REGISTRATION
            && frozen["source_bound_execution_admitted"] == true
            && frozen["target_complete"] == true
            && frozen["native_children_executed_by_auditor"] == 0
            && frozen["independent_replay_inside_online_interval"] == false,
        "original admission identity differs",
    )?;
    require(
        sha256(&read(&data.join("consumed.json"), 64 * 1024)?) == CLAIM,
        "original claim bytes differ",
    )?;
    let claim = load(&data.join("consumed.json"))?;
    let original_execution = claim["execution"]
        .as_str()
        .ok_or("missing original execution path")?;
    require(
        claim["registration_sha256"] == REGISTRATION && claim["status"] == "consumed-before-launch",
        "one-use claim differs",
    )?;
    let execution = data.join("execution");
    require(
        sha256(&read(&execution.join("terminal.json"), 64 * 1024)?) == TERMINAL,
        "original terminal bytes differ",
    )?;
    let terminal = load(&execution.join("terminal.json"))?;
    require(
        terminal["controller"]["receipt"] == load(&execution.join("controller.receipt.json"))?,
        "controller receipt differs",
    )?;
    require(
        terminal["registration_sha256"] == REGISTRATION
            && terminal["source_gate_passed"] == true
            && terminal["child_drain_passed"] == true
            && terminal["controller"]["exit_code"] == 0
            && terminal["controller"]["timed_out"] == false
            && terminal["controller"]["receipt"]["process_group_drain_confirmed"] == true,
        "terminal custody/drain differs",
    )?;
    let producer_bytes = read(&execution.join("producer.json"), 2 * 1024 * 1024)?;
    require(
        frozen["producer_sha256"] == sha256(&producer_bytes),
        "producer differs from frozen admission",
    )?;
    let producer: Value = serde_json::from_slice(&producer_bytes).map_err(|e| e.to_string())?;
    let (scalar, attempts, online) =
        sat_control::verify_report(&prep, &cfg, &producer, &execution, &reg)?;
    require(
        frozen["recovered_scalar"] == json!(scalar)
            && frozen["audited_attempts"] == json!(attempts)
            && frozen["online_wall_ns"] == online
            && frozen["online_phases_ns"] == producer["online_phases_ns"]
            && frozen["preparation_admission"] == preparation,
        "mathematics differs from original frozen audit",
    )?;
    // Also replay the exact role argv/environment and helper/output receipts.
    let capsule_path = terminal["controller"]["receipt"]["argv"][3]
        .as_str()
        .ok_or("missing original capsule argv")?;
    for trial in 0..attempts.len() {
        let dir = execution.join(format!("query-{trial:02}"));
        let original_dir = Path::new(original_execution).join(format!("query-{trial:02}"));
        let query = load(&dir.join("query.json"))?;
        for (role, binary, args) in [
            (
                "exporter",
                "exporter",
                vec![
                    "17".into(),
                    "6".into(),
                    "standard".into(),
                    cfg.export_nonce.to_string(),
                    cfg.conflict_budget.to_string(),
                    "instance".into(),
                    "1".into(),
                    "0".into(),
                    "--target-x".into(),
                    query["point"][0].to_string(),
                    "--target-y".into(),
                    query["point"][1].to_string(),
                    "--blind-instance-id".into(),
                    format!("native-control-{trial:02}"),
                    "--export-only".into(),
                ],
            ),
            (
                "cms",
                "cms",
                vec![
                    "--verb".into(),
                    "1".into(),
                    "--threads".into(),
                    "1".into(),
                    "--random".into(),
                    "1".into(),
                    "--maxsol".into(),
                    "1".into(),
                    "--maxconfl".into(),
                    cfg.conflict_budget.to_string(),
                    "instance/instance.xor.cnf".into(),
                ],
            ),
        ] {
            let launch = load(&dir.join(format!("{role}.launch.json")))?;
            let receipt = load(&dir.join(format!("{role}.receipt.json")))?;
            let helper_args = json!([
                Path::new(capsule_path).join("immutable/bin/prepared_sat_worker"),
                "child",
                "--capsule",
                capsule_path,
                "--attempt",
                original_dir,
                "--role",
                role
            ]);
            require(
                launch["program"]
                    == json!(Path::new(capsule_path).join(format!("immutable/assets/bin/{binary}")))
                    && launch["argv"] == json!(args)
                    && launch["cwd"] == json!(original_dir)
                    && receipt["argv"] == helper_args
                    && receipt["cwd"] == json!(original_dir)
                    && receipt["environment"] == json!({"LC_ALL":"C"})
                    && receipt["stdout_sha256"]
                        == sha256(&read(&dir.join(format!("{role}.stdout")), 8 * 1024 * 1024)?)
                    && receipt["stderr_sha256"]
                        == sha256(&read(&dir.join(format!("{role}.stderr")), 8 * 1024 * 1024)?),
                "role argv/environment/output differs",
            )?;
        }
    }
    require(
        online
            <= terminal["controller"]["receipt"]["child_wall_ns"]
                .as_u64()
                .ok_or("missing controller clock")?,
        "online exceeds controller",
    )?;
    Ok(
        json!({"schema_version":1,"status":"PASS_POSTEXECUTION_NATIVE_SAT_PUBLICATION_REPLAY",
        "registration_sha256":REGISTRATION,"original_frozen_audit_sha256":ORIGINAL_AUDIT,
        "recovered_scalar":scalar,"scalar_independently_verified":scalar.is_some(),"audited_attempts":attempts,
        "online_diagnostic_ns":online,"preparation_admission":preparation,"postexecution_context":true,
        "original_source_bound_admission_retained":true,"native_children_executed":0,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,"online_speedup":null}),
    )
}

pub fn replay_to(root: &Path, out: &Path) -> Result<String, String> {
    let result = replay(root)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    const RETAINED: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/result-v1";
    #[test]
    fn retained_result_replays_without_any_native_child() {
        let result = replay(&Path::new(env!("CARGO_MANIFEST_DIR")).join(RETAINED)).unwrap();
        assert_eq!(result["recovered_scalar"], 24886);
        assert_eq!(result["native_children_executed"], 0);
        assert_eq!(result["audited_attempts"].as_array().unwrap().len(), 3);
        assert_eq!(result["online_diagnostic_ns"], 56585480167u64);
        assert_eq!(result["fresh_paired_qualification"], false);
    }
    #[test]
    fn altered_role_policy_rejects_even_with_resealed_publication_inventory() {
        let dir =
            std::env::temp_dir().join(format!("native-sat-policy-control-{}", std::process::id()));
        native::copy_tree(&Path::new(env!("CARGO_MANIFEST_DIR")).join(RETAINED), &dir).unwrap();
        let launch_file = dir.join("data/execution/query-02/cms.launch.json");
        let mut launch = load(&launch_file).unwrap();
        launch["argv"][3] = json!("2");
        fs::write(&launch_file, serde_json::to_vec_pretty(&launch).unwrap()).unwrap();
        let mut manifest = load(&dir.join("manifest.json")).unwrap();
        manifest["files"] = native::inventory(&dir.join("data")).unwrap();
        fs::write(
            dir.join("manifest.json"),
            serde_json::to_vec_pretty(&manifest).unwrap(),
        )
        .unwrap();
        fs::write(
            dir.join("seal.json"),
            serde_json::to_vec_pretty(
                &json!({"manifest_sha256":canonical_sha(&manifest).unwrap()}),
            )
            .unwrap(),
        )
        .unwrap();
        assert_eq!(
            replay(&dir).unwrap_err(),
            "role argv/environment/output differs"
        );
        fs::remove_dir_all(&dir).unwrap();
    }
    #[test]
    fn modified_cnf_missing_attempt_and_promoted_manifest_reject() {
        for (i, path) in [
            "data/execution/query-02/instance/instance.xor.cnf",
            "data/execution/query-00/query.json",
            "manifest.json",
        ]
        .iter()
        .enumerate()
        {
            let dir = std::env::temp_dir().join(format!(
                "native-sat-data-control-{}-{i}",
                std::process::id()
            ));
            native::copy_tree(&Path::new(env!("CARGO_MANIFEST_DIR")).join(RETAINED), &dir).unwrap();
            if i == 0 {
                fs::write(dir.join(path), b"corrupted source").unwrap();
            } else if i == 1 {
                fs::remove_file(dir.join(path)).unwrap();
            } else {
                let mut manifest = load(&dir.join(path)).unwrap();
                manifest["fresh_paired_qualification"] = json!(true);
                fs::write(
                    dir.join(path),
                    serde_json::to_vec_pretty(&manifest).unwrap(),
                )
                .unwrap();
                fs::write(
                    dir.join("seal.json"),
                    serde_json::to_vec_pretty(
                        &json!({"manifest_sha256":canonical_sha(&manifest).unwrap()}),
                    )
                    .unwrap(),
                )
                .unwrap();
            }
            assert!(replay(&dir).is_err());
            fs::remove_dir_all(&dir).unwrap();
        }
    }
}
