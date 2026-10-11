//! Retain the sole consumed native F5 control and replay its data, never its worker.
//! This postexecution context does not replace the original frozen admission.
use super::{f5_control, sat_control::native};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path};

pub const REGISTRATION: &str = "5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060";
const ORIGINAL_AUDIT: &str = "994ab5bf1f3153eb6d22425eb15647a152cc1ebd14b68347147c279479a9ed69";
const CHECKER: &str = "5055542e85dc805c6a0c8ee9e09548d53f8bf2a84892e15e02695d580221e8b1";
const CLAIM: &str = "a6e8859aec821e5880107e493794edc456fb3328ee4c6a6b0c75158b91fc61aa";
const TERMINAL: &str = "84daae6c1f1a34d4da8751ce3357b388735c2d6ca1ea1962f466545f3adf97ed";
const ACCEPTANCE: &str = "96dec1470be30bfc8c1afbe2d7e0bd58c42270b726537b36b8f2234517df5c1f";
const SCOPE: &str = "closed-disclosed-n17-native-f5-control";
const ORIGINAL_CAPSULE: &str = "/private/tmp/ic-native-f5-control-20261003-preregister-v1";
const ORIGINAL_EXECUTION: &str = "/private/tmp/ic-native-f5-control-20261003-execution-v1";

pub fn publish(
    capsule: &Path,
    execution: &Path,
    audit: &Path,
    acceptance: &Path,
    out: &Path,
) -> Result<String, String> {
    let registration = f5_control::check_capsule(capsule)?;
    require(
        canonical_sha(&registration)? == REGISTRATION,
        "not the sole accepted F5 registration",
    )?;
    require(
        capsule.canonicalize().map_err(|e| e.to_string())? == Path::new(ORIGINAL_CAPSULE)
            && execution.canonicalize().map_err(|e| e.to_string())?
                == Path::new(ORIGINAL_EXECUTION),
        "publication is not from the original live F5 invocation",
    )?;
    require(
        sha256(&read(audit, 1_048_576)?) == ORIGINAL_AUDIT,
        "original frozen F5 audit differs",
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
        "job.json",
        "host-context.json",
        "consumed.json",
    ] {
        native::create(
            &data.join(name),
            &read(&capsule.join(name), 2 * 1024 * 1024)?,
        )?;
    }
    native::create(
        &data.join("independent-audit.json"),
        &read(audit, 1_048_576)?,
    )?;
    native::create(
        &data.join("source-acceptance.json"),
        &read(acceptance, 1_048_576)?,
    )?;
    let manifest = json!({"schema_version":1,"scope":SCOPE,"registration_sha256":REGISTRATION,
        "original_frozen_audit_sha256":ORIGINAL_AUDIT,"files":native::inventory(&data)?,
        "postexecution_publication_context":true,"fresh_paired_qualification":false,
        "headline_eligible":false,"promotion_eligible":false,"full_goal_complete":false,"online_speedup":null});
    save(&out.join("manifest.json"), &manifest)?;
    save(
        &out.join("seal.json"),
        &json!({"manifest_sha256":canonical_sha(&manifest)?}),
    )?;
    f5_control::check_capsule(capsule)?;
    let result = replay(out)?;
    save(&out.join("publication-replay.json"), &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

/// Portable result replay. Archived code and executables are never loaded or called.
pub fn replay(root: &Path) -> Result<Value, String> {
    let manifest = load(&root.join("manifest.json"))?;
    require(
        load(&root.join("seal.json"))?["manifest_sha256"] == canonical_sha(&manifest)?
            && manifest["schema_version"] == 1
            && manifest["scope"] == SCOPE
            && manifest["registration_sha256"] == REGISTRATION
            && manifest["original_frozen_audit_sha256"] == ORIGINAL_AUDIT
            && manifest["postexecution_publication_context"] == true
            && manifest["fresh_paired_qualification"] == false
            && manifest["headline_eligible"] == false
            && manifest["promotion_eligible"] == false
            && manifest["full_goal_complete"] == false
            && manifest["online_speedup"].is_null(),
        "F5 publication seal/scope differs",
    )?;
    let data = root.join("data");
    native::check_tree(&data, &manifest["files"])?;
    let reg = load(&data.join("registration.json"))?;
    require(
        canonical_sha(&reg)? == REGISTRATION
            && load(&data.join("seal.json"))?["registration_sha256"] == REGISTRATION
            && reg["validation_only"] == false
            && reg["auditor_sha256"] == CHECKER,
        "original executable F5 registration differs",
    )?;
    for name in [
        "config.json",
        "preparation.json",
        "job.json",
        "host-context.json",
    ] {
        require(
            reg[name] == sha256(&read(&data.join(name), 2 * 1024 * 1024)?),
            "frozen F5 input bytes differ",
        )?;
    }
    let audit_bytes = read(&data.join("independent-audit.json"), 1_048_576)?;
    require(
        sha256(&audit_bytes) == ORIGINAL_AUDIT,
        "original frozen F5 audit replaced",
    )?;
    let frozen: Value = serde_json::from_slice(&audit_bytes).map_err(|e| e.to_string())?;
    require(
        frozen["status"] == "PASS_NATIVE_SOURCE_BOUND_F5_CONTROL_AUDIT"
            && frozen["checker_binary_sha256"] == CHECKER
            && frozen["registration_sha256"] == REGISTRATION
            && frozen["source_bound_execution_admitted"] == true
            && frozen["target_complete"] == true
            && frozen["native_children_executed_by_auditor"] == 0
            && frozen["ordinary_queries_executed"] == 0
            && frozen["headline_eligible"] == false
            && frozen["fresh_paired_qualification"] == false
            && frozen["promotion_eligible"] == false
            && frozen["online_speedup"].is_null(),
        "original F5 admission identity/scope differs",
    )?;
    require(
        sha256(&read(&data.join("consumed.json"), 65_536)?) == CLAIM,
        "original F5 claim bytes differ",
    )?;
    require(
        sha256(&read(&data.join("execution/terminal.json"), 65_536)?) == TERMINAL,
        "original F5 terminal bytes differ",
    )?;
    // Original paths are retained provenance; no access to those paths occurs here.
    let mathematics = f5_control::verify_retained(
        &data,
        Path::new(ORIGINAL_CAPSULE),
        Path::new(ORIGINAL_EXECUTION),
    )?;
    require(
        mathematics == frozen["mathematics"],
        "mathematics differs from the original frozen F5 audit",
    )?;
    require(
        sha256(&read(&data.join("source-acceptance.json"), 1_048_576)?) == ACCEPTANCE,
        "original accepted-main gate bytes differ",
    )?;
    let accepted = load(&data.join("source-acceptance.json"))?;
    let checks = accepted["statusCheckRollup"]
        .as_array()
        .ok_or("missing acceptance checks")?;
    let merged_at = accepted["mergedAt"]
        .as_str()
        .ok_or("missing acceptance merge time")?;
    require(
        accepted["state"] == "MERGED"
            && accepted["headRefOid"] == "fd64db667d5be5c0e05aef9cdc55c2c462e415c0"
            && accepted["mergeCommit"]["oid"] == "367920ad035168fefa0965ddd9a957e7e4e664bb"
            && checks.len() == 19
            && checks
                .iter()
                .filter(|c| c["conclusion"] == "SUCCESS")
                .count()
                == 17
            && checks
                .iter()
                .filter(|c| c["conclusion"] == "NEUTRAL")
                .count()
                == 2
            && checks.iter().all(|c| c["status"] == "COMPLETED")
            && checks
                .iter()
                .filter(|c| c["conclusion"] == "SUCCESS")
                .all(|c| {
                    c["completedAt"]
                        .as_str()
                        .is_some_and(|done| done <= merged_at)
                }),
        "retained accepted-main F5 gate record differs",
    )?;
    Ok(
        json!({"schema_version":1,"status":"PASS_POSTEXECUTION_NATIVE_F5_PUBLICATION_REPLAY",
        "registration_sha256":REGISTRATION,"original_frozen_audit_sha256":ORIGINAL_AUDIT,
        "original_source_bound_admission_retained":true,"mathematics":mathematics,
        "postexecution_context":true,"native_children_executed":0,"ordinary_queries_executed":0,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_speedup":null}),
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
    const RETAINED: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-f5-control-registration-v1/result-v1";
    fn retained() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join(RETAINED)
    }
    #[test]
    fn original_result_replays_without_search_or_fresh_claim() {
        let result = replay(&retained()).unwrap();
        assert_eq!(result["native_children_executed"], 0);
        assert_eq!(result["mathematics"]["recovered_scalar"], 24886);
        assert_eq!(
            result["mathematics"]["audited_attempts"]
                .as_array()
                .unwrap()
                .len(),
            3
        );
        assert_eq!(result["fresh_paired_qualification"], false);
    }
    #[test]
    fn resealed_changes_do_not_replace_original_audit_claim_or_receipts() {
        for (i, path) in [
            "independent-audit.json",
            "consumed.json",
            "execution/worker.receipt.json",
            "execution/worker.stdout",
            "execution/worker-pids",
            "source-acceptance.json",
        ]
        .iter()
        .enumerate()
        {
            let dir = std::env::temp_dir().join(format!(
                "native-f5-publication-control-{}-{i}",
                std::process::id()
            ));
            native::copy_tree(&retained(), &dir).unwrap();
            let file = dir.join("data").join(path);
            let mut bytes = fs::read(&file).unwrap();
            if i == 4 {
                bytes.extend_from_slice(b"start 999999\n");
            } else if i == 2 {
                bytes = String::from_utf8(bytes)
                    .unwrap()
                    .replace(
                        "\"RAYON_NUM_THREADS\": \"1\"",
                        "\"RAYON_NUM_THREADS\": \"2\"",
                    )
                    .into_bytes();
            } else if i == 5 {
                bytes = String::from_utf8(bytes)
                    .unwrap()
                    .replace("SUCCESS", "FAILURE")
                    .into_bytes();
            } else {
                bytes.push(b' ');
            }
            fs::write(file, bytes).unwrap();
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
            assert!(replay(&dir).is_err(), "changed {path} passed replay");
            fs::remove_dir_all(dir).unwrap();
        }
    }
    #[test]
    fn resealed_fresh_promotion_and_overwrite_are_rejected() {
        let dir = std::env::temp_dir().join(format!(
            "native-f5-publication-promotion-{}",
            std::process::id()
        ));
        native::copy_tree(&retained(), &dir).unwrap();
        let mut manifest = load(&dir.join("manifest.json")).unwrap();
        manifest["fresh_paired_qualification"] = json!(true);
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
        assert!(replay(&dir).is_err());
        let out = dir.join("create-once.json");
        replay_to(&retained(), &out).unwrap();
        assert!(replay_to(&retained(), &out).is_err());
        fs::remove_dir_all(dir).unwrap();
    }
}
