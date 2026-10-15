//! Post-hoc, collect-only extraction from a fully verified immutable campaign.
//! This program can execute only sha256sum; it has no solver-launch operation.
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
fn err(s: impl Into<String>) -> Box<dyn std::error::Error> {
    s.into().into()
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(err(message))
    }
}
fn sha(bytes: &[u8]) -> Result<String> {
    let mut child = Command::new("/usr/bin/sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    child
        .stdin
        .take()
        .ok_or_else(|| err("missing hash stdin"))?
        .write_all(bytes)?;
    let out = child.wait_with_output()?;
    require(out.status.success(), "sha256sum failed")?;
    Ok(String::from_utf8(out.stdout)?
        .split_whitespace()
        .next()
        .ok_or_else(|| err("empty hash"))?
        .into())
}
fn digest_ok(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v[key]
        .as_str()
        .ok_or_else(|| err(format!("missing string {key}")))
}
fn parse(bytes: &[u8]) -> Result<Value> {
    Ok(serde_json::from_slice(bytes)?)
}
#[derive(Clone)]
struct Bound {
    path: PathBuf,
    digest: String,
    bytes: Vec<u8>,
}
fn pinned(v: &Value) -> Result<Bound> {
    let requested = PathBuf::from(text(v, "path")?);
    require(requested.is_absolute(), "pin path must be absolute")?;
    let path = fs::canonicalize(&requested)?;
    require(fs::metadata(&path)?.is_file(), "pin is not a regular file")?;
    let digest = text(v, "sha256")?.to_string();
    require(digest_ok(&digest), "malformed pinned SHA256")?;
    let bytes = fs::read(&path)?;
    require(
        sha(&bytes)? == digest,
        format!("pinned file changed: {}", path.display()),
    )?;
    Ok(Bound {
        path,
        digest,
        bytes,
    })
}
fn full_gate(r: &Value) -> Result<()> {
    require(
        r["kind"] == "native_n83_matched_answer_and_file_verification"
            && r["status"] == "PASS"
            && r["complete_panel_gate"] == true
            && r["aa_semantic_gate"] == true,
        "independent full-panel gate is not PASS",
    )?;
    for (key, n) in [
        ("planned_processes", 32),
        ("observed_processes", 32),
        ("completed_verified_processes", 32),
        ("planned_target_rows", 218),
        ("checked_target_rows", 218),
        ("independent_scalar_replays", 218),
        ("IC_scalar_replays", 104),
        ("rho_scalar_replays", 114),
        ("IC_four_point_witness_replays", 104),
        ("checked_aa_semantic_pairs", 5),
        ("invalid_rows", 0),
        ("missing_target_rows", 0),
        ("nonanswer_rows", 0),
        ("launch_failed_processes", 0),
        ("clean_environment_records_checked", 32),
        ("resource_envelopes_checked", 32),
        ("observation_execution_bindings_checked", 32),
        ("all_files_snapshotted", 336),
    ] {
        require(
            r[key].as_u64() == Some(n),
            format!("independent gate count mismatch: {key}"),
        )?;
    }
    require(
        r["gate_errors"].as_array().is_some_and(Vec::is_empty),
        "independent gate reports errors",
    )?;
    require(
        r["rows"].as_array().is_some_and(|a| a.len() == 218),
        "independent row audit incomplete",
    )?;
    Ok(())
}
struct Seal {
    bytes: BTreeMap<PathBuf, Vec<u8>>,
    hashes: BTreeMap<PathBuf, String>,
}
impl Seal {
    fn load(receipt: &Value) -> Result<Self> {
        let entries = receipt["files"]
            .as_array()
            .ok_or_else(|| err("receipt file snapshots missing"))?;
        require(
            entries.len() == 336,
            "receipt must seal all336 snapshot files",
        )?;
        let mut bytes = BTreeMap::new();
        let mut hashes = BTreeMap::new();
        for entry in entries {
            let bound =
                pinned(&json!({"path":text(entry,"path")?,"sha256":text(entry,"sha256")?}))?;
            require(
                !bytes.contains_key(&bound.path),
                format!("duplicate snapshot path: {}", bound.path.display()),
            )?;
            hashes.insert(bound.path.clone(), bound.digest);
            bytes.insert(bound.path, bound.bytes);
        }
        Ok(Self { bytes, hashes })
    }
    fn canonical(&self, p: &Path) -> Result<PathBuf> {
        let p = fs::canonicalize(p)?;
        require(
            self.bytes.contains_key(&p),
            format!("unsealed input: {}", p.display()),
        )?;
        Ok(p)
    }
    fn read(&self, p: &Path) -> Result<&[u8]> {
        let p = self.canonical(p)?;
        Ok(&self.bytes[&p])
    }
    fn json(&self, p: &Path) -> Result<Value> {
        parse(self.read(p)?)
    }
    fn binding(&self, p: &Path, digest: &str) -> Result<()> {
        let p = self.canonical(p)?;
        require(
            self.hashes[&p] == digest,
            format!("sealed binding mismatch: {}", p.display()),
        )
    }
    fn exact_tree(&self, root: &Path) -> Result<()> {
        fn visit(seal: &Seal, path: &Path) -> Result<()> {
            let meta = fs::symlink_metadata(path)?;
            require(
                !meta.file_type().is_symlink(),
                "symlink inside sealed campaign",
            )?;
            if meta.is_dir() {
                for entry in fs::read_dir(path)? {
                    visit(seal, &entry?.path())?;
                }
            } else {
                require(meta.is_file(), "nonregular campaign entry")?;
                seal.canonical(path)?;
            }
            Ok(())
        }
        visit(self, root)
    }
}
#[derive(Clone)]
struct Run {
    id: String,
    kind: String,
    arm: String,
    cohort: usize,
    pair: usize,
    leg: String,
    ids: Vec<String>,
    seeds: Vec<u64>,
}
fn ids() -> Vec<String> {
    (1..=32).map(|i| format!("target_{i:02}")).collect()
}
fn matrix() -> Vec<Run> {
    let mut out = Vec::new();
    for pair in 1..=5 {
        for leg in ["A", "B"] {
            out.push(Run {
                id: format!("AA{pair:02}-{leg}"),
                kind: "aa".into(),
                arm: "rho".into(),
                cohort: 0,
                pair,
                leg: leg.into(),
                ids: vec!["target_01".into()],
                seeds: vec![202610146410 + pair as u64],
            });
        }
    }
    for cohort in 1..=3 {
        let mut indices: Vec<_> = (0..32).map(|i| (i + (cohort - 1) * 11) % 32).collect();
        if cohort == 2 {
            indices.reverse();
        }
        let targets: Vec<_> = indices
            .iter()
            .map(|i| format!("target_{:02}", i + 1))
            .collect();
        let seeds: Vec<_> = indices
            .iter()
            .map(|i| 202610146501 + ((cohort - 1) * 32 + i) as u64)
            .collect();
        for arm in if cohort % 2 == 1 {
            ["rho", "ic"]
        } else {
            ["ic", "rho"]
        } {
            out.push(Run {
                id: format!("C{cohort:02}-{arm}"),
                kind: "cohort".into(),
                arm: arm.into(),
                cohort,
                pair: 0,
                leg: String::new(),
                ids: targets.clone(),
                seeds: seeds.clone(),
            });
        }
        for pair in match cohort {
            1 => 1..=3,
            2 => 4..=6,
            _ => 7..=8,
        } {
            for arm in if pair % 2 == 1 {
                ["ic", "rho"]
            } else {
                ["rho", "ic"]
            } {
                out.push(Run {
                    id: format!("K{pair:02}-{arm}"),
                    kind: "cold".into(),
                    arm: arm.into(),
                    cohort: 0,
                    pair,
                    leg: String::new(),
                    ids: vec![format!("target_{pair:02}")],
                    seeds: vec![202610146600 + pair as u64],
                });
            }
        }
    }
    out
}
fn matrix_entry(r: &Run, declared: &Value) -> Result<()> {
    require(
        declared["run_id"] == r.id
            && declared["kind"] == r.kind
            && declared["arm"] == r.arm
            && declared["cohort_id"].as_u64() == Some(r.cohort as u64)
            && declared["pair_id"].as_u64() == Some(r.pair as u64)
            && declared["leg"] == r.leg
            && declared["target_ids"] == json!(r.ids),
        "sealed run matrix identity mismatch",
    )?;
    if r.arm == "rho" {
        require(
            declared["walk_seeds"] == json!(r.seeds) && declared["rank_seed"].is_null(),
            "sealed rho seed matrix mismatch",
        )?;
    } else {
        require(
            declared["walk_seeds"].is_null() && declared["rank_seed"] == 29,
            "sealed IC seed matrix mismatch",
        )?;
    }
    Ok(())
}
fn output_path(path: &Path, campaign: &Path, protected: &BTreeSet<PathBuf>) -> Result<PathBuf> {
    require(
        !path.exists() && !fs::symlink_metadata(path).is_ok(),
        "refusing existing output",
    )?;
    let parent = fs::canonicalize(path.parent().ok_or_else(|| err("output parent missing"))?)?;
    let final_path = parent.join(
        path.file_name()
            .ok_or_else(|| err("output filename missing"))?,
    );
    require(
        !final_path.starts_with(campaign) && !protected.contains(&final_path),
        "output overlaps sealed campaign/input",
    )?;
    Ok(final_path)
}
fn source_provenance(config: &Path, config_bytes: &[u8]) -> Result<Value> {
    let exe = std::env::current_exe()?;
    Ok(
        json!({"kind":"separate_posthoc_collect_only_extractor","compiled_source_sha256":{"main":sha(include_bytes!("main.rs"))?,"manifest":sha(include_bytes!("../Cargo.toml"))?,"lock":sha(include_bytes!("../Cargo.lock"))?},"executable":exe,"executable_sha256":sha(&fs::read(&exe)?)?,"config":config,"config_sha256":sha(config_bytes)?,"numeric_copy_policy":"serde_json arbitrary_precision retains sealed observation number tokens; no observe reconstruction, float tolerance or clock recomputation","original_frozen_collect_succeeded":false,"solver_launch_operation_present":false}),
    )
}
fn extract(config_file: &Path, out: &Path, provenance_out: &Path) -> Result<()> {
    let config_path = fs::canonicalize(config_file)?;
    let config_bytes = fs::read(&config_path)?;
    let c = parse(&config_bytes)?;
    require(
        c["schema"] == "n83.matched-collect-only-config/v1",
        "collector config schema mismatch",
    )?;
    let campaign = fs::canonicalize(text(&c, "campaign_dir")?)?;
    let audit = pinned(&c["independent_receipt"])?;
    let receipt = parse(&audit.bytes)?;
    full_gate(&receipt)?;
    let protocol = pinned(&c["protocol"])?;
    let manifest = pinned(&c["manifest"])?;
    require(
        receipt["protocol_sha256"] == protocol.digest
            && receipt["campaign_manifest_sha256"] == manifest.digest,
        "independent receipt protocol/manifest mismatch",
    )?;
    require(
        manifest.path == campaign.join("campaign_manifest.json"),
        "manifest is outside declared campaign",
    )?;
    let p = parse(&protocol.bytes)?;
    let m = parse(&manifest.bytes)?;
    require(
        m["protocol_sha256"] == protocol.digest && m["protocol"] == p,
        "manifest protocol binding mismatch",
    )?;
    require(
        p["schema"] == "n83.matched-campaign/v1" && m["schema"] == "n83.matched-manifest/v1",
        "frozen schema mismatch",
    )?;
    let checker_source = pinned(&c["checker_source"])?;
    let checker_exe = pinned(&c["checker_executable"])?;
    let helper = pinned(&c["numeric_helper"])?;
    let v2_source = pinned(&c["original_v2_checker_source"])?;
    require(
        receipt["checker_source_sha256"] == checker_source.digest
            && receipt["checker_executable_sha256"] == checker_exe.digest
            && receipt["numeric_binding_helper_sha256"] == helper.digest
            && receipt["original_v2_checker_source_sha256"] == v2_source.digest,
        "checker source/executable/helper/v2 pins mismatch",
    )?;
    let runner = pinned(&c["original_runner_executable"])?;
    require(
        m["provenance"]["controller_executable_sha256"] == runner.digest
            && fs::canonicalize(text(&m["provenance"], "controller_executable")?)? == runner.path,
        "original frozen runner executable mismatch",
    )?;
    let failure_stderr = pinned(&c["original_collect_stderr"])?;
    let failure_stdout = pinned(&c["original_collect_stdout"])?;
    let failure_exit = pinned(&c["original_collect_exit"])?;
    require(
        std::str::from_utf8(&failure_exit.bytes)?.trim() == "2"
            && failure_stdout.bytes.is_empty()
            && std::str::from_utf8(&failure_stderr.bytes)?
                .contains("resume raw records do not reproduce the frozen observation"),
        "original frozen collect failure is not retained",
    )?;
    let seal = Seal::load(&receipt)?;
    seal.binding(&protocol.path, &protocol.digest)?;
    seal.binding(&manifest.path, &manifest.digest)?;
    seal.binding(&runner.path, &runner.digest)?;
    seal.exact_tree(&campaign)?;
    let cwd = PathBuf::from(text(&p, "cwd")?);
    let runner_source = cwd.join("research/n83_matched_20261014/native_runner");
    for (key, name) in [
        ("main", "src/main.rs"),
        ("monitor", "src/monitor.rs"),
        ("manifest", "Cargo.toml"),
        ("lock", "Cargo.lock"),
    ] {
        seal.binding(
            &runner_source.join(name),
            text(&m["provenance"]["compiled_source_sha256"], key)?,
        )?;
    }
    for (path, digest) in m["provenance"]["input_binary_bindings"]
        .as_object()
        .ok_or_else(|| err("frozen input/binary bindings missing"))?
    {
        seal.binding(
            Path::new(path),
            digest
                .as_str()
                .ok_or_else(|| err("invalid binary/input SHA"))?,
        )?;
    }
    for (path, digest) in m["generated_input_sha256"]
        .as_object()
        .ok_or_else(|| err("generated-input bindings missing"))?
    {
        seal.binding(
            Path::new(path),
            digest
                .as_str()
                .ok_or_else(|| err("invalid generated-input SHA"))?,
        )?;
    }
    let mut protected: BTreeSet<_> = seal.bytes.keys().cloned().collect();
    for b in [
        &audit,
        &protocol,
        &manifest,
        &checker_source,
        &checker_exe,
        &helper,
        &v2_source,
        &runner,
        &failure_stderr,
        &failure_stdout,
        &failure_exit,
    ] {
        protected.insert(b.path.clone());
    }
    protected.insert(config_path.clone());
    let output = output_path(out, &campaign, &protected)?;
    let provenance_output = output_path(provenance_out, &campaign, &protected)?;
    require(
        output != provenance_output,
        "output and provenance must differ",
    )?;
    let audit_rows = receipt["rows"].as_array().unwrap();
    let mut audit_map = BTreeMap::new();
    for row in audit_rows {
        require(
            row["process_status"] == "COMPLETED_VERIFIED"
                && row["row_present"] == true
                && row["answer_present"] == true
                && row["independent_scalar_replay"] == true
                && row["sidecar_checked_after_answer"] == true
                && row["complete_verified_row"] == true
                && row["error"].is_null(),
            "independent row gate incomplete",
        )?;
        if row["arm"] == "ic" {
            require(
                row["IC_four_point_witness_replay"] == true,
                "IC witness gate incomplete",
            )?;
        }
        let key = (
            text(row, "run_id")?.to_string(),
            text(row, "fixture_id")?.to_string(),
        );
        require(
            audit_map.insert(key, row).is_none(),
            "duplicate independently audited row",
        )?;
    }
    let expected = matrix();
    let declared = m["ordered_runs"]
        .as_array()
        .ok_or_else(|| err("frozen ordered runs missing"))?;
    require(declared.len() == 32, "frozen run count mismatch")?;
    let mut processes = Vec::new();
    let mut targets = Vec::new();
    let mut auxiliary = Vec::new();
    let mut extracted_audits = BTreeSet::new();
    let mut observation_pins = BTreeMap::new();
    for (r, declared) in expected.iter().zip(declared) {
        matrix_entry(r, declared)?;
        let dir = campaign.join("runs").join(&r.id);
        require(
            fs::canonicalize(text(declared, "out_dir")?)? == fs::canonicalize(&dir)?,
            "run output path mismatch",
        )?;
        for name in [
            "started.json",
            "execution.json",
            "observation.json",
            "stdout.jsonl",
            "stderr.txt",
            "heartbeat.jsonl",
            "targets.jsonl",
        ] {
            seal.canonical(&dir.join(name))?;
        }
        let observation_path = seal.canonical(&dir.join("observation.json"))?;
        let mut o = seal.json(&observation_path)?;
        observation_pins.insert(
            observation_path.to_string_lossy().to_string(),
            seal.hashes[&observation_path].clone(),
        );
        require(
            o["run_id"] == r.id
                && o["kind"] == r.kind
                && o["arm"] == r.arm
                && o["cohort_id"].as_u64() == Some(r.cohort as u64)
                && o["pair_id"].as_u64() == Some(r.pair as u64)
                && o["leg"] == r.leg
                && o["target_ids"] == json!(r.ids)
                && o["target_count"].as_u64() == Some(r.ids.len() as u64)
                && o["campaign_manifest_sha256"] == manifest.digest,
            "sealed observation identity/matrix mismatch",
        )?;
        require(
            o["status"] == "COMPLETED_VERIFIED"
                && o["all_verified"] == true
                && o["parent_known_answer_match_all"] == true
                && o["summary_ok"] == true
                && o["summary_counts_ok"] == true,
            "sealed observation is not already verified",
        )?;
        seal.binding(
            &dir.join("execution.json"),
            text(&o, "execution_receipt_sha256")?,
        )?;
        let mut rows = o["target_rows"]
            .as_array()
            .ok_or_else(|| err("sealed target rows missing"))?
            .clone();
        require(
            rows.len() == r.ids.len(),
            "sealed target-row count mismatch",
        )?;
        for (i, row) in rows.iter_mut().enumerate() {
            require(
                row["run_id"] == r.id
                    && row["target_id"] == r.ids[i]
                    && row["arm"] == r.arm
                    && row["kind"] == r.kind
                    && row["cohort_id"].as_u64() == Some(r.cohort as u64)
                    && row["pair_id"].as_u64() == Some(r.pair as u64)
                    && row["leg"] == r.leg
                    && row["status"] == "COMPLETED_VERIFIED"
                    && row["verified"] == true,
                "sealed target identity/status mismatch",
            )?;
            let key = (r.id.clone(), r.ids[i].clone());
            let audited = audit_map
                .get(&key)
                .ok_or_else(|| err("target absent from independent row seal"))?;
            require(audited["arm"] == r.arm, "independent target arm mismatch")?;
            require(extracted_audits.insert(key), "duplicate extracted target")?;
            if r.arm == "rho" {
                require(
                    row["seed"].as_u64() == Some(r.seeds[i])
                        && row["rawtargetrecord"]["walk_seed"].as_u64() == Some(r.seeds[i]),
                    "sealed target seed mismatch",
                )?;
            } else {
                require(row["rank_seed"] == 29, "sealed target rank seed mismatch")?;
            }
            row["independently_verified"] = json!(true);
        }
        if r.kind == "cohort" {
            targets.extend(rows);
        } else {
            auxiliary.extend(rows);
        }
        o.as_object_mut()
            .ok_or_else(|| err("observation object missing"))?
            .remove("target_rows");
        o["seed"] = if r.arm == "rho" && r.seeds.len() == 1 {
            json!(r.seeds[0])
        } else {
            Value::Null
        };
        o["rank_seed"] = if r.arm == "ic" {
            json!(29)
        } else {
            Value::Null
        };
        o["independently_verified"] = json!(true);
        processes.push(o);
    }
    require(
        processes.len() == 32
            && targets.len() == 192
            && auxiliary.len() == 26
            && extracted_audits.len() == 218,
        "normalized output matrix counts mismatch",
    )?;
    let provenance = source_provenance(&config_path, &config_bytes)?;
    let normalized = json!({"schema":"n83_matched_statistics_input_v1","process_rows":processes,"target_rows":targets,"auxiliary_target_rows":auxiliary,"independent_verification":receipt,"metadata":{"target_ids":ids(),"cold_target_ids":ids().into_iter().take(8).collect::<Vec<_>>(),"campaign_manifest_sha256":manifest.digest,"protocol_sha256":protocol.digest,"independent_replay_required":true,"independent_replay_receipt":audit.path,"independent_replay_receipt_sha256":audit.digest,"independent_verification_receipt_sha256":audit.digest,"independently_verified":true,"requested_processes":32,"requested_target_answers":218,"target_count":32,"cohort_repeats":3,"cold_subset_count":8,"aa_pairs":5,"ICrank_seed":29,"common_resource_envelope":{"cpu":3,"os_threads":1,"memory_bytes":12884901888u64,"timeout_s":2100},"IConline_excludes_input_validation_that_rho_includes":true,"ICfive_phase_timer_residuals_retained":true,"no_common_calibrated_operation_unit":true,"collection_method":"separate_posthoc_sealed_observation_extraction","original_frozen_collect_succeeded":false,"collector_provenance":provenance,"original_collect_failure":{"stderr":c["original_collect_stderr"],"stdout":c["original_collect_stdout"],"exit":c["original_collect_exit"]}}});
    let mut output_bytes = serde_json::to_vec_pretty(&normalized)?;
    output_bytes.push(b'\n');
    let extraction = json!({"schema":"n83.matched-collect-only-receipt/v1","status":"PASS","collector_provenance":provenance,"independent_verification_receipt":c["independent_receipt"],"protocol":c["protocol"],"manifest":c["manifest"],"validated_snapshot_files":seal.bytes.len(),"sealed_observation_sha256":observation_pins,"process_rows":32,"cohort_target_rows":192,"auxiliary_target_rows":26,"independent_target_rows":218,"original_frozen_collect_succeeded":false,"original_collect_failure":{"stderr":c["original_collect_stderr"],"stdout":c["original_collect_stdout"],"exit":c["original_collect_exit"]},"output":output,"output_sha256":sha(&output_bytes)?,"input_numbers_remeasured_or_recomputed":false});
    let mut receipt_bytes = serde_json::to_vec_pretty(&extraction)?;
    receipt_bytes.push(b'\n');
    let mut out_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)?;
    out_file.write_all(&output_bytes)?;
    out_file.sync_all()?;
    fs::set_permissions(&output, fs::Permissions::from_mode(0o444))?;
    let mut receipt_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&provenance_output)?;
    receipt_file.write_all(&receipt_bytes)?;
    receipt_file.sync_all()?;
    fs::set_permissions(&provenance_output, fs::Permissions::from_mode(0o444))?;
    println!(
        "{}",
        json!({"status":"PASS","output":output,"output_sha256":sha(&output_bytes)?,"collector_receipt":provenance_output,"validated_snapshot_files":336,"process_rows":32,"cohort_target_rows":192,"auxiliary_target_rows":26})
    );
    Ok(())
}
fn cli() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    require(
        args.len() == 6,
        "usage: n83-matched-collect-only --config CONFIG --out NEW_JSON --receipt NEW_RECEIPT",
    )?;
    let mut options = BTreeMap::new();
    for pair in args.chunks_exact(2) {
        require(
            options.insert(pair[0].clone(), pair[1].clone()).is_none(),
            "duplicate option",
        )?;
    }
    let option = |key: &str| {
        options
            .get(key)
            .map(PathBuf::from)
            .ok_or_else(|| err(format!("missing {key}")))
    };
    extract(
        &option("--config")?,
        &option("--out")?,
        &option("--receipt")?,
    )
}
fn main() {
    if let Err(e) = cli() {
        eprintln!("collect-only error: {e}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scratch() -> PathBuf {
        let p = std::env::temp_dir().join(format!(
            "n83-collect-only-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&p).unwrap();
        p
    }
    fn gate() -> Value {
        json!({"kind":"native_n83_matched_answer_and_file_verification","status":"PASS","complete_panel_gate":true,"aa_semantic_gate":true,"planned_processes":32,"observed_processes":32,"completed_verified_processes":32,"planned_target_rows":218,"checked_target_rows":218,"independent_scalar_replays":218,"IC_scalar_replays":104,"rho_scalar_replays":114,"IC_four_point_witness_replays":104,"checked_aa_semantic_pairs":5,"invalid_rows":0,"missing_target_rows":0,"nonanswer_rows":0,"launch_failed_processes":0,"clean_environment_records_checked":32,"resource_envelopes_checked":32,"observation_execution_bindings_checked":32,"all_files_snapshotted":336,"gate_errors":[],"rows":vec![json!({});218]})
    }
    #[test]
    fn partial_or_missing_aa_receipt_is_rejected() {
        full_gate(&gate()).unwrap();
        for (key, value) in [
            ("status", json!("FAIL")),
            ("checked_target_rows", json!(217)),
            ("checked_aa_semantic_pairs", json!(4)),
            ("aa_semantic_gate", json!(false)),
            ("IC_four_point_witness_replays", json!(103)),
            ("invalid_rows", json!(1)),
        ] {
            let mut bad = gate();
            bad[key] = value;
            assert!(full_gate(&bad).is_err());
        }
    }
    #[test]
    fn changed_or_wrongly_pinned_file_is_rejected() {
        let d = scratch();
        let p = d.join("input");
        fs::write(&p, b"original").unwrap();
        let pin = json!({"path":p,"sha256":sha(b"original").unwrap()});
        pinned(&pin).unwrap();
        fs::write(&p, b"changed").unwrap();
        assert!(pinned(&pin).is_err());
        let mut bad = pin;
        bad["sha256"] = json!("00");
        assert!(pinned(&bad).is_err());
    }
    #[test]
    fn unsealed_observation_is_rejected() {
        let d = scratch();
        let p = d.join("observation.json");
        fs::write(&p, b"{}").unwrap();
        let seal = Seal {
            bytes: BTreeMap::new(),
            hashes: BTreeMap::new(),
        };
        assert!(seal.json(&p).is_err());
        assert!(seal.exact_tree(&d).is_err());
    }
    #[test]
    fn existing_or_campaign_output_is_rejected() {
        let d = scratch();
        let campaign = d.join("campaign");
        fs::create_dir(&campaign).unwrap();
        let existing = d.join("existing");
        fs::write(&existing, b"keep").unwrap();
        assert!(output_path(&existing, &campaign, &BTreeSet::new()).is_err());
        assert!(output_path(&campaign.join("new"), &campaign, &BTreeSet::new()).is_err());
        assert!(output_path(&d.join("new"), &campaign, &BTreeSet::new()).is_ok());
    }
    #[test]
    fn number_tokens_are_preserved_without_another_float_roundtrip() {
        let bytes=br#"{"a":0.00019999999949504857,"b":3723.0182520000003,"n":340282366920938463463374607431768211455}"#;
        let parsed = parse(bytes).unwrap();
        assert_eq!(serde_json::to_vec(&parsed).unwrap(), bytes);
    }
    #[test]
    fn exact_matrix_remains32processes218answers() {
        let runs = matrix();
        assert_eq!(runs.len(), 32);
        assert_eq!(runs.iter().map(|r| r.ids.len()).sum::<usize>(), 218);
        let mut bad = json!({"run_id":runs[0].id,"kind":"aa","arm":"rho","cohort_id":0,"pair_id":1,"leg":"A","target_ids":["target_01"],"walk_seeds":[202610146411u64],"rank_seed":null});
        matrix_entry(&runs[0], &bad).unwrap();
        bad["walk_seeds"] = json!([202610146412u64]);
        assert!(matrix_entry(&runs[0], &bad).is_err());
    }
}
