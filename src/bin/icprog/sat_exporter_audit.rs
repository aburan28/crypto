//! Independent, data-only parity audit for the exporter built by a sealed SAT
//! target validation freeze. Disclosed controls never estimate natural yield.
use super::{sat_control::native, sat_target_build::capsule, sat_target_custody, target_math};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{read, require, save};
use serde_json::{json, Value};
use std::{collections::BTreeSet, fs, path::Path};

const PROBE_SHA: &str = "a9005ea75e03d35e251880166f0f69bb59b73b6726ae042205eff7b634b23f5e";
const CONTROLS: [([u64; 2], [&str; 3]); 3] = [
    (
        [40991, 73355],
        [
            "a7867583d28fe9593fd17dabdf6545d5bf79b0fda08888e5172f051a3d918eda",
            "05b0913d27231ae82601b44f042765cf9ba2727aeb66543467e92127b3a2d92f",
            "1e7609332d721678d763b1eb48658b505eae39506b0f6e4d0edff8261595eaab",
        ],
    ),
    (
        [73003, 104622],
        [
            "eab39a019e5298286e92169efcbe0eb381ce77af2d46749748d0ef202e1478ff",
            "03a017d1763f5be128f8069f78158723a456b0803757eaa03ffbb647d9c6d960",
            "c5f7baa8ef80be4d749fb3be5ac6c3d957978166d534826090455378ec6e5d7e",
        ],
    ),
    (
        [59775, 2910],
        [
            "1783698ff0fda61fc056b24b19d12e421f895c03ebbbc63fd399390b82def031",
            "2077d013463e6a58270a38f29752ec0395099d1db0abfcfaadb7f6b8970d784c",
            "958701bb9e154d2bea03b83e8243aa2b5dfc90278f421ddf3e21ac1d88159a0f",
        ],
    ),
];
const MARKER: &str = "c EXPORTER_PREPARED_STDIN_READY_v1";

fn load(path: &Path) -> Result<Value, String> {
    target_math::parse(&read(path, 16 * 1024 * 1024)?)
}
fn without_runtime(manifest: &Value) -> Result<Value, String> {
    let mut copy = manifest.clone();
    let object = copy
        .as_object_mut()
        .ok_or("exporter manifest is not an object")?;
    require(
        object.remove("timing_ns").is_some(),
        "exporter timing missing",
    )?;
    object.remove("transport");
    Ok(copy)
}
fn streams(dir: &Path, stem: &str, receipt: &Value) -> Result<(), String> {
    for stream in ["stdout", "stderr"] {
        let bytes = read(&dir.join(format!("{stem}.{stream}")), 16 * 1024 * 1024)?;
        require(
            receipt[format!("{stream}_bytes")] == bytes.len()
                && receipt[format!("{stream}_sha256")] == sha256(&bytes),
            "exporter control stream differs from retained receipt",
        )?;
    }
    require(
        load(&dir.join(format!("{stem}.receipt.json")))? == *receipt,
        "exporter control receipt differs from result",
    )
}
fn ready(dir: &Path, stem: &str, pin: &str) -> Result<Value, String> {
    let receipt = load(&dir.join(format!("{stem}.ready.json")))?;
    let stdout = read(&dir.join(format!("{stem}.stdout")), 16 * 1024 * 1024)?;
    let prefix = format!("{MARKER}\n");
    require(
        stdout.starts_with(prefix.as_bytes())
            && stdout
                .split(|&b| b == b'\n')
                .filter(|line| *line == MARKER.as_bytes())
                .count()
                == 1
            && receipt["schema_version"] == 1
            && receipt["state"] == "ready-without-input"
            && receipt["marker"] == MARKER
            && receipt["stdin_written_bytes"] == 0
            && receipt["executable_sha256"] == pin
            && receipt["stdout_at_ready_sha256"] == sha256(prefix.as_bytes()),
        "prepared exporter lacks an independently retained zero-input READY boundary",
    )?;
    Ok(receipt)
}
fn argv(dir: &Path, executable: &Path, role: &str, index: usize, point: [u64; 2]) -> Value {
    let mut out = vec![
        executable.to_string_lossy().into_owned(),
        "17".into(),
        "6".into(),
        "standard".into(),
        "2026100303".into(),
        "1000000".into(),
        dir.join(format!("{role}-instance"))
            .to_string_lossy()
            .into_owned(),
        "1".into(),
        "0".into(),
    ];
    if role == "accepted" {
        out.extend([
            "--target-x".into(),
            point[0].to_string(),
            "--target-y".into(),
            point[1].to_string(),
            "--blind-instance-id".into(),
            format!("native-control-{index:02}"),
            "--export-only".into(),
        ]);
    } else {
        out.extend(["--export-only".into(), "--prestart-stdin".into()]);
    }
    json!(out)
}
fn audit_control(
    dir: &Path,
    index: usize,
    accepted: &Path,
    prepared: &Path,
    prepared_pin: &str,
) -> Result<Value, String> {
    let (point, hashes) = CONTROLS[index];
    let row = load(&dir.join("result.json"))?;
    require(
        row["schema_version"] == 1
            && row["question"] == "disclosed-prepared-exporter-parity-v1"
            && row["control_index"] == index
            && row["public_point"] == json!(point)
            && row["accepted_exporter_sha256"] == native::EXPORTER_SHA
            && row["prepared_exporter_sha256"] == prepared_pin
            && row["source_bound_target_admitted"] == false
            && row["online_speedup"].is_null(),
        "exporter control identity or narrow claim differs",
    )?;
    let (accepted_manifest, accepted_anf, accepted_cnf, accepted_files) =
        native::validate_exports(&dir.join("accepted-instance"), point)?;
    let (prepared_manifest, prepared_anf, prepared_cnf, prepared_files) =
        native::validate_exports(&dir.join("prepared-instance"), point)?;
    for role in ["accepted", "prepared"] {
        load(&dir.join(format!("{role}-instance/manifest.json")))?;
    }
    require(
        accepted_anf == prepared_anf
            && accepted_cnf == prepared_cnf
            && accepted_files == prepared_files
            && row["accepted_files"] == accepted_files
            && row["prepared_files"] == prepared_files
            && without_runtime(&accepted_manifest)? == without_runtime(&prepared_manifest)?
            && accepted_manifest["source_instance"] == prepared_manifest["source_instance"],
        "exporter source files or deterministic manifests differ",
    )?;
    for (name, hash) in [
        ("instance.anf", hashes[0]),
        ("instance.xor.cnf", hashes[1]),
        ("instance.magma", hashes[2]),
    ] {
        require(
            accepted_files[name]["sha256"] == hash,
            "exporter fixed disclosed source hash differs",
        )?;
    }
    let mut pids = BTreeSet::new();
    for (role, program, pin) in [
        ("accepted", accepted, native::EXPORTER_SHA),
        ("prepared", prepared, prepared_pin),
    ] {
        let receipt = &row[format!("{role}_receipt")];
        streams(dir, role, receipt)?;
        require(
            receipt["argv"] == argv(dir, program, role, index, point)
                && receipt["cwd"] == json!(dir)
                && receipt["environment"] == json!({"LC_ALL":"C"})
                && receipt["executable_sha256_before"] == pin
                && receipt["executable_sha256_after"] == pin
                && receipt["exit_code"] == 0
                && receipt["timed_out"] == false
                && receipt["process_group_drain_confirmed"] == true,
            "exporter role execution receipt differs",
        )?;
        let pid = receipt["pid"].as_u64().ok_or("exporter PID missing")?;
        require(pids.insert(pid), "exporter role PID repeated")?;
    }
    let prepared_ready = ready(dir, "prepared", prepared_pin)?;
    require(
        row["prepared_ready_receipt"] == prepared_ready
            && prepared_ready["pid"] == row["prepared_receipt"]["pid"]
            && prepared_ready["argv"] == argv(dir, prepared, "prepared", index, point)
            && prepared_ready["cwd"] == json!(dir)
            && prepared_ready["environment"] == json!({"LC_ALL":"C"})
            && row["prepared_receipt"]["readiness_marker_count"] == 1
            && row["prepared_receipt"]["stdin_write_error"].is_null()
            && row["prepared_receipt"]["stdin_written_bytes"]
                .as_u64()
                .is_some_and(|n| n > 0),
        "exporter positive request lacks a prepared input boundary",
    )?;
    let request = format!(
        "{{\"blind_instance_id\":\"native-control-{index:02}\",\"target_x\":\"{}\",\"target_y\":\"{}\"}}\n",
        point[0], point[1]
    );
    require(
        row["prepared_receipt"]["stdin_bytes"] == request.len()
            && row["prepared_receipt"]["stdin_sha256"] == sha256(request.as_bytes()),
        "prepared exporter request differs from disclosed point",
    )?;
    let negatives = row["negative_controls"]
        .as_array()
        .ok_or("exporter negative controls absent")?;
    require(negatives.len() == 2, "exporter negative role count differs")?;
    for (i, (name, payload)) in [
        ("malformed", b"{bad}\n".to_vec()),
        ("oversized", {
            let mut bytes = vec![b' '; 511];
            bytes.push(b'\n');
            bytes
        }),
    ]
    .into_iter()
    .enumerate()
    {
        let entry = &negatives[i];
        let receipt = &entry["receipt"];
        let stem = format!("reject-{name}");
        streams(dir, &stem, receipt)?;
        let readiness = ready(dir, &stem, prepared_pin)?;
        let destination = dir.join(format!("reject-{name}-instance"));
        require(
            entry["name"] == name
                && entry["request_bytes"] == payload.len()
                && entry["request_sha256"] == sha256(&payload)
                && entry["rejected_as_expected"] == true
                && receipt["pid"] == readiness["pid"]
                && receipt["argv"] == argv(dir, prepared, &stem, index, point)
                && readiness["argv"] == receipt["argv"]
                && receipt["cwd"] == json!(dir)
                && readiness["cwd"] == json!(dir)
                && receipt["environment"] == json!({"LC_ALL":"C"})
                && readiness["environment"] == receipt["environment"]
                && receipt["stdin_bytes"] == payload.len()
                && receipt["stdin_sha256"] == sha256(&payload)
                && receipt["executable_sha256_before"] == prepared_pin
                && receipt["executable_sha256_after"] == prepared_pin
                && receipt["exit_code"].as_i64().is_some_and(|code| code != 0)
                && receipt["timed_out"] == false
                && receipt["process_group_drain_confirmed"] == true
                && !destination.join("manifest.json").exists(),
            "prepared exporter accepted a malformed or oversized request",
        )?;
        require(
            pids.insert(receipt["pid"].as_u64().ok_or("negative PID missing")?),
            "exporter negative PID repeated",
        )?;
    }
    let unused = &row["unused_cancel_receipt"];
    streams(dir, "cancel-unused", unused)?;
    let unused_ready = ready(dir, "cancel-unused", prepared_pin)?;
    require(
        unused["pid"] == unused_ready["pid"]
            && unused_ready["argv"] == argv(dir, prepared, "cancel-unused", index, point)
            && unused_ready["cwd"] == json!(dir)
            && unused_ready["environment"] == json!({"LC_ALL":"C"})
            && unused["state"] == "cancelled-unused"
            && unused["stdin_written_bytes"] == 0
            && unused["executable_sha256_before"] == prepared_pin
            && unused["executable_sha256_after"] == prepared_pin
            && unused["process_group_drain_confirmed"] == true
            && !dir.join("cancel-unused-instance/manifest.json").exists()
            && pids.insert(unused["pid"].as_u64().ok_or("unused PID missing")?),
        "never-fed exporter child was not cancelled and drained",
    )?;
    let ledger =
        String::from_utf8(read(&dir.join("child-pids"), 65536)?).map_err(|e| e.to_string())?;
    let lines = ledger.lines().collect::<Vec<_>>();
    require(
        lines.len() == pids.len() * 2,
        "exporter PID ledger is incomplete",
    )?;
    let mut seen = BTreeSet::new();
    for line in lines.chunks_exact(2) {
        let pid = line[0]
            .strip_prefix("start ")
            .ok_or("exporter PID ledger start missing")?
            .parse::<u64>()
            .map_err(|e| e.to_string())?;
        require(
            pids.contains(&pid) && seen.insert(pid) && line[1] == format!("done {pid}"),
            "exporter PID ledger lacks matched completion",
        )?;
    }
    require(seen == pids, "exporter PID ledger omitted a role")?;
    require(
        row["source_bytes_match"] == true
            && row["fixed_source_hashes_match"] == true
            && row["deterministic_manifest_match"] == true
            && row["source_instance_match"] == true
            && row["ready_before_delivery"] == true
            && row["negative_controls_pass"] == true
            && row["unused_cancel_pass"] == true
            && row["prepared_child_groups_drained"] == true
            && row["transport_control_pass"] == true,
        "exporter producer row does not claim full disclosed parity",
    )?;
    Ok(json!({"control_index":index,"public_point":point,
        "result_sha256":sha256(&read(&dir.join("result.json"),16*1024*1024)?),
        "fixed_source_hashes_verified":true,"prepared_roles_audited":4,
        "accepted_role_audited":true,"native_executables_called_by_auditor":0}))
}

pub(super) fn audit(
    publication: &Path,
    seal: &str,
    controls: &Path,
    probe: &Path,
    out: &Path,
) -> Result<String, String> {
    let result = verify(publication, seal, controls, probe)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

pub(super) fn verify(
    publication: &Path,
    seal: &str,
    controls: &Path,
    probe: &Path,
) -> Result<Value, String> {
    let custody =
        sat_target_custody::verify(publication, seal, sat_target_custody::Kind::Validation)?;
    let record = capsule::registration(publication, seal)?;
    let header = load(&publication.join("PUBLICATION.json"))?;
    let original = Path::new(
        header["capsule_root"]
            .as_str()
            .ok_or("SAT capsule path missing")?,
    );
    let accepted = record
        .preparation
        .capsule
        .join("immutable/assets/bin/exporter");
    let prepared = original.join("immutable/assets/bin/prepared-exporter");
    require(
        record.validation_only
            && sha256(&read(&accepted, 8 * 1024 * 1024)?) == native::EXPORTER_SHA
            && sha256(&read(&prepared, 8 * 1024 * 1024)?) == record.prepared_exporter_sha256
            && record.immutable_files["source/examples/ic_exporter_prestart_probe.rs"].is_object()
            && sha256(&read(probe, 8 * 1024 * 1024)?) == PROBE_SHA,
        "SAT exporter parity source or exact build differs",
    )?;
    let names = fs::read_dir(controls)
        .map_err(|e| e.to_string())?
        .map(|entry| {
            let entry = entry.map_err(|e| e.to_string())?;
            require(
                entry.file_type().map_err(|e| e.to_string())?.is_dir(),
                "exporter control root contains a non-directory",
            )?;
            entry
                .file_name()
                .into_string()
                .map_err(|_| "exporter control name is not UTF-8".to_string())
        })
        .collect::<Result<BTreeSet<_>, _>>()?;
    require(
        names
            == ["control-00", "control-01", "control-02"]
                .into_iter()
                .map(str::to_string)
                .collect(),
        "exporter disclosed control set differs",
    )?;
    let before = native::inventory(controls)?;
    let rows = (0..3)
        .map(|index| {
            audit_control(
                &controls.join(format!("control-{index:02}")),
                index,
                &accepted,
                &prepared,
                &record.prepared_exporter_sha256,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    require(
        native::inventory(controls)? == before,
        "exporter parity evidence changed during audit",
    )?;
    let result = json!({"schema_version":1,"status":"PASS_DISCLOSED_EXACT_BUILT_EXPORTER_PARITY",
        "registration_sha256":seal,"source_manifest_sha256":record.source_manifest_sha256,
        "publication_archive_sha256":custody["archive_sha256"],
        "prepared_exporter_sha256":record.prepared_exporter_sha256,
        "accepted_exporter_sha256":native::EXPORTER_SHA,
        "probe_sha256":PROBE_SHA,"evidence_inventory_sha256":canonical_sha(&before)?,
        "audited_controls":rows,"source_bound_target_admitted":false,
        "fresh_paired_qualification":false,"native_executables_called_by_auditor":0,
        "online_speedup":null});
    Ok(result)
}
