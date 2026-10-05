//! Disclosed-input parity for the one-request prepared n17 exporter.
//! Run only after the ordinary timing panels release the shared busy lock.
//! This is a transport control, not an IC target solve or speed measurement.
//!
//! cargo run --locked --release --example ic_exporter_prestart_probe -- \
//!   ACCEPTED_EXPORTER PREPARED_EXPORTER PREPARED_SHA256 CONTROL_INDEX NEW_OUTPUT_DIR

#[path = "../src/bin/prepared_sat_worker/native.rs"]
#[allow(dead_code)]
mod native;

use crypto_lib::cryptanalysis::prepared_sat_control::sha256;
use serde_json::{json, Value};
use std::{
    env, fs,
    path::{Path, PathBuf},
    process::ExitCode,
};

struct Control {
    point: [u64; 2],
    anf_sha256: &'static str,
    cnf_sha256: &'static str,
    magma_sha256: &'static str,
}

const CONTROLS: [Control; 3] = [
    Control {
        point: [40991, 73355],
        anf_sha256: "a7867583d28fe9593fd17dabdf6545d5bf79b0fda08888e5172f051a3d918eda",
        cnf_sha256: "05b0913d27231ae82601b44f042765cf9ba2727aeb66543467e92127b3a2d92f",
        magma_sha256: "1e7609332d721678d763b1eb48658b505eae39506b0f6e4d0edff8261595eaab",
    },
    Control {
        point: [73003, 104622],
        anf_sha256: "eab39a019e5298286e92169efcbe0eb381ce77af2d46749748d0ef202e1478ff",
        cnf_sha256: "03a017d1763f5be128f8069f78158723a456b0803757eaa03ffbb647d9c6d960",
        magma_sha256: "c5f7baa8ef80be4d749fb3be5ac6c3d957978166d534826090455378ec6e5d7e",
    },
    Control {
        point: [59775, 2910],
        anf_sha256: "1783698ff0fda61fc056b24b19d12e421f895c03ebbbc63fd399390b82def031",
        cnf_sha256: "2077d013463e6a58270a38f29752ec0395099d1db0abfcfaadb7f6b8970d784c",
        magma_sha256: "958701bb9e154d2bea03b83e8243aa2b5dfc90278f421ddf3e21ac1d88159a0f",
    },
];

fn expected_files(control: &Control, files: &Value) -> bool {
    files["instance.anf"]["sha256"] == control.anf_sha256
        && files["instance.xor.cnf"]["sha256"] == control.cnf_sha256
        && files["instance.magma"]["sha256"] == control.magma_sha256
}

fn deterministic_manifest(value: &Value) -> Result<Value, String> {
    let mut copy = value.clone();
    let object = copy.as_object_mut().ok_or("invalid exporter manifest")?;
    object
        .remove("timing_ns")
        .ok_or("missing exporter timing")?;
    object.remove("transport");
    Ok(copy)
}

fn run(
    accepted: &Path,
    prepared: &Path,
    prepared_sha256: &str,
    control_index: usize,
    out: &Path,
) -> Result<Value, String> {
    native::enforce_hardware()?;
    if ![accepted, prepared, out]
        .iter()
        .all(|path| path.is_absolute())
    {
        return Err("all executable and output paths must be absolute".into());
    }
    let control = CONTROLS
        .get(control_index)
        .ok_or("control index must be 0, 1 or 2")?;
    if sha256(&native::read(accepted, 2 * 1024 * 1024)?) != native::EXPORTER_SHA {
        return Err("accepted exporter binary pin differs".into());
    }
    if prepared_sha256.len() != 64
        || !prepared_sha256
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        || sha256(&native::read(prepared, 8 * 1024 * 1024)?) != prepared_sha256
    {
        return Err("prepared exporter binary pin differs".into());
    }
    let accepted_dir = out.join("accepted-instance");
    let prepared_dir = out.join("prepared-instance");
    let ledger = out.join("child-pids");
    let id = format!("native-control-{control_index:02}");
    let static_args = |directory: &Path| {
        vec![
            "17".into(),
            "6".into(),
            "standard".into(),
            "2026100303".into(),
            "1000000".into(),
            directory.to_string_lossy().into_owned(),
            "1".into(),
            "0".into(),
        ]
    };
    let mut accepted_args = static_args(&accepted_dir);
    accepted_args.extend([
        "--target-x".into(),
        control.point[0].to_string(),
        "--target-y".into(),
        control.point[1].to_string(),
        "--blind-instance-id".into(),
        id.clone(),
        "--export-only".into(),
    ]);
    let accepted_run = native::measured_child(
        accepted,
        &accepted_args,
        out,
        &out.join("accepted"),
        30_000,
        &ledger,
        false,
    )?;
    native::require(
        accepted_run.exit_code == Some(0) && !accepted_run.timed_out,
        "accepted exporter control failed",
    )?;
    let prepared_args = {
        let mut args = static_args(&prepared_dir);
        args.extend(["--export-only".into(), "--prestart-stdin".into()]);
        args
    };
    let environment = [("LC_ALL".to_string(), "C".to_string())];
    let child = native::PreparedChild::launch(native::PreparedLaunch {
        program: prepared,
        expected_sha256: prepared_sha256,
        args: &prepared_args,
        cwd: out,
        stem: &out.join("prepared"),
        ledger: &ledger,
        environment: &environment,
        marker: "c EXPORTER_PREPARED_STDIN_READY_v1",
        ready_deadline_ms: 30_000,
    })?;
    let request = serde_json::to_vec(&json!({
        "target_x":control.point[0].to_string(),
        "target_y":control.point[1].to_string(),
        "blind_instance_id":id,
    }))
    .map_err(|e| e.to_string())?;
    let mut request = request;
    request.push(b'\n');
    let prepared_run = child.deliver(&request, 30_000)?;
    native::require(
        prepared_run.exit_code == Some(0) && !prepared_run.timed_out,
        "prepared exporter control failed",
    )?;
    let (accepted_manifest, accepted_anf, accepted_cnf, accepted_files) =
        native::validate_exports(&accepted_dir, control.point)?;
    let (prepared_manifest, prepared_anf, prepared_cnf, prepared_files) =
        native::validate_exports(&prepared_dir, control.point)?;
    let source_bytes_match = accepted_anf == prepared_anf
        && accepted_cnf == prepared_cnf
        && accepted_files == prepared_files;
    let fixed_source_hashes_match =
        expected_files(control, &accepted_files) && expected_files(control, &prepared_files);
    let manifest_match =
        deterministic_manifest(&accepted_manifest)? == deterministic_manifest(&prepared_manifest)?;
    let source_instance_match =
        accepted_manifest["source_instance"] == prepared_manifest["source_instance"];
    let ready = native::load(&out.join("prepared.ready.json"))?;
    let ready_before_delivery = ready["state"] == "ready-without-input"
        && ready["stdin_written_bytes"] == 0
        && prepared_run.receipt["stdin_written_bytes"] == request.len();
    let mut negative_controls = Vec::new();
    let mut oversized = vec![b' '; 511];
    oversized.push(b'\n');
    for (name, payload) in [("malformed", b"{bad}\n".to_vec()), ("oversized", oversized)] {
        let rejected_dir = out.join(format!("reject-{name}-instance"));
        let mut args = static_args(&rejected_dir);
        args.extend(["--export-only".into(), "--prestart-stdin".into()]);
        let rejected = native::PreparedChild::launch(native::PreparedLaunch {
            program: prepared,
            expected_sha256: prepared_sha256,
            args: &args,
            cwd: out,
            stem: &out.join(format!("reject-{name}")),
            ledger: &ledger,
            environment: &environment,
            marker: "c EXPORTER_PREPARED_STDIN_READY_v1",
            ready_deadline_ms: 30_000,
        })?
        .deliver(&payload, 30_000)?;
        let rejected_as_expected = rejected.exit_code.is_some_and(|code| code != 0)
            && !rejected.timed_out
            && !rejected_dir.join("manifest.json").exists();
        negative_controls.push(json!({
            "name":name,
            "request_bytes":payload.len(),
            "request_sha256":sha256(&payload),
            "rejected_as_expected":rejected_as_expected,
            "receipt":rejected.receipt
        }));
    }
    let negative_controls_pass = negative_controls
        .iter()
        .all(|row| row["rejected_as_expected"] == true);
    let unused_dir = out.join("cancel-unused-instance");
    let mut unused_args = static_args(&unused_dir);
    unused_args.extend(["--export-only".into(), "--prestart-stdin".into()]);
    let unused = native::PreparedChild::launch(native::PreparedLaunch {
        program: prepared,
        expected_sha256: prepared_sha256,
        args: &unused_args,
        cwd: out,
        stem: &out.join("cancel-unused"),
        ledger: &ledger,
        environment: &environment,
        marker: "c EXPORTER_PREPARED_STDIN_READY_v1",
        ready_deadline_ms: 30_000,
    })?
    .cancel()?;
    let unused_cancel_pass = unused["state"] == "cancelled-unused"
        && unused["stdin_written_bytes"] == 0
        && unused["process_group_drain_confirmed"] == true
        && !unused_dir.join("manifest.json").exists();
    native::drain_ledger(&ledger)?;
    let pass = source_bytes_match
        && fixed_source_hashes_match
        && manifest_match
        && source_instance_match
        && ready_before_delivery
        && negative_controls_pass
        && unused_cancel_pass;
    Ok(json!({
        "schema_version":1,
        "question":"disclosed-prepared-exporter-parity-v1",
        "control_index":control_index,
        "public_point":control.point,
        "accepted_exporter_sha256":native::EXPORTER_SHA,
        "prepared_exporter_sha256":prepared_sha256,
        "accepted_files":accepted_files,"prepared_files":prepared_files,
        "accepted_receipt":accepted_run.receipt,
        "prepared_receipt":prepared_run.receipt,
        "prepared_ready_receipt":ready,
        "source_bytes_match":source_bytes_match,
        "fixed_source_hashes_match":fixed_source_hashes_match,
        "deterministic_manifest_match":manifest_match,
        "source_instance_match":source_instance_match,
        "ready_before_delivery":ready_before_delivery,
        "negative_controls":negative_controls,
        "negative_controls_pass":negative_controls_pass,
        "unused_cancel_receipt":unused,
        "unused_cancel_pass":unused_cancel_pass,
        "prepared_child_groups_drained":true,
        "transport_control_pass":pass,
        "source_bound_target_admitted":false,
        "online_speedup":null
    }))
}

fn main() -> ExitCode {
    let args = env::args_os().collect::<Vec<_>>();
    if args.len() != 6 {
        eprintln!("usage: ic_exporter_prestart_probe ACCEPTED_EXPORTER PREPARED_EXPORTER PREPARED_SHA256 CONTROL_INDEX NEW_OUTPUT_DIR");
        return ExitCode::FAILURE;
    }
    let out = PathBuf::from(&args[5]);
    if let Err(error) = fs::create_dir(&out) {
        eprintln!("new output directory required: {error}");
        return ExitCode::FAILURE;
    }
    let index = args[4].to_string_lossy().parse::<usize>();
    let result = index.map_err(|e| e.to_string()).and_then(|index| {
        run(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &args[3].to_string_lossy(),
            index,
            &out,
        )
    });
    let record = match &result {
        Ok(value) => value.clone(),
        Err(error) => json!({
            "schema_version":1,
            "question":"disclosed-prepared-exporter-parity-v1",
            "failure":error,
            "transport_control_pass":false,
            "source_bound_target_admitted":false,
            "online_speedup":null
        }),
    };
    if let Err(error) = native::save(&out.join("result.json"), &record) {
        eprintln!("cannot retain exporter control: {error}");
        return ExitCode::FAILURE;
    }
    println!("{}", out.join("result.json").display());
    if result
        .as_ref()
        .is_ok_and(|value| value["transport_control_pass"] == true)
    {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
