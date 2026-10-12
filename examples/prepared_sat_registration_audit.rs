//! Verify a published prepared SAT registration and archive as data; never execute it.
use crypto_lib::cryptanalysis::prepared_control_archive::{
    safe_path, verify_tar_inventory as inventory,
};
use crypto_lib::cryptanalysis::prepared_sat_control::{
    canonical_sha, sha256, ControlConfig, PREPARATION_SHA256,
};
use flate2::read::MultiGzDecoder;
use serde_json::{json, Value};
use std::{fs, io::Read, path::Path, process::ExitCode};

fn require(ok: bool, message: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}
fn read(path: &Path, limit: usize) -> Result<Vec<u8>, String> {
    let meta = fs::symlink_metadata(path).map_err(|e| e.to_string())?;
    require(
        meta.is_file() && meta.len() <= limit as u64,
        "not a bounded regular publication file",
    )?;
    let mut options = fs::OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NOFOLLOW);
    }
    let file = options.open(path).map_err(|e| e.to_string())?;
    require(
        file.metadata().map_err(|e| e.to_string())?.is_file(),
        "publication input is not regular",
    )?;
    let mut bytes = Vec::new();
    file.take(limit as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    require(
        bytes.len() <= limit && bytes.len() == meta.len() as usize,
        "publication file changed while reading",
    )?;
    Ok(bytes)
}
fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}

fn audit(dir: &Path, expected_seal: &str) -> Result<Value, String> {
    let publication: Value =
        serde_json::from_slice(&read(&dir.join("PUBLICATION.json"), 64 * 1024)?)
            .map_err(|e| e.to_string())?;
    require(
        publication["schema_version"] == 1
            && publication["stage"] == "preregistered-not-dispatched"
            && publication["registration_sha256"] == expected_seal,
        "publication scope or external seal differs",
    )?;
    let mut expected = serde_json::Map::new();
    let mut records = serde_json::Map::new();
    for name in [
        "registration.json",
        "seal.json",
        "config.json",
        "preparation.json",
    ] {
        let bytes = read(&dir.join(name), 2 * 1024 * 1024)?;
        let descriptor = desc(&bytes);
        require(
            publication["files"][name] == descriptor,
            "publication sidecar bytes differ",
        )?;
        expected.insert(name.into(), descriptor);
        records.insert(
            name.into(),
            serde_json::from_slice::<Value>(&bytes).map_err(|e| e.to_string())?,
        );
    }
    let reg = &records["registration.json"];
    require(
        canonical_sha(reg)? == expected_seal
            && records["seal.json"]["registration_sha256"] == expected_seal,
        "canonical registration seal differs",
    )?;
    require(
        reg["schema_version"] == 1
            && reg["scope"] == "disclosed-n17-native-prepared-control"
            && reg["status"] == "registered-not-dispatched"
            && reg["hardware"] == json!({"os":"macos","architecture":"aarch64"})
            && reg["fresh_paired_qualification"] == false
            && reg["headline_eligible"] == false
            && reg["online_speedup"].is_null()
            && reg["ordinary_queries_executed"] == 0
            && reg["source_commit"] == publication["source_commit"],
        "registration scope/claim differs",
    )?;
    let config: ControlConfig =
        serde_json::from_value(records["config.json"].clone()).map_err(|e| e.to_string())?;
    config.validate()?;
    require(
        config.descent_query_seed == 2026100302
            && config.export_nonce == 2026100303
            && config.conflict_budget == 1000000
            && config.max_queries == 8
            && config.exporter_timeout_ms == 30000
            && config.solver_timeout_ms == 60000
            && config.controller_timeout_ms == 900000,
        "configuration differs from the accepted control template",
    )?;
    require(
        canonical_sha(&records["preparation.json"])? == PREPARATION_SHA256,
        "whole preparation seal differs",
    )?;
    for name in ["config.json", "preparation.json"] {
        require(
            reg[name] == expected[name]["sha256"],
            "registration input hash differs",
        )?;
    }
    let files = reg["immutable_files"]
        .as_object()
        .ok_or("missing immutable inventory")?;
    require(
        publication["immutable_file_count"] == files.len(),
        "publication inventory count differs",
    )?;
    for (name, descriptor) in files {
        require(safe_path(name), "unsafe registered source path")?;
        expected.insert(format!("immutable/{name}"), descriptor.clone());
    }
    for (key, name) in [
        ("worker_sha256", "bin/prepared_sat_worker"),
        ("auditor_sha256", "bin/icprog"),
    ] {
        require(
            reg[key] == files[name]["sha256"] && publication[key] == reg[key],
            "registered executable pin differs",
        )?;
    }
    require(
        publication["archive"]["path"] == "capsule.tar.gz",
        "unexpected archive path",
    )?;
    let archive = read(&dir.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    require(
        desc(&archive)
            == json!({"bytes":publication["archive"]["bytes"],"sha256":publication["archive"]["sha256"]}),
        "archive bytes differ",
    )?;
    let mut tar = Vec::new();
    MultiGzDecoder::new(archive.as_slice())
        .take(320 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(
        tar.len() <= 320 * 1024 * 1024,
        "archive exceeds decompression limit",
    )?;
    let count = inventory(&tar, &Value::Object(expected))?;
    Ok(
        json!({"schema_version":1,"status":"PASS_NATIVE_REGISTRATION_PUBLICATION_CUSTODY",
        "registration_sha256":expected_seal,"source_commit":reg["source_commit"],"archive_sha256":sha256(&archive),
        "archive_regular_files":count,"immutable_files":files.len(),"worker_sha256":reg["worker_sha256"],
        "auditor_sha256":reg["auditor_sha256"],"native_solver_executions":0,"execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"online_speedup":null}),
    )
}
fn create(path: &Path, bytes: &[u8]) -> Result<(), String> {
    use std::io::Write;
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| e.to_string())?;
    out.write_all(bytes)
        .and_then(|()| out.sync_all())
        .map_err(|e| e.to_string())
}
/// Write sidecars for a thin-shell archive of an unconsumed capsule, then check it.
fn publish(capsule: &Path, dir: &Path) -> Result<Value, String> {
    require(
        !capsule.join("consumed.json").exists(),
        "registration already consumed; publication is audit-only",
    )?;
    let mut files = serde_json::Map::new();
    let reg_bytes = read(&capsule.join("registration.json"), 2 * 1024 * 1024)?;
    let reg: Value = serde_json::from_slice(&reg_bytes).map_err(|e| e.to_string())?;
    let seal = canonical_sha(&reg)?;
    for name in [
        "registration.json",
        "seal.json",
        "config.json",
        "preparation.json",
    ] {
        let bytes = read(&capsule.join(name), 2 * 1024 * 1024)?;
        files.insert(name.into(), desc(&bytes));
        create(&dir.join(name), &bytes)?;
    }
    let archive = read(&dir.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    let metadata = json!({"schema_version":1,"stage":"preregistered-not-dispatched",
        "registration_sha256":seal,"source_commit":reg["source_commit"],
        "immutable_file_count":reg["immutable_files"].as_object().ok_or("missing immutable inventory")?.len(),
        "worker_sha256":reg["worker_sha256"],"auditor_sha256":reg["auditor_sha256"],"files":files,
        "archive":{"path":"capsule.tar.gz","bytes":archive.len(),"sha256":sha256(&archive)},
        "native_solver_executions":0,"execution_admitted":false,"fresh_paired_qualification":false,
        "headline_eligible":false,"online_speedup":null});
    create(
        &dir.join("PUBLICATION.json"),
        &serde_json::to_vec_pretty(&metadata).map_err(|e| e.to_string())?,
    )?;
    audit(dir, &seal)
}
fn run() -> Result<(), String> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    require(
        args.len() == 6 && args[4] == "--out",
        "expected three explicit option/value pairs",
    )?;
    let result = match (args[0].as_str(), args[2].as_str()) {
        ("--publication", "--registration-sha256") => audit(Path::new(&args[1]), &args[3])?,
        ("--publish", "--publication") => publish(Path::new(&args[1]), Path::new(&args[3]))?,
        _ => return Err("use --publication DIR --registration-sha256 SHA or --publish CAPSULE --publication DIR; --out must be new".into()),
    };
    let bytes = serde_json::to_vec_pretty(&result).map_err(|e| e.to_string())?;
    create(Path::new(&args[5]), &bytes)?;
    println!("{}", String::from_utf8(bytes).map_err(|e| e.to_string())?);
    Ok(())
}
fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("registration custody: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn tar_file(name: &str, bytes: &[u8], kind: u8) -> Vec<u8> {
        let mut h = vec![0; 512];
        h[..name.len()].copy_from_slice(name.as_bytes());
        h[124..136].copy_from_slice(format!("{:011o}\0", bytes.len()).as_bytes());
        h[156] = kind;
        h[257..265].copy_from_slice(b"ustar\x0000");
        let sum: usize = h
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if (148..156).contains(&i) {
                    32
                } else {
                    b as usize
                }
            })
            .sum();
        h[148..156].copy_from_slice(format!("{sum:06o}\0 ").as_bytes());
        h.extend(bytes);
        h.resize(h.len().div_ceil(512) * 512, 0);
        h.extend(vec![0; 1024]);
        h
    }
    #[test]
    fn mutation_omission_duplicate_links_and_unsafe_paths_reject() {
        let expected = json!({"immutable/a":desc(b"frozen")});
        let good = tar_file("immutable/a", b"frozen", b'0');
        assert_eq!(inventory(&good, &expected).unwrap(), 1);
        let mut bad = good.clone();
        bad[512] ^= 1;
        assert!(inventory(&bad, &expected).is_err());
        assert!(inventory(&[0; 1024], &expected).is_err());
        let mut duplicate = good[..1024].to_vec();
        duplicate.extend(&good);
        assert!(inventory(&duplicate, &expected).is_err());
        for name in ["../escape", "/absolute", "immutable//a", "immutable/./a"] {
            assert!(inventory(&tar_file(name, b"frozen", b'0'), &expected).is_err());
        }
        assert!(inventory(&tar_file("immutable/a", b"frozen", b'2'), &expected).is_err());
        assert!(inventory(&good[..good.len() - 1], &expected).is_err());
        let mut trailer = good;
        trailer.push(1);
        assert!(inventory(&trailer, &expected).is_err());
    }
}
