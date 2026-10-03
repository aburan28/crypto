//! Bounded native transport. Scientific inputs are a sealed, one-use capsule.
use crypto_lib::cryptanalysis::prepared_sat_control::{
    canonical_sha, sha256, validate_manifest, ControlConfig, NativeOutput,
};
use serde_json::{json, Value};
#[cfg(unix)]
use std::os::unix::{
    fs::{OpenOptionsExt, PermissionsExt},
    process::CommandExt,
};
use std::{
    collections::BTreeMap,
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
    time::{Duration, Instant},
};

pub fn require(ok: bool, message: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}
pub fn read(path: &Path, limit: u64) -> Result<Vec<u8>, String> {
    let meta = fs::symlink_metadata(path).map_err(|e| format!("{}: {e}", path.display()))?;
    require(
        meta.is_file() && !meta.file_type().is_symlink() && meta.len() <= limit,
        "input is not a bounded regular file",
    )?;
    let mut opts = fs::OpenOptions::new();
    opts.read(true);
    #[cfg(unix)]
    opts.custom_flags(libc::O_NOFOLLOW);
    let file = opts.open(path).map_err(|e| e.to_string())?;
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    require(bytes.len() as u64 <= limit, "input exceeded limit")?;
    Ok(bytes)
}
pub fn load(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&read(path, 16 * 1024 * 1024)?).map_err(|e| e.to_string())
}
pub fn create(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    file.write_all(bytes)
        .and_then(|_| file.sync_all())
        .map_err(|e| e.to_string())?;
    #[cfg(unix)]
    fs::File::open(path.parent().ok_or("missing evidence parent")?)
        .and_then(|dir| dir.sync_all())
        .map_err(|e| e.to_string())?;
    Ok(())
}
pub fn save(path: &Path, value: &Value) -> Result<(), String> {
    create(
        path,
        &serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?,
    )
}
pub fn inventory(dir: &Path) -> Result<Value, String> {
    fn walk(root: &Path, dir: &Path, files: &mut BTreeMap<String, Value>) -> Result<(), String> {
        for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            let m = fs::symlink_metadata(&path).map_err(|e| e.to_string())?;
            require(
                !m.file_type().is_symlink(),
                "source tree contains a symlink",
            )?;
            if m.is_dir() {
                walk(root, &path, files)?
            } else {
                require(m.is_file(), "source tree contains a special file")?;
                let name = path
                    .strip_prefix(root)
                    .map_err(|e| e.to_string())?
                    .to_str()
                    .ok_or("source path is not UTF-8")?
                    .to_string();
                let bytes = read(&path, 256 * 1024 * 1024)?;
                files.insert(name, json!({"bytes":bytes.len(),"sha256":sha256(&bytes)}));
            }
        }
        Ok(())
    }
    require(
        fs::symlink_metadata(dir)
            .map_err(|e| e.to_string())?
            .is_dir(),
        "source root is not a directory",
    )?;
    let mut files = BTreeMap::new();
    walk(dir, dir, &mut files)?;
    Ok(json!(files))
}
pub fn check_tree(dir: &Path, expected: &Value) -> Result<(), String> {
    require(
        &inventory(dir)? == expected,
        "immutable source/binary tree differs",
    )
}
pub fn copy_tree(source: &Path, dest: &Path) -> Result<(), String> {
    fs::create_dir(dest).map_err(|e| e.to_string())?;
    for entry in fs::read_dir(source).map_err(|e| e.to_string())? {
        let path = entry.map_err(|e| e.to_string())?.path();
        let target = dest.join(path.file_name().ok_or("empty file name")?);
        let m = fs::symlink_metadata(&path).map_err(|e| e.to_string())?;
        require(
            !m.file_type().is_symlink(),
            "source copy contains a symlink",
        )?;
        if m.is_dir() {
            copy_tree(&path, &target)?
        } else {
            require(m.is_file(), "source copy contains a special file")?;
            create(&target, &read(&path, 256 * 1024 * 1024)?)?;
        }
    }
    Ok(())
}
pub const ARCHIVE_SHA: &str = "a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96";
const MANIFEST_SHA: &str = "e0c7a5edd01a6c64da810f2cf72a80b1b570d9031f96e8463cf1aaed2159e0f0";
pub const CMS_SHA: &str = "6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af";
pub const EXPORTER_SHA: &str = "4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a";
/// Extract only the already accepted archive. No ambient native binaries or download.
pub fn extract_assets(bundle: &Path, dest: &Path) -> Result<Value, String> {
    let bytes = read(&bundle.join("assets.tar.gz"), 8 * 1024 * 1024)?;
    require(
        bytes.len() == 7_292_102 && sha256(&bytes) == ARCHIVE_SHA,
        "accepted native archive differs",
    )?;
    let manifest = load(&bundle.join("manifest.json"))?;
    require(
        canonical_sha(&manifest)? == MANIFEST_SHA,
        "accepted native manifest differs",
    )?;
    // The complete archive itself is a fixed accepted pin. Reject all unsupported tar features.
    let mut tar = Vec::new();
    flate2::read::GzDecoder::new(bytes.as_slice())
        .take(40 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(tar.len() <= 40 * 1024 * 1024, "native tar exceeded limit")?;
    fs::create_dir(dest).map_err(|e| e.to_string())?;
    let mut offset = 0;
    let mut count = 0;
    while offset + 512 <= tar.len() {
        let h = &tar[offset..offset + 512];
        if h.iter().all(|&b| b == 0) {
            break;
        }
        let text = |b: &[u8]| {
            std::str::from_utf8(b)
                .map(|s| {
                    s.trim_matches(|c: char| c == '\0' || c.is_ascii_whitespace())
                        .to_string()
                })
                .map_err(|e| e.to_string())
        };
        let name = text(&h[..100])?;
        require(
            !name.is_empty()
                && name
                    .split('/')
                    .all(|p| !p.is_empty() && p != "." && p != "..")
                && h[345..500].iter().all(|&b| b == 0),
            "unsafe tar path/prefix",
        )?;
        require(
            h[156] == b'0' || h[156] == 0,
            "native tar member is not a regular file",
        )?;
        let expected = u64::from_str_radix(&text(&h[148..156])?, 8).map_err(|e| e.to_string())?;
        let checksum = h
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if (148..156).contains(&i) {
                    32
                } else {
                    b as u64
                }
            })
            .sum::<u64>();
        require(checksum == expected, "tar checksum differs")?;
        let size = usize::from_str_radix(&text(&h[124..136])?, 8).map_err(|e| e.to_string())?;
        offset += 512;
        require(size <= tar.len() - offset, "truncated tar member")?;
        let path = dest.join(&name);
        fs::create_dir_all(path.parent().ok_or("missing tar parent")?)
            .map_err(|e| e.to_string())?;
        create(&path, &tar[offset..offset + size])?;
        #[cfg(unix)]
        if name == "bin/cms" || name == "bin/exporter" {
            fs::set_permissions(&path, fs::Permissions::from_mode(0o755))
                .map_err(|e| e.to_string())?;
        }
        offset += size.div_ceil(512) * 512;
        count += 1;
    }
    require(
        count == 14 && tar[offset..].iter().all(|&b| b == 0),
        "native tar member count/trailer differs",
    )?;
    require(
        sha256(&read(&dest.join("bin/cms"), 4 * 1024 * 1024)?) == CMS_SHA
            && sha256(&read(&dest.join("bin/exporter"), 2 * 1024 * 1024)?) == EXPORTER_SHA,
        "native executable pin differs",
    )?;
    Ok(
        json!({"archive_sha256":ARCHIVE_SHA,"accepted_manifest_sha256":MANIFEST_SHA,"files":inventory(dest)?}),
    )
}
pub fn validate_exports(
    dir: &Path,
    point: [u64; 2],
) -> Result<(Value, String, String, Value), String> {
    let manifest = load(&dir.join("manifest.json"))?;
    validate_manifest(&manifest, point)?;
    let mut files = serde_json::Map::new();
    let mut contents = Vec::new();
    for (role, name) in [
        ("wdsat_anf", "instance.anf"),
        ("cryptominisat_xor_dimacs", "instance.xor.cnf"),
        ("magma_boolean_f4", "instance.magma"),
    ] {
        let bytes = read(&dir.join(name), 2 * 1024 * 1024)?;
        let descriptor = &manifest["exports"][role];
        require(
            descriptor["path"] == name
                && descriptor["bytes"] == bytes.len()
                && descriptor["blake3"] == blake3::hash(&bytes).to_hex().as_str(),
            "export metadata/hash differs",
        )?;
        files.insert(
            name.into(),
            json!({"bytes":bytes.len(),"sha256":sha256(&bytes)}),
        );
        contents.push(String::from_utf8(bytes).map_err(|e| e.to_string())?);
    }
    // Verify envelope/counts even for UNSAT/timeout: no unvalidated instance earns a status.
    for (i, n) in [
        51,
        manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
            .as_u64()
            .ok_or("missing model width")? as usize,
    ]
    .into_iter()
    .enumerate()
    {
        require((51..=4096).contains(&n), "source variable envelope differs")?;
        let text = &contents[i];
        let h = text
            .lines()
            .find(|l| !l.is_empty() && !l.starts_with('c'))
            .ok_or("missing source header")?;
        let tokens = h.split_whitespace().collect::<Vec<_>>();
        require(
            tokens.len() == 4
                && tokens[..2] == ["p", "cnf"]
                && tokens[2].parse::<usize>().ok() == Some(n),
            "source header differs from manifest",
        )?;
        let count = tokens[3].parse::<usize>().map_err(|e| e.to_string())?;
        let expected = if i == 0 {
            50
        } else {
            manifest["exports"]["cryptominisat_xor_dimacs"]["cnf_clauses"]
                .as_u64()
                .ok_or("missing clause count")? as usize
                + 50
        };
        require(
            count == expected
                && text
                    .lines()
                    .filter(|l| !l.is_empty() && !l.starts_with('c'))
                    .count()
                    == expected + 1,
            "source row count differs",
        )?;
    }
    crypto_lib::cryptanalysis::prepared_sat_control::validate_source_syntax(
        &contents[0],
        &contents[1],
        manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
            .as_u64()
            .ok_or("missing expanded width")? as usize,
    )?;
    Ok((
        manifest,
        contents.remove(0),
        contents.remove(0),
        json!(files),
    ))
}

/// The ledger is fsynced before a child leaves the controller's process group.
fn record_pid(ledger: &Path, pid: u32, event: &str) -> Result<(), String> {
    let mut file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(ledger)
        .map_err(|e| e.to_string())?;
    writeln!(file, "{event} {pid}")
        .and_then(|_| file.sync_all())
        .map_err(|e| e.to_string())
}
#[cfg(unix)]
fn kill_group(pid: u32) {
    if pid > 1 {
        unsafe {
            libc::kill(-(pid as i32), libc::SIGKILL);
        }
    }
}
#[cfg(not(unix))]
fn kill_group(_pid: u32) {}
#[cfg(unix)]
fn group_absent(pid: u32) -> bool {
    unsafe {
        libc::kill(-(pid as i32), 0) != 0
            && std::io::Error::last_os_error().raw_os_error() == Some(libc::ESRCH)
    }
}
#[cfg(not(unix))]
fn group_absent(_pid: u32) -> bool {
    false
}
fn confirm_drain(pid: u32) -> bool {
    kill_group(pid);
    let until = Instant::now() + Duration::from_secs(2);
    while !group_absent(pid) && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(5));
    }
    group_absent(pid)
}
pub fn drain_ledger(ledger: &Path) -> Result<(), String> {
    if ledger.exists() {
        let mut active = std::collections::BTreeSet::new();
        for line in std::str::from_utf8(&read(ledger, 64 * 1024)?)
            .map_err(|e| e.to_string())?
            .lines()
        {
            let (event, pid) = line.split_once(' ').ok_or("invalid PID ledger event")?;
            let pid = pid.parse::<u32>().map_err(|e| e.to_string())?;
            require(pid > 1, "unsafe child PID")?;
            match event {
                "start" => {
                    active.insert(pid);
                }
                "done" => {
                    active.remove(&pid);
                }
                _ => return Err("invalid PID ledger state".into()),
            }
        }
        for pid in active {
            require(
                confirm_drain(pid),
                "native child process group did not drain",
            )?;
        }
    }
    Ok(())
}
struct Guard {
    child: Child,
    ledger: PathBuf,
}
impl Drop for Guard {
    fn drop(&mut self) {
        kill_group(self.child.id());
        let _ = self.child.kill();
        let _ = self.child.wait();
        if group_absent(self.child.id()) {
            let _ = record_pid(&self.ledger, self.child.id(), "done");
        }
    }
}

/// Launch worker in its own group, or a helper which joins its own group only after GO.
pub fn measured_child(
    program: &Path,
    args: &[String],
    cwd: &Path,
    stem: &Path,
    deadline_ms: u64,
    ledger: &Path,
    helper: bool,
) -> Result<NativeOutput, String> {
    require(cfg!(unix), "native process watchdog requires Unix")?;
    let before = sha256(&read(program, 128 * 1024 * 1024)?);
    let stdout_path = stem.with_extension("stdout");
    let stderr_path = stem.with_extension("stderr");
    let stdout = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&stdout_path)
        .map_err(|e| e.to_string())?;
    let stderr = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&stderr_path)
        .map_err(|e| e.to_string())?;
    let mut command = Command::new(program);
    command
        .args(args)
        .current_dir(cwd)
        .env_clear()
        .env("LC_ALL", "C")
        .stdin(if helper {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(stdout)
        .stderr(stderr);
    #[cfg(unix)]
    if !helper {
        command.process_group(0);
    }
    let started = Instant::now();
    let mut guard = Guard {
        child: command.spawn().map_err(|e| e.to_string())?,
        ledger: ledger.to_path_buf(),
    };
    let pid = guard.child.id();
    record_pid(ledger, pid, "start")?;
    if helper {
        let mut stdin = guard.child.stdin.take().ok_or("missing helper handshake")?;
        stdin.write_all(b"GO\n").map_err(|e| e.to_string())?;
    }
    let mut timed_out = false;
    let mut output_limit = false;
    let status = loop {
        if let Some(status) = guard.child.try_wait().map_err(|e| e.to_string())? {
            break status;
        }
        output_limit = [&stdout_path, &stderr_path]
            .iter()
            .any(|p| fs::metadata(p).map_or(true, |m| m.len() > 8 * 1024 * 1024));
        if started.elapsed() >= Duration::from_millis(deadline_ms) || output_limit {
            timed_out = !output_limit;
            kill_group(pid);
            guard.child.kill().map_err(|e| e.to_string())?;
            break guard.child.wait().map_err(|e| e.to_string())?;
        }
        std::thread::sleep(Duration::from_millis(5));
    };
    let wall_ns = started.elapsed().as_nanos() as u64;
    // Always terminate descendants, including those retained after a normal parent exit.
    let drained = confirm_drain(pid);
    let after = sha256(&read(program, 128 * 1024 * 1024)?);
    let out = read(&stdout_path, 8 * 1024 * 1024)?;
    let err = read(&stderr_path, 8 * 1024 * 1024)?;
    let receipt = json!({"argv":std::iter::once(program.to_string_lossy().into_owned()).chain(args.iter().cloned()).collect::<Vec<_>>(),
        "cwd":cwd,"pid":pid,"deadline_ms":deadline_ms,"child_wall_ns":wall_ns,"exit_code":status.code(),
        "timed_out":timed_out,"output_limit":output_limit,"executable_sha256_before":before,"executable_sha256_after":after,
        "stdout_bytes":out.len(),"stdout_sha256":sha256(&out),"stderr_bytes":err.len(),"stderr_sha256":sha256(&err),
        "environment":{"LC_ALL":"C"},"process_group_drain_requested":true,"process_group_drain_confirmed":drained});
    save(&stem.with_extension("receipt.json"), &receipt)?;
    require(before == after, "child executable changed during execution")?;
    require(
        drained,
        "child process group did not drain; receipt retained",
    )?;
    Ok(NativeOutput {
        exit_code: status.code(),
        timed_out,
        stdout: String::from_utf8(out).map_err(|e| e.to_string())?,
        receipt,
    })
}
pub fn config(capsule: &Path) -> Result<ControlConfig, String> {
    let config: ControlConfig =
        serde_json::from_value(load(&capsule.join("config.json"))?).map_err(|e| e.to_string())?;
    config.validate()?;
    Ok(config)
}
pub fn check_capsule(capsule: &Path) -> Result<Value, String> {
    let seal = load(&capsule.join("seal.json"))?;
    let registration = load(&capsule.join("registration.json"))?;
    require(
        seal["registration_sha256"] == canonical_sha(&registration)?,
        "registration seal differs",
    )?;
    require(
        registration["schema_version"] == 1
            && registration["scope"] == "disclosed-n17-native-prepared-control"
            && registration["hardware"] == json!({"os":"macos","architecture":"aarch64"}),
        "capsule scope/hardware differs",
    )?;
    check_tree(&capsule.join("immutable"), &registration["immutable_files"])?;
    for name in ["config.json", "preparation.json"] {
        require(
            registration[name] == sha256(&read(&capsule.join(name), 16 * 1024 * 1024)?),
            "capsule input differs",
        )?;
    }
    config(capsule)?;
    Ok(registration)
}
pub fn enforce_hardware() -> Result<(), String> {
    require(
        std::env::consts::OS == "macos" && std::env::consts::ARCH == "aarch64",
        "retained native binaries are macOS ARM64 only; no fallback",
    )
}

/// Called inside a helper that cannot execute until its durable PID handshake completes.
pub fn helper(capsule: &Path, attempt: &Path, role: &str) -> Result<(), String> {
    enforce_hardware()?;
    let mut handshake = String::new();
    std::io::stdin()
        .take(4)
        .read_to_string(&mut handshake)
        .map_err(|e| e.to_string())?;
    require(handshake == "GO\n", "missing durable launch handshake")?;
    #[cfg(unix)]
    require(
        unsafe { libc::setpgid(0, 0) } == 0,
        "could not isolate native process group",
    )?;
    let cfg = config(capsule)?;
    let query = load(&attempt.join("query.json"))?;
    require(
        query["point"].as_array().is_some_and(|p| p.len() == 2),
        "missing registered query",
    )?;
    let coords = query["point"]
        .as_array()
        .ok_or("missing query point")?
        .iter()
        .map(|v| {
            v.as_u64()
                .filter(|&v| v < 1 << 17)
                .ok_or("invalid query coordinate")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let (program, args, pin) = match role {
        "exporter" => (
            capsule.join("immutable/assets/bin/exporter"),
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
                coords[0].to_string(),
                "--target-y".into(),
                coords[1].to_string(),
                "--blind-instance-id".into(),
                format!(
                    "native-control-{:02}",
                    query["trial"].as_u64().ok_or("invalid trial")?
                ),
                "--export-only".into(),
            ],
            EXPORTER_SHA,
        ),
        "cms" => (
            capsule.join("immutable/assets/bin/cms"),
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
            CMS_SHA,
        ),
        _ => return Err("unsupported native role".into()),
    };
    require(
        sha256(&read(&program, 4 * 1024 * 1024)?) == pin,
        "native role binary differs",
    )?;
    save(
        &attempt.join(format!("{role}.launch.json")),
        &json!({"role":role,"program":program,"argv":args,
        "binary_sha256_before_exec":pin,"config_sha256":canonical_sha(&json!(cfg))?,"query_sha256":canonical_sha(&query)?,
        "environment":{"LC_ALL":"C"},"cwd":attempt,"pid":std::process::id()}),
    )?;
    let mut command = Command::new(program);
    command
        .args(args)
        .current_dir(attempt)
        .env_clear()
        .env("LC_ALL", "C")
        .stdin(Stdio::null());
    #[cfg(unix)]
    {
        Err(command.exec().to_string())
    }
    #[cfg(not(unix))]
    Err("native helpers require Unix".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn temp(label: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("sat-native-{label}-{}", std::process::id()));
        fs::create_dir(&p).unwrap();
        p
    }
    #[test]
    fn create_only_and_tree_mutation_controls() {
        let p = temp("tree");
        create(&p.join("a"), b"accepted").unwrap();
        assert!(create(&p.join("a"), b"overwrite").is_err());
        let files = inventory(&p).unwrap();
        check_tree(&p, &files).unwrap();
        fs::write(p.join("a"), b"changed").unwrap();
        assert!(check_tree(&p, &files).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn retained_exports_have_exact_layout_and_hashes() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let dir=root.join("research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check");
        validate_exports(&dir, [62577, 27783]).unwrap();
        assert!(validate_exports(&dir, [62577, 27784]).is_err());
    }
    #[test]
    fn accepted_archive_custody_extracts_without_executing_binaries() {
        let p = temp("assets");
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let report=extract_assets(&root.join("research/ic_candidate_tournament_20260915/goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64"),&p.join("assets")).unwrap();
        assert_eq!(report["files"].as_object().unwrap().len(), 14);
        assert_eq!(report["archive_sha256"], ARCHIVE_SHA);
        fs::remove_dir_all(p).unwrap();
    }
    #[cfg(unix)]
    #[test]
    fn watchdog_preserves_exit_and_kills_a_timed_out_group() {
        let p = temp("watchdog");
        // Linux /bin/sh is often a symlink. Resolve only this non-scientific fixture;
        // registered executables must remain bounded regular files with their fixed pins.
        let shell = Path::new("/bin/sh").canonicalize().unwrap();
        let out = measured_child(
            &shell,
            &["-c".into(), "exit 3".into()],
            &p,
            &p.join("exit"),
            2000,
            &p.join("pids"),
            false,
        )
        .unwrap();
        assert_eq!(out.exit_code, Some(3));
        assert!(!out.timed_out);
        let out = measured_child(
            &shell,
            &["-c".into(), "sleep 30 & wait".into()],
            &p,
            &p.join("timeout"),
            30,
            &p.join("pids"),
            false,
        )
        .unwrap();
        assert!(out.timed_out);
        assert_eq!(out.exit_code, None);
        drain_ledger(&p.join("pids")).unwrap();
        fs::remove_dir_all(p).unwrap();
    }
}
