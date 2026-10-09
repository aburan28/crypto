//! Frozen-source one-target compact IC / strong-rho execution protocol.
//! Historical Python records remain archived; new runs use this native path.
#[allow(dead_code)] // Replay uses additional reference operations from this shared module.
mod reference;

use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use reference::{array, ensure, integer, point, sha256_file, BinaryCurve, Check};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const SUBDIR: &str = "experiments/koblitz-n53-compact-cold-20261009";
const SOURCE_FILES: &[&str] = &[
    "Cargo.toml",
    "Cargo.lock",
    "examples/koblitz_orbit_dlp_fast_online.rs",
    "examples/koblitz_rho_fixture.rs",
    "experiments/koblitz-n53-compact-cold-20261009/PROTOCOL.md",
    "experiments/koblitz-n53-compact-cold-20261009/native_run.rs",
    "experiments/koblitz-n53-compact-cold-20261009/native_replay.rs",
    "experiments/koblitz-n53-compact-cold-20261009/reference.rs",
    "experiments/koblitz-n53-compact-cold-20261009/target_points.jsonl",
    "experiments/koblitz-n53-compact-cold-20261009/workload.json",
];
const CAP_SECONDS: u64 = 900;
const CAP_RSS_BYTES: u64 = 16 * 1024 * 1024 * 1024;

fn path(name: &str) -> PathBuf {
    Path::new(ROOT).join(name)
}
fn protocol_path(name: &str) -> PathBuf {
    path(SUBDIR).join(name)
}

fn read_json(file: &Path) -> Check<Value> {
    serde_json::from_slice(&fs::read(file).map_err(|e| format!("{}: {e}", file.display()))?)
        .map_err(|e| format!("{}: {e}", file.display()))
}

fn write_new_json(file: &Path, value: &Value) -> Check<()> {
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(file)
        .map_err(|e| format!("{}: {e}", file.display()))?;
    let bytes = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    output.write_all(&bytes).map_err(|e| e.to_string())?;
    output.write_all(b"\n").map_err(|e| e.to_string())
}

fn git(args: &[&str]) -> Check<String> {
    let output = Command::new("git")
        .args(args)
        .current_dir(ROOT)
        .output()
        .map_err(|e| e.to_string())?;
    ensure(
        output.status.success(),
        format!(
            "git {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&output.stderr)
        ),
    )?;
    Ok(String::from_utf8(output.stdout)
        .map_err(|e| e.to_string())?
        .trim()
        .to_string())
}

fn assert_clean() -> Check<()> {
    ensure(
        git(&["status", "--porcelain", "--untracked-files=no"])?.is_empty(),
        "tracked source is dirty; commit or discard changes first",
    )
}

fn verify_workload() -> Check<Value> {
    let w = read_json(&protocol_path("workload.json"))?;
    let n = integer(&w, "n")? as u32;
    let a = integer(&w, "a")?;
    let curve = BinaryCurve::new(n, a, array(&w, "field_modulus_low_terms")?)?;
    let order = integer(&w, "subgroup_order")?;
    let seed = w["target_seed_ascii"]
        .as_str()
        .ok_or("missing seed ASCII")?;
    ensure(seed.is_ascii(), "target seed is not ASCII")?;
    let hash = BigUint::from_bytes_be(&sha256(seed.as_bytes()));
    let scalar = (hash % BigUint::from(order - 1) + BigUint::from(1u8))
        .to_u64()
        .ok_or("seed scalar overflow")?;
    let generator = point(&w["generator"])?;
    let target = curve.scale(scalar, generator)?;
    ensure(
        scalar == integer(&w, "verification_scalar")?,
        "seed scalar mismatch",
    )?;
    ensure(
        target == point(&w["primary_target"])?,
        "seed target mismatch",
    )?;
    ensure(
        curve.scale(order, generator)?.is_none(),
        "generator not in subgroup",
    )?;
    ensure(
        curve.on_curve(generator) && curve.on_curve(target),
        "workload point off curve",
    )?;
    ensure(
        integer(&w, "target_count")? == 1 && integer(&w, "warm_target_count")? == 0,
        "not a cold one-target workload",
    )?;
    let target_pair = target.ok_or("identity target")?;
    let text =
        fs::read_to_string(protocol_path("target_points.jsonl")).map_err(|e| e.to_string())?;
    ensure(
        text == format!("[{},{}]\n", target_pair.0, target_pair.1),
        "IC target input differs from workload",
    )?;
    Ok(w)
}

fn source_hashes() -> Check<BTreeMap<&'static str, String>> {
    SOURCE_FILES
        .iter()
        .map(|&name| Ok((name, sha256_file(&path(name))?)))
        .collect()
}

fn executable(file: &Path) -> Check<()> {
    ensure(
        file.is_file(),
        format!("executable missing: {}", file.display()),
    )?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        ensure(
            file.metadata()
                .map_err(|e| e.to_string())?
                .permissions()
                .mode()
                & 0o111
                != 0,
            format!("file not executable: {}", file.display()),
        )?;
    }
    Ok(())
}

fn freeze(ic: &Path, rho: &Path) -> Check<Value> {
    assert_clean()?;
    let w = verify_workload()?;
    executable(ic)?;
    executable(rho)?;
    Ok(json!({
        "schema_version":2, "kind":"n53_compact_cold_native_source_freeze",
        "source_commit":git(&["rev-parse","HEAD"])?,
        "source_sha256":source_hashes()?,
        "ic_binary_sha256":sha256_file(ic)?, "rho_binary_sha256":sha256_file(rho)?,
        "ic_argv_suffix":["construct:53:0:244","target_points.jsonl","20261009"],
        "rho_argv_suffix":["53","0","signed_frobenius","1","strong","20261009",integer(&w,"verification_scalar")?.to_string()],
        "wall_cap_seconds_per_arm":CAP_SECONDS,"rss_cap_bytes_per_arm":CAP_RSS_BYTES
    }))
}

fn verify_freeze(file: &Path, ic: &Path, rho: &Path) -> Check<Value> {
    let absolute = file.canonicalize().map_err(|e| e.to_string())?;
    let relative = absolute
        .strip_prefix(ROOT)
        .map_err(|_| "freeze file must be inside repository")?;
    let relative = relative.to_str().ok_or("non-UTF8 freeze path")?;
    let committed = Command::new("git")
        .args(["show", &format!("HEAD:{relative}")])
        .current_dir(ROOT)
        .output()
        .map_err(|e| e.to_string())?;
    ensure(committed.status.success(), "freeze file is not committed")?;
    ensure(
        committed.stdout == fs::read(&absolute).map_err(|e| e.to_string())?,
        "freeze differs from committed version",
    )?;
    let record: Value = serde_json::from_slice(&committed.stdout).map_err(|e| e.to_string())?;
    ensure(
        record["kind"] == "n53_compact_cold_native_source_freeze" && record["schema_version"] == 2,
        "unknown freeze schema",
    )?;
    assert_clean()?;
    git(&[
        "merge-base",
        "--is-ancestor",
        record["source_commit"]
            .as_str()
            .ok_or("missing source commit")?,
        "HEAD",
    ])?;
    ensure(
        record["source_sha256"] == json!(source_hashes()?),
        "source/input hash differs from freeze",
    )?;
    ensure(
        record["ic_binary_sha256"] == sha256_file(ic)?
            && record["rho_binary_sha256"] == sha256_file(rho)?,
        "executable hash differs from freeze",
    )?;
    ensure(
        record["wall_cap_seconds_per_arm"] == CAP_SECONDS
            && record["rss_cap_bytes_per_arm"] == CAP_RSS_BYTES,
        "resource cap differs from protocol",
    )?;
    let w = verify_workload()?;
    let frozen_scalar = record["rho_argv_suffix"][6]
        .as_str()
        .ok_or("missing frozen rho scalar")?
        .parse::<u64>()
        .map_err(|e| e.to_string())?;
    ensure(
        frozen_scalar == integer(&w, "verification_scalar")?,
        "rho scalar differs from freeze",
    )?;
    Ok(record)
}

fn rss_bytes(pid: u32) -> Check<u64> {
    let mut child = Command::new("/bin/ps")
        .args(["-o", "rss=", "-p", &pid.to_string()])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    let started = Instant::now();
    loop {
        if child.try_wait().map_err(|e| e.to_string())?.is_some() {
            break;
        }
        if started.elapsed() >= Duration::from_secs(2) {
            let _ = child.kill();
            let _ = child.wait();
            return Err("RSS monitor timed out after 2 seconds".into());
        }
        thread::sleep(Duration::from_millis(10));
    }
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    ensure(
        output.status.success(),
        format!(
            "RSS monitor failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ),
    )?;
    let kib = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse::<u64>()
        .map_err(|e| e.to_string())?;
    Ok(kib * 1024)
}

fn kill_group(pid: u32) {
    unsafe {
        libc::kill(-(pid as i32), libc::SIGKILL);
    }
}

fn wait4(pid: u32, options: i32) -> Check<Option<(i32, libc::rusage)>> {
    loop {
        let mut status = 0;
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
        let found = unsafe { libc::wait4(pid as i32, &mut status, options, usage.as_mut_ptr()) };
        if found < 0 {
            let error = std::io::Error::last_os_error();
            if error.kind() == std::io::ErrorKind::Interrupted {
                continue;
            }
            return Err(format!("wait4: {error}"));
        }
        if found == 0 {
            return Ok(None);
        }
        return Ok(Some((status, unsafe { usage.assume_init() })));
    }
}

fn run_arm(
    argv: &[String],
    stdout: &Path,
    stderr: &Path,
    environment: &BTreeMap<String, String>,
) -> Check<Value> {
    let started = Instant::now();
    let out = File::create(stdout).map_err(|e| e.to_string())?;
    let err = File::create(stderr).map_err(|e| e.to_string())?;
    let mut command = Command::new(&argv[0]);
    command
        .args(&argv[1..])
        .stdout(Stdio::from(out))
        .stderr(Stdio::from(err))
        .process_group(0);
    for (key, value) in environment {
        command.env(key, value);
    }
    let child = match command.spawn() {
        Ok(child) => child,
        Err(e) => {
            return Ok(
                json!({"argv":argv,"status":"launch_failure","exit_code":null,
            "wall_ms":started.elapsed().as_secs_f64()*1000.0,"sampled_peak_rss_bytes":0,
            "child_peak_rss_bytes":null,"monitor_error":e.to_string(),
            "stdout_sha256":sha256_file(stdout)?,"stderr_sha256":sha256_file(stderr)?}),
            )
        }
    };
    let pid = child.id();
    let mut peak = 0u64;
    let mut terminal = "exited";
    let mut monitor_error: Option<String> = None;
    let (status, usage) = loop {
        if let Some(done) = wait4(pid, libc::WNOHANG)? {
            break done;
        }
        if started.elapsed().as_secs_f64() >= CAP_SECONDS as f64 {
            terminal = "timeout";
            kill_group(pid);
            break wait4(pid, 0)?.ok_or("wait4 missing timed-out child")?;
        }
        match rss_bytes(pid) {
            Ok(value) => peak = peak.max(value),
            Err(e) => {
                if let Some(done) = wait4(pid, libc::WNOHANG)? {
                    break done;
                }
                terminal = "monitor_failure";
                monitor_error = Some(e);
                kill_group(pid);
                break wait4(pid, 0)?.ok_or("wait4 missing unmonitored child")?;
            }
        }
        if peak > CAP_RSS_BYTES {
            terminal = "memory_cap";
            kill_group(pid);
            break wait4(pid, 0)?.ok_or("wait4 missing capped child")?;
        }
        thread::sleep(Duration::from_millis(100));
    };
    let exit_code = if status & 0x7f == 0 {
        status >> 8
    } else {
        -(status & 0x7f)
    };
    let child_peak = (usage.ru_maxrss as u64) * (if cfg!(target_os = "macos") { 1 } else { 1024 });
    if terminal == "exited" && child_peak > CAP_RSS_BYTES {
        terminal = "memory_cap";
    }
    let outcome = if terminal == "exited" {
        if exit_code == 0 {
            "success"
        } else {
            "process_failure"
        }
    } else {
        terminal
    };
    Ok(json!({"argv":argv,"status":outcome,"exit_code":exit_code,
        "wall_ms":started.elapsed().as_secs_f64()*1000.0,"sampled_peak_rss_bytes":peak,
        "child_peak_rss_bytes":child_peak,"monitor_error":monitor_error,
        "stdout_sha256":sha256_file(stdout)?,"stderr_sha256":sha256_file(stderr)?}))
}

fn one_json_line(file: &Path, kind: &str) -> Check<Value> {
    reference::one(file, kind)
}

fn checked_phase_sum(values: &[f64], total: f64, label: &str) -> Check<()> {
    ensure(
        values.iter().all(|x| x.is_finite() && *x >= 0.0) && total.is_finite() && total >= 0.0,
        format!("nonfinite {label} time"),
    )?;
    ensure(
        (values.iter().sum::<f64>() - total).abs() <= (total * 1e-6).max(1e-6),
        format!("{label} phases do not sum to interval"),
    )
}

fn floating(v: &Value, key: &str) -> Check<f64> {
    v[key].as_f64().ok_or_else(|| format!("missing {key}"))
}

fn analyse(ic: &Value, rho: &Value, run_dir: &Path, w: &Value) -> Value {
    let mut result = json!({"status":"unverified","ic_online_ms":null,"rho_online_ms":null,
        "online_speedup":null,"independent_scalar_replay":false});
    if ic["status"] != "success" || rho["status"] != "success" {
        result["status"] = json!("arm_failure");
        return result;
    }
    let checked: Check<Value> = (|| {
        let target = point(&w["primary_target"])?;
        let scalar = integer(w, "verification_scalar")?;
        let ic_target =
            one_json_line(&run_dir.join("ic-target.jsonl"), "compact_orbit_dlp_target")?;
        let summary = one_json_line(
            &run_dir.join("ic.stdout.jsonl"),
            "compact_orbit_dlp_summary",
        )?;
        let rho_row = one_json_line(&run_dir.join("rho.stdout.jsonl"), "rho_public_fixture")?;
        let online = floating(&ic_target, "online_ms")?;
        let phases = [
            "target_query_ms",
            "target_pdp_ms",
            "target_relation_check_ms",
            "target_descent_ms",
            "target_recovery_check_ms",
        ]
        .iter()
        .map(|name| floating(&ic_target, name))
        .collect::<Check<Vec<_>>>()?;
        checked_phase_sum(&phases, online, "IC target")?;
        checked_phase_sum(
            &[floating(&ic_target, "target_phase_sum_ms")?],
            online,
            "IC target reported",
        )?;
        let cold = floating(&summary, "cold_in_process_ms")?;
        let cold_phases = summary["cold_phase_ms"]
            .as_object()
            .ok_or("missing cold phases")?;
        let costs: Vec<f64> = cold_phases
            .values()
            .map(|v| v.as_f64().ok_or("invalid cold phase".to_string()))
            .collect::<Check<_>>()?;
        checked_phase_sum(&costs, cold, "IC cold")?;
        ensure(
            point(&ic_target["published_q"])? == target
                && point(&rho_row["published_q"])? == target,
            "IC/rho point mismatch",
        )?;
        ensure(
            integer(&ic_target, "recovered_scalar")? == scalar
                && integer(&rho_row, "recovered_fixture_scalar")? == scalar,
            "IC/rho scalar mismatch",
        )?;
        ensure(
            ic_target["group_verified"] == true && rho_row["verified"] == true,
            "producer scalar verification absent",
        )?;
        let curve = BinaryCurve::new(
            integer(w, "n")? as u32,
            integer(w, "a")?,
            array(w, "field_modulus_low_terms")?,
        )?;
        ensure(
            curve.scale(scalar, point(&w["generator"])?)? == target,
            "independent scalar replay failed",
        )?;
        ensure(
            summary["rank"] == 244
                && summary["orbit_columns"] == 244
                && summary["factor_base_points"] == 25864,
            "IC base/rank differs from protocol",
        )?;
        ensure(
            rho_row["reference_grade"] == "strong"
                && rho_row["quotient_mode"] == "signed_frobenius",
            "rho policy differs from protocol",
        )?;
        Ok(
            json!({"status":"pending_independent_relation_replay","ic_online_ms":online,
            "rho_online_ms":floating(&rho_row,"walk_ms")?+floating(&rho_row,"validation_ms")?,
            "online_speedup":null,"independent_scalar_replay":true,
            "ic_cold_in_process_ms":cold,"ic_cold_phase_ms":cold_phases,
            "actual_factor_base_points":summary["factor_base_points"],"orbit_columns":summary["orbit_columns"],
            "rank_attempts":summary["rank_attempts"],"rank_relations":summary["rank_relations"],
            "rank_failures":summary["rank_failures"],"rho_walk_steps":rho_row["walk_steps"],
            "rho_automorphism_size":rho_row["automorphism_size"]}),
        )
    })();
    match checked {
        Ok(v) => v,
        Err(e) => {
            result["validation_error"] = json!(e);
            result
        }
    }
}

fn cpu_model() -> String {
    if cfg!(target_os = "macos") {
        for key in ["machdep.cpu.brand_string", "hw.model"] {
            if let Ok(out) = Command::new("/usr/sbin/sysctl").args(["-n", key]).output() {
                if out.status.success() {
                    return String::from_utf8_lossy(&out.stdout).trim().to_string();
                }
            }
        }
    } else if let Ok(cpuinfo) = fs::read_to_string("/proc/cpuinfo") {
        if let Some(line) = cpuinfo.lines().find(|s| s.starts_with("model name")) {
            return line
                .split_once(':')
                .map(|(_, s)| s.trim().to_string())
                .unwrap_or_default();
        }
    }
    "unknown".into()
}

fn execute(freeze_path: &Path, ic_binary: &Path, rho_binary: &Path, run_dir: &Path) -> Check<()> {
    let record = verify_freeze(freeze_path, ic_binary, rho_binary)?;
    let workload = verify_workload()?;
    fs::create_dir(run_dir).map_err(|e| format!("new run directory {}: {e}", run_dir.display()))?;
    let run_dir = run_dir.canonicalize().map_err(|e| e.to_string())?;
    let ambient: BTreeMap<_, _> = std::env::vars()
        .filter(|(k, _)| k.starts_with("KIC_"))
        .collect();
    let started = json!({"kind":"n53_compact_cold_attempt_started",
        "source_freeze_sha256":sha256_file(freeze_path)?,"source_commit":record["source_commit"],
        "workload_sha256":record["source_sha256"][format!("{SUBDIR}/workload.json")],
        "host":{"os":std::env::consts::OS,"arch":std::env::consts::ARCH,"cpu_model":cpu_model(),
            "logical_cpus":thread::available_parallelism().map(|n| n.get()).unwrap_or(0)},
        "environment_kic":ambient});
    write_new_json(&run_dir.join("attempt_started.json"), &started)?;
    let preflight = if !ambient.is_empty() {
        Err("unexpected ambient KIC_* variables".to_string())
    } else {
        rss_bytes(std::process::id()).map(|_| ())
    };
    if let Err(reason) = preflight {
        write_new_json(
            &run_dir.join("receipt.json"),
            &json!({"kind":"n53_compact_cold_preflight_failure",
            "reason":reason,"attempt_started_sha256":sha256_file(&run_dir.join("attempt_started.json"))?}),
        )?;
        println!(
            "{}",
            json!({"status":"preflight_failure","receipt":run_dir.join("receipt.json")})
        );
        std::process::exit(2);
    }
    let ic_args = vec![
        ic_binary
            .canonicalize()
            .map_err(|e| e.to_string())?
            .display()
            .to_string(),
        "construct:53:0:244".into(),
        protocol_path("target_points.jsonl").display().to_string(),
        "20261009".into(),
        run_dir.join("ic-target.jsonl").display().to_string(),
    ];
    let mut ic_environment = BTreeMap::new();
    ic_environment.insert(
        "KIC_DUMP_BASE".into(),
        run_dir.join("base.jsonl").display().to_string(),
    );
    ic_environment.insert(
        "KIC_DUMP_RANK".into(),
        run_dir.join("rank.jsonl").display().to_string(),
    );
    let ic = run_arm(
        &ic_args,
        &run_dir.join("ic.stdout.jsonl"),
        &run_dir.join("ic.stderr.txt"),
        &ic_environment,
    )?;
    let rho_args = vec![
        rho_binary
            .canonicalize()
            .map_err(|e| e.to_string())?
            .display()
            .to_string(),
        "53".into(),
        "0".into(),
        "signed_frobenius".into(),
        "1".into(),
        "strong".into(),
        "20261009".into(),
        integer(&workload, "verification_scalar")?.to_string(),
    ];
    let rho = run_arm(
        &rho_args,
        &run_dir.join("rho.stdout.jsonl"),
        &run_dir.join("rho.stderr.txt"),
        &BTreeMap::new(),
    )?;
    let mut files = BTreeMap::new();
    for name in [
        "attempt_started.json",
        "ic.stdout.jsonl",
        "ic.stderr.txt",
        "rho.stdout.jsonl",
        "rho.stderr.txt",
        "ic-target.jsonl",
        "base.jsonl",
        "rank.jsonl",
    ] {
        let file = run_dir.join(name);
        if file.exists() {
            files.insert(name, sha256_file(&file)?);
        }
    }
    let receipt = json!({"kind":"n53_compact_cold_attempt","ic":ic,"rho":rho,
        "analysis":analyse(&ic,&rho,&run_dir,&workload),"files_sha256":files,
        "source_freeze_sha256":sha256_file(freeze_path)?,"isolation_verified":false,
        "effective_environment_kic":{"ic":ic_environment,"rho":{}}});
    write_new_json(&run_dir.join("receipt.json"), &receipt)?;
    println!(
        "{}",
        json!({"status":receipt["analysis"]["status"],"receipt":run_dir.join("receipt.json")})
    );
    Ok(())
}

fn main() {
    if let Err(e) = entry() {
        eprintln!("native protocol failed: {e}");
        std::process::exit(1);
    }
}

fn entry() -> Check<()> {
    let mut args = std::env::args().skip(1);
    let mode = args
        .next()
        .ok_or("expected verify-workload, preflight, freeze, or run")?;
    let mut ic: Option<PathBuf> = None;
    let mut rho: Option<PathBuf> = None;
    let mut freeze_path: Option<PathBuf> = None;
    let mut run_dir: Option<PathBuf> = None;
    while let Some(arg) = args.next() {
        let value = args
            .next()
            .ok_or_else(|| format!("missing value for {arg}"))?;
        match arg.as_str() {
            "--ic-binary" => ic = Some(PathBuf::from(value)),
            "--rho-binary" => rho = Some(PathBuf::from(value)),
            "--freeze" => freeze_path = Some(PathBuf::from(value)),
            "--run-dir" => run_dir = Some(PathBuf::from(value)),
            _ => return Err(format!("unknown argument {arg}")),
        }
    }
    match mode.as_str() {
        "verify-workload" => {
            let w = verify_workload()?;
            println!(
                "{}",
                json!({"status":"PASS","target":w["primary_target"],"scalar":w["verification_scalar"]})
            );
        }
        "preflight" => {
            let rss = rss_bytes(std::process::id())?;
            println!("{}", json!({"status":"PASS","rss_bytes":rss}));
        }
        "freeze" => println!(
            "{}",
            serde_json::to_string_pretty(&freeze(
                &ic.ok_or("--ic-binary required")?,
                &rho.ok_or("--rho-binary required")?
            )?)
            .map_err(|e| e.to_string())?
        ),
        "run" => execute(
            &freeze_path.ok_or("--freeze required")?,
            &ic.ok_or("--ic-binary required")?,
            &rho.ok_or("--rho-binary required")?,
            &run_dir.ok_or("--run-dir required")?,
        )?,
        _ => return Err(format!("unknown mode {mode}")),
    }
    Ok(())
}
