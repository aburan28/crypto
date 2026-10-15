//! Native, sequential process receipts. Producer online clocks stay separate.
#![allow(dead_code)]
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::unix::fs::PermissionsExt;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicI32, Ordering};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

static INTERRUPTED: AtomicI32 = AtomicI32::new(0);
extern "C" fn interrupted(sig: i32) {
    INTERRUPTED.store(sig, Ordering::SeqCst);
}
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
fn err(s: impl Into<String>) -> Box<dyn std::error::Error> {
    s.into().into()
}
fn sha256(path: &Path) -> Result<String> {
    let out = Command::new("/usr/bin/sha256sum")
        .arg("--")
        .arg(path)
        .output()?;
    if !out.status.success() {
        return Err(err(format!("sha256sum failed: {}", path.display())));
    }
    let s = String::from_utf8(out.stdout)?;
    Ok(s.split_whitespace()
        .next()
        .ok_or_else(|| err("empty hash"))?
        .into())
}
fn epoch_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_millis()
}
fn write_new(path: &Path, value: &Value) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(path)?;
    serde_json::to_writer_pretty(&mut f, value)?;
    f.write_all(b"\n")?;
    f.sync_all()?;
    fs::set_permissions(path, fs::Permissions::from_mode(0o444))?;
    Ok(())
}
fn read_jsonl(path: &Path) -> Result<Vec<Value>> {
    let mut out = Vec::new();
    for (i, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        out.push(
            serde_json::from_str(&line)
                .map_err(|e| err(format!("{}:{}: {e}", path.display(), i + 1)))?,
        );
    }
    Ok(out)
}
fn allowed_cpus() -> Result<Vec<usize>> {
    let mut set: libc::cpu_set_t = unsafe { std::mem::zeroed() };
    if unsafe { libc::sched_getaffinity(0, std::mem::size_of_val(&set), &mut set) } != 0 {
        return Err(std::io::Error::last_os_error().into());
    }
    Ok((0..libc::CPU_SETSIZE as usize)
        .filter(|i| unsafe { libc::CPU_ISSET(*i, &set) })
        .collect())
}
fn point(v: &Value) -> Option<Vec<String>> {
    let a = v.as_array()?;
    if a.len() != 2 {
        return None;
    }
    a.iter()
        .map(|x| match x {
            Value::String(s) => Some(s.clone()),
            Value::Number(n) => Some(n.to_string()),
            _ => None,
        })
        .collect()
}
fn scalar(v: &Value) -> Option<String> {
    match v {
        Value::String(s) => Some(s.clone()),
        Value::Number(n) => Some(n.to_string()),
        _ => None,
    }
}
fn number(v: &Value, key: &str) -> Option<f64> {
    v.get(key)
        .and_then(Value::as_f64)
        .filter(|x| x.is_finite() && *x >= 0.)
}
pub fn rho_timing_valid(row: &Value) -> bool {
    let Some(phases) = row["online_phase_ns"].as_object() else {
        return false;
    };
    let mut sum = 0u64;
    for key in [
        "target_query_and_validation",
        "walk_and_collision",
        "generic_scalar_replay",
    ] {
        let Some(n) = phases.get(key).and_then(Value::as_u64) else {
            return false;
        };
        let Some(next) = sum.checked_add(n) else {
            return false;
        };
        sum = next;
    }
    let Some(total) = row["online_interval_ns"].as_u64() else {
        return false;
    };
    total > 0
        && total >= sum
        && row["online_phase_sum_ns"].as_u64() == Some(sum)
        && row["online_timer_boundary_residual_ns"].as_u64() == Some(total - sum)
        && row["target_ms"].as_f64() == Some(total as f64 / 1e6)
        && row["online_total_ms"].as_f64() == Some(total as f64 / 1e6)
}
fn preparation_ms(rows: &[Value]) -> Option<f64> {
    for row in rows {
        if let Some(t) = row.get("timing_ms") {
            let fs = ["setup_total", "rank_stage", "linear_algebra"];
            if fs
                .iter()
                .all(|k| t.get(k).and_then(Value::as_f64).is_some())
            {
                return Some(fs.iter().map(|k| t[k].as_f64().unwrap()).sum());
            }
        }
        if let Some(x) = number(row, "setup_ms") {
            return Some(x);
        }
        if let Some(x) = number(row, "preparation_ns_excluded") {
            return Some(x / 1e6);
        }
    }
    None
}
fn ic_summary_validate(path: &Path, expected_base_hash: Option<&str>) -> Value {
    ic_summary_validate_count(path, expected_base_hash, 1)
}
fn ic_summary_validate_count(
    path: &Path,
    expected_base_hash: Option<&str>,
    target_count: u64,
) -> Value {
    let rows = match read_jsonl(path) {
        Ok(r) => r,
        Err(e) => return json!({"valid":false,"reason":e.to_string()}),
    };
    let summaries: Vec<_> = rows
        .iter()
        .filter(|r| r["kind"] == "compact_orbit_dlp_summary")
        .collect();
    if summaries.len() != 1 {
        return json!({"valid":false,"reason":"expected exactly one compact_orbit_dlp_summary","summary_rows":summaries.len()});
    }
    let r = summaries[0];
    let rank_ok = r["rank"].as_u64() == r["orbit_columns"].as_u64()
        && matches!(r["orbit_columns"].as_u64(), Some(600 | 1271 | 1798));
    let target_ok = r["targets"].as_u64() == Some(target_count)
        && r["targets_solved"].as_u64() == Some(target_count)
        && r["targets_failed"].as_u64() == Some(0);
    let base_ok = expected_base_hash
        .map(|h| r["base_hash"].as_str() == Some(h))
        .unwrap_or(true);
    json!({"valid":rank_ok&&target_ok&&base_ok,"rank_ok":rank_ok,"target_counts_ok":target_ok,"base_hash_match":base_ok,"expected_base_hash":expected_base_hash,"selected_summary":r})
}
fn sem_validate(path: &Path, target: &Value, arm: &str) -> Value {
    let rows = match read_jsonl(path) {
        Ok(r) => r,
        Err(e) => return json!({"valid":false,"reason":e.to_string()}),
    };
    let target_point = point(target);
    let mut selected = Vec::new();
    for row in &rows {
        if row.get("recovered_scalar").is_some() || row.get("recovered_fixture_scalar").is_some() {
            selected.push(row);
        }
    }
    if selected.len() != 1 {
        return json!({"valid":false,"reason":"expected exactly one scalar-bearing target row", "target_rows":selected.len(),"jsonl_rows":rows.len()});
    }
    let r = selected[0];
    let observed = r
        .get("target")
        .or_else(|| r.get("published_q"))
        .and_then(point);
    let recovered = r
        .get("recovered_scalar")
        .or_else(|| r.get("recovered_fixture_scalar"))
        .and_then(scalar);
    let group = r
        .get("group_verified")
        .or_else(|| r.get("reference_group_validation"))
        .and_then(Value::as_bool)
        == Some(true);
    let verified = r
        .get("verified")
        .map(|x| x.as_bool() == Some(true))
        .unwrap_or(true);
    let no_explicit_failure = r
        .get("exit_code")
        .map(|x| x.as_i64() == Some(0))
        .unwrap_or(true);
    let scalar_present = recovered.as_ref().is_some_and(|s| !s.is_empty());
    let mut phases = json!({});
    let mut phase_sum: Option<f64> = None;
    let online = if arm == "ic" {
        let fields = [
            "target_query_ms",
            "target_pdp_ms",
            "target_relation_check_ms",
            "target_descent_ms",
            "target_recovery_check_ms",
        ];
        if fields.iter().all(|k| number(r, k).is_some()) {
            phase_sum = Some(fields.iter().map(|k| number(r, k).unwrap()).sum());
            for k in fields {
                phases[k] = r[k].clone();
            }
        }
        number(r, "target_ms").or_else(|| number(r, "online_total_ms"))
    } else {
        for k in [
            "walk_ms",
            "collision_ms",
            "recovery_check_ms",
            "validation_ms",
        ] {
            if let Some(x) = number(r, k) {
                phases[k] = json!(x);
            }
        }
        if let Some(ns) = r.get("online_phase_ns").and_then(Value::as_object) {
            for (k, v) in ns {
                if let Some(x) = v.as_f64() {
                    phases[k] = json!(x / 1e6);
                }
            }
            phase_sum = r
                .get("online_phase_sum_ns")
                .and_then(Value::as_f64)
                .map(|x| x / 1e6);
        }
        number(r, "online_total_ms").or_else(|| number(r, "target_ms"))
    };
    let phase_gap = online.zip(phase_sum).map(|(a, b)| a - b);
    let preparation = preparation_ms(&rows);
    json!({"valid":target_point.is_some()&&observed==target_point&&scalar_present&&group&&verified&&no_explicit_failure&&(!arm.starts_with("rho")||rho_timing_valid(r)),
        "public_target_match":observed==target_point,"scalar_present":scalar_present,"recovered_scalar":recovered,
        "group_verified":group,"verified_flag_ok":verified,"exit_code_field_ok":no_explicit_failure,
        "producer_online_ms":online,"producer_phase_sum_ms":phase_sum,"producer_phase_gap_ms":phase_gap,
        "producer_phases_ms":phases,"producer_preparation_ms":preparation,
        "online_boundary_qualified":online.is_some() && (arm!="ic" || phase_sum.is_some()),
        "scalar_is_verified_by_producer_group_replay_only":true,
        "independent_replay_required":true,"selected_row":r})
}
fn sem_validate_batch(path: &Path, targets: &[Value], arm: &str) -> Value {
    let rows = match read_jsonl(path) {
        Ok(r) => r,
        Err(e) => return json!({"valid":false,"reason":e.to_string()}),
    };
    let selected: Vec<_> = rows
        .iter()
        .filter(|r| {
            r.get("recovered_scalar").is_some() || r.get("recovered_fixture_scalar").is_some()
        })
        .collect();
    if selected.len() != targets.len() || targets.is_empty() {
        return json!({"valid":false,"reason":"expected one scalar-bearing row per declared public target","target_rows":selected.len(),"expected_targets":targets.len()});
    }
    let fields = [
        "target_query_ms",
        "target_pdp_ms",
        "target_relation_check_ms",
        "target_descent_ms",
        "target_recovery_check_ms",
    ];
    let mut results = Vec::new();
    let mut all_valid = true;
    for target in targets {
        let expected = point(target);
        let matching: Vec<_> = selected
            .iter()
            .filter(|r| {
                r.get("target")
                    .or_else(|| r.get("published_q"))
                    .and_then(point)
                    == expected
            })
            .collect();
        if matching.len() != 1 {
            all_valid = false;
            results.push(json!({"valid":false,"target":target,"matching_rows":matching.len()}));
            continue;
        }
        let r = *matching[0];
        let scalar_value = r.get("recovered_scalar").and_then(scalar);
        let ic_phases_ok = fields.iter().all(|k| number(r, k).is_some());
        let rho_phase_sum = r["online_phase_sum_ns"]
            .as_f64()
            .filter(|n| n.is_finite() && *n >= 0.)
            .map(|n| n / 1e6);
        let phases_ok = if arm == "ic_batch" {
            ic_phases_ok
        } else {
            rho_timing_valid(r)
        };
        let phase_sum = if ic_phases_ok {
            Some(fields.iter().map(|k| number(r, k).unwrap()).sum::<f64>())
        } else {
            rho_phase_sum
        };
        let online = number(r, "target_ms");
        let valid = scalar_value.is_some()
            && r["group_verified"] == true
            && phases_ok
            && online.is_some()
            && r["exit_code"].as_i64() == Some(0);
        all_valid &= valid;
        results.push(json!({"valid":valid,"target":target,"recovered_scalar":scalar_value,"group_verified":r["group_verified"],"target_ms":online,"producer_phase_sum_ms":phase_sum,
            "producer_phase_gap_ms":online.zip(phase_sum).map(|(a,b)|a-b),"selected_row":r}));
    }
    json!({"valid":all_valid,"target_count":targets.len(),"per_target_rows":results,"independent_replay_required":true,"shared_preparation_is_not_cold_single_target":true})
}
pub fn launch(args: &[String]) -> Result<i32> {
    run_arm(args)
}
pub fn was_interrupted() -> bool {
    INTERRUPTED.load(Ordering::SeqCst) != 0
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scratch(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root =
            std::env::temp_dir().join(format!("n83-matched-{name}-{}-{nonce}", std::process::id()));
        fs::create_dir(&root).unwrap();
        root
    }
    fn jsonl(root: &Path, name: &str, v: &Value) -> PathBuf {
        let p = root.join(name);
        fs::write(&p, format!("{v}\n")).unwrap();
        p
    }
    #[test]
    fn missing_rank_or_wrong_base_rejects_ic_summary() {
        let root = scratch("summary");
        let good = json!({"kind":"compact_orbit_dlp_summary","rank":600,"orbit_columns":600,"targets":32,"targets_solved":32,"targets_failed":0,"base_hash":"base"});
        let p = jsonl(&root, "good.jsonl", &good);
        assert_eq!(
            ic_summary_validate_count(&p, Some("base"), 32)["valid"],
            true
        );
        assert_eq!(
            ic_summary_validate_count(&p, Some("wrong"), 32)["valid"],
            false
        );
        let mut bad = good.clone();
        bad["rank"] = Value::Null;
        let p = jsonl(&root, "missing.jsonl", &bad);
        assert_eq!(
            ic_summary_validate_count(&p, Some("base"), 32)["valid"],
            false
        );
    }
    #[test]
    fn ic_batch_cannot_substitute_rho_phases_or_duplicate_points() {
        let root = scratch("batch");
        let p = jsonl(
            &root,
            "targets.jsonl",
            &json!({"target":["1","2"],"recovered_scalar":"3","group_verified":true,"exit_code":0,"target_ms":1.,"online_phase_sum_ns":1000000}),
        );
        assert_eq!(
            sem_validate_batch(&p, &[json!(["1", "2"])], "ic_batch")["valid"],
            false
        );
        assert_eq!(
            sem_validate_batch(&p, &[json!(["1", "2"]), json!(["1", "2"])], "rho_batch")["valid"],
            false
        );
    }
    fn launch_args(root: &Path, timeout: &str, memory: &str) -> Vec<String> {
        let target = jsonl(root, "point.jsonl", &json!(["1", "2"]));
        let cpu = allowed_cpus()
            .unwrap()
            .into_iter()
            .find(|c| *c == 3)
            .unwrap_or_else(|| allowed_cpus().unwrap()[0]);
        vec![
            "--arm".into(),
            "rho".into(),
            "--out-dir".into(),
            root.join("run").to_string_lossy().into(),
            "--target-file".into(),
            target.to_string_lossy().into(),
            "--cwd".into(),
            root.to_string_lossy().into(),
            "--timeout-s".into(),
            timeout.into(),
            "--memory-bytes".into(),
            memory.into(),
            "--cpu".into(),
            cpu.to_string(),
            "--".into(),
            "/usr/bin/sleep".into(),
            "10".into(),
        ]
    }
    #[test]
    fn zero_budget_or_memory_is_rejected_before_launch() {
        let root = scratch("budget");
        assert!(launch(&launch_args(&root, "0", "12884901888")).is_err());
        let root = scratch("memory");
        assert!(launch(&launch_args(&root, "1", "0")).is_err());
    }
    #[test]
    fn timeout_is_recorded_and_process_is_terminated() {
        let root = scratch("timeout");
        assert_eq!(
            launch(&launch_args(&root, "0.02", "12884901888")).unwrap(),
            1
        );
        let receipt: Value =
            serde_json::from_reader(File::open(root.join("run/execution.json")).unwrap()).unwrap();
        assert_eq!(receipt["status"], "TIMEOUT");
        assert_eq!(receipt["timeout"], true);
        assert!(receipt["elapsed_process_wall_s"].as_f64().unwrap() < 3.);
        assert_eq!(
            receipt["resource_envelope"]["timeout_s"].as_f64(),
            Some(0.02)
        );
        assert!(receipt["peak_rss_bytes"].is_number());
    }
    #[test]
    fn address_space_envelope_is_retained_on_allocation_failure() {
        let root = scratch("as-limit");
        assert_eq!(launch(&launch_args(&root, "1", "1048576")).unwrap(), 1);
        let receipt: Value =
            serde_json::from_reader(File::open(root.join("run/execution.json")).unwrap()).unwrap();
        assert_eq!(
            receipt["resource_envelope"]["address_space_limit_bytes"],
            1048576
        );
        assert_eq!(receipt["memory_limit_hit"], Value::Null);
        assert_eq!(receipt["status"], "PROCESS_FAILED");
    }
}

fn run_arm(args: &[String]) -> Result<i32> {
    let mut opts: BTreeMap<String, String> = BTreeMap::new();
    let mut supplied_env = BTreeMap::new();
    let mut i = 0;
    while i < args.len() && args[i] != "--" {
        let key = args[i].clone();
        i += 1;
        let val = args
            .get(i)
            .ok_or_else(|| err(format!("missing value for {key}")))?
            .clone();
        i += 1;
        if key == "--env" {
            let (k, v) = val
                .split_once('=')
                .ok_or_else(|| err("env needs KEY=VALUE"))?;
            supplied_env.insert(k.to_owned(), v.to_owned());
        } else if opts.insert(key.clone(), val).is_some() {
            return Err(err(format!("duplicate {key}")));
        }
    }
    if args.get(i).map(String::as_str) != Some("--") || i + 1 >= args.len() {
        return Err(err("run-arm needs -- COMMAND [ARGS]"));
    }
    let argv = &args[i + 1..];
    let req = |k: &str| {
        opts.get(k)
            .cloned()
            .ok_or_else(|| err(format!("missing {k}")))
    };
    let arm = req("--arm")?;
    let out = PathBuf::from(req("--out-dir")?);
    let target_file = PathBuf::from(req("--target-file")?);
    let cwd = PathBuf::from(req("--cwd")?);
    let timeout: f64 = req("--timeout-s")?.parse()?;
    if !timeout.is_finite() || timeout <= 0. {
        return Err(err("timeout must be finite and positive"));
    }
    let memory: u64 = opts
        .get("--memory-bytes")
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(12 * 1024 * 1024 * 1024);
    if memory == 0 {
        return Err(err("finite memory cap required"));
    }
    let cpus = allowed_cpus()?;
    let cpu: usize = opts
        .get("--cpu")
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(*cpus.first().ok_or_else(|| err("no allowed CPUs"))?);
    if !cpus.contains(&cpu) {
        return Err(err("selected CPU outside inherited allowed cpuset"));
    }
    let target_rows = read_jsonl(&target_file)?;
    let target_count: usize = opts
        .get("--target-count")
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(1);
    if ![1, 8].contains(&target_count)
        || target_rows.len() != target_count
        || target_rows.iter().any(|r| point(r).is_none())
    {
        return Err(err(
            "target-file needs the declared singleton or 8 two-coordinate public points",
        ));
    }
    // Public points only; no scalar-bearing environment reaches measured child.
    let allowed = [
        "KIC_RHO_PUBLIC_TARGET_POINT",
        "KIC_RHO_WALK_SEED",
        "KIC_DUMP_BASE",
        "KIC_DUMP_RANK",
        "KIC_PACKED_ROOT_TABLE",
        "KIC_PACKED_ROOT_TABLE_DENSE",
        "KIC_RANK_THREADS",
        "KIC_RHO_THREADS",
        "KIC_TARGET_THREADS",
        "KIC_THREADS",
        "RAYON_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
    ];
    for k in supplied_env.keys() {
        if !allowed.contains(&k.as_str()) {
            return Err(err(format!("environment key not whitelisted: {k}")));
        }
    }
    for k in [
        "KIC_RANK_THREADS",
        "KIC_RHO_THREADS",
        "KIC_TARGET_THREADS",
        "KIC_THREADS",
        "RAYON_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
    ] {
        if let Some(v) = supplied_env.get(k) {
            if v != "1" {
                return Err(err(format!("single-thread envelope requires {k}=1")));
            }
        }
        supplied_env.insert(k.into(), "1".into());
    }
    supplied_env.insert("PATH".into(), "/usr/bin:/bin".into());
    supplied_env.insert("LC_ALL".into(), "C".into());
    if supplied_env.contains_key("KIC_RHO_PUBLIC_TARGET_POINT")
        && (target_count != 1
            || serde_json::from_str::<Value>(&supplied_env["KIC_RHO_PUBLIC_TARGET_POINT"])?
                != target_rows[0])
    {
        return Err(err("env target differs from frozen public point"));
    }
    fs::create_dir(&out)?;
    let stdout_path = out.join("stdout.jsonl");
    let stderr_path = out.join("stderr.txt");
    let stdout = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&stdout_path)?;
    let stderr = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&stderr_path)?;
    let binary_hash = sha256(Path::new(&argv[0]))?;
    let target_hash = sha256(&target_file)?;
    let mut inputs = BTreeMap::new();
    for arg in &argv[1..] {
        let p = Path::new(arg);
        if p.is_file() {
            inputs.insert(arg.clone(), sha256(p)?);
        }
    }
    let started_epoch_ms = epoch_ms();
    let envelope = json!({"cpu_count":1,"cpu_affinity":[cpu],"threads":1,"address_space_limit_bytes":memory,"address_space_limit_mechanism":"RLIMIT_AS soft=hard","timeout_s":timeout,"rank_threads":1,"online_target_threads":1});
    write_new(
        &out.join("started.json"),
        &json!({"schema":"n83.process-start/v1","arm":arm,"argv":argv,"cwd":cwd,"environment":supplied_env,"resource_envelope":envelope,"target_file":target_file,"target_sha256":target_hash,"binary_sha256":binary_hash,"argv_existing_file_sha256":inputs,"started_epoch_ms":started_epoch_ms}),
    )?;
    let mut command = Command::new(&argv[0]);
    command
        .args(&argv[1..])
        .current_dir(&cwd)
        .env_clear()
        .envs(&supplied_env)
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(stderr);
    unsafe {
        command.pre_exec(move || {
            if libc::setpgid(0, 0) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            let lim = libc::rlimit {
                rlim_cur: memory as libc::rlim_t,
                rlim_max: memory as libc::rlim_t,
            };
            if libc::setrlimit(libc::RLIMIT_AS, &lim) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            let mut set: libc::cpu_set_t = std::mem::zeroed();
            libc::CPU_ZERO(&mut set);
            libc::CPU_SET(cpu, &mut set);
            if libc::sched_setaffinity(0, std::mem::size_of_val(&set), &set) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            Ok(())
        });
        libc::signal(libc::SIGINT, interrupted as *const () as libc::sighandler_t);
        libc::signal(
            libc::SIGTERM,
            interrupted as *const () as libc::sighandler_t,
        );
    }
    let clock = Instant::now();
    let child = match command.spawn() {
        Ok(c) => c,
        Err(e) => {
            let semantic_path = opts
                .get("--semantic-file")
                .map(PathBuf::from)
                .unwrap_or_else(|| stdout_path.clone());
            write_new(
                &out.join("execution.json"),
                &json!({"schema":"n83.process-receipt/v1","arm":arm,"status":"LAUNCH_FAILED","error":e.to_string(),"argv":argv,"cwd":cwd,"environment":supplied_env,"resource_envelope":envelope,"started_epoch_ms":started_epoch_ms,"finished_epoch_ms":epoch_ms(),"elapsed_process_wall_s":clock.elapsed().as_secs_f64(),"cpu_user_s":null,"cpu_system_s":null,"cpu_total_s":null,"peak_rss_bytes":null,"peak_rss_native_kib":null,"return_code":null,"termination_signal":null,"timeout":false,"memory_limit_hit":null,"memory_limit_note":"No child was successfully launched; resource usage is unavailable","binary_sha256":binary_hash,"target_sha256":target_hash,"argv_existing_file_sha256":inputs,"stdout_path":stdout_path,"stderr_path":stderr_path,"semantic_path":semantic_path,"stdout_sha256":sha256(&stdout_path)?,"stderr_sha256":sha256(&stderr_path)?,"semantic_sha256":if semantic_path.is_file(){Some(sha256(&semantic_path)?)}else{None},"semantic":{"valid":false,"reason":"child launch failed"},"known_answer_environment_present":false}),
            )?;
            for p in [&stdout_path, &stderr_path] {
                fs::set_permissions(p, fs::Permissions::from_mode(0o444))?;
            }
            return Ok(1);
        }
    };
    let pid = child.id() as libc::pid_t;
    let mut status = 0;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    let mut termination: Option<(Instant, &str)> = None;
    let mut last_beat = Instant::now();
    let mut heartbeats = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out.join("heartbeat.jsonl"))?;
    loop {
        let w = unsafe { libc::wait4(pid, &mut status, libc::WNOHANG, &mut usage) };
        if w == pid {
            break;
        }
        if w < 0 {
            return Err(std::io::Error::last_os_error().into());
        }
        if termination.is_none() {
            let sig = INTERRUPTED.load(Ordering::SeqCst);
            if sig != 0 {
                unsafe { libc::kill(-pid, libc::SIGTERM) };
                termination = Some((Instant::now(), "INTERRUPTED"));
            } else if clock.elapsed().as_secs_f64() >= timeout {
                unsafe { libc::kill(-pid, libc::SIGTERM) };
                termination = Some((Instant::now(), "TIMEOUT"));
            }
        } else if termination.as_ref().unwrap().0.elapsed() >= Duration::from_secs(2) {
            unsafe { libc::kill(-pid, libc::SIGKILL) };
        }
        if last_beat.elapsed() >= Duration::from_secs(30) {
            let beat = json!({"pid":pid,"epoch_ms":epoch_ms(),"elapsed_process_wall_s":clock.elapsed().as_secs_f64(),"termination":termination.as_ref().map(|x|x.1)});
            writeln!(heartbeats, "{}", beat)?;
            heartbeats.flush()?;
            eprintln!("{}", beat);
            last_beat = Instant::now();
        }
        std::thread::sleep(Duration::from_millis(5));
    }
    let elapsed = clock.elapsed().as_secs_f64();
    let exit_code = if libc::WIFEXITED(status) {
        Some(libc::WEXITSTATUS(status))
    } else {
        None
    };
    let signal = if libc::WIFSIGNALED(status) {
        Some(libc::WTERMSIG(status))
    } else {
        None
    };
    let semantic_path = opts
        .get("--semantic-file")
        .map(PathBuf::from)
        .unwrap_or_else(|| stdout_path.clone());
    let mut sem = if target_count == 1 {
        sem_validate(&semantic_path, &target_rows[0], &arm)
    } else {
        sem_validate_batch(&semantic_path, &target_rows, &arm)
    };
    if semantic_path != stdout_path {
        if let Ok(rows) = read_jsonl(&stdout_path) {
            if let Some(x) = preparation_ms(&rows) {
                sem["producer_preparation_ms"] = json!(x);
            }
        }
    }
    if arm == "ic" || arm == "ic_batch" {
        let summary = ic_summary_validate_count(
            &stdout_path,
            opts.get("--expected-base-hash").map(String::as_str),
            target_count as u64,
        );
        let phases = target_count == 8 || sem["producer_phase_sum_ms"].is_number();
        sem["valid"] = json!(sem["valid"] == true && summary["valid"] == true && phases);
        sem["summary_validation"] = summary;
        sem["five_phase_exact_conformance"] =
            json!(sem["producer_phase_gap_ms"].as_f64() == Some(0.));
        sem["five_phase_residual_reported"] = json!(sem["producer_phase_gap_ms"].is_number());
    }
    let success = exit_code == Some(0) && sem["valid"] == true && termination.is_none();
    let disposition = if let Some((_, s)) = termination.as_ref() {
        *s
    } else if success {
        "COMPLETED_VERIFIED"
    } else if exit_code != Some(0) {
        "PROCESS_FAILED"
    } else {
        "SEMANTIC_VALIDATION_FAILED"
    };
    let cpu_user = usage.ru_utime.tv_sec as f64 + usage.ru_utime.tv_usec as f64 / 1e6;
    let cpu_system = usage.ru_stime.tv_sec as f64 + usage.ru_stime.tv_usec as f64 / 1e6;
    let receipt = json!({"schema":"n83.process-receipt/v1","arm":arm,"status":disposition,"pid":pid,
        "argv":argv,"cwd":cwd,"environment":supplied_env,"resource_envelope":envelope,
        "started_epoch_ms":started_epoch_ms,"finished_epoch_ms":epoch_ms(),"elapsed_process_wall_s":elapsed,
        "cpu_user_s":cpu_user,"cpu_system_s":cpu_system,"cpu_total_s":cpu_user+cpu_system,
        "peak_rss_bytes":usage.ru_maxrss as u64*1024,"peak_rss_native_kib":usage.ru_maxrss,
        "peak_rss_source":"Linux wait4(2) ru_maxrss KiB multiplied by 1024; exact child process",
        "exit_observation_poll_ms":5,"exit_observation_delay_note":"Nominal 0-5 ms plus scheduling and poll-body delays; no strict bound under arbitrary scheduling",
        "return_code":exit_code,"termination_signal":signal,"raw_wait_status":status,
        "timeout":disposition=="TIMEOUT","interrupted_signal":INTERRUPTED.load(Ordering::SeqCst),
        "memory_limit_hit":null,"memory_limit_note":"RLIMIT_AS allocation refusal may appear as nonzero/abort; no inference from exit code alone",
        "stdout_path":stdout_path,"stderr_path":stderr_path,"semantic_path":semantic_path,
        "stdout_sha256":sha256(&stdout_path)?,"stderr_sha256":sha256(&stderr_path)?,
        "semantic_sha256":if semantic_path.is_file(){Some(sha256(&semantic_path)?)}else{None},
        "binary_sha256":binary_hash,"target_sha256":target_hash,"argv_existing_file_sha256":inputs,"semantic":sem,
        "accounting":"Whole-process wall/CPU include launch, load, preparation, online, replay. Producer online interval is separate.",
        "known_answer_environment_present":false});
    write_new(&out.join("execution.json"), &receipt)?;
    heartbeats.sync_all()?;
    drop(heartbeats);
    for p in [&stdout_path, &stderr_path, &out.join("heartbeat.jsonl")] {
        fs::set_permissions(p, fs::Permissions::from_mode(0o444))?;
    }
    if semantic_path != stdout_path && semantic_path.is_file() {
        fs::set_permissions(&semantic_path, fs::Permissions::from_mode(0o444))?;
    }
    println!(
        "{}",
        json!({"receipt":out.join("execution.json"),"status":disposition,"elapsed_process_wall_s":elapsed,"peak_rss_bytes":usage.ru_maxrss as u64*1024,"semantic_valid":success})
    );
    Ok(if success { 0 } else { 1 })
}
