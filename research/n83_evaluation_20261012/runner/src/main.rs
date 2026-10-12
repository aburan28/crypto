//! Native, sequential process receipts. Producer online clocks stay separate.
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
    let rank_ok = r["rank"].as_u64() == Some(600) && r["orbit_columns"].as_u64() == Some(600);
    let target_ok = r["targets"].as_u64() == Some(1)
        && r["targets_solved"].as_u64() == Some(1)
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
    json!({"valid":target_point.is_some()&&observed==target_point&&scalar_present&&group&&verified&&no_explicit_failure,
        "public_target_match":observed==target_point,"scalar_present":scalar_present,"recovered_scalar":recovered,
        "group_verified":group,"verified_flag_ok":verified,"exit_code_field_ok":no_explicit_failure,
        "producer_online_ms":online,"producer_phase_sum_ms":phase_sum,"producer_phase_gap_ms":phase_gap,
        "producer_phases_ms":phases,"producer_preparation_ms":preparation,
        "online_boundary_qualified":online.is_some() && (arm!="ic" || phase_sum.is_some()),
        "scalar_is_verified_by_producer_group_replay_only":true,
        "independent_replay_required":true,"selected_row":r})
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
    if target_rows.len() != 1 || point(&target_rows[0]).is_none() {
        return Err(err(
            "target-file needs exactly one two-coordinate public point",
        ));
    }
    // Public points only; no scalar-bearing environment reaches measured child.
    let allowed = [
        "KIC_RHO_PUBLIC_TARGET_POINT",
        "KIC_RHO_WALK_SEED",
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
        && serde_json::from_str::<Value>(&supplied_env["KIC_RHO_PUBLIC_TARGET_POINT"])?
            != target_rows[0]
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
            write_new(
                &out.join("execution.json"),
                &json!({"schema":"n83.process-receipt/v1","status":"LAUNCH_FAILED","error":e.to_string(),"argv":argv,"resource_envelope":envelope,"started_epoch_ms":started_epoch_ms,"elapsed_process_wall_s":clock.elapsed().as_secs_f64()}),
            )?;
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
        std::thread::sleep(Duration::from_millis(100));
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
    let mut sem = sem_validate(&semantic_path, &target_rows[0], &arm);
    if semantic_path != stdout_path {
        if let Ok(rows) = read_jsonl(&stdout_path) {
            if let Some(x) = preparation_ms(&rows) {
                sem["producer_preparation_ms"] = json!(x);
            }
        }
    }
    if arm == "ic" {
        let summary = ic_summary_validate(
            &stdout_path,
            opts.get("--expected-base-hash").map(String::as_str),
        );
        sem["valid"] = json!(
            sem["valid"] == true
                && summary["valid"] == true
                && sem["producer_phase_sum_ms"].is_number()
        );
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

#[derive(Clone)]
struct Sample {
    pair_id: String,
    kind: String,
    log_ratio: f64,
    cold_log_ratio: f64,
    baseline_ms: f64,
    candidate_ms: f64,
}
struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}
fn percentile(v: &[f64], p: f64) -> f64 {
    let pos = p * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    v[lo] * (1. - (pos - lo as f64)) + v[hi] * (pos - lo as f64)
}
fn stats(samples: &[Sample], seed: u64, draws: usize) -> Value {
    if samples.is_empty() {
        return json!({"n":0,"status":"no_qualified_pairs"});
    }
    let n = samples.len();
    let gm = (samples.iter().map(|s| s.log_ratio).sum::<f64>() / n as f64).exp();
    let mut ratios: Vec<f64> = samples.iter().map(|s| s.log_ratio.exp()).collect();
    ratios.sort_by(f64::total_cmp);
    let cold_gm = (samples.iter().map(|s| s.cold_log_ratio).sum::<f64>() / n as f64).exp();
    let mut rng = SplitMix(seed);
    let mut boot = Vec::with_capacity(draws);
    let mut cold_boot = Vec::with_capacity(draws);
    for _ in 0..draws {
        let mut sum = 0.;
        let mut cold = 0.;
        for _ in 0..n {
            let s = &samples[(rng.next() % n as u64) as usize];
            sum += s.log_ratio;
            cold += s.cold_log_ratio;
        }
        boot.push((sum / n as f64).exp());
        cold_boot.push((cold / n as f64).exp());
    }
    boot.sort_by(f64::total_cmp);
    cold_boot.sort_by(f64::total_cmp);
    json!({"n":n,"geometric_mean_paired_online_ratio":gm,"median_paired_online_ratio":percentile(&ratios,0.5),
        "paired_online_ratio_range":[ratios[0],ratios[n-1]],"paired_percentile_bootstrap_95_ci":if n>1{json!([percentile(&boot,0.025),percentile(&boot,0.975)])}else{Value::Null},
        "geometric_mean_paired_whole_process_ratio":cold_gm,"whole_process_bootstrap_95_ci":if n>1{json!([percentile(&cold_boot,0.025),percentile(&cold_boot,0.975)])}else{Value::Null},
        "bootstrap_resamples":draws,"bootstrap_seed":seed,"resampling_unit":"matched complete pair; arithmetic mean log-ratio then exponentiation",
        "confidence_scope":"seed/run variation on one frozen target and one host; no independent-target or host generalization",
        "insufficient_replication":n<5,
        "rows":samples.iter().map(|s|json!({"pair_id":s.pair_id,"kind":s.kind,"baseline_online_ms":s.baseline_ms,"candidate_online_ms":s.candidate_ms,"online_ratio":s.log_ratio.exp(),"whole_process_ratio":s.cold_log_ratio.exp()})).collect::<Vec<_>>()})
}
fn analyze(args: &[String]) -> Result<i32> {
    let mut opts = BTreeMap::new();
    if args.len() % 2 != 0 {
        return Err(err("analyze options need values"));
    }
    for p in args.chunks(2) {
        opts.insert(p[0].as_str(), p[1].as_str());
    }
    let pairs = PathBuf::from(*opts.get("--pairs").ok_or_else(|| err("missing --pairs"))?);
    let out = PathBuf::from(*opts.get("--out").ok_or_else(|| err("missing --out"))?);
    let draws: usize = opts.get("--draws").copied().unwrap_or("20000").parse()?;
    if draws < 1000 {
        return Err(err("bootstrap draws must be >=1000"));
    }
    let seed: u64 = opts
        .get("--seed")
        .copied()
        .unwrap_or("2026101283")
        .parse()?;
    let mut samples = Vec::new();
    let mut excluded = Vec::new();
    let rows = read_jsonl(&pairs)?;
    for p in &rows {
        let id = p["pair_id"]
            .as_str()
            .ok_or_else(|| err("pair_id missing"))?
            .to_string();
        let kind = p["kind"].as_str().unwrap_or("matched").to_string();
        let mut receipts = Vec::new();
        for k in ["baseline_receipt", "candidate_receipt"] {
            let path = p[k].as_str().ok_or_else(|| err(format!("{k} missing")))?;
            let r: Value = serde_json::from_reader(File::open(path)?)?;
            receipts.push(r);
        }
        let a = &receipts[0];
        let b = &receipts[1];
        let same_target =
            a["target_sha256"] == b["target_sha256"] && a["target_sha256"].is_string();
        let same_envelope = a["resource_envelope"] == b["resource_envelope"];
        let same_scalar = a["semantic"]["recovered_scalar"] == b["semantic"]["recovered_scalar"]
            && a["semantic"]["recovered_scalar"].is_string();
        let success = a["status"] == "COMPLETED_VERIFIED"
            && b["status"] == "COMPLETED_VERIFIED"
            && a["semantic"]["valid"] == true
            && b["semantic"]["valid"] == true
            && a["semantic"]["online_boundary_qualified"] == true
            && b["semantic"]["online_boundary_qualified"] == true
            && same_scalar;
        let a_ms = a["semantic"]["producer_online_ms"].as_f64();
        let b_ms = b["semantic"]["producer_online_ms"].as_f64();
        if !same_target
            || !same_envelope
            || !success
            || a_ms.is_none()
            || b_ms.is_none()
            || a_ms.unwrap_or(0.) <= 0.
            || b_ms.unwrap_or(0.) <= 0.
        {
            excluded.push(json!({"pair_id":id,"kind":kind,"same_target":same_target,"same_resource_envelope":same_envelope,"same_recovered_scalar":same_scalar,"both_successful":success,"baseline_status":a["status"],"candidate_status":b["status"],"baseline_online_ms":a_ms,"candidate_online_ms":b_ms}));
            continue;
        }
        let cold_a = a["elapsed_process_wall_s"]
            .as_f64()
            .ok_or_else(|| err("wall missing"))?;
        let cold_b = b["elapsed_process_wall_s"]
            .as_f64()
            .ok_or_else(|| err("wall missing"))?;
        samples.push(Sample {
            pair_id: id,
            kind,
            log_ratio: (a_ms.unwrap() / b_ms.unwrap()).ln(),
            cold_log_ratio: (cold_a / cold_b).ln(),
            baseline_ms: a_ms.unwrap(),
            candidate_ms: b_ms.unwrap(),
        });
    }
    let aa: Vec<_> = samples.iter().filter(|s| s.kind == "aa").cloned().collect();
    let matched: Vec<_> = samples
        .iter()
        .filter(|s| s.kind == "matched")
        .cloned()
        .collect();
    let summary = json!({"schema":"n83.paired-statistics/v1","pairs_sha256":sha256(&pairs)?,"input_pair_count":rows.len(),"qualified_pair_count":samples.len(),"excluded_pairs":excluded,"aa":stats(&aa,seed,draws),"matched":stats(&matched,seed+1,draws),"primary_measurement_gate":"Requires matching resource envelopes, verified scalar/point semantic rows and complete producer online boundaries. Independent replay and a second host still required for promotion.","cold_cost":"Whole-process wall includes preparation plus online and launch/load; measured separately from producer online."});
    write_new(&out, &summary)?;
    println!("{}", summary);
    Ok(0)
}
fn panel(args: &[String]) -> Result<i32> {
    let mut opts = BTreeMap::new();
    if args.len() % 2 != 0 {
        return Err(err("panel options need values"));
    }
    for p in args.chunks(2) {
        opts.insert(p[0].as_str(), p[1].as_str());
    }
    let protocol_path = PathBuf::from(
        *opts
            .get("--protocol")
            .ok_or_else(|| err("missing --protocol"))?,
    );
    let root = PathBuf::from(
        *opts
            .get("--out-dir")
            .ok_or_else(|| err("missing --out-dir"))?,
    );
    let protocol: Value = serde_json::from_reader(File::open(&protocol_path)?)?;
    let runs = protocol["runs"]
        .as_array()
        .ok_or_else(|| err("protocol.runs array missing"))?;
    let resume = opts.get("--resume").copied() == Some("true");
    if !resume {
        fs::create_dir(&root)?;
        write_new(&root.join("protocol.json"), &protocol)?;
        write_new(
            &root.join("panel_start.json"),
            &json!({"schema":"n83.panel-start/v1","protocol_path":protocol_path,"protocol_sha256":sha256(&protocol_path)?,"started_epoch_ms":epoch_ms(),"planned_arm_count":runs.len()}),
        )?;
        let mut pairs = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(root.join("pairs.jsonl"))?;
        for p in protocol["pairs"]
            .as_array()
            .ok_or_else(|| err("protocol.pairs missing"))?
        {
            let a = p["baseline_id"]
                .as_str()
                .ok_or_else(|| err("baseline_id missing"))?;
            let b = p["candidate_id"]
                .as_str()
                .ok_or_else(|| err("candidate_id missing"))?;
            writeln!(
                pairs,
                "{}",
                json!({"pair_id":p["pair_id"],"kind":p["kind"],"baseline_receipt":root.join(a).join("execution.json"),"candidate_receipt":root.join(b).join("execution.json")})
            )?;
        }
        pairs.sync_all()?;
        drop(pairs);
        fs::set_permissions(root.join("pairs.jsonl"), fs::Permissions::from_mode(0o444))?;
    } else {
        let frozen: Value = serde_json::from_reader(File::open(root.join("protocol.json"))?)?;
        if frozen != protocol {
            return Err(err("resume protocol differs from frozen panel protocol"));
        }
    }
    let mut dispositions = Vec::new();
    let mut failures = 0;
    for run in runs {
        let id = run["id"].as_str().ok_or_else(|| err("run id missing"))?;
        if id.is_empty() || id.contains('/') || id.contains("..") {
            return Err(err("run id must be a simple directory name"));
        }
        let out = root.join(id);
        if out.join("execution.json").exists() {
            let r: Value = serde_json::from_reader(File::open(out.join("execution.json"))?)?;
            let good = r["status"] == "COMPLETED_VERIFIED";
            if !good {
                failures += 1;
            }
            dispositions
                .push(json!({"id":id,"status":r["status"],"resumed_existing_receipt":true}));
            continue;
        }
        if out.exists() {
            return Err(err(format!(
                "incomplete preexisting arm {}; preserve it and use a new protocol/panel for rerun",
                out.display()
            )));
        }
        let subst = |s: &str| s.replace("{OUT}", &out.to_string_lossy());
        let mut arm_args = Vec::new();
        for (flag, key) in [
            ("--arm", "arm"),
            ("--target-file", "target_file"),
            ("--cwd", "cwd"),
            ("--timeout-s", "timeout_s"),
            ("--memory-bytes", "memory_bytes"),
            ("--cpu", "cpu"),
        ] {
            let v = run
                .get(key)
                .or_else(|| protocol.get(key))
                .ok_or_else(|| err(format!("panel missing {key}")))?;
            let s = v
                .as_str()
                .map(str::to_owned)
                .unwrap_or_else(|| v.to_string());
            arm_args.push(flag.into());
            arm_args.push(subst(&s));
        }
        arm_args.push("--out-dir".into());
        arm_args.push(out.to_string_lossy().to_string());
        if let Some(s) = run["semantic_file"].as_str() {
            arm_args.push("--semantic-file".into());
            arm_args.push(subst(s));
        }
        if let Some(s) = run
            .get("expected_base_hash")
            .or_else(|| protocol.get("expected_base_hash"))
            .and_then(Value::as_str)
        {
            arm_args.push("--expected-base-hash".into());
            arm_args.push(s.into());
        }
        if let Some(env) = run["env"].as_object() {
            for (k, v) in env {
                let s = v
                    .as_str()
                    .ok_or_else(|| err("env values must be strings"))?;
                arm_args.push("--env".into());
                arm_args.push(format!("{k}={s}"));
            }
        }
        arm_args.push("--".into());
        for a in run["argv"]
            .as_array()
            .ok_or_else(|| err("run argv missing"))?
        {
            arm_args.push(subst(
                a.as_str()
                    .ok_or_else(|| err("argv entries must be strings"))?,
            ));
        }
        eprintln!(
            "{}",
            json!({"panel_starting_arm":id,"epoch_ms":epoch_ms(),"arm":run["arm"]})
        );
        let code = run_arm(&arm_args)?;
        if code != 0 {
            failures += 1;
        }
        dispositions
            .push(json!({"id":id,"return_code":code,"excluded_warmup":run["excluded_warmup"]}));
        if INTERRUPTED.load(Ordering::SeqCst) != 0 {
            break;
        }
    }
    let summary = json!({"schema":"n83.panel-completion/v1","finished_epoch_ms":epoch_ms(),"completed_arm_count":dispositions.len(),"planned_arm_count":runs.len(),"failed_arms":failures,"interrupted_signal":INTERRUPTED.load(Ordering::SeqCst),"dispositions":dispositions});
    let final_path = root.join(format!("panel_completion_{}.json", epoch_ms()));
    write_new(&final_path, &summary)?;
    println!("{}", summary);
    Ok(
        if failures == 0 && INTERRUPTED.load(Ordering::SeqCst) == 0 {
            0
        } else {
            1
        },
    )
}
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result=match args.first().map(String::as_str){Some("run-arm")=>run_arm(&args[1..]),Some("analyze")=>analyze(&args[1..]),Some("panel")=>panel(&args[1..]),_=>Err(err("usage: n83-evaluation-runner run-arm --arm ARM --out-dir NEW_DIR --target-file JSONL --cwd DIR --timeout-s SECONDS [--memory-bytes BYTES] [--cpu CPU] [--semantic-file JSONL] [--env KEY=VALUE] -- COMMAND [ARGS]\n       n83-evaluation-runner analyze --pairs JSONL --out NEW_JSON [--draws 20000] [--seed 2026101283]\n       n83-evaluation-runner panel --protocol JSON --out-dir NEW_DIR [--resume true]"))};
    match result {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            eprintln!("runner error: {e}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(name: &str, text: &str) -> PathBuf {
        let root = PathBuf::from("/workspace/n83_evaluation/runner/smoke");
        fs::create_dir_all(&root).unwrap();
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let p = root.join(format!("{name}_{}_{}.jsonl", std::process::id(), nonce));
        fs::write(&p, text).unwrap();
        p
    }
    fn summary() -> Value {
        json!({"kind":"compact_orbit_dlp_summary","rank":600,"orbit_columns":600,"targets":1,"targets_solved":1,"targets_failed":0,"base_hash":"expected"})
    }
    #[test]
    fn truncated_json_is_rejected() {
        let p = fixture("truncated", "{\"kind\":");
        assert_eq!(ic_summary_validate(&p, Some("expected"))["valid"], false);
    }
    #[test]
    fn missing_rank_is_rejected() {
        let mut s = summary();
        s.as_object_mut().unwrap().remove("rank");
        let p = fixture("missing_rank", &s.to_string());
        assert_eq!(ic_summary_validate(&p, Some("expected"))["valid"], false);
    }
    #[test]
    fn exit_zero_with_failed_target_is_rejected() {
        let mut s = summary();
        s["targets_failed"] = json!(1);
        s["exit_code"] = json!(0);
        let p = fixture("failed_target", &s.to_string());
        assert_eq!(ic_summary_validate(&p, Some("expected"))["valid"], false);
    }
    #[test]
    fn base_digest_mismatch_is_rejected() {
        let p = fixture("base_mismatch", &summary().to_string());
        assert_eq!(ic_summary_validate(&p, Some("wrong"))["valid"], false);
    }
    #[test]
    fn exactly_one_complete_summary_is_required() {
        let p = fixture("good_summary", &summary().to_string());
        assert_eq!(ic_summary_validate(&p, Some("expected"))["valid"], true);
        let p = fixture(
            "duplicate_summary",
            &format!("{}\n{}\n", summary(), summary()),
        );
        assert_eq!(ic_summary_validate(&p, Some("expected"))["valid"], false);
    }
    #[test]
    fn mismatched_public_point_is_rejected() {
        let r = json!({"published_q":["3","4"],"recovered_scalar":"9","group_verified":true,"online_total_ms":1.});
        let p = fixture("point_mismatch", &r.to_string());
        assert_eq!(sem_validate(&p, &json!(["1", "2"]), "rho")["valid"], false);
    }
    #[test]
    fn zero_replicate_and_one_replicate_ci_are_not_fabricated() {
        assert_eq!(stats(&[], 1, 1000)["n"], 0);
        let s = Sample {
            pair_id: "one".into(),
            kind: "matched".into(),
            log_ratio: 2f64.ln(),
            cold_log_ratio: 0.,
            baseline_ms: 2.,
            candidate_ms: 1.,
        };
        let r = stats(&[s], 1, 1000);
        assert_eq!(r["geometric_mean_paired_online_ratio"], 2.);
        assert_eq!(r["paired_percentile_bootstrap_95_ci"], Value::Null);
    }
}
