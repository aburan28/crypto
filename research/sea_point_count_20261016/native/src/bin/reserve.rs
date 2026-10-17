//! Native equivalent of the repository isolation policy. No solver arithmetic.
//! Linux only: reserve one non-SMT CPU, evict movable threads, monitor /proc.
use isogeny_algos::json::Json;
use std::{
    collections::{BTreeMap, HashSet},
    fs::{self, OpenOptions},
    io::Write,
    os::fd::AsRawFd,
    process::{Command, Stdio},
    thread,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};
unsafe extern "C" {
    fn flock(fd: i32, op: i32) -> i32;
    fn sched_getaffinity(pid: i32, len: usize, mask: *mut u8) -> i32;
    fn sched_setaffinity(pid: i32, len: usize, mask: *const u8) -> i32;
    fn sysconf(name: i32) -> i64;
    fn getrusage(who: i32, usage: *mut i64) -> i32;
}
type Mask = [u8; 128];
fn affinity(pid: i32) -> Option<Mask> {
    let mut mask = [0u8; 128];
    (unsafe { sched_getaffinity(pid, 128, mask.as_mut_ptr()) } == 0).then_some(mask)
}
fn set(pid: i32, mask: &Mask) -> bool {
    (unsafe { sched_setaffinity(pid, 128, mask.as_ptr()) }) == 0
}
fn usage() -> [i64; 18] {
    let mut r = [0; 18];
    assert_eq!(unsafe { getrusage(-1, r.as_mut_ptr()) }, 0);
    r
}
fn psi_max() -> Option<f64> {
    let mut values = vec![];
    for name in ["cpu", "memory"] {
        if let Ok(s) = fs::read_to_string(format!("/proc/pressure/{name}")) {
            for line in s.lines().filter(|x| x.starts_with("some ")) {
                for v in line
                    .split_whitespace()
                    .filter_map(|x| x.strip_prefix("avg10="))
                {
                    if let Ok(v) = v.parse::<f64>() {
                        values.push(v);
                    }
                }
            }
        }
    }
    values.into_iter().reduce(f64::max)
}
fn snapshot() -> BTreeMap<i32, (String, u64)> {
    let mut out = BTreeMap::new();
    for e in fs::read_dir("/proc").unwrap().flatten() {
        let Ok(pid) = e.file_name().to_string_lossy().parse::<i32>() else {
            continue;
        };
        let Ok(s) = fs::read_to_string(e.path().join("stat")) else {
            continue;
        };
        let Some(end) = s.rfind(')') else { continue };
        let f: Vec<_> = s[end + 2..].split_whitespace().collect();
        if f.len() > 12 {
            out.insert(
                pid,
                (
                    s[s.find('(').unwrap() + 1..end].into(),
                    f[11].parse::<u64>().unwrap() + f[12].parse::<u64>().unwrap(),
                ),
            );
        }
    }
    out
}
fn other(
    before: &BTreeMap<i32, (String, u64)>,
    after: &BTreeMap<i32, (String, u64)>,
    exclude: &HashSet<i32>,
    hz: f64,
) -> (f64, Json) {
    let mut total = 0f64;
    let mut details = vec![];
    for (&pid, (name, ticks)) in after {
        if exclude.contains(&pid) {
            continue;
        }
        let delta = ticks.saturating_sub(before.get(&pid).map(|x| x.1).unwrap_or(0));
        if delta > 0 {
            let secs = delta as f64 / hz;
            total += secs;
            details.push((format!("{name}[{pid}]"), Json::Float(secs)));
        }
    }
    (total, Json::obj(details))
}
struct Reservation(Vec<(i32, Mask)>);
impl Drop for Reservation {
    fn drop(&mut self) {
        for (pid, mask) in &self.0 {
            set(*pid, mask);
        }
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: reserve BASELINE_BINARY CANDIDATE_BINARY OUTPUT_DIR".into());
    }
    let dir = std::path::Path::new(&args[3]);
    fs::create_dir_all(dir)?;
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .write(true)
        .open("/tmp/crypto-bench.lock")?;
    if unsafe { flock(lock.as_raw_fd(), 2 | 4) } != 0 {
        return Err("benchmark/busy lock held".into());
    }
    let me = std::process::id() as i32;
    let mask = affinity(0).ok_or("affinity unavailable")?;
    let cpus: Vec<usize> = (0..1024)
        .filter(|&c| mask[c / 8] & (1 << (c % 8)) != 0)
        .collect();
    if cpus.len() < 2 {
        return Err("need one reserved and one housekeeping CPU".into());
    }
    let cpu = *cpus.last().unwrap();
    let siblings = fs::read_to_string(format!(
        "/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"
    ))?;
    if siblings.trim() != cpu.to_string() {
        return Err("SMT reservation requires all siblings; select a non-SMT host".into());
    }
    let hz = unsafe { sysconf(2) } as f64;
    if hz <= 0. {
        return Err("clock tick rate unavailable".into());
    }
    let mut moved = Reservation(vec![]);
    let mut unmovable = vec![];
    for e in fs::read_dir("/proc")?.flatten() {
        let Ok(pid) = e.file_name().to_string_lossy().parse::<i32>() else {
            continue;
        };
        if pid == me {
            continue;
        }
        let Ok(tasks) = fs::read_dir(e.path().join("task")) else {
            continue;
        };
        for task in tasks.flatten() {
            let Ok(tid) = task.file_name().to_string_lossy().parse::<i32>() else {
                continue;
            };
            let Some(old) = affinity(tid) else { continue };
            if old[cpu / 8] & (1 << (cpu % 8)) == 0 {
                continue;
            }
            let mut target = old;
            target[cpu / 8] &= !(1 << (cpu % 8));
            if target.iter().all(|&v| v == 0) || !set(tid, &target) {
                unmovable.push(Json::Num(tid as i64));
            } else {
                moved.0.push((tid, old));
            }
        }
    }
    let mut pinned = [0u8; 128];
    pinned[cpu / 8] |= 1 << (cpu % 8);
    let mut housekeeping = mask;
    housekeeping[cpu / 8] &= !(1 << (cpu % 8));
    moved.0.push((me, mask));
    if !set(0, &housekeeping) {
        return Err("cannot move observer off benchmark CPU".into());
    }
    let psi_cpu = fs::read_to_string("/proc/pressure/cpu").ok();
    let psi_memory = fs::read_to_string("/proc/pressure/memory").ok();
    // Missing PSI is recorded, never silently treated as a successful observation.
    let host = Json::obj(vec![
        ("cpuinfo", Json::str(fs::read_to_string("/proc/cpuinfo")?)),
        ("meminfo", Json::str(fs::read_to_string("/proc/meminfo")?)),
        ("reserved_cpu", Json::Num(cpu as i64)),
        ("threads_moved", Json::Num(moved.0.len() as i64)),
        ("unmovable", Json::Arr(unmovable.clone())),
        ("psi_cpu", psi_cpu.map(Json::str).unwrap_or(Json::Null)),
        (
            "psi_memory",
            psi_memory.map(Json::str).unwrap_or(Json::Null),
        ),
        ("clock_ticks_per_second", Json::Float(hz)),
        ("max_other_cpu", Json::Float(0.10)),
        ("settle_seconds", Json::Num(2)),
        (
            "start_unix",
            Json::Float(SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64()),
        ),
    ]);
    fs::write(dir.join("host.json"), host.dump())?;
    let mut records = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(dir.join("runs.jsonl"))?;
    let cases = [
        "bsgs24",
        "bsgs32",
        "bsgs40",
        "sea40_cached",
        "sea192_cached",
        "count192_cold",
        "cm256_control",
    ];
    let deadline = Instant::now();
    for kind in ["aa", "ab"] {
        for case in cases {
            for round in 0..5 {
                for slot in 0..2 {
                    if deadline.elapsed() > Duration::from_secs(3600) {
                        return Err("suite budget exceeded".into());
                    }
                    let arm = if kind == "aa" || slot == 0 {
                        "baseline"
                    } else {
                        "candidate"
                    };
                    let binary = if arm == "baseline" {
                        &args[1]
                    } else {
                        &args[2]
                    };
                    let before_quiet = snapshot();
                    thread::sleep(Duration::from_secs(2));
                    let (quiet, _) = other(&before_quiet, &snapshot(), &HashSet::from([me]), hz);
                    let psi = psi_max();
                    if quiet > 0.20 || psi.is_some_and(|x| x > 5.) {
                        return Err(format!(
                            "preflight contention: {quiet} CPU seconds, PSI={psi:?}"
                        )
                        .into());
                    }
                    let before = snapshot();
                    let before_usage = usage();
                    let wall = Instant::now();
                    if !set(0, &pinned) {
                        return Err("cannot pin benchmark child".into());
                    }
                    let child = Command::new(binary)
                        .arg(case)
                        .env("RAYON_NUM_THREADS", "1")
                        .env("SEA_PHI_CACHE_DIR", dir.join("phi-cache"))
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn();
                    set(0, &housekeeping);
                    let mut child = child?;
                    let child_pid = child.id() as i32;
                    let mut killed = false;
                    loop {
                        if child.try_wait()?.is_some() {
                            break;
                        }
                        if wall.elapsed() > Duration::from_secs(120) {
                            child.kill()?;
                            killed = true;
                            break;
                        }
                        thread::sleep(Duration::from_millis(10));
                    }
                    let output = child.wait_with_output()?;
                    let after_usage = usage();
                    let wall_seconds = wall.elapsed().as_secs_f64();
                    let (used, details) =
                        other(&before, &snapshot(), &HashSet::from([me, child_pid]), hz);
                    let text = String::from_utf8_lossy(&output.stdout);
                    let result = Json::parse(text.trim()).unwrap_or(Json::Null);
                    let rec = Json::obj(vec![
                        ("kind", Json::str(kind)),
                        ("case", Json::str(case)),
                        ("round", Json::Num(round)),
                        ("slot", Json::Num(slot)),
                        ("arm", Json::str(arm)),
                        ("result", result),
                        ("wall_seconds", Json::Float(wall_seconds)),
                        ("other_cpu_seconds", Json::Float(used)),
                        ("other_processes", details),
                        ("quiet_cpu_seconds", Json::Float(quiet)),
                        ("psi_max_avg10", psi.map(Json::Float).unwrap_or(Json::Null)),
                        (
                            "user_seconds",
                            Json::Float(
                                (after_usage[0] - before_usage[0]) as f64
                                    + (after_usage[1] - before_usage[1]) as f64 / 1e6,
                            ),
                        ),
                        (
                            "system_seconds",
                            Json::Float(
                                (after_usage[2] - before_usage[2]) as f64
                                    + (after_usage[3] - before_usage[3]) as f64 / 1e6,
                            ),
                        ),
                        ("minor_faults", Json::Num(after_usage[8] - before_usage[8])),
                        ("major_faults", Json::Num(after_usage[9] - before_usage[9])),
                        (
                            "voluntary_switches",
                            Json::Num(after_usage[16] - before_usage[16]),
                        ),
                        (
                            "involuntary_switches",
                            Json::Num(after_usage[17] - before_usage[17]),
                        ),
                        ("contended", Json::Bool(used > 0.10 * wall_seconds)),
                        ("timeout", Json::Bool(killed)),
                        (
                            "exit_code",
                            Json::Num(output.status.code().unwrap_or(-1) as i64),
                        ),
                        ("stderr", Json::str(String::from_utf8_lossy(&output.stderr))),
                    ]);
                    writeln!(records, "{}", rec.dump())?;
                    records.flush()?;
                    println!(
                        "{kind} {case} round={round} slot={slot} arm={arm} exit={} contended={}",
                        output.status,
                        used > 0.10 * wall_seconds
                    );
                    if !output.status.success() || killed {
                        return Err("measurement failed; raw record retained".into());
                    }
                }
            }
        }
    }
    Ok(())
}
