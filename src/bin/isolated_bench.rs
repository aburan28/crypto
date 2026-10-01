//! `isolated-bench`: run a benchmark on a reserved, pinned CPU and record
//! what else was running (AGENTS.md §10).
//!
//! The native form of `tools/isolated_bench.py`, with the same modes,
//! options, lock and record (`isolated-bench/1`), so the two interoperate:
//! they take one lock, and a run tree may hold records from either.  The
//! record says which wrote it (`tool`), and `host.python` is `null` here.
//!
//! Wall time is only evidence when the run had the core to itself.
//! Pinning the benchmark is not enough on its own: it keeps the benchmark
//! on one core but does not keep anything else off that core, and it does
//! nothing about builds or other jobs on the neighbouring cores that share
//! the cache and the memory bus.  This tool does all three, and records
//! the evidence:
//!
//! 1. **One benchmark at a time.**  It takes an exclusive `flock` on
//!    `--lock` (default `/tmp/crypto-bench.lock`).  Heavy work that is not
//!    a benchmark runs through `isolated-bench busy -- CMD`, which takes the
//!    same lock, so a build cannot start in the middle of a timed stage.
//! 2. **A quiet machine before starting.**  It samples every other
//!    process's CPU use for `--settle` seconds and refuses to start if they
//!    use more than `--max-other-cpu` CPUs' worth between them, or if the
//!    10-second CPU or memory pressure (PSI) is above `--max-psi`.
//! 3. **A reserved core.**  It moves every other thread it may move off the
//!    benchmark CPUs (`--cpus`), pins the benchmark there, and restores the
//!    other threads' affinity afterwards.  On a host with SMT it refuses a
//!    CPU whose sibling is not reserved too.  A process started while a run
//!    had moved its parent off those CPUs inherits the narrowed mask, and no
//!    restore reaches it; so when this process's own mask lacks the CPUs
//!    asked for, it widens the mask to them where the system allows, and the
//!    record says so (`affinity_widened_from`).
//! 4. **A record of the conditions.**  Per run: wall time, user and system
//!    CPU, voluntary and involuntary context switches, page faults, load
//!    average and PSI before and after, and the CPU time every other
//!    process used while the benchmark ran.  A run where other processes
//!    used more than `--max-other-cpu` is marked `contended`.
//!
//! It cannot see or stop other tenants of a virtual machine's host, and it
//! cannot fix the CPU frequency; an A/A run measures that residual noise.
//!
//! ```text
//! isolated-bench run --cpus 3 --out rec.jsonl -- ./worker
//! isolated-bench reserve --cpus 3 --out cond.json -- ./harness --cpu 3
//! isolated-bench busy -- cargo build --release
//! ```
//!
//! `run` pins the command itself.  `reserve` is for a harness that pins its
//! own children: it takes the lock, reserves the CPUs and monitors the
//! machine for the whole command, but leaves pinning to it.  Linux only:
//! it reads `/proc` and PSI and sets other threads' affinity.

#[path = "icprog/json.rs"]
#[allow(dead_code)]
mod json;

#[cfg(target_os = "linux")]
mod linux {
    use std::collections::{BTreeSet, HashMap};
    use std::fs::{File, OpenOptions};
    use std::os::fd::AsRawFd;
    use std::process::{Command, Stdio};
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Arc, Mutex};
    use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

    use super::json::{self, J};

    pub const DEFAULT_LOCK: &str = "/tmp/crypto-bench.lock";

    type Cpus = BTreeSet<usize>;

    /// A refusal or failure, printed as the Python tool's `SystemExit`.
    pub struct Exit(pub String);

    fn exit<T>(msg: impl Into<String>) -> Result<T, Exit> {
        Err(Exit(msg.into()))
    }

    pub fn parse_cpus(text: &str) -> Result<Cpus, String> {
        let mut cpus = Cpus::new();
        for part in text.split(',').map(str::trim).filter(|p| !p.is_empty()) {
            let num = |s: &str| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| format!("bad CPU `{s}`"))
            };
            match part.split_once('-') {
                Some((lo, hi)) => cpus.extend(num(lo)?..=num(hi)?),
                None => {
                    cpus.insert(num(part)?);
                }
            }
        }
        if cpus.is_empty() {
            return Err("no CPUs given".into());
        }
        Ok(cpus)
    }

    fn smt_siblings(cpu: usize) -> Cpus {
        let path = format!("/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list");
        std::fs::read_to_string(path)
            .ok()
            .and_then(|t| parse_cpus(t.trim()).ok())
            .unwrap_or_else(|| Cpus::from([cpu]))
    }

    // ── affinity ────────────────────────────────────────────────────

    fn get_affinity(tid: i32) -> Option<Cpus> {
        // SAFETY: cpu_set_t is plain data; the kernel fills `set`.
        unsafe {
            let mut set: libc::cpu_set_t = std::mem::zeroed();
            if libc::sched_getaffinity(tid, std::mem::size_of::<libc::cpu_set_t>(), &mut set) != 0 {
                return None;
            }
            Some(
                (0..libc::CPU_SETSIZE as usize)
                    .filter(|&c| libc::CPU_ISSET(c, &set))
                    .collect(),
            )
        }
    }

    fn set_affinity(tid: i32, cpus: &Cpus) -> bool {
        // SAFETY: as above; `set` is initialised before the call.
        unsafe {
            let mut set: libc::cpu_set_t = std::mem::zeroed();
            libc::CPU_ZERO(&mut set);
            for &c in cpus {
                libc::CPU_SET(c, &mut set);
            }
            libc::sched_setaffinity(tid, std::mem::size_of::<libc::cpu_set_t>(), &set) == 0
        }
    }

    /// Refuse a reservation that cannot work.  Return the inherited mask if
    /// it lacked some of `cpus` and was widened to them.
    ///
    /// An isolated run moves every other thread off its CPUs and restores
    /// the threads it moved; a process forked meanwhile inherits the
    /// narrowed mask and keeps it.  A harness started that way would be
    /// refused on every later run, so the mask is widened instead, where
    /// the system allows it.
    fn check_cpus(cpus: &Cpus) -> Result<Option<Vec<usize>>, Exit> {
        let mut allowed = get_affinity(0).unwrap_or_default();
        let mut inherited = None;
        if !cpus.is_subset(&allowed) {
            set_affinity(0, &allowed.union(cpus).copied().collect());
            let now = get_affinity(0).unwrap_or_default();
            let missing: Vec<usize> = cpus.difference(&now).copied().collect();
            if !missing.is_empty() {
                let allowed: Vec<usize> = allowed.into_iter().collect();
                return exit(format!(
                    "CPUs {missing:?} are outside this process affinity {allowed:?}"
                ));
            }
            inherited = Some(allowed.into_iter().collect());
            allowed = now;
        }
        if *cpus == allowed {
            return exit(
                "reserving every CPU leaves nowhere to move other work; leave at least one free",
            );
        }
        for &cpu in cpus {
            let stray: Vec<usize> = smt_siblings(cpu).difference(cpus).copied().collect();
            if !stray.is_empty() {
                return exit(format!(
                    "CPU {cpu} shares a core with {stray:?}; reserve the siblings too"
                ));
            }
        }
        Ok(inherited)
    }

    // ── the machine's state ─────────────────────────────────────────

    fn psi(kind: &str) -> J {
        let Ok(text) = std::fs::read_to_string(format!("/proc/pressure/{kind}")) else {
            return J::Null;
        };
        let mut out = Vec::new();
        for line in text.split('\n').filter(|l| !l.is_empty()) {
            let mut it = line.split_whitespace();
            let name = it.next().unwrap_or("").to_string();
            let fields = it
                .filter_map(|f| f.split_once('='))
                .map(|(k, v)| (k.to_string(), J::Float(v.parse().unwrap_or(f64::NAN))))
                .collect();
            out.push((name, J::Obj(fields)));
        }
        J::Obj(out)
    }

    fn psi_some_avg10(p: &J) -> Option<f64> {
        p.get("some")?.get("avg10")?.as_f64()
    }

    fn loadavg() -> J {
        let text = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
        J::Arr(
            text.split_whitespace()
                .take(3)
                .map(|x| J::Float(x.parse().unwrap_or(f64::NAN)))
                .collect(),
        )
    }

    fn unix_time() -> f64 {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0.0, |d| d.as_secs_f64())
    }

    fn conditions() -> J {
        J::Obj(vec![
            ("unix_time".into(), J::Float(unix_time())),
            ("loadavg".into(), loadavg()),
            ("psi_cpu".into(), psi("cpu")),
            ("psi_memory".into(), psi("memory")),
        ])
    }

    /// Each process's command, CPU ticks and parent, from one `/proc` scan,
    /// in the scan's order.
    struct Snapshot(Vec<(i32, String, u64, i32)>);

    fn cpu_snapshot() -> Snapshot {
        let mut out = Vec::new();
        for entry in std::fs::read_dir("/proc").into_iter().flatten().flatten() {
            let name = entry.file_name();
            let Some(pid) = name.to_str().and_then(|s| s.parse::<i32>().ok()) else {
                continue;
            };
            let Ok(stat) = std::fs::read_to_string(entry.path().join("stat")) else {
                continue;
            };
            let (Some(open), Some(close)) = (stat.find('('), stat.rfind(')')) else {
                continue;
            };
            let comm = stat[open + 1..close].to_string();
            let fields: Vec<&str> = stat
                .get(close + 2..)
                .unwrap_or("")
                .split_whitespace()
                .collect();
            let (Some(u), Some(s), Some(pp)) = (fields.get(11), fields.get(12), fields.get(1))
            else {
                continue;
            };
            let ticks = u.parse::<u64>().unwrap_or(0) + s.parse::<u64>().unwrap_or(0);
            out.push((pid, comm, ticks, pp.parse().unwrap_or(0)));
        }
        Snapshot(out)
    }

    impl Snapshot {
        fn ticks(&self) -> HashMap<i32, u64> {
            self.0.iter().map(|(p, _, t, _)| (*p, *t)).collect()
        }

        /// Every descendant of `root` in this scan.
        fn descendants(&self, root: i32) -> BTreeSet<i32> {
            let mut children: HashMap<i32, Vec<i32>> = HashMap::new();
            for (pid, _, _, ppid) in &self.0 {
                children.entry(*ppid).or_default().push(*pid);
            }
            let (mut out, mut stack) = (BTreeSet::new(), vec![root]);
            while let Some(pid) = stack.pop() {
                for &c in children.get(&pid).into_iter().flatten() {
                    if out.insert(c) {
                        stack.push(c);
                    }
                }
            }
            out
        }
    }

    fn clock_ticks() -> f64 {
        // SAFETY: sysconf has no preconditions.
        unsafe { libc::sysconf(libc::_SC_CLK_TCK) as f64 }
    }

    /// CPU seconds used by processes outside `exclude` between two scans,
    /// the heaviest first.
    fn other_use(
        before: &HashMap<i32, u64>,
        after: &Snapshot,
        exclude: &BTreeSet<i32>,
    ) -> Vec<(String, f64)> {
        let tick = clock_ticks();
        let mut used: Vec<(String, f64)> = after
            .0
            .iter()
            .filter(|(pid, ..)| !exclude.contains(pid))
            .filter_map(|(pid, comm, ticks, _)| {
                let delta = *ticks as i64 - before.get(pid).copied().unwrap_or(0) as i64;
                (delta > 0).then(|| (format!("{comm}[{pid}]"), delta as f64 / tick))
            })
            .collect();
        used.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("finite"));
        used
    }

    /// Python's `sum` of the seconds: the integer 0 when there are none.
    fn total(used: &[(String, f64)]) -> J {
        if used.is_empty() {
            J::Int(0)
        } else {
            J::Float(used.iter().fold(0.0, |acc, (_, s)| acc + s))
        }
    }

    fn top(used: &[(String, f64)], k: usize) -> J {
        J::Obj(
            used.iter()
                .take(k)
                .map(|(n, s)| (n.clone(), J::Float(*s)))
                .collect(),
        )
    }

    fn settle(seconds: f64, exclude: &BTreeSet<i32>) -> J {
        let before = cpu_snapshot().ticks();
        std::thread::sleep(Duration::from_secs_f64(seconds));
        let used = other_use(&before, &cpu_snapshot(), exclude);
        J::Obj(vec![
            ("seconds".into(), J::Float(seconds)),
            ("other_cpu_seconds".into(), total(&used)),
            ("top".into(), top(&used, 8)),
        ])
    }

    /// Every thread on the machine, as (process, thread).
    fn threads() -> Vec<(i32, i32)> {
        let mut out = Vec::new();
        for proc_entry in std::fs::read_dir("/proc").into_iter().flatten().flatten() {
            let Some(pid) = proc_entry
                .file_name()
                .to_str()
                .and_then(|s| s.parse::<i32>().ok())
            else {
                continue;
            };
            for t in std::fs::read_dir(proc_entry.path().join("task"))
                .into_iter()
                .flatten()
                .flatten()
            {
                if let Some(tid) = t.file_name().to_str().and_then(|s| s.parse::<i32>().ok()) {
                    out.push((pid, tid));
                }
            }
        }
        out
    }

    /// Move every other thread off `cpus`; return the affinities to restore.
    fn evict(cpus: &Cpus, keep: &BTreeSet<i32>) -> Vec<(i32, Cpus)> {
        let mut moved = Vec::new();
        for (pid, tid) in threads() {
            if keep.contains(&pid) {
                continue;
            }
            let Some(current) = get_affinity(tid) else {
                continue;
            };
            if current.is_disjoint(cpus) {
                continue;
            }
            let target: Cpus = current.difference(cpus).copied().collect();
            if target.is_empty() {
                continue;
            }
            if set_affinity(tid, &target) {
                moved.push((tid, current));
            }
        }
        moved
    }

    fn restore(moved: &[(i32, Cpus)]) {
        for (tid, affinity) in moved {
            set_affinity(*tid, affinity);
        }
    }

    /// Threads still allowed on the reserved CPUs after eviction: user
    /// threads that refused the move by name, per-CPU kernel threads counted.
    fn unmovable_on(cpus: &Cpus, keep: &BTreeSet<i32>) -> J {
        let (mut user, mut kernel) = (Vec::new(), 0i128);
        for (pid, tid) in threads() {
            if keep.contains(&pid) {
                continue;
            }
            let Some(aff) = get_affinity(tid) else {
                continue;
            };
            if aff.is_disjoint(cpus) {
                continue;
            }
            let base = format!("/proc/{pid}");
            match std::fs::read(format!("{base}/cmdline")) {
                Ok(cmd) if !cmd.is_empty() => {
                    let comm = std::fs::read_to_string(format!("{base}/task/{tid}/comm"))
                        .unwrap_or_default();
                    user.push(J::Str(format!("{}[{tid}]", comm.trim())));
                }
                Ok(_) => kernel += 1,
                Err(_) => continue,
            }
        }
        J::Obj(vec![
            ("user_threads".into(), J::Arr(user)),
            ("kernel_threads".into(), J::Int(kernel)),
        ])
    }

    // ── the lock ────────────────────────────────────────────────────

    pub fn lock(path: &str, wait: bool) -> Result<File, Exit> {
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .read(true)
            .open(path)
            .map_err(|e| Exit(format!("{path}: {e}")))?;
        let flags = if wait {
            libc::LOCK_EX
        } else {
            libc::LOCK_EX | libc::LOCK_NB
        };
        loop {
            // SAFETY: the descriptor is open for the call's duration.
            if unsafe { libc::flock(file.as_raw_fd(), flags) } == 0 {
                return Ok(file);
            }
            let err = std::io::Error::last_os_error();
            match err.raw_os_error() {
                Some(libc::EINTR) => continue,
                Some(libc::EWOULDBLOCK) => {
                    return exit(format!("another benchmark or busy job holds {path}"))
                }
                _ => return exit(format!("{path}: {err}")),
            }
        }
    }

    fn host() -> J {
        let model = std::fs::read_to_string("/proc/cpuinfo")
            .ok()
            .and_then(|t| {
                t.split('\n')
                    .find(|l| l.starts_with("model name"))
                    .and_then(|l| l.split_once(':'))
                    .map(|(_, v)| v.trim().to_string())
            })
            .unwrap_or_default();
        let kernel = std::fs::read_to_string("/proc/sys/kernel/osrelease")
            .map(|s| s.trim().to_string())
            .unwrap_or_default();
        let cpus = std::thread::available_parallelism().map_or(0, |n| n.get());
        // available_parallelism honours the affinity mask; the Python tool
        // reports os.cpu_count(), every CPU online.
        // SAFETY: sysconf has no preconditions.
        let online = unsafe { libc::sysconf(libc::_SC_NPROCESSORS_ONLN) };
        let logical = if online > 0 {
            online as i128
        } else {
            cpus as i128
        };
        J::Obj(vec![
            ("cpu_model".into(), J::Str(model)),
            ("logical_cpus".into(), J::Int(logical)),
            ("kernel".into(), J::Str(kernel)),
            ("machine".into(), J::Str(std::env::consts::ARCH.into())),
            ("python".into(), J::Null),
        ])
    }

    // ── the modes ───────────────────────────────────────────────────

    pub struct Opts {
        pub lock: String,
        pub wait: bool,
        pub cpus: String,
        pub out: String,
        pub settle: f64,
        pub max_other_cpu: f64,
        pub max_psi: f64,
        pub label: String,
        pub stdin: Option<String>,
        pub period: f64,
    }

    fn preflight(o: &Opts, exclude: &BTreeSet<i32>) -> Result<J, Exit> {
        let quiet = settle(o.settle, exclude);
        let now = conditions();
        let worst = ["psi_cpu", "psi_memory"]
            .iter()
            .filter_map(|k| now.get(k).and_then(psi_some_avg10))
            .fold(0.0f64, f64::max);
        let budget = o.max_other_cpu * o.settle;
        let other = quiet
            .get("other_cpu_seconds")
            .and_then(J::as_f64)
            .unwrap_or(0.0);
        if other > budget || worst > o.max_psi {
            return exit(format!(
                "machine is busy: other processes used {other:.2} CPU s in {} s (limit {budget:.2}), \
                 PSI some avg10 {} (limit {}); top: {}",
                json::py_float(o.settle),
                json::py_float(worst),
                json::py_float(o.max_psi),
                json::dumps_line(quiet.get("top").unwrap_or(&J::Null), false)
            ));
        }
        Ok(J::Obj(vec![
            ("settle".into(), quiet),
            ("conditions".into(), now),
        ]))
    }

    fn exit_code(status: i32) -> i32 {
        if libc::WIFEXITED(status) {
            libc::WEXITSTATUS(status)
        } else if libc::WIFSIGNALED(status) {
            -libc::WTERMSIG(status)
        } else {
            status
        }
    }

    fn secs(t: libc::timeval) -> f64 {
        t.tv_sec as f64 + t.tv_usec as f64 / 1e6
    }

    fn run_pinned(o: &Opts, cpus: &Cpus, command: &[String]) -> Result<J, Exit> {
        let me = std::process::id() as i32;
        let before_ticks = cpu_snapshot().ticks();
        let before = conditions();
        let inherited = get_affinity(0).unwrap_or_default();
        set_affinity(0, cpus);
        let start = Instant::now();
        let mut cmd = Command::new(&command[0]);
        cmd.args(&command[1..]);
        if let Some(path) = &o.stdin {
            let f = File::open(path).map_err(|e| Exit(format!("{path}: {e}")))?;
            cmd.stdin(Stdio::from(f));
        }
        let spawned = cmd.spawn();
        set_affinity(0, &inherited);
        let child = spawned.map_err(|e| Exit(format!("{}: {e}", command[0])))?;
        let pid = child.id() as i32;
        let mut status = 0;
        // SAFETY: rusage is plain data that wait4 fills.
        let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
        loop {
            // SAFETY: `pid` is our child; we reap it exactly once, here.
            let r = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
            if r == pid {
                break;
            }
            if std::io::Error::last_os_error().raw_os_error() != Some(libc::EINTR) {
                return exit(format!("wait4: {}", std::io::Error::last_os_error()));
            }
        }
        let wall = start.elapsed().as_secs_f64();
        let after = conditions();
        let exclude = BTreeSet::from([me, pid]);
        let others = other_use(&before_ticks, &cpu_snapshot(), &exclude);
        let other_seconds = others.iter().fold(0.0, |acc, (_, s)| acc + s);
        Ok(J::Obj(vec![
            (
                "command".into(),
                J::Arr(command.iter().map(|s| J::Str(s.clone())).collect()),
            ),
            (
                "cpus".into(),
                J::Arr(cpus.iter().map(|&c| J::Int(c as i128)).collect()),
            ),
            ("exit_status".into(), J::Int(exit_code(status).into())),
            ("wall_seconds".into(), J::Float(wall)),
            ("user_seconds".into(), J::Float(secs(usage.ru_utime))),
            ("system_seconds".into(), J::Float(secs(usage.ru_stime))),
            ("voluntary_switches".into(), J::Int(usage.ru_nvcsw.into())),
            (
                "involuntary_switches".into(),
                J::Int(usage.ru_nivcsw.into()),
            ),
            ("minor_faults".into(), J::Int(usage.ru_minflt.into())),
            ("major_faults".into(), J::Int(usage.ru_majflt.into())),
            ("max_rss_kib".into(), J::Int(usage.ru_maxrss.into())),
            ("before".into(), before),
            ("after".into(), after),
            ("other_cpu_seconds".into(), total(&others)),
            ("other_top".into(), top(&others, 8)),
            (
                "contended".into(),
                J::Bool(other_seconds > o.max_other_cpu * wall),
            ),
        ]))
    }

    /// Sample other processes' CPU use every `period` while a harness runs.
    fn monitor(
        period: f64,
        root: i32,
        threshold: f64,
        stop: Arc<AtomicBool>,
        samples: Arc<Mutex<Vec<J>>>,
    ) -> std::thread::JoinHandle<()> {
        std::thread::spawn(move || {
            let me = std::process::id() as i32;
            let mut previous = cpu_snapshot().ticks();
            let step = Duration::from_millis(50);
            loop {
                let deadline = Instant::now() + Duration::from_secs_f64(period);
                while Instant::now() < deadline {
                    if stop.load(Ordering::SeqCst) {
                        return;
                    }
                    std::thread::sleep(
                        step.min(deadline.saturating_duration_since(Instant::now())),
                    );
                }
                let snapshot = cpu_snapshot();
                let mut exclude = snapshot.descendants(root);
                exclude.extend([me, root]);
                let used = other_use(&previous, &snapshot, &exclude);
                let sum = used.iter().fold(0.0, |acc, (_, s)| acc + s);
                samples.lock().expect("samples").push(J::Obj(vec![
                    ("unix_time".into(), J::Float(unix_time())),
                    ("other_cpu_seconds".into(), total(&used)),
                    ("contended".into(), J::Bool(sum > threshold * period)),
                    ("top".into(), top(&used, 5)),
                    ("loadavg".into(), loadavg()),
                ]));
                previous = snapshot.ticks();
            }
        })
    }

    pub fn busy(lock_path: &str, command: &[String]) -> Result<i32, Exit> {
        let _guard = lock(lock_path, true)?;
        let status = Command::new(&command[0])
            .args(&command[1..])
            .status()
            .map_err(|e| Exit(format!("{}: {e}", command[0])))?;
        Ok(status.code().unwrap_or_else(|| {
            use std::os::unix::process::ExitStatusExt;
            -status.signal().unwrap_or(1)
        }))
    }

    pub fn reserved(mode: &str, o: &Opts, command: &[String]) -> Result<i32, Exit> {
        let cpus = parse_cpus(&o.cpus).map_err(Exit)?;
        let inherited = check_cpus(&cpus)?;
        let record;
        let code;
        {
            let _guard = lock(&o.lock, o.wait)?;
            // Only this process is exempt.  The shell and agent that
            // launched it are moved off the reserved CPUs and charged as
            // contention like anything else.
            let mine = BTreeSet::from([std::process::id() as i32]);
            let pre = preflight(o, &mine)?;
            let moved = evict(&cpus, &mine);
            let mut kv = vec![
                ("schema".to_string(), J::Str("isolated-bench/1".into())),
                ("mode".into(), J::Str(mode.into())),
                ("label".into(), J::Str(o.label.clone())),
                ("host".into(), host()),
                (
                    "reserved_cpus".into(),
                    J::Arr(cpus.iter().map(|&c| J::Int(c as i128)).collect()),
                ),
                ("threads_moved".into(), J::Int(moved.len() as i128)),
                ("left_on_reserved".into(), unmovable_on(&cpus, &mine)),
                ("preflight".into(), pre),
                (
                    "tool".into(),
                    J::Obj(vec![
                        ("name".into(), J::Str("isolated-bench".into())),
                        ("implementation".into(), J::Str("native".into())),
                        ("version".into(), J::Str(env!("CARGO_PKG_VERSION").into())),
                    ]),
                ),
            ];
            if let Some(mask) = &inherited {
                kv.push((
                    "affinity_widened_from".into(),
                    J::Arr(mask.iter().map(|&c| J::Int(c as i128)).collect()),
                ));
            }
            let outcome = if mode == "run" {
                run_pinned(o, &cpus, command).map(|run| {
                    let c = run.get("exit_status").and_then(J::as_i128).unwrap_or(1) as i32;
                    kv.push(("run".into(), run));
                    c
                })
            } else {
                let child = Command::new(&command[0])
                    .args(&command[1..])
                    .spawn()
                    .map_err(|e| Exit(format!("{}: {e}", command[0])));
                child.and_then(|mut child| {
                    let stop = Arc::new(AtomicBool::new(false));
                    let samples = Arc::new(Mutex::new(Vec::new()));
                    let handle = monitor(
                        o.period,
                        child.id() as i32,
                        o.max_other_cpu,
                        stop.clone(),
                        samples.clone(),
                    );
                    let status = child.wait().map_err(|e| Exit(format!("wait: {e}")))?;
                    stop.store(true, Ordering::SeqCst);
                    handle.join().expect("monitor");
                    let c = status.code().unwrap_or_else(|| {
                        use std::os::unix::process::ExitStatusExt;
                        -status.signal().unwrap_or(1)
                    });
                    let samples = std::mem::take(&mut *samples.lock().expect("samples"));
                    let contended = samples
                        .iter()
                        .filter(|s| s.get("contended").is_some_and(J::truthy))
                        .count();
                    kv.push((
                        "command".into(),
                        J::Arr(command.iter().map(|s| J::Str(s.clone())).collect()),
                    ));
                    kv.push(("exit_status".into(), J::Int(c.into())));
                    kv.push(("samples".into(), J::Arr(samples)));
                    kv.push(("contended_samples".into(), J::Int(contended as i128)));
                    Ok(c)
                })
            };
            restore(&moved);
            code = outcome?;
            record = J::Obj(kv);
        }
        let line = json::dumps_line(&record, true);
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&o.out)
            .map_err(|e| Exit(format!("{}: {e}", o.out)))?;
        use std::io::Write;
        writeln!(f, "{line}").map_err(|e| Exit(format!("{}: {e}", o.out)))?;
        Ok(code)
    }
}

#[cfg(target_os = "linux")]
fn main() -> std::process::ExitCode {
    use clap::{Args, Parser, Subcommand};

    #[derive(Args)]
    struct Common {
        /// The lock every benchmark and busy job takes.
        #[arg(long, default_value = linux::DEFAULT_LOCK)]
        lock: String,
        /// Wait for the lock instead of failing.
        #[arg(long)]
        wait: bool,
    }

    #[derive(Args)]
    struct Reserve {
        /// Benchmark CPUs, e.g. 3 or 2-3.
        #[arg(long)]
        cpus: String,
        /// JSON lines file the record is appended to.
        #[arg(long)]
        out: String,
        #[arg(long, default_value_t = 2.0)]
        settle: f64,
        /// Other processes may use at most this many CPUs on average.
        #[arg(long, default_value_t = 0.10)]
        max_other_cpu: f64,
        #[arg(long, default_value_t = 5.0)]
        max_psi: f64,
        #[arg(long, default_value = "")]
        label: String,
    }

    #[derive(Subcommand)]
    enum Mode {
        /// Run one command pinned to the reserved CPUs.
        Run {
            #[command(flatten)]
            common: Common,
            #[command(flatten)]
            reserve: Reserve,
            /// File fed to the command on standard input.
            #[arg(long)]
            stdin: Option<String>,
            #[arg(last = true, required = true)]
            command: Vec<String>,
        },
        /// Reserve the CPUs and monitor the machine for a harness that pins its own children.
        Reserve {
            #[command(flatten)]
            common: Common,
            #[command(flatten)]
            reserve: Reserve,
            #[arg(long, default_value_t = 5.0)]
            period: f64,
            #[arg(last = true, required = true)]
            command: Vec<String>,
        },
        /// Run heavy work that is not a benchmark under the same lock.
        Busy {
            #[command(flatten)]
            common: Common,
            #[arg(last = true, required = true)]
            command: Vec<String>,
        },
    }

    #[derive(Parser)]
    #[command(
        name = "isolated-bench",
        version,
        about = "Run a benchmark on a reserved, pinned CPU and record what else was running"
    )]
    struct Cli {
        #[command(subcommand)]
        mode: Mode,
    }

    let opts = |common: Common, r: Reserve, stdin: Option<String>, period: f64| linux::Opts {
        lock: common.lock,
        wait: common.wait,
        cpus: r.cpus,
        out: r.out,
        settle: r.settle,
        max_other_cpu: r.max_other_cpu,
        max_psi: r.max_psi,
        label: r.label,
        stdin,
        period,
    };
    let result = match Cli::parse().mode {
        Mode::Busy { common, command } => linux::busy(&common.lock, &command),
        Mode::Run {
            common,
            reserve,
            stdin,
            command,
        } => linux::reserved("run", &opts(common, reserve, stdin, 5.0), &command),
        Mode::Reserve {
            common,
            reserve,
            period,
            command,
        } => linux::reserved("reserve", &opts(common, reserve, None, period), &command),
    };
    match result {
        Ok(code) => std::process::ExitCode::from((code & 0xff) as u8),
        Err(linux::Exit(msg)) => {
            eprintln!("{msg}");
            std::process::ExitCode::FAILURE
        }
    }
}

#[cfg(not(target_os = "linux"))]
fn main() -> std::process::ExitCode {
    eprintln!("isolated-bench runs on Linux only: it reads /proc and PSI and sets other threads' CPU affinity");
    std::process::ExitCode::FAILURE
}
