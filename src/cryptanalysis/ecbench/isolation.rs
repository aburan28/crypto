//! Isolation: reserve CPUs, keep everything else off them, pin the
//! measured process and its memory, and record what the kernel saw.
//!
//! This is AGENTS.md §10 in native code, for a harness that spawns its
//! own measured children.  It follows `tools/isolated_bench.py` (same
//! lock file, same preflight defaults, same refusals: no SMT sibling
//! left unreserved, never every CPU, never CPU 0's core by default) and
//! adds what `docs/ic/measurement/README.md` §7 found the `contended`
//! flag misses: hypervisor steal on the pinned CPU, foreign busy time on
//! it during the run, and memory placed on another NUMA node.
//!
//! Everything here reads `/proc` and `/sys` and calls
//! `sched_setaffinity`/`set_mempolicy`, so it is Linux only.  On any
//! other system each function says it cannot, and the runner records
//! that a run earned level L0 and why; it never pretends.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// One logical CPU's place in the machine.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CpuPlace {
    pub cpu: u32,
    pub core_id: Option<u32>,
    pub package: Option<u32>,
    pub node: Option<u32>,
    /// Logical CPUs sharing this CPU's core, itself included.
    pub siblings: Vec<u32>,
}

/// One NUMA node.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NumaNode {
    pub node: u32,
    pub cpus: Vec<u32>,
    pub mem_total_kib: Option<u64>,
    /// `/sys/devices/system/node/nodeN/distance`, one entry per node.
    pub distances: Vec<u32>,
}

/// The CPU and memory topology, as far as the kernel exposes it.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Topology {
    pub cpus: Vec<CpuPlace>,
    pub nodes: Vec<NumaNode>,
}

impl Topology {
    pub fn place(&self, cpu: u32) -> Option<&CpuPlace> {
        self.cpus.iter().find(|c| c.cpu == cpu)
    }
}

/// `"0-3,8,10-11"` → `[0, 1, 2, 3, 8, 10, 11]`.
pub fn parse_cpulist(s: &str) -> Vec<u32> {
    let mut out = Vec::new();
    for part in s.trim().split(',').filter(|p| !p.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => {
                if let (Ok(a), Ok(b)) = (a.trim().parse::<u32>(), b.trim().parse::<u32>()) {
                    out.extend(a..=b);
                }
            }
            None => {
                if let Ok(v) = part.trim().parse() {
                    out.push(v);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

#[cfg(target_os = "linux")]
fn read_trim(path: &str) -> Option<String> {
    std::fs::read_to_string(path)
        .ok()
        .map(|s| s.trim().to_string())
}

/// Read the topology from `/sys`.  Empty on systems without it.
pub fn topology() -> Topology {
    #[cfg(target_os = "linux")]
    {
        let mut nodes = Vec::new();
        let mut node_of: BTreeMap<u32, u32> = BTreeMap::new();
        if let Ok(rd) = std::fs::read_dir("/sys/devices/system/node") {
            let mut ids: Vec<u32> = rd
                .filter_map(|e| e.ok())
                .filter_map(|e| e.file_name().to_str()?.strip_prefix("node")?.parse().ok())
                .collect();
            ids.sort_unstable();
            for n in ids {
                let base = format!("/sys/devices/system/node/node{n}");
                let cpus = read_trim(&format!("{base}/cpulist"))
                    .map(|s| parse_cpulist(&s))
                    .unwrap_or_default();
                for &c in &cpus {
                    node_of.insert(c, n);
                }
                let mem_total_kib = read_trim(&format!("{base}/meminfo")).and_then(|m| {
                    m.lines()
                        .find(|l| l.contains("MemTotal:"))
                        .and_then(|l| l.split_whitespace().rev().nth(1)?.parse().ok())
                });
                let distances = read_trim(&format!("{base}/distance"))
                    .map(|d| {
                        d.split_whitespace()
                            .filter_map(|v| v.parse().ok())
                            .collect()
                    })
                    .unwrap_or_default();
                nodes.push(NumaNode {
                    node: n,
                    cpus,
                    mem_total_kib,
                    distances,
                });
            }
        }
        let online = read_trim("/sys/devices/system/cpu/online")
            .map(|s| parse_cpulist(&s))
            .unwrap_or_default();
        let cpus = online
            .into_iter()
            .map(|cpu| {
                let base = format!("/sys/devices/system/cpu/cpu{cpu}/topology");
                CpuPlace {
                    cpu,
                    core_id: read_trim(&format!("{base}/core_id")).and_then(|s| s.parse().ok()),
                    package: read_trim(&format!("{base}/physical_package_id"))
                        .and_then(|s| s.parse().ok()),
                    node: node_of.get(&cpu).copied(),
                    siblings: read_trim(&format!("{base}/thread_siblings_list"))
                        .map(|s| parse_cpulist(&s))
                        .unwrap_or_else(|| vec![cpu]),
                }
            })
            .collect();
        Topology { cpus, nodes }
    }
    #[cfg(not(target_os = "linux"))]
    {
        Topology::default()
    }
}

// ── The CPU request ────────────────────────────────────────────────

/// What a session asks for.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CpuRequest {
    /// The highest-numbered whole core not shared with CPU 0, all its
    /// SMT siblings reserved; the measured process runs on the first.
    Auto,
    /// The CPUs this process is already allowed (an isolab placement or
    /// a cgroup cpuset reserved them); run on the first.
    Inherit,
    /// These CPUs, reserved; run on the first.  Every SMT sibling of
    /// each must be in the list.
    Explicit(Vec<u32>),
    /// No reservation at all: the run is recorded and earns L0.
    None,
}

impl CpuRequest {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "auto" => Ok(CpuRequest::Auto),
            "inherit" => Ok(CpuRequest::Inherit),
            "none" => Ok(CpuRequest::None),
            list => {
                let cpus = parse_cpulist(list);
                if cpus.is_empty() {
                    Err(format!(
                        "--cpus takes auto, inherit, none or a CPU list, not `{list}`"
                    ))
                } else {
                    Ok(CpuRequest::Explicit(cpus))
                }
            }
        }
    }
}

/// The CPUs a session holds.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CpuPlan {
    /// The CPU the measured process is pinned to.
    pub run_cpu: u32,
    /// Every CPU held, the run CPU and its idle siblings included.
    pub reserved: Vec<u32>,
    /// The NUMA node the run CPU belongs to, and the one memory is bound
    /// to when the system has more than one node.
    pub node: Option<u32>,
    /// How the plan was chosen.
    pub request: CpuRequest,
}

/// The CPUs this process may run on.
pub fn own_affinity() -> Option<Vec<u32>> {
    #[cfg(target_os = "linux")]
    {
        // SAFETY: the set is zeroed and sized for the call.
        unsafe {
            let mut set: libc::cpu_set_t = std::mem::zeroed();
            if libc::sched_getaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &mut set) != 0 {
                return None;
            }
            Some(
                (0..libc::CPU_SETSIZE as u32)
                    .filter(|&c| libc::CPU_ISSET(c as usize, &set))
                    .collect(),
            )
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

/// Decide which CPUs to hold.  Refuses what would make the measurement
/// meaningless: a sibling left to a neighbour, every CPU, CPU 0's core
/// (where the kernel's housekeeping lands) unless explicitly listed.
pub fn plan_cpus(req: &CpuRequest, topo: &Topology) -> Result<Option<CpuPlan>, String> {
    if *req == CpuRequest::None {
        return Ok(None);
    }
    if topo.cpus.is_empty() {
        return Err(
            "this system exposes no CPU topology (not Linux?): run with --cpus none, which records every run at L0"
                .into(),
        );
    }
    let all: Vec<u32> = topo.cpus.iter().map(|c| c.cpu).collect();
    let siblings_of = |cpu: u32| -> Vec<u32> {
        topo.place(cpu)
            .map(|p| p.siblings.clone())
            .unwrap_or_else(|| vec![cpu])
    };
    let reserved = match req {
        CpuRequest::Auto => {
            let zero_core = siblings_of(0);
            let mut best: Option<Vec<u32>> = None;
            for c in all.iter().rev() {
                let sib = siblings_of(*c);
                if sib.iter().any(|s| zero_core.contains(s)) {
                    continue;
                }
                best = Some(sib);
                break;
            }
            best.ok_or("no core other than CPU 0's is online")?
        }
        CpuRequest::Inherit => {
            let mut own = own_affinity().ok_or("cannot read this process's affinity")?;
            own.sort_unstable();
            own
        }
        CpuRequest::Explicit(list) => {
            for c in list {
                if !all.contains(c) {
                    return Err(format!("CPU {c} is not online"));
                }
                for s in siblings_of(*c) {
                    if !list.contains(&s) {
                        return Err(format!(
                            "CPU {c}'s SMT sibling {s} is not reserved; a neighbour on it shares the core's pipeline and caches"
                        ));
                    }
                }
            }
            list.clone()
        }
        CpuRequest::None => unreachable!(),
    };
    if reserved.is_empty() {
        return Err("no CPU to reserve".into());
    }
    if *req != CpuRequest::Inherit && reserved.len() >= all.len() {
        return Err("reserving every CPU leaves nowhere to move the rest of the system".into());
    }
    let run_cpu = reserved[0];
    Ok(Some(CpuPlan {
        run_cpu,
        node: topo.place(run_cpu).and_then(|p| p.node),
        reserved,
        request: req.clone(),
    }))
}

// ── The lock ───────────────────────────────────────────────────────

/// The lock `tools/isolated_bench.py` and `src/bin/isolated_bench.rs`
/// take, so a session and those tools never time at once.
pub const DEFAULT_LOCK: &str = "/tmp/crypto-bench.lock";

/// An exclusive `flock`, held until dropped.
pub struct BenchLock {
    _file: std::fs::File,
    pub path: String,
}

impl BenchLock {
    pub fn acquire(path: &str, wait: bool) -> Result<Self, String> {
        let file = std::fs::OpenOptions::new()
            .create(true)
            .truncate(false)
            .write(true)
            .open(path)
            .map_err(|e| format!("cannot open lock {path}: {e}"))?;
        #[cfg(unix)]
        {
            use std::os::unix::io::AsRawFd;
            let op = libc::LOCK_EX | if wait { 0 } else { libc::LOCK_NB };
            // SAFETY: a valid descriptor owned by `file`.
            if unsafe { libc::flock(file.as_raw_fd(), op) } != 0 {
                return Err(format!(
                    "another benchmark holds {path}; pass --wait to queue behind it"
                ));
            }
        }
        Ok(Self {
            _file: file,
            path: path.to_string(),
        })
    }
}

// ── Kernel counters ────────────────────────────────────────────────

/// One CPU's line of `/proc/stat`, in clock ticks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CpuTicks {
    /// user + nice + system + irq + softirq: time something ran.
    pub busy: u64,
    pub idle: u64,
    /// Time the hypervisor ran something else while this vCPU wanted to run.
    pub steal: u64,
}

/// Ticks per CPU from `/proc/stat`.
pub fn cpu_ticks() -> BTreeMap<u32, CpuTicks> {
    #[allow(unused_mut)]
    let mut out = BTreeMap::new();
    #[cfg(target_os = "linux")]
    if let Ok(stat) = std::fs::read_to_string("/proc/stat") {
        for line in stat.lines() {
            let mut it = line.split_whitespace();
            let Some(name) = it.next() else { continue };
            let Some(cpu) = name.strip_prefix("cpu").and_then(|n| n.parse::<u32>().ok()) else {
                continue;
            };
            let v: Vec<u64> = it.filter_map(|x| x.parse().ok()).collect();
            if v.len() < 8 {
                continue;
            }
            // user nice system idle iowait irq softirq steal
            out.insert(
                cpu,
                CpuTicks {
                    busy: v[0] + v[1] + v[2] + v[5] + v[6],
                    idle: v[3] + v[4],
                    steal: v[7],
                },
            );
        }
    }
    out
}

/// Clock ticks per second (`_SC_CLK_TCK`).
pub fn ticks_per_second() -> u64 {
    #[cfg(unix)]
    {
        // SAFETY: sysconf has no preconditions.
        let v = unsafe { libc::sysconf(libc::_SC_CLK_TCK) };
        if v > 0 {
            return v as u64;
        }
    }
    100
}

/// `/proc/pressure/<resource>`'s `some` line: `(avg10, total µs)`.
pub fn psi_some(resource: &str) -> Option<(f64, u64)> {
    #[cfg(target_os = "linux")]
    {
        let text = std::fs::read_to_string(format!("/proc/pressure/{resource}")).ok()?;
        let line = text.lines().find(|l| l.starts_with("some"))?;
        let mut avg10 = None;
        let mut total = None;
        for kv in line.split_whitespace().skip(1) {
            match kv.split_once('=') {
                Some(("avg10", v)) => avg10 = v.parse().ok(),
                Some(("total", v)) => total = v.parse().ok(),
                _ => {}
            }
        }
        Some((avg10?, total?))
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = resource;
        None
    }
}

/// The 1-, 5- and 15-minute load averages.
pub fn loadavg() -> Option<[f64; 3]> {
    let mut v = [0f64; 3];
    #[cfg(unix)]
    {
        // SAFETY: three doubles, as getloadavg expects.
        let n = unsafe { libc::getloadavg(v.as_mut_ptr(), 3) };
        if n == 3 {
            return Some(v);
        }
    }
    None
}

/// The volatile state around a run: never hashed, always recorded.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Conditions {
    pub unix_ms: u128,
    pub loadavg: Option<[f64; 3]>,
    pub psi_cpu_some_avg10: Option<f64>,
    pub psi_cpu_some_total_us: Option<u64>,
    pub psi_memory_some_avg10: Option<f64>,
    pub psi_memory_some_total_us: Option<u64>,
    /// `/proc/stat` ticks of the reserved CPUs.
    pub reserved_ticks: BTreeMap<u32, CpuTicks>,
}

pub fn conditions(reserved: &[u32]) -> Conditions {
    let ticks = cpu_ticks();
    let cpu = psi_some("cpu");
    let mem = psi_some("memory");
    Conditions {
        unix_ms: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0),
        loadavg: loadavg(),
        psi_cpu_some_avg10: cpu.map(|c| c.0),
        psi_cpu_some_total_us: cpu.map(|c| c.1),
        psi_memory_some_avg10: mem.map(|m| m.0),
        psi_memory_some_total_us: mem.map(|m| m.1),
        reserved_ticks: reserved
            .iter()
            .filter_map(|c| ticks.get(c).map(|t| (*c, *t)))
            .collect(),
    }
}

// ── Preflight ──────────────────────────────────────────────────────

/// Thresholds, `isolated_bench`'s defaults.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thresholds {
    pub settle_seconds: f64,
    /// Busy fraction of a reserved CPU, per second of settle, above which
    /// the host is not quiet.
    pub max_other_cpu: f64,
    /// PSI `some avg10` ceiling, percent.
    pub max_psi: f64,
    /// Involuntary context switches per second of the measured process.
    pub max_preemptions_per_second: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            settle_seconds: 2.0,
            max_other_cpu: 0.10,
            max_psi: 5.0,
            max_preemptions_per_second: 50.0,
        }
    }
}

/// What the preflight saw.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Preflight {
    pub quiet: bool,
    pub reasons: Vec<String>,
    pub settle_seconds: f64,
    /// Busy seconds per reserved CPU during the settle window.
    pub reserved_busy_seconds: BTreeMap<u32, f64>,
    pub before: Conditions,
    pub after: Conditions,
}

/// Watch the reserved CPUs for the settle window; quiet when none was
/// busy past the threshold and PSI is below its ceiling.
pub fn preflight(plan: &CpuPlan, th: &Thresholds) -> Preflight {
    let before = conditions(&plan.reserved);
    std::thread::sleep(std::time::Duration::from_secs_f64(th.settle_seconds));
    let after = conditions(&plan.reserved);
    let hz = ticks_per_second() as f64;
    let mut reasons = Vec::new();
    let mut busy = BTreeMap::new();
    for c in &plan.reserved {
        match (before.reserved_ticks.get(c), after.reserved_ticks.get(c)) {
            (Some(a), Some(b)) => {
                let s = b.busy.saturating_sub(a.busy) as f64 / hz;
                busy.insert(*c, s);
                if s > th.max_other_cpu * th.settle_seconds {
                    reasons.push(format!(
                        "CPU {c} was busy {s:.3}s of the {:.1}s settle",
                        th.settle_seconds
                    ));
                }
            }
            _ => reasons.push(format!("no /proc/stat line for CPU {c}")),
        }
    }
    for (name, v) in [
        ("cpu", after.psi_cpu_some_avg10),
        ("memory", after.psi_memory_some_avg10),
    ] {
        match v {
            Some(v) if v > th.max_psi => {
                reasons.push(format!("PSI {name} some avg10 {v:.2} > {:.2}", th.max_psi))
            }
            Some(_) => {}
            None => reasons.push(format!("PSI {name} unavailable")),
        }
    }
    Preflight {
        quiet: reasons.is_empty(),
        reasons,
        settle_seconds: th.settle_seconds,
        reserved_busy_seconds: busy,
        before,
        after,
    }
}

// ── Eviction ───────────────────────────────────────────────────────

/// Threads moved off the reserved CPUs, to be put back.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Eviction {
    /// `(tid, original mask)`.
    #[serde(skip)]
    #[cfg_attr(not(target_os = "linux"), allow(dead_code))]
    moved: Vec<(i32, Vec<u32>)>,
    pub moved_count: u64,
    /// Threads whose mask could not be changed (another user's, a
    /// per-CPU kernel thread).
    pub failed_user_threads: u64,
    pub failed_kernel_threads: u64,
}

#[cfg(target_os = "linux")]
fn set_affinity(tid: i32, cpus: &[u32]) -> bool {
    // SAFETY: the set is zeroed and filled with valid indices.
    unsafe {
        let mut set: libc::cpu_set_t = std::mem::zeroed();
        for &c in cpus {
            libc::CPU_SET(c as usize, &mut set);
        }
        libc::sched_setaffinity(tid, std::mem::size_of::<libc::cpu_set_t>(), &set) == 0
    }
}

#[cfg(target_os = "linux")]
fn affinity_of(tid: i32) -> Option<Vec<u32>> {
    // SAFETY: as in own_affinity.
    unsafe {
        let mut set: libc::cpu_set_t = std::mem::zeroed();
        if libc::sched_getaffinity(tid, std::mem::size_of::<libc::cpu_set_t>(), &mut set) != 0 {
            return None;
        }
        Some(
            (0..libc::CPU_SETSIZE as u32)
                .filter(|&c| libc::CPU_ISSET(c as usize, &set))
                .collect(),
        )
    }
}

/// Move every other thread's affinity off `reserved`.  Kernel threads
/// (no `/proc/<pid>/exe`) bound to a CPU cannot move and are counted.
pub fn evict(reserved: &[u32]) -> Eviction {
    #[allow(unused_mut)]
    let mut ev = Eviction::default();
    #[cfg(target_os = "linux")]
    {
        let me = std::process::id() as i32;
        let Ok(procs) = std::fs::read_dir("/proc") else {
            return ev;
        };
        for p in procs.filter_map(|e| e.ok()) {
            let Some(pid) = p.file_name().to_str().and_then(|s| s.parse::<i32>().ok()) else {
                continue;
            };
            if pid == me {
                continue;
            }
            let kernel = std::fs::read_link(format!("/proc/{pid}/exe")).is_err()
                && std::fs::read_to_string(format!("/proc/{pid}/cmdline"))
                    .map(|c| c.is_empty())
                    .unwrap_or(true);
            let Ok(tasks) = std::fs::read_dir(format!("/proc/{pid}/task")) else {
                continue;
            };
            for t in tasks.filter_map(|e| e.ok()) {
                let Some(tid) = t.file_name().to_str().and_then(|s| s.parse::<i32>().ok()) else {
                    continue;
                };
                let Some(mask) = affinity_of(tid) else {
                    continue;
                };
                if !mask.iter().any(|c| reserved.contains(c)) {
                    continue;
                }
                let rest: Vec<u32> = mask
                    .iter()
                    .copied()
                    .filter(|c| !reserved.contains(c))
                    .collect();
                if !rest.is_empty() && set_affinity(tid, &rest) {
                    ev.moved.push((tid, mask));
                    ev.moved_count += 1;
                } else if kernel {
                    ev.failed_kernel_threads += 1;
                } else {
                    ev.failed_user_threads += 1;
                }
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = reserved;
    }
    ev
}

impl Eviction {
    /// Put every moved thread's mask back (threads that exited are gone).
    pub fn restore(&mut self) {
        #[cfg(target_os = "linux")]
        for (tid, mask) in self.moved.drain(..) {
            let _ = set_affinity(tid, &mask);
        }
    }
}

/// Keep this (the runner) process off the reserved CPUs, so polling and
/// record writing never land on the measured core.
pub fn leave_reserved(reserved: &[u32]) -> Result<(), String> {
    #[cfg(target_os = "linux")]
    {
        let own = own_affinity().ok_or("cannot read own affinity")?;
        let rest: Vec<u32> = own.into_iter().filter(|c| !reserved.contains(c)).collect();
        if rest.is_empty() {
            return Err("the runner would have no CPU left outside the reservation".into());
        }
        if !set_affinity(0, &rest) {
            return Err("sched_setaffinity on the runner failed".into());
        }
        Ok(())
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = reserved;
        Err("CPU affinity is not available on this OS".into())
    }
}

/// Pin the calling thread to `cpu` and bind its memory to `node`.  Called
/// in the child between fork and exec, so it uses system calls only and
/// allocates nothing.
///
/// # Safety
/// Must be async-signal-safe: no allocation, no locks.
#[cfg(target_os = "linux")]
pub unsafe fn pin_in_child(cpu: u32, node: Option<u32>) -> std::io::Result<()> {
    let mut set: libc::cpu_set_t = std::mem::zeroed();
    libc::CPU_SET(cpu as usize, &mut set);
    if libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &set) != 0 {
        return Err(std::io::Error::last_os_error());
    }
    if let Some(n) = node {
        if n < 128 {
            // MPOL_BIND = 2: allocate only on this node.
            let mut mask = [0u64; 2];
            mask[(n / 64) as usize] |= 1u64 << (n % 64);
            let rc = libc::syscall(libc::SYS_set_mempolicy, 2i64, mask.as_ptr(), 128u64 + 1);
            if rc != 0 {
                return Err(std::io::Error::last_os_error());
            }
        }
    }
    Ok(())
}

/// What the measured process reports about its own placement.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SelfPlacement {
    /// Its affinity mask as read back inside the process.
    pub affinity: Option<Vec<u32>>,
    pub cpu_at_start: Option<u32>,
    pub cpu_at_end: Option<u32>,
    /// `/proc/self/status` `Mems_allowed_list`.
    pub mems_allowed: Option<String>,
}

/// The CPU this thread is running on.
pub fn current_cpu() -> Option<u32> {
    #[cfg(target_os = "linux")]
    {
        // SAFETY: no preconditions.
        let c = unsafe { libc::sched_getcpu() };
        if c >= 0 {
            return Some(c as u32);
        }
    }
    None
}

pub fn mems_allowed() -> Option<String> {
    #[cfg(target_os = "linux")]
    {
        let s = std::fs::read_to_string("/proc/self/status").ok()?;
        return s
            .lines()
            .find(|l| l.starts_with("Mems_allowed_list:"))
            .map(|l| l["Mems_allowed_list:".len()..].trim().to_string());
    }
    #[allow(unreachable_code)]
    None
}

// ── Host measurement settings ──────────────────────────────────────

/// The settings L3 asks of the host, for one CPU.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CpuSettings {
    pub governor: Option<String>,
    /// `intel_pstate/no_turbo` (`1` is off).
    pub no_turbo: Option<String>,
    /// `cpufreq/boost` (`0` is off).
    pub boost: Option<String>,
    pub in_isolcpus: Option<bool>,
    pub in_nohz_full: Option<bool>,
    /// `smt/control`: `on`, `off`, `forceoff`, `notsupported`, ...
    pub smt_control: Option<String>,
}

pub fn cpu_settings(cpu: u32) -> CpuSettings {
    #[cfg(target_os = "linux")]
    {
        let member = |path: &str| read_trim(path).map(|s| parse_cpulist(&s).contains(&cpu));
        CpuSettings {
            governor: read_trim(&format!(
                "/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor"
            )),
            no_turbo: read_trim("/sys/devices/system/cpu/intel_pstate/no_turbo"),
            boost: read_trim("/sys/devices/system/cpu/cpufreq/boost"),
            in_isolcpus: member("/sys/devices/system/cpu/isolated"),
            in_nohz_full: member("/sys/devices/system/cpu/nohz_full"),
            smt_control: read_trim("/sys/devices/system/cpu/smt/control"),
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = cpu;
        CpuSettings::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpulists_parse() {
        assert_eq!(parse_cpulist("0-3,8,10-11"), vec![0, 1, 2, 3, 8, 10, 11]);
        assert_eq!(parse_cpulist(""), Vec::<u32>::new());
        assert_eq!(parse_cpulist("5\n"), vec![5]);
    }

    fn two_node_smt() -> Topology {
        // 4 cores × 2 threads, two nodes; siblings (c, c+4).
        let cpus = (0..8)
            .map(|c| CpuPlace {
                cpu: c,
                core_id: Some(c % 4),
                package: Some(0),
                node: Some(if c % 4 < 2 { 0 } else { 1 }),
                siblings: vec![c % 4, c % 4 + 4],
            })
            .collect();
        Topology {
            cpus,
            nodes: vec![],
        }
    }

    #[test]
    fn auto_takes_a_whole_core_away_from_cpu0() {
        let plan = plan_cpus(&CpuRequest::Auto, &two_node_smt())
            .unwrap()
            .unwrap();
        assert_eq!(plan.reserved, vec![3, 7]);
        assert_eq!(plan.run_cpu, 3);
        assert_eq!(plan.node, Some(1));
    }

    #[test]
    fn a_lone_sibling_and_every_cpu_are_refused() {
        let t = two_node_smt();
        assert!(plan_cpus(&CpuRequest::Explicit(vec![2]), &t).is_err());
        assert!(plan_cpus(&CpuRequest::Explicit((0..8).collect()), &t).is_err());
        assert!(plan_cpus(&CpuRequest::Explicit(vec![2, 6]), &t).is_ok());
        assert_eq!(plan_cpus(&CpuRequest::None, &t).unwrap(), None);
    }
}
