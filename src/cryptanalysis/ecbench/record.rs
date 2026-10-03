//! Records: one per execution, failures included.
//!
//! A record carries everything needed to audit its own figure: the
//! method and workload identities, the verified answer, the cost in the
//! unit with every phase, what was not priced, the process's resource
//! use, and the isolation level its own observations earned with the
//! reason for every level it did not.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ecbench::host::{is_virtual, HostCapsule};
use crate::cryptanalysis::ecbench::isolation::{
    anon_pages_by_node, current_cpu, mempolicy, mems_allowed, own_affinity, Conditions, CpuPlan,
    Eviction, Preflight, SelfPlacement, Thresholds,
};
use crate::cryptanalysis::ecbench::methods::{
    resolve, FactorBaseFacts, MethodSpec, PhaseRecord, ResolvedMethod, SolveReport,
};
use crate::cryptanalysis::ecbench::spec::Level;
use crate::cryptanalysis::ecbench::workload::{CurveSpec, Workload};

pub const RECORD_SCHEMA: &str = "ecbench.record/v1";
pub const CHILD_INPUT_SCHEMA: &str = "ecbench.child_input/v1";
pub const CHILD_OUTPUT_SCHEMA: &str = "ecbench.child_output/v1";

/// The unit every cost is in: group-addition equivalents (additions plus
/// doublings, plus native work at the repository's pinned ratios where
/// one exists), divided by `√r` for `S`.
pub const UNIT: &str = "ecbench.gae";

// ── The child protocol ─────────────────────────────────────────────

/// What the runner sends a measured child on stdin.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ChildInput {
    pub schema: String,
    pub curve: CurveSpec,
    pub target_seed: u64,
    pub target_index: u64,
    pub expected_workload_id: String,
    pub method: MethodSpec,
    pub expected_method_id: String,
    pub algorithm_seed: u64,
}

/// `/proc/self/schedstat`: on-CPU time, run-queue wait, timeslices.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SchedStat {
    pub cpu_ns: u64,
    pub run_delay_ns: u64,
    pub timeslices: u64,
}

pub fn schedstat() -> Option<SchedStat> {
    #[cfg(target_os = "linux")]
    {
        let s = std::fs::read_to_string("/proc/self/schedstat").ok()?;
        let v: Vec<u64> = s
            .split_whitespace()
            .filter_map(|x| x.parse().ok())
            .collect();
        if v.len() >= 3 {
            return Some(SchedStat {
                cpu_ns: v[0],
                run_delay_ns: v[1],
                timeslices: v[2],
            });
        }
    }
    None
}

/// What a measured child prints on stdout.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ChildOutput {
    pub schema: String,
    pub workload_id: Option<String>,
    pub method_id: Option<String>,
    pub report: Option<SolveReport>,
    pub error: Option<String>,
    pub placement: SelfPlacement,
    /// Scheduler statistics across the solve alone.
    pub schedstat_solve: Option<SchedStat>,
}

/// The measured child's whole job: rebuild the workload, check its
/// identity, solve, report.  It never sees the planted scalar's use: the
/// runner verifies the answer in its own process.
pub fn child_main(input: &str) -> ChildOutput {
    let mut placement = SelfPlacement {
        affinity: own_affinity(),
        cpu_at_start: current_cpu(),
        cpu_at_end: None,
        mems_allowed: mems_allowed(),
        mempolicy: mempolicy(),
        anon_pages_by_node: None,
    };
    let mut out = ChildOutput {
        schema: CHILD_OUTPUT_SCHEMA.into(),
        workload_id: None,
        method_id: None,
        report: None,
        error: None,
        placement: placement.clone(),
        schedstat_solve: None,
    };
    let job: ChildInput = match serde_json::from_str(input) {
        Ok(j) => j,
        Err(e) => {
            out.error = Some(format!("child input: {e}"));
            return out;
        }
    };
    let result = (|| -> Result<(SolveReport, Option<SchedStat>), String> {
        if job.schema != CHILD_INPUT_SCHEMA {
            return Err(format!("child input schema `{}`", job.schema));
        }
        let (w, inst) = Workload::build(&job.curve, job.target_seed, job.target_index)?;
        if w.workload_id != job.expected_workload_id {
            return Err(format!(
                "rebuilt workload {} is not the planned {}",
                w.workload_id, job.expected_workload_id
            ));
        }
        let m = resolve(&job.method)?;
        if m.method_id != job.expected_method_id {
            return Err(format!(
                "method resolved to {}, planned {}",
                m.method_id, job.expected_method_id
            ));
        }
        let before = schedstat();
        let report = crate::cryptanalysis::ecbench::methods::solve(
            &m,
            &inst,
            &w.curve,
            &w.target,
            job.algorithm_seed,
        )?;
        let after = schedstat();
        let delta = match (before, after) {
            (Some(a), Some(b)) => Some(SchedStat {
                cpu_ns: b.cpu_ns.saturating_sub(a.cpu_ns),
                run_delay_ns: b.run_delay_ns.saturating_sub(a.run_delay_ns),
                timeslices: b.timeslices.saturating_sub(a.timeslices),
            }),
            _ => None,
        };
        out.workload_id = Some(w.workload_id);
        out.method_id = Some(m.method_id);
        Ok((report, delta))
    })();
    placement.cpu_at_end = current_cpu();
    // Where the solve's memory landed, read after it ran.
    placement.anon_pages_by_node = anon_pages_by_node();
    out.placement = placement;
    match result {
        Ok((r, s)) => {
            out.report = Some(r);
            out.schedstat_solve = s;
        }
        Err(e) => out.error = Some(e),
    }
    out
}

// ── The record ─────────────────────────────────────────────────────

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Outcome {
    /// `verified`, `wrong_answer`, `exhausted`, `error`, `timeout` or
    /// `crashed`.  Only `verified` is a result.
    pub status: String,
    pub recovered: Option<String>,
    /// `[recovered]G = Q`, checked by the runner in its own process.
    pub matches_target: Option<bool>,
    pub matches_planted: Option<bool>,
    pub error: Option<String>,
    /// Exit code or signal of the measured process.
    pub exit: Option<String>,
    /// SHA-256 of the child's stderr, kept as `exec/<seq>.stderr` when
    /// nonempty.
    pub stderr_sha256: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Cost {
    pub unit: String,
    pub total_gae: Option<f64>,
    /// `total_gae / √r`: the column every row in this repository is read
    /// in, rho at about 1.25 for `A = 1`.
    pub s: Option<f64>,
    /// Some work was counted and not priced; the total is a floor.
    pub lower_bound: bool,
    pub unpriced: Vec<String>,
    pub deterministic: bool,
    pub nondeterminism: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Boundaries {
    /// Automorphisms a generic algorithm may use on this curve.
    pub automorphisms_available: u32,
    pub automorphisms_used: Option<u32>,
    /// `√(π / 2A)` with `A` available: the generic collision floor in `S`
    /// for this curve, which no generic method beats on average.
    pub floor_s: f64,
    pub ratio_to_floor: Option<f64>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Timing {
    /// In-process wall time of the algorithm alone.
    pub solve_wall_ns: Option<u64>,
    /// Fork to reap, as the runner timed it: a practicality figure.
    pub process_wall_ns: u64,
    pub user_ns: u64,
    pub sys_ns: u64,
    pub max_rss_kib: u64,
    pub minor_faults: u64,
    pub major_faults: u64,
    pub voluntary_switches: u64,
    pub involuntary_switches: u64,
    pub schedstat_solve: Option<SchedStat>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IsolationRecord {
    pub level: Level,
    /// Why the run did not earn each higher level.  Empty at L3.
    pub blockers: Vec<String>,
    pub run_cpu: Option<u32>,
    pub reserved: Vec<u32>,
    pub node: Option<u32>,
    pub placement: SelfPlacement,
    /// Busy time on the run CPU not accounted for by the child, ticks.
    pub foreign_busy_ticks: Option<i64>,
    pub steal_ticks: Option<u64>,
    pub idle_sibling_busy_ticks: Option<u64>,
    pub before: Conditions,
    pub after: Conditions,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Record {
    pub schema: String,
    /// `ECR1h` + 12 hex of this record with `record_id` empty.
    pub record_id: String,
    pub session_id: String,
    /// `<method_id><workload_id>R<n>`: the repository's
    /// `{candidate}W{workload}R{number}` convention.  `n` counts, from 1,
    /// every execution of this (method, workload) in the session in
    /// sequence order — warm-ups and A/A arms included — so a run id is
    /// unique within its session.
    pub run_id: String,
    pub seq: u64,
    pub round: u32,
    pub warmup: bool,
    pub arm: String,
    pub role: String,
    pub method: ResolvedMethod,
    pub workload: Workload,
    pub algorithm_seed: String,
    pub outcome: Outcome,
    pub cost: Cost,
    pub boundaries: Boundaries,
    pub phases: Vec<PhaseRecord>,
    pub counters: BTreeMap<String, u64>,
    pub factor_base: Option<FactorBaseFacts>,
    pub detail: Value,
    pub time: Timing,
    pub isolation: IsolationRecord,
    pub env_class_id: String,
    pub binary_sha256: Option<String>,
}

impl Record {
    /// Seal: compute `record_id` from the exact line the record is
    /// written as, with the id empty, and return that line with the id.
    ///
    /// The seal covers bytes, not a parsed value: `serde_json`'s default
    /// float parser is not correctly rounded, so a record parsed and
    /// re-serialised can differ in a last digit from what was written.
    /// [`check_line_seal`] therefore checks the raw line on disk.
    pub fn seal(&mut self) -> String {
        self.record_id = String::new();
        let body = serde_json::to_string(self).expect("a record serialises");
        self.record_id = format!("ECR1h{}", &sha256_hex(body.as_bytes())[..12]);
        serde_json::to_string(self).expect("a record serialises")
    }

    /// A result: verified, measured (not a warm-up).
    pub fn counts(&self) -> bool {
        !self.warmup && self.outcome.status == "verified"
    }
}

/// Whether a raw `records.jsonl` line hashes to the id it carries: the
/// line with its `"record_id":"…"` emptied is what [`Record::seal`]
/// hashed.
pub fn check_line_seal(line: &str, record_id: &str) -> bool {
    let field = format!("\"record_id\":\"{record_id}\"");
    if !line.contains(&field) || !record_id.starts_with("ECR1h") {
        return false;
    }
    let blank = line.replacen(&field, "\"record_id\":\"\"", 1);
    format!("ECR1h{}", &sha256_hex(blank.as_bytes())[..12]) == record_id
}

/// `x` after the serialise-parse round trip a record's floats take, so a
/// recomputed value compares with a parsed one bit for bit.
pub fn json_roundtrip(x: f64) -> f64 {
    serde_json::to_string(&x)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or(f64::NAN)
}

/// `√(π / 2A)`.
pub fn floor_s(automorphisms: u32) -> f64 {
    (std::f64::consts::PI / (2.0 * automorphisms.max(1) as f64)).sqrt()
}

// ── Grading ────────────────────────────────────────────────────────

/// Everything grading reads.  Each check is a blocker for the level it
/// belongs to; a run earns the highest level with no blocker at it or
/// below.
pub struct GradeInput<'a> {
    pub plan: Option<&'a CpuPlan>,
    pub preflight: Option<&'a Preflight>,
    pub eviction: Option<&'a Eviction>,
    pub runner_on_run_cpu: bool,
    pub capsule: &'a HostCapsule,
    pub thresholds: &'a Thresholds,
    pub placement: &'a SelfPlacement,
    pub before: &'a Conditions,
    pub after: &'a Conditions,
    pub timing: &'a Timing,
}

pub struct Grade {
    pub level: Level,
    pub blockers: Vec<String>,
    pub foreign_busy_ticks: Option<i64>,
    pub steal_ticks: Option<u64>,
    pub idle_sibling_busy_ticks: Option<u64>,
}

pub fn grade(g: &GradeInput) -> Grade {
    let mut l1: Vec<String> = Vec::new();
    let mut l2: Vec<String> = Vec::new();
    let mut l3: Vec<String> = Vec::new();
    let mut foreign = None;
    let mut steal = None;
    let mut sib_busy = None;
    let hz = crate::cryptanalysis::ecbench::isolation::ticks_per_second() as f64;
    let wall_s = g.timing.process_wall_ns as f64 / 1e9;
    let cpu_s = (g.timing.user_ns + g.timing.sys_ns) as f64 / 1e9;

    if g.capsule.build.debug_assertions {
        l1.push("debug build: its wall time measures the build, not the method".into());
    }
    match g.plan {
        None => {
            l1.push("no CPU reservation (--cpus none, or no affinity control on this OS)".into())
        }
        Some(plan) => {
            let run = plan.run_cpu;
            if g.placement.affinity.as_deref() != Some(&[run][..]) {
                l1.push(format!(
                    "affinity read back inside the child was {:?}, not [{run}]",
                    g.placement.affinity
                ));
            }
            if g.placement.cpu_at_start != Some(run) || g.placement.cpu_at_end != Some(run) {
                l1.push(format!(
                    "child ran on CPU {:?} → {:?}, not {run}",
                    g.placement.cpu_at_start, g.placement.cpu_at_end
                ));
            }
            let core0 = g
                .capsule
                .stable
                .topology
                .place(0)
                .map(|p| p.siblings.clone())
                .unwrap_or_else(|| vec![0]);
            if core0.contains(&run) {
                l1.push("the run CPU shares CPU 0's core, where housekeeping lands".into());
            }
            if g.runner_on_run_cpu {
                l1.push("the runner itself had no CPU outside the reservation".into());
            }
            match g.eviction {
                Some(ev) if ev.failed_user_threads > 0 => l1.push(format!(
                    "{} user threads could not be moved off the reserved CPUs (not root?)",
                    ev.failed_user_threads
                )),
                None => l1.push("no eviction was performed".into()),
                _ => {}
            }
            if wall_s > 0.0 && cpu_s / wall_s > 1.05 {
                l1.push(format!(
                    "CPU time {cpu_s:.3}s over wall {wall_s:.3}s: more than one thread ran"
                ));
            }
            if g.capsule.stable.numa_nodes > 1 {
                if let Some(node) = plan.node {
                    // The policy read back inside the child, and where its
                    // anonymous pages actually are.  `Mems_allowed` is the
                    // cpuset's permission, not the policy, and is not used.
                    let want = format!("bind:{node}");
                    if g.placement.mempolicy.as_deref() != Some(&want[..]) {
                        l1.push(format!(
                            "memory policy read back as {:?}, not {want}",
                            g.placement.mempolicy
                        ));
                    }
                    match &g.placement.anon_pages_by_node {
                        Some(pages) => {
                            let total: u64 = pages.values().sum();
                            let local = pages.get(&node).copied().unwrap_or(0);
                            if total > 0 && (local as f64) < 0.99 * total as f64 {
                                l1.push(format!(
                                    "{} of {total} anonymous pages on node {node}",
                                    local
                                ));
                            }
                        }
                        None => l1.push("no /proc/self/numa_maps to confirm page placement".into()),
                    }
                }
            }
            // L2: the window was quiet.
            match g.preflight {
                Some(p) if p.quiet => {}
                Some(p) => l2.push(format!("preflight not quiet: {}", p.reasons.join("; "))),
                None => l2.push("no preflight".into()),
            }
            let tick = |c: &Conditions, cpu: u32| c.reserved_ticks.get(&cpu).copied();
            match (tick(g.before, run), tick(g.after, run)) {
                (Some(a), Some(b)) => {
                    let busy = b.busy.saturating_sub(a.busy) as f64;
                    let child_ticks = cpu_s * hz;
                    let f = (busy - child_ticks).round() as i64;
                    foreign = Some(f);
                    let allowance = (2.0f64).max(0.01 * wall_s * hz);
                    if f as f64 > allowance {
                        l2.push(format!(
                            "{f} ticks of foreign work on the run CPU during the run"
                        ));
                    }
                    let st = b.steal.saturating_sub(a.steal);
                    steal = Some(st);
                    if st > 0 {
                        l2.push(format!("{st} ticks of hypervisor steal on the run CPU"));
                    }
                }
                _ => l2.push("no /proc/stat ticks for the run CPU".into()),
            }
            let mut sb = 0u64;
            for &c in plan.reserved.iter().filter(|&&c| c != run) {
                if let (Some(a), Some(b)) = (tick(g.before, c), tick(g.after, c)) {
                    sb += b.busy.saturating_sub(a.busy);
                }
            }
            sib_busy = Some(sb);
            if sb > 2 {
                l2.push(format!("{sb} busy ticks on the reserved idle siblings"));
            }
            match (g.timing.schedstat_solve, g.timing.solve_wall_ns) {
                (Some(s), Some(w)) if w > 0 => {
                    if s.run_delay_ns as f64 > 0.005 * w as f64 {
                        l2.push(format!(
                            "run-queue delay {:.3}% of the solve",
                            100.0 * s.run_delay_ns as f64 / w as f64
                        ));
                    }
                }
                _ => l2.push("no schedstat for the solve".into()),
            }
            if wall_s > 0.0 {
                let rate = g.timing.involuntary_switches as f64 / wall_s;
                if rate > g.thresholds.max_preemptions_per_second {
                    l2.push(format!("{rate:.1} preemptions per second"));
                }
            }
            match (
                g.before.psi_memory_some_total_us,
                g.after.psi_memory_some_total_us,
            ) {
                (Some(a), Some(b)) => {
                    let stall = b.saturating_sub(a) as f64 / 1e6;
                    if stall > 0.005 * wall_s.max(0.001) {
                        l2.push(format!("{stall:.4}s of memory stall system-wide"));
                    }
                }
                _ => l2.push("no memory PSI".into()),
            }
            // L3: the host is configured for measurement.
            let s = g.capsule.stable.cpu_settings.get(&run);
            let isolated = s.map(|s| s.in_isolcpus == Some(true) || s.in_nohz_full == Some(true));
            if isolated != Some(true) {
                l3.push("run CPU is in neither isolcpus nor nohz_full".into());
            }
            if s.and_then(|s| s.governor.as_deref()) != Some("performance") {
                l3.push(format!(
                    "governor is {:?}, not performance",
                    s.and_then(|s| s.governor.clone())
                ));
            }
            let turbo_off =
                s.map(|s| s.no_turbo.as_deref() == Some("1") || s.boost.as_deref() == Some("0"));
            if turbo_off != Some(true) {
                l3.push("turbo/boost is not known to be off".into());
            }
            if is_virtual(&g.capsule.stable) != Some(false) {
                l3.push("not known to be bare metal".into());
            }
        }
    }
    let level = if !l1.is_empty() {
        Level::L0
    } else if !l2.is_empty() {
        Level::L1
    } else if !l3.is_empty() {
        Level::L2
    } else {
        Level::L3
    };
    let mut blockers = Vec::new();
    for (lv, list) in [("L1", l1), ("L2", l2), ("L3", l3)] {
        blockers.extend(list.into_iter().map(|b| format!("{lv}: {b}")));
    }
    Grade {
        level,
        blockers,
        foreign_busy_ticks: foreign,
        steal_ticks: steal,
        idle_sibling_busy_ticks: sib_busy,
    }
}
