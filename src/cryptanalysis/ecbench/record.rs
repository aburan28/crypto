//! Records: one per execution, failures included.
//!
//! A record carries everything needed to audit its own figure: the
//! method and workload identities, the verified answer, the cost in the
//! unit with every phase, what was not priced, the process's resource
//! use, and the isolation level its own observations earned with the
//! reason for every level it did not.

use std::collections::BTreeMap;

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ecbench::host::{is_virtual, HostCapsule};
use crate::cryptanalysis::ecbench::isolation::{
    anon_pages_by_node, current_cpu, mempolicy, mems_allowed, own_affinity, Conditions, CpuPlan,
    EvictionSummary, HwCounters, HwCounts, Preflight, SelfPlacement, Thresholds,
};
use crate::cryptanalysis::ecbench::methods::{
    resolve, FactorBaseFacts, MethodSpec, OnlineWindow, PhaseRecord, ResolvedMethod, SolveReport,
    SolverStats,
};
use crate::cryptanalysis::ecbench::spec::Level;
use crate::cryptanalysis::ecbench::workload::{CurveSpec, TargetKind, Workload};
use crate::cryptanalysis::ic_boundary::FieldOps;
use crate::cryptanalysis::ic_measurement;

pub const RECORD_SCHEMA: &str = "ecbench.record/v1";

/// The rules [`grade`] applies.  Version 2 added the unreserved-sibling
/// check, the observed siblings, the policy-and-pages NUMA read-back and
/// the session's own tick rate.  An audit regrades a session only when it
/// was graded under the current version.
pub const GRADING_VERSION: u32 = 2;
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
    /// Planted or public; older inputs are planted.
    #[serde(default, skip_serializing_if = "TargetKind::is_planted")]
    pub target_kind: TargetKind,
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
    /// Hardware counters across the solve alone.
    #[serde(default)]
    pub hw_solve: Option<HwCounts>,
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
        hw_solve: None,
    };
    let job: ChildInput = match serde_json::from_str(input) {
        Ok(j) => j,
        Err(e) => {
            out.error = Some(format!("child input: {e}"));
            return out;
        }
    };
    let result = (|| -> Result<(SolveReport, Option<SchedStat>, HwCounts), String> {
        if job.schema != CHILD_INPUT_SCHEMA {
            return Err(format!("child input schema `{}`", job.schema));
        }
        let (w, inst) = Workload::build(
            &job.curve,
            job.target_seed,
            job.target_index,
            job.target_kind,
        )?;
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
        let counters = HwCounters::open();
        let before = schedstat();
        // Profiling is opt-in and the environment lookup is outside both
        // boundaries. Callgrind's first part is workload/setup prework;
        // sum every part after that marker through the solve marker, since
        // ic_measurement may dump additional internal phase parts.
        let profile_solve = std::env::var("ECBENCH_CALLGRIND_SOLVE").as_deref() == Ok("1");
        counters.start();
        if profile_solve {
            ic_measurement::callgrind_dump(b"ecbench_before_solve\0");
        }
        let report = crate::cryptanalysis::ecbench::methods::solve(
            &m,
            &inst,
            &w.curve,
            &w.target,
            job.algorithm_seed,
        );
        if profile_solve {
            ic_measurement::callgrind_dump(b"ecbench_solve\0");
        }
        let report = report?;
        let hw = counters.stop();
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
        Ok((report, delta, hw))
    })();
    placement.cpu_at_end = current_cpu();
    // Where the solve's memory landed, read after it ran.
    placement.anon_pages_by_node = anon_pages_by_node();
    out.placement = placement;
    match result {
        Ok((r, s, hw)) => {
            out.report = Some(r);
            out.schedstat_solve = s;
            out.hw_solve = Some(hw);
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
    /// User-space instructions and cycles across the solve, where the
    /// host has a hardware PMU.
    #[serde(default)]
    pub hw_solve: Option<HwCounts>,
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
    /// The one-target online window (AGENTS.md "IC measurements").
    #[serde(default)]
    pub online: Option<OnlineWindow>,
    /// The decomposition solver's statistics for an algebraic or SAT
    /// index-calculus run; `null` otherwise.  Informational, outside the
    /// replay comparison.
    #[serde(default)]
    pub solver: Option<SolverStats>,
    /// The field operations behind the group operations
    /// (`SolveReport::field_ops`): modular multiplications, squarings and
    /// inversions over the whole solve, when the method's group counted
    /// them.  Absent means unknown, never zero.  Written only when
    /// counted, so a record written before the block existed is byte for
    /// byte what it was and its seal still checks; the audit compares it
    /// on replay only when both the record and the replay carry it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub field_ops: Option<FieldOps>,
    pub time: Timing,
    pub isolation: IsolationRecord,
    pub env_class_id: String,
    pub binary_sha256: Option<String>,
}

impl Record {
    /// Seal: compute `record_id` from the exact line the record is
    /// written as, with the id empty, and return that line with the id.
    ///
    /// The seal covers bytes, not a parsed value: `ecbench.record/v1` was
    /// defined with serde_json's former, not-always-correctly-rounded float
    /// parser, so a record parsed and re-serialised can differ in a last digit
    /// from what was written. [`check_line_seal`] therefore checks the raw line
    /// on disk.
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
    legacy_json_float(x).unwrap_or(f64::NAN)
}

/// Reproduce the decimal parser used when `ecbench.record/v1` was defined.
///
/// serde_json formerly accumulated the significand as an integer and then
/// multiplied or divided by a power of ten. Enabling its `float_roundtrip`
/// feature repository-wide makes parsing correctly rounded, but sealed v1
/// records and all bounds derived from them retain the former result. Keep
/// that protocol behavior explicit instead of depending on a crate feature.
pub(crate) fn legacy_json_float(x: f64) -> Option<f64> {
    const POW10: [f64; 23] = [
        1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16,
        1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
    ];

    let encoded = serde_json::to_string(&x).ok()?;
    let (negative, unsigned) = encoded
        .strip_prefix('-')
        .map_or((false, encoded.as_str()), |rest| (true, rest));
    let (mantissa, explicit_exponent) =
        unsigned
            .split_once(['e', 'E'])
            .map_or((unsigned, 0), |(mantissa, exponent)| {
                exponent
                    .parse::<i32>()
                    .map(|value| (mantissa, value))
                    .unwrap_or((mantissa, i32::MAX))
            });
    if explicit_exponent == i32::MAX {
        return None;
    }

    let mut significand = 0u64;
    let mut fractional_digits = 0i32;
    let mut after_decimal = false;
    for byte in mantissa.bytes() {
        match byte {
            b'.' if !after_decimal => after_decimal = true,
            b'0'..=b'9' => {
                significand = significand
                    .checked_mul(10)?
                    .checked_add(u64::from(byte - b'0'))?;
                if after_decimal {
                    fractional_digits = fractional_digits.checked_add(1)?;
                }
            }
            _ => return None,
        }
    }
    let exponent = explicit_exponent.checked_sub(fractional_digits)?;
    let power = usize::try_from(exponent.unsigned_abs()).ok()?;
    let pow = *POW10.get(power)?;
    let mut value = significand as f64;
    if exponent >= 0 {
        value *= pow;
    } else {
        value /= pow;
    }
    Some(if negative { -value } else { value })
}

fn normalize_legacy_json_floats(value: &mut Value) -> Result<(), String> {
    match value {
        Value::Array(values) => {
            for value in values {
                normalize_legacy_json_floats(value)?;
            }
        }
        Value::Object(values) => {
            for value in values.values_mut() {
                normalize_legacy_json_floats(value)?;
            }
        }
        Value::Number(number) if number.is_f64() => {
            let current = number
                .as_f64()
                .ok_or_else(|| format!("record float {number} does not fit f64"))?;
            let legacy = legacy_json_float(current).ok_or_else(|| {
                format!("record float {number} is outside the v1 decimal-parser domain")
            })?;
            let number = serde_json::Number::from_f64(legacy)
                .ok_or_else(|| format!("legacy parse of record float {number} is not finite"))?;
            *value = Value::Number(number);
        }
        _ => {}
    }
    Ok(())
}

/// Deserialize a sealed `ecbench.record/v1` JSON value with the decimal
/// semantics under which the schema and its derived bounds were frozen.
pub(crate) fn from_str_legacy_floats<T: DeserializeOwned>(text: &str) -> Result<T, String> {
    let mut value: Value = serde_json::from_str(text).map_err(|error| error.to_string())?;
    normalize_legacy_json_floats(&mut value)?;
    serde_json::from_value(value).map_err(|error| error.to_string())
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
    pub eviction: Option<&'a EvictionSummary>,
    pub runner_on_run_cpu: bool,
    pub capsule: &'a HostCapsule,
    pub thresholds: &'a Thresholds,
    pub placement: &'a SelfPlacement,
    pub before: &'a Conditions,
    pub after: &'a Conditions,
    pub timing: &'a Timing,
    /// Clock ticks per second on the host that recorded the ticks.
    pub hz: f64,
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
    let hz = g.hz;
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
            if !plan.unreserved_siblings.is_empty() {
                l1.push(format!(
                    "the run CPU's SMT siblings {:?} are outside the reservation",
                    plan.unreserved_siblings
                ));
            }
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
            // Every other observed CPU: the reserved idle siblings, and under
            // an inherited placement the run CPU's siblings that are not ours.
            let mut sb = 0u64;
            for c in plan.observed_cpus().into_iter().filter(|&c| c != run) {
                if let (Some(a), Some(b)) = (tick(g.before, c), tick(g.after, c)) {
                    sb += b.busy.saturating_sub(a.busy);
                }
            }
            sib_busy = Some(sb);
            if sb > 2 {
                l2.push(format!("{sb} busy ticks on the run CPU's idle siblings and the rest of the reservation"));
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

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Deserialize)]
    struct LegacyFloatFixture {
        gae: f64,
        nested: Vec<f64>,
        integer: u64,
    }

    #[test]
    fn legacy_record_deserializer_reproduces_v1_decimal_semantics() {
        let text = r#"{"gae":3934.0137668849675,"nested":[0.37222665730402227],"integer":17}"#;
        let current: LegacyFloatFixture = serde_json::from_str(text).unwrap();
        let legacy: LegacyFloatFixture = from_str_legacy_floats(text).unwrap();

        assert_eq!(
            legacy.gae.to_bits(),
            legacy_json_float(current.gae).unwrap().to_bits()
        );
        assert_eq!(
            legacy.nested[0].to_bits(),
            legacy_json_float(current.nested[0]).unwrap().to_bits()
        );
        // A newer serde_json parser may already yield the historical value.
        // Compatibility is the archived bit pattern, independent of whether
        // the current parser happens to differ on this decimal fixture.
        assert_eq!(legacy.gae.to_bits(), 3934.013766884967_f64.to_bits());
        assert_eq!(legacy.nested[0].to_bits(), 0.3722266573040223_f64.to_bits());
        assert_eq!(legacy.integer, 17);
    }

    #[test]
    fn legacy_record_roundtrip_is_not_a_numeric_tolerance() {
        let recorded: f64 = 3934.013766884967;
        let replayed: f64 = 3934.0137668849675;
        assert_ne!(recorded.to_bits(), replayed.to_bits());
        assert_eq!(
            json_roundtrip(replayed).to_bits(),
            json_roundtrip(recorded).to_bits()
        );

        let two_ulp_mutation = f64::from_bits(recorded.to_bits() + 2);
        assert_ne!(
            json_roundtrip(two_ulp_mutation).to_bits(),
            json_roundtrip(recorded).to_bits()
        );
    }
}
