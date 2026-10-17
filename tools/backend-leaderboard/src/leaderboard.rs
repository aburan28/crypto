//! Backend frontiers on exact single-target workloads, derived from sealed sessions.
//!
//! Sessions, targets and binaries are never pooled. Native counters retain their
//! units; the pricing band is a sensitivity analysis, not a calibrated total.
//! Failed attempts and interrupted plans remain visible and cannot win.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Component, Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::band;
use super::canonical::sha256_hex;
use super::evidence::{self as audit, AuditReport};
use super::evidence::{
    read_records, read_session, Cost, FactorBaseFacts, HostCapsule, Level, OnlineWindow, Outcome,
    PhaseRecord, Record, ResolvedMethod, SolverStats, Spec, Timing, Workload,
};
use super::stats::{bootstrap_ci, median};

pub const SCHEMA: &str = "ecbench.backend_leaderboard/v1";
pub const MANIFEST_SCHEMA: &str = "ecbench.backend_leaderboard_manifest/v1";

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema: String,
    pub description: String,
    /// Extension degrees, not subgroup sizes or numbers of polynomial variables.
    pub field_degrees: Vec<u32>,
    pub sessions: Vec<Input>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub dir: PathBuf,
    pub replay_receipt: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
pub struct Source {
    pub dir: String,
    pub session_id: String,
    pub status: String,
    pub git_commit: Option<String>,
    pub git_dirty: Option<bool>,
    pub binary_sha256: Option<String>,
    pub env_class_id: String,
    pub cpu: Option<String>,
    pub architecture: String,
    pub threads: BTreeMap<String, String>,
    pub required_level: Level,
    pub resource_envelope: serde_json::Value,
    pub executions_planned: u64,
    pub records_written: u64,
    pub integrity_ok: bool,
    pub integrity_problems: Vec<String>,
    pub file_hashes: BTreeMap<String, String>,
    pub receipt: Option<String>,
    pub receipt_sha256: Option<String>,
    pub receipt_bound: bool,
    pub receipt_problems: Vec<String>,
    /// Historical replay evidence, not a fresh replay by this report builder.
    pub reproduced_record_ids: BTreeSet<String>,
}

#[derive(Debug, Serialize)]
pub struct Sample {
    pub record_id: String,
    pub run_id: String,
    pub seq: u64,
    pub round: u32,
    pub algorithm_seed: String,
    pub outcome: Outcome,
    pub cost: Cost,
    pub recorded_s: Option<f64>,
    pub band_s: Vec<Option<f64>>,
    pub unpriced_counts: BTreeMap<String, Option<u64>>,
    pub phases: Vec<PhaseRecord>,
    pub counters: BTreeMap<String, u64>,
    pub solver: Option<SolverStats>,
    pub solver_origin: Option<String>,
    pub solver_context: serde_json::Value,
    pub fbm_id: Option<String>,
    pub field_ops: Option<serde_json::Value>,
    pub online: Option<OnlineWindow>,
    pub time: Timing,
    pub isolation: Level,
    pub isolation_blockers: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Entry {
    pub entry_id: String,
    pub source: usize,
    pub arm: String,
    pub role: String,
    pub backend: String,
    pub method: ResolvedMethod,
    pub workload: Workload,
    pub factor_base: Option<FactorBaseFacts>,
    pub attempted: usize,
    pub expected: u32,
    pub verified: usize,
    pub statuses: BTreeMap<String, usize>,
    pub lower_bound: bool,
    pub band_mean_s: Vec<Option<f64>>,
    pub cold_median_ms: Option<f64>,
    pub cold_range_ms: Option<(f64, f64)>,
    pub process_median_ms: Option<f64>,
    pub online_median_ms: Option<f64>,
    pub peak_rss_kib: u64,
    pub solver_median_ops: Option<f64>,
    pub solver_op_unit: Option<String>,
    pub sat_conflicts_median: Option<f64>,
    pub sat_learnt_clauses_median: Option<f64>,
    pub count_blockers: Vec<String>,
    pub full_historical_replay: bool,
    pub count_dominated_by: Vec<String>,
    /// The current IC recipe boundary, keeping rho/BSGS outside the incumbent
    /// set while retaining their separate cross-method comparison.
    pub ic_count_dominated_by: Vec<String>,
    pub ic_count_frontier: bool,
    /// A nondominated *observed* count-band result, not a ledger promotion.
    pub count_frontier: bool,
    /// Memory versus count-band tradeoffs, only for complete measured vectors.
    pub resource_frontier: bool,
    pub samples: Vec<Sample>,
}

#[derive(Debug, Serialize)]
pub struct WallComparison {
    pub workload_id: String,
    pub a: String,
    pub b: String,
    pub pairs: usize,
    pub median_b_over_a: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    pub aa_arms: Option<(String, String)>,
    pub aa_ci95: Option<(f64, f64)>,
    pub winner: Option<String>,
    pub blockers: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Coverage {
    pub field_degree: u32,
    pub field_size: String,
    pub note: String,
    /// Observations / count-eligible exact-target entries, by backend.
    pub backends: BTreeMap<String, (usize, usize)>,
    pub planned_without_measurements: BTreeMap<String, usize>,
}

#[derive(Debug, Serialize)]
pub struct Unattempted {
    pub source: usize,
    pub arm: String,
    pub backend: String,
    pub workload: Workload,
    pub expected_rounds: u32,
    pub reason: &'static str,
}

#[derive(Debug, Serialize)]
pub struct Leaderboard {
    pub schema: &'static str,
    pub description: String,
    pub manifest_sha256: String,
    pub producer_input_hashes: BTreeMap<String, String>,
    pub scope: &'static str,
    pub weights_gae_per_native_op: Vec<f64>,
    pub sources: Vec<Source>,
    pub coverage: Vec<Coverage>,
    pub entries: Vec<Entry>,
    pub unattempted: Vec<Unattempted>,
    pub wall_comparisons: Vec<WallComparison>,
}

fn relative(path: &Path) -> Result<String, String> {
    if path.as_os_str().is_empty()
        || path
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
    {
        return Err(format!(
            "expected a repository-relative path: {}",
            path.display()
        ));
    }
    Ok(path.to_string_lossy().replace('\\', "/"))
}

/// Categorise the configured engine, never an arm's user-chosen label.
pub fn backend(method: &ResolvedMethod) -> String {
    let oracle = method
        .params
        .get("oracle")
        .map(String::as_str)
        .unwrap_or("");
    let solver = method
        .params
        .get("solver")
        .map(String::as_str)
        .unwrap_or("");
    let engine = oracle
        .split([',', ':'])
        .find_map(|p| p.strip_prefix("engine="));
    match engine {
        Some("matrix-f4") => "F4 matrix",
        Some("inherited-f4") => "F4 inherited",
        Some("matrix-f5") => "F5 matrix",
        Some("f6-ic") => "F6-IC",
        Some(e) if e.starts_with("rootfind") => "Root finding",
        _ if oracle.contains("rootfind") || oracle.contains("subtract") => "Root finding",
        _ if solver.contains("sat") => "SAT",
        _ if solver.contains("f4") => "F4 matrix",
        _ if method.family == "rho" => "Pollard rho",
        _ if oracle.contains("mitm") || oracle.contains("pair") => "Pair / MITM",
        _ => {
            return if solver.is_empty() {
                method.id.clone()
            } else {
                solver.into()
            }
        }
    }
    .into()
}

fn all_mean(values: impl Iterator<Item = Option<f64>>) -> Option<f64> {
    let values: Option<Vec<f64>> = values.collect();
    let values = values?;
    (!values.is_empty() && values.iter().all(|v| v.is_finite() && *v >= 0.0))
        .then(|| values.iter().sum::<f64>() / values.len() as f64)
}

fn times(records: &[&Record], f: impl Fn(&Record) -> Option<u64>) -> Vec<f64> {
    records
        .iter()
        .filter_map(|r| f(r).map(|n| n as f64 / 1e6))
        .collect()
}

fn complete_online(o: &OnlineWindow, family: &str) -> bool {
    if o.start_event.is_empty()
        || o.stop_event.is_empty()
        || o.phases_ns
            .values()
            .try_fold(0u64, |a, b| a.checked_add(*b))
            != Some(o.wall_ns)
    {
        return false;
    }
    if family == "ic" {
        [
            "target_query",
            "target_PDP",
            "target_relation_check",
            "target_descent",
            "target_recovery_check",
        ]
        .iter()
        .all(|name| o.phases_ns.contains_key(*name))
    } else {
        !o.phases_ns.is_empty()
    }
}

fn receipt(root: &Path, input: &Input, source: &mut Source) -> Result<(), String> {
    let Some(path) = &input.replay_receipt else {
        return Ok(());
    };
    source.receipt = Some(relative(path)?);
    let bytes = std::fs::read(root.join(path)).map_err(|e| format!("{}: {e}", path.display()))?;
    source.receipt_sha256 = Some(sha256_hex(&bytes));
    let r: AuditReport = serde_json::from_slice(&bytes).map_err(|e| format!("receipt: {e}"))?;
    if r.schema != "ecbench.audit/v1"
        || r.session_id != source.session_id
        || !r.ok
        || r.records_bound == Some(false)
    {
        source
            .receipt_problems
            .push("receipt schema, session or audit result does not match".into());
    }
    for name in [
        "session.json",
        "host.json",
        "spec.json",
        "plan.json",
        "records.jsonl",
    ] {
        if !r.files.contains_key(name) || r.files.get(name) != source.file_hashes.get(name) {
            source
                .receipt_problems
                .push(format!("receipt does not bind {name}"));
        }
    }
    source.receipt_bound = source.receipt_problems.is_empty();
    if source.receipt_bound {
        source.reproduced_record_ids = r
            .replays
            .into_iter()
            .filter(|r| r.reproduced)
            .map(|r| r.record_id)
            .collect();
    }
    Ok(())
}

fn summarize(source: usize, s: &Source, records: &[&Record], expected: u32) -> Entry {
    let first = records[0];
    let mut blockers = BTreeSet::new();
    if !s.integrity_ok {
        blockers.insert("session integrity audit failed".to_string());
    }
    if !s.receipt_bound {
        blockers.insert("no hash-bound historical native arithmetic audit".to_string());
    }
    if s.status != "complete" {
        blockers.insert(format!("session is {}", s.status));
    }
    if records.len() != expected as usize {
        blockers.insert("measured rounds are incomplete".into());
    }
    let verified = records
        .iter()
        .filter(|r| r.counts() && r.outcome.matches_target == Some(true))
        .count();
    if verified != records.len() {
        blockers.insert("failed or unverified attempts".into());
    }
    if records.iter().any(|r| !r.cost.deterministic) {
        blockers.insert("nondeterministic charged work".into());
    }
    let mut statuses = BTreeMap::new();
    let samples: Vec<Sample> = records
        .iter()
        .map(|r| {
            *statuses.entry(r.outcome.status.clone()).or_default() += 1;
            let unpriced_counts = band::unpriced_counts(r);
            for (name, count) in &unpriced_counts {
                if count.is_none() {
                    blockers.insert(format!("missing native counter: {name}"));
                }
            }
            let band_s = band::DEFAULT_WEIGHTS
                .iter()
                .map(|&w| match band::priced_s(r, w) {
                    Ok(value) if value.is_some_and(|v| v.is_finite() && v >= 0.0) => value,
                    Ok(_) => {
                        blockers.insert("missing or invalid priced count".into());
                        None
                    }
                    Err(e) => {
                        blockers.insert(e);
                        None
                    }
                })
                .collect();
            Sample {
                record_id: r.record_id.clone(),
                run_id: r.run_id.clone(),
                seq: r.seq,
                round: r.round,
                algorithm_seed: r.algorithm_seed.clone(),
                outcome: r.outcome.clone(),
                cost: r.cost.clone(),
                recorded_s: r.cost.s,
                band_s,
                unpriced_counts,
                phases: r.phases.clone(),
                counters: r.counters.clone(),
                solver: r.solver.clone(),
                solver_origin: r.solver_origin.clone(),
                solver_context: r
                    .detail
                    .as_object()
                    .map(|d| {
                        serde_json::Value::Object(
                            d.iter()
                                .filter(|(k, _)| {
                                    [
                                        "decomposition",
                                        "linear_algebra",
                                        "label",
                                        "calibration_pinned_units",
                                        "calibration_fallback_units",
                                        "framework_total_gae_before_solver_removal",
                                    ]
                                    .contains(&k.as_str())
                                })
                                .map(|(k, v)| (k.clone(), v.clone()))
                                .collect(),
                        )
                    })
                    .unwrap_or(serde_json::Value::Null),
                fbm_id: r.fbm_id.clone(),
                field_ops: r.field_ops.clone(),
                online: r.online.clone(),
                time: r.time.clone(),
                isolation: r.isolation.level,
                isolation_blockers: r.isolation.blockers.clone(),
            }
        })
        .collect();
    let band_mean_s = (0..band::DEFAULT_WEIGHTS.len())
        .map(|i| all_mean(samples.iter().map(|r| r.band_s[i])))
        .collect();
    let cold = times(records, |r| r.time.solve_wall_ns);
    let solver_units: BTreeSet<_> = records
        .iter()
        .filter_map(|r| r.solver.as_ref().map(|s| s.op_unit.clone()))
        .collect();
    let full_median = |values: Vec<f64>| {
        (values.len() == records.len())
            .then(|| median(&values))
            .flatten()
    };
    let solver_median_ops = full_median(
        records
            .iter()
            .filter_map(|r| r.solver.as_ref().map(|s| s.ops as f64))
            .collect(),
    );
    let sat_conflicts_median = full_median(
        records
            .iter()
            .filter_map(|r| {
                let s = r.solver.as_ref()?;
                s.rest
                    .get("sat_conflicts")
                    .and_then(serde_json::Value::as_f64)
                    .or_else(|| (s.op_unit == "conflicts").then_some(s.ops as f64))
            })
            .collect(),
    );
    let sat_learnt_clauses_median = full_median(
        records
            .iter()
            .filter_map(|r| {
                let s = r.solver.as_ref()?;
                s.rest
                    .get("sat_learnt_clauses")
                    .and_then(serde_json::Value::as_f64)
                    .or_else(|| s.rest.get("extra")?.get("learnt_clauses")?.as_f64())
            })
            .collect(),
    );
    let process = times(records, |r| Some(r.time.process_wall_ns));
    // Missing windows and incomplete clocks stay unknown; do not average a subset.
    let online = times(records, |r| {
        r.online
            .as_ref()
            .filter(|o| complete_online(o, &r.method.family))
            .map(|o| o.wall_ns)
    });
    Entry {
        entry_id: format!(
            "{}:{}:{}",
            s.session_id, first.workload.workload_id, first.arm
        ),
        source,
        arm: first.arm.clone(),
        role: first.role.clone(),
        backend: backend(&first.method),
        method: first.method.clone(),
        workload: first.workload.clone(),
        factor_base: first.factor_base.clone(),
        attempted: records.len(),
        expected,
        verified,
        statuses,
        lower_bound: records.iter().any(|r| r.cost.lower_bound),
        band_mean_s,
        cold_median_ms: (cold.len() == records.len())
            .then(|| median(&cold))
            .flatten(),
        cold_range_ms: (cold.len() == records.len() && !cold.is_empty()).then(|| {
            (
                *cold.iter().min_by(|a, b| a.total_cmp(b)).unwrap(),
                *cold.iter().max_by(|a, b| a.total_cmp(b)).unwrap(),
            )
        }),
        process_median_ms: median(&process),
        online_median_ms: (online.len() == records.len())
            .then(|| median(&online))
            .flatten(),
        peak_rss_kib: records
            .iter()
            .map(|r| r.time.max_rss_kib)
            .max()
            .unwrap_or(0),
        count_blockers: blockers.into_iter().collect(),
        solver_median_ops,
        solver_op_unit: (solver_units.len() == 1)
            .then(|| solver_units.into_iter().next())
            .flatten(),
        sat_conflicts_median,
        sat_learnt_clauses_median,
        full_historical_replay: s.receipt_bound
            && records
                .iter()
                .all(|r| s.reproduced_record_ids.contains(&r.record_id)),
        count_dominated_by: vec![],
        ic_count_dominated_by: vec![],
        ic_count_frontier: false,
        count_frontier: false,
        resource_frontier: false,
        samples,
    }
}

fn matched(a: &Entry, b: &Entry) -> bool {
    a.source == b.source
        && a.workload.workload_id == b.workload.workload_id
        && a.workload.workload_sha256 == b.workload.workload_sha256
        && a.samples
            .iter()
            .map(|r| (r.round, &r.algorithm_seed))
            .collect::<BTreeSet<_>>()
            == b.samples
                .iter()
                .map(|r| (r.round, &r.algorithm_seed))
                .collect::<BTreeSet<_>>()
}

fn dominates(a: &Entry, b: &Entry, memory: bool) -> bool {
    if !a.count_blockers.is_empty() || !b.count_blockers.is_empty() || !matched(a, b) {
        return false;
    }
    let costs: Option<Vec<(f64, f64)>> = a
        .band_mean_s
        .iter()
        .zip(&b.band_mean_s)
        .map(|(a, b)| Some(((*a)?, (*b)?)))
        .collect();
    let Some(costs) = costs else {
        return false;
    };
    if memory {
        a.peak_rss_kib > 0
            && b.peak_rss_kib > 0
            && a.peak_rss_kib <= b.peak_rss_kib
            && costs.iter().all(|(a, b)| a <= b)
            && (costs.iter().any(|(a, b)| a < b) || a.peak_rss_kib < b.peak_rss_kib)
    } else {
        costs.iter().all(|(a, b)| a < b)
    }
}

fn wall_ratios(a: &Entry, b: &Entry, required: Level) -> Vec<f64> {
    let by_round: BTreeMap<_, _> = b.samples.iter().map(|r| (r.round, r)).collect();
    a.samples
        .iter()
        .filter_map(|a| {
            let b = by_round.get(&a.round)?;
            if a.isolation < required
                || b.isolation < required
                || a.outcome.status != "verified"
                || b.outcome.status != "verified"
            {
                return None;
            }
            let (x, y) = (a.time.solve_wall_ns?, b.time.solve_wall_ns?);
            (x > 0 && y > 0).then_some(y as f64 / x as f64)
        })
        .collect()
}

fn interval(values: &[f64]) -> Option<(f64, f64)> {
    bootstrap_ci(&[values.to_vec()], 2000, 20261017, |s| median(&s[0]))
}

fn wall_comparison(a: &Entry, b: &Entry, entries: &[Entry], required: Level) -> WallComparison {
    let required = required.max(Level::L2);
    let mut blockers = Vec::new();
    let pairs = wall_ratios(a, b, required);
    let all = wall_ratios(a, b, Level::L0);
    if !matched(a, b) {
        blockers.push("rounds or algorithm seeds do not match".into());
    }
    for e in [a, b] {
        if e.attempted != e.expected as usize || e.verified != e.attempted {
            blockers.push(format!("{} has incomplete or failed attempts", e.arm));
        }
    }
    if pairs.len() != a.expected as usize || pairs.len() < super::evidence::MIN_WALL_PAIRS {
        blockers.push(format!(
            "{} paired rounds at {}; need every round and at least five",
            pairs.len(),
            required.name()
        ));
    }
    let peers: Vec<_> = entries.iter().filter(|e| matched(a, e)).collect();
    let mut aa = None;
    for (i, x) in peers.iter().enumerate() {
        for y in &peers[i + 1..] {
            if x.method.method_id == y.method.method_id
                && x.verified == x.expected as usize
                && y.verified == y.expected as usize
            {
                let values = wall_ratios(x, y, required);
                if values.len() >= super::evidence::MIN_WALL_PAIRS
                    && values.len() == x.expected as usize
                {
                    aa = interval(&values).map(|ci| ((x.arm.clone(), y.arm.clone()), ci));
                    break;
                }
            }
        }
        if aa.is_some() {
            break;
        }
    }
    if aa.is_none() {
        blockers.push("no complete same-target A/A noise interval at L2 or higher".into());
    }
    let ci95 = (pairs.len() >= super::evidence::MIN_WALL_PAIRS)
        .then(|| interval(&pairs))
        .flatten();
    let mut winner = None;
    if blockers.is_empty() {
        if let (Some((lo, hi)), Some((_, (alo, ahi)))) = (ci95, &aa) {
            if hi < 1.0 && hi < *alo {
                winner = Some(b.entry_id.clone());
            } else if lo > 1.0 && lo > *ahi {
                winner = Some(a.entry_id.clone());
            } else {
                blockers.push("paired 95% interval overlaps parity or A/A noise".into());
            }
        }
    }
    WallComparison {
        workload_id: a.workload.workload_id.clone(),
        a: a.entry_id.clone(),
        b: b.entry_id.clone(),
        pairs: pairs.len(),
        median_b_over_a: median(&all),
        ci95,
        aa_arms: aa.as_ref().map(|(arms, _)| arms.clone()),
        aa_ci95: aa.map(|(_, ci)| ci),
        winner,
        blockers,
    }
}

/// Build without running or retiming a solver. Every input is integrity-audited;
/// a failed audit is displayed but excludes the source from every frontier.
pub fn build(root: &Path, manifest_path: &Path) -> Result<Leaderboard, String> {
    let bytes =
        std::fs::read(manifest_path).map_err(|e| format!("{}: {e}", manifest_path.display()))?;
    let manifest: Manifest =
        serde_json::from_slice(&bytes).map_err(|e| format!("manifest: {e}"))?;
    if manifest.schema != MANIFEST_SCHEMA || manifest.sessions.is_empty() {
        return Err("expected a backend leaderboard manifest with sessions".into());
    }
    let mut sources = Vec::new();
    let mut entries = Vec::new();
    let mut wall_comparisons = Vec::new();
    let mut unattempted = Vec::new();
    let mut session_ids = BTreeSet::new();
    for input in &manifest.sessions {
        let dir_text = relative(&input.dir)?;
        let dir = root.join(&input.dir);
        let session = read_session(&dir)?;
        if !session_ids.insert(session.session_id.clone()) {
            return Err("duplicate session in manifest".into());
        }
        let spec_text =
            std::fs::read_to_string(dir.join("spec.json")).map_err(|e| e.to_string())?;
        let spec = Spec::from_json(&spec_text)?;
        let host: HostCapsule = serde_json::from_str(
            &std::fs::read_to_string(dir.join("host.json")).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        let audit = audit::audit_with(&dir, 0, true)?;
        let records = read_records(&dir)?;
        let mut source = Source {
            dir: dir_text,
            session_id: session.session_id,
            status: session.status,
            git_commit: session.git_commit,
            git_dirty: session.git_dirty,
            binary_sha256: session.binary_sha256,
            env_class_id: session.env_class_id,
            cpu: host.stable.cpu_model,
            architecture: host.stable.arch,
            threads: session
                .child_env
                .into_iter()
                .filter(|(k, _)| k.ends_with("NUM_THREADS"))
                .collect(),
            required_level: spec.measurement.isolation_required,
            resource_envelope: serde_json::json!({"options":session.options,"cpu_plan":session.cpu_plan}),
            executions_planned: session.executions_planned,
            records_written: session.records_written,
            integrity_ok: audit.ok && audit.records_bound != Some(false),
            integrity_problems: audit.problems,
            file_hashes: audit.files,
            receipt: None,
            receipt_sha256: None,
            receipt_bound: false,
            receipt_problems: vec![],
            reproduced_record_ids: BTreeSet::new(),
        };
        receipt(root, input, &mut source)?;
        let mut by = BTreeMap::<(String, String), Vec<&Record>>::new();
        for r in records.iter().filter(|r| !r.warmup) {
            by.entry((r.workload.workload_id.clone(), r.arm.clone()))
                .or_default()
                .push(r);
        }
        let plan: super::evidence::PlanDoc = serde_json::from_slice(
            &std::fs::read(dir.join("plan.json")).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        for w in &plan.workloads {
            for a in &plan.arms {
                if !by.contains_key(&(w.workload_id.clone(), a.name.clone())) {
                    unattempted.push(Unattempted {
                        source: sources.len(),
                        arm: a.name.clone(),
                        backend: backend(&a.method),
                        workload: w.clone(),
                        expected_rounds: spec.measurement.rounds,
                        reason: "planned arm has no measured records; warmups are not measurements",
                    });
                }
            }
        }
        let start = entries.len();
        for rs in by.values() {
            entries.push(summarize(
                sources.len(),
                &source,
                rs,
                spec.measurement.rounds,
            ));
        }
        for i in start..entries.len() {
            for j in i + 1..entries.len() {
                if entries[i].workload.workload_id != entries[j].workload.workload_id {
                    continue;
                }
                let mut c = wall_comparison(
                    &entries[i],
                    &entries[j],
                    &entries[start..],
                    spec.measurement.isolation_required,
                );
                if !source.integrity_ok || !source.receipt_bound || source.status != "complete" {
                    c.winner = None;
                    c.blockers.push("session is interrupted, failed integrity, or lacks a bound native arithmetic audit".into());
                }
                wall_comparisons.push(c);
            }
        }
        sources.push(source);
    }
    for i in 0..entries.len() {
        let dominated_by = entries
            .iter()
            .filter(|e| dominates(e, &entries[i], false))
            .map(|e| e.entry_id.clone())
            .collect::<Vec<_>>();
        let eligible = entries[i].count_blockers.is_empty();
        let resource_frontier = eligible
            && entries[i].peak_rss_kib > 0
            && !entries.iter().any(|e| dominates(e, &entries[i], true));
        entries[i].count_frontier = eligible && dominated_by.is_empty();
        entries[i].ic_count_dominated_by = entries
            .iter()
            .filter(|e| e.method.family == "ic" && dominates(e, &entries[i], false))
            .map(|e| e.entry_id.clone())
            .collect();
        entries[i].ic_count_frontier = eligible
            && entries[i].method.family == "ic"
            && entries[i].ic_count_dominated_by.is_empty();
        entries[i].resource_frontier = resource_frontier;
        entries[i].count_dominated_by = dominated_by;
    }
    let degrees: BTreeSet<u32> = manifest
        .field_degrees
        .into_iter()
        .chain(entries.iter().filter_map(|e| e.workload.curve.field_degree))
        .collect();
    let coverage = degrees.into_iter().map(|m| {
        let mut backends = BTreeMap::<String,(usize,usize)>::new();
        let mut planned_without_measurements = BTreeMap::<String,usize>::new();
        for e in entries.iter().filter(|e| e.workload.curve.field_degree == Some(m)) {
            let v = backends.entry(e.backend.clone()).or_default();
            v.0 += 1; v.1 += usize::from(e.count_blockers.is_empty());
        }
        for gap in unattempted.iter().filter(|e| e.workload.curve.field_degree == Some(m)) {
            *planned_without_measurements.entry(gap.backend.clone()).or_default() += 1;
        }
        Coverage { field_degree: m, field_size: format!("2^{m}"), backends, planned_without_measurements,
            note: match m {
                51 => "Composite degree: GF(2^3) and GF(2^17) are proper intermediate subfields.",
                83 => "Required smaller-curve scaling gate for ECC2K-130; exact E_0 model and subgroup must also match.",
                31 => "Exploratory degree; its Frobenius module differs from degrees 53, 83 and 131.",
                _ => "Coverage describes only the sessions in this manifest; absence is not an impossibility result.",
            }.into() }
    }).collect();
    let producer_input_hashes = [
        (
            "tools/backend-leaderboard/src/leaderboard.rs",
            include_str!("leaderboard.rs"),
        ),
        (
            "tools/backend-leaderboard/src/evidence.rs",
            include_str!("evidence.rs"),
        ),
        (
            "src/cryptanalysis/ecbench/band.rs",
            include_str!("../../../src/cryptanalysis/ecbench/band.rs"),
        ),
        (
            "src/cryptanalysis/ecbench/stats.rs",
            include_str!("../../../src/cryptanalysis/ecbench/stats.rs"),
        ),
        (
            "src/cryptanalysis/ecbench/canonical.rs",
            include_str!("../../../src/cryptanalysis/ecbench/canonical.rs"),
        ),
        (
            "src/cryptanalysis/ecbench/methods.rs",
            include_str!("../../../src/cryptanalysis/ecbench/methods.rs"),
        ),
        (
            "docs/ic/calibration.json",
            include_str!("../../../docs/ic/calibration.json"),
        ),
        (
            "src/cryptanalysis/curve_aliases.json",
            include_str!("../../../src/cryptanalysis/curve_aliases.json"),
        ),
    ]
    .into_iter()
    .map(|(path, text)| (path.into(), sha256_hex(text.as_bytes())))
    .collect();
    Ok(Leaderboard { schema: SCHEMA, description: manifest.description, manifest_sha256: sha256_hex(&bytes), producer_input_hashes,
        scope: "observed single-target, cold solve count-band frontiers within one session and exact workload; not ledger promotions, global records, or cross-host wall rankings",
        weights_gae_per_native_op: band::DEFAULT_WEIGHTS.to_vec(), sources, coverage, entries, unattempted, wall_comparisons })
}

fn number(v: Option<f64>) -> String {
    v.map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "unknown".into())
}
fn cell(s: &str) -> String {
    s.replace('|', "\\|").replace(['\n', '\r'], " ")
}

pub fn markdown(board: &Leaderboard) -> String {
    let mut out = format!(
        "# Backend and algorithm frontier\n\n{}\n\n{}\n\n",
        board.description, board.scope
    );
    out.push_str("Generated by the native `backend-leaderboard` tool. Each row is one exact target, backend configuration and session. Field degree m, field size 2^m and subgroup order r are separate. Unknown costs remain unknown.\n\n");
    out.push_str("Cold solve ms includes the method's setup, search, extraction and internal verification. Process ms is fork-to-reap and includes startup. The runner's independent scalar check is outside both clocks. Online ms excludes reusable setup. All clocks are descriptive unless a same-target L2-or-higher paired comparison passes its A/A and 95% interval gates.\n\n");
    out.push_str("## Coverage\n\n| degree m | field size | backend | observed / count eligible |\n|---:|---|---|---:|\n");
    for c in &board.coverage {
        if c.backends.is_empty() {
            out.push_str(&format!(
                "| {} | {} | no imported measurements | 0 / 0 |\n",
                c.field_degree, c.field_size
            ));
        }
        for (b, (n, k)) in &c.backends {
            out.push_str(&format!(
                "| {} | {} | {} | {} / {} |\n",
                c.field_degree,
                c.field_size,
                cell(b),
                n,
                k
            ));
        }
        for (b, n) in &c.planned_without_measurements {
            out.push_str(&format!(
                "| {} | {} | {}: {} planned targets without measurements | 0 / 0 |\n",
                c.field_degree,
                c.field_size,
                cell(b),
                n
            ));
        }
    }
    out.push_str("\n## Exact-target results\n\nS(w) charges all arms alike at w = 0, 0.1 and 1 GAE per native operation. These are sensitivity scenarios, not calibrated full costs. An observed count frontier has no strictly cheaper eligible arm at all three weights in its own cohort. Equal or crossing bands remain unresolved. The resource frontier also includes peak process RSS. Counts are finite-workload observations; one round provides no scaling inference.\n\n| session / target | backend / arm | verified / attempted / expected | S(0) / S(0.1) / S(1) | cold ms | online ms | peak KiB | status |\n|---|---|---:|---|---:|---:|---:|---|\n");
    for e in &board.entries {
        let status = if e.count_blockers.is_empty() {
            if e.ic_count_frontier {
                "observed IC recipe frontier".into()
            } else if e.count_frontier {
                "observed count frontier".into()
            } else {
                "dominated count band".into()
            }
        } else {
            e.count_blockers.join("; ")
        };
        out.push_str(&format!(
            "| {} / {} | {} / {} | {} / {} / {} | {} | {} | {} | {} | {} |\n",
            board.sources[e.source].session_id,
            e.workload.workload_id,
            cell(&e.backend),
            cell(&e.arm),
            e.verified,
            e.attempted,
            e.expected,
            e.band_mean_s
                .iter()
                .map(|v| number(*v))
                .collect::<Vec<_>>()
                .join(" / "),
            number(e.cold_median_ms),
            number(e.online_median_ms),
            e.peak_rss_kib,
            cell(&status)
        ));
    }
    out.push_str("\nComplete method parameters, curve/EC1/UID, subgroup r, target, factor-base identity, per-run phases, native counters, SAT conflicts/clauses/learned clauses, solving degrees, statuses, memory and evidence hashes are in [leaderboard.json](leaderboard.json). The interactive view is [index.html](index.html).\n\n## Evidence\n\n");
    for s in &board.sources {
        out.push_str(&format!("- `{}`: `{}`, {}, integrity {}; host {} / {}; planned {}, written {}; bound historical replay receipt: {}. [Session](https://github.com/aburan28/crypto/tree/main/{}).\n",
            s.session_id,s.dir,s.status,s.integrity_ok,cell(s.cpu.as_deref().unwrap_or("unknown")),s.architecture,s.executions_planned,s.records_written,s.receipt_bound,s.dir));
        for p in s.integrity_problems.iter().chain(&s.receipt_problems) {
            out.push_str(&format!("  - {}\n", cell(p)));
        }
    }
    out.push_str("\n## Move the frontier\n\n1. Select an exact workload and observed incumbent from the JSON. Freeze its source/binary, configuration, input hashes, seeds, output checks, resource envelope and accounting interval before timing.\n2. Change one solver, conflict/learned-clause policy, F4/F5/F6 engine, factor base, root-finding strategy or combination. Keep a matching baseline and fresh holdouts. Native counters are diagnostic; fewer conflicts alone does not establish an end-to-end improvement.\n3. Use `ecbench run` for public single-target workloads. Retain every timeout, OOM, budget exhaustion and wrong answer. Run A/A, then at least five interleaved paired rounds; wall claims require L2 or higher and an interval outside parity and A/A noise.\n4. Run `ecbench verify --replay-all` on an independent host. Compare the fastest available external implementation on the same inputs and host; name unavailable comparisons. A historical replay badge here is not independent-host promotion.\n5. Add the sealed session and receipt to the manifest, rebuild all three views, and submit the evidence and benchmark report by PR. Promote the baseline only after full-cost, correctness, independent replay and holdout gates pass. Growing-size comparisons, including degree 83 on the exact E_0 family, remain required for scaling claims.\n");
    out
}

pub fn html(board: &Leaderboard) -> Result<String, String> {
    let mut value = serde_json::to_value(board).map_err(|e| e.to_string())?;
    // Keep u64 counters exact in the browser too; JS numbers cannot hold them all.
    fn exact(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Number(n)
                if n.as_u64().is_some_and(|n| n > 9_007_199_254_740_991) =>
            {
                *v = serde_json::Value::String(n.to_string())
            }
            serde_json::Value::Array(a) => a.iter_mut().for_each(exact),
            serde_json::Value::Object(o) => o.values_mut().for_each(exact),
            _ => (),
        }
    }
    exact(&mut value);
    let data = serde_json::to_string(&value)
        .map_err(|e| e.to_string())?
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('&', "\\u0026");
    Ok(include_str!("leaderboard.html").replace("__LEADERBOARD_DATA__", &data))
}

/// All views share one native producer. `--check` is read-only and fails closed.
pub fn write(board: &Leaderboard, dir: &Path, check: bool) -> Result<(), String> {
    let json = serde_json::to_string(board).map_err(|e| e.to_string())? + "\n";
    let outputs = [
        ("leaderboard.json", json),
        ("LEADERBOARD.md", markdown(board)),
        ("index.html", html(board)?),
    ];
    if !check {
        std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    }
    for (name, text) in outputs {
        let p = dir.join(name);
        if check {
            if std::fs::read_to_string(&p).ok().as_deref() != Some(&text) {
                return Err(format!("stale backend leaderboard: {}", p.display()));
            }
        } else {
            std::fs::write(&p, text).map_err(|e| format!("{}: {e}", p.display()))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::evidence::Curve;
    use super::*;

    fn repository() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    #[test]
    fn historical_inputs_keep_backend_coverage_and_unfinished_attempts() {
        let root = repository();
        let board = build(&root, &root.join("docs/backend-leaderboard/manifest.json")).unwrap();
        for name in [
            "F4 matrix",
            "F4 inherited",
            "F5 matrix",
            "F6-IC",
            "Root finding",
            "SAT",
        ] {
            assert!(
                board.entries.iter().any(|e| e.backend == name),
                "missing backend {name}"
            );
        }
        assert!(board
            .coverage
            .iter()
            .any(|c| c.field_degree == 51 && c.backends.is_empty()));
        assert!(board
            .unattempted
            .iter()
            .any(|e| e.workload.curve.field_degree == Some(31) && e.backend == "F6-IC"));
        assert!(board
            .entries
            .iter()
            .any(|e| e.samples.iter().any(|r| r.outcome.status != "verified")));
        assert!(board.entries.iter().flat_map(|e| &e.samples).any(|r| r
            .solver
            .as_ref()
            .is_some_and(|s| s.op_unit == "conflicts" && s.ops > 0)));
        assert!(
            board.wall_comparisons.iter().all(|c| c.winner.is_none()),
            "L0/L1 historical observations cannot become timing wins"
        );
        assert!(board
            .entries
            .iter()
            .filter(|e| !board.sources[e.source].receipt_bound)
            .all(|e| !e.count_frontier));
        assert!(
            board
                .entries
                .iter()
                .any(|e| e.ic_count_frontier && !e.count_frontier),
            "the IC incumbent remains visible while a reference has a lower band"
        );
    }

    #[test]
    fn changed_raw_records_and_misbound_receipts_fail_closed() {
        let root = repository();
        let dir = std::env::temp_dir().join(format!("backend-import-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("s")).unwrap();
        for name in [
            "session.json",
            "host.json",
            "spec.json",
            "plan.json",
            "records.jsonl",
        ] {
            std::fs::copy(
                root.join("research/f4_f5_f6_cross_family_20261011/sessions/K13")
                    .join(name),
                dir.join("s").join(name),
            )
            .unwrap();
        }
        let receipt_path =
            root.join("research/f4_f5_f6_cross_family_20261011/audits/K13-replay-all.json");
        std::fs::copy(receipt_path, dir.join("audit.json")).unwrap();
        std::fs::write(
            dir.join("manifest.json"),
            serde_json::to_vec(&serde_json::json!({
                "schema":MANIFEST_SCHEMA,"description":"test","field_degrees":[51],
                "sessions":[{"dir":"s","replay_receipt":"audit.json"}]
            }))
            .unwrap(),
        )
        .unwrap();
        let good = build(&dir, &dir.join("manifest.json")).unwrap();
        assert!(good.sources[0].integrity_ok && good.sources[0].receipt_bound);
        let mut receipt: serde_json::Value =
            super::super::evidence::read_json(&dir.join("audit.json")).unwrap();
        receipt["files"]["records.jsonl"] = "wrong hash".into();
        std::fs::write(
            dir.join("audit.json"),
            serde_json::to_vec(&receipt).unwrap(),
        )
        .unwrap();
        let bad = build(&dir, &dir.join("manifest.json")).unwrap();
        assert!(!bad.sources[0].receipt_bound);
        assert!(bad.entries.iter().all(|e| !e.count_frontier));
        let records = std::fs::read_to_string(dir.join("s/records.jsonl")).unwrap();
        std::fs::write(
            dir.join("s/records.jsonl"),
            records.replacen("\"status\":\"verified\"", "\"status\":\"wrong_answer\"", 1),
        )
        .unwrap();
        let bad = build(&dir, &dir.join("manifest.json")).unwrap();
        assert!(!bad.sources[0].integrity_ok);
        assert!(bad.entries.iter().all(|e| !e.count_frontier));
        std::fs::remove_dir_all(dir).unwrap();
    }

    fn entry(arm: &str, band: [f64; 3], wall: u64) -> Entry {
        let workload = Workload {
            workload_id: "target".into(),
            workload_sha256: "target-sha".into(),
            curve: Curve {
                slug: "synthetic-test-model".into(),
                family: "koblitz".into(),
                field_degree: Some(13),
                field_bits: None,
                r: 2003,
                rest: BTreeMap::new(),
            },
            target: ["1".into(), "2".into()],
            rest: BTreeMap::new(),
        };
        let method = ResolvedMethod {
            id: "rho.negation".into(),
            family: "rho".into(),
            params: BTreeMap::new(),
            method_id: "method".into(),
            method_sha256: "method-sha".into(),
            entry: "synthetic-test".into(),
        };
        let samples = (0..5)
            .map(|round| Sample {
                record_id: format!("record-{arm}-{round}"),
                run_id: format!("run-{arm}-{round}"),
                seq: round as u64,
                round,
                algorithm_seed: round.to_string(),
                outcome: Outcome {
                    status: "verified".into(),
                    recovered: Some("1".into()),
                    matches_target: Some(true),
                    matches_planted: None,
                    error: None,
                    exit: None,
                    stderr_sha256: None,
                },
                cost: Cost {
                    unit: "ecbench.gae".into(),
                    total_gae: Some(band[0] * 2003f64.sqrt()),
                    s: Some(band[0]),
                    lower_bound: false,
                    unpriced: vec![],
                    deterministic: true,
                    nondeterminism: vec![],
                },
                recorded_s: Some(band[0]),
                band_s: band.into_iter().map(Some).collect(),
                unpriced_counts: BTreeMap::new(),
                phases: vec![],
                counters: BTreeMap::new(),
                solver: None,
                solver_origin: None,
                solver_context: serde_json::Value::Null,
                fbm_id: None,
                field_ops: None,
                online: None,
                time: Timing {
                    solve_wall_ns: Some(wall),
                    ..Default::default()
                },
                isolation: Level::L2,
                isolation_blockers: vec![],
            })
            .collect();
        Entry {
            entry_id: arm.into(),
            source: 0,
            arm: arm.into(),
            role: "candidate".into(),
            backend: "test".into(),
            method,
            workload,
            factor_base: None,
            attempted: 5,
            expected: 5,
            verified: 5,
            statuses: BTreeMap::from([("verified".into(), 5)]),
            lower_bound: false,
            band_mean_s: band.into_iter().map(Some).collect(),
            cold_median_ms: Some(wall as f64 / 1e6),
            cold_range_ms: None,
            process_median_ms: None,
            online_median_ms: None,
            peak_rss_kib: 100,
            solver_median_ops: None,
            solver_op_unit: None,
            sat_conflicts_median: None,
            sat_learnt_clauses_median: None,
            count_blockers: vec![],
            full_historical_replay: false,
            count_dominated_by: vec![],
            ic_count_dominated_by: vec![],
            ic_count_frontier: false,
            count_frontier: false,
            resource_frontier: false,
            samples,
        }
    }

    #[test]
    fn incomplete_online_phase_sets_remain_unknown_even_if_the_sum_matches() {
        let mut window = OnlineWindow {
            wall_ns: 1,
            start_event: "start".into(),
            stop_event: "verified".into(),
            rest: BTreeMap::new(),
            phases_ns: BTreeMap::from([
                ("target_query".into(), 1),
                ("target_PDP".into(), 0),
                ("target_relation_check".into(), 0),
                ("target_descent".into(), 0),
                ("target_recovery_check".into(), 0),
            ]),
        };
        assert!(complete_online(&window, "ic"));
        window.phases_ns.remove("target_recovery_check");
        assert!(!complete_online(&window, "ic"));
        window.phases_ns.insert("target_recovery_check".into(), 0);
        window.wall_ns = 2;
        assert!(!complete_online(&window, "ic"));
    }

    #[test]
    fn crossing_bands_unknowns_and_failures_cannot_dominate() {
        let a = entry("a", [1., 50., 100.], 10);
        let mut b = entry("b", [2., 40., 80.], 10);
        assert!(!dominates(&a, &b, false));
        assert!(!dominates(&b, &a, false));
        b.band_mean_s = vec![Some(2.), Some(60.), None];
        assert!(!dominates(&a, &b, false));
        b.band_mean_s = vec![Some(2.), Some(60.), Some(120.)];
        assert!(dominates(&a, &b, false));
        b.count_blockers.push("timeout".into());
        assert!(!dominates(&a, &b, false));
    }

    #[test]
    fn no_cross_session_target_seed_or_memory_shortcuts() {
        let a = entry("a", [1., 1., 1.], 10);
        let mut b = entry("b", [2., 2., 2.], 10);
        assert!(dominates(&a, &b, false));
        b.source = 1;
        assert!(!dominates(&a, &b, false));
        b.source = 0;
        b.workload.workload_sha256.push('x');
        assert!(!dominates(&a, &b, false));
        b.workload.workload_sha256.pop();
        b.samples[0].algorithm_seed.push('x');
        assert!(!dominates(&a, &b, false));
        b.samples[0].algorithm_seed.pop();
        b.peak_rss_kib = 50;
        assert!(!dominates(&a, &b, true));
        b.peak_rss_kib = 0;
        assert!(!dominates(&a, &b, true));
    }

    #[test]
    fn timing_requires_complete_l2_pairs_and_same_target_aa() {
        let a = entry("a", [1.; 3], 10_000);
        let mut b = entry("b", [1.; 3], 5_000);
        b.method.method_id = "different-method".into();
        let aa = entry("aa", [1.; 3], 10_000);
        let mut entries = vec![a, b, aa];
        let c = wall_comparison(&entries[0], &entries[1], &entries, Level::L0);
        assert_eq!(c.winner.as_deref(), Some("b"));
        entries[1].samples[0].isolation = Level::L0;
        assert!(
            wall_comparison(&entries[0], &entries[1], &entries, Level::L0)
                .winner
                .is_none()
        );
        entries[1].samples[0].isolation = Level::L2;
        entries[2].workload.workload_id = "another-target".into();
        assert!(
            wall_comparison(&entries[0], &entries[1], &entries, Level::L2)
                .winner
                .is_none()
        );
        entries[2].workload.workload_id = entries[0].workload.workload_id.clone();
        entries[1].verified = 4;
        assert!(
            wall_comparison(&entries[0], &entries[1], &entries, Level::L2)
                .winner
                .is_none()
        );
    }

    #[test]
    fn categories_read_configuration_and_preserve_root_finding() {
        let mut e = entry("misleading-f6-label", [1.; 3], 1);
        e.method
            .params
            .insert("oracle".into(), "pdp3-binary:m=3,engine=matrix-f5".into());
        assert_eq!(backend(&e.method), "F5 matrix");
        e.method.params.insert(
            "oracle".into(),
            "pdp3-koblitz:m=3,engine=rootfind-window".into(),
        );
        assert_eq!(backend(&e.method), "Root finding");
        e.method
            .params
            .insert("oracle".into(), "descent-algebraic:m=2".into());
        e.method
            .params
            .insert("solver".into(), "sat-cdcl:sat_xor_encoding=native".into());
        assert_eq!(backend(&e.method), "SAT");
    }

    #[test]
    fn output_escapes_embedded_json_and_check_never_repairs_staleness() {
        let mut e = entry("a", [1.; 3], 1);
        e.samples[0].counters.insert("large".into(), u64::MAX);
        let mut b = Leaderboard {
            schema: SCHEMA,
            description: "</script><img onerror='bad'>".into(),
            manifest_sha256: "test".into(),
            producer_input_hashes: BTreeMap::new(),
            scope: "test",
            weights_gae_per_native_op: vec![0., 0.1, 1.],
            sources: vec![],
            coverage: vec![],
            entries: vec![],
            unattempted: vec![],
            wall_comparisons: vec![],
        };
        b.entries.push(e);
        let h = html(&b).unwrap();
        assert!(!h.contains("</script><img"));
        assert!(h.contains("\\u003c/script\\u003e"));
        assert!(h.contains("\"large\":\"18446744073709551615\""));
        b.entries.clear();
        assert!(relative(Path::new("../outside")).is_err());
        assert!(relative(Path::new("/absolute")).is_err());
        let dir = std::env::temp_dir().join(format!("backend-board-check-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        write(&b, &dir, false).unwrap();
        write(&b, &dir, true).unwrap();
        std::fs::write(dir.join("index.html"), "stale").unwrap();
        assert!(write(&b, &dir, true).is_err());
        assert_eq!(
            std::fs::read_to_string(dir.join("index.html")).unwrap(),
            "stale"
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
