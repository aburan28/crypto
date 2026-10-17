//! Read-only ecbench schema projection. Unknown fields in source records are
//! retained where exposed. This is a seal/file/plan integrity check, not a fresh
//! arithmetic replay: count admission additionally needs the hash-bound native
//! audit receipt from the historical experiment.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::canonical::sha256_hex;

pub const MIN_WALL_PAIRS: usize = 5;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Level {
    L0,
    L1,
    L2,
    L3,
}
impl Level {
    pub fn name(self) -> &'static str {
        match self {
            Self::L0 => "L0",
            Self::L1 => "L1",
            Self::L2 => "L2",
            Self::L3 => "L3",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResolvedMethod {
    pub id: String,
    pub family: String,
    pub params: BTreeMap<String, String>,
    pub method_id: String,
    pub method_sha256: String,
    pub entry: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Curve {
    pub slug: String,
    pub family: String,
    pub field_degree: Option<u32>,
    pub field_bits: Option<u32>,
    #[serde(with = "crate::canonical::compat_u128")]
    pub r: u128,
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Workload {
    pub workload_id: String,
    pub workload_sha256: String,
    pub curve: Curve,
    pub target: [String; 2],
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhaseRecord {
    pub name: String,
    pub gae: f64,
    pub wall_ns: Option<u64>,
    #[serde(default)]
    pub native: BTreeMap<String, u64>,
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

pub type FactorBaseFacts = Value;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SolverStats {
    pub ops: u64,
    pub op_unit: String,
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OnlineWindow {
    pub wall_ns: u64,
    pub phases_ns: BTreeMap<String, u64>,
    pub start_event: String,
    pub stop_event: String,
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Timing {
    pub solve_wall_ns: Option<u64>,
    pub process_wall_ns: u64,
    pub max_rss_kib: u64,
    #[serde(flatten)]
    pub rest: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Outcome {
    pub status: String,
    pub recovered: Option<String>,
    pub matches_target: Option<bool>,
    pub matches_planted: Option<bool>,
    pub error: Option<String>,
    pub exit: Option<String>,
    pub stderr_sha256: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Cost {
    pub unit: String,
    pub total_gae: Option<f64>,
    pub s: Option<f64>,
    pub lower_bound: bool,
    pub unpriced: Vec<String>,
    pub deterministic: bool,
    #[serde(default)]
    pub nondeterminism: Vec<String>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Isolation {
    pub level: Level,
    pub blockers: Vec<String>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Record {
    pub schema: String,
    pub record_id: String,
    pub run_id: String,
    pub session_id: String,
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
    pub phases: Vec<PhaseRecord>,
    pub counters: BTreeMap<String, u64>,
    pub factor_base: Option<FactorBaseFacts>,
    pub fbm_id: Option<String>,
    pub field_ops: Option<Value>,
    pub detail: Value,
    pub online: Option<OnlineWindow>,
    pub solver: Option<SolverStats>,
    #[serde(skip)]
    pub solver_origin: Option<String>,
    pub time: Timing,
    pub isolation: Isolation,
    pub env_class_id: String,
    pub binary_sha256: Option<String>,
}

impl Record {
    pub fn counts(&self) -> bool {
        !self.warmup && self.outcome.status == "verified"
    }
}

#[derive(Debug, Deserialize)]
pub struct Session {
    pub schema: String,
    pub session_id: String,
    pub status: String,
    pub git_commit: Option<String>,
    pub git_dirty: Option<bool>,
    pub binary_sha256: Option<String>,
    pub env_class_id: String,
    pub child_env: BTreeMap<String, String>,
    pub options: Value,
    pub cpu_plan: Value,
    pub executions_planned: u64,
    pub records_written: u64,
    pub records_sha256: Option<String>,
    pub host_sha256: String,
    pub plan_sha256: String,
}

#[derive(Debug, Deserialize)]
pub struct HostCapsule {
    pub stable: HostFacts,
}
#[derive(Debug, Deserialize)]
pub struct HostFacts {
    pub cpu_model: Option<String>,
    pub arch: String,
}
#[derive(Debug, Deserialize)]
pub struct Spec {
    pub schema: String,
    pub measurement: Measurement,
}
#[derive(Debug, Deserialize)]
pub struct Measurement {
    pub rounds: u32,
    pub isolation_required: Level,
}
impl Spec {
    pub fn from_json(s: &str) -> Result<Self, String> {
        let spec: Self = serde_json::from_str(s).map_err(|e| e.to_string())?;
        if spec.schema != "ecbench.spec/v1" || spec.measurement.rounds == 0 {
            return Err("invalid ecbench spec".into());
        }
        Ok(spec)
    }
}
#[derive(Debug, Deserialize)]
pub struct PlanArm {
    pub name: String,
    pub method: ResolvedMethod,
}
#[derive(Debug, Deserialize)]
pub struct Execution {
    pub seq: u64,
    pub round: u32,
    pub warmup: bool,
    pub arm: usize,
    pub workload: usize,
    pub algorithm_seed: u64,
}
#[derive(Debug, Deserialize)]
pub struct PlanDoc {
    pub arms: Vec<PlanArm>,
    pub workloads: Vec<Workload>,
    pub executions: Vec<Execution>,
}

#[derive(Debug, Deserialize)]
pub struct Replay {
    pub record_id: String,
    pub reproduced: bool,
}
#[derive(Debug, Deserialize)]
pub struct AuditReport {
    pub schema: String,
    pub session_id: String,
    pub ok: bool,
    #[serde(default)]
    pub records_bound: Option<bool>,
    #[serde(default)]
    pub problems: Vec<String>,
    #[serde(default)]
    pub replays: Vec<Replay>,
    pub files: BTreeMap<String, String>,
}

pub fn read_session(dir: &Path) -> Result<Session, String> {
    read_json(&dir.join("session.json"))
}
pub fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, String> {
    serde_json::from_slice(&std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?)
        .map_err(|e| format!("{}: {e}", path.display()))
}
pub fn read_records(dir: &Path) -> Result<Vec<Record>, String> {
    let text = std::fs::read_to_string(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
    text.lines()
        .enumerate()
        .map(|(i, line)| {
            let mut r: Record = serde_json::from_str(line)
                .map_err(|e| format!("{} records line {}: {e}", dir.display(), i + 1))?;
            if r.solver.is_some() {
                r.solver_origin = Some("/solver".into());
            } else if r
                .cost
                .unpriced
                .iter()
                .any(|name| name.starts_with("solver_"))
            {
                // Before the top-level block existed, counted SAT records
                // emitted this same observation under detail. Adapt explicitly
                // unpriced solver work only; never invent a calibration.
                if let Some(value) = r.detail.pointer("/decomposition/solver") {
                    if !value.is_null() {
                        r.solver = Some(serde_json::from_value(value.clone()).map_err(|e| {
                            format!("{} record {} legacy solver: {e}", dir.display(), r.seq)
                        })?);
                        r.solver_origin = Some(
                            "/detail/decomposition/solver (historical unpriced observation)".into(),
                        );
                    }
                }
            }
            Ok(r)
        })
        .collect()
}

pub fn audit_with(dir: &Path, _: usize, _: bool) -> Result<AuditReport, String> {
    let session = read_session(dir)?;
    let plan: PlanDoc = read_json(&dir.join("plan.json"))?;
    let raw = std::fs::read_to_string(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
    let records = read_records(dir)?;
    let mut files = BTreeMap::new();
    let mut problems = Vec::new();
    for name in [
        "session.json",
        "spec.json",
        "plan.json",
        "host.json",
        "records.jsonl",
    ] {
        let hash = sha256_hex(&std::fs::read(dir.join(name)).map_err(|e| e.to_string())?);
        files.insert(name.to_string(), hash);
    }
    // spec_sha256 is the harness's *normalised identity*, not the file hash.
    // The historical native audit binds spec.json's exact bytes below.
    for (name, hash) in [
        ("plan.json", &session.plan_sha256),
        ("host.json", &session.host_sha256),
    ] {
        if files.get(name) != Some(hash) {
            problems.push(format!("file hash mismatch: {name}"));
        }
    }
    let bound = session.records_sha256.as_ref() == files.get("records.jsonl");
    if !bound {
        problems.push("records file is unbound or its hash mismatches".into());
    }
    if session.schema != "ecbench.session/v1" {
        problems.push("unsupported session schema".into());
    }
    if session.executions_planned != plan.executions.len() as u64
        || (session.status == "complete" && records.len() != plan.executions.len())
        || (session.status != "running" && records.len() as u64 != session.records_written)
    {
        problems.push("record or plan count mismatch".into());
    }
    for (i, (r, line)) in records.iter().zip(raw.lines()).enumerate() {
        let field = format!("\"record_id\":\"{}\"", r.record_id);
        let empty = line.replacen(&field, "\"record_id\":\"\"", 1);
        let valid_seal = line.matches(&field).count() == 1
            && r.record_id == format!("ECR1h{}", &sha256_hex(empty.as_bytes())[..12]);
        let valid_plan = plan.executions.get(i).is_some_and(|x| {
            x.seq == r.seq
                && x.round == r.round
                && x.warmup == r.warmup
                && x.algorithm_seed.to_string() == r.algorithm_seed
                && plan
                    .arms
                    .get(x.arm)
                    .is_some_and(|a| a.name == r.arm && a.method == r.method)
                && plan
                    .workloads
                    .get(x.workload)
                    .is_some_and(|w| *w == r.workload)
        });
        if !valid_seal
            || !valid_plan
            || r.seq != i as u64
            || r.session_id != session.session_id
            || r.schema != "ecbench.record/v1"
            || r.env_class_id != session.env_class_id
            || r.binary_sha256 != session.binary_sha256
            || r.cost.unit != "ecbench.gae"
        {
            problems.push(format!(
                "record {}: seal, plan, identity or unit mismatch",
                r.seq
            ));
        }
        if r.counts() && (r.outcome.matches_target != Some(true) || r.outcome.recovered.is_none()) {
            problems.push(format!(
                "record {}: verified status without a checked answer",
                r.seq
            ));
        }
    }
    Ok(AuditReport {
        schema: "backend-leaderboard.input-integrity/v1".into(),
        session_id: session.session_id,
        ok: problems.is_empty(),
        records_bound: Some(bound),
        problems,
        replays: vec![],
        files,
    })
}

pub fn pinned_calibration_ratios(
    regime: &str,
    instance: &str,
) -> Option<serde_json::Map<String, Value>> {
    static TABLE: OnceLock<Value> = OnceLock::new();
    static ALIASES: OnceLock<Value> = OnceLock::new();
    let table = TABLE.get_or_init(|| {
        serde_json::from_str(include_str!("../../../docs/ic/calibration.json"))
            .expect("calibration JSON")
    });
    let aliases = ALIASES.get_or_init(|| {
        serde_json::from_str(include_str!(
            "../../../src/cryptanalysis/curve_aliases.json"
        ))
        .expect("curve aliases")
    });
    let resolve = |s: &str| {
        let key: String = s
            .trim()
            .trim_matches('`')
            .chars()
            .filter(|c| !c.is_whitespace() && !matches!(c, '{' | '}' | '_'))
            .map(|c| match c {
                '₀'..='₉' => char::from(b'0' + (c as u32 - '₀' as u32) as u8),
                _ => c.to_ascii_lowercase(),
            })
            .collect();
        aliases
            .get("aliases")
            .and_then(|a| a.get(&key))
            .and_then(Value::as_str)
            .unwrap_or(s)
            .to_string()
    };
    let prefix = format!("{regime}/");
    let rows = table.get("instances")?.as_object()?;
    let row = rows.get(&format!("{prefix}{instance}")).or_else(|| {
        rows.iter().find_map(|(key, row)| {
            let name = key.strip_prefix(&prefix)?;
            (resolve(name) == resolve(instance)).then_some(row)
        })
    })?;
    row.as_object().cloned()
}
