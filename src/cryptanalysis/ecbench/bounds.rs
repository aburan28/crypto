//! Bounds: what a method costs, as a sealed record another session can
//! cite, re-derive and challenge.
//!
//! A *bound* takes one arm of one or more sessions and states, in the
//! harness's unit, the cost of that method on that problem: the constant
//! at the declared exponent (`S = gae / √r`, and the dimensionless ratio
//! to the curve's generic floor `√(π / 2A)`), the exponent fitted across
//! sizes (`gae = C · r^α`, least squares in log–log space, as
//! `ca_bench complexity` fits it in aburan28/cryptanalysis), both with
//! two-stage bootstrap intervals, the share and exponent of every phase,
//! the memory the method stored and the work it counted but did not
//! charge, and the sessions and audit receipts it was read from.
//!
//! A bound is scoped by its *domain* — problem, curve family, target
//! kind, unit, tier, resource envelope — and compares with nothing outside
//! it.  Its id is the SHA-256 of its own bytes (`seal_document`), so a
//! changed bound is a new bound, two sessions fitting the same sessions
//! write the same file, and `bound check` re-derives a committed bound
//! from the sessions it names.  `docs/bounds/README.md` is the protocol;
//! [`frontier`](super::frontier) ranks bounds and
//! [`challenge`](super::challenge) replaces them.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::cryptanalysis::ecbench::audit::AuditReport;
use crate::cryptanalysis::ecbench::canonical::{compat_u128, sha256_hex, short_id};
use crate::cryptanalysis::ecbench::methods::expected_s;
use crate::cryptanalysis::ecbench::record::{floor_s, Record};
use crate::cryptanalysis::ecbench::runner::{read_records, read_session, PlanDoc};
use crate::cryptanalysis::ecbench::stats::{bootstrap_ci, cluster_bootstrap_ci, mean};

pub const BOUND_SCHEMA: &str = "ecbench.bound/v1";
pub const BOUND_PREFIX: &str = "ECBND1h";
pub const DOMAIN_PREFIX: &str = "ECDOM1";
/// The problem every ecbench workload poses: one cold target, one
/// subgroup, no precomputation.
pub const PROBLEM: &str = "ecdlp.single_target";
/// Bootstrap resamples and seed, fixed so a bound re-derives bit for bit.
pub const RESAMPLES: usize = 2000;
pub const SEED: u64 = 20261005;
/// A scaling claim (an exponent) needs this many sizes (AGENTS.md §5).
pub const MIN_SIZES_FOR_SCALING: usize = 4;

// ── Documents sealed by their bytes ────────────────────────────────

/// Seal a pretty-printed JSON document: `field` is written empty, the
/// text is hashed, and the id (`prefix` + 12 hex) is written in.  Returns
/// the id and the text.  The seal covers bytes, not a parsed value, for
/// the reason `Record::seal` gives: a float parsed and re-serialised can
/// differ in a last digit from what was written.
pub fn seal_document<T: Serialize>(
    doc: &T,
    field: &str,
    prefix: &str,
) -> Result<(String, String), String> {
    let mut v = serde_json::to_value(doc).map_err(|e| e.to_string())?;
    let obj = v
        .as_object_mut()
        .ok_or("a sealed document is a JSON object")?;
    if !obj.contains_key(field) {
        return Err(format!("document has no `{field}` field to seal"));
    }
    obj.insert(field.to_string(), Value::String(String::new()));
    // Key order: the struct's own order is lost in a `Value`; the text is
    // what is hashed, so whatever order it has is the order that counts.
    let blank = serde_json::to_string_pretty(&v).map_err(|e| e.to_string())? + "\n";
    let needle = format!("\"{field}\": \"\"");
    if blank.matches(&needle).count() != 1 {
        return Err(format!(
            "the `{field}` field must appear exactly once in the document"
        ));
    }
    let id = format!("{prefix}{}", &sha256_hex(blank.as_bytes())[..12]);
    let text = blank.replacen(&needle, &format!("\"{field}\": \"{id}\""), 1);
    Ok((id, text))
}

/// Whether `text` hashes to the id it carries in `field`; the id itself.
pub fn check_document_seal(text: &str, field: &str, prefix: &str) -> Result<String, String> {
    let v: Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
    let id = v
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("no `{field}`"))?
        .to_string();
    if !id.starts_with(prefix) {
        return Err(format!("`{field}` {id} does not start with {prefix}"));
    }
    let filled = format!("\"{field}\": \"{id}\"");
    if text.matches(&filled).count() != 1 {
        return Err(format!("`{field}` must appear exactly once as written"));
    }
    let blank = text.replacen(&filled, &format!("\"{field}\": \"\""), 1);
    let expect = format!("{prefix}{}", &sha256_hex(blank.as_bytes())[..12]);
    if expect != id {
        return Err(format!(
            "seal mismatch: bytes hash to {expect}, document says {id}"
        ));
    }
    Ok(id)
}

/// Read and parse a sealed document, checking its seal first.
pub fn read_sealed<T: DeserializeOwned>(
    path: &Path,
    field: &str,
    prefix: &str,
) -> Result<T, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    check_document_seal(&text, field, prefix).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))
}

// ── The domain ─────────────────────────────────────────────────────

/// The resource envelope a bound was measured in.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Envelope {
    pub targets: u32,
    /// `none`: every run starts cold.
    pub precomputation: String,
    pub threads: u32,
}

impl Default for Envelope {
    fn default() -> Self {
        Self {
            targets: 1,
            precomputation: "none".into(),
            threads: 1,
        }
    }
}

/// What a bound is about.  Bounds compare only inside one domain.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Domain {
    pub problem: String,
    /// `prime`, `koblitz` or `binary`: the automorphism group, and so the
    /// floor, differs.
    pub family: String,
    /// `planted` or `public`.
    pub target_kind: String,
    pub unit: String,
    /// `toy`, `medium` or `crypto`, from the field size.
    pub tier: String,
    pub envelope: Envelope,
}

impl Domain {
    /// `ECDOM1h` + 12 hex of the canonical domain.
    pub fn id(&self) -> Result<String, String> {
        let v = serde_json::to_value(self).map_err(|e| e.to_string())?;
        Ok(short_id(DOMAIN_PREFIX, &v)?.0)
    }
}

/// The claim tier a field size contributes to (crypto-autoresearcher,
/// `docs/claims-and-verification.md`): `toy` up to 32 bits, `medium` up
/// to 96, `crypto` above.
pub fn tier_of_bits(field_bits: u32) -> &'static str {
    if field_bits <= 32 {
        "toy"
    } else if field_bits <= 96 {
        "medium"
    } else {
        "crypto"
    }
}

/// The field size of a record's curve: bits of `p`, or the degree over
/// `F_2`.
pub fn field_size(r: &Record) -> Option<u32> {
    let c = &r.workload.curve;
    if c.family == "prime" {
        c.field_bits
    } else {
        c.field_degree
    }
}

// ── The record ─────────────────────────────────────────────────────

/// A statistic with its interval.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Estimate {
    pub value: f64,
    pub ci95: Option<(f64, f64)>,
    /// `cluster` (sizes, then runs), `within` (runs of one size) or `none`.
    pub ci_method: String,
}

/// The method a bound is about: the configuration (`ECM1h…`), never the
/// code, whose hash is in the provenance.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MethodRef {
    pub id: String,
    pub params: BTreeMap<String, String>,
    pub method_id: String,
    pub family: String,
    pub entry: String,
}

/// One curve size.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SizeRow {
    pub slug: String,
    pub log2_r: f64,
    #[serde(with = "compat_u128")]
    pub r: u128,
    pub field_bits: u32,
    pub automorphisms_available: u32,
    pub floor_s: f64,
    pub workloads: u64,
    /// Measured runs of the arm on this curve, and those that verified.
    pub runs: u64,
    pub verified: u64,
    pub mean_gae: Option<f64>,
    pub mean_s: Option<f64>,
    pub s_ci95: Option<(f64, f64)>,
    pub ratio_to_floor: Option<f64>,
    /// `expected_s` of the method on this curve, and its ratio to the
    /// floor, where the registry derives one.
    pub declared_s: Option<f64>,
    pub declared_ratio_to_floor: Option<f64>,
    pub memory_entries_per_sqrt_r: Option<f64>,
    pub uncharged_per_sqrt_r: Option<f64>,
    pub levels: BTreeMap<String, u64>,
}

/// `gae = C · r^α`, least squares on every verified run's `(ln r, ln gae)`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fit {
    pub size_parameter: String,
    pub model: String,
    pub points: u64,
    pub sizes: u64,
    pub alpha: Option<f64>,
    pub alpha_ci95: Option<(f64, f64)>,
    pub log2_c: Option<f64>,
    pub r_squared: Option<f64>,
    pub ci_method: String,
    /// `0.5` for every method the registry derives an expectation for.
    pub declared_alpha: Option<f64>,
    /// The declared exponent lies inside the fitted interval.
    pub alpha_agrees_with_declared: Option<bool>,
    /// At least [`MIN_SIZES_FOR_SCALING`] sizes: the exponent is a result,
    /// not a description.
    pub scaling_claim: bool,
}

/// The constant at the declared exponent.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Constant {
    /// Mean `S = gae / √r` over verified runs.
    pub s: Estimate,
    /// Mean `S / √(π / 2A)`: comparable across curve families.
    pub ratio_to_floor: Estimate,
    /// The method's derived ratio to the floor when it is the same on
    /// every size (rho.negation: 1; bsgs.textbook: 1.5 / √(π/4)).
    pub declared_ratio_to_floor: Option<f64>,
}

/// One cost axis.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Dimension {
    pub statistic: String,
    pub source: String,
    pub lower_is_better: bool,
    /// Every verified run reported it.  Unknown is not zero.
    pub known: bool,
    pub value: Option<f64>,
    pub ci95: Option<(f64, f64)>,
}

/// One phase of the method: its share of the total and its own exponent.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StageRow {
    pub name: String,
    pub mean_gae: f64,
    /// `Σ phase gae / Σ total gae` over verified runs.
    pub share_of_gae: f64,
    pub sizes_with_work: u64,
    pub alpha: Option<f64>,
    pub alpha_ci95: Option<(f64, f64)>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SessionRef {
    pub dir: String,
    pub session_id: String,
    pub spec_id: String,
    pub status: String,
    pub arm: String,
    pub role: String,
    pub binary_sha256: Option<String>,
    pub git_commit: Option<String>,
    pub env_class_id: String,
    pub records_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AuditRef {
    pub path: String,
    pub sha256: String,
    pub session_id: String,
    pub ok: bool,
    pub replays: u64,
    pub replays_reproduced: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Provenance {
    pub sessions: Vec<SessionRef>,
    pub audits: Vec<AuditRef>,
    pub records: u64,
    pub verified: u64,
    pub statuses: BTreeMap<String, u64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Admissibility {
    /// `admissible`, or `inadmissible` (a measured run did not verify).
    pub status: String,
    /// Some work was counted and not priced: every total is a floor.
    pub bounded: bool,
    pub unpriced: Vec<String>,
    pub deterministic: bool,
    pub reasons: Vec<String>,
}

/// How the bound was fitted, so `bound check` can fit it again.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FitOptions {
    pub tier: Option<String>,
    /// Curve slugs kept; empty keeps every curve of the arm.
    pub curves: Vec<String>,
    pub resamples: usize,
    pub seed: u64,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            tier: None,
            curves: vec![],
            resamples: RESAMPLES,
            seed: SEED,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bound {
    pub schema: String,
    /// `ECBND1h` + 12 hex of this document with the field empty.
    pub bound_id: String,
    pub label: String,
    pub domain: Domain,
    pub domain_id: String,
    pub method: MethodRef,
    /// `exponent` (a scaling claim) or `constant`.
    pub level: String,
    pub sizes: Vec<SizeRow>,
    pub fit: Fit,
    pub constant: Constant,
    pub dimensions: BTreeMap<String, Dimension>,
    pub stages: Vec<StageRow>,
    pub provenance: Provenance,
    pub admissibility: Admissibility,
    pub fit_options: FitOptions,
    /// Bounds this one was measured against and beat or traded with, set
    /// by the verdict that produced it.
    pub improves_on: Vec<String>,
    pub verdict_id: Option<String>,
    pub notes: String,
}

impl Bound {
    /// Seal: compute `bound_id` and return the text to write.
    pub fn seal(&mut self) -> Result<String, String> {
        self.bound_id = String::new();
        let (id, text) = seal_document(self, "bound_id", BOUND_PREFIX)?;
        self.bound_id = id;
        Ok(text)
    }

    pub fn read(path: &Path) -> Result<Self, String> {
        let b: Bound = read_sealed(path, "bound_id", BOUND_PREFIX)?;
        if b.schema != BOUND_SCHEMA {
            return Err(format!(
                "{}: schema `{}`, expected {BOUND_SCHEMA}",
                path.display(),
                b.schema
            ));
        }
        Ok(b)
    }

    /// The ops axis: mean ratio to the floor.
    pub fn ops(&self) -> Option<&Dimension> {
        self.dimensions.get("ops")
    }

    pub fn admissible(&self) -> bool {
        self.admissibility.status == "admissible"
    }
}

// ── Fitting ────────────────────────────────────────────────────────

/// Least squares of `y = alpha · x + b` over `(x, y)`: `(alpha, b, R²)`.
/// `None` below two distinct `x`.
pub fn ols(points: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    let n = points.len() as f64;
    if points.len() < 2 {
        return None;
    }
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for (x, y) in points {
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return None;
    }
    let alpha = (n * sxy - sx * sy) / denom;
    let b = (sy - alpha * sx) / n;
    let ybar = sy / n;
    let (mut ss_res, mut ss_tot) = (0.0, 0.0);
    for (x, y) in points {
        let yh = alpha * x + b;
        ss_res += (y - yh) * (y - yh);
        ss_tot += (y - ybar) * (y - ybar);
    }
    let r2 = if ss_tot > 0.0 {
        1.0 - ss_res / ss_tot
    } else {
        1.0
    };
    Some((alpha, b, r2))
}

/// `gae = C · r^α` over strata of `(ln r, ln gae)` points, one stratum per
/// size: the slope, its two-stage interval (sizes, then runs), the
/// intercept as `log2 C`, and R².
pub fn power_law(
    strata: &[Vec<(f64, f64)>],
    resamples: usize,
    seed: u64,
) -> (
    Option<f64>,
    Option<(f64, f64)>,
    Option<f64>,
    Option<f64>,
    String,
) {
    let all: Vec<(f64, f64)> = strata.iter().flatten().copied().collect();
    let Some((alpha, b, r2)) = ols(&all) else {
        return (None, None, None, None, "none".into());
    };
    let slope = |s: &[Vec<(f64, f64)>]| {
        let pts: Vec<(f64, f64)> = s.iter().flatten().copied().collect();
        ols(&pts).map(|(a, _, _)| a)
    };
    let nonempty: Vec<Vec<(f64, f64)>> = strata.iter().filter(|s| !s.is_empty()).cloned().collect();
    let (ci, method) = if nonempty.len() >= 2 {
        (
            cluster_bootstrap_ci(&nonempty, resamples, seed, slope),
            "cluster".to_string(),
        )
    } else {
        (None, "none".to_string())
    };
    (
        Some(alpha),
        ci,
        Some(b / std::f64::consts::LN_2),
        Some(r2),
        method,
    )
}

/// Mean over strata with its two-stage interval; within-stratum when
/// there is one stratum.
fn estimate(strata: &[Vec<f64>], resamples: usize, seed: u64) -> Option<Estimate> {
    let all: Vec<f64> = strata.iter().flatten().copied().collect();
    let value = mean(&all)?;
    let stat = |s: &[Vec<f64>]| mean(&s.concat());
    let nonempty: Vec<Vec<f64>> = strata.iter().filter(|s| !s.is_empty()).cloned().collect();
    let (ci95, ci_method) = if nonempty.len() >= 2 {
        (
            cluster_bootstrap_ci(&nonempty, resamples, seed, stat),
            "cluster".to_string(),
        )
    } else {
        (
            bootstrap_ci(&nonempty, resamples, seed, stat),
            "within".to_string(),
        )
    };
    Some(Estimate {
        value,
        ci95,
        ci_method,
    })
}

/// The counters a method reports its stored table through, in the order
/// they are read: BSGS, the kangaroo and the claw count `inserts_uncharged`;
/// the strong Koblitz rho counts `table_inserts_uncharged` and reports
/// `table_entries`; the tuned rho walks store one distinguished point per
/// walk (`distinguished_points`).
pub const MEMORY_COUNTERS: &[&str] = &[
    "inserts_uncharged",
    "table_inserts_uncharged",
    "table_entries",
    "distinguished_points",
];

/// Table entries a run stored, as a multiple of `√r`, from the first of
/// [`MEMORY_COUNTERS`] the record carries.  `None` when it carries none:
/// unknown, not zero.
pub fn memory_entries(r: &Record) -> Option<f64> {
    let entries = MEMORY_COUNTERS
        .iter()
        .find_map(|k| r.counters.get(*k))
        .copied()?;
    Some(entries as f64 / (r.workload.curve.r as f64).sqrt())
}

/// Every counter the unit does not charge, summed, per `√r`.
pub fn uncharged_work(r: &Record) -> f64 {
    let total: u64 = r
        .counters
        .iter()
        .filter(|(k, _)| k.ends_with("_uncharged"))
        .map(|(_, v)| *v)
        .sum();
    total as f64 / (r.workload.curve.r as f64).sqrt()
}

/// The primitive-level axes: modular multiplications, squarings and
/// inversions behind the group operations, per `√r`, read from a record's
/// `field_ops` block.  A bound carries them only when every verified run
/// of the arm reports the block (prime-field curves through the generic
/// walks and tables); a verdict reports them always and lets them decide
/// only when the challenge names them.
pub const FIELD_AXES: &[&str] = &["field_muls", "field_sqrs", "field_invs"];

/// What a field axis counts, for the record's `statistic` and `source`.
pub fn field_axis_describes(axis: &str) -> Option<(&'static str, &'static str)> {
    match axis {
        "field_muls" => Some(("field multiplications", "field_ops.muls")),
        "field_sqrs" => Some(("field squarings", "field_ops.sqrs")),
        "field_invs" => Some(("field inversions", "field_ops.invs")),
        _ => None,
    }
}

/// One field axis of a run as a multiple of `√r`, or `None` when the
/// record carries no field-operation block (unknown, not zero) or `axis`
/// is not one of [`FIELD_AXES`].
pub fn field_axis(r: &Record, axis: &str) -> Option<f64> {
    let f = r.field_ops.as_ref()?;
    let count = match axis {
        "field_muls" => f.muls,
        "field_sqrs" => f.sqrs,
        "field_invs" => f.invs,
        _ => return None,
    };
    Some(count as f64 / (r.workload.curve.r as f64).sqrt())
}

/// What `bound fit` reads besides the sessions.
pub struct FitInputs<'a> {
    pub dirs: &'a [PathBuf],
    pub arm: &'a str,
    pub audits: &'a [PathBuf],
    pub label: Option<&'a str>,
    pub notes: &'a str,
    /// Paths in the record are written relative to this root.
    pub root: &'a Path,
    pub options: FitOptions,
}

fn relative(root: &Path, p: &Path) -> String {
    let abs = std::fs::canonicalize(p).unwrap_or_else(|_| p.to_path_buf());
    let root = std::fs::canonicalize(root).unwrap_or_else(|_| root.to_path_buf());
    abs.strip_prefix(&root)
        .map(|r| r.display().to_string())
        .unwrap_or_else(|_| p.display().to_string())
}

fn sha_of_file(p: &Path) -> Result<String, String> {
    let bytes = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
    Ok(sha256_hex(&bytes))
}

/// Fit the bound of `arm` from the sessions in `dirs`.
pub fn bound_from_sessions(inputs: &FitInputs) -> Result<Bound, String> {
    if inputs.dirs.is_empty() {
        return Err("a bound needs at least one session".into());
    }
    let opts = &inputs.options;
    let mut sessions = Vec::new();
    let mut recs: Vec<Record> = Vec::new();
    let mut method: Option<MethodRef> = None;
    for dir in inputs.dirs {
        let s = read_session(dir)?;
        let plan: PlanDoc = serde_json::from_str(
            &std::fs::read_to_string(dir.join("plan.json"))
                .map_err(|e| format!("{}: {e}", dir.join("plan.json").display()))?,
        )
        .map_err(|e| format!("plan.json: {e}"))?;
        let arm = plan
            .arms
            .iter()
            .find(|a| a.name == inputs.arm)
            .ok_or_else(|| {
                format!(
                    "session {} has no arm `{}` (arms: {})",
                    s.session_id,
                    inputs.arm,
                    plan.arms
                        .iter()
                        .map(|a| a.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            })?;
        let m = MethodRef {
            id: arm.method.id.clone(),
            params: arm.method.params.clone(),
            method_id: arm.method.method_id.clone(),
            family: arm.method.family.clone(),
            entry: arm.method.entry.clone(),
        };
        match &method {
            None => method = Some(m),
            Some(prev) if prev.method_id != m.method_id => {
                return Err(format!(
                    "arm `{}` is {} in one session and {} in another; a bound is about one method",
                    inputs.arm, prev.method_id, m.method_id
                ))
            }
            _ => {}
        }
        let role = format!("{:?}", arm.role).to_lowercase();
        sessions.push(SessionRef {
            dir: relative(inputs.root, dir),
            session_id: s.session_id.clone(),
            spec_id: s.spec_id.clone(),
            status: s.status.clone(),
            arm: inputs.arm.to_string(),
            role,
            binary_sha256: s.binary_sha256.clone(),
            git_commit: s.git_commit.clone(),
            env_class_id: s.env_class_id.clone(),
            records_sha256: sha_of_file(&dir.join("records.jsonl"))?,
        });
        recs.extend(
            read_records(dir)?
                .into_iter()
                .filter(|r| r.arm == inputs.arm && !r.warmup),
        );
    }
    let method = method.expect("at least one session");

    // Filters: named curves, then tier.
    if !opts.curves.is_empty() {
        recs.retain(|r| opts.curves.contains(&r.workload.curve.slug));
    }
    if let Some(t) = &opts.tier {
        recs.retain(|r| field_size(r).map(tier_of_bits) == Some(t.as_str()));
    }
    if recs.is_empty() {
        return Err(format!(
            "no measured runs of arm `{}` after filtering",
            inputs.arm
        ));
    }

    // One domain.
    let families: BTreeSet<&str> = recs
        .iter()
        .map(|r| r.workload.curve.family.as_str())
        .collect();
    if families.len() != 1 {
        return Err(format!(
            "the arm ran on {} curve families ({}); a bound is about one, so pass --curve to pick",
            families.len(),
            families.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    let kinds: BTreeSet<String> = recs
        .iter()
        .map(|r| format!("{:?}", r.workload.kind()).to_lowercase())
        .collect();
    if kinds.len() != 1 {
        return Err(
            "the arm ran on planted and public targets; a bound is about one target kind".into(),
        );
    }
    let units: BTreeSet<&str> = recs.iter().map(|r| r.cost.unit.as_str()).collect();
    if units.len() != 1 {
        return Err("the runs report different units".into());
    }
    let mut tiers: BTreeSet<&str> = BTreeSet::new();
    for r in &recs {
        let bits = field_size(r)
            .ok_or_else(|| format!("{}: curve facts carry no field size", r.record_id))?;
        tiers.insert(tier_of_bits(bits));
    }
    if tiers.len() != 1 {
        return Err(format!(
            "the sizes span tiers {}; tiers are not fungible, so fit one bound per tier with --tier",
            tiers.into_iter().collect::<Vec<_>>().join(" and ")
        ));
    }
    let domain = Domain {
        problem: PROBLEM.into(),
        family: families.into_iter().next().unwrap().to_string(),
        target_kind: kinds.into_iter().next().unwrap(),
        unit: units.into_iter().next().unwrap().to_string(),
        tier: tiers.into_iter().next().unwrap().to_string(),
        envelope: Envelope::default(),
    };
    let domain_id = domain.id()?;

    // Sizes.
    let mut by_slug: BTreeMap<String, Vec<&Record>> = BTreeMap::new();
    for r in &recs {
        by_slug
            .entry(r.workload.curve.slug.clone())
            .or_default()
            .push(r);
    }
    let mut sizes: Vec<SizeRow> = Vec::new();
    for (slug, rs) in &by_slug {
        let any = rs[0];
        let ok: Vec<&&Record> = rs.iter().filter(|r| r.counts()).collect();
        let a = any.workload.curve.automorphisms_available;
        let floor = floor_s(a);
        let mut by_w: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
        for r in &ok {
            if let Some(s) = r.cost.s {
                by_w.entry(r.workload.workload_id.as_str())
                    .or_default()
                    .push(s);
            }
        }
        let strata: Vec<Vec<f64>> = by_w.into_values().collect();
        let s_est = estimate(&strata, opts.resamples, opts.seed);
        let mut levels = BTreeMap::new();
        for r in rs {
            *levels
                .entry(r.isolation.level.name().to_string())
                .or_insert(0) += 1;
        }
        let gaes: Vec<f64> = ok.iter().filter_map(|r| r.cost.total_gae).collect();
        let mems: Vec<f64> = ok.iter().filter_map(|r| memory_entries(r)).collect();
        let unch: Vec<f64> = ok.iter().map(|r| uncharged_work(r)).collect();
        let declared = expected_s(&method.id, a);
        sizes.push(SizeRow {
            slug: slug.clone(),
            log2_r: (any.workload.curve.r as f64).log2(),
            r: any.workload.curve.r,
            field_bits: field_size(any).unwrap_or(0),
            automorphisms_available: a,
            floor_s: floor,
            workloads: rs
                .iter()
                .map(|r| r.workload.workload_id.as_str())
                .collect::<BTreeSet<_>>()
                .len() as u64,
            runs: rs.len() as u64,
            verified: ok.len() as u64,
            mean_gae: mean(&gaes),
            mean_s: s_est.as_ref().map(|e| e.value),
            s_ci95: s_est.as_ref().and_then(|e| e.ci95),
            ratio_to_floor: s_est.as_ref().map(|e| e.value / floor),
            declared_s: declared,
            declared_ratio_to_floor: declared.map(|d| d / floor),
            memory_entries_per_sqrt_r: if mems.len() == ok.len() {
                mean(&mems)
            } else {
                None
            },
            uncharged_per_sqrt_r: mean(&unch),
            levels,
        });
    }
    sizes.sort_by_key(|a| a.r);

    // The fit, over verified runs, one stratum per size.
    let verified: Vec<&Record> = recs.iter().filter(|r| r.counts()).collect();
    let mut strata: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut s_strata: Vec<Vec<f64>> = Vec::new();
    let mut f_strata: Vec<Vec<f64>> = Vec::new();
    let mut m_strata: Vec<Vec<f64>> = Vec::new();
    let mut u_strata: Vec<Vec<f64>> = Vec::new();
    let mut memory_known = true;
    for row in &sizes {
        let mine: Vec<&&Record> = verified
            .iter()
            .filter(|r| r.workload.curve.slug == row.slug)
            .collect();
        let ln_r = (row.r as f64).ln();
        strata.push(
            mine.iter()
                .filter_map(|r| r.cost.total_gae)
                .filter(|g| *g > 0.0)
                .map(|g| (ln_r, g.ln()))
                .collect(),
        );
        s_strata.push(mine.iter().filter_map(|r| r.cost.s).collect());
        f_strata.push(
            mine.iter()
                .filter_map(|r| r.cost.s)
                .map(|s| s / row.floor_s)
                .collect(),
        );
        let mems: Vec<f64> = mine.iter().filter_map(|r| memory_entries(r)).collect();
        if mems.len() != mine.len() {
            memory_known = false;
        }
        m_strata.push(mems);
        u_strata.push(mine.iter().map(|r| uncharged_work(r)).collect());
    }
    let n_sizes = sizes.iter().filter(|s| s.verified > 0).count();
    let (alpha, alpha_ci, log2_c, r2, ci_method) = power_law(&strata, opts.resamples, opts.seed);
    let declared_alpha = expected_s(&method.id, sizes[0].automorphisms_available).map(|_| 0.5);
    let scaling_claim = n_sizes >= MIN_SIZES_FOR_SCALING && alpha.is_some();
    let fit = Fit {
        size_parameter: "r".into(),
        model: "gae = C * r^alpha".into(),
        points: strata.iter().map(|s| s.len() as u64).sum(),
        sizes: n_sizes as u64,
        alpha,
        alpha_ci95: alpha_ci,
        log2_c,
        r_squared: r2,
        ci_method,
        declared_alpha,
        alpha_agrees_with_declared: match (declared_alpha, alpha_ci) {
            (Some(d), Some((lo, hi))) => Some(lo <= d && d <= hi),
            _ => None,
        },
        scaling_claim,
    };
    let none = || Estimate {
        value: f64::NAN,
        ci95: None,
        ci_method: "none".into(),
    };
    let s_est = estimate(&s_strata, opts.resamples, opts.seed ^ 0x5).unwrap_or_else(none);
    let f_est = estimate(&f_strata, opts.resamples, opts.seed ^ 0xF).unwrap_or_else(none);
    let declared_ratios: Vec<f64> = sizes
        .iter()
        .filter_map(|s| s.declared_ratio_to_floor)
        .collect();
    let declared_ratio = if declared_ratios.len() == sizes.len()
        && declared_ratios
            .iter()
            .all(|d| (d - declared_ratios[0]).abs() < 1e-9)
    {
        declared_ratios.first().copied()
    } else {
        None
    };
    let constant = Constant {
        s: s_est,
        ratio_to_floor: f_est.clone(),
        declared_ratio_to_floor: declared_ratio,
    };

    let mut dimensions = BTreeMap::new();
    dimensions.insert(
        "ops".to_string(),
        Dimension {
            statistic: "mean S / sqrt(pi / 2A) over verified runs".into(),
            source: "cost.total_gae, workload.curve.r, automorphisms_available".into(),
            lower_is_better: true,
            known: f_est.value.is_finite(),
            value: f_est.value.is_finite().then_some(f_est.value),
            ci95: f_est.ci95,
        },
    );
    let m_est = if memory_known {
        estimate(&m_strata, opts.resamples, opts.seed ^ 0x3)
    } else {
        None
    };
    dimensions.insert(
        "memory".to_string(),
        Dimension {
            statistic: "mean table entries / sqrt(r) over verified runs".into(),
            source: format!("first of counters.{}", MEMORY_COUNTERS.join(", counters.")),
            lower_is_better: true,
            known: m_est.is_some(),
            value: m_est.as_ref().map(|e| e.value),
            ci95: m_est.as_ref().and_then(|e| e.ci95),
        },
    );
    let u_est = estimate(&u_strata, opts.resamples, opts.seed ^ 0x7);
    dimensions.insert(
        "uncharged".to_string(),
        Dimension {
            statistic: "mean sum of *_uncharged counters / sqrt(r) over verified runs".into(),
            source: "counters.*_uncharged".into(),
            lower_is_better: true,
            known: u_est.is_some(),
            value: u_est.as_ref().map(|e| e.value),
            ci95: u_est.as_ref().and_then(|e| e.ci95),
        },
    );
    // The field-operation axes, with the same two-stage interval, only
    // when there are verified runs and every one of them carries the
    // block.  Otherwise nothing is added — not an unknown axis — so a
    // bound fitted from records written before the block existed is
    // byte for byte the bound it was.
    for (i, axis) in FIELD_AXES.iter().enumerate() {
        let mut strata: Vec<Vec<f64>> = Vec::new();
        let mut known = !verified.is_empty();
        for row in &sizes {
            let mine: Vec<&&Record> = verified
                .iter()
                .filter(|r| r.workload.curve.slug == row.slug)
                .collect();
            let vals: Vec<f64> = mine.iter().filter_map(|r| field_axis(r, axis)).collect();
            if vals.len() != mine.len() {
                known = false;
            }
            strata.push(vals);
        }
        if !known {
            continue;
        }
        let (what, source) = field_axis_describes(axis).expect("a field axis");
        let est = estimate(&strata, opts.resamples, opts.seed ^ (0x11 + 2 * i as u64));
        dimensions.insert(
            axis.to_string(),
            Dimension {
                statistic: format!("mean {what} / sqrt(r) over verified runs"),
                source: format!("{source}, workload.curve.r"),
                lower_is_better: true,
                known: est.is_some(),
                value: est.as_ref().map(|e| e.value),
                ci95: est.as_ref().and_then(|e| e.ci95),
            },
        );
    }

    // Stages.
    let mut names: Vec<String> = Vec::new();
    for r in &verified {
        for p in &r.phases {
            if !names.contains(&p.name) {
                names.push(p.name.clone());
            }
        }
    }
    let total_gae: f64 = verified.iter().filter_map(|r| r.cost.total_gae).sum();
    let mut stages = Vec::new();
    for name in names {
        let mut sum = 0.0;
        let mut st: Vec<Vec<(f64, f64)>> = Vec::new();
        for row in &sizes {
            let mut pts = Vec::new();
            for r in verified
                .iter()
                .filter(|r| r.workload.curve.slug == row.slug)
            {
                for p in r.phases.iter().filter(|p| p.name == name) {
                    sum += p.gae;
                    if p.gae > 0.0 {
                        pts.push(((row.r as f64).ln(), p.gae.ln()));
                    }
                }
            }
            st.push(pts);
        }
        let with_work = st.iter().filter(|s| !s.is_empty()).count();
        let (a, ci) = if with_work >= MIN_SIZES_FOR_SCALING {
            let (a, ci, _, _, _) = power_law(&st, opts.resamples, opts.seed ^ 0x51);
            (a, ci)
        } else {
            (None, None)
        };
        stages.push(StageRow {
            name,
            mean_gae: if verified.is_empty() {
                0.0
            } else {
                sum / verified.len() as f64
            },
            share_of_gae: if total_gae > 0.0 {
                sum / total_gae
            } else {
                0.0
            },
            sizes_with_work: with_work as u64,
            alpha: a,
            alpha_ci95: ci,
        });
    }

    // Provenance and admissibility.
    let mut audits = Vec::new();
    let session_ids: BTreeSet<&str> = sessions.iter().map(|s| s.session_id.as_str()).collect();
    for p in inputs.audits {
        let bytes = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
        let rep: AuditReport =
            serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", p.display()))?;
        if !session_ids.contains(rep.session_id.as_str()) {
            return Err(format!(
                "{} audits session {}, which this bound does not read",
                p.display(),
                rep.session_id
            ));
        }
        audits.push(AuditRef {
            path: relative(inputs.root, p),
            sha256: sha256_hex(&bytes),
            session_id: rep.session_id.clone(),
            ok: rep.ok,
            replays: rep.replays.len() as u64,
            replays_reproduced: rep.replays.iter().filter(|r| r.reproduced).count() as u64,
        });
    }
    let mut statuses = BTreeMap::new();
    for r in &recs {
        *statuses.entry(r.outcome.status.clone()).or_insert(0u64) += 1;
    }
    let mut reasons = Vec::new();
    let unverified = recs.len() - verified.len();
    if unverified > 0 {
        reasons.push(format!(
            "{unverified} of {} measured runs did not verify; only a wholly verified arm is admissible",
            recs.len()
        ));
    }
    if !scaling_claim {
        reasons.push(format!(
            "{n_sizes} size(s): below {MIN_SIZES_FOR_SCALING}, so the exponent is descriptive and this is a constant bound"
        ));
    }
    for s in &sessions {
        if s.status != "complete" {
            reasons.push(format!("session {} is {}", s.session_id, s.status));
        }
    }
    if audits.iter().any(|a| !a.ok) {
        reasons.push("an audit receipt reports problems".into());
    }
    let bounded = recs.iter().any(|r| r.cost.lower_bound);
    let mut unpriced: Vec<String> = recs.iter().flat_map(|r| r.cost.unpriced.clone()).collect();
    unpriced.sort();
    unpriced.dedup();
    let admissible = unverified == 0 && audits.iter().all(|a| a.ok);
    let admissibility = Admissibility {
        status: if admissible {
            "admissible"
        } else {
            "inadmissible"
        }
        .into(),
        bounded,
        unpriced,
        deterministic: recs.iter().all(|r| r.cost.deterministic),
        reasons,
    };
    let label = inputs.label.map(String::from).unwrap_or_else(|| {
        format!(
            "{} on {} {} curves, {} tier",
            method.id, n_sizes, domain.family, domain.tier
        )
    });
    let mut b = Bound {
        schema: BOUND_SCHEMA.into(),
        bound_id: String::new(),
        label,
        domain,
        domain_id,
        method,
        level: if scaling_claim {
            "exponent"
        } else {
            "constant"
        }
        .into(),
        sizes,
        fit,
        constant,
        dimensions,
        stages,
        provenance: Provenance {
            sessions,
            audits,
            records: recs.len() as u64,
            verified: verified.len() as u64,
            statuses,
        },
        admissibility,
        fit_options: opts.clone(),
        improves_on: vec![],
        verdict_id: None,
        notes: inputs.notes.to_string(),
    };
    b.seal()?;
    Ok(b)
}

/// Re-derive a committed bound from the sessions it names and require the
/// same id: the bound's analogue of a replay.  Returns the re-derived
/// bound.
pub fn check(path: &Path, root: &Path) -> Result<Bound, String> {
    let recorded = Bound::read(path)?;
    let dirs: Vec<PathBuf> = recorded
        .provenance
        .sessions
        .iter()
        .map(|s| root.join(&s.dir))
        .collect();
    let audits: Vec<PathBuf> = recorded
        .provenance
        .audits
        .iter()
        .map(|a| root.join(&a.path))
        .collect();
    let arm = recorded
        .provenance
        .sessions
        .first()
        .map(|s| s.arm.clone())
        .ok_or("the bound names no session")?;
    let mut fresh = bound_from_sessions(&FitInputs {
        dirs: &dirs,
        arm: &arm,
        audits: &audits,
        label: Some(&recorded.label),
        notes: &recorded.notes,
        root,
        options: recorded.fit_options.clone(),
    })?;
    // What a verdict wrote onto the bound is not re-derivable from the
    // sessions; carry it over, then compare the whole document.
    fresh.improves_on = recorded.improves_on.clone();
    fresh.verdict_id = recorded.verdict_id.clone();
    fresh.seal()?;
    if fresh.bound_id != recorded.bound_id {
        let mut diffs = Vec::new();
        if fresh.sizes != recorded.sizes {
            diffs.push("sizes");
        }
        if fresh.fit != recorded.fit {
            diffs.push("fit");
        }
        if fresh.constant != recorded.constant {
            diffs.push("constant");
        }
        if fresh.dimensions != recorded.dimensions {
            diffs.push("dimensions");
        }
        if fresh.provenance != recorded.provenance {
            diffs.push("provenance");
        }
        if fresh.admissibility != recorded.admissibility {
            diffs.push("admissibility");
        }
        return Err(format!(
            "{}: re-derived bound is {}, recorded {}; differs in: {}",
            path.display(),
            fresh.bound_id,
            recorded.bound_id,
            if diffs.is_empty() {
                "other fields".to_string()
            } else {
                diffs.join(", ")
            }
        ));
    }
    Ok(fresh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ols_recovers_an_exact_power_law() {
        // cost = 3 · r^0.5 at four sizes.
        let pts: Vec<(f64, f64)> = [16.0f64, 18.0, 20.0, 22.0]
            .iter()
            .map(|b| {
                let r = 2f64.powf(*b);
                (r.ln(), (3.0 * r.sqrt()).ln())
            })
            .collect();
        let (alpha, b, r2) = ols(&pts).unwrap();
        assert!((alpha - 0.5).abs() < 1e-9, "{alpha}");
        assert!((b.exp() - 3.0).abs() < 1e-9);
        assert!((r2 - 1.0).abs() < 1e-12);
        assert!(ols(&pts[..1]).is_none());
        assert!(ols(&[(1.0, 1.0), (1.0, 2.0)]).is_none());
    }

    #[test]
    fn power_law_interval_brackets_the_slope_and_needs_two_sizes() {
        let mut strata = Vec::new();
        for (i, b) in [16.0f64, 18.0, 20.0, 22.0].iter().enumerate() {
            let r = 2f64.powf(*b);
            // Spread the runs around the law, deterministically.
            strata.push(
                (0..8)
                    .map(|k| {
                        let noise = 1.0 + 0.3 * (((k * 7 + i * 3) % 5) as f64 - 2.0) / 2.0;
                        (r.ln(), (1.25 * r.sqrt() * noise).ln())
                    })
                    .collect::<Vec<_>>(),
            );
        }
        let (alpha, ci, log2_c, r2, method) = power_law(&strata, 500, 1);
        let a = alpha.unwrap();
        assert!((a - 0.5).abs() < 0.05, "{a}");
        let (lo, hi) = ci.unwrap();
        assert!(lo <= a && a <= hi);
        assert_eq!(method, "cluster");
        assert!(log2_c.unwrap().is_finite());
        assert!(r2.unwrap() > 0.9);
        let (_, ci1, _, _, m1) = power_law(&strata[..1], 500, 1);
        assert!(ci1.is_none());
        assert_eq!(m1, "none");
    }

    #[test]
    fn tiers_follow_the_field_size() {
        assert_eq!(tier_of_bits(13), "toy");
        assert_eq!(tier_of_bits(32), "toy");
        assert_eq!(tier_of_bits(33), "medium");
        assert_eq!(tier_of_bits(96), "medium");
        assert_eq!(tier_of_bits(131), "crypto");
    }

    #[test]
    fn the_domain_id_ignores_nothing_and_is_stable() {
        let d = Domain {
            problem: PROBLEM.into(),
            family: "prime".into(),
            target_kind: "planted".into(),
            unit: "ecbench.gae".into(),
            tier: "toy".into(),
            envelope: Envelope::default(),
        };
        let id = d.id().unwrap();
        assert!(id.starts_with("ECDOM1h"), "{id}");
        let mut k = d.clone();
        k.family = "koblitz".into();
        assert_ne!(k.id().unwrap(), id);
    }

    #[derive(Serialize, Deserialize)]
    struct Doc {
        schema: String,
        doc_id: String,
        value: f64,
        inner: BTreeMap<String, String>,
    }

    #[test]
    fn a_sealed_document_checks_and_a_changed_byte_does_not() {
        let mut inner = BTreeMap::new();
        inner.insert(
            "doc_id".to_string(),
            "nested keys of the same name are fine".to_string(),
        );
        let d = Doc {
            schema: "t/v1".into(),
            doc_id: String::new(),
            value: 0.8862269254527579,
            inner,
        };
        let (id, text) = seal_document(&d, "doc_id", "T1h").unwrap();
        assert!(id.starts_with("T1h") && id.len() == 3 + 12);
        assert_eq!(check_document_seal(&text, "doc_id", "T1h").unwrap(), id);
        let forged = text.replace("0.8862269254527579", "0.5862269254527579");
        assert!(check_document_seal(&forged, "doc_id", "T1h").is_err());
        // The nested `doc_id` is not the sealed field: the top-level one is
        // the first occurrence in the pretty text because keys sort.
        assert!(text.find("\"doc_id\": \"T1h").unwrap() < text.find("nested keys").unwrap());
    }
}
