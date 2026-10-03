//! `vs_rho` claims: an IC row and a rho row of one session, assembled into
//! the report `docs/ic/boundary_targets.json` defines and checked against
//! it natively.
//!
//! A claim names its candidate by the tournament's IC1 identity
//! (`research/ic_candidate_tournament_20260915/identity.py`), which this
//! module ports: the curve record and its EC1 id, the factor-base
//! inventory computed from the base's own points (cofactor projection,
//! sign-and-Frobenius orbits), the method record with every stage
//! resolved and bound to source hashes, and the label
//! `IC1N…C…fb…PDP…RC…LA…TD…ISO0h<12 hex>`.  The workload record is
//! `identity.py`'s too, and like it this refuses a planted target: a claim
//! is about a previously unseen public target.  All hashing is
//! `identity.py`'s canonical JSON (sorted keys, compact, UTF-8, no
//! floats), pinned against values computed by `identity.py` in the tests.
//!
//! The candidate's implementation is the source of the IC stages, hashed
//! into this binary at compile time.  A claim is therefore built only by
//! the binary that measured the session (`binary_sha256` must match), so
//! the identity names the code that ran.

use std::collections::BTreeSet;
use std::path::Path;

use serde_json::{json, Map, Value};

use crate::cryptanalysis::ecbench::audit::{audit, AuditReport};
use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ecbench::methods::{dump_factor_base, FactorBaseFacts, ResolvedMethod};
use crate::cryptanalysis::ecbench::record::Record;
use crate::cryptanalysis::ecbench::runner::{read_records, read_session, PlanDoc, Session};
use crate::cryptanalysis::ecbench::spec::Spec;
use crate::cryptanalysis::ecbench::workload::{Instance, PUBLIC_TARGET_LAW};
use crate::cryptanalysis::ic_boundary::{BinaryGroup, BinaryInstance, CountedGroup, GroupOps};
use crate::cryptanalysis::koblitz_fast::FastPoint;

pub const CLAIM_LEDGER: &str = include_str!("../../../docs/ic/boundary_targets.json");

/// The `schema` field of a claim file, by which `ecbench db` knows it.
/// The report's shape is the ledger's `vs_rho` schema; this names the
/// producer.
pub const CLAIM_SCHEMA: &str = "ecbench.claim/v1";

// ── identity.py's canonical JSON ───────────────────────────────────

/// `json.dumps(v, sort_keys=True, separators=(',', ':'), ensure_ascii=False,
/// allow_nan=False)`: like [`super::canonical::canonical`] but with
/// non-ASCII characters written as UTF-8, as `identity.py` hashes them.
pub fn canonical_utf8(v: &Value) -> Result<String, String> {
    fn write(v: &Value, out: &mut String) -> Result<(), String> {
        match v {
            Value::Null => out.push_str("null"),
            Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Value::Number(n) => {
                if n.is_f64() {
                    return Err(format!("float {n} in an identity record"));
                }
                out.push_str(&n.to_string());
            }
            Value::String(s) => {
                out.push('"');
                for c in s.chars() {
                    match c {
                        '"' => out.push_str("\\\""),
                        '\\' => out.push_str("\\\\"),
                        '\n' => out.push_str("\\n"),
                        '\r' => out.push_str("\\r"),
                        '\t' => out.push_str("\\t"),
                        '\u{08}' => out.push_str("\\b"),
                        '\u{0c}' => out.push_str("\\f"),
                        c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
                        c => out.push(c),
                    }
                }
                out.push('"');
            }
            Value::Array(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    write(item, out)?;
                }
                out.push(']');
            }
            Value::Object(map) => {
                let mut keys: Vec<&String> = map.keys().collect();
                keys.sort();
                out.push('{');
                for (i, k) in keys.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    write(&Value::String((*k).clone()), out)?;
                    out.push(':');
                    write(&map[*k], out)?;
                }
                out.push('}');
            }
        }
        Ok(())
    }
    let mut out = String::new();
    write(v, &mut out)?;
    Ok(out)
}

/// `identity.sha256`: SHA-256 of the canonical bytes.
pub fn id_sha256(v: &Value) -> Result<String, String> {
    Ok(sha256_hex(canonical_utf8(v)?.as_bytes()))
}

// ── The curve record ───────────────────────────────────────────────

/// The Koblitz facts `identity.py`'s `Curve` checks, from an instance.
struct Kb<'a> {
    inst: &'a BinaryInstance,
    n: u32,
    a: u64,
    modulus: u64,
    r: u64,
    h: u64,
}

impl<'a> Kb<'a> {
    fn of(inst: &'a BinaryInstance) -> Result<Self, String> {
        let kc = inst
            .koblitz
            .as_ref()
            .ok_or("a claim's IC1 identity covers Koblitz curves only")?;
        let n = inst.n;
        // identity.py's adapter (oracle.Curve) admits odd degrees 5..=61.
        if !(5..=61).contains(&n) || n.is_multiple_of(2) {
            return Err(format!(
                "the IC1 identity admits odd field degrees 5 to 61, not {n}"
            ));
        }
        let modulus = (1u64 << n)
            | inst
                .irreducible
                .low_terms
                .iter()
                .fold(0u64, |m, &t| m | (1u64 << t));
        let _ = kc;
        Ok(Self {
            inst,
            n,
            a: inst.a,
            modulus,
            r: inst.r,
            h: inst.cofactor,
        })
    }

    fn order(&self) -> u64 {
        self.h * self.r
    }

    fn frob(&self, p: FastPoint) -> FastPoint {
        if p.infinity {
            p
        } else {
            FastPoint::affine(self.inst.gf.sqr(p.x), self.inst.gf.sqr(p.y))
        }
    }

    fn neg(&self, p: FastPoint) -> FastPoint {
        if p.infinity {
            p
        } else {
            FastPoint::affine(p.x, p.x ^ p.y)
        }
    }

    /// Remove the signed Frobenius orbit of `start` from `set`, including
    /// the part of it `set` does not hold.
    fn remove_orbit(&self, set: &mut BTreeSet<(u64, u64)>, start: (u64, u64)) {
        let mut p = FastPoint::affine(start.0, start.1);
        for _ in 0..self.n {
            set.remove(&(p.x, p.y));
            let q = self.neg(p);
            set.remove(&(q.x, q.y));
            p = self.frob(p);
        }
    }
}

/// `identity.curve_record`.
fn curve_record(kb: &Kb) -> Result<Value, String> {
    let g = kb.inst.generator;
    let trace = (1i128 << kb.n) + 1 - kb.order() as i128;
    let mut record = json!({
        "field": {"p": 2, "n": kb.n, "basis": "polynomial", "modulus": kb.modulus,
                  "element_encoding": "unsigned-integer-polynomial-bits"},
        "curve": {"tag": format!("kb{}", kb.a), "model": "y^2+x*y=x^3+a2*x^2+a6",
                  "coefficients": {"a1": 1, "a2": kb.a, "a3": 0, "a4": 0, "a6": 1},
                  "order": kb.order(), "trace": trace as i64,
                  "r": kb.r, "cofactor": kb.h, "generator": [g.x, g.y],
                  "target_group": "prime-order-subgroup-generated-by-G"},
    });
    let id = format!("EC1N{}Ckb{}h{}", kb.n, kb.a, &id_sha256(&record)?[..12]);
    record["curve"]["curve_id"] = json!(id);
    Ok(record)
}

// ── The factor-base inventory ──────────────────────────────────────

fn encode(points: &[(u64, u64)]) -> Value {
    Value::Array(points.iter().map(|(x, y)| json!([x, y])).collect())
}

/// How a run's columns relate to the usable base.  ecbench's Koblitz
/// bases fold the *raw* lifted points by sign and Frobenius, so a raw
/// orbit whose cofactor image is the identity (the 2-torsion point
/// `(0, 1)` above `x = 0`, which every subspace contains) is a column
/// whose logarithm is identically zero, and two raw orbits with one
/// image orbit carry one unknown twice.  `identity.py` counts neither;
/// the claim discloses both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ColumnMap {
    /// The signed Frobenius orbits of the raw points: the run's columns.
    raw: u64,
    /// Raw orbits whose cofactor image is the identity.
    identity_image: u64,
    /// Raw orbits whose image orbit another raw orbit already covers.
    merged: u64,
}

impl ColumnMap {
    fn json(&self) -> Value {
        json!({
            "fold": "signed Frobenius orbits of the raw lifted points",
            "raw_columns": self.raw,
            "identity_image_columns": self.identity_image,
            "merged_columns": self.merged,
        })
    }
}

/// `identity.factor_base_inventory` under the `cofactor` convention: the
/// usable base is the image `[h]B` without the identity, its columns the
/// sign-and-Frobenius orbits of that image.  `run_columns` is the column
/// count the run reported; it must be the signed Frobenius orbits of the
/// raw points, recounted here, and the inventory's `effective_columns`
/// is the image's orbit count, as `identity.py` requires.
fn factor_base_inventory(
    kb: &Kb,
    points: &[(u64, u64)],
    run_columns: u64,
) -> Result<(Value, ColumnMap), String> {
    let set: BTreeSet<(u64, u64)> = points.iter().copied().collect();
    if points.is_empty() || set.len() != points.len() {
        return Err("empty or duplicate geometric base".into());
    }
    let g = BinaryGroup(&kb.inst.fast);
    let mut ops = GroupOps::default();
    let mut image = |(x, y): (u64, u64)| {
        let p = g.mul(&mut ops, FastPoint::affine(x, y), kb.h);
        (!p.infinity).then_some((p.x, p.y))
    };
    let projected: Vec<Option<(u64, u64)>> = points.iter().map(|&p| image(p)).collect();
    // The run's fold, recounted on the raw points.  Projection commutes
    // with sign and Frobenius, so one image per raw orbit decides it.
    let (raw, identity_image) = {
        let (mut raw, mut dead) = (0u64, 0u64);
        let mut remaining = set.clone();
        while let Some(&start) = remaining.iter().next() {
            kb.remove_orbit(&mut remaining, start);
            raw += 1;
            dead += u64::from(image(start).is_none());
        }
        (raw, dead)
    };
    if run_columns != raw {
        return Err(format!(
            "the run reports {run_columns} columns, but its points have {raw} signed Frobenius \
             orbits: the IC1 identity's quotient is sign-and-Frobenius, so this base's fold \
             cannot be named by it"
        ));
    }
    let identity_images = projected.iter().filter(|p| p.is_none()).count();
    let usable: BTreeSet<(u64, u64)> = projected.iter().flatten().copied().collect();
    if usable.is_empty() {
        return Err("base has no usable points".into());
    }
    // subgroup_orbits: partition by Frobenius and negation; check one
    // point of each orbit lies in the r-subgroup.
    let mut remaining = usable.clone();
    let mut orbits = 0u64;
    while let Some(&start) = remaining.iter().next() {
        if !g
            .mul(&mut ops, FastPoint::affine(start.0, start.1), kb.r)
            .infinity
        {
            return Err("non-subgroup factor-base point".into());
        }
        kb.remove_orbit(&mut remaining, start);
        orbits += 1;
    }
    let map = ColumnMap {
        raw,
        identity_image,
        merged: raw
            .checked_sub(identity_image + orbits)
            .ok_or("more image orbits than raw orbits")?,
    };
    let sorted: Vec<(u64, u64)> = set.iter().copied().collect();
    let usable_sorted: Vec<(u64, u64)> = usable.iter().copied().collect();
    let inventory = json!({
        "geometric_point_count": points.len(),
        "geometric_set_sha256": id_sha256(&encode(&sorted))?,
        "geometric_order_sha256": id_sha256(&encode(points))?,
        "usable_point_count": usable.len(),
        "usable_set_sha256": id_sha256(&encode(&usable_sorted))?,
        "subgroup_map": "cofactor-multiplication",
        "identity_images": identity_images,
        "duplicate_nonidentity_images": projected.len() - identity_images - usable.len(),
        "quotient": "sign-and-Frobenius",
        "effective_columns": orbits,
    });
    Ok((inventory, map))
}

// ── The candidate record ───────────────────────────────────────────

/// The source files the IC stages run, hashed at compile time.
pub struct Implementation {
    pub components: Vec<(String, String)>,
}

impl Implementation {
    /// This binary's IC sources.
    pub fn this_binary() -> Self {
        let h = |bytes: &[u8]| sha256_hex(bytes);
        Self {
            components: vec![
                (
                    "ic_pipeline".into(),
                    h(include_bytes!("../ic_framework/mod.rs")),
                ),
                (
                    "ic_plugins".into(),
                    h(include_bytes!("../ic_framework/plugins.rs")),
                ),
                (
                    "ic_solvers".into(),
                    h(include_bytes!("../ic_framework/solvers.rs")),
                ),
                (
                    "ic_linalg".into(),
                    h(include_bytes!("../ic_framework/linalg.rs")),
                ),
                (
                    "relation_loop".into(),
                    h(include_bytes!("../ic_boundary.rs")),
                ),
                ("harness".into(), h(include_bytes!("methods.rs"))),
            ],
        }
    }

    fn sha(&self, role: &str) -> String {
        self.components
            .iter()
            .find(|(r, _)| r == role)
            .map(|(_, s)| s.clone())
            .unwrap_or_default()
    }

    fn json(&self) -> Value {
        Value::Array(
            self.components
                .iter()
                .map(|(role, sha)| json!({"role": role, "sha256": sha}))
                .collect(),
        )
    }
}

/// A compact stage code: lower-case letters and digits only.
fn code(s: &str) -> String {
    s.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .flat_map(|c| c.to_lowercase())
        .collect()
}

/// The method record of an `ic.pipeline` configuration, every stage
/// resolved as `identity.candidate_record` requires.
fn method_record(
    m: &ResolvedMethod,
    kb: &Kb,
    fb: &FactorBaseFacts,
    columns: &ColumnMap,
    imp: &Implementation,
) -> Result<Value, String> {
    let p = &m.params;
    let (oracle, oparams) = p
        .get("oracle")
        .map(|s| s.split_once(':').unwrap_or((s.as_str(), "")))
        .ok_or("ic.pipeline without an oracle")?;
    let summands: u64 = oparams
        .split(',')
        .find_map(|kv| kv.strip_prefix("m=")?.parse().ok())
        .unwrap_or(2);
    let solver = p.get("solver").cloned().unwrap_or_default();
    let solver_code = if solver.is_empty() {
        code(oracle)
    } else {
        code(&format!(
            "{oracle}{}",
            solver.split(':').next().unwrap_or("")
        ))
    };
    let collector = code(p.get("targets").map(String::as_str).unwrap_or("walk"));
    let la = code(
        p.get("linalg")
            .map(String::as_str)
            .unwrap_or("incremental-gauss"),
    );
    let algebraic = oracle == "descent-algebraic";
    let manifest = id_sha256(&imp.json())?;
    Ok(json!({
        "isogeny": "none",
        "endomorphism": {"order_conductor": null, "frobenius_order_conductor": null, "volcano_levels": []},
        "factor_base": {
            "construction": {
                "family": fb.family,
                "params": fb.params,
                "builder": "ic_framework::plugins",
                "points_sha256": fb.points_sha256,
                "solver_columns": columns.json(),
            },
            "nominal_bound": fb.signed_points,
        },
        "point_decomposition": {
            "summands": summands,
            "solver": solver_code,
            "summation_polynomial": if algebraic { "Semaev summation polynomial, Weil-descended over the base subspace" } else { "none: exact table lookup of point sums" },
            "encoding": if algebraic { format!("boolean system for solver {solver}") } else { format!("pair table ({oracle})") },
            "equation_order": if algebraic { "solver default" } else { "not applicable" },
            "monomial_order": if algebraic { "solver default" } else { "not applicable" },
            "internal_matrix_kernel": if algebraic { "solver default" } else { "not applicable" },
            "limits": format!("max_trials={}, solver_budget_seconds={}", p.get("max_trials").map(String::as_str).unwrap_or("?"), p.get("solver_budget_seconds").map(String::as_str).unwrap_or("0")),
            "cache_policy": "tables built per run (cold)",
            "source_sha256": imp.sha("ic_plugins"),
        },
        "relation_collection": {
            "collector": collector,
            "query_distribution": "[a]G + [b]Q from the r-adding walk or uniform draws of ic_boundary::collect_and_solve_with",
            "query_rule": format!("targets={}", p.get("targets").map(String::as_str).unwrap_or("walk")),
            "filtering": "target guard: a repeated query is skipped",
            "verification": "decompositions are exact by construction; the recovered logarithm is checked by [d]G = Q",
            "duplicates": "guarded by key",
            "dependencies": "incremental elimination drops dependent rows",
            "stop_rule": "stop when the target column is pinned",
            "source_sha256": imp.sha("relation_loop"),
        },
        "relation_linear_algebra": {
            "solver": la,
            "modulus": kb.r,
            "matrix_construction": "one row per relation over the factor-base columns plus the target column",
            "orbit_quotient": "sign-and-Frobenius",
            "rank_criterion": "the target column is pinned",
            "block_parameters": "none",
            "preconditioner": "none",
            "source_sha256": imp.sha("ic_linalg"),
        },
        "target_descent": {
            "method": "joint",
            "policy": "classic index calculus: the target is a column of the relation matrix",
            "recursive_solvers": "none",
            "success_rule": "[d]G = Q",
            "stop_rule": "first pinned target column",
            "source_sha256": imp.sha("ic_pipeline"),
        },
        "implementation": {
            "source_manifest_sha256": manifest,
            "components": imp.json(),
            "flags": {"ecbench_method_id": m.method_id, "params": p},
        },
    }))
}

/// `identity.candidate_manifest`: `(candidate_id, record_sha256, record)`.
fn candidate_manifest(
    kb: &Kb,
    method: Value,
    inventory: Value,
) -> Result<(String, String, Value), String> {
    let mut record = curve_record(kb)?;
    let usable = inventory["usable_point_count"].as_u64().unwrap_or(0);
    let mut method = method;
    method["factor_base"]["inventory"] = inventory;
    let obj = record.as_object_mut().expect("an object");
    obj.insert("schema_version".into(), json!(1));
    for (k, v) in method.as_object().expect("an object") {
        obj.insert(k.clone(), v.clone());
    }
    let sha = id_sha256(&record)?;
    let pd = &record["point_decomposition"];
    let label = format!(
        "IC1N{}C{}fb{}PDP{}{}RC{}LA{}TD{}ISO0h{}",
        kb.n,
        record["curve"]["tag"].as_str().unwrap_or(""),
        usable,
        pd["summands"],
        pd["solver"].as_str().unwrap_or(""),
        record["relation_collection"]["collector"]
            .as_str()
            .unwrap_or(""),
        record["relation_linear_algebra"]["solver"]
            .as_str()
            .unwrap_or(""),
        record["target_descent"]["method"].as_str().unwrap_or(""),
        &sha[..12]
    );
    Ok((label, sha, record))
}

/// `identity.workload_manifest`: `(workload_id, record_sha256, record)`.
fn workload_manifest(
    curve_id: &str,
    target: (u64, u64),
    target_seed: Value,
    algorithm_seed: u64,
    envelope: &Value,
) -> Result<(String, String, Value), String> {
    let record = json!({
        "curve_id": curve_id,
        "targets": [[target.0, target.1]],
        "input_law": format!("ecbench.{PUBLIC_TARGET_LAW}"),
        "target_seeds": [target_seed],
        "algorithm_seed": algorithm_seed,
        "cache_policy": "cold",
        "target_count": 1,
        "resource_envelope": envelope,
    });
    let h = id_sha256(&record)?;
    Ok((h[..12].to_string(), h, record))
}

// ── The claim ──────────────────────────────────────────────────────

/// An independent replay, as a claim cites it.
pub struct Independent {
    pub receipt: AuditReport,
    pub receipt_sha256: String,
    pub pointer: String,
}

/// Why `i` is not an independent replay of `runs` in the session at
/// `dir`: empty when it is one.  It must be a passing audit of this
/// session's exact files, run on another host class, that reproduced
/// every one of the runs.
fn independent_problems(
    i: &Independent,
    dir: &Path,
    session: &Session,
    runs: &[&Record],
) -> Vec<String> {
    let r = &i.receipt;
    let mut problems = Vec::new();
    if !r.ok {
        problems.push(format!("the audit failed ({} problems)", r.problems.len()));
    }
    if r.session_id != session.session_id {
        problems.push(format!(
            "it audits session {}, not {}",
            r.session_id, session.session_id
        ));
    }
    match r.auditor_env_class_id.as_deref() {
        None => problems.push("it records no auditor host class".into()),
        Some(c) if c == session.env_class_id => {
            problems.push(format!("it ran on the session's own host class {c}"))
        }
        Some(_) => {}
    }
    for name in [
        "session.json",
        "spec.json",
        "host.json",
        "plan.json",
        "records.jsonl",
    ] {
        let here = std::fs::read(dir.join(name)).ok().map(|b| sha256_hex(&b));
        if r.files.get(name) != here.as_ref() {
            problems.push(format!("its {name} is not this session's"));
        }
    }
    let replayed: BTreeSet<&str> = r
        .replays
        .iter()
        .filter(|x| x.reproduced)
        .map(|x| x.record_id.as_str())
        .collect();
    for run in runs {
        if !replayed.contains(run.record_id.as_str()) {
            problems.push(format!("it did not reproduce run {}", run.run_id));
        }
    }
    if i.pointer.trim().is_empty() {
        problems.push("no durable pointer to the receipt".into());
    }
    problems
}

fn hex(s: &str) -> Result<u64, String> {
    u64::from_str_radix(s.trim_start_matches("0x"), 16).map_err(|_| format!("bad hex `{s}`"))
}

fn ms(ns: u64) -> f64 {
    ns as f64 / 1e6
}

/// Assemble the `vs_rho` claim for the IC arm and rho arm on one workload
/// and round of the session in `dir` (rounds count from 0, warm-ups
/// included; `None` takes the first measured round both arms ran).
pub fn build(
    dir: &Path,
    ic_arm: &str,
    rho_arm: &str,
    workload_id: &str,
    round: Option<u32>,
    independent: Option<&Independent>,
) -> Result<Value, String> {
    let session = read_session(dir)?;
    let records = read_records(dir)?;
    let read =
        |name: &str| std::fs::read_to_string(dir.join(name)).map_err(|e| format!("{name}: {e}"));
    let plan: PlanDoc =
        serde_json::from_str(&read("plan.json")?).map_err(|e| format!("plan.json: {e}"))?;
    let spec = Spec::from_json(&read("spec.json")?)?;
    if session.status != "complete" {
        return Err(format!(
            "session {} is {}, not complete",
            session.session_id, session.status
        ));
    }
    let local = audit(dir, 0)?;
    if !local.ok {
        return Err(format!(
            "the session does not audit clean: {}",
            local.problems.join("; ")
        ));
    }
    let this = std::env::current_exe()
        .ok()
        .and_then(|p| std::fs::read(p).ok())
        .map(|b| sha256_hex(&b));
    if this.is_none() || this != session.binary_sha256 {
        return Err(
            "a claim is built by the binary that measured the session: the candidate identity binds this binary's sources"
                .into(),
        );
    }
    let measured = |arm: &str, round: u32| {
        records.iter().find(|r| {
            r.arm == arm && !r.warmup && r.round == round && r.workload.workload_id == workload_id
        })
    };
    // Rounds count from 0, warm-ups included; by default the first
    // measured round both arms ran.
    let round = match round {
        Some(n) => n,
        None => records
            .iter()
            .filter(|r| !r.warmup && r.arm == ic_arm && r.workload.workload_id == workload_id)
            .map(|r| r.round)
            .filter(|&n| measured(rho_arm, n).is_some())
            .min()
            .ok_or_else(|| {
                format!("no measured round of `{ic_arm}` and `{rho_arm}` on {workload_id}")
            })?,
    };
    let find = |arm: &str| -> Result<&Record, String> {
        measured(arm, round)
            .ok_or_else(|| format!("no measured run of `{arm}` on {workload_id} in round {round}"))
    };
    let (ic, rho) = (find(ic_arm)?, find(rho_arm)?);
    if ic.method.family != "ic" {
        return Err(format!("`{ic_arm}` is not an index-calculus arm"));
    }
    if rho.method.id != "rho.signed_frobenius_strong" {
        return Err(format!(
            "the rho arm must be the strong single-target reference (rho.signed_frobenius_strong), not {}",
            rho.method.id
        ));
    }
    for r in [ic, rho] {
        if r.outcome.status != "verified" || r.outcome.matches_target != Some(true) {
            return Err(format!("run {} did not verify", r.run_id));
        }
        if r.online.is_none() {
            return Err(format!("run {} has no online window", r.run_id));
        }
    }
    let w = &ic.workload;
    if w.planted.is_some() {
        return Err(
            "the workload's target is planted; a claim needs a public target (target_kind: public)"
                .into(),
        );
    }
    let pw = plan
        .workloads
        .iter()
        .find(|x| x.workload_id == workload_id)
        .ok_or("the workload is not in the plan")?;
    let inst = pw.curve_spec.build()?;
    let Instance::Binary(bi) = &inst else {
        return Err("a claim's IC1 identity covers Koblitz curves only".into());
    };
    let kb = Kb::of(bi)?;
    let curve = curve_record(&kb)?;
    let curve_id = curve["curve"]["curve_id"]
        .as_str()
        .unwrap_or("")
        .to_string();

    // The factor base, rebuilt, must be the one the run used.
    let fbf = ic
        .factor_base
        .as_ref()
        .ok_or("the IC run reports no factor base")?;
    let fb_spec = ic
        .method
        .params
        .get("factor_base")
        .ok_or("ic.pipeline without a factor base")?;
    let dump = dump_factor_base(&pw.curve_spec, fb_spec)?;
    if dump.factor_base.fb_id != fbf.fb_id {
        return Err(format!(
            "the rebuilt factor base {} is not the run's {}",
            dump.factor_base.fb_id, fbf.fb_id
        ));
    }
    let points: Vec<(u64, u64)> = dump
        .points
        .iter()
        .map(|p| Ok((hex(&p.x)?, hex(&p.y)?)))
        .collect::<Result<_, String>>()?;
    let (inventory, columns) = factor_base_inventory(&kb, &points, fbf.columns)?;
    let imp = Implementation::this_binary();
    let method = method_record(&ic.method, &kb, fbf, &columns, &imp)?;
    let (candidate_id, candidate_sha, _) = candidate_manifest(&kb, method, inventory)?;

    // One resource envelope: the session's, identical for both arms.
    let plan_cpu = session.cpu_plan.as_ref();
    let envelope = json!({
        "binary_sha256": session.binary_sha256,
        "env_class_id": session.env_class_id,
        "cpus": plan_cpu.map(|p| json!([p.run_cpu])).unwrap_or(json!("none")),
        "numa_node": plan_cpu.and_then(|p| p.node),
        "threads": 1,
        "timeout_seconds": spec.measurement.timeout_seconds,
    });
    let target = (hex(&w.target[0])?, hex(&w.target[1])?);
    let target_seed = json!({"law": PUBLIC_TARGET_LAW, "seed": w.target_seed,
                             "index": w.target_index, "counter": w.target_counter});
    let algorithm_seed: u64 = ic
        .algorithm_seed
        .parse()
        .map_err(|_| "bad algorithm seed")?;
    let (wid, wsha, _) = workload_manifest(
        &curve_id,
        target,
        target_seed.clone(),
        algorithm_seed,
        &envelope,
    )?;
    let n_run = ic.run_id.rsplit('R').next().unwrap_or("1").to_string();
    let run_id = format!("{candidate_id}W{wid}R{n_run}");

    let (icw, rhow) = (ic.online.as_ref().unwrap(), rho.online.as_ref().unwrap());
    let phase = |k: &str| ms(icw.phases_ns.get(k).copied().unwrap_or(0));
    let ic_ms = ms(icw.wall_ns);
    let rho_ms = ms(rhow.wall_ns);
    let target_hash = id_sha256(&json!([target.0, target.1]))?;
    let required = spec.measurement.isolation_required;
    let levels = (ic.isolation.level, rho.isolation.level);
    let admissible = levels.0 >= required && levels.1 >= required;
    let speedup = rho_ms / ic_ms;
    // A receipt that is not an independent replay of these two runs is
    // refused, not cited: without one the certificates stay null and the
    // checker says so.
    if let Some(i) = independent {
        let problems = independent_problems(i, dir, &session, &[ic, rho]);
        if !problems.is_empty() {
            return Err(format!(
                "the receipt is not an independent replay of these runs: {}",
                problems.join("; ")
            ));
        }
    }
    let (cert, pointer) = match independent {
        Some(i) => (json!(i.receipt_sha256), json!(i.pointer)),
        None => (Value::Null, json!("none: no independent replay supplied")),
    };
    let setup_ns: u64 = ic
        .phases
        .iter()
        .filter(|p| p.name == "factor_base" || p.name == "oracle_setup")
        .filter_map(|p| p.wall_ns)
        .sum();
    let table_entries = rho.counters.get("table_entries").copied().unwrap_or(0);
    let lanes = rho.method.params.get("lanes").cloned().unwrap_or_default();
    let dp_bits = rho
        .method
        .params
        .get("dp_bits")
        .cloned()
        .unwrap_or_default();
    let mut claim = Map::new();
    let mut put = |k: &str, v: Value| {
        claim.insert(k.to_string(), v);
    };
    put("schema", json!(CLAIM_SCHEMA));
    put("n_or_bits", json!(kb.n));
    put("candidate_id", json!(candidate_id));
    put("candidate_manifest_sha256", json!(candidate_sha));
    put("workload_id", json!(wid));
    put("workload_manifest_sha256", json!(wsha));
    put("run_id", json!(run_id));
    put("target_count", json!(1));
    put("ic_target_hash", json!(target_hash));
    put("rho_target_hash", json!(target_hash));
    put("timing_class", json!("single_target_online_wall"));
    put("ic_online_wall_ms", json!(ic_ms));
    put("rho_online_wall_ms", json!(rho_ms));
    put("online_speedup", json!(speedup));
    put(
        "ic_online_phase_ms",
        json!({
            "T_target_query_ms": phase("target_query"),
            "T_target_PDP_ms": phase("target_PDP"),
            "T_target_relation_check_ms": phase("target_relation_check"),
            "T_target_descent_ms": phase("target_descent"),
            "T_target_recovery_check_ms": phase("target_recovery_check"),
        }),
    );
    put(
        "online_interval",
        json!({
            "ic_start_event": icw.start_event,
            "ic_stop_event": icw.stop_event,
            "rho_start_event": rhow.start_event,
            "rho_stop_event": rhow.stop_event,
            "ic_included_stages": icw.included_stages,
            "rho_included_stages": rhow.included_stages,
            "ic_zero_phases": icw.zero_phases,
        }),
    );
    put("same_resource_envelope", json!(true));
    put("independent_validation", json!(independent.is_some()));
    put("ic_replay_certificate_sha256", cert.clone());
    put("rho_replay_certificate_sha256", cert);
    put("ic_resource_envelope", envelope.clone());
    put("rho_resource_envelope", envelope.clone());
    put("ic_scalar_verified", json!(true));
    put("rho_scalar_verified", json!(true));
    put(
        "rho_policy",
        json!({
            "worker_count": 1,
            "distinguished_point_memory_bytes": table_entries * 40,
            "walk_policy": format!("signed-Frobenius lockstep walks with batched inversion (koblitz_strong_rho), lanes={lanes}"),
            "collision_policy": format!("distinguished points, dp_bits={dp_bits}; candidates checked by [d]G = Q"),
        }),
    );
    put(
        "verdict",
        json!(match (admissible, independent.is_some()) {
            (true, true) => format!(
                "online speedup {speedup:.4} (rho / IC) on one public target, both runs at {} or above, replayed independently",
                required.name()
            ),
            (true, false) => format!(
                "online speedup {speedup:.4} (rho / IC) on one public target, both runs at {} or above; not a claim until an independent replay is cited",
                required.name()
            ),
            (false, _) => format!(
                "descriptive only, not admissible: the runs earned {} and {}, below the spec's {}",
                levels.0.name(),
                levels.1.name(),
                required.name()
            ),
        }),
    );
    put(
        "claim_boundary",
        json!(format!(
            "one public target on {}; one-target online windows as AGENTS.md 'IC measurements' defines them; the IC factor base and oracle tables are reusable set-up outside the window ({:.3} ms, reported as ic_reusable_setup_ms); host class {}",
            w.curve.slug,
            ms(setup_ns),
            session.env_class_id
        )),
    );
    put("independent_replay_pointer", pointer);
    put("ic_reusable_setup_ms", json!(ms(setup_ns)));
    put("fixture_hash", json!(wsha));
    put("executable_or_source_hash", json!(session.binary_sha256));
    put("host_id", json!(session.env_class_id));
    put("resource_caps", envelope);
    put(
        "seeds",
        json!({"algorithm_seed": ic.algorithm_seed, "target": target_seed}),
    );
    put(
        "claim_boundary_non_claims",
        json!([
            "no claim about other targets, curves or sizes",
            "no claim about cold end-to-end cost, which the records' S reports",
            "no multi-target or batch claim",
            "no transfer to m = 83 or m = 131",
        ]),
    );
    put(
        "isolation_levels",
        json!({"ic": levels.0.name(), "rho": levels.1.name(), "required": required.name()}),
    );
    put(
        "ecbench_records",
        json!({"ic": ic.record_id, "rho": rho.record_id,
                                 "session": session.session_id}),
    );
    put(
        "independent_replay",
        independent.map_or(Value::Null, |i| {
            json!({
                "receipt_sha256": i.receipt_sha256,
                "auditor_env_class_id": i.receipt.auditor_env_class_id,
                "auditor_binary_sha256": i.receipt.auditor_binary_sha256,
                "replays": i.receipt.replays.len(),
            })
        }),
    );
    Ok(Value::Object(claim))
}

// ── The checker ────────────────────────────────────────────────────

/// `boundary_autolab.FIELD_ALIASES`: the report keys that satisfy an
/// abstract required field.
const FIELD_ALIASES: &[(&str, &[&str])] = &[
    ("n_or_bits", &["n_or_bits", "n", "bits"]),
    (
        "factor_base_size_F",
        &["factor_base_size_F", "factor_base_size", "F", "|F|"],
    ),
    (
        "orbit_count_K",
        &["orbit_count_K", "orbit_columns", "K", "orbit_columns_K"],
    ),
    (
        "dimension_l_or_dim",
        &["dimension_l_or_dim", "dimension", "l", "ell", "dim"],
    ),
    ("m_summands", &["m_summands", "m", "summands"]),
    (
        "ffd_or_degree_of_regularity",
        &[
            "ffd_or_degree_of_regularity",
            "ffd",
            "degree_of_regularity",
            "DoR",
            "dor",
        ],
    ),
    (
        "base_id_or_hash",
        &["base_id_or_hash", "base_id", "base_hash"],
    ),
    (
        "eta_or_coverage_policy",
        &["eta_or_coverage_policy", "eta", "coverage_policy"],
    ),
    (
        "pr_decomposition_or_hit_rate_with_ci",
        &[
            "pr_decomposition_or_hit_rate_with_ci",
            "pr_decomposition",
            "hit_rate",
            "hit_rate_with_ci",
        ],
    ),
    (
        "orbit_columns_K",
        &["orbit_columns_K", "orbit_columns", "K"],
    ),
    (
        "la_wall_ms_or_la_charged_ms",
        &["la_wall_ms_or_la_charged_ms", "la_wall_ms", "la_charged_ms"],
    ),
    (
        "recovered_d_verified",
        &["recovered_d_verified", "recovered_d", "d_verified"],
    ),
    (
        "executable_or_source_hash",
        &[
            "executable_or_source_hash",
            "executable_hash",
            "source_hash",
        ],
    ),
    ("seeds", &["seeds", "seed"]),
    (
        "claim_boundary_non_claims",
        &["claim_boundary_non_claims", "non_claims", "claim_boundary"],
    ),
];

/// `boundary_autolab.ONLINE_REQUIRED_STAGES`.
const ONLINE_REQUIRED_STAGES: &[(&str, &[&str])] = &[
    (
        "ic_included_stages",
        &[
            "target_query",
            "target_PDP",
            "target_relation_check",
            "target_descent",
            "target_recovery_check",
        ],
    ),
    (
        "rho_included_stages",
        &["walk", "collision", "recovery_check"],
    ),
];

/// `boundary_autolab.field_present`: present under some alias, not null,
/// not a blank string.
fn field_present(report: &Value, key: &str) -> bool {
    let aliases = FIELD_ALIASES
        .iter()
        .find(|(k, _)| *k == key)
        .map(|(_, a)| *a)
        .unwrap_or(std::slice::from_ref(&key));
    aliases.iter().any(|a| match report.get(*a) {
        None | Some(Value::Null) => false,
        Some(Value::String(s)) => !s.trim().is_empty(),
        Some(_) => true,
    })
}

/// Python's `math.isclose`.
fn isclose(a: f64, b: f64, rel_tol: f64, abs_tol: f64) -> bool {
    a == b || (a - b).abs() <= (rel_tol * a.abs().max(b.abs())).max(abs_tol)
}

/// `boundary_autolab.positive_cost`: a finite positive number (not a bool).
fn positive_cost(v: Option<&Value>) -> Option<f64> {
    v.and_then(Value::as_f64)
        .filter(|x| x.is_finite() && *x > 0.0)
}

fn is_lower_hex(s: &str) -> bool {
    s.chars().all(|c| matches!(c, '0'..='9' | 'a'..='f'))
}

/// A nonempty, non-blank string.
fn text(v: Option<&Value>) -> Option<&str> {
    v.and_then(Value::as_str).filter(|s| !s.trim().is_empty())
}

/// `boundary_autolab.validate_claim(stage="vs_rho")`, natively, with its
/// output: `status` is `PASS` exactly when the three lists are empty.
pub fn check_vs_rho(report: &Value) -> Result<Value, String> {
    let ledger: Value = serde_json::from_str(CLAIM_LEDGER).map_err(|e| e.to_string())?;
    let schema = &ledger["measurement_schema"];
    if schema["fail_closed"] != Value::Bool(true) {
        return Err("measurement_schema.fail_closed must be true".into());
    }
    let stage = &schema["vs_rho"];
    let strs = |v: &Value| -> Vec<String> {
        v.as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default()
    };
    let required = strs(&stage["required"]);
    let global = strs(&schema["global_provenance_required"]);
    let missing = |keys: &[String]| -> Vec<String> {
        keys.iter()
            .filter(|k| !field_present(report, k))
            .cloned()
            .collect()
    };
    let (missing_stage, missing_global) = (missing(&required), missing(&global));
    let mut errors: Vec<String> = Vec::new();
    let mut err = |s: String| errors.push(s);
    let get = |k: &str| report.get(k);

    for key in ["candidate_id", "workload_id", "run_id"] {
        if text(get(key)).is_none() {
            err(format!("{key} must be a nonempty string"));
        }
    }
    let candidate_id = get("candidate_id").and_then(Value::as_str);
    let workload_id = get("workload_id").and_then(Value::as_str);
    let run_id = get("run_id").and_then(Value::as_str);
    let csha = get("candidate_manifest_sha256").and_then(Value::as_str);
    let wsha = get("workload_manifest_sha256").and_then(Value::as_str);
    for (key, v) in [
        ("candidate_manifest_sha256", csha),
        ("workload_manifest_sha256", wsha),
    ] {
        match v {
            Some(h) if h.chars().count() == 64 => {
                if !is_lower_hex(h) {
                    err(format!("{key} must be lowercase hex"));
                }
            }
            _ => err(format!("{key} must be a full SHA-256 hex digest")),
        }
    }
    if let (Some(cid), Some(csha)) = (candidate_id, csha) {
        // re.search(r"h([0-9a-f]{12,64})$", cid): `$` also matches
        // before one trailing newline.
        let body = cid.strip_suffix('\n').unwrap_or(cid);
        let digest = body
            .rsplit_once('h')
            .map(|(_, d)| d)
            .filter(|d| (12..=64).contains(&d.len()) && is_lower_hex(d));
        match digest {
            Some(d) if cid.starts_with("IC1") => {
                if !csha.starts_with(d) {
                    err("candidate_id digest must match candidate_manifest_sha256".into());
                }
            }
            _ => err("candidate_id must use the IC1 identity format with a digest suffix".into()),
        }
    }
    if let (Some(wid), Some(wsha)) = (workload_id, wsha) {
        if wid.chars().count() != 12 || !is_lower_hex(wid) {
            err("workload_id must be 12 lowercase hex digits".into());
        } else if !wsha.starts_with(wid) {
            err("workload_id must match workload_manifest_sha256".into());
        }
    }
    if let (Some(cid), Some(wid), Some(rid)) = (candidate_id, workload_id, run_id) {
        if !cid.is_empty() && !wid.is_empty() && !rid.is_empty() {
            let n = rid.strip_prefix(&format!("{cid}W{wid}R"));
            let ok = n.is_some_and(|n| {
                n.starts_with(|c: char| matches!(c, '1'..='9'))
                    && n.chars().all(|c| c.is_ascii_digit())
            });
            if !ok {
                err("run_id must be <candidate_id>W<workload_id>R<run-number>".into());
            }
        }
    }
    if !get("target_count").is_some_and(|v| v.is_i64() || v.is_u64())
        || get("target_count").and_then(Value::as_i64) != Some(1)
    {
        err("target_count must equal 1".into());
    }
    match (text(get("ic_target_hash")), text(get("rho_target_hash"))) {
        (Some(a), Some(b)) => {
            if a != b {
                err("IC and rho target hashes must match".into());
            }
        }
        _ => err("IC and rho target hashes must be nonempty strings".into()),
    }
    if !stage["timing_class_enum"]
        .as_array()
        .is_some_and(|e| get("timing_class").is_some_and(|t| e.contains(t)))
    {
        err("timing_class must be single_target_online_wall".into());
    }
    for key in [
        "same_resource_envelope",
        "ic_scalar_verified",
        "rho_scalar_verified",
        "independent_validation",
    ] {
        if get(key) != Some(&Value::Bool(true)) {
            err(format!("{key} must be true"));
        }
    }
    for key in [
        "ic_replay_certificate_sha256",
        "rho_replay_certificate_sha256",
    ] {
        let ok = get(key)
            .and_then(Value::as_str)
            .is_some_and(|h| h.chars().count() == 64 && is_lower_hex(h));
        if !ok {
            err(format!("{key} must be a full lowercase SHA-256 digest"));
        }
    }
    let ic_env = get("ic_resource_envelope").and_then(Value::as_object);
    let rho_env = get("rho_resource_envelope").and_then(Value::as_object);
    for (key, v) in [
        ("ic_resource_envelope", ic_env),
        ("rho_resource_envelope", rho_env),
    ] {
        if !v.is_some_and(|o| !o.is_empty()) {
            err(format!("{key} must be a nonempty object"));
        }
    }
    if let (Some(a), Some(b)) = (ic_env, rho_env) {
        if a != b {
            err("IC and rho resource envelopes must match exactly".into());
        }
    }
    let ic_ms = positive_cost(get("ic_online_wall_ms"));
    let rho_ms = positive_cost(get("rho_online_wall_ms"));
    let speedup = positive_cost(get("online_speedup"));
    match (ic_ms, rho_ms, speedup) {
        (Some(ic), Some(rho), Some(sp)) => {
            if !isclose(sp, rho / ic, 1e-9, 1e-12) {
                err("online_speedup must equal rho_online_wall_ms / ic_online_wall_ms".into());
            }
        }
        _ => err("online times and speedup must be positive finite numbers".into()),
    }
    match get("ic_online_phase_ms").and_then(Value::as_object) {
        None => err("ic_online_phase_ms must be an object".into()),
        Some(phases) => {
            let fields = strs(&stage["ic_online_phase_fields"]);
            let mut values = Vec::new();
            for key in &fields {
                match phases
                    .get(key)
                    .and_then(Value::as_f64)
                    .filter(|v| v.is_finite() && *v >= 0.0)
                {
                    Some(v) => values.push(v),
                    None => err(format!(
                        "ic_online_phase_ms.{key} must be a finite nonnegative number"
                    )),
                }
            }
            if let Some(ic) = ic_ms.filter(|_| values.len() == fields.len()) {
                // math.fsum: an exactly rounded sum; five terms of this
                // size sum exactly enough in f64 for a 1e-6 tolerance.
                let sum: f64 = values.iter().sum();
                if !isclose(sum, ic, 1e-6, 1e-3) {
                    err("IC exclusive phase costs must sum to ic_online_wall_ms".into());
                }
            }
        }
    }
    match get("online_interval").and_then(Value::as_object) {
        None => err("online_interval must be an object".into()),
        Some(interval) => {
            let keys = strs(&stage["online_interval_event_fields"])
                .into_iter()
                .chain(strs(&stage["online_interval_stage_fields"]));
            for key in keys {
                let v = interval.get(&key);
                let valid = if key.ends_with("_event") {
                    text(v).is_some()
                } else {
                    v.and_then(Value::as_array)
                        .is_some_and(|a| !a.is_empty() && a.iter().all(|x| text(Some(x)).is_some()))
                };
                if !valid {
                    err(format!("online_interval.{key} is missing or invalid"));
                }
            }
            for (key, need) in ONLINE_REQUIRED_STAGES {
                if let Some(have) = interval.get(*key).and_then(Value::as_array) {
                    let mut missing: Vec<&str> = need
                        .iter()
                        .copied()
                        .filter(|s| !have.iter().any(|h| h.as_str() == Some(s)))
                        .collect();
                    missing.sort_unstable();
                    if !missing.is_empty() {
                        err(format!(
                            "online_interval.{key} is missing required stages: {}",
                            missing.join(", ")
                        ));
                    }
                }
            }
        }
    }
    match get("rho_policy").and_then(Value::as_object) {
        None => err("rho_policy must be an object".into()),
        Some(policy) => {
            for key in strs(&stage["rho_policy_integer_fields"]) {
                let min = stage["rho_policy_minimums"]
                    .get(&key)
                    .and_then(Value::as_i64);
                let v = policy.get(&key);
                let int = v.filter(|v| v.is_i64() || v.is_u64());
                let ok = int.is_some()
                    && min.is_none_or(|m| int.and_then(Value::as_i64).is_none_or(|x| x >= m));
                if !ok {
                    err(format!(
                        "rho_policy.{key} must be an integer >= {}",
                        min.map_or("None".to_string(), |m| m.to_string())
                    ));
                }
            }
            for key in strs(&stage["rho_policy_string_fields"]) {
                if text(policy.get(&key)).is_none() {
                    err(format!("rho_policy.{key} must be a nonempty string"));
                }
            }
        }
    }
    let ok = missing_stage.is_empty() && missing_global.is_empty() && errors.is_empty();
    Ok(json!({
        "schema_version": ledger.get("schema_version"),
        "stage": "vs_rho",
        "fail_closed": true,
        "status": if ok { "PASS" } else { "FAIL" },
        "missing_stage_fields": missing_stage,
        "missing_global_provenance": missing_global,
        "validation_errors": errors,
        "required_stage_fields": required,
        "required_global_provenance": global,
    }))
}

/// Every problem a [`check_vs_rho`] result lists, one line each.
pub fn check_problems(check: &Value) -> Vec<String> {
    let list = |k: &str| -> Vec<String> {
        check[k]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default()
    };
    list("missing_stage_fields")
        .into_iter()
        .map(|k| format!("missing field {k}"))
        .chain(
            list("missing_global_provenance")
                .into_iter()
                .map(|k| format!("missing provenance field {k}")),
        )
        .chain(list("validation_errors"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;

    #[test]
    fn canonical_utf8_writes_non_ascii_raw() {
        assert_eq!(
            canonical_utf8(&json!({"b": "é", "a": [1, -2]})).unwrap(),
            "{\"a\":[1,-2],\"b\":\"é\"}"
        );
        assert!(canonical_utf8(&json!(1.5)).is_err());
    }

    /// The K_1, n = 17 cross-check fixture: the curve, the
    /// `koblitz-orbit:divisor=0;1` base, an `ic.pipeline` method with
    /// placeholder source hashes, and the generator as a public target.
    struct Fixture {
        fixture: Value,
        report: Value,
        method: Value,
        inventory: Value,
        columns: ColumnMap,
        curve_id: String,
        candidate: (String, String),
        workload: (String, String),
        envelope: Value,
    }

    fn k1n17() -> Fixture {
        let inst = koblitz_instance(1, 17).unwrap();
        let kb = Kb::of(&inst).unwrap();
        let kc = inst.koblitz.as_ref().unwrap();
        use num_traits::ToPrimitive;
        let lambda = (&kc.lambda % num_bigint::BigUint::from(inst.r))
            .to_u64()
            .unwrap();
        let fixture = json!({
            "degree": inst.n, "curve_a": inst.a, "subgroup_order": inst.r,
            "irreducible": {"degree": inst.n, "low_terms": inst.irreducible.low_terms},
            "group_order": inst.group_order, "cofactor": inst.cofactor,
            "generator": [inst.generator.x, inst.generator.y], "lambda": lambda,
            "targets": [[inst.generator.x, inst.generator.y]],
            "target_seeds": [{"law": "hash_to_subgroup_v1", "seed": 1, "index": 0, "counter": 0}],
            "target_scalar_constructed": false,
        });
        let spec = crate::cryptanalysis::ecbench::workload::CurveSpec::Koblitz { a: 1, n: 17 };
        let dump = dump_factor_base(&spec, "koblitz-orbit:divisor=0;1").unwrap();
        let points: Vec<(u64, u64)> = dump
            .points
            .iter()
            .map(|p| (hex(&p.x).unwrap(), hex(&p.y).unwrap()))
            .collect();
        let (inventory, columns) =
            factor_base_inventory(&kb, &points, dump.factor_base.columns).unwrap();
        // identity.py is told the usable base's orbit count; the dead and
        // merged raw columns ride in the method record.
        let report = json!({"factor_base": encode(&points),
                            "columns": inventory["effective_columns"],
                            "column_convention": "cofactor"});
        let imp = Implementation {
            components: vec![
                ("ic_pipeline".into(), "1".repeat(64)),
                ("ic_plugins".into(), "2".repeat(64)),
                ("ic_solvers".into(), "3".repeat(64)),
                ("ic_linalg".into(), "4".repeat(64)),
                ("relation_loop".into(), "5".repeat(64)),
                ("harness".into(), "6".repeat(64)),
            ],
        };
        let m = crate::cryptanalysis::ecbench::methods::resolve(
            &crate::cryptanalysis::ecbench::methods::MethodSpec {
                id: "ic.pipeline".into(),
                params: [
                    (
                        "factor_base".to_string(),
                        "koblitz-orbit:divisor=0;1".to_string(),
                    ),
                    ("oracle".to_string(), "mitm-frobenius:m=2".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        )
        .unwrap();
        let method = method_record(&m, &kb, &dump.factor_base, &columns, &imp).unwrap();
        let (cid, csha, _) = candidate_manifest(&kb, method.clone(), inventory.clone()).unwrap();
        let curve = curve_record(&kb).unwrap();
        let curve_id = curve["curve"]["curve_id"].as_str().unwrap().to_string();
        let envelope = json!({"threads": 1});
        let (wid, wsha, _) = workload_manifest(
            &curve_id,
            (inst.generator.x, inst.generator.y),
            fixture["target_seeds"][0].clone(),
            7,
            &envelope,
        )
        .unwrap();
        Fixture {
            fixture,
            report,
            method,
            inventory,
            columns,
            curve_id,
            candidate: (cid, csha),
            workload: (wid, wsha),
            envelope,
        }
    }

    /// Every identity this module computes on the K_1, n = 17 fixture,
    /// against `identity.py`'s values for the same inputs (computed with
    /// research/ic_candidate_tournament_20260915/identity.py from the
    /// output of `print_identity_inputs`, Python 3.13.1, 2026-10-02).  A
    /// change to the IC sources' identity fields, the factor-base plug-in
    /// or `ic.pipeline`'s defaults moves these: rerun that cross-check
    /// before repinning.
    #[test]
    fn identities_match_identity_py() {
        let f = k1n17();
        assert_eq!(f.curve_id, "EC1N17Ckb1hbbe2b5b6b1e6");
        // The base holds (0, 1), the 2-torsion point above x = 0, whose
        // image under [2] is the identity: one dead raw column.
        assert_eq!(
            f.columns,
            ColumnMap {
                raw: 14,
                identity_image: 1,
                merged: 0
            }
        );
        let inv = &f.inventory;
        assert_eq!(inv["geometric_point_count"], 443);
        assert_eq!(inv["usable_point_count"], 442);
        assert_eq!(inv["identity_images"], 1);
        assert_eq!(inv["duplicate_nonidentity_images"], 0);
        assert_eq!(inv["effective_columns"], 13);
        assert_eq!(
            inv["geometric_set_sha256"],
            "0c9cb58e5214b4f1ff6fcdae1661225858eb1a4d5ba8555cc253ff0cf21366ee"
        );
        assert_eq!(
            inv["geometric_order_sha256"],
            "3245ed0371655245736234b0ccf7db5d1a277194339a8c8ae743e04fe8f96ef9"
        );
        assert_eq!(
            inv["usable_set_sha256"],
            "8ce9f0959ce27c4dec3dd54ac6e2c77603d6a0e886408ac5ec5e48a6393952fa"
        );
        assert_eq!(
            f.candidate.0,
            "IC1N17Ckb1fb442PDP2mitmfrobeniusRCwalkLAincrementalgaussTDjointISO0h9b819c32f0ca"
        );
        assert_eq!(
            f.candidate.1,
            "9b819c32f0caa3bd830af437089fd30c658c608ebcd065c49a0d8b218114b7cd"
        );
        assert_eq!(f.workload.0, "daaac62e5a1e");
        assert_eq!(
            f.workload.1,
            "daaac62e5a1eb463f185f441fe5a59f4ba8c442fffdb839be2dc161056b76341"
        );
    }

    #[test]
    fn a_fold_other_than_sign_and_frobenius_is_refused() {
        let inst = koblitz_instance(1, 17).unwrap();
        let kb = Kb::of(&inst).unwrap();
        let spec = crate::cryptanalysis::ecbench::workload::CurveSpec::Koblitz { a: 1, n: 17 };
        let dump = dump_factor_base(&spec, "koblitz-orbit:divisor=0;1,no_fold=1").unwrap();
        let points: Vec<(u64, u64)> = dump
            .points
            .iter()
            .map(|p| (hex(&p.x).unwrap(), hex(&p.y).unwrap()))
            .collect();
        let err = factor_base_inventory(&kb, &points, dump.factor_base.columns).unwrap_err();
        assert!(err.contains("signed Frobenius"), "{err}");
    }

    /// The inputs identity.py needs for the cross-check, and this port's
    /// answers:
    /// `cargo test --release --lib print_identity_inputs -- --ignored --nocapture`,
    /// then feed the `IDENTITY_INPUTS` JSON to `identity.curve_record`,
    /// `factor_base_inventory`, `candidate_manifest` and `workload_manifest`.
    #[test]
    #[ignore]
    fn print_identity_inputs() {
        let f = k1n17();
        println!(
            "IDENTITY_INPUTS {}",
            json!({"fixture": f.fixture, "report": f.report, "method": f.method,
                   "input_law": format!("ecbench.{PUBLIC_TARGET_LAW}"), "algorithm_seed": 7,
                   "resource_envelope": f.envelope,
                   "rust": {"curve_id": f.curve_id, "inventory": f.inventory,
                            "candidate_id": f.candidate.0, "record_sha256": f.candidate.1,
                            "workload_id": f.workload.0, "workload_sha256": f.workload.1}})
        );
    }
}
