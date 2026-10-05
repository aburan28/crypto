//! Queued trait detection: shard a walk's curves into `taskq` jobs
//! (`taskq/README.md`), run the detectors one shard per job, and merge the
//! results.
//!
//! - [`plan`] writes one `taskq.task-spec/v1` per shard.  The walk is
//!   deterministic, so each job rebuilds it in its setup step from a
//!   pinned commit and needs no shared storage.  The rebuilt routes file
//!   is hashed into the shard's output, so [`collect`] can prove every
//!   shard saw the same walk.
//! - [`run_shard`] runs [`super::traits::default_detectors`] on the
//!   curves with `index ≡ shard (mod of)`, the strided ownership of
//!   PR #1330's contract.  Shard 0 also runs the class audits.
//! - [`collect`] merges shard directories.  It refuses a missing or
//!   duplicated shard, a curve seen twice or not at all, and shards from
//!   different walks or shard counts.
//!
//! Any queue works if it can run `argv` with `TASKQ_OUTPUT_DIR` set, or
//! with `--out`; without a queue the same shards run locally in a loop.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use num_bigint::BigUint;
use rayon::prelude::*;

use super::curve::Model;
use super::field::Field;
use super::record::{sha256_hex, V};
use super::traits::{self, CurveCtx};
use super::walk::StartCurve;
use crate::cryptanalysis::curve_id;

pub const SHARD_SCHEMA: &str = "isogeny-walk-trait-shard/v1";
pub const TRAITS_FILE: &str = "traits.jsonl";
pub const METRICS_FILE: &str = "metrics.json";
pub const CLASS_FILE: &str = "class_audits.json";

/// What one shard writes.
pub struct ShardOutput {
    /// One canonical JSON object per curve, in index order.
    pub jsonl: String,
    /// `metrics.json`: taskq parses it into the run's metrics.
    pub metrics: V,
    /// Shard 0 only.
    pub class_audits: Option<V>,
}

fn dec(v: &serde_json::Value, what: &str) -> Result<BigUint, String> {
    v.as_str()
        .and_then(|s| BigUint::parse_bytes(s.as_bytes(), 10))
        .ok_or_else(|| format!("{what}: expected a decimal string"))
}

/// Run the detectors on shard `shard` of `of` of a walk's
/// `isogeny_routes.json` (as text).  `audit_sample` walked curves of this
/// shard are re-audited by the class audits when `shard == 0`.
pub fn run_shard(
    routes_text: &str,
    start: &StartCurve,
    shard: usize,
    of: usize,
    audit_sample: usize,
) -> Result<ShardOutput, String> {
    if of == 0 || shard >= of {
        return Err(format!("shard {shard} of {of}: need 0 <= shard < of"));
    }
    let routes: serde_json::Value = serde_json::from_str(routes_text).map_err(|e| e.to_string())?;
    let nodes = routes["curve_nodes"].as_array().ok_or("curve_nodes")?;
    let f = Field::new(&start.p).ok_or("bad field")?;
    let detectors = traits::default_detectors();
    let mine: Vec<(usize, &serde_json::Value)> = nodes
        .iter()
        .enumerate()
        .filter(|(i, _)| i % of == shard)
        .collect();
    type Row = (usize, String, (BigUint, BigUint, BigUint, BigUint));
    let rows: Vec<Result<Row, String>> = mine
        .par_iter()
        .map(|(i, node)| {
            let r = node["ref"].as_str().ok_or("ref")?;
            let co = node["coefficients_a1_a2_a3_a4_a6"]
                .as_array()
                .ok_or("coefficients")?;
            let (a, b) = (dec(&co[3], "a4")?, dec(&co[4], "a6")?);
            let g = node["generator_G"].as_array().ok_or("generator")?;
            let (gx, gy) = (dec(&g[0], "gx")?, dec(&g[1], "gy")?);
            let id = curve_id::prime(&start.p, &a, &b, &start.order).ok_or("identity")?;
            if node["icv1_slug"].as_str() != Some(id.slug.as_str()) {
                return Err(format!("{r}: slug differs from the recomputed {}", id.slug));
            }
            let model = Model {
                a: f.from_big(&a),
                b: f.from_big(&b),
            };
            let ctx = CurveCtx {
                field: &f,
                model: &model,
                generator: (f.from_big(&gx), f.from_big(&gy)),
                start,
            };
            let traits = V::Map(
                detectors
                    .iter()
                    .map(|d| {
                        (
                            d.name().to_string(),
                            V::map(vec![
                                ("value", d.detect(&ctx)),
                                ("status", V::s(d.status())),
                            ]),
                        )
                    })
                    .collect(),
            );
            let line = V::map(vec![
                ("index", V::int(i)),
                ("ref", V::s(r)),
                ("icv1_slug", V::s(id.slug)),
                (
                    "curve_uid",
                    V::s(node["curve_uid"].as_str().unwrap_or("").to_string()),
                ),
                ("traits", traits),
            ])
            .canonical_json();
            Ok((*i, line, (a, b, gx, gy)))
        })
        .collect();
    let mut jsonl = String::new();
    let mut sample = Vec::new();
    for row in rows {
        let (i, line, params) = row?;
        jsonl.push_str(&line);
        jsonl.push('\n');
        if i > 0 && sample.len() < audit_sample {
            sample.push(params);
        }
    }
    let class_audits = (shard == 0).then(|| traits::class_audits(start, &sample));
    let metrics = V::map(vec![
        ("schema", V::s(SHARD_SCHEMA)),
        ("isogeny_routes_sha256", V::s(sha256_hex(routes_text))),
        ("root", V::s(start.name.clone())),
        ("shard", V::int(shard)),
        ("of", V::int(of)),
        ("curves_total", V::int(nodes.len())),
        ("curves_in_shard", V::int(mine.len())),
        ("traits_sha256", V::s(sha256_hex(&jsonl))),
        (
            "detectors",
            V::Seq(detectors.iter().map(|d| V::s(d.name())).collect()),
        ),
    ]);
    Ok(ShardOutput {
        jsonl,
        metrics,
        class_audits,
    })
}

/// Write a shard's files into `dir`.
pub fn write_shard(dir: &Path, out: &ShardOutput) -> Result<(), String> {
    std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    let w = |name: &str, text: &str| {
        std::fs::write(dir.join(name), text).map_err(|e| format!("{name}: {e}"))
    };
    w(TRAITS_FILE, &out.jsonl)?;
    w(METRICS_FILE, &out.metrics.json())?;
    if let Some(c) = &out.class_audits {
        w(CLASS_FILE, &c.json())?;
    }
    Ok(())
}

/// Merge shard directories into one `traits.jsonl` (index order) and a
/// summary.  Every shard must come from the same walk and shard count,
/// and every curve must appear exactly once.
pub fn collect(dirs: &[PathBuf]) -> Result<(String, V, Option<String>), String> {
    let mut walk: Option<(String, u64, u64)> = None;
    let mut shards = BTreeMap::new();
    let mut lines: BTreeMap<u64, String> = BTreeMap::new();
    let mut class = None;
    for d in dirs {
        let read = |n: &str| {
            std::fs::read_to_string(d.join(n)).map_err(|e| format!("{}: {e}", d.join(n).display()))
        };
        let m: serde_json::Value =
            serde_json::from_str(&read(METRICS_FILE)?).map_err(|e| e.to_string())?;
        if m["schema"].as_str() != Some(SHARD_SCHEMA) {
            return Err(format!("{}: not a trait shard", d.display()));
        }
        let key = (
            m["isogeny_routes_sha256"]
                .as_str()
                .unwrap_or("")
                .to_string(),
            m["of"].as_u64().ok_or("of")?,
            m["curves_total"].as_u64().ok_or("curves_total")?,
        );
        match &walk {
            None => walk = Some(key.clone()),
            Some(w) if *w != key => {
                return Err(format!("{}: from another walk or shard count", d.display()))
            }
            _ => {}
        }
        let shard = m["shard"].as_u64().ok_or("shard")?;
        if shards.insert(shard, d.display().to_string()).is_some() {
            return Err(format!("shard {shard} given twice"));
        }
        let text = read(TRAITS_FILE)?;
        if m["traits_sha256"].as_str() != Some(sha256_hex(&text).as_str()) {
            return Err(format!(
                "{}: traits.jsonl does not match its metrics",
                d.display()
            ));
        }
        let mut count = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            let i = v["index"].as_u64().ok_or("index")?;
            if i % key.1 != shard {
                return Err(format!("curve {i} is not owned by shard {shard}"));
            }
            if lines.insert(i, line.to_string()).is_some() {
                return Err(format!("curve {i} appears twice"));
            }
            count += 1;
        }
        if Some(count) != m["curves_in_shard"].as_u64() {
            return Err(format!(
                "{}: line count differs from its metrics",
                d.display()
            ));
        }
        if shard == 0 {
            class = std::fs::read_to_string(d.join(CLASS_FILE)).ok();
        }
    }
    let (sha, of, total) = walk.ok_or("no shards")?;
    let missing: Vec<u64> = (0..of).filter(|s| !shards.contains_key(s)).collect();
    if !missing.is_empty() {
        return Err(format!("missing shards {missing:?} of {of}"));
    }
    if lines.len() as u64 != total || lines.keys().next_back().is_some_and(|&k| k + 1 != total) {
        return Err(format!("{} of {total} curves covered", lines.len()));
    }
    let mut merged = String::new();
    // trait (flattened, `a.b` for nested values) -> value -> curves.
    let mut dist: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    for l in lines.values() {
        merged.push_str(l);
        merged.push('\n');
        let v: serde_json::Value = serde_json::from_str(l).map_err(|e| e.to_string())?;
        if let Some(traits) = v["traits"].as_object() {
            for (name, t) in traits {
                tally(&mut dist, name, &t["value"]);
            }
        }
    }
    let summary = V::map(vec![
        ("schema", V::s("isogeny-walk-traits-collected/v1")),
        ("isogeny_routes_sha256", V::s(sha)),
        ("shards", V::int(of)),
        ("curves", V::int(total)),
        ("traits_sha256", V::s(sha256_hex(&merged))),
        (
            "class_audits",
            V::s(if class.is_some() {
                CLASS_FILE
            } else {
                "absent"
            }),
        ),
        (
            "distributions",
            V::Map(
                dist.into_iter()
                    .map(|(k, m)| {
                        (
                            k,
                            V::Map(m.into_iter().map(|(val, c)| (val, V::int(c))).collect()),
                        )
                    })
                    .collect(),
            ),
        ),
    ]);
    Ok((merged, summary, class))
}

/// Count a trait value; objects are flattened to `name.key`.
fn tally(dist: &mut BTreeMap<String, BTreeMap<String, u64>>, name: &str, v: &serde_json::Value) {
    match v {
        serde_json::Value::Object(m) => {
            for (k, x) in m {
                tally(dist, &format!("{name}.{k}"), x);
            }
        }
        other => {
            let key = match other {
                serde_json::Value::String(s) => s.clone(),
                x => x.to_string(),
            };
            *dist
                .entry(name.to_string())
                .or_default()
                .entry(key)
                .or_default() += 1;
        }
    }
}

/// The walk a plan rebuilds, as `isogeny_walk walk` arguments (without
/// `--out`).
#[derive(Clone, Debug)]
pub struct PlanWalk {
    pub curve_args: Vec<String>,
    pub walk_args: Vec<String>,
}

/// The `taskq.task-spec/v1` documents for a sharded trait run: one per
/// shard, each rebuilding the walk at `commit` and detecting one shard.
impl PlanWalk {
    /// Sixteen hex digits naming this walk's arguments.
    pub fn key(&self) -> String {
        sha256_hex(&format!("{:?}|{:?}", self.curve_args, self.walk_args))[..16].to_string()
    }

    /// `<curve>-<key>`: the run id a stored walk is filed under.
    pub fn run_id(&self) -> String {
        let curve = self
            .curve_args
            .iter()
            .position(|a| a == "--curve")
            .and_then(|i| self.curve_args.get(i + 1))
            .map_or("curve".to_string(), |c| {
                c.to_ascii_lowercase()
                    .chars()
                    .filter(|ch| ch.is_ascii_alphanumeric())
                    .collect()
            });
        format!("{curve}-{}", self.key())
    }
}

/// Where planned jobs read the walk from and write their shards to.
#[derive(Clone, Debug, Default)]
pub struct PlanStore {
    /// `s3://bucket/prefix`: each shard publishes its output there.
    pub store: Option<String>,
    /// Fetch the walk from `store` instead of rebuilding it.
    pub walk_from_store: bool,
}

pub fn plan(
    walk: &PlanWalk,
    commit: &str,
    queue: &str,
    shards: usize,
    timeout_seconds: u64,
    storage: &PlanStore,
) -> Result<Vec<V>, String> {
    if storage.walk_from_store && storage.store.is_none() {
        return Err("--walk-from-store needs --store".into());
    }
    if let Some(s) = &storage.store {
        super::store::S3Loc::parse(s)?;
    }
    if commit.len() != 40
        || !commit
            .bytes()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
    {
        return Err(
            "--commit must be a full lower-case 40-hex sha (taskq refuses branch names)".into(),
        );
    }
    if shards == 0 {
        return Err("--shards must be at least 1".into());
    }
    let key = walk.key();
    let run_id = walk.run_id();
    let dir = format!("target/isogeny-walk/{key}");
    let bin = "target/release/isogeny_walk".to_string();
    let mut walk_cmd = vec![bin.clone(), "walk".into()];
    walk_cmd.extend(walk.curve_args.iter().cloned());
    walk_cmd.extend(walk.walk_args.iter().cloned());
    // Shard 0 runs the class audits; the rebuilt walk skips them.
    walk_cmd.extend([
        "--no-class-audits".to_string(),
        "--out".to_string(),
        dir.clone(),
    ]);
    if let (true, Some(store)) = (storage.walk_from_store, &storage.store) {
        walk_cmd = vec![
            bin.clone(),
            "fetch".into(),
            "--from".into(),
            store.clone(),
            "--run".into(),
            run_id.clone(),
            "--out".into(),
            dir.clone(),
        ];
    }
    let strs = |v: &[String]| V::Seq(v.iter().map(|s| V::s(s.clone())).collect());
    Ok((0..shards)
        .map(|i| {
            let mut argv = vec![bin.clone(), "traits".into()];
            argv.extend(walk.curve_args.iter().cloned());
            argv.extend([
                "--dir".to_string(),
                dir.clone(),
                "--shard".into(),
                i.to_string(),
                "--of".into(),
                shards.to_string(),
            ]);
            if let Some(store) = &storage.store {
                argv.extend([
                    "--store".to_string(),
                    store.clone(),
                    "--run".into(),
                    run_id.clone(),
                ]);
            }
            V::map(vec![
                ("schema", V::s("taskq.task-spec/v1")),
                ("queue", V::s(queue)),
                ("kind", V::s("command")),
                (
                    "source",
                    V::map(vec![("repo", V::s("crypto")), ("commit", V::s(commit))]),
                ),
                (
                    "command",
                    V::map(vec![
                        ("argv", strs(&argv)),
                        (
                            "setup",
                            V::Seq(vec![
                                strs(&[
                                    "cargo".into(),
                                    "build".into(),
                                    "--release".into(),
                                    "--bin".into(),
                                    "isogeny_walk".into(),
                                ]),
                                strs(&walk_cmd),
                            ]),
                        ),
                        ("cwd", V::s(".")),
                    ]),
                ),
                (
                    "limits",
                    V::map(vec![
                        ("timeout_seconds", V::int(timeout_seconds)),
                        ("setup_timeout_seconds", V::int(timeout_seconds.max(3600))),
                    ]),
                ),
                ("retry", V::map(vec![("max_attempts", V::int(3))])),
                (
                    "idempotency_key",
                    V::s(format!(
                        "isogeny-walk-traits/{}/{}/{i}-of-{shards}",
                        key,
                        &commit[..12]
                    )),
                ),
                (
                    "labels",
                    V::map(vec![
                        ("experiment", V::s("isogeny-walk-traits")),
                        ("walk", V::s(&key)),
                        ("shard", V::s(format!("{i}/{shards}"))),
                    ]),
                ),
            ])
        })
        .collect())
}
