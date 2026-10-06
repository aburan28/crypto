//! The single-target rule's comparison (ledger §23; IC_TOOL_PROGRAM.md,
//! "The rule's comparison"), natively: the claim rows §23's `claims.py`
//! built and checked, and the figures its `analyse.py` read from them.
//!
//! - `claims`: for each row the run produced (the first clean, complete
//!   attempt of each `T<i>-R<run>`), the canonical IC1 candidate and
//!   workload and the run ID ([`identity`]); both arms' certificates
//!   replayed outside `ic`, in [`oracle`]'s arithmetic; rho's EC1 reference
//!   identity; and the `vs_rho` claim, checked by [`autolab`].
//! - `analyse`: per size, the online speedup, rho at the canonical step,
//!   the cold ratio and the break-even count, each with a bootstrap
//!   interval drawn as `random.Random` drew it ([`pyrandom`]); the
//!   precomputation boundary, the A/A and the diagnostics; across sizes,
//!   the slopes of `ln S` on `ln r`.  Every float is computed as Python
//!   3.11 computed it (`fsum` for `fmean`, a left-to-right `sum`, libm's
//!   `pow`), so §23's `analysis.json` is reproduced byte for byte.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use crypto_lib::cryptanalysis::curve_id;

use super::autolab;
use super::bench;
use super::identity;
use super::json::{self, J};
use super::oracle::Curve;
use super::pin;
use super::pyrandom::Random;
use super::runs;
use super::stats;
use super::suite::{self, Row};

/// §23's sizes, `(a, n)` in the order of `r`.
pub const SIZES: [(u8, u32); 6] = [(1, 47), (0, 57), (0, 41), (0, 53), (1, 59), (0, 61)];
const TARGET_SEED_BASE: i128 = 23000;
const BOOT: usize = 10_000;
const BOOT_SEED: u128 = 2_309_300;
const BL_PRODUCT: f64 = 1.93 * 1.21;
/// Student's t, 97.5%, by degrees of freedom.
const T975: [f64; 8] = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306];

/// Who replayed the certificates, named in every replay record.
pub const CHECKER: &str = "icprog's oracle (its own field arithmetic, not ic's)";

/// What a comparison's analysis sets beside each size's figures.
pub enum Against {
    /// §23: the figures its `predict.py` declared before the run, in the
    /// comparison's `prediction.json`.
    Prediction,
    /// A later comparison: the same figures from an earlier one's analysis
    /// (a path under the repository), under `key`.
    Earlier {
        analysis: &'static str,
        key: &'static str,
    },
}

/// What a comparison's records say about how its rows were made.
pub struct Constants {
    /// The sources each arm executes, hashed at the binary's commit.
    pub ic_sources: &'static [&'static str],
    pub rho_sources: &'static [&'static str],
    pub input_law: &'static str,
    pub non_claims: &'static [&'static str],
    /// The resource envelope's `isolation` and `process`.
    pub isolation: &'static str,
    pub process: &'static str,
    pub what_this_is: &'static str,
    /// The six curves' EC1 records, under the repository.
    pub curve_ids: &'static str,
    pub against: Against,
    /// An earlier comparison's claims, under the repository, whose
    /// logarithms every row must repeat; and its run tree, whose parameter
    /// files this one prices.
    pub earlier: Option<(&'static str, &'static str)>,
    /// Whether the comparison checks its rho against the strong reference
    /// (`docs/ic/BOUNDARY_TARGETS.md`, 2026-10-01), which postdates §23.
    pub reference_check: bool,
    /// Whether the comparison's runs are frozen.
    pub frozen: bool,
}

const S23_IC_SOURCES: &[&str] = &[
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/koblitz_sparse_la.rs",
    "src/cryptanalysis/semaev_decomp.rs",
    "src/cryptanalysis/ic_measurement.rs",
    "src/bin/ic/price.rs",
    "src/bin/ic/workflow.rs",
    "src/bin/ic/experiment.rs",
];
const S23_RHO_SOURCES: &[&str] = &[
    "src/cryptanalysis/ic_boundary.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/semaev_decomp.rs",
    "src/cryptanalysis/ic_measurement.rs",
    "src/bin/ic/price.rs",
];
const S23_INPUT_LAW: &str =
    "public hash-to-curve, ic workflow's domain ic-workflow-public-target-v1, seed 23000 + i";
const S23_NON_CLAIMS: &[&str] = &[
    "one public target per row; no multi-target or batch result",
    "the online interval excludes the reusable set-up, which is reported separately and dominates the cold cost",
    "rho here has no precomputed table of distinguished points; a generic walk with one (Bernstein and Lange) \
     is a different reference and is not measured",
    "no statement about ECC2K-130 or any curve past n = 61; m = 83 is not run (AGENTS.md section 8a)",
    "wall time on one x86-64 container; no claim for other hardware classes",
];
const S23_PROCESS: &str = "one process: the index calculus, then rho, both rebuilt each repetition";
const S23_CURVE_IDS: &str = "research/ic_single_target_20260930/curve_ids.json";

/// Ledger §23's.
pub const S23: Constants = Constants {
    ic_sources: S23_IC_SOURCES,
    rho_sources: S23_RHO_SOURCES,
    input_law: S23_INPUT_LAW,
    non_claims: S23_NON_CLAIMS,
    isolation: "tools/isolated_bench.py run --wait --cpus 2, uncontended",
    process: S23_PROCESS,
    what_this_is: "Ledger §23's figures, from the checked rows.",
    curve_ids: S23_CURVE_IDS,
    against: Against::Prediction,
    earlier: None,
    reference_check: false,
    frozen: true,
};

/// The rule's comparison at baseline v2 (`research/ic_tool_program/rule/v2`):
/// §23's rows, sources, law and boundary, on the native runner.  The arms'
/// sources are §23's list: v1 and v2 changed `koblitz_index_calculus.rs`
/// only, and nothing the arms run has been added since.
pub const V2: Constants = Constants {
    ic_sources: S23_IC_SOURCES,
    rho_sources: S23_RHO_SOURCES,
    input_law: S23_INPUT_LAW,
    non_claims: S23_NON_CLAIMS,
    isolation: "isolated_bench (src/bin/isolated_bench.rs) run --wait --cpus 2, uncontended",
    process: S23_PROCESS,
    what_this_is: "The rule's comparison at baseline v2: §23's figures, from the checked rows.",
    curve_ids: S23_CURVE_IDS,
    against: Against::Earlier {
        analysis: "research/ic_single_target_20260930/analysis.json",
        key: "s23_then_v2",
    },
    earlier: Some((
        "research/ic_single_target_20260930/claims",
        "research/ic_single_target_20260930/runs",
    )),
    reference_check: true,
    frozen: false,
};

/// The rule's comparison at baseline v3 (`research/ic_tool_program/rule/v3`):
/// rule v2's protocol on v3's binary, main's `995ea207` (R07).  The arms'
/// sources are still §23's list: main changed five of them after v2, and
/// added no file the arms run.
pub const V3: Constants = Constants {
    ic_sources: S23_IC_SOURCES,
    rho_sources: S23_RHO_SOURCES,
    input_law: S23_INPUT_LAW,
    non_claims: S23_NON_CLAIMS,
    isolation: "isolated_bench (src/bin/isolated_bench.rs) run --wait --cpus 2, uncontended",
    process: S23_PROCESS,
    what_this_is: "The rule's comparison at baseline v3: §23's figures, from the checked rows.",
    curve_ids: S23_CURVE_IDS,
    against: Against::Earlier {
        analysis: "research/ic_single_target_20260930/analysis.json",
        key: "s23_then_v3",
    },
    earlier: Some((
        "research/ic_single_target_20260930/claims",
        "research/ic_single_target_20260930/runs",
    )),
    reference_check: true,
    frozen: false,
};

/// Where a comparison lives and what it writes.
pub struct Comparison {
    /// The repository: the ledger, and what record paths are relative to.
    pub root: PathBuf,
    /// The repository whose objects hold the binary's build commit.
    pub git: PathBuf,
    /// The comparison's directory (`curve_ids.json`, `prediction.json`).
    pub here: PathBuf,
    pub runs: PathBuf,
    /// Where `claims/` and `manifests/` go.
    pub out: PathBuf,
    pub constants: &'static Constants,
}

fn read(path: &Path) -> Result<J, String> {
    json::read(path)
}

fn write(path: &Path, text: String) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// `str(v)` for the values the records interpolate.
fn py_str(v: &J) -> String {
    match v {
        J::Null => "None".into(),
        J::Bool(b) => if *b { "True" } else { "False" }.into(),
        J::Int(i) => i.to_string(),
        J::Float(x) => json::py_float(*x),
        J::Str(s) => s.clone(),
        other => json::dumps_line(other, false),
    }
}

fn num(v: &J, what: &str) -> Result<f64, String> {
    v.as_f64().ok_or_else(|| format!("{what} is not a number"))
}

/// The names in `dir`, sorted as Python sorts paths in one directory.
fn names(dir: &Path) -> Vec<String> {
    let mut out: Vec<String> = std::fs::read_dir(dir)
        .map(|it| {
            it.filter_map(Result::ok)
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| !n.starts_with('.'))
                .collect()
        })
        .unwrap_or_default();
    out.sort();
    out
}

/// `prefix*middle*suffix` as `glob` reads it, for one `*` on each side of
/// `middle` (or none, when `middle` is empty).
fn glob2(name: &str, prefix: &str, middle: &str, suffix: &str) -> bool {
    name.len() >= prefix.len() + middle.len() + suffix.len()
        && name.starts_with(prefix)
        && name.ends_with(suffix)
        && name[prefix.len()..name.len() - suffix.len()].contains(middle)
}

/// `T(\d+)-R(\d+)<suffix>`, matched whole.
fn target_run(name: &str, suffix: &str) -> Option<(u32, u32)> {
    let (t, r) = name
        .strip_suffix(suffix)?
        .strip_prefix('T')?
        .split_once("-R")?;
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    (digits(t) && digits(r)).then(|| Some((t.parse().ok()?, r.parse().ok()?)))?
}

struct Found {
    i: u32,
    run: u32,
    chosen: Option<(PathBuf, J, J)>,
    tried: Vec<J>,
}

impl Comparison {
    fn rel(&self, path: &Path) -> String {
        path.strip_prefix(&self.root)
            .unwrap_or(path)
            .to_string_lossy()
            .into_owned()
    }

    fn git_bytes(&self, commit: &str, path: &str) -> Result<Vec<u8>, String> {
        let out = Command::new("git")
            .arg("-C")
            .arg(&self.git)
            .args(["show", &format!("{commit}:{path}")])
            .output()
            .map_err(|e| format!("git: {e}"))?;
        if !out.status.success() {
            return Err(format!(
                "git show {commit}:{path}: {}",
                String::from_utf8_lossy(&out.stderr).trim()
            ));
        }
        Ok(out.stdout)
    }

    /// `{path: sha256}` of each source at the binary's commit.
    fn source_hashes(&self, commit: &str, paths: &[&str]) -> Result<J, String> {
        let mut kv = Vec::new();
        for p in paths {
            let h = identity::sha256_hex(&self.git_bytes(commit, p)?);
            kv.push((p.to_string(), J::Str(h)));
        }
        Ok(J::Obj(kv))
    }

    /// Each `(target, run)` of a size: every attempt, and the first clean,
    /// complete one.
    fn rows(&self, a: u8, n: u32) -> Result<Vec<Found>, String> {
        let d = self.runs.join(format!("k{a}n{n}"));
        let all = names(&d);
        let mut rows = Vec::new();
        for name in all.iter().filter(|x| glob2(x, "T", "-R", ".price.json")) {
            if name.contains("-retry") {
                continue;
            }
            let Some((i, run)) = target_run(name, ".price.json") else {
                continue;
            };
            let stem = &name[..name.len() - ".price.json".len()];
            let mut attempts = vec![name.clone()];
            attempts.extend(
                all.iter()
                    .filter(|x| glob2(x, &format!("{stem}-retry"), "", ".price.json"))
                    .cloned(),
            );
            let (mut chosen, mut tried) = (None, Vec::new());
            for att in attempts {
                let path = d.join(&att);
                let rec = d.join(format!("{}.isolation.jsonl", &att[..att.len() - 11]));
                let iso = if rec.exists() {
                    let text = std::fs::read_to_string(&rec)
                        .map_err(|e| format!("{}: {e}", rec.display()))?;
                    let last = text
                        .lines()
                        .last()
                        .ok_or_else(|| format!("{}: empty", rec.display()))?;
                    Some(json::parse(last)?.at("run")?.clone())
                } else {
                    None
                };
                let rep = if path.exists() {
                    read(&path)?
                } else {
                    json::obj([("status", J::Str("no report".into()))])
                };
                let status = rep.get("status").cloned().unwrap_or(J::Null);
                let ok = iso.as_ref().is_some_and(|r| {
                    r.get("exit_status")
                        .is_some_and(|s| json::py_eq(s, &J::Int(0)))
                        && !r.get("contended").is_some_and(J::truthy)
                }) && status.as_str() == Some("complete");
                let field = |k: &str| {
                    iso.as_ref()
                        .map_or(J::Null, |r| r.get(k).cloned().unwrap_or(J::Null))
                };
                tried.push(J::Obj(vec![
                    ("file".into(), J::Str(att.clone())),
                    ("status".into(), status),
                    ("contended".into(), field("contended")),
                    ("exit_status".into(), field("exit_status")),
                ]));
                if ok && chosen.is_none() {
                    chosen = Some((path, rep, iso.expect("checked")));
                }
            }
            rows.push(Found {
                i,
                run,
                chosen,
                tried,
            });
        }
        Ok(rows)
    }

    /// The commit the `ic` binary was built from, and its digest.
    fn build(host: &J) -> Result<(String, J), String> {
        // The native runner's manifest names the binary under `binaries`.
        if let Some(ic) = host.get("binaries").and_then(|b| b.get("ic")) {
            let commit = ic
                .at("built_from")?
                .as_str()
                .ok_or("host.json names no commit")?;
            return Ok((commit.to_string(), ic.at("sha256")?.clone()));
        }
        let built = host.get("ic_built_from").filter(|v| v.truthy());
        let commit = built
            .or_else(|| host.get("commit"))
            .and_then(J::as_str)
            .ok_or("host.json names no commit")?;
        Ok((commit.to_string(), host.at("ic_binary_sha256")?.clone()))
    }

    pub fn claims(&self) -> Result<String, String> {
        let k = self.constants;
        let host = read(&self.runs.join("host.json"))?;
        let (commit, binary_sha) = Self::build(&host)?;
        let ic_hashes = self.source_hashes(&commit, k.ic_sources)?;
        let rho_hashes = self.source_hashes(&commit, k.rho_sources)?;
        let curve_ids = read(&self.root.join(self.constants.curve_ids))?;
        let ledger = autolab::load_ledger(&self.root)?;
        let claims_dir = self.out.join("claims");
        let manifests = self.out.join("manifests");
        let mut curves: HashMap<String, Curve> = HashMap::new();
        let mut candidates: HashMap<String, J> = HashMap::new();
        let mut summary = Vec::new();
        for (a, n) in SIZES {
            for row in self.rows(a, n)? {
                let label = format!("k{a}n{n}/T{:02}-R{}", row.i, row.run);
                let Some((path, rep, iso)) = row.chosen else {
                    summary.push(J::Obj(vec![
                        ("row".into(), J::Str(label)),
                        ("status".into(), J::Str("no clean complete attempt".into())),
                        ("attempts".into(), J::Arr(row.tried)),
                    ]));
                    continue;
                };
                let params = read(
                    &self
                        .runs
                        .join(format!("k{a}n{n}/T{:02}.params.json", row.i)),
                )?;
                let ii = rep.at("identity_inputs")?;
                let fixture = ii.at("fixture")?;
                let curve_key = curve_key(fixture)?;
                if !curves.contains_key(&curve_key) {
                    curves.insert(curve_key.clone(), Curve::new(fixture)?);
                }
                let curve = &curves[&curve_key];
                let report = J::Obj(vec![
                    (
                        "factor_base_orbits".into(),
                        ii.at("factor_base_orbits")?.clone(),
                    ),
                    ("columns".into(), ii.at("columns")?.clone()),
                    (
                        "column_convention".into(),
                        ii.at("column_convention")?.clone(),
                    ),
                ]);
                let method = method_record(&params, fixture, &ic_hashes)?;
                let cand_key = json::dumps_compact(
                    &J::Arr(vec![J::Str(curve_key), report.clone(), method.clone()]),
                    true,
                    true,
                );
                if !candidates.contains_key(&cand_key) {
                    let c = identity::candidate_manifest(curve, fixture, &report, &method)?;
                    candidates.insert(cand_key.clone(), c);
                }
                let cand = &candidates[&cand_key];
                let envelope = J::Obj(vec![
                    ("worker_count".into(), J::Int(1)),
                    ("rayon_threads".into(), J::Int(1)),
                    ("isolation".into(), J::Str(k.isolation.into())),
                    (
                        "cpus_allowed_list".into(),
                        rep.at("resource_envelope")?
                            .at("cpus_allowed_list")?
                            .clone(),
                    ),
                    ("process".into(), J::Str(k.process.into())),
                ]);
                let seed = params
                    .at("seed")?
                    .as_i128()
                    .ok_or("invalid algorithm seed")?;
                let work = identity::workload_manifest(
                    curve,
                    fixture,
                    k.input_law,
                    seed,
                    &envelope,
                    "cold",
                )?;
                let cand_id = cand.at("candidate_id")?.as_str().unwrap_or("");
                let work_id = work.at("workload_id")?.as_str().unwrap_or("");
                let run_id = identity::run_id(cand_id, work_id, row.run.into())?;
                let rho_ref = rho_manifest(&curve_ids, a, n, &rep, &rho_hashes)?;
                let rho_uid = identity::ec1_candidate_identity(&rho_ref)?;
                let rho_sha = rho_uid.at("candidate_sha256")?.as_str().unwrap_or("");
                identity::write_immutable(
                    &manifests.join("candidates").join(format!("{cand_id}.json")),
                    cand,
                )?;
                identity::write_immutable(
                    &manifests.join("workloads").join(format!("{work_id}.json")),
                    &work,
                )?;
                let J::Obj(mut rho_doc) = rho_uid.clone() else {
                    unreachable!("an identity is an object")
                };
                rho_doc.push(("manifest".into(), rho_ref));
                identity::write_immutable(
                    &manifests.join("rho").join(format!("{rho_sha}.json")),
                    &J::Obj(rho_doc),
                )?;
                let replays = J::Obj(vec![
                    ("ic".into(), replay(curve, fixture, &rep, "ic")?),
                    ("rho".into(), replay(curve, fixture, &rep, "rho")?),
                ]);
                let stem = format!("k{a}n{n}/T{:02}-R{}", row.i, row.run);
                let replay_path = claims_dir.join(format!("{stem}.replay.json"));
                write(&replay_path, json::dumps(&replays, 1) + "\n")?;
                let claim = self.claim(
                    (a, n, row.i),
                    &rep,
                    &iso,
                    &params,
                    [cand, &work],
                    &run_id,
                    &rho_uid,
                    &envelope,
                    &replays,
                    &replay_path,
                    &path,
                    row.tried,
                    &curve_ids,
                    [&host, &ic_hashes, &rho_hashes],
                    (&commit, &binary_sha),
                    fixture,
                )?;
                let result = autolab::validate_claim(&claim, "vs_rho", &ledger)?;
                let doc = J::Obj(vec![
                    ("claim".into(), claim.clone()),
                    ("validation".into(), result.clone()),
                ]);
                write(
                    &claims_dir.join(format!("{stem}.claim.json")),
                    json::dumps(&doc, 1) + "\n",
                )?;
                let mut missing = result
                    .at("missing_stage_fields")?
                    .as_arr()
                    .unwrap_or(&[])
                    .to_vec();
                missing.extend(
                    result
                        .at("missing_global_provenance")?
                        .as_arr()
                        .unwrap_or(&[])
                        .iter()
                        .cloned(),
                );
                summary.push(J::Obj(vec![
                    ("row".into(), J::Str(label)),
                    ("run_id".into(), J::Str(run_id)),
                    ("status".into(), result.at("status")?.clone()),
                    ("errors".into(), result.at("validation_errors")?.clone()),
                    ("missing".into(), J::Arr(missing)),
                    (
                        "independent_validation".into(),
                        claim.at("independent_validation")?.clone(),
                    ),
                ]));
            }
        }
        let passed = summary
            .iter()
            .filter(|s| s.get("status").and_then(J::as_str) == Some("PASS"))
            .count();
        let total = summary.len();
        write(
            &claims_dir.join("summary.json"),
            json::dumps(&J::Arr(summary), 1) + "\n",
        )?;
        Ok(format!("{passed} of {total} rows pass the vs_rho check"))
    }

    #[allow(clippy::too_many_arguments)]
    fn claim(
        &self,
        (a, n, i): (u8, u32, u32),
        rep: &J,
        iso: &J,
        params: &J,
        [cand, work]: [&J; 2],
        run_id: &str,
        rho_uid: &J,
        envelope: &J,
        replays: &J,
        replay_path: &Path,
        report_path: &Path,
        tried: Vec<J>,
        curve_ids: &J,
        [host, ic_hashes, rho_hashes]: [&J; 3],
        (commit, binary_sha): (&str, &J),
        fixture: &J,
    ) -> Result<J, String> {
        let k = self.constants;
        let m = rep.at("median")?;
        let all_verified = |arm: &str| -> Result<bool, String> {
            let reps = rep.at("repetitions")?.as_arr().ok_or("no repetitions")?;
            reps.iter().try_fold(
                true,
                |acc, r| Ok(acc && r.at(arm)?.at("verified")?.truthy()),
            )
        };
        let ic_ok = all_verified("ic_online")?;
        let rho_ok = all_verified("rho_online")?;
        let ic_ms = num(m.at("ic_online_wall_ms")?, "ic_online_wall_ms")?;
        let rho_ms = num(m.at("rho_online_wall_ms")?, "rho_online_wall_ms")?;
        let speedup = rho_ms / ic_ms;
        let independent = ["ic", "rho"].iter().try_fold(true, |acc, arm| {
            let r = replays.at(arm)?;
            Ok::<bool, String>(
                acc && r.at("digest_matches")?.truthy() && r.at("statement_holds")?.truthy(),
            )
        })?;
        let curve_id = curve_ids
            .at("curves")?
            .as_arr()
            .unwrap_or(&[])
            .iter()
            .find(|c| is_size(c, a, n))
            .ok_or("curve_ids.json has no such curve")?
            .at("curve_id")?
            .clone();
        let certs = rep.at("certificates")?;
        let policy = rep.at("rho_policy")?;
        let pick = |src: &J, keys: &[&str]| -> Result<J, String> {
            Ok(J::Obj(
                keys.iter()
                    .map(|key| Ok((key.to_string(), src.at(key)?.clone())))
                    .collect::<Result<_, String>>()?,
            ))
        };
        let target_hash = rep.at("target")?.at("sha256")?.clone();
        let curve = rep.at("curve")?;
        Ok(J::Obj(vec![
            ("n".into(), J::Int(n.into())),
            ("a".into(), J::Int(a.into())),
            ("log2_r".into(), rep.at("log2_r")?.clone()),
            ("curve".into(), curve.clone()),
            ("curve_id".into(), curve_id),
            ("candidate_id".into(), cand.at("candidate_id")?.clone()),
            (
                "candidate_manifest_sha256".into(),
                cand.at("record_sha256")?.clone(),
            ),
            ("workload_id".into(), work.at("workload_id")?.clone()),
            (
                "workload_manifest_sha256".into(),
                work.at("record_sha256")?.clone(),
            ),
            ("run_id".into(), J::Str(run_id.into())),
            (
                "rho_reference_uid".into(),
                rho_uid.at("candidate_uid")?.clone(),
            ),
            ("target_count".into(), J::Int(1)),
            ("ic_target_hash".into(), target_hash.clone()),
            ("rho_target_hash".into(), target_hash),
            (
                "timing_class".into(),
                J::Str("single_target_online_wall".into()),
            ),
            (
                "ic_online_wall_ms".into(),
                m.at("ic_online_wall_ms")?.clone(),
            ),
            (
                "rho_online_wall_ms".into(),
                m.at("rho_online_wall_ms")?.clone(),
            ),
            ("online_speedup".into(), J::Float(speedup)),
            (
                "ic_online_phase_ms".into(),
                m.at("ic_online_phase_ms")?.clone(),
            ),
            (
                "ic_online_phases_not_entered".into(),
                m.at("ic_online_phases_not_entered")?.clone(),
            ),
            ("online_interval".into(), rep.at("online_interval")?.clone()),
            ("same_resource_envelope".into(), J::Bool(true)),
            ("independent_validation".into(), J::Bool(independent)),
            (
                "ic_replay_certificate_sha256".into(),
                certs.at("ic_sha256")?.clone(),
            ),
            (
                "rho_replay_certificate_sha256".into(),
                certs.at("rho_sha256")?.clone(),
            ),
            ("ic_resource_envelope".into(), envelope.clone()),
            ("rho_resource_envelope".into(), envelope.clone()),
            ("ic_scalar_verified".into(), J::Bool(ic_ok)),
            ("rho_scalar_verified".into(), J::Bool(rho_ok)),
            (
                "rho_policy".into(),
                pick(
                    policy,
                    &[
                        "worker_count",
                        "walk_policy",
                        "collision_policy",
                        "distinguished_point_memory_bytes",
                        "lanes",
                        "distinguished_point_bits",
                    ],
                )?,
            ),
            (
                "verdict".into(),
                J::Str(
                    if speedup > 1.0 {
                        "index calculus online faster"
                    } else {
                        "rho online faster"
                    }
                    .into(),
                ),
            ),
            (
                "claim_boundary".into(),
                J::Str(format!(
                    "one public target on {}; the index calculus's reusable set-up is outside its \
                     online interval and is reported beside it",
                    py_str(curve)
                )),
            ),
            (
                "independent_replay_pointer".into(),
                J::Str(self.rel(replay_path)),
            ),
            ("fixture_hash".into(), J::Str(identity::sha256(fixture)?)),
            (
                "executable_or_source_hash".into(),
                J::Obj(vec![
                    ("ic_binary_sha256".into(), binary_sha.clone()),
                    ("source_commit".into(), J::Str(commit.into())),
                    (
                        "ic_source_manifest_sha256".into(),
                        J::Str(identity::sha256(ic_hashes)?),
                    ),
                    (
                        "rho_source_manifest_sha256".into(),
                        J::Str(identity::sha256(rho_hashes)?),
                    ),
                ]),
            ),
            (
                "host_id".into(),
                J::Obj(vec![
                    ("cpu_model".into(), host.at("cpu_model")?.clone()),
                    ("logical_cores".into(), host.at("logical_cores")?.clone()),
                    ("os".into(), host.at("os")?.clone()),
                    ("host_manifest".into(), J::Str("runs/host.json".into())),
                ]),
            ),
            (
                "resource_caps".into(),
                J::Obj(vec![
                    ("rayon_threads".into(), J::Int(1)),
                    ("cpus".into(), J::Str("2".into())),
                    ("contended".into(), iso.at("contended")?.clone()),
                ]),
            ),
            (
                "seeds".into(),
                J::Obj(vec![
                    ("recipe_seed".into(), params.at("seed")?.clone()),
                    (
                        "public_hash_seed".into(),
                        J::Int(TARGET_SEED_BASE + i128::from(i)),
                    ),
                    ("rho_seed".into(), policy.at("seed")?.clone()),
                ]),
            ),
            (
                "claim_boundary_non_claims".into(),
                J::Arr(k.non_claims.iter().map(|s| J::Str(s.to_string())).collect()),
            ),
            (
                "units".into(),
                pick(
                    m,
                    &[
                        "ic_online_units",
                        "rho_online_units",
                        "setup_units",
                        "rho_setup_units",
                        "rho_model_units",
                        "s_ic_online",
                        "s_rho_online",
                        "s_setup",
                        "s_rho_setup",
                        "s_ic_cold",
                        "s_rho_cold",
                        "s_rho_model",
                        "cold_ratio_ic_over_rho",
                        "online_speedup_rho_model",
                        "rho_units_per_step",
                        "rho_step_over_model",
                        "ic_replay_ns",
                        "rho_replay_ns",
                    ],
                )?,
            ),
            ("report".into(), J::Str(self.rel(report_path))),
            ("attempts".into(), J::Arr(tried)),
        ]))
    }
}

fn is_size(c: &J, a: u8, n: u32) -> bool {
    c.get("a")
        .is_some_and(|v| json::py_eq(v, &J::Int(a.into())))
        && c.get("n")
            .is_some_and(|v| json::py_eq(v, &J::Int(n.into())))
}

/// The fields the checker's curve is a function of, as one key.
fn curve_key(fixture: &J) -> Result<String, String> {
    const FIELDS: [&str; 8] = [
        "degree",
        "curve_a",
        "subgroup_order",
        "irreducible",
        "group_order",
        "cofactor",
        "generator",
        "lambda",
    ];
    let kv = FIELDS
        .iter()
        .map(|k| Ok((k.to_string(), fixture.at(k)?.clone())))
        .collect::<Result<Vec<_>, String>>()?;
    Ok(json::dumps_compact(&J::Obj(kv), true, true))
}

/// The method the index calculus ran: §20's recipe, stage by stage.
fn method_record(params: &J, fixture: &J, hashes: &J) -> Result<J, String> {
    let lib = hashes
        .at("src/cryptanalysis/koblitz_index_calculus.rs")?
        .clone();
    let spec = params.at("factor_base")?.at("spec")?;
    let m = params
        .get("descent_summands")
        .unwrap_or(params.at("summands")?);
    let max_trials = params.at("max_trials")?;
    let order = match fixture.at("subgroup_order")? {
        J::Int(i) => *i,
        J::Str(s) => s.trim().parse().map_err(|_| "subgroup_order")?,
        _ => return Err("subgroup_order".into()),
    };
    let s = |x: &str| J::Str(x.into());
    let components = match hashes {
        J::Obj(kv) => {
            let mut kv = kv.clone();
            kv.sort_by(|x, y| x.0.cmp(&y.0));
            kv.into_iter()
                .map(|(p, h)| J::Obj(vec![("role".into(), J::Str(p)), ("sha256".into(), h)]))
                .collect()
        }
        _ => Vec::new(),
    };
    Ok(J::Obj(vec![
        ("isogeny".into(), s("none")),
        (
            "endomorphism".into(),
            json::obj([
                ("order_conductor", J::Null),
                ("frobenius_order_conductor", J::Null),
                ("volcano_levels", J::Arr(vec![])),
            ]),
        ),
        (
            "factor_base".into(),
            json::obj([
                ("construction", spec.clone()),
                ("nominal_bound", spec.at("points")?.clone()),
            ]),
        ),
        (
            "point_decomposition".into(),
            json::obj([
                ("summands", params.at("summands")?.clone()),
                ("solver", s("pairtable")),
                (
                    "summation_polynomial",
                    s("none: the pair-sum table is looked up"),
                ),
                (
                    "encoding",
                    s(
                        "pair sums folded by sign and Frobenius, keyed by the normal-basis least \
                       rotation of the abscissa",
                    ),
                ),
                ("equation_order", s("not applicable: no polynomial system")),
                ("monomial_order", s("not applicable: no polynomial system")),
                (
                    "internal_matrix_kernel",
                    s("not applicable: no polynomial system"),
                ),
                (
                    "limits",
                    json::obj([
                        ("max_trials", max_trials.clone()),
                        ("table_tier", s("auto, the workflow's probe budget")),
                    ]),
                ),
                ("cache_policy", s("rebuilt from nothing in every process")),
                ("source_sha256", lib.clone()),
            ]),
        ),
        (
            "relation_collection".into(),
            json::obj([
                ("collector", s("aimed")),
                (
                    "query_distribution",
                    s("[a]G + [b]Q-free probes from seeded work units"),
                ),
                (
                    "query_rule",
                    s("aimed at the least-mentioned projected columns"),
                ),
                (
                    "filtering",
                    J::Str(format!(
                        "collection window {}",
                        py_str(params.at("collection_window")?)
                    )),
                ),
                (
                    "verification",
                    s("every relation checked in the group before it is pushed"),
                ),
                ("duplicates", s("dropped by the log solver")),
                ("dependencies", s("rank of the sparse system")),
                (
                    "stop_rule",
                    s("extend by units until every column is determined"),
                ),
                ("source_sha256", lib.clone()),
            ]),
        ),
        (
            "relation_linear_algebra".into(),
            json::obj([
                ("solver", s("sparse")),
                ("modulus", J::Int(order)),
                (
                    "matrix_construction",
                    s("one column per projected signed Frobenius orbit"),
                ),
                ("orbit_quotient", s("sign-and-Frobenius")),
                ("rank_criterion", s("every column determined and verified")),
                ("block_parameters", s("library defaults")),
                ("preconditioner", s("library defaults")),
                ("source_sha256", lib.clone()),
            ]),
        ),
        (
            "target_descent".into(),
            json::obj([
                ("method", s("walk")),
                (
                    "policy",
                    J::Str(format!(
                        "{}-summand pair-table walk of 64 lanes stepped by G, one start by two \
                         scalar multiplications",
                        py_str(m)
                    )),
                ),
                ("recursive_solvers", s("none")),
                ("success_rule", s("[d]G = Q on the single-word ladder")),
                (
                    "stop_rule",
                    J::Str(format!("max_trials {}", py_str(max_trials))),
                ),
                ("source_sha256", lib),
            ]),
        ),
        (
            "implementation".into(),
            json::obj([
                ("source_manifest_sha256", J::Str(identity::sha256(hashes)?)),
                ("components", J::Arr(components)),
                (
                    "flags",
                    json::obj([
                        ("rayon_threads", J::Int(1)),
                        ("single_target", J::Bool(true)),
                    ]),
                ),
            ]),
        ),
    ]))
}

/// Rho's reference identity: the curve's EC1 record and the walk it ran.
fn rho_manifest(curve_ids: &J, a: u8, n: u32, rep: &J, hashes: &J) -> Result<J, String> {
    let rec = curve_ids
        .at("curves")?
        .as_arr()
        .unwrap_or(&[])
        .iter()
        .find(|c| is_size(c, a, n))
        .ok_or("curve_ids.json has no such curve")?;
    let policy = rep.at("rho_policy")?;
    let J::Obj(mut curve) = rec.at("curve")?.clone() else {
        return Err("curve_ids.json: a curve record is not an object".into());
    };
    // `{**rec["curve"], "curve_id": ...}`: a present key keeps its place.
    let id = rec.at("curve_id")?.clone();
    match curve.iter_mut().find(|(k, _)| k == "curve_id") {
        Some((_, v)) => *v = id,
        None => curve.push(("curve_id".into(), id)),
    }
    let components = match hashes {
        J::Obj(kv) => {
            let mut kv = kv.clone();
            kv.sort_by(|x, y| x.0.cmp(&y.0));
            kv.into_iter()
                .map(|(p, h)| J::Obj(vec![("role".into(), J::Str(p)), ("sha256".into(), h)]))
                .collect()
        }
        _ => Vec::new(),
    };
    Ok(J::Obj(vec![
        ("field".into(), rec.at("field")?.clone()),
        ("curve".into(), J::Obj(curve)),
        ("method".into(), J::Str("pollard-rho".into())),
        ("factor_base".into(), J::Str("none".into())),
        ("isogeny".into(), J::Str("none".into())),
        (
            "configuration".into(),
            J::Obj(vec![
                ("walk".into(), policy.at("walk_policy")?.clone()),
                ("collision".into(), policy.at("collision_policy")?.clone()),
                ("lanes".into(), policy.at("lanes")?.clone()),
                (
                    "distinguished_point_bits".into(),
                    policy.at("distinguished_point_bits")?.clone(),
                ),
                ("automorphisms".into(), J::Int(2 * i128::from(n))),
                ("jumps".into(), policy.at("counters")?.at("jumps")?.clone()),
            ]),
        ),
        (
            "implementation".into(),
            J::Obj(vec![
                (
                    "source_manifest_sha256".into(),
                    J::Str(identity::sha256(hashes)?),
                ),
                ("components".into(), J::Arr(components)),
            ]),
        ),
    ]))
}

/// One arm's certificate, replayed: its digest recomputed from its canonical
/// JSON, and `[scalar]G = target` checked in [`oracle`]'s arithmetic.
fn replay(curve: &Curve, fixture: &J, rep: &J, arm: &str) -> Result<J, String> {
    let certs = rep.at("certificates")?;
    let cert = certs.at(arm)?;
    let digest = identity::sha256_hex(json::dumps_compact(cert, true, false).as_bytes());
    let target = curve.decode(cert.at("target")?)?;
    let scalar = match cert.at("scalar")? {
        J::Int(i) => *i,
        J::Str(s) => s
            .trim()
            .parse::<i128>()
            .map_err(|_| format!("{arm} certificate: the scalar is not an integer"))?,
        _ => return Err(format!("{arm} certificate: the scalar is not an integer")),
    };
    let first = fixture
        .at("targets")?
        .as_arr()
        .and_then(<[J]>::first)
        .ok_or("the fixture has no target")?;
    let holds = 0 <= scalar
        && scalar < i128::from(curve.r)
        && curve.mul(Some(curve.g), scalar as u128) == target
        && json::py_eq(cert.at("target")?, first);
    let recorded = certs.at(&format!("{arm}_sha256"))?;
    Ok(J::Obj(vec![
        ("arm".into(), J::Str(arm.into())),
        ("digest_recomputed".into(), J::Str(digest.clone())),
        (
            "digest_matches".into(),
            J::Bool(recorded.as_str() == Some(digest.as_str())),
        ),
        ("statement_holds".into(), J::Bool(holds)),
        ("scalar".into(), cert.at("scalar")?.clone()),
        ("checker".into(), J::Str(CHECKER.into())),
    ]))
}

// ── the analysis ────────────────────────────────────────────────────

/// `x ** 2` as CPython's `float_pow` computes it: the platform's `pow`,
/// which need not be the product `x * x` to the last bit.
fn pow2(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    // The exponent is hidden from the optimiser, which would otherwise
    // fold `pow(x, 2)` into `x * x`.
    x.abs().powf(std::hint::black_box(2.0))
}

/// Python's built-in `sum` of floats before 3.12: left to right from 0.
fn py_sum(xs: impl IntoIterator<Item = f64>) -> f64 {
    xs.into_iter().fold(0.0, |acc, x| acc + x)
}

/// `statistics.median` of JSON numbers: the middle one as it is, or the
/// mean of the two middle ones as a float.
fn py_median(mut v: Vec<J>) -> Result<J, String> {
    if v.is_empty() {
        return Err("no median for empty data".into());
    }
    for x in &v {
        num(x, "a median's datum")?;
    }
    v.sort_by(|x, y| {
        x.as_f64()
            .unwrap()
            .partial_cmp(&y.as_f64().unwrap())
            .expect("no NaN in the records")
    });
    let n = v.len();
    if n % 2 == 1 {
        return Ok(v[n / 2].clone());
    }
    Ok(J::Float(match (&v[n / 2 - 1], &v[n / 2]) {
        (J::Int(a), J::Int(b)) => (a + b) as f64 / 2.0,
        (a, b) => (a.as_f64().unwrap() + b.as_f64().unwrap()) / 2.0,
    }))
}

/// `round(x, 2)`: the correctly rounded decimal, read back.
fn round2(x: f64) -> f64 {
    format!("{x:.2}").parse().expect("a decimal")
}

/// The slope of `y` on `x` with its 95% `t` interval; `null` below two
/// points.
fn fit(points: &[(f64, f64)]) -> J {
    let k = points.len();
    if k < 2 {
        return J::Null;
    }
    let xs: Vec<f64> = points.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.1).collect();
    let (mx, my) = (stats::fmean(&xs), stats::fmean(&ys));
    let sxx = py_sum(xs.iter().map(|x| pow2(x - mx)));
    let slope = py_sum(points.iter().map(|(x, y)| (x - mx) * (y - my))) / sxx;
    let icept = my - slope * mx;
    if k < 3 {
        return J::Obj(vec![
            ("slope".into(), J::Float(slope)),
            ("points".into(), J::Int(k as i128)),
        ]);
    }
    let resid = py_sum(points.iter().map(|(x, y)| pow2(y - icept - slope * x)));
    let se = (resid / (k - 2) as f64 / sxx).sqrt();
    let t = T975.get(k - 3).copied().unwrap_or(1.96);
    J::Obj(vec![
        ("slope".into(), J::Float(slope)),
        (
            "ci95".into(),
            J::Arr(vec![J::Float(slope - t * se), J::Float(slope + t * se)]),
        ),
        ("points".into(), J::Int(k as i128)),
    ])
}

/// A checked row, with its target and run.
struct Checked {
    claim: J,
    target: u32,
    run: u32,
}

impl Checked {
    fn unit(&self, key: &str) -> f64 {
        self.claim
            .get("units")
            .and_then(|u| u.get(key))
            .and_then(J::as_f64)
            .unwrap_or(f64::NAN)
    }
}

fn mean(rows: &[&Checked], key: &str) -> f64 {
    let v: Vec<f64> = rows.iter().map(|r| r.unit(key)).collect();
    stats::fmean(&v)
}

/// The 2.5% and 97.5% points of `stat` over `BOOT` resamples of `rows`.
fn boot_ci(
    rows: &[Checked],
    stat: impl Fn(&[&Checked]) -> Option<f64>,
    seed: u128,
) -> Result<J, String> {
    let mut rng = Random::new(seed);
    let len = u32::try_from(rows.len()).map_err(|_| "too many rows")?;
    let mut vals = Vec::with_capacity(BOOT);
    let mut sample: Vec<&Checked> = Vec::with_capacity(rows.len());
    for _ in 0..BOOT {
        sample.clear();
        for _ in 0..rows.len() {
            sample.push(&rows[rng.randbelow(len) as usize]);
        }
        if let Some(v) = stat(&sample) {
            if v.is_finite() {
                vals.push(v);
            }
        }
    }
    if vals.is_empty() {
        return Err("no finite resample".into());
    }
    vals.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let lo = (0.025 * vals.len() as f64) as usize;
    let hi = (0.975 * vals.len() as f64) as usize - 1;
    Ok(J::Arr(vec![J::Float(vals[lo]), J::Float(vals[hi])]))
}

impl Comparison {
    /// `(R1 rows, R2 rows, rows that did not pass)`.
    fn load_rows(&self, a: u8, n: u32) -> Result<(Vec<Checked>, Vec<Checked>, Vec<J>), String> {
        let d = self.out.join("claims").join(format!("k{a}n{n}"));
        let (mut r1, mut r2, mut bad) = (Vec::new(), Vec::new(), Vec::new());
        for name in names(&d)
            .into_iter()
            .filter(|x| glob2(x, "T", "-R", ".claim.json"))
        {
            let doc = read(&d.join(&name))?;
            let validation = doc.at("validation")?;
            let ok = validation.at("status")?.as_str() == Some("PASS");
            let (target, run) = target_run(&name, ".claim.json")
                .ok_or_else(|| format!("{name} is not a claim row"))?;
            let claim = doc.at("claim")?.clone();
            if !(ok && claim.at("independent_validation")?.truthy()) {
                bad.push(J::Obj(vec![
                    ("row".into(), J::Str(name)),
                    ("status".into(), validation.at("status")?.clone()),
                    ("errors".into(), validation.at("validation_errors")?.clone()),
                ]));
                continue;
            }
            let row = Checked { claim, target, run };
            if run == 1 {
                r1.push(row);
            } else {
                r2.push(row);
            }
        }
        Ok((r1, r2, bad))
    }

    fn report_of(&self, claim: &J) -> Result<J, String> {
        let p = PathBuf::from(claim.at("report")?.as_str().ok_or("report")?);
        read(&if p.is_absolute() {
            p
        } else {
            self.root.join(p)
        })
    }

    fn isolation(&self, a: u8, n: u32) -> Result<J, String> {
        let d = self.runs.join(format!("k{a}n{n}"));
        let (mut runs, mut contended, mut failed, mut retries) = (0i128, 0i128, 0i128, 0i128);
        for name in names(&d) {
            if glob2(&name, "", "-retry", ".price.json") {
                retries += 1;
            }
            if !name.ends_with(".isolation.jsonl") {
                continue;
            }
            let path = d.join(&name);
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            for line in text.lines() {
                let run = json::parse(line)?.at("run")?.clone();
                runs += 1;
                contended += i128::from(run.at("contended")?.truthy());
                failed += i128::from(!json::py_eq(run.at("exit_status")?, &J::Int(0)));
            }
        }
        Ok(J::Obj(vec![
            ("processes".into(), J::Int(runs)),
            ("contended".into(), J::Int(contended)),
            ("failed".into(), J::Int(failed)),
            ("retries".into(), J::Int(retries)),
        ]))
    }

    fn size(&self, a: u8, n: u32, seed: u128) -> Result<Option<J>, String> {
        let (r1, r2, bad) = self.load_rows(a, n)?;
        if r1.is_empty() {
            return Ok(None);
        }
        let reports: Vec<J> = r1
            .iter()
            .map(|c| self.report_of(&c.claim))
            .collect::<Result<_, _>>()?;
        let r_json = reports[0].at("r")?.clone();
        let r = num(&r_json, "r")?;
        let sqrt_r = r.sqrt();
        let canonical_json = py_median(
            reports
                .iter()
                .map(|rep| Ok(rep.at("step_prices")?.at("canonical_step_units")?.clone()))
                .collect::<Result<_, String>>()?,
        )?;
        let canonical = num(&canonical_json, "canonical_step_units")?;
        let first_rep = |rep: &J| -> Result<J, String> {
            Ok(rep
                .at("repetitions")?
                .as_arr()
                .and_then(<[J]>::first)
                .ok_or("no repetitions")?
                .clone())
        };
        let mut steps = Vec::new();
        for rep in &reports {
            steps.push(num(
                first_rep(rep)?.at("rho_online")?.at("steps")?,
                "rho steps",
            )?);
        }
        let nf = f64::from(n);
        let floor_steps = (std::f64::consts::PI * r / (4.0 * nf)).sqrt();
        let unit_median = |key: &str| -> Result<J, String> {
            py_median(
                r1.iter()
                    .map(|c| Ok(c.claim.at("units")?.at(key)?.clone()))
                    .collect::<Result<_, String>>()?,
            )
        };
        let setup_json = unit_median("setup_units")?;
        let setup = num(&setup_json, "setup_units")?;
        let rho_setup = num(&unit_median("rho_setup_units")?, "rho_setup_units")?;

        let speedup = |rows: &[&Checked]| {
            Some(mean(rows, "rho_online_units") / mean(rows, "ic_online_units"))
        };
        let speedup_model =
            |rows: &[&Checked]| Some(mean(rows, "rho_model_units") / mean(rows, "ic_online_units"));
        let cold = |rows: &[&Checked]| {
            Some(
                (setup + mean(rows, "ic_online_units"))
                    / (rho_setup + mean(rows, "rho_online_units")),
            )
        };
        let break_even = |rows: &[&Checked]| {
            let saving = mean(rows, "rho_online_units") - mean(rows, "ic_online_units");
            (saving > 0.0).then(|| setup / saving)
        };
        let all: Vec<&Checked> = r1.iter().collect();
        let bl_online = BL_PRODUCT * r * pow2(canonical) / (2.0 * nf * setup);
        let per_row: Vec<J> = r1
            .iter()
            .map(|c| Ok(c.claim.at("online_speedup")?.clone()))
            .collect::<Result<_, String>>()?;
        let above_one = per_row
            .iter()
            .filter(|s| s.as_f64().is_some_and(|x| x > 1.0))
            .count();
        let mut aa = Vec::new();
        for c2 in &r2 {
            if let Some(c1) = r1.iter().find(|c| c.target == c2.target) {
                let ratio = |key: &str| -> Result<f64, String> {
                    Ok(num(c2.claim.at(key)?, key)? / num(c1.claim.at(key)?, key)?)
                };
                aa.push(J::Obj(vec![
                    ("target".into(), J::Int(c2.target.into())),
                    (
                        "online_speedup_R2_over_R1".into(),
                        J::Float(ratio("online_speedup")?),
                    ),
                    (
                        "ic_online_R2_over_R1".into(),
                        J::Float(c2.unit("ic_online_units") / c1.unit("ic_online_units")),
                    ),
                    (
                        "rho_online_R2_over_R1".into(),
                        J::Float(c2.unit("rho_online_units") / c1.unit("rho_online_units")),
                    ),
                ]));
            }
        }
        let break_even_value = break_even(&all).map_or(J::Null, J::Float);
        let f = J::Float;
        let first = &r1[0].claim;
        let diagnostics = self.diagnostics(a, n, &all, &reports, setup, rho_setup)?;
        let decisions = self.decisions(a, n, &r1, &r2)?;
        let reference = self.reference_check(a, n, &r1)?;
        let mut out = vec![
            ("curve".into(), first.at("curve")?.clone()),
            ("a".into(), J::Int(a.into())),
            ("n".into(), J::Int(n.into())),
            ("r".into(), r_json),
            ("log2_r".into(), f(round2(r.log2()))),
            ("curve_id".into(), first.at("curve_id")?.clone()),
            ("candidate_id".into(), first.at("candidate_id")?.clone()),
            (
                "rho_reference_uid".into(),
                first.at("rho_reference_uid")?.clone(),
            ),
            (
                "rows".into(),
                J::Obj(vec![
                    ("R1".into(), J::Int(r1.len() as i128)),
                    ("R2".into(), J::Int(r2.len() as i128)),
                    ("not_passing".into(), J::Arr(bad.clone())),
                ]),
            ),
            ("all_rows_pass_the_check".into(), J::Bool(bad.is_empty())),
            (
                "online_speedup".into(),
                J::Obj(vec![
                    ("mean_ratio".into(), f(speedup(&all).expect("a ratio"))),
                    ("ci95".into(), boot_ci(&r1, speedup, seed)?),
                    ("median_of_rows".into(), py_median(per_row)?),
                    ("rows_above_one".into(), J::Int(above_one as i128)),
                ]),
            ),
            (
                "online_speedup_rho_model".into(),
                J::Obj(vec![
                    (
                        "mean_ratio".into(),
                        f(speedup_model(&all).expect("a ratio")),
                    ),
                    ("ci95".into(), boot_ci(&r1, speedup_model, seed + 1)?),
                ]),
            ),
            (
                "cold_ratio_ic_over_rho".into(),
                J::Obj(vec![
                    ("value".into(), f(cold(&all).expect("a ratio"))),
                    ("ci95".into(), boot_ci(&r1, cold, seed + 2)?),
                ]),
            ),
            (
                "break_even_targets".into(),
                J::Obj(vec![
                    ("value".into(), break_even_value),
                    ("ci95".into(), boot_ci(&r1, break_even, seed + 3)?),
                ]),
            ),
            (
                "s_ic_online_mean".into(),
                f(mean(&all, "ic_online_units") / sqrt_r),
            ),
            (
                "s_rho_online_mean".into(),
                f(mean(&all, "rho_online_units") / sqrt_r),
            ),
            (
                "s_rho_model_mean".into(),
                f(mean(&all, "rho_model_units") / sqrt_r),
            ),
            ("s_setup".into(), f(setup / sqrt_r)),
            ("s_rho_setup".into(), f(rho_setup / sqrt_r)),
            (
                "s_floor_steps".into(),
                f((std::f64::consts::PI / (4.0 * nf)).sqrt()),
            ),
            (
                "rho_steps_over_floor_mean".into(),
                f(stats::fmean(&steps) / floor_steps),
            ),
            (
                "rho_units_per_step_median".into(),
                unit_median("rho_units_per_step")?,
            ),
            ("canonical_step_units".into(), canonical_json),
            (
                "rho_step_over_model_median".into(),
                unit_median("rho_step_over_model")?,
            ),
            (
                "precomputation_boundary_model".into(),
                J::Obj(vec![
                    (
                        "what".into(),
                        J::Str(
                            "Bernstein-Lange: a generic walk given the index calculus's set-up \
                             as its budget, online ~ 1.93*1.21 r c^2 / (2n P) units; a model"
                                .into(),
                        ),
                    ),
                    ("bl_online_units".into(), f(bl_online)),
                    ("s_bl_online".into(), f(bl_online / sqrt_r)),
                    (
                        "ic_online_over_bl".into(),
                        f(mean(&all, "ic_online_units") / bl_online),
                    ),
                ]),
            ),
            ("ic_replay_ns_median".into(), unit_median("ic_replay_ns")?),
            ("rho_replay_ns_median".into(), unit_median("rho_replay_ns")?),
            ("aa".into(), J::Arr(aa)),
            ("isolation".into(), self.isolation(a, n)?),
            ("diagnostics".into(), diagnostics),
        ];
        // A later comparison also says whether each row recovered the earlier
        // one's logarithms; §23's analysis has no such key.
        if let Some(d) = decisions {
            out.insert(10, ("logarithms_as_earlier".into(), d));
        }
        if let Some(d) = reference {
            out.insert(11, ("reference_check".into(), d));
        }
        Ok(Some(J::Obj(out)))
    }

    /// Read-outs beside the declared figures: where the online time goes,
    /// the curve's construction, the cold ratio without it and with rho at
    /// the canonical step, and the online interval against §22's
    /// per-target descent.
    fn diagnostics(
        &self,
        a: u8,
        n: u32,
        r1: &[&Checked],
        reports: &[J],
        setup: f64,
        rho_setup: f64,
    ) -> Result<J, String> {
        let mut phases: Vec<(String, Vec<f64>)> = Vec::new();
        for rep in reports {
            let at = rep
                .at("median")?
                .at("ic_online_repetition")?
                .as_i128()
                .ok_or("ic_online_repetition")?;
            let r = rep
                .at("repetitions")?
                .as_arr()
                .and_then(|v| v.get(at as usize))
                .ok_or("the median repetition")?;
            let unit = num(r.at("unit_ns")?, "unit_ns")?;
            for (k, v) in r
                .at("ic_online")?
                .at("phases_ns")?
                .as_obj()
                .ok_or("phases_ns")?
            {
                let x = (if v.truthy() { num(v, k)? } else { 0.0 }) / unit;
                match phases.iter_mut().find(|(name, _)| name == k) {
                    Some((_, list)) => list.push(x),
                    None => phases.push((k.clone(), vec![x])),
                }
            }
        }
        let online = mean(r1, "ic_online_units");
        let mut per_report = Vec::new();
        for rep in reports {
            let mut inner = Vec::new();
            for x in rep.at("repetitions")?.as_arr().ok_or("repetitions")? {
                let setup_ns = num(x.at("setup_phases_ns")?.at("setup")?, "setup")?;
                inner.push(J::Float(setup_ns / num(x.at("unit_ns")?, "unit_ns")?));
            }
            per_report.push(py_median(inner)?);
        }
        let curve = num(&py_median(per_report)?, "curve construction")?;
        let (rho_online, rho_model) = (mean(r1, "rho_online_units"), mean(r1, "rho_model_units"));
        let first_rep = |rep: &J| -> Result<J, String> {
            Ok(rep
                .at("repetitions")?
                .as_arr()
                .and_then(<[J]>::first)
                .ok_or("no repetitions")?
                .clone())
        };
        let mut trials_ours = Vec::new();
        for rep in reports {
            trials_ours.push(num(
                first_rep(rep)?.at("ic_online")?.at("trials")?,
                "trials",
            )?);
        }
        let online_trials_mean = stats::fmean(&trials_ours);
        let r = num(reports[0].at("r")?, "r")?;
        let mut out = vec![
            (
                "online_phase_shares".to_string(),
                J::Obj(
                    phases
                        .iter()
                        .map(|(k, v)| (k.clone(), J::Float(stats::fmean(v) / online)))
                        .collect(),
                ),
            ),
            ("curve_construction_units".into(), J::Float(curve)),
            ("s_curve_construction".into(), J::Float(curve / r.sqrt())),
            (
                "cold_ratio_without_curve_construction".into(),
                J::Float((setup - curve + online) / rho_online),
            ),
            (
                "cold_ratio_rho_at_canonical_step".into(),
                J::Float((setup + online) / (rho_setup + rho_model)),
            ),
            ("online_trials_mean".into(), J::Float(online_trials_mean)),
        ];
        // §22's isolated batch descent at the same size.
        let s22 = self
            .root
            .join("research/ic_descent_20260930/runs-isolated/main")
            .join(format!("k{a}n{n}"));
        let mut files = Vec::new();
        for m in names(&s22).into_iter().filter(|x| x.starts_with('M')) {
            let d = s22.join(&m);
            if !d.is_dir() {
                continue;
            }
            for f in names(&d)
                .into_iter()
                .filter(|x| glob2(x, "r", "", "-candidate.price.json"))
            {
                files.push(d.join(f));
            }
        }
        files.sort();
        let (mut trials, mut per_trial, mut per_target) = (Vec::new(), Vec::new(), Vec::new());
        for path in files {
            let d = read(&path)?;
            if d.get("status").and_then(J::as_str) != Some("complete") {
                continue;
            }
            let t: Vec<f64> = d
                .at("counts")?
                .at("descent")?
                .at("trials_per_target")?
                .as_arr()
                .ok_or("trials_per_target")?
                .iter()
                .map(|x| num(x, "trials"))
                .collect::<Result<_, _>>()?;
            let units = d.at("median")?.at("phases_units")?;
            let descent = num(units.at("descent")?, "descent")?;
            let verify = num(units.at("verify_final")?, "verify_final")?;
            let targets = num(d.at("targets")?, "targets")?;
            let t_mean = stats::fmean(&t);
            trials.push(t_mean);
            per_trial.push(J::Float(descent / targets / t_mean));
            per_target.push(J::Float((descent + verify) / targets));
        }
        if !trials.is_empty() {
            let mut units_sum = Vec::new();
            let mut trials_sum = 0i128;
            for rep in reports {
                units_sum.push(num(
                    rep.at("median")?.at("ic_online_units")?,
                    "ic_online_units",
                )?);
                trials_sum += first_rep(rep)?
                    .at("ic_online")?
                    .at("trials")?
                    .as_i128()
                    .ok_or("trials")?;
            }
            let ours = py_sum(units_sum) / trials_sum as f64;
            out.push((
                "against_s22_batch_descent".into(),
                J::Obj(vec![
                    ("s22_files".into(), J::Int(trials.len() as i128)),
                    (
                        "online_over_s22_per_target_descent".into(),
                        J::Float(online / num(&py_median(per_target)?, "median")?),
                    ),
                    (
                        "trials_ratio".into(),
                        J::Float(online_trials_mean / stats::fmean(&trials)),
                    ),
                    (
                        "units_per_trial_ratio".into(),
                        J::Float(ours / num(&py_median(per_trial)?, "median")?),
                    ),
                ]),
            ));
        }
        Ok(J::Obj(out))
    }

    pub fn analyse(&self) -> Result<J, String> {
        let mut sizes = Vec::new();
        for (k, (a, n)) in SIZES.iter().enumerate() {
            if let Some(s) = self.size(*a, *n, BOOT_SEED + 10 * k as u128)? {
                sizes.push(s);
            }
        }
        let ln = |s: &J, path: &[&str]| -> Result<(f64, f64), String> {
            let mut v = s;
            for p in path {
                v = v.at(p)?;
            }
            Ok((num(s.at("r")?, "r")?.ln(), num(v, path[0])?.ln()))
        };
        let series = |path: &[&str]| -> Result<J, String> {
            let points = sizes
                .iter()
                .map(|s| ln(s, path))
                .collect::<Result<Vec<_>, _>>()?;
            Ok(fit(&points))
        };
        let mut fits = Vec::new();
        for key in [
            "s_ic_online_mean",
            "s_rho_online_mean",
            "s_setup",
            "s_rho_model_mean",
        ] {
            fits.push((key.to_string(), series(&[key])?));
        }
        fits.push((
            "online_speedup".into(),
            series(&["online_speedup", "mean_ratio"])?,
        ));
        fits.push((
            "cold_ratio".into(),
            series(&["cold_ratio_ic_over_rho", "value"])?,
        ));
        let (against_key, against) = self.against(&sizes)?;
        Ok(J::Obj(vec![
            (
                "what_this_is".into(),
                J::Str(self.constants.what_this_is.into()),
            ),
            ("sizes".into(), J::Arr(sizes)),
            ("fits_ln_s_on_ln_r".into(), J::Obj(fits)),
            (against_key, J::Arr(against)),
        ]))
    }
}

/// The four figures a size is read against: its own path, and the key
/// §23's `prediction.json` gives the predicted value under.
const FIGURES: [(&str, [&str; 2], &str); 4] = [
    (
        "online_speedup",
        ["online_speedup", "mean_ratio"],
        "online_speedup_probe_step",
    ),
    (
        "online_speedup_model",
        ["online_speedup_rho_model", "mean_ratio"],
        "online_speedup_canonical_step",
    ),
    (
        "cold_ratio",
        ["cold_ratio_ic_over_rho", "value"],
        "cold_ratio_probe_step",
    ),
    (
        "ic_online_over_bl",
        ["precomputation_boundary_model", "ic_online_over_bl"],
        "ic_online_over_bl_model",
    ),
];

fn same_size(row: &J, s: &J) -> Result<bool, String> {
    let (a, n) = (s.at("a")?, s.at("n")?);
    Ok(row.get("a").is_some_and(|v| json::py_eq(v, a))
        && row.get("n").is_some_and(|v| json::py_eq(v, n)))
}

impl Comparison {
    /// Each size's figures beside what they are read against: §23's
    /// prediction, or an earlier comparison's measurement of the same
    /// figure.
    fn against(&self, sizes: &[J]) -> Result<(String, Vec<J>), String> {
        let (key, then, rows_key) = match &self.constants.against {
            Against::Prediction => (
                "prediction_then_measurement".to_string(),
                read(&self.here.join("prediction.json"))?,
                "rows",
            ),
            Against::Earlier { analysis, key } => {
                (key.to_string(), read(&self.root.join(analysis))?, "sizes")
            }
        };
        let mut out = Vec::new();
        for s in sizes {
            let mut earlier = None;
            for row in then.at(rows_key)?.as_arr().unwrap_or(&[]) {
                if same_size(row, s)? {
                    earlier = Some(row);
                    break;
                }
            }
            let earlier = earlier.ok_or("the earlier figures have no row for a size")?;
            let mut row = vec![("curve".to_string(), s.at("curve")?.clone())];
            for (name, [k1, k2], predicted) in FIGURES {
                let was = match self.constants.against {
                    Against::Prediction => earlier.at(predicted)?,
                    Against::Earlier { .. } => earlier.at(k1)?.at(k2)?,
                };
                row.push((
                    name.to_string(),
                    J::Arr(vec![was.clone(), s.at(k1)?.at(k2)?.clone()]),
                ));
            }
            out.push(J::Obj(row));
        }
        Ok((key, out))
    }

    /// Whether every checked row recovered the earlier comparison's
    /// logarithms, for both arms, on the same target.
    fn decisions(
        &self,
        a: u8,
        n: u32,
        r1: &[Checked],
        r2: &[Checked],
    ) -> Result<Option<J>, String> {
        let Some((claims, _)) = self.constants.earlier else {
            return Ok(None);
        };
        let replay_of = |claim: &J| -> Result<J, String> {
            let p = claim
                .at("independent_replay_pointer")?
                .as_str()
                .ok_or("replay pointer")?;
            read(&self.root.join(p))
        };
        let (mut equal, mut differing) = (0i128, Vec::new());
        for row in r1.iter().chain(r2) {
            let name = format!("T{:02}-R{}", row.target, row.run);
            let then_path = self
                .root
                .join(claims)
                .join(format!("k{a}n{n}/{name}.claim.json"));
            let then = read(&then_path)?;
            let then = then.at("claim")?;
            let (now_replay, then_replay) = (replay_of(&row.claim)?, replay_of(then)?);
            let mut same = row.claim.at("ic_target_hash")? == then.at("ic_target_hash")?;
            for arm in ["ic", "rho"] {
                same &= now_replay.at(arm)?.at("scalar")? == then_replay.at(arm)?.at("scalar")?;
            }
            if same {
                equal += 1;
            } else {
                differing.push(J::Str(name));
            }
        }
        Ok(Some(J::Obj(vec![
            ("rows".into(), J::Int((r1.len() + r2.len()) as i128)),
            ("equal".into(), J::Int(equal)),
            ("differing".into(), J::Arr(differing)),
        ])))
    }
}

// ── the runner (§23's `run.py`, ported) ─────────────────────────────

const RHO_SEED_BASE: i128 = 0x23_0000;
const TARGETS: u32 = 64;
const AA_TARGETS: u32 = 4;

/// What a comparison's runs price with.
pub struct Tools {
    pub ic: PathBuf,
    /// The commit `ic` was built from, for a fresh manifest.
    pub ic_commit: Option<String>,
    pub isolate: PathBuf,
}

fn uptime() -> String {
    Command::new("uptime")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

impl Comparison {
    fn refuse_frozen(&self) -> Result<(), String> {
        if self.constants.frozen {
            Err("this comparison's runs are frozen".into())
        } else {
            Ok(())
        }
    }

    /// The host and the binary the runs ran on, written once; a resumed
    /// run must price with the binary its manifest names.
    pub fn manifest(&self, tools: &Tools) -> Result<J, String> {
        self.refuse_frozen()?;
        let path = self.runs.join("host.json");
        let sha = bench::sha256_file(&tools.ic)?;
        if path.exists() {
            let host = read(&path)?;
            let (_, named) = Self::build(&host)?;
            if named.as_str() != Some(sha.as_str()) {
                return Err(format!(
                    "{} is not the binary host.json names",
                    tools.ic.display()
                ));
            }
            return Ok(host);
        }
        let commit = tools
            .ic_commit
            .as_ref()
            .ok_or("a fresh manifest needs the commit ic was built from")?;
        let ic = J::Obj(vec![
            (
                "path_basename".into(),
                J::Str(
                    tools
                        .ic
                        .file_name()
                        .map(|f| f.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                ),
            ),
            ("sha256".into(), J::Str(sha)),
            ("built_from".into(), J::Str(commit.clone())),
        ]);
        bench::host_manifest(
            &path,
            &self.root,
            J::Obj(vec![("ic".into(), ic)]),
            &tools.isolate,
            "the host the rule's comparison ran on, at its start",
        )
    }

    /// The parameter file of `T<i>` at a size: the earlier comparison's
    /// frozen file, copied once and the same bytes after.
    fn params(&self, a: u8, n: u32, i: u32) -> Result<PathBuf, String> {
        let (_, runs) = self
            .constants
            .earlier
            .ok_or("no earlier comparison to take the parameter files from")?;
        let name = format!("k{a}n{n}/T{i:02}.params.json");
        let (from, to) = (self.root.join(runs).join(&name), self.runs.join(&name));
        let bytes = std::fs::read(&from).map_err(|e| format!("{}: {e}", from.display()))?;
        if to.exists() {
            let have = std::fs::read(&to).map_err(|e| format!("{}: {e}", to.display()))?;
            if have != bytes {
                return Err(format!("{} differs from {}", to.display(), from.display()));
            }
        } else {
            if let Some(parent) = to.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| format!("{}: {e}", parent.display()))?;
            }
            std::fs::write(&to, &bytes).map_err(|e| format!("{}: {e}", to.display()))?;
        }
        Ok(to)
    }

    fn row(&self, a: u8, n: u32, i: u32) -> Result<Row, String> {
        Ok(Row {
            id: format!("k{a}n{n}/T{i:02}"),
            a,
            n,
            r: 0.0,
            recipe_seed: None,
            params: Some(self.params(a, n, i)?),
            rho_seed: Some(RHO_SEED_BASE + i128::from(i)),
            suite_id: None,
        })
    }

    /// The pin: `T01` at every size, untimed under `taskset -c 2`, decides
    /// as the earlier comparison's `T01-R1` did, in each of the outputs a
    /// pin compares, and every name is the curve's slug.
    pub fn pin(&self, tools: &Tools) -> Result<J, String> {
        self.refuse_frozen()?;
        let out = self.runs.join("pin").join("pin.json");
        if out.exists() {
            return read(&out);
        }
        let (_, runs) = self
            .constants
            .earlier
            .ok_or("no earlier comparison to pin against")?;
        let mut rows = Vec::new();
        for (a, n) in SIZES {
            let slug = suite::curve_slug(a, n)?;
            let path = self.runs.join(format!("pin/k{a}n{n}-T01.price.json"));
            let new = pin::untimed(&tools.ic, &self.row(a, n, 1)?, &path)?;
            let old = read(
                &self
                    .root
                    .join(runs)
                    .join(format!("k{a}n{n}/T01-R1.price.json")),
            )?;
            let (x, y) = (pin::outputs(&new), pin::outputs(&old));
            let mut differs: Vec<&str> = x
                .iter()
                .zip(&y)
                .filter(|((_, p), (_, q))| !json::py_eq(p, q))
                .map(|((k, _), _)| *k)
                .collect();
            differs.sort_unstable();
            let then_name = curve_id::resolve(old.get("curve").and_then(J::as_str).unwrap_or(""));
            rows.push(json::obj([
                ("row", J::Str(format!("{slug}/T01"))),
                ("slug", J::Str(slug.clone())),
                ("equal", J::Bool(differs.is_empty())),
                (
                    "differs_in",
                    J::Arr(differs.iter().map(|k| J::Str((*k).into())).collect()),
                ),
                (
                    "names_agree",
                    J::Bool(
                        then_name == Some(slug.as_str())
                            && new.get("curve").and_then(J::as_str) == Some(slug.as_str()),
                    ),
                ),
            ]));
        }
        let all = |key: &str| rows.iter().all(|r| r.get(key) == Some(&J::Bool(true)));
        let (held, names_agree) = (all("equal"), all("names_agree"));
        let doc = json::obj([
            ("rows", J::Arr(rows)),
            ("held", J::Bool(held)),
            ("names_agree", J::Bool(names_agree)),
        ]);
        write(&out, json::dumps(&doc, 1) + "\n")?;
        Ok(doc)
    }

    /// One size: `T01`–`T64` once, then `T01`–`T04` again (the A/A), each
    /// the first clean run of up to three attempts through the isolation
    /// tool.  A row that does not complete stops the size; nothing is
    /// overwritten, so a rerun resumes where the last one stopped.
    pub fn size_runs(&self, a: u8, n: u32, tools: &Tools) -> Result<(), String> {
        self.refuse_frozen()?;
        self.manifest(tools)?;
        let d = self.runs.join(format!("k{a}n{n}"));
        std::fs::create_dir_all(&d).map_err(|e| format!("{}: {e}", d.display()))?;
        let stop = d.join("STOPPED");
        if stop.exists() {
            let why = std::fs::read_to_string(&stop).unwrap_or_default();
            println!(
                "{}: stopped earlier ({}); not resumed",
                suite::curve_slug(a, n)?,
                why.trim()
            );
            return Ok(());
        }
        let before = d.join("uptime-before.txt");
        if !before.exists() {
            write(&before, uptime() + "\n")?;
        }
        let bench = bench::Bench {
            isolate: tools.isolate.clone(),
        };
        let order = (1..=TARGETS)
            .map(|i| (i, 1))
            .chain((1..=AA_TARGETS).map(|i| (i, 2)));
        for (i, run) in order {
            let row = self.row(a, n, i)?;
            let out = d.join(format!("T{i:02}-R{run}.price.json"));
            let rep = bench.price(&tools.ic, &row, &out, &self.runs, &[])?;
            let status = rep.get("status").cloned().unwrap_or(J::Null);
            if status.as_str() != Some("complete") {
                write(
                    &stop,
                    format!("T{i:02}-R{run}: status {}\n", py_str(&status)),
                )?;
                println!(
                    "{}: T{i:02}-R{run} is {}; size stopped",
                    suite::curve_slug(a, n)?,
                    py_str(&status)
                );
                return Ok(());
            }
            let speedup = rep
                .get("median")
                .and_then(|m| m.get("online_speedup"))
                .and_then(J::as_f64)
                .unwrap_or(f64::NAN);
            println!(
                "{}: T{i:02}-R{run} online speedup {speedup:.2}",
                suite::curve_slug(a, n)?
            );
        }
        let after = d.join("uptime-after.txt");
        if !after.exists() {
            write(&after, uptime() + "\n")?;
        }
        Ok(())
    }

    /// The manifest, the pin, then every size in the order of `r`.
    pub fn all(&self, tools: &Tools) -> Result<(), String> {
        self.manifest(tools)?;
        let pinned = self.pin(tools)?;
        if pinned.get("held") != Some(&J::Bool(true))
            || pinned.get("names_agree") != Some(&J::Bool(true))
        {
            return Err("the pin failed: the binary does not decide as the earlier comparison did, or a name disagrees; stopped".into());
        }
        for (a, n) in SIZES {
            println!("{}: start", suite::curve_slug(a, n)?);
            self.size_runs(a, n, tools)?;
            println!("{}: done", suite::curve_slug(a, n)?);
        }
        Ok(())
    }
}

// ── the reference check (docs/ic/BOUNDARY_TARGETS.md, 2026-10-01) ────

/// The targets the strong walk prices at each size.
const REFERENCE_TARGETS: u32 = 8;

/// The strong fixture's output for `T<i>`'s attempt `k`.
fn strong_attempt(d: &Path, i: u32, k: usize) -> PathBuf {
    if k == 0 {
        d.join(format!("T{i:02}.strong.json"))
    } else {
        d.join(format!("T{i:02}-retry{k}.strong.json"))
    }
}

/// The fixture's one JSON row, if the attempt wrote one.
fn strong_row(path: &Path) -> Option<J> {
    let text = std::fs::read_to_string(path).ok()?;
    json::parse(text.lines().next()?).ok()
}

/// A point's coordinates as integers, whether written as numbers or as
/// decimal strings.
fn coordinates(p: &J) -> Option<Vec<i128>> {
    p.as_arr()?
        .iter()
        .map(|c| match c {
            J::Int(i) => Some(*i),
            J::Str(s) => s.parse().ok(),
            _ => None,
        })
        .collect()
}

impl Comparison {
    /// The reference check's runs: the strong rho fixture
    /// (`koblitz_rho_fixture <n> <a> signed_frobenius 1 strong`, its
    /// defaults) on `T01`–`T08`'s public targets (`hash:<23000 + i>`, the
    /// domain `ic workflow` hashes with) at every size, each the first clean,
    /// verified run of up to three through the isolation tool.
    pub fn reference(
        &self,
        fixture: &Path,
        fixture_commit: Option<&str>,
        isolate: &Path,
    ) -> Result<(), String> {
        self.refuse_frozen()?;
        let d = self.runs.join("reference");
        let manifest = d.join("fixture.json");
        let sha = bench::sha256_file(fixture)?;
        if manifest.exists() {
            if read(&manifest)?.at("sha256")?.as_str() != Some(sha.as_str()) {
                return Err(format!(
                    "{} is not the fixture {} names",
                    fixture.display(),
                    manifest.display()
                ));
            }
        } else {
            let commit =
                fixture_commit.ok_or("a fresh reference check needs the fixture's commit")?;
            let doc = J::Obj(vec![
                (
                    "path_basename".into(),
                    J::Str(
                        fixture
                            .file_name()
                            .map(|f| f.to_string_lossy().into_owned())
                            .unwrap_or_default(),
                    ),
                ),
                ("sha256".into(), J::Str(sha)),
                ("built_from".into(), J::Str(commit.into())),
                (
                    "command".into(),
                    J::Str(
                        "koblitz_rho_fixture <n> <a> signed_frobenius 1 strong <0x230000 + i> \
                         hash:<23000 + i>, its default lanes and distinguished-point bits"
                            .into(),
                    ),
                ),
            ]);
            write(&manifest, json::dumps(&doc, 1) + "\n")?;
        }
        let bench = bench::Bench {
            isolate: isolate.to_path_buf(),
        };
        for (a, n) in SIZES {
            let sd = d.join(format!("k{a}n{n}"));
            std::fs::create_dir_all(&sd).map_err(|e| format!("{}: {e}", sd.display()))?;
            for i in 1..=REFERENCE_TARGETS {
                let cmd: Vec<String> = vec![
                    fixture.to_string_lossy().into_owned(),
                    n.to_string(),
                    a.to_string(),
                    "signed_frobenius".into(),
                    "1".into(),
                    "strong".into(),
                    (RHO_SEED_BASE + i128::from(i)).to_string(),
                    format!("hash:{}", TARGET_SEED_BASE + i128::from(i)),
                ];
                let mut note = String::from("no clean verified run");
                for k in 0..=runs::RETRIES {
                    let attempt = strong_attempt(&sd, i, k);
                    if !attempt.exists() {
                        bench.launch_to(&cmd, &attempt, &self.runs, true)?;
                    }
                    let row = strong_row(&attempt);
                    let verified = row
                        .as_ref()
                        .is_some_and(|r| r.get("verified") == Some(&J::Bool(true)));
                    if runs::clean(&attempt)? && verified {
                        let r = row.expect("verified");
                        let ns = r.get("walk_ms").and_then(J::as_f64).unwrap_or(f64::NAN) * 1e6
                            / r.get("walk_steps").and_then(J::as_f64).unwrap_or(f64::NAN);
                        note = format!("{ns:.1} ns a step");
                        break;
                    }
                }
                println!("reference {}: T{i:02} {note}", suite::curve_slug(a, n)?);
            }
        }
        Ok(())
    }

    /// The comparison's rho against the strong walk, per step, on the same
    /// targets: the median over `T01`–`T08` of each.  The comparison's rho
    /// is admissible at a size when it is no slower per step, every strong
    /// walk verified, and every target is the comparison's own.
    fn reference_check(&self, a: u8, n: u32, r1: &[Checked]) -> Result<Option<J>, String> {
        if !self.constants.reference_check {
            return Ok(None);
        }
        let d = self.runs.join(format!("reference/k{a}n{n}"));
        let (mut rows, mut strong, mut ours) = (Vec::new(), Vec::new(), Vec::new());
        let mut complete = true;
        for i in 1..=REFERENCE_TARGETS {
            let claim = r1.iter().find(|c| c.target == i);
            let found = (0..=runs::RETRIES)
                .map(|k| strong_attempt(&d, i, k))
                .find_map(|p| {
                    let row = strong_row(&p)?;
                    (runs::clean(&p).ok()? && row.get("verified") == Some(&J::Bool(true)))
                        .then_some(row)
                });
            let (Some(claim), Some(row)) = (claim, found) else {
                complete = false;
                rows.push(J::Obj(vec![
                    ("target".into(), J::Int(i.into())),
                    ("measured".into(), J::Bool(false)),
                ]));
                continue;
            };
            let rep = self.report_of(&claim.claim)?;
            let wall = num(
                rep.at("median")?.at("rho_online_wall_ns")?,
                "rho_online_wall_ns",
            )?;
            let steps = num(
                rep.at("repetitions")?
                    .as_arr()
                    .and_then(<[J]>::first)
                    .ok_or("no repetitions")?
                    .at("rho_online")?
                    .at("steps")?,
                "rho steps",
            )?;
            let s = num(row.at("walk_ms")?, "walk_ms")? * 1e6
                / num(row.at("walk_steps")?, "walk_steps")?;
            let same =
                coordinates(row.at("published_q")?) == coordinates(rep.at("target")?.at("point")?);
            complete &= same;
            strong.push(s);
            ours.push(wall / steps);
            rows.push(J::Obj(vec![
                ("target".into(), J::Int(i.into())),
                ("measured".into(), J::Bool(true)),
                ("same_target".into(), J::Bool(same)),
                ("strong_ns_per_step".into(), J::Float(s)),
                ("rho_ns_per_step".into(), J::Float(wall / steps)),
            ]));
        }
        let (ms, mo) = if strong.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            (stats::median(&strong), stats::median(&ours))
        };
        let ratio = mo / ms;
        Ok(Some(J::Obj(vec![
            (
                "what".into(),
                J::Str(
                    "the comparison's rho against the strong walk (koblitz_rho_fixture strong), ns \
                     per step on the same public targets T01-T08; admissible when no slower"
                        .into(),
                ),
            ),
            ("targets".into(), J::Arr(rows)),
            ("strong_ns_per_step_median".into(), opt(ms)),
            ("rho_ns_per_step_median".into(), opt(mo)),
            ("rho_over_strong".into(), opt(ratio)),
            (
                "admissible".into(),
                J::Bool(complete && ratio.is_finite() && ratio <= 1.0),
            ),
        ])))
    }
}

fn opt(x: f64) -> J {
    if x.is_finite() {
        J::Float(x)
    } else {
        J::Null
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_are_matched_as_glob_and_the_pattern_match_them() {
        assert!(glob2("T01-R1.price.json", "T", "-R", ".price.json"));
        assert!(glob2("T01-R1-retry1.price.json", "T", "-R", ".price.json"));
        assert!(!glob2("T01.params.json", "T", "-R", ".price.json"));
        assert!(glob2(
            "T01-R1-retry2.price.json",
            "T01-R1-retry",
            "",
            ".price.json"
        ));
        assert_eq!(target_run("T07-R2.price.json", ".price.json"), Some((7, 2)));
        assert_eq!(target_run("T07-R2-retry1.price.json", ".price.json"), None);
    }

    #[test]
    fn medians_keep_pythons_types() {
        let ints = |v: &[i128]| v.iter().map(|&i| J::Int(i)).collect::<Vec<_>>();
        assert_eq!(py_median(ints(&[3, 1, 2])).unwrap(), J::Int(2));
        assert_eq!(py_median(ints(&[4, 1, 2, 3])).unwrap(), J::Float(2.5));
        assert_eq!(round2(39.000_001_514_334_05), 39.0);
        assert_eq!(round2(45.678), 45.68);
    }

    #[test]
    fn a_fit_through_a_line_has_its_slope() {
        let pts: Vec<(f64, f64)> = (0..5).map(|i| (i as f64, 0.5 * i as f64 + 1.0)).collect();
        let f = fit(&pts);
        assert_eq!(f.at("slope").unwrap(), &J::Float(0.5));
        assert_eq!(fit(&pts[..1]), J::Null);
    }
}
