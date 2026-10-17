//! Two deliberately separate checks: reproduce a report's bytes, or compare
//! fresh native arithmetic audits and ledgers. Neither establishes a timing win
//! or authenticates a third party's claimed machine or executable provenance.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::canonical::sha256_hex;
use crate::evidence::{self, Record};
use crate::leaderboard::{self, Leaderboard, Manifest};

const LOCK_SCHEMA: &str = "ecbench.backend_reproduction_lock/v1";
const TOOLCHAIN: &str = "1.99.1";
const SESSION_FILES: [&str; 5] = [
    "session.json",
    "host.json",
    "spec.json",
    "plan.json",
    "records.jsonl",
];

/// All source and build inputs to this standalone producer. Compare the bytes
/// compiled into the executable too, so an old executable cannot seal new code.
fn producer_inputs() -> BTreeMap<String, String> {
    macro_rules! inputs {
        ($($path:literal),+ $(,)?) => {
            [$(($path.into(), sha256_hex(include_bytes!(concat!("../../../", $path))))),+]
                .into_iter().collect()
        };
    }
    inputs![
        "tools/backend-leaderboard/Cargo.toml",
        "tools/backend-leaderboard/Cargo.lock",
        "tools/backend-leaderboard/rust-toolchain.toml",
        "tools/backend-leaderboard/build.rs",
        "tools/backend-leaderboard/src/main.rs",
        "tools/backend-leaderboard/src/evidence.rs",
        "tools/backend-leaderboard/src/leaderboard.rs",
        "tools/backend-leaderboard/src/leaderboard.html",
        "tools/backend-leaderboard/src/reproduction.rs",
        "src/cryptanalysis/ecbench/band.rs",
        "src/cryptanalysis/ecbench/canonical.rs",
        "src/cryptanalysis/ecbench/stats.rs",
        "src/cryptanalysis/ecbench/methods.rs",
        "src/cryptanalysis/curve_aliases.json",
        "docs/ic/calibration.json",
        "docs/ecbench/BACKEND_REPRODUCTION.md",
    ]
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct Summary {
    sessions: usize,
    entries: usize,
    planned_without_measurements: usize,
    count_eligible: usize,
    ic_frontiers: usize,
    cross_method_frontiers: usize,
    failed_samples: usize,
    timing_winners: usize,
}

fn summary(b: &Leaderboard) -> Summary {
    Summary {
        sessions: b.sources.len(),
        entries: b.entries.len(),
        planned_without_measurements: b.unattempted.len(),
        count_eligible: b
            .entries
            .iter()
            .filter(|e| e.count_blockers.is_empty())
            .count(),
        ic_frontiers: b.entries.iter().filter(|e| e.ic_count_frontier).count(),
        cross_method_frontiers: b.entries.iter().filter(|e| e.count_frontier).count(),
        failed_samples: b
            .entries
            .iter()
            .flat_map(|e| &e.samples)
            .filter(|s| s.outcome.status != "verified")
            .count(),
        timing_winners: b
            .wall_comparisons
            .iter()
            .filter(|c| c.winner.is_some())
            .count(),
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Lock {
    schema: String,
    toolchain: String,
    manifest: String,
    /// Includes the manifest, every session file, historical receipt and producer input.
    files: BTreeMap<String, String>,
    outputs: BTreeMap<String, String>,
    expected: Summary,
}

fn file_hash(path: &Path) -> Result<String, String> {
    std::fs::read(path)
        .map(|b| sha256_hex(&b))
        .map_err(|e| format!("{}: {e}", path.display()))
}

fn input_hash(root: &Path, relative: &str) -> Result<String, String> {
    leaderboard::relative(Path::new(relative))?;
    let base = root.canonicalize().map_err(|e| e.to_string())?;
    let path = root
        .join(relative)
        .canonicalize()
        .map_err(|e| format!("{relative}: {e}"))?;
    if !path.starts_with(base) {
        return Err(format!("input escapes repository: {relative}"));
    }
    file_hash(&path)
}

fn inventory(root: &Path, manifest: &Path) -> Result<BTreeSet<String>, String> {
    let manifest_path = leaderboard::relative(manifest)?;
    let m: Manifest = evidence::read_json(&root.join(manifest))?;
    if m.schema != leaderboard::MANIFEST_SCHEMA || m.sessions.is_empty() {
        return Err("invalid manifest for reproduction lock".into());
    }
    let mut files: BTreeSet<String> = producer_inputs().into_keys().collect();
    files.insert(manifest_path);
    for input in m.sessions {
        let dir = leaderboard::relative(&input.dir)?;
        for name in SESSION_FILES {
            files.insert(format!("{dir}/{name}"));
        }
        if let Some(path) = input.replay_receipt {
            files.insert(leaderboard::relative(&path)?);
        }
    }
    Ok(files)
}

fn compiler_check() -> Result<(), String> {
    let actual = env!("BACKEND_RUSTC");
    if actual.split_whitespace().nth(1) != Some(TOOLCHAIN) {
        return Err(format!(
            "report compiler must be rustc {TOOLCHAIN}, got {actual}"
        ));
    }
    Ok(())
}

pub(crate) fn compiled_check(root: &Path) -> Result<(), String> {
    compiler_check()?;
    for (path, compiled) in producer_inputs() {
        if input_hash(root, &path)? != compiled {
            return Err(format!(
                "executable was built from different bytes: {path}; rebuild it"
            ));
        }
    }
    Ok(())
}

pub fn verify_inputs(
    root: &Path,
    manifest: &Path,
    path: &Path,
    trusted_hash: Option<&str>,
) -> Result<Lock, String> {
    let bytes = std::fs::read(path).map_err(|e| {
        format!(
            "{}: {e}; a reviewed reproduction lock is required",
            path.display()
        )
    })?;
    if trusted_hash.is_some_and(|h| h != sha256_hex(&bytes)) {
        return Err("reproduction lock does not match the trusted SHA-256".into());
    }
    let lock: Lock =
        serde_json::from_slice(&bytes).map_err(|e| format!("reproduction lock: {e}"))?;
    if lock.schema != LOCK_SCHEMA
        || lock.toolchain != TOOLCHAIN
        || lock.manifest != leaderboard::relative(manifest)?
    {
        return Err("reproduction lock schema, compiler or manifest differs".into());
    }
    if inventory(root, manifest)? != lock.files.keys().cloned().collect() {
        return Err(
            "reproduction lock does not cover exactly the manifest and producer inputs".into(),
        );
    }
    for (file, expected) in &lock.files {
        if input_hash(root, file)? != *expected {
            return Err(format!("locked input SHA-256 mismatch: {file}"));
        }
    }
    compiled_check(root)?;
    Ok(lock)
}

fn output_hashes(board: &Leaderboard) -> Result<BTreeMap<String, String>, String> {
    Ok(leaderboard::render(board)?
        .into_iter()
        .map(|(p, s)| (p, sha256_hex(s.as_bytes())))
        .collect())
}

pub fn verify_outputs(lock: &Lock, board: &Leaderboard) -> Result<(), String> {
    if summary(board) != lock.expected || output_hashes(board)? != lock.outputs {
        return Err(
            "recomputed report differs from the pinned output hashes or expected results".into(),
        );
    }
    Ok(())
}

pub fn seal(root: &Path, manifest: &Path, path: &Path, board: &Leaderboard) -> Result<(), String> {
    compiled_check(root)?;
    let files = inventory(root, manifest)?
        .into_iter()
        .map(|p| Ok((p.clone(), input_hash(root, &p)?)))
        .collect::<Result<_, String>>()?;
    let lock = Lock {
        schema: LOCK_SCHEMA.into(),
        toolchain: TOOLCHAIN.into(),
        manifest: leaderboard::relative(manifest)?,
        files,
        outputs: output_hashes(board)?,
        expected: summary(board),
    };
    let bytes = serde_json::to_vec_pretty(&lock).map_err(|e| e.to_string())?;
    let mut bytes = bytes;
    bytes.push(b'\n');
    std::fs::write(path, &bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    eprintln!(
        "review and pin reproduction lock SHA-256 {}",
        sha256_hex(&bytes)
    );
    Ok(())
}

pub fn report_receipt(path: &Path, lock: &Lock, board: &Leaderboard) -> Result<Value, String> {
    Ok(json!({
        "schema":"ecbench.backend_report_reproduction/v1", "ok":true,
        "lock_sha256":file_hash(path)?, "compiler":env!("BACKEND_RUSTC"),
        "files_checked":lock.files.len(), "outputs":lock.outputs, "expected":summary(board),
        "solver_executed":false, "fresh_arithmetic_replay":false,
        "independent_host":"not_established", "timing_reproduction":"not_tested"
    }))
}

/// Command cards use argv arrays. Placeholders are literal parameters to fill,
/// not executable shell. The entire source session is the acceptance unit.
pub fn plan(board: &Leaderboard) -> Result<String, String> {
    let sessions: Vec<Value> = board.sources.iter().enumerate().map(|(i, s)| {
        let samples: Vec<_> = board.entries.iter().filter(|e| e.source == i).flat_map(|e| &e.samples).collect();
        let mut expected: Vec<_> = samples.iter().filter(|r| r.outcome.status == "verified" && r.cost.deterministic).map(|r| json!({"seq":r.seq,"record_id":r.record_id})).collect();
        expected.sort_by_key(|r| r["seq"].as_u64());
        let mut blockers = Vec::new();
        if s.status != "complete" { blockers.push("original session is incomplete".to_string()); }
        if !s.integrity_ok { blockers.push("original integrity is not complete/bound".into()); }
        if s.git_commit.is_none() || s.git_dirty != Some(false)
            || !matches!(s.git_commit_source.as_deref(), Some("CRYPTO_BUILD_GIT_COMMIT" | "GITHUB_SHA")) {
            blockers.push("clean build-time source is not bound; a runtime checkout is insufficient: obtain the exact original binary or a bound complete source snapshot".into());
        }
        if samples.iter().any(|r| r.outcome.status != "verified" || !r.cost.deterministic) { blockers.push("source contains failed or nondeterministic attempts; strict whole-session ledger reproduction cannot pass".into()); }
        json!({"source":s.dir,"session_id":s.session_id,"git_commit":s.git_commit,"git_dirty":s.git_dirty,"git_commit_source":s.git_commit_source,"original_binary_sha256":s.binary_sha256,
            "original_environment":s.env_class_id,"file_hashes":s.file_hashes,"expected_verified_deterministic_replays":expected,
            "required_measured_records":samples.len(),"known_blockers":blockers,
            "execution_order":["fresh_run_argv","original_replay_argv","fresh_audit_argv","check_argv"],
            "original_replay_argv":["<PINNED_ECBENCH>","verify","--dir",s.dir,"--replay-all","--exit-code","--out","<NEW_EVIDENCE>/original.audit.json"],
            "fresh_run_argv":["<PINNED_ECBENCH>","run","--spec",format!("{}/spec.json",s.dir),"--out","<NEW_EVIDENCE>/session"],
            "fresh_audit_argv":["<PINNED_ECBENCH>","verify","--dir","<NEW_EVIDENCE>/session","--replay-all","--exit-code","--out","<NEW_EVIDENCE>/rerun.audit.json"],
            "check_argv":["<BACKEND_LEADERBOARD>","check-rerun","--source",s.dir,"--rerun","<NEW_EVIDENCE>/session","--original-audit","<NEW_EVIDENCE>/original.audit.json","--rerun-audit","<NEW_EVIDENCE>/rerun.audit.json","--binary","<PINNED_ECBENCH>","--expected-lock-sha256","<TRUSTED_LOCK_SHA256>"]})
    }).collect();
    serde_json::to_string_pretty(&json!({"schema":"ecbench.backend_replay_plan/v1",
        "acceptance":"docs/ecbench/BACKEND_REPRODUCTION.md", "sessions":sessions,
        "gates": {
            "R0":{"checker":"reproduce","toolchain":TOOLCHAIN,"agreement":"SHA-256 of every pinned input and generated output; exact summary and frontier memberships","solver_execution":false},
            "R1":{"checker":"check-rerun","acceptance_unit":"entire session","input_agreement":"exact spec.json and plan.json bytes","audits":"original and fresh session: every deterministic verified measurement, no duplicates, all file bindings and checks pass","ledger_agreement":"exact; includes both solver observations and native counters; only specified direct clock fields excluded","failed_or_nondeterministic_measured_attempts_allowed":0,"independent_host_automatically_established":false},
            "R2":{"checker":"external preregistered statistical receipt; not implemented by this report tool","minimum_matched_rounds_per_target":20,"minimum_aa_pairs_per_target":20,"minimum_isolation":"L2","same_environment_class_required":true,"window":"declare cold, process or online before running; no substitution","paired_ratio":"rerun/original","equivalence_interval":[0.95,1.05],"confidence":0.95,"resamples":10000,"bootstrap_seed":20261017,"statistic":"median paired ratio","pass":"entire paired CI and A/A CIs contained in equivalence interval","performance_promotion":"also require full-cost gain beyond parity and A/A noise, fresh holdouts and matched external references"}
        },
        "notice":"Commands are a protocol, not evidence they ran. Exact result/counter agreement is separate from independent-host evidence and statistical wall-time agreement. No session is freshly reproduced by generating this plan."}))
        .map(|s| s + "\n").map_err(|e| e.to_string())
}

fn raw_records(dir: &Path) -> Result<Vec<Value>, String> {
    std::fs::read_to_string(dir.join("records.jsonl"))
        .map_err(|e| e.to_string())?
        .lines()
        .map(|l| serde_json::from_str(l).map_err(|e| e.to_string()))
        .collect()
}

fn expected_replays(records: &[Record]) -> BTreeMap<u64, String> {
    records
        .iter()
        .filter(|r| r.counts() && r.cost.deterministic)
        .map(|r| (r.seq, r.record_id.clone()))
        .collect()
}

fn check_audit(
    dir: &Path,
    path: &Path,
    records: &[Record],
    binary_sha: &str,
) -> Result<Value, String> {
    let a: Value = evidence::read_json(path)?;
    let session = evidence::read_session(dir)?;
    let verified = records.iter().filter(|r| r.counts()).count();
    if a["schema"] != "ecbench.audit/v1"
        || a["session_id"] != session.session_id
        || a["session_status"] != "complete"
        || a["ok"] != true
        || a["records_bound"] == false
        || (!a["records_bound"].is_null() && a["records_bound"] != true)
        || a["problems"].as_array().is_none_or(|p| !p.is_empty())
        || a["records"].as_u64() != Some(records.len() as u64)
        || a["verified_records"].as_u64() != Some(verified as u64)
        || a["auditor_binary_sha256"] != binary_sha
    {
        return Err(format!(
            "{}: incomplete, failed or mismatched native audit/binary",
            path.display()
        ));
    }
    let checks = a["checks"].as_object().ok_or("missing audit checks")?;
    if checks.is_empty()
        || checks.values().any(|v| v != true)
        || ["replays", "records_hash", "records", "session_complete"]
            .iter()
            .any(|k| checks.get(*k) != Some(&Value::Bool(true)))
    {
        return Err("native audit checks did not all pass including replays".into());
    }
    for name in SESSION_FILES {
        if a["files"][name] != file_hash(&dir.join(name))? {
            return Err(format!("native audit does not bind {name}"));
        }
    }
    let expected = expected_replays(records);
    if expected.is_empty() {
        return Err("no deterministic verified measurements to reproduce".into());
    }
    let replays = a["replays"].as_array().ok_or("missing replay list")?;
    let mut actual = BTreeMap::new();
    for r in replays {
        let seq = r["seq"].as_u64().ok_or("replay lacks sequence")?;
        let id = r["record_id"].as_str().ok_or("replay lacks record id")?;
        if r["reproduced"] != true || actual.insert(seq, id.to_string()).is_some() {
            return Err("failed or duplicate replay".into());
        }
    }
    if actual != expected {
        return Err(
            "native audit must replay exactly every deterministic verified measured record".into(),
        );
    }
    Ok(a)
}

/// Include every solver observation. Native ecbench audit's legacy equality
/// contract does not cover this block, so fresh record streams are necessary.
/// Only clocks are stripped. New/missing observations are differences, not zero.
fn deterministic_ledger(record: &Value) -> Value {
    fn without_clocks(v: &Value) -> Value {
        match v {
            Value::Object(o) => Value::Object(
                o.iter()
                    .filter(|(k, _)| {
                        !matches!(
                            k.as_str(),
                            "wall_ns"
                                | "wall_seconds"
                                | "wall_ms"
                                | "cpu_ns"
                                | "solve_wall_ns"
                                | "elapsed_ns"
                                | "elapsed_ms"
                                | "elapsed_seconds"
                        )
                    })
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect(),
            ),
            _ => v.clone(),
        }
    }
    let keys = [
        "schema",
        "seq",
        "round",
        "warmup",
        "arm",
        "role",
        "method",
        "workload",
        "algorithm_seed",
        "cost",
        "counters",
        "factor_base",
        "fbm_id",
        "field_ops",
        "boundaries",
    ];
    let mut out = serde_json::Map::new();
    for k in keys {
        out.insert(k.into(), record.get(k).cloned().unwrap_or(Value::Null));
    }
    out.insert("outcome".into(), json!({"status":record["outcome"]["status"],"recovered":record["outcome"]["recovered"],"matches_target":record["outcome"]["matches_target"],"matches_planted":record["outcome"]["matches_planted"]}));
    out.insert(
        "phases".into(),
        record["phases"]
            .as_array()
            .map(|p| Value::Array(p.iter().map(without_clocks).collect()))
            .unwrap_or(Value::Null),
    );
    out.insert("solver".into(), without_clocks(&record["solver"]));
    out.insert(
        "decomposition_solver".into(),
        without_clocks(&record["detail"]["decomposition"]["solver"]),
    );
    out.insert(
        "decomposition_system".into(),
        record["detail"]["decomposition"]["system"].clone(),
    );
    Value::Object(out)
}

fn clean_build_identity(host: &Value, commit: &str, binary: &str) -> bool {
    host["build"]["git_commit"] == commit
        && host["build"]["git_dirty"] == false
        && host["build"]["binary_sha256"] == binary
        && matches!(
            host["build"]["git_commit_source"].as_str(),
            Some("CRYPTO_BUILD_GIT_COMMIT" | "GITHUB_SHA")
        )
}

pub fn check_rerun(
    root: &Path,
    board: &Leaderboard,
    source: &Path,
    rerun: &Path,
    original_audit: &Path,
    rerun_audit: &Path,
    binary: &Path,
) -> Result<Value, String> {
    let source_name = leaderboard::relative(source)?;
    let s = board
        .sources
        .iter()
        .find(|s| s.dir == source_name)
        .ok_or("source is not in the locked manifest")?;
    let original = root.join(source);
    if original.canonicalize().map_err(|e| e.to_string())?
        == rerun.canonicalize().map_err(|e| e.to_string())?
    {
        return Err("rerun must be a new session, not the original directory".into());
    }
    for dir in [&original, &rerun.to_path_buf()] {
        let a = evidence::audit_with(dir, 0, false)?;
        if !a.ok
            || a.records_bound != Some(true)
            || evidence::read_session(dir)?.status != "complete"
        {
            return Err(format!(
                "{}: complete bound session integrity is required",
                dir.display()
            ));
        }
    }
    let rs = evidence::read_records(&original)?;
    let fresh = evidence::read_records(rerun)?;
    if rs
        .iter()
        .chain(&fresh)
        .any(|r| !r.warmup && (!r.counts() || !r.cost.deterministic))
    {
        return Err("strict reproduction requires all measured attempts verified and deterministic; failures remain evidence".into());
    }
    let old_session = evidence::read_session(&original)?;
    let new_session = evidence::read_session(rerun)?;
    if old_session.session_id == new_session.session_id {
        return Err("rerun must have a fresh session identity".into());
    }
    for name in ["spec.json", "plan.json"] {
        if file_hash(&original.join(name))? != file_hash(&rerun.join(name))? {
            return Err(format!("rerun changed exact {name}"));
        }
    }
    let binary_sha = file_hash(binary)?;
    if new_session.binary_sha256.as_deref() != Some(&binary_sha) {
        return Err("rerun binary does not match the supplied executable".into());
    }
    let exact_binary = old_session.binary_sha256.as_deref() == Some(&binary_sha);
    let old_host: Value = evidence::read_json(&original.join("host.json"))?;
    let new_host: Value = evidence::read_json(&rerun.join("host.json"))?;
    if old_host["build"]["binary_sha256"].as_str() != old_session.binary_sha256.as_deref()
        || new_host["build"]["binary_sha256"] != binary_sha
    {
        return Err("host build capsules do not bind the measured executables".into());
    }
    if !exact_binary
        && (old_session.git_commit.is_none()
            || old_session.git_dirty != Some(false)
            || old_session.git_commit != new_session.git_commit
            || new_session.git_dirty != Some(false)
            || !clean_build_identity(
                &old_host,
                old_session.git_commit.as_deref().unwrap_or(""),
                old_session.binary_sha256.as_deref().unwrap_or(""),
            )
            || !clean_build_identity(
                &new_host,
                old_session.git_commit.as_deref().unwrap_or(""),
                &binary_sha,
            ))
    {
        return Err(
            "need the exact original binary or matching clean build-time source provenance; runtime checkout identities do not bind code"
                .into(),
        );
    }
    let a = check_audit(&original, original_audit, &rs, &binary_sha)?;
    let b = check_audit(rerun, rerun_audit, &fresh, &binary_sha)?;
    let new_meta: Value = evidence::read_json(&rerun.join("session.json"))?;
    let old_meta: Value = evidence::read_json(&original.join("session.json"))?;
    let start = new_meta["started_unix_ms"]
        .as_u64()
        .ok_or("rerun lacks start timestamp")?;
    let finish = new_meta["finished_unix_ms"]
        .as_u64()
        .ok_or("rerun lacks finish timestamp")?;
    if finish < start
        || old_meta["finished_unix_ms"]
            .as_u64()
            .is_none_or(|old| start <= old)
        || [&a, &b]
            .iter()
            .any(|audit| audit["audited_unix_ms"].as_u64().is_none_or(|t| t < finish))
    {
        return Err(
            "require a later fresh session and both audits generated after it finished".into(),
        );
    }
    let raw_old = raw_records(&original)?;
    let raw_new = raw_records(rerun)?;
    if raw_old.len() != raw_new.len() {
        return Err("rerun omitted or added attempts".into());
    }
    let mut ledger_bytes = Vec::new();
    for (old, new) in raw_old.iter().zip(&raw_new) {
        let expected = deterministic_ledger(old);
        let actual = deterministic_ledger(new);
        if expected != actual {
            let fields: Vec<_> = expected
                .as_object()
                .unwrap()
                .keys()
                .filter(|k| expected[*k] != actual[*k])
                .cloned()
                .collect();
            return Err(format!(
                "rerun seq {} differs in {}",
                old["seq"],
                fields.join(", ")
            ));
        }
        serde_json::to_writer(&mut ledger_bytes, &expected).map_err(|e| e.to_string())?;
        ledger_bytes.push(b'\n');
    }
    Ok(
        json!({"schema":"ecbench.backend_ledger_reproduction/v1", "ok":true,
        "source":s.dir,"source_session_id":s.session_id,"rerun_session_id":new_session.session_id,
        "binary_sha256":binary_sha,"exact_original_binary":exact_binary,"declared_source_commit":new_session.git_commit,
        "original_audit_sha256":file_hash(original_audit)?,"rerun_audit_sha256":file_hash(rerun_audit)?,
        "original_files":a["files"],"rerun_files":b["files"],
        "original_auditor_environment":a["auditor_env_class_id"],"rerun_auditor_environment":b["auditor_env_class_id"],
        "records_compared":raw_old.len(),"deterministic_ledger_sha256":sha256_hex(&ledger_bytes),
        "solver_statistics_compared":true,"fresh_arithmetic_audits_checked":true,
        "independent_host":"requires external runner provenance; environment IDs alone do not establish physical independence",
        "timing_reproduction":"not_tested","eligible_for_performance_promotion":false,
        "trust_boundary":"checks supplied bytes and audit assertions; not a signature, remote attestation, or proof that a third party executed them"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_dir(label: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("backend-{label}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn save(path: &Path, value: &Value) {
        std::fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
    }

    #[test]
    fn reproduction_lock_rejects_changed_inputs_outputs_inventory_and_stale_executables() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let dir = test_dir("lock-adversarial");
        let source = "research/f4_f5_f6_cross_family_20261011/sessions/K13";
        let receipt = "research/f4_f5_f6_cross_family_20261011/audits/K13-replay-all.json";
        let mut paths: Vec<_> = producer_inputs().into_keys().collect();
        paths.extend(SESSION_FILES.map(|name| format!("{source}/{name}")));
        paths.push(receipt.into());
        for path in paths {
            let target = dir.join(&path);
            std::fs::create_dir_all(target.parent().unwrap()).unwrap();
            std::fs::copy(root.join(path), target).unwrap();
        }
        let manifest = Path::new("manifest.json");
        save(
            &dir.join(manifest),
            &json!({"schema":leaderboard::MANIFEST_SCHEMA,"description":"test fixture only", "field_degrees":[13],"sessions":[{"dir":source,"replay_receipt":receipt}]}),
        );
        let board = leaderboard::build(&dir, &dir.join(manifest)).unwrap();
        let lock_path = dir.join("lock.json");
        seal(&dir, manifest, &lock_path, &board).unwrap();
        let trusted = file_hash(&lock_path).unwrap();
        let lock = verify_inputs(&dir, manifest, &lock_path, Some(&trusted)).unwrap();
        verify_outputs(&lock, &board).unwrap();
        assert!(verify_inputs(&dir, manifest, &lock_path, Some(&"0".repeat(64))).is_err());
        let original_lock_bytes = std::fs::read(&lock_path).unwrap();
        let original_lock: Value = serde_json::from_slice(&original_lock_bytes).unwrap();
        let mut changed = original_lock.clone();
        changed["files"]
            .as_object_mut()
            .unwrap()
            .remove(&format!("{source}/host.json"));
        save(&lock_path, &changed);
        assert!(verify_inputs(&dir, manifest, &lock_path, None)
            .unwrap_err()
            .contains("cover exactly"));
        std::fs::write(&lock_path, &original_lock_bytes).unwrap();
        let record_file = dir.join(source).join("records.jsonl");
        let original_records = std::fs::read(&record_file).unwrap();
        std::fs::write(&record_file, b"changed evidence").unwrap();
        assert!(verify_inputs(&dir, manifest, &lock_path, Some(&trusted))
            .unwrap_err()
            .contains("SHA-256 mismatch"));
        std::fs::write(record_file, original_records).unwrap();
        let mut changed = original_lock.clone();
        changed["outputs"]["index.html"] = json!("wrong");
        let wrong: Lock = serde_json::from_value(changed).unwrap();
        assert!(verify_outputs(&wrong, &board).is_err());
        let source_path = "tools/backend-leaderboard/src/main.rs";
        std::fs::write(dir.join(source_path), b"changed producer").unwrap();
        let mut changed = original_lock;
        changed["files"][source_path] = json!(input_hash(&dir, source_path).unwrap());
        save(&lock_path, &changed);
        assert!(verify_inputs(&dir, manifest, &lock_path, None)
            .unwrap_err()
            .contains("executable was built from different bytes"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    // This is a synthetic receipt/record fixture, never a claimed solver run.
    // Reseal it so semantic rejection cannot be attributed to a broken hash.
    fn reseal_fixture(dir: &Path, rows: &mut [Value], audit: &mut Value) {
        let mut text = String::new();
        for r in rows {
            r["record_id"] = json!("");
            let empty = serde_json::to_string(r).unwrap();
            r["record_id"] = json!(format!("ECR1h{}", &sha256_hex(empty.as_bytes())[..12]));
            text.push_str(&serde_json::to_string(r).unwrap());
            text.push('\n');
        }
        std::fs::write(dir.join("records.jsonl"), &text).unwrap();
        let mut session: Value = evidence::read_json(&dir.join("session.json")).unwrap();
        session["records_sha256"] = json!(sha256_hex(text.as_bytes()));
        session["host_sha256"] = json!(file_hash(&dir.join("host.json")).unwrap());
        save(&dir.join("session.json"), &session);
        for name in SESSION_FILES {
            audit["files"][name] = json!(file_hash(&dir.join(name)).unwrap());
        }
        let rs = evidence::read_records(dir).unwrap();
        audit["replays"] = Value::Array(expected_replays(&rs).into_iter().map(|(seq, id)| json!({"seq":seq,"record_id":id,"reproduced":true,"detail":"synthetic test fixture"})).collect());
    }

    #[test]
    fn fresh_session_comparison_accepts_clock_changes_but_rejects_resealed_solver_changes() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let dir = test_dir("fresh-ledger-fixture");
        let fresh = dir.join("session");
        std::fs::create_dir(&fresh).unwrap();
        let source = Path::new("research/f4_f5_f6_cross_family_20261011/sessions/K13");
        for name in SESSION_FILES {
            std::fs::copy(root.join(source).join(name), fresh.join(name)).unwrap();
        }
        let binary = dir.join("fixture-binary");
        std::fs::write(&binary, b"synthetic test fixture; never executed").unwrap();
        let hash = file_hash(&binary).unwrap();
        let original_fixture = dir.join("original");
        std::fs::create_dir(&original_fixture).unwrap();
        for name in SESSION_FILES {
            std::fs::copy(root.join(source).join(name), original_fixture.join(name)).unwrap();
        }
        let mut original_session: Value =
            evidence::read_json(&original_fixture.join("session.json")).unwrap();
        original_session["binary_sha256"] = json!(hash);
        save(&original_fixture.join("session.json"), &original_session);
        for capsule_dir in [&original_fixture, &fresh] {
            let mut host: Value = evidence::read_json(&capsule_dir.join("host.json")).unwrap();
            host["build"]["binary_sha256"] = json!(hash);
            save(&capsule_dir.join("host.json"), &host);
        }
        let mut session: Value = evidence::read_json(&fresh.join("session.json")).unwrap();
        session["session_id"] = json!("synthetic-new-session");
        session["binary_sha256"] = json!(hash);
        let timestamp = session["finished_unix_ms"].as_u64().unwrap();
        session["started_unix_ms"] = json!(timestamp + 1);
        session["finished_unix_ms"] = json!(timestamp + 2);
        save(&fresh.join("session.json"), &session);
        let mut rows = raw_records(&fresh).unwrap();
        for r in &mut rows {
            r["session_id"] = session["session_id"].clone();
            r["binary_sha256"] = json!(hash);
            r["time"]["solve_wall_ns"] = json!(123);
            if let Some(phases) = r["phases"].as_array_mut() {
                for p in phases {
                    p["wall_ns"] = json!(777);
                }
            }
        }
        let original_path =
            root.join("research/f4_f5_f6_cross_family_20261011/audits/K13-replay-all.json");
        let mut original_audit: Value = evidence::read_json(&original_path).unwrap();
        original_audit["auditor_binary_sha256"] = json!(hash);
        original_audit["audited_unix_ms"] = json!(timestamp + 3);
        let mut original_rows = raw_records(&original_fixture).unwrap();
        for r in &mut original_rows {
            r["binary_sha256"] = json!(hash);
        }
        reseal_fixture(&original_fixture, &mut original_rows, &mut original_audit);
        save(&dir.join("original.json"), &original_audit);
        let mut fresh_audit = original_audit;
        fresh_audit["session_id"] = session["session_id"].clone();
        reseal_fixture(&fresh, &mut rows, &mut fresh_audit);
        save(&dir.join("fresh.json"), &fresh_audit);
        save(
            &dir.join("manifest.json"),
            &json!({"schema":leaderboard::MANIFEST_SCHEMA,"description":"synthetic fixture, never solver evidence","field_degrees":[13],"sessions":[{"dir":"original","replay_receipt":null}]}),
        );
        let board = leaderboard::build(&dir, &dir.join("manifest.json")).unwrap();
        let check = || {
            check_rerun(
                &dir,
                &board,
                Path::new("original"),
                &fresh,
                &dir.join("original.json"),
                &dir.join("fresh.json"),
                &binary,
            )
        };
        let receipt = check().unwrap();
        assert_eq!(receipt["records_compared"], 48);
        assert_eq!(receipt["eligible_for_performance_promotion"], false);
        rows[0]["detail"]["decomposition"] =
            json!({"solver":{"ops":1,"op_unit":"conflicts","conflicts":1}});
        reseal_fixture(&fresh, &mut rows, &mut fresh_audit);
        save(&dir.join("fresh.json"), &fresh_audit);
        assert!(check().unwrap_err().contains("decomposition_solver"));
        // A record can be well sealed and its receipt updated, yet disagree.
        assert!(evidence::audit_with(&fresh, 0, false).unwrap().ok);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn runtime_checkout_cannot_substitute_for_clean_build_provenance() {
        let mut host = json!({"build":{"git_commit":"revision","git_dirty":false,"binary_sha256":"binary","git_commit_source":"runtime_worktree"}});
        assert!(!clean_build_identity(&host, "revision", "binary"));
        host["build"]["git_commit_source"] = json!("CRYPTO_BUILD_GIT_COMMIT");
        assert!(clean_build_identity(&host, "revision", "binary"));
        host["build"]["git_dirty"] = json!(true);
        assert!(!clean_build_identity(&host, "revision", "binary"));
    }

    #[test]
    fn solver_counters_are_exact_but_solver_clocks_are_not_reproducibility_inputs() {
        let mut a = json!({"schema":"ecbench.record/v1","seq":1,"cost":{"deterministic":true},"solver":{"ops":90,"conflicts":17,"learnt_clauses":6,"wall_seconds":1.2},"phases":[{"name":"relations","wall_ns":100,"native":{"xors":91}}],"detail":{"decomposition":{"solver":{"conflicts":17,"wall_seconds":1.2}}}});
        let mut b = a.clone();
        b["solver"]["wall_seconds"] = json!(8.7);
        b["detail"]["decomposition"]["solver"]["wall_seconds"] = json!(8.7);
        b["phases"][0]["wall_ns"] = json!(900);
        assert_eq!(deterministic_ledger(&a), deterministic_ledger(&b));
        b["solver"]["conflicts"] = json!(18);
        assert_ne!(deterministic_ledger(&a), deterministic_ledger(&b));
        b = a.clone();
        b["detail"]["decomposition"]["solver"]["conflicts"] = json!(18);
        assert_ne!(deterministic_ledger(&a), deterministic_ledger(&b));
        b = a.clone();
        b["solver"]["learnt_clauses"] = json!(7);
        assert_ne!(deterministic_ledger(&a), deterministic_ledger(&b));
        a["solver"]["counters"] = json!({"wall_ns":7});
        b = a.clone();
        b["solver"]["counters"]["wall_ns"] = json!(8);
        assert_ne!(
            deterministic_ledger(&a),
            deterministic_ledger(&b),
            "clock exclusions cannot remove a nested counter"
        );
        a["solver"]["ops"] = json!(9_007_199_254_740_993u64);
        b = a.clone();
        b["solver"]["ops"] = json!(9_007_199_254_740_994u64);
        assert_ne!(deterministic_ledger(&a), deterministic_ledger(&b));
    }

    #[test]
    fn native_receipt_must_cover_all_records_without_duplicates_or_wrong_files() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let dir = root.join("research/f4_f5_f6_cross_family_20261011/sessions/K13");
        let original =
            root.join("research/f4_f5_f6_cross_family_20261011/audits/K13-replay-all.json");
        let records = evidence::read_records(&dir).unwrap();
        let audit: Value = evidence::read_json(&original).unwrap();
        let binary = audit["auditor_binary_sha256"].as_str().unwrap().to_string();
        let path = std::env::temp_dir().join(format!(
            "backend-replay-receipt-{}.json",
            std::process::id()
        ));
        let check = |v: &Value| {
            std::fs::write(&path, serde_json::to_vec(v).unwrap()).unwrap();
            check_audit(&dir, &path, &records, &binary)
        };
        assert!(check(&audit).is_ok());
        let mut bad = audit.clone();
        bad["replays"].as_array_mut().unwrap().pop();
        assert!(check(&bad).is_err());
        let mut bad = audit.clone();
        let duplicate = bad["replays"][0].clone();
        bad["replays"].as_array_mut().unwrap().push(duplicate);
        assert!(check(&bad).is_err());
        let mut bad = audit.clone();
        bad["files"]["plan.json"] = json!("wrong");
        assert!(check(&bad).is_err());
        let mut bad = audit.clone();
        bad["replays"][0]["reproduced"] = json!(false);
        assert!(check(&bad).is_err());
        let mut bad = audit.clone();
        bad["records_bound"] = json!(false);
        assert!(check(&bad).is_err());
        let mut bad = audit.clone();
        bad["checks"]
            .as_object_mut()
            .unwrap()
            .remove("records_hash");
        assert!(check(&bad).is_err());
        std::fs::remove_file(path).unwrap();
    }
}
