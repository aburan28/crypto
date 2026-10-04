//! New target-free capsule contract; old target-control capsules are never admitted.
use super::native;
use crypto_lib::cryptanalysis::{
    koblitz_fast::FastCurve,
    koblitz_index_calculus::{probe_scalar, KoblitzCurve},
    prepared_ordinary::{Family, Observer, Plan},
    prepared_sat_control::{canonical_sha, encoded, sha256},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

pub const SCOPE: &str = "disclosed-n17-native-ordinary-preparation-v1";
pub const WORKER: &str = "prepared_ordinary_worker";
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub schema_version: u32,
    pub question: String,
    pub plan: Plan,
    pub export_nonce: u64,
    pub conflict_budget: u64,
    pub exporter_timeout_ms: u64,
    pub solver_timeout_ms: u64,
    pub worker_timeout_ms: u64,
}
impl Config {
    pub fn validate(&self) -> Result<(), String> {
        self.plan.validate()?;
        native::require(
            self.schema_version == 1
                && self.question == SCOPE
                && (1..=1_000_000).contains(&self.conflict_budget)
                && (1..=60_000).contains(&self.exporter_timeout_ms)
                && (1..=120_000).contains(&self.solver_timeout_ms)
                && (1..=43_200_000).contains(&self.worker_timeout_ms),
            "outside bounded ordinary-worker configuration",
        )
    }
}
pub fn config(capsule: &Path) -> Result<Config, String> {
    let cfg: Config =
        serde_json::from_slice(&native::read(&capsule.join("config.json"), 64 * 1024)?)
            .map_err(|e| e.to_string())?;
    cfg.validate()?;
    Ok(cfg)
}
pub fn environment() -> BTreeMap<String, String> {
    [
        ("LC_ALL".into(), "C".into()),
        ("RAYON_NUM_THREADS".into(), "1".into()),
    ]
    .into()
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
pub fn identity() -> Value {
    json!({"schema_version":1,"source_manifest_sha256":option_env!("IC_ORDINARY_SOURCE_MANIFEST_SHA256"),
        "build_sha256":option_env!("IC_ORDINARY_BUILD_SHA256"),
        "target_os":std::env::consts::OS,"target_arch":std::env::consts::ARCH})
}
pub fn require_identity(compiled: &Value, expected: &Value) -> Result<(), String> {
    native::require(
        compiled == expected
            && compiled["schema_version"] == 1
            && compiled["source_manifest_sha256"]
                .as_str()
                .is_some_and(digest)
            && compiled["build_sha256"].as_str().is_some_and(digest)
            && compiled["target_os"] == "macos"
            && compiled["target_arch"] == "aarch64",
        "worker lacks matching compiled source/build identity",
    )
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Registration {
    pub schema_version: u32,
    pub scope: String,
    pub source_commit: String,
    pub hardware: Value,
    pub immutable_files: Value,
    pub source_manifest_sha256: String,
    pub config_sha256: String,
    pub host_context_sha256: String,
    pub worker_sha256: String,
    pub auditor_sha256: String,
    pub worker_build_identity: Value,
    pub runtime_environment: BTreeMap<String, String>,
    pub validation_only: bool,
}
pub fn registration(capsule: &Path) -> Result<(Registration, String), String> {
    let bytes = native::read(&capsule.join("registration.json"), 16 * 1024 * 1024)?;
    let record: Registration = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let hash = canonical_sha(&serde_json::from_slice(&bytes).map_err(|e| e.to_string())?)?;
    native::require(
        native::load(&capsule.join("seal.json"))?["registration_sha256"] == hash,
        "ordinary registration seal differs",
    )?;
    native::require(
        record.schema_version == 1
            && record.scope == SCOPE
            && record.hardware == json!({"os":"macos","architecture":"aarch64"})
            && record.source_commit.len() == 40
            && record
                .source_commit
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            && [
                &record.source_manifest_sha256,
                &record.config_sha256,
                &record.host_context_sha256,
                &record.worker_sha256,
                &record.auditor_sha256,
            ]
            .into_iter()
            .all(|s| digest(s))
            && record.runtime_environment == environment(),
        "ordinary registration contract differs",
    )?;
    require_identity(&record.worker_build_identity, &record.worker_build_identity)?;
    native::require(
        record.worker_build_identity["source_manifest_sha256"] == record.source_manifest_sha256,
        "ordinary compiled source binding differs",
    )?;
    native::require(
        record.config_sha256 == sha256(&native::read(&capsule.join("config.json"), 64 * 1024)?)
            && record.host_context_sha256
                == sha256(&native::read(
                    &capsule.join("host-context.json"),
                    64 * 1024,
                )?),
        "ordinary configuration/host binding differs",
    )?;
    config(capsule)?;
    native::require(
        !capsule.join("preparation.json").exists() && !capsule.join("job.json").exists(),
        "ordinary capsule must not import a target-control job or known logs",
    )?;
    Ok((record, hash))
}
pub fn check_capsule(capsule: &Path) -> Result<(Registration, String), String> {
    let (record, hash) = registration(capsule)?;
    native::check_tree(&capsule.join("immutable"), &record.immutable_files)?;
    native::require(
        canonical_sha(&native::inventory(&capsule.join("immutable/source"))?)?
            == record.source_manifest_sha256
            && sha256(&native::read(
                &capsule.join("immutable/bin").join(WORKER),
                128 * 1024 * 1024,
            )?) == record.worker_sha256
            && sha256(&native::read(
                &capsule.join("immutable/bin/icprog"),
                128 * 1024 * 1024,
            )?) == record.auditor_sha256,
        "ordinary source or executable binding differs",
    )?;
    let cfg = config(capsule)?;
    if cfg.plan.family == Family::Cryptominisat {
        check_roles(capsule)?;
    }
    Ok((record, hash))
}
pub fn check_roles(capsule: &Path) -> Result<(), String> {
    native::require(
        sha256(&native::read(
            &capsule.join("immutable/assets/bin/cms"),
            4 * 1024 * 1024,
        )?) == native::CMS_SHA
            && sha256(&native::read(
                &capsule.join("immutable/assets/bin/exporter"),
                2 * 1024 * 1024,
            )?) == native::EXPORTER_SHA,
        "accepted ordinary native role binary differs",
    )
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Claim {
    pub schema_version: u32,
    pub scope: String,
    pub registration_sha256: String,
    pub execution: PathBuf,
    pub worker_sha256: String,
}
pub fn claim(
    capsule: &Path,
    execution: &Path,
    record: &Registration,
    hash: &str,
) -> Result<(), String> {
    let claimed: Claim =
        serde_json::from_slice(&native::read(&capsule.join("consumed.json"), 64 * 1024)?)
            .map_err(|e| e.to_string())?;
    native::require(
        !record.validation_only
            && claimed.schema_version == 1
            && claimed.scope == SCOPE
            && claimed.registration_sha256 == hash
            && claimed.worker_sha256 == record.worker_sha256
            && execution.is_absolute()
            && claimed.execution == execution,
        "ordinary one-use execution claim differs or registration is validation-only",
    )
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Query {
    pub trial: usize,
    pub scalar: u64,
    pub point: [u64; 2],
}
pub fn query(attempt: &Path, cfg: &Config) -> Result<Query, String> {
    let q: Query = serde_json::from_slice(&native::read(&attempt.join("query.json"), 64 * 1024)?)
        .map_err(|e| e.to_string())?;
    native::require(
        cfg.plan.family == Family::Cryptominisat
            && q.trial < cfg.plan.planned_queries
            && q.scalar == probe_scalar(cfg.plan.algorithm_seed, q.trial as u64, 65587),
        "ordinary native query differs from bounded frozen law",
    )?;
    let curve = KoblitzCurve::new(1, 17).ok_or("ordinary n17 curve unavailable")?;
    let fast = FastCurve::new(&curve.curve).ok_or("ordinary portable curve unavailable")?;
    native::require(
        encoded(fast.mul_u64(fast.lift(curve.generator()), q.scalar)) == Some(q.point),
        "ordinary public point differs from its frozen scalar",
    )?;
    Ok(q)
}
pub fn helper_context(
    capsule: &Path,
    attempt: &Path,
    parent: u32,
) -> Result<(Config, Query), String> {
    let (record, hash) = registration(capsule)?;
    let claimed: Claim =
        serde_json::from_slice(&native::read(&capsule.join("consumed.json"), 64 * 1024)?)
            .map_err(|e| e.to_string())?;
    let execution = claimed
        .execution
        .canonicalize()
        .map_err(|e| e.to_string())?;
    claim(capsule, &execution, &record, &hash)?;
    let started = native::load(&execution.join("worker-started.json"))?;
    native::require(
        parent > 1
            && started["pid"] == parent
            && started["registration_sha256"] == hash
            && started["worker_sha256"] == record.worker_sha256,
        "ordinary native helper is not a child of the claimed worker",
    )?;
    let cfg = config(capsule)?;
    let q = query(attempt, &cfg)?;
    native::require(
        attempt.canonicalize().map_err(|e| e.to_string())?
            == execution
                .join(format!("query-{:03}", q.trial))
                .canonicalize()
                .map_err(|e| e.to_string())?,
        "ordinary native query is outside claimed execution",
    )?;
    let start = native::load(&execution.join(format!("start-{:03}.json", q.trial)))?;
    native::require(
        start["trial"] == q.trial && start["scalar"] == q.scalar,
        "ordinary native query lacks its durable start",
    )?;
    Ok((cfg, q))
}
pub fn role(
    capsule: &Path,
    cfg: &Config,
    q: &Query,
    role: &str,
) -> Result<(PathBuf, Vec<String>, &'static str), String> {
    cfg.validate()?;
    native::require(
        cfg.plan.family == Family::Cryptominisat && q.trial < cfg.plan.planned_queries,
        "ordinary native role is outside external CMS panel",
    )?;
    match role {
        "exporter" => Ok((
            capsule.join("immutable/assets/bin/exporter"),
            vec![
                "17".into(),
                "6".into(),
                "standard".into(),
                cfg.export_nonce.to_string(),
                cfg.conflict_budget.to_string(),
                "instance".into(),
                "1".into(),
                "0".into(),
                "--target-x".into(),
                q.point[0].to_string(),
                "--target-y".into(),
                q.point[1].to_string(),
                "--blind-instance-id".into(),
                format!("native-ordinary-{:03}", q.trial),
                "--export-only".into(),
            ],
            native::EXPORTER_SHA,
        )),
        "cms" => Ok((
            capsule.join("immutable/assets/bin/cms"),
            vec![
                "--verb".into(),
                "1".into(),
                "--threads".into(),
                "1".into(),
                "--random".into(),
                "1".into(),
                "--maxsol".into(),
                "1".into(),
                "--maxconfl".into(),
                cfg.conflict_budget.to_string(),
                "instance/instance.xor.cnf".into(),
            ],
            native::CMS_SHA,
        )),
        _ => Err("unsupported ordinary native role".into()),
    }
}

pub struct DurableObserver {
    execution: PathBuf,
    next: usize,
    pending: Option<(usize, u64)>,
}
impl DurableObserver {
    pub fn new(execution: &Path) -> Self {
        Self {
            execution: execution.into(),
            next: 0,
            pending: None,
        }
    }
}
impl Observer for DurableObserver {
    fn started(&mut self, trial: usize, scalar: u64) -> Result<(), String> {
        native::require(
            self.pending.is_none() && trial == self.next,
            "ordinary start chronology differs",
        )?;
        native::save(
            &self.execution.join(format!("start-{trial:03}.json")),
            &json!({"trial":trial,"scalar":scalar}),
        )?;
        self.pending = Some((trial, scalar));
        Ok(())
    }
    fn completed(&mut self, record: &Value) -> Result<(), String> {
        let (trial, scalar) = self.pending.ok_or("ordinary completion has no start")?;
        native::require(
            record["attempt"]["trial"] == trial && record["attempt"]["scalar"] == scalar,
            "ordinary completion differs from durable start",
        )?;
        native::save(
            &self.execution.join(format!("progress-{trial:03}.json")),
            record,
        )?;
        self.pending = None;
        self.next += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    fn cfg() -> Config {
        Config {
            schema_version: 1,
            question: SCOPE.into(),
            plan: Plan {
                schema_version: 1,
                question: "native-target-free-preparation-n17-v1".into(),
                family: Family::Cryptominisat,
                algorithm_seed: 2026100311,
                planned_queries: 2,
            },
            export_nonce: 2026100311,
            conflict_budget: 100000,
            exporter_timeout_ms: 30000,
            solver_timeout_ms: 30000,
            worker_timeout_ms: 120000,
        }
    }
    fn temp(label: &str) -> PathBuf {
        let p =
            std::env::temp_dir().join(format!("ordinary-contract-{label}-{}", std::process::id()));
        fs::create_dir(&p).unwrap();
        p
    }
    #[test]
    fn strict_configuration_and_compiled_identity_are_required() {
        let mut cfg = cfg();
        cfg.validate().unwrap();
        let mut v = serde_json::to_value(&cfg).unwrap();
        v["target"] = json!([52411, 72106]);
        assert!(serde_json::from_value::<Config>(v).is_err());
        cfg.plan.planned_queries = 513;
        assert!(cfg.validate().is_err());
        cfg = tests::cfg();
        cfg.solver_timeout_ms = 120001;
        assert!(cfg.validate().is_err());
        let good = json!({"schema_version":1,"source_manifest_sha256":"a".repeat(64),"build_sha256":"b".repeat(64),
            "target_os":"macos","target_arch":"aarch64"});
        require_identity(&good, &good).unwrap();
        let mut wrong = good.clone();
        wrong["build_sha256"] = json!(null);
        assert!(require_identity(&wrong, &wrong).is_err());
        wrong = good.clone();
        wrong["source_manifest_sha256"] = json!("c".repeat(64));
        assert!(require_identity(&good, &wrong).is_err());
    }
    #[test]
    fn durable_records_preserve_pending_prefix_and_reject_overwrites() {
        let p = temp("progress");
        let mut observer = DurableObserver::new(&p);
        assert!(observer.started(1, 10).is_err());
        observer.started(0, 10).unwrap();
        assert!(observer.started(1, 20).is_err());
        assert!(observer
            .completed(&json!({"attempt":{"trial":0,"scalar":11}}))
            .is_err());
        assert!(p.join("start-000.json").exists());
        assert!(!p.join("progress-000.json").exists());
        let row = json!({"attempt":{"trial":0,"scalar":10,"outcome":"timeout","indices":null}});
        observer.completed(&row).unwrap();
        assert!(observer.completed(&row).is_err());
        observer.started(1, 20).unwrap();
        let prefix = native::load(&p.join("progress-000.json")).unwrap();
        assert_eq!(prefix, row);
        assert!(DurableObserver::new(&p).started(0, 10).is_err());
        assert!(p.join("start-001.json").exists());
        assert!(!p.join("progress-001.json").exists());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn helper_requires_new_claim_parent_durable_start_and_exact_public_query() {
        let p = temp("helper");
        let execution = p.join("execution");
        fs::create_dir(&execution).unwrap();
        let cfg = cfg();
        native::save(&p.join("config.json"), &json!(cfg)).unwrap();
        native::save(
            &p.join("host-context.json"),
            &json!({"mock_contract_control":true}),
        )
        .unwrap();
        let source = "a".repeat(64);
        let record = Registration {
            schema_version: 1,
            scope: SCOPE.into(),
            source_commit: "c".repeat(40),
            hardware: json!({"os":"macos","architecture":"aarch64"}),
            immutable_files: json!({}),
            source_manifest_sha256: source.clone(),
            config_sha256: sha256(&native::read(&p.join("config.json"), 65536).unwrap()),
            host_context_sha256: sha256(
                &native::read(&p.join("host-context.json"), 65536).unwrap(),
            ),
            worker_sha256: "d".repeat(64),
            auditor_sha256: "e".repeat(64),
            worker_build_identity: json!({"schema_version":1,"source_manifest_sha256":source,
                "build_sha256":"b".repeat(64),"target_os":"macos","target_arch":"aarch64"}),
            runtime_environment: environment(),
            validation_only: false,
        };
        let hash = canonical_sha(&json!(record)).unwrap();
        native::save(&p.join("registration.json"), &json!(record)).unwrap();
        native::save(&p.join("seal.json"), &json!({"registration_sha256":hash})).unwrap();
        let claimed = Claim {
            schema_version: 1,
            scope: SCOPE.into(),
            registration_sha256: hash.clone(),
            execution: execution.clone(),
            worker_sha256: record.worker_sha256.clone(),
        };
        native::save(&p.join("consumed.json"), &json!(claimed)).unwrap();
        assert!(native::save(&p.join("consumed.json"), &json!(claimed)).is_err());
        claim(&p, &execution, &record, &hash).unwrap();
        assert!(claim(&p, &p, &record, &hash).is_err());
        native::save(
            &execution.join("worker-started.json"),
            &json!({"pid":1337,
            "registration_sha256":hash,"worker_sha256":record.worker_sha256}),
        )
        .unwrap();
        let attempt = execution.join("query-000");
        fs::create_dir(&attempt).unwrap();
        let scalar = probe_scalar(cfg.plan.algorithm_seed, 0, 65587);
        let curve = KoblitzCurve::new(1, 17).unwrap();
        let fast = FastCurve::new(&curve.curve).unwrap();
        let q = Query {
            trial: 0,
            scalar,
            point: encoded(fast.mul_u64(fast.lift(curve.generator()), scalar)).unwrap(),
        };
        native::save(&attempt.join("query.json"), &json!(q)).unwrap();
        assert!(helper_context(&p, &attempt, 1337).is_err());
        native::save(
            &execution.join("start-000.json"),
            &json!({"trial":0,"scalar":scalar}),
        )
        .unwrap();
        helper_context(&p, &attempt, 1337).unwrap();
        assert!(helper_context(&p, &attempt, 1338).is_err());
        let mut wrong = json!(q);
        wrong["point"][0] = json!(q.point[0] ^ 1);
        fs::write(
            attempt.join("query.json"),
            serde_json::to_vec(&wrong).unwrap(),
        )
        .unwrap();
        assert!(helper_context(&p, &attempt, 1337).is_err());
        let mut old = json!(record);
        old["scope"] = json!("disclosed-n17-native-prepared-control");
        fs::write(
            p.join("registration.json"),
            serde_json::to_vec(&old).unwrap(),
        )
        .unwrap();
        fs::write(
            p.join("seal.json"),
            serde_json::to_vec(&json!({"registration_sha256":canonical_sha(&old).unwrap()}))
                .unwrap(),
        )
        .unwrap();
        assert!(registration(&p).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn native_roles_use_only_accepted_paths_and_frozen_single_thread_argv() {
        let cfg = cfg();
        let q = Query {
            trial: 1,
            scalar: 1,
            point: [5, 6],
        };
        let p = Path::new("/sealed");
        let (program, args, pin) = role(p, &cfg, &q, "cms").unwrap();
        assert_eq!(program, Path::new("/sealed/immutable/assets/bin/cms"));
        assert_eq!(pin, native::CMS_SHA);
        assert_eq!(
            args,
            vec![
                "--verb",
                "1",
                "--threads",
                "1",
                "--random",
                "1",
                "--maxsol",
                "1",
                "--maxconfl",
                "100000",
                "instance/instance.xor.cnf"
            ]
        );
        let (_, args, pin) = role(p, &cfg, &q, "exporter").unwrap();
        assert_eq!(pin, native::EXPORTER_SHA);
        assert!(args.ends_with(&["native-ordinary-001".into(), "--export-only".into()]));
        assert!(role(p, &cfg, &q, "rust-cdcl").is_err());
        assert!(role(p, &cfg, &Query { trial: 2, ..q }, "cms").is_err());
    }
}
