//! New SAT target capsule: no target point enters the frozen source/config.
use super::{journal, native, preparation};
use crypto_lib::cryptanalysis::{
    koblitz_fast::{FastCurve, FastPoint},
    koblitz_index_calculus::KoblitzCurve,
    prepared_n17_target::SatTargetPlan,
    prepared_sat_control::{canonical_sha, sha256},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

pub const SCOPE: &str = "disclosed-n17-native-one-target-cms-v1";
pub const WORKER: &str = "prepared_sat_target_worker";
pub const CARD_QUESTION: &str = "fresh-paired-n17-public-point-card-v1";
const CARD_CURVE: &str = "EC1N17Ce1hdfbf24105ef5";
const CARD_LAW: &str = "sha256-seed-first-lift-min-y-cofactor2-v1";

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub schema_version: u32,
    pub question: String,
    pub plan: SatTargetPlan,
    pub worker_timeout_ms: u64,
}
impl Config {
    pub fn validate(&self) -> Result<(), String> {
        self.plan.validate()?;
        native::require(
            self.schema_version == 1
                && self.question == SCOPE
                && (30_000..=900_000).contains(&self.worker_timeout_ms)
                && self.plan.controller_timeout_ms <= self.worker_timeout_ms,
            "outside bounded SAT target worker configuration",
        )
    }
    pub fn journal_binding(&self, seal: &str, worker: &str, target: [u64; 2]) -> journal::Binding {
        journal::Binding {
            schema_version: 1,
            question: journal::QUESTION.into(),
            registration_sha256: seal.into(),
            worker_sha256: worker.into(),
            family: journal::Family::Cryptominisat,
            target,
            algorithm_seed: self.plan.algorithm_seed,
            max_queries: self.plan.max_queries,
        }
    }
}
pub fn config(capsule: &Path) -> Result<Config, String> {
    let cfg: Config = serde_json::from_slice(&native::read(&capsule.join("config.json"), 65536)?)
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
pub fn identity() -> Value {
    json!({"schema_version":1,
        "source_manifest_sha256":option_env!("IC_SAT_TARGET_SOURCE_MANIFEST_SHA256"),
        "build_sha256":option_env!("IC_SAT_TARGET_BUILD_SHA256"),
        "target_os":std::env::consts::OS,"target_arch":std::env::consts::ARCH})
}
pub fn require_identity(compiled: &Value, expected: &Value) -> Result<(), String> {
    native::require(
        compiled == expected
            && compiled["schema_version"] == 1
            && compiled["source_manifest_sha256"]
                .as_str()
                .is_some_and(journal::digest)
            && compiled["build_sha256"]
                .as_str()
                .is_some_and(journal::digest)
            && compiled["target_os"] == "macos"
            && compiled["target_arch"] == "aarch64",
        "SAT target worker lacks matching compiled source/build identity",
    )
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PreparationBinding {
    pub capsule: PathBuf,
    pub execution: PathBuf,
    pub registration_sha256: String,
    pub auditor_sha256: String,
    pub producer_sha256: String,
    pub mathematics_sha256: String,
}
#[derive(Deserialize, Serialize)]
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
    pub prepared_exporter_sha256: String,
    pub prepared_cms_sha256: String,
    pub worker_build_identity: Value,
    pub runtime_environment: BTreeMap<String, String>,
    pub preparation: PreparationBinding,
    pub validation_only: bool,
}
pub fn registration(capsule: &Path, external: &str) -> Result<Registration, String> {
    native::require(
        journal::digest(external),
        "external SAT target seal required",
    )?;
    let bytes = native::read(&capsule.join("registration.json"), 16 * 1024 * 1024)?;
    let record: Registration = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let p = &record.preparation;
    native::require(
        record.schema_version == 1
            && record.scope == SCOPE
            && record.hardware == json!({"os":"macos","architecture":"aarch64"})
            && record.source_commit.len() == 40
            && record
                .source_commit
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            && record.runtime_environment == environment()
            && p.capsule.is_absolute()
            && p.execution.is_absolute()
            && [
                &record.source_manifest_sha256,
                &record.config_sha256,
                &record.host_context_sha256,
                &record.worker_sha256,
                &record.auditor_sha256,
                &record.prepared_exporter_sha256,
                &record.prepared_cms_sha256,
                &p.registration_sha256,
                &p.auditor_sha256,
                &p.producer_sha256,
                &p.mathematics_sha256,
            ]
            .into_iter()
            .all(|value| journal::digest(value)),
        "SAT target registration/source/preparation pins differ",
    )?;
    require_identity(&record.worker_build_identity, &record.worker_build_identity)?;
    native::require(
        record.worker_build_identity["source_manifest_sha256"] == record.source_manifest_sha256
            && canonical_sha(&json!(record))? == external
            && native::load(&capsule.join("seal.json"))?["registration_sha256"] == external
            && sha256(&native::read(&capsule.join("config.json"), 65536)?) == record.config_sha256
            && sha256(&native::read(&capsule.join("host-context.json"), 65536)?)
                == record.host_context_sha256,
        "SAT target registration or sealed sidecars differ",
    )?;
    config(capsule)?;
    Ok(record)
}
pub fn check_capsule(capsule: &Path, external: &str) -> Result<Registration, String> {
    let record = registration(capsule, external)?;
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
        "SAT target source or executable binding differs",
    )?;
    check_roles(capsule, &record)?;
    Ok(record)
}
pub fn check_roles(capsule: &Path, record: &Registration) -> Result<(), String> {
    native::require(
        sha256(&native::read(
            &capsule.join("immutable/assets/bin/prepared-cms"),
            16 * 1024 * 1024,
        )?) == record.prepared_cms_sha256
            && sha256(&native::read(
                &capsule.join("immutable/assets/bin/prepared-exporter"),
                8 * 1024 * 1024,
            )?) == record.prepared_exporter_sha256
            && record.prepared_cms_sha256 != native::CMS_SHA
            && record.prepared_exporter_sha256 != native::EXPORTER_SHA,
        "SAT target prepared role binary differs from frozen pin",
    )
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Card {
    pub schema_version: u32,
    pub question: String,
    pub curve_id: String,
    pub point_law: String,
    pub seed_hex: String,
    pub start_x: u64,
    pub offset: u64,
    pub lift_x: u64,
    pub target: [u64; 2],
    pub source_publications: BTreeMap<String, String>,
    pub scalar_generated: bool,
}
pub fn card(path: &Path) -> Result<(Card, String), String> {
    let bytes = native::read(path, 65536)?;
    let value: Card = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let curve = KoblitzCurve::new(1, 17).ok_or("n17 SAT target curve unavailable")?;
    let fast = FastCurve::new(&curve.curve).ok_or("native SAT target curve unavailable")?;
    let point = FastPoint::affine(value.target[0], value.target[1]);
    native::require(
        value.schema_version == 1
            && value.question == CARD_QUESTION
            && value.curve_id == CARD_CURVE
            && value.point_law == CARD_LAW
            && !value.scalar_generated
            && value.seed_hex.len() == 64
            && value
                .seed_hex
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            && value.start_x < 1 << 17
            && value.offset < 1 << 17
            && value.lift_x == (value.start_x + value.offset) % (1 << 17)
            && value.target.iter().all(|&v| v < 1 << 17)
            && value.source_publications.len() == 4
            && ["cms", "f5", "incumbent", "rho"].into_iter().all(|name| {
                value
                    .source_publications
                    .get(name)
                    .is_some_and(|digest| journal::digest(digest))
            })
            && fast.is_on_curve(point)
            && fast.mul_u64(point, 65587).infinity,
        "SAT target card is malformed or outside exact nonidentity subgroup",
    )?;
    native::require(
        native::read(path, 65536)? == bytes,
        "SAT target card changed during read",
    )?;
    Ok((value, sha256(&bytes)))
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Claim {
    pub schema_version: u32,
    pub scope: String,
    pub registration_sha256: String,
    pub execution: PathBuf,
    pub worker_sha256: String,
    pub card: PathBuf,
    pub card_sha256: String,
    pub publication_sha256: String,
}
pub fn claim(
    capsule: &Path,
    execution: &Path,
    card_path: &Path,
    record: &Registration,
    seal: &str,
) -> Result<(Card, String), String> {
    let claimed: Claim =
        serde_json::from_slice(&native::read(&capsule.join("consumed.json"), 65536)?)
            .map_err(|e| e.to_string())?;
    let card_path = card_path.canonicalize().map_err(|e| e.to_string())?;
    let (card, hash) = card(&card_path)?;
    native::require(
        !record.validation_only
            && claimed.schema_version == 1
            && claimed.scope == SCOPE
            && claimed.registration_sha256 == seal
            && claimed.worker_sha256 == record.worker_sha256
            && claimed.execution == execution
            && claimed.card == card_path
            && claimed.card_sha256 == hash
            && card.source_publications["cms"] == claimed.publication_sha256,
        "SAT target one-use claim/card/publication differs",
    )?;
    Ok((card, hash))
}
pub fn preparation_paths(
    binding: &PreparationBinding,
) -> Result<(PathBuf, PathBuf, PathBuf), String> {
    let root = binding.capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = binding
        .execution
        .canonicalize()
        .map_err(|e| e.to_string())?;
    native::require(
        root == binding.capsule && execution == binding.execution,
        "SAT original preparation paths are not canonical",
    )?;
    let (record, hash) = preparation::check_capsule(&root)?;
    preparation::claim(&root, &execution, &record, &hash)?;
    native::require(
        hash == binding.registration_sha256
            && record.auditor_sha256 == binding.auditor_sha256
            && preparation::config(&root)?.plan.family
                == crypto_lib::cryptanalysis::prepared_ordinary::Family::Cryptominisat,
        "SAT target preparation family/source/auditor differs",
    )?;
    let auditor = root.join("immutable/bin/icprog");
    native::require(
        sha256(&native::read(&auditor, 128 * 1024 * 1024)?) == binding.auditor_sha256,
        "SAT original preparation auditor changed",
    )?;
    Ok((root, execution, auditor))
}
pub fn preparation_receipt(
    binding: &PreparationBinding,
    receipt: &Value,
    producer: &Value,
    producer_bytes: &[u8],
) -> Result<Value, String> {
    let math = &producer["mathematical_input"];
    native::require(
        receipt["schema_version"] == 1
            && receipt["status"] == "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT"
            && receipt["registration_sha256"] == binding.registration_sha256
            && receipt["checker_sha256"] == binding.auditor_sha256
            && receipt["source_bound_execution_admitted"] == true
            && receipt["mathematical_preparation_complete"] == true
            && receipt["producer_sha256"] == binding.producer_sha256
            && sha256(producer_bytes) == binding.producer_sha256
            && receipt["mathematics"]["input_sha256"] == binding.mathematics_sha256
            && canonical_sha(math)? == binding.mathematics_sha256
            && receipt["mathematics"]["declared_family"] == "cryptominisat"
            && receipt["mathematics"]["rank"] == 29
            && receipt["mathematics"]["folded_columns"] == 29
            && receipt["mathematics"]["usable_points"] == 62
            && receipt["mathematics"]["panel_complete"] == true,
        "SAT preparation recheck lacks complete source-bound rank/log evidence",
    )?;
    Ok(math.clone())
}

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Query {
    pub trial: usize,
    pub point: [u64; 2],
}
pub struct PreparedRole {
    pub program: PathBuf,
    pub args: Vec<String>,
    pub pin: String,
    pub marker: &'static str,
}
pub fn prepared_role(
    capsule: &Path,
    cfg: &Config,
    trial: usize,
    role: &str,
    exporter_pin: &str,
    cms_pin: &str,
) -> Result<PreparedRole, String> {
    cfg.validate()?;
    native::require(
        trial < cfg.plan.max_queries && journal::digest(exporter_pin) && journal::digest(cms_pin),
        "SAT prepared role outside frozen plan or pin",
    )?;
    match role {
        "exporter" => Ok(PreparedRole {
            program: capsule.join("immutable/assets/bin/prepared-exporter"),
            args: vec![
                "17".into(),
                "6".into(),
                "standard".into(),
                cfg.plan.export_nonce.to_string(),
                cfg.plan.conflict_budget.to_string(),
                "instance".into(),
                "1".into(),
                "0".into(),
                "--export-only".into(),
                "--prestart-stdin".into(),
            ],
            pin: exporter_pin.into(),
            marker: "c EXPORTER_PREPARED_STDIN_READY_v1",
        }),
        "cms" => Ok(PreparedRole {
            program: capsule.join("immutable/assets/bin/prepared-cms"),
            args: vec![
                "--verb".into(),
                "1".into(),
                "--threads".into(),
                "1".into(),
                "--random".into(),
                "1".into(),
                "--maxsol".into(),
                "1".into(),
                "--maxconfl".into(),
                cfg.plan.conflict_budget.to_string(),
            ],
            pin: cms_pin.into(),
            marker: "c PREPARED_STDIN_READY_v1",
        }),
        _ => Err("unsupported SAT target native role".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn cfg() -> Config {
        Config {
            schema_version: 1,
            question: SCOPE.into(),
            plan: SatTargetPlan {
                schema_version: 1,
                question: "prepared-public-target-n17-cms-v1".into(),
                algorithm_seed: 7,
                max_queries: 8,
                export_nonce: 11,
                conflict_budget: 1_000_000,
                exporter_timeout_ms: 30_000,
                solver_timeout_ms: 60_000,
                controller_timeout_ms: 900_000,
            },
            worker_timeout_ms: 900_000,
        }
    }

    #[test]
    fn no_target_enters_source_config_and_roles_use_public_query() {
        let mut config = cfg();
        config.validate().unwrap();
        let mut source = json!(config);
        source["target"] = json!([61889, 74818]);
        assert!(serde_json::from_value::<Config>(source).is_err());
        config.worker_timeout_ms = 899_999;
        assert!(config.validate().is_err());
        let config = cfg();
        let exporter = prepared_role(
            Path::new("/sealed"),
            &config,
            0,
            "exporter",
            &"a".repeat(64),
            &"b".repeat(64),
        )
        .unwrap();
        assert_eq!(exporter.pin, "a".repeat(64));
        assert_eq!(exporter.args.last().unwrap(), "--prestart-stdin");
        assert!(!exporter
            .args
            .iter()
            .any(|arg| arg.contains("61889") || arg.contains("74818")));
        let cms = prepared_role(
            Path::new("/sealed"),
            &config,
            0,
            "cms",
            &"a".repeat(64),
            &"b".repeat(64),
        )
        .unwrap();
        assert_eq!(cms.pin, "b".repeat(64));
        assert_eq!(cms.marker, "c PREPARED_STDIN_READY_v1");
        assert!(prepared_role(
            Path::new("/sealed"),
            &config,
            0,
            "f5",
            &"a".repeat(64),
            &"b".repeat(64)
        )
        .is_err());
    }

    #[test]
    fn card_rejects_nonpoint_and_missing_source_publication() {
        let root = std::env::temp_dir().join(format!(
            "sat-target-card-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let path = root.join("card-valid.json");
        let sources = ["cms", "f5", "incumbent", "rho"]
            .into_iter()
            .map(|name| (name.to_string(), "a".repeat(64)))
            .collect();
        let mut value = Card {
            schema_version: 1,
            question: CARD_QUESTION.into(),
            curve_id: CARD_CURVE.into(),
            point_law: CARD_LAW.into(),
            seed_hex: "b".repeat(64),
            start_x: 0,
            offset: 0,
            lift_x: 0,
            target: [61889, 74818],
            source_publications: sources,
            scalar_generated: false,
        };
        native::save(&path, &json!(value)).unwrap();
        card(&path).unwrap();
        value.target = [0, 0];
        let invalid_point = root.join("card-invalid-point.json");
        native::save(&invalid_point, &json!(value)).unwrap();
        assert!(card(&invalid_point).is_err());
        value.target = [61889, 74818];
        value.source_publications.remove("rho");
        let missing_source = root.join("card-missing-source.json");
        native::save(&missing_source, &json!(value)).unwrap();
        assert!(card(&missing_source).is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
