//! Native, scalar-free point cards for a later fresh paired n17 holdout.
//! Source-publication bytes and point arithmetic are replayed as data only.
//! An outer frozen controller must establish publication order and runtime
//! custody before a card can qualify a one-target comparison.
use super::{json, oracle::Curve, sat_control::native, target_math};
use crypto_lib::cryptanalysis::{
    ecbench::workload::{CurveSpec, TargetKind, Workload, PUBLIC_TARGET_LAW},
    prepared_sat_control::{canonical_sha, sha256},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

const QUESTION: &str = "fresh-paired-n17-ecbench-public-point-card-v1";
const SOURCE_QUESTION: &str = "fresh-paired-n17-source-publication-v1";
// The tournament's canonical decimal polynomial-bit record for this fixture.
// The older curve registry also has an e1/hex-encoding alias; it is not the
// identity used by ecbench's paired IC/rho candidate and workload records.
const CURVE_ID: &str = "EC1N17Ckb1hbbe2b5b6b1e6";
const ROLES: [&str; 4] = ["cms", "f5", "incumbent", "rho"];

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SourceDescriptor {
    schema_version: u32,
    question: String,
    role: String,
    curve_id: String,
    registration_sha256: String,
    source_manifest_sha256: String,
    archive_sha256: String,
    target_free: bool,
    source_bound_execution_admitted: bool,
}

fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn source_pins(paths: &BTreeMap<String, PathBuf>) -> Result<BTreeMap<String, String>, String> {
    native::require(
        paths.len() == ROLES.len()
            && ROLES.iter().all(|role| paths.contains_key(*role))
            && paths.values().all(|path| path.is_absolute()),
        "point card requires four absolute published source descriptors",
    )?;
    paths
        .iter()
        .map(|(role, path)| {
            let bytes = native::read(path, 65536)?;
            let descriptor: SourceDescriptor =
                serde_json::from_value(target_math::parse(&bytes)?).map_err(|e| e.to_string())?;
            native::require(
                descriptor.schema_version == 1
                    && descriptor.question == SOURCE_QUESTION
                    && descriptor.role.as_str() == role.as_str()
                    && descriptor.curve_id == CURVE_ID
                    && digest(&descriptor.registration_sha256)
                    && digest(&descriptor.source_manifest_sha256)
                    && digest(&descriptor.archive_sha256)
                    && descriptor.target_free
                    && !descriptor.source_bound_execution_admitted,
                "source descriptor differs from target-free sealed role and curve",
            )?;
            Ok((role.clone(), sha256(&bytes)))
        })
        .collect()
}
fn paths(path: &Path) -> Result<BTreeMap<String, PathBuf>, String> {
    let bytes = native::read(path, 65536)?;
    let map: BTreeMap<String, PathBuf> =
        serde_json::from_value(target_math::parse(&bytes)?).map_err(|e| e.to_string())?;
    source_pins(&map)?;
    native::require(
        native::read(path, 65536)? == bytes,
        "point-card publication list changed during read",
    )?;
    Ok(map)
}

#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Card {
    schema_version: u32,
    question: String,
    curve_id: String,
    point_law: String,
    target_seed: u64,
    target_index: u64,
    target_counter: u64,
    ecbench_workload_id: String,
    ecbench_workload_sha256: String,
    target: [u64; 2],
    source_publications: BTreeMap<String, String>,
    scalar_generated: bool,
}
fn fixture() -> Value {
    json!({"degree":17,"curve_a":1,"subgroup_order":65587,"group_order":131174,
        "cofactor":2,"generator":[43693,23339],"lambda":17184,
        "irreducible":{"degree":17,"low_terms":[0,3]},"targets":[],
        "target_seeds":[],"target_scalar_constructed":false})
}
fn independently_check(point: [u64; 2]) -> Result<(), String> {
    let curve = Curve::new(&json::parse(&fixture().to_string())?)?;
    let decoded = curve.decode(&json::parse(&json!(point).to_string())?)?;
    native::require(
        decoded.is_some() && curve.mul(decoded, 65587).is_none(),
        "derived public point is outside the nonidentity subgroup",
    )
}
fn native_workload(seed: u64) -> Result<(Workload, [u64; 2]), String> {
    let spec = CurveSpec::Koblitz { a: 1, n: 17 };
    let instance = spec.build()?;
    let workload = Workload::on(&spec, &instance, seed, 0, TargetKind::Public)?;
    native::require(
        workload.curve.family == "koblitz"
            && workload.curve.field_degree == Some(17)
            && workload.curve.group_order == 131174
            && workload.curve.r == 65587
            && workload.curve.cofactor == 2
            && workload.curve.generator[0] == "0xaaad"
            && workload.curve.generator[1] == "0x5b2b"
            && workload.target_law == PUBLIC_TARGET_LAW
            && workload.planted.is_none()
            && workload.target_counter.is_some(),
        "ecbench public workload differs from exact n17 subgroup",
    )?;
    let parse = |hex: &str| {
        u64::from_str_radix(
            hex.strip_prefix("0x")
                .ok_or_else(|| "noncanonical point hex".to_string())?,
            16,
        )
        .map_err(|_| "invalid point coordinate hex".to_string())
    };
    let target = [parse(&workload.target[0])?, parse(&workload.target[1])?];
    independently_check(target)?;
    Ok((workload, target))
}
fn derive(seed: u64, sources: BTreeMap<String, String>) -> Result<Card, String> {
    native::require(
        sources.len() == ROLES.len()
            && ROLES
                .iter()
                .all(|role| sources.get(*role).is_some_and(|value| digest(value))),
        "point card source digests are incomplete",
    )?;
    let (workload, target) = native_workload(seed)?;
    Ok(Card {
        schema_version: 1,
        question: QUESTION.into(),
        curve_id: CURVE_ID.into(),
        point_law: PUBLIC_TARGET_LAW.into(),
        target_seed: seed,
        target_index: 0,
        target_counter: workload
            .target_counter
            .ok_or("public target counter absent")?,
        ecbench_workload_id: workload.workload_id,
        ecbench_workload_sha256: workload.workload_sha256,
        target,
        source_publications: sources,
        scalar_generated: false,
    })
}
fn replay(card: &Card, sources: BTreeMap<String, String>) -> Result<(), String> {
    native::require(
        card.schema_version == 1
            && card.question == QUESTION
            && card.curve_id == CURVE_ID
            && card.point_law == PUBLIC_TARGET_LAW
            && !card.scalar_generated
            && card.source_publications == sources
            && *card == derive(card.target_seed, sources)?,
        "public point or published source pins differ from the frozen card",
    )
}
pub(super) fn generate(publications: &Path, out: &Path) -> Result<String, String> {
    let paths = paths(publications)?;
    let before = source_pins(&paths)?;
    let mut seed = [0u8; 8];
    getrandom::getrandom(&mut seed).map_err(|e| e.to_string())?;
    let card = derive(u64::from_le_bytes(seed), before.clone())?;
    native::require(
        source_pins(&paths)? == before,
        "source publication changed during point generation",
    )?;
    native::save(out, &json!(card))?;
    let bytes = native::read(out, 65536)?;
    serde_json::to_string_pretty(
        &json!({"schema_version":1,"status":"PUBLIC_POINT_CARD_CREATED",
        "card_sha256":sha256(&bytes),"card_identity_sha256":canonical_sha(&json!(card))?,
        "target":card.target,"ecbench_workload_id":card.ecbench_workload_id,
        "ecbench_workload_sha256":card.ecbench_workload_sha256,
        "source_publications":before,"scalar_generated":false,
        "publication_order_and_freshness_certified":false,"solver_calls":0,
        "source_bound_execution_admitted":false,"campaign_workload_id":null,
        "online_speedup":null}),
    )
    .map_err(|e| e.to_string())
}
pub(super) fn audit(card_path: &Path, publications: &Path, out: &Path) -> Result<String, String> {
    let paths = paths(publications)?;
    let sources = source_pins(&paths)?;
    let bytes = native::read(card_path, 65536)?;
    let card: Card =
        serde_json::from_value(target_math::parse(&bytes)?).map_err(|e| e.to_string())?;
    replay(&card, sources.clone())?;
    native::require(
        native::read(card_path, 65536)? == bytes && source_pins(&paths)? == sources,
        "public card or publication descriptors changed during audit",
    )?;
    let result = json!({"schema_version":1,"status":"PASS_PUBLIC_POINT_CARD_DATA_ONLY_AUDIT",
        "card_sha256":sha256(&bytes),"card_identity_sha256":canonical_sha(&json!(card))?,
        "target":card.target,"curve_id":CURVE_ID,
        "ecbench_workload_id":card.ecbench_workload_id,
        "ecbench_workload_sha256":card.ecbench_workload_sha256,
        "source_publications":sources,"scalar_generated":false,"solver_calls_by_auditor":0,
        "publication_order_and_freshness_certified":false,
        "source_bound_execution_admitted":false,"campaign_workload_id":null,
        "online_speedup":null});
    native::save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};
    #[test]
    fn native_point_card_is_scalar_free_and_detects_point_or_source_changes() {
        let sources = ROLES
            .iter()
            .enumerate()
            .map(|(i, role)| (role.to_string(), format!("{i:064x}")))
            .collect::<BTreeMap<_, _>>();
        let card = derive(7, sources.clone()).unwrap();
        assert_eq!(card.curve_id, CURVE_ID);
        assert_eq!(card.point_law, PUBLIC_TARGET_LAW);
        assert_eq!(card.target_index, 0);
        assert!(card.ecbench_workload_id.starts_with('W'));
        assert_eq!(card.ecbench_workload_id.len(), 13);
        assert!(digest(&card.ecbench_workload_sha256));
        assert!(!card.scalar_generated);
        replay(&card, sources.clone()).unwrap();
        let mut changed = card.clone();
        changed.target[0] ^= 1;
        assert!(replay(&changed, sources.clone()).is_err());
        let mut changed = card.clone();
        changed.ecbench_workload_id = "W000000000000".into();
        assert!(replay(&changed, sources.clone()).is_err());
        let mut changed = sources;
        changed.insert("cms".into(), "f".repeat(64));
        assert!(replay(&card, changed).is_err());
    }

    #[test]
    fn card_binds_original_publication_bytes_without_a_solver_or_known_scalar() {
        let root = std::env::temp_dir().join(format!(
            "ic-target-card-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&root).unwrap();
        let root = root.canonicalize().unwrap();
        let mut paths = BTreeMap::new();
        for (index, role) in ROLES.iter().enumerate() {
            let path = root.join(format!("{role}-publication.json"));
            native::save(
                &path,
                &json!({"schema_version":1,"question":SOURCE_QUESTION,
                    "role":role,"curve_id":CURVE_ID,
                    "registration_sha256":format!("{index:064x}"),
                    "source_manifest_sha256":format!("{:064x}",index+10),
                    "archive_sha256":format!("{:064x}",index+20),
                    "target_free":true,"source_bound_execution_admitted":false}),
            )
            .unwrap();
            paths.insert((*role).to_string(), path);
        }
        let listing = root.join("publications.json");
        native::save(&listing, &json!(paths)).unwrap();
        let card = root.join("card.json");
        let created: Value = serde_json::from_str(&generate(&listing, &card).unwrap()).unwrap();
        let audited: Value =
            serde_json::from_str(&audit(&card, &listing, &root.join("audit.json")).unwrap())
                .unwrap();
        assert_eq!(created["card_sha256"], audited["card_sha256"]);
        assert_eq!(
            created["ecbench_workload_id"],
            audited["ecbench_workload_id"]
        );
        assert_eq!(audited["scalar_generated"], false);
        assert_eq!(audited["solver_calls_by_auditor"], 0);
        assert_eq!(audited["publication_order_and_freshness_certified"], false);
        std::fs::write(
            paths.get("cms").unwrap(),
            serde_json::to_vec(&json!({"schema_version":1,"question":SOURCE_QUESTION,
                "role":"cms","curve_id":CURVE_ID,
                "registration_sha256":"f".repeat(64),
                "source_manifest_sha256":format!("{:064x}",10),
                "archive_sha256":format!("{:064x}",20),
                "target_free":true,"source_bound_execution_admitted":false}))
            .unwrap(),
        )
        .unwrap();
        assert!(audit(&card, &listing, &root.join("tampered-audit.json")).is_err());
        let wrong_curve = root.join("wrong-curve-publication.json");
        native::save(
            &wrong_curve,
            &json!({"schema_version":1,"question":SOURCE_QUESTION,
                "role":"cms","curve_id":"EC1N17Ce1hdfbf24105ef5",
                "registration_sha256":"f".repeat(64),
                "source_manifest_sha256":format!("{:064x}",10),
                "archive_sha256":format!("{:064x}",20),
                "target_free":true,"source_bound_execution_admitted":false}),
        )
        .unwrap();
        paths.insert("cms".into(), wrong_curve);
        native::save(&root.join("wrong-curve-list.json"), &json!(paths)).unwrap();
        assert!(super::paths(&root.join("wrong-curve-list.json")).is_err());
        paths.insert("cms".into(), root.join("cms-publication.json"));
        let mut leaked: Value =
            target_math::parse(&native::read(paths.get("f5").unwrap(), 65536).unwrap()).unwrap();
        leaked["target"] = json!([43693, 23339]);
        native::save(&root.join("target-leak.json"), &leaked).unwrap();
        paths.insert("f5".into(), root.join("target-leak.json"));
        native::save(&root.join("target-leak-list.json"), &json!(paths)).unwrap();
        assert!(super::paths(&root.join("target-leak-list.json")).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }
}
