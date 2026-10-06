//! Native, scalar-free point cards for a later fresh paired n17 holdout.
//! Source-publication bytes and point arithmetic are replayed as data only.
//! An outer frozen controller must establish publication order and runtime
//! custody before a card can qualify a one-target comparison.
use super::{json, oracle::Curve, sat_control::native, target_math};
use crypto_lib::{
    binary_ecc::{BinaryPoint, F2mElement},
    cryptanalysis::{
        koblitz_fast::FastCurve,
        koblitz_index_calculus::{points_with_x_fast, KoblitzCurve},
        prepared_sat_control::{canonical_sha, sha256},
    },
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

const QUESTION: &str = "fresh-paired-n17-public-point-card-v1";
const CURVE_ID: &str = "EC1N17Ce1hdfbf24105ef5";
const POINT_LAW: &str = "sha256-seed-first-lift-min-y-cofactor2-v1";
const DOMAIN: &[u8] = b"IC-fresh-paired-n17-point-v1\0";
const ROLES: [&str; 4] = ["cms", "f5", "incumbent", "rho"];

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
            let bytes = native::read(path, 16 * 1024 * 1024)?;
            let descriptor = target_math::parse(&bytes)?;
            native::require(
                descriptor["registration_sha256"]
                    .as_str()
                    .is_some_and(digest)
                    && descriptor["source_bound_execution_admitted"] == false,
                "source descriptor is not a sealed unexecuted registration",
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
    seed_hex: String,
    start_x: u64,
    offset: u64,
    lift_x: u64,
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
fn first_point(seed: &[u8; 32]) -> Result<(u64, u64, [u64; 2]), String> {
    let mut input = DOMAIN.to_vec();
    input.extend_from_slice(seed);
    let digest = sha256(&input);
    let start_x = u64::from_str_radix(&digest[..8], 16)
        .map_err(|_| "point-seed digest is invalid")?
        % (1 << 17);
    let curve = KoblitzCurve::new(1, 17).ok_or("n17 point-card curve unavailable")?;
    native::require(
        curve.n == 17
            && curve.a == 1
            && curve.k == 1
            && curve.subgroup_order.to_u64() == Some(65587)
            && curve.cofactor.to_u64() == Some(2),
        "point-card curve differs from declared n17 subgroup",
    )?;
    let fast = FastCurve::new(&curve.curve).ok_or("n17 native curve unavailable")?;
    for offset in 0..(1u64 << 17) {
        let x = (start_x + offset) % (1 << 17);
        let element = F2mElement::from_biguint(&BigUint::from(x), 17);
        let chosen = points_with_x_fast(&fast, &element)
            .into_iter()
            .filter_map(|point| match &point {
                BinaryPoint::Affine { y, .. } => {
                    y.to_biguint().to_u64().map(|value| (value, point))
                }
                BinaryPoint::Infinity => None,
            })
            .min_by_key(|(y, _)| *y);
        let Some((_, point)) = chosen else {
            continue;
        };
        let subgroup = fast.mul_u64(fast.lift(&point), 2);
        if subgroup.infinity {
            continue;
        }
        let target = [subgroup.x, subgroup.y];
        independently_check(target)?;
        return Ok((start_x, offset, target));
    }
    Err("no nonidentity subgroup point in the declared x search".into())
}
fn derive(seed: &[u8; 32], sources: BTreeMap<String, String>) -> Result<Card, String> {
    native::require(
        sources.len() == ROLES.len()
            && ROLES
                .iter()
                .all(|role| sources.get(*role).is_some_and(|value| digest(value))),
        "point card source digests are incomplete",
    )?;
    let (start_x, offset, target) = first_point(seed)?;
    Ok(Card {
        schema_version: 1,
        question: QUESTION.into(),
        curve_id: CURVE_ID.into(),
        point_law: POINT_LAW.into(),
        seed_hex: hex::encode(seed),
        start_x,
        offset,
        lift_x: (start_x + offset) % (1 << 17),
        target,
        source_publications: sources,
        scalar_generated: false,
    })
}
fn replay(card: &Card, sources: BTreeMap<String, String>) -> Result<(), String> {
    let seed = hex::decode(&card.seed_hex).map_err(|e| e.to_string())?;
    let seed: [u8; 32] = seed.try_into().map_err(|_| "card seed is not 32 bytes")?;
    native::require(
        card.schema_version == 1
            && card.question == QUESTION
            && card.curve_id == CURVE_ID
            && card.point_law == POINT_LAW
            && !card.scalar_generated
            && card.source_publications == sources
            && *card == derive(&seed, sources)?,
        "public point or published source pins differ from the frozen card",
    )
}
pub(super) fn generate(publications: &Path, out: &Path) -> Result<String, String> {
    let paths = paths(publications)?;
    let before = source_pins(&paths)?;
    let mut seed = [0u8; 32];
    getrandom::getrandom(&mut seed).map_err(|e| e.to_string())?;
    let card = derive(&seed, before.clone())?;
    native::require(
        source_pins(&paths)? == before,
        "source publication changed during point generation",
    )?;
    native::save(out, &json!(card))?;
    let bytes = native::read(out, 65536)?;
    serde_json::to_string_pretty(
        &json!({"schema_version":1,"status":"PUBLIC_POINT_CARD_CREATED",
        "card_sha256":sha256(&bytes),"card_identity_sha256":canonical_sha(&json!(card))?,
        "target":card.target,
        "source_publications":before,"scalar_generated":false,
        "publication_order_and_freshness_certified":false,"solver_calls":0,
        "source_bound_execution_admitted":false,"online_speedup":null}),
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
        "source_publications":sources,"scalar_generated":false,"solver_calls_by_auditor":0,
        "publication_order_and_freshness_certified":false,
        "source_bound_execution_admitted":false,"online_speedup":null});
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
        let card = derive(&[7u8; 32], sources.clone()).unwrap();
        assert_eq!(card.curve_id, CURVE_ID);
        assert!(!card.scalar_generated);
        replay(&card, sources.clone()).unwrap();
        let mut changed = card.clone();
        changed.target[0] ^= 1;
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
                &json!({"registration_sha256":format!("{index:064x}"),
                    "source_bound_execution_admitted":false}),
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
        assert_eq!(audited["scalar_generated"], false);
        assert_eq!(audited["solver_calls_by_auditor"], 0);
        assert_eq!(audited["publication_order_and_freshness_certified"], false);
        std::fs::write(
            paths.get("cms").unwrap(),
            serde_json::to_vec(&json!({"registration_sha256":"f".repeat(64),
                "source_bound_execution_admitted":false}))
            .unwrap(),
        )
        .unwrap();
        assert!(audit(&card, &listing, &root.join("tampered-audit.json")).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }
}
