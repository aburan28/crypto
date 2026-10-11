//! Frozen producer-side Role-1 controls.

use std::collections::{BTreeMap, BTreeSet};

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde_json::{json, Value};

use super::arithmetic::{
    d23_control, hnf_from_vectors, integer, multiply_hnf, norm, p192_parameters,
    p192_sqrt_discriminant_control, pow_hnf, principal_hnf, Hnf,
};
use super::certificate::{encoded_sint, residual_prime_admissible};
use super::digest::sha256_hex;
use super::evidence::curve_object;
use super::factor_base::{characteristic_roots, is_characteristic_root, Entry};
use super::primality;
use super::provenance::BuildProvenance;
use super::ref0::{
    build_reference_segmented, build_reference_with_order, ReferenceOutput, RECORD_COUNT,
};
use super::schema::{canonical_json_bytes, CURVE_UID, EXPERIMENT_ID, PROTOCOL_VERSION};
use super::Result;

fn assertion(assertion_id: &str, passed: bool) -> Value {
    json!({"assertion_id":assertion_id,"passed":passed})
}

fn control(control_id: &str, values: &[(&str, bool)]) -> Value {
    let assertions = values
        .iter()
        .map(|(assertion_id, passed)| assertion(assertion_id, *passed))
        .collect::<Vec<_>>();
    json!({
        "assertions": assertions,
        "control_id": control_id,
        "status": if values.iter().all(|(_, passed)| *passed) { "PASS" } else { "FAIL" },
    })
}

fn split_dispositions(bytes: &[u8]) -> Result<Vec<&[u8]>> {
    let mut output = Vec::new();
    let mut offset = 0usize;
    for expected_index in 0..RECORD_COUNT {
        if bytes.len().saturating_sub(offset) < 10 {
            return Err("truncated REF-0 disposition control input".to_owned());
        }
        let status = bytes[offset];
        let index = u64::from_be_bytes(
            bytes[offset + 1..offset + 9]
                .try_into()
                .map_err(|_| "disposition index width".to_owned())?,
        );
        let flag = bytes[offset + 9];
        if index != expected_index || status > 6 || flag > 1 {
            return Err("noncanonical REF-0 disposition control input".to_owned());
        }
        let retained = matches!(status, 1..=3);
        if retained != (flag == 1) {
            return Err("REF-0 disposition certificate flag mismatch".to_owned());
        }
        let length = if retained { 74 } else { 10 };
        if bytes.len().saturating_sub(offset) < length {
            return Err("truncated REF-0 disposition hashes".to_owned());
        }
        output.push(&bytes[offset..offset + length]);
        offset += length;
    }
    if offset != bytes.len() {
        return Err("trailing REF-0 disposition bytes".to_owned());
    }
    Ok(output)
}

fn certificate_store_count(bytes: &[u8]) -> Result<u64> {
    const PREFIX: &[u8] = b"P192-WCM-REF-CERT-STORE-v1\0";
    if bytes.len() < PREFIX.len() + 8 || &bytes[..PREFIX.len()] != PREFIX {
        return Err("invalid REF-0 certificate-store prefix".to_owned());
    }
    Ok(u64::from_be_bytes(
        bytes[PREFIX.len()..PREFIX.len() + 8]
            .try_into()
            .map_err(|_| "certificate-store count width".to_owned())?,
    ))
}

fn same_reference(left: &ReferenceOutput, right: &ReferenceOutput) -> bool {
    left.candidate_records == right.candidate_records
        && left.candidate_shard_sha256 == right.candidate_shard_sha256
        && left.candidate_shell_sha256 == right.candidate_shell_sha256
        && left.disposition_records == right.disposition_records
        && left.disposition_shard_sha256 == right.disposition_shard_sha256
        && left.disposition_shell_sha256 == right.disposition_shell_sha256
        && left.certificates == right.certificates
        && left.status_counts == right.status_counts
}

fn endpoint_allowed(
    candidate: u64,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> bool {
    candidate <= bound
        && primality::is_prime(candidate)
        && residual_prime_admissible(candidate, algebraic_bound, trace, field_norm, discriminant)
}

fn residual_status(
    residual: &BigUint,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<u8> {
    let bound_squared = BigUint::from(bound) * bound;
    if residual > &bound_squared {
        return Ok(4);
    }
    let value = residual
        .to_u64()
        .ok_or_else(|| "bounded LP control residual does not fit u64".to_owned())?;
    let factors = primality::factor(value)?;
    let multiplicity = factors.iter().map(|(_, exponent)| *exponent).sum::<u32>();
    if factors.iter().any(|(prime, _)| {
        !endpoint_allowed(
            *prime,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )
    }) {
        return Ok(4);
    }
    Ok(match multiplicity {
        1 => 2,
        2 => 3,
        _ => 4,
    })
}

#[derive(Clone)]
struct D23Fixture {
    selected_root: u64,
    signed_coordinate: BigInt,
    alpha_u: BigInt,
    alpha_v: BigInt,
    claimed_norm: BigInt,
    rational_factor_exponent: u32,
    ideal_power_top_right: BigInt,
    signed_coordinate_encoding: Vec<u8>,
}

fn d23_fixture() -> D23Fixture {
    D23Fixture {
        selected_root: 1,
        signed_coordinate: BigInt::from(-3),
        alpha_u: BigInt::one(),
        alpha_v: BigInt::one(),
        claimed_norm: BigInt::from(8u8),
        rational_factor_exponent: 3,
        ideal_power_top_right: BigInt::from(-7),
        signed_coordinate_encoding: vec![0x02, 0, 0, 0, 1, 3],
    }
}

fn validate_d23_fixture(fixture: &D23Fixture) -> Result<bool> {
    let p = BigInt::from(6u8);
    let t = BigInt::one();
    if !is_characteristic_root(fixture.selected_root, 2, &t, &p)
        || (&fixture.alpha_u + &fixture.alpha_v * BigInt::from(fixture.selected_root))
            .mod_floor(&BigInt::from(2u8))
            != BigInt::zero()
    {
        return Ok(false);
    }
    let expected_coordinate = if fixture.selected_root == 0 {
        BigInt::from(fixture.rational_factor_exponent)
    } else {
        -BigInt::from(fixture.rational_factor_exponent)
    };
    if fixture.signed_coordinate != expected_coordinate
        || encoded_sint(&fixture.signed_coordinate)? != fixture.signed_coordinate_encoding
        || norm(&fixture.alpha_u, &fixture.alpha_v, &t, &p) != fixture.claimed_norm
        || BigInt::from(2u8).pow(fixture.rational_factor_exponent) != fixture.claimed_norm
    {
        return Ok(false);
    }
    let ideal_power = pow_hnf(
        &Hnf::prime(2, fixture.selected_root)?,
        fixture.rational_factor_exponent,
        &p,
        &t,
    )?;
    let recorded_hnf = hnf_from_vectors(&[
        (BigInt::from(8u8), BigInt::zero()),
        (fixture.ideal_power_top_right.clone(), BigInt::one()),
    ])?;
    Ok(ideal_power == recorded_hnf
        && ideal_power == principal_hnf(&fixture.alpha_u, &fixture.alpha_v, &p, &t)?)
}

#[derive(Clone)]
struct LpFixture {
    q: u64,
    smaller_root: u64,
    larger_root: u64,
    sign: u8,
    multiplicity: u32,
}

fn validate_lp_fixture(
    fixture: &LpFixture,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> bool {
    fixture.smaller_root < fixture.larger_root
        && fixture.sign == 1
        && fixture.multiplicity == 1
        && is_characteristic_root(fixture.smaller_root, fixture.q, trace, field_norm)
        && is_characteristic_root(fixture.larger_root, fixture.q, trace, field_norm)
        && endpoint_allowed(
            fixture.q,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )
}

const MUTATION_FILE_SCHEMA: &str = "p192-wcm-implementation-mutation-file-v1";
const MUTATION_RECORD_DOMAIN: &[u8] = b"P192-WCM-MUTATION-RECORD-v1\0";
const MUTATION_FILE_DOMAIN: &[u8] = b"P192-WCM-MUTATION-FILE-v1\0";

#[derive(Clone, Debug)]
struct MutationCase {
    id: &'static str,
    fixture_kind: &'static str,
    target_field_path: &'static str,
    carrier_field_path: &'static str,
    before: Value,
    after: Value,
    baseline: Value,
    mutated: Value,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MutationHashFailure {
    Envelope,
    RecordDigest,
    FileDigest,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct MutationReplayEvidence {
    cases: u64,
    mutation_identities: Vec<(&'static str, &'static str)>,
    exact_single_scalar_changes: u64,
    stale_record_hash_rejections: u64,
    stale_file_hash_rejections: u64,
    recomputed_hash_chain_acceptances: u64,
    semantic_rejections_after_hash_acceptance: u64,
}

impl MutationReplayEvidence {
    fn complete(&self) -> bool {
        self.cases == 9
            && self.mutation_identities.as_slice() == EXPECTED_MUTATION_IDENTITIES
            && self.exact_single_scalar_changes == 9
            && self.stale_record_hash_rejections == 9
            && self.stale_file_hash_rejections == 9
            && self.recomputed_hash_chain_acceptances == 9
            && self.semantic_rejections_after_hash_acceptance == 9
    }
}

const EXPECTED_MUTATION_IDENTITIES: &[(&str, &str)] = &[
    (
        "root",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/selected_root",
    ),
    (
        "orientation_sign",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
    ),
    (
        "exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
    ),
    (
        "alpha_u",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/alpha/u",
    ),
    (
        "norm",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/norm",
    ),
    (
        "rational_factor_exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/rational_factors/0/exponent",
    ),
    (
        "lp_smaller_root",
        "/role1_interface_addendum/role1_control_fixtures/LP-BOUNDARY/endpoint_fixtures/1/smaller_root",
    ),
    (
        "hnf_matrix_entry",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/ideal_power_column_matrix_rows/0/1",
    ),
    (
        "sint_sign_byte",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate_sint_hex",
    ),
];

#[derive(Debug, PartialEq, Eq)]
struct FrozenMutationResults {
    root: bool,
    orientation_sign: bool,
    exponent: bool,
    alpha_u: bool,
    norm: bool,
    rational_factor_exponent: bool,
    lp_smaller_root: bool,
    hnf_matrix_entry: bool,
    sint_sign_byte: bool,
    evidence: MutationReplayEvidence,
}

impl FrozenMutationResults {
    fn all_rejected(&self) -> bool {
        self.root
            && self.orientation_sign
            && self.exponent
            && self.alpha_u
            && self.norm
            && self.rational_factor_exponent
            && self.lp_smaller_root
            && self.hnf_matrix_entry
            && self.sint_sign_byte
            && self.evidence.complete()
    }
}

fn d23_mutation_carrier() -> Value {
    json!({
        "alpha": {"u":"1","v":"1"},
        "ideal_power_column_matrix_rows": [[8,-7],[0,1]],
        "norm": "8",
        "rational_factors": [{"ell":"2","exponent":3}],
        "selected_root": 1,
        "signed_coordinate": "-3",
        "signed_coordinate_sint_hex": "020000000103",
    })
}

fn lp_mutation_carrier(bound: u64) -> Value {
    json!({
        "larger_root": 2_025_854_371u64,
        "multiplicity": 1,
        "q": bound,
        "sign": 1,
        "smaller_root": 1_109_020_142u64,
    })
}

fn mutation_case(
    id: &'static str,
    fixture_kind: &'static str,
    target_field_path: &'static str,
    carrier_field_path: &'static str,
    before: Value,
    after: Value,
    baseline: &Value,
) -> Result<MutationCase> {
    if baseline.pointer(carrier_field_path) != Some(&before) {
        return Err(format!(
            "mutation {id} before-value differs from its frozen carrier"
        ));
    }
    let mut mutated = baseline.clone();
    let target = mutated
        .pointer_mut(carrier_field_path)
        .ok_or_else(|| format!("mutation {id} carrier pointer is absent"))?;
    *target = after.clone();
    Ok(MutationCase {
        id,
        fixture_kind,
        target_field_path,
        carrier_field_path,
        before,
        after,
        baseline: baseline.clone(),
        mutated,
    })
}

fn mutation_cases(d23: &Value, lp: &Value) -> Result<Vec<MutationCase>> {
    Ok(vec![
        mutation_case(
            "root",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/selected_root",
            "/selected_root",
            json!(1),
            json!(0),
            d23,
        )?,
        mutation_case(
            "orientation_sign",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
            "/signed_coordinate",
            json!("-3"),
            json!("3"),
            d23,
        )?,
        mutation_case(
            "exponent",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
            "/signed_coordinate",
            json!("-3"),
            json!("-2"),
            d23,
        )?,
        mutation_case(
            "alpha_u",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/alpha/u",
            "/alpha/u",
            json!("1"),
            json!("2"),
            d23,
        )?,
        mutation_case(
            "norm",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/norm",
            "/norm",
            json!("8"),
            json!("4"),
            d23,
        )?,
        mutation_case(
            "rational_factor_exponent",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/rational_factors/0/exponent",
            "/rational_factors/0/exponent",
            json!(3),
            json!(2),
            d23,
        )?,
        mutation_case(
            "lp_smaller_root",
            "LP-BOUNDARY",
            "/role1_interface_addendum/role1_control_fixtures/LP-BOUNDARY/endpoint_fixtures/1/smaller_root",
            "/smaller_root",
            json!(1_109_020_142u64),
            json!(1_109_020_143u64),
            lp,
        )?,
        mutation_case(
            "hnf_matrix_entry",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/ideal_power_column_matrix_rows/0/1",
            "/ideal_power_column_matrix_rows/0/1",
            json!(-7),
            json!(-6),
            d23,
        )?,
        mutation_case(
            "sint_sign_byte",
            "D23-ORDER3",
            "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate_sint_hex",
            "/signed_coordinate_sint_hex",
            json!("020000000103"),
            json!("010000000103"),
            d23,
        )?,
    ])
}

fn collect_scalar_bytes(
    value: &Value,
    path: &str,
    output: &mut BTreeMap<String, Vec<u8>>,
) -> Result<()> {
    match value {
        Value::Object(fields) => {
            for (name, child) in fields {
                let escaped = name.replace('~', "~0").replace('/', "~1");
                collect_scalar_bytes(child, &format!("{path}/{escaped}"), output)?;
            }
        }
        Value::Array(items) => {
            for (index, child) in items.iter().enumerate() {
                collect_scalar_bytes(child, &format!("{path}/{index}"), output)?;
            }
        }
        _ => {
            if output
                .insert(path.to_owned(), canonical_json_bytes(value)?)
                .is_some()
            {
                return Err("duplicate scalar path in mutation carrier".to_owned());
            }
        }
    }
    Ok(())
}

fn exact_single_scalar_change(case: &MutationCase) -> Result<bool> {
    let mut baseline = BTreeMap::new();
    let mut mutated = BTreeMap::new();
    collect_scalar_bytes(&case.baseline, "", &mut baseline)?;
    collect_scalar_bytes(&case.mutated, "", &mut mutated)?;
    if baseline.keys().collect::<Vec<_>>() != mutated.keys().collect::<Vec<_>>() {
        return Ok(false);
    }
    let changed = baseline
        .iter()
        .filter(|(path, bytes)| mutated.get(*path) != Some(*bytes))
        .map(|(path, bytes)| (path.as_str(), bytes.as_slice(), mutated[path].as_slice()))
        .collect::<Vec<_>>();
    let before = canonical_json_bytes(&case.before)?;
    let after = canonical_json_bytes(&case.after)?;
    Ok(changed == vec![(case.carrier_field_path, before.as_slice(), after.as_slice())])
}

fn exact_object_keys(value: &Value, expected: &[&str]) -> bool {
    let Some(object) = value.as_object() else {
        return false;
    };
    let actual = object.keys().map(String::as_str).collect::<BTreeSet<_>>();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    actual == expected
}

fn mutation_record_digest(record: &Value) -> Result<String> {
    let projection = json!({
        "fixture": record.get("fixture").cloned().ok_or_else(|| "mutation record fixture absent".to_owned())?,
        "fixture_kind": record.get("fixture_kind").cloned().ok_or_else(|| "mutation record fixture_kind absent".to_owned())?,
        "mutation_id": record.get("mutation_id").cloned().ok_or_else(|| "mutation record mutation_id absent".to_owned())?,
        "target_field_path": record.get("target_field_path").cloned().ok_or_else(|| "mutation record target_field_path absent".to_owned())?,
    });
    let mut bytes = MUTATION_RECORD_DOMAIN.to_vec();
    bytes.extend_from_slice(&canonical_json_bytes(&projection)?);
    Ok(sha256_hex(&bytes))
}

fn mutation_file_digest(file: &Value) -> Result<String> {
    let projection = json!({
        "records": file.get("records").cloned().ok_or_else(|| "mutation file records absent".to_owned())?,
        "schema": file.get("schema").cloned().ok_or_else(|| "mutation file schema absent".to_owned())?,
    });
    let mut bytes = MUTATION_FILE_DOMAIN.to_vec();
    bytes.extend_from_slice(&canonical_json_bytes(&projection)?);
    Ok(sha256_hex(&bytes))
}

fn sealed_mutation_file(case: &MutationCase, fixture: Value) -> Result<Vec<u8>> {
    let mut record = json!({
        "fixture": fixture,
        "fixture_kind": case.fixture_kind,
        "mutation_id": case.id,
        "record_sha256": "",
        "target_field_path": case.target_field_path,
    });
    record["record_sha256"] = Value::String(mutation_record_digest(&record)?);
    let mut file = json!({
        "file_sha256": "",
        "records": [record],
        "schema": MUTATION_FILE_SCHEMA,
    });
    file["file_sha256"] = Value::String(mutation_file_digest(&file)?);
    canonical_json_bytes(&file)
}

fn mutation_hash_gate(bytes: &[u8]) -> std::result::Result<Value, MutationHashFailure> {
    let file: Value = serde_json::from_slice(bytes).map_err(|_| MutationHashFailure::Envelope)?;
    if canonical_json_bytes(&file).map_err(|_| MutationHashFailure::Envelope)? != bytes
        || !exact_object_keys(&file, &["schema", "records", "file_sha256"])
        || file.get("schema").and_then(Value::as_str) != Some(MUTATION_FILE_SCHEMA)
    {
        return Err(MutationHashFailure::Envelope);
    }
    let records = file
        .get("records")
        .and_then(Value::as_array)
        .ok_or(MutationHashFailure::Envelope)?;
    if records.len() != 1
        || !exact_object_keys(
            &records[0],
            &[
                "mutation_id",
                "fixture_kind",
                "target_field_path",
                "fixture",
                "record_sha256",
            ],
        )
    {
        return Err(MutationHashFailure::Envelope);
    }
    let recorded_record_digest = records[0]
        .get("record_sha256")
        .and_then(Value::as_str)
        .ok_or(MutationHashFailure::Envelope)?;
    let computed_record_digest =
        mutation_record_digest(&records[0]).map_err(|_| MutationHashFailure::Envelope)?;
    if recorded_record_digest != computed_record_digest {
        return Err(MutationHashFailure::RecordDigest);
    }
    let recorded_file_digest = file
        .get("file_sha256")
        .and_then(Value::as_str)
        .ok_or(MutationHashFailure::Envelope)?;
    let computed_file_digest =
        mutation_file_digest(&file).map_err(|_| MutationHashFailure::Envelope)?;
    if recorded_file_digest != computed_file_digest {
        return Err(MutationHashFailure::FileDigest);
    }
    Ok(file)
}

fn required_string<'a>(fixture: &'a Value, pointer: &str) -> Result<&'a str> {
    fixture
        .pointer(pointer)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("mutation carrier {pointer} is not a string"))
}

fn required_u64(fixture: &Value, pointer: &str) -> Result<u64> {
    fixture
        .pointer(pointer)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("mutation carrier {pointer} is not a u64"))
}

fn required_i64(fixture: &Value, pointer: &str) -> Result<i64> {
    fixture
        .pointer(pointer)
        .and_then(Value::as_i64)
        .ok_or_else(|| format!("mutation carrier {pointer} is not an i64"))
}

fn lower_hex_bytes(text: &str) -> Result<Vec<u8>> {
    fn nibble(byte: u8) -> Option<u8> {
        match byte {
            b'0'..=b'9' => Some(byte - b'0'),
            b'a'..=b'f' => Some(byte - b'a' + 10),
            _ => None,
        }
    }
    if !text.len().is_multiple_of(2) {
        return Err("mutation carrier hex has odd length".to_owned());
    }
    text.as_bytes()
        .as_chunks::<2>()
        .0
        .iter()
        .map(|pair| {
            let high =
                nibble(pair[0]).ok_or_else(|| "mutation carrier hex is invalid".to_owned())?;
            let low =
                nibble(pair[1]).ok_or_else(|| "mutation carrier hex is invalid".to_owned())?;
            Ok((high << 4) | low)
        })
        .collect()
}

fn d23_from_carrier(fixture: &Value) -> Result<D23Fixture> {
    if !exact_object_keys(
        fixture,
        &[
            "selected_root",
            "signed_coordinate",
            "alpha",
            "norm",
            "rational_factors",
            "ideal_power_column_matrix_rows",
            "signed_coordinate_sint_hex",
        ],
    ) || !exact_object_keys(
        fixture
            .get("alpha")
            .ok_or_else(|| "mutation alpha absent".to_owned())?,
        &["u", "v"],
    ) {
        return Err("D=-23 mutation carrier keys differ from the frozen target".to_owned());
    }
    let factors = fixture
        .get("rational_factors")
        .and_then(Value::as_array)
        .ok_or_else(|| "D=-23 rational_factors is not an array".to_owned())?;
    if factors.len() != 1
        || !exact_object_keys(&factors[0], &["ell", "exponent"])
        || required_string(&factors[0], "/ell")? != "2"
    {
        return Err("D=-23 rational factor carrier is noncanonical".to_owned());
    }
    let matrix = fixture
        .get("ideal_power_column_matrix_rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "D=-23 HNF carrier is not an array".to_owned())?;
    if matrix.len() != 2
        || matrix[0].as_array().map(Vec::len) != Some(2)
        || matrix[1].as_array().map(Vec::len) != Some(2)
        || required_i64(fixture, "/ideal_power_column_matrix_rows/0/0")? != 8
        || required_i64(fixture, "/ideal_power_column_matrix_rows/1/0")? != 0
        || required_i64(fixture, "/ideal_power_column_matrix_rows/1/1")? != 1
    {
        return Err("D=-23 HNF carrier has nonfrozen columns".to_owned());
    }
    let exponent = u32::try_from(required_u64(&factors[0], "/exponent")?)
        .map_err(|_| "D=-23 exponent does not fit u32".to_owned())?;
    Ok(D23Fixture {
        selected_root: required_u64(fixture, "/selected_root")?,
        signed_coordinate: integer(required_string(fixture, "/signed_coordinate")?)?,
        alpha_u: integer(required_string(fixture, "/alpha/u")?)?,
        alpha_v: integer(required_string(fixture, "/alpha/v")?)?,
        claimed_norm: integer(required_string(fixture, "/norm")?)?,
        rational_factor_exponent: exponent,
        ideal_power_top_right: BigInt::from(required_i64(
            fixture,
            "/ideal_power_column_matrix_rows/0/1",
        )?),
        signed_coordinate_encoding: lower_hex_bytes(required_string(
            fixture,
            "/signed_coordinate_sint_hex",
        )?)?,
    })
}

fn lp_from_carrier(fixture: &Value) -> Result<LpFixture> {
    if !exact_object_keys(
        fixture,
        &["q", "smaller_root", "larger_root", "sign", "multiplicity"],
    ) {
        return Err("LP mutation carrier keys differ from the frozen target".to_owned());
    }
    Ok(LpFixture {
        q: required_u64(fixture, "/q")?,
        smaller_root: required_u64(fixture, "/smaller_root")?,
        larger_root: required_u64(fixture, "/larger_root")?,
        sign: u8::try_from(required_u64(fixture, "/sign")?)
            .map_err(|_| "LP sign does not fit u8".to_owned())?,
        multiplicity: u32::try_from(required_u64(fixture, "/multiplicity")?)
            .map_err(|_| "LP multiplicity does not fit u32".to_owned())?,
    })
}

fn semantically_rejected(
    record: &Value,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<bool> {
    let kind = record
        .get("fixture_kind")
        .and_then(Value::as_str)
        .ok_or_else(|| "mutation record fixture_kind is absent".to_owned())?;
    let fixture = record
        .get("fixture")
        .ok_or_else(|| "mutation record fixture is absent".to_owned())?;
    match kind {
        "D23-ORDER3" => Ok(!validate_d23_fixture(&d23_from_carrier(fixture)?)?),
        "LP-BOUNDARY" => Ok(!validate_lp_fixture(
            &lp_from_carrier(fixture)?,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )),
        _ => Err("unknown mutation fixture_kind".to_owned()),
    }
}

fn run_frozen_mutations(
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<FrozenMutationResults> {
    run_frozen_mutations_with_bases(
        bound,
        algebraic_bound,
        trace,
        field_norm,
        discriminant,
        d23_mutation_carrier(),
        lp_mutation_carrier(bound),
    )
}

fn run_frozen_mutations_with_bases(
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
    d23_base: Value,
    lp_base: Value,
) -> Result<FrozenMutationResults> {
    if !validate_d23_fixture(&d23_from_carrier(&d23_base)?)? {
        return Err("frozen D23 mutation base fixture failed replay".to_owned());
    }
    if !validate_lp_fixture(
        &lp_from_carrier(&lp_base)?,
        bound,
        algebraic_bound,
        trace,
        field_norm,
        discriminant,
    ) {
        return Err("frozen LP mutation base fixture failed replay".to_owned());
    }
    let mut evidence = MutationReplayEvidence::default();
    let mut outcomes = BTreeMap::new();
    for case in mutation_cases(&d23_base, &lp_base)? {
        evidence.cases += 1;
        evidence
            .mutation_identities
            .push((case.id, case.target_field_path));
        if !exact_single_scalar_change(&case)? {
            return Err(format!(
                "mutation {} changed more than its frozen scalar or altered held-fixed bytes",
                case.id
            ));
        }
        evidence.exact_single_scalar_changes += 1;

        let baseline_bytes = sealed_mutation_file(&case, case.baseline.clone())?;
        let mut stale_record: Value = serde_json::from_slice(&baseline_bytes)
            .map_err(|error| format!("parse sealed mutation carrier: {error}"))?;
        stale_record["records"][0]["fixture"] = case.mutated.clone();
        let stale_record_bytes = canonical_json_bytes(&stale_record)?;
        if mutation_hash_gate(&stale_record_bytes) != Err(MutationHashFailure::RecordDigest) {
            return Err(format!(
                "mutation {} stale record hash did not fail at the record gate",
                case.id
            ));
        }
        evidence.stale_record_hash_rejections += 1;

        let mut stale_file: Value = serde_json::from_slice(&baseline_bytes)
            .map_err(|error| format!("parse sealed mutation carrier: {error}"))?;
        stale_file["records"][0]["fixture"] = case.mutated.clone();
        let record_digest = mutation_record_digest(&stale_file["records"][0])?;
        stale_file["records"][0]["record_sha256"] = Value::String(record_digest);
        let stale_file_bytes = canonical_json_bytes(&stale_file)?;
        if mutation_hash_gate(&stale_file_bytes) != Err(MutationHashFailure::FileDigest) {
            return Err(format!(
                "mutation {} stale file hash did not fail at the file gate",
                case.id
            ));
        }
        evidence.stale_file_hash_rejections += 1;

        let mutated_bytes = sealed_mutation_file(&case, case.mutated.clone())?;
        let admitted = mutation_hash_gate(&mutated_bytes).map_err(|failure| {
            format!(
                "mutation {} recomputed hash chain was rejected: {failure:?}",
                case.id
            )
        })?;
        evidence.recomputed_hash_chain_acceptances += 1;
        let record = &admitted["records"][0];
        if record["mutation_id"] != case.id
            || record["fixture_kind"] != case.fixture_kind
            || record["target_field_path"] != case.target_field_path
            || record["fixture"] != case.mutated
        {
            return Err(format!(
                "mutation {} changed identity after hashing",
                case.id
            ));
        }
        let rejected = semantically_rejected(
            record,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?;
        if !rejected {
            return Err(format!(
                "mutation {} passed semantic replay after its hash chain passed",
                case.id
            ));
        }
        evidence.semantic_rejections_after_hash_acceptance += 1;
        if outcomes.insert(case.id, rejected).is_some() {
            return Err("duplicate frozen mutation id".to_owned());
        }
    }
    if !evidence.complete() || outcomes.len() != 9 {
        return Err("frozen mutation evidence counters are incomplete".to_owned());
    }
    let passed = |id| outcomes.get(id).copied().unwrap_or(false);
    Ok(FrozenMutationResults {
        root: passed("root"),
        orientation_sign: passed("orientation_sign"),
        exponent: passed("exponent"),
        alpha_u: passed("alpha_u"),
        norm: passed("norm"),
        rational_factor_exponent: passed("rational_factor_exponent"),
        lp_smaller_root: passed("lp_smaller_root"),
        hnf_matrix_entry: passed("hnf_matrix_entry"),
        sint_sign_byte: passed("sint_sign_byte"),
        evidence,
    })
}

const SQRT_D_ACTION_SCALAR_DEC: &str = "6277101735386680763835789423144451611450480845975359606884";
const SQRT_D_MAPPABLE_DEGREE: u64 = 113;
const SQRT_D_LARGE_PRIME_BOUND: u64 = 2_147_483_647;

#[derive(Clone, Debug, PartialEq, Eq)]
struct SqrtDActionFixture {
    alpha_u: BigInt,
    alpha_v: BigInt,
    subgroup_order: BigInt,
    frobenius_scalar: BigInt,
    residual_prime: BigInt,
    ramified_orientations: [(u64, u64); 3],
    residual_degree_map_available: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SqrtDActionReplay {
    frobenius_scalar: BigInt,
    alpha_scalar: BigInt,
    square_congruence: bool,
    canonical_frobenius_root: bool,
    residual_degree_map_unavailable: bool,
    action_admitted: bool,
}

fn frozen_sqrt_d_action_fixture(alpha_u: BigInt, alpha_v: BigInt) -> Result<SqrtDActionFixture> {
    Ok(SqrtDActionFixture {
        alpha_u,
        alpha_v,
        subgroup_order: integer(super::arithmetic::N_DEC)?,
        frobenius_scalar: BigInt::one(),
        residual_prime: integer(super::arithmetic::C_DEC)?,
        ramified_orientations: [(5, 2), (11, 3), (31, 7)],
        residual_degree_map_available: false,
    })
}

fn curve_decimal(curve: &Value, field: &str) -> Result<BigInt> {
    let encoded = curve
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("curve field {field} is not a decimal string"))?;
    integer(encoded)
}

fn curve_unsigned_decimal(curve: &Value, field: &str) -> Result<BigUint> {
    let encoded = curve
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("curve field {field} is not an unsigned decimal string"))?;
    BigUint::parse_bytes(encoded.as_bytes(), 10)
        .ok_or_else(|| format!("curve field {field} is not an unsigned decimal integer"))
}

fn curve_hex(curve: &Value, field: &str) -> Result<BigUint> {
    let encoded = curve
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("curve field {field} is not a hexadecimal string"))?;
    BigUint::parse_bytes(encoded.as_bytes(), 16)
        .ok_or_else(|| format!("curve field {field} is not hexadecimal"))
}

fn replay_sqrt_d_subgroup_action(
    fixture: &SqrtDActionFixture,
    factor_base: &[Entry],
    curve: &Value,
) -> Result<SqrtDActionReplay> {
    let field_norm = integer(super::arithmetic::P_DEC)?;
    let trace = integer(super::arithmetic::T_DEC)?;
    let discriminant = integer(super::arithmetic::D_DEC)?;
    let field_modulus = curve_unsigned_decimal(curve, "p")?;
    let curve_a = curve_unsigned_decimal(curve, "a")?;
    let curve_b = curve_unsigned_decimal(curve, "b")?;
    let generator_x = curve_hex(curve, "generator_x_hex")?;
    let generator_y = curve_hex(curve, "generator_y_hex")?;

    let nonsingular = (BigUint::from(4u8) * curve_a.modpow(&BigUint::from(3u8), &field_modulus)
        + BigUint::from(27u8) * curve_b.modpow(&BigUint::from(2u8), &field_modulus))
        % &field_modulus
        != BigUint::zero();
    let generator_on_curve = generator_x < field_modulus
        && generator_y < field_modulus
        && generator_y.modpow(&BigUint::from(2u8), &field_modulus)
            == (generator_x.modpow(&BigUint::from(3u8), &field_modulus)
                + &curve_a * &generator_x
                + &curve_b)
                % &field_modulus;
    let generator_frobenius_fixed = generator_on_curve
        && generator_x.modpow(&field_modulus, &field_modulus) == generator_x
        && generator_y.modpow(&field_modulus, &field_modulus) == generator_y;

    let curve_obligations = [
        (
            "named curve uid",
            curve.get("curve_uid").and_then(Value::as_str) == Some(CURVE_UID),
        ),
        ("serialized field", curve_decimal(curve, "p")? == field_norm),
        (
            "serialized trace",
            curve_decimal(curve, "trace_t")? == trace,
        ),
        (
            "serialized subgroup order",
            curve_decimal(curve, "n")? == fixture.subgroup_order,
        ),
        (
            "cofactor one",
            curve.get("cofactor").and_then(Value::as_u64) == Some(1),
        ),
        ("nonsingular named curve", nonsingular),
        ("named generator on curve", generator_on_curve),
        (
            "named generator fixed by p-Frobenius",
            generator_frobenius_fixed,
        ),
    ];
    if let Some((name, _)) = curve_obligations.iter().find(|(_, passed)| !passed) {
        return Err(format!("sqrt(D) subgroup curve obligation failed: {name}"));
    }

    let order_identity = &field_norm + BigInt::one() - &trace;
    if fixture.subgroup_order <= BigInt::one() || fixture.subgroup_order != order_identity {
        return Err("sqrt(D) subgroup order is not p+1-t".to_owned());
    }
    if fixture.alpha_u != -&trace || fixture.alpha_v != BigInt::from(2u8) {
        return Err("sqrt(D) coordinates differ from (-t,2)".to_owned());
    }
    let alpha_norm = &fixture.alpha_u * &fixture.alpha_u
        + &trace * &fixture.alpha_u * &fixture.alpha_v
        + &field_norm * &fixture.alpha_v * &fixture.alpha_v;
    if alpha_norm != discriminant.abs() {
        return Err("sqrt(D) alpha norm differs from |D|".to_owned());
    }

    // The characteristic equation alone has two roots modulo n.  The named
    // generator is F_p-rational, so its actual Frobenius action selects 1,
    // not merely whichever polynomial root was supplied by the fixture.
    let characteristic_at_scalar = &fixture.frobenius_scalar * &fixture.frobenius_scalar
        - &trace * &fixture.frobenius_scalar
        + &field_norm;
    let canonical_frobenius_root = fixture.frobenius_scalar == BigInt::one()
        && generator_frobenius_fixed
        && characteristic_at_scalar
            .mod_floor(&fixture.subgroup_order)
            .is_zero();
    if !canonical_frobenius_root {
        return Err("sqrt(D) fixture did not select the generator's Frobenius scalar".to_owned());
    }

    let expected_orientations = [(5u64, 2u64), (11, 3), (31, 7)];
    if fixture.ramified_orientations != expected_orientations {
        return Err("sqrt(D) ramified orientation table is noncanonical".to_owned());
    }
    for (ell, root) in fixture.ramified_orientations {
        let matching = factor_base
            .iter()
            .filter(|entry| entry.ell == ell)
            .collect::<Vec<_>>();
        if matching.len() != 1
            || matching[0].kind != 1
            || matching[0].kronecker != 0
            || matching[0].roots != [root]
            || characteristic_roots(ell, &trace, &field_norm, &discriminant) != [root]
            || (&fixture.alpha_u + &fixture.alpha_v * BigInt::from(root))
                .mod_floor(&BigInt::from(ell))
                != BigInt::zero()
        {
            return Err(format!(
                "sqrt(D) canonical ramified root replay failed at ell={ell}"
            ));
        }
    }

    let alpha_scalar = (&fixture.alpha_u + &fixture.alpha_v * &fixture.frobenius_scalar)
        .mod_floor(&fixture.subgroup_order);
    let expected_scalar = integer(SQRT_D_ACTION_SCALAR_DEC)?;
    let square_congruence = (&alpha_scalar * &alpha_scalar - &discriminant)
        .mod_floor(&fixture.subgroup_order)
        .is_zero();
    if alpha_scalar <= BigInt::zero()
        || alpha_scalar >= fixture.subgroup_order
        || alpha_scalar != expected_scalar
        || !square_congruence
    {
        return Err("sqrt(D) derived subgroup scalar failed its frozen identities".to_owned());
    }

    let large_prime_bound = BigInt::from(SQRT_D_LARGE_PRIME_BOUND);
    let residual_absent = !factor_base
        .iter()
        .any(|entry| BigInt::from(entry.ell) == fixture.residual_prime);
    let residual_degree_map_unavailable = fixture.residual_prime
        == integer(super::arithmetic::C_DEC)?
        && -&discriminant
            == BigInt::from(5u8)
                * BigInt::from(11u8)
                * BigInt::from(31u8)
                * &fixture.residual_prime
        && fixture.residual_prime > BigInt::from(SQRT_D_MAPPABLE_DEGREE)
        && fixture.residual_prime > &large_prime_bound * &large_prime_bound
        && residual_absent
        && !fixture.residual_degree_map_available;
    if !residual_degree_map_unavailable {
        return Err("sqrt(D) residual-degree action unexpectedly became available".to_owned());
    }

    Ok(SqrtDActionReplay {
        frobenius_scalar: fixture.frobenius_scalar.clone(),
        alpha_scalar,
        square_congruence,
        canonical_frobenius_root,
        residual_degree_map_unavailable,
        action_admitted: false,
    })
}

pub fn build(
    provenance: &BuildProvenance,
    factor_base_sha256: &str,
    factor_base_digest: &[u8; 32],
    factor_base: &[Entry],
    reference: &ReferenceOutput,
) -> Result<Value> {
    let d23 = d23_control()?;
    let d23_fixture = d23_fixture();
    let d23_fixture_valid = validate_d23_fixture(&d23_fixture)?;
    let d23_principal = d23.prime_ideal_cubed == d23.principal_alpha;
    let form_a = BigInt::from(2u8);
    let form_b = BigInt::one();
    let form_c = BigInt::from(3u8);
    let form_discriminant = &form_b * &form_b - BigInt::from(4u8) * &form_a * &form_c;
    let canonical_form = form_discriminant == BigInt::from(-23)
        && form_b.abs() <= form_a
        && form_a <= form_c
        && form_a.gcd(&form_b).gcd(&form_c).is_one();
    let d23_ideal_matches = is_characteristic_root(1, 2, &BigInt::one(), &BigInt::from(6u8))
        && Hnf::prime(2, 1)?
            == Hnf {
                a: BigInt::from(2u8),
                b: BigInt::one(),
                d: BigInt::one(),
            };
    let d23_x = BigInt::from(2u8) * &d23_fixture.alpha_u + &d23_fixture.alpha_v;
    let d23_alpha_valid =
        d23_fixture_valid && d23.norm == d23_fixture.claimed_norm && d23_x == BigInt::from(3u8);
    let d23_values = [
        ("canonical_reduced_form_norm_2", canonical_form),
        ("ideal_matches_form", d23_ideal_matches),
        ("alpha_certificate_valid", d23_alpha_valid),
        (
            "encoding_replays",
            encoded_sint(&d23_fixture.signed_coordinate)? == d23_fixture.signed_coordinate_encoding,
        ),
        ("relation_g2_cubed_identity", d23_principal),
    ];

    let sqrt_d = p192_sqrt_discriminant_control()?;
    let (p, t, d) = p192_parameters()?;
    let c = integer(super::arithmetic::C_DEC)?;
    let inverse_two_mod_c = (&c + BigInt::one()) / BigInt::from(2u8);
    let root_c = (&t * inverse_two_mod_c).mod_floor(&c);
    let prime_c = Hnf {
        a: c.clone(),
        b: (-&root_c).mod_floor(&c),
        d: BigInt::one(),
    };
    let mut sqrt_ideal = Hnf::identity();
    for (ell, root) in [(5u64, 2u64), (11, 3), (31, 7)] {
        sqrt_ideal = multiply_hnf(&sqrt_ideal, &Hnf::prime(ell, root)?, &p, &t)?;
    }
    sqrt_ideal = multiply_hnf(&sqrt_ideal, &prime_c, &p, &t)?;
    let sqrt_principal = principal_hnf(&sqrt_d.u, &sqrt_d.v, &p, &t)?;
    let p192_algebra_passes =
        sqrt_ideal == sqrt_principal && sqrt_ideal.determinant() == sqrt_d.norm;
    let action_fixture = frozen_sqrt_d_action_fixture(sqrt_d.u.clone(), sqrt_d.v.clone())?;
    let action_replay =
        replay_sqrt_d_subgroup_action(&action_fixture, factor_base, &curve_object())?;
    let c_degree_unavailable = action_replay.residual_degree_map_unavailable
        && !action_replay.action_admitted
        && sqrt_d.residual_exceeds_bound_squared;
    let scalar_action_rejected = action_replay.frobenius_scalar == BigInt::one()
        && action_replay.alpha_scalar == integer(SQRT_D_ACTION_SCALAR_DEC)?
        && action_replay.square_congruence
        && action_replay.canonical_frobenius_root
        && !action_replay.action_admitted;
    let p192_values = [
        (
            "alpha_coordinates_minus_t_2",
            sqrt_d.u == -integer(super::arithmetic::T_DEC)? && sqrt_d.v == BigInt::from(2u8),
        ),
        (
            "norm_equals_abs_D",
            sqrt_d.norm == integer(super::arithmetic::D_DEC)?.abs(),
        ),
        (
            "factorization_equals_5_11_31_C",
            sqrt_d.factorization.factors == vec![(5, 1), (11, 1), (31, 1)]
                && sqrt_d.factorization.residual.to_string() == super::arithmetic::C_DEC,
        ),
        ("algebra_passes", p192_algebra_passes),
        ("C_degree_map_unavailable", c_degree_unavailable),
        ("C_scalar_action_rejected", scalar_action_rejected),
    ];

    let mut ramified = Vec::new();
    for (ell, root, id) in [
        (5u64, 2u64, "ell_5_scalar"),
        (11, 3, "ell_11_scalar"),
        (31, 7, "ell_31_scalar"),
    ] {
        let square = pow_hnf(&Hnf::prime(ell, root)?, 2, &p, &t)?;
        let scalar = Hnf {
            a: BigInt::from(ell),
            b: BigInt::zero(),
            d: BigInt::from(ell),
        };
        ramified.push((id, square == scalar));
    }

    let inert_absent = !factor_base.iter().any(|entry| entry.ell == 2)
        && characteristic_roots(2, &t, &p, &d).is_empty();
    let inert_residual_prime = 65_537u64;
    let inert_residual_rejected = inert_residual_prime > 65_521
        && primality::is_prime(inert_residual_prime)
        && !residual_prime_admissible(inert_residual_prime, 65_521, &t, &p, &d)
        && characteristic_roots(inert_residual_prime, &t, &p, &d).is_empty();
    let inert_values = [
        ("base_injection_rejected", inert_absent),
        ("residual_injection_rejected", inert_residual_rejected),
    ];

    let large = 2_147_483_647u64;
    let mutations = run_frozen_mutations(large, 65_521, &t, &p, &d)?;
    if !mutations.all_rejected() {
        return Err("one or more frozen Role-1 mutations were accepted".to_owned());
    }
    let mutation_values = [
        ("root_mutation_rejected", mutations.root),
        (
            "orientation_sign_mutation_rejected",
            mutations.orientation_sign && mutations.sint_sign_byte,
        ),
        ("exponent_mutation_rejected", mutations.exponent),
        ("alpha_coordinate_mutation_rejected", mutations.alpha_u),
        (
            "norm_factor_mutation_rejected",
            mutations.norm && mutations.rational_factor_exponent,
        ),
        ("lp_endpoint_mutation_rejected", mutations.lp_smaller_root),
        ("hnf_column_mutation_rejected", mutations.hnf_matrix_entry),
    ];
    let lp_roots = characteristic_roots(large, &t, &p, &d);
    if lp_roots.len() != 2 {
        return Err("LP cancellation fixture is not split".to_owned());
    }
    let lp_positive = Hnf::prime(large, lp_roots[0])?;
    let lp_negative = Hnf::prime(large, lp_roots[1])?;
    let lp_scalar = Hnf {
        a: BigInt::from(large),
        b: BigInt::zero(),
        d: BigInt::from(large),
    };
    let cancellation_values = [
        ("unmatched_orientation_rejected", lp_positive != lp_scalar),
        (
            "exact_conjugate_cancels",
            multiply_hnf(&lp_positive, &lp_negative, &p, &t)? == lp_scalar,
        ),
    ];

    let canonical_order = (0..RECORD_COUNT).collect::<Vec<_>>();
    let rebuilt = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &canonical_order,
    )?;
    let segmented = build_reference_segmented(factor_base, factor_base_digest, 65_521, large)?;
    let disposition_records = split_dispositions(&segmented.disposition_records)?;
    let retained_record_count = disposition_records
        .iter()
        .filter(|record| matches!(record[0], 1..=3))
        .count() as u64;
    let scalar_segmented_values = [
        (
            "candidate_bytes_equal",
            segmented.candidate_records == reference.candidate_records,
        ),
        (
            "candidate_shard_root_equal",
            segmented.candidate_shard_sha256 == reference.candidate_shard_sha256,
        ),
        (
            "candidate_shell_root_equal",
            segmented.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "disposition_bytes_equal",
            segmented.disposition_records == reference.disposition_records,
        ),
        (
            "disposition_shard_root_equal",
            segmented.disposition_shard_sha256 == reference.disposition_shard_sha256,
        ),
        (
            "disposition_shell_root_equal",
            segmented.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "retained_records_equal",
            retained_record_count == reference.status_counts.retained()
                && segmented.certificates == reference.certificates
                && segmented.status_counts == reference.status_counts
                && certificate_store_count(&segmented.certificates)?
                    == reference.status_counts.retained(),
        ),
    ];
    let reverse_order = (0..RECORD_COUNT).rev().collect::<Vec<_>>();
    let reverse = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &reverse_order,
    )?;
    let mut worker_order = Vec::with_capacity(RECORD_COUNT as usize);
    for lane in 0..8u64 {
        worker_order.extend((lane..RECORD_COUNT).step_by(8));
    }
    let workers = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &worker_order,
    )?;
    let direction_values = [
        (
            "forward_reverse_candidate_root_equal",
            reverse.candidate_shard_sha256 == reference.candidate_shard_sha256
                && reverse.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "forward_reverse_disposition_root_equal",
            reverse.disposition_shard_sha256 == reference.disposition_shard_sha256
                && reverse.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "one_eight_worker_candidate_root_equal",
            workers.candidate_shard_sha256 == reference.candidate_shard_sha256
                && workers.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "one_eight_worker_disposition_root_equal",
            workers.disposition_shard_sha256 == reference.disposition_shard_sha256
                && workers.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "independent_in_memory_ref0_regenerations_streams_roots_counters_equal",
            same_reference(&rebuilt, reference),
        ),
    ];

    let large_roots = characteristic_roots(large, &t, &p, &d);
    let large_minus_one = large - 1;
    let large_plus_one = large + 1;
    let large_squared = BigUint::from(large) * large;
    let above_large_squared = &large_squared + BigUint::one();
    let lp_values = [
        (
            "q_L_minus_1_not_prime_rejected_as_lp",
            !endpoint_allowed(large_minus_one, large, 65_521, &t, &p, &d),
        ),
        (
            "q_L_prime_split_accepted_as_one_lp",
            endpoint_allowed(large, large, 65_521, &t, &p, &d)
                && large_roots == vec![1_109_020_142, 2_025_854_371],
        ),
        (
            "q_L_plus_1_out_of_range_rejected_as_lp",
            !endpoint_allowed(large_plus_one, large, 65_521, &t, &p, &d),
        ),
        (
            "residual_L_is_one_large_prime",
            residual_status(&BigUint::from(large), large, 65_521, &t, &p, &d)? == 2,
        ),
        (
            "residual_L_squared_is_two_large_prime_with_multiplicity_two",
            residual_status(&large_squared, large, 65_521, &t, &p, &d)? == 3
                && primality::factor(large * large)? == vec![(large, 2)],
        ),
        (
            "residual_L_squared_plus_one_is_rejected_by_inequality",
            residual_status(&above_large_squared, large, 65_521, &t, &p, &d)? == 4,
        ),
    ];

    let controls = vec![
        control("D23-ORDER3", &d23_values),
        control("P192-SQRT-D", &p192_values),
        control("P192-RAMIFIED-SQUARES", &ramified),
        control("INERT-FACTOR-REJECTION", &inert_values),
        control("MUTATION-REJECTION", &mutation_values),
        control("LP-CANCELLATION", &cancellation_values),
        control("REF0-SCALAR-SEGMENTED", &scalar_segmented_values),
        control("REF0-DIRECTION-WORKERS", &direction_values),
        control("LP-BOUNDARY", &lp_values),
    ];
    if controls.iter().any(|value| value["status"] != "PASS") {
        return Err("one or more frozen Role-1 controls failed".to_owned());
    }
    let mut digest_preimage = b"P192-WCM-CONTROL-RESULTS-v1\0".to_vec();
    digest_preimage.extend_from_slice(&canonical_json_bytes(&Value::Array(controls.clone()))?);
    let control_results_sha256 = sha256_hex(&digest_preimage);
    Ok(json!({
        "control_results_sha256": control_results_sha256,
        "controls": controls,
        "curve_uid": CURVE_UID,
        "experiment_id": EXPERIMENT_ID,
        "factor_base_sha256": factor_base_sha256,
        "overall_status": "PASS",
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": "p192-wcm-controls-v1",
        "source_commit": provenance.source_commit,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disposition_splitter_rejects_wrong_certificate_flag() {
        let mut bytes = Vec::new();
        for index in 0..RECORD_COUNT {
            bytes.push(4);
            bytes.extend_from_slice(&index.to_be_bytes());
            bytes.push(0);
        }
        assert_eq!(
            split_dispositions(&bytes).unwrap().len() as u64,
            RECORD_COUNT
        );
        bytes[9] = 1;
        assert!(split_dispositions(&bytes).is_err());
    }

    #[test]
    fn all_nine_frozen_mutations_replay_through_common_validator() {
        let (field_norm, trace, discriminant) = p192_parameters().unwrap();
        let results =
            run_frozen_mutations(2_147_483_647, 65_521, &trace, &field_norm, &discriminant)
                .unwrap();
        assert_eq!(
            [
                results.root,
                results.orientation_sign,
                results.exponent,
                results.alpha_u,
                results.norm,
                results.rational_factor_exponent,
                results.lp_smaller_root,
                results.hnf_matrix_entry,
                results.sint_sign_byte,
            ],
            [true; 9]
        );
        assert!(results.all_rejected());
        assert_eq!(results.evidence.cases, 9);
        assert_eq!(
            results.evidence.mutation_identities.as_slice(),
            EXPECTED_MUTATION_IDENTITIES
        );
        assert_eq!(results.evidence.exact_single_scalar_changes, 9);
        assert_eq!(results.evidence.stale_record_hash_rejections, 9);
        assert_eq!(results.evidence.stale_file_hash_rejections, 9);
        assert_eq!(results.evidence.recomputed_hash_chain_acceptances, 9);
        assert_eq!(
            results.evidence.semantic_rejections_after_hash_acceptance,
            9
        );
    }

    #[test]
    fn mutation_replay_rejects_invalid_held_fixed_bases() {
        let (field_norm, trace, discriminant) = p192_parameters().unwrap();
        let bound = 2_147_483_647;

        let mut invalid_d23 = d23_mutation_carrier();
        invalid_d23["alpha"]["v"] = json!("2");
        assert_eq!(
            run_frozen_mutations_with_bases(
                bound,
                65_521,
                &trace,
                &field_norm,
                &discriminant,
                invalid_d23,
                lp_mutation_carrier(bound),
            )
            .unwrap_err(),
            "frozen D23 mutation base fixture failed replay"
        );

        let mut invalid_lp = lp_mutation_carrier(bound);
        invalid_lp["larger_root"] = json!(2_025_854_370u64);
        assert_eq!(
            run_frozen_mutations_with_bases(
                bound,
                65_521,
                &trace,
                &field_norm,
                &discriminant,
                d23_mutation_carrier(),
                invalid_lp,
            )
            .unwrap_err(),
            "frozen LP mutation base fixture failed replay"
        );
    }

    #[test]
    fn sqrt_d_subgroup_action_derives_the_frozen_scalar() {
        let factor_base = super::super::factor_base::build(65_521).unwrap();
        let sqrt_d = p192_sqrt_discriminant_control().unwrap();
        let fixture = frozen_sqrt_d_action_fixture(sqrt_d.u.clone(), sqrt_d.v.clone()).unwrap();
        let replay =
            replay_sqrt_d_subgroup_action(&fixture, &factor_base, &curve_object()).unwrap();
        let discriminant = integer(super::super::arithmetic::D_DEC).unwrap();

        assert_eq!(replay.frobenius_scalar, BigInt::one());
        assert_eq!(replay.alpha_scalar.to_string(), SQRT_D_ACTION_SCALAR_DEC);
        assert!(replay.alpha_scalar > BigInt::zero());
        assert!(replay.alpha_scalar < fixture.subgroup_order);
        assert!((&replay.alpha_scalar * &replay.alpha_scalar - discriminant)
            .mod_floor(&fixture.subgroup_order)
            .is_zero());
        assert!(replay.square_congruence);
        assert!(replay.canonical_frobenius_root);
        assert!(replay.residual_degree_map_unavailable);
        assert!(!replay.action_admitted);
    }

    #[test]
    fn sqrt_d_subgroup_action_rejects_noncanonical_and_unavailable_paths() {
        let factor_base = super::super::factor_base::build(65_521).unwrap();
        let sqrt_d = p192_sqrt_discriminant_control().unwrap();
        let fixture = frozen_sqrt_d_action_fixture(sqrt_d.u, sqrt_d.v).unwrap();
        let curve = curve_object();

        let mut wrong_u = fixture.clone();
        wrong_u.alpha_u += BigInt::one();
        assert!(replay_sqrt_d_subgroup_action(&wrong_u, &factor_base, &curve).is_err());

        let mut wrong_order = fixture.clone();
        wrong_order.subgroup_order -= BigInt::one();
        assert!(replay_sqrt_d_subgroup_action(&wrong_order, &factor_base, &curve).is_err());

        // p mod n is the other root of X^2-tX+p modulo n.  It must still be
        // rejected because it is not the action of Frobenius on the named
        // F_p-rational generator.
        let field_norm = integer(super::super::arithmetic::P_DEC).unwrap();
        let trace = integer(super::super::arithmetic::T_DEC).unwrap();
        let alternative_root = field_norm.mod_floor(&fixture.subgroup_order);
        let alternative_characteristic =
            &alternative_root * &alternative_root - &trace * &alternative_root + &field_norm;
        assert_ne!(alternative_root, BigInt::one());
        assert!(alternative_characteristic
            .mod_floor(&fixture.subgroup_order)
            .is_zero());
        let mut noncanonical_root = fixture.clone();
        noncanonical_root.frobenius_scalar = alternative_root;
        assert!(replay_sqrt_d_subgroup_action(&noncanonical_root, &factor_base, &curve).is_err());

        let mut wrong_orientation = fixture.clone();
        wrong_orientation.ramified_orientations[0].1 = 3;
        assert!(replay_sqrt_d_subgroup_action(&wrong_orientation, &factor_base, &curve).is_err());

        let mut mapped_residual = fixture.clone();
        mapped_residual.residual_degree_map_available = true;
        assert!(replay_sqrt_d_subgroup_action(&mapped_residual, &factor_base, &curve).is_err());

        let mut nonrational_generator = curve;
        nonrational_generator["generator_y_hex"] = Value::String("0".to_owned());
        assert!(
            replay_sqrt_d_subgroup_action(&fixture, &factor_base, &nonrational_generator).is_err()
        );
    }

    #[test]
    fn complete_role1_control_set_recomputes() {
        let factor_base = super::super::factor_base::build(65_521).unwrap();
        let digest = [0x5au8; 32];
        let reference =
            super::super::ref0::build_reference(&factor_base, &digest, 65_521, 2_147_483_647)
                .unwrap();
        let provenance = BuildProvenance {
            protocol_commit: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            source_commit: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            source_dirty: false,
            git_metadata_present: true,
            build_profile: "release",
        };
        let document = build(
            &provenance,
            &super::super::ref0::hex(&digest),
            &digest,
            &factor_base,
            &reference,
        )
        .unwrap();
        assert_eq!(document["overall_status"], "PASS");
        assert_eq!(document["controls"].as_array().unwrap().len(), 9);
        assert_eq!(
            document["controls"][7]["assertions"][4]["assertion_id"],
            "independent_in_memory_ref0_regenerations_streams_roots_counters_equal"
        );
    }
}
