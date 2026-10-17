//! Independent finite-field point-count replay for one fixed GF(8) fixture.
//!
//! This module deliberately does not import the production binary-curve,
//! Mumford-divisor, Jacobian, or transfer-map implementations. Its only shared
//! mathematical input is the public field and curve equation recorded in the
//! certificate.

use crypto_lib::hash::sha256;
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use serde_json::{json, Map, Number, Value};
use std::collections::BTreeSet;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

const SCHEMA: &str = "ecdlp.fixed-f8-genus3-order-certificate/v1";
const UID_PREFIX: &str = "urn:ecdlp-order-certificate:1:sha256:";
const MAX_CERTIFICATE_BYTES: u64 = 1 << 20;
const VERIFIER_PATH: &str = "src/bin/ecdlp_audit/fixed_f8_order.rs";
const VERIFIER_SOURCE: &[u8] = include_bytes!("fixed_f8_order.rs");

#[derive(Debug)]
struct StrictValue(Value);

impl<'de> Deserialize<'de> for StrictValue {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_any(StrictValueVisitor)
    }
}

struct StrictValueVisitor;

impl<'de> Visitor<'de> for StrictValueVisitor {
    type Value = StrictValue;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("strict JSON containing only finite integer numbers")
    }

    fn visit_bool<E>(self, value: bool) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::Bool(value)))
    }

    fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::Number(Number::from(value))))
    }

    fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::Number(Number::from(value))))
    }

    fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Err(E::custom(format!(
            "non-integer or non-standard JSON number: {value}"
        )))
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(StrictValue(Value::String(value.to_owned())))
    }

    fn visit_string<E>(self, value: String) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::String(value)))
    }

    fn visit_none<E>(self) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::Null))
    }

    fn visit_unit<E>(self) -> Result<Self::Value, E> {
        Ok(StrictValue(Value::Null))
    }

    fn visit_seq<A>(self, mut sequence: A) -> Result<Self::Value, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let mut values = Vec::new();
        while let Some(StrictValue(value)) = sequence.next_element::<StrictValue>()? {
            values.push(value);
        }
        Ok(StrictValue(Value::Array(values)))
    }

    fn visit_map<A>(self, mut entries: A) -> Result<Self::Value, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut keys = BTreeSet::new();
        let mut object = Map::new();
        while let Some(key) = entries.next_key::<String>()? {
            if !keys.insert(key.clone()) {
                return Err(de::Error::custom(format!("duplicate JSON key: {key}")));
            }
            let StrictValue(value) = entries.next_value::<StrictValue>()?;
            object.insert(key, value);
        }
        Ok(StrictValue(Value::Object(object)))
    }
}

fn require(condition: bool, message: impl Into<String>) -> Result<(), String> {
    if condition {
        Ok(())
    } else {
        Err(message.into())
    }
}

fn parse_strict_json(bytes: &[u8]) -> Result<Value, String> {
    let mut deserializer = serde_json::Deserializer::from_slice(bytes);
    let StrictValue(value) = StrictValue::deserialize(&mut deserializer)
        .map_err(|error| format!("cannot parse certificate: {error}"))?;
    deserializer
        .end()
        .map_err(|error| format!("cannot parse certificate: {error}"))?;
    Ok(value)
}

fn canonical_json(value: &Value) -> Result<String, String> {
    fn write(value: &Value, output: &mut String) -> Result<(), String> {
        match value {
            Value::Null => output.push_str("null"),
            Value::Bool(flag) => output.push_str(if *flag { "true" } else { "false" }),
            Value::Number(number) => {
                require(
                    !number.is_f64(),
                    "floats are forbidden in certificate identities",
                )?;
                output.push_str(&number.to_string());
            }
            Value::String(text) => output.push_str(
                &serde_json::to_string(text)
                    .map_err(|error| format!("cannot encode certificate string: {error}"))?,
            ),
            Value::Array(items) => {
                output.push('[');
                for (index, item) in items.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    write(item, output)?;
                }
                output.push(']');
            }
            Value::Object(object) => {
                output.push('{');
                let mut keys: Vec<&String> = object.keys().collect();
                keys.sort();
                for (index, key) in keys.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    output.push_str(
                        &serde_json::to_string(key)
                            .map_err(|error| format!("cannot encode certificate key: {error}"))?,
                    );
                    output.push(':');
                    write(&object[*key], output)?;
                }
                output.push('}');
            }
        }
        Ok(())
    }

    let mut output = String::new();
    write(value, &mut output)?;
    Ok(output)
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn expected_record(verifier_sha256: &str) -> Value {
    json!({
        "claim": {
            "statement": "For the stated smooth genus-3 curve H over GF(8), #Jac(H)(GF(8)) = 486.",
            "status": "supported",
            "scope": "exact fixed polynomial-basis GF(8) source curve only",
            "evidence_kind": "implementation-separated exhaustive extension point count and zeta reconstruction"
        },
        "fixture": {
            "fixture_id": "binary-cover-fixed-f8-a1-b3-d7",
            "base_field": {
                "characteristic": 2,
                "degree": 3,
                "cardinality": 8,
                "binary_modulus": "0xb",
                "polynomial": "z^3 + z + 1",
                "element_encoding": "polynomial-basis bits, least-significant bit is the constant coefficient"
            },
            "construction_parameters": {"a": "0x1", "b": "0x3", "d": "0x7"},
            "source_curve": {
                "model": "v^2 + h(u)*v = f(u)",
                "equation": "v^2 + u^2*v = u^7 + u^4 + (z^2+z+1)*u",
                "h_coefficients_ascending_hex": ["0x0", "0x0", "0x1"],
                "f_coefficients_ascending_hex": ["0x0", "0x7", "0x0", "0x0", "0x1", "0x0", "0x0", "0x1"],
                "genus": 3,
                "infinity_configuration": "one rational point"
            },
            "cover_context": {
                "target_equation": "y^2 + x*y = x^3 + x^2 + 3",
                "map": "(u,v) -> (u^3,u*v+7)",
                "degree": 3,
                "certified_by_this_record": false
            }
        },
        "geometric_preconditions": {
            "smooth": true,
            "genus": 3,
            "rational_points_at_infinity": 1,
            "checks": [
                "d^2=b and d is nonzero in GF(8)",
                "the affine v-derivative is u^2; at u=0 the u-derivative is d, so there is no affine singularity",
                "in the infinity chart t=1/u and w=v/u^4, w^2+t^2*w=t+t^4+d*t^7 is smooth at t=w=0",
                "s=v/u^2 gives s^2+s=u^3+1+d/u^3 with two order-3 poles; the Artin-Schreier genus formula gives genus 3"
            ]
        },
        "method": {
            "point_count": "For n=1,2,3, enumerate every (u,v) in GF(8^n)^2 satisfying the exact source equation and add its unique rational point at infinity.",
            "extension_construction": "Use independently checked irreducible binary polynomials and an explicitly checked root of z^3+z+1 in each extension.",
            "zeta_reconstruction": "Set S_n=8^n+1-N_n, apply Newton identities through n=3, then apply the genus-3 functional equation c[6-i]=8^(3-i)c[i].",
            "jacobian_order": "For a smooth projective curve over GF(8), #Jac(H)(GF(8))=P(1), where P is the zeta numerator.",
            "implementation_separation": {
                "language": "Rust 1.90.0 native verifier",
                "production_curve_arithmetic_imported": false,
                "mumford_arithmetic_used": false,
                "production_486_class_table_used": false,
                "shared_input": "the public field and curve equations"
            }
        },
        "results": {
            "extension_point_counts": [
                {"extension_degree_over_f8": 1, "binary_degree": 3, "field_cardinality": 8, "binary_modulus": "0xb", "base_generator_embedding": "0x2", "d_embedding": "0x7", "points": 6},
                {"extension_degree_over_f8": 2, "binary_degree": 6, "field_cardinality": 64, "binary_modulus": "0x43", "base_generator_embedding": "0xe", "d_embedding": "0x18", "points": 104},
                {"extension_degree_over_f8": 3, "binary_degree": 9, "field_cardinality": 512, "binary_modulus": "0x203", "base_generator_embedding": "0xfc", "d_embedding": "0x1d3", "points": 558}
            ],
            "frobenius_power_sums": [3, -39, -45],
            "zeta_numerator": {
                "definition": "Z(H,T)=P(T)/((1-T)*(1-8*T))",
                "variable": "T",
                "coefficient_order": "ascending",
                "coefficients": [1, -3, 24, -48, 192, -192, 512],
                "polynomial": "1 - 3*T + 24*T^2 - 48*T^3 + 192*T^4 - 192*T^5 + 512*T^6",
                "functional_equation": "c[6-i] = 8^(3-i) * c[i] for 0 <= i <= 3"
            },
            "jacobian_order_over_f8": {"derivation": "P(1)", "value": 486},
            "derived_extension_point_counts_n1_to_n6": [6, 104, 558, 3920, 32646, 263192]
        },
        "provenance": {
            "created_utc": "2026-10-14",
            "repository": "aburan28/crypto",
            "repository_base_commit": "1e447956d3616092a223264f1c6c31eff8ac9979",
            "verifier": {
                "path": VERIFIER_PATH,
                "sha256": verifier_sha256,
                "tested_with": "rustc 1.90.0",
                "imports_production_arithmetic": false
            },
            "canonical_uid_rule": "SHA-256 of UTF-8 JSON(record) with recursively sorted keys, no insignificant whitespace, ensure_ascii=false"
        },
        "limitations": [
            "This certifies only the cardinality 486 of the stated Jacobian over GF(8).",
            "It does not certify the Jacobian group structure, the A/B presentation, divisor reduction, addition, scalar multiplication, or enumeration of reduced representatives.",
            "It does not certify the cover map, arbitrary-class pushforward, norm homomorphism, pull-push composition, kernel fibers, or subgroup preservation.",
            "It does not generalize to another field representation, source curve, cover, or catalog row.",
            "The replay is a separate implementation but has not received an external CAS replay or independent reviewer signature."
        ]
    })
}

fn polynomial_degree(value: u64) -> i32 {
    if value == 0 {
        -1
    } else {
        63 - value.leading_zeros() as i32
    }
}

fn polynomial_mod(mut value: u64, modulus: u64) -> Result<u64, String> {
    require(modulus != 0, "polynomial modulus must be nonzero")?;
    let modulus_degree = polynomial_degree(modulus);
    if modulus_degree == 0 {
        return Ok(0);
    }
    while polynomial_degree(value) >= modulus_degree {
        value ^= modulus << (polynomial_degree(value) - modulus_degree);
    }
    Ok(value)
}

fn polynomial_gcd(mut left: u64, mut right: u64) -> Result<u64, String> {
    while right != 0 {
        (left, right) = (right, polynomial_mod(left, right)?);
    }
    Ok(left)
}

fn field_mul(mut left: u64, mut right: u64, modulus: u64) -> Result<u64, String> {
    let mut product = 0;
    while right != 0 {
        if right & 1 != 0 {
            product ^= left;
        }
        right >>= 1;
        left <<= 1;
    }
    polynomial_mod(product, modulus)
}

fn field_square(value: u64, modulus: u64) -> Result<u64, String> {
    field_mul(value, value, modulus)
}

fn field_pow(mut value: u64, mut exponent: u64, modulus: u64) -> Result<u64, String> {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = field_mul(result, value, modulus)?;
        }
        value = field_square(value, modulus)?;
        exponent >>= 1;
    }
    Ok(result)
}

fn is_irreducible_binary(modulus: u64, degree: u32) -> Result<bool, String> {
    if polynomial_degree(modulus) != degree as i32 || modulus & 1 == 0 {
        return Ok(false);
    }
    let x = 0b10;
    let mut power = x;
    for _ in 1..=degree / 2 {
        power = field_square(power, modulus)?;
        if polynomial_gcd(power ^ x, modulus)? != 1 {
            return Ok(false);
        }
    }
    power = x;
    for _ in 0..degree {
        power = field_square(power, modulus)?;
    }
    Ok(power == x)
}

fn integer(record: &Value, key: &str) -> Result<u64, String> {
    record
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("{key} must be a nonnegative integer"))
}

fn parse_hex(record: &Value, key: &str) -> Result<u64, String> {
    let encoded = record
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("{key} must be hex"))?;
    require(encoded.starts_with("0x"), format!("{key} must be hex"))?;
    u64::from_str_radix(&encoded[2..], 16).map_err(|_| format!("{key} is not valid hex"))
}

fn count_curve_points(field_record: &Value) -> Result<u64, String> {
    let extension_degree = integer(field_record, "extension_degree_over_f8")?;
    require(
        (1..=3).contains(&extension_degree),
        "invalid extension degree",
    )?;
    let binary_degree = integer(field_record, "binary_degree")?;
    require(
        binary_degree == 3 * extension_degree,
        "binary degree mismatch",
    )?;
    let field_cardinality = integer(field_record, "field_cardinality")?;
    require(
        field_cardinality == 1u64 << binary_degree,
        "field cardinality mismatch",
    )?;

    let modulus = parse_hex(field_record, "binary_modulus")?;
    let z = parse_hex(field_record, "base_generator_embedding")?;
    let declared_d = parse_hex(field_record, "d_embedding")?;
    require(
        is_irreducible_binary(modulus, binary_degree as u32)?,
        "reducible field modulus",
    )?;
    require(
        z < field_cardinality,
        "base generator embedding is outside the field",
    )?;
    require(
        field_pow(z, 3, modulus)? ^ z ^ 1 == 0,
        "base generator embedding does not satisfy z^3+z+1",
    )?;
    let d = field_square(z, modulus)? ^ z ^ 1;
    require(d == declared_d, "d embedding mismatch")?;
    require(field_square(d, modulus)? == z ^ 1, "d^2 != 3")?;

    let mut points = 1;
    for u in 0..field_cardinality {
        let u2 = field_square(u, modulus)?;
        let u4 = field_square(u2, modulus)?;
        let u7 = field_mul(field_mul(u4, u2, modulus)?, u, modulus)?;
        let right = u7 ^ u4 ^ field_mul(d, u, modulus)?;
        for v in 0..field_cardinality {
            let left = field_square(v, modulus)? ^ field_mul(u2, v, modulus)?;
            if left == right {
                points += 1;
            }
        }
    }
    Ok(points)
}

fn reconstructed_zeta_coefficients(q: i64, counts: &[i64]) -> Result<(Vec<i64>, Vec<i64>), String> {
    require(
        counts.len() == 3,
        "a genus-3 replay requires N1, N2, and N3",
    )?;
    let power_sums: Vec<i64> = counts
        .iter()
        .enumerate()
        .map(|(index, count)| q.pow((index + 1) as u32) + 1 - count)
        .collect();
    let (s1, s2, s3) = (power_sums[0], power_sums[1], power_sums[2]);
    let e1 = s1;
    require((e1 * s1 - s2) % 2 == 0, "nonintegral second elementary sum")?;
    let e2 = (e1 * s1 - s2) / 2;
    require(
        (s3 - e1 * s2 + e2 * s1) % 3 == 0,
        "nonintegral third elementary sum",
    )?;
    let e3 = (s3 - e1 * s2 + e2 * s1) / 3;
    let coefficients = vec![1, -e1, e2, -e3, q * e2, -(q * q) * e1, q.pow(3)];
    Ok((power_sums, coefficients))
}

fn derived_counts(q: i64, coefficients: &[i64]) -> Vec<i64> {
    let mut power_sums = Vec::new();
    for n in 1..=6 {
        let subtotal: i64 = (1..n)
            .map(|j| coefficients[j] * power_sums[n - j - 1])
            .sum();
        power_sums.push(-(subtotal + n as i64 * coefficients[n]));
    }
    (1..=6)
        .map(|n| q.pow(n as u32) + 1 - power_sums[n - 1])
        .collect()
}

fn replay_bytes(bytes: &[u8]) -> Result<Value, String> {
    require(
        bytes.len() as u64 <= MAX_CERTIFICATE_BYTES,
        "certificate exceeds the 1 MiB replay limit",
    )?;
    let certificate = parse_strict_json(bytes)?;
    let object = certificate
        .as_object()
        .ok_or_else(|| "certificate must be a JSON object".to_owned())?;
    let actual_keys: BTreeSet<&str> = object.keys().map(String::as_str).collect();
    let expected_keys = BTreeSet::from(["schema", "certificate_uid", "record"]);
    require(actual_keys == expected_keys, "unexpected top-level keys")?;
    require(
        certificate["schema"] == SCHEMA,
        "unsupported certificate schema",
    )?;
    let record = certificate
        .get("record")
        .ok_or_else(|| "record is missing".to_owned())?;
    require(record.is_object(), "record must be an object")?;
    let expected_uid = format!(
        "{UID_PREFIX}{}",
        sha256_hex(canonical_json(record)?.as_bytes())
    );
    require(
        certificate["certificate_uid"] == expected_uid,
        "certificate UID mismatch",
    )?;

    let verifier_sha256 = sha256_hex(VERIFIER_SOURCE);
    require(
        record == &expected_record(&verifier_sha256),
        "claim-bearing certificate record mismatch",
    )?;

    let results = &record["results"];
    let extension_counts = results["extension_point_counts"]
        .as_array()
        .ok_or_else(|| "extension point counts must be an array".to_owned())?;
    let mut replayed_counts = Vec::new();
    for (index, field_record) in extension_counts.iter().enumerate() {
        let expected_degree = index as u64 + 1;
        require(
            integer(field_record, "extension_degree_over_f8")? == expected_degree,
            "point-count records must be ordered N1, N2, N3",
        )?;
        let actual_count = count_curve_points(field_record)?;
        require(
            integer(field_record, "points")? == actual_count,
            format!("N{expected_degree} mismatch"),
        )?;
        replayed_counts.push(actual_count as i64);
    }

    let (power_sums, zeta_coefficients) = reconstructed_zeta_coefficients(8, &replayed_counts)?;
    require(
        results["frobenius_power_sums"] == json!(power_sums),
        "Frobenius power sums mismatch",
    )?;
    require(
        results["zeta_numerator"]["coefficients"] == json!(zeta_coefficients),
        "zeta coefficients mismatch",
    )?;
    let jacobian_order: i64 = zeta_coefficients.iter().sum();
    require(
        results["jacobian_order_over_f8"] == json!({"derivation": "P(1)", "value": jacobian_order}),
        "Jacobian-order derivation mismatch",
    )?;
    let replayed_derived_counts = derived_counts(8, &zeta_coefficients);
    require(
        results["derived_extension_point_counts_n1_to_n6"] == json!(replayed_derived_counts),
        "derived counts mismatch",
    )?;

    for (index, count) in replayed_counts.iter().enumerate() {
        let n = index as u32 + 1;
        let delta = count - (8i64.pow(n) + 1);
        require(
            delta * delta <= 36 * 8i64.pow(n),
            format!("N{n} violates the genus-3 Hasse-Weil bound"),
        )?;
    }
    require(
        jacobian_order == 486,
        "unexpected fixed-fixture Jacobian order",
    )?;

    Ok(json!({
        "schema": SCHEMA,
        "certificate_uid": expected_uid,
        "extension_point_counts": replayed_counts,
        "zeta_numerator_coefficients_ascending": zeta_coefficients,
        "jacobian_order_over_f8": jacobian_order,
        "status": "verified"
    }))
}

pub(super) fn verify_certificate_file(input: &Path, repo_root: &Path) -> Result<Value, String> {
    let root = fs::canonicalize(repo_root)
        .map_err(|error| format!("cannot resolve repository root: {error}"))?;
    let unresolved: PathBuf = if input.is_absolute() {
        input.to_owned()
    } else {
        root.join(input)
    };
    let path = fs::canonicalize(&unresolved).map_err(|error| {
        format!(
            "cannot resolve certificate {}: {error}",
            unresolved.display()
        )
    })?;
    require(
        path.starts_with(&root),
        format!(
            "certificate resolves outside repository: {}",
            path.display()
        ),
    )?;
    require(
        path.is_file(),
        format!("certificate is not a file: {}", path.display()),
    )?;
    let metadata = fs::metadata(&path)
        .map_err(|error| format!("cannot stat certificate {}: {error}", path.display()))?;
    require(
        metadata.len() <= MAX_CERTIFICATE_BYTES,
        "certificate exceeds the 1 MiB replay limit",
    )?;
    let bytes = fs::read(&path)
        .map_err(|error| format!("cannot read certificate {}: {error}", path.display()))?;
    replay_bytes(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CERTIFICATE: &[u8] =
        include_bytes!("../../../docs/curves/fixed-f8-genus3-order-certificate.json");

    fn fixture() -> Value {
        parse_strict_json(CERTIFICATE).expect("checked-in certificate must parse")
    }

    fn refresh_uid(certificate: &mut Value) {
        let canonical = canonical_json(&certificate["record"]).expect("record must canonicalize");
        certificate["certificate_uid"] =
            Value::String(format!("{UID_PREFIX}{}", sha256_hex(canonical.as_bytes())));
    }

    fn assert_mutation_rejected(mutate: impl FnOnce(&mut Value), refresh: bool) {
        let mut certificate = fixture();
        mutate(&mut certificate);
        if refresh {
            refresh_uid(&mut certificate);
        }
        let encoded = serde_json::to_vec(&certificate).expect("mutant must encode");
        assert!(replay_bytes(&encoded).is_err());
    }

    #[test]
    fn checked_in_certificate_replays() {
        let result = replay_bytes(CERTIFICATE).expect("certificate must replay");
        assert_eq!(result["extension_point_counts"], json!([6, 104, 558]));
        assert_eq!(result["jacobian_order_over_f8"], 486);
    }

    #[test]
    fn scope_expansion_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| {
                value["record"]["claim"]["statement"] = "All catalog transfers are correct.".into();
                value["record"]["claim"]["scope"] =
                    "all fields, curves, arithmetic, and transfer maps".into();
            },
            true,
        );
    }

    #[test]
    fn cover_certification_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| {
                value["record"]["fixture"]["cover_context"]["certified_by_this_record"] =
                    true.into()
            },
            true,
        );
    }

    #[test]
    fn displayed_equation_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| value["record"]["fixture"]["source_curve"]["equation"] = "v^2 = u".into(),
            true,
        );
    }

    #[test]
    fn displayed_zeta_polynomial_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| value["record"]["results"]["zeta_numerator"]["polynomial"] = "1".into(),
            true,
        );
    }

    #[test]
    fn removed_limitations_are_rejected_after_uid_refresh() {
        assert_mutation_rejected(|value| value["record"]["limitations"] = json!([]), true);
    }

    #[test]
    fn extra_nested_key_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| value["record"]["method"]["catalog_transfer_certified"] = true.into(),
            true,
        );
    }

    #[test]
    fn bad_embedding_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| {
                value["record"]["results"]["extension_point_counts"][1]
                    ["base_generator_embedding"] = "0xf".into()
            },
            true,
        );
    }

    #[test]
    fn bad_count_is_rejected_after_uid_refresh() {
        assert_mutation_rejected(
            |value| value["record"]["results"]["extension_point_counts"][2]["points"] = 559.into(),
            true,
        );
    }

    #[test]
    fn uid_mismatch_is_rejected() {
        assert_mutation_rejected(
            |value| value["record"]["claim"]["status"] = "unknown".into(),
            false,
        );
    }

    #[test]
    fn duplicate_key_is_rejected() {
        let original = std::str::from_utf8(CERTIFICATE).expect("certificate is UTF-8");
        let needle = format!("\"schema\": \"{SCHEMA}\",");
        let replacement = format!("{needle}\n  {needle}");
        let mutant = original.replacen(&needle, &replacement, 1);
        assert!(replay_bytes(mutant.as_bytes()).is_err());
    }

    #[test]
    fn nan_is_rejected() {
        let original = std::str::from_utf8(CERTIFICATE).expect("certificate is UTF-8");
        let mutant = original.replacen("\"status\": \"supported\"", "\"status\": NaN", 1);
        assert!(replay_bytes(mutant.as_bytes()).is_err());
    }
}
