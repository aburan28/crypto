//! Data-only producer schema helpers.
//!
//! The additive experiment contract is intentionally localized here.  Math
//! modules return exact integers and byte strings; changes to artifact field
//! names or packaging should not alter their arithmetic.

use serde_json::Value;

use super::Result;

pub const SCHEMA_REFERENCE_BOX: &str = "p192-wcm-reference-box-v1";
pub const EXPERIMENT_ID: &str = "EXP-SCURVE-1a8daf";
pub const PROTOCOL_VERSION: u64 = 2;
pub const CURVE_UID: &str =
    "urn:ec-record:1:sha256:5531c4a08bdb64b6e86a6e30e9a08aa57edef7af15ac5f6d4d2a83a53bf2f646";
pub const ICV1: &str = "icv1-fp192-t31607402316713927207482677199-52e4af59";
pub const EC1: &str = "EC1P192Cp192h5531c4a08bdb";

pub fn verdicts(check_ids: &[&str]) -> Value {
    Value::Array(
        check_ids
            .iter()
            .map(|check_id| serde_json::json!({"check_id":check_id,"status":"PASS"}))
            .collect(),
    )
}

fn reject_noncanonical_numbers(value: &Value) -> Result<()> {
    match value {
        Value::Null | Value::Bool(_) | Value::String(_) => Ok(()),
        Value::Number(number) => {
            const INTEROPERABLE_MAX: u64 = 9_007_199_254_740_991;
            let interoperable = number
                .as_u64()
                .map(|value| value <= INTEROPERABLE_MAX)
                .or_else(|| {
                    number
                        .as_i64()
                        .map(|value| value >= -(INTEROPERABLE_MAX as i64))
                })
                .unwrap_or(false);
            if interoperable {
                Ok(())
            } else {
                Err(
                    "canonical producer JSON requires interoperable-range integer numbers"
                        .to_owned(),
                )
            }
        }
        Value::Array(values) => {
            for value in values {
                reject_noncanonical_numbers(value)?;
            }
            Ok(())
        }
        Value::Object(values) => {
            for value in values.values() {
                reject_noncanonical_numbers(value)?;
            }
            Ok(())
        }
    }
}

/// Serialize the ASCII/integer-only schema subset as RFC 8785 bytes.
///
/// `serde_json` uses lexicographically ordered maps unless its
/// `preserve_order` feature is enabled (it is not enabled here).  The schema
/// excludes floats and non-ASCII strings, avoiding the two representation
/// cases where generic JSON serialization is not by itself a JCS contract.
pub fn canonical_json_bytes(value: &Value) -> Result<Vec<u8>> {
    reject_noncanonical_numbers(value)?;
    serde_json::to_vec(value).map_err(|error| format!("serialize canonical JSON: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn canonical_json_is_compact_sorted_and_has_no_newline() {
        assert_eq!(
            canonical_json_bytes(&json!({"z": 2, "a": "x", "m": null})).unwrap(),
            br#"{"a":"x","m":null,"z":2}"#
        );
        assert!(canonical_json_bytes(&json!({"not_exact": 0.5})).is_err());
        assert!(canonical_json_bytes(&json!({"too_large": 9_007_199_254_740_992u64})).is_err());
    }
}
