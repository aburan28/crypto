//! Frozen recursive Pocklington certificate for the discriminant cofactor C.

use std::collections::{BTreeMap, BTreeSet};

use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};

use super::arithmetic::C_DEC;
use super::Result;

const Q_DEC: &str = "3654106800343397140285403541412623476381";

const C_MINUS_ONE: [(&str, u32); 8] = [
    ("2", 1),
    ("3", 2),
    ("7", 1),
    ("17", 1),
    ("139", 1),
    ("11471", 1),
    ("1133039", 1),
    (Q_DEC, 1),
];
const Q_MINUS_ONE: [(&str, u32); 5] = [
    ("2", 2),
    ("5", 1),
    ("28929853", 1),
    ("386058915559", 1),
    ("16358799425486714897", 1),
];
const P28929853_MINUS_ONE: [(&str, u32); 5] =
    [("2", 2), ("3", 3), ("7", 1), ("17", 1), ("2251", 1)];
const P386058915559_MINUS_ONE: [(&str, u32); 4] = [("2", 1), ("3", 3), ("17", 1), ("420543481", 1)];
const P16358799425486714897_MINUS_ONE: [(&str, u32); 4] =
    [("2", 4), ("7", 1), ("449", 1), ("325302247563767", 1)];
const P1133039_MINUS_ONE: [(&str, u32); 3] = [("2", 1), ("397", 1), ("1427", 1)];
const P420543481_MINUS_ONE: [(&str, u32); 6] = [
    ("2", 3),
    ("3", 1),
    ("5", 1),
    ("7", 2),
    ("37", 1),
    ("1933", 1),
];
const P325302247563767_MINUS_ONE: [(&str, u32); 5] =
    [("2", 1), ("11", 1), ("839", 1), ("39371", 1), ("447637", 1)];
const P447637_MINUS_ONE: [(&str, u32); 4] = [("2", 2), ("3", 1), ("7", 1), ("73", 2)];

fn unsigned(text: &str) -> BigUint {
    BigUint::parse_bytes(text.as_bytes(), 10).expect("frozen Pocklington integer")
}

fn gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }
    left
}

fn trial_prime(value: u32) -> bool {
    if value < 2 {
        return false;
    }
    if value.is_multiple_of(2) {
        return value == 2;
    }
    let mut divisor = 3u32;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }
    true
}

fn witness(value: &str) -> Result<&'static str> {
    match value {
        C_DEC | Q_DEC | "447637" => Ok("2"),
        "28929853" => Ok("15"),
        "386058915559" => Ok("6"),
        "16358799425486714897" => Ok("3"),
        "1133039" => Ok("13"),
        "420543481" => Ok("11"),
        "325302247563767" => Ok("5"),
        _ => Err(format!("no frozen Pocklington witness for {value}")),
    }
}

fn node(value_text: &str, factors: &[(&str, u32)], children: Vec<Value>) -> Result<Value> {
    let value = unsigned(value_text);
    let mut child_by_value = BTreeMap::new();
    for child in children {
        let key = child["value"]
            .as_str()
            .ok_or_else(|| "recursive Pocklington child has no value".to_owned())?
            .to_owned();
        if child_by_value.insert(key.clone(), child).is_some() {
            return Err(format!("duplicate recursive Pocklington child value {key}"));
        }
    }
    let mut product = BigUint::one();
    let mut factor_records = Vec::new();
    let mut witness_records = Vec::new();
    let common_base = witness(value_text)?;
    let mut seen = BTreeSet::new();
    let mut previous = BigUint::zero();
    for &(prime_text, exponent) in factors {
        if !seen.insert(prime_text) {
            return Err(format!("duplicate frozen Pocklington factor {prime_text}"));
        }
        let prime = unsigned(prime_text);
        if prime <= previous || exponent == 0 {
            return Err(format!(
                "non-increasing or zero-exponent factor {prime_text}"
            ));
        }
        previous = prime.clone();
        product *= prime.pow(exponent);
        let proof = if let Some(child) = child_by_value.get(prime_text) {
            child.clone()
        } else {
            let leaf = prime_text
                .parse::<u32>()
                .map_err(|_| format!("Pocklington leaf does not fit u32: {prime_text}"))?;
            if !trial_prime(leaf) {
                return Err(format!("Pocklington leaf is composite: {prime_text}"));
            }
            json!({"method":"trial_division_u32","value":prime_text})
        };
        factor_records.push(json!({
            "exponent": exponent,
            "prime": prime_text,
            "proof": proof,
        }));
        witness_records.push(json!({"base":common_base,"prime":prime_text}));

        let base = unsigned(common_base);
        let value_minus_one = &value - 1u8;
        if base.modpow(&value_minus_one, &value) != BigUint::one() {
            return Err(format!(
                "Pocklington Fermat equality failed for {value_text}"
            ));
        }
        let probe = base.modpow(&(&value_minus_one / &prime), &value);
        let delta = if probe.is_zero() {
            &value - 1u8
        } else {
            probe - 1u8
        };
        if !gcd(delta, value.clone()).is_one() {
            return Err(format!(
                "Pocklington gcd equality failed for {value_text}/{prime_text}"
            ));
        }
    }
    if product != &value - 1u8 {
        return Err(format!("frozen factors do not multiply to {value_text}-1"));
    }
    if &product * &product <= value {
        return Err(format!(
            "Pocklington known-factor product does not satisfy F^2 > {value_text}"
        ));
    }
    if child_by_value
        .keys()
        .any(|child| !seen.contains(child.as_str()))
    {
        return Err(format!(
            "unused recursive Pocklington child for {value_text}"
        ));
    }
    Ok(json!({
        "cofactor": "1",
        "factors": factor_records,
        "method": "pocklington_v1",
        "value": value_text,
        "witnesses": witness_records,
    }))
}

pub fn c_proof() -> Result<Value> {
    let p447637 = node("447637", &P447637_MINUS_ONE, vec![])?;
    let p325302247563767 = node(
        "325302247563767",
        &P325302247563767_MINUS_ONE,
        vec![p447637],
    )?;
    let p420543481 = node("420543481", &P420543481_MINUS_ONE, vec![])?;
    let p1133039 = node("1133039", &P1133039_MINUS_ONE, vec![])?;
    let p28929853 = node("28929853", &P28929853_MINUS_ONE, vec![])?;
    let p386058915559 = node("386058915559", &P386058915559_MINUS_ONE, vec![p420543481])?;
    let p16358799425486714897 = node(
        "16358799425486714897",
        &P16358799425486714897_MINUS_ONE,
        vec![p325302247563767],
    )?;
    let q = node(
        Q_DEC,
        &Q_MINUS_ONE,
        vec![p28929853, p386058915559, p16358799425486714897],
    )?;
    node(C_DEC, &C_MINUS_ONE, vec![p1133039, q])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_c_certificate_builds_and_has_exact_variant_keys() {
        let proof = c_proof().unwrap();
        assert_eq!(proof["method"], "pocklington_v1");
        assert_eq!(proof["value"], C_DEC);
        assert_eq!(proof["cofactor"], "1");
        assert_eq!(proof["factors"].as_array().unwrap().len(), 8);
        assert_eq!(proof["witnesses"].as_array().unwrap().len(), 8);
    }

    #[test]
    fn duplicate_recursive_child_values_are_rejected() {
        let child = json!({"value":"397"});
        let error = node("1133039", &P1133039_MINUS_ONE, vec![child.clone(), child]).unwrap_err();
        assert!(error.contains("duplicate recursive Pocklington child"));
    }
}
