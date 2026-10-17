//! Native finite-support census for the P-192 Type I character-guided law.
//! The 16 quotient points are frozen exports of the existing checked base.

mod ec;
use ec::{curve_zoo, FieldElement, Point};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::path::Path;

fn dec(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 10).expect("decimal integer")
}

fn sha(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn coefficients(input: &Value, key: &str) -> Vec<BigUint> {
    input[key]
        .as_array()
        .expect("coefficient array")
        .iter()
        .map(|x| dec(x.as_str().expect("decimal coefficient")))
        .collect()
}

fn eval(coeff: &[BigUint], x: &FieldElement, p: &BigUint) -> FieldElement {
    let mut result = FieldElement::zero(p.clone());
    for c in coeff.iter().rev() {
        result = result.mul(x).add(&FieldElement::new(c.clone(), p.clone()));
    }
    result
}

fn key(point: &Point) -> String {
    match point {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => format!("{}:{}", x.value, y.value),
    }
}

fn character(point: &Point, sextic: &[BigUint], factors: &[Vec<BigUint>], p: &BigUint) -> char {
    let Point::Affine { x, .. } = point else {
        return 'I';
    };
    let value = eval(sextic, x, p);
    if value.is_zero() {
        return '0';
    }
    let exponent = (p - BigUint::one()) >> 1usize;
    let legendre = value.pow(&exponent);
    let factor_product = factors.iter().fold(
        FieldElement::new(sextic[6].clone(), p.clone()),
        |acc, factor| acc.mul(&eval(factor, x, p)),
    );
    assert_eq!(value, factor_product, "independent cubic-factor evaluation");
    if legendre.value == BigUint::one() {
        '+'
    } else {
        assert_eq!(legendre.value, p - BigUint::one());
        '-'
    }
}

fn audit(fixture: &Value, density: &Value, collector: &Value) -> Value {
    let curve = curve_zoo::CurveParams::p192();
    let p = curve.p.clone();
    let r = curve.n.clone();
    assert_eq!(dec(fixture["field_prime"].as_str().unwrap()), p);
    assert_eq!(dec(fixture["subgroup_order"].as_str().unwrap()), r);
    assert_eq!(dec(density["p"].as_str().unwrap()), p);
    assert_eq!(dec(density["r"].as_str().unwrap()), r);
    let sextic = coefficients(density, "sextic_coefficients_decimal");
    assert_eq!(sextic.len(), 7);
    let factors: Vec<Vec<BigUint>> = density["sextic_factor_coefficients_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|c| dec(c.as_str().unwrap()))
                .collect()
        })
        .collect();
    assert_eq!(factors.len(), 2);
    assert!(factors.iter().all(|factor| factor.len() == 4));

    let a = curve.a_fe();
    let points: Vec<Point> = fixture["points"]
        .as_array()
        .unwrap()
        .iter()
        .map(|entry| {
            let point = Point::Affine {
                x: curve.fe(dec(entry["x"].as_str().unwrap())),
                y: curve.fe(dec(entry["y"].as_str().unwrap())),
            };
            assert!(curve.is_on_curve(&point));
            assert_eq!(point.scalar_mul(&r, &a), Point::Infinity);
            assert_ne!(point, Point::Infinity);
            point
        })
        .collect();
    assert_eq!(points.len(), 16);
    let half = (&r + BigUint::one()) >> 1usize;
    assert_eq!((&half << 1usize) % &r, BigUint::one());
    let negative_half: Vec<Point> = points
        .iter()
        .map(|point| point.scalar_mul(&half, &a).neg())
        .collect();
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut targets = BTreeSet::new();
    let mut cell_keys = Vec::new();
    for i in 0..16 {
        for j in i..16 {
            for k in j..16 {
                let q = points[i]
                    .add_vartime(&points[j], &a)
                    .add_vartime(&points[k], &a);
                assert!(curve.is_on_curve(&q));
                let q_key = key(&q);
                assert!(targets.insert(q_key.clone()), "three-sum collision");
                let p_query = negative_half[i]
                    .add_vartime(&negative_half[j], &a)
                    .add_vartime(&negative_half[k], &a);
                assert_eq!(p_query.double_vartime(&a).neg(), q);
                let first = character(&p_query, &sextic, &factors, &p);
                let twice = p_query.double_vartime(&a);
                let second = character(&twice, &sextic, &factors, &p);
                let cell = format!("{first}{second}");
                *counts.entry(cell.clone()).or_default() += 1;
                cell_keys.push(format!("{q_key}|{cell}"));
            }
        }
    }
    assert_eq!(targets.len(), 816);
    assert_eq!(counts.get("++"), Some(&224));
    assert_eq!(counts.get("+-"), Some(&211));
    assert_eq!(counts.get("-+"), Some(&203));
    assert_eq!(counts.get("--"), Some(&178));
    assert_eq!(counts.values().sum::<usize>(), 816);
    let accepted_scalars = collector["accepted_scalars_decimal"].as_array().unwrap();
    assert_eq!(accepted_scalars.len(), 64);
    assert_eq!(collector["accepted_queries"], 64);
    let generator = curve.generator();
    assert!(curve.is_on_curve(&generator));
    let mut collector_hits = 0usize;
    for scalar in accepted_scalars {
        let scalar = dec(scalar.as_str().unwrap());
        assert!(!scalar.is_zero() && scalar < r);
        let p_query = generator.scalar_mul(&scalar, &a);
        assert_eq!(character(&p_query, &sextic, &factors, &p), '+');
        let twice = p_query.double_vartime(&a);
        assert_eq!(character(&twice, &sextic, &factors, &p), '+');
        let q = twice.neg();
        assert_eq!(
            q,
            generator.double_vartime(&a).neg().scalar_mul(&scalar, &a)
        );
        if targets.contains(&key(&q)) {
            collector_hits += 1;
        }
    }
    assert_eq!(collector_hits, 0);
    assert_eq!(collector["exact_three_sum_matches"], 0);
    cell_keys.sort();
    let lower = dec(density["accepted_nonzero_scalar_lower_bound"]
        .as_str()
        .unwrap());
    assert!(BigUint::from(224u32) * (BigUint::one() << 182usize) < lower);
    json!({
        "status": "native_exact_type_i_fixed_base_character_support_partition",
        "field_prime": p.to_string(),
        "subgroup_order": r.to_string(),
        "base_classes": points.len(),
        "unordered_triples": 816,
        "distinct_three_sum_targets": targets.len(),
        "counts_by_character_pair": counts,
        "accepted_support_targets": 224,
        "frozen_known_scalar_queries_replayed": accepted_scalars.len(),
        "frozen_known_scalar_three_sum_hits": collector_hits,
        "accepted_nonzero_scalar_lower_bound": lower.to_string(),
        "conditioned_hit_probability_upper_bound": format!("224/{lower}"),
        "target_cell_sha256": sha(cell_keys.join("\n").as_bytes()),
        "query_equation": "q=-2P, P=aG for seeded known-scalar queries",
        "boundary": "Native arithmetic on frozen, Sage-verified quotient-point inputs; cubic residuals and original Fp3-character derivation have separate certificates."
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err(
            "usage: p192-type-i-conditioned-support BASE.json DENSITY.json COLLECTOR.json OUTPUT.json".into(),
        );
    }
    let fixture_bytes = fs::read(&args[1])?;
    let density_bytes = fs::read(&args[2])?;
    let collector_bytes = fs::read(&args[3])?;
    let fixture: Value = serde_json::from_slice(&fixture_bytes)?;
    let density: Value = serde_json::from_slice(&density_bytes)?;
    let collector: Value = serde_json::from_slice(&collector_bytes)?;
    let mut result = audit(&fixture, &density, &collector);
    result["base_fixture_sha256"] = json!(sha(&fixture_bytes));
    result["density_receipt_sha256"] = json!(sha(&density_bytes));
    result["collector_receipt_sha256"] = json!(sha(&collector_bytes));
    result["source_sha256"] = json!(sha(include_bytes!("main.rs")));
    result["ec_source_sha256"] = json!(sha(include_bytes!("ec.rs")));
    result["base_source_receipt_sha256"] = fixture["base_receipt_sha256"].clone();
    result["base_source_verifier_sha256"] = fixture["base_verifier_source_sha256"].clone();
    let output = serde_json::to_vec_pretty(&result)?;
    fs::write(Path::new(&args[4]), [&output[..], b"\n"].concat())?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_base_and_character_partition() {
        let fixture: Value =
            serde_json::from_str(include_str!("../../type_i_base_points_legacy_fixture.json"))
                .unwrap();
        let density: Value =
            serde_json::from_str(include_str!("../../type_i_character_density.json")).unwrap();
        let collector: Value =
            serde_json::from_str(include_str!("../../type_i_geom_character_collector.json"))
                .unwrap();
        let receipt = audit(&fixture, &density, &collector);
        assert_eq!(receipt["accepted_support_targets"], 224);
        assert_eq!(receipt["distinct_three_sum_targets"], 816);
        assert_eq!(receipt["frozen_known_scalar_queries_replayed"], 64);
        assert_eq!(receipt["frozen_known_scalar_three_sum_hits"], 0);
    }
}
