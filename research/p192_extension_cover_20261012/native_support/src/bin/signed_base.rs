//! Exact sign-closed three-sum census for the fixed P-192 Type I quotient base.
//! The 32 signed points still have only 16 independent logarithm columns.

#[path = "../ec.rs"]
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

fn dec(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 10).expect("decimal integer")
}

fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn point_key(point: &Point) -> String {
    match point {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => format!("{}:{}", x.value, y.value),
    }
}

fn eval(coeff: &[BigUint], x: &FieldElement, p: &BigUint) -> FieldElement {
    let mut result = FieldElement::zero(p.clone());
    for c in coeff.iter().rev() {
        result = result.mul(x).add(&FieldElement::new(c.clone(), p.clone()));
    }
    result
}

fn character(point: &Point, coeff: &[BigUint], p: &BigUint) -> char {
    let Point::Affine { x, .. } = point else {
        return 'I';
    };
    let value = eval(coeff, x, p);
    if value.is_zero() {
        return '0';
    }
    let legendre = value.pow(&((p - BigUint::one()) >> 1usize));
    if legendre.value == BigUint::one() {
        '+'
    } else {
        assert_eq!(legendre.value, p - BigUint::one());
        '-'
    }
}

fn audit(base: &Value, density: &Value, collector: &Value) -> Value {
    let curve = curve_zoo::CurveParams::p192();
    let p = &curve.p;
    let r = &curve.n;
    assert_eq!(dec(base["field_prime"].as_str().unwrap()), *p);
    assert_eq!(dec(density["p"].as_str().unwrap()), *p);
    assert_eq!(dec(density["r"].as_str().unwrap()), *r);
    let a = curve.a_fe();
    let sextic: Vec<BigUint> = density["sextic_coefficients_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| dec(v.as_str().unwrap()))
        .collect();
    assert_eq!(sextic.len(), 7);
    let unsigned: Vec<Point> = base["points"]
        .as_array()
        .unwrap()
        .iter()
        .map(|entry| {
            let point = Point::Affine {
                x: curve.fe(dec(entry["x"].as_str().unwrap())),
                y: curve.fe(dec(entry["y"].as_str().unwrap())),
            };
            assert!(curve.is_on_curve(&point));
            assert_eq!(point.scalar_mul(r, &a), Point::Infinity);
            point
        })
        .collect();
    assert_eq!(unsigned.len(), 16);
    let mut signed = unsigned.clone();
    signed.extend(unsigned.iter().map(Point::neg));
    assert_eq!(
        signed.iter().map(point_key).collect::<BTreeSet<_>>().len(),
        32
    );
    let negative_half = (r + BigUint::one()) >> 1usize;
    let lifted: Vec<Point> = signed
        .iter()
        .map(|q| q.scalar_mul(&negative_half, &a).neg())
        .collect();
    for (q, p_query) in signed.iter().zip(&lifted) {
        assert_eq!(p_query.double_vartime(&a).neg(), *q);
    }

    let mut target_rows: BTreeMap<String, Vec<[i8; 16]>> = BTreeMap::new();
    let mut target_points: BTreeMap<String, Point> = BTreeMap::new();
    let mut target_preimages: BTreeMap<String, Point> = BTreeMap::new();
    let mut cancellation_triples = 0usize;
    let mut cancellation_targets = BTreeSet::new();
    let mut noncancellation_targets = BTreeSet::new();
    for i in 0..32 {
        for j in i..32 {
            for k in j..32 {
                let mut row = [0i8; 16];
                for idx in [i, j, k] {
                    row[idx % 16] += if idx < 16 { 1 } else { -1 };
                }
                let cancellation = row.iter().map(|v| v.abs() as usize).sum::<usize>() < 3;
                let q = signed[i]
                    .add_vartime(&signed[j], &a)
                    .add_vartime(&signed[k], &a);
                let preimage = lifted[i]
                    .add_vartime(&lifted[j], &a)
                    .add_vartime(&lifted[k], &a);
                assert_eq!(preimage.double_vartime(&a).neg(), q);
                let key = point_key(&q);
                if cancellation {
                    cancellation_triples += 1;
                    cancellation_targets.insert(key.clone());
                } else {
                    noncancellation_targets.insert(key.clone());
                }
                target_points.entry(key.clone()).or_insert(q);
                target_preimages.entry(key.clone()).or_insert(preimage);
                target_rows.entry(key).or_default().push(row);
            }
        }
    }
    assert_eq!(target_rows.values().map(Vec::len).sum::<usize>(), 5984);
    assert_eq!(cancellation_triples, 512);
    assert_eq!(cancellation_targets.len(), 32);
    assert_eq!(noncancellation_targets.len(), 5472);
    assert!(cancellation_targets.is_disjoint(&noncancellation_targets));
    assert_eq!(target_rows.len(), target_points.len());
    assert_eq!(target_rows.len(), 5504);
    let mut cells: BTreeMap<String, usize> = BTreeMap::new();
    let mut accepted_keys = BTreeSet::new();
    for (key, q) in &target_points {
        let preimage = &target_preimages[key];
        let first = character(preimage, &sextic, p);
        let second = character(&preimage.double_vartime(&a), &sextic, p);
        let cell = format!("{first}{second}");
        *cells.entry(cell.clone()).or_default() += 1;
        if cell == "++" {
            accepted_keys.insert(key.clone());
        }
        assert!(curve.is_on_curve(q));
    }

    let generator = curve.generator();
    let mut hits = Vec::new();
    for (index, scalar) in collector["accepted_scalars_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
    {
        let scalar = dec(scalar.as_str().unwrap());
        assert!(!scalar.is_zero() && scalar < *r);
        let preimage = generator.scalar_mul(&scalar, &a);
        assert_eq!(character(&preimage, &sextic, p), '+');
        assert_eq!(character(&preimage.double_vartime(&a), &sextic, p), '+');
        let q = preimage.double_vartime(&a).neg();
        if let Some(rows) = target_rows.get(&point_key(&q)) {
            hits.push(json!({"query_index": index, "scalar": scalar.to_string(), "rows": rows}));
        }
    }
    let lower = dec(density["accepted_nonzero_scalar_lower_bound"]
        .as_str()
        .unwrap());
    let accepted = accepted_keys.len();
    json!({
        "status": "exact_signed_type_i_fixed_base_support_census",
        "field_prime": p.to_string(), "subgroup_order": r.to_string(),
        "unsigned_base_classes": 16, "signed_base_classes": 32,
        "logarithm_columns_before_orbit_folding": 16,
        "unordered_signed_triples": 5984,
        "cancellation_triples": cancellation_triples,
        "cancellation_targets": cancellation_targets.len(),
        "noncancellation_targets": noncancellation_targets.len(),
        "distinct_three_sum_targets": target_rows.len(),
        "counts_by_character_pair": cells,
        "accepted_support_targets": accepted,
        "accepted_nonzero_scalar_lower_bound": lower.to_string(),
        "conditional_hit_probability_upper_bound": format!("{accepted}/{lower}"),
        "frozen_known_scalar_queries": collector["accepted_queries"],
        "frozen_query_hits": hits,
        "target_keys_sha256": hash(target_rows.keys().cloned().collect::<Vec<_>>().join("\n").as_bytes()),
        "accepted_keys_sha256": hash(accepted_keys.into_iter().collect::<Vec<_>>().join("\n").as_bytes()),
        "query_equation": "q=a(-2G), with signed base entries ±B_i and coefficient columns B_i",
        "geometry_status": "The 16 negative quotient classes are exact subgroup points; fixed-origin genus-three divisor representatives are a separate construction obligation."
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_signed_support_and_collector() {
        let base: Value = serde_json::from_str(include_str!(
            "../../../type_i_base_points_legacy_fixture.json"
        ))
        .unwrap();
        let density: Value =
            serde_json::from_str(include_str!("../../../type_i_character_density.json")).unwrap();
        let collector: Value = serde_json::from_str(include_str!(
            "../../../type_i_geom_character_collector.json"
        ))
        .unwrap();
        let result = audit(&base, &density, &collector);
        assert_eq!(result["distinct_three_sum_targets"], 5504);
        assert_eq!(result["accepted_support_targets"], 1344);
        assert_eq!(result["cancellation_targets"], 32);
        assert_eq!(result["frozen_query_hits"].as_array().unwrap().len(), 0);
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err("usage: signed_base BASE.json DENSITY.json COLLECTOR.json OUTPUT.json".into());
    }
    let input: Vec<Vec<u8>> = args[1..4].iter().map(fs::read).collect::<Result<_, _>>()?;
    let values: Vec<Value> = input
        .iter()
        .map(|v| serde_json::from_slice(v))
        .collect::<Result<_, _>>()?;
    let mut result = audit(&values[0], &values[1], &values[2]);
    result["input_sha256"] = json!(input.iter().map(|v| hash(v)).collect::<Vec<_>>());
    result["source_sha256"] = json!(hash(include_bytes!("signed_base.rs")));
    result["ec_source_sha256"] = json!(hash(include_bytes!("../ec.rs")));
    fs::write(
        &args[4],
        [serde_json::to_vec_pretty(&result)?, b"\n".to_vec()].concat(),
    )?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
