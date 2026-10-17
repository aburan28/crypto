//! Frozen adaptive-base audit: one new quotient point per known-scalar query.
//! The new relation has a private logarithm column unless the point already
//! lies in the old signed base or has a second old-base decomposition.
#![recursion_limit = "256"]
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
use std::time::Instant;

fn dec(v: &Value) -> BigUint {
    BigUint::parse_bytes(v.as_str().unwrap().as_bytes(), 10).unwrap()
}
fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn key(q: &Point) -> String {
    match q {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => format!("{}:{}", x.value, y.value),
    }
}
fn eval(coeff: &[BigUint], x: &FieldElement, p: &BigUint) -> FieldElement {
    coeff
        .iter()
        .rev()
        .fold(FieldElement::zero(p.clone()), |z, c| {
            z.mul(x).add(&FieldElement::new(c.clone(), p.clone()))
        })
}
fn character(q: &Point, coeff: &[BigUint], p: &BigUint) -> char {
    let Point::Affine { x, .. } = q else {
        return 'I';
    };
    let v = eval(coeff, x, p);
    if v.is_zero() {
        return '0';
    }
    if v.pow(&((p - BigUint::one()) >> 1usize)).value == BigUint::one() {
        '+'
    } else {
        '-'
    }
}
fn cell(pre: &Point, coeff: &[BigUint], p: &BigUint, a: &FieldElement) -> String {
    format!(
        "{}{}",
        character(pre, coeff, p),
        character(&pre.double_vartime(a), coeff, p)
    )
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        return Err(
            "usage: adaptive_private_columns BASE.json DENSITY.json GUIDE.json MIXED.json OUT.json"
                .into(),
        );
    }
    let start = Instant::now();
    let inputs: Vec<Vec<u8>> = args[1..5].iter().map(fs::read).collect::<Result<_, _>>()?;
    let docs: Vec<Value> = inputs
        .iter()
        .map(|b| serde_json::from_slice(b))
        .collect::<Result<_, _>>()?;
    let (base, density, guide, mixed) = (&docs[0], &docs[1], &docs[2], &docs[3]);
    let curve = curve_zoo::CurveParams::p192();
    let p = &curve.p;
    let r = &curve.n;
    let a = curve.a_fe();
    let g = curve.generator();
    assert_eq!(dec(&base["field_prime"]), *p);
    assert_eq!(dec(&density["p"]), *p);
    assert_eq!(dec(&density["r"]), *r);
    assert_eq!(dec(&mixed["field_prime"]), *p);
    assert_eq!(dec(&mixed["subgroup_order"]), *r);
    assert_eq!(mixed["input_sha256"][0], hash(&inputs[0]));
    assert_eq!(mixed["input_sha256"][1], hash(&inputs[1]));
    assert_eq!(mixed["input_sha256"][2], hash(&inputs[2]));
    let coeff: Vec<BigUint> = density["sextic_coefficients_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(dec)
        .collect();
    let unsigned: Vec<Point> = base["points"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            let q = Point::Affine {
                x: curve.fe(dec(&v["x"])),
                y: curve.fe(dec(&v["y"])),
            };
            assert!(curve.is_on_curve(&q));
            q
        })
        .collect();
    assert_eq!(unsigned.len(), 16);
    let guide_scalars: Vec<BigUint> = guide["accepted_scalars_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(dec)
        .collect();
    assert_eq!(guide_scalars.len(), 64);
    let mut old_points = unsigned.clone();
    old_points.extend(unsigned.iter().map(Point::neg));
    for s in &guide_scalars {
        old_points.push(g.scalar_mul(s, &a).double_vartime(&a).neg());
    }
    let guided_positive = old_points[32..96].to_vec();
    old_points.extend(guided_positive.iter().map(Point::neg));
    assert_eq!(old_points.len(), 160);
    let old_keys: BTreeSet<String> = old_points.iter().map(key).collect();
    assert_eq!(old_keys.len(), 160);
    let minus_half = (r - BigUint::one()) >> 1usize;
    let fresh = mixed["independent_accepted_queries"].as_array().unwrap();
    assert_eq!(fresh.len(), 64);
    let mut candidates = Vec::new();
    let mut candidate_keys = BTreeMap::new();
    let mut cell_counts = BTreeMap::<String, usize>::new();
    let mut duplicate_new_points = 0usize;
    let mut old_base_points = 0usize;
    let mut identity_new_points = 0usize;
    let mut query_points = Vec::new();
    for (j, record) in fresh.iter().enumerate() {
        let s = dec(&record["scalar"]);
        assert!(!s.is_zero() && s < *r);
        let pre = g.scalar_mul(&s, &a);
        assert_eq!(cell(&pre, &coeff, p, &a), "++");
        let q = pre.double_vartime(&a).neg();
        query_points.push(q.clone());
        let i = j % 16;
        let k = (7 * j + 3) % 16;
        assert_ne!(i, k);
        let c = q
            .add_vartime(&unsigned[i].neg(), &a)
            .add_vartime(&unsigned[k].neg(), &a);
        assert_eq!(
            unsigned[i]
                .add_vartime(&unsigned[k], &a)
                .add_vartime(&c, &a),
            q
        );
        assert!(curve.is_on_curve(&c));
        assert_eq!(c.scalar_mul(r, &a), Point::Infinity);
        let ckey = key(&c);
        if c == Point::Infinity {
            identity_new_points += 1;
        }
        let cpre = c.scalar_mul(&minus_half, &a);
        let ccell = cell(&cpre, &coeff, p, &a);
        *cell_counts.entry(ccell.clone()).or_default() += 1;
        if old_keys.contains(&ckey) {
            old_base_points += 1;
        }
        if candidate_keys.insert(ckey.clone(), j).is_some() {
            duplicate_new_points += 1;
        }
        candidates.push(json!({"index":j,"query_scalar":s.to_string(),
            "query_point_key":key(&q),"base_indices":[i,k],
            "new_point_key":ckey,"new_point_character_cell":ccell}));
    }
    let mut old_unknown_pairs = Vec::new();
    for i in 0..16 {
        for k in i..16 {
            old_unknown_pairs.push((i, k, unsigned[i].add_vartime(&unsigned[k], &a)));
        }
    }
    assert_eq!(old_unknown_pairs.len(), 136);
    let mut full_pair_keys = BTreeMap::<String, (usize, usize, usize)>::new();
    let mut full_pair_duplicates = 0usize;
    let mut full_pair_old_base_points = 0usize;
    let mut full_pair_identity_points = 0usize;
    let mut first_full_pair_collisions = Vec::new();
    for (j, q) in query_points.iter().enumerate() {
        for (i, k, pair) in &old_unknown_pairs {
            let c = q.add_vartime(&pair.neg(), &a);
            let ckey = key(&c);
            if c == Point::Infinity {
                full_pair_identity_points += 1;
            }
            if old_keys.contains(&ckey) {
                full_pair_old_base_points += 1;
            }
            if let Some(previous) = full_pair_keys.insert(ckey.clone(), (j, *i, *k)) {
                full_pair_duplicates += 1;
                if first_full_pair_collisions.len() < 20 {
                    first_full_pair_collisions.push(json!({"target":ckey,
                        "first":previous,"second":[j,i,k]}));
                }
            }
        }
    }
    assert!(candidate_keys
        .keys()
        .all(|k| full_pair_keys.contains_key(k)));
    let full_pair_keys_sha256 = hash(
        full_pair_keys
            .keys()
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
            .as_bytes(),
    );
    let setup_ms = start.elapsed().as_millis();
    let pairs: Vec<Vec<Point>> = (0..160)
        .map(|i| {
            (i..160)
                .map(|j| old_points[i].add_vartime(&old_points[j], &a))
                .collect()
        })
        .collect();
    let pair_index_ms = start.elapsed().as_millis() - setup_ms;
    let mut old_support_hits = 0usize;
    let mut old_support_candidates = BTreeSet::new();
    let mut first_witness = Vec::new();
    let mut full_pair_support_hits = 0usize;
    let mut full_pair_support_targets = BTreeSet::new();
    let mut first_full_pair_support_witnesses = Vec::new();
    for i in 0..160 {
        for j in i..160 {
            for k in j..160 {
                let target = pairs[i][j - i].add_vartime(&old_points[k], &a);
                let target_key = key(&target);
                if let Some(candidate_index) = candidate_keys.get(&target_key) {
                    old_support_hits += 1;
                    old_support_candidates.insert(*candidate_index);
                    if first_witness.len() < 20 {
                        first_witness.push(json!({"candidate_index":candidate_index,
                            "old_base_indices":[i,j,k]}));
                    }
                }
                if let Some(source) = full_pair_keys.get(&target_key) {
                    full_pair_support_hits += 1;
                    full_pair_support_targets.insert(target_key);
                    if first_full_pair_support_witnesses.len() < 20 {
                        first_full_pair_support_witnesses.push(json!({"adaptive_source":source,
                            "old_base_indices":[i,j,k]}));
                    }
                }
            }
        }
    }
    let private_identity =
        duplicate_new_points == 0 && old_base_points == 0 && identity_new_points == 0;
    let full_pair_private_identity = full_pair_duplicates == 0
        && full_pair_old_base_points == 0
        && full_pair_identity_points == 0;
    let out = json!({
        "status":"exact_adaptive_private_column_audit",
        "field_prime":p.to_string(),"subgroup_order":r.to_string(),
        "query_law":"first 64 frozen independent accepted query scalars; i=j mod 16, k=(7j+3) mod 16; C_j=Q_j-B_i-B_k",
        "old_unknown_columns":16,"old_known_signed_points":128,
        "new_relations":64,"new_candidate_points":64,
        "distinct_new_points":candidate_keys.len(),
        "new_point_duplicates":duplicate_new_points,
        "identity_new_points":identity_new_points,
        "new_points_already_in_old_base":old_base_points,
        "new_point_character_cells":cell_counts,
        "new_points_with_old_three_sum_decomposition":old_support_candidates.len(),
        "old_three_sum_witnesses":old_support_hits,
        "first_old_three_sum_witnesses":first_witness,
        "private_identity_block":private_identity,
        "augmented_relation_rank_if_private":if private_identity {json!(64)} else {Value::Null},
        "augmented_unknown_columns_if_private":if private_identity {json!(80)} else {Value::Null},
        "old_unknown_nullity_if_private":if private_identity {json!(16)} else {Value::Null},
        "full_pair_policy":"all 64 accepted known-scalar targets times all 136 unordered old-unknown pairs, with repetition",
        "full_pair_candidate_rows":64*136,
        "full_pair_distinct_points":full_pair_keys.len(),
        "full_pair_duplicate_points":full_pair_duplicates,
        "full_pair_identity_points":full_pair_identity_points,
        "full_pair_points_already_in_old_base":full_pair_old_base_points,
        "full_pair_keys_sha256":full_pair_keys_sha256,
        "full_pair_first_collisions":first_full_pair_collisions,
        "full_pair_points_with_old_three_sum_decomposition":full_pair_support_targets.len(),
        "full_pair_old_three_sum_witnesses":full_pair_support_hits,
        "full_pair_first_old_three_sum_witnesses":first_full_pair_support_witnesses,
        "full_pair_private_identity_block":full_pair_private_identity,
        "full_pair_rank_if_private":if full_pair_private_identity {json!(64*136)} else {Value::Null},
        "full_pair_unknown_columns_if_private":if full_pair_private_identity {json!(16+64*136)} else {Value::Null},
        "full_pair_old_unknown_nullity_if_private":if full_pair_private_identity {json!(16)} else {Value::Null},
        "candidates":candidates,
        "input_sha256":inputs.iter().map(|b|hash(b)).collect::<Vec<_>>(),
        "source_sha256":hash(include_bytes!("adaptive_private_columns.rs")),
        "ec_source_sha256":hash(include_bytes!("../ec.rs")),
        "setup_ms":setup_ms,"pair_index_ms":pair_index_ms,
        "support_scan_ms":start.elapsed().as_millis()-setup_ms-pair_index_ms
    });
    fs::write(
        &args[5],
        [serde_json::to_vec_pretty(&out)?, b"\n".to_vec()].concat(),
    )?;
    println!(
        "ADAPTIVE_PRIVATE_COLUMNS_PASS distinct={} old_support={} old_nullity={}",
        candidate_keys.len(),
        old_support_candidates.len(),
        out["old_unknown_nullity_if_private"]
    );
    Ok(())
}
