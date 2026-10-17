//! Exact collision audit for the signed P-192 Type I mixed quotient base.
//! A repeated elliptic target contributes a new relation only if its formal
//! signed base-point vector differs from a prior representation.
#[path = "../ec.rs"]
mod ec;

use ec::{curve_zoo, Point};
use num_bigint::BigUint;
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
fn formal(indices: [usize; 3]) -> [i8; 80] {
    let mut out = [0; 80];
    for i in indices {
        let (column, sign) = match i {
            0..=15 => (i, 1),
            16..=31 => (i - 16, -1),
            32..=95 => (i - 16, 1),
            96..=159 => (i - 80, -1),
            _ => unreachable!(),
        };
        out[column] += sign;
    }
    out
}
fn nonzero_unknown(row: &[i8; 80]) -> bool {
    row[..16].iter().any(|v| *v != 0)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        return Err("usage: formal_collision BASE.json COLLECTOR.json OUT.json".into());
    }
    let start = Instant::now();
    let base_bytes = fs::read(&args[1])?;
    let guide_bytes = fs::read(&args[2])?;
    let base: Value = serde_json::from_slice(&base_bytes)?;
    let guide: Value = serde_json::from_slice(&guide_bytes)?;
    let curve = curve_zoo::CurveParams::p192();
    let a = curve.a_fe();
    let g = curve.generator();
    assert_eq!(dec(&base["field_prime"]), curve.p);
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
    let scalars: Vec<BigUint> = guide["accepted_scalars_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(dec)
        .collect();
    assert_eq!(scalars.len(), 64);
    let mut points = unsigned.clone();
    points.extend(unsigned.iter().map(Point::neg));
    for s in &scalars {
        let q = g.scalar_mul(s, &a).double_vartime(&a).neg();
        points.push(q);
    }
    let positive_known = points[32..96].to_vec();
    points.extend(positive_known.iter().map(Point::neg));
    assert_eq!(points.len(), 160);
    assert_eq!(points.iter().map(key).collect::<BTreeSet<_>>().len(), 160);
    let setup_ms = start.elapsed().as_millis();
    let pairs: Vec<Vec<Point>> = (0..160)
        .map(|i| {
            (i..160)
                .map(|j| points[i].add_vartime(&points[j], &a))
                .collect()
        })
        .collect();
    let pair_index_ms = start.elapsed().as_millis() - setup_ms;

    let mut bearing = BTreeMap::<String, ([i8; 80], [usize; 3])>::new();
    let mut anchored = BTreeMap::<String, ([i8; 80], [usize; 3])>::new();
    let mut bearing_triples = 0usize;
    let mut known_only_triples = 0usize;
    let mut zero_unknown_triples = 0usize;
    let mut repeated_bearing = 0usize;
    let mut repeated_anchored = 0usize;
    let mut distinct_formal_bearing = BTreeSet::new();
    let mut distinct_formal_anchored = BTreeSet::new();
    let mut bearing_examples = Vec::new();
    let mut anchored_examples = Vec::new();
    let mut zero_unknown_targets = BTreeSet::new();
    let mut known_only_targets = BTreeSet::new();
    for i in 0..160 {
        for j in i..160 {
            for k in j..160 {
                let indices = [i, j, k];
                let qkey = key(&pairs[i][j - i].add_vartime(&points[k], &a));
                let vector = formal(indices);
                let is_bearing = nonzero_unknown(&vector);
                if is_bearing {
                    bearing_triples += 1;
                } else if i < 32 {
                    zero_unknown_triples += 1;
                    zero_unknown_targets.insert(qkey.clone());
                } else {
                    known_only_triples += 1;
                    known_only_targets.insert(qkey.clone());
                }
                let (map, repeats, novel, examples) = if is_bearing {
                    (
                        &mut bearing,
                        &mut repeated_bearing,
                        &mut distinct_formal_bearing,
                        &mut bearing_examples,
                    )
                } else {
                    (
                        &mut anchored,
                        &mut repeated_anchored,
                        &mut distinct_formal_anchored,
                        &mut anchored_examples,
                    )
                };
                if let Some((old_vector, old_indices)) = map.get(&qkey) {
                    *repeats += 1;
                    if old_vector != &vector && novel.insert(qkey.clone()) && examples.len() < 20 {
                        examples.push(json!({"target":qkey,"first":old_indices,"second":indices,
                            "first_vector":old_vector.to_vec(),"second_vector":vector.to_vec()}));
                    }
                } else {
                    map.insert(qkey, (vector, indices));
                }
            }
        }
    }
    let bearing_anchored_overlap = bearing.keys().filter(|k| anchored.contains_key(*k)).count();
    let zero_unknown_in_known_only = zero_unknown_targets
        .iter()
        .filter(|k| known_only_targets.contains(*k))
        .count();
    let out = json!({
        "status":"exact_formal_collision_census",
        "field_prime":curve.p.to_string(),"subgroup_order":curve.n.to_string(),
        "all_unordered_triples":bearing_triples+known_only_triples+zero_unknown_triples,
        "bearing_triples":bearing_triples,"known_only_triples":known_only_triples,
        "zero_unknown_triples":zero_unknown_triples,
        "bearing_targets":bearing.len(),"known_only_targets":known_only_targets.len(),
        "anchored_targets":anchored.len(),
        "repeated_bearing_triples":repeated_bearing,
        "repeated_anchored_triples":repeated_anchored,
        "bearing_targets_with_distinct_formal_vectors":distinct_formal_bearing.len(),
        "anchored_targets_with_distinct_formal_vectors":distinct_formal_anchored.len(),
        "bearing_anchored_target_overlap":bearing_anchored_overlap,
        "zero_unknown_distinct_targets":zero_unknown_targets.len(),
        "zero_unknown_targets_in_anchored_support":zero_unknown_in_known_only,
        "bearing_distinct_formal_examples":bearing_examples,
        "anchored_distinct_formal_examples":anchored_examples,
        "setup_ms":setup_ms,"pair_index_ms":pair_index_ms,
        "triple_census_ms":start.elapsed().as_millis()-setup_ms-pair_index_ms,
        "input_sha256":[hash(&base_bytes),hash(&guide_bytes)],
        "source_sha256":hash(include_bytes!("formal_collision.rs")),
        "ec_source_sha256":hash(include_bytes!("../ec.rs"))
    });
    fs::write(
        &args[3],
        [serde_json::to_vec_pretty(&out)?, b"\n".to_vec()].concat(),
    )?;
    println!(
        "FORMAL_COLLISION_PASS bearing={} anchored={} distinct_formal_bearing={} anchor_overlap={}",
        bearing.len(),
        anchored.len(),
        distinct_formal_bearing.len(),
        bearing_anchored_overlap
    );
    Ok(())
}
