//! Exact two- and four-sum membership test for frozen adaptive Type I points.
#![recursion_limit = "256"]
#[path = "../ec.rs"]
mod ec;

use ec::{curve_zoo, Point};
use num_bigint::BigUint;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
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
fn parse_point(s: &str, curve: &curve_zoo::CurveParams) -> Point {
    if s == "O" {
        return Point::Infinity;
    }
    let (x, y) = s.split_once(':').unwrap();
    Point::Affine {
        x: curve.fe(BigUint::parse_bytes(x.as_bytes(), 10).unwrap()),
        y: curve.fe(BigUint::parse_bytes(y.as_bytes(), 10).unwrap()),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err(
            "usage: adaptive_short_sums BASE.json GUIDE.json ADAPTIVE.json OUT.json".into(),
        );
    }
    let start = Instant::now();
    let raw: Vec<Vec<u8>> = args[1..4].iter().map(fs::read).collect::<Result<_, _>>()?;
    let doc: Vec<Value> = raw
        .iter()
        .map(|b| serde_json::from_slice(b))
        .collect::<Result<_, _>>()?;
    let (base, guide, adaptive) = (&doc[0], &doc[1], &doc[2]);
    let curve = curve_zoo::CurveParams::p192();
    let a = curve.a_fe();
    let g = curve.generator();
    assert_eq!(dec(&base["field_prime"]), curve.p);
    assert_eq!(dec(&adaptive["field_prime"]), curve.p);
    assert_eq!(dec(&adaptive["subgroup_order"]), curve.n);
    assert_eq!(adaptive["input_sha256"][0], hash(&raw[0]));
    assert_eq!(adaptive["input_sha256"][2], hash(&raw[1]));
    let unknown: Vec<Point> = base["points"]
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
    assert_eq!(unknown.len(), 16);
    let mut old = unknown.clone();
    old.extend(unknown.iter().map(Point::neg));
    let guided: Vec<Point> = guide["accepted_scalars_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| g.scalar_mul(&dec(v), &a).double_vartime(&a).neg())
        .collect();
    assert_eq!(guided.len(), 64);
    old.extend(guided.iter().cloned());
    old.extend(guided.iter().map(Point::neg));
    assert_eq!(old.len(), 160);
    let mut pairs = Vec::with_capacity(160 * 161 / 2);
    let mut index = BTreeMap::<String, (usize, usize)>::new();
    for i in 0..160 {
        for j in i..160 {
            let sum = old[i].add_vartime(&old[j], &a);
            index.entry(key(&sum)).or_insert((i, j));
            pairs.push((i, j, sum));
        }
    }
    let setup_ms = start.elapsed().as_millis();
    let mut two_hits = Vec::new();
    let mut four_hits = Vec::new();
    let candidates = adaptive["candidates"].as_array().unwrap();
    assert_eq!(candidates.len(), 64);
    for (k, rec) in candidates.iter().enumerate() {
        assert_eq!(rec["index"].as_u64().unwrap() as usize, k);
        let c = parse_point(rec["new_point_key"].as_str().unwrap(), &curve);
        assert!(curve.is_on_curve(&c));
        assert_eq!(c.scalar_mul(&curve.n, &a), Point::Infinity);
        if let Some(&(i, j)) = index.get(&key(&c)) {
            two_hits.push(json!({"candidate_index":k,"old_indices":[i,j]}));
        }
        for (i, j, sum) in &pairs {
            let complement = c.add_vartime(&sum.neg(), &a);
            if let Some(&(u, v)) = index.get(&key(&complement)) {
                four_hits.push(json!({"candidate_index":k,"old_indices":[i,j,u,v]}));
                break;
            }
        }
    }
    let out = json!({
        "status":"exact_adaptive_short_sum_audit",
        "field_prime":curve.p.to_string(),"subgroup_order":curve.n.to_string(),
        "old_signed_points":160,"old_unordered_pairs_with_repetition":pairs.len(),
        "distinct_pair_sum_points":index.len(),
        "selected_new_points":candidates.len(),
        "two_sum_membership_tests":candidates.len(),
        "four_sum_complement_tests_upper_bound":candidates.len()*pairs.len(),
        "two_sum_hits":two_hits,"four_sum_hits":four_hits,
        "input_sha256":raw.iter().map(|b| hash(b)).collect::<Vec<_>>(),
        "source_sha256":hash(include_bytes!("adaptive_short_sums.rs")),
        "ec_source_sha256":hash(include_bytes!("../ec.rs")),
        "setup_ms":setup_ms,"search_ms":start.elapsed().as_millis()-setup_ms
    });
    fs::write(
        &args[4],
        [serde_json::to_vec_pretty(&out)?, b"\n".to_vec()].concat(),
    )?;
    println!(
        "ADAPTIVE_SHORT_SUMS_PASS selected={} two_hits={} four_hits={}",
        candidates.len(),
        out["two_sum_hits"].as_array().unwrap().len(),
        out["four_sum_hits"].as_array().unwrap().len()
    );
    Ok(())
}
