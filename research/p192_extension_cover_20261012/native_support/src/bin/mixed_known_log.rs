//! Exact P-192 Type I mixed-base support and two known-scalar query laws.
//! Geometry of all 160 fixed-origin classes is checked separately by GP.
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
fn row(indices: [usize; 3]) -> [i8; 16] {
    let mut out = [0; 16];
    for i in indices {
        if i < 32 {
            out[i % 16] += if i < 16 { 1 } else { -1 };
        }
    }
    out
}
fn nonzero(v: &[i8; 16]) -> bool {
    v.iter().any(|x| *x != 0)
}
fn known_sum(indices: [usize; 3], scalars: &[BigUint], r: &BigUint) -> BigUint {
    indices
        .iter()
        .fold(BigUint::zero(), |z, i| (z + &scalars[*i]) % r)
}
fn candidate_scalar(counter: u64, r: &BigUint) -> Option<BigUint> {
    let input = format!("p192-type-i-mixed-independent-v1:{counter}");
    let digest = Sha256::digest(input.as_bytes());
    let n = BigUint::from_bytes_be(&digest[..24]);
    if n.is_zero() || n >= *r {
        None
    } else {
        Some(n)
    }
}
fn audit(base: &Value, density: &Value, collector: &Value) -> Value {
    let start = Instant::now();
    let curve = curve_zoo::CurveParams::p192();
    let p = &curve.p;
    let r = &curve.n;
    let a = curve.a_fe();
    let g = curve.generator();
    assert_eq!(dec(&base["field_prime"]), *p);
    assert_eq!(dec(&density["p"]), *p);
    assert_eq!(dec(&density["r"]), *r);
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
            assert_eq!(q.scalar_mul(r, &a), Point::Infinity);
            q
        })
        .collect();
    assert_eq!(unsigned.len(), 16);
    let minus_half = (r - BigUint::one()) >> 1usize;
    let mut points = unsigned.clone();
    points.extend(unsigned.iter().map(Point::neg));
    let mut preimages: Vec<Point> = points
        .iter()
        .map(|q| q.scalar_mul(&minus_half, &a))
        .collect();
    let mut scalars = vec![BigUint::zero(); 32];
    let guide: Vec<BigUint> = collector["accepted_scalars_decimal"]
        .as_array()
        .unwrap()
        .iter()
        .map(dec)
        .collect();
    assert_eq!(guide.len(), 64);
    for s in &guide {
        assert!(!s.is_zero() && s < r);
        let pre = g.scalar_mul(s, &a);
        assert_eq!(cell(&pre, &coeff, p, &a), "++");
        let q = pre.double_vartime(&a).neg();
        points.push(q.clone());
        preimages.push(pre);
        scalars.push(s.clone());
    }
    for s in &guide {
        let pre = g.scalar_mul(s, &a);
        points.push(pre.double_vartime(&a));
        preimages.push(pre.neg());
        scalars.push(r - s);
    }
    assert_eq!(points.len(), 160);
    for i in 0..160 {
        assert_eq!(preimages[i].double_vartime(&a).neg(), points[i]);
    }
    assert_eq!(points.iter().map(key).collect::<BTreeSet<_>>().len(), 160);
    for i in 32..160 {
        assert_eq!(
            g.scalar_mul(&scalars[i], &a).double_vartime(&a).neg(),
            points[i]
        );
    }
    let setup_ms = start.elapsed().as_millis();
    let pairs: Vec<Vec<Point>> = (0..160)
        .map(|i| {
            (i..160)
                .map(|j| points[i].add_vartime(&points[j], &a))
                .collect()
        })
        .collect();
    let pre_pairs: Vec<Vec<Point>> = (0..160)
        .map(|i| {
            (i..160)
                .map(|j| preimages[i].add_vartime(&preimages[j], &a))
                .collect()
        })
        .collect();
    let index_ms = start.elapsed().as_millis() - setup_ms;

    // Keep one nonzero equation per target, while counting every triple.
    let mut bearing: BTreeMap<String, (Point, [i8; 16], BigUint)> = BTreeMap::new();
    let mut known_only: BTreeMap<String, (Point, BigUint)> = BTreeMap::new();
    let mut bearing_triples = 0usize;
    let mut zero_row_triples = 0usize;
    let mut known_triples = 0usize;
    for i in 0..160 {
        for j in i..160 {
            for k in j..160 {
                let ix = [i, j, k];
                let q = pairs[i][j - i].add_vartime(&points[k], &a);
                let pre = pre_pairs[i][j - i].add_vartime(&preimages[k], &a);
                let s = known_sum(ix, &scalars, r);
                if i >= 32 {
                    known_triples += 1;
                    known_only.entry(key(&q)).or_insert((pre, s));
                } else {
                    let v = row(ix);
                    if nonzero(&v) {
                        bearing_triples += 1;
                        bearing.entry(key(&q)).or_insert((pre, v, s));
                    } else {
                        zero_row_triples += 1;
                    }
                }
            }
        }
    }
    assert_eq!(bearing_triples + zero_row_triples + known_triples, 695520);
    assert_eq!(known_triples, 357760);
    let enumeration_ms = start.elapsed().as_millis() - setup_ms - index_ms;

    let mut bearing_cells = BTreeMap::<String, usize>::new();
    let mut accepted_bearing = BTreeSet::new();
    for (k, (pre, _, _)) in &bearing {
        let c = cell(pre, &coeff, p, &a);
        *bearing_cells.entry(c.clone()).or_default() += 1;
        if c == "++" {
            accepted_bearing.insert(k.clone());
        }
    }
    let mut known_cells = BTreeMap::<String, usize>::new();
    let mut steered = Vec::new();
    let mut overlap = 0usize;
    for (k, (pre, scalar)) in &known_only {
        let c = cell(pre, &coeff, p, &a);
        *known_cells.entry(c.clone()).or_default() += 1;
        if bearing.contains_key(k) {
            overlap += 1;
        }
        if c == "++" && !scalar.is_zero() && steered.len() < 64 {
            let witness = bearing.get(k);
            steered.push(json!({"scalar":scalar.to_string(),"point_key":k,"rank_bearing_hit":witness.is_some(),
                "row":witness.map(|(_,v,_)|v.to_vec())}));
        }
    }
    assert_eq!(steered.len(), 64);
    let census_ms = start.elapsed().as_millis() - setup_ms - index_ms - enumeration_ms;

    let frozen: BTreeSet<String> = guide.iter().map(ToString::to_string).collect();
    let mut sampled = Vec::new();
    let mut attempts = 0u64;
    let mut character_rejects = 0u64;
    while sampled.len() < 64 && attempts < 4096 {
        let Some(s) = candidate_scalar(attempts, r) else {
            attempts += 1;
            continue;
        };
        attempts += 1;
        if frozen.contains(&s.to_string()) {
            continue;
        }
        let pre = g.scalar_mul(&s, &a);
        if cell(&pre, &coeff, p, &a) != "++" {
            character_rejects += 1;
            continue;
        }
        let q = pre.double_vartime(&a).neg();
        let witness = bearing.get(&key(&q));
        sampled.push(
            json!({"scalar":s.to_string(),"rank_bearing_hit":witness.is_some(),
            "row":witness.map(|(_,v,_)|v.to_vec())}),
        );
    }
    assert_eq!(sampled.len(), 64);
    let query_ms = start.elapsed().as_millis() - setup_ms - index_ms - enumeration_ms - census_ms;
    let steered_hits = steered
        .iter()
        .filter(|v| v["rank_bearing_hit"] == true)
        .count();
    let independent_hits = sampled
        .iter()
        .filter(|v| v["rank_bearing_hit"] == true)
        .count();
    let lower = dec(&density["accepted_nonzero_scalar_lower_bound"]);
    json!({
        "status":"exact_mixed_known_log_type_i_support_census",
        "field_prime":p.to_string(),"subgroup_order":r.to_string(),
        "unknown_unsigned_classes":16,"unknown_signed_classes":32,
        "known_unsigned_classes":64,"known_signed_classes":128,"logarithm_columns":16,
        "all_unordered_triples":695520,"known_only_triples":known_triples,
        "zero_unknown_row_triples":zero_row_triples,"rank_bearing_triples":bearing_triples,
        "distinct_known_only_targets":known_only.len(),"distinct_rank_bearing_targets":bearing.len(),
        "known_and_bearing_target_overlap":overlap,
        "known_character_cells":known_cells,"bearing_character_cells":bearing_cells,
        "accepted_rank_bearing_support":accepted_bearing.len(),
        "accepted_nonzero_scalar_lower_bound":lower.to_string(),
        "uniform_accepted_hit_upper_bound":format!("{}/{}",accepted_bearing.len(),lower),
        "steered_known_only_queries":steered,"steered_rank_bearing_hits":steered_hits,
        "independent_accepted_queries":sampled,"independent_rank_bearing_hits":independent_hits,
        "independent_attempts":attempts,"independent_character_rejects":character_rejects,
        "accepted_bearing_keys_sha256":hash(accepted_bearing.iter().cloned().collect::<Vec<_>>().join("\n").as_bytes()),
        "setup_ms":setup_ms,"pair_index_ms":index_ms,"triple_enumeration_ms":enumeration_ms,
        "character_census_ms":census_ms,"fresh_queries_ms":query_ms,
        "query_laws":["first 64 lexicographically indexed distinct ++ known-only three-sums with known scalar", "first 64 SHA256-counter nonzero scalars below r, disjoint from guide, passing ++"],
        "geometry_receipts":["type_i_fixed_origin_reduction.json", "type_i_negative_fixed_origin_candidates.gp", "type_i_guided_negative_fixed_origin_candidates.gp"]
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err("usage: mixed_known_log BASE.json DENSITY.json COLLECTOR.json OUT.json".into());
    }
    let inputs: Vec<Vec<u8>> = args[1..4].iter().map(fs::read).collect::<Result<_, _>>()?;
    let values: Vec<Value> = inputs
        .iter()
        .map(|s| serde_json::from_slice(s))
        .collect::<Result<_, _>>()?;
    let mut out = audit(&values[0], &values[1], &values[2]);
    out["input_sha256"] = json!(inputs.iter().map(|v| hash(v)).collect::<Vec<_>>());
    out["source_sha256"] = json!(hash(include_bytes!("mixed_known_log.rs")));
    out["ec_source_sha256"] = json!(hash(include_bytes!("../ec.rs")));
    fs::write(
        &args[4],
        [serde_json::to_vec_pretty(&out)?, b"\n".to_vec()].concat(),
    )?;
    println!(
        "MIXED_KNOWN_LOG_PASS bearing={} accepted={} steered_hits={} independent_hits={}",
        out["distinct_rank_bearing_targets"],
        out["accepted_rank_bearing_support"],
        out["steered_rank_bearing_hits"],
        out["independent_rank_bearing_hits"]
    );
    Ok(())
}
