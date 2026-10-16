//! Frozen different-j inverse-Velu factor-base comparison.
//! Solver modes read only the public input. Scalars exist only in checker input.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use std::sync::OnceLock;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex};
use crypto_lib::cryptanalysis::ext_curve::{
    take_field_counters, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::ic_boundary::{rho_reference_negation, CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::isogeny_walk::record::{ec1_identity, V};
use crypto_lib::cryptanalysis::residual_walk::pow_mod;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;
static CELL: OnceLock<(u64, u64, u64)> = OnceLock::new();
fn set_cell(p: u64, r: u64, h: u64) {
    assert!([
        (409, 433, 387),
        (421, 457, 387),
        (457, 421, 495),
        (1237, 1201, 1275)
    ]
    .contains(&(p, r, h)));
    assert!(CELL.set((p, r, h)).is_ok());
}
fn p() -> u64 {
    CELL.get().unwrap().0
}
fn r() -> u64 {
    CELL.get().unwrap().1
}
fn h() -> u64 {
    CELL.get().unwrap().2
}
fn rel_cap() -> u64 {
    8 * r()
}
fn rho_cap() -> u64 {
    20 * r()
}
fn baby_size() -> u64 {
    (r() as f64).sqrt().ceil() as u64
}
fn j_pair() -> (u64, u64) {
    match p() {
        409 => (84, 368),
        421 => (55, 176),
        457 => (119, 57),
        1237 => (1066, 727),
        _ => unreachable!(),
    }
}
const N: usize = 64;
const COLUMNS: usize = 8;
const ARMS: [&str; 5] = [
    "inverse-velu",
    "full-image",
    "source-bsgs",
    "rho-negation",
    "rho-negation-aa",
];

fn source() -> Curve {
    let f = Fp2::new(p()).unwrap();
    Curve {
        f,
        a: f.from_fp(p() - 21),
        b: f.from_fp(37),
    }
}

fn kernel(src: &Curve) -> Point {
    Point::affine(src.f.from_fp(3), src.f.one())
}

fn destination(src: &Curve, t: Point) -> Curve {
    let f = &src.f;
    let a = f.sub(f.neg(f.scale(src.a, 9)), f.scale(f.sqr(t.x), 30));
    let x3 = f.mul(f.sqr(t.x), t.x);
    let b = f.sub(
        f.sub(f.neg(f.scale(src.b, 27)), f.scale(x3, 70)),
        f.scale(f.mul(src.a, t.x), 42),
    );
    let dst = Curve { f: *f, a, b };
    assert_eq!((dst.a, dst.b), ([p() - 81, 0], [p() - 243, 0]));
    dst
}

fn velu(src: &Curve, t: Point, p: Point, ops: &mut GroupOps) -> Point {
    if p.infinity || p == t || p == src.neg(t) {
        return src.identity();
    }
    let plus = src.add(ops, p, t);
    let minus = src.add(ops, p, src.neg(t));
    assert!(!plus.infinity && !minus.infinity);
    let f = &src.f;
    Point::affine(
        f.sub(f.add(f.add(p.x, plus.x), minus.x), f.scale(t.x, 2)),
        f.add(f.add(p.y, plus.y), minus.y),
    )
}

fn parse_point(v: &Value) -> Point {
    if v["infinity"].as_bool().unwrap() {
        return Point::infinity([0, 0]);
    }
    Point::affine(
        [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
    )
}

fn points(v: &Value) -> Vec<Point> {
    v.as_array().unwrap().iter().map(parse_point).collect()
}

fn tuple(p: Point) -> String {
    format!(
        "{},{},{},{},{}",
        u64::from(p.infinity),
        p.x[0],
        p.x[1],
        p.y[0],
        p.y[1]
    )
}

fn identities(curve: &Curve, g: Point, j: u64) -> Value {
    let f = &curve.f;
    let c0 = p() - f.nu;
    let field = format!(
        "fpk-{}-2-{}",
        p(),
        &sha256_hex(format!("fpk-modulus:{}:{c0},0", p()).as_bytes())[..8]
    );
    let model = canonical(&json!({
        "a":[curve.a[0].to_string(),curve.a[1].to_string()],
        "b":[curve.b[0].to_string(),curve.b[1].to_string()],
        "field":field,"form":"y^2=x^3+a*x+b","k":"2",
        "modulus":[c0.to_string(),"0"],"p":p().to_string(),"v":"1"
    }))
    .unwrap();
    let model_sha256 = sha256_hex(model.as_bytes());
    let n2 = h() * r();
    let trace = (p() * p() + 1) as i64 - n2 as i64;
    let icv1 = format!(
        "ICV1:{field}:{trace}:{n2}:{j},0:unk:unk:r:{}",
        &model_sha256[..12]
    );
    let coords = |x: [u64; 2]| V::Seq(vec![V::int(x[0]), V::int(x[1])]);
    let field_record = V::map(vec![
        ("characteristic", V::int(p())),
        ("degree", V::int(2)),
        ("representation", V::s("polynomial basis")),
        ("modulus", V::Seq(vec![V::int(c0), V::int(0), V::int(1)])),
    ]);
    let curve_record = V::map(vec![
        ("model", V::s("short Weierstrass")),
        ("a", coords(curve.a)),
        ("b", coords(curve.b)),
        ("subgroup_order", V::s(r().to_string())),
        ("cofactor", V::s(h().to_string())),
        ("generator", V::Seq(vec![coords(g.x), coords(g.y)])),
        ("target_group", V::s("prime-order subgroup")),
    ]);
    let (ec1, uid) = ec1_identity(&BigUint::from(p()), &field_record, &curve_record, "fpk2");
    json!({"icv1":icv1,"model":model,"model_sha256":model_sha256,
        "ec1":ec1,"ec1_uid":uid,"ec1_field_record":field_record.canonical_json(),
        "ec1_curve_record":curve_record.canonical_json(),
        "ec1_status":"self-contained extension record, not checked against registry"})
}

fn scan_base(curve: &Curve, ops: &mut GroupOps) -> ([u64; 2], Vec<Point>, Vec<Point>, u64) {
    let f = &curve.f;
    let mut scanned = 0u64;
    for a in 0..p() {
        let line = [a, 1];
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        for t in 1..p() {
            let mut lifts = curve.lift_x(f.scale(line, t));
            lifts.sort_by_key(|point| f.pack(point.y));
            for point in lifts {
                scanned += 1;
                let q = curve.mul(ops, point, h());
                if q.infinity {
                    continue;
                }
                let orbit = curve.key(&q).min(curve.key(&curve.neg(q)));
                if seen.insert(orbit) {
                    raw.push(point);
                    projected.push(q);
                    if raw.len() == COLUMNS {
                        break;
                    }
                }
            }
            if raw.len() == COLUMNS {
                break;
            }
        }
        if raw.len() == COLUMNS {
            return (line, raw, projected, scanned);
        }
    }
    panic!("eight projected columns absent")
}

type Poly = Vec<[u64; 2]>;

fn trim_poly(f: &Fp2, a: &mut Poly) {
    while a.last().is_some_and(|&x| f.is_zero(x)) {
        a.pop();
    }
}

fn sub_poly(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![f.zero(); a.len().max(b.len())];
    for (i, &v) in a.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    for (i, &v) in b.iter().enumerate() {
        out[i] = f.sub(out[i], v);
    }
    trim_poly(f, &mut out);
    out
}

fn div_rem_poly(f: &Fp2, a: &Poly, b: &Poly) -> (Poly, Poly) {
    assert!(!b.is_empty());
    let mut rem = a.clone();
    let mut quo = vec![f.zero(); a.len().saturating_sub(b.len()) + 1];
    let inv_lead = f.inv(*b.last().unwrap()).unwrap();
    while rem.len() >= b.len() {
        let shift = rem.len() - b.len();
        let coeff = f.mul(*rem.last().unwrap(), inv_lead);
        quo[shift] = coeff;
        for (j, &term) in b.iter().enumerate() {
            rem[shift + j] = f.sub(rem[shift + j], f.mul(coeff, term));
        }
        trim_poly(f, &mut rem);
    }
    trim_poly(f, &mut quo);
    (quo, rem)
}

fn mul_mod_poly(f: &Fp2, a: &Poly, b: &Poly, modulus: &Poly) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut product = vec![f.zero(); a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            product[i + j] = f.add(product[i + j], f.mul(x, y));
        }
    }
    div_rem_poly(f, &product, modulus).1
}

fn pow_mod_poly(f: &Fp2, a: &Poly, mut exponent: u64, modulus: &Poly) -> Poly {
    let mut out = vec![f.one()];
    let mut base = div_rem_poly(f, a, modulus).1;
    while exponent > 0 {
        if exponent & 1 != 0 {
            out = mul_mod_poly(f, &out, &base, modulus);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = mul_mod_poly(f, &base, &base, modulus);
        }
    }
    out
}

fn gcd_poly(f: &Fp2, mut a: Poly, mut b: Poly) -> Poly {
    while !b.is_empty() {
        let rem = div_rem_poly(f, &a, &b).1;
        a = b;
        b = rem;
    }
    if let Some(&lead) = a.last() {
        let inverse = f.inv(lead).unwrap();
        for coeff in &mut a {
            *coeff = f.mul(*coeff, inverse);
        }
    }
    a
}

fn inverse_velu_poly(f: &Fp2, u: [u64; 2]) -> Poly {
    vec![
        f.neg(f.add(f.scale(u, 9), f.from_fp(32))),
        f.add(f.scale(u, 6), f.from_fp(21)),
        f.neg(f.add(u, f.from_fp(6))),
        f.one(),
    ]
}

fn inverse_velu_root(f: &Fp2, u: [u64; 2], split_attempts: &mut u64) -> Option<[u64; 2]> {
    let cubic = inverse_velu_poly(f, u);
    let x = vec![f.zero(), f.one()];
    let rational = gcd_poly(
        f,
        cubic.clone(),
        sub_poly(f, &pow_mod_poly(f, &x, f.order(), &cubic), &x),
    );
    if rational.len() <= 1 {
        return None;
    }
    let linear = if rational.len() == 2 {
        rational
    } else {
        let mut factor = None;
        for shift in 0..32u64 {
            *split_attempts += 1;
            let shifted_x = vec![f.from_fp(shift), f.one()];
            let character = pow_mod_poly(f, &shifted_x, (f.order() - 1) / 2, &rational);
            for test in [sub_poly(f, &character, &vec![f.one()]), character] {
                let divisor = gcd_poly(f, rational.clone(), test);
                if divisor.len() > 1 && divisor.len() < rational.len() {
                    let linear = if divisor.len() == 2 {
                        divisor
                    } else {
                        let (quotient, remainder) = div_rem_poly(f, &rational, &divisor);
                        assert!(remainder.is_empty());
                        quotient
                    };
                    assert_eq!(linear.len(), 2);
                    factor = Some(linear);
                    break;
                }
            }
            if factor.is_some() {
                break;
            }
        }
        factor.expect("inverse-Velu polynomial split cap")
    };
    let root = f.neg(f.mul(linear[0], f.inv(linear[1]).unwrap()));
    assert!(poly_eval(f, &cubic, root) == f.zero());
    Some(root)
}

fn poly_eval(f: &Fp2, poly: &Poly, x: [u64; 2]) -> [u64; 2] {
    poly.iter()
        .rev()
        .fold(f.zero(), |acc, &term| f.add(f.mul(acc, x), term))
}

#[derive(Default)]
struct OracleStats {
    destination_points: u64,
    polynomial_calls: u64,
    split_attempts: u64,
    rational_x_roots: u64,
    rejected_no_root: u64,
    rejected_no_lift: u64,
    exact_map_checks: u64,
    rejected_map_mismatch: u64,
    accepted_image_points: u64,
}

fn inverse_velu_base(
    src: &Curve,
    dst: &Curve,
    t: Point,
    ops: &mut GroupOps,
) -> (Vec<Point>, Vec<Point>, Vec<Point>, OracleStats) {
    let f = &src.f;
    let mut raw = Vec::new();
    let mut projected = Vec::new();
    let mut witnesses = Vec::new();
    let mut seen = BTreeSet::new();
    let mut stats = OracleStats::default();
    for packed_x in 0..f.order() {
        let x = [packed_x % p(), packed_x / p()];
        let mut lifts = dst.lift_x(x);
        lifts.sort_by_key(|point| dst.key(point));
        for mapped in lifts {
            stats.destination_points += 1;
            stats.polynomial_calls += 1;
            let Some(root) = inverse_velu_root(f, x, &mut stats.split_attempts) else {
                stats.rejected_no_root += 1;
                continue;
            };
            stats.rational_x_roots += 1;
            let source_lifts = src.lift_x(root);
            if source_lifts.is_empty() {
                stats.rejected_no_lift += 1;
                continue;
            }
            let mut found = None;
            for point in source_lifts {
                stats.exact_map_checks += 1;
                if velu(src, t, point, ops) == mapped {
                    found = Some(point);
                    break;
                }
            }
            let Some(point) = found else {
                stats.rejected_map_mismatch += 1;
                continue;
            };
            stats.accepted_image_points += 1;
            let plus = src.add(ops, point, t);
            let minus = src.add(ops, point, src.neg(t));
            let mut preimages = [point, plus, minus];
            assert!(preimages.iter().all(|&q| !q.infinity && velu(src, t, q, ops) == mapped));
            preimages.sort_by_key(|q| (q.x[0], q.x[1], q.y[0], q.y[1]));
            let q = dst.mul(ops, mapped, h());
            if q.infinity {
                continue;
            }
            let orbit = dst.key(&q).min(dst.key(&dst.neg(q)));
            if seen.insert(orbit) {
                raw.push(mapped);
                projected.push(q);
                witnesses.push(preimages[0]);
                if raw.len() == COLUMNS {
                    return (raw, projected, witnesses, stats);
                }
            }
        }
    }
    panic!("inverse-Velu oracle did not find eight columns")
}

fn full_image_base(
    src: &Curve,
    dst: &Curve,
    t: Point,
    ops: &mut GroupOps,
) -> (Vec<Point>, Vec<Point>, Vec<Point>, u64, u64, u64) {
    let f = &src.f;
    let mut squares = BTreeMap::<[u64; 2], Vec<[u64; 2]>>::new();
    for y0 in 0..p() {
        for y1 in 0..p() {
            let y = [y0, y1];
            squares.entry(f.sqr(y)).or_default().push(y);
        }
    }
    let mut image = BTreeMap::<u64, (Point, Point, u64)>::new();
    let mut source_points = 0u64;
    let mut visit = |point: Point, ops: &mut GroupOps| {
        let mapped = velu(src, t, point, ops);
        assert!(dst.is_on_curve(mapped));
        let slot = image.entry(dst.key(&mapped)).or_insert((mapped, point, 0));
        assert_eq!(slot.0, mapped);
        slot.2 += 1;
        source_points += 1;
    };
    visit(src.identity(), ops);
    for x0 in 0..p() {
        for x1 in 0..p() {
            let x = [x0, x1];
            let rhs = f.add(f.add(f.mul(f.sqr(x), x), f.mul(src.a, x)), src.b);
            if let Some(ys) = squares.get(&rhs) {
                for &y in ys {
                    visit(Point::affine(x, y), ops);
                }
            }
        }
    }
    assert_eq!(source_points, h() * r());
    assert_eq!(image.len() as u64, h() * r() / 3);
    assert!(image.values().all(|slot| slot.2 == 3));
    let mut raw = Vec::new();
    let mut projected = Vec::new();
    let mut witness = Vec::new();
    let mut seen = BTreeSet::new();
    let mut candidate_count = 0u64;
    for (_, (mapped, preimage, _)) in image {
        if mapped.infinity {
            continue;
        }
        candidate_count += 1;
        let q = dst.mul(ops, mapped, h());
        if q.infinity {
            continue;
        }
        let orbit = dst.key(&q).min(dst.key(&dst.neg(q)));
        if seen.insert(orbit) {
            raw.push(mapped);
            projected.push(q);
            witness.push(preimage);
            if raw.len() == COLUMNS {
                break;
            }
        }
    }
    assert_eq!(raw.len(), COLUMNS);
    (
        raw,
        projected,
        witness,
        source_points,
        h() * r() / 3,
        candidate_count,
    )
}

fn assert_public(v: &Value) {
    fn walk(v: &Value) {
        match v {
            Value::Object(map) => {
                for (k, child) in map {
                    assert!(
                        !["answer_scalar", "scalars", "answers", "private_answers"]
                            .contains(&k.as_str()),
                        "answer-bearing public input"
                    );
                    walk(child);
                }
            }
            Value::Array(items) => {
                for item in items {
                    walk(item);
                }
            }
            _ => {}
        }
    }
    walk(v);
    assert_eq!(v["schema"], "different_j_trend_batch_public/v1");
    assert_eq!((v["p"].as_u64(), v["r"].as_u64()), (Some(p()), Some(r())));
    assert_eq!(v["targets"].as_array().unwrap().len(), N);
}

fn freeze(fiber_path: &str, public_path: &str, checker_path: &str) {
    assert!(
        !Path::new(public_path).exists() && !Path::new(checker_path).exists(),
        "refuse overwrite"
    );
    let fiber_bytes = fs::read(fiber_path).unwrap();
    let fiber: Value = serde_json::from_slice(&fiber_bytes).unwrap();
    assert_eq!(fiber["schema"], "different_j_trend_fiber/v1");
    assert_eq!(fiber["all_checks_passed"], true);
    set_cell(
        fiber["p"].as_u64().unwrap(),
        fiber["r"].as_u64().unwrap(),
        fiber["h"].as_u64().unwrap(),
    );
    let src = source();
    let t = kernel(&src);
    let dst = destination(&src, t);
    let g = parse_point(&fiber["source_generator"]);
    let mapped_g = parse_point(&fiber["image_generator"]);
    let mut ops = GroupOps::default();
    assert!(src.is_on_curve(g) && dst.is_on_curve(mapped_g));
    assert_eq!(velu(&src, t, g, &mut ops), mapped_g);
    assert!(src.mul(&mut ops, g, r()).infinity && dst.mul(&mut ops, mapped_g, r()).infinity);
    let (source_line, source_raw, source_base, source_scan) = scan_base(&src, &mut ops);
    let (dest_line, dest_raw, dest_base, dest_scan) = scan_base(&dst, &mut ops);
    let (full_raw, full_base, full_witness, full_scan, image_points, full_candidates) =
        full_image_base(&src, &dst, t, &mut ops);
    let mut entries = Vec::new();
    for scalar in 1..r() {
        let point = src.mul(&mut ops, g, scalar);
        let mapped = dst.mul(&mut ops, mapped_g, scalar);
        assert!(!point.infinity && !mapped.infinity);
        assert_eq!(velu(&src, t, point, &mut ops), mapped);
        let key = tuple(point);
        let digest = sha256_hex(format!("different-j-prime-to-char-target-v2:13:{key}").as_bytes());
        entries.push((digest, key, point, mapped, scalar));
    }
    entries.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    let targets: Vec<Value> = entries
        .iter()
        .take(N)
        .map(|e| json!({"source":e.2,"destination":e.3}))
        .collect();
    let scalars: Vec<u64> = entries.iter().take(N).map(|e| e.4).collect();
    let (source_j, destination_j) = j_pair();
    let public = json!({
        "schema":"different_j_trend_batch_public/v1","p":p(),"r":r(),"cofactor":h(),
        "field":format!("F_{}[s]/(s^2-{})",p(),src.f.nu),"source_a":src.a,"source_b":src.b,
        "destination_a":dst.a,"destination_b":dst.b,"source_j":source_j,"destination_j":destination_j,
        "source_identity":identities(&src,g,source_j),"destination_identity":identities(&dst,mapped_g,destination_j),
        "fiber_receipt_sha256":sha256_hex(&fiber_bytes),"kernel":t,
        "generator":g,"mapped_generator":mapped_g,
        "source_line":source_line,"source_base_raw":source_raw,"source_base_projected":source_base,
        "source_scan_count":source_scan,
        "destination_line":dest_line,"destination_base_raw":dest_raw,
        "destination_base_projected":dest_base,"destination_scan_count":dest_scan,
        "full_base_raw":full_raw,"full_base_projected":full_base,
        "full_base_witness_source":full_witness,"full_source_points":full_scan,
        "full_image_points":image_points,"full_candidate_count":full_candidates,
        "target_count":N,"target_order_domain":"different-j-prime-to-char-target-v2:13",
        "targets":targets
    });
    assert_public(&public);
    fs::write(
        public_path,
        format!("{}\n", serde_json::to_string_pretty(&public).unwrap()),
    )
    .unwrap();
    let public_hash = sha256_hex(&fs::read(public_path).unwrap());
    let checker = json!({"schema":"different_j_trend_batch_checker/v1",
        "public_sha256":public_hash,"target_count":N,"scalars":scalars});
    fs::write(
        checker_path,
        format!("{}\n", serde_json::to_string_pretty(&checker).unwrap()),
    )
    .unwrap();
    println!(
        "{}",
        json!({"schema":"different_j_trend_batch_freeze/v1",
        "public_sha256":public_hash,"checker_sha256":sha256_hex(&fs::read(checker_path).unwrap()),
        "target_count":N,"source_scan":source_scan,"destination_scan":dest_scan,
        "full_source_points":full_scan,"full_image_points":image_points,"full_candidates":full_candidates,"group_ops":ops})
    );
}

fn modmul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % r() as u128) as u64
}

fn insert_relation(pivots: &mut [Option<Vec<u64>>; 9], mut row: Vec<u64>) -> bool {
    for col in 0..9 {
        if row[col] == 0 {
            continue;
        }
        if let Some(pivot) = &pivots[col] {
            let scale = row[col];
            for j in col..10 {
                row[j] = (row[j] + r() - modmul(scale, pivot[j])) % r();
            }
        } else {
            let inv = pow_mod(row[col], r() - 2, r());
            for item in row.iter_mut().skip(col) {
                *item = modmul(*item, inv);
            }
            pivots[col] = Some(row);
            return true;
        }
    }
    assert_eq!(row[9], 0, "contradictory relation");
    false
}

fn recover(pivots: &[Option<Vec<u64>>; 9]) -> Option<[u64; 9]> {
    if pivots.iter().any(Option::is_none) {
        return None;
    }
    let mut x = [0u64; 9];
    for i in (0..9).rev() {
        let row = pivots[i].as_ref().unwrap();
        let mut rhs = row[9];
        for j in i + 1..9 {
            rhs = (rhs + r() - modmul(row[j], x[j])) % r();
        }
        x[i] = rhs;
    }
    Some(x)
}

struct RelationSetup {
    signed: Vec<Point>,
    table_points: Vec<Point>,
    table: BTreeMap<u64, Vec<(usize, usize)>>,
    pair_entries: u64,
    scan_count: u64,
    image_points: u64,
    candidate_count: u64,
    oracle: Option<OracleStats>,
}

fn prepare_relation(
    public: &Value,
    src: &Curve,
    dst: &Curve,
    t: Point,
    arm: &str,
    ops: &mut GroupOps,
) -> RelationSetup {
    let (raw, base, scan_count, image_points, candidate_count, oracle) = match arm {
        "inverse-velu" => {
            let (raw, base, witnesses, stats) = inverse_velu_base(src, dst, t, ops);
            assert_eq!(raw, points(&public["full_base_raw"]));
            assert_eq!(base, points(&public["full_base_projected"]));
            assert_eq!(witnesses, points(&public["full_base_witness_source"]));
            assert_eq!(
                stats.accepted_image_points,
                public["full_candidate_count"].as_u64().unwrap()
            );
            let scan = stats.destination_points;
            let candidates = stats.accepted_image_points;
            (raw, base, scan, candidates, candidates, Some(stats))
        }
        "full-image" => {
            let (raw, base, witnesses, scan, images, candidates) =
                full_image_base(src, dst, t, ops);
            assert_eq!(raw, points(&public["full_base_raw"]));
            assert_eq!(base, points(&public["full_base_projected"]));
            assert_eq!(witnesses, points(&public["full_base_witness_source"]));
            (raw, base, scan, images, candidates, None)
        }
        "destination-direct" => {
            let (line, raw, base, scan) = scan_base(dst, ops);
            assert_eq!(public["destination_line"], json!(line));
            assert_eq!(raw, points(&public["destination_base_raw"]));
            assert_eq!(base, points(&public["destination_base_projected"]));
            (raw, base, scan, 0, 0, None)
        }
        "source" | "isogeny-image" => {
            let (line, raw, base, scan) = scan_base(src, ops);
            assert_eq!(public["source_line"], json!(line));
            assert_eq!(raw, points(&public["source_base_raw"]));
            assert_eq!(base, points(&public["source_base_projected"]));
            (raw, base, scan, 0, 0, None)
        }
        _ => panic!("unknown relation arm"),
    };
    assert_eq!(raw.len(), COLUMNS);
    let sign_curve = if arm == "source" || arm == "isogeny-image" {
        src
    } else {
        dst
    };
    let signed: Vec<Point> = base.iter().flat_map(|&p| [p, sign_curve.neg(p)]).collect();
    let table_curve = if arm == "source" { src } else { dst };
    let table_points: Vec<Point> = if arm == "isogeny-image" {
        signed.iter().map(|&p| velu(src, t, p, ops)).collect()
    } else {
        signed.clone()
    };
    let mut table = BTreeMap::<u64, Vec<(usize, usize)>>::new();
    let mut pair_entries = 0u64;
    for i in 0..table_points.len() {
        for j in i..table_points.len() {
            let sum = table_curve.add(ops, table_points[i], table_points[j]);
            if !sum.infinity {
                table.entry(table_curve.key(&sum)).or_default().push((i, j));
                pair_entries += 1;
            }
        }
    }
    RelationSetup {
        signed,
        table_points,
        table,
        pair_entries,
        scan_count,
        image_points,
        candidate_count,
        oracle,
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_relation(
    src: &Curve,
    dst: &Curve,
    t: Point,
    g: Point,
    q: Point,
    mapped_g: Point,
    mapped_q: Point,
    arm: &str,
    setup: &RelationSetup,
    seed: u64,
) -> Value {
    let direct = ["destination-direct", "full-image", "inverse-velu"].contains(&arm);
    let query_curve = if direct { dst } else { src };
    let table_curve = if arm == "source" { src } else { dst };
    let query_g = if direct { mapped_g } else { g };
    let query_q = if direct { mapped_q } else { q };
    let mut rng = StdRng::seed_from_u64(seed);
    let mut pivots: [Option<Vec<u64>>; 9] = std::array::from_fn(|_| None);
    let mut witnesses = Vec::new();
    let mut trials = 0u64;
    let mut relations = 0u64;
    let mut rank = 0u64;
    let mut invalid = 0u64;
    let mut map_queries = 0u64;
    let mut la_cpu_ns = 0u64;
    let mut ops = GroupOps::default();
    for _ in 0..rel_cap() {
        trials += 1;
        let a = rng.gen_range(0..r());
        let b = rng.gen_range(0..r());
        let a_g = query_curve.mul(&mut ops, query_g, a);
        let b_q = query_curve.mul(&mut ops, query_q, b);
        let candidate = query_curve.add(&mut ops, a_g, b_q);
        let projected = query_curve.mul(&mut ops, candidate, h());
        let searched = if arm == "isogeny-image" {
            map_queries += 1;
            velu(src, t, projected, &mut ops)
        } else {
            projected
        };
        if searched.infinity {
            continue;
        }
        let Some(pairs) = setup.table.get(&table_curve.key(&searched)) else {
            continue;
        };
        let &(i, j) = pairs.first().unwrap();
        let sum = table_curve.add(&mut ops, setup.table_points[i], setup.table_points[j]);
        if sum != searched {
            invalid += 1;
            continue;
        }
        if arm == "isogeny-image"
            && src.add(&mut ops, setup.signed[i], setup.signed[j]) != projected
        {
            invalid += 1;
            continue;
        }
        let mut row = vec![0u64; 10];
        for index in [i, j] {
            let col = index / 2;
            let coeff = if index % 2 == 0 { 1 } else { r() - 1 };
            row[col] = (row[col] + coeff) % r();
        }
        row[8] = (r() - modmul(h() % r(), b)) % r();
        row[9] = modmul(h() % r(), a);
        witnesses.push(json!({"a":a,"b":b,"i":i,"j":j}));
        relations += 1;
        let la_start = thread_cpu_ns();
        if insert_relation(&mut pivots, row) {
            rank += 1;
        }
        la_cpu_ns += thread_cpu_ns() - la_start;
        if rank == 9 {
            break;
        }
    }
    let recovered = recover(&pivots).map(|x| x[8]);
    let mut verify_ops = GroupOps::default();
    let verified = recovered.is_some_and(|k| {
        let source_ok = src.mul(&mut verify_ops, g, k) == q;
        source_ok && (arm == "source" || dst.mul(&mut verify_ops, mapped_g, k) == mapped_q)
    });
    json!({"recovered":recovered,"verified":verified,"seed":seed,"trials":trials,
        "relations":relations,"rank":rank,"failed_queries":trials-relations,
        "invalid_candidates":invalid,"map_queries":map_queries,"la_cpu_ns_included":la_cpu_ns,
        "relation_witnesses":witnesses,"query_group_ops":ops,"verify_group_ops":verify_ops})
}

#[allow(clippy::too_many_arguments)]
fn solve_bsgs(
    src: &Curve,
    g: Point,
    q: Point,
    baby: &BTreeMap<u64, u64>,
    giant_step: Point,
) -> Value {
    let mut ops = GroupOps::default();
    let mut current = q;
    let mut recovered = None;
    let mut lookups = 0u64;
    for giant in 0..r().div_ceil(baby_size()) {
        lookups += 1;
        if let Some(&small) = baby.get(&src.key(&current)) {
            recovered = Some((giant * baby_size() + small) % r());
            break;
        }
        current = src.add(&mut ops, current, giant_step);
    }
    let mut verify_ops = GroupOps::default();
    let verified = recovered.is_some_and(|k| src.mul(&mut verify_ops, g, k) == q);
    json!({"recovered":recovered,"verified":verified,"giant_lookups":lookups,"query_group_ops":ops,"verify_group_ops":verify_ops})
}

#[allow(clippy::too_many_arguments)]
fn solve_rho(src: &Curve, g: Point, q: Point, seed: u64) -> Value {
    let mut attempts = Vec::new();
    let mut recovered = None;
    let mut rho_ops = GroupOps::default();
    for retry in 0..32u64 {
        let try_seed = seed.wrapping_add(retry * 0x0100_0000);
        let result = rho_reference_negation(src, g, q, r(), try_seed, rho_cap());
        rho_ops.merge(result.group_ops);
        attempts.push(json!({"retry":retry,"seed":try_seed,"steps":result.steps,
            "walks":result.walks,"group_ops":result.group_ops,"recovered":result.recovered,"verified":result.verified}));
        if result.verified {
            recovered = result.recovered;
            break;
        }
    }
    let mut verify_ops = GroupOps::default();
    let verified = recovered.is_some_and(|k| src.mul(&mut verify_ops, g, k) == q);
    json!({"recovered":recovered,"verified":verified,"seed":seed,"attempts":attempts,
        "rho_group_ops":rho_ops,"verify_group_ops":verify_ops})
}

fn thread_cpu_ns() -> u64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    assert_eq!(
        unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts) },
        0
    );
    ts.tv_sec as u64 * 1_000_000_000 + ts.tv_nsec as u64
}

fn peak_rss_bytes() -> u64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    assert_eq!(
        unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) },
        0
    );
    let raw = unsafe { usage.assume_init() }.ru_maxrss as u64;
    if cfg!(target_os = "macos") {
        raw
    } else {
        raw * 1024
    }
}

fn run(path: &str, round: u64, arm: &str, prefix_length: usize) {
    assert!(round < 3 && ARMS.contains(&arm));
    assert!([1, 4, 16, 64].contains(&prefix_length));
    let transfer = ["isogeny-image", "destination-direct", "full-image", "inverse-velu"].contains(&arm);
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    take_field_counters();
    let public: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    set_cell(
        public["p"].as_u64().unwrap(),
        public["r"].as_u64().unwrap(),
        public["cofactor"].as_u64().unwrap(),
    );
    assert_public(&public);
    let src = source();
    let t = kernel(&src);
    let dst = destination(&src, t);
    let g = parse_point(&public["generator"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let mut common_ops = GroupOps::default();
    assert!(src.is_on_curve(g));
    assert!(src.mul(&mut common_ops, g, r()).infinity);
    if transfer {
        assert!(dst.is_on_curve(mapped_g));
        assert!(dst.mul(&mut common_ops, mapped_g, r()).infinity);
        assert_eq!(velu(&src, t, g, &mut common_ops), mapped_g);
    }
    let mut source_keys = BTreeSet::new();
    let mut dest_keys = BTreeSet::new();
    let targets: Vec<(Point, Point)> = public["targets"]
        .as_array()
        .unwrap()
        .iter()
        .take(prefix_length)
        .map(|row| {
            let q = parse_point(&row["source"]);
            let mq = parse_point(&row["destination"]);
            assert!(!q.infinity && src.is_on_curve(q));
            assert!(src.mul(&mut common_ops, q, r()).infinity);
            assert!(source_keys.insert(src.key(&q)));
            if transfer {
                assert!(!mq.infinity && dst.is_on_curve(mq));
                assert!(dst.mul(&mut common_ops, mq, r()).infinity);
                assert_eq!(velu(&src, t, q, &mut common_ops), mq);
                assert!(dest_keys.insert(dst.key(&mq)));
            }
            (q, mq)
        })
        .collect();
    assert_eq!(targets.len(), prefix_length);
    assert_eq!(source_keys.len(), prefix_length);
    if transfer {
        assert_eq!(dest_keys.len(), prefix_length);
    }
    let common_cpu_ns = thread_cpu_ns() - cpu;
    let common_fp = take_field_counters();
    let setup_start = thread_cpu_ns();
    let mut setup_ops = GroupOps::default();
    let relation = if [
        "source",
        "isogeny-image",
        "destination-direct",
        "full-image",
        "inverse-velu",
    ]
    .contains(&arm)
    {
        Some(prepare_relation(
            &public,
            &src,
            &dst,
            t,
            arm,
            &mut setup_ops,
        ))
    } else {
        None
    };
    let mut baby = BTreeMap::new();
    let mut giant_step = src.identity();
    if arm == "source-bsgs" {
        let mut current = src.identity();
        for j in 0..baby_size() {
            assert!(baby.insert(src.key(&current), j).is_none());
            current = src.add(&mut setup_ops, current, g);
        }
        giant_step = src.neg(src.mul(&mut setup_ops, g, baby_size()));
    }
    let setup_cpu_ns = thread_cpu_ns() - setup_start;
    let setup_fp = take_field_counters();
    let mut total_fp_muls = common_fp.0 + setup_fp.0;
    let mut total_fp_invs = common_fp.1 + setup_fp.1;
    let mut target_records = Vec::new();
    let mut prefixes = Vec::new();
    let mut all_verified = true;
    for (index, &(q, mq)) in targets.iter().enumerate() {
        let target_start = thread_cpu_ns();
        let seed = (0x2026_1015u64 ^ p())
            .wrapping_add(round.wrapping_mul(0x9e37_79b9))
            .wrapping_add((index as u64).wrapping_mul(0x0123_4567));
        let mut result = if let Some(setup) = &relation {
            solve_relation(&src, &dst, t, g, q, mapped_g, mq, arm, setup, seed)
        } else if arm == "source-bsgs" {
            solve_bsgs(&src, g, q, &baby, giant_step)
        } else {
            solve_rho(&src, g, q, seed)
        };
        let fp = take_field_counters();
        total_fp_muls += fp.0;
        total_fp_invs += fp.1;
        all_verified &= result["verified"] == true;
        result["position"] = json!(index);
        result["target_cpu_ns"] = json!(thread_cpu_ns() - target_start);
        result["target_fp_muls"] = json!(fp.0);
        result["target_fp_invs"] = json!(fp.1);
        target_records.push(result);
        let count = index + 1;
        if [1, 4, 16, 64].contains(&count) {
            prefixes.push(json!({"targets":count,"total_cpu_ns":thread_cpu_ns()-cpu,
                "total_wall_ns":wall.elapsed().as_nanos() as u64,
                "total_fp_muls":total_fp_muls,"total_fp_invs":total_fp_invs,
                "peak_rss_bytes":peak_rss_bytes(),
                "verified":target_records.iter().all(|row| row["verified"]==true)}));
        }
    }
    let total_cpu_ns = thread_cpu_ns() - cpu;
    let total_wall_ns = wall.elapsed().as_nanos() as u64;
    let peak_rss = peak_rss_bytes();
    let relation_meta = relation.as_ref().map(|setup| {
        json!({
            "pair_entries":setup.pair_entries,"pair_keys":setup.table.len(),
            "support_keys":setup.table.keys().copied().collect::<Vec<_>>(),
            "scan_count":setup.scan_count,"image_points":setup.image_points,
            "candidate_count":setup.candidate_count,
            "oracle":setup.oracle.as_ref().map(|stats| json!({
                "destination_points":stats.destination_points,
                "polynomial_calls":stats.polynomial_calls,
                "split_attempts":stats.split_attempts,
                "rational_x_roots":stats.rational_x_roots,
                "rejected_no_root":stats.rejected_no_root,
                "rejected_no_lift":stats.rejected_no_lift,
                "exact_map_checks":stats.exact_map_checks,
                "rejected_map_mismatch":stats.rejected_map_mismatch,
                "accepted_image_points":stats.accepted_image_points,
            }))})
    });
    println!(
        "{}",
        json!({"schema":"different_j_trend_batch_run/v1","p":p(),"r":r(),"h":h(),
        "relation_cap":rel_cap(),"rho_step_cap":rho_cap(),"round":round,"arm":arm,
        "prefix_length":prefix_length,"target_count":targets.len(),"public_sha256":sha256_hex(&fs::read(path).unwrap()),
        "common_cpu_ns":common_cpu_ns,"common_fp_muls":common_fp.0,"common_fp_invs":common_fp.1,
        "common_group_ops":common_ops,"setup_cpu_ns":setup_cpu_ns,"setup_fp_muls":setup_fp.0,
        "setup_fp_invs":setup_fp.1,"setup_group_ops":setup_ops,"relation_setup":relation_meta,
        "baby_entries":baby.len(),"prefixes":prefixes,"target_records":target_records,
        "construction":"inverse-velu-v1",
        "total_cpu_ns":total_cpu_ns,"total_wall_ns":total_wall_ns,
        "total_fp_muls":total_fp_muls,"total_fp_invs":total_fp_invs,
        "peak_rss_bytes":peak_rss,"all_verified":all_verified})
    );
    assert!(all_verified, "unrecovered or invalid target");
    assert!(peak_rss <= 512 * 1024 * 1024, "RSS ceiling");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("freeze") if args.len()==5 => freeze(&args[2],&args[3],&args[4]),
        Some("run" | "batch-run") if args.len()==6 => run(&args[2],args[3].parse().unwrap(),&args[4],args[5].parse().unwrap()),
        _ => panic!("usage: different_j_trend_batch freeze FIBER.json PUBLIC.json CHECKER.json | run PUBLIC.json ROUND ARM PREFIX")
    }
}
