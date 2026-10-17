//! Native screen and frozen instance for an explicit different-j degree-three map over F_(p^2).
//! The selector and stop conditions are frozen in the dated study protocol.

#![recursion_limit = "256"]

use std::cell::Cell;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex};
use crypto_lib::cryptanalysis::ext_curve::{
    jacobi, random_point, take_field_counters, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{rho_reference_negation, CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::isogeny_walk::record::{ec1_identity, V};
use crypto_lib::cryptanalysis::residual_walk::{is_prime_u64, pow_mod};
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;

fn source(p: u64) -> Curve {
    let f = Fp2::new(p).unwrap();
    Curve {
        f,
        a: f.from_fp(p - 21 % p),
        b: f.from_fp(37),
    }
}

fn destination(src: &Curve, t: Point) -> Curve {
    let f = &src.f;
    let a = f.sub(f.neg(f.scale(src.a, 9)), f.scale(f.sqr(t.x), 30));
    let x3 = f.mul(f.sqr(t.x), t.x);
    let b = f.sub(
        f.sub(f.neg(f.scale(src.b, 27)), f.scale(x3, 70)),
        f.scale(f.mul(src.a, t.x), 42),
    );
    Curve { f: *f, a, b }
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

fn count_fp(curve: &Curve) -> u64 {
    let p = curve.f.p;
    let mut n = 1u64;
    for x in 0..p {
        let xx = curve.f.from_fp(x);
        let x3 = curve.f.mul(curve.f.sqr(xx), xx);
        let rhs = curve
            .f
            .add(curve.f.add(x3, curve.f.mul(curve.a, xx)), curve.b)[0];
        n += if rhs == 0 {
            1
        } else if jacobi(rhs, p) == 1 {
            2
        } else {
            0
        };
    }
    n
}

fn j(curve: &Curve) -> Option<[u64; 2]> {
    let f = &curve.f;
    let a3 = f.mul(f.sqr(curve.a), curve.a);
    let four_a3 = f.scale(a3, 4);
    let den = f.add(four_a3, f.scale(f.sqr(curve.b), 27));
    f.inv(den).map(|inv| f.mul(f.scale(four_a3, 1728), inv))
}

fn screen() {
    for p in (61..=271).filter(|p| p % 3 == 1 && is_prime_u64(*p)) {
        let src = source(p);
        let t = Point::affine(src.f.from_fp(3), src.f.one());
        let mut ops = GroupOps::default();
        let torsion = src.is_on_curve(t) && src.mul(&mut ops, t, 3).infinity;
        if !torsion {
            println!(
                "{}",
                json!({"p":p,"status":"reject","reason":"kernel_check"})
            );
            continue;
        }
        let dst = destination(&src, t);
        let js = j(&src);
        let jd = j(&dst);
        if js.is_none() || jd.is_none() {
            println!("{}", json!({"p":p,"status":"reject","reason":"singular"}));
            continue;
        }
        if js == jd {
            println!(
                "{}",
                json!({"p":p,"status":"reject","reason":"equal_j","j":js})
            );
            continue;
        }
        let n1 = count_fp(&src);
        let n1_dst = count_fp(&dst);
        if n1 != n1_dst {
            println!(
                "{}",
                json!({"p":p,"status":"reject","reason":"base_order_mismatch","n1":n1,"n1_dst":n1_dst})
            );
            continue;
        }
        let trace = p as i64 + 1 - n1 as i64;
        let n2 = ((p + 1) * (p + 1)) as i64 - trace * trace;
        let n2 = n2 as u64;
        let factors = factor_u64(n2);
        let r = factors
            .iter()
            .rev()
            .find(|(prime, power)| *prime >= 19 && *prime != 3 && *power == 1)
            .map(|x| x.0);
        let Some(r) = r else {
            println!(
                "{}",
                json!({"p":p,"status":"reject","reason":"subgroup","n1":n1,"n2":n2,"factors":factors})
            );
            continue;
        };
        let mut rng = StdRng::seed_from_u64(0x2026_1015 ^ p);
        let mut sample_checks = 0u64;
        let mut good = true;
        for _ in 0..64 {
            let a = random_point(&src, &mut rng);
            let b = random_point(&src, &mut rng);
            let ma = velu(&src, t, a, &mut ops);
            let mb = velu(&src, t, b, &mut ops);
            let mab = velu(&src, t, src.add(&mut ops, a, b), &mut ops);
            if !dst.is_on_curve(ma) || !dst.is_on_curve(mb) || mab != dst.add(&mut ops, ma, mb) {
                good = false;
                break;
            }
            sample_checks += 1;
        }
        if !good {
            println!(
                "{}",
                json!({"p":p,"status":"reject","reason":"map_check","sample_checks":sample_checks})
            );
            continue;
        }
        println!(
            "{}",
            json!({"schema":"different_j_screen/v1","p":p,"status":"accept","source":{"a":src.a,"b":src.b,"j":js},"destination":{"a":dst.a,"b":dst.b,"j":jd},"kernel":[t,src.neg(t)],"n1":n1,"n2":n2,"r":r,"cofactor":n2/r,"factors":factors,"sample_checks":sample_checks,"ops":ops})
        );
        break;
    }
}

fn identities(curve: &Curve, n2: u64, r: u64, g: Point) -> serde_json::Value {
    let f = &curve.f;
    let p = f.p;
    let c0 = p - f.nu;
    let field = format!(
        "fpk-{p}-2-{}",
        &sha256_hex(format!("fpk-modulus:{p}:{c0},0").as_bytes())[..8]
    );
    let model = canonical(&json!({
        "a":[curve.a[0].to_string(),curve.a[1].to_string()],
        "b":[curve.b[0].to_string(),curve.b[1].to_string()],
        "field":field,"form":"y^2=x^3+a*x+b","k":"2",
        "modulus":[c0.to_string(),"0"],"p":p.to_string(),"v":"1"
    }))
    .unwrap();
    let model_sha256 = sha256_hex(model.as_bytes());
    let jv = j(curve).unwrap();
    let trace = p as i64 * p as i64 + 1 - n2 as i64;
    let icv1 = format!(
        "ICV1:{field}:{trace}:{n2}:{},{}:unk:unk:r:{}",
        jv[0],
        jv[1],
        &model_sha256[..12]
    );
    let coords = |x: [u64; 2]| V::Seq(vec![V::int(x[0]), V::int(x[1])]);
    let field_record = V::map(vec![
        ("characteristic", V::int(p)),
        ("degree", V::int(2)),
        ("representation", V::s("polynomial basis")),
        ("modulus", V::Seq(vec![V::int(c0), V::int(0), V::int(1)])),
    ]);
    let curve_record = V::map(vec![
        ("model", V::s("short Weierstrass")),
        ("a", coords(curve.a)),
        ("b", coords(curve.b)),
        ("subgroup_order", V::s(r.to_string())),
        ("cofactor", V::s((n2 / r).to_string())),
        ("generator", V::Seq(vec![coords(g.x), coords(g.y)])),
        ("target_group", V::s("prime-order subgroup")),
    ]);
    let (ec1, uid) = ec1_identity(&BigUint::from(p), &field_record, &curve_record, "fpk2");
    json!({"icv1":icv1,"model":model,"model_sha256":model_sha256,
        "ec1":ec1,"ec1_uid":uid,"ec1_field_record":field_record.canonical_json(),
        "ec1_curve_record":curve_record.canonical_json(),
        "ec1_status":"self-contained extension record, not checked against registry"})
}

fn freeze() {
    let p = 67u64;
    let r = 61u64;
    let n2 = 4575u64;
    let h = n2 / r;
    let src = source(p);
    let f = &src.f;
    let t = Point::affine(f.from_fp(3), f.one());
    let dst = destination(&src, t);
    assert_eq!(count_fp(&src), 75);
    assert_eq!(count_fp(&dst), 75);
    assert_eq!(j(&src), Some([7, 0]));
    assert_eq!(j(&dst), Some([11, 0]));
    let mut ops = GroupOps::default();
    assert!(src.mul(&mut ops, t, 3).infinity);
    assert!(velu(&src, t, t, &mut ops).infinity);
    assert!(velu(&src, t, src.neg(t), &mut ops).infinity);
    let mut map_checks = 0u64;
    for x0 in 0..p {
        for x1 in 0..p {
            for point in src.lift_x([x0, x1]) {
                let mapped = velu(&src, t, point, &mut ops);
                assert!(dst.is_on_curve(mapped));
                map_checks += 1;
            }
        }
    }
    assert_eq!(map_checks, n2 - 1);
    let target_seed = format!("different-j/2026-10-15/p={p}/r={r}");
    let seed_bytes = sha256(target_seed.as_bytes());
    let mut rng = StdRng::from_seed(seed_bytes);
    let g = loop {
        let candidate = src.mul(&mut ops, random_point(&src, &mut rng), h);
        if !candidate.infinity {
            break candidate;
        }
    };
    assert!(src.mul(&mut ops, g, r).infinity);
    let k = 1 + u64::from_be_bytes(seed_bytes[0..8].try_into().unwrap()) % (r - 1);
    let q = src.mul(&mut ops, g, k);
    let mapped_g = velu(&src, t, g, &mut ops);
    let mapped_q = velu(&src, t, q, &mut ops);
    assert!(!mapped_g.infinity && !mapped_q.infinity);
    assert!(dst.mul(&mut ops, mapped_g, r).infinity);
    assert_eq!(dst.mul(&mut ops, mapped_g, k), mapped_q);
    let mut selected = None;
    for a in 0..p {
        let line = [a, 1];
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        let mut scanned = 0u64;
        for t0 in 1..p {
            let mut lifts = src.lift_x(f.scale(line, t0));
            lifts.sort_by_key(|point| f.pack(point.y));
            for point in lifts {
                scanned += 1;
                let projected_point = src.mul(&mut ops, point, h);
                if projected_point.infinity {
                    continue;
                }
                let orbit = src
                    .key(&projected_point)
                    .min(src.key(&src.neg(projected_point)));
                if seen.insert(orbit) {
                    raw.push(point);
                    projected.push(projected_point);
                    if raw.len() == 8 {
                        break;
                    }
                }
            }
            if raw.len() == 8 {
                break;
            }
        }
        if raw.len() == 8 {
            selected = Some((line, raw, projected, scanned));
            break;
        }
    }
    let (line, raw, projected, scanned) = selected.expect("eight projected columns absent");
    let mapped_raw: Vec<_> = raw
        .iter()
        .map(|&point| velu(&src, t, point, &mut ops))
        .collect();
    let mapped_projected: Vec<_> = projected
        .iter()
        .map(|&point| velu(&src, t, point, &mut ops))
        .collect();
    for (&a, &b) in mapped_raw.iter().zip(&mapped_projected) {
        assert_eq!(dst.mul(&mut ops, a, h), b);
    }
    println!(
        "{}",
        json!({"schema":"different_j_fixture/v1","p":p,"nu":f.nu,
        "r":r,"n1":75,"n2":n2,"cofactor":h,"kernel":[t,src.neg(t)],
        "source":identities(&src, n2, r, g),"destination":identities(&dst, n2, r, mapped_g),
        "source_a":src.a,"source_b":src.b,"destination_a":dst.a,"destination_b":dst.b,
        "source_j":j(&src),"destination_j":j(&dst),"generator":g,"target":q,
        "mapped_generator":mapped_g,"mapped_target":mapped_q,"answer_scalar":k,
        "target_seed":target_seed,"line":line,"base_raw":raw,"base_projected":projected,
        "mapped_base_raw":mapped_raw,"mapped_base_projected":mapped_projected,
        "raw_points_scanned":scanned,"map_points_checked":map_checks,"ops":ops})
    );
}

fn parse_point(v: &Value) -> Point {
    if v["infinity"].as_bool().unwrap() {
        return Point::infinity([0, 0]);
    }
    let pair = |name: &str| -> [u64; 2] {
        let a = v[name].as_array().unwrap();
        [a[0].as_u64().unwrap(), a[1].as_u64().unwrap()]
    };
    Point::affine(pair("x"), pair("y"))
}

fn points(v: &Value) -> Vec<Point> {
    v.as_array().unwrap().iter().map(parse_point).collect()
}

fn scan_base(curve: &Curve, h: u64, columns: usize) -> ([u64; 2], Vec<Point>, Vec<Point>, u64) {
    let f = &curve.f;
    let mut scanned = 0;
    for a in 0..f.p {
        let line = [a, 1];
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        let mut ops = GroupOps::default();
        for t in 1..f.p {
            let mut lifts = curve.lift_x(f.scale(line, t));
            lifts.sort_by_key(|point| f.pack(point.y));
            for point in lifts {
                scanned += 1;
                let q = curve.mul(&mut ops, point, h);
                if q.infinity {
                    continue;
                }
                let orbit = curve.key(&q).min(curve.key(&curve.neg(q)));
                if seen.insert(orbit) {
                    raw.push(point);
                    projected.push(q);
                    if raw.len() == columns {
                        break;
                    }
                }
            }
            if raw.len() == columns {
                break;
            }
        }
        if raw.len() == columns {
            return (line, raw, projected, scanned);
        }
    }
    panic!("eight distinct projected columns absent");
}

/// Enumerate the entire rational image before choosing any subgroup columns.
/// The first source preimage is retained as an independently checkable witness.
fn scan_full_image(
    src: &Curve,
    dst: &Curve,
    t: Point,
    h: u64,
    columns: usize,
    ops: &mut GroupOps,
) -> (Vec<Point>, Vec<Point>, Vec<Point>, u64, u64, u64) {
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
    for x1 in 0..src.f.p {
        for x0 in 0..src.f.p {
            let mut lifts = src.lift_x([x0, x1]);
            lifts.sort_by_key(|point| src.f.pack(point.y));
            for point in lifts {
                visit(point, ops);
            }
        }
    }
    assert_eq!(source_points, 4575);
    assert_eq!(image.len(), 1525);
    assert!(image.values().all(|entry| entry.2 == 3));
    let mut raw = Vec::new();
    let mut projected = Vec::new();
    let mut witness = Vec::new();
    let mut orbits = BTreeSet::new();
    let mut image_candidates = 0u64;
    for (_, (mapped, preimage, _)) in image {
        if mapped.infinity {
            continue;
        }
        image_candidates += 1;
        let q = dst.mul(ops, mapped, h);
        if q.infinity {
            continue;
        }
        let orbit = dst.key(&q).min(dst.key(&dst.neg(q)));
        if orbits.insert(orbit) {
            raw.push(mapped);
            projected.push(q);
            witness.push(preimage);
            if raw.len() == columns {
                break;
            }
        }
    }
    assert_eq!(raw.len(), columns);
    (
        raw,
        projected,
        witness,
        source_points,
        1525,
        image_candidates,
    )
}

fn freeze_public(fixture_path: &str, public_path: &str) {
    let mut v: Value = serde_json::from_str(&fs::read_to_string(fixture_path).unwrap()).unwrap();
    assert_eq!(v["schema"], "different_j_fixture/v1");
    let answer = v.as_object_mut().unwrap().remove("answer_scalar");
    assert!(answer.is_some());
    let src = source(67);
    let t = parse_point(&v["kernel"][0]);
    let dst = destination(&src, t);
    let (line, raw, projected, scanned) = scan_base(&dst, 75, 8);
    assert_ne!(j(&src), j(&dst));
    v["schema"] = json!("different_j_public/v1");
    v["direct_line"] = json!(line);
    v["direct_base_raw"] = json!(raw);
    v["direct_base_projected"] = json!(projected);
    v["direct_raw_points_scanned"] = json!(scanned);
    assert!(v.get("answer_scalar").is_none());
    fs::write(public_path, serde_json::to_vec_pretty(&v).unwrap()).unwrap();
    println!(
        "{}",
        json!({"schema":"different_j_public_freeze/v1",
        "fixture_sha256":sha256_hex(&fs::read(fixture_path).unwrap()),
        "public_sha256":sha256_hex(&fs::read(public_path).unwrap()),
        "direct_line":line,"direct_points_scanned":scanned,
        "direct_columns":projected.len()})
    );
}

fn freeze_full_image(public_v1_path: &str, public_v2_path: &str) {
    let mut v: Value = serde_json::from_str(&fs::read_to_string(public_v1_path).unwrap()).unwrap();
    assert_eq!(v["schema"], "different_j_public/v1");
    assert!(v.get("answer_scalar").is_none());
    let src = source(67);
    let t = parse_point(&v["kernel"][0]);
    let dst = destination(&src, t);
    let mut ops = GroupOps::default();
    let (raw, projected, witnesses, source_points, image_points, image_candidates) =
        scan_full_image(&src, &dst, t, 75, 8, &mut ops);
    v["schema"] = json!("different_j_public/v2");
    v["full_base_raw"] = json!(raw);
    v["full_base_projected"] = json!(projected);
    v["full_base_witness_source"] = json!(witnesses);
    v["full_image_source_points"] = json!(source_points);
    v["full_image_points"] = json!(image_points);
    v["full_image_selection_candidates"] = json!(image_candidates);
    fs::write(public_v2_path, serde_json::to_vec_pretty(&v).unwrap()).unwrap();
    println!(
        "{}",
        json!({"schema":"different_j_full_image_freeze/v1",
        "source_public_sha256":sha256_hex(&fs::read(public_v1_path).unwrap()),
        "full_image_public_sha256":sha256_hex(&fs::read(public_v2_path).unwrap()),
        "source_points":source_points,"image_points":image_points,
        "selection_candidates":image_candidates,"columns":projected.len(),"group_ops":ops})
    );
}

fn read_public(path: &str) -> Value {
    let v: Value = serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert!(
        ["different_j_public/v1", "different_j_public/v2"].contains(&v["schema"].as_str().unwrap())
    );
    assert!(v.get("answer_scalar").is_none());
    assert_eq!((v["p"].as_u64(), v["r"].as_u64()), (Some(67), Some(61)));
    v
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

fn modmul(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn insert_relation(pivots: &mut [Option<Vec<u64>>; 9], mut row: Vec<u64>, r: u64) -> bool {
    for col in 0..9 {
        if row[col] == 0 {
            continue;
        }
        if let Some(pivot) = &pivots[col] {
            let scale = row[col];
            for j in col..10 {
                row[j] = (row[j] + r - modmul(scale, pivot[j], r)) % r;
            }
        } else {
            let inv = pow_mod(row[col], r - 2, r);
            for item in row.iter_mut().skip(col) {
                *item = modmul(*item, inv, r);
            }
            pivots[col] = Some(row);
            return true;
        }
    }
    assert_eq!(row[9], 0, "contradictory relation");
    false
}

fn recover(pivots: &[Option<Vec<u64>>; 9], r: u64) -> Option<[u64; 9]> {
    if pivots.iter().any(Option::is_none) {
        return None;
    }
    let mut x = [0u64; 9];
    for i in (0..9).rev() {
        let row = pivots[i].as_ref().unwrap();
        let mut rhs = row[9];
        for j in i + 1..9 {
            rhs = (rhs + r - modmul(row[j], x[j], r)) % r;
        }
        x[i] = rhs;
    }
    Some(x)
}

fn signed_base(curve: &Curve, base: &[Point]) -> Vec<Point> {
    base.iter().flat_map(|&p| [p, curve.neg(p)]).collect()
}

fn relation_run(path: &str, round: u64, arm: &str) {
    assert!(round < 3);
    assert!([
        "source",
        "isogeny-image",
        "destination-direct",
        "full-image"
    ]
    .contains(&arm));
    let public = read_public(path);
    let is_v2 = public["schema"] == "different_j_public/v2";
    assert!(arm != "full-image" || is_v2);
    let src = source(67);
    let t = parse_point(&public["kernel"][0]);
    let dst = destination(&src, t);
    assert_eq!(dst.a, [53, 0]);
    assert_eq!(dst.b, [25, 0]);
    let g = parse_point(&public["generator"]);
    let q = parse_point(&public["target"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let mapped_q = parse_point(&public["mapped_target"]);
    let r = 61u64;
    let h = 75u64;
    let seed = (0x2026_1015u64 ^ 67).wrapping_add(round * 0x9e37_79b9);
    let cap = 20 * r;
    take_field_counters();
    let start_wall = Instant::now();
    let start_cpu = thread_cpu_ns();
    let setup_cpu_start = thread_cpu_ns();
    let setup_wall_start = Instant::now();
    let mut setup_ops = GroupOps::default();
    let (curve, raw, base, scan_count, image_points, selection_candidates) = if arm == "full-image"
    {
        let (raw, projected, witnesses, scanned, image_count, candidates) =
            scan_full_image(&src, &dst, t, h, 8, &mut setup_ops);
        assert_eq!(raw, points(&public["full_base_raw"]));
        assert_eq!(projected, points(&public["full_base_projected"]));
        assert_eq!(witnesses, points(&public["full_base_witness_source"]));
        assert_eq!(public["full_image_source_points"], scanned);
        assert_eq!(public["full_image_points"], image_count);
        assert_eq!(public["full_image_selection_candidates"], candidates);
        (&dst, raw, projected, scanned, image_count, candidates)
    } else if arm == "destination-direct" {
        let (line, raw, projected, scanned) = scan_base(&dst, h, 8);
        assert_eq!(public["direct_line"], json!(line));
        assert_eq!(raw, points(&public["direct_base_raw"]));
        assert_eq!(projected, points(&public["direct_base_projected"]));
        (&dst, raw, projected, scanned, 0, 0)
    } else {
        let (line, raw, projected, scanned) = scan_base(&src, h, 8);
        assert_eq!(public["line"], json!(line));
        assert_eq!(raw, points(&public["base_raw"]));
        assert_eq!(projected, points(&public["base_projected"]));
        (&src, raw, projected, scanned, 0, 0)
    };
    let signed = signed_base(curve, &base);
    let table_curve = if arm == "isogeny-image" { &dst } else { curve };
    let table_points: Vec<Point> = if arm == "isogeny-image" {
        signed
            .iter()
            .map(|&p| velu(&src, t, p, &mut setup_ops))
            .collect()
    } else {
        signed.clone()
    };
    let mut table = BTreeMap::<u64, Vec<(usize, usize)>>::new();
    let mut pair_entries = 0u64;
    for i in 0..table_points.len() {
        for j in i..table_points.len() {
            let sum = table_curve.add(&mut setup_ops, table_points[i], table_points[j]);
            if !sum.infinity {
                table.entry(table_curve.key(&sum)).or_default().push((i, j));
                pair_entries += 1;
            }
        }
    }
    let setup_fp = take_field_counters();
    let setup_cpu_ns = thread_cpu_ns() - setup_cpu_start;
    let setup_wall_ns = setup_wall_start.elapsed().as_nanos() as u64;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut pivots: [Option<Vec<u64>>; 9] = std::array::from_fn(|_| None);
    let mut trials = 0u64;
    let mut relations = 0u64;
    let mut relation_witnesses = Vec::<Value>::new();
    let mut rank = 0u64;
    let mut map_queries = 0u64;
    let mut invalid_candidates = 0u64;
    let mut query_ops = GroupOps::default();
    let mut la_cpu_ns = 0u64;
    take_field_counters();
    let collect_cpu_start = thread_cpu_ns();
    let collect_wall_start = Instant::now();
    for _ in 0..cap {
        trials += 1;
        let a = rand::Rng::gen_range(&mut rng, 0..r);
        let b = rand::Rng::gen_range(&mut rng, 0..r);
        let (query_curve, query_g, query_q) = if arm == "destination-direct" || arm == "full-image"
        {
            (&dst, mapped_g, mapped_q)
        } else {
            (&src, g, q)
        };
        let a_g = query_curve.mul(&mut query_ops, query_g, a);
        let b_q = query_curve.mul(&mut query_ops, query_q, b);
        let candidate = query_curve.add(&mut query_ops, a_g, b_q);
        let projected = query_curve.mul(&mut query_ops, candidate, h);
        let searched = if arm == "isogeny-image" {
            map_queries += 1;
            velu(&src, t, projected, &mut query_ops)
        } else {
            projected
        };
        if searched.infinity {
            continue;
        }
        let Some(pairs) = table.get(&table_curve.key(&searched)) else {
            continue;
        };
        let &(i, j) = pairs.first().unwrap();
        let actual = table_curve.add(&mut query_ops, table_points[i], table_points[j]);
        if actual != searched {
            invalid_candidates += 1;
            continue;
        }
        if arm == "isogeny-image" {
            let source_sum = src.add(&mut query_ops, signed[i], signed[j]);
            if source_sum != projected {
                invalid_candidates += 1;
                continue;
            }
        }
        let mut row = vec![0u64; 10];
        for index in [i, j] {
            let col = index / 2;
            let coeff = if index % 2 == 0 { 1 } else { r - 1 };
            row[col] = (row[col] + coeff) % r;
        }
        row[8] = (r - modmul(h % r, b, r)) % r;
        row[9] = modmul(h % r, a, r);
        relations += 1;
        if is_v2 {
            relation_witnesses.push(json!({"a":a,"b":b,"i":i,"j":j}));
        }
        let la_start = thread_cpu_ns();
        if insert_relation(&mut pivots, row, r) {
            rank += 1;
        }
        la_cpu_ns += thread_cpu_ns() - la_start;
        if rank == 9 {
            break;
        }
    }
    let collect_cpu_ns = thread_cpu_ns() - collect_cpu_start;
    let collect_wall_ns = collect_wall_start.elapsed().as_nanos() as u64;
    let collect_fp = take_field_counters();
    let verify_start = thread_cpu_ns();
    let solve_start = thread_cpu_ns();
    let solved = recover(&pivots, r);
    let solve_cpu_ns = thread_cpu_ns() - solve_start;
    let recovered = solved.map(|x| x[8]);
    let mut verify_ops = GroupOps::default();
    let verified = recovered.is_some_and(|k| {
        src.mul(&mut verify_ops, g, k) == q && dst.mul(&mut verify_ops, mapped_g, k) == mapped_q
    });
    let verify_cpu_ns = thread_cpu_ns() - verify_start - solve_cpu_ns;
    let total_cpu_ns = thread_cpu_ns() - start_cpu;
    let total_wall_ns = start_wall.elapsed().as_nanos() as u64;
    let final_fp = take_field_counters();
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    println!(
        "{}",
        json!({"schema":if is_v2 {"different_j_relation_run/v2"} else {"different_j_relation_run/v1"},"arm":arm,
        "p":67,"r":r,"round":round,"seed":seed,"trial_cap":cap,
        "source_icv1":public["source"]["icv1"],
        "destination_icv1":public["destination"]["icv1"],
        "raw_base":raw.len(),"points_scanned":scan_count,
        "full_image_points":image_points,"selection_candidates":selection_candidates,
        "pair_entries":pair_entries,"pair_keys":table.len(),
        "setup_cpu_ns":setup_cpu_ns,"setup_wall_ns":setup_wall_ns,
        "setup_fp_muls":setup_fp.0,"setup_fp_invs":setup_fp.1,
        "setup_group_ops":setup_ops,"trials":trials,"relations":relations,
        "rank":rank,"failed_queries":trials-relations,"map_queries":map_queries,
        "relation_witnesses":relation_witnesses,
        "invalid_candidates":invalid_candidates,"query_group_ops":query_ops,
        "collect_cpu_ns":collect_cpu_ns,"collect_wall_ns":collect_wall_ns,
        "collect_fp_muls":collect_fp.0,"collect_fp_invs":collect_fp.1,
        "la_cpu_ns_included_in_collection":la_cpu_ns,"solve_cpu_ns":solve_cpu_ns,
        "verify_cpu_ns":verify_cpu_ns,"verify_group_ops":verify_ops,
        "recovered":recovered,"verified":verified,
        "total_cpu_ns":total_cpu_ns,"total_wall_ns":total_wall_ns,
        "final_fp_muls":final_fp.0,"final_fp_invs":final_fp.1,
        "peak_rss_bytes":peak_rss_bytes()})
    );
}

fn rho_run(path: &str, round: u64, aa: bool) {
    assert!(round < 3);
    let public = read_public(path);
    let src = source(67);
    let g = parse_point(&public["generator"]);
    let q = parse_point(&public["target"]);
    let seed = (0x2026_1015u64 ^ 67)
        .wrapping_add(round * 0x9e37_79b9)
        .wrapping_add(if aa { 0x00aa_0000 } else { 0 });
    let cap = 20 * 61;
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    take_field_counters();
    let mut attempts = Vec::new();
    let mut answer = None;
    for retry in 0..32u64 {
        let try_seed = seed.wrapping_add(retry * 0x0100_0000);
        let attempt_cpu = thread_cpu_ns();
        let attempt_wall = Instant::now();
        take_field_counters();
        let result = rho_reference_negation(&src, g, q, 61, try_seed, cap);
        let fp = take_field_counters();
        attempts.push(json!({"retry":retry,"seed":try_seed,"rho":result,
            "cpu_ns":thread_cpu_ns()-attempt_cpu,
            "wall_ns":attempt_wall.elapsed().as_nanos() as u64,
            "fp_muls":fp.0,"fp_invs":fp.1}));
        if result.verified {
            answer = result.recovered;
            break;
        }
    }
    let total_fp = take_field_counters();
    println!(
        "{}",
        json!({"schema":"different_j_rho_run/v1",
        "arm":if aa {"rho-negation-aa"} else {"rho-negation"},
        "round":round,"seed":seed,"step_cap":cap,"attempts":attempts,
        "recovered":answer,"verified":answer.is_some(),
        "total_cpu_ns":thread_cpu_ns()-cpu,
        "total_wall_ns":wall.elapsed().as_nanos() as u64,
        "final_fp_muls":total_fp.0,"final_fp_invs":total_fp.1,
        "peak_rss_bytes":peak_rss_bytes()})
    );
}

// Independent replay arithmetic: separate field and point-add formulas from ExtCurve.
struct Replay {
    p: u64,
    nu: u64,
    muls: Cell<u64>,
    inversions: Cell<u64>,
    group_adds: Cell<u64>,
}

impl Replay {
    fn new() -> Self {
        Self {
            p: 67,
            nu: 2,
            muls: Cell::new(0),
            inversions: Cell::new(0),
            group_adds: Cell::new(0),
        }
    }
    fn m(&self, a: u64, b: u64) -> u64 {
        self.muls.set(self.muls.get() + 1);
        (a * b) % self.p
    }
    fn addm(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }
    fn subm(&self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }
    fn add(&self, a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
        [self.addm(a[0], b[0]), self.addm(a[1], b[1])]
    }
    fn sub(&self, a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
        [self.subm(a[0], b[0]), self.subm(a[1], b[1])]
    }
    fn mul(&self, a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
        [
            self.addm(self.m(a[0], b[0]), self.m(self.nu, self.m(a[1], b[1]))),
            self.addm(self.m(a[0], b[1]), self.m(a[1], b[0])),
        ]
    }
    fn pow_fp(&self, mut a: u64, mut e: u64) -> u64 {
        let mut out = 1;
        while e > 0 {
            if e & 1 == 1 {
                out = self.m(out, a);
            }
            e >>= 1;
            if e > 0 {
                a = self.m(a, a);
            }
        }
        out
    }
    fn inv(&self, a: [u64; 2]) -> [u64; 2] {
        let norm = self.subm(self.m(a[0], a[0]), self.m(self.nu, self.m(a[1], a[1])));
        assert_ne!(norm, 0);
        self.inversions.set(self.inversions.get() + 1);
        let ni = self.pow_fp(norm, self.p - 2);
        let result = [self.m(a[0], ni), self.m(self.subm(0, a[1]), ni)];
        assert_eq!(self.mul(a, result), [1, 0]);
        result
    }
    fn on_curve(&self, point: Point, a: u64, b: u64) -> bool {
        point.infinity
            || self.mul(point.y, point.y)
                == self.add(
                    self.add(
                        self.mul(self.mul(point.x, point.x), point.x),
                        self.mul([a, 0], point.x),
                    ),
                    [b, 0],
                )
    }
    fn add_points(&self, lhs: Point, rhs: Point, a: u64) -> Point {
        self.group_adds.set(self.group_adds.get() + 1);
        if lhs.infinity {
            return rhs;
        }
        if rhs.infinity {
            return lhs;
        }
        if lhs.x == rhs.x && self.add(lhs.y, rhs.y) == [0, 0] {
            return Point::infinity([0, 0]);
        }
        let slope = if lhs.x == rhs.x {
            self.mul(
                self.add(self.mul([3, 0], self.mul(lhs.x, lhs.x)), [a, 0]),
                self.inv(self.mul([2, 0], lhs.y)),
            )
        } else {
            self.mul(self.sub(rhs.y, lhs.y), self.inv(self.sub(rhs.x, lhs.x)))
        };
        let x = self.sub(self.sub(self.mul(slope, slope), lhs.x), rhs.x);
        let y = self.sub(self.mul(slope, self.sub(lhs.x, x)), lhs.y);
        Point::affine(x, y)
    }
    fn scalar_mul(&self, point: Point, mut k: u64, a: u64) -> Point {
        let mut out = Point::infinity([0, 0]);
        let mut base = point;
        while k > 0 {
            if k & 1 == 1 {
                out = self.add_points(out, base, a);
            }
            k >>= 1;
            if k > 0 {
                base = self.add_points(base, base, a);
            }
        }
        out
    }
    fn neg_point(&self, point: Point) -> Point {
        if point.infinity {
            point
        } else {
            Point::affine(point.x, self.sub([0, 0], point.y))
        }
    }
    fn key(&self, point: Point) -> u64 {
        if point.infinity {
            return 0;
        }
        let pack = |x: [u64; 2]| x[0] + self.p * x[1];
        let sign = u64::from(pack(point.y) > pack(self.neg_point(point).y));
        ((pack(point.x) + 1) << 1) | sign
    }
    fn velu(&self, point: Point, t: Point) -> Point {
        let neg_t = self.neg_point(t);
        if point.infinity || point == t || point == neg_t {
            return Point::infinity([0, 0]);
        }
        let plus = self.add_points(point, t, 46);
        let minus = self.add_points(point, neg_t, 46);
        assert!(!plus.infinity && !minus.infinity);
        Point::affine(
            self.sub(
                self.add(self.add(point.x, plus.x), minus.x),
                self.mul([2, 0], t.x),
            ),
            self.add(self.add(point.y, plus.y), minus.y),
        )
    }
}

fn replay_relation_row(
    checker: &Replay,
    public: &Value,
    row: &Value,
    g: Point,
    q: Point,
    mapped_g: Point,
    mapped_q: Point,
    t: Point,
) -> u64 {
    let arm = row["arm"].as_str().unwrap();
    let direct = arm == "destination-direct" || arm == "full-image";
    let (base, curve_a) = match arm {
        "source" | "isogeny-image" => (points(&public["base_projected"]), 46),
        "destination-direct" => (points(&public["direct_base_projected"]), 53),
        "full-image" => (points(&public["full_base_projected"]), 53),
        _ => panic!("unexpected relation arm"),
    };
    assert_eq!(base.len(), 8);
    let mut sign_orbits = BTreeSet::new();
    for &point in &base {
        assert!(!point.infinity);
        assert!(checker.scalar_mul(point, 61, curve_a).infinity);
        let orbit = checker
            .key(point)
            .min(checker.key(checker.neg_point(point)));
        assert!(sign_orbits.insert(orbit));
    }
    let signed: Vec<Point> = base
        .iter()
        .flat_map(|&point| [point, checker.neg_point(point)])
        .collect();
    let table_points: Vec<Point> = if arm == "isogeny-image" {
        signed.iter().map(|&point| checker.velu(point, t)).collect()
    } else {
        signed.clone()
    };
    let table_a = if arm == "isogeny-image" { 53 } else { curve_a };
    let mut support = BTreeSet::new();
    let mut entries = 0u64;
    for i in 0..table_points.len() {
        for j in i..table_points.len() {
            let sum = checker.add_points(table_points[i], table_points[j], table_a);
            if !sum.infinity {
                support.insert(checker.key(sum));
                entries += 1;
            }
        }
    }
    assert_eq!(row["pair_entries"], entries);
    assert_eq!(row["pair_keys"], support.len());
    let witnesses = row["relation_witnesses"].as_array().unwrap();
    assert_eq!(witnesses.len(), row["relations"].as_u64().unwrap() as usize);
    for witness in witnesses {
        let a = witness["a"].as_u64().unwrap();
        let b = witness["b"].as_u64().unwrap();
        let i = witness["i"].as_u64().unwrap() as usize;
        let j = witness["j"].as_u64().unwrap() as usize;
        assert!(i <= j && j < table_points.len());
        let query_a = if direct { 53 } else { 46 };
        let query_g = if direct { mapped_g } else { g };
        let query_q = if direct { mapped_q } else { q };
        let candidate = checker.add_points(
            checker.scalar_mul(query_g, a, query_a),
            checker.scalar_mul(query_q, b, query_a),
            query_a,
        );
        let projected = checker.scalar_mul(candidate, 75, query_a);
        let searched = if arm == "isogeny-image" {
            checker.velu(projected, t)
        } else {
            projected
        };
        assert!(!searched.infinity);
        assert_eq!(
            checker.add_points(table_points[i], table_points[j], table_a),
            searched
        );
        if arm == "isogeny-image" {
            assert_eq!(checker.add_points(signed[i], signed[j], 46), projected);
        }
    }
    witnesses.len() as u64
}

fn replay(public_path: &str, fixture_path: &str, records_path: &str, old_public: Option<&str>) {
    let public = read_public(public_path);
    let full_image = public["schema"] == "different_j_public/v2";
    assert_eq!(full_image, old_public.is_some());
    if let Some(old_path) = old_public {
        let old = read_public(old_path);
        assert_eq!(old["schema"], "different_j_public/v1");
        let mut normalized = public.clone();
        let object = normalized.as_object_mut().unwrap();
        for name in [
            "full_base_raw",
            "full_base_projected",
            "full_base_witness_source",
            "full_image_source_points",
            "full_image_points",
            "full_image_selection_candidates",
        ] {
            assert!(object.remove(name).is_some());
        }
        object.insert("schema".to_string(), json!("different_j_public/v1"));
        assert_eq!(normalized, old);
    }
    let fixture: Value = serde_json::from_str(&fs::read_to_string(fixture_path).unwrap()).unwrap();
    let answer = fixture["answer_scalar"].as_u64().unwrap();
    let g = parse_point(&public["generator"]);
    let q = parse_point(&public["target"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let mapped_q = parse_point(&public["mapped_target"]);
    let checker = Replay::new();
    assert!(checker.on_curve(g, 46, 37) && checker.on_curve(q, 46, 37));
    assert!(checker.on_curve(mapped_g, 53, 25) && checker.on_curve(mapped_q, 53, 25));
    assert!(checker.scalar_mul(g, 61, 46).infinity);
    assert!(checker.scalar_mul(mapped_g, 61, 53).infinity);
    for scalar in 1..61 {
        assert_eq!(
            checker.scalar_mul(g, scalar, 46).infinity,
            checker.scalar_mul(mapped_g, scalar, 53).infinity
        );
    }
    let t = parse_point(&public["kernel"][0]);
    let mut base_witness_checks = 0u64;
    if full_image {
        let raw = points(&public["full_base_raw"]);
        let projected = points(&public["full_base_projected"]);
        let witnesses = points(&public["full_base_witness_source"]);
        assert_eq!((raw.len(), projected.len(), witnesses.len()), (8, 8, 8));
        for ((&point, &projection), &preimage) in raw.iter().zip(&projected).zip(&witnesses) {
            assert!(checker.on_curve(preimage, 46, 37));
            assert!(checker.on_curve(point, 53, 25));
            assert_eq!(checker.velu(preimage, t), point);
            assert_eq!(checker.scalar_mul(point, 75, 53), projection);
            assert!(checker.scalar_mul(projection, 61, 53).infinity);
            base_witness_checks += 1;
        }
    }
    let mut checks = BTreeMap::<(u64, String), u64>::new();
    let mut relation_witness_checks = 0u64;
    for line in fs::read_to_string(records_path).unwrap().lines() {
        let row: Value = serde_json::from_str(line).unwrap();
        let schema = row["schema"].as_str().unwrap();
        assert!([
            "different_j_relation_run/v1",
            "different_j_relation_run/v2",
            "different_j_rho_run/v1"
        ]
        .contains(&schema));
        let arm = row["arm"].as_str().unwrap().to_string();
        let round = row["round"].as_u64().unwrap();
        assert!(round < 3);
        let scalar = row["recovered"].as_u64().expect("complete scalar");
        assert_eq!(row["verified"], true);
        assert_eq!(scalar, answer);
        assert_eq!(checker.scalar_mul(g, scalar, 46), q);
        assert_eq!(checker.scalar_mul(mapped_g, scalar, 53), mapped_q);
        if full_image && schema == "different_j_relation_run/v2" {
            relation_witness_checks +=
                replay_relation_row(&checker, &public, &row, g, q, mapped_g, mapped_q, t);
        }
        *checks.entry((round, arm)).or_default() += 1;
    }
    assert_eq!(checks.len(), if full_image { 18 } else { 15 });
    assert!(checks.values().all(|&v| v == 1));
    if full_image {
        for round in 0..3 {
            assert_eq!(checks.get(&(round, "full-image".to_string())), Some(&1));
        }
    }
    println!(
        "{}",
        json!({"schema":if full_image {"different_j_replay/v2"} else {"different_j_replay/v1"},
        "checked":checks.len(),"base_witness_checks":base_witness_checks,
        "relation_witness_checks":relation_witness_checks,
        "all_rounds_and_arms":true,"independent_fp_muls":checker.muls.get(),
        "independent_fp_invs":checker.inversions.get(),
        "independent_group_adds":checker.group_adds.get(),
        "fixture_sha256":sha256_hex(&fs::read(fixture_path).unwrap()),
        "source_public_sha256":old_public.map(|path| sha256_hex(&fs::read(path).unwrap())),
        "public_sha256":sha256_hex(&fs::read(public_path).unwrap()),
        "records_sha256":sha256_hex(&fs::read(records_path).unwrap())})
    );
}

fn audit_fibers(public_path: &str) {
    let public = read_public(public_path);
    let src = source(67);
    let t = parse_point(&public["kernel"][0]);
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    take_field_counters();
    let dst = destination(&src, t);
    let mut ops = GroupOps::default();
    let mut fibers = BTreeMap::<u64, u64>::new();
    let mut source_points = 0u64;
    let mut insert = |point: Point, ops: &mut GroupOps| {
        let mapped = velu(&src, t, point, ops);
        assert!(dst.is_on_curve(mapped));
        *fibers.entry(dst.key(&mapped)).or_default() += 1;
        source_points += 1;
    };
    insert(src.identity(), &mut ops);
    for x0 in 0..67 {
        for x1 in 0..67 {
            for point in src.lift_x([x0, x1]) {
                insert(point, &mut ops);
            }
        }
    }
    assert_eq!(source_points, 4575);
    let mut distribution = BTreeMap::<u64, u64>::new();
    for &size in fibers.values() {
        *distribution.entry(size).or_default() += 1;
    }
    assert_eq!(fibers.len(), 1525);
    assert_eq!(distribution.get(&3), Some(&1525));
    let g = parse_point(&public["generator"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let subgroup: Vec<Point> = (0..61).map(|i| src.mul(&mut ops, g, i)).collect();
    let image: Vec<Point> = subgroup
        .iter()
        .map(|&point| velu(&src, t, point, &mut ops))
        .collect();
    assert_eq!(image[1], mapped_g);
    let unique_images: BTreeSet<_> = image.iter().map(|point| dst.key(point)).collect();
    assert_eq!(unique_images.len(), 61);
    let mut checked_pairs = 0u64;
    for a in 0..61usize {
        for b in 0..61usize {
            let k = (a + b) % 61;
            assert_eq!(src.add(&mut ops, subgroup[a], subgroup[b]), subgroup[k]);
            assert_eq!(dst.add(&mut ops, image[a], image[b]), image[k]);
            assert_eq!(velu(&src, t, subgroup[k], &mut ops), image[k]);
            checked_pairs += 1;
        }
    }
    let fp = take_field_counters();
    println!(
        "{}",
        json!({"schema":"different_j_fiber_audit/v1",
        "p":67,"r":61,"source_points":source_points,
        "image_points":fibers.len(),"fiber_size_distribution":distribution,
        "kernel":[t,src.neg(t)],"all_images_on_destination":true,
        "subgroup_images":unique_images.len(),
        "subgroup_pair_homomorphism_checks":checked_pairs,
        "source_icv1":public["source"]["icv1"],
        "destination_icv1":public["destination"]["icv1"],
        "group_ops":ops,"fp_muls":fp.0,"fp_invs":fp.1,
        "cpu_ns":thread_cpu_ns()-cpu,"wall_ns":wall.elapsed().as_nanos() as u64,
        "peak_rss_bytes":peak_rss_bytes()})
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("screen") => screen(),
        Some("freeze") => freeze(),
        Some("freeze-public") => freeze_public(&args[2], &args[3]),
        Some("freeze-full-image") => freeze_full_image(&args[2], &args[3]),
        Some("run") => relation_run(&args[2], args[3].parse().unwrap(), &args[4]),
        Some("rho") => rho_run(&args[2], args[3].parse().unwrap(), args[4] == "aa"),
        Some("replay") => replay(&args[2], &args[3], &args[4], None),
        Some("replay-full-image") => replay(&args[2], &args[4], &args[5], Some(&args[3])),
        Some("audit-fibers") => audit_fibers(&args[2]),
        _ => panic!(
            "usage: nonisomorphic_extension_isogeny screen|freeze|freeze-public|freeze-full-image|run|rho|replay|replay-full-image|audit-fibers"
        ),
    }
}
