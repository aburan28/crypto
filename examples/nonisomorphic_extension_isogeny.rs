//! Native screen and frozen instance for an explicit different-j degree-three map over F_(p^2).
//! The selector and stop conditions are frozen in the dated study protocol.

use std::collections::BTreeSet;

use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex};
use crypto_lib::cryptanalysis::ext_curve::{
    jacobi, random_point, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::isogeny_walk::record::{ec1_identity, V};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

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
    })).unwrap();
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

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("screen") => screen(),
        Some("freeze") => freeze(),
        _ => panic!("usage: nonisomorphic_extension_isogeny screen|freeze"),
    }
}
