//! Native instance screen for the frozen quadratic-extension comparison.
//! `screen` prints candidate parameters only; the paired DLP run is a
//! separate, frozen stage in the study protocol.

use std::collections::BTreeSet;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ext_curve::jacobi;
use crypto_lib::cryptanalysis::ext_curve::{random_point, ExtCurve, ExtField, ExtPoint, Fp2};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::hash::sha256::sha256;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

type Point = ExtPoint<[u64; 2]>;

fn digest(data: &str) -> String {
    sha256_hex(data.as_bytes())
}

fn icv1(p: u64, nu: u64, n2: u64, b: u64) -> serde_json::Value {
    let c0 = (p - nu) % p;
    let field = format!(
        "fpk-{p}-2-{}",
        &digest(&format!("fpk-modulus:{p}:{c0},0"))[..8]
    );
    let model = format!("{{\"a\":[\"0\",\"0\"],\"b\":[\"{b}\",\"0\"],\"field\":\"{field}\",\"form\":\"y^2=x^3+a*x+b\",\"k\":\"2\",\"modulus\":[\"{c0}\",\"0\"],\"p\":\"{p}\",\"v\":\"1\"}}");
    let model_hash = digest(&model);
    let trace = p as i64 * p as i64 + 1 - n2 as i64;
    let id = format!(
        "ICV1:{field}:{trace}:{n2}:0,0:unk:unk:r:{}",
        &model_hash[..12]
    );
    json!({"id":id,"model":model,"field":field,"model_sha256":model_hash})
}

fn chi(curve: &ExtCurve<Fp2>, point: Point) -> Option<[u64; 2]> {
    let f = &curve.f;
    if point.infinity {
        return Some(f.one());
    }
    let den = f.sub(point.y, f.one());
    f.inv(den).map(|inv| f.mul(f.add(point.y, f.one()), inv))
}

fn mobius(f: &Fp2, z: Option<[u64; 2]>) -> Option<[u64; 2]> {
    match z {
        None => Some(f.zero()),
        Some(z) => f.inv(f.sub(f.one(), z)),
    }
}

fn phi(curve: &ExtCurve<Fp2>, point: Point) -> Point {
    let f = &curve.f;
    if point.infinity || point.x == f.zero() {
        return curve.identity();
    }
    let x2 = f.sqr(point.x);
    let x3 = f.mul(x2, point.x);
    let xp = f.mul(f.add(x3, f.from_fp(4)), f.inv(x2).unwrap());
    let yp = f.mul(f.mul(point.y, f.sub(x3, f.from_fp(8))), f.inv(x3).unwrap());
    Point::affine(xp, yp)
}

fn freeze() {
    // Chosen by selector v2, SCREEN_V2.jsonl. The solver reads these
    // parameters, never the answer scalar supplied only to verification.
    for (index, (p, expected_r)) in [(61, 19), (271, 61), (787, 211), (3319, 811)]
        .into_iter()
        .enumerate()
    {
        let n1 = base_order(p);
        let trace = p as i64 + 1 - n1 as i64;
        let n2 = n1 * (p as i64 + 1 + trace) as u64;
        let r = factor_u64(n2 / n1).last().unwrap().0;
        assert_eq!(r, expected_r);
        let h = n2 / r;
        assert_ne!(h % r, 0);
        let f = Fp2::new(p).unwrap();
        let curve = ExtCurve {
            f,
            a: f.zero(),
            b: f.one(),
        };
        let dest = ExtCurve {
            f,
            a: f.zero(),
            b: f.from_fp(p - 27),
        };
        let line = f.s();
        let columns = if index == 0 {
            8
        } else {
            ((r as f64).sqrt().ceil() as usize).min(32)
        };
        let mut raw = Vec::new();
        let mut projected = Vec::new();
        let mut seen = BTreeSet::new();
        let mut ops = GroupOps::default();
        let mut scanned = 0u64;
        for t in 1..p {
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
        assert_eq!(raw.len(), columns, "insufficient projected base");
        let mut rng = StdRng::seed_from_u64(0x2026_1014 ^ p);
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        assert!(curve.mul(&mut ops, generator, r).infinity);
        assert!(generator.x[1] != 0 || generator.y[1] != 0);
        let target_seed = format!("torsion-extension/2026-10-14/p={p}/r={r}");
        let hash = sha256(target_seed.as_bytes());
        let k = 1 + u64::from_be_bytes(hash[0..8].try_into().unwrap()) % (r - 1);
        let target = curve.mul(&mut ops, generator, k);
        assert!(curve.is_on_curve(target));
        let torsion = Point::affine(f.zero(), f.one());
        assert!(curve.mul(&mut ops, torsion, 3).infinity);
        let mut chi_checks = 0;
        let mut map_checks = 0;
        for _ in 0..64 {
            let a = random_point(&curve, &mut rng);
            let b = random_point(&curve, &mut rng);
            let sum = curve.add(&mut ops, a, b);
            assert_eq!(
                chi(&curve, curve.add(&mut ops, a, torsion)),
                mobius(&f, chi(&curve, a))
            );
            chi_checks += 1;
            let mapped = phi(&curve, a);
            assert!(dest.is_on_curve(mapped));
            assert_eq!(phi(&curve, sum), dest.add(&mut ops, mapped, phi(&curve, b)));
            map_checks += 1;
        }
        assert_eq!(phi(&curve, torsion), dest.identity());
        assert_eq!(phi(&curve, curve.neg(torsion)), dest.identity());
        println!(
            "{}",
            json!({
                "schema":"torsion_extension_fixture/v1", "p":p,"nu":f.nu,
                "source":icv1(p,f.nu,n2,1),"destination":icv1(p,f.nu,n2,p-27),
                "n1":n1,"trace_fp":trace,"n2":n2,"r":r,"cofactor":h,
                "line":line,"generator":generator,"target":target,"answer_scalar":k,
                "target_seed":target_seed,"base_raw":raw,"base_projected":projected,
                "columns":columns,"raw_points_scanned":scanned,
                "chi_checks":chi_checks,"map_checks":map_checks,"setup_group_ops":ops,
                "embedding_degree":embedding_degree(p,r),
            })
        );
    }
}

fn base_order(p: u64) -> u64 {
    let mut n = 1;
    for x in 0..p {
        let rhs = ((x as u128 * x as u128 % p as u128 * x as u128 + 1) % p as u128) as u64;
        n += (1 + jacobi(rhs, p)) as u64;
    }
    n
}

fn embedding_degree(p: u64, r: u64) -> u64 {
    let mut power = p % r;
    for k in 1..=r {
        if power == 1 {
            return k;
        }
        power = ((power as u128 * p as u128) % r as u128) as u64;
    }
    unreachable!("r must be prime and coprime to p")
}

fn screen() {
    let intervals = [
        (31, 127, 17),
        (128, 511, 53),
        (512, 2047, 211),
        (2048, 4093, 811),
    ];
    for (lo, hi, min_r) in intervals {
        let mut selected = None;
        for p in lo..=hi {
            if p % 3 != 1 || !is_prime_u64(p) {
                continue;
            }
            let n1 = base_order(p);
            let t = p as i64 + 1 - n1 as i64;
            let anti = (p as i64 + 1 + t) as u64;
            let n2 = n1 * anti;
            let r = factor_u64(anti).last().map(|&(v, _)| v).unwrap_or(0);
            if r < min_r || r == 3 || n1 % r == 0 || (n2 / r) % r == 0 {
                continue;
            }
            let f = Fp2::new(p).unwrap();
            let curve = ExtCurve {
                f,
                a: f.zero(),
                b: f.one(),
            };
            let h = n2 / r;
            let mut rng = StdRng::seed_from_u64(0x2026_1014 ^ p);
            let mut ops = GroupOps::default();
            let mut generator = None;
            for _ in 0..64 {
                let raw = random_point(&curve, &mut rng);
                let g = curve.mul(&mut ops, raw, h);
                if !g.infinity && (g.x[1] != 0 || g.y[1] != 0) && curve.mul(&mut ops, g, r).infinity
                {
                    generator = Some(g);
                    break;
                }
            }
            let Some(g) = generator else {
                continue;
            };
            let mut line = None;
            for a in 0..p {
                let c = [a, 1];
                let mut projected = BTreeSet::new();
                let mut raw_points = 0u64;
                for s in 1..p {
                    let x = f.scale(c, s);
                    for point in curve.lift_x(x) {
                        raw_points += 1;
                        let q = curve.mul(&mut ops, point, h);
                        if !q.infinity {
                            projected.insert(curve.key(&q));
                        }
                    }
                    if projected.len() >= 16 {
                        break;
                    }
                }
                if projected.len() >= 16 {
                    line = Some((c, raw_points, projected.len()));
                    break;
                }
            }
            if let Some((c, raw_points, projected_points)) = line {
                selected = Some(json!({
                    "p":p,"nu":f.nu,"n1":n1,"trace":t,"n2":n2,
                    "r":r,"h":h,"embedding_degree":embedding_degree(p,r),
                    "generator":{"x":g.x,"y":g.y},
                    "line":c,"raw_points_inspected":raw_points,
                    "projected_signed_points_seen":projected_points,
                    "group_ops_screen":ops,
                }));
                break;
            }
        }
        println!(
            "{}",
            json!({"interval":[lo,hi],"min_r":min_r,"selected":selected})
        );
    }
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("screen") => screen(),
        Some("freeze") => freeze(),
        _ => panic!("usage: isogeny_torsion_extension screen"),
    }
}
