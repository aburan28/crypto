//! Native instance screen for the frozen quadratic-extension comparison.
//! `screen` prints candidate parameters only; the paired DLP run is a
//! separate, frozen stage in the study protocol.

use std::collections::BTreeSet;

use crypto_lib::cryptanalysis::ext_curve::jacobi;
use crypto_lib::cryptanalysis::ext_curve::{random_point, ExtCurve, ExtField, Fp2};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

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
        _ => panic!("usage: isogeny_torsion_extension screen"),
    }
}
