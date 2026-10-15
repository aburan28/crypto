//! Native screen for an explicit different-j degree-three map over F_(p^2).
//! The selector and stop conditions are frozen in the dated study protocol.

use crypto_lib::cryptanalysis::ext_curve::{
    jacobi, random_point, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
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

fn main() {
    assert_eq!(std::env::args().nth(1).as_deref(), Some("screen"));
    screen();
}
