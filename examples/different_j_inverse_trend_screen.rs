//! Four-bin direct inverse-Velu order screen for the explicit different-j family.
//! See research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu/FOUR_GATE_PROTOCOL_V1.md.

use std::time::Instant;

use crypto_lib::cryptanalysis::ext_curve::{
    jacobi, random_point, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::glv_invariant_base::factor_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;
const BINS: [(u64, u64); 4] = [(2400, 4800), (4800, 9600), (9600, 19200), (19200, 38400)];

fn source(p: u64) -> Curve {
    let f = Fp2::new(p).expect("eligible prime field");
    Curve {
        f,
        a: f.from_fp((p - 21 % p) % p),
        b: f.from_fp(37 % p),
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

fn velu(src: &Curve, t: Point, point: Point, ops: &mut GroupOps) -> Point {
    if point.infinity || point == t || point == src.neg(t) {
        return src.identity();
    }
    let plus = src.add(ops, point, t);
    let minus = src.add(ops, point, src.neg(t));
    assert!(!plus.infinity && !minus.infinity);
    let f = &src.f;
    Point::affine(
        f.sub(f.add(f.add(point.x, plus.x), minus.x), f.scale(t.x, 2)),
        f.add(f.add(point.y, plus.y), minus.y),
    )
}

fn count_fp(curve: &Curve) -> u64 {
    let p = curve.f.p;
    let mut count = 1u64;
    for x in 0..p {
        let xx = curve.f.from_fp(x);
        let x3 = curve.f.mul(curve.f.sqr(xx), xx);
        let rhs = curve
            .f
            .add(curve.f.add(x3, curve.f.mul(curve.a, xx)), curve.b)[0];
        count += if rhs == 0 {
            1
        } else if jacobi(rhs, p) == 1 {
            2
        } else {
            0
        };
    }
    count
}

fn j(curve: &Curve) -> Option<[u64; 2]> {
    let f = &curve.f;
    let a3 = f.mul(f.sqr(curve.a), curve.a);
    let four_a3 = f.scale(a3, 4);
    let den = f.add(four_a3, f.scale(f.sqr(curve.b), 27));
    f.inv(den).map(|inv| f.mul(f.scale(four_a3, 1728), inv))
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

fn row(p: u64) -> Value {
    let src = source(p);
    let t = Point::affine(src.f.from_fp(3), src.f.one());
    let mut ops = GroupOps::default();
    let js = j(&src);
    if js.is_none() {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"source_singular"});
    }
    if !src.is_on_curve(t) || t.infinity || !src.mul(&mut ops, t, 3).infinity {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"kernel_check"});
    }
    let dst = destination(&src, t);
    let jd = j(&dst);
    if jd.is_none() {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"destination_singular"});
    }
    if js == jd {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"equal_j","j":js});
    }
    let n1 = count_fp(&src);
    let n1_dst = count_fp(&dst);
    if n1 != n1_dst {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"base_order_mismatch","n1":n1,"n1_dst":n1_dst});
    }
    let trace = p as i64 + 1 - n1 as i64;
    let n2 = ((p + 1) * (p + 1)) as i64 - trace * trace;
    assert!(n2 > 0);
    let n2 = n2 as u64;
    let factors = factor_u64(n2);
    let eligible_orders: Vec<u64> = factors
        .iter()
        .filter(|(prime, power)| {
            *prime != p && *power == 1 && BINS.iter().any(|&(lo, hi)| lo <= *prime && *prime < hi)
        })
        .map(|x| x.0)
        .collect();
    if eligible_orders.is_empty() {
        return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"subgroup","n1":n1,"n2":n2,"factors":factors});
    }
    assert!(velu(&src, t, t, &mut ops).infinity);
    assert!(velu(&src, t, src.neg(t), &mut ops).infinity);
    let mut rng = StdRng::seed_from_u64(0x2026_1015 ^ p);
    let mut samples = Vec::with_capacity(64);
    for _ in 0..64 {
        let a = random_point(&src, &mut rng);
        let b = random_point(&src, &mut rng);
        let sum = src.add(&mut ops, a, b);
        let ma = velu(&src, t, a, &mut ops);
        let mb = velu(&src, t, b, &mut ops);
        let msum = velu(&src, t, sum, &mut ops);
        if !dst.is_on_curve(ma) || !dst.is_on_curve(mb) || msum != dst.add(&mut ops, ma, mb) {
            return json!({"schema":"different_j_inverse_trend_screen/v1","p":p,"status":"reject","reason":"map_check","checked_samples":samples.len()});
        }
        samples.push(json!({"a":a,"b":b,"sum":sum,"ma":ma,"mb":mb,"msum":msum}));
    }
    json!({
        "schema":"different_j_inverse_trend_screen/v1","p":p,"status":"accept",
        "field_nu":src.f.nu,"source":{"a":src.a,"b":src.b,"j":js},
        "destination":{"a":dst.a,"b":dst.b,"j":jd},
        "kernel":[t,src.neg(t)],"n1":n1,"n1_dst":n1_dst,
        "trace":trace,"n2":n2,"factors":factors,"eligible_orders":eligible_orders,
        "sample_seed":0x2026_1015 ^ p,"samples":samples,"group_ops":ops
    })
}

fn main() {
    let extension = match std::env::args().nth(1).as_deref() {
        Some("screen-inverse-trend") => false,
        Some("screen-inverse-trend-extension") => true,
        _ => panic!("expected screen-inverse-trend or screen-inverse-trend-extension"),
    };
    let (min_p, max_p) = if extension {
        (10009, 30011)
    } else {
        (1238, 10007)
    };
    let start = Instant::now();
    let mut selected: [Option<Value>; 4] = std::array::from_fn(|_| None);
    let mut accepted = 0usize;
    let mut candidates = 0usize;
    for p in (min_p..=max_p).filter(|p| p % 3 == 1 && is_prime_u64(*p)) {
        assert!(
            start.elapsed().as_secs() < 3600,
            "whole-screen wall cap exhausted"
        );
        assert!(peak_rss_bytes() <= 512 * 1024 * 1024, "RSS cap exhausted");
        let record = row(p);
        if record["status"] == "accept" {
            accepted += 1;
            for r in record["eligible_orders"]
                .as_array()
                .unwrap()
                .iter()
                .filter_map(Value::as_u64)
            {
                for (bin, &(lo, hi)) in BINS.iter().enumerate() {
                    if extension && bin != 3 {
                        continue;
                    }
                    if selected[bin].is_none() && lo <= r && r < hi {
                        selected[bin] = Some(
                            json!({"p":p,"r":r,"cofactor":record["n2"].as_u64().unwrap()/r,"bin":[lo,hi]}),
                        );
                    }
                }
            }
        }
        candidates += 1;
        println!("{record}");
    }
    assert!(
        start.elapsed().as_secs() < 3600,
        "whole-screen wall cap exhausted"
    );
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024, "RSS cap exhausted");
    println!(
        "{}",
        json!({"schema":if extension {"different_j_inverse_trend_screen_summary/v2"}
                 else {"different_j_inverse_trend_screen_summary/v1"},
        "candidates":candidates,"accepted":accepted,"selected":selected,
        "min_p":min_p,"max_p":max_p,"order_bins":BINS,"prime_to_characteristic":true,
        "v1_raw_sha256":if extension {Some("ef232a54b22938239819fd6b14182d1e49de5778e9c00542ce257c54097608a5")} else {None},
        "rss_ceiling_bytes":512*1024*1024,"wall_ceiling_seconds":3600})
    );
}
