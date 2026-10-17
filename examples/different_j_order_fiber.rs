//! Full rational-fiber and order-229 subgroup audit for the selected p=211 pair.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ext_curve::{ExtCurve, ExtField, ExtPoint, Fp2};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use serde_json::json;

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;

fn key(p: Point) -> [u64; 5] {
    [p.infinity as u64, p.x[0], p.x[1], p.y[0], p.y[1]]
}

fn key_text(p: Point) -> String {
    let k = key(p);
    format!("{},{},{},{},{}", k[0], k[1], k[2], k[3], k[4])
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

fn all_points(curve: &Curve) -> Vec<Point> {
    let f = &curve.f;
    let mut roots = BTreeMap::<[u64; 2], Vec<[u64; 2]>>::new();
    for y0 in 0..f.p {
        for y1 in 0..f.p {
            let y = [y0, y1];
            roots.entry(f.sqr(y)).or_default().push(y);
        }
    }
    let mut points = vec![curve.identity()];
    for x0 in 0..f.p {
        for x1 in 0..f.p {
            let x = [x0, x1];
            let rhs = f.add(f.add(f.mul(f.sqr(x), x), f.mul(curve.a, x)), curve.b);
            if let Some(ys) = roots.get(&rhs) {
                for &y in ys {
                    let point = Point::affine(x, y);
                    assert!(curve.is_on_curve(point));
                    points.push(point);
                }
            }
        }
    }
    points.sort_by_key(|&p| key(p));
    assert!(points.windows(2).all(|w| key(w[0]) != key(w[1])));
    points
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

fn main() {
    assert_eq!(std::env::args().nth(1).as_deref(), Some("audit"));
    let wall = Instant::now();
    let f = Fp2::new(211).unwrap();
    assert_eq!(f.nu, 2);
    let src = Curve {
        f,
        a: f.from_fp(190),
        b: f.from_fp(37),
    };
    let dst = Curve {
        f,
        a: f.from_fp(130),
        b: f.from_fp(179),
    };
    let t = Point::affine(f.from_fp(3), f.one());
    let mut ops = GroupOps::default();
    assert!(src.is_on_curve(t));
    assert!(src.mul(&mut ops, t, 3).infinity);
    let points = all_points(&src);
    assert_eq!(points.len(), 44_655);
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    assert!(wall.elapsed().as_secs() < 120);

    let mut fibers = BTreeMap::<[u64; 5], Vec<Point>>::new();
    let mut stream = String::new();
    for &point in &points {
        let mapped = velu(&src, t, point, &mut ops);
        assert!(dst.is_on_curve(mapped));
        fibers.entry(key(mapped)).or_default().push(point);
        stream.push_str(&key_text(point));
        stream.push('|');
        stream.push_str(&key_text(mapped));
        stream.push('\n');
    }
    assert_eq!(fibers.len(), 14_885);
    assert!(fibers.values().all(|fiber| fiber.len() == 3));
    let kernel = fibers.get(&key(dst.identity())).unwrap();
    let expected_kernel = BTreeSet::from([key(src.identity()), key(t), key(src.neg(t))]);
    let actual_kernel: BTreeSet<[u64; 5]> = kernel.iter().copied().map(key).collect();
    assert_eq!(actual_kernel, expected_kernel);

    let mut g = None;
    for &point in &points {
        let candidate = src.mul(&mut ops, point, 195);
        if !candidate.infinity && src.mul(&mut ops, candidate, 229).infinity {
            g = Some(candidate);
            break;
        }
    }
    let g = g.expect("order-229 point");
    let mg = velu(&src, t, g, &mut ops);
    assert!(!mg.infinity && dst.mul(&mut ops, mg, 229).infinity);
    let mut source_subgroup = Vec::with_capacity(229);
    let mut image_subgroup = Vec::with_capacity(229);
    let mut subgroup_stream = String::new();
    let mut current = src.identity();
    for _ in 0..229 {
        let mapped = velu(&src, t, current, &mut ops);
        assert!(dst.is_on_curve(mapped));
        source_subgroup.push(current);
        image_subgroup.push(mapped);
        subgroup_stream.push_str(&key_text(current));
        subgroup_stream.push('|');
        subgroup_stream.push_str(&key_text(mapped));
        subgroup_stream.push('\n');
        current = src.add(&mut ops, current, g);
    }
    assert!(current.infinity);
    assert_eq!(
        source_subgroup
            .iter()
            .copied()
            .map(key)
            .collect::<BTreeSet<_>>()
            .len(),
        229
    );
    assert_eq!(
        image_subgroup
            .iter()
            .copied()
            .map(key)
            .collect::<BTreeSet<_>>()
            .len(),
        229
    );
    let mut checks = 0usize;
    for &a in &source_subgroup {
        for &b in &source_subgroup {
            let lhs = velu(&src, t, src.add(&mut ops, a, b), &mut ops);
            let ma = velu(&src, t, a, &mut ops);
            let mb = velu(&src, t, b, &mut ops);
            let rhs = dst.add(&mut ops, ma, mb);
            assert_eq!(lhs, rhs);
            checks += 1;
        }
    }
    assert_eq!(checks, 229 * 229);
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    assert!(wall.elapsed().as_secs() < 120);
    println!(
        "{}",
        json!({"schema":"different_j_order_fiber/v1","p":211,"r":229,
            "source_points":points.len(),"destination_points":44655,
            "image_points":fibers.len(),"fiber_histogram":{"3":fibers.len()},
            "identity_fiber":actual_kernel,"full_map_sha256":sha256_hex(stream.as_bytes()),
            "source_generator":g,"image_generator":mg,
            "subgroup_points":source_subgroup.len(),"subgroup_map_sha256":sha256_hex(subgroup_stream.as_bytes()),
            "subgroup_pair_checks":checks,"peak_rss_bytes":peak_rss_bytes(),
            "group_ops":ops,"all_checks_passed":true})
    );
}
