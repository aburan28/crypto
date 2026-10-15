//! Independent field, curve and Vélu replay of the frozen order-trend fibers.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

type E = [u64; 2];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Point {
    x: E,
    y: E,
    infinity: bool,
}

#[derive(Clone, Copy)]
struct Field {
    p: u64,
    nu: u64,
}

impl Field {
    fn add(self, a: E, b: E) -> E {
        [(a[0] + b[0]) % self.p, (a[1] + b[1]) % self.p]
    }

    fn sub(self, a: E, b: E) -> E {
        [
            (a[0] + self.p - b[0]) % self.p,
            (a[1] + self.p - b[1]) % self.p,
        ]
    }

    fn mul(self, a: E, b: E) -> E {
        let p = self.p as u128;
        let real =
            (a[0] as u128 * b[0] as u128 + self.nu as u128 * a[1] as u128 * b[1] as u128) % p;
        let imag = (a[0] as u128 * b[1] as u128 + a[1] as u128 * b[0] as u128) % p;
        [real as u64, imag as u64]
    }

    fn inv(self, a: E) -> E {
        let norm = (a[0] * a[0] + self.p - self.nu * a[1] * a[1] % self.p) % self.p;
        assert_ne!(norm, 0);
        let inv_norm = self.pow_fp(norm, self.p - 2);
        self.mul([a[0], (self.p - a[1]) % self.p], [inv_norm, 0])
    }

    fn pow_fp(self, mut a: u64, mut e: u64) -> u64 {
        let mut result = 1u64;
        while e != 0 {
            if e & 1 != 0 {
                result = result * a % self.p;
            }
            a = a * a % self.p;
            e >>= 1;
        }
        result
    }

    fn scale(self, a: E, n: u64) -> E {
        self.mul(a, [n % self.p, 0])
    }

    fn neg(self, a: E) -> E {
        self.sub([0, 0], a)
    }
}

fn identity() -> Point {
    Point {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    }
}

fn key(p: Point) -> [u64; 5] {
    [p.infinity as u64, p.x[0], p.x[1], p.y[0], p.y[1]]
}

fn key_text(p: Point) -> String {
    let k = key(p);
    format!("{},{},{},{},{}", k[0], k[1], k[2], k[3], k[4])
}

fn parse_point(v: &Value) -> Point {
    Point {
        x: [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        y: [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
        infinity: v["infinity"].as_bool().unwrap(),
    }
}

fn neg(f: Field, p: Point) -> Point {
    Point { y: f.neg(p.y), ..p }
}

fn add(f: Field, a: E, p: Point, q: Point) -> Point {
    if p.infinity {
        return q;
    }
    if q.infinity {
        return p;
    }
    if p.x == q.x && p.y == f.neg(q.y) {
        return identity();
    }
    let slope = if p == q {
        let numer = f.add(f.scale(f.mul(p.x, p.x), 3), a);
        f.mul(numer, f.inv(f.scale(p.y, 2)))
    } else {
        f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)))
    };
    let x = f.sub(f.sub(f.mul(slope, slope), p.x), q.x);
    let y = f.sub(f.mul(slope, f.sub(p.x, x)), p.y);
    Point {
        x,
        y,
        infinity: false,
    }
}

fn scalar(f: Field, a: E, mut p: Point, mut k: u64) -> Point {
    let mut result = identity();
    while k != 0 {
        if k & 1 != 0 {
            result = add(f, a, result, p);
        }
        p = add(f, a, p, p);
        k >>= 1;
    }
    result
}

fn on_curve(f: Field, a: E, b: E, p: Point) -> bool {
    p.infinity || f.mul(p.y, p.y) == f.add(f.add(f.mul(f.mul(p.x, p.x), p.x), f.mul(a, p.x)), b)
}

fn velu(f: Field, a: E, t: Point, p: Point) -> Point {
    if p.infinity || p == t || p == neg(f, t) {
        return identity();
    }
    let plus = add(f, a, p, t);
    let minus = add(f, a, p, neg(f, t));
    assert!(!plus.infinity && !minus.infinity);
    Point {
        x: f.sub(f.add(f.add(p.x, plus.x), minus.x), f.scale(t.x, 2)),
        y: f.add(f.add(p.y, plus.y), minus.y),
        infinity: false,
    }
}

fn all_points(f: Field, a: E, b: E) -> Vec<Point> {
    let mut square_roots = BTreeMap::<E, Vec<E>>::new();
    for y0 in 0..f.p {
        for y1 in 0..f.p {
            let y = [y0, y1];
            square_roots.entry(f.mul(y, y)).or_default().push(y);
        }
    }
    let mut points = vec![identity()];
    for x0 in 0..f.p {
        for x1 in 0..f.p {
            let x = [x0, x1];
            let rhs = f.add(f.add(f.mul(f.mul(x, x), x), f.mul(a, x)), b);
            if let Some(ys) = square_roots.get(&rhs) {
                for &y in ys {
                    let p = Point {
                        x,
                        y,
                        infinity: false,
                    };
                    assert!(on_curve(f, a, b, p));
                    points.push(p);
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
    let receipt_path = std::env::args().nth(1).expect("runner receipt JSON");
    let receipt: Value = serde_json::from_slice(&fs::read(receipt_path).unwrap()).unwrap();
    assert_eq!(receipt["schema"], "different_j_trend_fiber/v1");
    assert_eq!(receipt["all_checks_passed"], true);
    let wall = Instant::now();
    let p = receipt["p"].as_u64().unwrap();
    let r = receipt["r"].as_u64().unwrap();
    let h = receipt["h"].as_u64().unwrap();
    assert!([(409, 433, 387), (421, 457, 387), (457, 421, 495)].contains(&(p, r, h)));
    let nu = (2..p)
        .find(|&n| {
            let candidate = Field { p, nu: n };
            candidate.pow_fp(n, (p - 1) / 2) == p - 1
        })
        .unwrap();
    let f = Field { p, nu };
    assert_eq!(receipt["field_nu"].as_u64(), Some(nu));
    let sa = [p - 21, 0];
    let sb = [37, 0];
    let da = [p - 81, 0];
    let db = [p - 243, 0];
    let t = Point {
        x: [3, 0],
        y: [1, 0],
        infinity: false,
    };
    assert!(on_curve(f, sa, sb, t));
    assert!(scalar(f, sa, t, 3).infinity);
    let points = all_points(f, sa, sb);
    assert_eq!(points.len() as u64, h * r);
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    let mut fibers = BTreeMap::<[u64; 5], Vec<Point>>::new();
    let mut map_stream = String::new();
    for &point in &points {
        let image = velu(f, sa, t, point);
        assert!(on_curve(f, da, db, image));
        fibers.entry(key(image)).or_default().push(point);
        map_stream.push_str(&key_text(point));
        map_stream.push('|');
        map_stream.push_str(&key_text(image));
        map_stream.push('\n');
    }
    assert_eq!(fibers.len() as u64, h * r / 3);
    assert!(fibers.values().all(|v| v.len() == 3));
    let kernel: BTreeSet<[u64; 5]> = fibers
        .get(&key(identity()))
        .unwrap()
        .iter()
        .copied()
        .map(key)
        .collect();
    assert_eq!(
        kernel,
        BTreeSet::from([key(identity()), key(t), key(neg(f, t))])
    );
    let map_hash = sha256_hex(map_stream.as_bytes());
    assert_eq!(receipt["source_points"].as_u64(), Some(points.len() as u64));
    assert_eq!(receipt["image_points"].as_u64(), Some(fibers.len() as u64));
    assert_eq!(receipt["fiber_histogram"], json!({"3":fibers.len()}));
    assert_eq!(receipt["identity_fiber"], json!(kernel));
    assert_eq!(receipt["full_map_sha256"], map_hash);

    let mut generator = None;
    for &point in &points {
        let projected = scalar(f, sa, point, h);
        if !projected.infinity && scalar(f, sa, projected, r).infinity {
            generator = Some(projected);
            break;
        }
    }
    let g = generator.expect("order-r point");
    let mg = velu(f, sa, t, g);
    assert_eq!(g, parse_point(&receipt["source_generator"]));
    assert_eq!(mg, parse_point(&receipt["image_generator"]));
    assert!(!mg.infinity && scalar(f, da, mg, r).infinity);
    let mut source_subgroup = Vec::with_capacity(r as usize);
    let mut image_subgroup = Vec::with_capacity(r as usize);
    let mut subgroup_stream = String::new();
    let mut current = identity();
    for _ in 0..r {
        let image = velu(f, sa, t, current);
        assert!(on_curve(f, da, db, image));
        source_subgroup.push(current);
        image_subgroup.push(image);
        subgroup_stream.push_str(&key_text(current));
        subgroup_stream.push('|');
        subgroup_stream.push_str(&key_text(image));
        subgroup_stream.push('\n');
        current = add(f, sa, current, g);
    }
    assert!(current.infinity);
    assert_eq!(
        source_subgroup
            .iter()
            .copied()
            .map(key)
            .collect::<BTreeSet<_>>()
            .len(),
        r as usize
    );
    assert_eq!(
        image_subgroup
            .iter()
            .copied()
            .map(key)
            .collect::<BTreeSet<_>>()
            .len(),
        r as usize
    );
    let subgroup_hash = sha256_hex(subgroup_stream.as_bytes());
    assert_eq!(receipt["subgroup_map_sha256"], subgroup_hash);
    assert_eq!(receipt["subgroup_points"].as_u64(), Some(r));
    let mut checks = 0usize;
    for (i, &p) in source_subgroup.iter().enumerate() {
        for (j, &q) in source_subgroup.iter().enumerate() {
            let lhs = velu(f, sa, t, add(f, sa, p, q));
            let rhs = add(f, da, image_subgroup[i], image_subgroup[j]);
            assert_eq!(lhs, rhs);
            checks += 1;
        }
    }
    assert_eq!(checks as u64, r * r);
    assert_eq!(
        receipt["subgroup_pair_checks"].as_u64(),
        Some(checks as u64)
    );
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    assert!(wall.elapsed().as_secs() < 300, "checker wall cap exhausted");
    println!(
        "{}",
        json!({"schema":"different_j_trend_fiber_check/v1","p":p,"r":r,"h":h,
        "status":"pass","source_points":points.len(),"image_points":fibers.len(),
        "three_point_fibers":fibers.len(),"subgroup_points":r,
        "subgroup_pair_checks":checks,"peak_rss_bytes":peak_rss_bytes(),
        "full_map_sha256":map_hash,
        "subgroup_map_sha256":subgroup_hash})
    );
}
