//! Independent arithmetic replay for the bounded different-j candidate screen.

use std::fs;

use serde_json::{json, Value};

#[derive(Clone, Copy, Debug)]
struct Field {
    p: u64,
    nu: u64,
}

type E = [u64; 2];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Point {
    x: E,
    y: E,
    infinity: bool,
}

impl Field {
    fn scalar(self, n: i64) -> u64 {
        n.rem_euclid(self.p as i64) as u64
    }

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

    fn mul_scalar(self, a: E, n: u64) -> E {
        self.mul(a, [n % self.p, 0])
    }

    fn inv(self, a: E) -> E {
        let p = self.p;
        let norm = (a[0] * a[0] + p - self.nu * a[1] * a[1] % p) % p;
        assert_ne!(norm, 0);
        let inv_norm = pow_mod(norm, p - 2, p);
        self.mul([a[0], (p - a[1]) % p], [inv_norm, 0])
    }

    fn neg(self, a: E) -> E {
        self.sub([0, 0], a)
    }
}

fn pow_mod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut z = 1u64;
    while e != 0 {
        if e & 1 != 0 {
            z = ((z as u128 * a as u128) % p as u128) as u64;
        }
        a = ((a as u128 * a as u128) % p as u128) as u64;
        e >>= 1;
    }
    z
}

fn prime(n: u64) -> bool {
    n >= 2
        && (2..)
            .take_while(|d| d * d <= n)
            .all(|d| !n.is_multiple_of(d))
}

fn factors(mut n: u64) -> Vec<(u64, u32)> {
    let mut result = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            let mut power = 0;
            while n.is_multiple_of(d) {
                n /= d;
                power += 1;
            }
            result.push((d, power));
        }
        d += 1;
    }
    if n > 1 {
        result.push((n, 1));
    }
    result
}

fn first_nonsquare(p: u64) -> u64 {
    (2..p)
        .find(|&x| pow_mod(x, (p - 1) / 2, p) == p - 1)
        .unwrap()
}

fn j(a: u64, b: u64, p: u64) -> Option<u64> {
    let a3 = pow_mod(a, 3, p);
    let b2 = pow_mod(b, 2, p);
    let v = 4 * a3 % p;
    let den = (v + 27 * b2) % p;
    (den != 0).then(|| 1728 * v % p * pow_mod(den, p - 2, p) % p)
}

fn count_fp(a: u64, b: u64, p: u64) -> u64 {
    1 + (0..p)
        .map(|x| {
            let rhs = (pow_mod(x, 3, p) + a * x + b) % p;
            if rhs == 0 {
                1
            } else if pow_mod(rhs, (p - 1) / 2, p) == 1 {
                2
            } else {
                0
            }
        })
        .sum::<u64>()
}

fn identity() -> Point {
    Point {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    }
}

fn negate(f: Field, p: Point) -> Point {
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
    let lambda = if p == q {
        let numer = f.add(f.mul_scalar(f.mul(p.x, p.x), 3), a);
        f.mul(numer, f.inv(f.mul_scalar(p.y, 2)))
    } else {
        f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)))
    };
    let x = f.sub(f.sub(f.mul(lambda, lambda), p.x), q.x);
    let y = f.sub(f.mul(lambda, f.sub(p.x, x)), p.y);
    Point {
        x,
        y,
        infinity: false,
    }
}

fn on_curve(f: Field, a: E, b: E, p: Point) -> bool {
    p.infinity || f.mul(p.y, p.y) == f.add(f.add(f.mul(f.mul(p.x, p.x), p.x), f.mul(a, p.x)), b)
}

fn velu(f: Field, a: E, t: Point, p: Point) -> Point {
    if p.infinity || p == t || p == negate(f, t) {
        return identity();
    }
    let plus = add(f, a, p, t);
    let minus = add(f, a, p, negate(f, t));
    assert!(!plus.infinity && !minus.infinity);
    Point {
        x: f.sub(f.add(f.add(p.x, plus.x), minus.x), f.mul_scalar(t.x, 2)),
        y: f.add(f.add(p.y, plus.y), minus.y),
        infinity: false,
    }
}

fn point(v: &Value) -> Point {
    let pair =
        |field: &str| -> E { [v[field][0].as_u64().unwrap(), v[field][1].as_u64().unwrap()] };
    Point {
        x: pair("x"),
        y: pair("y"),
        infinity: v["infinity"].as_bool().unwrap(),
    }
}

fn verify_row(p: u64, row: &Value) -> (bool, usize) {
    assert_eq!(row["schema"], "different_j_order_screen/v1");
    assert_eq!(row["p"].as_u64(), Some(p));
    let f = Field {
        p,
        nu: first_nonsquare(p),
    };
    let a = f.scalar(-21);
    let b = f.scalar(37);
    let da = f.scalar(-81);
    let db = f.scalar(-243);
    let t = Point {
        x: [3, 0],
        y: [1, 0],
        infinity: false,
    };
    let js = j(a, b, p);
    let reason = if js.is_none() {
        Some("source_singular")
    } else if !on_curve(f, [a, 0], [b, 0], t) || !add(f, [a, 0], add(f, [a, 0], t, t), t).infinity {
        Some("kernel_check")
    } else if j(da, db, p).is_none() {
        Some("destination_singular")
    } else if js == j(da, db, p) {
        Some("equal_j")
    } else {
        None
    };
    if let Some(reason) = reason {
        assert_eq!(row["status"], "reject");
        assert_eq!(row["reason"], reason);
        return (false, 0);
    }
    let n1 = count_fp(a, b, p);
    let n1_dst = count_fp(da, db, p);
    if n1 != n1_dst {
        assert_eq!(row["status"], "reject");
        assert_eq!(row["reason"], "base_order_mismatch");
        assert_eq!(row["n1"].as_u64(), Some(n1));
        assert_eq!(row["n1_dst"].as_u64(), Some(n1_dst));
        return (false, 0);
    }
    let trace = p as i64 + 1 - n1 as i64;
    let n2 = ((p + 1) * (p + 1)) as i64 - trace * trace;
    assert!(n2 > 0);
    let n2 = n2 as u64;
    let factors = factors(n2);
    let r = factors
        .iter()
        .rev()
        .find(|(prime, power)| *prime >= 211 && *prime != 3 && *power == 1)
        .map(|x| x.0);
    if r.is_none() {
        assert_eq!(row["status"], "reject");
        assert_eq!(row["reason"], "subgroup");
        assert_eq!(row["n1"].as_u64(), Some(n1));
        assert_eq!(row["n2"].as_u64(), Some(n2));
        assert_eq!(row["factors"], json!(factors));
        return (false, 0);
    }
    let r = r.unwrap();
    assert_eq!(row["status"], "accept", "unexpected map rejection at p={p}");
    assert_eq!(row["field_nu"].as_u64(), Some(f.nu));
    assert_eq!(
        row["source"],
        json!({"a":[a,0],"b":[b,0],"j":[js.unwrap(),0]})
    );
    assert_eq!(
        row["destination"],
        json!({"a":[da,0],"b":[db,0],"j":[j(da,db,p).unwrap(),0]})
    );
    assert_eq!(point(&row["kernel"][0]), t);
    assert_eq!(point(&row["kernel"][1]), negate(f, t));
    assert_eq!(row["n1"].as_u64(), Some(n1));
    assert_eq!(row["n1_dst"].as_u64(), Some(n1_dst));
    assert_eq!(row["trace"].as_i64(), Some(trace));
    assert_eq!(row["n2"].as_u64(), Some(n2));
    assert_eq!(row["factors"], json!(factors));
    assert_eq!(row["r"].as_u64(), Some(r));
    assert_eq!(row["cofactor"].as_u64(), Some(n2 / r));
    assert_eq!(row["sample_seed"].as_u64(), Some(0x2026_1015 ^ p));
    assert_eq!(velu(f, [a, 0], t, t), identity());
    assert_eq!(velu(f, [a, 0], t, negate(f, t)), identity());
    let samples = row["samples"].as_array().unwrap();
    assert_eq!(samples.len(), 64);
    for sample in samples {
        let u = point(&sample["a"]);
        let v = point(&sample["b"]);
        let sum = point(&sample["sum"]);
        let mu = point(&sample["ma"]);
        let mv = point(&sample["mb"]);
        let msum = point(&sample["msum"]);
        assert!(on_curve(f, [a, 0], [b, 0], u));
        assert!(on_curve(f, [a, 0], [b, 0], v));
        assert_eq!(add(f, [a, 0], u, v), sum);
        assert_eq!(velu(f, [a, 0], t, u), mu);
        assert_eq!(velu(f, [a, 0], t, v), mv);
        assert_eq!(velu(f, [a, 0], t, sum), msum);
        assert!(on_curve(f, [da, 0], [db, 0], mu));
        assert!(on_curve(f, [da, 0], [db, 0], mv));
        assert_eq!(add(f, [da, 0], mu, mv), msum);
    }
    (true, samples.len())
}

fn main() {
    let path = std::env::args().nth(1).expect("screen JSONL path");
    let raw = fs::read_to_string(&path).expect("screen JSONL");
    let mut lines = raw
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap());
    let mut accepted = 0usize;
    let mut witnesses = 0usize;
    let mut selected = None;
    let mut candidates = 0usize;
    for p in (61..=271).filter(|p| p % 3 == 1 && prime(*p)) {
        let row = lines.next().expect("missing candidate row");
        let (yes, checked) = verify_row(p, &row);
        candidates += 1;
        witnesses += checked;
        if yes {
            accepted += 1;
            if selected.is_none() {
                selected = Some(json!({"p":p,"r":row["r"]}));
            }
        }
    }
    let summary = lines.next().expect("missing summary row");
    assert!(lines.next().is_none(), "extra screen output rows");
    assert_eq!(summary["schema"], "different_j_order_screen_summary/v1");
    assert_eq!(summary["candidates"].as_u64(), Some(candidates as u64));
    assert_eq!(summary["accepted"].as_u64(), Some(accepted as u64));
    assert_eq!(summary["selected"], selected.clone().unwrap_or(Value::Null));
    assert_eq!(summary["max_p"].as_u64(), Some(271));
    assert_eq!(summary["min_r"].as_u64(), Some(211));
    println!(
        "{}",
        json!({"schema":"different_j_order_check/v1","status":"pass",
        "candidates":candidates,"accepted":accepted,"selected":selected,
        "map_witnesses":witnesses,"bytes":raw.len()})
    );
}
