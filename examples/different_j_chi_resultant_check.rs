//! Independent prime-field and curve replay of the χ-resultant pair sets.

use std::collections::BTreeSet;
use std::fs;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

type El = [u64; 2];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Point {
    x: El,
    y: El,
    infinity: bool,
}

#[derive(Clone, Copy)]
struct Field {
    p: u64,
    nu: u64,
}

impl Field {
    fn add(self, x: El, y: El) -> El {
        [(x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p]
    }
    fn sub(self, x: El, y: El) -> El {
        [
            (x[0] + self.p - y[0]) % self.p,
            (x[1] + self.p - y[1]) % self.p,
        ]
    }
    fn mul(self, x: El, y: El) -> El {
        let p = self.p as u128;
        let real =
            (x[0] as u128 * y[0] as u128 + self.nu as u128 * x[1] as u128 * y[1] as u128) % p;
        let imag = (x[0] as u128 * y[1] as u128 + x[1] as u128 * y[0] as u128) % p;
        [real as u64, imag as u64]
    }
    fn scale(self, x: El, scalar: u64) -> El {
        self.mul(x, [scalar % self.p, 0])
    }
    fn inv(self, x: El) -> El {
        let norm = self.sub(
            self.mul([x[0], 0], [x[0], 0]),
            self.scale([x[1] * x[1] % self.p, 0], self.nu),
        )[0];
        assert_ne!(norm, 0);
        let inv = pow_mod(norm, self.p - 2, self.p);
        self.mul([x[0], (self.p - x[1]) % self.p], [inv, 0])
    }
    fn neg(self, x: El) -> El {
        self.sub([0, 0], x)
    }
}

fn pow_mod(mut x: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = ((result as u128 * x as u128) % modulus as u128) as u64;
        }
        x = ((x as u128 * x as u128) % modulus as u128) as u64;
        exponent >>= 1;
    }
    result
}

fn identity() -> Point {
    Point {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    }
}

fn neg(f: Field, point: Point) -> Point {
    Point {
        y: f.neg(point.y),
        ..point
    }
}

fn add(f: Field, a: El, p: Point, q: Point) -> Point {
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
        let numerator = f.add(f.scale(f.mul(p.x, p.x), 3), a);
        f.mul(numerator, f.inv(f.scale(p.y, 2)))
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

fn scalar_mul(f: Field, a: El, point: Point, mut scalar: u64) -> Point {
    let mut result = identity();
    let mut base = point;
    while scalar != 0 {
        if scalar & 1 != 0 {
            result = add(f, a, result, base);
        }
        base = add(f, a, base, base);
        scalar >>= 1;
    }
    result
}

fn on_curve(f: Field, a: El, b: El, p: Point) -> bool {
    p.infinity || f.mul(p.y, p.y) == f.add(f.add(f.mul(f.mul(p.x, p.x), p.x), f.mul(a, p.x)), b)
}

fn point(v: &Value) -> Point {
    if v["infinity"] == true {
        return identity();
    }
    Point {
        x: [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        y: [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
        infinity: false,
    }
}

fn expected_pairs(f: Field, a: El, signed: &[Point], target: Point) -> BTreeSet<(usize, usize)> {
    let mut out = BTreeSet::new();
    for (i, &p) in signed.iter().enumerate() {
        for (j, &q) in signed.iter().enumerate().skip(i) {
            if add(f, a, p, q) == target {
                out.insert((i, j));
            }
        }
    }
    out
}

fn check_case(f: Field, a: El, signed: &[Point], target: Point, label: Value, row: Value) -> usize {
    assert_eq!(row["schema"], "different_j_chi_resultant_case/v1");
    assert_eq!(row["scalar_label"], label);
    assert_eq!(point(&row["target"]), target);
    assert_eq!(row["status"], "pass");
    let reported: BTreeSet<_> = row["pairs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| {
            let i = p[0].as_u64().unwrap() as usize;
            let j = p[1].as_u64().unwrap() as usize;
            assert!(i <= j && j < signed.len());
            (i, j)
        })
        .collect();
    assert_eq!(reported.len(), row["pairs"].as_array().unwrap().len());
    assert_eq!(reported, expected_pairs(f, a, signed, target));
    reported.len()
}

fn check_s3(public_path: &str, raw_path: &str) {
    let public_bytes = fs::read(public_path).unwrap();
    let public: Value = serde_json::from_slice(&public_bytes).unwrap();
    assert_eq!(public["p"], 211);
    assert_eq!(public["r"], 229);
    let f = Field { p: 211, nu: 2 };
    let source_a = [190, 0];
    let quotient_a = [130, 0];
    let g = point(&public["generator"]);
    let mapped_g = point(&public["mapped_generator"]);
    let t = Point {
        x: [3, 0],
        y: [1, 0],
        infinity: false,
    };
    let mut signed = Vec::new();
    for witness in public["full_base_witness_source"].as_array().unwrap() {
        let projected = scalar_mul(f, source_a, point(witness), 195);
        signed.push(projected);
        signed.push(neg(f, projected));
    }
    let mapped: Vec<Point> = signed
        .iter()
        .map(|&q| {
            if q.infinity || q == t || q == neg(f, t) {
                return identity();
            }
            let plus = add(f, source_a, q, t);
            let minus = add(f, source_a, q, neg(f, t));
            Point {
                x: f.sub(f.add(f.add(q.x, plus.x), minus.x), f.scale(t.x, 2)),
                y: f.add(f.add(q.y, plus.y), minus.y),
                infinity: false,
            }
        })
        .collect();
    let bytes = fs::read(raw_path).unwrap();
    let mut rows = bytes
        .split(|&byte| byte == b'\n')
        .filter(|row| !row.is_empty())
        .map(|row| serde_json::from_slice::<Value>(row).unwrap());
    let mut source_target = identity();
    let mut quotient_target = identity();
    let mut source_pairs = 0u64;
    let mut quotient_pairs = 0u64;
    let mut zero_source = 0u64;
    let mut zero_quotient = 0u64;
    for scalar in 0..229 {
        let row = rows.next().unwrap();
        assert_eq!(row["schema"], "different_j_chi_s3_case/v2");
        assert_eq!(row["scalar_label"], scalar);
        assert_eq!(point(&row["source_target"]), source_target);
        assert_eq!(point(&row["quotient_target"]), quotient_target);
        assert_eq!(row["status"], "pass");
        let expected_source = expected_pairs(f, source_a, &signed, source_target);
        let expected_quotient = expected_pairs(f, quotient_a, &mapped, quotient_target);
        let parse_pairs = |v: &Value| -> BTreeSet<(usize, usize)> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|pair| {
                    (
                        pair[0].as_u64().unwrap() as usize,
                        pair[1].as_u64().unwrap() as usize,
                    )
                })
                .collect()
        };
        let actual_source = parse_pairs(&row["source_pairs"]);
        let actual_quotient = parse_pairs(&row["quotient_pairs"]);
        assert_eq!(actual_source, expected_source);
        assert_eq!(actual_quotient, expected_quotient);
        assert_eq!(actual_source, actual_quotient);
        source_pairs += actual_source.len() as u64;
        quotient_pairs += actual_quotient.len() as u64;
        zero_source += u64::from(actual_source.is_empty());
        zero_quotient += u64::from(actual_quotient.is_empty());
        source_target = add(f, source_a, source_target, g);
        quotient_target = add(f, quotient_a, quotient_target, mapped_g);
    }
    assert_eq!(source_target, identity());
    assert_eq!(quotient_target, identity());
    let summary = rows.next().unwrap();
    assert!(rows.next().is_none());
    assert_eq!(summary["schema"], "different_j_chi_s3_exhaustive/v2");
    assert_eq!(summary["public_sha256"], sha256_hex(&public_bytes));
    assert_eq!(summary["source_cases"], 229);
    assert_eq!(summary["quotient_cases"], 229);
    assert_eq!(summary["source_pairs"], source_pairs);
    assert_eq!(summary["quotient_pairs"], quotient_pairs);
    assert_eq!(summary["source_zero_cases"], zero_source);
    assert_eq!(summary["quotient_zero_cases"], zero_quotient);
    assert_eq!(summary["status"], "pass");
    println!(
        "{}",
        json!({"schema":"different_j_chi_s3_check/v2",
        "status":"pass","source_cases":229,"quotient_cases":229,
        "source_pairs":source_pairs,"quotient_pairs":quotient_pairs,
        "source_zero_cases":zero_source,"quotient_zero_cases":zero_quotient,
        "public_sha256":sha256_hex(&public_bytes),
        "raw_sha256":sha256_hex(&bytes)})
    );
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("check-s3") {
        assert_eq!(args.len(), 4);
        check_s3(&args[2], &args[3]);
        return;
    }
    assert_eq!(args.len(), 3, "usage: checker PUBLIC.json EXHAUSTIVE.jsonl");
    let public_bytes = fs::read(&args[1]).unwrap();
    let public: Value = serde_json::from_slice(&public_bytes).unwrap();
    assert_eq!(public["p"], 211);
    assert_eq!(public["r"], 229);
    assert_eq!(public["cofactor"], 195);
    let f = Field { p: 211, nu: 2 };
    let a = [190, 0];
    let b = [37, 0];
    let g = point(&public["generator"]);
    assert!(on_curve(f, a, b, g));
    assert_eq!(scalar_mul(f, a, g, 229), identity());
    let mut signed = Vec::new();
    for row in public["full_base_witness_source"].as_array().unwrap() {
        let raw = point(row);
        assert!(on_curve(f, a, b, raw));
        let projected = scalar_mul(f, a, raw, 195);
        assert!(!projected.infinity && on_curve(f, a, b, projected));
        assert_eq!(scalar_mul(f, a, projected, 229), identity());
        signed.push(projected);
        signed.push(neg(f, projected));
    }
    assert_eq!(signed.len(), 16);
    let raw = fs::read_to_string(&args[2]).unwrap();
    let mut rows = raw
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap());
    let mut target = identity();
    let mut pair_count = 0;
    let mut zero_cases = 0;
    for scalar in 0..229 {
        let count = check_case(f, a, &signed, target, json!(scalar), rows.next().unwrap());
        pair_count += count;
        zero_cases += usize::from(count == 0);
        target = add(f, a, target, g);
    }
    assert_eq!(target, identity());
    let t = Point {
        x: [3, 0],
        y: [1, 0],
        infinity: false,
    };
    for (label, exceptional) in [("torsion_plus", t), ("torsion_minus", neg(f, t))] {
        let count = check_case(
            f,
            a,
            &signed,
            exceptional,
            json!(label),
            rows.next().unwrap(),
        );
        pair_count += count;
        zero_cases += usize::from(count == 0);
    }
    let summary = rows.next().unwrap();
    assert!(rows.next().is_none());
    assert_eq!(summary["schema"], "different_j_chi_resultant_summary/v1");
    assert_eq!(summary["public_sha256"], sha256_hex(&public_bytes));
    assert_eq!(summary["cases"], 231);
    assert_eq!(summary["pair_count"], pair_count);
    assert_eq!(summary["zero_cases"], zero_cases);
    assert_eq!(summary["signed_points"], 16);
    assert_eq!(summary["chi_polynomial_degree"], 16);
    assert_eq!(summary["status"], "pass");
    println!(
        "{}",
        json!({"schema":"different_j_chi_resultant_check/v1","status":"pass",
        "cases":231,"pairs":pair_count,"zero_cases":zero_cases,
        "public_sha256":sha256_hex(&public_bytes),"raw_sha256":sha256_hex(raw.as_bytes())})
    );
}
