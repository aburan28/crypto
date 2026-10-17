//! Direct F_(p^2) point census for each accepted different-j screen pair.
//! Uses independent small-field arithmetic and quadratic characters.

use std::fs;
use std::time::Instant;

use serde_json::{json, Value};

type E = [u64; 2];

#[derive(Clone, Copy)]
struct F {
    p: u64,
    nu: u64,
}

impl F {
    fn add(self, a: E, b: E) -> E {
        [(a[0] + b[0]) % self.p, (a[1] + b[1]) % self.p]
    }

    fn mul(self, a: E, b: E) -> E {
        let p = self.p as u128;
        let x = (a[0] as u128 * b[0] as u128 + self.nu as u128 * a[1] as u128 * b[1] as u128) % p;
        let y = (a[0] as u128 * b[1] as u128 + a[1] as u128 * b[0] as u128) % p;
        [x as u64, y as u64]
    }

    fn pow(self, mut base: E, mut exp: u64) -> E {
        let mut result = [1, 0];
        while exp != 0 {
            if exp & 1 != 0 {
                result = self.mul(result, base);
            }
            base = self.mul(base, base);
            exp >>= 1;
        }
        result
    }

    fn count(self, a: u64, b: u64) -> u64 {
        let q = self.p * self.p;
        let mut count = 1u64;
        for x0 in 0..self.p {
            for x1 in 0..self.p {
                let x = [x0, x1];
                let x3 = self.mul(self.mul(x, x), x);
                let rhs = self.add(self.add(x3, self.mul([a, 0], x)), [b, 0]);
                if rhs == [0, 0] {
                    count += 1;
                } else {
                    let sign = self.pow(rhs, (q - 1) / 2);
                    assert!(sign == [1, 0] || sign == [self.p - 1, 0]);
                    if sign == [1, 0] {
                        count += 2;
                    }
                }
            }
        }
        count
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("frozen screen JSONL");
    let raw = fs::read_to_string(path).expect("screen JSONL");
    let start = Instant::now();
    let mut accepted = 0usize;
    for line in raw.lines() {
        let row: Value = serde_json::from_str(line).unwrap();
        if row["schema"] != "different_j_order_screen/v1" || row["status"] != "accept" {
            continue;
        }
        let p = row["p"].as_u64().unwrap();
        let nu = row["field_nu"].as_u64().unwrap();
        let f = F { p, nu };
        assert_eq!(f.pow([nu, 0], (p - 1) / 2), [p - 1, 0]);
        let src_a = row["source"]["a"][0].as_u64().unwrap();
        let src_b = row["source"]["b"][0].as_u64().unwrap();
        let dst_a = row["destination"]["a"][0].as_u64().unwrap();
        let dst_b = row["destination"]["b"][0].as_u64().unwrap();
        let source_points = f.count(src_a, src_b);
        let destination_points = f.count(dst_a, dst_b);
        let expected = row["n2"].as_u64().unwrap();
        assert_eq!(source_points, expected, "source order at p={p}");
        assert_eq!(destination_points, expected, "destination order at p={p}");
        accepted += 1;
        println!(
            "{}",
            json!({"schema":"different_j_order_census/v1",
            "p":p,"field_nu":nu,"source_points":source_points,
            "destination_points":destination_points,"expected_n2":expected,
            "r":row["r"],"status":"pass"})
        );
        assert!(start.elapsed().as_secs() < 120, "census wall cap exhausted");
    }
    assert_eq!(accepted, 3, "accepted pair count changed");
    println!(
        "{}",
        json!({"schema":"different_j_order_census_summary/v1",
        "accepted":accepted,"status":"pass"})
    );
}
