//! Post-hoc replay reconstructs the base from public parameters, not backend files.
use super::{read, write, ROOT};
use crypto_lib::cryptanalysis::{
    ecbench::canonical::sha256_hex,
    ic_boundary::{CountedGroup, GroupOps, IncrementalGauss, PrimeCurve, PrimePoint, RowStatus},
};
use serde_json::{json, Value};
use std::{collections::BTreeSet, fs, path::Path};
fn u(v: &Value, k: &str) -> u64 {
    v[k].as_u64().unwrap()
}
fn mul(a: u64, b: u64, r: u64) -> u64 {
    (a as u128 * b as u128 % r as u128) as u64
}
pub fn run() {
    let root = Path::new(ROOT);
    let snapshot = read(root.join("results/campaign-progress.json"));
    let mut receipts = Vec::new();
    for summary in snapshot["records"].as_array().unwrap() {
        if summary["arm"] == "rho" || summary["status"] == "error" {
            continue;
        }
        let text = fs::read_to_string(
            root.join("sessions")
                .join(summary["session"].as_str().unwrap())
                .join("records.jsonl"),
        )
        .unwrap();
        let v: Value = text
            .lines()
            .map(|l| serde_json::from_str::<Value>(l).unwrap())
            .find(|v| v["record_id"] == summary["record_id"])
            .unwrap();
        let spec = &v["workload"]["curve_spec"];
        let c = PrimeCurve {
            p: u(spec, "p"),
            a: u(spec, "a"),
            b: u(spec, "b"),
        };
        let r = u(spec, "r");
        let g = PrimePoint::affine(u(spec, "gx"), u(spec, "gy"));
        let parse = |i: usize| {
            u64::from_str_radix(
                v["workload"]["target"][i]
                    .as_str()
                    .unwrap()
                    .trim_start_matches("0x"),
                16,
            )
            .unwrap()
        };
        let q = PrimePoint::affine(parse(0), parse(1));
        let d = &v["detail"];
        let lambda = u(d, "lambda");
        let degree = if v["method"]["id"] == "ic.cm_folded_native_d43_v1" {
            11
        } else {
            5
        };
        assert_eq!(
            (mul(lambda, lambda, r) as u128 + r as u128 - lambda as u128 + degree) % r as u128,
            0
        );
        let size = u(&d["native"], "B") as usize;
        let length = size / 4;
        let mut ops = GroupOps::default();
        let mut base = Vec::new();
        let mut x = 1;
        let mut csv = String::from("point,column,coefficient,x,y\n");
        let mut abscissae = BTreeSet::new();
        for column in 0..4 {
            let seed = loop {
                let cand = x;
                x += 1;
                if cand == g.x {
                    continue;
                }
                if let Some(y) = c.sqrt(c.rhs(cand)) {
                    if y != 0 {
                        break PrimePoint::affine(cand, y.min(c.p - y));
                    }
                }
            };
            let mut point = seed;
            let mut coefficient = 1;
            for level in 0..length {
                assert!(c.is_on_curve(point) && abscissae.insert(point.x));
                csv += &format!(
                    "{},{column},{coefficient},{},{}\n",
                    base.len(),
                    point.x,
                    point.y
                );
                base.push((point, column, coefficient));
                if level + 1 < length {
                    point = c.mul(&mut ops, point, lambda);
                    coefficient = mul(coefficient, lambda, r);
                    if point.y > c.p - point.y {
                        point = c.neg(point);
                        coefficient = r - coefficient;
                    }
                }
            }
        }
        assert_eq!(sha256_hex(csv.as_bytes()), d["replay"]["base_sha256"]);
        let mut matrix = IncrementalGauss::new(5, r);
        let mut normalized = BTreeSet::new();
        let mut mutation_rejected = false;
        let seed: u64 = v["algorithm_seed"].as_str().unwrap().parse().unwrap();
        let queries = u(&d["native"], "queries");
        for rel in d["relations"].as_array().unwrap() {
            let mut row = vec![0; 5];
            let mut sum = PrimePoint::INFINITY;
            for (key, signkey) in [("i", "si"), ("j", "sj")] {
                let (point, col, coef) = base[u(rel, key) as usize];
                let sign = rel[signkey].as_i64().unwrap();
                assert!(sign == 1 || sign == -1);
                let signed = if sign == 1 { point } else { c.neg(point) };
                sum = c.add(&mut ops, sum, signed);
                row[col] = ((row[col] as u128 + if sign == 1 { coef } else { r - coef } as u128)
                    % r as u128) as u64;
            }
            let alpha = u(rel, "alpha");
            assert!(
                (alpha as u128 + r as u128 - (seed % r) as u128) % (r as u128) < queries as u128
            );
            let ag = c.mul(&mut ops, g, alpha);
            let expected = c.add(&mut ops, ag, q);
            assert_eq!(sum, expected);
            assert_eq!(sum, PrimePoint::affine(u(rel, "Rx"), u(rel, "Ry")));
            if !mutation_rejected {
                let p = base[u(rel, "i") as usize].0;
                let twice = c.mul(&mut ops, p, 2);
                let changed = c.add(
                    &mut ops,
                    sum,
                    if rel["si"] == 1 { c.neg(twice) } else { twice },
                );
                assert_ne!(changed, expected);
                mutation_rejected = true;
            }
            let first = *row[..4]
                .iter()
                .find(|&&x| x != 0)
                .expect("zero factor-base column part");
            let mut inverse = 1;
            let mut b = first;
            let mut e = r - 2;
            while e > 0 {
                if e & 1 != 0 {
                    inverse = mul(inverse, b, r);
                }
                b = mul(b, b, r);
                e >>= 1;
            }
            let profile: Vec<_> = row[..4].iter().map(|&x| mul(x, inverse, r)).collect();
            assert!(normalized.insert(profile));
            row[4] = r - 1;
            assert!(!matches!(
                matrix.add_row(row, alpha),
                RowStatus::Inconsistent
            ));
            assert_eq!(matrix.rank() as u64, u(rel, "rank"));
        }
        if v["outcome"]["status"] == "verified" {
            let recovered: u64 = v["outcome"]["recovered"].as_str().unwrap().parse().unwrap();
            assert_eq!(matrix.pinned(4), Some(recovered));
            assert_eq!(c.mul(&mut ops, g, recovered), q);
        }
        receipts.push(json!({"record_id":v["record_id"],"base_points":size,"base_sha256":sha256_hex(csv.as_bytes()),"relations":d["relations"].as_array().unwrap().len(),"rank":matrix.rank(),"mutated_sign_rejected":mutation_rejected,"status":"independent-public-base-and-relation-replay-passed"}));
    }
    write(
        root.join("results/posthoc-replay.json"),
        &json!({"schema":"conductor-gap.posthoc/v1","snapshot_sha256":sha256_hex(&fs::read(root.join("results/campaign-progress.json")).unwrap()),"records":receipts}),
    );
}
