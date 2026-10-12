//! Validate census totals and independent PARI trace controls.
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

#[path = "../../../../src/hash/sha256.rs"]
mod sha256_source;

#[derive(Debug, Clone)]
struct Row {
    p: u64,
    trace: i128,
    count: u64,
    depth: Option<u32>,
}

fn parse(input: &str) -> Vec<Row> {
    let mut lines = input.lines();
    assert_eq!(lines.next().unwrap(), "p,q,trace,weak_representatives,fundamental_discriminant,frobenius_2_depth,two_split,v2_delta_quarter,class_number_odd,trace_status");
    lines
        .map(|line| {
            let c: Vec<_> = line.split(',').collect();
            assert_eq!(c.len(), 10);
            let p: u64 = c[0].parse().unwrap();
            assert_eq!(c[1].parse::<u64>().unwrap(), p * p);
            let t: i128 = c[2].parse().unwrap();
            let depth = if c[9] == "ordinary" {
                assert_ne!(t % p as i128, 0);
                let dk: i128 = c[4].parse().unwrap();
                assert!(dk < 0);
                let delta = t * t - 4 * (p as i128).pow(6);
                let quotient = delta / dk;
                assert_eq!(delta % dk, 0);
                let mut f = (quotient as f64).sqrt() as i128;
                while f * f < quotient {
                    f += 1;
                }
                while f * f > quotient {
                    f -= 1;
                }
                assert_eq!(f * f, quotient);
                let d: u32 = c[5].parse().unwrap();
                assert_eq!(f.trailing_zeros(), d);
                assert_eq!(
                    (-(delta / 4) as u128).trailing_zeros(),
                    c[7].parse::<u32>().unwrap()
                );
                let split = if dk % 2 == 0 {
                    0
                } else if dk.rem_euclid(8) == 1 {
                    1
                } else {
                    -1
                };
                assert_eq!(split, c[6].parse::<i8>().unwrap());
                Some(d)
            } else {
                assert_eq!(c[9], "nonordinary_or_unrealized");
                None
            };
            Row {
                p,
                trace: t,
                count: c[3].parse().unwrap(),
                depth,
            }
        })
        .collect()
}

fn audit(rows: &[Row]) -> (u64, u64, u64, u64, u64) {
    let p = rows[0].p;
    let lim = 2 * (p as i128).pow(3);
    assert_eq!(rows.len() as i128, lim / 2 + 1);
    let labels: BTreeMap<_, _> = rows.iter().map(|r| (r.trace, r.count)).collect();
    assert_eq!(labels.len(), rows.len());
    let total: u64 = rows.iter().map(|r| r.count).sum();
    assert_eq!(total, 2 * p.pow(4) + 2 * p * p);
    for (i, r) in rows.iter().enumerate() {
        assert_eq!(r.p, p);
        assert_eq!(r.trace, -lim + 4 * i as i128);
        assert_eq!(r.count, *labels.get(&-r.trace).unwrap());
        if r.count > 0 {
            let v = ((p as i128).pow(6) + 1).rem_euclid(16);
            assert!(r.trace.rem_euclid(16) == v || (-r.trace).rem_euclid(16) == v);
            if let Some(d) = r.depth {
                assert!(d >= 2);
            }
        }
    }
    let ordinary = rows.iter().filter(|r| r.depth.is_some()).count() as u64;
    let weak = rows
        .iter()
        .filter(|r| r.depth.is_some() && r.count > 0)
        .count() as u64;
    let low = rows.iter().filter(|r| r.depth == Some(1)).count() as u64;
    let high_zero = rows
        .iter()
        .filter(|r| r.depth.is_some_and(|d| d >= 2) && r.count == 0)
        .count() as u64;
    (total, ordinary, weak, low, high_zero)
}

fn digest(path: &Path) -> String {
    sha256_source::sha256(&std::fs::read(path).unwrap())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        args.len() == 1 || args.len() == 3,
        "usage: iso1_census_audit CENSUS.csv [samples|orbits GP.csv]"
    );
    let path = Path::new(&args[0]);
    let rows = parse(&std::fs::read_to_string(path).unwrap());
    let (total, ordinary, weak, low, high_zero) = audit(&rows);
    println!("p={} rows={} normalized_representatives={} ordinary={} weak={} depth1_weak=0 depth1_rows={} high_depth_zero={} high_depth_rows={}",rows[0].p,rows.len(),total,ordinary,weak,low,high_zero,ordinary-low);
    println!("census_sha256={}", digest(path));
    if args.len() == 3 {
        let control = Path::new(&args[2]);
        let input = std::fs::read_to_string(control).unwrap();
        let positive: BTreeMap<_, _> = rows
            .iter()
            .filter(|r| r.trace > 0)
            .map(|r| (r.trace, r.count))
            .collect();
        let mut observed = BTreeMap::<i128, u64>::new();
        let mut histogram = BTreeMap::<u64, u64>::new();
        let mut distinct = BTreeSet::new();
        let mut samples = 0;
        for line in input.lines() {
            let c: Vec<_> = line.split(',').collect();
            assert_eq!(c.len(), 2);
            let (trace, weight) = match args[1].as_str() {
                "samples" => {
                    assert_eq!(c[0].parse::<u64>().unwrap(), rows[0].p);
                    (c[1].parse::<i128>().unwrap().abs(), 1)
                }
                "orbits" => (c[0].parse::<i128>().unwrap().abs(), c[1].parse().unwrap()),
                _ => panic!("unknown control mode"),
            };
            assert!(
                positive.get(&trace).is_some_and(|v| *v > 0),
                "GP trace {trace} has no native positive label"
            );
            *observed.entry(trace).or_default() += weight;
            *histogram.entry(weight).or_default() += 1;
            distinct.insert(trace);
            samples += 1;
        }
        assert!(samples > 0);
        if args[1] == "orbits" {
            let p = rows[0].p;
            assert_eq!(samples, (p.pow(4) + 3 * p * p + 8) / 12);
            assert_eq!(histogram.get(&2), Some(&1));
            assert_eq!(histogram.get(&6), Some(&((p * p - 1) / 3)));
            assert_eq!(histogram.get(&12), Some(&((p.pow(4) - p * p) / 12)));
            for (&trace, &weight) in &positive {
                assert_eq!(
                    weight,
                    *observed.get(&trace).unwrap_or(&0),
                    "trace-weight mismatch at {trace}"
                );
            }
        }
        println!(
            "GP_control={} samples={} distinct_abs_traces={} missing_positive=0 histogram={:?}",
            args[1],
            samples,
            distinct.len(),
            histogram
        );
        println!("control_sha256={}", digest(control));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_p13_census_and_independent_orbits_agree() {
        let study = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let rows = parse(&std::fs::read_to_string(study.join("p13_twist_derived.csv")).unwrap());
        assert_eq!(audit(&rows), (57460, 2028, 928, 1014, 86));
    }
    #[test]
    fn shifted_trace_weight_is_rejected() {
        let study = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let mut rows =
            parse(&std::fs::read_to_string(study.join("p13_twist_derived.csv")).unwrap());
        let i = rows.iter().position(|r| r.count > 0).unwrap();
        rows[i].count += 1;
        assert!(std::panic::catch_unwind(|| audit(&rows)).is_err());
    }
}
