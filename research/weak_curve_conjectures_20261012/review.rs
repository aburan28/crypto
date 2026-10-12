//! Fixed published-census metadata review; no external curve input.
#[allow(dead_code)]
#[cfg(not(test))]
#[path = "../../src/hash/sha256.rs"]
mod sha256;
#[cfg(not(test))]
use std::{collections::BTreeSet, fs};
#[cfg(not(test))]
const PRIMES: [i128; 7] = [7, 11, 13, 17, 19, 29, 37];
#[derive(Debug)]
struct Row {
    p: i128,
    t: i128,
    c: bool,
    q: bool,
    u: bool,
}
fn bit(s: &str) -> bool {
    match s {
        "0" => false,
        "1" => true,
        _ => panic!("invalid label"),
    }
}
fn row(line: &str) -> Row {
    let fields: Vec<_> = line.split(',').collect();
    assert_eq!(fields.len(), 5);
    let r = Row {
        p: fields[0].parse().unwrap(),
        t: fields[1].parse().unwrap(),
        c: bit(fields[2]),
        q: bit(fields[3]),
        u: bit(fields[4]),
    };
    assert_eq!(r.u, r.c || r.q, "union label mismatch");
    r
}
fn depth(p: i128, t: i128) -> u32 {
    let delta = t * t - 4 * p.pow(6);
    assert!(delta < 0);
    let s = delta.unsigned_abs().trailing_zeros();
    let odd = delta / (1_i128 << s);
    let fundamental_depth = if s % 2 == 1 {
        3
    } else if odd.rem_euclid(4) == 1 {
        0
    } else {
        2
    };
    assert!(s >= fundamental_depth && (s - fundamental_depth) % 2 == 0);
    (s - fundamental_depth) / 2
}
#[cfg(not(test))]
fn main() {
    assert_eq!(std::env::args().count(), 1, "fixed inputs only");
    let mut table = String::from("p\tclasses\tcubic\tquadratic\tboth\tcubic_only\tquadratic_only\tunion_zero\tcubic_depth_admissible\tcubic_high_depth_zero\tmax_zero_boundary_gap\tzeros_outside_64p\n");
    let mut failures = String::from("claim\tp\ttrace\tboundary_gap\tlimit\n");
    let mut hashes = String::new();
    let mut reviewed = 0_usize;
    let mut total_zero = 0_usize;
    let mut total_outside = 0_usize;
    for p in PRIMES {
        let path = format!("../iso1_weak_classes_20261007/two_branch_20261009/p{p}_classes.csv");
        let bytes = fs::read(&path).unwrap();
        let hash: String = sha256::sha256(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        hashes.push_str(&format!("{hash}\t{path}\n"));
        let input = String::from_utf8(bytes).unwrap();
        let mut lines = input.lines();
        assert_eq!(lines.next(), Some("p,trace,cubic,quadratic,combined"));
        let bound = 2 * p.pow(3);
        let mut seen = BTreeSet::new();
        let (mut c, mut q, mut both, mut zero, mut admissible, mut high_zero, mut max_gap, mut out) = (
            0_usize, 0_usize, 0_usize, 0_usize, 0_usize, 0_usize, 0_i128, 0_usize,
        );
        for line in lines {
            let r = row(line);
            assert_eq!(r.p, p);
            assert!(r.t.abs() < bound && r.t.rem_euclid(4) == 2 && r.t % p != 0);
            assert!(seen.insert(r.t), "duplicate signed trace");
            let d = depth(p, r.t);
            assert!(!r.c || d >= 2, "cubic depth obstruction violated");
            c += usize::from(r.c);
            q += usize::from(r.q);
            both += usize::from(r.c && r.q);
            admissible += usize::from(d >= 2);
            high_zero += usize::from(d >= 2 && !r.c);
            if !r.u {
                zero += 1;
                let gap = bound - r.t.abs();
                max_gap = max_gap.max(gap);
                if gap > 64 * p {
                    out += 1;
                    failures.push_str(&format!(
                        "boundary_strip_64p\t{p}\t{}\t{gap}\t{}\n",
                        r.t,
                        64 * p
                    ));
                }
            }
        }
        let mut expected = BTreeSet::new();
        let mut t = -bound;
        while t <= bound {
            if t % p != 0 {
                expected.insert(t);
            }
            t += 4;
        }
        assert_eq!(seen, expected, "incomplete ordinary signed trace domain");
        let n = seen.len();
        assert_eq!(n as i128, p * p * (p - 1));
        assert_eq!(c + q, n - zero + both);
        reviewed += n;
        total_zero += zero;
        total_outside += out;
        table.push_str(&format!("{p}\t{n}\t{c}\t{q}\t{both}\t{}\t{}\t{zero}\t{admissible}\t{high_zero}\t{max_gap}\t{out}\n",c-both,q-both));
    }
    fs::write("REVIEW.tsv", table).unwrap();
    fs::write("FALSIFIERS.tsv", failures).unwrap();
    fs::write("INPUTS_SHA256.tsv", hashes).unwrap();
    println!("FIXED_CENSUS_REVIEW_COMPLETE|primes=7|rows={reviewed}|union_zero={total_zero}|strip_falsifiers={total_outside}");
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[should_panic(expected = "union label mismatch")]
    fn rejects_inconsistent_union() {
        row("7,-682,0,0,1");
    }
    #[test]
    #[should_panic(expected = "invalid label")]
    fn rejects_unknown_label_as_boolean() {
        row("7,-682,?,0,0");
    }
    #[test]
    fn counterexample_depth_agrees_with_discriminant() {
        assert_eq!(depth(7, 22), 1);
    }
}
