//! A direct chi-coordinate chord resultant on the frozen different-j cell.
//! This correctness gate does not use a pair-sum table or scan base pairs
//! inside the polynomial oracle. It tests every point of the order-229 group.

use std::collections::BTreeSet;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ext_curve::{
    take_field_counters, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use serde_json::{json, Value};

type El = [u64; 2];
type Point = ExtPoint<El>;
type Curve = ExtCurve<Fp2>;
type Poly = Vec<El>;

fn trim(f: &Fp2, p: &mut Poly) {
    while p.last().is_some_and(|&c| c == f.zero()) {
        p.pop();
    }
}

fn add(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![f.zero(); a.len().max(b.len())];
    for (i, &v) in a.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    for (i, &v) in b.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    trim(f, &mut out);
    out
}

fn sub(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![f.zero(); a.len().max(b.len())];
    for (i, &v) in a.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    for (i, &v) in b.iter().enumerate() {
        out[i] = f.sub(out[i], v);
    }
    trim(f, &mut out);
    out
}

fn scale(f: &Fp2, a: &Poly, k: El) -> Poly {
    let mut out: Poly = a.iter().map(|&v| f.mul(v, k)).collect();
    trim(f, &mut out);
    out
}

fn mul(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![f.zero(); a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = f.add(out[i + j], f.mul(x, y));
        }
    }
    trim(f, &mut out);
    out
}

fn div_rem(f: &Fp2, a: &Poly, b: &Poly) -> (Poly, Poly) {
    assert!(!b.is_empty());
    let mut rem = a.clone();
    let mut quo = vec![f.zero(); a.len().saturating_sub(b.len()) + 1];
    let inv_lead = f.inv(*b.last().unwrap()).unwrap();
    while rem.len() >= b.len() {
        let shift = rem.len() - b.len();
        let coeff = f.mul(*rem.last().unwrap(), inv_lead);
        quo[shift] = coeff;
        for (j, &term) in b.iter().enumerate() {
            rem[shift + j] = f.sub(rem[shift + j], f.mul(coeff, term));
        }
        trim(f, &mut rem);
    }
    trim(f, &mut quo);
    (quo, rem)
}

fn gcd(f: &Fp2, mut a: Poly, mut b: Poly) -> Poly {
    while !b.is_empty() {
        let rem = div_rem(f, &a, &b).1;
        a = b;
        b = rem;
    }
    if let Some(&lead) = a.last() {
        let inverse = f.inv(lead).unwrap();
        for coeff in &mut a {
            *coeff = f.mul(*coeff, inverse);
        }
    }
    a
}

fn mul_mod(f: &Fp2, a: &Poly, b: &Poly, modulus: &Poly) -> Poly {
    div_rem(f, &mul(f, a, b), modulus).1
}

fn pow_mod(f: &Fp2, a: &Poly, mut exponent: u64, modulus: &Poly) -> Poly {
    let mut out = vec![f.one()];
    let mut base = div_rem(f, a, modulus).1;
    while exponent > 0 {
        if exponent & 1 != 0 {
            out = mul_mod(f, &out, &base, modulus);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = mul_mod(f, &base, &base, modulus);
        }
    }
    out
}

fn eval(f: &Fp2, a: &Poly, x: El) -> El {
    a.iter()
        .rev()
        .fold(f.zero(), |acc, &c| f.add(f.mul(acc, x), c))
}

fn root_multiplicity(f: &Fp2, a: &Poly, root: El) -> u64 {
    let mut q = a.clone();
    let mut count = 0;
    let linear = vec![f.neg(root), f.one()];
    while q.len() > 1 {
        let (next, rem) = div_rem(f, &q, &linear);
        if !rem.is_empty() {
            break;
        }
        q = next;
        count += 1;
    }
    count
}

#[derive(Default)]
struct Stats {
    resultants: u64,
    gcds: u64,
    split_attempts: u64,
    exhaustive_fallback_evaluations: u64,
    root_candidates: u64,
    chi_lift_scans: u64,
    exact_checks: u64,
    false_candidates: u64,
    vertical_cases: u64,
    tangent_cases: u64,
    horizontal_cases: u64,
    infinity_cases: u64,
    exceptional_chi_cases: u64,
    multiple_resultant_roots: u64,
}

fn all_roots(f: &Fp2, polynomial: &Poly, stats: &mut Stats) -> Vec<El> {
    if polynomial.len() <= 1 {
        return Vec::new();
    }
    let x = vec![f.zero(), f.one()];
    let rational = gcd(
        f,
        polynomial.clone(),
        sub(f, &pow_mod(f, &x, f.order(), polynomial), &x),
    );
    stats.gcds += 1;
    fn split(f: &Fp2, a: Poly, stats: &mut Stats, output: &mut Vec<El>) {
        if a.len() <= 1 {
            return;
        }
        if a.len() == 2 {
            output.push(f.neg(f.mul(a[0], f.inv(a[1]).unwrap())));
            return;
        }
        if a.len() == 3 {
            let four_ac = f.scale(f.mul(a[2], a[0]), 4);
            let disc = f.sub(f.sqr(a[1]), four_ac);
            if let Some(root_disc) = f.sqrt(disc) {
                let inv = f.inv(f.scale(a[2], 2)).unwrap();
                output.push(f.mul(f.add(f.neg(a[1]), root_disc), inv));
                output.push(f.mul(f.sub(f.neg(a[1]), root_disc), inv));
                return;
            }
        }
        if a.len() == 4 {
            for shift in 0..128u64 {
                stats.split_attempts += 1;
                let linear = vec![[shift % f.p, 1], f.one()];
                let character = pow_mod(f, &linear, (f.order() - 1) / 2, &a);
                for test in [
                    sub(f, &character, &vec![f.one()]),
                    add(f, &character, &vec![f.one()]),
                    character,
                ] {
                    let divisor = gcd(f, a.clone(), test);
                    stats.gcds += 1;
                    if divisor.len() > 1 && divisor.len() < a.len() {
                        let (quotient, rem) = div_rem(f, &a, &divisor);
                        assert!(rem.is_empty());
                        split(f, divisor, stats, output);
                        split(f, quotient, stats, output);
                        return;
                    }
                }
            }
        }
        // Complete, explicitly charged fallback for a failed deterministic split.
        for packed in 0..f.order() {
            let z = [packed % f.p, packed / f.p];
            stats.exhaustive_fallback_evaluations += 1;
            if eval(f, &a, z) == f.zero() {
                output.push(z);
            }
        }
    }
    let mut roots = Vec::new();
    split(f, rational, stats, &mut roots);
    roots.sort_unstable();
    roots.dedup();
    assert!(roots.iter().all(|&z| eval(f, polynomial, z) == f.zero()));
    roots
}

fn chi(curve: &Curve, point: Point) -> Option<El> {
    if point.infinity {
        return None;
    }
    let f = &curve.f;
    f.inv(f.sub(point.y, f.one()))
        .map(|inv| f.mul(f.add(point.y, f.one()), inv))
}

fn chi_base(curve: &Curve, signed: &[Point]) -> Poly {
    let f = &curve.f;
    let mut base = vec![f.one()];
    for &point in signed {
        if let Some(z) = chi(curve, point) {
            base = mul(f, &base, &vec![f.neg(z), f.one()]);
        }
    }
    base
}

fn chord_resultant(curve: &Curve, point: Point, slope: El) -> Poly {
    let f = &curve.f;
    let d = vec![f.neg(f.one()), f.one()];
    let n = vec![f.one(), f.one()];
    let m3 = f.mul(f.sqr(slope), slope);
    let l = add(f, &scale(f, &d, f.sub(f.mul(slope, point.x), point.y)), &n);
    let d2 = mul(f, &d, &d);
    let d3 = mul(f, &d2, &d);
    let l3 = mul(f, &mul(f, &l, &l), &l);
    let ax = scale(f, &mul(f, &d2, &l), f.mul(curve.a, f.sqr(slope)));
    let b = scale(f, &d3, f.mul(curve.b, m3));
    let y2 = scale(f, &mul(f, &d, &mul(f, &n, &n)), m3);
    sub(f, &add(f, &add(f, &l3, &ax), &b), &y2)
}

fn parse_point(v: &Value) -> Point {
    if v["infinity"] == true {
        return Point::infinity([0, 0]);
    }
    Point::affine(
        [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
    )
}

fn direct_pairs(
    curve: &Curve,
    signed: &[Point],
    target: Point,
    ops: &mut GroupOps,
) -> BTreeSet<(usize, usize)> {
    let mut out = BTreeSet::new();
    for (i, &a) in signed.iter().enumerate() {
        for (j, &b) in signed.iter().enumerate().skip(i) {
            if curve.add(ops, a, b) == target {
                out.insert((i, j));
            }
        }
    }
    out
}

fn polynomial_pairs(
    curve: &Curve,
    signed: &[Point],
    base: &Poly,
    target: Point,
    ops: &mut GroupOps,
    stats: &mut Stats,
) -> BTreeSet<(usize, usize)> {
    let f = &curve.f;
    let mut out = BTreeSet::new();
    for (i, &point) in signed.iter().enumerate() {
        if target.infinity {
            stats.infinity_cases += 1;
            for (j, &candidate) in signed.iter().enumerate() {
                if candidate == curve.neg(point) {
                    out.insert((i.min(j), i.max(j)));
                }
            }
            continue;
        }
        if point == target {
            stats.vertical_cases += 1; // The second summand would be infinity.
            continue;
        }
        let minus_target = curve.neg(target);
        let slope = if point == minus_target {
            stats.tangent_cases += 1;
            let numerator = f.add(f.scale(f.sqr(point.x), 3), curve.a);
            let Some(denominator) = f.inv(f.scale(point.y, 2)) else {
                stats.vertical_cases += 1;
                continue;
            };
            f.mul(numerator, denominator)
        } else {
            let Some(denominator) = f.inv(f.sub(minus_target.x, point.x)) else {
                stats.vertical_cases += 1;
                continue;
            };
            f.mul(f.sub(minus_target.y, point.y), denominator)
        };
        if slope == f.zero() {
            stats.horizontal_cases += 1;
            let z = chi(curve, point);
            for (j, &candidate) in signed.iter().enumerate() {
                stats.chi_lift_scans += 1;
                if chi(curve, candidate) == z {
                    stats.exact_checks += 1;
                    if curve.add(ops, point, candidate) == target {
                        out.insert((i.min(j), i.max(j)));
                    } else {
                        stats.false_candidates += 1;
                    }
                }
            }
            continue;
        }
        let resultant = chord_resultant(curve, point, slope);
        stats.resultants += 1;
        let intersection = gcd(f, resultant.clone(), base.clone());
        stats.gcds += 1;
        for root in all_roots(f, &intersection, stats) {
            stats.root_candidates += 1;
            if root_multiplicity(f, &resultant, root) > 1 {
                stats.multiple_resultant_roots += 1;
            }
            for (j, &candidate) in signed.iter().enumerate() {
                stats.chi_lift_scans += 1;
                if chi(curve, candidate) == Some(root) {
                    stats.exact_checks += 1;
                    if curve.add(ops, point, candidate) == target {
                        out.insert((i.min(j), i.max(j)));
                    } else {
                        stats.false_candidates += 1;
                    }
                }
            }
        }
        for (j, &candidate) in signed.iter().enumerate() {
            if chi(curve, candidate).is_none() {
                stats.exceptional_chi_cases += 1;
                stats.exact_checks += 1;
                if curve.add(ops, point, candidate) == target {
                    out.insert((i.min(j), i.max(j)));
                }
            }
        }
    }
    out
}

fn main() {
    assert_eq!(std::env::args().nth(1).as_deref(), Some("exhaustive"));
    let path = std::env::args().nth(2).expect("public input path");
    let bytes = fs::read(&path).unwrap();
    let public: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(public["p"], 211);
    assert_eq!(public["r"], 229);
    assert_eq!(public["cofactor"], 195);
    let f = Fp2::new(211).unwrap();
    let curve = Curve {
        f,
        a: [190, 0],
        b: [37, 0],
    };
    let generator = parse_point(&public["generator"]);
    let mut ops = GroupOps::default();
    assert!(curve.is_on_curve(generator));
    assert!(curve.mul(&mut ops, generator, 229).infinity);
    let mut signed = Vec::new();
    for witness in public["full_base_witness_source"].as_array().unwrap() {
        let raw = parse_point(witness);
        let projected = curve.mul(&mut ops, raw, 195);
        assert!(!projected.infinity && curve.mul(&mut ops, projected, 229).infinity);
        signed.push(projected);
        signed.push(curve.neg(projected));
    }
    assert_eq!(signed.len(), 16);
    let base = chi_base(&curve, &signed);
    take_field_counters();
    let start = Instant::now();
    let mut stats = Stats::default();
    let mut target = curve.identity();
    let mut checked = 0;
    let mut pair_count = 0;
    let mut zero_cases = 0;
    for scalar in 0..229 {
        let expected = direct_pairs(&curve, &signed, target, &mut ops);
        let actual = polynomial_pairs(&curve, &signed, &base, target, &mut ops, &mut stats);
        assert_eq!(
            actual, expected,
            "chi resultant mismatch at scalar {scalar}"
        );
        pair_count += actual.len();
        zero_cases += usize::from(actual.is_empty());
        checked += 1;
        println!(
            "{}",
            json!({"schema":"different_j_chi_resultant_case/v1","scalar_label":scalar,
            "target":target,"pairs":actual.iter().copied().collect::<Vec<_>>(),"status":"pass"})
        );
        target = curve.add(&mut ops, target, generator);
    }
    assert!(target.infinity);
    let torsion = Point::affine(f.from_fp(3), f.one());
    for (label, exceptional) in [
        ("torsion_plus", torsion),
        ("torsion_minus", curve.neg(torsion)),
    ] {
        let expected = direct_pairs(&curve, &signed, exceptional, &mut ops);
        let actual = polynomial_pairs(&curve, &signed, &base, exceptional, &mut ops, &mut stats);
        assert_eq!(actual, expected, "chi resultant mismatch at {label}");
        pair_count += actual.len();
        zero_cases += usize::from(actual.is_empty());
        checked += 1;
        println!(
            "{}",
            json!({"schema":"different_j_chi_resultant_case/v1","scalar_label":label,
            "target":exceptional,"pairs":actual.iter().copied().collect::<Vec<_>>(),"status":"pass"})
        );
    }
    let field = take_field_counters();
    println!(
        "{}",
        json!({"schema":"different_j_chi_resultant_summary/v1",
        "public_sha256":sha256_hex(&bytes),"p":211,"r":229,"factor_base_columns":8,
        "signed_points":signed.len(),"chi_polynomial_degree":base.len()-1,
        "cases":checked,"pair_count":pair_count,"zero_cases":zero_cases,
        "resultants":stats.resultants,"gcds":stats.gcds,"split_attempts":stats.split_attempts,
        "exhaustive_fallback_evaluations":stats.exhaustive_fallback_evaluations,
        "root_candidates":stats.root_candidates,"chi_lift_scans":stats.chi_lift_scans,
        "exact_checks":stats.exact_checks,"false_candidates":stats.false_candidates,
        "vertical_cases":stats.vertical_cases,"tangent_cases":stats.tangent_cases,
        "horizontal_cases":stats.horizontal_cases,"infinity_cases":stats.infinity_cases,
        "exceptional_chi_cases":stats.exceptional_chi_cases,
        "multiple_resultant_roots":stats.multiple_resultant_roots,
        "group_ops":ops,"fp_muls":field.0,"fp_invs":field.1,
        "wall_ns":start.elapsed().as_nanos(),"status":"pass"})
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chord_elimination_vanishes_on_group_sum() {
        let f = Fp2::new(211).unwrap();
        let curve = Curve {
            f,
            a: [190, 0],
            b: [37, 0],
        };
        let generator = Point::affine([170, 0], [0, 109]);
        let mut ops = GroupOps::default();
        let p1 = curve.mul(&mut ops, generator, 17);
        let p2 = curve.mul(&mut ops, generator, 73);
        let target = curve.add(&mut ops, p1, p2);
        let minus_target = curve.neg(target);
        let slope = f.mul(
            f.sub(minus_target.y, p1.y),
            f.inv(f.sub(minus_target.x, p1.x)).unwrap(),
        );
        assert_ne!(slope, f.zero());
        let z = chi(&curve, p2).unwrap();
        let h = chord_resultant(&curve, p1, slope);
        assert_eq!(eval(&f, &h, z), f.zero());
    }

    #[test]
    fn rational_root_split_returns_every_distinct_root() {
        let f = Fp2::new(211).unwrap();
        let roots = [[5, 1], [17, 3], [201, 40]];
        let polynomial = roots.iter().fold(vec![f.one()], |acc, &root| {
            mul(&f, &acc, &vec![f.neg(root), f.one()])
        });
        let mut stats = Stats::default();
        let actual = all_roots(&f, &polynomial, &mut stats);
        let mut expected = roots.to_vec();
        expected.sort_unstable();
        assert_eq!(actual, expected);
    }
}
