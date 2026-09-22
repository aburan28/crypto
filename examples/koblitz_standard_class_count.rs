//! Exact cofactor-class admissibility count for a standard algebraic factor base.
//! Usage: koblitz_standard_class_count <n> <a> <ell>
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, point_key, KoblitzCurve,
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::json;
use std::collections::HashMap;
use std::time::Instant;

fn factors(mut n: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            out.push(p);
            while n % p == 0 { n /= p; }
        }
        p += if p == 2 { 1 } else { 2 };
    }
    if n > 1 { out.push(n); }
    out
}

fn admissible_counts(frequency: &[u64]) -> (u128, u128, u128, u128) {
    let h = frequency.len();
    let mut ordered = 0u128;
    for (x, &fx) in frequency.iter().enumerate() {
        if fx == 0 { continue; }
        for (y, &fy) in frequency.iter().enumerate() {
            if fy == 0 { continue; }
            let z = (h - (x + y) % h) % h;
            ordered += u128::from(fx) * u128::from(fy) * u128::from(frequency[z]);
        }
    }
    let transposition: u128 = frequency
        .iter()
        .enumerate()
        .map(|(x, &fx)| {
            let y = (h - (2 * x) % h) % h;
            u128::from(fx) * u128::from(frequency[y])
        })
        .sum();
    let three_cycle: u128 = frequency
        .iter()
        .enumerate()
        .filter(|(x, _)| (3 * *x) % h == 0)
        .map(|(_, &fx)| u128::from(fx))
        .sum();
    let unordered = (ordered + 3 * transposition + 2 * three_cycle) / 6;
    (ordered, transposition, three_cycle, unordered)
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 4, "usage: <n> <a> <ell>");
    let n: u32 = args[1].parse().unwrap();
    let a: u8 = args[2].parse().unwrap();
    let ell: u32 = args[3].parse().unwrap();
    let begin = Instant::now();
    let kc = KoblitzCurve::new(a, n).expect("curve");
    let h = kc.cofactor.to_u64().expect("cofactor fits u64");
    assert!(h <= 100_000, "bounded public torsion quotient only");
    let fb = build_standard_subspace_factor_base(&kc, ell).expect("standard base");
    let classes = fb.cofactor_classes(&kc);
    let primes = factors(h);
    let generator = classes
        .iter()
        .find(|p| {
            **p != BinaryPoint::Infinity
                && primes.iter().all(|q| {
                    kc.mul(p, &BigUint::from(h / q)) != BinaryPoint::Infinity
                })
        })
        .expect("factor-base classes contain a generator of the cofactor group")
        .clone();
    let mut labels = HashMap::new();
    let mut current = BinaryPoint::Infinity;
    for k in 0..h {
        assert!(labels.insert(point_key(&current), k).is_none(), "short generator orbit");
        current = kc.add(&current, &generator);
    }
    assert_eq!(current, BinaryPoint::Infinity, "cofactor generator did not close");
    assert_eq!(labels.len(), h as usize);
    let mut frequency = vec![0u64; h as usize];
    for class in &classes {
        let label = *labels.get(&point_key(class)).expect("class outside cyclic quotient");
        frequency[label as usize] += 1;
    }
    let support = frequency.iter().filter(|&&x| x != 0).count();
    let (ordered, transposition, three_cycle, unordered) = admissible_counts(&frequency);
    let all = (fb.points.len() as u128)
        * (fb.points.len() as u128 + 1)
        * (fb.points.len() as u128 + 2)
        / 6;
    let generator_key = point_key(&generator);
    println!("{}", json!({
        "schema":"koblitz_standard_class_count.v1",
        "n":n,"a":a,"ell":ell,"subgroup_order":kc.subgroup_order.to_string(),
        "cofactor":h,"factor_points":fb.points.len(),"projected_columns":
            crypto_lib::cryptanalysis::koblitz_index_calculus::projected_signed_orbit_count(&kc,&fb),
        "cofactor_group_cyclic":true,"cofactor_generator":[generator_key.0.to_string(),generator_key.1.to_string()],
        "class_support":support,"class_frequency":frequency,
        "all_unordered_triples":all.to_string(),
        "cofactor_admissible_ordered_triples":ordered.to_string(),
        "cofactor_admissible_unordered_triples":unordered.to_string(),
        "transposition_fixed":transposition.to_string(),"three_cycle_fixed":three_cycle.to_string(),
        "admissible_fraction_of_triples":unordered as f64 / all as f64,
        "uniform_subgroup_coverage_ceiling":unordered as f64 / kc.subgroup_order.to_f64().unwrap(),
        "elapsed_seconds":begin.elapsed().as_secs_f64(),
        "scope":"public cofactor-torsion classes only; no prime-subgroup discrete logarithms or target enumeration"
    }));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn burnside_count_matches_explicit_unordered_triples() {
        let frequency = [2u64, 1, 3, 0, 2];
        let labels: Vec<usize> = frequency
            .iter()
            .enumerate()
            .flat_map(|(label, &count)| std::iter::repeat_n(label, count as usize))
            .collect();
        let mut brute = 0u128;
        for i in 0..labels.len() {
            for j in i..labels.len() {
                for k in j..labels.len() {
                    brute += u128::from((labels[i] + labels[j] + labels[k]) % frequency.len() == 0);
                }
            }
        }
        assert_eq!(admissible_counts(&frequency).3, brute);
    }
}
