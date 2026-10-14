//! **Candidate C4 of `research/pdp_root_finder_requirements_20261012` §8:
//! does the decomposition condition on PKM's Kummer tower have additive
//! structure in the exponents?**
//!
//! On the tower `V = μ_N`, `N = 2^t`, every base abscissa is `ζ^e`.  A
//! representation-technique (Wagner / Howgrave-Graham–Joux) search for
//! `R = P_a + P_b + P_c` needs a homomorphism to a small quotient: some
//! linear form of the exponents `(e_a, e_b, e_c)` whose value is predictable
//! from `R`.  The simplest test: is `e_a + e_b + e_c mod 2^k` concentrated,
//! for any `k`, among the decompositions of one target?  If it is not, no
//! sum-of-exponents split exists and C4 is falsified on this base.
//!
//! The detector is checked on a positive control in which the "group law"
//! is multiplication in `F_p^*` (`x_a·x_b·x_c = x_R`), where the exponent
//! sum is fixed by the target: there the concentration must be 1.
//!
//! For each target, every decomposition over the signed tower base is
//! enumerated exactly (pair table of `P_a + P_b`, then `R − P_c` looked up).
//! Reported per `k`: the mean over targets of the largest residue-class
//! share of `e_a + e_b + e_c mod 2^k` (1/2^k for no structure, 1 for full
//! structure), with the same statistic for random triples of the same size
//! as the null.  Deterministic in the seed.
//!
//! ```text
//! cargo run --release --example kummer_exponent_bias -- --p 946177 --a A --b B --t 8 --targets 2000 --seed 20261013
//! ```

use std::collections::HashMap;
use std::env;

use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use crypto_lib::cryptanalysis::ic_framework::prime_tower::primitive_two_power_root;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

/// Largest share of one residue class of `v mod 2^k` among `vals`.
fn top_share(vals: &[u64], k: u32) -> f64 {
    let mut c: HashMap<u64, usize> = HashMap::new();
    for &v in vals {
        *c.entry(v & ((1 << k) - 1)).or_default() += 1;
    }
    *c.values().max().unwrap_or(&0) as f64 / vals.len().max(1) as f64
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .map(|i| args[i + 1].parse::<u64>().expect(k))
    };
    let p = get("--p").expect("--p");
    let a = get("--a").expect("--a");
    let b = get("--b").expect("--b");
    let t = get("--t").expect("--t") as u32;
    let targets = get("--targets").unwrap_or(2000) as usize;
    let seed = get("--seed").unwrap_or(20261013);
    let curve = PrimeCurve { p, a, b };
    let n = 1u64 << t;
    let zeta = primitive_two_power_root(p, t).expect("2^t | p - 1");
    // Signed tower points with their exponents.
    let mut pts: Vec<(PrimePoint, u64)> = Vec::new();
    let mut z = 1u64;
    for e in 0..n {
        if let Some(y) = curve.sqrt(curve.rhs(z)) {
            let q = PrimePoint::affine(z, y);
            pts.push((q, e));
            if y != 0 {
                pts.push((curve.neg(q), e));
            }
        }
        z = mulmod(z, zeta, p);
    }
    let mut ops = GroupOps::default();
    // Pair table: key of P_a + P_b → exponent sums e_a + e_b (all pairs a ≤ b).
    let mut pairs: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    for i in 0..pts.len() {
        for j in i..pts.len() {
            let s = curve.add(&mut ops, pts[i].0, pts[j].0);
            pairs.entry(curve.key(&s)).or_default().push((i, j));
        }
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let kmax = t.min(6);
    let (mut sum_curve, mut sum_null, mut sum_mult) = (
        vec![0f64; kmax as usize + 1],
        vec![0f64; kmax as usize + 1],
        vec![0f64; kmax as usize + 1],
    );
    let (mut used, mut sols_total) = (0usize, 0usize);
    let mut tries = 0usize;
    while used < targets && tries < targets * 50 {
        tries += 1;
        // A uniform point of the curve as target.
        let x = rng.gen_range(0..p);
        let Some(y) = curve.sqrt(curve.rhs(x)) else {
            continue;
        };
        let r = PrimePoint::affine(x, y);
        // Each decomposition once: the sorted index triple.
        let mut triples: std::collections::BTreeSet<[usize; 3]> = Default::default();
        for (c, &(pc, _)) in pts.iter().enumerate() {
            let rest = curve.add(&mut ops, r, curve.neg(pc));
            if let Some(list) = pairs.get(&curve.key(&rest)) {
                for &(i, j) in list {
                    let mut tr = [i, j, c];
                    tr.sort_unstable();
                    triples.insert(tr);
                }
            }
        }
        let esums: Vec<u64> = triples
            .iter()
            .map(|tr| tr.iter().map(|&i| pts[i].1).sum())
            .collect();
        if esums.len() < 4 {
            continue;
        }
        used += 1;
        sols_total += esums.len();
        // Null: random exponent triples, as many as this target had.
        let null: Vec<u64> = (0..esums.len())
            .map(|_| rng.gen_range(0..n) + rng.gen_range(0..n) + rng.gen_range(0..n))
            .collect();
        // Positive control: x_a·x_b·x_c = x_R in μ_N·g, exponents sum to a
        // fixed residue; take as many random solutions of e_a + e_b + e_c ≡ s.
        let s = rng.gen_range(0..n);
        let mult: Vec<u64> = (0..esums.len())
            .map(|_| {
                let (ea, eb) = (rng.gen_range(0..n), rng.gen_range(0..n));
                let ec = (s + 2 * n - ea - eb) % n;
                ea + eb + ec
            })
            .collect();
        for k in 1..=kmax {
            sum_curve[k as usize] += top_share(&esums, k);
            sum_null[k as usize] += top_share(&null, k);
            sum_mult[k as usize] += top_share(&mult, k);
        }
    }
    println!(
        "{{\"p\":{p},\"a\":{a},\"b\":{b},\"t\":{t},\"base_points\":{},\"targets\":{used},\"decompositions\":{sols_total},\"rows\":[",
        pts.len()
    );
    for k in 1..=kmax {
        let u = used as f64;
        println!(
            "  {{\"k\":{k},\"uniform\":{:.5},\"curve\":{:.5},\"null_same_size\":{:.5},\"positive_control\":{:.5}}}{}",
            1.0 / (1u64 << k) as f64,
            sum_curve[k as usize] / u,
            sum_null[k as usize] / u,
            sum_mult[k as usize] / u,
            if k < kmax { "," } else { "" }
        );
    }
    println!("]}}");
}
