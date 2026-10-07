//! CSIDH-style class-group action (Castryck–Lange–Martindale–Panny–Renes 2018), the supersingular
//! F_p instantiation of the Couveignes–Rostovtsev–Stolbunov hard homogeneous space, on Montgomery
//! curves with x-only Vélu. p = 4 l_1 ... l_n - 1; curves E_A: y^2 = x^3 + A x^2 + x with A reachable
//! from A = 0. The ideal class of l_i acts by the l_i-isogeny whose kernel lies in E_A(F_p) (positive
//! direction) or in the twist (negative): the direction is read off by whether x^3+Ax^2+x is a
//! square, so no eigenvalue computation is needed (compare `couveignes.rs` for ordinary curves).
//!
//! Also: the meet-in-the-middle solver for the group-action inversion problem on exponent boxes,
//! the order of an ideal class from its isogeny cycle (Couveignes 2006, Stolbunov 2010), and an
//! independent class-number count to validate cycle lengths.
use crate::field::{is_prime, Field, Rng, Zp};
use crate::kernel::montgomery::*;
use std::collections::HashMap;

pub const ODD_PRIMES: [u64; 16] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59];

pub struct Csidh {
    pub fp: Zp,
    pub primes: Vec<u64>,
}

impl Csidh {
    /// None unless p = 4 prod(l_i) - 1 is a prime below 2^62.
    pub fn new(primes: &[u64]) -> Option<Self> {
        let prod: u128 = primes.iter().map(|&l| l as u128).product();
        let p = 4 * prod - 1;
        if p >= (1u128 << 62) || !is_prime(p as u64) {
            return None;
        }
        Some(Csidh {
            fp: Zp::new(p as u64),
            primes: primes.to_vec(),
        })
    }

    /// n primes: the first n-1 odd primes plus the smallest further odd prime making p prime.
    pub fn with_n_primes(n: usize) -> Option<Self> {
        let mut base: Vec<u64> = ODD_PRIMES[..n - 1].to_vec();
        for &last in &ODD_PRIMES[n - 1..] {
            base.push(last);
            if let Some(c) = Csidh::new(&base) {
                return Some(c);
            }
            base.pop();
        }
        None
    }

    pub fn p(&self) -> u64 {
        self.fp.p
    }

    fn rhs(&self, a: u64, x: u64) -> u64 {
        let fp = &self.fp;
        fp.add(fp.mul(fp.add(fp.mul(x, x), fp.mul(a, x)), x), x)
    }

    /// One l_i-isogeny from E_A in the given direction; returns the new Montgomery coefficient.
    pub fn step(&self, a: u64, idx: usize, positive: bool, rng: &mut Rng) -> u64 {
        let fp = &self.fp;
        let ell = self.primes[idx];
        let e24 = a24(fp, a);
        loop {
            let x = fp.random(rng);
            let r = self.rhs(a, x);
            if r == 0 || fp.sqrt(r).is_some() != positive {
                continue;
            }
            let k = ladder(fp, e24, x, ((fp.p + 1) / ell) as u128);
            if k.1 == 0 {
                continue;
            }
            let d = ((ell - 1) / 2) as usize;
            let ms = multiples(fp, e24, k, d);
            return velu_codomain(fp, a, &ms, ell);
        }
    }

    /// [prod l_i^{e_i}] E_A, one isogeny at a time (Algorithm 1 of CLMPR, no batching).
    pub fn action(&self, mut a: u64, exps: &[i32], rng: &mut Rng) -> u64 {
        for (i, &e) in exps.iter().enumerate() {
            for _ in 0..e.unsigned_abs() {
                a = self.step(a, i, e > 0, rng);
            }
        }
        a
    }

    /// Order of the ideal class of l_i in Cl(Z[sqrt(-p)]): length of the cycle E_0 -> ... -> E_0.
    pub fn ideal_order(&self, idx: usize, rng: &mut Rng, cap: usize) -> Option<usize> {
        let mut a = 0u64;
        for n in 1..=cap {
            a = self.step(a, idx, true, rng);
            if a == 0 {
                return Some(n);
            }
        }
        None
    }

    /// Meet-in-the-middle for E_target = [e] E_start with e in [-m, m]^n: boxes [0,m]^n from both
    /// ends, collision [a]E_start = [b]E_target gives e = a - b. Returns (e, nodes visited).
    pub fn mitm(
        &self,
        start: u64,
        target: u64,
        m: usize,
        rng: &mut Rng,
    ) -> Option<(Vec<i32>, usize)> {
        let n = self.primes.len();
        let base = m + 1;
        let total = base.pow(n as u32);
        let build = |origin: u64, rng: &mut Rng| -> Vec<u64> {
            let mut nodes = vec![origin; total];
            for k in 1..total {
                let (mut digit, mut i) = (k, 0);
                while digit % base == 0 {
                    digit /= base;
                    i += 1;
                }
                nodes[k] = self.step(nodes[k - base.pow(i as u32)], i, true, rng);
            }
            nodes
        };
        let from_start = build(start, rng);
        let mut index: HashMap<u64, usize> = HashMap::with_capacity(total);
        for (i, &a) in from_start.iter().enumerate() {
            index.entry(a).or_insert(i);
        }
        let from_target = build(target, rng);
        let digits = |mut k: usize| -> Vec<i32> {
            (0..n)
                .map(|_| {
                    let d = (k % base) as i32;
                    k /= base;
                    d
                })
                .collect()
        };
        for (j, a) in from_target.iter().enumerate() {
            if let Some(&i) = index.get(a) {
                let (da, db) = (digits(i), digits(j));
                let e: Vec<i32> = da.iter().zip(&db).map(|(x, y)| x - y).collect();
                return Some((e, 2 * total));
            }
        }
        None
    }
}

/// h(D) for a negative discriminant D: number of reduced primitive forms ax^2+bxy+cy^2.
/// O(|D|); for validating cycle lengths on toy parameters only.
pub fn class_number(d: i64) -> u64 {
    assert!(d < 0 && (d.rem_euclid(4) == 0 || d.rem_euclid(4) == 1));
    let dd = -d;
    let gcd = |mut a: i64, mut b: i64| {
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a.abs()
    };
    let mut h = 0u64;
    let amax = ((dd as f64) / 3.0).sqrt() as i64 + 1;
    for a in 1..=amax {
        let mut b = -a + 1;
        while b <= a {
            if (b - d).rem_euclid(2) == 0 {
                let num = b * b - d; // 4ac
                if num % (4 * a) == 0 {
                    let c = num / (4 * a);
                    if c >= a && !(a == c && b < 0) && gcd(gcd(a, b), c) == 1 {
                        h += 1;
                    }
                }
            }
            b += 1;
        }
        // b = -a only when a == c is excluded above unless b == a; handled by range -a+1..=a
    }
    h
}
