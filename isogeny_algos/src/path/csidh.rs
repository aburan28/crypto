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
use crate::bigint::Big;
use crate::field::{is_prime, Field, Rng, Zp};
use crate::fpm::FpM;
use crate::kernel::montgomery::*;
use std::collections::HashMap;

pub const ODD_PRIMES: [u64; 16] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59];

pub struct Csidh<F: Field> {
    pub fp: F,
    pub primes: Vec<u64>,
    /// p + 1
    pub p1: Big,
}

impl Csidh<Zp> {
    /// None unless p = 4 prod(l_i) - 1 is a prime below 2^62.
    pub fn new(primes: &[u64]) -> Option<Self> {
        let prod: u128 = primes.iter().map(|&l| l as u128).product();
        let p = 4 * prod - 1;
        if p >= (1u128 << 62) || !is_prime(p as u64) {
            return None;
        }
        Some(Csidh { fp: Zp::new(p as u64), primes: primes.to_vec(), p1: Big::from_u128(p + 1) })
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
}

impl<const N: usize> Csidh<FpM<N>> {
    /// p = 4 prod(primes) - 1 over N limbs (the caller vouches that p is prime; checked by Fermat).
    pub fn new_big(primes: &[u64]) -> Option<Self> {
        let mut prod = Big::from_u64(4);
        for &l in primes {
            prod = prod.mul_small(l);
        }
        let p = prod.sub_small(1);
        let f = FpM::<N>::new(&p);
        let a = f.from_u64(3);
        if f.pow_big(a, &p.sub_small(1)) != f.one() {
            return None;
        }
        Some(Csidh { fp: f, primes: primes.to_vec(), p1: prod })
    }

    /// CSIDH-512: primes 3, 5, ..., 373 (the first 73 odd primes) and 587.
    pub fn csidh512() -> Csidh<FpM<8>> {
        let mut primes = vec![];
        let mut l = 3u64;
        while primes.len() < 73 {
            if is_prime(l) {
                primes.push(l);
            }
            l += 2;
        }
        primes.push(587);
        Csidh::<FpM<8>>::new_big(&primes).expect("CSIDH-512 prime")
    }
}

impl<F: Field> Csidh<F> {
    fn rhs(&self, a: F::E, x: F::E) -> F::E {
        let f = &self.fp;
        f.add(f.mul(f.add(f.mul(x, x), f.mul(a, x)), x), x)
    }

    fn cofactor(&self, divide_by: &[u64]) -> Big {
        let mut k = self.p1.clone();
        for &l in divide_by {
            k = k.divrem_small(l).0;
        }
        k
    }

    /// One l_i-isogeny from E_A in the given direction; returns the new Montgomery coefficient.
    pub fn step(&self, a: F::E, idx: usize, positive: bool, rng: &mut Rng) -> F::E {
        let f = &self.fp;
        let ell = self.primes[idx];
        let e24 = a24(f, a);
        let cof = self.cofactor(&[ell]);
        loop {
            let x = f.random(rng);
            let r = self.rhs(a, x);
            if f.is_zero(r) || f.sqrt(r).is_some() != positive {
                continue;
            }
            let k = ladder_xz(f, e24, (x, f.one()), &cof);
            if f.is_zero(k.1) {
                continue;
            }
            let d = ((ell - 1) / 2) as usize;
            let ms = multiples(f, e24, k, d);
            return velu_codomain(f, a, &ms, ell);
        }
    }

    /// [prod l_i^{e_i}] E_A, one isogeny at a time (Algorithm 1 of CLMPR, no batching).
    pub fn action(&self, mut a: F::E, exps: &[i32], rng: &mut Rng) -> F::E {
        for (i, &e) in exps.iter().enumerate() {
            for _ in 0..e.unsigned_abs() {
                a = self.step(a, i, e > 0, rng);
            }
        }
        a
    }

    /// [prod l_i^{e_i}] E_A by CLMPR's algorithm: one sampled point serves every prime whose exponent
    /// has the sign of the point (E or its twist): clear the cofactor once, then peel off one prime at
    /// a time, pushing the point through each isogeny.
    pub fn action_batched(&self, mut a: F::E, exps: &[i32], rng: &mut Rng) -> F::E {
        let f = &self.fp;
        let mut e = exps.to_vec();
        let is_square = |v: F::E| f.pow_big(v, &self.p1.sub_small(2).shr(1)) == f.one();
        while e.iter().any(|&x| x != 0) {
            let x = f.random(rng);
            let r = self.rhs(a, x);
            if f.is_zero(r) {
                continue;
            }
            let positive = is_square(r);
            let set: Vec<usize> = (0..e.len()).filter(|&i| e[i] != 0 && (e[i] > 0) == positive).collect();
            if set.is_empty() {
                continue;
            }
            let set_primes: Vec<u64> = set.iter().map(|&i| self.primes[i]).collect();
            let mut q = ladder_xz(f, a24(f, a), (x, f.one()), &self.cofactor(&set_primes));
            let mut remaining = set_primes.clone();
            for &i in set.iter().rev() {
                let ell = self.primes[i];
                remaining.retain(|&l| l != ell);
                if f.is_zero(q.1) {
                    break;
                }
                let e24 = a24(f, a);
                let mut k = Big::from_u64(1);
                for &l in &remaining {
                    k = k.mul_small(l);
                }
                let kp = ladder_xz(f, e24, q, &k);
                if f.is_zero(kp.1) {
                    continue; // q has no l-part: this prime waits for another point
                }
                let d = ((ell - 1) / 2) as usize;
                let ms = multiples(f, e24, kp, d);
                let a_new = velu_codomain(f, a, &ms, ell);
                q = isog_xz(f, &ms, q);
                a = a_new;
                e[i] -= if positive { 1 } else { -1 };
            }
        }
        a
    }

    /// [prod l_i^{e_i}] E_A with projective coefficients (no inversion per isogeny) and a tree
    /// strategy for the kernel points: for the primes S served by one sampled point Q (order
    /// dividing prod S), split S = L u R, compute [prod R] Q for the L-subtree while Q waits on a
    /// stack and is pushed through every L-isogeny, then recurse on R with the pushed Q. Ladder
    /// work drops from O(|S|^2) prime-sized multiplications (CLMPR's loop in `action_batched`)
    /// to O(|S| log |S|) for a balanced tree; the shape is chosen per round by `plan`, which
    /// trades ladder length against pushes (in the spirit of the SIDH strategies of De Feo–Jao–
    /// Plût and their CSIDH adaptations, e.g. Hutchinson–LeGrow–Koziel–Azarderakhsh 2020).
    pub fn action_fast(&self, a: F::E, exps: &[i32], rng: &mut Rng) -> F::E {
        let f = &self.fp;
        let mut e = exps.to_vec();
        let mut k = proj24(f, a);
        let half = self.p1.sub_small(2).shr(1);
        let mut stack = Vec::new();
        while e.iter().any(|&x| x != 0) {
            let a_aff = affine_a(f, k);
            let x = f.random(rng);
            let r = self.rhs(a_aff, x);
            if f.is_zero(r) {
                continue;
            }
            let positive = f.pow_big(r, &half) == f.one();
            let set: Vec<usize> = (0..e.len()).filter(|&i| e[i] != 0 && (e[i] > 0) == positive).collect();
            if set.is_empty() {
                continue;
            }
            let set_primes: Vec<u64> = set.iter().map(|&i| self.primes[i]).collect();
            let q = ladder_p(f, k, (x, f.one()), &self.cofactor(&set_primes));
            let sign = if positive { 1 } else { -1 };
            let split = self.plan(&set);
            self.tree(&mut k, q, &set, 0, set.len() - 1, &split, &mut e, sign, &mut stack);
        }
        affine_a(f, k)
    }

    /// Cheapest tree shape for the primes `set` (ascending), as split[i][j] = last index of the
    /// left subtree of the node covering set[i..=j]. A leaf at stack depth s costs its isogeny
    /// plus s pushes, and a left edge adds one stacked point to every leaf below it, so
    /// C(i,j) = min_k ladder(set[k+1..=j]) + C(i,k) + sum_{t<=k} push(t) + C(k+1,j):
    /// O(|set|^3) with costs in field multiplications (xADD = xDBL = 6, push = 4d + 4,
    /// isogeny = 8d + 3 log l + 10 for l = 2d + 1).
    fn plan(&self, set: &[usize]) -> Vec<Vec<usize>> {
        let n = set.len();
        let lg: Vec<f64> = set.iter().map(|&i| (self.primes[i] as f64).log2()).collect();
        let d: Vec<f64> = set.iter().map(|&i| ((self.primes[i] - 1) / 2) as f64).collect();
        let mut pre_lg = vec![0.0; n + 1];
        let mut pre_push = vec![0.0; n + 1];
        for t in 0..n {
            pre_lg[t + 1] = pre_lg[t] + lg[t];
            pre_push[t + 1] = pre_push[t] + 4.0 * d[t] + 4.0;
        }
        let mut c = vec![vec![0.0f64; n]; n];
        let mut split = vec![vec![0usize; n]; n];
        for t in 0..n {
            c[t][t] = 8.0 * d[t] + 3.0 * lg[t] + 10.0;
        }
        for len in 2..=n {
            for i in 0..=n - len {
                let j = i + len - 1;
                let mut best = f64::INFINITY;
                for k in i..j {
                    let v = 12.0 * (pre_lg[j + 1] - pre_lg[k + 1])
                        + c[i][k]
                        + (pre_push[k + 1] - pre_push[i])
                        + c[k + 1][j];
                    if v < best {
                        best = v;
                        split[i][j] = k;
                    }
                }
                c[i][j] = best;
            }
        }
        split
    }

    #[allow(clippy::too_many_arguments)]
    fn tree(
        &self,
        k: &mut Proj24<F::E>,
        q: XZ<F::E>,
        set: &[usize],
        lo: usize,
        hi: usize,
        split: &[Vec<usize>],
        e: &mut [i32],
        sign: i32,
        stack: &mut Vec<XZ<F::E>>,
    ) {
        let f = &self.fp;
        if f.is_zero(q.1) {
            return; // no component of order in `set`: these primes wait for another point
        }
        if lo == hi {
            let i = set[lo];
            let ell = self.primes[i];
            let ms = multiples_p(f, *k, q, ((ell - 1) / 2) as usize);
            let pre = kernel_pre(f, &ms);
            *k = velu_codomain_p(f, *k, &pre, ell);
            for p in stack.iter_mut() {
                *p = isog_xz_pre(f, &pre, *p);
            }
            e[i] -= sign;
            return;
        }
        // left subtree first while Q waits on the stack; Q, pushed through every left isogeny,
        // then has order dividing the right subtree's primes
        let mid = split[lo][hi];
        let mut m = Big::from_u64(1);
        for &i in &set[mid + 1..=hi] {
            m = m.mul_small(self.primes[i]);
        }
        let ql = ladder_p(f, *k, q, &m);
        stack.push(q);
        self.tree(k, ql, set, lo, mid, split, e, sign, stack);
        let q = stack.pop().unwrap();
        self.tree(k, q, set, mid + 1, hi, split, e, sign, stack);
    }

    /// Order of the ideal class of l_i in Cl(Z[sqrt(-p)]): length of the cycle E_0 -> ... -> E_0.
    pub fn ideal_order(&self, idx: usize, rng: &mut Rng, cap: usize) -> Option<usize> {
        let mut a = self.fp.zero();
        for n in 1..=cap {
            a = self.step(a, idx, true, rng);
            if self.fp.is_zero(a) {
                return Some(n);
            }
        }
        None
    }

    /// Meet-in-the-middle for E_target = [e] E_start with e in [-m, m]^n: boxes [0,m]^n from both
    /// ends, collision [a]E_start = [b]E_target gives e = a - b. Returns (e, nodes visited).
    pub fn mitm(&self, start: F::E, target: F::E, m: usize, rng: &mut Rng) -> Option<(Vec<i32>, usize)> {
        let n = self.primes.len();
        let base = m + 1;
        let total = base.pow(n as u32);
        let build = |origin: F::E, rng: &mut Rng| -> Vec<F::E> {
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
        let mut index: HashMap<F::E, usize> = HashMap::with_capacity(total);
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
