//! Strassen's modular factorial and deterministic factorisation (Bernstein–De Feo–Leroux–Smith
//! 2020, "Faster computation of isogenies of large prime degree", section 2).
//!
//! * §2.3: l! mod n from b = floor(sqrt l), f(X) = prod_{i=1}^{b} (X + i) (product tree), the
//!   values f(0), f(b), ..., f((b-1) b) by multipoint evaluation (remainder tree) — their product
//!   is (b^2)! — and the tail (b^2 + 1) ... l naively. O~(sqrt l) ring operations in Z/nZ.
//! * §2.1: gcd(n, l! mod n) > 1 iff n has a prime factor <= l, so a binary search over
//!   l <= sqrt n finds the smallest prime factor; repeating factors n completely, O~(n^{1/4}).
//!
//! n may be composite, so the arithmetic is the ring Z/nZ ([`ZnRing`]). It implements
//! [`Field`] only so the generic polynomial code runs over it: every polynomial that is divided
//! by here is monic (products of X - r), so the only inversions ever requested are of 1 (the
//! leading coefficient, and the constant term of a reversed monic polynomial in the Newton
//! reciprocal). [`ZnRing::inv`] panics on a non-unit rather than return garbage.
//!
//! The remainder tree is [`RemTree`]: every node's reversed polynomial has its power-series
//! reciprocal precomputed (Newton), and each remainder is two Karatsuba (low) products, so a
//! tree over m points costs O(M(m) log m) with M the multiplication cost. (`poly::rem` divides
//! schoolbook whenever the dividend is longer than 2 deg - 1, which is the normal case in a
//! remainder tree, and `poly::EvalTree`'s middle products are schoolbook dot products; both make
//! a tree over m points cost Θ(m^2) — fine for √élu's sizes, but it would hide the
//! square-root behaviour.) With Karatsuba M(m) = m^1.58 the cost of the factorial is
//! O(l^0.79 log l), not the paper's O~(l^{1/2}) (which assumes quasi-linear multiplication). At
//! opt-level 1 on 61-bit n it was slower than the l - 1 naive multiplications at l = 10^4 and
//! faster from l = 10^5 on (2.5x at 10^6; `strassen_factorial_timing` prints the table).
use crate::field::{lazy_conv, lazy_dot_rev, Field, Red128, Rng};
use crate::poly::{self, Poly};
use crate::series;

/// The ring Z/nZ for 2 <= n < 2^62, any n (prime, composite, even). Elements are residues in
/// [0, n). Not a field when n is composite: see the module documentation for why the generic
/// [`Field`] interface is still sound for the monic-tree operations used here.
#[derive(Clone, Copy, Debug)]
pub struct ZnRing {
    pub n: u64,
    /// division-free reduction for the lazy kernels; Montgomery needs n odd
    red: Option<Red128>,
}

impl ZnRing {
    pub fn new(n: u64) -> Self {
        assert!((2..1u64 << 62).contains(&n), "ZnRing: need 2 <= n < 2^62");
        ZnRing {
            n,
            red: (n % 2 == 1).then(|| Red128::new(n)),
        }
    }
    /// Inverse of a unit, None if gcd(a, n) > 1.
    pub fn try_inv(&self, a: u64) -> Option<u64> {
        let (mut r0, mut r1) = (self.n as i128, (a % self.n) as i128);
        let (mut s0, mut s1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        if r0 != 1 {
            return None;
        }
        Some(s0.rem_euclid(self.n as i128) as u64)
    }
}

impl Field for ZnRing {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.n {
            s - self.n
        } else {
            s
        }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.n - b
        }
    }
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.n - a
        }
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.n as u128) as u64
    }
    /// Inverse of a unit; panics on a non-unit (never reached by the monic-tree code).
    fn inv(&self, a: u64) -> u64 {
        self.try_inv(a)
            .expect("ZnRing: inverse of a non-unit (the monic-tree code never asks for one)")
    }
    fn from_u64(&self, v: u64) -> u64 {
        v % self.n
    }
    /// The modulus n (a characteristic only when n is prime).
    fn char(&self) -> u64 {
        self.n
    }
    fn size(&self) -> u128 {
        self.n as u128
    }
    fn random(&self, rng: &mut Rng) -> u64 {
        rng.below(self.n)
    }
    /// Square roots are not defined for this use; always None.
    fn sqrt(&self, _a: u64) -> Option<u64> {
        None
    }
    // lazy accumulation is plain modular arithmetic: valid for any odd modulus below 2^62
    // (the kernels reduce by Montgomery steps); even n accumulates in u128 directly
    fn conv(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        self.conv_trunc(a, b, (a.len() + b.len()).saturating_sub(1))
    }
    fn conv_trunc(&self, a: &[u64], b: &[u64], k: usize) -> Vec<u64> {
        match &self.red {
            Some(red) => lazy_conv(a, b, red, k, |x| red.reduce(x)),
            None if a.is_empty() || b.is_empty() => vec![0; k],
            None => (0..k)
                .map(|t| {
                    let lo = (t + 1).saturating_sub(b.len());
                    let hi = t.min(a.len().saturating_sub(1));
                    (lo..=hi).fold(0u64, |acc, i| self.add(acc, self.mul(a[i], b[t - i])))
                })
                .collect(),
        }
    }
    fn dot_rev(&self, a: &[u64], b: &[u64]) -> u64 {
        match &self.red {
            Some(red) => lazy_dot_rev(a, b, red, |x| red.reduce(x)),
            None => {
                let (n, bl) = (a.len().min(b.len()), b.len());
                (0..n).fold(0, |acc, i| self.add(acc, self.mul(a[i], b[bl - 1 - i])))
            }
        }
    }
}

/// floor(sqrt(v)).
pub fn isqrt(v: u64) -> u64 {
    if v < 2 {
        return v;
    }
    let mut x = (v as f64).sqrt() as u64;
    while (x as u128) * (x as u128) > v as u128 {
        x -= 1;
    }
    while ((x + 1) as u128) * ((x + 1) as u128) <= v as u128 {
        x += 1;
    }
    x
}

pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Low `k` coefficients of a * b (a "low product"; Karatsuba above the series threshold).
fn low_mul<F: Field>(f: &F, a: &[F::E], b: &[F::E], k: usize) -> Vec<F::E> {
    let (la, lb) = (a.len().min(k), b.len().min(k));
    if la == 0 || lb == 0 {
        return vec![f.zero(); k];
    }
    if la.min(lb) >= f.karatsuba_threshold() {
        let mut r = poly::mul_raw(f, &a[..la], &b[..lb]);
        r.resize(k.max(r.len()), f.zero());
        r.truncate(k);
        return r;
    }
    f.conv_trunc(&a[..la], &b[..lb], k)
}

/// a mod m for monic m of degree n, given `rinv` = 1 / rev(m) as a power series with at least
/// a.len() - n coefficients: quotient from one low product with rinv, remainder from one low
/// product with m (A.5: "precomputing reciprocals of tree nodes", "low products").
pub fn rem_pre<F: Field>(f: &F, a: &[F::E], m: &[F::E], rinv: &[F::E]) -> Poly<F> {
    let n = m.len() - 1;
    if a.len() <= n {
        let mut r = a.to_vec();
        poly::trim(f, &mut r);
        return r;
    }
    let k = a.len() - n;
    assert!(rinv.len() >= k, "rem_pre: reciprocal too short");
    let rev_a: Vec<F::E> = a.iter().rev().take(k).copied().collect();
    let qr = low_mul(f, &rev_a, &rinv[..k], k);
    let q: Vec<F::E> = qr.into_iter().rev().collect();
    let qm = low_mul(f, &q, m, n);
    let mut r: Vec<F::E> = (0..n).map(|i| f.sub(a[i], qm[i])).collect();
    poly::trim(f, &mut r);
    r
}

/// Product tree over monic leaves with the reciprocal of every reversed node precomputed;
/// remainders of one polynomial modulo every leaf (the remainder tree of §2.3 and §3.2).
pub struct RemTree<F: Field> {
    /// level 0 = the leaves, last level = the root (a single node)
    pub levels: Vec<Vec<Poly<F>>>,
    /// rinv[l][i] = 1 / rev(levels[l][i]) to deg(parent) - deg(node) + 1 coefficients (enough
    /// for the parent's remainder); the root's is computed per call to the needed length
    rinv: Vec<Vec<Vec<F::E>>>,
}

impl<F: Field> RemTree<F> {
    /// Tree over monic, non-constant leaves (not empty).
    pub fn new(f: &F, leaves: Vec<Poly<F>>) -> Self {
        assert!(!leaves.is_empty(), "RemTree: no leaves");
        let mut levels = vec![leaves];
        while levels.last().unwrap().len() > 1 {
            let prev = levels.last().unwrap();
            let next: Vec<Poly<F>> = prev
                .chunks(2)
                .map(|c| {
                    if c.len() == 2 {
                        poly::mul_raw(f, &c[0], &c[1])
                    } else {
                        c[0].clone()
                    }
                })
                .collect();
            levels.push(next);
        }
        let top = levels.len() - 1;
        let mut rinv: Vec<Vec<Vec<F::E>>> = Vec::with_capacity(top);
        for l in 0..top {
            let mut row = Vec::with_capacity(levels[l].len());
            for (i, node) in levels[l].iter().enumerate() {
                let parent = &levels[l + 1][i / 2];
                let k = (parent.len() - node.len()).max(1) + 1;
                row.push(rev_inv(f, node, k));
            }
            rinv.push(row);
        }
        RemTree { levels, rinv }
    }

    /// Tree over the linear leaves X - r.
    pub fn from_roots(f: &F, roots: &[F::E]) -> Self {
        Self::new(f, roots.iter().map(|&r| vec![f.neg(r), f.one()]).collect())
    }

    pub fn root(&self) -> &Poly<F> {
        &self.levels.last().unwrap()[0]
    }

    pub fn len(&self) -> usize {
        self.levels[0].len()
    }

    pub fn is_empty(&self) -> bool {
        self.levels[0].is_empty()
    }

    /// a mod leaf, for every leaf, in leaf order.
    pub fn rems(&self, f: &F, a: &Poly<F>) -> Vec<Poly<F>> {
        let top = self.levels.len() - 1;
        let root = self.root();
        let r0 = if a.len() >= root.len() {
            let k = a.len() - (root.len() - 1);
            rem_pre(f, a, root, &rev_inv(f, root, k))
        } else {
            let mut r = a.clone();
            poly::trim(f, &mut r);
            r
        };
        let mut cur = vec![r0];
        for l in (0..top).rev() {
            let nodes = &self.levels[l];
            let mut next = Vec::with_capacity(nodes.len());
            for (i, r) in cur.iter().enumerate() {
                for c in [2 * i, 2 * i + 1] {
                    if c < nodes.len() {
                        next.push(rem_pre(f, r, &nodes[c], &self.rinv[l][c]));
                    }
                }
            }
            cur = next;
        }
        cur
    }

    /// a(r_i) for linear leaves X - r_i.
    pub fn eval(&self, f: &F, a: &Poly<F>) -> Vec<F::E> {
        self.rems(f, a)
            .into_iter()
            .map(|r| if r.is_empty() { f.zero() } else { r[0] })
            .collect()
    }
}

/// Root of a product tree over the given polynomials (Karatsuba at every level).
pub fn product_tree_root<F: Field>(f: &F, mut layer: Vec<Poly<F>>) -> Poly<F> {
    while layer.len() > 1 {
        layer = layer
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    poly::mul_raw(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    layer.pop().unwrap_or_else(|| vec![f.one()])
}

/// 1 / rev(m) to k coefficients (m monic: rev(m) has constant term 1).
fn rev_inv<F: Field>(f: &F, m: &[F::E], k: usize) -> Vec<F::E> {
    let rev: Vec<F::E> = m.iter().rev().copied().collect();
    series::inv(f, &rev, k)
}

/// l! mod n by l - 1 multiplications (reference).
pub fn factorial_naive(l: u64, n: u64) -> u64 {
    let nn = n as u128;
    let mut r = 1u128 % nn;
    for i in 2..=l {
        r = r * (i as u128 % nn) % nn;
    }
    r as u64
}

/// Below this b = floor(sqrt l) the naive product is used: measured crossover between l = 10^4
/// and 10^5 (b = 100..316) with Karatsuba trees over 61-bit n.
pub const STRASSEN_MIN_B: u64 = 200;

/// l! mod n by Strassen's algorithm (§2.3), 1 <= n < 2^62.
pub fn factorial_mod(l: u64, n: u64) -> u64 {
    factorial_mod_with(l, n, STRASSEN_MIN_B)
}

/// As [`factorial_mod`], with the naive cut-off as a parameter (tests force the trees).
pub fn factorial_mod_with(l: u64, n: u64, min_b: u64) -> u64 {
    if n == 1 {
        return 0;
    }
    let b = isqrt(l);
    if b < min_b.max(2) {
        return factorial_naive(l, n);
    }
    let r = ZnRing::new(n);
    // f(X) = prod_{i=1}^{b} (X + i) and its values at X = k b, k = 0..b-1:
    // f(k b) = (k b + 1) ... (k b + b), so the product of all values is (b^2)!
    let f_poly = product_tree_root(&r, (1..=b).map(|i| vec![r.from_u64(i), r.one()]).collect());
    let pts: Vec<u64> = (0..b).map(|k| r.from_u64(k * b)).collect();
    let vals = RemTree::from_roots(&r, &pts).eval(&r, &f_poly);
    let mut acc = vals.into_iter().fold(r.one(), |a, v| r.mul(a, v));
    for i in b * b + 1..=l {
        acc = r.mul(acc, r.from_u64(i));
    }
    acc
}

/// Smallest prime factor of n >= 2 by §2.1: binary search for the least l with
/// gcd(n, l! mod n) > 1, over lo < l <= floor(sqrt n) (`lo`: no prime factor <= lo is known).
pub fn smallest_factor_from(n: u64, lo: u64) -> u64 {
    assert!(n >= 2);
    let s = isqrt(n);
    if s <= lo || gcd(n, factorial_mod(s, n)) == 1 {
        return n; // no factor <= sqrt n: n is prime
    }
    // invariant: gcd(n, lo!) = 1, gcd(n, hi!) > 1
    let (mut lo, mut hi) = (lo.max(1), s);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if gcd(n, factorial_mod(mid, n)) > 1 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

pub fn smallest_factor(n: u64) -> u64 {
    smallest_factor_from(n, 1)
}

/// Complete factorisation (prime, exponent) in increasing order by repeated smallest-factor
/// searches (§2.1), deterministic, O~(n^{1/4}) ring operations per search.
pub fn factor(n: u64) -> Vec<(u64, u32)> {
    let mut out = vec![];
    let mut m = n;
    let mut lo = 1u64;
    while m > 1 {
        let p = smallest_factor_from(m, lo);
        let mut e = 0;
        while m.is_multiple_of(p) {
            m /= p;
            e += 1;
        }
        out.push((p, e));
        lo = p; // every remaining prime factor exceeds p
    }
    out
}
