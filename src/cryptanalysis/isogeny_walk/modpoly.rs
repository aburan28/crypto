//! The classical modular polynomial `Φ_ℓ(X, Y) mod p`, computed natively
//! from the `q`-expansion of `j`.
//!
//! `Φ_ℓ(X, Y) = X^{ℓ+1} + Y^{ℓ+1} − X^ℓY^ℓ + Σ c_{ab} X^a Y^b` with
//! `c_{ab} = c_{ba}` and `a, b ≤ ℓ`, and `Φ_ℓ(j(q), j(q^ℓ)) = 0`.  The
//! symmetric unknown `c_{ab}` (`a ≤ b`) first appears at `q^{−(a + ℓb)}`
//! with coefficient one, and those pole orders are distinct, so the
//! system is triangular: walking the exponents from `q^{−ℓ(ℓ+1)}` upward
//! determines each unknown at its own exponent without a division.  Every
//! other exponent, and a band of positive ones, must then vanish; they are
//! checked, so a wrong assumption about the form of `Φ_ℓ` or a wrong `j`
//! series fails the construction instead of yielding a wrong polynomial.
//! Because no division happens, the result is the reduction of the
//! integral `Φ_ℓ` for every `p`.
//!
//! `j = 1728·E4³ / (E4³ − E6²)` is computed mod `p` from the divisor sums
//! in `E4 = 1 + 240 Σ σ₃(n) qⁿ` and `E6 = 1 − 504 Σ σ₅(n) qⁿ`.
//!
//! # Cost
//!
//! Series products use Kronecker substitution into one big-integer
//! product, so each of the `ℓ + 2` powers of `j` costs `M(ℓ²)` word
//! operations instead of `ℓ⁴` field multiplications.  The triangular solve
//! is organised by blocks `b = ℓ + 1, ℓ, …, 0`: inside a block the
//! unknowns `c_{ab}`, `a ≤ b`, are determined from `O(1)` lazily computed
//! terms each, and the block's whole row `Σ_a c_{ab} X^a Y^b` and its
//! mirror `Σ_a c_{ab} X^b Y^a` are then applied to the residual as two
//! products, `O(ℓ³)` per block.  The solve is `O(ℓ⁴)` field
//! multiplications in all, and memory is dominated by the `ℓ + 2` stored
//! powers, `(ℓ + 2) · ℓ² · 32` bytes.

use super::field::{Fe, Field};
use super::poly::{self, Poly};
use num_bigint::BigUint;

/// Extra positive exponents checked beyond those that determine `Φ_ℓ`.
const CHECK_BAND: usize = 8;

/// `Φ_ℓ mod p` as a dense `(ℓ+2) × (ℓ+2)` coefficient table.
#[derive(Clone, Debug)]
pub struct ModPoly {
    pub ell: u64,
    /// `c[a][b]`, the coefficient of `X^a Y^b`.
    pub c: Vec<Vec<Fe>>,
    /// Exponents whose coefficient was required to vanish and did.
    pub checks_passed: usize,
}

/// `J_k` with `j(q) = q^{-1} Σ_{k≥0} J_k q^k`, for `k < len`.
pub fn j_series(f: &Field, len: usize) -> Vec<Fe> {
    let n = len + 1;
    // σ₃ and σ₅ by a divisor sieve, O(n log n) instead of O(n²).
    let mut s3 = vec![0u128; n];
    let mut s5 = vec![0u128; n];
    for d in 1..n {
        let (d3, d5) = ((d as u128).pow(3), (d as u128).pow(5));
        let mut m = d;
        while m < n {
            s3[m] += d3;
            s5[m] += d5;
            m += d;
        }
    }
    let mut e4 = vec![f.zero(); n];
    let mut e6 = vec![f.zero(); n];
    e4[0] = f.one();
    e6[0] = f.one();
    for m in 1..n {
        e4[m] = f.mul(&f.from_u64(240), &f.from_u128(s3[m]));
        e6[m] = f.neg(&f.mul(&f.from_u64(504), &f.from_u128(s5[m])));
    }
    let e4_3 = mul_trunc(f, &mul_trunc(f, &e4, &e4, n), &e4, n);
    let e6_2 = mul_trunc(f, &e6, &e6, n);
    // E4³ − E6² = 1728 q Π(1 − qⁿ)^24 has no constant term.
    let diff: Vec<Fe> = (0..n).map(|i| f.sub(&e4_3[i], &e6_2[i])).collect();
    debug_assert!(f.is_zero(&diff[0]));
    let den: Vec<Fe> = diff[1..].to_vec(); // (E4³ − E6²)/q, constant 1728
    let inv = inv_series(f, &den, len);
    let num = mul_trunc(f, &e4_3[..len], &inv, len);
    let k1728 = f.from_u64(1728);
    (0..len)
        .map(|k| f.mul(&k1728, num.get(k).unwrap_or(&f.zero())))
        .collect()
}

/// Series inverse of `a` (unit constant term) to `len` terms, by Newton
/// iteration `g ← g(2 − ag)` with doubling precision.
fn inv_series(f: &Field, a: &[Fe], len: usize) -> Vec<Fe> {
    let mut g = vec![f.inv(&a[0]).expect("unit constant term")];
    let two = f.from_u64(2);
    while g.len() < len {
        let m = (2 * g.len()).min(len);
        let ag = mul_trunc(f, &a[..a.len().min(m)], &g, m);
        let mut t: Vec<Fe> = ag.iter().map(|x| f.neg(x)).collect();
        t[0] = f.add(&t[0], &two);
        g = mul_trunc(f, &g, &t, m);
    }
    g
}

/// Bytes per Kronecker slot: a product coefficient is below `len · p² <
/// 2^(512 + 64)` for every `len < 2^64`.
const KRONECKER_SLOT: usize = 72;

/// An element `x + y·i` of `F_{p²} = F_p[i]/(i² − d)`, `d` the least
/// non-residue.  Used only for transforms: `p + 1` carries the large power
/// of two for P-256 and P-192, where `p − 1` does not.
#[derive(Clone, Copy)]
struct Fe2(Fe, Fe);

struct Ext<'a> {
    f: &'a Field,
    d: Fe,
}

impl Ext<'_> {
    fn mul(&self, a: &Fe2, b: &Fe2) -> Fe2 {
        let f = self.f;
        let re = f.add(&f.mul(&a.0, &b.0), &f.mul(&self.d, &f.mul(&a.1, &b.1)));
        let im = f.add(&f.mul(&a.0, &b.1), &f.mul(&a.1, &b.0));
        Fe2(re, im)
    }
    fn mul_re(&self, a: &Fe2, c: &Fe) -> Fe2 {
        Fe2(self.f.mul(&a.0, c), self.f.mul(&a.1, c))
    }
    fn add(&self, a: &Fe2, b: &Fe2) -> Fe2 {
        Fe2(self.f.add(&a.0, &b.0), self.f.add(&a.1, &b.1))
    }
    fn sub(&self, a: &Fe2, b: &Fe2) -> Fe2 {
        Fe2(self.f.sub(&a.0, &b.0), self.f.sub(&a.1, &b.1))
    }
    fn one(&self) -> Fe2 {
        Fe2(self.f.one(), self.f.zero())
    }
    fn is_one(&self, a: &Fe2) -> bool {
        a.0 == self.f.one() && self.f.is_zero(&a.1)
    }
    fn pow(&self, a: &Fe2, e: &BigUint) -> Fe2 {
        let mut acc = self.one();
        for bit in (0..e.bits()).rev() {
            acc = self.mul(&acc, &acc);
            if e.bit(bit) {
                acc = self.mul(&acc, a);
            }
        }
        acc
    }
    /// `a^{-1}` through the norm: `a^{-1} = ā / N(a)`.
    fn inv(&self, a: &Fe2) -> Option<Fe2> {
        let f = self.f;
        let norm = f.sub(&f.mul(&a.0, &a.0), &f.mul(&self.d, &f.mul(&a.1, &a.1)));
        let ni = f.inv(&norm)?;
        Some(Fe2(f.mul(&a.0, &ni), f.neg(&f.mul(&a.1, &ni))))
    }
    /// An element of exact order `n` (a power of two dividing `p² − 1`).
    fn root_of_unity(&self, n: usize) -> Option<Fe2> {
        let p = self.f.modulus();
        let order = p * p - 1u32;
        if n == 0 || !n.is_power_of_two() || (&order % n as u64) != BigUint::from(0u32) {
            return None;
        }
        let e = &order / n as u64;
        for c in 1u64..64 {
            let z = Fe2(self.f.from_u64(c), self.f.one());
            let w = self.pow(&z, &e);
            if n == 1 || !self.is_one(&self.pow(&w, &BigUint::from((n / 2) as u64))) {
                return Some(w);
            }
        }
        None
    }
}

/// In-place transform of length `n = a.len()`, a power of two, with `root`
/// of exact order `n`.
fn ntt_in_place(x: &Ext, a: &mut [Fe2], root: &Fe2) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let w_len = x.pow(root, &BigUint::from((n / len) as u64));
        let half = len / 2;
        let mut tw = Vec::with_capacity(half);
        let mut w = x.one();
        for _ in 0..half {
            tw.push(w);
            w = x.mul(&w, &w_len);
        }
        for start in (0..n).step_by(len) {
            for (k, t) in tw.iter().enumerate() {
                let u = a[start + k];
                let v = x.mul(&a[start + k + half], t);
                a[start + k] = x.add(&u, &v);
                a[start + k + half] = x.sub(&u, &v);
            }
        }
        len <<= 1;
    }
}

/// Truncated product of two `F_p` series by a transform over `F_{p²}`,
/// when that field has the roots of unity; otherwise `None`.
fn mul_trunc_ntt(f: &Field, a: &[Fe], b: &[Fe], len: usize) -> Option<Vec<Fe>> {
    let n = (a.len() + b.len() - 1).next_power_of_two();
    let x = Ext {
        f,
        d: f.from_u64(f.least_nonresidue()),
    };
    let root = x.root_of_unity(n)?;
    let inv_root = x.inv(&root)?;
    let inv_n = f.inv(&f.from_u64(n as u64))?;
    let lift = |s: &[Fe]| -> Vec<Fe2> {
        let mut v = vec![Fe2(f.zero(), f.zero()); n];
        for (i, c) in s.iter().enumerate() {
            v[i] = Fe2(*c, f.zero());
        }
        v
    };
    let mut fa = lift(a);
    let mut fb = lift(b);
    ntt_in_place(&x, &mut fa, &root);
    ntt_in_place(&x, &mut fb, &root);
    for (u, v) in fa.iter_mut().zip(fb.iter()) {
        *u = x.mul(u, v);
    }
    ntt_in_place(&x, &mut fa, &inv_root);
    let mut out = vec![f.zero(); len];
    for (o, c) in out.iter_mut().zip(fa.iter()) {
        let c = x.mul_re(c, &inv_n);
        debug_assert!(f.is_zero(&c.1), "real input must give a real product");
        *o = c.0;
    }
    Some(out)
}

/// Truncated product of two series: schoolbook when short, an NTT when
/// `F_p` has the roots of unity, else Kronecker substitution into one
/// big-integer product.
fn mul_trunc(f: &Field, a: &[Fe], b: &[Fe], len: usize) -> Vec<Fe> {
    let a = &a[..a.len().min(len)];
    let b = &b[..b.len().min(len)];
    if a.is_empty() || b.is_empty() {
        return vec![f.zero(); len];
    }
    if a.len().min(b.len()) < 32 {
        let mut out = vec![f.zero(); len];
        for (i, ai) in a.iter().enumerate() {
            if f.is_zero(ai) {
                continue;
            }
            for (j, bj) in b.iter().enumerate().take(len - i) {
                out[i + j] = f.add(&out[i + j], &f.mul(ai, bj));
            }
        }
        return out;
    }
    if let Some(out) = mul_trunc_ntt(f, a, b, len) {
        return out;
    }
    let pack = |s: &[Fe]| -> BigUint {
        let mut bytes = vec![0u8; s.len() * KRONECKER_SLOT];
        for (i, c) in s.iter().enumerate() {
            let v = f.to_big(c).to_bytes_le();
            bytes[i * KRONECKER_SLOT..i * KRONECKER_SLOT + v.len()].copy_from_slice(&v);
        }
        BigUint::from_bytes_le(&bytes)
    };
    let prod = pack(a) * pack(b);
    let bytes = prod.to_bytes_le();
    let p = f.modulus();
    (0..len)
        .map(|i| {
            let lo = i * KRONECKER_SLOT;
            if lo >= bytes.len() {
                return f.zero();
            }
            let hi = (lo + KRONECKER_SLOT).min(bytes.len());
            f.from_big(&(BigUint::from_bytes_le(&bytes[lo..hi]) % p))
        })
        .collect()
}

/// `Φ_ℓ mod p` for a prime `ℓ < p`.  Returns `Err` if any check fails.
pub fn modular_polynomial(f: &Field, ell: u64) -> Result<ModPoly, String> {
    let l = ell as usize;
    let top = l * (l + 1); // the deepest pole, of Y^{ℓ+1} and X^ℓY^ℓ
    let len = top + CHECK_BAND + 1; // exponents −top ..= CHECK_BAND
    let lenq = len.div_ceil(l); // terms of a series in Q = q^ℓ
    let trace = std::env::var_os("MODPOLY_TRACE").is_some();
    let t0 = std::time::Instant::now();
    let jser = j_series(f, len);
    if trace {
        eprintln!("Φ_{ell}: j-series {} ms", t0.elapsed().as_millis());
    }
    // A[a] = (Σ J_k q^k)^a to `len` terms; Bq[b] = (Σ J_k Q^k)^b to `lenq` terms.
    // With a transform available, the series of j is transformed once and
    // each power costs one pointwise product, one inverse transform and
    // one forward transform of the truncated result.
    let mut apow: Vec<Vec<Fe>> = Vec::with_capacity(l + 2);
    let mut one = vec![f.zero(); len];
    one[0] = f.one();
    apow.push(one);
    let n = (2 * len - 1).next_power_of_two();
    let x = Ext {
        f,
        d: f.from_u64(f.least_nonresidue()),
    };
    match x.root_of_unity(n) {
        Some(root) if len >= 32 => {
            let inv_root = x.inv(&root).expect("root is a unit");
            let inv_n = f.inv(&f.from_u64(n as u64)).expect("n is a unit");
            let lift = |s: &[Fe]| -> Vec<Fe2> {
                let mut v = vec![Fe2(f.zero(), f.zero()); n];
                for (i, c) in s.iter().enumerate() {
                    v[i] = Fe2(*c, f.zero());
                }
                v
            };
            let mut jhat = lift(&jser);
            ntt_in_place(&x, &mut jhat, &root);
            let mut prev = lift(&apow[0]);
            ntt_in_place(&x, &mut prev, &root);
            for _ in 1..=l + 1 {
                let mut prod: Vec<Fe2> = prev.iter().zip(jhat.iter()).map(|(u, v)| x.mul(u, v)).collect();
                ntt_in_place(&x, &mut prod, &inv_root);
                let next: Vec<Fe> = prod[..len].iter().map(|c| f.mul(&c.0, &inv_n)).collect();
                prev = lift(&next);
                ntt_in_place(&x, &mut prev, &root);
                apow.push(next);
            }
        }
        _ => {
            for a in 1..=l + 1 {
                let next = mul_trunc(f, &apow[a - 1], &jser, len);
                apow.push(next);
            }
        }
    }
    let jq: Vec<Fe> = jser[..lenq.min(len)].to_vec();
    let mut bpow: Vec<Vec<Fe>> = Vec::with_capacity(l + 2);
    let mut oneq = vec![f.zero(); lenq];
    oneq[0] = f.one();
    bpow.push(oneq);
    for b in 1..=l + 1 {
        let next = mul_trunc(f, &bpow[b - 1], &jq, lenq);
        bpow.push(next);
    }
    if trace {
        eprintln!("Φ_{ell}: powers {} ms", t0.elapsed().as_millis());
    }
    // Residual R[e + top] = coefficient of q^e of Φ_ℓ(j(q), j(q^ℓ)) so far.
    let mut r = vec![f.zero(); len];
    // Add q^{−pole} · S(q) · T(q^ℓ) to the residual.
    // Add q^{−pole} · S(q) · T(q^ℓ) to the residual.  `T` has at most
    // `lenq` terms, so this is `O(len · lenq)`; a transform only wins
    // beyond ℓ ≈ 300.  The row pole may exceed `top` (block ℓ + 1), so the
    // index is bounded directly: a term lands at r[top − pole + ℓkq + k]
    // and is kept exactly when that index lies in 0..len.
    let add_product = |r: &mut Vec<Fe>, pole: usize, s: &[Fe], t: &[Fe]| {
        let base = top as isize - pole as isize;
        for (kq, tk) in t.iter().enumerate() {
            if f.is_zero(tk) {
                continue;
            }
            let start = base + (kq * l) as isize;
            if start >= len as isize {
                break;
            }
            let k_lo = (-start).max(0) as usize;
            let k_hi = ((len as isize - start).max(0) as usize).min(s.len());
            for (k, sk) in s.iter().enumerate().take(k_hi).skip(k_lo) {
                let idx = (start + k as isize) as usize;
                r[idx] = f.add(&r[idx], &f.mul(tk, sk));
            }
        }
    };
    // Coefficient of q^e in X^a Y^b, for e ≥ −(a + ℓb), from the stored powers.
    let coef_at = |a: usize, b: usize, e: isize| -> Fe {
        let shift = e + (a + l * b) as isize;
        if shift < 0 {
            return f.zero();
        }
        let shift = shift as usize;
        let mut acc = f.zero();
        let mut kq = 0;
        while kq * l <= shift {
            let k = shift - kq * l;
            if k < len && kq < lenq {
                acc = f.add(&acc, &f.mul(&apow[a][k], &bpow[b][kq]));
            }
            kq += 1;
        }
        acc
    };
    let one = f.one();
    let minus_one = f.neg(&one);
    let mut c = vec![vec![f.zero(); l + 2]; l + 2];
    c[l + 1][0] = one;
    c[0][l + 1] = one;
    c[l][l] = minus_one;
    // Block b = ℓ + 1 holds only the known Y^{ℓ+1} and its mirror X^{ℓ+1}.
    // Blocks b = ℓ, …, 0 determine c_{ab}, a ≤ b, in decreasing pole order
    // a + ℓb; c_{ℓℓ} is preset and only applied.
    let mut unknowns = 0usize;
    for b in (0..=l + 1).rev() {
        if b <= l {
            for a in (0..=b).rev() {
                if a == l && b == l {
                    continue;
                }
                let e = -((a + l * b) as isize);
                // Residual at this exponent: completed blocks (in r) plus the
                // already determined members of this block and their mirrors.
                let mut acc = r[(e + top as isize) as usize];
                for a2 in a + 1..=b {
                    let coef = c[a2][b];
                    if f.is_zero(&coef) {
                        continue;
                    }
                    acc = f.add(&acc, &f.mul(&coef, &coef_at(a2, b, e)));
                    if a2 != b {
                        acc = f.add(&acc, &f.mul(&coef, &coef_at(b, a2, e)));
                    }
                }
                // X^a Y^b has leading coefficient one at its own pole.
                let v = f.neg(&acc);
                c[a][b] = v;
                c[b][a] = v;
                unknowns += 1;
            }
        }
        // Row: Σ_{a ≤ b} c_{ab} X^a Y^b = q^{−(b + ℓb)} · (Σ_a c_{ab} q^{b−a} A[a]) · Bq[b].
        let pole = b + l * b;
        let mut row = vec![f.zero(); len];
        for a in 0..=b.min(l + 1) {
            let coef = c[a][b];
            if f.is_zero(&coef) {
                continue;
            }
            let sh = b - a;
            for (k, ak) in apow[a].iter().enumerate().take(len - sh) {
                row[k + sh] = f.add(&row[k + sh], &f.mul(&coef, ak));
            }
        }
        add_product(&mut r, pole, &row, &bpow[b]);
        // Mirror: Σ_{a < b} c_{ab} X^b Y^a = q^{−(b + ℓb)} · A[b] · (Σ_a c_{ab} Q^{b−a} Bq[a]).
        if b >= 1 {
            let mut mirror = vec![f.zero(); lenq];
            for a in 0..b {
                let coef = c[b][a];
                if f.is_zero(&coef) {
                    continue;
                }
                let sh = b - a;
                for (kq, bk) in bpow[a].iter().enumerate().take(lenq.saturating_sub(sh)) {
                    mirror[kq + sh] = f.add(&mirror[kq + sh], &f.mul(&coef, bk));
                }
            }
            add_product(&mut r, pole, &apow[b], &mirror);
        }
    }
    if trace {
        eprintln!("Φ_{ell}: solve {} ms", t0.elapsed().as_millis());
    }
    // Every exponent in range must now vanish.
    for (idx, v) in r.iter().enumerate() {
        if !f.is_zero(v) {
            return Err(format!(
                "Φ_{ell}: q^{} does not vanish after the solve",
                idx as isize - top as isize
            ));
        }
    }
    Ok(ModPoly {
        ell,
        c,
        checks_passed: len - unknowns,
    })
}

impl ModPoly {
    /// `Φ_ℓ(x, Y)` as a polynomial in `Y`.
    pub fn at_x(&self, f: &Field, x: &Fe) -> Poly {
        let n = self.c.len();
        let mut xp = vec![f.one(); n];
        for i in 1..n {
            xp[i] = f.mul(&xp[i - 1], x);
        }
        let out = (0..n)
            .map(|b| {
                (0..n).fold(f.zero(), |acc, a| {
                    f.add(&acc, &f.mul(&self.c[a][b], &xp[a]))
                })
            })
            .collect();
        poly::trim(f, out)
    }

    /// `∂^{dx+dy} Φ / ∂X^{dx} ∂Y^{dy}` at `(x, y)`, for `dx + dy ≤ 2`.
    pub fn partial(&self, f: &Field, x: &Fe, y: &Fe, dx: u32, dy: u32) -> Fe {
        let n = self.c.len();
        let falling =
            |i: usize, d: u32| -> u64 { (0..d as usize).map(|t| (i - t) as u64).product::<u64>() };
        let mut acc = f.zero();
        for a in (dx as usize)..n {
            let xa = f.mul(
                &f.from_u64(falling(a, dx)),
                &f.pow(x, &((a - dx as usize) as u64).into()),
            );
            for b in (dy as usize)..n {
                if f.is_zero(&self.c[a][b]) {
                    continue;
                }
                let yb = f.mul(
                    &f.from_u64(falling(b, dy)),
                    &f.pow(y, &((b - dy as usize) as u64).into()),
                );
                acc = f.add(&acc, &f.mul(&self.c[a][b], &f.mul(&xa, &yb)));
            }
        }
        acc
    }
}

#[cfg(test)]
mod reference {
    //! The original schoolbook construction, kept only as the oracle the
    //! fast construction is compared against.
    use super::*;

    /// Truncated product of two series.
    fn mul_trunc_schoolbook(f: &Field, a: &[Fe], b: &[Fe], len: usize) -> Vec<Fe> {
        let mut out = vec![f.zero(); len];
        for (i, ai) in a.iter().enumerate().take(len) {
            if f.is_zero(ai) {
                continue;
            }
            for (j, bj) in b.iter().enumerate().take(len - i) {
                out[i + j] = f.add(&out[i + j], &f.mul(ai, bj));
            }
        }
        out
    }

    /// The original schoolbook construction, kept as the reference the fast one is tested against.
    pub fn modular_polynomial_schoolbook(f: &Field, ell: u64) -> Result<ModPoly, String> {
        let l = ell as usize;
        let top = l * (l + 1); // the deepest pole, of Y^{ℓ+1} and X^ℓY^ℓ
        let len = top + CHECK_BAND + 1; // exponents −top ..= CHECK_BAND
        let jser = j_series(f, len);
        // A[a] = (Σ J_k q^k)^a, B[b] = (Σ J_k q^{ℓk})^b, both truncated.
        let mut apow = vec![{
            let mut one = vec![f.zero(); len];
            one[0] = f.one();
            one
        }];
        for a in 1..=l + 1 {
            let next = mul_trunc_schoolbook(f, &apow[a - 1], &jser, len);
            apow.push(next);
        }
        let mut jl = vec![f.zero(); len];
        for k in 0..len {
            if k * l < len {
                jl[k * l] = jser[k];
            }
        }
        let mut bpow = vec![apow[0].clone()];
        for b in 1..=l + 1 {
            let next = mul_trunc_schoolbook(f, &bpow[b - 1], &jl, len);
            bpow.push(next);
        }
        // R[e + top] accumulates the coefficient of q^e, e ∈ [−top, CHECK_BAND].
        let mut r = vec![f.zero(); len];
        // X^aY^b = q^{−pole} · A[a]·B[b]; its q^e lands at index e + top, so
        // only the first `pole + CHECK_BAND + 1` terms of the product matter.
        // B[b] is a series in q^ℓ: iterate its nonzero terms only.
        let add_monomial = |r: &mut Vec<Fe>, a: usize, b: usize, coef: &Fe| {
            let pole = a + l * b;
            let need = (pole + CHECK_BAND + 1).min(len);
            let base = top - pole;
            for i in (0..need).step_by(l) {
                let bi = bpow[b][i];
                if f.is_zero(&bi) {
                    continue;
                }
                let cb = f.mul(coef, &bi);
                for (k, ak) in apow[a].iter().enumerate().take(need - i) {
                    let idx = base + i + k;
                    r[idx] = f.add(&r[idx], &f.mul(&cb, ak));
                }
            }
        };
        let one = f.one();
        let minus_one = f.neg(&one);
        add_monomial(&mut r, l + 1, 0, &one);
        add_monomial(&mut r, 0, l + 1, &one);
        add_monomial(&mut r, l, l, &minus_one);
        // The unknown whose leading pole is P, if any.
        let mut by_pole: Vec<Option<(usize, usize)>> = vec![None; top + 1];
        for b in 0..=l {
            for a in 0..=b {
                if a == l && b == l {
                    continue;
                }
                let pole = a + l * b;
                assert!(by_pole[pole].is_none(), "pole orders are distinct");
                by_pole[pole] = Some((a, b));
            }
        }
        let mut c = vec![vec![f.zero(); l + 2]; l + 2];
        c[l + 1][0] = one;
        c[0][l + 1] = one;
        c[l][l] = minus_one;
        let mut checks = 0usize;
        for idx in 0..len {
            let pole = top as isize - idx as isize; // e = −pole
            let slot = if pole >= 0 {
                by_pole[pole as usize]
            } else {
                None
            };
            match slot {
                Some((a, b)) => {
                    let coef = f.neg(&r[idx]);
                    c[a][b] = coef;
                    c[b][a] = coef;
                    add_monomial(&mut r, a, b, &coef);
                    if a != b {
                        add_monomial(&mut r, b, a, &coef);
                    }
                    if !f.is_zero(&r[idx]) {
                        return Err(format!(
                            "Φ_{ell}: unknown ({a},{b}) left q^{} nonzero",
                            -pole
                        ));
                    }
                }
                None => {
                    if !f.is_zero(&r[idx]) {
                        return Err(format!("Φ_{ell}: q^{} does not vanish", -pole));
                    }
                    checks += 1;
                }
            }
        }
        Ok(ModPoly {
            ell,
            c,
            checks_passed: checks,
        })
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fast_construction_matches_the_schoolbook_reference() {
        for p in [BigUint::from(1009u64), BigUint::parse_bytes(
            b"ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap()] {
            let f = Field::new(&p).unwrap();
            for ell in [3u64, 5, 7, 11, 13] {
                let fast = modular_polynomial(&f, ell).unwrap();
                let slow = reference::modular_polynomial_schoolbook(&f, ell).unwrap();
                assert_eq!(fast.c, slow.c, "Φ_{ell} mod {p}");
                assert_eq!(fast.checks_passed, slow.checks_passed, "checks for Φ_{ell}");
            }
        }
    }
    use crate::cryptanalysis::modular_polynomial::{phi_2, phi_3};
    use num_bigint::{BigInt, BigUint};
    use num_traits::Signed;

    fn reduce(f: &Field, v: &BigInt) -> Fe {
        let m = f.from_big(&v.abs().to_biguint().unwrap());
        if v.is_negative() {
            f.neg(&m)
        } else {
            m
        }
    }

    #[test]
    fn j_series_matches_the_integer_expansion() {
        let f = Field::new(&BigUint::from(1_000_000_007u64)).unwrap();
        let j = j_series(&f, 4);
        // j = q^-1 + 744 + 196884 q + 21493760 q^2 + 864299970 q^3.
        let want = [1u64, 744, 196_884, 21_493_760, 864_299_970];
        for (k, w) in want.iter().take(4).enumerate() {
            assert_eq!(f.to_big(&j[k]), BigUint::from(*w) % 1_000_000_007u64);
        }
    }

    #[test]
    fn phi2_and_phi3_match_the_literature_tables() {
        let p = BigUint::parse_bytes(
            b"ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let f = Field::new(&p).unwrap();
        for table in [phi_2(), phi_3()] {
            let m = modular_polynomial(&f, table.l).unwrap();
            let l = table.l as u32;
            for a in 0..=l + 1 {
                for b in 0..=l + 1 {
                    assert_eq!(
                        m.c[a as usize][b as usize],
                        reduce(&f, &table.coeff(a, b)),
                        "Φ_{l} coefficient ({a},{b})"
                    );
                }
            }
        }
    }
}
