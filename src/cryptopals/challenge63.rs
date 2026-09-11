//! # Challenge 63 — GCM Authentication-Key Recovery via Nonce Reuse
//!
//! GMAC computes
//!
//! ```text
//!   t  =  AAD₁·H^N + AAD₂·H^(N-1) + … + len·H + s
//! ```
//!
//! where `H = AES_K(0)` is the authentication key and
//! `s = AES_K(nonce ‖ 1)` is the nonce-derived mask.  Reuse the
//! nonce across two messages and `s` cancels when you XOR the tags:
//!
//! ```text
//!   t₁ ⊕ t₂  =  Σ (Δaᵢ)·H^i
//! ```
//!
//! H is now a root of a known polynomial over GF(2¹²⁸).  Factor the
//! polynomial, harvest the degree-1 factors, and you've leaked the
//! authentication key in one shot — game over for forging on that
//! GCM key.
//!
//! ## What this module provides
//!
//! * `GF2_128` — element type with add/mul/inv operating over
//!   `x^128 + x^7 + x^2 + x + 1`.
//! * `Poly` — polynomials with `GF2_128` coefficients with the
//!   standard ring operations + Euclidean div/mod and gcd.
//! * `factor` — distinct-degree + Cantor–Zassenhaus equal-degree
//!   factoring, returning the list of irreducible factors.  We only
//!   actually need the degree-1 ones (the roots), so a thin
//!   `roots` wrapper extracts those.
//! * `attack` — given two messages encrypted under the same
//!   `(K, nonce)`, recover candidate H values and verify.

use crate::cryptopals::Report;
use crate::symmetric::aes::{aes_ctr, encrypt_block, AesKey};
use num_bigint::BigUint;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

// ── GF(2^128) element ─────────────────────────────────────────────

/// 128-bit field element packed as a little-endian `u128`.  Cryptopals
/// GHASH actually uses a *big-endian, bit-reversed-per-byte* layout;
/// for educational clarity we operate in plain GF(2)\[x]/(p(x)) where
/// `x^k` is bit `k`.  The conversion to GCM byte order happens at the
/// edges of the attack (functions `pack`/`unpack`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Gf128(pub u128);

impl Gf128 {
    pub const fn zero() -> Self {
        Gf128(0)
    }
    pub const fn one() -> Self {
        Gf128(1)
    }
    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Reduction polynomial: `x^128 + x^7 + x^2 + x + 1`.
    /// Represented in u128 form: low 128 bits = x^0..x^127.
    const REDUCTION_LOW: u128 = (1u128 << 7) | (1u128 << 2) | (1u128 << 1) | 1u128;

    /// GF(2^128) multiplication with carryless propagation.  O(128).
    pub fn mul(self, rhs: Self) -> Self {
        // Schoolbook with reduction folded in.  Standard
        // "two-limb shift-XOR" routine.
        let mut a = self.0;
        let mut b = rhs.0;
        let mut p: u128 = 0;
        for _ in 0..128 {
            if b & 1 != 0 {
                p ^= a;
            }
            let hi = (a >> 127) & 1;
            a <<= 1;
            if hi != 0 {
                a ^= Self::REDUCTION_LOW;
            }
            b >>= 1;
        }
        Gf128(p)
    }

    pub fn add(self, rhs: Self) -> Self {
        Gf128(self.0 ^ rhs.0)
    }

    pub fn sub(self, rhs: Self) -> Self {
        self.add(rhs)
    }

    /// Inverse via Fermat: `a^(2^128 − 2)` since the multiplicative
    /// group has order `2^128 − 1`.
    pub fn inv(self) -> Self {
        assert!(!self.is_zero(), "Gf128::inv(0)");
        // a^(2^128 - 2) = a · (a^(2^128 - 1) / a) but use simpler:
        // square-and-multiply over the binary expansion of 2^128 - 2.
        // 2^128 - 2 = 0b1...10 (127 ones, then a zero), so the
        // routine is: square 128 times, multiplying when the bit is 1.
        let mut acc = Gf128::one();
        let mut base = self;
        let mut e: u128 = (u128::MAX) ^ 1; // 2^128 - 2 = 1...10
        for _ in 0..128 {
            if e & 1 != 0 {
                acc = acc.mul(base);
            }
            base = base.mul(base);
            e >>= 1;
        }
        // We covered 128 of the 128 bits; the leading `1`s at top
        // of e weren't shifted yet, but actually they were — the
        // loop runs 128 times and `e` is 128 bits wide.  Done.
        acc
    }

    pub fn pow(self, exp: &BigUint) -> Self {
        let bytes = exp.to_bytes_le();
        let mut acc = Gf128::one();
        let mut base = self;
        for byte in bytes {
            let mut b = byte;
            for _ in 0..8 {
                if b & 1 != 0 {
                    acc = acc.mul(base);
                }
                base = base.mul(base);
                b >>= 1;
            }
        }
        acc
    }
}

/// Convert a 16-byte GCM block to a `Gf128` element.  GCM byte order:
/// the first bit (most significant) of the first byte is the
/// coefficient of x^127.  We convert by reversing each byte's bits
/// then interpreting little-endian.
pub fn pack(block: &[u8; 16]) -> Gf128 {
    let mut v: u128 = 0;
    for (i, &b) in block.iter().enumerate() {
        v |= (reverse_bits(b) as u128) << (i * 8);
    }
    Gf128(v)
}

pub fn unpack(g: Gf128) -> [u8; 16] {
    let mut out = [0u8; 16];
    let mut v = g.0;
    for i in 0..16 {
        out[i] = reverse_bits((v & 0xff) as u8);
        v >>= 8;
    }
    out
}

fn reverse_bits(b: u8) -> u8 {
    let mut r: u8 = 0;
    for i in 0..8 {
        if b & (1 << i) != 0 {
            r |= 1 << (7 - i);
        }
    }
    r
}

// ── Polynomials over GF(2^128) ────────────────────────────────────

/// `coeffs[i]` is the coefficient of `x^i`.  Degree is `len - 1`
/// when the trailing coefficient is non-zero (or the poly is zero).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    pub coeffs: Vec<Gf128>,
}

impl Poly {
    pub fn new(mut coeffs: Vec<Gf128>) -> Self {
        while coeffs.len() > 1 && coeffs.last().unwrap().is_zero() {
            coeffs.pop();
        }
        Poly { coeffs }
    }

    pub fn zero() -> Self {
        Poly { coeffs: vec![Gf128::zero()] }
    }
    pub fn one() -> Self {
        Poly { coeffs: vec![Gf128::one()] }
    }

    pub fn deg(&self) -> i32 {
        if self.coeffs.len() == 1 && self.coeffs[0].is_zero() {
            return -1;
        }
        (self.coeffs.len() - 1) as i32
    }

    pub fn leading(&self) -> Gf128 {
        *self.coeffs.last().unwrap()
    }

    pub fn add(&self, rhs: &Self) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = vec![Gf128::zero(); n];
        for i in 0..self.coeffs.len() {
            out[i] = out[i].add(self.coeffs[i]);
        }
        for i in 0..rhs.coeffs.len() {
            out[i] = out[i].add(rhs.coeffs[i]);
        }
        Poly::new(out)
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(rhs)
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = vec![Gf128::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for i in 0..self.coeffs.len() {
            for j in 0..rhs.coeffs.len() {
                out[i + j] = out[i + j].add(self.coeffs[i].mul(rhs.coeffs[j]));
            }
        }
        Poly::new(out)
    }

    pub fn scale(&self, c: Gf128) -> Self {
        Poly::new(self.coeffs.iter().map(|x| x.mul(c)).collect())
    }

    /// Polynomial division: returns `(q, r)` with `self = q·div + r`,
    /// `deg(r) < deg(div)`.
    pub fn divmod(&self, div: &Self) -> (Self, Self) {
        let mut r = self.clone();
        let mut q_coeffs = vec![Gf128::zero(); (self.deg() - div.deg() + 1).max(0) as usize];
        let div_lead_inv = div.leading().inv();
        while r.deg() >= div.deg() && !(r.deg() == -1) {
            let shift = (r.deg() - div.deg()) as usize;
            let factor = r.leading().mul(div_lead_inv);
            q_coeffs[shift] = factor;
            // r -= factor · x^shift · div
            for i in 0..div.coeffs.len() {
                let idx = i + shift;
                r.coeffs[idx] = r.coeffs[idx].add(factor.mul(div.coeffs[i]));
            }
            // Trim leading zero.
            while r.coeffs.len() > 1 && r.coeffs.last().unwrap().is_zero() {
                r.coeffs.pop();
            }
        }
        (Poly::new(q_coeffs), r)
    }

    /// Make monic (leading coefficient = 1).
    pub fn monic(&self) -> Self {
        if self.deg() < 0 {
            return self.clone();
        }
        let inv = self.leading().inv();
        self.scale(inv)
    }

    pub fn gcd(a: &Self, b: &Self) -> Self {
        if b.deg() < 0 {
            return a.monic();
        }
        let (_, r) = a.divmod(b);
        Self::gcd(b, &r)
    }

    /// Compute `self^e mod m` for `BigUint` exponent.
    pub fn modpow(&self, e: &BigUint, m: &Self) -> Self {
        let mut acc = Poly::one();
        let mut base = self.divmod(m).1;
        let bytes = e.to_bytes_le();
        for byte in bytes {
            let mut b = byte;
            for _ in 0..8 {
                if b & 1 != 0 {
                    acc = acc.mul(&base).divmod(m).1;
                }
                base = base.mul(&base).divmod(m).1;
                b >>= 1;
            }
        }
        acc
    }

    /// Evaluate at `x`.
    pub fn eval(&self, x: Gf128) -> Gf128 {
        let mut acc = Gf128::zero();
        for &c in self.coeffs.iter().rev() {
            acc = acc.mul(x).add(c);
        }
        acc
    }
}

// ── Factoring: distinct-degree + Cantor–Zassenhaus ────────────────

/// Distinct-degree factorisation.  Returns `[(g_d, d)]` where each
/// `g_d` is the product of all irreducible factors of degree `d`.
pub fn distinct_degree_factor(f: &Poly) -> Vec<(Poly, usize)> {
    let mut out = Vec::new();
    let mut f = f.clone();
    let mut d = 1usize;
    // x^(2^128) - x mod f
    let mut h = Poly::new(vec![Gf128::zero(), Gf128::one()]); // x
    while f.deg() >= 2 * (d as i32) {
        // h := h^(2^128) mod f  ≡ x^(2^(128·d)) mod f after d iterations
        h = h.modpow(&(BigUint::from(1u32) << 128), &f);
        // g_d := gcd(h - x, f)
        let h_minus_x = h.sub(&Poly::new(vec![Gf128::zero(), Gf128::one()]));
        let g = Poly::gcd(&h_minus_x, &f).monic();
        if g.deg() > 0 {
            out.push((g.clone(), d));
            let (q, _) = f.divmod(&g);
            f = q;
            // Restart h relative to the new f.
            let (_q, r) = h.divmod(&f);
            h = r;
        }
        d += 1;
        if d > 4 {
            // Higher-degree factors don't matter for the cryptopals
            // attack — we only want degree-1 roots.
            break;
        }
    }
    if f.deg() > 0 {
        out.push((f.clone(), f.deg() as usize));
    }
    out
}

/// Equal-degree factorisation: given `g_d` of degree `r·d`, split it
/// into `r` irreducible degree-`d` factors using Cantor–Zassenhaus
/// over GF(2^128).  Because we operate in characteristic 2, we use
/// the trace-map trick instead of the `(q^d − 1)/3` exponent
/// described in the cryptopals spec (which only works in odd
/// characteristic).
pub fn equal_degree_factor(g: &Poly, d: usize, seed: u64) -> Vec<Poly> {
    let mut s = vec![g.clone()];
    let mut rng = StdRng::seed_from_u64(seed);
    let r = (g.deg() as usize) / d.max(1);
    while s.len() < r {
        let mut h = random_poly(g.deg() as usize, &mut rng);
        // Trace polynomial: T(h) = h + h^2 + h^4 + … + h^(2^(128·d − 1)) mod g.
        let mut t = Poly::zero();
        let mut hh = h.clone();
        for _ in 0..(128 * d) {
            t = t.add(&hh);
            hh = hh.mul(&hh).divmod(g).1;
        }
        t = t.divmod(g).1;
        h = t;
        for f in s.clone() {
            if f.deg() <= d as i32 {
                continue;
            }
            let gcd = Poly::gcd(&h, &f);
            if gcd.deg() > 0 && gcd.deg() < f.deg() {
                let (q, _) = f.divmod(&gcd);
                s.retain(|p| p != &f);
                s.push(gcd.monic());
                s.push(q.monic());
            }
        }
    }
    s
}

/// Roots of `f` (degree-1 factors).  Returns the set of `r` such
/// that `f(r) = 0`.
pub fn roots(f: &Poly) -> Vec<Gf128> {
    let monic = f.monic();
    let dds = distinct_degree_factor(&monic);
    let mut out = Vec::new();
    for (g, d) in dds {
        if d != 1 {
            continue;
        }
        // Equal-degree split into linear factors.
        if g.deg() == 1 {
            // y + c → root c
            out.push(g.coeffs[0]);
            continue;
        }
        let factors = equal_degree_factor(&g, 1, 42);
        for f in factors {
            if f.deg() == 1 {
                out.push(f.coeffs[0]);
            }
        }
    }
    out
}

fn random_poly(deg: usize, rng: &mut StdRng) -> Poly {
    let mut c = vec![Gf128::zero(); deg];
    for x in c.iter_mut() {
        *x = Gf128(rng.gen::<u128>());
    }
    c.push(Gf128::one());
    Poly::new(c)
}

// ── The attack ────────────────────────────────────────────────────

/// Build the GCM-style polynomial whose root is the authentication
/// key, given two ciphertext / tag pairs encrypted under the same
/// `(K, nonce)` and (for simplicity) no AAD.
///
/// Polynomial form (with `ct` reversed so block n is x^1, block 1 is x^n):
///
/// ```text
///   f(y) = Σ (ct1[i] ⊕ ct2[i]) · y^(n - i)  ⊕  (len ⊕ len) · y  ⊕  (t1 ⊕ t2)
/// ```
///
/// With both ciphertexts the *same length* the length-block delta
/// is zero, so the linear term drops out.
pub fn build_nonce_reuse_poly(
    ct1: &[u8],
    tag1: &[u8; 16],
    ct2: &[u8],
    tag2: &[u8; 16],
) -> Poly {
    assert_eq!(ct1.len(), ct2.len(), "messages must be same length");
    let n_blocks = (ct1.len() + 15) / 16;
    let mut coeffs = vec![Gf128::zero(); n_blocks + 2];
    // y^0 = t1 ⊕ t2
    let mut t_delta = [0u8; 16];
    for i in 0..16 {
        t_delta[i] = tag1[i] ^ tag2[i];
    }
    coeffs[0] = pack(&t_delta);
    // y^1 = length block diff = 0 (same length).
    // y^2.. = (ct1[i] ⊕ ct2[i]) for i descending.
    for i in 0..n_blocks {
        let mut blk = [0u8; 16];
        let start = i * 16;
        let end = (start + 16).min(ct1.len());
        let len = end - start;
        for j in 0..len {
            blk[j] = ct1[start + j] ^ ct2[start + j];
        }
        let exp = n_blocks + 1 - i;
        coeffs[exp] = pack(&blk);
    }
    Poly::new(coeffs)
}

pub fn run() -> Report {
    let mut r = Report::new(63, "GCM Authentication-Key Recovery via Nonce Reuse");

    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let nonce = [0u8; 12];
    let msg1 = b"the magic words are squeamish ossifrage!!";
    let msg2 = b"we attack at dawn, repeat, attack at dawn";
    // Same length ensures the length-block delta is zero, simplifying.
    assert_eq!(msg1.len(), msg2.len());
    let aad = b"";
    let ct1_tag = crate::symmetric::aes::aes_gcm_encrypt(msg1, &key, &nonce, aad);
    let ct2_tag = crate::symmetric::aes::aes_gcm_encrypt(msg2, &key, &nonce, aad);
    let (ct1, t1) = ct1_tag.split_at(ct1_tag.len() - 16);
    let (ct2, t2) = ct2_tag.split_at(ct2_tag.len() - 16);
    let mut tag1 = [0u8; 16];
    let mut tag2 = [0u8; 16];
    tag1.copy_from_slice(t1);
    tag2.copy_from_slice(t2);

    let f = build_nonce_reuse_poly(ct1, &tag1, ct2, &tag2);
    r.line(format!("Polynomial degree: {}", f.deg()));
    let candidates = roots(&f);
    r.line(format!("Candidate roots: {}", candidates.len()));

    // The true H is AES_K(0).
    let h_truth = pack(&encrypt_block(&[0u8; 16], &key));
    r.line(format!(
        "True H = AES_K(0)  : {:032x}",
        u128::from_le_bytes(unpack(h_truth))
    ));
    let mut found = false;
    for c in &candidates {
        if c == &h_truth {
            found = true;
        }
    }
    r.line(format!("True H among roots: {}", found));
    let _ = aes_ctr;
    if found {
        r.succeed()
    } else {
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gf128_inverse() {
        let a = Gf128(0x123456789abcdef0_fedcba9876543210);
        let inv = a.inv();
        assert_eq!(a.mul(inv), Gf128::one());
    }

    #[test]
    fn pack_round_trip() {
        let block = *b"YELLOW SUBMARINE";
        let g = pack(&block);
        let back = unpack(g);
        assert_eq!(back, block);
    }

    #[test]
    fn poly_eval_at_root() {
        // f(y) = (y - 3)(y - 5)  with coefficients in Gf128
        let a = Gf128(3);
        let b = Gf128(5);
        let f1 = Poly::new(vec![a, Gf128::one()]); // (y + 3)
        let f2 = Poly::new(vec![b, Gf128::one()]); // (y + 5)
        let prod = f1.mul(&f2);
        assert_eq!(prod.eval(a), Gf128::zero());
        assert_eq!(prod.eval(b), Gf128::zero());
    }

    #[test]
    fn factor_finds_linear_factors() {
        let a = Gf128(0xDEAD_BEEF_CAFE_BABE_0001_0002_0003_0004);
        let b = Gf128(0x1111_2222_3333_4444_5555_6666_7777_8888);
        let f1 = Poly::new(vec![a, Gf128::one()]);
        let f2 = Poly::new(vec![b, Gf128::one()]);
        let prod = f1.mul(&f2);
        let rs = roots(&prod);
        assert!(rs.contains(&a) && rs.contains(&b), "got roots: {:?}", rs);
    }
}
