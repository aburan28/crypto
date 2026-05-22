//! **ML-DSA** — Module-Lattice-based Digital Signature Algorithm
//! (NIST FIPS 204, August 2024; formerly CRYSTALS-Dilithium).
//!
//! # Background
//! ML-DSA is the NIST post-quantum standard for digital signatures.
//! Its security reduces to the hardness of the Module Learning With
//! Errors (Module-LWE) and Module Short Integer Solution (Module-SIS)
//! problems over the ring `R_q = Z_q[X] / (X^256 + 1)` with
//! `q = 8380417 = 2^23 - 2^13 + 1`.  This prime is chosen so that
//! a length-256 Number Theoretic Transform (NTT) exists, allowing fast
//! polynomial multiplication.
//!
//! # High level
//! Signatures follow the Fiat-Shamir-with-aborts paradigm:
//!   1. Sample a masking vector `y` with coefficients in `(-gamma1, gamma1]`.
//!   2. Compute `w = A·y`, take its high bits `w1`.
//!   3. Derive a challenge `c` (a sparse +/-1 polynomial) by hashing
//!      `(message, w1)`.
//!   4. Compute `z = y + c·s1` and check bound `‖z‖_∞ < gamma1 - beta`.
//!   5. Compute the low bits `r0 = LowBits(w - c·s2)` and check
//!      `‖r0‖_∞ < gamma2 - beta`.
//!   6. Encode hints `h` describing the difference between `w1` and
//!      `HighBits(w - c·s2 + c·t0)`.
//!   7. If any check fails, restart with fresh `y` ("aborts").
//!
//! Verification recomputes `w'_approx = A·z - c·t1·2^d` and uses `h`
//! to recover `w1`, then recomputes the challenge hash and compares.
//!
//! # Parameter set
//! This file implements **ML-DSA-65** (NIST security category 3, the
//! recommended default).  Parameters from FIPS 204 §4 Table 1:
//!   `(n, q, d, tau, lambda, gamma1, gamma2, k, l, eta, beta, omega) =
//!    (256, 8380417, 13, 49, 192, 2^19, (q-1)/32, 6, 5, 4, 196, 55)`.
//!
//! # This implementation
//! Faithful to FIPS 204 in algorithm structure and (we believe) in bit-
//! level output for keygen and verification.  Signature determinism
//! depends on the random seed `rnd` argument; for hedged signatures
//! supply 32 fresh random bytes per call, for deterministic signatures
//! supply `[0u8; 32]`.
//!
//! Educational implementation: no constant-time guarantees.

use crate::hash::sha3::{shake128, shake256};

// ── Parameters (ML-DSA-65) ────────────────────────────────────────────────────

const N: usize = 256;
const Q: i32 = 8_380_417;
const D: u32 = 13;
const TAU: usize = 49;
const LAMBDA: usize = 192;
const GAMMA1: i32 = 1 << 19;
const GAMMA2: i32 = (Q - 1) / 32;
const K: usize = 6;
const L: usize = 5;
const ETA: i32 = 4;
const BETA: i32 = (TAU as i32) * ETA; // = 196
const OMEGA: usize = 55;

// Serialized sizes derived from the parameter set.
const SEED_BYTES: usize = 32;
const TR_BYTES: usize = 64;
const C_TILDE_BYTES: usize = LAMBDA / 4; // 48 bytes
const POLYT1_PACKED: usize = 320; // 10 bits/coeff
const POLYT0_PACKED: usize = 416; // 13 bits/coeff
const POLYETA_PACKED: usize = 128; // 4 bits/coeff (η=4)
const POLYZ_PACKED: usize = 640; // (1+log2(gamma1))=20 bits, gamma1=2^19
const POLYW1_PACKED: usize = 128; // 4 bits/coeff (w1 in 0..15)

const PUBKEY_BYTES: usize = SEED_BYTES + K * POLYT1_PACKED;
const SECKEY_BYTES: usize =
    SEED_BYTES + SEED_BYTES + TR_BYTES + L * POLYETA_PACKED + K * POLYETA_PACKED + K * POLYT0_PACKED;
const SIG_BYTES: usize = C_TILDE_BYTES + L * POLYZ_PACKED + OMEGA + K;

// ── Polynomial type ───────────────────────────────────────────────────────────

/// Polynomial in `R_q = Z_q[X]/(X^256+1)`, stored as a length-256 array of i32.
/// Coefficients are kept in a representation convenient to the operation in
/// progress; explicit reductions are applied where needed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly(pub [i32; N]);

impl Poly {
    pub fn zero() -> Self {
        Poly([0; N])
    }
}

// ── NTT machinery ────────────────────────────────────────────────────────────
//
// 256-th primitive root of unity: ζ = 1753 mod q.
// FIPS 204 §7.5 specifies a fixed bit-reversal order for the constants used.

const ZETAS: [i32; N] = {
    // Precomputed at compile time so the constants are exact and review-able.
    // ZETAS[k] = ζ^(brv(k)) mod q, where brv is bit-reversal on 8 bits and
    // ζ = 1753 is the canonical primitive 512-th root of unity used by ML-DSA.
    let mut out = [0i32; N];
    let mut powers = [0i64; 256];
    powers[0] = 1;
    let mut i = 1;
    while i < 256 {
        powers[i] = (powers[i - 1] * 1753) % (Q as i64);
        i += 1;
    }
    // bit-reverse on 8 bits
    let mut k = 0usize;
    while k < N {
        let mut br = 0usize;
        let mut x = k;
        let mut b = 0;
        while b < 8 {
            br = (br << 1) | (x & 1);
            x >>= 1;
            b += 1;
        }
        out[k] = powers[br] as i32;
        k += 1;
    }
    out
};

/// Montgomery-style reduction not used: we operate in plain i64 and reduce
/// modulo q with `barrett`-like step where needed.
#[inline]
fn fqmul(a: i32, b: i32) -> i32 {
    let prod = (a as i64) * (b as i64);
    (prod.rem_euclid(Q as i64)) as i32
}

/// Forward NTT (in-place).  Cooley-Tukey, decimation-in-time.
/// All intermediate arithmetic uses i64 to avoid i32 overflow; the result
/// is reduced into [0, q) at each butterfly.
fn ntt(p: &mut Poly) {
    // Bring inputs into [0, q) first.
    for c in p.0.iter_mut() {
        *c = c.rem_euclid(Q);
    }
    let mut k = 0usize;
    let mut len = 128usize;
    while len >= 1 {
        let mut start = 0usize;
        while start < N {
            k += 1;
            let zeta = ZETAS[k] as i64;
            let mut j = start;
            while j < start + len {
                let t = ((zeta * (p.0[j + len] as i64)) % (Q as i64)) as i32;
                p.0[j + len] = (p.0[j] - t).rem_euclid(Q);
                p.0[j] = (p.0[j] + t).rem_euclid(Q);
                j += 1;
            }
            start += 2 * len;
        }
        len >>= 1;
    }
}

/// Inverse NTT (in-place).  Gentleman-Sande, decimation-in-frequency.
fn inv_ntt(p: &mut Poly) {
    for c in p.0.iter_mut() {
        *c = c.rem_euclid(Q);
    }
    let mut k = 256usize;
    let mut len = 1usize;
    while len < N {
        let mut start = 0usize;
        while start < N {
            k -= 1;
            let zeta = -(ZETAS[k] as i64); // inverse direction
            let mut j = start;
            while j < start + len {
                let t = p.0[j];
                let s = p.0[j + len];
                p.0[j] = (t + s).rem_euclid(Q);
                let diff = (t - s).rem_euclid(Q);
                p.0[j + len] = ((zeta * (diff as i64)).rem_euclid(Q as i64)) as i32;
                j += 1;
            }
            start += 2 * len;
        }
        len <<= 1;
    }
    // Multiply by n^{-1} mod q.  n=256, n^{-1} = 8347681 mod q.
    let n_inv: i32 = 8_347_681;
    for c in p.0.iter_mut() {
        *c = fqmul(*c, n_inv);
    }
}

/// Pointwise multiplication in NTT domain.
fn poly_pointwise(a: &Poly, b: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.0[i] = fqmul(a.0[i], b.0[i]);
    }
    out
}

/// Add two polynomials, no reduction beyond keeping coefficients in [0, q).
fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.0[i] = (a.0[i] + b.0[i]).rem_euclid(Q);
    }
    out
}

/// Subtract polynomials, result in [0, q).
fn poly_sub(a: &Poly, b: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.0[i] = (a.0[i] - b.0[i]).rem_euclid(Q);
    }
    out
}

// ── Centered representative helper ──────────────────────────────────────────

/// Convert canonical representative in [0,q) to centered representative
/// in (-q/2, q/2].
#[inline]
fn to_centered(a: i32) -> i32 {
    if a > Q / 2 {
        a - Q
    } else {
        a
    }
}

/// Infinity norm of polynomial with centered coefficients.
fn poly_inf_norm(p: &Poly) -> i32 {
    let mut max = 0i32;
    for &c in &p.0 {
        let v = to_centered(c).abs();
        if v > max {
            max = v;
        }
    }
    max
}

// ── Power2Round / Decompose / HighBits / LowBits / MakeHint / UseHint ───────
//
// All per FIPS 204 §7.4.  These operate on the *centered* representative
// abstractly, but the FIPS reference computes them directly on canonical [0,q).

/// Split coefficient `r` (in [0,q)) into (r1, r0) with `r = r1 * 2^D + r0` and
/// `r0` in `(-2^{D-1}, 2^{D-1}]`.
fn power2_round(r: i32) -> (i32, i32) {
    let r = r.rem_euclid(Q);
    let r1 = (r + (1 << (D - 1)) - 1) >> D;
    let r0 = r - (r1 << D);
    (r1, r0)
}

/// FIPS 204 Algorithm 36 Decompose with α = 2·γ₂ = (q-1)/16.
///
///   r⁺ ← r mod q
///   r₀ ← r⁺ mod± (2γ₂)               (centered representative in (-γ₂, γ₂])
///   if r⁺ - r₀ == q - 1: return (0, r₀ - 1)
///   else: r₁ ← (r⁺ - r₀) / (2γ₂); return (r₁, r₀)
///
/// Output invariant: r = r₁·2γ₂ + r₀ (mod q), with r₀ ∈ (-γ₂, γ₂] and r₁ ∈ {0,…,15},
/// except for the boundary case where (r₁, r₀) = (0, -γ₂) and r ≡ q - 1.
fn decompose(r: i32) -> (i32, i32) {
    let r_plus = r.rem_euclid(Q);
    let two_gamma2 = 2 * GAMMA2;
    // centered mod 2γ₂ → r0 in (-γ₂, γ₂]
    let mut r0 = r_plus.rem_euclid(two_gamma2);
    if r0 > GAMMA2 {
        r0 -= two_gamma2;
    }
    if r_plus - r0 == Q - 1 {
        // Boundary case: collapse r1 to 0, decrement r0.
        (0, r0 - 1)
    } else {
        ((r_plus - r0) / two_gamma2, r0)
    }
}

fn high_bits(r: i32) -> i32 {
    decompose(r).0
}

fn low_bits(r: i32) -> i32 {
    decompose(r).1
}

/// MakeHint per FIPS 204 Algorithm 39.  Returns 1 if the high bits of
/// `z + r` differ from those of `r`.
fn make_hint(z: i32, r: i32) -> u8 {
    let r1 = high_bits(r);
    let v1 = high_bits((r + z).rem_euclid(Q));
    if r1 != v1 {
        1
    } else {
        0
    }
}

/// UseHint per FIPS 204 Algorithm 40.
fn use_hint(h: u8, r: i32) -> i32 {
    let (r1, r0) = decompose(r);
    if h == 0 {
        return r1;
    }
    if r0 > 0 {
        (r1 + 1).rem_euclid(16)
    } else {
        (r1 - 1).rem_euclid(16)
    }
}

// ── Vector / matrix types ────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct PolyVecL([Poly; L]);
#[derive(Clone, Debug)]
struct PolyVecK([Poly; K]);

impl PolyVecL {
    fn zero() -> Self {
        PolyVecL(core::array::from_fn(|_| Poly::zero()))
    }
}
impl PolyVecK {
    fn zero() -> Self {
        PolyVecK(core::array::from_fn(|_| Poly::zero()))
    }
}

fn polyvecl_ntt(v: &mut PolyVecL) {
    for p in v.0.iter_mut() {
        ntt(p);
    }
}
fn polyveck_ntt(v: &mut PolyVecK) {
    for p in v.0.iter_mut() {
        ntt(p);
    }
}
fn polyveck_inv_ntt(v: &mut PolyVecK) {
    for p in v.0.iter_mut() {
        inv_ntt(p);
    }
}
fn polyvecl_inv_ntt(v: &mut PolyVecL) {
    for p in v.0.iter_mut() {
        inv_ntt(p);
    }
}

fn polyveck_add(a: &PolyVecK, b: &PolyVecK) -> PolyVecK {
    PolyVecK(core::array::from_fn(|i| poly_add(&a.0[i], &b.0[i])))
}
fn polyveck_sub(a: &PolyVecK, b: &PolyVecK) -> PolyVecK {
    PolyVecK(core::array::from_fn(|i| poly_sub(&a.0[i], &b.0[i])))
}
fn polyvecl_add(a: &PolyVecL, b: &PolyVecL) -> PolyVecL {
    PolyVecL(core::array::from_fn(|i| poly_add(&a.0[i], &b.0[i])))
}

/// Matrix `A` has shape K×L (K rows, L columns).  Stored as `[K][L]`.
type Matrix = [[Poly; L]; K];

fn matrix_vector_mul(a: &Matrix, v: &PolyVecL) -> PolyVecK {
    let mut out = PolyVecK::zero();
    for i in 0..K {
        let mut acc = Poly::zero();
        for j in 0..L {
            let prod = poly_pointwise(&a[i][j], &v.0[j]);
            acc = poly_add(&acc, &prod);
        }
        out.0[i] = acc;
    }
    out
}

// ── Sampling ─────────────────────────────────────────────────────────────────

/// Rejection sampling from SHAKE128 stream into uniform [0,q) coefficients.
/// Used by ExpandA to build matrix entries.  FIPS 204 §7.3 Algorithm 30
/// `RejNTTPoly` — samples coefficients directly in the NTT domain (each
/// coefficient is independent so this is fine).
fn rej_uniform(seed: &[u8], nonce: u16) -> Poly {
    // SHAKE128(seed || nonce_le16) — output stream long enough.
    let mut input = Vec::with_capacity(seed.len() + 2);
    input.extend_from_slice(seed);
    input.push((nonce & 0xff) as u8);
    input.push((nonce >> 8) as u8);

    let mut len = 5 * 168; // plenty for ~256 coeffs even with ~1/4 rejection
    let mut p = Poly::zero();
    loop {
        let buf = shake128(&input, len);
        let mut ctr = 0usize;
        let mut pos = 0usize;
        while pos + 3 <= buf.len() && ctr < N {
            let b0 = buf[pos] as u32;
            let b1 = buf[pos + 1] as u32;
            let b2 = (buf[pos + 2] as u32) & 0x7f;
            let t = b0 | (b1 << 8) | (b2 << 16);
            pos += 3;
            if (t as i32) < Q {
                p.0[ctr] = t as i32;
                ctr += 1;
            }
        }
        if ctr == N {
            return p;
        }
        len *= 2;
    }
}

/// Rejection sampling from SHAKE256 stream into [-η, η] (η=4) coefficients.
/// FIPS 204 Algorithm 31 `RejBoundedPoly`.
fn rej_eta(seed: &[u8], nonce: u16) -> Poly {
    let mut input = Vec::with_capacity(seed.len() + 2);
    input.extend_from_slice(seed);
    input.push((nonce & 0xff) as u8);
    input.push((nonce >> 8) as u8);

    let mut len = 3 * 136;
    let mut p = Poly::zero();
    loop {
        let buf = shake256(&input, len);
        let mut ctr = 0usize;
        let mut pos = 0usize;
        while pos < buf.len() && ctr < N {
            let byte = buf[pos];
            pos += 1;
            let nib_lo = byte & 0x0f;
            let nib_hi = byte >> 4;
            // For η=4 we accept nibbles in [0,8] and map to [η, -η].
            if nib_lo < 9 && ctr < N {
                p.0[ctr] = ETA - (nib_lo as i32);
                ctr += 1;
            }
            if nib_hi < 9 && ctr < N {
                p.0[ctr] = ETA - (nib_hi as i32);
                ctr += 1;
            }
        }
        if ctr == N {
            return p;
        }
        len *= 2;
    }
}

/// ExpandA: derive the K×L matrix A from the public seed `rho`.
/// Each `A[i][j]` is sampled via `rej_uniform(rho, (i<<8)|j)`.
fn expand_a(rho: &[u8; SEED_BYTES]) -> Matrix {
    let mut a: Matrix = core::array::from_fn(|_| core::array::from_fn(|_| Poly::zero()));
    for i in 0..K {
        for j in 0..L {
            let nonce = ((i as u16) << 8) | (j as u16);
            a[i][j] = rej_uniform(rho, nonce);
        }
    }
    a
}

/// ExpandS: derive (s1, s2) with coefficients in [-η, η] from a seed.
/// Per FIPS 204, s1 has L polys with nonces 0..L-1, s2 has K polys with
/// nonces L..L+K-1.
fn expand_s(rho_prime: &[u8]) -> (PolyVecL, PolyVecK) {
    let mut s1 = PolyVecL::zero();
    let mut s2 = PolyVecK::zero();
    for i in 0..L {
        s1.0[i] = rej_eta(rho_prime, i as u16);
    }
    for i in 0..K {
        s2.0[i] = rej_eta(rho_prime, (L + i) as u16);
    }
    (s1, s2)
}

/// Sample a polynomial with coefficients in (-γ₁, γ₁].  Used for the mask y.
/// FIPS 204 Algorithm 34 `ExpandMask` per-polynomial step.
fn expand_mask_poly(seed: &[u8], nonce: u16) -> Poly {
    // γ₁ = 2^19, so coefficients fit in 20 bits.  Each coeff uses 20 bits =
    // 2.5 bytes; we read in 5-byte groups producing 2 coefficients.
    let mut input = Vec::with_capacity(seed.len() + 2);
    input.extend_from_slice(seed);
    input.push((nonce & 0xff) as u8);
    input.push((nonce >> 8) as u8);
    // Each polynomial needs ceil(256 * 20 / 8) = 640 bytes of stream.
    let buf = shake256(&input, POLYZ_PACKED);
    let mut p = Poly::zero();
    for i in 0..(N / 2) {
        // Two 20-bit values per 5-byte group.
        let b0 = buf[5 * i] as u32;
        let b1 = buf[5 * i + 1] as u32;
        let b2 = buf[5 * i + 2] as u32;
        let b3 = buf[5 * i + 3] as u32;
        let b4 = buf[5 * i + 4] as u32;
        let a = b0 | (b1 << 8) | ((b2 & 0x0f) << 16);
        let b = (b2 >> 4) | (b3 << 4) | (b4 << 12);
        p.0[2 * i] = GAMMA1 - (a as i32);
        p.0[2 * i + 1] = GAMMA1 - (b as i32);
    }
    p
}

/// ExpandMask producing a length-L vector with nonces `kappa, kappa+1, …`.
fn expand_mask(seed: &[u8], kappa: u16) -> PolyVecL {
    let mut y = PolyVecL::zero();
    for i in 0..L {
        y.0[i] = expand_mask_poly(seed, kappa + i as u16);
    }
    y
}

/// SampleInBall: hash `c_tilde` into a sparse polynomial with `tau`
/// coefficients in {-1,+1} and the rest 0.  FIPS 204 Algorithm 29.
fn sample_in_ball(c_tilde: &[u8]) -> Poly {
    // SHAKE256 with `c_tilde` as seed; first 8 bytes are sign bits.
    let mut stream_len = 8 + 4 * TAU;
    loop {
        let buf = shake256(c_tilde, stream_len);
        let signs = u64::from_le_bytes(buf[0..8].try_into().unwrap());
        let mut c = Poly::zero();
        let mut pos = 8usize;
        let mut sign_idx = 0u32;
        let mut ok = true;
        for i in (N - TAU)..N {
            // Find a byte j ≤ i via rejection sampling.
            loop {
                if pos >= buf.len() {
                    ok = false;
                    break;
                }
                let j = buf[pos] as usize;
                pos += 1;
                if j <= i {
                    c.0[i] = c.0[j];
                    let s = ((signs >> sign_idx) & 1) as i32;
                    c.0[j] = 1 - 2 * s;
                    sign_idx += 1;
                    break;
                }
            }
            if !ok {
                break;
            }
        }
        if ok {
            return c;
        }
        stream_len *= 2;
    }
}

// ── Packing helpers ──────────────────────────────────────────────────────────

/// Pack a polynomial whose coefficients are in [0, 2^10) (i.e. t1 high bits)
/// into 320 bytes (10 bits per coefficient).
fn pack_t1(p: &Poly) -> [u8; POLYT1_PACKED] {
    let mut out = [0u8; POLYT1_PACKED];
    for i in 0..(N / 4) {
        let a0 = p.0[4 * i] as u32 & 0x3ff;
        let a1 = p.0[4 * i + 1] as u32 & 0x3ff;
        let a2 = p.0[4 * i + 2] as u32 & 0x3ff;
        let a3 = p.0[4 * i + 3] as u32 & 0x3ff;
        out[5 * i] = a0 as u8;
        out[5 * i + 1] = ((a0 >> 8) | (a1 << 2)) as u8;
        out[5 * i + 2] = ((a1 >> 6) | (a2 << 4)) as u8;
        out[5 * i + 3] = ((a2 >> 4) | (a3 << 6)) as u8;
        out[5 * i + 4] = (a3 >> 2) as u8;
    }
    out
}

fn unpack_t1(bytes: &[u8]) -> Poly {
    let mut p = Poly::zero();
    for i in 0..(N / 4) {
        let b0 = bytes[5 * i] as u32;
        let b1 = bytes[5 * i + 1] as u32;
        let b2 = bytes[5 * i + 2] as u32;
        let b3 = bytes[5 * i + 3] as u32;
        let b4 = bytes[5 * i + 4] as u32;
        p.0[4 * i] = (b0 | (b1 << 8)) as i32 & 0x3ff;
        p.0[4 * i + 1] = ((b1 >> 2) | (b2 << 6)) as i32 & 0x3ff;
        p.0[4 * i + 2] = ((b2 >> 4) | (b3 << 4)) as i32 & 0x3ff;
        p.0[4 * i + 3] = ((b3 >> 6) | (b4 << 2)) as i32 & 0x3ff;
    }
    p
}

/// Pack t0 — centered representative in (-2^{D-1}, 2^{D-1}] (D=13) — into
/// 416 bytes (13 bits/coeff).  Convention: store `(1<<(D-1)) - coeff` so the
/// stored value is in [0, 2^D).
fn pack_t0(p: &Poly) -> [u8; POLYT0_PACKED] {
    let mut out = [0u8; POLYT0_PACKED];
    let bias = 1i32 << (D - 1);
    for i in 0..(N / 8) {
        let t: [u32; 8] = core::array::from_fn(|j| ((bias - p.0[8 * i + j]) as u32) & 0x1fff);
        // 8 × 13 bits = 104 bits = 13 bytes
        out[13 * i] = t[0] as u8;
        out[13 * i + 1] = ((t[0] >> 8) | (t[1] << 5)) as u8;
        out[13 * i + 2] = (t[1] >> 3) as u8;
        out[13 * i + 3] = ((t[1] >> 11) | (t[2] << 2)) as u8;
        out[13 * i + 4] = ((t[2] >> 6) | (t[3] << 7)) as u8;
        out[13 * i + 5] = (t[3] >> 1) as u8;
        out[13 * i + 6] = ((t[3] >> 9) | (t[4] << 4)) as u8;
        out[13 * i + 7] = (t[4] >> 4) as u8;
        out[13 * i + 8] = ((t[4] >> 12) | (t[5] << 1)) as u8;
        out[13 * i + 9] = ((t[5] >> 7) | (t[6] << 6)) as u8;
        out[13 * i + 10] = (t[6] >> 2) as u8;
        out[13 * i + 11] = ((t[6] >> 10) | (t[7] << 3)) as u8;
        out[13 * i + 12] = (t[7] >> 5) as u8;
    }
    out
}

fn unpack_t0(bytes: &[u8]) -> Poly {
    let mut p = Poly::zero();
    let bias = 1i32 << (D - 1);
    for i in 0..(N / 8) {
        let b: [u32; 13] = core::array::from_fn(|j| bytes[13 * i + j] as u32);
        let mut t = [0u32; 8];
        t[0] = (b[0] | (b[1] << 8)) & 0x1fff;
        t[1] = ((b[1] >> 5) | (b[2] << 3) | (b[3] << 11)) & 0x1fff;
        t[2] = ((b[3] >> 2) | (b[4] << 6)) & 0x1fff;
        t[3] = ((b[4] >> 7) | (b[5] << 1) | (b[6] << 9)) & 0x1fff;
        t[4] = ((b[6] >> 4) | (b[7] << 4) | (b[8] << 12)) & 0x1fff;
        t[5] = ((b[8] >> 1) | (b[9] << 7)) & 0x1fff;
        t[6] = ((b[9] >> 6) | (b[10] << 2) | (b[11] << 10)) & 0x1fff;
        t[7] = ((b[11] >> 3) | (b[12] << 5)) & 0x1fff;
        for j in 0..8 {
            p.0[8 * i + j] = bias - t[j] as i32;
        }
    }
    p
}

/// Pack a polynomial with coefficients in [-η, η] (η=4) into 128 bytes
/// (4 bits/coeff).
fn pack_eta(p: &Poly) -> [u8; POLYETA_PACKED] {
    let mut out = [0u8; POLYETA_PACKED];
    for i in 0..(N / 2) {
        let a = (ETA - p.0[2 * i]) as u32 & 0x0f;
        let b = (ETA - p.0[2 * i + 1]) as u32 & 0x0f;
        out[i] = (a | (b << 4)) as u8;
    }
    out
}

fn unpack_eta(bytes: &[u8]) -> Poly {
    let mut p = Poly::zero();
    for i in 0..(N / 2) {
        let a = bytes[i] & 0x0f;
        let b = bytes[i] >> 4;
        p.0[2 * i] = ETA - (a as i32);
        p.0[2 * i + 1] = ETA - (b as i32);
    }
    p
}

/// Pack a polynomial with coefficients in (-γ₁, γ₁] into 640 bytes
/// (20 bits/coeff).  Accepts coefficients in either canonical [0, q) form or
/// centered form — `to_centered` normalises them before packing.
fn pack_z(p: &Poly) -> [u8; POLYZ_PACKED] {
    let mut out = [0u8; POLYZ_PACKED];
    for i in 0..(N / 2) {
        let a = (GAMMA1 - to_centered(p.0[2 * i])) as u32 & 0xfffff;
        let b = (GAMMA1 - to_centered(p.0[2 * i + 1])) as u32 & 0xfffff;
        out[5 * i] = a as u8;
        out[5 * i + 1] = (a >> 8) as u8;
        out[5 * i + 2] = ((a >> 16) | (b << 4)) as u8;
        out[5 * i + 3] = (b >> 4) as u8;
        out[5 * i + 4] = (b >> 12) as u8;
    }
    out
}

fn unpack_z(bytes: &[u8]) -> Poly {
    let mut p = Poly::zero();
    for i in 0..(N / 2) {
        let b0 = bytes[5 * i] as u32;
        let b1 = bytes[5 * i + 1] as u32;
        let b2 = bytes[5 * i + 2] as u32;
        let b3 = bytes[5 * i + 3] as u32;
        let b4 = bytes[5 * i + 4] as u32;
        let a = b0 | (b1 << 8) | ((b2 & 0x0f) << 16);
        let b = (b2 >> 4) | (b3 << 4) | (b4 << 12);
        p.0[2 * i] = GAMMA1 - (a as i32);
        p.0[2 * i + 1] = GAMMA1 - (b as i32);
    }
    p
}

/// Pack w1 with 4 bits/coeff (γ₂ small case: w1 ∈ {0,…,15}).
fn pack_w1(p: &Poly) -> [u8; POLYW1_PACKED] {
    let mut out = [0u8; POLYW1_PACKED];
    for i in 0..(N / 2) {
        out[i] = ((p.0[2 * i] & 0x0f) | ((p.0[2 * i + 1] & 0x0f) << 4)) as u8;
    }
    out
}

// ── Public types ─────────────────────────────────────────────────────────────

/// ML-DSA public key — opaque byte blob conforming to FIPS 204 §7.2 pkEncode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MlDsaPublicKey(pub Vec<u8>);

/// ML-DSA secret key — opaque byte blob conforming to FIPS 204 §7.2 skEncode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MlDsaSecretKey(pub Vec<u8>);

impl Drop for MlDsaSecretKey {
    fn drop(&mut self) {
        for b in self.0.iter_mut() {
            *b = 0;
        }
    }
}

// ── Key generation ───────────────────────────────────────────────────────────

/// ML-DSA-65 key generation from a 32-byte seed.  Deterministic in `seed`.
pub fn ml_dsa_65_keygen(seed: &[u8; SEED_BYTES]) -> (MlDsaPublicKey, MlDsaSecretKey) {
    // FIPS 204 §6 Algorithm 1.
    // Expand seed: SHAKE256(seed || k_byte || l_byte, 128) → rho || rho' || K
    let mut shake_input = Vec::with_capacity(SEED_BYTES + 2);
    shake_input.extend_from_slice(seed);
    shake_input.push(K as u8);
    shake_input.push(L as u8);
    let h = shake256(&shake_input, 128);

    let mut rho = [0u8; SEED_BYTES];
    rho.copy_from_slice(&h[0..32]);
    let mut rho_prime = [0u8; 64];
    rho_prime.copy_from_slice(&h[32..96]);
    let mut big_k = [0u8; SEED_BYTES];
    big_k.copy_from_slice(&h[96..128]);

    // A in NTT domain
    let a_hat = expand_a(&rho);
    let (s1, s2) = expand_s(&rho_prime);

    // t = A·s1 + s2  (computed in NTT domain, then inverse NTT)
    let mut s1_hat = s1.clone();
    polyvecl_ntt(&mut s1_hat);
    let mut t = matrix_vector_mul(&a_hat, &s1_hat);
    polyveck_inv_ntt(&mut t);
    let t = polyveck_add(&t, &s2);

    // Decompose t into (t1, t0): t1 holds the high D=13 bits.
    let mut t1 = PolyVecK::zero();
    let mut t0 = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            let (h1, h0) = power2_round(t.0[i].0[j]);
            t1.0[i].0[j] = h1;
            t0.0[i].0[j] = h0;
        }
    }

    // Pack public key: rho || t1
    let mut pk = Vec::with_capacity(PUBKEY_BYTES);
    pk.extend_from_slice(&rho);
    for i in 0..K {
        pk.extend_from_slice(&pack_t1(&t1.0[i]));
    }

    // tr = SHAKE256(pk, 64)
    let tr = shake256(&pk, TR_BYTES);

    // Pack secret key: rho || K || tr || s1 || s2 || t0
    let mut sk = Vec::with_capacity(SECKEY_BYTES);
    sk.extend_from_slice(&rho);
    sk.extend_from_slice(&big_k);
    sk.extend_from_slice(&tr);
    for i in 0..L {
        sk.extend_from_slice(&pack_eta(&s1.0[i]));
    }
    for i in 0..K {
        sk.extend_from_slice(&pack_eta(&s2.0[i]));
    }
    for i in 0..K {
        sk.extend_from_slice(&pack_t0(&t0.0[i]));
    }

    (MlDsaPublicKey(pk), MlDsaSecretKey(sk))
}

// ── Secret-key unpacking ─────────────────────────────────────────────────────

struct UnpackedSk {
    rho: [u8; SEED_BYTES],
    big_k: [u8; SEED_BYTES],
    tr: [u8; TR_BYTES],
    s1: PolyVecL,
    s2: PolyVecK,
    t0: PolyVecK,
}

fn unpack_sk(sk: &[u8]) -> UnpackedSk {
    let mut rho = [0u8; SEED_BYTES];
    rho.copy_from_slice(&sk[0..SEED_BYTES]);
    let mut big_k = [0u8; SEED_BYTES];
    big_k.copy_from_slice(&sk[SEED_BYTES..2 * SEED_BYTES]);
    let mut tr = [0u8; TR_BYTES];
    tr.copy_from_slice(&sk[2 * SEED_BYTES..2 * SEED_BYTES + TR_BYTES]);

    let mut off = 2 * SEED_BYTES + TR_BYTES;
    let mut s1 = PolyVecL::zero();
    for i in 0..L {
        s1.0[i] = unpack_eta(&sk[off..off + POLYETA_PACKED]);
        off += POLYETA_PACKED;
    }
    let mut s2 = PolyVecK::zero();
    for i in 0..K {
        s2.0[i] = unpack_eta(&sk[off..off + POLYETA_PACKED]);
        off += POLYETA_PACKED;
    }
    let mut t0 = PolyVecK::zero();
    for i in 0..K {
        t0.0[i] = unpack_t0(&sk[off..off + POLYT0_PACKED]);
        off += POLYT0_PACKED;
    }
    UnpackedSk {
        rho,
        big_k,
        tr,
        s1,
        s2,
        t0,
    }
}

// ── Signing ──────────────────────────────────────────────────────────────────

/// ML-DSA-65 signing (FIPS 204 §6 Algorithm 2).
///
/// `rnd` is the "hedge" randomness: for deterministic signatures pass
/// `[0u8; 32]`, for randomized/hedged signatures pass 32 fresh random bytes.
pub fn ml_dsa_65_sign(sk: &MlDsaSecretKey, msg: &[u8], rnd: &[u8; SEED_BYTES]) -> Vec<u8> {
    let unp = unpack_sk(&sk.0);

    // mu = SHAKE256(tr || msg, 64)
    let mut mu_in = Vec::with_capacity(TR_BYTES + msg.len());
    mu_in.extend_from_slice(&unp.tr);
    mu_in.extend_from_slice(msg);
    let mu = shake256(&mu_in, 64);

    // rho'' = SHAKE256(K || rnd || mu, 64)
    let mut rho2_in = Vec::with_capacity(SEED_BYTES + SEED_BYTES + 64);
    rho2_in.extend_from_slice(&unp.big_k);
    rho2_in.extend_from_slice(rnd);
    rho2_in.extend_from_slice(&mu);
    let rho_prime_prime = shake256(&rho2_in, 64);

    // Re-derive A in NTT form, and NTT-transform s1, s2, t0.
    let a_hat = expand_a(&unp.rho);
    let mut s1_hat = unp.s1.clone();
    polyvecl_ntt(&mut s1_hat);
    let mut s2_hat = unp.s2.clone();
    polyveck_ntt(&mut s2_hat);
    let mut t0_hat = unp.t0.clone();
    polyveck_ntt(&mut t0_hat);

    let mut kappa: u16 = 0;
    loop {
        // Sample mask y with coefficients in (-γ₁, γ₁].
        let y = expand_mask(&rho_prime_prime, kappa);
        kappa += L as u16;

        // w = A · NTT(y), then inverse NTT.
        let mut y_hat = y.clone();
        polyvecl_ntt(&mut y_hat);
        let mut w = matrix_vector_mul(&a_hat, &y_hat);
        polyveck_inv_ntt(&mut w);

        // w1 = HighBits(w)
        let mut w1 = PolyVecK::zero();
        for i in 0..K {
            for j in 0..N {
                w1.0[i].0[j] = high_bits(w.0[i].0[j]);
            }
        }

        // c_tilde = SHAKE256(mu || w1_packed, λ/4)
        let mut ch_in = Vec::with_capacity(64 + K * POLYW1_PACKED);
        ch_in.extend_from_slice(&mu);
        for i in 0..K {
            ch_in.extend_from_slice(&pack_w1(&w1.0[i]));
        }
        let c_tilde = shake256(&ch_in, C_TILDE_BYTES);

        // c ← SampleInBall(c_tilde); compute c_hat in NTT form.
        let c = sample_in_ball(&c_tilde);
        let mut c_hat = c.clone();
        ntt(&mut c_hat);

        // z = y + c·s1
        let mut cs1 = PolyVecL::zero();
        for i in 0..L {
            cs1.0[i] = poly_pointwise(&c_hat, &s1_hat.0[i]);
        }
        polyvecl_inv_ntt(&mut cs1);
        let z = polyvecl_add(&y, &cs1);

        // Check ‖z‖_∞ < γ₁ - β
        let mut z_ok = true;
        for i in 0..L {
            if poly_inf_norm(&z.0[i]) >= GAMMA1 - BETA {
                z_ok = false;
                break;
            }
        }
        if !z_ok {
            continue;
        }

        // r0 = LowBits(w - c·s2)
        let mut cs2 = PolyVecK::zero();
        for i in 0..K {
            cs2.0[i] = poly_pointwise(&c_hat, &s2_hat.0[i]);
        }
        polyveck_inv_ntt(&mut cs2);
        let w_minus_cs2 = polyveck_sub(&w, &cs2);

        let mut r0_ok = true;
        for i in 0..K {
            for j in 0..N {
                let r0 = low_bits(w_minus_cs2.0[i].0[j]);
                if to_centered(r0).abs() >= GAMMA2 - BETA {
                    r0_ok = false;
                    break;
                }
            }
            if !r0_ok {
                break;
            }
        }
        if !r0_ok {
            continue;
        }

        // c·t0
        let mut ct0 = PolyVecK::zero();
        for i in 0..K {
            ct0.0[i] = poly_pointwise(&c_hat, &t0_hat.0[i]);
        }
        polyveck_inv_ntt(&mut ct0);

        // Check ‖c·t0‖_∞ < γ₂
        let mut ct0_ok = true;
        for i in 0..K {
            if poly_inf_norm(&ct0.0[i]) >= GAMMA2 {
                ct0_ok = false;
                break;
            }
        }
        if !ct0_ok {
            continue;
        }

        // h = MakeHint(-c·t0, w - c·s2 + c·t0)
        let arg = polyveck_add(&w_minus_cs2, &ct0);
        let mut hints = vec![[0u8; N]; K];
        let mut total_hints = 0usize;
        for i in 0..K {
            for j in 0..N {
                let neg_ct0 = (-ct0.0[i].0[j]).rem_euclid(Q);
                let h = make_hint(neg_ct0, arg.0[i].0[j].rem_euclid(Q));
                hints[i][j] = h;
                total_hints += h as usize;
            }
        }
        if total_hints > OMEGA {
            continue;
        }

        // Encode signature: c_tilde || z_packed || hints_packed
        let mut sig = Vec::with_capacity(SIG_BYTES);
        sig.extend_from_slice(&c_tilde);
        for i in 0..L {
            sig.extend_from_slice(&pack_z(&z.0[i]));
        }
        // Pack hints: list of OMEGA byte positions followed by K running counts.
        let mut hint_bytes = vec![0u8; OMEGA + K];
        let mut k_ptr = 0usize;
        for i in 0..K {
            for j in 0..N {
                if hints[i][j] == 1 {
                    hint_bytes[k_ptr] = j as u8;
                    k_ptr += 1;
                }
            }
            hint_bytes[OMEGA + i] = k_ptr as u8;
        }
        sig.extend_from_slice(&hint_bytes);
        return sig;
    }
}

// ── Verification ─────────────────────────────────────────────────────────────

pub fn ml_dsa_65_verify(pk: &MlDsaPublicKey, msg: &[u8], sig: &[u8]) -> bool {
    if pk.0.len() != PUBKEY_BYTES || sig.len() != SIG_BYTES {
        return false;
    }
    // Unpack public key: rho || t1
    let mut rho = [0u8; SEED_BYTES];
    rho.copy_from_slice(&pk.0[0..SEED_BYTES]);
    let mut t1 = PolyVecK::zero();
    for i in 0..K {
        let off = SEED_BYTES + i * POLYT1_PACKED;
        t1.0[i] = unpack_t1(&pk.0[off..off + POLYT1_PACKED]);
    }

    // Unpack signature: c_tilde || z || hints
    let c_tilde = &sig[0..C_TILDE_BYTES];
    let mut z = PolyVecL::zero();
    for i in 0..L {
        let off = C_TILDE_BYTES + i * POLYZ_PACKED;
        z.0[i] = unpack_z(&sig[off..off + POLYZ_PACKED]);
    }

    // Hint-packed format: positions (OMEGA bytes), then K running totals.
    let hint_off = C_TILDE_BYTES + L * POLYZ_PACKED;
    let hint_bytes = &sig[hint_off..hint_off + OMEGA + K];
    let mut hints = vec![[0u8; N]; K];
    let mut k_ptr = 0usize;
    for i in 0..K {
        let end_idx = hint_bytes[OMEGA + i] as usize;
        if end_idx < k_ptr || end_idx > OMEGA {
            return false;
        }
        for j in k_ptr..end_idx {
            // positions must be strictly increasing within a polynomial
            if j > k_ptr && hint_bytes[j] <= hint_bytes[j - 1] {
                return false;
            }
            hints[i][hint_bytes[j] as usize] = 1;
        }
        k_ptr = end_idx;
    }
    // Trailing bytes in the omega region must be zero (no extra positions).
    for j in k_ptr..OMEGA {
        if hint_bytes[j] != 0 {
            return false;
        }
    }

    // ‖z‖_∞ < γ₁ - β
    for i in 0..L {
        if poly_inf_norm(&z.0[i]) >= GAMMA1 - BETA {
            return false;
        }
    }

    // tr = SHAKE256(pk, 64); mu = SHAKE256(tr || msg, 64)
    let tr = shake256(&pk.0, TR_BYTES);
    let mut mu_in = Vec::with_capacity(TR_BYTES + msg.len());
    mu_in.extend_from_slice(&tr);
    mu_in.extend_from_slice(msg);
    let mu = shake256(&mu_in, 64);

    // c ← SampleInBall(c_tilde); c_hat
    let c = sample_in_ball(c_tilde);
    let mut c_hat = c.clone();
    ntt(&mut c_hat);

    // Compute A·z - c·t1·2^D (in NTT, then inverse).
    let a_hat = expand_a(&rho);
    let mut z_hat = z.clone();
    polyvecl_ntt(&mut z_hat);
    let mut az = matrix_vector_mul(&a_hat, &z_hat);

    // Shift t1 by 2^D, transform, multiply by c, subtract.
    let mut t1_shift = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            t1_shift.0[i].0[j] = (t1.0[i].0[j] << D).rem_euclid(Q);
        }
        ntt(&mut t1_shift.0[i]);
    }
    for i in 0..K {
        let ct1 = poly_pointwise(&c_hat, &t1_shift.0[i]);
        az.0[i] = poly_sub(&az.0[i], &ct1);
    }
    polyveck_inv_ntt(&mut az);

    // w1' = UseHint(h, A·z - c·t1·2^D)
    let mut w1 = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            w1.0[i].0[j] = use_hint(hints[i][j], az.0[i].0[j]);
        }
    }

    // c_tilde' = SHAKE256(mu || w1_packed, λ/4)
    let mut ch_in = Vec::with_capacity(64 + K * POLYW1_PACKED);
    ch_in.extend_from_slice(&mu);
    for i in 0..K {
        ch_in.extend_from_slice(&pack_w1(&w1.0[i]));
    }
    let c_tilde_prime = shake256(&ch_in, C_TILDE_BYTES);

    c_tilde_prime.as_slice() == c_tilde
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn fixed_seed(byte: u8) -> [u8; SEED_BYTES] {
        [byte; SEED_BYTES]
    }

    #[test]
    fn ntt_roundtrip_identity() {
        // For random-ish input, inv_ntt(ntt(p)) == p.
        let mut p = Poly::zero();
        for i in 0..N {
            p.0[i] = ((i as i64 * 12345 + 7) as i32).rem_euclid(Q);
        }
        let original = p.clone();
        ntt(&mut p);
        inv_ntt(&mut p);
        assert_eq!(p, original);
    }

    #[test]
    fn ntt_multiplication_matches_schoolbook() {
        // For two small polynomials, NTT-domain pointwise multiply then
        // inverse NTT should match the schoolbook product in R_q.
        let mut a = Poly::zero();
        let mut b = Poly::zero();
        for i in 0..16 {
            a.0[i] = (i as i32 * 7 + 3).rem_euclid(Q);
            b.0[i] = (i as i32 * 5 + 11).rem_euclid(Q);
        }
        // Schoolbook reference
        let mut expected = Poly::zero();
        for i in 0..N {
            for j in 0..N {
                let idx = (i + j) % N;
                let sign: i64 = if i + j >= N { -1 } else { 1 };
                let term = sign * (a.0[i] as i64) * (b.0[j] as i64);
                expected.0[idx] = ((expected.0[idx] as i64 + term).rem_euclid(Q as i64)) as i32;
            }
        }
        // NTT path
        let mut ah = a.clone();
        let mut bh = b.clone();
        ntt(&mut ah);
        ntt(&mut bh);
        let mut prod = poly_pointwise(&ah, &bh);
        inv_ntt(&mut prod);
        assert_eq!(prod, expected);
    }

    #[test]
    fn power2_round_reconstructs() {
        for r in [0i32, 1, 12345, Q - 1, (1 << D) - 1, 1 << D, Q / 2] {
            let (r1, r0) = power2_round(r);
            let recon = (r1 * (1 << D) + r0).rem_euclid(Q);
            assert_eq!(recon, r, "power2_round failed for r={}", r);
            assert!(r0.abs() <= (1 << (D - 1)));
        }
    }

    #[test]
    fn decompose_reconstructs() {
        // r = r1 * 2γ₂ + r0 (mod q), with r0 in (-γ₂, γ₂].
        for r in [0i32, 1, 100_000, Q - 1, Q / 2, 2 * GAMMA2 - 1, 2 * GAMMA2] {
            let (r1, r0) = decompose(r);
            let recon = (r1 * 2 * GAMMA2 + r0).rem_euclid(Q);
            assert_eq!(recon, r, "decompose failed for r={}", r);
            assert!(r0 > -GAMMA2 && r0 <= GAMMA2);
            assert!((0..16).contains(&r1));
        }
    }

    #[test]
    fn keygen_is_deterministic() {
        let seed = fixed_seed(0x42);
        let (pk1, sk1) = ml_dsa_65_keygen(&seed);
        let (pk2, sk2) = ml_dsa_65_keygen(&seed);
        assert_eq!(pk1.0, pk2.0);
        assert_eq!(sk1.0, sk2.0);
        assert_eq!(pk1.0.len(), PUBKEY_BYTES);
        assert_eq!(sk1.0.len(), SECKEY_BYTES);
    }

    #[test]
    fn different_seeds_produce_different_keys() {
        let (pk1, _) = ml_dsa_65_keygen(&fixed_seed(0));
        let (pk2, _) = ml_dsa_65_keygen(&fixed_seed(1));
        assert_ne!(pk1.0, pk2.0);
    }

    #[test]
    fn sign_verify_roundtrip_deterministic() {
        let seed = fixed_seed(7);
        let (pk, sk) = ml_dsa_65_keygen(&seed);
        let msg = b"hello, ML-DSA-65";
        let rnd = [0u8; SEED_BYTES];
        let sig = ml_dsa_65_sign(&sk, msg, &rnd);
        assert_eq!(sig.len(), SIG_BYTES);
        assert!(ml_dsa_65_verify(&pk, msg, &sig));
    }

    #[test]
    fn sign_verify_roundtrip_hedged() {
        let seed = fixed_seed(9);
        let (pk, sk) = ml_dsa_65_keygen(&seed);
        let msg = b"hedged signature payload";
        let rnd1 = [0xa5u8; SEED_BYTES];
        let rnd2 = [0x5au8; SEED_BYTES];
        let sig1 = ml_dsa_65_sign(&sk, msg, &rnd1);
        let sig2 = ml_dsa_65_sign(&sk, msg, &rnd2);
        // Different rnd → different signatures (overwhelming probability).
        assert_ne!(sig1, sig2);
        assert!(ml_dsa_65_verify(&pk, msg, &sig1));
        assert!(ml_dsa_65_verify(&pk, msg, &sig2));
    }

    #[test]
    fn wrong_message_rejected() {
        let (pk, sk) = ml_dsa_65_keygen(&fixed_seed(3));
        let sig = ml_dsa_65_sign(&sk, b"correct", &[0u8; SEED_BYTES]);
        assert!(!ml_dsa_65_verify(&pk, b"tampered", &sig));
    }

    #[test]
    fn flipped_bit_rejected() {
        let (pk, sk) = ml_dsa_65_keygen(&fixed_seed(4));
        let msg = b"sign me";
        let mut sig = ml_dsa_65_sign(&sk, msg, &[0u8; SEED_BYTES]);
        // Flip a bit deep in the z portion (well past c_tilde).
        sig[200] ^= 0x01;
        assert!(!ml_dsa_65_verify(&pk, msg, &sig));
    }

    #[test]
    fn wrong_public_key_rejected() {
        let (pk_a, sk_a) = ml_dsa_65_keygen(&fixed_seed(11));
        let (pk_b, _sk_b) = ml_dsa_65_keygen(&fixed_seed(22));
        let msg = b"audience of one";
        let sig = ml_dsa_65_sign(&sk_a, msg, &[0u8; SEED_BYTES]);
        assert!(ml_dsa_65_verify(&pk_a, msg, &sig));
        assert!(!ml_dsa_65_verify(&pk_b, msg, &sig));
    }

    #[test]
    fn malformed_sig_length_rejected() {
        let (pk, _) = ml_dsa_65_keygen(&fixed_seed(13));
        assert!(!ml_dsa_65_verify(&pk, b"x", b""));
        assert!(!ml_dsa_65_verify(&pk, b"x", &vec![0u8; SIG_BYTES - 1]));
    }

    #[test]
    fn pubkey_size_matches_spec() {
        // FIPS 204 ML-DSA-65: pk = 1952 bytes, sk = 4032 bytes, sig = 3309 bytes.
        assert_eq!(PUBKEY_BYTES, 1952);
        assert_eq!(SECKEY_BYTES, 4032);
        assert_eq!(SIG_BYTES, 3309);
    }
}
