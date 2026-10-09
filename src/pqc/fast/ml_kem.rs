//! **ML-KEM** (FIPS 203) written for speed.
//!
//! Byte-for-byte interchangeable with [`crate::pqc::ml_kem`], which stays the
//! readable reference and is validated against the NIST ACVP vectors. The tests
//! at the bottom of this file hold the two implementations against each other
//! on every parameter set, so the ACVP validation carries over.
//!
//! # Where the time was going
//!
//! Measured on the reference implementation, ML-KEM-768 key generation cost
//! 337 kilocycles, of which hashing was 40. The other 88% was ring arithmetic,
//! and almost all of that came from two lines rather than from the modular
//! reductions one would suspect:
//!
//! * `ntt` and `ntt_inv` each began by calling `zetas()`, which rebuilds all
//!   128 twiddle factors by square-and-multiply. That is roughly 1800 modular
//!   multiplications to set up a transform that only performs 1024.
//! * `multiply_ntts` called `zeta_pow(2·BitRev7(i) + 1)` inside its loop, so a
//!   single pointwise product ran 128 modular exponentiations.
//!
//! Both are pure setup that does not depend on the input, so this module keeps
//! them in `const` tables computed at compile time.
//!
//! # What else changed
//!
//! * **Signed coefficients and Montgomery multiplication.** Coefficients are
//!   `i16` in `(-q, q)` rather than `u16` in `[0, q)`, and products go through
//!   `montgomery_reduce`, which is two multiplications and a shift. The
//!   reference reduced with `% q` on `u32`, which the compiler turns into a
//!   multiply-shift pair as well, so this is worth less than it looks; the
//!   gain is that intermediate values need reducing far less often.
//! * **Streaming rejection sampling.** `sample_ntt` needs an unpredictable
//!   number of SHAKE128 bytes. The reference asked for a fixed length and, when
//!   that ran short, re-ran the whole sponge with a doubled request, repeating
//!   every permutation it had already done. This one squeezes another block
//!   from the same state.
//! * **No allocation in the hot paths.** Polynomials and vectors are fixed
//!   arrays. The reference returned `Vec` from `byte_encode`, `prf` and the
//!   polynomial operators, so a single key generation made hundreds of
//!   short-lived heap allocations.
//!
//! # Second round: system-level work and SIMD
//!
//! After the changes above, ML-KEM-768 encapsulation measured ~94 kilocycles
//! and the profile had moved: Keccak was ~45% of instructions, expanding `Â`
//! ~42%, and hashing the encapsulation key ~10%. The second round follows
//! Sayed, Taha and Nijjer, *System-Level Optimization Beyond Cryptographic
//! Kernels* (arXiv:2610.01960), which finds the same split on a Cortex-M7 and
//! argues for looking past the arithmetic kernels at work repeated across
//! calls. In the order it was done, with what each step measured:
//!
//! * **Public-data reuse.** [`MlKemPreparedEncapsKey`] and
//!   [`MlKemPreparedDecapsKey`] keep `Â`, `t̂`, `H(ek)` and `ŝ`, all derived
//!   from public data or (for `ŝ`) from the key itself. Key generation
//!   computes all of it anyway, so [`ml_kem_keygen_prepared_internal`]
//!   captures it for the price of copies. Wire formats are untouched.
//! * **Four-way Keccak** ([`super::keccak4`]): the k² matrix entries and the
//!   noise PRFs are independent sponges, run four to a 256-bit register.
//! * **AVX2 ring arithmetic** ([`avx2`]): NTT, inverse NTT, pointwise product,
//!   Barrett and Montgomery steps, bit-identical to the scalar code.
//! * **Grouped packing, vector compression, vector rejection sampling**, and
//!   single sponges (`H`, `G`, `J`) through one lane of the AVX-512 Keccak.
//!
//! Every SIMD path is chosen at run time and has the scalar code as its
//! fallback; `docs/pqc-speed.md` has the per-step table and the per-class
//! (AVX-512, AVX2-only, scalar) measurements.
//!
//! # Domains
//!
//! The convention is the one the pq-crystals reference uses, because it is the
//! one whose constant factors are known to work out:
//!
//! * twiddle factors are stored in Montgomery form, so `fqmul(x, zeta)` is the
//!   true product `x·ζ` rather than `x·ζ/2^16`;
//! * `ntt` therefore computes the true transform;
//! * `basemul` multiplies two ordinary values and so produces the true product
//!   divided by `2^16`;
//! * `invntt` folds in `2^16/128`, so it produces its result multiplied by
//!   `2^16`, which exactly cancels what `basemul` left behind.
//!
//! The one place this does not cancel is `t̂ = Â∘ŝ + ê`, where the matrix
//! product stays in the NTT domain and never meets an `invntt`. There the
//! product is scaled back up with [`poly_tomont`] before `ê` is added. Getting
//! this wrong produces output that is wrong by a constant factor, which the
//! differential tests catch immediately.
//!
//! # Not constant-time
//!
//! As with the rest of this library: the rejection samplers branch on their
//! output and the compressions are not written to avoid data-dependent timing.
//! Do not use this for anything real. See `SECURITY.md`.

use super::keccak::{sha3_256, sha3_512_2, shake256_2_into, Sponge, SHAKE128_RATE};
use super::keccak4;
use crate::pqc::ml_kem::{MlKemDecapsKey, MlKemEncapsKey, MlKemParams, SHARED_SECRET_BYTES};
use crate::utils::random::random_bytes;
use subtle::{ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroize;

// ── ring parameters ──────────────────────────────────────────────────────────

/// Polynomial degree: `R_q = Z_q[X]/(X^256 + 1)`.
const N: usize = 256;
/// Coefficient modulus.
const Q: i16 = 3329;
const QI: i32 = Q as i32;
/// `q^-1 mod 2^16`, as the `i16` it is reinterpreted as (62209 − 65536).
const QINV: i16 = -3327;
/// `ζ = 17` generates the 256th roots of unity mod q.
const ZETA: i32 = 17;
/// The largest module rank over the three parameter sets, so that vectors can
/// be fixed-size arrays instead of `Vec`.
const KMAX: usize = 4;

// ── Montgomery and Barrett reduction ─────────────────────────────────────────

/// `a · 2^-16 mod q`, as a value in `(-q, q)`.
///
/// Requires `|a| < q · 2^15`, which every caller here satisfies because both
/// factors of the product are bounded by `q` in absolute value.
#[inline(always)]
const fn montgomery_reduce(a: i32) -> i16 {
    let t = (a as i16).wrapping_mul(QINV);
    ((a - (t as i32) * QI) >> 16) as i16
}

/// `a mod q` as a value in `[-(q-1)/2, (q-1)/2]`, for any `i16`.
#[inline(always)]
const fn barrett_reduce(a: i16) -> i16 {
    // v = round(2^26 / q); the +2^25 makes the shift round to nearest.
    const V: i32 = ((1i32 << 26) / QI) + 1;
    let t = ((V * a as i32 + (1 << 25)) >> 26) as i16;
    a - t.wrapping_mul(Q)
}

/// The true product of two ordinary values, divided by `2^16`. With a
/// Montgomery-form second argument this is the true product.
#[inline(always)]
const fn fqmul(a: i16, b: i16) -> i16 {
    montgomery_reduce(a as i32 * b as i32)
}

/// Map a value in `(-q, q)` to its representative in `[0, q)`.
#[inline(always)]
const fn to_canonical(a: i16) -> u16 {
    (a + ((a >> 15) & Q)) as u16
}

// ── compile-time twiddle tables ──────────────────────────────────────────────

/// `base^exp mod q`, for table construction only.
const fn pow_mod(base: i32, exp: u32) -> i32 {
    let mut r: i64 = 1;
    let mut b: i64 = (base % QI) as i64;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % QI as i64;
        }
        b = b * b % QI as i64;
        e >>= 1;
    }
    r as i32
}

/// `v · 2^16 mod q`, as a signed value in `[-(q-1)/2, (q-1)/2]`.
const fn to_mont_signed(v: i32) -> i16 {
    let m = ((v as i64) * 65536 % QI as i64) as i32;
    let m = if m > QI / 2 { m - QI } else { m };
    m as i16
}

/// Reverse the low seven bits (FIPS 203 §2.3 `BitRev7`).
const fn bitrev7(x: usize) -> usize {
    let mut r = 0;
    let mut b = 0;
    while b < 7 {
        r |= ((x >> b) & 1) << (6 - b);
        b += 1;
    }
    r
}

/// `ζ^BitRev7(i)` in Montgomery form: the twiddles the seven NTT layers use.
const fn zetas_table() -> [i16; 128] {
    let mut z = [0i16; 128];
    let mut i = 0;
    while i < 128 {
        z[i] = to_mont_signed(pow_mod(ZETA, bitrev7(i) as u32));
        i += 1;
    }
    z
}

/// `ζ^(2·BitRev7(i)+1)` in Montgomery form: the modulus of the i-th quadratic
/// residue class, needed once per `basemul` pair.
///
/// The reference recomputed this by modular exponentiation inside the pointwise
/// multiplication loop, 128 times per product.
const fn basemul_zetas_table() -> [i16; 128] {
    let mut z = [0i16; 128];
    let mut i = 0;
    while i < 128 {
        z[i] = to_mont_signed(pow_mod(ZETA, 2 * bitrev7(i) as u32 + 1));
        i += 1;
    }
    z
}

static ZETAS: [i16; 128] = zetas_table();
static BASEMUL_ZETAS: [i16; 128] = basemul_zetas_table();

/// `2^32 mod q`: multiplying by this through `fqmul` scales by `2^16`.
const MONT2: i16 = ((1i64 << 32) % QI as i64) as i16;
/// `2^32 / 128 mod q`: the inverse NTT's final scaling, which folds the `1/128`
/// together with the move into the Montgomery domain.
const INVNTT_F: i16 = 1441;

// ── polynomials ──────────────────────────────────────────────────────────────

/// A polynomial of `R_q`, coefficients in `(-q, q)` unless stated otherwise.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Poly([i16; N]);

impl Poly {
    const fn zero() -> Self {
        Poly([0; N])
    }
}

#[inline]
fn poly_add(r: &mut Poly, a: &Poly) {
    for i in 0..N {
        r.0[i] = r.0[i].wrapping_add(a.0[i]);
    }
}

#[inline]
fn poly_sub(r: &mut Poly, a: &Poly) {
    for i in 0..N {
        r.0[i] = r.0[i].wrapping_sub(a.0[i]);
    }
}

#[inline]
fn poly_reduce_scalar(r: &mut Poly) {
    for c in r.0.iter_mut() {
        *c = barrett_reduce(*c);
    }
}

/// Scale by `2^16`, cancelling the factor `basemul` divides out.
#[inline]
fn poly_tomont_scalar(r: &mut Poly) {
    for c in r.0.iter_mut() {
        *c = fqmul(*c, MONT2);
    }
}

/// Forward NTT, in place. Input coefficients bounded by `q` in absolute value;
/// output bounded by `8q` before the final reduction, which is applied here.
fn ntt_scalar(r: &mut Poly) {
    let mut k = 1usize;
    let mut len = 128usize;
    while len >= 2 {
        let mut start = 0usize;
        while start < N {
            let zeta = ZETAS[k];
            k += 1;
            for j in start..start + len {
                let t = fqmul(zeta, r.0[j + len]);
                r.0[j + len] = r.0[j] - t;
                r.0[j] += t;
            }
            start += 2 * len;
        }
        len >>= 1;
    }
    poly_reduce_scalar(r);
}

/// Inverse NTT, in place, including the `1/128` and the move into the
/// Montgomery domain (see the domain note in the module docs).
fn invntt_scalar(r: &mut Poly) {
    let mut k = 127usize;
    let mut len = 2usize;
    while len <= 128 {
        let mut start = 0usize;
        while start < N {
            let zeta = ZETAS[k];
            k = k.wrapping_sub(1);
            for j in start..start + len {
                let t = r.0[j];
                // Gentleman–Sande: the sum needs reducing, the difference does
                // not because fqmul brings it back into (-q, q).
                r.0[j] = barrett_reduce(t + r.0[j + len]);
                r.0[j + len] -= t;
                r.0[j + len] = fqmul(zeta, r.0[j + len]);
            }
            start += 2 * len;
        }
        len <<= 1;
    }
    for c in r.0.iter_mut() {
        *c = fqmul(*c, INVNTT_F);
    }
}

/// One pair of the pointwise product: multiply two linear polynomials modulo
/// `X² − γ`.
#[inline(always)]
fn basemul(r: &mut [i16], a: &[i16], b: &[i16], gamma: i16) {
    r[0] = fqmul(fqmul(a[1], b[1]), gamma) + fqmul(a[0], b[0]);
    r[1] = fqmul(a[0], b[1]) + fqmul(a[1], b[0]);
}

/// `r += a ∘ b` in the NTT domain.
fn poly_basemul_acc_scalar(r: &mut Poly, a: &Poly, b: &Poly) {
    let mut t = [0i16; 2];
    for i in 0..128 {
        basemul(
            &mut t,
            &a.0[2 * i..2 * i + 2],
            &b.0[2 * i..2 * i + 2],
            BASEMUL_ZETAS[i],
        );
        r.0[2 * i] = r.0[2 * i].wrapping_add(t[0]);
        r.0[2 * i + 1] = r.0[2 * i + 1].wrapping_add(t[1]);
    }
}

// ── dispatch ─────────────────────────────────────────────────────────────────
//
// The vector kernels in [`avx2`] return exactly what the scalar ones do (see
// that module's docs), so the choice is invisible except in time. It is made
// per call from the cached CPU feature bit; there is no global `target-cpu`
// setting, so one binary runs everywhere.

#[cfg(target_arch = "x86_64")]
mod avx2;

#[inline(always)]
#[cfg(target_arch = "x86_64")]
fn use_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        // Tests force the portable path through the hashing override, so one
        // switch exercises the fully scalar configuration.
        #[cfg(test)]
        if keccak4::FORCE.with(|f| f.get()) == Some(keccak4::Backend::Portable) {
            return false;
        }
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

macro_rules! dispatch {
    ($name:ident, $scalar:ident, $vector:ident, ($($arg:ident: $ty:ty),*)) => {
        #[inline]
        fn $name($($arg: $ty),*) {
            #[cfg(target_arch = "x86_64")]
            if use_avx2() {
                // SAFETY: AVX2 was detected on the running CPU.
                unsafe { avx2::$vector($($arg),*) };
                return;
            }
            $scalar($($arg),*)
        }
    };
}

dispatch!(ntt, ntt_scalar, ntt, (r: &mut Poly));
dispatch!(invntt, invntt_scalar, invntt, (r: &mut Poly));
dispatch!(poly_reduce, poly_reduce_scalar, reduce, (r: &mut Poly));
dispatch!(poly_tomont, poly_tomont_scalar, tomont, (r: &mut Poly));
dispatch!(
    poly_basemul_acc,
    poly_basemul_acc_scalar,
    basemul_acc,
    (r: &mut Poly, a: &Poly, b: &Poly)
);

// ── sampling ─────────────────────────────────────────────────────────────────

/// Rejection-sampling scratch, with spare slots past N so neither parse loop
/// needs a bounds branch: the scalar one writes at most one value past the
/// end, the vector one stores eight lanes at a time.
struct RejBuf {
    c: [i16; N + 16],
    n: usize,
}

impl RejBuf {
    const fn new() -> Self {
        RejBuf {
            c: [0; N + 16],
            n: 0,
        }
    }

    #[inline]
    fn done(&self) -> bool {
        self.n >= N
    }

    /// Consume a whole number of 3-byte groups from `buf`, appending each
    /// 12-bit value below q, until N have been accepted.
    ///
    /// Branchless on the accept test: every candidate is written and the
    /// count advances by the comparison. The second value of a group whose
    /// first filled the polynomial lands in the spare slot and is ignored,
    /// which is what `SampleNTT` does by not looking at it.
    #[inline]
    fn parse(&mut self, buf: &[u8]) {
        debug_assert_eq!(buf.len() % 3, 0);
        let mut n = self.n;
        #[cfg(target_arch = "x86_64")]
        let mut buf = buf;
        #[cfg(target_arch = "x86_64")]
        if use_avx2() {
            // SAFETY: AVX2 was detected on the running CPU.
            let (m, used) = unsafe { avx2::rej_uniform(&mut self.c, n, buf) };
            n = m;
            buf = &buf[used..];
        }
        for g in buf.chunks_exact(3) {
            if n >= N {
                break;
            }
            let c0 = g[0] as u16;
            let c1 = g[1] as u16;
            let c2 = g[2] as u16;
            let d1 = c0 | ((c1 & 0x0f) << 8);
            let d2 = (c1 >> 4) | (c2 << 4);
            self.c[n] = d1 as i16;
            n += (d1 < Q as u16) as usize;
            self.c[n] = d2 as i16;
            n += (d2 < Q as u16) as usize;
        }
        self.n = n;
    }

    fn poly(&self) -> Poly {
        let mut f = Poly::zero();
        f.0.copy_from_slice(&self.c[..N]);
        f
    }
}

/// `SampleNTT`: rejection-sample a uniform NTT-domain polynomial from
/// SHAKE128(ρ ‖ j ‖ i).
///
/// The sponge is squeezed a block at a time and continues where it left off,
/// so an unlucky draw costs one more permutation rather than a restart.
fn sample_ntt(rho: &[u8; 32], j: u8, i: u8) -> Poly {
    let mut xof = Sponge::<SHAKE128_RATE>::new();
    xof.absorb(rho);
    xof.absorb(&[j, i]);
    xof.finalize(0x1f);

    let mut buf = [0u8; SHAKE128_RATE];
    let mut r = RejBuf::new();
    while !r.done() {
        xof.squeeze(&mut buf);
        r.parse(&buf);
    }
    r.poly()
}

/// Four `SampleNTT`s in four SIMD lanes. `idx[l] = (j, i)` names the matrix
/// entry lane `l` produces.
///
/// Three blocks are squeezed up front, since 504 bytes give 336 candidates
/// for 256 slots at an acceptance rate of 3329/4096 and almost always
/// suffice; any lane still short gets further single blocks, all four lanes
/// permuted together because that costs no more than one.
fn sample_ntt_x4(rho: &[u8; 32], idx: [(u8, u8); 4]) -> [Poly; 4] {
    const R: usize = SHAKE128_RATE;
    let mut inp = [[0u8; 34]; 4];
    for l in 0..4 {
        inp[l][..32].copy_from_slice(rho);
        inp[l][32] = idx[l].0;
        inp[l][33] = idx[l].1;
    }
    let mut s = keccak4::absorb_short_x4::<R>([&inp[0], &inp[1], &inp[2], &inp[3]], 0x1f);
    let mut r = [RejBuf::new(), RejBuf::new(), RejBuf::new(), RejBuf::new()];
    let mut buf = [0u8; R];
    let mut first = true;
    loop {
        // Three blocks on the first pass, one on each later pass.
        for blk in 0..(if first { 3 } else { 1 }) {
            if blk > 0 || !first {
                keccak4::keccak_f1600_x4(&mut s);
            }
            for l in 0..4 {
                if !r[l].done() {
                    keccak4::extract_block::<R>(&s, l, &mut buf);
                    r[l].parse(&buf);
                }
            }
        }
        first = false;
        if r.iter().all(|x| x.done()) {
            break;
        }
    }
    [r[0].poly(), r[1].poly(), r[2].poly(), r[3].poly()]
}

/// `SamplePolyCBD_2` from 128 bytes: each coefficient is a difference of two
/// two-bit popcounts, computed four coefficients at a time out of one word.
fn cbd2(buf: &[u8; 128]) -> Poly {
    let mut f = Poly::zero();
    for i in 0..N / 8 {
        let t = u32::from_le_bytes([buf[4 * i], buf[4 * i + 1], buf[4 * i + 2], buf[4 * i + 3]]);
        // Sum adjacent bit pairs in parallel across the word.
        let d = (t & 0x5555_5555) + ((t >> 1) & 0x5555_5555);
        for j in 0..8 {
            let a = ((d >> (4 * j)) & 0x3) as i16;
            let b = ((d >> (4 * j + 2)) & 0x3) as i16;
            f.0[8 * i + j] = a - b;
        }
    }
    f
}

/// `SamplePolyCBD_3` from 192 bytes: three-bit popcounts, four coefficients per
/// three bytes.
fn cbd3(buf: &[u8; 192]) -> Poly {
    let mut f = Poly::zero();
    for i in 0..N / 4 {
        let t = u32::from_le_bytes([buf[3 * i], buf[3 * i + 1], buf[3 * i + 2], 0]);
        let d = (t & 0x0024_9249) + ((t >> 1) & 0x0024_9249) + ((t >> 2) & 0x0024_9249);
        for j in 0..4 {
            let a = ((d >> (6 * j)) & 0x7) as i16;
            let b = ((d >> (6 * j + 3)) & 0x7) as i16;
            f.0[4 * i + j] = a - b;
        }
    }
    f
}

/// `SamplePolyCBD_η(PRF_η(s, b))`.
fn sample_cbd_prf(s: &[u8; 32], b: u8, eta: usize) -> Poly {
    if eta == 2 {
        let mut buf = [0u8; 128];
        shake256_2_into(&mut buf, s, &[b]);
        cbd2(&buf)
    } else {
        debug_assert_eq!(eta, 3);
        let mut buf = [0u8; 192];
        shake256_2_into(&mut buf, s, &[b]);
        cbd3(&buf)
    }
}

/// `out[t] = SamplePolyCBD_{etas[t]}(PRF(s, first + t))` for every `t`.
///
/// With four-way hashing the calls go four at a time. A group squeezes the
/// longest output any member needs (two blocks if any is η = 3), and members
/// with η = 2 read the first 128 bytes of it, which is their PRF output since
/// SHAKE is a stream. A lone leftover call goes through the vector unit only
/// when that is cheaper than one scalar permutation, which it is with
/// AVX-512 and is not with AVX2.
fn sample_noise(s: &[u8; 32], first: u8, etas: &[usize], out: &mut [Poly]) {
    debug_assert_eq!(etas.len(), out.len());
    let n = etas.len();
    let backend = keccak4::backend();
    if backend == keccak4::Backend::Portable {
        for t in 0..n {
            out[t] = sample_cbd_prf(s, first + t as u8, etas[t]);
        }
        return;
    }
    let mut t = 0;
    while t < n {
        let left = n - t;
        if left == 1 && backend != keccak4::Backend::Avx512 {
            out[t] = sample_cbd_prf(s, first + t as u8, etas[t]);
            break;
        }
        let used = left.min(4);
        let mut inp = [[0u8; 33]; 4];
        for l in 0..4 {
            inp[l][..32].copy_from_slice(s);
            // Unused lanes repeat the last real nonce and are discarded.
            inp[l][32] = first + (t + l.min(used - 1)) as u8;
        }
        let eta_max = etas[t..t + used].iter().copied().max().unwrap();
        let len = 64 * eta_max;
        let mut b = [[0u8; 192]; 4];
        {
            let [b0, b1, b2, b3] = &mut b;
            keccak4::shake256_x4(
                [
                    &mut b0[..len],
                    &mut b1[..len],
                    &mut b2[..len],
                    &mut b3[..len],
                ],
                [&inp[0], &inp[1], &inp[2], &inp[3]],
            );
        }
        for l in 0..used {
            out[t + l] = if etas[t + l] == 2 {
                cbd2(b[l][..128].try_into().unwrap())
            } else {
                debug_assert_eq!(etas[t + l], 3);
                cbd3(&b[l])
            };
        }
        t += used;
    }
}

// ── compression and byte encoding ────────────────────────────────────────────

/// `Compress_d(x)` for canonical `x`, exactly as FIPS 203 §4.2.1 defines it.
#[inline(always)]
fn compress(x: u16, d: usize) -> u16 {
    ((((x as u32) << d) + (Q as u32) / 2) / (Q as u32)) as u16 & ((1u16 << d) - 1)
}

/// `Decompress_d(y)`.
#[inline(always)]
fn decompress(y: u16, d: usize) -> i16 {
    (((y as u32) * (Q as u32) + (1 << (d - 1))) >> d) as i16
}

/// Pack 256 values of `D` bits each, eight at a time: eight `D`-bit values
/// are exactly `D` bytes, for every `D`, so each group is assembled in one
/// register and written with one copy.
///
/// This replaced a byte-at-a-time accumulator whose inner `while` loop and
/// per-byte store were 18% of a prepared ML-KEM-768 encapsulation once the
/// hashing and the NTT had been vectorised. `D` is a const parameter so each
/// width gets its own fully unrolled body.
#[inline(always)]
fn pack8<const D: usize>(out: &mut [u8], mut val: impl FnMut(usize) -> u16) {
    debug_assert_eq!(out.len(), 32 * D);
    for (g, chunk) in out.chunks_exact_mut(D).enumerate() {
        if D <= 8 {
            let mut acc = 0u64;
            for j in 0..8 {
                acc |= (val(8 * g + j) as u64) << (D * j);
            }
            chunk.copy_from_slice(&acc.to_le_bytes()[..D]);
        } else {
            let mut acc = 0u128;
            for j in 0..8 {
                acc |= (val(8 * g + j) as u128) << (D * j);
            }
            chunk.copy_from_slice(&acc.to_le_bytes()[..D]);
        }
    }
}

/// The inverse of [`pack8`]: hand each `D`-bit field to `put(index, value)`.
#[inline(always)]
fn unpack8<const D: usize>(bytes: &[u8], mut put: impl FnMut(usize, u16)) {
    debug_assert_eq!(bytes.len(), 32 * D);
    let mask = ((1u32 << D) - 1) as u16;
    for (g, chunk) in bytes.chunks_exact(D).enumerate() {
        if D <= 8 {
            let mut b = [0u8; 8];
            b[..D].copy_from_slice(chunk);
            let acc = u64::from_le_bytes(b);
            for j in 0..8 {
                put(8 * g + j, (acc >> (D * j)) as u16 & mask);
            }
        } else {
            let mut b = [0u8; 16];
            b[..D].copy_from_slice(chunk);
            let acc = u128::from_le_bytes(b);
            for j in 0..8 {
                put(8 * g + j, (acc >> (D * j)) as u16 & mask);
            }
        }
    }
}

/// Instantiate `$body` once per encoding width ML-KEM uses, with `$D` bound
/// to it as a constant.
macro_rules! with_width {
    ($d:expr, $D:ident => $body:expr) => {
        match $d {
            1 => {
                const $D: usize = 1;
                $body
            }
            4 => {
                const $D: usize = 4;
                $body
            }
            5 => {
                const $D: usize = 5;
                $body
            }
            10 => {
                const $D: usize = 10;
                $body
            }
            11 => {
                const $D: usize = 11;
                $body
            }
            12 => {
                const $D: usize = 12;
                $body
            }
            _ => unreachable!("ML-KEM encodes only at d in {{1, 4, 5, 10, 11, 12}}"),
        }
    };
}

/// `ByteEncode_d` of a polynomial whose coefficients are already canonical.
fn byte_encode_into(out: &mut [u8], f: &Poly, d: usize) {
    with_width!(d, D => pack8::<D>(out, |i| f.0[i] as u16 & ((1 << D) - 1) as u16))
}

/// `ByteDecode_d`. For `d = 12` the values are reduced mod q, which makes the
/// decode lossy exactly as the reference documents; callers must not treat a
/// successful decode as proof that the input was canonical. A 12-bit value is
/// below `2q`, so the reduction is one conditional subtraction.
fn byte_decode(bytes: &[u8], d: usize) -> Poly {
    let mut f = Poly::zero();
    with_width!(d, D => unpack8::<D>(bytes, |i, v| {
        f.0[i] = if D == 12 {
            (v as i16) - (Q & -((v >= Q as u16) as i16))
        } else {
            v as i16
        };
    }));
    f
}

/// `ByteEncode_d(Compress_d(f))`, with `f` in any representative.
fn poly_compress_encode(out: &mut [u8], f: &Poly, d: usize) {
    #[cfg(target_arch = "x86_64")]
    if use_avx2() {
        let mut v = [0u16; N];
        // SAFETY: AVX2 was detected on the running CPU.
        unsafe { avx2::compress(f, d as u32, &mut v) };
        return with_width!(d, D => pack8::<D>(out, |i| v[i]));
    }
    poly_compress_encode_scalar(out, f, d)
}

fn poly_compress_encode_scalar(out: &mut [u8], f: &Poly, d: usize) {
    with_width!(d, D => pack8::<D>(out, |i| compress(to_canonical(barrett_reduce(f.0[i])), D)))
}

/// `Decompress_d(ByteDecode_d(bytes))`.
fn poly_decode_decompress(bytes: &[u8], d: usize) -> Poly {
    let mut r = Poly::zero();
    with_width!(d, D => unpack8::<D>(bytes, |i, v| r.0[i] = decompress(v, D)));
    r
}

// ── K-PKE ────────────────────────────────────────────────────────────────────

/// The k×k matrix `Â`, row-major with stride `KMAX`.
type Matrix = [[Poly; KMAX]; KMAX];

/// Expand ρ into the k×k matrix `Â` in NTT form. Entry `(i, j)` comes from
/// `XOF(ρ ‖ j ‖ i)`: the column index goes first, as in the final standard.
///
/// With four-way hashing the k² entries go four at a time: one batch for
/// ML-KEM-512, two and a remainder of one for ML-KEM-768, four for
/// ML-KEM-1024. The remainder goes through the vector unit only when that is
/// cheaper than one scalar sponge (AVX-512, not AVX2).
fn expand_matrix(rho: &[u8; 32], k: usize, a: &mut Matrix) {
    let backend = keccak4::backend();
    if backend == keccak4::Backend::Portable {
        for i in 0..k {
            for j in 0..k {
                a[i][j] = sample_ntt(rho, j as u8, i as u8);
            }
        }
        return;
    }
    let total = k * k;
    let mut t = 0;
    while t < total {
        let left = total - t;
        if left == 1 && backend != keccak4::Backend::Avx512 {
            let (i, j) = (t / k, t % k);
            a[i][j] = sample_ntt(rho, j as u8, i as u8);
            break;
        }
        let used = left.min(4);
        let mut idx = [(0u8, 0u8); 4];
        for (l, item) in idx.iter_mut().enumerate() {
            // Unused lanes repeat the last real entry and are discarded.
            let e = t + l.min(used - 1);
            *item = ((e % k) as u8, (e / k) as u8);
        }
        let polys = sample_ntt_x4(rho, idx);
        for (l, f) in polys.into_iter().enumerate().take(used) {
            let e = t + l;
            a[e / k][e % k] = f;
        }
        t += used;
    }
}

/// Decode `t̂` from the first `384·k` bytes of an encapsulation key.
fn decode_t_hat(k: usize, ek: &[u8], t_hat: &mut [Poly; KMAX]) {
    for (i, item) in t_hat.iter_mut().enumerate().take(k) {
        *item = byte_decode(&ek[384 * i..384 * (i + 1)], 12);
    }
}

/// Canonicalise in place, so that a value cached from key generation and the
/// same value decoded from bytes are equal as arrays, not just mod q.
fn poly_canonical(r: &mut Poly) {
    for c in r.0.iter_mut() {
        *c = to_canonical(barrett_reduce(*c)) as i16;
    }
}

/// `K-PKE.KeyGen`, also handing back `Â`, `t̂` and `ŝ`, which it computes
/// anyway, so that a caller who wants the public-data cache gets it without
/// redoing the work (the paper's "capture `A` during key generation").
fn kpke_keygen(
    p: &MlKemParams,
    d: &[u8; 32],
    ek: &mut [u8],
    dk: &mut [u8],
    a: &mut Matrix,
    t_out: &mut [Poly; KMAX],
    s_out: &mut [Poly; KMAX],
) {
    let k = p.k;
    let g = sha3_512_2(d, &[k as u8]);
    let mut rho = [0u8; 32];
    let mut sigma = [0u8; 32];
    rho.copy_from_slice(&g[..32]);
    sigma.copy_from_slice(&g[32..]);

    expand_matrix(&rho, k, a);

    // ŝ and ê: 2k noise polynomials under σ with nonces 0..2k.
    let mut noise = [Poly::zero(); 2 * KMAX];
    sample_noise(&sigma, 0, &[p.eta1; 2 * KMAX][..2 * k], &mut noise[..2 * k]);
    for f in noise[..2 * k].iter_mut() {
        ntt(f);
    }
    let s_hat = s_out;
    s_hat[..k].copy_from_slice(&noise[..k]);
    let mut e_hat = [Poly::zero(); KMAX];
    e_hat[..k].copy_from_slice(&noise[k..2 * k]);
    for f in noise.iter_mut() {
        f.0.zeroize();
    }

    // t̂ = Â∘ŝ + ê. The accumulated product is short by 2^16 (see the domain
    // note), so it is scaled back before ê joins it.
    for i in 0..k {
        let mut acc = Poly::zero();
        for j in 0..k {
            poly_basemul_acc(&mut acc, &a[i][j], &s_hat[j]);
        }
        poly_reduce(&mut acc);
        poly_tomont(&mut acc);
        poly_add(&mut acc, &e_hat[i]);
        poly_canonical(&mut acc);
        byte_encode_into(&mut ek[384 * i..384 * (i + 1)], &acc, 12);
        t_out[i] = acc;
    }
    ek[384 * k..384 * k + 32].copy_from_slice(&rho);

    for i in 0..k {
        poly_canonical(&mut s_hat[i]);
        byte_encode_into(&mut dk[384 * i..384 * (i + 1)], &s_hat[i], 12);
    }
}

/// `K-PKE.Encrypt` given an already expanded `Â` and decoded `t̂`.
///
/// This is the whole of encryption except the two pieces of work that depend
/// only on the public key, which is what makes them cacheable.
fn kpke_encrypt_with(
    p: &MlKemParams,
    a: &Matrix,
    t_hat: &[Poly; KMAX],
    m: &[u8; 32],
    r: &[u8; 32],
    c: &mut [u8],
) {
    let k = p.k;
    // ŷ, e1 and e2: 2k + 1 noise polynomials under r with nonces 0..=2k.
    let mut etas = [0usize; 2 * KMAX + 1];
    for (t, e) in etas.iter_mut().enumerate().take(2 * k + 1) {
        *e = if t < k { p.eta1 } else { p.eta2 };
    }
    let mut noise = [Poly::zero(); 2 * KMAX + 1];
    sample_noise(r, 0, &etas[..2 * k + 1], &mut noise[..2 * k + 1]);
    let mut y_hat = [Poly::zero(); KMAX];
    let mut e1 = [Poly::zero(); KMAX];
    for i in 0..k {
        y_hat[i] = noise[i];
        ntt(&mut y_hat[i]);
        e1[i] = noise[k + i];
    }
    let e2 = noise[2 * k];

    // u = NTT^-1(Â^T ∘ ŷ) + e1
    let du_bytes = 32 * p.du;
    for i in 0..k {
        let mut acc = Poly::zero();
        for j in 0..k {
            poly_basemul_acc(&mut acc, &a[j][i], &y_hat[j]);
        }
        poly_reduce(&mut acc);
        invntt(&mut acc);
        poly_add(&mut acc, &e1[i]);
        poly_reduce(&mut acc);
        poly_compress_encode(&mut c[du_bytes * i..du_bytes * (i + 1)], &acc, p.du);
    }

    // v = NTT^-1(t̂^T ∘ ŷ) + e2 + Decompress_1(m)
    let mut v = Poly::zero();
    for j in 0..k {
        poly_basemul_acc(&mut v, &t_hat[j], &y_hat[j]);
    }
    poly_reduce(&mut v);
    invntt(&mut v);
    poly_add(&mut v, &e2);
    let mu = poly_decode_decompress(m, 1);
    poly_add(&mut v, &mu);
    poly_reduce(&mut v);
    poly_compress_encode(&mut c[du_bytes * k..], &v, p.dv);
}

/// `K-PKE.Encrypt` from the encoded key: decode `t̂`, expand `Â`, encrypt.
fn kpke_encrypt(p: &MlKemParams, ek: &[u8], m: &[u8; 32], r: &[u8; 32], c: &mut [u8]) {
    let k = p.k;
    let mut t_hat = [Poly::zero(); KMAX];
    decode_t_hat(k, ek, &mut t_hat);
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&ek[384 * k..384 * k + 32]);
    let mut a = [[Poly::zero(); KMAX]; KMAX];
    expand_matrix(&rho, k, &mut a);
    kpke_encrypt_with(p, &a, &t_hat, m, r, c);
}

/// `K-PKE.Decrypt` given a decoded `ŝ`.
fn kpke_decrypt_with(p: &MlKemParams, s_hat: &[Poly; KMAX], c: &[u8]) -> [u8; 32] {
    let k = p.k;
    let du_bytes = 32 * p.du;
    let mut w = Poly::zero();
    for i in 0..k {
        let mut u = poly_decode_decompress(&c[du_bytes * i..du_bytes * (i + 1)], p.du);
        ntt(&mut u);
        poly_basemul_acc(&mut w, &s_hat[i], &u);
    }
    poly_reduce(&mut w);
    invntt(&mut w);

    let mut v = poly_decode_decompress(&c[du_bytes * k..], p.dv);
    poly_sub(&mut v, &w);
    poly_reduce(&mut v);

    let mut out = [0u8; 32];
    poly_compress_encode(&mut out, &v, 1);
    out
}

fn kpke_decrypt(p: &MlKemParams, dk: &[u8], c: &[u8]) -> [u8; 32] {
    let mut s_hat = [Poly::zero(); KMAX];
    decode_t_hat(p.k, dk, &mut s_hat);
    kpke_decrypt_with(p, &s_hat, c)
}

// ── ML-KEM ───────────────────────────────────────────────────────────────────

/// `ML-KEM.KeyGen_internal(d, z)`. Deterministic; the seeds are the caller's.
pub fn ml_kem_keygen_internal(
    p: &MlKemParams,
    d: &[u8; 32],
    z: &[u8; 32],
) -> (MlKemEncapsKey, MlKemDecapsKey) {
    let mut a = [[Poly::zero(); KMAX]; KMAX];
    let mut t = [Poly::zero(); KMAX];
    let mut s = [Poly::zero(); KMAX];
    let (ek, dk, _) = keygen_parts(p, d, z, &mut a, &mut t, &mut s);
    s.iter_mut().for_each(|x| x.0.zeroize());
    (ek, dk)
}

/// The key generation body, returning `H(ek)` alongside the keys so that the
/// preparing variant does not hash the key twice.
fn keygen_parts(
    p: &MlKemParams,
    d: &[u8; 32],
    z: &[u8; 32],
    a: &mut Matrix,
    t: &mut [Poly; KMAX],
    s: &mut [Poly; KMAX],
) -> (MlKemEncapsKey, MlKemDecapsKey, [u8; 32]) {
    let k = p.k;
    let mut ek = vec![0u8; p.ek_len()];
    let mut dk = vec![0u8; p.dk_len()];
    kpke_keygen(p, d, &mut ek, &mut dk[..384 * k], a, t, s);
    let h = sha3_256(&ek);
    dk[384 * k..768 * k + 32].copy_from_slice(&ek);
    dk[768 * k + 32..768 * k + 64].copy_from_slice(&h);
    dk[768 * k + 64..].copy_from_slice(z);
    (MlKemEncapsKey(ek), MlKemDecapsKey(dk), h)
}

/// `ML-KEM.KeyGen`.
pub fn ml_kem_keygen(p: &MlKemParams) -> (MlKemEncapsKey, MlKemDecapsKey) {
    let mut d = [0u8; 32];
    let mut z = [0u8; 32];
    random_bytes(&mut d);
    random_bytes(&mut z);
    ml_kem_keygen_internal(p, &d, &z)
}

/// `ML-KEM.Encaps_internal(ek, m)`. Deterministic; the message is the caller's.
pub fn ml_kem_encaps_internal(
    p: &MlKemParams,
    ek: &MlKemEncapsKey,
    m: &[u8; 32],
) -> Option<(Vec<u8>, [u8; SHARED_SECRET_BYTES])> {
    if ek.0.len() != p.ek_len() {
        return None;
    }
    let g = sha3_512_2(m, &sha3_256(&ek.0));
    let mut key = [0u8; 32];
    let mut r = [0u8; 32];
    key.copy_from_slice(&g[..32]);
    r.copy_from_slice(&g[32..]);

    let mut c = vec![0u8; p.ct_len()];
    kpke_encrypt(p, &ek.0, m, &r, &mut c);
    Some((c, key))
}

/// `ML-KEM.Encaps`.
pub fn ml_kem_encaps(
    p: &MlKemParams,
    ek: &MlKemEncapsKey,
) -> Option<(Vec<u8>, [u8; SHARED_SECRET_BYTES])> {
    let mut m = [0u8; 32];
    random_bytes(&mut m);
    ml_kem_encaps_internal(p, ek, &m)
}

/// The FO tail shared by both decapsulation paths: derive `(K', r')`, compute
/// the rejection key, re-encrypt, and select without branching.
///
/// `J(z ‖ c)` is computed on every call, valid ciphertext or not. Skipping it
/// on success would be the cheapest saving left in decapsulation and would
/// also hand an attacker a timing oracle for ciphertext validity, which is the
/// plaintext-checking oracle that chosen-ciphertext attacks on Kyber are built
/// from. So it stays.
fn decaps_tail(
    p: &MlKemParams,
    m_prime: &[u8; 32],
    h: &[u8],
    z: &[u8],
    c: &[u8],
    reencrypt: impl FnOnce(&[u8; 32], &[u8; 32], &mut [u8]),
) -> [u8; SHARED_SECRET_BYTES] {
    let g = sha3_512_2(m_prime, h);
    let mut k_prime = [0u8; 32];
    let mut r_prime = [0u8; 32];
    k_prime.copy_from_slice(&g[..32]);
    r_prime.copy_from_slice(&g[32..]);

    let mut k_bar = [0u8; 32];
    shake256_2_into(&mut k_bar, z, c);

    let mut c_prime = vec![0u8; p.ct_len()];
    reencrypt(m_prime, &r_prime, &mut c_prime);

    // Implicit rejection: select without branching on the comparison.
    let same = c_prime.ct_eq(c);
    let mut out = [0u8; 32];
    for i in 0..32 {
        out[i] = u8::conditional_select(&k_bar[i], &k_prime[i], same);
    }
    k_prime.zeroize();
    out
}

/// `ML-KEM.Decaps`. Returns `None` only for malformed input lengths; a
/// well-formed but forged ciphertext yields the implicit-rejection secret
/// `J(z ‖ c)`, never an error.
pub fn ml_kem_decaps(
    p: &MlKemParams,
    dk: &MlKemDecapsKey,
    c: &[u8],
) -> Option<[u8; SHARED_SECRET_BYTES]> {
    let k = p.k;
    if dk.0.len() != p.dk_len() || c.len() != p.ct_len() {
        return None;
    }
    let dk_pke = &dk.0[..384 * k];
    let ek = &dk.0[384 * k..768 * k + 32];
    let h = &dk.0[768 * k + 32..768 * k + 64];
    let z = &dk.0[768 * k + 64..];

    let m_prime = kpke_decrypt(p, dk_pke, c);
    Some(decaps_tail(p, &m_prime, h, z, c, |m, r, out| {
        kpke_encrypt(p, ek, m, r, out)
    }))
}

/// `ML-KEM.EncapsKeyCheck`: the §7.2 modulus check, by decode/re-encode.
pub fn ml_kem_check_ek(p: &MlKemParams, ek: &[u8]) -> bool {
    if ek.len() != p.ek_len() {
        return false;
    }
    let mut buf = [0u8; 384];
    for i in 0..p.k {
        let f = byte_decode(&ek[384 * i..384 * (i + 1)], 12);
        byte_encode_into(&mut buf, &f, 12);
        if buf[..] != ek[384 * i..384 * (i + 1)] {
            return false;
        }
    }
    true
}

// ── public-data reuse ────────────────────────────────────────────────────────
//
// Sayed, Taha and Nijjer (arXiv:2610.01960) measure, on a Cortex-M7, that
// once the arithmetic kernels are tuned the largest remaining cost in ML-KEM
// is work that depends only on public data and is redone on every call:
// expanding `Â` from ρ, hashing the encapsulation key, and decoding `t̂`.
// None of it involves a secret, so it can be computed once per key and kept.
//
// The types below are that cache. The standardised byte formats of keys,
// ciphertexts and shared secrets are untouched: a prepared key is built
// *from* those bytes (or captured during key generation, which computes all
// of it anyway) and every output is byte-identical to the unprepared path,
// which the tests check against the reference. The cache is bound to its key
// by construction, since it owns a copy of the key it was derived from; there
// is no way to pair a cache with a different key through this API.
//
// Cost, per key: k² + k polynomials of 512 bytes, plus the key bytes. That is
// 6 KiB for ML-KEM-512, 12 KiB for ML-KEM-768, 20 KiB for ML-KEM-1024.

/// An encapsulation key with its public-only derived data precomputed:
/// `Â`, `t̂` and `H(ek)`.
#[derive(Clone)]
pub struct MlKemPreparedEncapsKey {
    p: MlKemParams,
    a: Box<Matrix>,
    t_hat: [Poly; KMAX],
    h_ek: [u8; 32],
    ek: MlKemEncapsKey,
}

impl MlKemPreparedEncapsKey {
    /// Prepare from the standard encoding. `None` on a wrong length, exactly
    /// where [`ml_kem_encaps_internal`] would return `None`.
    pub fn new(p: &MlKemParams, ek: &MlKemEncapsKey) -> Option<Self> {
        if ek.0.len() != p.ek_len() {
            return None;
        }
        let k = p.k;
        let mut t_hat = [Poly::zero(); KMAX];
        decode_t_hat(k, &ek.0, &mut t_hat);
        let mut rho = [0u8; 32];
        rho.copy_from_slice(&ek.0[384 * k..384 * k + 32]);
        let mut a = Box::new([[Poly::zero(); KMAX]; KMAX]);
        expand_matrix(&rho, k, &mut a);
        Some(MlKemPreparedEncapsKey {
            p: *p,
            a,
            t_hat,
            h_ek: sha3_256(&ek.0),
            ek: ek.clone(),
        })
    }

    /// The key this cache was derived from.
    pub fn encaps_key(&self) -> &MlKemEncapsKey {
        &self.ek
    }

    /// `H(ek)`, which ML-KEM also stores inside the decapsulation key.
    pub fn h_ek(&self) -> &[u8; 32] {
        &self.h_ek
    }

    pub fn params(&self) -> &MlKemParams {
        &self.p
    }

    /// `ML-KEM.Encaps_internal(ek, m)` against the cache. Byte-identical to
    /// [`ml_kem_encaps_internal`] on the same key and message.
    pub fn encaps_internal(&self, m: &[u8; 32]) -> (Vec<u8>, [u8; SHARED_SECRET_BYTES]) {
        let p = &self.p;
        let g = sha3_512_2(m, &self.h_ek);
        let mut key = [0u8; 32];
        let mut r = [0u8; 32];
        key.copy_from_slice(&g[..32]);
        r.copy_from_slice(&g[32..]);
        let mut c = vec![0u8; p.ct_len()];
        kpke_encrypt_with(p, &self.a, &self.t_hat, m, &r, &mut c);
        (c, key)
    }

    /// `ML-KEM.Encaps` against the cache, with fresh randomness.
    pub fn encaps(&self) -> (Vec<u8>, [u8; SHARED_SECRET_BYTES]) {
        let mut m = [0u8; 32];
        random_bytes(&mut m);
        self.encaps_internal(&m)
    }
}

/// A decapsulation key with `ŝ` decoded and the public cache of its embedded
/// encapsulation key. Re-encryption inside decapsulation is where the cache
/// pays off: it removes the matrix expansion from every call.
#[derive(Clone)]
pub struct MlKemPreparedDecapsKey {
    s_hat: [Poly; KMAX],
    z: [u8; 32],
    ek: MlKemPreparedEncapsKey,
}

impl Drop for MlKemPreparedDecapsKey {
    fn drop(&mut self) {
        for s in self.s_hat.iter_mut() {
            s.0.zeroize();
        }
        self.z.zeroize();
    }
}

impl MlKemPreparedDecapsKey {
    /// Prepare from the standard encoding. `None` on a wrong length.
    ///
    /// The `H(ek)` stored in `dk` is used as stored, as [`ml_kem_decaps`] uses
    /// it, rather than recomputed; the encapsulation side of this cache
    /// therefore carries whatever `dk` says, which for a well-formed key is
    /// the true hash.
    pub fn new(p: &MlKemParams, dk: &MlKemDecapsKey) -> Option<Self> {
        if dk.0.len() != p.dk_len() {
            return None;
        }
        let k = p.k;
        let mut s_hat = [Poly::zero(); KMAX];
        decode_t_hat(k, &dk.0[..384 * k], &mut s_hat);
        let ek = MlKemEncapsKey(dk.0[384 * k..768 * k + 32].to_vec());
        let mut t_hat = [Poly::zero(); KMAX];
        decode_t_hat(k, &ek.0, &mut t_hat);
        let mut rho = [0u8; 32];
        rho.copy_from_slice(&ek.0[384 * k..384 * k + 32]);
        let mut a = Box::new([[Poly::zero(); KMAX]; KMAX]);
        expand_matrix(&rho, k, &mut a);
        let mut h_ek = [0u8; 32];
        h_ek.copy_from_slice(&dk.0[768 * k + 32..768 * k + 64]);
        let mut z = [0u8; 32];
        z.copy_from_slice(&dk.0[768 * k + 64..]);
        Some(MlKemPreparedDecapsKey {
            s_hat,
            z,
            ek: MlKemPreparedEncapsKey {
                p: *p,
                a,
                t_hat,
                h_ek,
                ek,
            },
        })
    }

    /// The public half, for handing to an encapsulator that wants it.
    pub fn encaps_key(&self) -> &MlKemPreparedEncapsKey {
        &self.ek
    }

    /// `ML-KEM.Decaps` against the cache. Byte-identical to [`ml_kem_decaps`]
    /// on the same key and ciphertext, including the implicit-rejection key.
    pub fn decaps(&self, c: &[u8]) -> Option<[u8; SHARED_SECRET_BYTES]> {
        let p = &self.ek.p;
        if c.len() != p.ct_len() {
            return None;
        }
        let m_prime = kpke_decrypt_with(p, &self.s_hat, c);
        let ek = &self.ek;
        Some(decaps_tail(
            p,
            &m_prime,
            &ek.h_ek,
            &self.z,
            c,
            |m, r, out| kpke_encrypt_with(p, &ek.a, &ek.t_hat, m, r, out),
        ))
    }
}

/// `ML-KEM.KeyGen_internal(d, z)`, also returning the prepared decapsulation
/// key. Key generation computes `Â`, `t̂`, `ŝ` and `H(ek)` anyway, so this
/// costs only the copies; the returned keys are byte-identical to
/// [`ml_kem_keygen_internal`]'s.
pub fn ml_kem_keygen_prepared_internal(
    p: &MlKemParams,
    d: &[u8; 32],
    z: &[u8; 32],
) -> (MlKemEncapsKey, MlKemDecapsKey, MlKemPreparedDecapsKey) {
    let mut a = Box::new([[Poly::zero(); KMAX]; KMAX]);
    let mut t_hat = [Poly::zero(); KMAX];
    let mut s_hat = [Poly::zero(); KMAX];
    let (ek, dk, h_ek) = keygen_parts(p, d, z, &mut a, &mut t_hat, &mut s_hat);
    let prepared = MlKemPreparedDecapsKey {
        s_hat,
        z: *z,
        ek: MlKemPreparedEncapsKey {
            p: *p,
            a,
            t_hat,
            h_ek,
            ek: ek.clone(),
        },
    };
    s_hat.iter_mut().for_each(|x| x.0.zeroize());
    (ek, dk, prepared)
}

/// `ML-KEM.KeyGen` with the prepared decapsulation key.
pub fn ml_kem_keygen_prepared(
    p: &MlKemParams,
) -> (MlKemEncapsKey, MlKemDecapsKey, MlKemPreparedDecapsKey) {
    let mut d = [0u8; 32];
    let mut z = [0u8; 32];
    random_bytes(&mut d);
    random_bytes(&mut z);
    let out = ml_kem_keygen_prepared_internal(p, &d, &z);
    d.zeroize();
    z.zeroize();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pqc::ml_kem as slow;

    /// A cheap deterministic generator, so a failure names a seed.
    fn seeded(n: usize, tag: u8) -> Vec<u8> {
        let mut out = vec![0u8; n];
        let mut x = 0x9e37_79b9_7f4a_7c15u64 ^ ((tag as u64) << 32);
        for b in out.iter_mut() {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *b = (x >> 24) as u8;
        }
        out
    }

    fn sets() -> [MlKemParams; 3] {
        [slow::ML_KEM_512, slow::ML_KEM_768, slow::ML_KEM_1024]
    }

    /// The constants the domain argument rests on. If any of these is wrong the
    /// scheme is wrong by a constant factor, which is exactly the bug class the
    /// Montgomery rewrite invites.
    #[test]
    fn montgomery_constants_are_what_they_claim() {
        assert_eq!((Q as i32 * QINV as i32) as i16, 1, "q * q^-1 = 1 mod 2^16");
        assert_eq!(MONT2 as i32, (1i64 << 32).rem_euclid(QI as i64) as i32);
        assert_eq!(
            (INVNTT_F as i64 * 128) % QI as i64,
            MONT2 as i64,
            "invntt scaling folds 1/128 into the Montgomery factor"
        );
        // fqmul with a Montgomery-form constant is the true product.
        for v in [1i32, 2, 17, 1000, 3328] {
            let m = to_mont_signed(v);
            assert_eq!(
                to_canonical(barrett_reduce(fqmul(5, m))) as i32,
                (5 * v).rem_euclid(QI),
                "fqmul(5, mont({v}))"
            );
        }
    }

    /// The twiddle tables, against a schoolbook exponentiation written here so
    /// the check does not share any code with the `const fn` it is checking.
    #[test]
    fn twiddle_tables_hold_the_right_powers() {
        fn pow(mut e: u32) -> u32 {
            let (mut r, mut b) = (1u32, 17u32);
            while e > 0 {
                if e & 1 == 1 {
                    r = r * b % 3329;
                }
                b = b * b % 3329;
                e >>= 1;
            }
            r
        }
        fn brv(x: u32) -> u32 {
            (0..7).fold(0, |r, b| r | (((x >> b) & 1) << (6 - b)))
        }
        for i in 0..128u32 {
            let k = i as usize;
            // fqmul(1, mont(v)) = v, so this reads the table back out of the
            // Montgomery domain without using to_mont_signed again.
            assert_eq!(
                to_canonical(barrett_reduce(fqmul(1, ZETAS[k]))) as u32,
                pow(brv(i)),
                "ZETAS[{i}]"
            );
            assert_eq!(
                to_canonical(barrett_reduce(fqmul(1, BASEMUL_ZETAS[k]))) as u32,
                pow(2 * brv(i) + 1),
                "BASEMUL_ZETAS[{i}]"
            );
        }
    }

    /// NTT then inverse NTT is the identity, which pins the scaling constants
    /// independently of everything above.
    #[test]
    fn ntt_round_trips() {
        for tag in 0..8u8 {
            let bytes = seeded(N * 2, tag);
            let mut f = Poly::zero();
            for i in 0..N {
                f.0[i] = (u16::from_le_bytes([bytes[2 * i], bytes[2 * i + 1]]) % Q as u16) as i16;
            }
            let orig = f;
            ntt(&mut f);
            invntt(&mut f);
            for i in 0..N {
                // invntt leaves the result multiplied by 2^16, so dividing by
                // 2^16 is one fqmul against 1.
                let got = to_canonical(barrett_reduce(fqmul(f.0[i], 1)));
                assert_eq!(got as i16, orig.0[i], "coefficient {i}, tag {tag}");
            }
        }
    }

    /// Keys must be byte-identical to the reference, over many seeds and all
    /// three parameter sets. This is what carries the ACVP validation across.
    #[test]
    fn keygen_matches_the_reference_byte_for_byte() {
        for p in sets() {
            for tag in 0..12u8 {
                let ds = seeded(32, tag);
                let zs = seeded(32, tag ^ 0x5a);
                let d: [u8; 32] = ds.try_into().unwrap();
                let z: [u8; 32] = zs.try_into().unwrap();

                let (ek_f, dk_f) = ml_kem_keygen_internal(&p, &d, &z);
                let (ek_s, dk_s) = slow::ml_kem_keygen_internal(&p, &d, &z);
                assert_eq!(ek_f.0, ek_s.0, "{} ek, tag {tag}", p.name);
                assert_eq!(dk_f.0, dk_s.0, "{} dk, tag {tag}", p.name);
            }
        }
    }

    #[test]
    fn encaps_matches_the_reference_byte_for_byte() {
        for p in sets() {
            for tag in 0..12u8 {
                let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                let m: [u8; 32] = seeded(32, tag ^ 0xa5).try_into().unwrap();
                let (ek, _) = ml_kem_keygen_internal(&p, &d, &z);

                let (c_f, k_f) = ml_kem_encaps_internal(&p, &ek, &m).unwrap();
                let (c_s, k_s) = slow::ml_kem_encaps_internal(&p, &ek, &m).unwrap();
                assert_eq!(c_f, c_s, "{} ciphertext, tag {tag}", p.name);
                assert_eq!(k_f, k_s, "{} shared secret, tag {tag}", p.name);
            }
        }
    }

    #[test]
    fn decaps_matches_the_reference_including_implicit_rejection() {
        for p in sets() {
            for tag in 0..8u8 {
                let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                let m: [u8; 32] = seeded(32, tag ^ 0xa5).try_into().unwrap();
                let (ek, dk) = ml_kem_keygen_internal(&p, &d, &z);
                let (c, k) = ml_kem_encaps_internal(&p, &ek, &m).unwrap();

                // Honest ciphertext: both agree, and agree with encapsulation.
                let kf = ml_kem_decaps(&p, &dk, &c).unwrap();
                let ks = slow::ml_kem_decaps(&p, &dk, &c).unwrap();
                assert_eq!(kf, ks, "{} decaps, tag {tag}", p.name);
                assert_eq!(kf, k, "{} round trip, tag {tag}", p.name);

                // Forged ciphertext: the rejection key must match too, which is
                // the part an implementation can get wrong without any test
                // that only checks round trips ever noticing.
                for bit in [0usize, 7, 100, c.len() * 8 - 1] {
                    let mut bad = c.clone();
                    bad[bit / 8] ^= 1 << (bit % 8);
                    let kf = ml_kem_decaps(&p, &dk, &bad).unwrap();
                    let ks = slow::ml_kem_decaps(&p, &dk, &bad).unwrap();
                    assert_eq!(kf, ks, "{} rejection key, tag {tag}, bit {bit}", p.name);
                    assert_ne!(kf, k, "{} rejection differs from real key", p.name);
                }
            }
        }
    }

    #[test]
    fn malformed_lengths_are_rejected() {
        let p = slow::ML_KEM_768;
        let d = [1u8; 32];
        let z = [2u8; 32];
        let (ek, dk) = ml_kem_keygen_internal(&p, &d, &z);

        assert!(ml_kem_decaps(&p, &dk, &vec![0u8; p.ct_len() - 1]).is_none());
        assert!(ml_kem_decaps(&p, &dk, &vec![0u8; p.ct_len() + 1]).is_none());
        let short = MlKemEncapsKey(ek.0[..ek.0.len() - 1].to_vec());
        assert!(ml_kem_encaps_internal(&p, &short, &[0u8; 32]).is_none());
    }

    #[test]
    fn encaps_key_check_matches_the_reference() {
        let p = slow::ML_KEM_768;
        let (ek, _) = ml_kem_keygen_internal(&p, &[3u8; 32], &[4u8; 32]);
        assert!(ml_kem_check_ek(&p, &ek.0));
        assert_eq!(ml_kem_check_ek(&p, &ek.0), slow::ml_kem_check_ek(&p, &ek.0));

        // A non-canonical coefficient: the first 12-bit field set to q itself,
        // which decodes to 0 and so re-encodes differently.  The packing is
        // little-endian by bit, so coefficient 0 is byte 0 plus the LOW nibble
        // of byte 1, and q = 0xd01 splits as 0x01 and 0xd.
        let mut bad = ek.0.clone();
        bad[0] = 0x01;
        bad[1] = (bad[1] & 0xf0) | 0x0d;
        assert!(!ml_kem_check_ek(&p, &bad));
        assert_eq!(ml_kem_check_ek(&p, &bad), slow::ml_kem_check_ek(&p, &bad));

        assert!(!ml_kem_check_ek(&p, &ek.0[..ek.0.len() - 1]));
    }

    /// Randomised end-to-end agreement, without fixed seeds, so the suite is
    /// not only exercising one path through the rejection samplers.
    #[test]
    fn random_round_trips_agree_with_the_reference() {
        for p in sets() {
            for _ in 0..4 {
                let (ek, dk) = ml_kem_keygen(&p);
                let (c, k) = ml_kem_encaps(&p, &ek).unwrap();
                assert_eq!(ml_kem_decaps(&p, &dk, &c).unwrap(), k);
                assert_eq!(slow::ml_kem_decaps(&p, &dk, &c).unwrap(), k);
                assert!(slow::ml_kem_check_ek(&p, &ek.0));
            }
        }
    }

    // ── public-data reuse ────────────────────────────────────────────────────

    /// The cache captured during key generation must be the cache one would
    /// build from the key bytes: same `Â`, `t̂`, `ŝ`, `H(ek)`, `z`. If keygen
    /// left `t̂` or `ŝ` in a non-canonical representative this would still be
    /// right mod q, but two paths to "the same" cache would disagree, so the
    /// test pins them as arrays.
    #[test]
    fn keygen_capture_equals_preparing_from_bytes() {
        for p in sets() {
            for tag in 0..6u8 {
                let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                let (ek, dk, kp) = ml_kem_keygen_prepared_internal(&p, &d, &z);
                let (ek2, dk2) = slow::ml_kem_keygen_internal(&p, &d, &z);
                assert_eq!(ek.0, ek2.0, "{} ek", p.name);
                assert_eq!(dk.0, dk2.0, "{} dk", p.name);

                let fp = MlKemPreparedDecapsKey::new(&p, &dk).unwrap();
                let pe = MlKemPreparedEncapsKey::new(&p, &ek).unwrap();
                for i in 0..p.k {
                    assert_eq!(kp.s_hat[i], fp.s_hat[i], "{} s_hat[{i}]", p.name);
                    assert_eq!(kp.ek.t_hat[i], fp.ek.t_hat[i], "{} t_hat[{i}]", p.name);
                    assert_eq!(kp.ek.t_hat[i], pe.t_hat[i], "{} t_hat[{i}] ek", p.name);
                    for j in 0..p.k {
                        assert_eq!(kp.ek.a[i][j], fp.ek.a[i][j], "{} a[{i}][{j}]", p.name);
                        assert_eq!(kp.ek.a[i][j], pe.a[i][j], "{} a[{i}][{j}] ek", p.name);
                    }
                }
                assert_eq!(kp.ek.h_ek, pe.h_ek);
                assert_eq!(fp.ek.h_ek, pe.h_ek);
                assert_eq!(kp.z, fp.z);
                assert_eq!(kp.ek.ek.0, ek.0);
            }
        }
    }

    /// Prepared encapsulation is byte-identical to the reference.
    #[test]
    fn prepared_encaps_matches_the_reference() {
        for p in sets() {
            for tag in 0..12u8 {
                let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                let (ek, _) = slow::ml_kem_keygen_internal(&p, &d, &z);
                let pe = MlKemPreparedEncapsKey::new(&p, &ek).unwrap();
                for mtag in 0..4u8 {
                    let m: [u8; 32] = seeded(32, tag ^ 0xa5 ^ (mtag << 4)).try_into().unwrap();
                    let (c, k) = pe.encaps_internal(&m);
                    let (c_s, k_s) = slow::ml_kem_encaps_internal(&p, &ek, &m).unwrap();
                    assert_eq!(c, c_s, "{} ciphertext, tag {tag}/{mtag}", p.name);
                    assert_eq!(k, k_s, "{} key, tag {tag}/{mtag}", p.name);
                }
            }
        }
    }

    /// Prepared decapsulation agrees with the reference on honest and forged
    /// ciphertexts, so the cache cannot have changed the rejection path.
    #[test]
    fn prepared_decaps_matches_the_reference_including_rejection() {
        for p in sets() {
            for tag in 0..8u8 {
                let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                let m: [u8; 32] = seeded(32, tag ^ 0xa5).try_into().unwrap();
                let (ek, dk, kp) = ml_kem_keygen_prepared_internal(&p, &d, &z);
                let fp = MlKemPreparedDecapsKey::new(&p, &dk).unwrap();
                let (c, k) = slow::ml_kem_encaps_internal(&p, &ek, &m).unwrap();
                assert_eq!(kp.decaps(&c).unwrap(), k, "{} round trip", p.name);
                assert_eq!(fp.decaps(&c).unwrap(), k, "{} round trip (bytes)", p.name);
                for bit in [0usize, 9, 300, c.len() * 8 - 1] {
                    let mut bad = c.clone();
                    bad[bit / 8] ^= 1 << (bit % 8);
                    let ks = slow::ml_kem_decaps(&p, &dk, &bad).unwrap();
                    assert_eq!(
                        kp.decaps(&bad).unwrap(),
                        ks,
                        "{} rejection bit {bit}",
                        p.name
                    );
                    assert_ne!(ks, k);
                }
                assert!(kp.decaps(&c[..c.len() - 1]).is_none());
            }
        }
    }

    #[test]
    fn prepared_rejects_wrong_lengths() {
        let p = slow::ML_KEM_768;
        let (ek, dk) = ml_kem_keygen_internal(&p, &[1u8; 32], &[2u8; 32]);
        assert!(MlKemPreparedEncapsKey::new(&p, &MlKemEncapsKey(ek.0[1..].to_vec())).is_none());
        assert!(MlKemPreparedDecapsKey::new(&p, &MlKemDecapsKey(dk.0[1..].to_vec())).is_none());
        // A 512 key is the wrong length for 768 parameters.
        let (ek5, _) = ml_kem_keygen_internal(&slow::ML_KEM_512, &[1u8; 32], &[2u8; 32]);
        assert!(MlKemPreparedEncapsKey::new(&p, &ek5).is_none());
    }

    // ── four-way hashing ─────────────────────────────────────────────────────

    /// Run `f` once per Keccak backend the host supports, with that backend
    /// forced on this thread.
    fn for_each_backend(mut f: impl FnMut(keccak4::Backend)) {
        for b in keccak4::tests::backends() {
            keccak4::FORCE.with(|c| c.set(Some(b)));
            f(b);
            keccak4::FORCE.with(|c| c.set(None));
        }
    }

    /// Four-lane `SampleNTT` against the one-lane version over enough seeds
    /// that some lanes need a fourth block (about 0.8% of polynomials do),
    /// which is the branch a fixed three-block implementation gets wrong.
    #[test]
    fn sample_ntt_x4_matches_sample_ntt_including_extra_blocks() {
        let mut extra = 0;
        for tag in 0..=255u8 {
            let rho: [u8; 32] = seeded(32, tag).try_into().unwrap();
            let idx = [(0, 0), (1, 0), (tag % 4, 3), (2, tag % 3)];
            let got = sample_ntt_x4(&rho, idx);
            for l in 0..4 {
                assert_eq!(
                    got[l],
                    sample_ntt(&rho, idx[l].0, idx[l].1),
                    "tag {tag} lane {l}"
                );
                // Count lanes whose first 504 bytes were not enough.
                let mut buf = [0u8; 504];
                let mut inp = rho.to_vec();
                inp.extend_from_slice(&[idx[l].0, idx[l].1]);
                crate::pqc::fast::keccak::shake128_into(&mut buf, &inp);
                let mut r = RejBuf::new();
                r.parse(&buf);
                extra += (!r.done()) as usize;
            }
        }
        assert!(extra > 0, "no seed exercised the extra-block path");
    }

    #[test]
    fn sample_noise_matches_one_call_at_a_time_on_every_backend() {
        let s: [u8; 32] = seeded(32, 77).try_into().unwrap();
        // Mixed η in one group is the ML-KEM-512 encryption shape.
        let shapes: [&[usize]; 5] = [&[3, 3, 2, 2, 2], &[2; 7], &[2; 9], &[3; 4], &[2]];
        for_each_backend(|b| {
            for etas in shapes {
                let mut out = vec![Poly::zero(); etas.len()];
                sample_noise(&s, 5, etas, &mut out);
                for (t, &eta) in etas.iter().enumerate() {
                    assert_eq!(
                        out[t],
                        sample_cbd_prf(&s, 5 + t as u8, eta),
                        "{b:?} {etas:?} {t}"
                    );
                }
            }
        });
    }

    /// The whole scheme on every backend, against the reference.
    #[test]
    fn every_backend_matches_the_reference() {
        for_each_backend(|b| {
            for p in sets() {
                for tag in 0..6u8 {
                    let d: [u8; 32] = seeded(32, tag).try_into().unwrap();
                    let z: [u8; 32] = seeded(32, tag ^ 0x5a).try_into().unwrap();
                    let m: [u8; 32] = seeded(32, tag ^ 0xa5).try_into().unwrap();
                    let (ek, dk) = ml_kem_keygen_internal(&p, &d, &z);
                    let (ek_s, dk_s) = slow::ml_kem_keygen_internal(&p, &d, &z);
                    assert_eq!((&ek.0, &dk.0), (&ek_s.0, &dk_s.0), "{b:?} {} keys", p.name);
                    let (c, k) = ml_kem_encaps_internal(&p, &ek, &m).unwrap();
                    assert_eq!(
                        (c.clone(), k),
                        slow::ml_kem_encaps_internal(&p, &ek, &m).unwrap(),
                        "{b:?} {} encaps",
                        p.name
                    );
                    let pd = MlKemPreparedDecapsKey::new(&p, &dk).unwrap();
                    assert_eq!(pd.encaps_key().encaps_internal(&m), (c.clone(), k));
                    let mut bad = c.clone();
                    bad[3] ^= 0x10;
                    for ct in [&c, &bad] {
                        let want = slow::ml_kem_decaps(&p, &dk, ct).unwrap();
                        assert_eq!(ml_kem_decaps(&p, &dk, ct).unwrap(), want, "{b:?}");
                        assert_eq!(pd.decaps(ct).unwrap(), want, "{b:?} prepared");
                    }
                }
            }
        });
    }

    // ── vector arithmetic ────────────────────────────────────────────────────

    #[cfg(target_arch = "x86_64")]
    fn random_poly(x: &mut u64, bound: i16) -> Poly {
        let mut f = Poly::zero();
        for c in f.0.iter_mut() {
            *x ^= *x << 13;
            *x ^= *x >> 7;
            *x ^= *x << 17;
            *c = ((*x >> 20) % (2 * bound as u64 - 1)) as i16 - (bound - 1);
        }
        f
    }

    /// The AVX2 kernels against the scalar ones, for exact array equality:
    /// the module promises bit-identical output, not just congruence, and
    /// that is what lets the two be mixed.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_kernels_are_bit_identical_to_scalar() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut x = 0xfeed_beef_dead_cafeu64;
        for _ in 0..2000 {
            let f = random_poly(&mut x, Q);
            let (mut a, mut b) = (f, f);
            ntt_scalar(&mut a);
            unsafe { avx2::ntt(&mut b) };
            assert_eq!(a, b, "ntt");

            let (mut a, mut b) = (f, f);
            invntt_scalar(&mut a);
            unsafe { avx2::invntt(&mut b) };
            assert_eq!(a, b, "invntt");

            let g = random_poly(&mut x, Q);
            let acc = random_poly(&mut x, 4 * Q);
            let (mut a, mut b) = (acc, acc);
            poly_basemul_acc_scalar(&mut a, &f, &g);
            unsafe { avx2::basemul_acc(&mut b, &f, &g) };
            assert_eq!(a, b, "basemul_acc");
        }
        // Reduction and the Montgomery scaling over every input up to 9q in
        // absolute value, past the 8q the transforms can produce. Beyond
        // about 9.5q the scalar Barrett step overflows i16 in its product
        // (harmlessly in release, where it wraps, but a debug build panics),
        // so that range is outside both functions' contract.
        for base in (-9 * Q as i32..9 * Q as i32 - N as i32).step_by(N) {
            let mut f = Poly::zero();
            for (i, c) in f.0.iter_mut().enumerate() {
                *c = (base + i as i32) as i16;
            }
            let (mut a, mut b) = (f, f);
            poly_reduce_scalar(&mut a);
            unsafe { avx2::reduce(&mut b) };
            assert_eq!(a, b, "reduce from {base}");
            let (mut a, mut b) = (f, f);
            poly_tomont_scalar(&mut a);
            unsafe { avx2::tomont(&mut b) };
            assert_eq!(a, b, "tomont from {base}");
        }
    }

    /// The grouped packers against `ByteEncode`/`ByteDecode` written bit by
    /// bit from FIPS 203 Algorithms 5 and 6, at every width ML-KEM uses,
    /// including the d = 12 decode of every 12-bit value (the only width
    /// where decoding reduces).
    #[test]
    fn packers_match_bitwise_encode_and_decode() {
        fn encode_bits(f: &[u16; N], d: usize) -> Vec<u8> {
            let mut out = vec![0u8; 32 * d];
            for i in 0..N {
                for j in 0..d {
                    let bit = (f[i] >> j) & 1;
                    out[(i * d + j) / 8] |= (bit as u8) << ((i * d + j) % 8);
                }
            }
            out
        }
        let mut x = 0x1357_9bdf_2468_ace0u64;
        for d in [1usize, 4, 5, 10, 11, 12] {
            for _ in 0..50 {
                let mut vals = [0u16; N];
                let mut f = Poly::zero();
                for i in 0..N {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    let m = if d == 12 { Q as u64 } else { 1u64 << d };
                    vals[i] = (x % m) as u16;
                    f.0[i] = vals[i] as i16;
                }
                let want = encode_bits(&vals, d);
                let mut got = vec![0u8; 32 * d];
                byte_encode_into(&mut got, &f, d);
                assert_eq!(got, want, "encode d={d}");
                assert_eq!(byte_decode(&want, d), f, "decode d={d}");
            }
        }
        for v in 0..4096u16 {
            let mut vals = [0u16; N];
            vals[17] = v;
            let got = byte_decode(&encode_bits(&vals, 12), 12);
            assert_eq!(got.0[17] as u16, v % Q as u16, "d=12 decode of {v}");
        }
    }

    /// Vector compression against scalar, over every representative a
    /// coefficient can arrive in, at every width ML-KEM compresses to.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_compress_matches_scalar_for_every_input() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        for d in [1usize, 4, 5, 10, 11] {
            // Every value in (-q, q), a polynomial at a time.
            let mut start = -(Q as i32) + 1;
            while start < Q as i32 {
                let mut f = Poly::zero();
                for (i, c) in f.0.iter_mut().enumerate() {
                    *c = (start + i as i32).min(Q as i32 - 1) as i16;
                }
                let mut v = [0u16; N];
                unsafe { avx2::compress(&f, d as u32, &mut v) };
                for i in 0..N {
                    let want = compress(to_canonical(barrett_reduce(f.0[i])), d);
                    assert_eq!(v[i], want, "d={d} x={}", f.0[i]);
                }
                start += N as i32;
            }
        }
    }

    /// Vector rejection sampling against the scalar loop, from many starting
    /// counts (including ones that leave fewer than 16 free slots, where the
    /// vector loop must hand over) and on bytes biased so that both accept
    /// and reject are common.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_rejection_sampling_matches_scalar() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let mut x = 0x2545_f491_4f6c_dd1du64;
        for trial in 0..3000 {
            let mut buf = [0u8; 168];
            for b in buf.iter_mut() {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                // Every third trial skews toward 0xff, so most candidates
                // are rejected and the compaction tables see sparse masks.
                *b = if trial % 3 == 0 {
                    (x >> 56) as u8 | 0xc0
                } else {
                    (x >> 56) as u8
                };
            }
            let start = (x as usize) % N;
            let mut a = RejBuf::new();
            a.n = start;
            a.parse(&buf);
            // The same input through the scalar loop alone.
            keccak4::FORCE.with(|f| f.set(Some(keccak4::Backend::Portable)));
            let mut b = RejBuf::new();
            b.n = start;
            b.parse(&buf);
            keccak4::FORCE.with(|f| f.set(None));
            assert_eq!(a.n, b.n, "trial {trial}");
            assert_eq!(a.c[..N.min(a.n)], b.c[..N.min(b.n)], "trial {trial}");
        }
    }
}
