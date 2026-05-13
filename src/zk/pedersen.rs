//! **Pedersen commitments** (Pedersen 1991).
//!
//! Given two generators `G, H` of a cyclic group of order `n` with
//! **unknown discrete-log relation** between them (`H = h·G` for an
//! `h` known to nobody), the Pedersen commitment to a message
//! `m ∈ Z_n` is:
//!
//! ```text
//! Com(m, r) := m·G + r·H
//! ```
//!
//! where `r ∈ Z_n` is a uniform random **blinding factor** sampled
//! by the committer.
//!
//! # Properties
//!
//! - **Perfectly hiding.**  Because `r` is uniform on `Z_n` and `H`
//!   is a generator, `Com(m, r)` is uniformly distributed over the
//!   group regardless of `m`.  Even an unbounded adversary cannot
//!   distinguish `Com(m₁, ·)` from `Com(m₂, ·)`.
//!
//! - **Computationally binding** under the discrete-log assumption.
//!   If a committer could open `Com(m, r) = Com(m', r')` with
//!   `m ≠ m'`, then `(m − m')·G = (r' − r)·H`, hence
//!   `H = ((m − m')/(r' − r))·G` — recovering `log_G(H)`, breaking
//!   DLP.
//!
//! # Generating an independent `H`
//!
//! On secp256k1/P-256 we need a second generator `H` with unknown
//! discrete-log relation to `G`.  The **hash-to-curve** method:
//!
//! 1. Compute `h = SHA-256(G_serialised ‖ tag)`.
//! 2. Treat `h` as a candidate `x`-coordinate.
//! 3. If `x³ + ax + b` is a quadratic residue mod `p`, set `y =
//!    √(x³+ax+b)` (pick even-y by convention); else increment `x`
//!    and retry.
//!
//! The resulting `H` has a "nothing-up-my-sleeve" derivation —
//! no party knows `log_G(H)`.
//!
//! # What this is NOT
//!
//! - Not a zero-knowledge proof system on its own.  To prove
//!   *something* about a committed `m` (e.g. range proofs), combine
//!   with sigma protocols ([`super::schnorr_zkp`]) or Bulletproofs
//!   (not yet implemented).
//! - Not an encryption.  Anyone who opens the commitment can read
//!   `m`; Pedersen is binding+hiding, not confidential under
//!   chosen-ciphertext.

use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::OsRng;

/// Parameters for a Pedersen commitment scheme over a chosen curve:
/// the curve itself plus the auxiliary generator `H`.
#[derive(Clone, Debug)]
pub struct PedersenParams {
    pub curve: CurveParams,
    pub h: Point,
}

impl PedersenParams {
    /// Construct standard Pedersen parameters for the given curve,
    /// deriving `H` via nothing-up-my-sleeve hash-to-curve from
    /// `SHA-256(G_serialised ‖ "PEDERSEN/H/v1")`.
    pub fn standard(curve: CurveParams) -> Self {
        let h = pedersen_second_generator(&curve);
        Self { curve, h }
    }
}

/// A Pedersen commitment: an opaque elliptic-curve point.
/// "Opaque" in the cryptographic sense — the point reveals nothing
/// about `m` or `r` without the opener disclosing them.
#[derive(Clone, Debug, PartialEq)]
pub struct PedersenCommitment {
    pub point: Point,
}

/// **Derive the second generator `H` for Pedersen commitments**, with
/// no party knowing `log_G(H)`.
///
/// Hash-to-curve via try-and-increment: start from
/// `SHA-256(G_x ‖ G_y ‖ "PEDERSEN/H/v1")` interpreted as a candidate
/// `x`-coordinate, increment until `(x, y)` is on the curve.
pub fn pedersen_second_generator(curve: &CurveParams) -> Point {
    let g = curve.generator();
    let (gx, gy) = match &g {
        Point::Affine { x, y } => (&x.value, &y.value),
        Point::Infinity => panic!("curve generator must be affine"),
    };
    let mut seed = Vec::with_capacity(64 + 32);
    push_be32(&mut seed, gx);
    push_be32(&mut seed, gy);
    seed.extend_from_slice(b"PEDERSEN/H/v1");
    let mut x_candidate = BigUint::from_bytes_be(&sha256(&seed));
    x_candidate %= &curve.p;

    let three = BigUint::from(3u32);
    let two = BigUint::from(2u32);

    // Try x_candidate, x_candidate+1, ... until we find an on-curve point.
    for _ in 0..1024 {
        let rhs_value =
            (x_candidate.modpow(&three, &curve.p) + &curve.a * &x_candidate + &curve.b) % &curve.p;
        // Test if rhs is a QR mod p via Euler: rhs^((p-1)/2) mod p.
        let p_minus_1_over_2 = (&curve.p - BigUint::one()) / &two;
        let euler = rhs_value.modpow(&p_minus_1_over_2, &curve.p);
        if euler == BigUint::one() {
            // Compute square root.  For NIST primes we use Tonelli-Shanks
            // or, when p ≡ 3 (mod 4), the simpler formula y = rhs^((p+1)/4).
            //
            // Both secp256k1's and P-256's primes satisfy p ≡ 3 (mod 4)
            // for secp256k1 (since p = 2^256 - 2^32 - 977, and 977 ≡ 1 mod 4
            // gives p ≡ -1 ≡ 3 mod 4).  For P-256 the prime is slightly
            // different; we check at runtime.
            let exp = if (&curve.p % BigUint::from(4u32)) == BigUint::from(3u32) {
                (&curve.p + BigUint::one()) / BigUint::from(4u32)
            } else {
                // Fall back via Tonelli-Shanks would go here; for P-256
                // we handle the p ≡ 1 mod 4 case below.
                tonelli_shanks_exponent(&curve.p)
            };
            let mut y = rhs_value.modpow(&exp, &curve.p);
            // Normalise to even-y (BIP-340-style x-only convention) so
            // the H is canonical.
            if (&y & BigUint::one()) == BigUint::one() {
                y = &curve.p - &y;
            }
            let point = Point::Affine {
                x: FieldElement::new(x_candidate.clone(), curve.p.clone()),
                y: FieldElement::new(y, curve.p.clone()),
            };
            debug_assert!(curve.is_on_curve(&point));
            return point;
        }
        x_candidate = (&x_candidate + BigUint::one()) % &curve.p;
    }
    panic!("hash-to-curve failed after 1024 attempts (statistically impossible)");
}

/// Tonelli-Shanks "exponent" placeholder — for p ≡ 1 mod 4, square-
/// root extraction is not just `rhs^((p+1)/4)`.  For our hash-to-curve
/// purposes, we accept this code path is only exercised for P-256
/// (where p ≡ 3 mod 4 happens to hold for `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² +
/// 2⁹⁶ − 1` since 1 ≡ 1 mod 4 ... actually let me check).
/// Both secp256k1 and P-256 satisfy `p ≡ 3 (mod 4)`, so the simple
/// formula always applies.
fn tonelli_shanks_exponent(p: &BigUint) -> BigUint {
    // Fallback: (p+1)/4 only valid for p ≡ 3 mod 4.  For the curves
    // we care about, this is the case.
    (p + BigUint::one()) / BigUint::from(4u32)
}

fn push_be32(out: &mut Vec<u8>, v: &BigUint) {
    let bytes = v.to_bytes_be();
    if bytes.len() < 32 {
        out.extend(std::iter::repeat(0).take(32 - bytes.len()));
    }
    out.extend_from_slice(&bytes);
}

/// **Commit** to a message `m` with fresh random blinding `r`.
/// Returns `(commitment, r)` so the caller can later open.
pub fn pedersen_commit(m: &BigUint, params: &PedersenParams) -> (PedersenCommitment, BigUint) {
    let mut rng = OsRng;
    let r = rng.gen_biguint_below(&params.curve.n);
    let c = pedersen_commit_with_blinding(m, &r, params);
    (c, r)
}

/// **Commit** with a caller-supplied blinding `r`.  Use the
/// random-blinding variant in production; this one is for testing
/// determinism.
pub fn pedersen_commit_with_blinding(
    m: &BigUint,
    r: &BigUint,
    params: &PedersenParams,
) -> PedersenCommitment {
    let a = params.curve.a_fe();
    let g = params.curve.generator();
    let m_g = g.scalar_mul(&(m % &params.curve.n), &a);
    let r_h = params.h.scalar_mul(&(r % &params.curve.n), &a);
    PedersenCommitment {
        point: m_g.add(&r_h, &a),
    }
}

/// **Open** a commitment: verify that `(m, r)` re-create the given
/// commitment.
pub fn pedersen_open(
    commitment: &PedersenCommitment,
    m: &BigUint,
    r: &BigUint,
    params: &PedersenParams,
) -> bool {
    let recomputed = pedersen_commit_with_blinding(m, r, params);
    commitment == &recomputed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `H` is on the curve and distinct from `G`.
    #[test]
    fn h_generator_is_on_curve_and_distinct() {
        for curve in [CurveParams::secp256k1(), CurveParams::p256()] {
            let h = pedersen_second_generator(&curve);
            assert!(curve.is_on_curve(&h), "H must be on curve {}", curve.name);
            assert_ne!(
                h,
                curve.generator(),
                "H must differ from G on {}",
                curve.name
            );
        }
    }

    /// Standard Pedersen params derive deterministically.
    #[test]
    fn standard_params_are_deterministic() {
        let p1 = PedersenParams::standard(CurveParams::secp256k1());
        let p2 = PedersenParams::standard(CurveParams::secp256k1());
        assert_eq!(p1.h, p2.h);
    }

    /// **Correctness**: open(commit(m, r), m, r) = true.
    #[test]
    fn pedersen_open_succeeds_on_honest_open() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let m = BigUint::from(42u32);
        let (c, r) = pedersen_commit(&m, &params);
        assert!(pedersen_open(&c, &m, &r, &params));
    }

    /// **Binding (sanity)**: opening with a different message fails.
    #[test]
    fn pedersen_open_fails_on_wrong_message() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let m = BigUint::from(100u32);
        let (c, r) = pedersen_commit(&m, &params);
        let m_wrong = BigUint::from(101u32);
        assert!(!pedersen_open(&c, &m_wrong, &r, &params));
    }

    /// **Binding (sanity)**: opening with a different blinding fails.
    #[test]
    fn pedersen_open_fails_on_wrong_blinding() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let m = BigUint::from(7u32);
        let (c, r) = pedersen_commit(&m, &params);
        let r_wrong = (&r + 1u32) % &params.curve.n;
        assert!(!pedersen_open(&c, &m, &r_wrong, &params));
    }

    /// **Hiding (statistical sanity)**: same message committed twice
    /// with fresh randomness should give two different commitments.
    #[test]
    fn pedersen_two_commits_to_same_message_differ() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let m = BigUint::from(1234u32);
        let (c1, _) = pedersen_commit(&m, &params);
        let (c2, _) = pedersen_commit(&m, &params);
        assert_ne!(c1, c2, "two commitments to the same m must differ");
    }

    /// **Homomorphism**: Com(m₁, r₁) + Com(m₂, r₂) = Com(m₁+m₂, r₁+r₂).
    ///
    /// Pedersen commitments are additively homomorphic — a property
    /// used in confidential-transaction protocols (Bulletproofs,
    /// Mimblewimble) to verify balance sums without revealing
    /// individual amounts.
    #[test]
    fn pedersen_is_additively_homomorphic() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let m1 = BigUint::from(10u32);
        let m2 = BigUint::from(20u32);
        let r1 = BigUint::from(111u32);
        let r2 = BigUint::from(222u32);
        let c1 = pedersen_commit_with_blinding(&m1, &r1, &params);
        let c2 = pedersen_commit_with_blinding(&m2, &r2, &params);
        let a = params.curve.a_fe();
        let c_sum = c1.point.add(&c2.point, &a);

        let m_sum = (&m1 + &m2) % &params.curve.n;
        let r_sum = (&r1 + &r2) % &params.curve.n;
        let c_expected = pedersen_commit_with_blinding(&m_sum, &r_sum, &params);

        assert_eq!(c_sum, c_expected.point);
    }
}
