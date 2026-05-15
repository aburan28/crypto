//! **Chaum-Pedersen** equality of discrete logarithms (Chaum-Pedersen
//! 1992).
//!
//! Given two pairs `(G, Y)` and `(H, Z)` with `Y = x·G` and `Z = x·H`
//! for the **same** witness `x`, prove knowledge of `x` and equality
//! of the two discrete logs **without revealing `x`**.
//!
//! This composes the Schnorr sigma protocol with itself, sharing the
//! nonce `r` across the two equations.  Used widely:
//!
//! - **ElGamal re-encryption proofs** (mix-nets, voting).
//! - **Ring signatures / OR proofs** (Camenisch-Stadler 1997).
//! - **Threshold signature schemes** (FROST, Pedersen DKG).
//! - **VDF-based randomness beacons** (where one proves the same
//!   exponent was applied to two bases).
//!
//! # The protocol
//!
//! Statement: `(G, Y, H, Z)` with `Y = x·G` and `Z = x·H` (claimed).
//! Witness:   `x ∈ Z_n` known only to prover.
//!
//! ```text
//! Prover                                       Verifier
//! ──────                                       ────────
//! r ← Z_n
//! R₁ = r·G,  R₂ = r·H            ─── R₁, R₂ ──>
//!                                              c ← Z_n
//!                               <── c ───
//! s = r + c·x  (mod n)          ─── s ──>
//!                                              accept iff
//!                                                s·G = R₁ + c·Y
//!                                                s·H = R₂ + c·Z
//! ```
//!
//! Non-interactive (Fiat-Shamir): `c = H(tag ‖ G ‖ Y ‖ H ‖ Z ‖ R₁ ‖ R₂)`.
//!
//! # Properties
//!
//! - **Completeness**: trivial.
//! - **Soundness**: from two transcripts `(R₁, R₂, c, s)` and
//!   `(R₁, R₂, c', s')` with `c ≠ c'`, extract
//!   `x = (s − s')/(c − c')`.  Witness consistency across both
//!   equations is enforced by the shared `r`.
//! - **HVZK**: same simulator as Schnorr but emitting `(R₁, R₂, c, s)`.
//!
//! # Security caveat
//!
//! Sharing the nonce `r` across the two equations is essential to
//! the protocol — but it **must not** be reused across different
//! Chaum-Pedersen proofs, just like Schnorr.

use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::rngs::OsRng;

/// A Chaum-Pedersen proof of equality of two discrete logs.
#[derive(Clone, Debug, PartialEq)]
pub struct ChaumPedersenProof {
    /// `R₁ = r·G`.
    pub r1: Point,
    /// `R₂ = r·H`.
    pub r2: Point,
    /// `s = r + c·x (mod n)`.
    pub s: BigUint,
}

const FS_TAG: &str = "ZK-CHAUM-PEDERSEN/v1";

/// Compute Fiat-Shamir challenge `c = H(tag ‖ G ‖ Y ‖ H ‖ Z ‖ R₁ ‖ R₂) mod n`.
fn fs_challenge(
    curve: &CurveParams,
    h_pt: &Point,
    y: &Point,
    z: &Point,
    r1: &Point,
    r2: &Point,
) -> BigUint {
    let g = curve.generator();
    for counter in 0u32.. {
        let mut buf = Vec::with_capacity(FS_TAG.len() + 65 * 6 + 4);
        buf.extend_from_slice(FS_TAG.as_bytes());
        buf.extend_from_slice(&counter.to_be_bytes());
        buf.extend_from_slice(&encode_point(&g));
        buf.extend_from_slice(&encode_point(y));
        buf.extend_from_slice(&encode_point(h_pt));
        buf.extend_from_slice(&encode_point(z));
        buf.extend_from_slice(&encode_point(r1));
        buf.extend_from_slice(&encode_point(r2));
        let mut c = BigUint::from_bytes_be(&sha256(&buf));
        let q_bits = curve.n.bits() as usize;
        if q_bits < 256 {
            c >>= 256 - q_bits;
        }
        if !c.is_zero() && c < curve.n {
            return c;
        }
    }
    unreachable!("u32 counter space is enough for Fiat-Shamir rejection sampling")
}

fn encode_point(p: &Point) -> Vec<u8> {
    match p {
        Point::Infinity => vec![0x00; 65],
        Point::Affine { x, y } => {
            let mut out = Vec::with_capacity(65);
            out.push(0x04);
            push_be32(&mut out, &x.value);
            push_be32(&mut out, &y.value);
            out
        }
    }
}

fn push_be32(out: &mut Vec<u8>, v: &BigUint) {
    let bytes = v.to_bytes_be();
    if bytes.len() < 32 {
        out.extend(std::iter::repeat(0).take(32 - bytes.len()));
    }
    out.extend_from_slice(&bytes);
}

/// **Prove** that `log_G(Y) = log_H(Z) = x`, with fresh OS randomness.
pub fn chaum_pedersen_prove(x: &BigUint, h_pt: &Point, curve: &CurveParams) -> ChaumPedersenProof {
    let mut rng = OsRng;
    let r = rng.gen_biguint_below(&curve.n);
    chaum_pedersen_prove_with_nonce(x, &r, h_pt, curve)
}

/// **Prove** with caller-supplied nonce.  Testing only.
pub fn chaum_pedersen_prove_with_nonce(
    x: &BigUint,
    r: &BigUint,
    h_pt: &Point,
    curve: &CurveParams,
) -> ChaumPedersenProof {
    let a = curve.a_fe();
    let g = curve.generator();
    let r1 = g.scalar_mul(r, &a);
    let r2 = h_pt.scalar_mul(r, &a);
    let y = g.scalar_mul(x, &a);
    let z = h_pt.scalar_mul(x, &a);
    let c = fs_challenge(curve, h_pt, &y, &z, &r1, &r2);
    let s = (r + &c * x) % &curve.n;
    ChaumPedersenProof { r1, r2, s }
}

/// **Verify** a Chaum-Pedersen proof.  Checks both equations
/// `s·G = R₁ + c·Y` and `s·H = R₂ + c·Z`.
pub fn chaum_pedersen_verify(
    y: &Point,
    z: &Point,
    h_pt: &Point,
    proof: &ChaumPedersenProof,
    curve: &CurveParams,
) -> bool {
    // Reject identity commitments.
    match (&proof.r1, &proof.r2) {
        (Point::Infinity, _) | (_, Point::Infinity) => return false,
        _ => {}
    }
    if proof.s >= curve.n {
        return false;
    }
    if !curve.is_valid_public_point(y)
        || !curve.is_valid_public_point(z)
        || !curve.is_valid_public_point(h_pt)
    {
        return false;
    }
    if !curve.is_valid_public_point(&proof.r1) || !curve.is_valid_public_point(&proof.r2) {
        return false;
    }

    let a = curve.a_fe();
    let g = curve.generator();
    let c = fs_challenge(curve, h_pt, y, z, &proof.r1, &proof.r2);

    // First check: s·G = R₁ + c·Y.
    let s_g = g.scalar_mul(&proof.s, &a);
    let cy = y.scalar_mul(&c, &a);
    let r1_plus_cy = proof.r1.add(&cy, &a);
    if s_g != r1_plus_cy {
        return false;
    }

    // Second check: s·H = R₂ + c·Z.
    let s_h = h_pt.scalar_mul(&proof.s, &a);
    let cz = z.scalar_mul(&c, &a);
    let r2_plus_cz = proof.r2.add(&cz, &a);
    s_h == r2_plus_cz
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zk::pedersen::pedersen_second_generator;

    /// **Completeness**: honest prover convinces verifier.
    #[test]
    fn chaum_pedersen_completeness() {
        let curve = CurveParams::secp256k1();
        let h_pt = pedersen_second_generator(&curve);
        let x = BigUint::from(0xCAFEBABEu64);
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);
        let z = h_pt.scalar_mul(&x, &a);

        let proof = chaum_pedersen_prove(&x, &h_pt, &curve);
        assert!(chaum_pedersen_verify(&y, &z, &h_pt, &proof, &curve));
    }

    /// **Soundness**: tampering with `s` fails.
    #[test]
    fn chaum_pedersen_rejects_tampered_s() {
        let curve = CurveParams::secp256k1();
        let h_pt = pedersen_second_generator(&curve);
        let x = BigUint::from(42u32);
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);
        let z = h_pt.scalar_mul(&x, &a);

        let mut proof = chaum_pedersen_prove(&x, &h_pt, &curve);
        proof.s = (&proof.s + 1u32) % &curve.n;
        assert!(!chaum_pedersen_verify(&y, &z, &h_pt, &proof, &curve));
    }

    /// **Soundness**: mismatched DLs are rejected.
    /// `Y = x·G` and `Z = x'·H` with `x ≠ x'` should not have a
    /// valid Chaum-Pedersen proof (no prover can construct one).
    #[test]
    fn chaum_pedersen_rejects_mismatched_dls() {
        let curve = CurveParams::secp256k1();
        let h_pt = pedersen_second_generator(&curve);
        let x = BigUint::from(11u32);
        let x_other = BigUint::from(22u32);
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);
        let z_wrong = h_pt.scalar_mul(&x_other, &a); // Z uses DIFFERENT x

        // Prover *attempts* to prove equality using x (the DL of Y).
        let proof = chaum_pedersen_prove(&x, &h_pt, &curve);
        assert!(!chaum_pedersen_verify(&y, &z_wrong, &h_pt, &proof, &curve));
    }

    /// **Determinism with fixed nonce**.
    #[test]
    fn chaum_pedersen_deterministic_with_nonce() {
        let curve = CurveParams::secp256k1();
        let h_pt = pedersen_second_generator(&curve);
        let x = BigUint::from(99u32);
        let r = BigUint::from(7u32);
        let p1 = chaum_pedersen_prove_with_nonce(&x, &r, &h_pt, &curve);
        let p2 = chaum_pedersen_prove_with_nonce(&x, &r, &h_pt, &curve);
        assert_eq!(p1, p2);
    }
}
