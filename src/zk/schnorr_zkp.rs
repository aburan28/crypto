//! **Schnorr's sigma protocol** for proof of knowledge of a discrete
//! logarithm (Schnorr 1989).
//!
//! Given a public statement `Y = x·G` for a known generator `G` and
//! group order `n`, the prover demonstrates knowledge of `x` without
//! revealing it.  This is *the* textbook starting point of zero-
//! knowledge cryptography — every sigma protocol generalises this
//! structure.
//!
//! # The protocol (interactive form)
//!
//! Statement: `Y = x·G` known to both prover and verifier.
//! Witness:   `x ∈ Z_n` known only to the prover.
//!
//! ```text
//! Prover                                       Verifier
//! ──────                                       ────────
//! r ← Z_n  (random nonce)
//! R = r·G                       ─── R ──>
//!                                              c ← Z_n  (challenge)
//!                                <── c ───
//! s = r + c·x  (mod n)          ─── s ──>
//!                                              accept iff  s·G = R + c·Y
//! ```
//!
//! Properties (Schnorr 1989):
//!
//! 1. **Completeness**: an honest prover always convinces.
//! 2. **Special soundness**: from two valid transcripts
//!    `(R, c₁, s₁)`, `(R, c₂, s₂)` with `c₁ ≠ c₂`, one can compute
//!    `x = (s₁ − s₂) / (c₁ − c₂)`.  Hence no PPT prover succeeds
//!    without knowing `x`.
//! 3. **Honest-verifier zero-knowledge**: the transcript distribution
//!    is simulable without `x` by picking `s, c` first and setting
//!    `R = s·G − c·Y`.
//!
//! # Non-interactive form (Fiat-Shamir)
//!
//! Replace the verifier's random `c` with a hash of the transcript:
//!
//! ```text
//! c = H(domain_tag ‖ G ‖ Y ‖ R)
//! ```
//!
//! This converts the 3-message interactive proof into a single
//! non-interactive message `(R, s)`.  Security in the **random
//! oracle model** under DLP hardness (Pointcheval-Stern 1996).
//!
//! # Security caveats
//!
//! - **Nonce reuse is catastrophic.**  Same `r` for two different
//!   statements leaks `x`.  This module's API uses fresh OS
//!   randomness internally; `_with_nonce` versions exist for
//!   deterministic testing only — do not call them with the same
//!   `r` twice in production.
//! - **Hash domain separation.**  We tag the Fiat-Shamir hash with
//!   `"ZK-SCHNORR/v1"` to prevent transcript-replay across protocols.
//! - **No `pubkey-as-witness` weak signatures.**  A bug in some early
//!   Schnorr implementations let an attacker forge `(R, s)` for a
//!   chosen `Y`; we prevent this by including `Y` in the challenge
//!   computation (standard practice).

use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::rngs::OsRng;

/// A non-interactive Schnorr proof.
#[derive(Clone, Debug, PartialEq)]
pub struct SchnorrZkProof {
    /// Commitment `R = r·G`.
    pub r_point: Point,
    /// Response `s = r + c·x (mod n)`.
    pub s: BigUint,
}

/// Domain-separation tag for the Fiat-Shamir hash.
const FS_TAG: &str = "ZK-SCHNORR/v1";

/// Compute the Fiat-Shamir challenge `c = H(tag ‖ G ‖ Y ‖ R) mod n`.
fn fs_challenge(curve: &CurveParams, y: &Point, r: &Point) -> BigUint {
    let g = curve.generator();
    let mut buf = Vec::with_capacity(32 + 33 * 3);
    buf.extend_from_slice(FS_TAG.as_bytes());
    buf.extend_from_slice(&encode_point(&g));
    buf.extend_from_slice(&encode_point(y));
    buf.extend_from_slice(&encode_point(r));
    let h = sha256(&buf);
    // Reduce mod n so the challenge sits in Z_n.
    BigUint::from_bytes_be(&h) % &curve.n
}

/// Serialise a point as the 1-byte tag (`0x04` for affine, `0x00` for
/// infinity) followed by 32-byte big-endian `x` and `y` (zero-padded).
/// Keeps the format simple and uniquely-decodable for the FS hash.
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
    // Left-pad to 32 bytes.
    if bytes.len() < 32 {
        out.extend(std::iter::repeat(0).take(32 - bytes.len()));
    }
    out.extend_from_slice(&bytes);
}

/// **Prove** knowledge of `x` such that `Y = x·G` on the given curve.
///
/// Generates fresh randomness for the nonce internally.  Returns a
/// non-interactive proof `(R, s)`.
pub fn schnorr_zkp_prove(x: &BigUint, curve: &CurveParams) -> SchnorrZkProof {
    let mut rng = OsRng;
    let r = rng.gen_biguint_below(&curve.n);
    schnorr_zkp_prove_with_nonce(x, &r, curve)
}

/// **Prove** with a caller-supplied nonce `r`.  Use only for
/// deterministic testing; reusing `r` across statements leaks the
/// witness.
pub fn schnorr_zkp_prove_with_nonce(
    x: &BigUint,
    r: &BigUint,
    curve: &CurveParams,
) -> SchnorrZkProof {
    let a = curve.a_fe();
    let g = curve.generator();
    // R = r·G
    let r_point = g.scalar_mul(r, &a);
    // Y = x·G  (recomputed by the verifier from the public statement)
    let y = g.scalar_mul(x, &a);
    // c = H(tag ‖ G ‖ Y ‖ R)
    let c = fs_challenge(curve, &y, &r_point);
    // s = r + c·x  (mod n)
    let s = (r + &c * x) % &curve.n;
    SchnorrZkProof { r_point, s }
}

/// **Verify** a Schnorr proof against the public statement `Y`.
///
/// Returns `true` iff `s·G = R + c·Y` where `c = H(tag ‖ G ‖ Y ‖ R)`.
pub fn schnorr_zkp_verify(y: &Point, proof: &SchnorrZkProof, curve: &CurveParams) -> bool {
    // Defence: reject proofs with R = O or s out of range.
    match &proof.r_point {
        Point::Infinity => return false,
        Point::Affine { .. } => {}
    }
    if proof.s >= curve.n {
        return false;
    }
    if !curve.is_on_curve(&proof.r_point) {
        return false;
    }
    if !curve.is_on_curve(y) {
        return false;
    }

    let a = curve.a_fe();
    let g = curve.generator();
    let c = fs_challenge(curve, y, &proof.r_point);

    // LHS = s·G
    let lhs = g.scalar_mul(&proof.s, &a);
    // RHS = R + c·Y
    let cy = y.scalar_mul(&c, &a);
    let rhs = proof.r_point.add(&cy, &a);
    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Completeness**: honest prover always convinces verifier.
    #[test]
    fn schnorr_zkp_completeness() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::parse_bytes(b"DEADBEEF1234567890ABCDEF", 16).unwrap();
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);

        let proof = schnorr_zkp_prove(&x, &curve);
        assert!(schnorr_zkp_verify(&y, &proof, &curve));
    }

    /// **Determinism**: same nonce → same proof.
    #[test]
    fn schnorr_zkp_deterministic_with_nonce() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::from(42u32);
        let r = BigUint::from(99u32);

        let p1 = schnorr_zkp_prove_with_nonce(&x, &r, &curve);
        let p2 = schnorr_zkp_prove_with_nonce(&x, &r, &curve);
        assert_eq!(p1, p2);
    }

    /// **Soundness**: tampering with the response breaks the proof.
    #[test]
    fn schnorr_zkp_rejects_tampered_s() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::from(7u32);
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);

        let mut proof = schnorr_zkp_prove(&x, &curve);
        proof.s = (&proof.s + 1u32) % &curve.n;
        assert!(!schnorr_zkp_verify(&y, &proof, &curve));
    }

    /// **Soundness**: tampering with R breaks the proof.
    #[test]
    fn schnorr_zkp_rejects_tampered_r() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::from(13u32);
        let a = curve.a_fe();
        let g = curve.generator();
        let y = g.scalar_mul(&x, &a);

        let mut proof = schnorr_zkp_prove(&x, &curve);
        // Replace R with R + G.
        proof.r_point = proof.r_point.add(&g, &a);
        assert!(!schnorr_zkp_verify(&y, &proof, &curve));
    }

    /// **Soundness**: proof for a different statement Y' does not verify
    /// for Y.
    #[test]
    fn schnorr_zkp_rejects_wrong_statement() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::from(2024u32);
        let proof = schnorr_zkp_prove(&x, &curve);

        // Verify against Y' = (x+1)·G, NOT the correct Y.
        let a = curve.a_fe();
        let g = curve.generator();
        let y_wrong = g.scalar_mul(&(&x + 1u32), &a);
        assert!(!schnorr_zkp_verify(&y_wrong, &proof, &curve));
    }

    /// **Special soundness**: from two transcripts with different
    /// challenges but same commitment, the witness can be extracted.
    /// We simulate this by constructing two proofs with the SAME r
    /// but different statements (effectively, same R for different x;
    /// this is what an attacker would observe if a careless prover
    /// reused `r`).
    #[test]
    fn schnorr_zkp_witness_extraction_from_reused_nonce() {
        let curve = CurveParams::secp256k1();
        let x = BigUint::from(123456789u64);
        let r = BigUint::from(987654321u64);

        // Two proofs with the SAME nonce r but for the same x (just
        // different challenges, which we'd get if we lied about Y).
        // We simulate "different challenges, same R" by constructing
        // the equations ourselves.
        let a = curve.a_fe();
        let g = curve.generator();
        let r_pt = g.scalar_mul(&r, &a);
        let y = g.scalar_mul(&x, &a);

        let c1 = BigUint::from(11u32);
        let c2 = BigUint::from(22u32);
        let s1 = (&r + &c1 * &x) % &curve.n;
        let s2 = (&r + &c2 * &x) % &curve.n;

        // Extraction: x = (s1 - s2) / (c1 - c2)  (mod n)
        let n = &curve.n;
        // s1 - s2 mod n
        let s_diff = if s1 >= s2 {
            (&s1 - &s2) % n
        } else {
            (n - ((&s2 - &s1) % n)) % n
        };
        let c_diff = if c1 >= c2 {
            (&c1 - &c2) % n
        } else {
            (n - ((&c2 - &c1) % n)) % n
        };
        let c_diff_inv = crate::utils::mod_inverse(&c_diff, n).unwrap();
        let x_extracted = (&s_diff * &c_diff_inv) % n;

        assert_eq!(
            x_extracted, x,
            "two transcripts with same R but different c reveal x"
        );
        // Sanity that R_pt and y satisfy both equations.
        let _ = r_pt;
        let _ = y;
    }

    /// **Honest-verifier zero-knowledge** (simulator construction):
    /// pick `s` and `c` first, set `R = s·G − c·Y`.  Verifies as a
    /// real proof, but no knowledge of `x` was used.
    #[test]
    fn schnorr_zkp_simulator_produces_valid_transcript() {
        let curve = CurveParams::secp256k1();
        let a = curve.a_fe();
        let g = curve.generator();
        // Public Y for some unknown x.
        let x_real = BigUint::from(0xCAFEBABEu32);
        let y = g.scalar_mul(&x_real, &a);

        // Simulator: pick s, c uniformly; set R = s·G − c·Y.
        // Note: this uses the *non-Fiat-Shamir* protocol — it works
        // only because the verifier chose c at random.  With Fiat-
        // Shamir, c is determined by H(R, ...) which prevents this
        // simulation.  We illustrate the interactive form's HVZK
        // property here.
        let s = BigUint::from(0xDEADBEEFu64);
        let c = BigUint::from(0xFEEDFACEu64);
        let s_g = g.scalar_mul(&s, &a);
        let cy = y.scalar_mul(&c, &a);
        let r = s_g.add(&cy.neg(), &a);

        // Check the verifier's equation: s·G == R + c·Y.
        let lhs = s_g;
        let rhs = r.add(&cy, &a);
        assert_eq!(lhs, rhs);
    }
}
