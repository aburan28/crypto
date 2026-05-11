//! **KZG (Kate-Zaverucha-Goldberg 2010) polynomial commitments.**
//!
//! A KZG commitment to a polynomial `p(X) ∈ F_r[X]` is a single
//! group element `C = p(τ)·G₁`, where `τ ∈ F_r` is a **trusted-setup
//! secret** that nobody knows after the setup ceremony.
//!
//! Given `C`, the prover can later open `p(z) = y` to any chosen
//! point `z` using a **single G₁ element** as proof.  Verification
//! is one pairing equation.
//!
//! ## The protocol
//!
//! **Setup** (one-time, requires trusted ceremony):
//!
//! ```text
//! Sample τ ←$ F_r
//! Powers in G₁: [G₁, τ·G₁, τ²·G₁, …, τⁿ·G₁]
//! Powers in G₂: [G₂, τ·G₂]
//! Public:       SRS = (G1_powers, G2_powers); τ is discarded.
//! ```
//!
//! **Commit** to `p(X)`:
//!
//! ```text
//! C = Σᵢ pᵢ · τⁱ·G₁    =   p(τ)·G₁
//! ```
//!
//! **Open** at `z` with claimed value `y = p(z)`:
//!
//! ```text
//! Compute quotient q(X) = (p(X) − y) / (X − z)
//! Proof: π = q(τ)·G₁
//! ```
//!
//! **Verify**:
//!
//! ```text
//! Check  e(C − y·G₁,  G₂)  ==  e(π,  τ·G₂ − z·G₂)
//! ```
//!
//! The equation says `(p(τ) − y) = q(τ) · (τ − z)`, which (by
//! polynomial division uniqueness) implies `p(z) = y`.
//!
//! ## Limitations of this implementation
//!
//! - **Insecure setup** — we generate `τ` in-process for testing.
//!   A real deployment uses MPC ceremonies (Powers of Tau, etc.).
//!   Production code MUST consume an externally-generated SRS.
//! - **Slow** — uses `BigUint` BLS12-381 ops; commits and opens
//!   in seconds, not microseconds.
//! - **Single-point opening only** — batch openings and multi-point
//!   openings (used in PLONK) are not implemented.

use super::polynomial::{fr_mul, fr_neg, fr_sub, Poly};
use crate::bls12_381::fq::scalar_modulus;
use crate::bls12_381::fq2::Fq2;
use crate::bls12_381::fq12::Fq12;
use crate::bls12_381::g1::G1Point;
use crate::bls12_381::g2::G2Point;
use crate::bls12_381::pairing::pairing;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::OsRng;

/// Structured reference string (SRS) for KZG.  Powers of `τ` in
/// both groups.
#[derive(Clone, Debug)]
pub struct KzgSrs {
    /// `[G₁, τ·G₁, …, τⁿ·G₁]`.
    pub g1_powers: Vec<G1Point>,
    /// `[G₂, τ·G₂]`.
    pub g2_powers: Vec<G2Point>,
}

impl KzgSrs {
    /// Maximum polynomial degree this SRS supports.
    pub fn max_degree(&self) -> usize {
        self.g1_powers.len().saturating_sub(1)
    }
}

/// **INSECURE**: generate a fresh `τ` and SRS in-process.  For
/// tests / educational use only.  A production deployment must
/// consume an externally-generated SRS where `τ` was destroyed.
pub fn kzg_insecure_setup(max_degree: usize) -> KzgSrs {
    let mut rng = OsRng;
    let r = scalar_modulus();
    let tau = rng.gen_biguint_below(&r);

    let g1 = G1Point::generator();
    let g2 = G2Point::generator();

    let mut g1_powers = Vec::with_capacity(max_degree + 1);
    let mut tau_pow = BigUint::one();
    for _ in 0..=max_degree {
        g1_powers.push(g1.scalar_mul(&tau_pow));
        tau_pow = (&tau_pow * &tau) % &r;
    }
    let g2_powers = vec![g2.clone(), g2.scalar_mul(&tau)];
    KzgSrs { g1_powers, g2_powers }
}

/// A KZG commitment to a polynomial.
#[derive(Clone, Debug, PartialEq)]
pub struct KzgCommitment {
    pub point: G1Point,
}

/// A KZG opening proof: `q(τ)·G₁`.
#[derive(Clone, Debug, PartialEq)]
pub struct KzgProof {
    pub point: G1Point,
}

/// **Commit** to a polynomial: `C = p(τ)·G₁ = Σᵢ pᵢ · (τⁱ·G₁)`.
pub fn kzg_commit(p: &Poly, srs: &KzgSrs) -> KzgCommitment {
    assert!(
        p.coeffs.len() <= srs.g1_powers.len(),
        "polynomial degree exceeds SRS capacity"
    );
    let mut acc = G1Point::Infinity;
    for (i, c) in p.coeffs.iter().enumerate() {
        if c == &BigUint::from(0u32) { continue; }
        acc = acc.add(&srs.g1_powers[i].scalar_mul(c));
    }
    KzgCommitment { point: acc }
}

/// **Open** the commitment at `z`: returns `(y, π)` with `y = p(z)`
/// and `π = q(τ)·G₁` where `q(X) = (p(X) − y) / (X − z)`.
pub fn kzg_open(p: &Poly, z: &BigUint, srs: &KzgSrs) -> (BigUint, KzgProof) {
    let y = p.evaluate(z);
    // p(X) − y
    let mut shifted = p.clone();
    if shifted.coeffs.is_empty() {
        shifted.coeffs.push(BigUint::from(0u32));
    }
    shifted.coeffs[0] = fr_sub(&shifted.coeffs[0], &y);
    // X − z = (−z, 1)
    let divisor = Poly::from_coeffs(vec![fr_neg(z), BigUint::one()]);
    let (q, rem) = shifted.divmod(&divisor);
    debug_assert!(rem.coeffs.is_empty() || rem.coeffs.iter().all(|c| c.is_zero()));
    let _ = rem;

    // Proof: q(τ)·G₁
    let proof_point = commit_poly_g1(&q, srs);
    (y, KzgProof { point: proof_point })
}

fn commit_poly_g1(p: &Poly, srs: &KzgSrs) -> G1Point {
    let mut acc = G1Point::Infinity;
    for (i, c) in p.coeffs.iter().enumerate() {
        if c.is_zero_fr() { continue; }
        acc = acc.add(&srs.g1_powers[i].scalar_mul(c));
    }
    acc
}

trait FrZero {
    fn is_zero_fr(&self) -> bool;
}
impl FrZero for BigUint {
    fn is_zero_fr(&self) -> bool {
        *self == BigUint::from(0u32)
    }
}

/// **Verify** an opening: check `e(C − y·G₁, G₂) == e(π, τ·G₂ − z·G₂)`.
pub fn kzg_verify(
    commit: &KzgCommitment,
    z: &BigUint,
    y: &BigUint,
    proof: &KzgProof,
    srs: &KzgSrs,
) -> bool {
    let g1 = &srs.g1_powers[0];
    let g2 = &srs.g2_powers[0];
    let tau_g2 = &srs.g2_powers[1];

    // LHS point in G₁: C − y·G₁
    let y_g1 = g1.scalar_mul(y);
    let lhs_g1 = commit.point.add(&y_g1.neg());

    // RHS point in G₂: τ·G₂ − z·G₂
    let z_g2 = g2.scalar_mul(z);
    let rhs_g2 = tau_g2.add(&z_g2.neg());

    // Pairing check.
    let lhs = pairing(&lhs_g1, g2);
    let rhs = pairing(&proof.point, &rhs_g2);
    let _ = Fq2::one(); // suppress unused-import warning if Fq2 wasn't used otherwise
    let _ = Fq12::one();
    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Commitments are additively homomorphic: Com(p) + Com(q) = Com(p + q).
    #[test]
    fn kzg_homomorphic() {
        let srs = kzg_insecure_setup(4);
        let p = Poly::from_coeffs(vec![
            BigUint::from(1u32), BigUint::from(2u32), BigUint::from(3u32),
        ]);
        let q = Poly::from_coeffs(vec![
            BigUint::from(4u32), BigUint::from(5u32),
        ]);
        let c_p = kzg_commit(&p, &srs);
        let c_q = kzg_commit(&q, &srs);
        let c_sum_lhs = c_p.point.add(&c_q.point);
        let c_sum_rhs = kzg_commit(&p.add(&q), &srs).point;
        assert_eq!(c_sum_lhs, c_sum_rhs);
    }

    /// **Completeness**: open & verify roundtrip.
    /// This is the headline test; it exercises the full pairing path
    /// and is slow (~minutes due to BigUint Miller loop).
    #[test]
    #[ignore]
    fn kzg_open_verify_roundtrip() {
        let srs = kzg_insecure_setup(8);
        let p = Poly::from_coeffs(
            (1u32..=5).map(BigUint::from).collect(),
        );
        let z = BigUint::from(3u32);
        let (y, proof) = kzg_open(&p, &z, &srs);
        let commit = kzg_commit(&p, &srs);
        // Sanity: y = p(z) (no pairing).
        assert_eq!(y, p.evaluate(&z));
        // Full verify (uses pairing).
        assert!(kzg_verify(&commit, &z, &y, &proof, &srs));
    }
}
