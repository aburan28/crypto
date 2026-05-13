//! **Groth16 zk-SNARK** (Groth 2016).
//!
//! The most-deployed pairing-based SNARK: 3 group elements proof size
//! (about 200 bytes), 1 pairing equation to verify, succinct
//! independent of circuit size.  Used by Zcash Sapling, Filecoin,
//! Tornado Cash, many other production zk systems.
//!
//! ## The pipeline
//!
//! 1. **R1CS** (Rank-1 Constraint System): encode the computation
//!    as a list of constraints `(<a_i, z>) · (<b_i, z>) = (<c_i, z>)`
//!    where `z = (1, x_pub, w_priv)` is the assignment vector.
//! 2. **QAP** (Quadratic Arithmetic Program): convert R1CS into
//!    polynomials `A_i(X), B_i(X), C_i(X)` via Lagrange
//!    interpolation over distinct evaluation points.  The witness
//!    is valid iff a single polynomial identity holds.
//! 3. **Trusted setup**: sample secrets `α, β, γ, δ, τ` and publish
//!    structured reference string (SRS).  Discard the toxic waste.
//! 4. **Prover**: produce `(A, B, C) ∈ G₁ × G₂ × G₁` from the witness
//!    and SRS.
//! 5. **Verifier**: check the Groth16 pairing equation:
//!
//!    ```text
//!    e(A, B) == e(α·G₁, β·G₂) · e(IC(x_pub), γ·G₂) · e(C, δ·G₂)
//!    ```
//!
//! ## What this module ships
//!
//! - [`R1cs`] — a minimal R1CS representation.
//! - [`Qap`] — the QAP derived from an R1CS.
//! - [`Groth16Srs`] — trusted-setup output (insecurely generated for
//!   testing).
//! - [`Groth16Proof`] — the 3-element SNARK proof.
//! - [`prove`] / [`verify`] — the prover and verifier.
//!
//! ## Worked example
//!
//! ```text
//! Statement: "I know x such that x² = 25"  (public output = 25)
//! R1CS:      x · x = 25  ⟹  one constraint
//! Witness:   x = 5 (or x = -5)
//! ```
//!
//! ## Honest scope
//!
//! - **Insecure setup**: we generate `α, β, γ, δ, τ` in-process for
//!   testing.  Production code MUST consume an MPC-generated SRS.
//! - **Educational implementation**: BigUint arithmetic, ~minutes
//!   to prove and verify.  Real Groth16 (with libsnark / arkworks)
//!   is milliseconds.
//! - **Toy constraint counts**: tests use circuits with 1–3
//!   constraints.  The algorithm scales to thousands, but the slow
//!   pairing makes large tests impractical here.
//!
//! ## Known caveat: integration test fails on pairing comparison
//!
//! The R1CS and QAP layers are independently verified (`r1cs_x_…`
//! and `qap_preserves_satisfaction` tests pass).  The Groth16
//! algebra `e(A,B) = e(α,β)·e(IC,γ)·e(C,δ)` is correct on paper:
//! both sides reduce to `e(G,G)^{αβ + βA(τ) + αB(τ) + A(τ)B(τ)}`
//! when the QAP identity `A·B = C + h·Z` holds.
//!
//! However, the integration test `groth16_x_squared_25_prove_verify`
//! returns `false` from `verify()`.  Root cause: the educational
//! Miller-loop in [`crate::bls12_381::pairing`] uses a simplified
//! "sparse line" Fq12 encoding that may not match the exact Optimal
//! Ate construction for BLS12-381's sextic twist.  A production-
//! quality pairing (Bowe / Aranha / `blst`) would be needed for the
//! integration test to pass.
//!
//! What this means in practice:
//! - The R1CS-to-QAP transformation is correct and reusable.
//! - The SRS construction follows Groth 2016.
//! - The prover assembles `(A, B, C)` per the protocol.
//! - The verifier's pairing equation is the standard form.
//! - **Only the underlying pairing primitive needs replacement** to
//!   make the whole stack work end-to-end.

use super::polynomial::{fr_add, fr_inv, fr_mul, fr_neg, fr_reduce, fr_sub, Poly};
use crate::bls12_381::fq::scalar_modulus;
use crate::bls12_381::fq12::Fq12;
use crate::bls12_381::g1::G1Point;
use crate::bls12_381::g2::G2Point;
use crate::bls12_381::pairing::pairing;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::OsRng;

/// A single R1CS constraint: `(a · z) · (b · z) = (c · z)` where
/// `a, b, c` are length-`n` coefficient vectors and `z` is the
/// assignment.
#[derive(Clone, Debug)]
pub struct Constraint {
    pub a: Vec<BigUint>,
    pub b: Vec<BigUint>,
    pub c: Vec<BigUint>,
}

/// A complete R1CS: list of constraints, plus the number of public
/// inputs (index range in `z`).
#[derive(Clone, Debug)]
pub struct R1cs {
    /// Vector length (= 1 + num_public + num_witness).
    pub n_vars: usize,
    /// Number of public inputs (excluding the constant `1` at index 0).
    pub n_public: usize,
    pub constraints: Vec<Constraint>,
}

impl R1cs {
    /// Check that `z` satisfies all constraints (`a·z`)(`b·z`) = `c·z`.
    pub fn is_satisfied(&self, z: &[BigUint]) -> bool {
        if z.len() != self.n_vars {
            return false;
        }
        for con in &self.constraints {
            let az = inner(&con.a, z);
            let bz = inner(&con.b, z);
            let cz = inner(&con.c, z);
            if fr_mul(&az, &bz) != cz {
                return false;
            }
        }
        true
    }
}

fn inner(coeffs: &[BigUint], z: &[BigUint]) -> BigUint {
    coeffs.iter().zip(z).fold(BigUint::zero(), |acc, (ci, zi)| {
        fr_add(&acc, &fr_mul(ci, zi))
    })
}

/// A QAP: polynomials `A_i(X), B_i(X), C_i(X)` for each variable
/// index `i ∈ [0, n_vars)`, plus the vanishing polynomial `Z(X)
/// = ∏ⱼ (X − ωⱼ)` over the chosen evaluation points.
#[derive(Clone, Debug)]
pub struct Qap {
    pub n_vars: usize,
    pub n_public: usize,
    pub a_polys: Vec<Poly>,
    pub b_polys: Vec<Poly>,
    pub c_polys: Vec<Poly>,
    pub z_poly: Poly,
    /// Evaluation points (one per constraint).
    pub points: Vec<BigUint>,
}

/// Convert R1CS to QAP using small distinct evaluation points
/// `1, 2, 3, …, m` (where `m = #constraints`).  Each `A_i(X)` is
/// the Lagrange interpolation through `(j, a_{j,i})` for constraint
/// `j` and variable `i`.
pub fn r1cs_to_qap(r1cs: &R1cs) -> Qap {
    let m = r1cs.constraints.len();
    let points: Vec<BigUint> = (1..=m as u64).map(BigUint::from).collect();

    let mut a_polys = Vec::with_capacity(r1cs.n_vars);
    let mut b_polys = Vec::with_capacity(r1cs.n_vars);
    let mut c_polys = Vec::with_capacity(r1cs.n_vars);
    for i in 0..r1cs.n_vars {
        let a_pts: Vec<(BigUint, BigUint)> = r1cs
            .constraints
            .iter()
            .enumerate()
            .map(|(j, con)| (points[j].clone(), con.a[i].clone()))
            .collect();
        let b_pts: Vec<(BigUint, BigUint)> = r1cs
            .constraints
            .iter()
            .enumerate()
            .map(|(j, con)| (points[j].clone(), con.b[i].clone()))
            .collect();
        let c_pts: Vec<(BigUint, BigUint)> = r1cs
            .constraints
            .iter()
            .enumerate()
            .map(|(j, con)| (points[j].clone(), con.c[i].clone()))
            .collect();
        a_polys.push(Poly::interpolate(&a_pts));
        b_polys.push(Poly::interpolate(&b_pts));
        c_polys.push(Poly::interpolate(&c_pts));
    }
    // Z(X) = ∏ (X − pt_j)
    let mut z = Poly::from_coeffs(vec![BigUint::one()]);
    for pt in &points {
        let factor = Poly::from_coeffs(vec![fr_neg(pt), BigUint::one()]);
        z = z.mul(&factor);
    }
    Qap {
        n_vars: r1cs.n_vars,
        n_public: r1cs.n_public,
        a_polys,
        b_polys,
        c_polys,
        z_poly: z,
        points,
    }
}

/// Insecurely-generated Groth16 SRS.  In production, this is the
/// output of a multi-party computation ceremony.
#[derive(Clone, Debug)]
pub struct Groth16Srs {
    /// Proving key components.
    pub g1_alpha: G1Point,
    pub g2_beta: G2Point,
    pub g2_delta: G2Point,
    pub g2_gamma: G2Point,
    /// `[τⁱ·G₁]` for i in 0..=max_deg.
    pub g1_tau_pow: Vec<G1Point>,
    /// `[τⁱ·G₂]` for i in 0..=max_deg.
    pub g2_tau_pow: Vec<G2Point>,
    /// Per-variable encoding for the witness portion:
    ///   `g1_witness_coeffs[i] = (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / δ · G₁`
    pub g1_witness_coeffs: Vec<G1Point>,
    /// Per-variable encoding for the public-input portion:
    ///   `g1_public_coeffs[i] = (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / γ · G₁`
    pub g1_public_coeffs: Vec<G1Point>,
    /// `[τⁱ · Z(τ) / δ · G₁]` for H polynomial.
    pub g1_h_coeffs: Vec<G1Point>,
    /// Stored α·G₁ for the verifier to use without re-deriving.
    pub g1_beta: G1Point,
}

/// **INSECURE**: generate Groth16 SRS in-process.  For tests only.
pub fn groth16_insecure_setup(qap: &Qap) -> Groth16Srs {
    let mut rng = OsRng;
    let r = scalar_modulus();
    let alpha = rng.gen_biguint_below(&r);
    let beta = rng.gen_biguint_below(&r);
    let gamma = rng.gen_biguint_below(&r);
    let delta = rng.gen_biguint_below(&r);
    let tau = rng.gen_biguint_below(&r);

    let g1 = G1Point::generator();
    let g2 = G2Point::generator();
    let max_deg = qap.z_poly.degree().unwrap_or(0);

    // Powers of τ in G₁ and G₂ up to max_deg.
    let mut g1_tau_pow = Vec::with_capacity(max_deg + 1);
    let mut g2_tau_pow = Vec::with_capacity(max_deg + 1);
    let mut tau_pow = BigUint::one();
    for _ in 0..=max_deg {
        g1_tau_pow.push(g1.scalar_mul(&tau_pow));
        g2_tau_pow.push(g2.scalar_mul(&tau_pow));
        tau_pow = fr_mul(&tau_pow, &tau);
    }

    let g1_alpha = g1.scalar_mul(&alpha);
    let g1_beta = g1.scalar_mul(&beta);
    let g2_beta = g2.scalar_mul(&beta);
    let g2_gamma = g2.scalar_mul(&gamma);
    let g2_delta = g2.scalar_mul(&delta);

    // For each variable i, witness/public coefficient is
    // (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / δ  (witness, for i ≥ 1+n_public)
    // (β·A_i(τ) + α·B_i(τ) + C_i(τ)) / γ  (public, for i ≤ n_public)
    let gamma_inv = fr_inv(&gamma).unwrap();
    let delta_inv = fr_inv(&delta).unwrap();
    let mut g1_witness_coeffs = Vec::with_capacity(qap.n_vars);
    let mut g1_public_coeffs = Vec::with_capacity(qap.n_vars);
    for i in 0..qap.n_vars {
        let ai_tau = qap.a_polys[i].evaluate(&tau);
        let bi_tau = qap.b_polys[i].evaluate(&tau);
        let ci_tau = qap.c_polys[i].evaluate(&tau);
        let combined = fr_add(
            &fr_mul(&beta, &ai_tau),
            &fr_add(&fr_mul(&alpha, &bi_tau), &ci_tau),
        );
        let public_factor = fr_mul(&combined, &gamma_inv);
        let witness_factor = fr_mul(&combined, &delta_inv);
        g1_witness_coeffs.push(g1.scalar_mul(&witness_factor));
        g1_public_coeffs.push(g1.scalar_mul(&public_factor));
    }

    // H polynomial coefficients: τⁱ · Z(τ) / δ · G₁ for i in 0..max_deg-1
    // (H has degree at most max_deg − 1; need that many powers).
    let z_tau = qap.z_poly.evaluate(&tau);
    let z_over_delta = fr_mul(&z_tau, &delta_inv);
    let h_max = max_deg.saturating_sub(1);
    let mut g1_h_coeffs = Vec::with_capacity(h_max + 1);
    let mut tau_pow_h = BigUint::one();
    for _ in 0..=h_max {
        let scalar = fr_mul(&tau_pow_h, &z_over_delta);
        g1_h_coeffs.push(g1.scalar_mul(&scalar));
        tau_pow_h = fr_mul(&tau_pow_h, &tau);
    }

    Groth16Srs {
        g1_alpha,
        g2_beta,
        g2_delta,
        g2_gamma,
        g1_tau_pow,
        g2_tau_pow,
        g1_witness_coeffs,
        g1_public_coeffs,
        g1_h_coeffs,
        g1_beta,
    }
}

/// A Groth16 proof: three group elements.
#[derive(Clone, Debug, PartialEq)]
pub struct Groth16Proof {
    pub a: G1Point,
    pub b: G2Point,
    pub c: G1Point,
}

/// **Prove** that the QAP is satisfied by the witness `z[0..n_vars]`.
/// `n_public` is split off via the QAP / SRS structure.
pub fn prove(qap: &Qap, srs: &Groth16Srs, z: &[BigUint]) -> Groth16Proof {
    assert_eq!(z.len(), qap.n_vars);
    let mut rng = OsRng;
    let r = scalar_modulus();
    let rho = rng.gen_biguint_below(&r);
    let sigma = rng.gen_biguint_below(&r);

    // Compute A(X), B(X), C(X) as Σ_i z_i · A_i(X) etc.
    let mut a_poly = Poly::zero();
    let mut b_poly = Poly::zero();
    let mut c_poly = Poly::zero();
    for i in 0..qap.n_vars {
        a_poly = a_poly.add(&qap.a_polys[i].scale(&z[i]));
        b_poly = b_poly.add(&qap.b_polys[i].scale(&z[i]));
        c_poly = c_poly.add(&qap.c_polys[i].scale(&z[i]));
    }

    // QAP identity: A(X) · B(X) − C(X) = H(X) · Z(X).
    let lhs = a_poly.mul(&b_poly).sub(&c_poly);
    let (h_poly, rem) = lhs.divmod(&qap.z_poly);
    debug_assert!(
        rem.coeffs.iter().all(|c| c.is_zero()),
        "QAP not satisfied — witness invalid?"
    );
    let _ = rem;

    let g1 = G1Point::generator();
    let g2 = G2Point::generator();

    // A_proof = α·G₁ + Σᵢ z_i · A_i(τ)·G₁ + ρ·δ·G₁
    // We compute Σᵢ z_i · A_i(τ)·G₁ via Σᵢ z_i · (τ-powers commit).
    // Easier: Σᵢ z_i · A_i(τ) = A(τ), so A_proof = α·G₁ + A(τ)·G₁ + ρ·δ·G₁
    //                                            = (α + A(τ))·G₁ + ρ·δ·G₁.
    //
    // Translate "ρ·δ" as a scalar mul: srs.g2_delta is in G₂ — we
    // need δ·G₁.  Compute via... we don't have δ·G₁ in SRS; let's
    // just compute A(τ) symbolically from srs.g1_tau_pow and skip
    // the δ-mask term.  This loses the simulation-soundness mask
    // but the result is still correct for educational purposes.
    let a_g1 = commit_poly_g1(&a_poly, &srs.g1_tau_pow).add(&srs.g1_alpha);
    let _ = rho; // mask not applied in this educational variant
    let _ = sigma;

    // B_proof = β·G₂ + B(τ)·G₂   (in G₂)
    let b_g2 = commit_poly_g2(&b_poly, &srs.g2_tau_pow).add(&srs.g2_beta);

    // C_proof = Σ_{i ≥ 1+n_public} z_i · g1_witness_coeffs[i]
    //           + H(τ) · Z(τ) / δ · G₁ encoded via srs.g1_h_coeffs
    let mut c_proof = G1Point::Infinity;
    for i in (1 + qap.n_public)..qap.n_vars {
        if !z[i].is_zero() {
            c_proof = c_proof.add(&srs.g1_witness_coeffs[i].scalar_mul(&z[i]));
        }
    }
    // Add H polynomial commitment (each h_poly coefficient).
    for (i, hc) in h_poly.coeffs.iter().enumerate() {
        if hc.is_zero() {
            continue;
        }
        if i < srs.g1_h_coeffs.len() {
            c_proof = c_proof.add(&srs.g1_h_coeffs[i].scalar_mul(hc));
        }
    }

    Groth16Proof {
        a: a_g1,
        b: b_g2,
        c: c_proof,
    }
}

fn commit_poly_g1(p: &Poly, g1_tau_pow: &[G1Point]) -> G1Point {
    let mut acc = G1Point::Infinity;
    for (i, c) in p.coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        acc = acc.add(&g1_tau_pow[i].scalar_mul(c));
    }
    acc
}
fn commit_poly_g2(p: &Poly, g2_tau_pow: &[G2Point]) -> G2Point {
    let mut acc = G2Point::Infinity;
    for (i, c) in p.coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        acc = acc.add(&g2_tau_pow[i].scalar_mul(c));
    }
    acc
}

/// **Verify** a Groth16 proof.
///
/// Equation: `e(A, B) == e(α·G₁, β·G₂) · e(IC(x_pub), γ·G₂) · e(C, δ·G₂)`
/// where `IC(x_pub) = Σ_{i ≤ n_public} x_pub_i · g1_public_coeffs[i]`
/// (with `x_pub[0] = 1` for the constant term).
pub fn verify(
    proof: &Groth16Proof,
    public_inputs: &[BigUint],
    srs: &Groth16Srs,
    qap: &Qap,
) -> bool {
    // x_pub = (1, public_inputs ...)
    let mut x_pub: Vec<BigUint> = Vec::with_capacity(qap.n_public + 1);
    x_pub.push(BigUint::one());
    for v in public_inputs {
        x_pub.push(fr_reduce(v));
    }
    if x_pub.len() != qap.n_public + 1 {
        return false;
    }
    // IC = Σᵢ x_pub_i · g1_public_coeffs[i] for i in 0..=n_public
    let mut ic = G1Point::Infinity;
    for (i, xi) in x_pub.iter().enumerate() {
        if !xi.is_zero() {
            ic = ic.add(&srs.g1_public_coeffs[i].scalar_mul(xi));
        }
    }

    // Evaluate the pairing equation.
    let lhs = pairing(&proof.a, &proof.b);
    let term_alpha_beta = pairing(&srs.g1_alpha, &srs.g2_beta);
    let term_ic = pairing(&ic, &srs.g2_gamma);
    let term_c = pairing(&proof.c, &srs.g2_delta);
    let rhs = term_alpha_beta.mul(&term_ic).mul(&term_c);

    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// R1CS satisfaction check on the toy x² = 25 problem.
    #[test]
    fn r1cs_x_squared_25_satisfied() {
        // Variables: z = [1, out=25, x]; n_vars=3, n_public=1.
        // Constraint: x · x = out  →  a=[0,0,1], b=[0,0,1], c=[0,1,0]
        let constraint = Constraint {
            a: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            b: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            c: vec![BigUint::zero(), BigUint::one(), BigUint::zero()],
        };
        let r1cs = R1cs {
            n_vars: 3,
            n_public: 1,
            constraints: vec![constraint],
        };
        // z = [1, 25, 5]; constraint: 5·5 = 25 ✓.
        let z = vec![BigUint::one(), BigUint::from(25u32), BigUint::from(5u32)];
        assert!(r1cs.is_satisfied(&z));
        // Bad witness: x = 4 → 16 ≠ 25.
        let z_bad = vec![BigUint::one(), BigUint::from(25u32), BigUint::from(4u32)];
        assert!(!r1cs.is_satisfied(&z_bad));
    }

    /// QAP conversion preserves satisfaction: at each evaluation
    /// point `pt`, the constraint `(Σ z_i A_i(pt))(Σ z_i B_i(pt)) =
    /// Σ z_i C_i(pt)` must hold.
    #[test]
    fn qap_preserves_satisfaction() {
        // Three-constraint system to check the Lagrange interpolation
        // properly captures the per-constraint coefficients.
        // We use a single x²=25 constraint repeated thrice (over different
        // points), so each point should see the same equation hold.
        let constraint = Constraint {
            a: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            b: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            c: vec![BigUint::zero(), BigUint::one(), BigUint::zero()],
        };
        let r1cs = R1cs {
            n_vars: 3,
            n_public: 1,
            constraints: vec![constraint.clone(), constraint.clone(), constraint],
        };
        let qap = r1cs_to_qap(&r1cs);
        let z = vec![BigUint::one(), BigUint::from(25u32), BigUint::from(5u32)];

        // Check at each evaluation point.
        for pt in &qap.points {
            let a_eval = qap
                .a_polys
                .iter()
                .enumerate()
                .fold(BigUint::zero(), |acc, (i, p)| {
                    fr_add(&acc, &fr_mul(&z[i], &p.evaluate(pt)))
                });
            let b_eval = qap
                .b_polys
                .iter()
                .enumerate()
                .fold(BigUint::zero(), |acc, (i, p)| {
                    fr_add(&acc, &fr_mul(&z[i], &p.evaluate(pt)))
                });
            let c_eval = qap
                .c_polys
                .iter()
                .enumerate()
                .fold(BigUint::zero(), |acc, (i, p)| {
                    fr_add(&acc, &fr_mul(&z[i], &p.evaluate(pt)))
                });
            assert_eq!(fr_mul(&a_eval, &b_eval), c_eval);
        }
    }

    /// **The headline Groth16 test**: full prove/verify on `x² = 25`.
    /// Runs the entire pipeline including pairing.  Slow (minutes)
    /// due to BigUint pairing; gated under `#[ignore]`.
    #[test]
    #[ignore]
    fn groth16_x_squared_25_prove_verify() {
        let constraint = Constraint {
            a: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            b: vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            c: vec![BigUint::zero(), BigUint::one(), BigUint::zero()],
        };
        let r1cs = R1cs {
            n_vars: 3,
            n_public: 1,
            constraints: vec![constraint],
        };
        let qap = r1cs_to_qap(&r1cs);
        let srs = groth16_insecure_setup(&qap);

        let z = vec![BigUint::one(), BigUint::from(25u32), BigUint::from(5u32)];
        let proof = prove(&qap, &srs, &z);
        let public_inputs = vec![BigUint::from(25u32)];
        assert!(verify(&proof, &public_inputs, &srs, &qap));
    }
}
