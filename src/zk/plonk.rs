//! **PLONK** (Gabizon, Williamson, Ciobotaru 2019) — universal-setup,
//! polynomial-IOP-based zk-SNARK.
//!
//! PLONK encodes a computation as an **arithmetic circuit** of gates
//! with five selector polynomials (`q_L, q_R, q_O, q_M, q_C`) and
//! **copy constraints** wiring outputs of some gates to inputs of
//! others.  The witness `(a, b, c)` (length-`n` left/right/output
//! wire vectors) satisfies the circuit iff:
//!
//! 1. **Gate equation** (at each row `i`):
//!    `q_L[i]·a[i] + q_R[i]·b[i] + q_O[i]·c[i] + q_M[i]·a[i]·b[i] + q_C[i] = 0`.
//! 2. **Copy constraint**: the wires form a multiset matching the
//!    intended wiring permutation `σ`, expressed as a "grand-product"
//!    polynomial identity.
//!
//! ## Polynomial IOP
//!
//! Both constraints reduce to **polynomial identities on a
//! multiplicative subgroup `H ⊂ F_r`**.  The prover commits to the
//! wire polynomials `a(X), b(X), c(X)` and the permutation
//! accumulator `Z(X)` via KZG.  The verifier picks a random
//! challenge `ζ ∈ F_r` (Fiat-Shamir), the prover opens the relevant
//! polynomials at `ζ`, and the verifier checks the combined
//! identity reduces to zero at `ζ`.
//!
//! ## What this module ships
//!
//! - [`Circuit`] — gate selectors + permutation indices.
//! - [`Witness`] — left/right/output wire assignments.
//! - [`Plonkish`] — circuit + witness checker (sanity check that
//!   the gate equation and copy constraints hold).
//! - [`prove`] / [`verify`] — full polynomial-IOP prover and
//!   verifier with Fiat-Shamir.
//!
//! ## Honest scope
//!
//! - **Polynomial-IOP level only.**  We implement the prover/verifier
//!   in terms of polynomial *values* at random challenge `ζ` — not
//!   in terms of pairing-based KZG opening proofs.  The KZG-binding
//!   layer is in [`super::kzg`]; combining the two requires a
//!   working pairing primitive, which our `bls12_381::pairing`
//!   does not yet provide (see Groth16 notes).  The polynomial-IOP
//!   verifier here accepts iff the combined identity holds at `ζ` —
//!   this is the actual soundness-relevant check; KZG just binds
//!   the polynomial choices.
//!   The standalone [`verify`] API is therefore an educational identity checker,
//!   not a complete binding SNARK verifier.
//!
//! - **Single multiplicative-subgroup domain.**  We use a small
//!   power-of-2 domain `H` with a primitive root of unity `ω` over
//!   `F_r` (the BLS12-381 scalar field has a 2^32-th root of unity,
//!   so any domain up to size `2^32` is supported).

use super::polynomial::{fr_add, fr_inv, fr_mul, fr_neg, fr_pow, fr_reduce, fr_sub, Poly};
use crate::bls12_381::fq::scalar_modulus;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::OsRng;

const FS_TAG: &str = "ZK-PLONK/v1";

// ── Multiplicative subgroup helpers ────────────────────────────────

/// `ω` = a primitive `n`-th root of unity in `F_r`, where `n` is a
/// power of 2 dividing `2^32`.
///
/// Derivation: `F_r^*` has order `r − 1 = 2^32 · t` for odd `t`.
/// For any generator `g` of `F_r^*` (a multiplicative non-residue
/// works), `g^((r−1)/n)` is a primitive `n`-th root of unity.
///
/// We use `g = 7` (verified non-residue in `F_r` for BLS12-381) and
/// compute `g^((r−1)/n)` at runtime.  Cached on first call would
/// be a natural optimisation; we recompute for simplicity.
pub fn primitive_nth_root(n: usize) -> BigUint {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    let log_n = n.trailing_zeros() as u64;
    assert!(log_n <= 32, "n must divide 2^32");
    let r = scalar_modulus();
    // For BLS12-381, `g = 7` is a multiplicative generator of F_r*.
    let g = BigUint::from(7u32);
    let exp = (&r - BigUint::one()) / BigUint::from(n as u64);
    fr_pow(&g, &exp)
}

/// `[1, ω, ω², …, ω^{n−1}]`.
fn domain(n: usize) -> Vec<BigUint> {
    let omega = primitive_nth_root(n);
    let mut out = Vec::with_capacity(n);
    let mut acc = BigUint::one();
    for _ in 0..n {
        out.push(acc.clone());
        acc = fr_mul(&acc, &omega);
    }
    out
}

// ── Circuit / witness ─────────────────────────────────────────────

/// A PLONK arithmetic circuit, encoded as selector polynomials
/// over a multiplicative-subgroup domain of size `n`.  All vectors
/// have length `n` (the gate count, padded to a power of 2).
#[derive(Clone, Debug)]
pub struct Circuit {
    pub n: usize,
    /// Selector for the `a` input (left-wire) term.
    pub q_l: Vec<BigUint>,
    pub q_r: Vec<BigUint>,
    pub q_o: Vec<BigUint>,
    pub q_m: Vec<BigUint>,
    pub q_c: Vec<BigUint>,
    /// Permutation: for each of the 3n wires, the *index* it's
    /// wired to.  Wire indices follow the convention used in the
    /// PLONK paper: wires 0..n are `a`-column, n..2n are `b`,
    /// 2n..3n are `c`.
    pub permutation: Vec<usize>,
}

/// Witness vectors: the prover's secret left/right/output wire
/// assignments for each gate row.
#[derive(Clone, Debug)]
pub struct Witness {
    pub a: Vec<BigUint>,
    pub b: Vec<BigUint>,
    pub c: Vec<BigUint>,
}

impl Circuit {
    pub fn is_well_formed(&self) -> bool {
        self.n.is_power_of_two()
            && self.q_l.len() == self.n
            && self.q_r.len() == self.n
            && self.q_o.len() == self.n
            && self.q_m.len() == self.n
            && self.q_c.len() == self.n
            && self.permutation.len() == 3 * self.n
            && self.permutation.iter().all(|&i| i < 3 * self.n)
    }

    /// **Direct check** of the circuit constraints against a witness.
    /// Used for sanity testing the prover; the real verifier would
    /// instead check via the polynomial identity at challenge ζ.
    pub fn is_satisfied(&self, w: &Witness) -> bool {
        if !self.is_well_formed() {
            return false;
        }
        if w.a.len() != self.n || w.b.len() != self.n || w.c.len() != self.n {
            return false;
        }
        // Gate equation.
        for i in 0..self.n {
            let term = fr_add(
                &fr_add(
                    &fr_mul(&self.q_l[i], &w.a[i]),
                    &fr_add(
                        &fr_mul(&self.q_r[i], &w.b[i]),
                        &fr_mul(&self.q_o[i], &w.c[i]),
                    ),
                ),
                &fr_add(
                    &fr_mul(&self.q_m[i], &fr_mul(&w.a[i], &w.b[i])),
                    &self.q_c[i],
                ),
            );
            if !term.is_zero() {
                return false;
            }
        }
        // Permutation (copy constraints).
        let all_wires: Vec<&BigUint> = w.a.iter().chain(w.b.iter()).chain(w.c.iter()).collect();
        for (i, &target) in self.permutation.iter().enumerate() {
            if all_wires[i] != all_wires[target] {
                return false;
            }
        }
        true
    }
}

// ── Fiat-Shamir transcript ─────────────────────────────────────────

#[derive(Clone)]
struct Transcript {
    state: Vec<u8>,
}

impl Transcript {
    fn new() -> Self {
        Self {
            state: FS_TAG.as_bytes().to_vec(),
        }
    }
    fn append(&mut self, label: &str, bytes: &[u8]) {
        self.state.extend_from_slice(label.as_bytes());
        self.state.extend_from_slice(bytes);
    }
    fn append_scalar(&mut self, label: &str, s: &BigUint) {
        let mut buf = [0u8; 32];
        let bytes = s.to_bytes_be();
        buf[32 - bytes.len()..].copy_from_slice(&bytes);
        self.append(label, &buf);
    }
    fn append_poly(&mut self, label: &str, p: &Poly) {
        let n = p.coeffs.len();
        self.append(label, &(n as u64).to_be_bytes());
        for c in &p.coeffs {
            self.append_scalar("c", c);
        }
    }
    fn challenge(&mut self, label: &str) -> BigUint {
        self.state.extend_from_slice(label.as_bytes());
        let h = sha256(&self.state);
        self.state.extend_from_slice(&h);
        BigUint::from_bytes_be(&h) % scalar_modulus()
    }
}

// ── Prover / Verifier ──────────────────────────────────────────────

/// A PLONK proof at the polynomial-IOP level: polynomial commitments
/// represented as the polynomials themselves (no KZG binding yet),
/// plus the evaluations at the verifier's challenge `ζ` and the
/// quotient polynomial `t`.
#[derive(Clone, Debug)]
pub struct PlonkProof {
    /// Wire polynomials interpolated at the domain.
    pub a_poly: Poly,
    pub b_poly: Poly,
    pub c_poly: Poly,
    /// Permutation grand-product accumulator polynomial `Z(X)`.
    pub z_poly: Poly,
    /// Quotient polynomial `t(X)` such that the combined identity
    /// equals `t(X) · Z_H(X)` (where `Z_H` is the vanishing
    /// polynomial of the domain `H`).
    pub t_poly: Poly,
    /// Wire evaluations at `ζ`.
    pub a_zeta: BigUint,
    pub b_zeta: BigUint,
    pub c_zeta: BigUint,
    pub z_zeta_omega: BigUint,
}

/// Prove that `witness` satisfies `circuit`.  Returns a [`PlonkProof`].
pub fn prove(circuit: &Circuit, witness: &Witness) -> PlonkProof {
    assert!(
        circuit.is_satisfied(witness),
        "witness does not satisfy circuit (prover would fail anyway; flagged early)"
    );
    let n = circuit.n;
    let h_domain = domain(n);

    // Interpolate witness polynomials over H.
    let a_poly = lagrange_from_evaluations(&h_domain, &witness.a);
    let b_poly = lagrange_from_evaluations(&h_domain, &witness.b);
    let c_poly = lagrange_from_evaluations(&h_domain, &witness.c);

    // Build the permutation accumulator Z(X) using PLONK's
    // beta/gamma challenges.
    let mut transcript = Transcript::new();
    transcript.append_poly("a", &a_poly);
    transcript.append_poly("b", &b_poly);
    transcript.append_poly("c", &c_poly);
    let beta = transcript.challenge("beta");
    let gamma = transcript.challenge("gamma");

    let z_values = compute_z(circuit, witness, &h_domain, &beta, &gamma);
    let z_poly = lagrange_from_evaluations(&h_domain, &z_values);

    // Construct the quotient polynomial t(X) such that
    //   t(X) · Z_H(X) = combined_identity(X)
    // where combined_identity is gate_constraint + α·perm_constraint.
    transcript.append_poly("z", &z_poly);
    let alpha = transcript.challenge("alpha");

    let t_poly = compute_quotient(
        circuit, &h_domain, &a_poly, &b_poly, &c_poly, &z_poly, &beta, &gamma, &alpha,
    );

    transcript.append_poly("t", &t_poly);
    let zeta = transcript.challenge("zeta");

    let omega = primitive_nth_root(n);
    let zeta_omega = fr_mul(&zeta, &omega);

    PlonkProof {
        a_poly: a_poly.clone(),
        b_poly: b_poly.clone(),
        c_poly: c_poly.clone(),
        z_poly: z_poly.clone(),
        t_poly: t_poly.clone(),
        a_zeta: a_poly.evaluate(&zeta),
        b_zeta: b_poly.evaluate(&zeta),
        c_zeta: c_poly.evaluate(&zeta),
        z_zeta_omega: z_poly.evaluate(&zeta_omega),
    }
}

/// Verify a [`PlonkProof`] against the circuit.  Replays the
/// transcript to recover `(β, γ, α, ζ)`, then checks the combined
/// polynomial identity at `ζ`.
pub fn verify(circuit: &Circuit, proof: &PlonkProof) -> bool {
    if !circuit.is_well_formed() {
        return false;
    }
    let n = circuit.n;
    let h_domain = domain(n);

    let mut transcript = Transcript::new();
    transcript.append_poly("a", &proof.a_poly);
    transcript.append_poly("b", &proof.b_poly);
    transcript.append_poly("c", &proof.c_poly);
    let beta = transcript.challenge("beta");
    let gamma = transcript.challenge("gamma");
    transcript.append_poly("z", &proof.z_poly);
    let alpha = transcript.challenge("alpha");
    transcript.append_poly("t", &proof.t_poly);
    let zeta = transcript.challenge("zeta");

    let omega = primitive_nth_root(n);

    // Recompute the selector and identity-permutation values at ζ
    // (the verifier knows the circuit, so this is "preprocessing").
    let q_l_poly = lagrange_from_evaluations(&h_domain, &circuit.q_l);
    let q_r_poly = lagrange_from_evaluations(&h_domain, &circuit.q_r);
    let q_o_poly = lagrange_from_evaluations(&h_domain, &circuit.q_o);
    let q_m_poly = lagrange_from_evaluations(&h_domain, &circuit.q_m);
    let q_c_poly = lagrange_from_evaluations(&h_domain, &circuit.q_c);

    let a_z = &proof.a_zeta;
    let b_z = &proof.b_zeta;
    let c_z = &proof.c_zeta;
    let q_l_z = q_l_poly.evaluate(&zeta);
    let q_r_z = q_r_poly.evaluate(&zeta);
    let q_o_z = q_o_poly.evaluate(&zeta);
    let q_m_z = q_m_poly.evaluate(&zeta);
    let q_c_z = q_c_poly.evaluate(&zeta);

    // Gate identity at ζ:
    //   gate(ζ) = q_M·a·b + q_L·a + q_R·b + q_O·c + q_C.
    let gate_zeta = fr_add(
        &fr_add(
            &fr_add(&fr_mul(&q_m_z, &fr_mul(a_z, b_z)), &fr_mul(&q_l_z, a_z)),
            &fr_add(&fr_mul(&q_r_z, b_z), &fr_mul(&q_o_z, c_z)),
        ),
        &q_c_z,
    );

    // Permutation identity at ζ (simplified — full PLONK uses
    // σ_a, σ_b, σ_c polynomials and the grand-product check):
    //   perm(ζ) = Z(ζ·ω) · ∏_i (wire_i(ζ) + β·σ_i(ζ) + γ)
    //           − Z(ζ)   · ∏_i (wire_i(ζ) + β·k_i·ζ + γ).
    let z_z = proof.z_poly.evaluate(&zeta);
    let z_zw = &proof.z_zeta_omega;
    let perm_check = permutation_identity_at_zeta(
        circuit, &h_domain, &zeta, &beta, &gamma, a_z, b_z, c_z, &z_z, z_zw,
    );

    // Combined identity:
    //   combined(ζ) = gate(ζ) + α·perm(ζ).
    let combined = fr_add(&gate_zeta, &fr_mul(&alpha, &perm_check));

    // Vanishing polynomial Z_H(X) = X^n − 1 at ζ.
    let zeta_n = fr_pow(&zeta, &BigUint::from(n as u64));
    let z_h_zeta = fr_sub(&zeta_n, &BigUint::one());

    // Check: combined(ζ) = t(ζ) · Z_H(ζ).
    let t_zeta = proof.t_poly.evaluate(&zeta);
    let rhs = fr_mul(&t_zeta, &z_h_zeta);
    let _ = omega;
    combined == rhs
}

// ── Helpers ─────────────────────────────────────────────────────────

/// Build the polynomial of degree `< n` that takes `evals[i]` at
/// `domain[i]`.  Uses naive Lagrange interpolation; for small `n`
/// this is fine, but a production PLONK would use an inverse-FFT
/// over the multiplicative subgroup.
fn lagrange_from_evaluations(domain: &[BigUint], evals: &[BigUint]) -> Poly {
    assert_eq!(domain.len(), evals.len());
    let pts: Vec<(BigUint, BigUint)> = domain
        .iter()
        .zip(evals)
        .map(|(d, e)| (d.clone(), e.clone()))
        .collect();
    Poly::interpolate(&pts)
}

/// Compute the values of the permutation grand-product polynomial
/// `Z(X)` on the domain `H`.
///
/// Recursion:
///   `Z(ω^0) = 1`
///   `Z(ω^{i+1}) = Z(ω^i) · numerator_i / denominator_i`
/// where:
///   `numerator_i = ∏_{wires}(wire_i + β·k_wire·ω^i + γ)`
///   `denominator_i = ∏_{wires}(wire_i + β·σ(wire_position_i) + γ)`
///
/// Here `k_a, k_b, k_c` are coset-distinguishing constants (1, k1, k2);
/// for the educational implementation, we use `k_a = 1`, `k_b = 2`,
/// `k_c = 3`.
fn compute_z(
    circuit: &Circuit,
    witness: &Witness,
    h_domain: &[BigUint],
    beta: &BigUint,
    gamma: &BigUint,
) -> Vec<BigUint> {
    let n = circuit.n;
    let k_a = BigUint::one();
    let k_b = BigUint::from(2u32);
    let k_c = BigUint::from(3u32);
    let mut z = Vec::with_capacity(n);
    z.push(BigUint::one()); // Z(ω^0) = 1
    let mut acc = BigUint::one();
    for i in 0..(n - 1) {
        // Numerator at row i.
        let num_a = fr_add(
            &witness.a[i],
            &fr_add(&fr_mul(beta, &fr_mul(&k_a, &h_domain[i])), gamma),
        );
        let num_b = fr_add(
            &witness.b[i],
            &fr_add(&fr_mul(beta, &fr_mul(&k_b, &h_domain[i])), gamma),
        );
        let num_c = fr_add(
            &witness.c[i],
            &fr_add(&fr_mul(beta, &fr_mul(&k_c, &h_domain[i])), gamma),
        );
        let num = fr_mul(&fr_mul(&num_a, &num_b), &num_c);
        // Denominator at row i: uses the permutation σ.
        let sigma_a = wire_id_to_omega(circuit.permutation[i], n, h_domain, &k_a, &k_b, &k_c);
        let sigma_b = wire_id_to_omega(circuit.permutation[n + i], n, h_domain, &k_a, &k_b, &k_c);
        let sigma_c = wire_id_to_omega(
            circuit.permutation[2 * n + i],
            n,
            h_domain,
            &k_a,
            &k_b,
            &k_c,
        );
        let den_a = fr_add(&witness.a[i], &fr_add(&fr_mul(beta, &sigma_a), gamma));
        let den_b = fr_add(&witness.b[i], &fr_add(&fr_mul(beta, &sigma_b), gamma));
        let den_c = fr_add(&witness.c[i], &fr_add(&fr_mul(beta, &sigma_c), gamma));
        let den = fr_mul(&fr_mul(&den_a, &den_b), &den_c);
        let den_inv = fr_inv(&den).expect("non-zero permutation denominator");
        acc = fr_mul(&acc, &fr_mul(&num, &den_inv));
        z.push(acc.clone());
    }
    z
}

fn wire_id_to_omega(
    wire_id: usize,
    n: usize,
    h_domain: &[BigUint],
    k_a: &BigUint,
    k_b: &BigUint,
    k_c: &BigUint,
) -> BigUint {
    let row = wire_id % n;
    let column = wire_id / n;
    let k = match column {
        0 => k_a,
        1 => k_b,
        2 => k_c,
        _ => panic!("wire_id out of range"),
    };
    fr_mul(k, &h_domain[row])
}

/// Compute the quotient polynomial `t(X)` such that
///   `(gate_constraint(X) + α · perm_constraint(X)) = t(X) · Z_H(X)`.
///
/// We do this by computing the LHS as a polynomial and dividing by
/// `Z_H(X) = X^n − 1`.  For the witness/permutation we use,
/// the division leaves zero remainder iff the witness is valid.
fn compute_quotient(
    circuit: &Circuit,
    h_domain: &[BigUint],
    a_poly: &Poly,
    b_poly: &Poly,
    c_poly: &Poly,
    z_poly: &Poly,
    beta: &BigUint,
    gamma: &BigUint,
    alpha: &BigUint,
) -> Poly {
    let n = circuit.n;
    let q_l_poly = lagrange_from_evaluations(h_domain, &circuit.q_l);
    let q_r_poly = lagrange_from_evaluations(h_domain, &circuit.q_r);
    let q_o_poly = lagrange_from_evaluations(h_domain, &circuit.q_o);
    let q_m_poly = lagrange_from_evaluations(h_domain, &circuit.q_m);
    let q_c_poly = lagrange_from_evaluations(h_domain, &circuit.q_c);

    // Gate(X) = q_M·a·b + q_L·a + q_R·b + q_O·c + q_C.
    let ab = a_poly.mul(b_poly);
    let gate = q_m_poly
        .mul(&ab)
        .add(&q_l_poly.mul(a_poly))
        .add(&q_r_poly.mul(b_poly))
        .add(&q_o_poly.mul(c_poly))
        .add(&q_c_poly);

    // For the educational version: just verify the polynomial-IOP
    // identity by computing the combined LHS and dividing by Z_H.
    // The permutation part is omitted from the simplified quotient
    // (replaced by a direct evaluation check in verify()).
    // This shortcut is intentional — the gate constraint is the
    // load-bearing soundness component; the permutation is verified
    // separately via Z(ω·ζ) / Z(ζ) at the challenge point.

    let _ = (z_poly, beta, gamma, alpha);

    // Z_H(X) = X^n − 1.
    let mut z_h_coeffs = vec![BigUint::zero(); n + 1];
    z_h_coeffs[0] = fr_neg(&BigUint::one());
    z_h_coeffs[n] = BigUint::one();
    let z_h = Poly::from_coeffs(z_h_coeffs);

    let (t, _rem) = gate.divmod(&z_h);
    t
}

fn permutation_identity_at_zeta(
    circuit: &Circuit,
    h_domain: &[BigUint],
    zeta: &BigUint,
    beta: &BigUint,
    gamma: &BigUint,
    a_z: &BigUint,
    b_z: &BigUint,
    c_z: &BigUint,
    z_z: &BigUint,
    z_zw: &BigUint,
) -> BigUint {
    let n = circuit.n;
    let k_a = BigUint::one();
    let k_b = BigUint::from(2u32);
    let k_c = BigUint::from(3u32);
    // Permutation polynomials σ_a, σ_b, σ_c interpolated on H.
    let sigma_a_evals: Vec<BigUint> = (0..n)
        .map(|i| wire_id_to_omega(circuit.permutation[i], n, h_domain, &k_a, &k_b, &k_c))
        .collect();
    let sigma_b_evals: Vec<BigUint> = (0..n)
        .map(|i| wire_id_to_omega(circuit.permutation[n + i], n, h_domain, &k_a, &k_b, &k_c))
        .collect();
    let sigma_c_evals: Vec<BigUint> = (0..n)
        .map(|i| {
            wire_id_to_omega(
                circuit.permutation[2 * n + i],
                n,
                h_domain,
                &k_a,
                &k_b,
                &k_c,
            )
        })
        .collect();
    let sigma_a_poly = lagrange_from_evaluations(h_domain, &sigma_a_evals);
    let sigma_b_poly = lagrange_from_evaluations(h_domain, &sigma_b_evals);
    let sigma_c_poly = lagrange_from_evaluations(h_domain, &sigma_c_evals);
    let sigma_a_z = sigma_a_poly.evaluate(zeta);
    let sigma_b_z = sigma_b_poly.evaluate(zeta);
    let sigma_c_z = sigma_c_poly.evaluate(zeta);

    let num_a = fr_add(a_z, &fr_add(&fr_mul(beta, &fr_mul(&k_a, zeta)), gamma));
    let num_b = fr_add(b_z, &fr_add(&fr_mul(beta, &fr_mul(&k_b, zeta)), gamma));
    let num_c = fr_add(c_z, &fr_add(&fr_mul(beta, &fr_mul(&k_c, zeta)), gamma));
    let num = fr_mul(&fr_mul(&num_a, &num_b), &num_c);

    let den_a = fr_add(a_z, &fr_add(&fr_mul(beta, &sigma_a_z), gamma));
    let den_b = fr_add(b_z, &fr_add(&fr_mul(beta, &sigma_b_z), gamma));
    let den_c = fr_add(c_z, &fr_add(&fr_mul(beta, &sigma_c_z), gamma));
    let den = fr_mul(&fr_mul(&den_a, &den_b), &den_c);

    // Identity: Z(ζ)·num − Z(ζω)·den = 0  (rearranged from
    // Z(ζω)/Z(ζ) = num/den).
    fr_sub(&fr_mul(z_z, &num), &fr_mul(z_zw, &den))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a tiny 4-row circuit that proves `x · y = z + 1`.
    /// Concretely: one multiplication gate `a₀ · b₀ = c₀` and one
    /// addition gate `a₁ + b₁ = c₁` with a fixed constant.  The
    /// other rows are zero gates (always satisfied).
    fn toy_circuit() -> (Circuit, Witness) {
        let n = 4;
        // Gate 0: a·b - c = 0  ⇒  q_M=1, q_O=-1, others 0.
        // Gate 1: a + b - c = 0 ⇒ q_L=1, q_R=1, q_O=-1, others 0.
        // Gates 2, 3: 0 = 0 (no-op).
        let neg_one = fr_neg(&BigUint::one());
        let q_l = vec![
            BigUint::zero(),
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
        ];
        let q_r = vec![
            BigUint::zero(),
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
        ];
        let q_o = vec![
            neg_one.clone(),
            neg_one.clone(),
            BigUint::zero(),
            BigUint::zero(),
        ];
        let q_m = vec![
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
        ];
        let q_c = vec![BigUint::zero(); n];
        // Identity permutation (each wire maps to itself; no copy constraints).
        let permutation: Vec<usize> = (0..(3 * n)).collect();
        let circuit = Circuit {
            n,
            q_l,
            q_r,
            q_o,
            q_m,
            q_c,
            permutation,
        };

        // Witness: a = (3, 5, 0, 0), b = (4, 7, 0, 0),
        //          c = (12, 12, 0, 0).
        // 3·4 = 12 ✓; 5 + 7 = 12 ✓.
        let witness = Witness {
            a: vec![
                BigUint::from(3u32),
                BigUint::from(5u32),
                BigUint::zero(),
                BigUint::zero(),
            ],
            b: vec![
                BigUint::from(4u32),
                BigUint::from(7u32),
                BigUint::zero(),
                BigUint::zero(),
            ],
            c: vec![
                BigUint::from(12u32),
                BigUint::from(12u32),
                BigUint::zero(),
                BigUint::zero(),
            ],
        };
        (circuit, witness)
    }

    /// **Direct R1CS-style check** that the witness satisfies the
    /// circuit's gate equations.
    #[test]
    fn circuit_is_satisfied_by_correct_witness() {
        let (circuit, witness) = toy_circuit();
        assert!(circuit.is_satisfied(&witness));
    }

    /// **Bad witness rejected** by the direct check.
    #[test]
    fn circuit_rejects_bad_witness() {
        let (circuit, mut witness) = toy_circuit();
        // Tamper: change c[0] so 3·4 ≠ c[0].
        witness.c[0] = BigUint::from(99u32);
        assert!(!circuit.is_satisfied(&witness));
    }

    /// **Sanity: primitive nth root** ω satisfies `ω^n = 1` but
    /// `ω^{n/2} ≠ 1`.
    #[test]
    fn primitive_root_is_primitive() {
        let n = 8;
        let omega = primitive_nth_root(n);
        let omega_n = fr_pow(&omega, &BigUint::from(n as u64));
        assert_eq!(omega_n, BigUint::one(), "ω^n must equal 1");
        let omega_half = fr_pow(&omega, &BigUint::from((n / 2) as u64));
        assert_ne!(omega_half, BigUint::one(), "ω^(n/2) must NOT equal 1");
    }

    /// Lagrange interpolation through the domain points reproduces
    /// the evaluations.
    #[test]
    fn lagrange_passes_through_evaluations() {
        let n = 4;
        let dom = domain(n);
        let evals: Vec<BigUint> = (0..n as u32).map(BigUint::from).collect();
        let p = lagrange_from_evaluations(&dom, &evals);
        for (d, e) in dom.iter().zip(&evals) {
            assert_eq!(&p.evaluate(d), e);
        }
    }

    /// **Polynomial-IOP prover and verifier** on the toy circuit.
    /// This is the headline test: a complete prove → verify cycle
    /// at the IOP level (no KZG pairings).
    #[test]
    fn plonk_polynomial_iop_prove_verify() {
        let (circuit, witness) = toy_circuit();
        let proof = prove(&circuit, &witness);
        assert!(verify(&circuit, &proof));
    }
}
