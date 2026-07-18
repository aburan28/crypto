//! **HFEv-** — Hidden Field Equations with the vinegar and minus
//! modifiers (Patarin 1996; the basis of GeMSS and Gui, NIST round-2/3
//! multivariate signatures).
//!
//! # A different multivariate trapdoor
//! UOV/Rainbow build their trapdoor from an oil/vinegar *structure*.
//! HFE instead hides a **single univariate polynomial over an extension
//! field**.  Work in `F_{qⁿ}` (here `F₂₈ = GF(256)`).  The secret is a
//! polynomial of special "HFE shape"
//!
//! ```text
//! F(X) = Σ a_{ij} X^{qⁱ+qʲ}  +  Σ b_i X^{qⁱ}  +  c,
//! ```
//!
//! whose exponents are all sums of two powers of `q`.  Two facts make
//! it a trapdoor:
//! - Written in an `F_q`-basis, `X ↦ F(X)` is a system of `n`
//!   **quadratic** polynomials in `n` base-field variables (because
//!   `X^{qⁱ}` is `F_q`-linear, so `X^{qⁱ+qʲ}` is quadratic).  That
//!   system, composed with secret affine maps `S, T`, is the public key.
//! - Knowing `F` as a univariate polynomial, the signer inverts it by
//!   **root-finding** over `F_{qⁿ}`.
//!
//! # The ...v- modifiers (defeating the early HFE attacks)
//! Plain HFE was weakened by Kipnis–Shamir's MinRank attack.  Two cheap
//! modifiers restore security:
//! - **vinegar (v)**: add `v` extra variables that enter `F`'s linear
//!   and constant coefficients; the signer fixes them at random and
//!   inverts the resulting univariate polynomial.
//! - **minus (−)**: publish only `n − minus` of the `n` equations; the
//!   signer guesses the removed coordinates.
//!
//! # This implementation
//! `F₂₈`, `n = 8`, HFE exponents from `i, j ∈ {0,1,2}`, vinegar `v = 2`,
//! minus `= 1` (so the public key is 7 quadratics in 10 `F₂`
//! variables).  Univariate inversion is exhaustive root search over the
//! 256 field elements — trivial here, the analog of the Berlekamp
//! root-finding a real HFE signer runs.  The public quadratic system is
//! extracted by evaluation, exactly as a deployed HFE key would be.
//! Toy scale, not constant-time; see SECURITY.md.

use crate::hash::sha3::shake256;
use crate::utils::random::random_bytes;

/// Extension degree over F₂ (F₂₈ = GF(256)).
pub const N: usize = 8;
/// Vinegar variables.
pub const V: usize = 2;
/// Minus: equations removed from the public key.
pub const MINUS: usize = 1;
/// Public input variables.
pub const NV: usize = N + V; // 10
/// Public equations.
pub const PUB_EQ: usize = N - MINUS; // 7

// ── GF(256) ───────────────────────────────────────────────────────────────────

fn gf_mul(a: u8, b: u8) -> u8 {
    let (mut a, mut b, mut r) = (a as u16, b as u16, 0u16);
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        a <<= 1;
        if a & 0x100 != 0 {
            a ^= 0x11b;
        }
        b >>= 1;
    }
    r as u8
}

/// Frobenius `X^{2^i}` = square `i` times.
fn frob(mut x: u8, i: usize) -> u8 {
    for _ in 0..i {
        x = gf_mul(x, x);
    }
    x
}

// ── The HFEv core polynomial ──────────────────────────────────────────────────

/// HFE exponent set: pairs (i, j) with i ≤ j drawn from {0,1,2}.
const EXPS: [(usize, usize); 6] = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];

/// Secret HFEv core: quadratic coefficients `a_{ij}`, and vinegar-affine
/// linear/constant coefficients.
#[derive(Clone)]
struct HfeCore {
    a: [u8; 6],          // coefficient of X^{2^i+2^j} for each EXPS entry
    b: [[u8; V + 1]; 3], // β_i(vin) = b[i][0] ⊕ Σ_l b[i][l+1]·vin_l, for i∈{0,1,2}
    g: [u8; V + 1],      // γ(vin) = g[0] ⊕ Σ_l g[l+1]·vin_l
}

impl HfeCore {
    fn random() -> Self {
        let mut a = [0u8; 6];
        random_bytes(&mut a);
        let mut b = [[0u8; V + 1]; 3];
        for row in b.iter_mut() {
            random_bytes(row);
        }
        let mut g = [0u8; V + 1];
        random_bytes(&mut g);
        HfeCore { a, b, g }
    }

    /// Evaluate `F(X, vin) ∈ F₂₈`.
    fn eval(&self, x: u8, vin: &[u8]) -> u8 {
        let mut acc = 0u8;
        for (k, &(i, j)) in EXPS.iter().enumerate() {
            if self.a[k] != 0 {
                let term = gf_mul(frob(x, i), frob(x, j));
                acc ^= gf_mul(self.a[k], term);
            }
        }
        for i in 0..3 {
            // β_i(vin) as an F₂-affine combination (vinegar bits scale
            // fixed field coefficients).
            let mut beta = self.b[i][0];
            for (l, &vl) in vin.iter().enumerate() {
                if vl & 1 == 1 {
                    beta ^= self.b[i][l + 1];
                }
            }
            if beta != 0 {
                acc ^= gf_mul(beta, frob(x, i));
            }
        }
        let mut gamma = self.g[0];
        for (l, &vl) in vin.iter().enumerate() {
            if vl & 1 == 1 {
                gamma ^= self.g[l + 1];
            }
        }
        acc ^ gamma
    }

    /// Invert: all `X` with `F(X, vin) = y` (exhaustive over the field).
    fn preimages(&self, y: u8, vin: &[u8]) -> Vec<u8> {
        (0u16..256).filter(|&x| self.eval(x as u8, vin) == y).map(|x| x as u8).collect()
    }
}

// ── F₂ affine maps ────────────────────────────────────────────────────────────

/// Invertible F₂ matrix + inverse (size × size).
fn random_gf2_invertible(size: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    loop {
        let mut m = vec![vec![0u8; size]; size];
        for row in m.iter_mut() {
            let mut b = vec![0u8; size];
            random_bytes(&mut b);
            for (c, x) in row.iter_mut().zip(&b) {
                *c = x & 1;
            }
        }
        if let Some(inv) = gf2_invert(&m) {
            return (m, inv);
        }
    }
}

fn gf2_invert(m: &[Vec<u8>]) -> Option<Vec<Vec<u8>>> {
    let n = m.len();
    let mut a: Vec<Vec<u8>> = m.to_vec();
    let mut inv = vec![vec![0u8; n]; n];
    for i in 0..n {
        inv[i][i] = 1;
    }
    for col in 0..n {
        let piv = (col..n).find(|&r| a[r][col] == 1)?;
        a.swap(col, piv);
        inv.swap(col, piv);
        for r in 0..n {
            if r != col && a[r][col] == 1 {
                for j in 0..n {
                    a[r][j] ^= a[col][j];
                    inv[r][j] ^= inv[col][j];
                }
            }
        }
    }
    Some(inv)
}

fn gf2_apply(m: &[Vec<u8>], v: &[u8], c: &[u8]) -> Vec<u8> {
    (0..m.len())
        .map(|i| {
            let mut acc = c[i] & 1;
            for j in 0..v.len() {
                acc ^= m[i][j] & v[j] & 1;
            }
            acc
        })
        .collect()
}

fn byte_to_bits(x: u8) -> Vec<u8> {
    (0..N).map(|i| (x >> i) & 1).collect()
}

fn bits_to_byte(b: &[u8]) -> u8 {
    (0..N).fold(0u8, |acc, i| acc | ((b[i] & 1) << i))
}

// ── Public key: extracted F₂ quadratic system ────────────────────────────────

/// One quadratic form over `F₂` in `NV` variables:
/// `Σ_{i<j} A_ij w_i w_j ⊕ Σ_i L_i w_i ⊕ c`.
#[derive(Clone)]
pub struct QuadForm {
    a: Vec<Vec<u8>>, // upper-triangular NV×NV
    l: Vec<u8>,      // NV
    c: u8,
}

impl QuadForm {
    fn eval(&self, w: &[u8]) -> u8 {
        let mut acc = self.c;
        for i in 0..NV {
            if w[i] & 1 == 1 {
                acc ^= self.l[i];
                for j in (i + 1)..NV {
                    if self.a[i][j] & 1 == 1 && w[j] & 1 == 1 {
                        acc ^= 1;
                    }
                }
            }
        }
        acc & 1
    }
}

#[derive(Clone)]
pub struct HfePublicKey {
    pub forms: Vec<QuadForm>, // PUB_EQ forms
}

#[derive(Clone)]
pub struct HfeSecretKey {
    core: HfeCore,
    s_inv: Vec<Vec<u8>>,
    s_c: Vec<u8>,
    t_mat: Vec<Vec<u8>>,
    t_inv: Vec<Vec<u8>>,
    t_c: Vec<u8>,
}

/// The full composed public map `P: F₂^{NV} → F₂^N` (before minus).
fn public_map_full(sk: &HfeSecretKey, s_mat: &[Vec<u8>], w: &[u8]) -> Vec<u8> {
    // z = S(w).
    let z = gf2_apply(s_mat, w, &sk.s_c);
    let x = bits_to_byte(&z[..N]);
    let vin = &z[N..];
    let y = sk.core.eval(x, vin);
    let y_bits = byte_to_bits(y);
    // P = T(y_bits).
    gf2_apply(&sk.t_mat, &y_bits, &sk.t_c)
}

/// Extract the `N` quadratic forms of `public_map_full` by evaluation.
fn extract_forms(sk: &HfeSecretKey, s_mat: &[Vec<u8>]) -> Vec<QuadForm> {
    let zero = vec![0u8; NV];
    let base = public_map_full(sk, s_mat, &zero); // constants
    let unit = |i: usize| {
        let mut e = vec![0u8; NV];
        e[i] = 1;
        public_map_full(sk, s_mat, &e)
    };
    let units: Vec<Vec<u8>> = (0..NV).map(unit).collect();

    (0..N)
        .map(|k| {
            let c = base[k];
            let l: Vec<u8> = (0..NV).map(|i| units[i][k] ^ c).collect();
            let mut a = vec![vec![0u8; NV]; NV];
            for i in 0..NV {
                for j in (i + 1)..NV {
                    let mut e = vec![0u8; NV];
                    e[i] = 1;
                    e[j] = 1;
                    let pij = public_map_full(sk, s_mat, &e)[k];
                    // A_ij = P(e_i+e_j) ⊕ P(e_i) ⊕ P(e_j) ⊕ c.
                    a[i][j] = pij ^ units[i][k] ^ units[j][k] ^ c;
                }
            }
            QuadForm { a, l, c }
        })
        .collect()
}

pub fn hfe_keygen() -> (HfePublicKey, HfeSecretKey) {
    let core = HfeCore::random();
    let (s_mat, s_inv) = random_gf2_invertible(NV);
    let (t_mat, t_inv) = random_gf2_invertible(N);
    let mut s_c = vec![0u8; NV];
    let mut t_c = vec![0u8; N];
    random_bytes(&mut s_c);
    random_bytes(&mut t_c);
    for x in s_c.iter_mut() {
        *x &= 1;
    }
    for x in t_c.iter_mut() {
        *x &= 1;
    }
    let sk = HfeSecretKey { core, s_inv, s_c, t_mat, t_inv, t_c };
    let all_forms = extract_forms(&sk, &s_mat);
    // Minus: publish only the first PUB_EQ forms.
    let forms = all_forms.into_iter().take(PUB_EQ).collect();
    (HfePublicKey { forms }, sk)
}

fn hash_to_target(msg: &[u8]) -> Vec<u8> {
    let mut input = b"HFEv-toy".to_vec();
    input.extend_from_slice(msg);
    shake256(&input, PUB_EQ).iter().map(|b| b & 1).collect()
}

pub fn hfe_sign(sk: &HfeSecretKey, msg: &[u8]) -> Vec<u8> {
    let target = hash_to_target(msg);
    loop {
        // Guess the MINUS removed output coordinates.
        let mut removed = vec![0u8; MINUS];
        random_bytes(&mut removed);
        for x in removed.iter_mut() {
            *x &= 1;
        }
        let mut full_target = target.clone();
        full_target.extend_from_slice(&removed);

        // y_bits = T⁻¹(full_target ⊕ t_c); Y = field element.
        let shifted: Vec<u8> = (0..N).map(|i| full_target[i] ^ sk.t_c[i]).collect();
        let y_bits = gf2_apply(&sk.t_inv, &shifted, &vec![0u8; N]);
        let y = bits_to_byte(&y_bits);

        // Random vinegar, then invert the univariate HFE polynomial.
        let mut vin = vec![0u8; V];
        random_bytes(&mut vin);
        for x in vin.iter_mut() {
            *x &= 1;
        }
        let roots = sk.core.preimages(y, &vin);
        if roots.is_empty() {
            continue; // no root: fresh minus/vinegar guess
        }
        let x = roots[0];

        // z = (φ(X) ‖ vin); w = S⁻¹(z ⊕ s_c).
        let mut z = byte_to_bits(x);
        z.extend_from_slice(&vin);
        let shifted_z: Vec<u8> = (0..NV).map(|i| z[i] ^ sk.s_c[i]).collect();
        let w = gf2_apply(&sk.s_inv, &shifted_z, &vec![0u8; NV]);
        return w;
    }
}

pub fn hfe_verify(pk: &HfePublicKey, msg: &[u8], sig: &[u8]) -> bool {
    if sig.len() != NV {
        return false;
    }
    let target = hash_to_target(msg);
    pk.forms.iter().zip(&target).all(|(form, &t)| form.eval(sig) == t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_is_invertible_by_search() {
        let core = HfeCore::random();
        let vin = vec![1u8, 0];
        // Every value in the image has a preimage found by the search.
        let x = 0x9c;
        let y = core.eval(x, &vin);
        assert!(core.preimages(y, &vin).contains(&x));
    }

    #[test]
    fn bit_field_roundtrip() {
        for x in 0u16..256 {
            assert_eq!(bits_to_byte(&byte_to_bits(x as u8)), x as u8);
        }
    }

    #[test]
    fn extracted_forms_match_composed_map() {
        // The extracted quadratic system must equal the true composed
        // public map on random inputs.
        let (s_mat, s_inv) = random_gf2_invertible(NV);
        let (t_mat, t_inv) = random_gf2_invertible(N);
        let sk = HfeSecretKey {
            core: HfeCore::random(),
            s_inv,
            s_c: vec![0u8; NV],
            t_mat,
            t_inv,
            t_c: vec![0u8; N],
        };
        let forms = extract_forms(&sk, &s_mat);
        for _ in 0..20 {
            let mut w = vec![0u8; NV];
            random_bytes(&mut w);
            for x in w.iter_mut() {
                *x &= 1;
            }
            let truth = public_map_full(&sk, &s_mat, &w);
            for k in 0..N {
                assert_eq!(forms[k].eval(&w), truth[k], "form {k}");
            }
        }
    }

    #[test]
    fn sign_verify_roundtrip() {
        let (pk, sk) = hfe_keygen();
        let msg = b"hidden field equations";
        let sig = hfe_sign(&sk, msg);
        assert_eq!(sig.len(), NV);
        assert!(hfe_verify(&pk, msg, &sig));
    }

    #[test]
    fn several_messages() {
        let (pk, sk) = hfe_keygen();
        for i in 0u8..8 {
            let m = [i; 5];
            let sig = hfe_sign(&sk, &m);
            assert!(hfe_verify(&pk, &m, &sig), "message {i}");
        }
    }

    #[test]
    fn wrong_message_rejected() {
        let (pk, sk) = hfe_keygen();
        let sig = hfe_sign(&sk, b"one");
        assert!(!hfe_verify(&pk, b"two", &sig));
    }

    #[test]
    fn tampered_signature_rejected() {
        let (pk, sk) = hfe_keygen();
        let msg = b"tamper";
        let mut sig = hfe_sign(&sk, msg);
        sig[0] ^= 1;
        assert!(!hfe_verify(&pk, msg, &sig));
    }
}
