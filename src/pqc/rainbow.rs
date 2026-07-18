//! **Rainbow** — layered Unbalanced Oil-and-Vinegar signatures
//! (Ding–Schmidt 2005; NIST round-3 *finalist*, broken by Beullens
//! 2022).
//!
//! # Rainbow = UOV in layers
//! Plain UOV (`pqc::uov`) has one oil/vinegar split.  Rainbow stacks
//! **two** layers to shrink keys.  Partition the `n` variables into
//! nested sets `V₁ ⊂ V₂ ⊂ V₃ = {1..n}`:
//!
//! - **Layer 1** treats `V₁` as vinegar and `O₁ = V₂ ∖ V₁` as oil:
//!   its `|O₁|` central equations contain no oil×oil term among `O₁`.
//! - **Layer 2** treats *all of* `V₂` as vinegar and `O₂ = V₃ ∖ V₂` as
//!   oil: its `|O₂|` equations contain no oil×oil term among `O₂`.
//!
//! Signing inverts the layers in turn: fix the `V₁` vinegars at random,
//! solve layer 1's linear system for `O₁`; now all of `V₂` is known, so
//! layer 2 is linear in `O₂` — solve it.  As in UOV, the public map is
//! `P = S ∘ F ∘ T` for secret affine `S, T`, hiding the layered
//! structure so `P` looks like a random quadratic system.
//!
//! # The 2022 break
//! Beullens' "Breaking Rainbow takes a weekend on a laptop" recovered
//! the secret layer structure via a **MinRank / differential** attack
//! exploiting exactly the band structure that makes signing efficient —
//! the layers leave a low-rank signature in the differential of the
//! public map.  This knocked Rainbow out of NIST standardisation and is
//! why the surviving multivariate candidates (UOV, MAYO, …) avoid the
//! extra layer.  A cautionary tale about structure enabling attacks.
//!
//! # This implementation
//! Two layers over GF(256), toy params `v₁ = 6, o₁ = 4, o₂ = 4`
//! (`n = 14`, `m = 8`).  Faithful to the layered central map and the
//! layer-by-layer inversion; toy sizes offer no security (and the
//! scheme is broken even at real sizes).  Not constant-time; see
//! SECURITY.md.

use crate::hash::sha3::shake256;
use crate::utils::random::random_bytes;

/// Layer-1 vinegar variables.
pub const V1: usize = 6;
/// Layer-1 oil variables.
pub const O1: usize = 4;
/// Layer-2 oil variables.
pub const O2: usize = 4;
/// Total variables.
pub const N: usize = V1 + O1 + O2; // 14
/// Total equations.
pub const M: usize = O1 + O2; // 8
/// End of layer-1 vinegar block (= start of O1) — also |V₂| once O1 added.
const V2: usize = V1 + O1; // 10

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

fn gf_inv(a: u8) -> u8 {
    let mut r = 1u8;
    let mut base = a;
    let mut e = 254u32;
    while e > 0 {
        if e & 1 == 1 {
            r = gf_mul(r, base);
        }
        base = gf_mul(base, base);
        e >>= 1;
    }
    r
}

type Matrix = Vec<Vec<u8>>;

fn zero(r: usize, c: usize) -> Matrix {
    vec![vec![0u8; c]; r]
}

fn solve(a: &Matrix, b: &[u8]) -> Option<Vec<u8>> {
    let n = a.len();
    if n == 0 {
        return Some(vec![]);
    }
    let mut m: Matrix = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let piv = (col..n).find(|&r| m[r][col] != 0)?;
        m.swap(col, piv);
        rhs.swap(col, piv);
        let inv = gf_inv(m[col][col]);
        for x in m[col].iter_mut() {
            *x = gf_mul(*x, inv);
        }
        rhs[col] = gf_mul(rhs[col], inv);
        for r in 0..n {
            if r != col && m[r][col] != 0 {
                let f = m[r][col];
                for j in 0..n {
                    m[r][j] ^= gf_mul(f, m[col][j]);
                }
                rhs[r] ^= gf_mul(f, rhs[col]);
            }
        }
    }
    Some(rhs)
}

fn invert(a: &Matrix) -> Option<Matrix> {
    let n = a.len();
    let mut m = a.to_vec();
    let mut inv = zero(n, n);
    for i in 0..n {
        inv[i][i] = 1;
    }
    for col in 0..n {
        let piv = (col..n).find(|&r| m[r][col] != 0)?;
        m.swap(col, piv);
        inv.swap(col, piv);
        let iv = gf_inv(m[col][col]);
        for j in 0..n {
            m[col][j] = gf_mul(m[col][j], iv);
            inv[col][j] = gf_mul(inv[col][j], iv);
        }
        for r in 0..n {
            if r != col && m[r][col] != 0 {
                let f = m[r][col];
                for j in 0..n {
                    m[r][j] ^= gf_mul(f, m[col][j]);
                    inv[r][j] ^= gf_mul(f, inv[col][j]);
                }
            }
        }
    }
    Some(inv)
}

fn random_invertible(n: usize) -> (Matrix, Matrix) {
    loop {
        let mut m = zero(n, n);
        for row in m.iter_mut() {
            random_bytes(row);
        }
        if let Some(inv) = invert(&m) {
            return (m, inv);
        }
    }
}

fn apply_matrix(m: &Matrix, v: &[u8]) -> Vec<u8> {
    m.iter().map(|row| row.iter().zip(v).fold(0u8, |a, (&x, &y)| a ^ gf_mul(x, y))).collect()
}

fn eval_form(q: &Matrix, x: &[u8]) -> u8 {
    let mut acc = 0u8;
    for (i, row) in q.iter().enumerate() {
        if x[i] == 0 {
            continue;
        }
        let mut inner = 0u8;
        for (j, &qij) in row.iter().enumerate() {
            if qij != 0 {
                inner ^= gf_mul(qij, x[j]);
            }
        }
        acc ^= gf_mul(x[i], inner);
    }
    acc
}

// ── Keys ──────────────────────────────────────────────────────────────────────

#[derive(Clone)]
pub struct RainbowPublicKey {
    pub forms: Vec<Matrix>, // m public quadratic forms
}

#[derive(Clone)]
pub struct RainbowSecretKey {
    central: Vec<Matrix>,
    s_inv: Matrix,
    t_inv: Matrix,
}

/// Build a central-map form for equation `k` respecting the layer rule:
/// upper-triangular entry `(i, j)` (i ≤ j) is allowed unless it is an
/// oil×oil term of that equation's layer.
fn random_central() -> Vec<Matrix> {
    let mut central = Vec::with_capacity(M);
    for k in 0..M {
        // Layer 1: equations 0..O1, oil block is [V1, V2); vinegar [0, V1).
        // Layer 2: equations O1..M, oil block is [V2, N); vinegar [0, V2).
        let (vinegar_end, oil_start, oil_end) =
            if k < O1 { (V1, V1, V2) } else { (V2, V2, N) };
        let mut q = zero(N, N);
        for i in 0..N {
            for j in i..N {
                // Skip variables the equation must not see: layer-2 oil
                // is invisible to layer-1 equations.
                let max_var = oil_end;
                if i >= max_var || j >= max_var {
                    continue;
                }
                // Forbid oil×oil in this layer.
                let i_oil = i >= oil_start;
                let j_oil = j >= oil_start;
                if i_oil && j_oil {
                    continue;
                }
                // Also forbid vinegar-only quadratics from exceeding the
                // vinegar/oil blocks incoherently — allowed set is
                // vinegar×vinegar and vinegar×oil.
                let _ = vinegar_end;
                let mut b = [0u8; 1];
                random_bytes(&mut b);
                q[i][j] = b[0];
            }
        }
        central.push(q);
    }
    central
}

fn hash_to_target(msg: &[u8]) -> Vec<u8> {
    let mut input = b"Rainbow-toy".to_vec();
    input.extend_from_slice(msg);
    shake256(&input, M)
}

pub fn rainbow_keygen() -> (RainbowPublicKey, RainbowSecretKey) {
    let central = random_central();
    let (s, s_inv) = random_invertible(M);
    let (t, t_inv) = random_invertible(N);
    // Public forms P_k = Σ_l S[k][l] · (Tᵀ · F_l · T).
    let t_t: Matrix = (0..N).map(|i| (0..N).map(|j| t[j][i]).collect()).collect();
    let composed: Vec<Matrix> = central
        .iter()
        .map(|f| {
            // Tᵀ F T.
            let ft: Matrix = (0..N)
                .map(|i| {
                    (0..N)
                        .map(|j| (0..N).fold(0u8, |a, l| a ^ gf_mul(f[i][l], t[l][j])))
                        .collect()
                })
                .collect();
            (0..N)
                .map(|i| {
                    (0..N)
                        .map(|j| (0..N).fold(0u8, |a, l| a ^ gf_mul(t_t[i][l], ft[l][j])))
                        .collect()
                })
                .collect()
        })
        .collect();
    // Mix equations by S.
    let forms: Vec<Matrix> = (0..M)
        .map(|k| {
            let mut acc = zero(N, N);
            for l in 0..M {
                if s[k][l] != 0 {
                    for i in 0..N {
                        for j in 0..N {
                            acc[i][j] ^= gf_mul(s[k][l], composed[l][i][j]);
                        }
                    }
                }
            }
            acc
        })
        .collect();
    (RainbowPublicKey { forms }, RainbowSecretKey { central, s_inv, t_inv })
}

pub fn rainbow_sign(sk: &RainbowSecretKey, msg: &[u8]) -> Vec<u8> {
    let target = hash_to_target(msg);
    loop {
        // u = S⁻¹ · target: invert the output mixing.
        let u = apply_matrix(&sk.s_inv, &target);

        // Fix layer-1 vinegar at random.
        let mut x = vec![0u8; N];
        random_bytes(&mut x[..V1]);

        // Layer 1: equations 0..O1 are linear in O1 = [V1, V2).
        if let Some(o1) = solve_layer(&sk.central, &u, &x, 0, O1, V1, V2) {
            for (i, val) in o1.iter().enumerate() {
                x[V1 + i] = *val;
            }
        } else {
            continue; // singular: fresh vinegar
        }

        // Layer 2: equations O1..M are linear in O2 = [V2, N).
        if let Some(o2) = solve_layer(&sk.central, &u, &x, O1, M, V2, N) {
            for (i, val) in o2.iter().enumerate() {
                x[V2 + i] = *val;
            }
        } else {
            continue;
        }

        // s = T⁻¹ · x.
        return apply_matrix(&sk.t_inv, &x);
    }
}

/// Solve one layer: equations `k ∈ [eq_lo, eq_hi)` are affine in the oil
/// block `[oil_lo, oil_hi)` once earlier variables are fixed in `x`.
fn solve_layer(
    central: &[Matrix],
    u: &[u8],
    x: &[u8],
    eq_lo: usize,
    eq_hi: usize,
    oil_lo: usize,
    oil_hi: usize,
) -> Option<Vec<u8>> {
    let rows = eq_hi - eq_lo;
    let cols = oil_hi - oil_lo;
    debug_assert_eq!(rows, cols);
    let mut a = zero(rows, cols);
    let mut b = vec![0u8; rows];
    for (r, k) in (eq_lo..eq_hi).enumerate() {
        let q = &central[k];
        // Constant part: evaluate with oil = 0 (they are 0 in x already).
        let mut xt = x.to_vec();
        for i in oil_lo..oil_hi {
            xt[i] = 0;
        }
        let c = eval_form(q, &xt);
        // Linear coefficient of each oil var: difference when that oil is 1.
        for (col, ovar) in (oil_lo..oil_hi).enumerate() {
            let mut xu = xt.clone();
            xu[ovar] = 1;
            a[r][col] = eval_form(q, &xu) ^ c;
        }
        b[r] = u[k] ^ c;
    }
    solve(&a, &b)
}

pub fn rainbow_verify(pk: &RainbowPublicKey, msg: &[u8], sig: &[u8]) -> bool {
    if sig.len() != N {
        return false;
    }
    let target = hash_to_target(msg);
    pk.forms.iter().zip(&target).all(|(q, &t)| eval_form(q, sig) == t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layer_structure_forbids_oil_squares() {
        // Layer-1 equations must have zero O1×O1 block and never touch O2;
        // layer-2 equations must have zero O2×O2 block.
        let central = random_central();
        for (k, q) in central.iter().enumerate() {
            if k < O1 {
                for i in V1..V2 {
                    for j in V1..V2 {
                        assert_eq!(q[i][j], 0, "layer-1 O1×O1 at ({i},{j})");
                    }
                }
                for i in V2..N {
                    for j in 0..N {
                        assert_eq!(q[i][j], 0, "layer-1 must not see O2");
                        assert_eq!(q[j][i], 0);
                    }
                }
            } else {
                for i in V2..N {
                    for j in V2..N {
                        assert_eq!(q[i][j], 0, "layer-2 O2×O2 at ({i},{j})");
                    }
                }
            }
        }
    }

    #[test]
    fn sign_verify_roundtrip() {
        let (pk, sk) = rainbow_keygen();
        let msg = b"layered oil and vinegar";
        let sig = rainbow_sign(&sk, msg);
        assert_eq!(sig.len(), N);
        assert!(rainbow_verify(&pk, msg, &sig));
    }

    #[test]
    fn wrong_message_rejected() {
        let (pk, sk) = rainbow_keygen();
        let sig = rainbow_sign(&sk, b"one");
        assert!(!rainbow_verify(&pk, b"two", &sig));
    }

    #[test]
    fn tampered_signature_rejected() {
        let (pk, sk) = rainbow_keygen();
        let msg = b"tamper";
        let mut sig = rainbow_sign(&sk, msg);
        sig[0] ^= 1;
        assert!(!rainbow_verify(&pk, msg, &sig));
    }

    #[test]
    fn signatures_randomized() {
        let (pk, sk) = rainbow_keygen();
        let msg = b"same";
        let s1 = rainbow_sign(&sk, msg);
        let s2 = rainbow_sign(&sk, msg);
        assert!(rainbow_verify(&pk, msg, &s1));
        assert!(rainbow_verify(&pk, msg, &s2));
        assert_ne!(s1, s2);
    }
}
