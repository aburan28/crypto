//! **GGH** — the Goldreich–Goldwasser–Halevi lattice cryptosystem
//! (1997), and a demonstration of *why it is broken*.
//!
//! # The idea
//! GGH was one of the first lattice encryption schemes and the direct
//! ancestor of NTRUSign/Falcon-style "hash and sign".  A lattice `Λ`
//! has many bases; some are **good** (short, near-orthogonal vectors)
//! and some are **bad** (long, skewed).  With a *good* basis you can
//! solve the closest-vector problem (CVP) approximately via Babai's
//! rounding; with a *bad* basis you cannot.
//!
//! - **Private key**: a good basis `R` (here `R = large·I + small`).
//! - **Public key**: a bad basis `B = U·R` for a random unimodular `U`.
//! - **Encrypt** a message vector `m` (small integers): the lattice
//!   point `m·B` plus a small error `e`: `c = m·B + e`.
//! - **Decrypt**: Babai-round `c` against the good basis `R` to recover
//!   the lattice point, then multiply by `B⁻¹` to read off `m`.
//!
//! # Why it's broken (Nguyen 1999)
//! The error `e` is chosen from a *fixed small box* (e.g. each
//! coordinate `±σ`).  Nguyen observed that `c ≡ m·B + e (mod 2σ)`
//! leaks `e mod 2σ`, and that subtracting the known error-centre
//! massively shrinks the effective error, reducing the CVP instance to
//! one an *attacker's* lattice reduction can solve.  The scheme's error
//! distribution simply isn't wide enough to hide the message.  This is
//! the cautionary tale that motivates Regev's Gaussian errors (LWE) and
//! Falcon's Gaussian sampler.
//!
//! # This implementation
//! Toy dimension `n = 8` with integer arithmetic.  Provides keygen,
//! encrypt, decrypt, **and** `nguyen_reduce_error`, which demonstrates
//! the `mod 2σ` leak that guts the scheme's security.  Not
//! constant-time, and — as the module name warns — **not secure by
//! design**.  Educational only; see SECURITY.md.

use crate::utils::random::random_bytes;

/// Lattice dimension.
pub const N: usize = 8;
/// Diagonal magnitude of the good (private) basis.
const DIAG: i64 = 100;
/// Off-diagonal perturbation bound of the good basis.
const PERT: i64 = 4;
/// Error box half-width σ (the fatal fixed error distribution).
pub const SIGMA: i64 = 3;
/// Number of unimodular shear rounds building the public basis.
const SHEAR_ROUNDS: usize = 3;

type Vector = Vec<i64>;
type Matrix = Vec<Vector>;

fn rand_range(lo: i64, hi: i64) -> i64 {
    let mut b = [0u8; 8];
    random_bytes(&mut b);
    let span = (hi - lo + 1) as u64;
    lo + (u64::from_le_bytes(b) % span) as i64
}

fn identity(n: usize) -> Matrix {
    let mut m = vec![vec![0i64; n]; n];
    for i in 0..n {
        m[i][i] = 1;
    }
    m
}

fn mat_vec(m: &Matrix, v: &[i64]) -> Vector {
    (0..m.len()).map(|i| (0..v.len()).map(|j| m[i][j] * v[j]).sum()).collect()
}

/// Row-vector times matrix: `v·M` (used for `m·B`).
fn vec_mat(v: &[i64], m: &Matrix) -> Vector {
    let cols = m[0].len();
    (0..cols).map(|j| (0..v.len()).map(|i| v[i] * m[i][j]).sum()).collect()
}

fn mat_mul(a: &Matrix, b: &Matrix) -> Matrix {
    let n = a.len();
    let mut out = vec![vec![0i64; n]; n];
    for i in 0..n {
        for k in 0..n {
            if a[i][k] == 0 {
                continue;
            }
            for j in 0..n {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

/// Exact integer matrix inverse via Bareiss/adjugate — feasible at the
/// toy dimension.  Returns `(inverse_scaled, det)` with
/// `inverse_scaled = det · M⁻¹` integral.
fn inverse_scaled(m: &Matrix) -> (Matrix, i64) {
    // Gaussian elimination over rationals tracked as (num, den) is
    // overkill here; use cofactor expansion for n = 8.
    let n = m.len();
    let det = determinant(m);
    let mut adj = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let minor = minor_matrix(m, j, i); // transpose of cofactor
            let sign = if (i + j) % 2 == 0 { 1 } else { -1 };
            adj[i][j] = sign * determinant(&minor);
        }
    }
    (adj, det)
}

fn minor_matrix(m: &Matrix, skip_r: usize, skip_c: usize) -> Matrix {
    let n = m.len();
    let mut out = Vec::with_capacity(n - 1);
    for i in 0..n {
        if i == skip_r {
            continue;
        }
        let mut row = Vec::with_capacity(n - 1);
        for j in 0..n {
            if j == skip_c {
                continue;
            }
            row.push(m[i][j]);
        }
        out.push(row);
    }
    out
}

fn determinant(m: &Matrix) -> i64 {
    let n = m.len();
    if n == 1 {
        return m[0][0];
    }
    if n == 2 {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }
    let mut det = 0i64;
    for j in 0..n {
        if m[0][j] == 0 {
            continue;
        }
        let sign = if j % 2 == 0 { 1 } else { -1 };
        det += sign * m[0][j] * determinant(&minor_matrix(m, 0, j));
    }
    det
}

/// Public key: the bad basis `B`.
#[derive(Clone)]
pub struct GghPublicKey {
    pub b: Matrix,
}

/// Secret key: the good basis `R`, plus `det·B⁻¹` and `det` for reading
/// off the message coordinates after decryption.
#[derive(Clone)]
pub struct GghSecretKey {
    pub r: Matrix,
    b_inv_scaled: Matrix,
    b_det: i64,
}

pub fn ggh_keygen() -> (GghPublicKey, GghSecretKey) {
    // Good basis: dominant diagonal + small perturbation ⇒ short,
    // near-orthogonal rows.
    let mut r = vec![vec![0i64; N]; N];
    for i in 0..N {
        for j in 0..N {
            r[i][j] = if i == j { DIAG } else { rand_range(-PERT, PERT) };
        }
    }
    // Unimodular U (product of ±1 shears) ⇒ same lattice, bad basis.
    let mut u = identity(N);
    for _ in 0..SHEAR_ROUNDS {
        for i in 0..N {
            let mut j = rand_range(0, N as i64 - 1) as usize;
            if j == i {
                j = (j + 1) % N;
            }
            let c = if rand_range(0, 1) == 0 { 1 } else { -1 };
            for k in 0..N {
                u[i][k] += c * u[j][k];
            }
        }
    }
    let b = mat_mul(&u, &r);
    let (b_inv_scaled, b_det) = inverse_scaled(&b);
    (GghPublicKey { b: b.clone() }, GghSecretKey { r, b_inv_scaled, b_det })
}

/// Fixed-box error: each coordinate ±σ (the design flaw).
fn error_vector() -> Vector {
    (0..N).map(|_| if rand_range(0, 1) == 0 { SIGMA } else { -SIGMA }).collect()
}

/// Encrypt a small integer message vector `m` (length n).
pub fn ggh_encrypt(pk: &GghPublicKey, m: &[i64]) -> (Vector, Vector) {
    let e = error_vector();
    let mb = vec_mat(m, &pk.b);
    let c: Vector = (0..N).map(|i| mb[i] + e[i]).collect();
    (c, e)
}

/// Babai rounding against the good basis `R`: nearest lattice point to
/// `c`, returned as its `R`-coordinates rounded to integers.
fn babai_round(r: &Matrix, c: &[i64]) -> Vector {
    // Solve c ≈ y·R for real y, round y, return the lattice point y·R.
    // With R = DIAG·I + small, R⁻¹ ≈ (1/DIAG)·I; do it exactly over f64.
    let (r_inv_scaled, det) = inverse_scaled(r);
    // y = c · R⁻¹ = c · (r_inv_scaled / det).
    let y: Vec<f64> = (0..N)
        .map(|j| {
            let s: i64 = (0..N).map(|i| c[i] * r_inv_scaled[i][j]).sum();
            (s as f64) / (det as f64)
        })
        .collect();
    let y_round: Vec<i64> = y.iter().map(|v| v.round() as i64).collect();
    vec_mat(&y_round, r)
}

/// Decrypt to the message vector.
pub fn ggh_decrypt(sk: &GghSecretKey, c: &[i64]) -> Vector {
    // Nearest lattice point ≈ m·B; recover m = (m·B)·B⁻¹.
    let lattice_pt = babai_round(&sk.r, c);
    (0..N)
        .map(|j| {
            let s: i64 = (0..N).map(|i| lattice_pt[i] * sk.b_inv_scaled[i][j]).sum();
            // s = det · m_j, and det | s exactly for a true lattice point.
            s / sk.b_det
        })
        .collect()
}

/// Nguyen's 1999 break, the key-free part.  Every error coordinate is
/// `±σ`, and `+σ ≡ −σ ≡ σ (mod 2σ)`, so an attacker knows the *entire*
/// error residue vector `e ≡ σ·1 (mod 2σ)` with no secret key — the
/// fixed error box is far too narrow to hide anything modulo `2σ`.
/// Stripping that known residue reveals `(m·B) mod 2σ`:
///
/// ```text
/// (c − σ·1) mod 2σ = (m·B + e − σ·1) mod 2σ = (m·B) mod 2σ,
/// ```
///
/// because `e − σ·1 ∈ {0, −2σ}ⁿ ≡ 0`.  That is message information
/// leaking straight out of the ciphertext; iterating the idea (Nguyen
/// solved for `m` outright) breaks GGH entirely.  Returns
/// `(known_error_residue = σ·1, leaked_mB_mod_2σ)`.
pub fn nguyen_leak(c: &[i64]) -> (Vector, Vector) {
    let two_sigma = 2 * SIGMA;
    let known_error_residue = vec![SIGMA; c.len()];
    let leaked_mb_mod = c.iter().map(|&ci| (ci - SIGMA).rem_euclid(two_sigma)).collect();
    (known_error_residue, leaked_mb_mod)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_message() -> Vec<i64> {
        (0..N).map(|_| rand_range(-3, 3)).collect()
    }

    #[test]
    fn good_basis_is_near_orthogonal() {
        // The private basis is diagonally dominant: |R_ii| ≫ Σ|R_ij|.
        let (_, sk) = ggh_keygen();
        for i in 0..N {
            let off: i64 = (0..N).filter(|&j| j != i).map(|j| sk.r[i][j].abs()).sum();
            assert!(sk.r[i][i].abs() > off, "row {i} not diagonally dominant");
        }
    }

    #[test]
    fn encrypt_decrypt_roundtrip() {
        let (pk, sk) = ggh_keygen();
        for _ in 0..20 {
            let m = random_message();
            let (c, _e) = ggh_encrypt(&pk, &m);
            assert_eq!(ggh_decrypt(&sk, &c), m);
        }
    }

    #[test]
    fn public_and_private_span_same_lattice() {
        // det(B) = ±det(R): U is unimodular, so both bases share Λ.
        let (pk, sk) = ggh_keygen();
        assert_eq!(determinant(&pk.b).abs(), determinant(&sk.r).abs());
    }

    #[test]
    fn nguyen_leak_exposes_message_information() {
        // The break: from a real ciphertext, with no secret key, the
        // attacker (a) knows every error coordinate ≡ σ (mod 2σ) and
        // (b) thereby strips it to learn m·B mod 2σ — information the
        // ciphertext was supposed to hide.
        let (pk, _sk) = ggh_keygen();
        let m = random_message();
        let (c, e) = ggh_encrypt(&pk, &m);

        // (a) The fixed error box collapses mod 2σ to the known σ·1.
        for &ei in &e {
            assert_eq!(ei.rem_euclid(2 * SIGMA), SIGMA);
        }

        // (b) The key-free strip reveals exactly (m·B) mod 2σ.
        let (known_residue, leaked_mb_mod) = nguyen_leak(&c);
        assert_eq!(known_residue, vec![SIGMA; N]);
        let mb = vec_mat(&m, &pk.b);
        let mb_mod_true: Vec<i64> = mb.iter().map(|&x| x.rem_euclid(2 * SIGMA)).collect();
        assert_eq!(leaked_mb_mod, mb_mod_true, "attacker learns m·B mod 2σ without the key");
    }
}
