//! **Rank-metric cryptography** — Gabidulin-code public-key encryption,
//! the core of RQC/ROLLO (NIST round-2 rank-based candidates).
//!
//! # A different hard problem
//! McEliece/Niederreiter (`pqc::mceliece`, `pqc::niederreiter`) use the
//! **Hamming** metric — weight = number of nonzero coordinates.  Rank-
//! metric schemes instead measure a vector `x ∈ F_{qᵐ}ⁿ` by its **rank
//! weight**: expand each coordinate over an `F_q`-basis to get an `m×n`
//! matrix, and take that matrix's rank.  The hard problem is **rank
//! syndrome decoding** — finding a low-*rank* error with a given
//! syndrome.  It has no known efficient (even quantum) attack and gives
//! much smaller keys than Hamming-metric codes, which is why NIST ran
//! RQC and ROLLO as round-2 candidates.
//!
//! # Gabidulin codes (the rank-metric analog of Reed–Solomon)
//! A codeword is `(f(g₁), …, f(gₙ))` where `g₁…gₙ ∈ F_{qᵐ}` are
//! `F_q`-linearly independent and `f` is a **linearized polynomial**
//! `f(x) = Σ_{i<k} fᵢ·x^{qⁱ}` (additive: `f(a+b)=f(a)+f(b)`).  Like
//! Reed–Solomon it corrects up to `⌊(n−k)/2⌋` errors, but in the rank
//! metric.  This is the trapdoor code behind RQC/GPT.
//!
//! # This implementation (GPT-style encryption)
//! Toy parameters `q = 2, m = 8` (so `F_{qᵐ} = F₂₈ = GF(256)`),
//! `n = 4, k = 2`, correcting `t = 1` rank error.
//!
//! - **KeyGen**: secret `F_q`-independent basis `g` and invertible
//!   scrambler `S`; public generator `G_pub = S·G` where `G` is the
//!   Gabidulin (Moore) generator.  `G_pub` spans the same code as `G`
//!   but hides the evaluation-point structure.
//! - **Encrypt** `m ∈ F₂₈²`: `c = m·G_pub ⊕ e` with a fresh rank-1
//!   error `e`.
//! - **Decrypt**: rank-1 decode `c` to the codeword, read off `m·S`
//!   with the secret basis, then unscramble by `S⁻¹`.
//!
//! Decoding here is an exhaustive rank-1 search (feasible only at these
//! toy sizes — its infeasibility at real parameters is the security
//! claim).  Not constant-time; see SECURITY.md.

use crate::utils::random::random_bytes;

/// Extension degree: F_{2^M} = GF(256).
pub const M: usize = 8;
/// Code length.
pub const N: usize = 4;
/// Code dimension.
pub const K: usize = 2;
/// Correctable rank-error weight `⌊(n−k)/2⌋`.
pub const T: usize = (N - K) / 2;

// ── GF(256) (AES polynomial x⁸+x⁴+x³+x+1) ────────────────────────────────────

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

/// Frobenius: `x^{2^i}` = square `i` times.
fn frob(mut x: u8, i: usize) -> u8 {
    for _ in 0..i {
        x = gf_mul(x, x);
    }
    x
}

// ── Rank weight ───────────────────────────────────────────────────────────────

/// Rank weight of `x ∈ GF(256)ⁿ`: the F₂-rank of the `m×n` bit matrix
/// whose columns are the bit-expansions of the coordinates.
pub fn rank_weight(x: &[u8]) -> usize {
    // Build m×n bit matrix (row r, col j = bit r of x[j]) and row-reduce
    // over F₂.
    let mut rows: Vec<u16> = (0..M)
        .map(|r| {
            let mut bits = 0u16;
            for (j, &xj) in x.iter().enumerate() {
                if (xj >> r) & 1 == 1 {
                    bits |= 1 << j;
                }
            }
            bits
        })
        .collect();
    let mut rank = 0;
    let mut pivot_col = 0;
    while pivot_col < x.len() && rank < rows.len() {
        if let Some(sel) = (rank..rows.len()).find(|&i| (rows[i] >> pivot_col) & 1 == 1) {
            rows.swap(rank, sel);
            for i in 0..rows.len() {
                if i != rank && (rows[i] >> pivot_col) & 1 == 1 {
                    rows[i] ^= rows[rank];
                }
            }
            rank += 1;
        }
        pivot_col += 1;
    }
    rank
}

// ── Gabidulin encoding ────────────────────────────────────────────────────────

/// Encode `f ∈ GF(256)^k` (linearized-polynomial coefficients) to a
/// codeword `c_j = Σ_i f_i · g_j^{2^i}`.
fn gabidulin_encode(f: &[u8], g: &[u8]) -> Vec<u8> {
    (0..N)
        .map(|j| {
            let mut acc = 0u8;
            for i in 0..K {
                acc ^= gf_mul(f[i], frob(g[j], i));
            }
            acc
        })
        .collect()
}

/// Recover `f` from a *clean* Gabidulin codeword using its first `k`
/// coordinates: solve the `k×k` Moore system `M·f = c[..k]`.
fn gabidulin_decode_clean(c: &[u8], g: &[u8]) -> Option<Vec<u8>> {
    // Moore matrix rows: M[j][i] = g_j^{2^i}, for j, i < k.
    let mut mat = vec![vec![0u8; K]; K];
    let mut rhs = vec![0u8; K];
    for j in 0..K {
        for i in 0..K {
            mat[j][i] = frob(g[j], i);
        }
        rhs[j] = c[j];
    }
    // Gaussian elimination over GF(256).
    for col in 0..K {
        let piv = (col..K).find(|&r| mat[r][col] != 0)?;
        mat.swap(col, piv);
        rhs.swap(col, piv);
        let inv = gf_inv(mat[col][col]);
        for x in mat[col].iter_mut() {
            *x = gf_mul(*x, inv);
        }
        rhs[col] = gf_mul(rhs[col], inv);
        for r in 0..K {
            if r != col && mat[r][col] != 0 {
                let f = mat[r][col];
                for i in 0..K {
                    mat[r][i] ^= gf_mul(f, mat[col][i]);
                }
                rhs[r] ^= gf_mul(f, rhs[col]);
            }
        }
    }
    let f = rhs;
    // Verify the recovered f reproduces the whole codeword.
    if gabidulin_encode(&f, g) == c {
        Some(f)
    } else {
        None
    }
}

/// Rank-`t` bounded-distance decode: strip a rank-≤1 error by exhaustive
/// search, then decode the clean codeword.
fn gabidulin_decode(r: &[u8], g: &[u8]) -> Option<Vec<u8>> {
    // Zero-error case.
    if let Some(f) = gabidulin_decode_clean(r, g) {
        return Some(f);
    }
    // Rank-1 errors: e_j = b_j · a for a ∈ GF(256)*, b ∈ F₂ⁿ \ {0}.
    for a in 1u16..256 {
        for b in 1u8..(1 << N) {
            let e: Vec<u8> = (0..N).map(|j| if (b >> j) & 1 == 1 { a as u8 } else { 0 }).collect();
            let candidate: Vec<u8> = (0..N).map(|j| r[j] ^ e[j]).collect();
            if let Some(f) = gabidulin_decode_clean(&candidate, g) {
                return Some(f);
            }
        }
    }
    None
}

// ── Keys / encrypt / decrypt ──────────────────────────────────────────────────

#[derive(Clone)]
pub struct RqcPublicKey {
    /// G_pub = S·G, a k×n generator over GF(256).
    pub g_pub: Vec<Vec<u8>>,
}

#[derive(Clone)]
pub struct RqcSecretKey {
    g: Vec<u8>,          // secret evaluation basis (length n)
    s_inv: Vec<Vec<u8>>, // k×k inverse scrambler
}

fn rand_nonzero() -> u8 {
    loop {
        let mut b = [0u8; 1];
        random_bytes(&mut b);
        if b[0] != 0 {
            return b[0];
        }
    }
}

/// Random F₂-linearly-independent basis of length `n` in GF(256).
fn random_independent_basis() -> Vec<u8> {
    loop {
        let mut g = vec![0u8; N];
        random_bytes(&mut g);
        if g.iter().all(|&x| x != 0) && rank_weight(&g) == N {
            return g;
        }
    }
}

/// Random invertible k×k matrix over GF(256) with its inverse.
fn random_invertible() -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    loop {
        let mut m = vec![vec![0u8; K]; K];
        for row in m.iter_mut() {
            for x in row.iter_mut() {
                *x = rand_nonzero();
            }
        }
        if let Some(inv) = invert(&m) {
            return (m, inv);
        }
    }
}

fn invert(m: &[Vec<u8>]) -> Option<Vec<Vec<u8>>> {
    let n = m.len();
    let mut a: Vec<Vec<u8>> = m.to_vec();
    let mut inv = vec![vec![0u8; n]; n];
    for i in 0..n {
        inv[i][i] = 1;
    }
    for col in 0..n {
        let piv = (col..n).find(|&r| a[r][col] != 0)?;
        a.swap(col, piv);
        inv.swap(col, piv);
        let iv = gf_inv(a[col][col]);
        for j in 0..n {
            a[col][j] = gf_mul(a[col][j], iv);
            inv[col][j] = gf_mul(inv[col][j], iv);
        }
        for r in 0..n {
            if r != col && a[r][col] != 0 {
                let f = a[r][col];
                for j in 0..n {
                    a[r][j] ^= gf_mul(f, a[col][j]);
                    inv[r][j] ^= gf_mul(f, inv[col][j]);
                }
            }
        }
    }
    Some(inv)
}

fn vec_mat(v: &[u8], m: &[Vec<u8>]) -> Vec<u8> {
    let cols = m[0].len();
    (0..cols)
        .map(|j| (0..v.len()).fold(0u8, |acc, i| acc ^ gf_mul(v[i], m[i][j])))
        .collect()
}

pub fn rqc_keygen() -> (RqcPublicKey, RqcSecretKey) {
    let g = random_independent_basis();
    let (s, s_inv) = random_invertible();
    // Gabidulin generator G (k×n): G[i][j] = g_j^{2^i}.
    let gen: Vec<Vec<u8>> =
        (0..K).map(|i| (0..N).map(|j| frob(g[j], i)).collect()).collect();
    // Public generator G_pub = S·G.
    let g_pub: Vec<Vec<u8>> = (0..K)
        .map(|i| {
            (0..N)
                .map(|j| (0..K).fold(0u8, |acc, l| acc ^ gf_mul(s[i][l], gen[l][j])))
                .collect()
        })
        .collect();
    (RqcPublicKey { g_pub }, RqcSecretKey { g, s_inv })
}

/// A fresh rank-1 error vector.
fn rank1_error() -> Vec<u8> {
    let a = rand_nonzero();
    loop {
        let mut b = [0u8; 1];
        random_bytes(&mut b);
        let pattern = b[0] & ((1 << N) - 1);
        if pattern != 0 {
            return (0..N).map(|j| if (pattern >> j) & 1 == 1 { a } else { 0 }).collect();
        }
    }
}

pub fn rqc_encrypt(pk: &RqcPublicKey, msg: &[u8; K]) -> Vec<u8> {
    let codeword = vec_mat(msg, &pk.g_pub);
    let e = rank1_error();
    (0..N).map(|j| codeword[j] ^ e[j]).collect()
}

pub fn rqc_decrypt(sk: &RqcSecretKey, ct: &[u8]) -> Option<[u8; K]> {
    // Decode to the Gabidulin message f = m·S, then unscramble.
    let f = gabidulin_decode(ct, &sk.g)?;
    let m = vec_mat(&f, &sk.s_inv);
    m.try_into().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_weight_basics() {
        // A single nonzero coordinate has rank 1.
        assert_eq!(rank_weight(&[5, 0, 0, 0]), 1);
        // Two coordinates that are F₂-multiples of each other: still rank 1.
        assert_eq!(rank_weight(&[5, 5, 0, 0]), 1);
        // The zero vector has rank 0.
        assert_eq!(rank_weight(&[0, 0, 0, 0]), 0);
        // An independent basis has full rank n.
        let g = random_independent_basis();
        assert_eq!(rank_weight(&g), N);
    }

    #[test]
    fn rank1_errors_have_rank_one() {
        for _ in 0..20 {
            let e = rank1_error();
            assert_eq!(rank_weight(&e), 1);
        }
    }

    #[test]
    fn clean_codeword_decodes() {
        let g = random_independent_basis();
        let f = [0x1f, 0xa3];
        let c = gabidulin_encode(&f, &g);
        assert_eq!(gabidulin_decode_clean(&c, &g), Some(f.to_vec()));
    }

    #[test]
    fn encrypt_decrypt_roundtrip() {
        let (pk, sk) = rqc_keygen();
        for _ in 0..30 {
            let mut m = [0u8; K];
            random_bytes(&mut m);
            let ct = rqc_encrypt(&pk, &m);
            assert_eq!(rqc_decrypt(&sk, &ct), Some(m));
        }
    }

    #[test]
    fn ciphertext_carries_a_rank1_error() {
        let (pk, sk) = rqc_keygen();
        let m = [0x42, 0x99];
        let ct = rqc_encrypt(&pk, &m);
        // Ciphertext minus the true codeword is exactly a rank-1 vector.
        let codeword = vec_mat(&m, &pk.g_pub);
        let e: Vec<u8> = (0..N).map(|j| ct[j] ^ codeword[j]).collect();
        assert_eq!(rank_weight(&e), 1);
        assert_eq!(rqc_decrypt(&sk, &ct), Some(m));
    }

    #[test]
    fn wrong_key_fails() {
        let (pk, _) = rqc_keygen();
        let (_, sk_other) = rqc_keygen();
        let m = [0x11, 0x22];
        let ct = rqc_encrypt(&pk, &m);
        assert_ne!(rqc_decrypt(&sk_other, &ct), Some(m));
    }
}
