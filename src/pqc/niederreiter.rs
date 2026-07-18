//! **Niederreiter** — the syndrome-based dual of McEliece
//! (Niederreiter 1986).
//!
//! # McEliece vs. Niederreiter
//! Both rest on the hardness of decoding a random-looking linear code,
//! and they are provably equivalent in security.  They differ in what
//! plays the role of the message:
//!
//! | | McEliece (`pqc::mceliece`) | Niederreiter (this module) |
//! |---|---|---|
//! | public key | scrambled **generator** `G'` | scrambled **parity check** `H'` |
//! | message | `k` bits, placed in a codeword | a **weight-`t` error pattern** |
//! | ciphertext | `m·G' + e` (length `n`) | syndrome `H'·eᵀ` (length `n−k`) |
//! | decrypt | strip scramble, decode, read `m` | syndrome-decode to recover `e` |
//!
//! Niederreiter's ciphertext is just the syndrome — shorter than
//! McEliece's — and the message is encoded as the *choice* of which `t`
//! of the `n` coordinates carry the error.  This is exactly the
//! encoding Classic McEliece (the NIST KEM) uses internally.
//!
//! # This implementation
//! Reuses the Goppa-code machinery of `pqc::mceliece` (parity-check
//! construction, systematic form, and Patterson decoding) with the same
//! toy parameters `m = 6, n = 32, t = 3`.  The message is a `t`-subset
//! of the `n` coordinates, encoded via the combinatorial number system
//! (so the plaintext space is `C(32, 3) = 4960`).  Decryption rebuilds
//! a word with the ciphertext syndrome, runs Patterson decoding to peel
//! off the weight-`t` error, and reads back the subset.
//!
//! As with the sibling `pqc::mceliece`, these toy parameters use a fixed
//! Goppa support (`0…n−1`), so the code structure is **not** hidden and
//! the scheme offers no real security — the value here is the exact
//! syndrome-encryption structure and its equivalence to McEliece, not
//! confidentiality.  Not constant-time; see SECURITY.md.

use super::mceliece::{
    expand_to_binary, parity_check_gf, patterson_decode, random_irreducible_poly, systematic_form,
    GfTables, M, N, T,
};

/// Parity-check rows = `n − k = m·t`.
pub const SYND_LEN: usize = M * T; // 18

/// Public key: the systematic parity-check matrix `H = [I | B]`
/// (`(n−k) × n` over F₂), with the code's structure hidden.
#[derive(Clone)]
pub struct NiederreiterPublicKey {
    pub h: Vec<Vec<u8>>,
}

/// Secret key: the Goppa trapdoor (support order + polynomial + field).
#[derive(Clone)]
pub struct NiederreiterSecretKey {
    support: Vec<u32>,
    g: Vec<u32>,
    gf: GfTables,
}

pub fn niederreiter_keygen() -> (NiederreiterPublicKey, NiederreiterSecretKey) {
    let gf = GfTables::new();
    loop {
        let g = random_irreducible_poly(T, &gf);
        let support: Vec<u32> = (0..N as u32).collect();
        let h_gf = parity_check_gf(&support, &g, &gf);
        let h_bin = expand_to_binary(&h_gf);
        let (h_sys, perm) = match systematic_form(&h_bin) {
            Some(p) => p,
            None => continue, // singular parity check: retry
        };
        // Reorder the support to match the systematic column order, so
        // Patterson decoding later sees columns in the public order.
        let support_sys: Vec<u32> = perm.iter().map(|&i| support[i]).collect();
        return (
            NiederreiterPublicKey { h: h_sys },
            NiederreiterSecretKey { support: support_sys, g, gf },
        );
    }
}

// ── Combinatorial number system: rank/unrank of t-subsets of [0, n) ──────────

fn binom(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut num = 1u64;
    for i in 0..k {
        num = num * (n - i) as u64 / (i as u64 + 1);
    }
    num
}

/// Number of distinct plaintexts = C(n, t).
pub fn message_space() -> u64 {
    binom(N, T)
}

/// Unrank: message index → sorted `t`-subset of `[0, n)`, using the
/// combinatorial number system (combinadics).  Elements are selected in
/// decreasing order: for `i = t … 1`, take the largest `c` with
/// `C(c, i) ≤ idx`, then subtract.
fn unrank(mut idx: u64) -> Vec<usize> {
    let mut subset = vec![0usize; T];
    for i in (1..=T).rev() {
        let mut c = i - 1; // smallest c with C(c, i) potentially > 0
        while binom(c + 1, i) <= idx {
            c += 1;
        }
        idx -= binom(c, i);
        subset[i - 1] = c; // fills largest→smallest as i descends
    }
    subset // already ascending by construction
}

/// Rank: sorted `t`-subset → message index (inverse of `unrank`).
/// The `j`-th smallest element contributes `C(s[j], j+1)`.
fn rank(subset: &[usize]) -> u64 {
    let mut sorted = subset.to_vec();
    sorted.sort_unstable();
    sorted.iter().enumerate().map(|(j, &c)| binom(c, j + 1)).sum()
}

/// Encode a message index as a weight-`t` error vector in F₂ⁿ.
fn message_to_error(idx: u64) -> Vec<u8> {
    let subset = unrank(idx);
    let mut e = vec![0u8; N];
    for &p in &subset {
        e[p] = 1;
    }
    e
}

/// Matrix-vector product `H·eᵀ` over F₂.
fn syndrome(h: &[Vec<u8>], e: &[u8]) -> Vec<u8> {
    h.iter()
        .map(|row| row.iter().zip(e).fold(0u8, |acc, (&hij, &ej)| acc ^ (hij & ej)))
        .collect()
}

/// Encrypt a message index (`< message_space()`), returning its
/// syndrome (length `n−k`).
pub fn niederreiter_encrypt(pk: &NiederreiterPublicKey, idx: u64) -> Vec<u8> {
    assert!(idx < message_space(), "message out of range");
    syndrome(&pk.h, &message_to_error(idx))
}

/// Decrypt a syndrome back to the message index, or `None` if it does
/// not decode to a weight-`t` error.
pub fn niederreiter_decrypt(sk: &NiederreiterSecretKey, synd: &[u8]) -> Option<u64> {
    if synd.len() != SYND_LEN {
        return None;
    }
    // Build any word r0 with H·r0ᵀ = synd: since H = [I | B] is
    // systematic, r0 = (synd ‖ 0) works.  r0 = e + codeword, so it has
    // the same syndrome as the planted error e.
    let mut r0 = vec![0u8; N];
    r0[..SYND_LEN].copy_from_slice(synd);
    // Patterson decoding recovers the weight-t error e.
    let (_codeword, error) = patterson_decode(&r0, &sk.support, &sk.g, &sk.gf);
    let positions: Vec<usize> = (0..N).filter(|&i| error[i] == 1).collect();
    if positions.len() != T {
        return None; // not a valid weight-t ciphertext
    }
    Some(rank(&positions))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combinatorial_rank_unrank_roundtrip() {
        for idx in 0..message_space() {
            let subset = unrank(idx);
            assert_eq!(subset.len(), T);
            assert_eq!(rank(&subset), idx, "rank∘unrank must be identity");
        }
    }

    #[test]
    fn ciphertext_is_a_syndrome() {
        let (pk, _) = niederreiter_keygen();
        let ct = niederreiter_encrypt(&pk, 42);
        assert_eq!(ct.len(), SYND_LEN);
    }

    #[test]
    fn encrypt_decrypt_roundtrip() {
        let (pk, sk) = niederreiter_keygen();
        // Sample messages across the space.
        for idx in [0u64, 1, 7, 100, 1000, message_space() - 1] {
            let ct = niederreiter_encrypt(&pk, idx);
            assert_eq!(niederreiter_decrypt(&sk, &ct), Some(idx), "idx {idx}");
        }
    }

    #[test]
    fn many_random_messages_roundtrip() {
        use crate::utils::random::random_bytes;
        let (pk, sk) = niederreiter_keygen();
        for _ in 0..30 {
            let mut b = [0u8; 8];
            random_bytes(&mut b);
            let idx = u64::from_le_bytes(b) % message_space();
            let ct = niederreiter_encrypt(&pk, idx);
            assert_eq!(niederreiter_decrypt(&sk, &ct), Some(idx));
        }
    }

    #[test]
    fn tampered_syndrome_changes_message() {
        // Flipping a syndrome bit yields a different coset, so decoding
        // (when it still succeeds) returns a different message — never
        // the original.  This exercises that decryption genuinely tracks
        // the syndrome rather than ignoring it.
        let (pk, sk) = niederreiter_keygen();
        let idx = 321u64;
        let mut ct = niederreiter_encrypt(&pk, idx);
        ct[0] ^= 1;
        assert_ne!(niederreiter_decrypt(&sk, &ct), Some(idx));
    }

    #[test]
    fn malformed_syndrome_length_rejected() {
        let (_, sk) = niederreiter_keygen();
        assert_eq!(niederreiter_decrypt(&sk, &vec![0u8; SYND_LEN - 1]), None);
    }
}
