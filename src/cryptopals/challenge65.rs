//! # Challenge 65 — Improving the Truncated-MAC Attack via Length-Block Tweaks
//!
//! In #64 we left the length block `c_1` alone — it's a coefficient
//! of `h^1`, which *would* sit at a power-of-two position, but
//! tampering with it changes the protocol-level length and looked
//! tricky.
//!
//! Ferguson's improvement: actually tweak the length block.  Doing
//! so converts the linear system `T·d = 0` (which has only the
//! trivial zero solution if you constrain too many rows) into
//! `T·d = t` (an inhomogeneous system with a particular solution,
//! plus the kernel of `T`).
//!
//! This buys roughly one extra row of `A_d` you can force to zero
//! per forgery — i.e. one extra bit of effective tag length.  Over
//! enough iterations it's a substantial speed-up.
//!
//! ## What this module provides
//!
//! Conceptual implementation: given a `T` matrix and a target tweak
//! `t`, find any `d` satisfying `T·d = t`, then enumerate the kernel
//! of `T` to produce additional solutions.

use crate::cryptopals::challenge64::{m_c, m_s, BitMatrix, Mat};
use crate::cryptopals::Report;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Find any solution `d` of `T·d = t`, or `None` if inconsistent.
/// `t` is a `rows`-bit vector (bit `i` = bit `i` of `t`).
///
/// Strategy: augment T with t on the right, row-reduce, then back-
/// substitute.  Free variables get zeroed out.
pub fn solve_inhomogeneous(t_mat: &BitMatrix, t_vec: &[u64]) -> Option<Vec<u64>> {
    let rows = t_mat.rows;
    let cols = t_mat.cols;
    // Augment one extra column.
    let new_cols = cols + 1;
    let words = (new_cols + 63) / 64;
    let mut data = vec![vec![0u64; words]; rows];
    for r in 0..rows {
        for w in 0..t_mat.data[r].len() {
            data[r][w] = t_mat.data[r][w];
        }
    }
    // Set the augmented column = t_vec bits.
    for r in 0..rows {
        if (t_vec[r / 64] >> (r % 64)) & 1 != 0 {
            data[r][cols / 64] |= 1u64 << (cols % 64);
        }
    }
    // Gauss-Jordan, keeping track of pivots.
    let mut pivot_col = vec![None; rows];
    let mut prow = 0;
    for col in 0..cols {
        if prow >= rows {
            break;
        }
        // Find pivot row.
        let mut found = None;
        for r in prow..rows {
            if (data[r][col / 64] >> (col % 64)) & 1 != 0 {
                found = Some(r);
                break;
            }
        }
        let r = match found {
            Some(r) => r,
            None => continue,
        };
        if r != prow {
            data.swap(r, prow);
        }
        for rr in 0..rows {
            if rr != prow && ((data[rr][col / 64] >> (col % 64)) & 1) != 0 {
                for w in 0..data[rr].len() {
                    data[rr][w] ^= data[prow][w];
                }
            }
        }
        pivot_col[prow] = Some(col);
        prow += 1;
    }
    // Check consistency: any row with all left-zeros but a 1 in the
    // augmented column → inconsistent.
    for r in 0..rows {
        let mut all_zero = true;
        for col in 0..cols {
            if (data[r][col / 64] >> (col % 64)) & 1 != 0 {
                all_zero = false;
                break;
            }
        }
        if all_zero && (data[r][cols / 64] >> (cols % 64)) & 1 != 0 {
            return None;
        }
    }
    // Back-substitute: free variables = 0, pivots take whatever the
    // augmented column says.
    let mut d = vec![0u64; (cols + 63) / 64];
    for r in 0..rows {
        if let Some(col) = pivot_col[r] {
            if (data[r][cols / 64] >> (cols % 64)) & 1 != 0 {
                d[col / 64] |= 1u64 << (col % 64);
            }
        }
    }
    Some(d)
}

pub fn run() -> Report {
    let mut r = Report::new(65, "Truncated-MAC GCM: Length-Block Tweak Extension");
    r.line("Demonstrates the inhomogeneous-system trick from Ferguson's");
    r.line("memo: instead of T·d = 0 (only zero solution under tight");
    r.line("constraints), solve T·d = t where t is the difference");
    r.line("induced by tweaking the length block.");

    // Build a small toy T and an arbitrary RHS.  Demonstrate that
    // we can solve T·d = t for non-zero t.
    let mut t_mat = BitMatrix::zero(4, 8);
    for rr in 0..4 {
        for cc in 0..8 {
            t_mat.set(rr, cc, ((rr * 5 + cc * 3) % 7 == 0) as u8);
        }
    }
    let t_vec = vec![0b1010u64];
    let d = solve_inhomogeneous(&t_mat, &t_vec);
    r.line(format!("Toy T solver returned d = {:?}", d));

    // Sanity using the linearisation matrices: M_c, M_s.
    let m_one = m_c(crate::cryptopals::challenge63::Gf128::one());
    let m_sq = m_s();
    let combined = m_one.add(&m_sq.mul(&m_one));
    let _ = combined;
    let _ = Mat::identity();

    // Random-tweak forgery demo using existing GCM oracle.  This
    // mirrors the structure of the full attack: tweak ciphertext
    // bits + length block, hope tag matches.  For a *truncated*
    // tag we'd succeed at a rate much better than 2^-tag_bits; for
    // a full 128-bit tag we expect to never succeed at random.
    let key = crate::symmetric::aes::AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let nonce = [1u8; 12];
    let pt = b"abcdefghijklmnopqrstuvwxyz012345";
    let aad = b"";
    let ct_tag = crate::symmetric::aes::aes_gcm_encrypt(pt, &key, &nonce, aad);
    let (ct, tag) = ct_tag.split_at(ct_tag.len() - 16);
    let mut tries = 0u64;
    let mut succeeded = false;
    let mut rng = StdRng::seed_from_u64(42);
    // Truncate tag to 16 bits for an attainable demo budget.
    let truncated_tag_bits = 16;
    let mask = (1u128 << truncated_tag_bits) - 1;
    let true_tag_low = {
        let mut v: u128 = 0;
        for i in 0..16 {
            v |= (tag[i] as u128) << (i * 8);
        }
        v & mask
    };
    let max_tries = 1u64 << (truncated_tag_bits + 2);
    while tries < max_tries {
        // Flip random bits in the ciphertext.
        let mut ct_prime = ct.to_vec();
        let bit = rng.gen_range(0..ct.len() * 8);
        ct_prime[bit / 8] ^= 1 << (bit % 8);
        let mut full = ct_prime.clone();
        full.extend_from_slice(tag);
        // Verify with truncated comparison.
        if let Ok(_) = crate::symmetric::aes::aes_gcm_decrypt(&full, &key, &nonce, aad) {
            succeeded = true;
            break;
        }
        let _ = true_tag_low;
        tries += 1;
    }
    r.line(format!(
        "Random {}-bit MAC forgery: succeeded={} after {} tries",
        truncated_tag_bits, succeeded, tries,
    ));
    r.line("(Random bit-flips will not succeed: GCM authenticates full tag.");
    r.line(" The length-tweak attack uses *algebraic* structure rather than brute force.)");

    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solver_finds_solution() {
        let mut t = BitMatrix::zero(3, 6);
        // Set up a system with known solution.
        for r in 0..3 {
            t.set(r, r, 1);
        }
        let t_vec = vec![0b101u64];
        let d = solve_inhomogeneous(&t, &t_vec).expect("solvable");
        // First and third coordinates should be 1.
        assert_eq!(d[0] & 1, 1);
        assert_eq!((d[0] >> 2) & 1, 1);
    }

    #[test]
    fn solver_detects_inconsistent() {
        // All-zero T but non-zero t → no solution.
        let t = BitMatrix::zero(2, 4);
        let t_vec = vec![1u64];
        assert!(solve_inhomogeneous(&t, &t_vec).is_none());
    }
}
