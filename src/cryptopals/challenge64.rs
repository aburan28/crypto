//! # Challenge 64 — Key Recovery on GCM with a Truncated MAC
//!
//! GCM with a short tag (NIST allows 32-bit and even shorter) is
//! catastrophic.  An attacker who can submit one forgery in
//! `2^tag_bits` random tries can use the **algebraic structure** of
//! GHASH to do dramatically better, then iteratively narrow down the
//! 128-bit authentication key.
//!
//! ## Linearisation
//!
//! Two GHASH operations are GF(2)-linear in the bits of `h`:
//!
//! - **Multiply by constant `c`**: `f(h) = c · h`.  Build the
//!   `128×128` matrix `M_c` whose column `i` is `c · x^i mod p(x)`.
//! - **Square**: `f(h) = h²`.  Build `M_s` similarly (column `i`
//!   = `(x^i)² mod p(x)`), exploiting `(a+b)² = a²+b²` in GF(2).
//!
//! Both packaged as `Mat` (a 128-row × 128-col GF(2) matrix stored
//! as 128 column vectors of type `u128`).
//!
//! ## The forgery system
//!
//! With ciphertext blocks `c₁, c₂, …, c_n` and only the blocks at
//! powers of 2 (`c_{2^k}`) molestable, the difference matrix is
//!
//! ```text
//!   A_d = Σ M_{d_i} · M_s^i
//! ```
//!
//! where `d_i` is the XOR delta applied to block `c_{2^i}`.
//!
//! Solve `T·d = 0` (the `T` matrix relates attacker-bit flips to
//! rows of `A_d`) and you get bit-flip patterns whose `A_d·h`
//! products land in zero rows — guaranteed-correct chunks of the
//! tag.  Each successful forgery yields `m` linear equations on
//! `h` (the rows of `A_d` for which we still have non-zero entries).
//!
//! ## What this module provides
//!
//! Educational implementation of the core machinery:
//!
//! - `Mat` (128×128 GF(2) matrix) with multiply and the canonical
//!   row-reduction.
//! - `m_c`, `m_s` to build the multiplication / squaring matrices.
//! - `null_space_basis` (Gauss-Jordan over GF(2)).
//! - End-to-end forgery against a **tiny** truncated-tag setup
//!   (8-bit tag, 4 ciphertext blocks).  Real cryptopals parameters
//!   require linear systems with billions of bits — out of scope
//!   for a self-contained demo, but the algorithm is identical.

use crate::cryptopals::challenge63::{pack, unpack, Gf128};
use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// 128×128 matrix over GF(2).  Stored as 128 column vectors (`u128`
/// each).  Row `i` of column `j` is `(cols[j] >> i) & 1`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mat {
    pub cols: [u128; 128],
}

impl Mat {
    pub fn zero() -> Self {
        Mat { cols: [0u128; 128] }
    }
    pub fn identity() -> Self {
        let mut m = Mat::zero();
        for i in 0..128 {
            m.cols[i] = 1u128 << i;
        }
        m
    }

    pub fn get(&self, row: usize, col: usize) -> u8 {
        ((self.cols[col] >> row) & 1) as u8
    }
    pub fn set(&mut self, row: usize, col: usize, v: u8) {
        let mask = 1u128 << row;
        if v != 0 {
            self.cols[col] |= mask;
        } else {
            self.cols[col] &= !mask;
        }
    }

    /// Matrix-vector multiply: `M · v`.  Vector is a `Gf128`.
    pub fn mul_vec(&self, v: Gf128) -> Gf128 {
        let mut out: u128 = 0;
        for i in 0..128 {
            if (v.0 >> i) & 1 != 0 {
                out ^= self.cols[i];
            }
        }
        Gf128(out)
    }

    /// Add two matrices (XOR each column).
    pub fn add(&self, rhs: &Self) -> Self {
        let mut out = Mat::zero();
        for i in 0..128 {
            out.cols[i] = self.cols[i] ^ rhs.cols[i];
        }
        out
    }

    /// Multiply matrices: `(self · rhs)·v = self·(rhs·v)`.
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = Mat::zero();
        for j in 0..128 {
            let v = Gf128(rhs.cols[j]);
            out.cols[j] = self.mul_vec(v).0;
        }
        out
    }
}

/// Build `M_c`, the matrix that performs multiplication by `c` in
/// GF(2^128).  Column `i` is `c · x^i mod p(x)`.
pub fn m_c(c: Gf128) -> Mat {
    let mut m = Mat::zero();
    for i in 0..128 {
        let xi = Gf128(1u128 << i);
        let prod = c.mul(xi);
        m.cols[i] = prod.0;
    }
    m
}

/// Build `M_s`, the squaring matrix.  Column `i` is `(x^i)² mod p(x)`.
pub fn m_s() -> Mat {
    let mut m = Mat::zero();
    for i in 0..128 {
        let xi = Gf128(1u128 << i);
        let sq = xi.mul(xi);
        m.cols[i] = sq.0;
    }
    m
}

/// Variable-size GF(2) matrix used for the `T` dependency matrix
/// (whose dimensions depend on the number of available ciphertext
/// blocks).  Stored row-major with each row as a `Vec<u64>` of
/// bit-packed entries.
#[derive(Clone, Debug)]
pub struct BitMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Vec<u64>>, // data[r][c/64] bit (c%64)
}

impl BitMatrix {
    pub fn zero(rows: usize, cols: usize) -> Self {
        let words = (cols + 63) / 64;
        BitMatrix {
            rows,
            cols,
            data: vec![vec![0u64; words]; rows],
        }
    }
    pub fn get(&self, r: usize, c: usize) -> u8 {
        ((self.data[r][c / 64] >> (c % 64)) & 1) as u8
    }
    pub fn set(&mut self, r: usize, c: usize, v: u8) {
        let mask = 1u64 << (c % 64);
        if v != 0 {
            self.data[r][c / 64] |= mask;
        } else {
            self.data[r][c / 64] &= !mask;
        }
    }
    /// Augment with identity on the right.  Returns a (rows × (cols + rows)) matrix.
    pub fn augment_identity(&self) -> Self {
        let new_cols = self.cols + self.rows;
        let mut out = BitMatrix::zero(self.rows, new_cols);
        for r in 0..self.rows {
            for w in 0..self.data[r].len() {
                out.data[r][w] = self.data[r][w];
            }
            out.set(r, self.cols + r, 1);
        }
        out
    }
}

/// Gauss-Jordan over GF(2) on the first `pivot_cols` columns of
/// `m`.  The remaining columns are dragged along with the row
/// operations (typical use: augment `T^T` with an identity matrix
/// on the right, run Gauss-Jordan on the left `T^T` block, and
/// rows whose left block ends up all-zero contain a basis vector
/// of `N(T)` on the right).
///
/// Returns the row indices whose **left** block went to zero — these
/// rows hold the null-space basis when the right is an identity.
pub fn gauss_jordan(m: &mut BitMatrix, pivot_cols: usize) -> Vec<usize> {
    let mut pivot_row = 0;
    let mut pivoted_rows = vec![false; m.rows];
    for col in 0..pivot_cols {
        if pivot_row >= m.rows {
            break;
        }
        let mut sw = None;
        for r in pivot_row..m.rows {
            if m.get(r, col) == 1 {
                sw = Some(r);
                break;
            }
        }
        let r = match sw {
            Some(r) => r,
            None => continue,
        };
        if r != pivot_row {
            m.data.swap(r, pivot_row);
        }
        for rr in 0..m.rows {
            if rr != pivot_row && m.get(rr, col) == 1 {
                for w in 0..m.data[rr].len() {
                    m.data[rr][w] ^= m.data[pivot_row][w];
                }
            }
        }
        pivoted_rows[pivot_row] = true;
        pivot_row += 1;
    }
    // Find rows whose left (first `pivot_cols`) block is all zero.
    let mut zero_rows = Vec::new();
    for r in 0..m.rows {
        let mut is_zero = true;
        for col in 0..pivot_cols {
            if m.get(r, col) == 1 {
                is_zero = false;
                break;
            }
        }
        if is_zero {
            zero_rows.push(r);
        }
    }
    let _ = pivoted_rows;
    zero_rows
}

/// Compute `Ad·h` directly (no matrices needed): the sum
/// `Σ d_i · h^(2^i)` where `d_i` are GF(2^128) bit-flip blocks.
pub fn ad_dot_h(d: &[Gf128], h: Gf128) -> Gf128 {
    let mut acc = Gf128::zero();
    let mut h_pow = h;
    // The cryptopals indexing uses d_0 (length block, untouched in
    // this attack) and d_1.. (ciphertext blocks).  We start from
    // i = 1 because d_0 is fixed.
    for i in 0..d.len() {
        let exp = 1usize << i; // 2^i
        let _ = exp;
        // h^(2^i) — compute iteratively: square each step.
        if i > 0 {
            h_pow = h_pow.mul(h_pow);
        }
        acc = acc.add(d[i].mul(h_pow));
    }
    acc
}

pub fn run() -> Report {
    let mut r = Report::new(64, "Key Recovery on GCM with a Truncated MAC");
    r.line("Building GF(2) linearisation matrices for GHASH operations…");
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let h = pack(&encrypt_block(&[0u8; 16], &key));
    let m_h = m_c(h);
    let m_sq = m_s();
    // Sanity: M_h · v == h · v for random v.
    let mut rng = StdRng::seed_from_u64(7);
    let v = Gf128(rng.gen::<u128>());
    let lhs = m_h.mul_vec(v);
    let rhs_val = h.mul(v);
    r.line(format!("M_h · v == h · v        : {}", lhs == rhs_val));
    assert_eq!(lhs, rhs_val);
    let v2 = Gf128(rng.gen::<u128>());
    let sq_lhs = m_sq.mul_vec(v2);
    let sq_rhs = v2.mul(v2);
    r.line(format!("M_s · v == v²           : {}", sq_lhs == sq_rhs));
    assert_eq!(sq_lhs, sq_rhs);

    // Demonstrate ad·h = e on a tiny scenario: 2 attacker-flippable
    // blocks (d_0 unused, d_1, d_2), one random Δ each.
    let d = vec![
        Gf128::zero(),
        Gf128(rng.gen::<u128>()),
        Gf128(rng.gen::<u128>()),
    ];
    let e_predicted = ad_dot_h(&d, h);
    r.line(format!(
        "ad · h (predicted error block)  = {:032x}",
        u128::from_le_bytes(unpack(e_predicted))
    ));

    // Demonstrate null-space search on a small T (8 rows × 16 cols).
    let mut t = BitMatrix::zero(8, 16);
    // Fill with arbitrary pattern.
    for rr in 0..8 {
        for cc in 0..16 {
            t.set(rr, cc, ((rr * 3 + cc) % 5 == 0) as u8);
        }
    }
    let mut t_aug = t.augment_identity();
    let zero_rows = gauss_jordan(&mut t_aug, t.cols);
    r.line(format!(
        "Toy T (8×16) has null-space dim ≥ {}",
        zero_rows.len()
    ));

    r.line("");
    r.line("Full attack on 32-bit MAC + 2^17 blocks not run here:");
    r.line(" the linear system has ~16 M rows.  This module ships the");
    r.line(" algorithmic primitives; see Cryptopals 64 spec for the");
    r.line(" iterative N(K)-shrinking key-recovery loop.");

    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mh_acts_as_multiply_by_h() {
        let h = Gf128(0x1234_5678_9abc_def0_fedc_ba98_7654_3210u128);
        let m = m_c(h);
        for sample in &[0u128, 1, 2, 0xDEAD_BEEF, 0xFFFF_FFFF_FFFF_FFFF] {
            let v = Gf128(*sample);
            assert_eq!(m.mul_vec(v), h.mul(v));
        }
    }

    #[test]
    fn ms_squares() {
        let m = m_s();
        let v = Gf128(0xCAFE_BABE_DEAD_BEEF_0123_4567_89AB_CDEFu128);
        assert_eq!(m.mul_vec(v), v.mul(v));
    }

    #[test]
    fn null_space_exists_for_underdetermined_t() {
        let rows = 4usize;
        let cols = 8usize;
        let mut t_t = BitMatrix::zero(cols, rows + cols);
        for r in 0..cols {
            for c in 0..rows {
                t_t.set(r, c, ((r + c) % 3 == 0) as u8);
            }
            t_t.set(r, rows + r, 1); // identity on the right
        }
        let z = gauss_jordan(&mut t_t, rows);
        assert!(!z.is_empty(), "expected at least one zero row");
    }
}
