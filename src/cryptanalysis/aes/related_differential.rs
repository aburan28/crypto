//! Related differentials — grouping structurally similar differential
//! trails into a single, stronger distinguishing event.
//!
//! This module takes its terminology from Yan, Tan, and Qi,
//! "Related-Differential Distinguishers on up to 7 Rounds of AES"
//! (2024), which classifies pairs of differences over the same cipher
//! state as:
//!
//! - **exchanged** — the same multiset of nonzero byte values, at the
//!   same diagonal, just permuted among its four positions (the
//!   "exchange" [`super::mixture`] already implements between two
//!   whole plaintexts);
//! - **shifted** — literally the same difference, rotated onto a
//!   different diagonal — the symmetry ShiftRows itself gives AES,
//!   already used structurally by [`super::yoyo`];
//! - **mixed** — a difference that needs *both* a rotation and a
//!   permutation to match another; and
//! - **independent** — no shared structure at all.
//!
//! # What this module is and isn't
//!
//! The paper's headline results — distinguishers reaching 5, 6, and 7
//! rounds at complexities up to `2^116` — come from an automated trail
//! search over the full related-differential space, the same kind of
//! tooling this crate's AES subtree already declines to build (see
//! `DEFERRED.md`: MILP/SAT trail search is explicitly out of scope
//! here).
//!
//! What *is* small enough to verify directly, exactly, and by running
//! real code rather than by trusting a paper's numbers, is the
//! mechanism the classification is built on: AES's round function
//! treats every diagonal identically, so a **shifted family** of
//! one-round differentials — the same post-SubBytes difference,
//! landing on each of the four rows MixColumns mixes — has an *exactly
//! computable* aggregate probability, straight from the S-box DDT in
//! [`super::differential`]. This module:
//!
//! 1. Implements the exchanged/shifted/mixed/independent taxonomy
//!    itself ([`classify`]), as a general, reusable comparison between
//!    two diagonal-confined differences.
//! 2. Derives a genuine shifted family from AES's own MixColumns
//!    matrix ([`mix_columns_shifted_family`]) and checks with
//!    [`classify`] that its four members really do classify as
//!    pairwise `Shifted` — not asserted, computed.
//! 3. Computes the exact aggregate probability of "some member of the
//!    family fires" from the real DDT ([`shifted_family_report`]).
//! 4. Confirms, by running the real reduced-round cipher, the boundary
//!    this whole family of techniques lives at: a diagonal-confined
//!    difference is *always* confined to a single column after
//!    exactly 1 round (branch number 5), and *never* confined to any
//!    single diagonal or column after 2 ([`diffusion_confinement`]) —
//!    both deterministic structural facts, not empirical tendencies.

use super::differential::AesDdt;
use super::reduced::{bytes_to_state, state_to_bytes, ReducedAes128, RoundOps};

/// Number of diagonals (and columns) in the AES state.
pub const DIAGONALS: usize = 4;

/// The four `(column, row)` positions of diagonal `d` (`d` in `0..4`):
/// exactly the positions ShiftRows collects into column `d`.
pub fn diagonal_positions(d: usize) -> [(usize, usize); 4] {
    assert!(d < DIAGONALS);
    std::array::from_fn(|c| (c, (c + DIAGONALS - d) % DIAGONALS))
}

/// A state difference confined to one diagonal: one byte per column,
/// held in column order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiagonalDifference {
    pub diagonal: usize,
    pub delta: [u8; 4],
}

impl DiagonalDifference {
    pub fn new(diagonal: usize, delta: [u8; 4]) -> Self {
        assert!(diagonal < DIAGONALS);
        Self { diagonal, delta }
    }

    pub fn single_active(diagonal: usize, column: usize, value: u8) -> Self {
        let mut delta = [0u8; 4];
        delta[column] = value;
        Self::new(diagonal, delta)
    }

    /// Embed into a full 16-byte state difference (zero elsewhere).
    pub fn to_state_diff(self) -> [u8; 16] {
        let mut out = [0u8; 16];
        for (col, &(c, r)) in diagonal_positions(self.diagonal).iter().enumerate() {
            out[4 * c + r] = self.delta[col];
        }
        out
    }

    /// The same four values, rotated onto `to_diagonal` — the relation
    /// [`classify`] calls `Shifted`.
    fn rotated(self, to_diagonal: usize) -> Self {
        let offset = (to_diagonal + DIAGONALS - self.diagonal) % DIAGONALS;
        let mut delta = [0u8; 4];
        for i in 0..DIAGONALS {
            delta[(i + offset) % DIAGONALS] = self.delta[i];
        }
        Self::new(to_diagonal, delta)
    }
}

/// How `b` relates to `a`, per the paper's taxonomy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RelatedDifferenceClass {
    /// Same diagonal, same multiset of values — `b` is `a` with its
    /// four values permuted across positions.
    Exchanged,
    /// Different diagonal, values in exactly the rotated order — `b`
    /// is `a` moved by ShiftRows' own symmetry.
    Shifted,
    /// Different diagonal, same multiset, but not a pure rotation —
    /// both a rotation and a permutation are needed.
    Mixed,
    /// Different multiset entirely: no shared structure.
    Independent,
}

/// Classify how `b` relates to `a`.
pub fn classify(a: DiagonalDifference, b: DiagonalDifference) -> RelatedDifferenceClass {
    let mut ma = a.delta;
    ma.sort_unstable();
    let mut mb = b.delta;
    mb.sort_unstable();
    if ma != mb {
        return RelatedDifferenceClass::Independent;
    }
    if a.diagonal == b.diagonal {
        return RelatedDifferenceClass::Exchanged;
    }
    if a.rotated(b.diagonal) == b {
        return RelatedDifferenceClass::Shifted;
    }
    RelatedDifferenceClass::Mixed
}

/// The exact linear half of MixColumns applied to a single active
/// byte at row `row` of a column, given its post-SubBytes difference
/// `beta`. Delegates to [`RoundOps::mix_columns`] so the coefficients
/// can never drift from the cipher's own MixColumns.
fn mix_single_row(row: usize, beta: u8) -> [u8; 4] {
    let mut block = [0u8; 16];
    block[row] = beta; // column 0, row `row`.
    let mut state = bytes_to_state(&block);
    RoundOps::mix_columns(&mut state);
    let out = state_to_bytes(&state);
    [out[0], out[1], out[2], out[3]]
}

/// The shifted family MixColumns itself produces: the same
/// post-SubBytes difference `beta`, landing on each of the four rows
/// of a column in turn. Indexed here by row number so [`classify`]
/// can be run directly on the members.
pub fn mix_columns_shifted_family(beta: u8) -> [DiagonalDifference; DIAGONALS] {
    std::array::from_fn(|row| DiagonalDifference::new(row, mix_single_row(row, beta)))
}

/// Aggregate probability of a shifted family of one-round
/// differentials, computed exactly from the real S-box DDT.
#[derive(Debug, Clone, Copy)]
pub struct ShiftedFamilyReport {
    pub alpha: u8,
    pub beta: u8,
    /// `Pr[SubBytes(alpha) = beta] = DDT[alpha][beta] / 256`, the
    /// probability of any *one* family member firing.
    pub single_trail_probability: f64,
    /// `1 - (1 - p)^4`: the probability that *some* member of the
    /// 4-way shifted family fires, treating each row's trail as an
    /// independent event (a standard assumption for this kind of
    /// aggregate bound — each row's differential is driven by an
    /// independent S-box call on an independent pair of bytes).
    pub any_member_probability: f64,
}

pub fn shifted_family_report(ddt: &AesDdt, alpha: u8, beta: u8) -> ShiftedFamilyReport {
    let p = ddt[alpha as usize][beta as usize] as f64 / 256.0;
    ShiftedFamilyReport {
        alpha,
        beta,
        single_trail_probability: p,
        any_member_probability: 1.0 - (1.0 - p).powi(DIAGONALS as i32),
    }
}

fn xor16(a: [u8; 16], b: [u8; 16]) -> [u8; 16] {
    std::array::from_fn(|i| a[i] ^ b[i])
}

/// If `state_diff` is nonzero only within a single diagonal or a
/// single column, return `(is_column, index)`. Otherwise `None`.
fn confined_line(state_diff: &[u8; 16]) -> Option<(bool, usize)> {
    for is_column in [false, true] {
        for idx in 0..DIAGONALS {
            let positions: [(usize, usize); 4] = if is_column {
                std::array::from_fn(|r| (idx, r))
            } else {
                diagonal_positions(idx)
            };
            let mut mask = [false; 16];
            for &(c, r) in &positions {
                mask[4 * c + r] = true;
            }
            let outside_zero = (0..16).all(|i| mask[i] || state_diff[i] == 0);
            let inside_nonzero = positions.iter().any(|&(c, r)| state_diff[4 * c + r] != 0);
            if outside_zero && inside_nonzero {
                return Some((is_column, idx));
            }
        }
    }
    None
}

/// Result of [`diffusion_confinement`].
#[derive(Debug, Clone, Copy)]
pub struct DiffusionReport {
    pub rounds: usize,
    pub trials: usize,
    /// Number of trials whose ciphertext difference stayed confined
    /// to a single diagonal or column.
    pub confined_count: usize,
}

/// Run `trials` random (key, plaintext) samples of `input` through
/// `rounds` rounds of AES and count how often the ciphertext
/// difference is still confined to a single diagonal or column.
///
/// At `rounds = 1` this is always `trials` (branch number 5: a single
/// active byte forces its whole column active, nothing else). At
/// `rounds = 2` it is always `0`: round 2's ShiftRows spreads that
/// column across all four columns (one active byte in each), and its
/// MixColumns then forces every one of those columns fully active —
/// full diffusion, deterministically, not just in practice.
pub fn diffusion_confinement(
    rounds: usize,
    trials: usize,
    seed: u64,
    input: DiagonalDifference,
) -> DiffusionReport {
    let mut state = seed;
    let mut byte = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u8
    };
    let diff = input.to_state_diff();
    let mut confined = 0usize;
    for _ in 0..trials {
        let mut key = [0u8; 16];
        let mut plaintext = [0u8; 16];
        for b in key.iter_mut() {
            *b = byte();
        }
        for b in plaintext.iter_mut() {
            *b = byte();
        }
        let cipher = ReducedAes128::new(&key, rounds, true);
        let p2 = xor16(plaintext, diff);
        let c1 = cipher.encrypt(&plaintext);
        let c2 = cipher.encrypt(&p2);
        if confined_line(&xor16(c1, c2)).is_some() {
            confined += 1;
        }
    }
    DiffusionReport {
        rounds,
        trials,
        confined_count: confined,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::aes::differential::aes_sbox_ddt;

    #[test]
    fn classify_same_diagonal_permutation_is_exchanged() {
        let a = DiagonalDifference::new(0, [1, 2, 3, 4]);
        let b = DiagonalDifference::new(0, [4, 3, 2, 1]);
        assert_eq!(classify(a, b), RelatedDifferenceClass::Exchanged);
        // The identity permutation counts too.
        assert_eq!(classify(a, a), RelatedDifferenceClass::Exchanged);
    }

    #[test]
    fn classify_rotation_is_shifted() {
        let a = DiagonalDifference::new(0, [10, 20, 30, 40]);
        for d in 1..DIAGONALS {
            let b = a.rotated(d);
            assert_eq!(
                classify(a, b),
                RelatedDifferenceClass::Shifted,
                "offset {d}"
            );
        }
    }

    #[test]
    fn classify_permutation_at_another_diagonal_is_mixed() {
        let a = DiagonalDifference::new(0, [1, 2, 3, 4]);
        // Same multiset, different diagonal, but not the rotation of `a`.
        let b = DiagonalDifference::new(1, [1, 2, 4, 3]);
        assert_ne!(
            a.rotated(1).delta,
            b.delta,
            "test setup: must not accidentally be a rotation"
        );
        assert_eq!(classify(a, b), RelatedDifferenceClass::Mixed);
    }

    #[test]
    fn classify_disjoint_values_is_independent() {
        let a = DiagonalDifference::new(0, [1, 2, 3, 4]);
        let b = DiagonalDifference::new(0, [5, 6, 7, 8]);
        assert_eq!(classify(a, b), RelatedDifferenceClass::Independent);
    }

    /// MixColumns' own shifted family really is pairwise `Shifted` —
    /// computed, not assumed.
    #[test]
    fn mix_columns_family_is_pairwise_shifted() {
        for beta in [1u8, 2, 0x57, 0xff] {
            let family = mix_columns_shifted_family(beta);
            for i in 0..DIAGONALS {
                for j in 0..DIAGONALS {
                    if i == j {
                        continue;
                    }
                    assert_eq!(
                        classify(family[i], family[j]),
                        RelatedDifferenceClass::Shifted,
                        "beta={beta:#x} i={i} j={j}"
                    );
                }
            }
        }
    }

    #[test]
    fn shifted_family_probability_matches_ddt() {
        let ddt = aes_sbox_ddt();
        // The AES S-box's maximum DDT entry is 4, giving p = 4/256.
        let (alpha, beta, count) = super::super::differential::max_differential_probability(&ddt);
        assert_eq!(count, 4);
        let report = shifted_family_report(&ddt, alpha, beta);
        assert!((report.single_trail_probability - 4.0 / 256.0).abs() < 1e-12);
        let p = report.single_trail_probability;
        let expected_any = 1.0 - (1.0 - p).powi(4);
        assert!((report.any_member_probability - expected_any).abs() < 1e-12);
        // The family is a real (if small) amplification over any one trail.
        assert!(report.any_member_probability > report.single_trail_probability);
    }

    #[test]
    fn one_round_is_always_column_confined() {
        let input = DiagonalDifference::single_active(0, 0, 0x7);
        let report = diffusion_confinement(1, 200, 0xc0ffee, input);
        assert_eq!(report.confined_count, report.trials);
    }

    #[test]
    fn two_rounds_never_stays_confined() {
        let input = DiagonalDifference::single_active(0, 0, 0x7);
        let report = diffusion_confinement(2, 200, 0xdecaf, input);
        assert_eq!(report.confined_count, 0);
    }
}
