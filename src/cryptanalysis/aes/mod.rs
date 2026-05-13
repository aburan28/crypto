//! AES cryptanalysis — reduced-round attacks for teaching.
//!
//! Full AES has never been broken in any practical sense. The biclique
//! attack on AES-256 (Bogdanov-Khovratovich-Rechberger, ASIACRYPT 2011)
//! reaches the entire cipher at `2^254.4` operations — better than brute
//! force, useless in practice. The interesting cryptanalysis happens on
//! *reduced-round* variants, where the techniques are visible without the
//! work factor exploding past anything we can run on a laptop.
//!
//! This module ships two canonical reduced-round attacks against a
//! configurable AES-128:
//!
//! - [`square`] — the **Square / integral attack**
//!   (Daemen, Knudsen, Rijmen, FSE 1997). It exploits the fact that if
//!   you fix 15 bytes of the AES-128 plaintext and let the 16th cycle
//!   through all 256 values, every byte of the state after 3 full
//!   rounds XORs to zero. Extending one round at the end gives a
//!   4-round key recovery: guess one byte of the last round key,
//!   partially decrypt, check the XOR sum. Cost: 256 chosen plaintexts
//!   and 2¹⁶ partial decryptions to recover the full last round key.
//!
//! - [`differential`] — **differential cryptanalysis** of reduced AES
//!   (Biham-Shamir style, adapted by the AES designers themselves in
//!   their proposal document). The AES S-box has maximum differential
//!   probability `4/256 = 2⁻⁶`, and the MixColumns matrix has branch
//!   number 5, so every two-round differential has at least 5 active
//!   S-boxes — already enough to give a working 2-round key recovery
//!   from a 1-round differential characteristic. We show the
//!   Differential Distribution Table, propagation through one round,
//!   and a key-recovery attack against 2 rounds of AES.
//!
//! The shared substrate is [`reduced::ReducedAes128`], a from-scratch
//! AES-128 with the round count `Nr` exposed as a parameter (1 ≤ Nr ≤ 10)
//! and a flag that controls whether the *last* round drops MixColumns,
//! as FIPS 197 specifies. The implementation tracks the production AES at
//! `Nr = 10` (verified by a test against the known NIST vector).
//!
//! ## What this is and isn't
//!
//! **Is**: a working, runnable demonstration of how 4-round AES gives
//! up its key to 256 chosen plaintexts, and how 2-round AES gives up
//! key bytes to a differential trail. Both attacks complete in
//! milliseconds. The reader can step through them, modify Nr, and see
//! how the cost grows.
//!
//! **Isn't**: an MILP/SAT-based trail search, a 7-round
//! Demirci-Selçuk meet-in-the-middle, a mixture-differential
//! distinguisher, or the biclique attack. Those rely on automated
//! tooling and large-state structures that would dwarf the rest of the
//! cryptanalysis module. They are tracked in `DEFERRED.md`.

pub mod differential;
pub mod reduced;
pub mod square;

pub use differential::{
    aes_sbox_ddt, key_recovery_two_round, max_differential_probability, propagate_one_round,
    AesDdt, DifferentialAttackReport, RoundDifference,
};
pub use reduced::{ReducedAes128, RoundOps};
pub use square::{
    is_balanced, key_recovery_four_round, lambda_set, square_distinguisher_three_round,
    SquareAttackReport,
};
