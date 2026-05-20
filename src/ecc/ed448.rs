//! **Ed448-Goldilocks signatures** — RFC 8032 §5.2.
//!
//! ## Status: not implemented
//!
//! The agent assigned to this module hit an organisational usage
//! limit before producing meaningful code.  This file currently
//! contains no functional signature operations.
//!
//! For 128-bit-security signatures, see [`super::ed25519`].  Ed448
//! is required only when 224-bit security is the bar.
//!
//! ## What an implementation would need
//!
//! - 448-bit field arithmetic over `p = 2^448 - 2^224 - 1`.
//! - Twisted Edwards curve `x² + y² = 1 + d·x²·y²` with `d = -39081 mod p`.
//! - Scalar multiplication and point encoding (57-byte compressed form).
//! - SHAKE256 (available at [`crate::hash::sha3::shake256`]).
//! - The "SigEd448" context prefix and 57-byte signature scalar `S` mod `L`.
//!
//! RFC 8032 §7.4 has the test vectors to verify against.
