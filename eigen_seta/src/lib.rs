//! Eigen-Séta: research implementation of an isogeny-based KEM whose
//! decryption key is a secret endomorphism of the public curve (see
//! `RESEARCH_ENDRING_PKE_STRATEGY.md` at the repository root).
//!
//! **Research code. Not reviewed, not standardised, security not established.**
//!
//! Constant-time discipline (see `SIDE_CHANNELS.md`): every routine reachable
//! from decapsulation with secret-dependent inputs is written without
//! secret-dependent branches, memory indices or early exits; field inversion
//! and square roots use fixed public exponents.  Key generation is
//! variable-time by design and is documented as such.

#![forbid(unsafe_code)]

pub mod ct;
pub mod fp;
pub mod safegcd;
pub mod fp2;
pub mod field;
pub mod scalar;
pub mod mont;
pub mod weier;
pub mod isog3;
pub mod strategy;
pub mod params;
pub mod kani;
pub mod ntheory;
pub mod chain;
pub mod cm;
pub mod pairing;
pub mod kem;
