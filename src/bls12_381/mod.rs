//! # BLS12-381 — pairing-friendly elliptic curve.
//!
//! Barreto-Lynn-Scott family, `k = 12` embedding degree.  The
//! standard SNARK / BLS-signature curve adopted by Zcash, Ethereum
//! 2.0, IETF draft-irtf-cfrg-bls-signature, and most modern
//! pairing-based cryptography deployments.
//!
//! ## Parameters
//!
//! - `u = -0xd201000000010000` (BLS parameter; defines `r`, `p`)
//! - `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
//!   (381-bit prime, base field characteristic)
//! - `r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`
//!   (255-bit prime, group order on G1 and G2)
//! - `G1` over `F_p`, `G2` over `F_{p²}`, target group `F_{p¹²}`.
//!
//! ## What this module provides
//!
//! - [`fq`] — base field `F_p` (381-bit arithmetic via BigUint).
//! - [`fq2`] — quadratic extension `F_p² = F_p[u]/(u² + 1)`.
//! - [`fq6`] — cubic extension `F_p⁶ = F_p²[v]/(v³ − (u + 1))`.
//! - [`fq12`] — quadratic extension `F_p¹² = F_p⁶[w]/(w² − v)`.
//! - [`g1`] — curve `E(F_p): y² = x³ + 4`.
//! - [`g2`] — twist `E'(F_p²): y² = x³ + 4(u + 1)`.
//! - [`pairing`] — Miller loop + final exponentiation, computing
//!   `e: G1 × G2 → F_{p¹²}`.
//!
//! ## Limitations
//!
//! - **Educational, not production.**  Uses `BigUint` throughout —
//!   slow (multi-second pairings).  A production library would
//!   use Montgomery-form 64-bit-limb arithmetic with constant-time
//!   ops.  See the `blst` or `bls12_381` crates for real
//!   performance.
//! - **Variable-time.**  Suitable for prover-side code only.
//! - **No subgroup-check optimisations.**  We perform full scalar
//!   multiplication for subgroup membership tests; BLS exposes
//!   faster methods via the GLS endomorphism that we don't yet
//!   implement.
//!
//! ## Reference
//!
//! - Barreto, Lynn, Scott 2003: "Constructing elliptic curves with
//!   prescribed embedding degrees."
//! - Bowe 2017: "BLS12-381: New zk-SNARK Elliptic Curve Construction."
//! - Sean Bowe / Zcash protocol specification §5.4.9.

pub mod fq;
pub mod fq2;
pub mod fq6;
pub mod fq12;
pub mod g1;
pub mod g2;
pub mod pairing;

pub use fq::Fq;
pub use fq2::Fq2;
pub use fq6::Fq6;
pub use fq12::Fq12;
pub use g1::G1Point;
pub use g2::G2Point;
pub use pairing::pairing;
