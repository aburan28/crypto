//! # Zero-knowledge proofs and commitments
//!
//! This module collects the **foundational primitives** for building
//! zero-knowledge protocols over the elliptic-curve infrastructure
//! in [`crate::ecc`].
//!
//! ## What's included (Phase 1)
//!
//! - [`schnorr_zkp`] — Schnorr's **sigma protocol** for proof of
//!   knowledge of a discrete logarithm: given a public key
//!   `Y = x·G`, the prover demonstrates knowledge of `x` without
//!   revealing it.  Both interactive and **Fiat-Shamir** (non-
//!   interactive) variants.  This is the textbook starting point
//!   for almost every ZK system.
//!
//! - [`pedersen`] — Pedersen **commitments**: `Com(m, r) = m·G + r·H`
//!   where `G, H` are independent generators (with unknown discrete-
//!   log relation between them).  **Perfectly hiding** (uniform
//!   `r` makes `Com` uniform), **computationally binding** under
//!   the discrete-log assumption.
//!
//! - [`chaum_pedersen`] — Chaum-Pedersen **equality of discrete
//!   logs**: given public `Y = x·G` and `Z = x·H`, prove that both
//!   share the same secret `x` *without* revealing it.  Composes
//!   into more complex sigma protocols (e.g. ElGamal re-encryption
//!   proofs, ring signatures).
//!
//! - [`merkle`] — Binary **Merkle tree** over SHA-256.  Used as a
//!   commitment to a list of leaves; supports `O(log n)` inclusion
//!   proofs.  Foundation for STARKs, blockchain transaction
//!   commitments, certificate transparency, etc.
//!
//! ## Conventions
//!
//! - **Curve choice**: all sigma protocols default to **secp256k1**
//!   for compatibility with the existing Schnorr-signature code.
//!   Generalising to P-256 is trivial — pass a different
//!   [`crate::ecc::curve::CurveParams`].
//!
//! - **Field encoding**: scalars are `BigUint` reduced mod `n`
//!   (the group order).  Points are [`crate::ecc::point::Point`].
//!
//! - **Hash domain separation**: Fiat-Shamir challenges use a
//!   protocol-specific tag prefix to prevent transcript-replay
//!   attacks across different proof systems.
//!
//! ## What this is NOT
//!
//! - **Not constant-time.**  These primitives use `BigUint`
//!   arithmetic via the existing `Point::scalar_mul` (not the CT
//!   ladder).  Suitable for prover-side computation; **not** for
//!   verifier-side validation of attacker-controlled inputs in a
//!   timing-sensitive context.
//!
//! - **Not SNARKs.**  Sigma protocols are *interactive* (or
//!   non-interactive via Fiat-Shamir) but **not succinct**.  A
//!   Schnorr proof for a 256-bit DLP is ~64 bytes; a SNARK proof
//!   for the same statement could be ~200 bytes but with a much
//!   more complex prover.  See `RESEARCH_P256.md`'s ZK section for
//!   the roadmap toward Bulletproofs, KZG, Groth16, PLONK, STARKs.
//!
//! - **Not pairing-based.**  No pairings (BLS12-381, BN254) yet —
//!   those are required for Groth16/PLONK SNARKs, and are a
//!   substantial separate undertaking.

pub mod chaum_pedersen;
pub mod merkle;
pub mod pedersen;
pub mod schnorr_zkp;

pub use chaum_pedersen::{
    chaum_pedersen_prove, chaum_pedersen_prove_with_nonce, chaum_pedersen_verify,
    ChaumPedersenProof,
};
pub use merkle::{
    merkle_root, merkle_proof, merkle_verify, MerkleProof, MerkleTree,
};
pub use pedersen::{
    pedersen_commit, pedersen_commit_with_blinding, pedersen_open,
    pedersen_second_generator, PedersenCommitment, PedersenParams,
};
pub use schnorr_zkp::{
    schnorr_zkp_prove, schnorr_zkp_prove_with_nonce, schnorr_zkp_verify,
    SchnorrZkProof,
};
