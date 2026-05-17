//! # Cryptanalysis — toolkit for analysing novel ciphers and hashes.
//!
//! This module gives security researchers a set of building blocks for
//! evaluating proposed symmetric primitives **before** they are deployed.
//! The intent is to make the kinds of analysis that used to require
//! ad-hoc Python scripts (DDT / LAT printouts, Walsh transforms, SAC
//! matrices, chi-squared distinguishers) into a single library that
//! plugs cleanly into a cipher's existing Rust implementation.
//!
//! ## What's included
//!
//! - [`sbox::Sbox`] — generic n-in / m-out S-box with the full modern
//!   distinguishing-table battery:
//!   [`Sbox::ddt`] (Differential Distribution Table),
//!   [`Sbox::lat`] (Linear Approximation Table),
//!   [`Sbox::bct`] (Boomerang Connectivity Table — Cid et al.,
//!   EUROCRYPT 2018),
//!   [`Sbox::dlct`] (Differential-Linear Connectivity Table —
//!   Bar-On et al., EUROCRYPT 2019),
//!   [`Sbox::truncated_ddt`] (Knudsen, FSE 1994),
//!   plus the derived metrics
//!   [`Sbox::differential_uniformity`],
//!   [`Sbox::max_differential_probability`],
//!   [`Sbox::max_linear_bias`],
//!   [`Sbox::nonlinearity`],
//!   [`Sbox::boomerang_uniformity`],
//!   [`Sbox::max_dlct_bias`],
//!   [`Sbox::is_balanced`],
//!   [`Sbox::is_bijective`],
//!   and [`Sbox::algebraic_degree`].
//!
//! - [`boolean`] — Boolean-function helpers: Walsh–Hadamard transform,
//!   algebraic normal form, algebraic degree.  These are the
//!   primitives the S-box analysis is layered on top of, exposed
//!   independently because designers of stream ciphers and
//!   non-table-based round functions often want them.
//!
//! - [`avalanche`] — diffusion measurements over arbitrary functions
//!   (`fn(&[u8]) -> Vec<u8>`).  Full avalanche matrix, Strict
//!   Avalanche Criterion (SAC) score, Bit-Independence Criterion
//!   (BIC) score.  Works on any function you can call from Rust —
//!   useful for end-to-end testing of your full cipher, not just
//!   individual components.
//!
//! - [`statistical`] — chi-squared, monobit, runs, and byte-frequency
//!   distinguishers over the output of an arbitrary cipher / hash /
//!   PRF.  These are the classic "is the output statistically
//!   distinguishable from random?" tests.
//!
//! ## How to use it
//!
//! 1. Wrap your S-box(es) in [`Sbox::new`] and call [`Sbox::report`]
//!    to get a one-page summary of differential / linear strength.
//! 2. For your full round function or full cipher, expose it as a
//!    closure `|input: &[u8]| -> Vec<u8>` and pass it to
//!    [`avalanche::full_avalanche`] or
//!    [`statistical::chi_squared_byte_test`].
//! 3. For a Boolean function (say, a tap from your nonlinear filter),
//!    convert it to a truth table (`Vec<u8>` of 0s/1s) and call
//!    [`boolean::walsh_hadamard`] or [`boolean::algebraic_degree`].
//!
//! ## What this is NOT
//!
//! - **Not a full automated trail-search engine.**  Searching for the
//!   best `r`-round differential trail in a cipher with state larger
//!   than ~16 bits is genuinely hard (MILP / SAT).  We give you the
//!   per-round building blocks (DDT, LAT) so you can plug them into
//!   your own search; we do not ship a CP/SAT/MILP backend.
//! - **Not a substitute for academic peer review.**  Passing every
//!   test in this module does NOT mean a cipher is secure.  Lots of
//!   broken ciphers had clean DDT/LAT and high SAC scores.  These
//!   tools are necessary, not sufficient.
//! - **Not constant-time.**  This module is *analytical*, not
//!   operational — it is meant to be run on a designer's workstation
//!   against a candidate cipher, not on a production server.  The
//!   code uses `Vec` allocations and `match`-based dispatch
//!   throughout.
//!
//! ## Example: analyse a 4-bit S-box
//!
//! ```
//! use crypto::cryptanalysis::sbox::Sbox;
//!
//! // Serpent S0 — Anderson/Biham/Knudsen 1998.
//! let s0 = Sbox::new(4, 4, vec![3, 8, 15, 1, 10, 6, 5, 11, 14, 13, 4, 2, 7, 0, 9, 12]).unwrap();
//!
//! assert!(s0.is_bijective());
//! assert!(s0.is_balanced());
//! // Serpent's S-boxes have differential uniformity 4 (i.e. max DDT entry = 4).
//! assert_eq!(s0.differential_uniformity(), 4);
//! // ... and max DP = 4/16 = 0.25.
//! assert!((s0.max_differential_probability() - 0.25).abs() < 1e-9);
//! ```

pub mod aes;
pub mod ai_schoof;
pub mod aut_folded_rho;
pub mod auto_attack;
pub mod avalanche;
pub mod b_seed_profile;
pub mod binary_isogeny;
pub mod binary_semaev;
pub mod bleichenbacher;
pub mod boolean;
pub mod boomerang;
pub mod canonical_lift;
pub mod cga_hnc;
pub mod cipher_registry;
pub mod cm_canonical_lift;
pub mod coleman_integration;
pub mod ec_index_calculus;
pub mod ec_index_calculus_j0;
pub mod ec_trapdoor;
pub mod ecdsa_audit;
pub mod ecm;
pub mod fght_snfs;
pub mod ghs_descent;
pub mod ghs_full_attack;
pub mod hash_attacks;
pub mod hilbert_class_poly;
pub mod hnp_ecdsa;
pub mod j0_twists;
pub mod lattice;
pub mod legacy_curve_attacks;
pub mod mazur_tate_sigma;
pub mod md5_differential;
pub mod ml_rho_walks;
pub mod modular_polynomial;
pub mod multi_key_hnp;
pub mod cheon_attack;
pub mod invalid_curve_attack;
pub mod mov_attack;
pub mod nonanom_formal_log;
pub mod orbit_homology;
pub mod p256_attacks;
pub mod p256_isogeny_cover;
pub mod p256_structural;
pub mod pohlig_hellman;
pub mod pollard_rho;
pub mod preprocessing_rho;
pub mod quantum_estimator;
pub mod research_bench;
pub mod sbox;
pub mod sha1_differential;
pub mod shor;
pub mod signature_corpus;
pub mod solinas_correlations;
pub mod statistical;
pub mod signal_ratchet;
pub mod tls12_kdf;
pub mod tls13_kdf;
pub mod visual_demos;
pub mod visualize;

pub use aut_folded_rho::{
    apply_aut, aut_folded_rho_dlp, canonical_form, AutElt, FoldedRhoOptions, FoldedRhoSolution,
    J0CurveAut,
};
pub use avalanche::{bit_independence_score, full_avalanche, sac_score, AvalancheReport};
pub use bleichenbacher::{
    bias_magnitude, bleichenbacher_direct, signature_to_sample, BleichenbacherPeak,
    BleichenbacherSample,
};
pub use boolean::{algebraic_degree, anf_coefficients, walsh_hadamard};
pub use boomerang::{
    boomerang_distinguisher, boomerang_trail_search, differential_trail_search, rectangle_attack,
    sandwich_distinguisher, BlockCipher, BoomerangResult, BoomerangTrailPair, DifferentialTrail,
    RectangleResult, SandwichResult, SpnTrailModel, ToySpn,
};
pub use canonical_lift::{
    find_anomalous_curve, hensel_lift_point, smart_attack_anomalous, ZpCurve, ZpInt, ZpPoint,
};
pub use ec_index_calculus::{
    build_factor_base, ec_index_calculus_dlp, find_one_relation, find_roots_fp,
    gaussian_eliminate_mod_n, pollard_rho_ecdlp, semaev_s3, semaev_s3_in_x3, semaev_s4_in_x4,
    sqrt_mod_p, FactorBaseEntry, Relation,
};
pub use ec_index_calculus_j0::{
    build_eisenstein_factor_base, eisenstein_smooth_ic_dlp, j0_index_calculus_dlp,
};
pub use ecdsa_audit::{
    audit_ecdsa_transcript, quick_bias_score, AuditOptions, AuditResult, EcdsaSample,
};
pub use hnp_ecdsa::{
    hnp_recover_key, hnp_recover_key_with_reduction, BiasedSignature, HnpReduction,
};
pub use j0_twists::{
    enumerate_twists, factorise_small, format_twist_table, max_prime_factor, naive_point_count,
    primitive_root, twist_coefficients, TwistInfo,
};
pub use lattice::{bkz_reduce, lll_reduce};
pub use legacy_curve_attacks::{
    bounded_bsgs_binary, bounded_bsgs_prime, legacy_curve_attack_report,
    run_legacy_curve_attack_demos, BoundedDlpSolution, LegacyCurveAttackDemo,
};
pub use multi_key_hnp::{build_transcript, multi_key_hnp_recover_master, ChildKeySignature};
pub use pollard_rho::{
    pollard_rho_dlp, pollard_rho_dlp_zp, pollard_rho_dlp_zp_multi, pollard_rho_dp_dlp_zp,
    pollard_rho_dp_dlp_zp_multi, DpRhoOptions, RhoOptions, RhoSolution,
};
pub use preprocessing_rho::{
    build_preprocessing_table, expected_online_cost, online_solve, preprocessing_rho_dlp,
    PreprocessingOptions, PreprocessingTable,
};
pub use sbox::{Sbox, SboxReport};
pub use sha1_differential::{
    estimate_differential, find_near_collision, round_function_truth_table, sha1, sha1_avalanche,
    sha1_compress, DifferentialEstimate, NearCollision,
};
pub use shor::{shor_factor, shor_order_find};
pub use signature_corpus::{
    CorpusAnalyzer, CorpusReport, Finding, ReportRow, Severity, SignatureRecord,
};
pub use statistical::{chi_squared_byte_test, monobit_test, runs_test, ChiSquaredReport};
