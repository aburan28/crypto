//! Supersingular endomorphism-ring engine at toy parameters.
//!
//! Layers, bottom up:
//! - [`fp`]: `F_p` and `F_{p²}` for `p ≡ 3 (mod 4)`, `p < 2⁶³`.
//! - [`curve`]: short Weierstrass curves over `F_{p²}`, torsion bases of
//!   `E[N]` for `N | p + 1`, smooth two-dimensional discrete logarithms.
//! - [`isogeny`]: Vélu's formulas with full `(x, y)` images, prime-degree
//!   steps and smooth-degree chains, evaluation on points.
//!
//! The quaternion side (the order, its left ideals, the Deuring
//! correspondence, KLPT) is not in this crate.
//!
//! Correctness-first research code; not constant-time.  See the README.

pub mod curve;
pub mod fp;
pub mod isogeny;
