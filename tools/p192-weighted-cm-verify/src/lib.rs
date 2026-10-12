//! Independent verification for the P-192 weighted-CM protocol.
//!
//! This package deliberately has no dependency on the repository's root
//! `crypto_lib` crate and does not include source files from it.  It regenerates
//! the small REF-0 admission fixture with separately implemented arithmetic.

pub mod controls;
pub mod encoding;
pub mod factor_base;
pub mod ideal;
pub mod identity;
pub mod params;
pub mod ref0;
pub mod verify;

pub type Result<T> = std::result::Result<T, String>;
