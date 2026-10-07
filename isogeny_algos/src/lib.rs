//! Baseline implementations of isogeny-computation algorithms.
//!
//! Kernel -> isogeny:        `kernel::{velu, kohel, sqrt_velu}`
//! Find l-isogenies from E:  `find::{divpoly, elkies (Phi_l + BMSS)}`
//! Find an isogeny E1 -> E2: `path::{galbraith, ghs, volcano (Kohel), couveignes, delfs_galbraith}`
pub mod bigint;
pub mod binary;
pub mod curve;
pub mod ext;
pub mod field;
pub mod fp2;
pub mod find;
pub mod fpm;
pub mod genus2;
pub mod gf2n;
pub mod int;
pub mod kernel;
pub mod path;
pub mod poly;
pub mod quat;
pub mod series;
pub mod testdata;
pub mod theta;
pub mod theta_g;
