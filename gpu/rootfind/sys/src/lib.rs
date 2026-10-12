//! Raw FFI for the GPU root-finding 2-decomposition oracle
//! (`gpu/rootfind/rootfind_api.h`).  Safe wrappers live in the root crate
//! at `crypto_lib::cryptanalysis::gpu_rootfind`; this crate only declares
//! the C ABI and (via `build.rs`) compiles or links the native library.
//!
//! Layouts mirror `rootfind_api.h` exactly: `rf_relation` is 16 bytes,
//! `rf_options` / `rf_result` are plain structs of `int` / `uint32_t` /
//! `double`.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int};

/// `rf_relation`: one decomposition `R_target = eps·F_i + tau·F_j`
/// (`tau == 0` marks the degenerate `R = eps·F_i`, with `j == i`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RfRelation {
    pub target: u32,
    pub i: u32,
    pub j: u32,
    pub eps: i8,
    pub tau: i8,
    pub pad: [u8; 2],
}

/// `rf_options`.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct RfOptions {
    /// `enum rf_backend`: 0 auto, 1 emulate, 2 cuda.
    pub backend: c_int,
    /// `enum rf_lookup`: 0 auto, 1 sorted, 2 hash.
    pub lookup: c_int,
    /// Factor-base slice per thread: 16, 32 or 64 (0 = 32).
    pub chunk: u32,
}

/// `rf_result`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct RfResult {
    pub n_written: u32,
    pub n_found: u32,
    pub backend_used: c_int,
    pub lookup_used: c_int,
    pub chunk_used: u32,
    pub seconds: f64,
}

pub const RF_OK: c_int = 0;
pub const RF_ERR_OVERFLOW: c_int = 1;
pub const RF_ERR_BAD_ARG: c_int = -1;
pub const RF_ERR_NO_CUDA: c_int = -2;
pub const RF_ERR_CUDA: c_int = -3;
pub const RF_ERR_BAD_CHUNK: c_int = -4;
pub const RF_ERR_NOMEM: c_int = -5;

extern "C" {
    /// 64-bit field, odd prime `3 < p < 2^62`; canonical coordinates in,
    /// up to `cap` relations out.
    pub fn rf64_decompose(
        p: u64,
        fb_x: *const u64,
        fb_y: *const u64,
        n_fb: u32,
        tx: *const u64,
        ty: *const u64,
        n_targets: u32,
        out: *mut RfRelation,
        cap: u32,
        opt: *const RfOptions,
        res: *mut RfResult,
    ) -> c_int;
    /// NIST P-256; every coordinate is 4 little-endian `u64` limbs.
    pub fn rf256_decompose(
        fb_x: *const u64,
        fb_y: *const u64,
        n_fb: u32,
        tx: *const u64,
        ty: *const u64,
        n_targets: u32,
        out: *mut RfRelation,
        cap: u32,
        opt: *const RfOptions,
        res: *mut RfResult,
    ) -> c_int;
    pub fn rf_cuda_available() -> c_int;
    pub fn rf_backend_info() -> *const c_char;
    pub fn rf_strerror(status: c_int) -> *const c_char;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layouts_match_the_c_header() {
        assert_eq!(std::mem::size_of::<RfRelation>(), 16);
        assert_eq!(std::mem::size_of::<RfOptions>(), 12);
        assert_eq!(std::mem::size_of::<RfResult>(), 32);
    }

    #[test]
    fn library_links_and_reports_a_backend() {
        // SAFETY: both functions take no arguments and return static
        // storage / an int.
        let info = unsafe { std::ffi::CStr::from_ptr(rf_backend_info()) };
        assert!(!info.to_bytes().is_empty());
        let ok = unsafe { std::ffi::CStr::from_ptr(rf_strerror(RF_OK)) };
        assert_eq!(ok.to_str().unwrap(), "ok");
    }
}
