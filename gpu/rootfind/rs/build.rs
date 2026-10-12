//! Build script for the `gpu-rootfind` crate: compiles or links the
//! oracle's native library in `gpu/rootfind/` (the parent directory).
//!
//! Two modes, selected by the environment:
//!
//! - default: compile `../rootfind_api.cpp` with `-DRF_HAVE_CUDA=0` through
//!   the `cc` crate.  No CUDA toolchain needed; the oracle runs the
//!   kernel's per-thread function on the host (emulation).
//! - `GPU_ROOTFIND_CUDA_LIB_DIR=<dir>`: link the prebuilt
//!   `librootfind_cuda.a` from that directory (built with
//!   `make -C gpu/rootfind lib-cuda ARCH=sm_120`) plus `cudart`.  Set
//!   `CUDA_PATH` (or `CUDA_HOME`) if the runtime library is not in the
//!   default `/usr/local/cuda/lib64`.
//!
//! See RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md §5 for the exact commands.

use std::path::PathBuf;

fn main() {
    let src = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"))
        .join("..");
    println!("cargo:rerun-if-env-changed=GPU_ROOTFIND_CUDA_LIB_DIR");
    println!("cargo:rerun-if-env-changed=CUDA_PATH");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    for f in [
        "rootfind_api.cpp",
        "rootfind_api.h",
        "rootfind.cuh",
        "rootfind_host.hpp",
        "fp64.cuh",
        "fp256x64.cuh",
    ] {
        println!("cargo:rerun-if-changed={}", src.join(f).display());
    }

    if let Some(dir) = std::env::var_os("GPU_ROOTFIND_CUDA_LIB_DIR") {
        let dir = PathBuf::from(dir);
        println!("cargo:rustc-link-search=native={}", dir.display());
        println!("cargo:rustc-link-lib=static=rootfind_cuda");
        let cuda = std::env::var_os("CUDA_PATH")
            .or_else(|| std::env::var_os("CUDA_HOME"))
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("/usr/local/cuda"));
        println!(
            "cargo:rustc-link-search=native={}",
            cuda.join("lib64").display()
        );
        println!(
            "cargo:rustc-link-search=native={}",
            cuda.join("lib").display()
        );
        println!("cargo:rustc-link-lib=dylib=cudart");
        println!("cargo:rustc-link-lib=dylib=stdc++");
        return;
    }

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .file(src.join("rootfind_api.cpp"))
        .include(&src)
        .define("RF_HAVE_CUDA", "0")
        .flag_if_supported("-std=c++17")
        .flag_if_supported("-pthread")
        .flag_if_supported("-Wno-unknown-pragmas")
        .flag_if_supported("-Wno-unused-parameter")
        .opt_level(3);
    build.compile("rootfind_oracle");
}
