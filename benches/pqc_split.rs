//! How much of an ML-KEM-768 key generation is hashing.
//!
//! This exists because the received wisdom is that Keccak dominates the lattice
//! schemes, and acting on that without checking would have sent the
//! optimisation work to the wrong place. Run it alongside the `keygen` row of
//! `cargo bench --bench pqc_speed -- ml-kem` and compare: the difference is
//! what the ring arithmetic costs. The measured split is recorded in
//! `docs/pqc-speed.md`.
use crypto_lib::hash::sha3::{shake128, shake256};
use std::hint::black_box;

#[cfg(target_arch = "x86_64")]
#[inline]
fn cycles() -> u64 { unsafe { core::arch::x86_64::_rdtsc() } }

fn measure<T>(reps: usize, mut f: impl FnMut() -> T) -> f64 {
    for _ in 0..reps / 4 { black_box(f()); }
    let mut v = Vec::with_capacity(reps);
    for _ in 0..reps {
        let a = cycles();
        black_box(f());
        v.push(cycles() - a);
    }
    v.sort_unstable();
    v[v.len() / 2] as f64
}

fn main() {
    let seed = [3u8; 34];
    // An ML-KEM-768 key generation hashes k*k = 9 matrix polynomials out of
    // SHAKE128 (about three blocks each) and 2k = 6 CBD samples out of SHAKE256
    // (128*eta bytes each), plus one SHA3-512 that rounds to nothing here.
    let mat = measure(200, || {
        let mut acc = 0u64;
        for _ in 0..9 { acc += shake128(&seed, 3 * 168)[0] as u64; }
        acc
    });
    let cbd = measure(200, || {
        let mut acc = 0u64;
        for _ in 0..6 { acc += shake256(&seed[..33], 128 * 2)[0] as u64; }
        acc
    });
    println!("ML-KEM-768 keygen hashing (9x shake128 3-block): {:.1} kc", mat / 1000.0);
    println!("ML-KEM-768 keygen hashing (6x shake256 256B)   : {:.1} kc", cbd / 1000.0);
    println!("hash subtotal                                  : {:.1} kc", (mat + cbd) / 1000.0);
    println!("compare against the ML-KEM-768 keygen row of `cargo bench --bench pqc_speed`");
}
