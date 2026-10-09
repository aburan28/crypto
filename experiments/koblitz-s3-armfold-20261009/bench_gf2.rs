//! Fixed-operand field-kernel diagnostic; startup and input generation excluded.

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use std::hint::black_box;
use std::time::Instant;

fn next(seed: &mut u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed
}

fn main() {
    assert_eq!(std::env::consts::ARCH, "aarch64");
    for (n, low_terms) in [(41, vec![0, 3]), (53, vec![0, 1, 2, 6])] {
        let gf = Gf2::new(&IrreduciblePoly { degree: n, low_terms });
        let mask = (1u64 << n) - 1;
        let mut seed = 0xA24B_2026_1009_B3E5u64;
        let inputs: Vec<_> = (0..65_536)
            .map(|_| (next(&mut seed) & mask, next(&mut seed) & mask))
            .collect();
        for repetition in 0..9 {
            let mut checksum = 0u64;
            let started = Instant::now();
            for _ in 0..16 {
                for &(a, b) in &inputs {
                    checksum = checksum.wrapping_add(gf.mul(black_box(a), black_box(b)));
                }
            }
            let elapsed_ns = started.elapsed().as_nanos();
            println!(
                "n={n} fold={} repetition={repetition} ns_per_mul={:.6} checksum={checksum}",
                gf.folds(), elapsed_ns as f64 / (16 * inputs.len()) as f64
            );
            black_box(checksum);
        }
    }
}
