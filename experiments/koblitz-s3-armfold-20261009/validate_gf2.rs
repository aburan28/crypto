//! Standalone differential check of the installed ARM PMULL folding path.
//! Compile against the release `crypto_lib` rlib from this exact source tree.

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

fn next(seed: &mut u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed
}

fn reference_mul(a: u64, b: u64, n: u32, modulus: u64) -> u64 {
    let mut wide = 0u128;
    for bit in 0..n {
        if (b >> bit) & 1 != 0 {
            wide ^= (a as u128) << bit;
        }
    }
    for bit in (n..2 * n).rev() {
        if (wide >> bit) & 1 != 0 {
            wide ^= (modulus as u128) << (bit - n);
        }
    }
    (wide as u64) & ((1u64 << n) - 1)
}

fn main() {
    assert_eq!(std::env::consts::ARCH, "aarch64");
    let fold_requested = std::env::var("GF2_ARM_PMULL_FOLD").ok().as_deref() == Some("1");
    let mut seed = 0xA24B_2026_1009_F01Du64;
    let mut checks = 0usize;
    for n in 2u32..=63 {
        let tail_degrees = [1u32, n / 4, (n + 1) / 2, n - 1];
        for d in tail_degrees {
            if d == 0 || d >= n {
                continue;
            }
            let irr = IrreduciblePoly {
                degree: n,
                low_terms: vec![0, d],
            };
            let gf = Gf2::new(&irr);
            assert_eq!(gf.kernel_name(), "pmull");
            assert_eq!(gf.folds(), fold_requested && 2 * d <= n + 1,
                       "fold dispatch n={n} d={d}");
            let modulus = (1u64 << n) | (1u64 << d) | 1;
            let mask = (1u64 << n) - 1;
            for _ in 0..1000 {
                let a = next(&mut seed) & mask;
                let b = next(&mut seed) & mask;
                assert_eq!(gf.mul(a, b), reference_mul(a, b, n, modulus),
                           "n={n} d={d} a={a:#x} b={b:#x}");
                assert_eq!(gf.sqr(a), reference_mul(a, a, n, modulus),
                           "square n={n} d={d} a={a:#x}");
                checks += 2;
            }
        }
    }
    println!("ARM PMULL field products checked: {checks}; fold requested: {fold_requested}");
}
