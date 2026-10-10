//! Compare the table and ARM PMULL fold kernels on the exact N41/N53 fields.

use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

fn next(seed: &mut u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed
}

fn main() {
    for n in [41, 53] {
        let curve = KoblitzCurve::new(0, n).expect("frozen Koblitz curve");
        std::env::set_var("GF2_ARM_PMULL_FOLD", "0");
        let table = Gf2::new(&curve.curve.irreducible);
        std::env::set_var("GF2_ARM_PMULL_FOLD", "1");
        let folded = Gf2::new(&curve.curve.irreducible);
        assert!(!table.folds() && folded.folds());
        let mask = (1u64 << n) - 1;
        let mut seed = 0x91e1_2026_1010_5eed ^ n as u64;
        for i in 0..4096 {
            let a = next(&mut seed) & mask;
            let b = next(&mut seed) & mask;
            assert_eq!(table.mul(a, b), folded.mul(a, b), "mul n={n} i={i}");
            assert_eq!(table.sqr(a), folded.sqr(a), "sqr n={n} i={i}");
            let k = (i % 53) as u32;
            assert_eq!(table.sqr_k(a, k), folded.sqr_k(a, k), "sqr_k n={n} i={i}");
            let nonzero = a.max(1);
            let x = table.inv(nonzero);
            let y = folded.inv(nonzero);
            assert_eq!(x, y, "inv n={n} i={i}");
            assert_eq!(table.mul(nonzero, x), 1, "inverse n={n} i={i}");
        }
        for _ in 0..32 {
            let original = (0..32)
                .map(|_| (next(&mut seed) & mask).max(1))
                .collect::<Vec<_>>();
            let mut left = original.clone();
            let mut right = original.clone();
            table.batch_inv(&mut left, &mut Vec::new());
            folded.batch_inv(&mut right, &mut Vec::new());
            assert_eq!(left, right, "batch inverse n={n}");
            for (value, inverse) in original.iter().zip(right) {
                assert_eq!(folded.mul(*value, inverse), 1);
            }
        }
        println!("n={n} pmull_fold_equivalence_passed 4096 scalar cases 32 batches");
    }
}
