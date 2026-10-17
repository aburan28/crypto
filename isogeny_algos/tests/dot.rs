//! The four-row dot-product kernel (`Field::dot_rev_rows`, behind the table reducer and the
//! Brent–Kung composer) equals one `dot_rev` per row: word-size and one-limb Montgomery fields,
//! lengths around the lazy-accumulation cap, row counts not a multiple of four, entries p - 1;
//! and the word-size convolution, truncated product and square kernels against schoolbook.
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;

fn check<F: Field>(f: &F, top: F::E, seed: u64) {
    let mut rng = Rng::new(seed);
    for n in [1usize, 7, 14, 15, 16, 31, 32, 33, 63, 64, 65, 100, 200] {
        for nrows in [1usize, 3, 4, 5, 9] {
            for maxed in [false, true] {
                let mut e = |_: usize| if maxed { top } else { f.random(&mut rng) };
                let a: Vec<F::E> = (0..n).map(&mut e).collect();
                let rows: Vec<Vec<F::E>> =
                    (0..nrows).map(|_| (0..n).map(&mut e).collect()).collect();
                let refs: Vec<&[F::E]> = rows.iter().map(|r| &r[..]).collect();
                let want: Vec<F::E> = rows.iter().map(|r| f.dot_rev(&a, r)).collect();
                assert_eq!(
                    f.dot_rev_rows(&a, &refs),
                    want,
                    "n {n} rows {nrows} maxed {maxed}"
                );
            }
        }
    }
}

#[test]
fn rows_equal_single_dots() {
    for p in [
        1_000_003u64,
        1_099_511_627_791,
        2_305_843_009_213_693_951,
        4_611_686_018_427_387_847,
    ] {
        check(&Zp::new(p), p - 1, p);
    }
    let f = FpM::<1>::from_dec("2305843009213693951");
    let mut pm1 = f.p;
    pm1[0] -= 1;
    check(&f, pm1, 1);
}

/// Schoolbook product with one field multiplication and addition per term.
fn naive<F: Field>(f: &F, a: &[F::E], b: &[F::E], nout: usize) -> Vec<F::E> {
    let mut r = vec![f.zero(); nout];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            if i + j < nout {
                r[i + j] = f.add(r[i + j], f.mul(x, y));
            }
        }
    }
    r
}

fn check_conv<F: Field>(f: &F, top: F::E, seed: u64) {
    let mut rng = Rng::new(seed);
    let lens = [1usize, 3, 8, 9, 15, 16, 17, 31, 32, 33, 40, 64, 70];
    for &la in &lens {
        for &lb in &lens {
            for maxed in [false, true] {
                let mut e = |_: usize| if maxed { top } else { f.random(&mut rng) };
                let a: Vec<F::E> = (0..la).map(&mut e).collect();
                let b: Vec<F::E> = (0..lb).map(&mut e).collect();
                let full = la + lb - 1;
                assert_eq!(
                    f.conv(&a, &b),
                    naive(f, &a, &b, full),
                    "{la}x{lb} maxed {maxed}"
                );
                let k = full / 2 + 1;
                assert_eq!(
                    f.conv_trunc(&a, &b, k),
                    naive(f, &a, &b, k),
                    "{la}x{lb} trunc"
                );
                if la == lb {
                    assert_eq!(f.conv_sq(&a), naive(f, &a, &a, full), "{la} square");
                }
            }
        }
    }
}

#[test]
fn word_kernels_equal_naive_products() {
    for p in [
        1_000_003u64,
        1_099_511_627_791,
        2_305_843_009_213_693_951,
        4_611_686_018_427_387_847,
    ] {
        check_conv(&Zp::new(p), p - 1, p);
    }
    let f = FpM::<1>::from_dec("2305843009213693951");
    let mut pm1 = f.p;
    pm1[0] -= 1;
    check_conv(&f, pm1, 2);
}
