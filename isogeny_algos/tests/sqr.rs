//! Squarings (`poly::sqr_raw`: Karatsuba on squares over `Field::conv_sq`) equal the general
//! product of two equal vectors, for every limb count, at lengths around the Karatsuba and NTT
//! thresholds, with the NTT on and off, and with the largest operands (the integer p - 1).
use isogeny_algos::field::*;
use isogeny_algos::fpm::{FpM, NTT_N};
use isogeny_algos::poly;

fn check<const N: usize>(f: &FpM<N>, seed: u64) {
    let mut rng = Rng::new(seed);
    let mut pm1 = f.p;
    pm1[0] -= 1;
    for n in [
        1usize, 2, 3, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 129, 255, 300, 513,
    ] {
        for top in [false, true] {
            let a: Vec<[u64; N]> = (0..n)
                .map(|_| if top { pm1 } else { f.random(&mut rng) })
                .collect();
            let b = a.clone();
            for ntt in [usize::MAX, 1] {
                NTT_N.with(|k| k.set(Some(ntt)));
                let want = poly::mul_raw(f, &a, &b);
                let got = poly::sqr_raw(f, &a);
                let via_mul = poly::mul_raw(f, &a, &a);
                NTT_N.with(|k| k.set(None));
                assert_eq!(got, want, "N={N} n={n} top={top} ntt={ntt}");
                assert_eq!(
                    via_mul, want,
                    "N={N} n={n} top={top} ntt={ntt} (mul_raw(a, a))"
                );
            }
            assert_eq!(f.conv_sq(&a), f.conv(&a, &b), "N={N} n={n} conv_sq");
        }
    }
}

#[test]
fn squaring_one_limb() {
    check(&FpM::<1>::from_dec("2305843009213693951"), 1);
    check(&FpM::<1>::from_dec("18446744073709551557"), 2);
}

#[test]
fn squaring_two_limbs() {
    check(
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
        3,
    );
    check(
        &FpM::<2>::from_dec("340282366920938463463374607431768211297"),
        4,
    );
}

#[test]
fn squaring_four_and_eight_limbs() {
    check(
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
        5,
    );
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    check(&FpM::<8>::from_dec(csidh), 6);
}
