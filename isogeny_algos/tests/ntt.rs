//! The multi-prime NTT product (`FpM::conv_long`) equals the Karatsuba/schoolbook product
//! exactly: every limb count, balanced and unbalanced lengths, squarings, lengths around the
//! thresholds and just past powers of two, and operands at the top of the range (p - 1).
use isogeny_algos::field::*;
use isogeny_algos::fpm::{FpM, NTT_N};
use isogeny_algos::poly;

fn check<const N: usize>(f: &FpM<N>, seed: u64) {
    // both transform families: AVX-512 IFMA where the CPU has it, and the portable one
    for simd in [true, false] {
        isogeny_algos::ntt::SIMD.with(|s| s.set(simd));
        check_mode(f, seed);
    }
    isogeny_algos::ntt::SIMD.with(|s| s.set(true));
}

fn check_mode<const N: usize>(f: &FpM<N>, seed: u64) {
    let mut rng = Rng::new(seed);
    let thr = f.ntt_threshold();
    let lens: Vec<(usize, usize)> = vec![
        (1, 1),
        (2, 3),
        (17, 9),
        (thr - 1, thr - 1),
        (thr, thr),
        (thr + 1, 2 * thr + 5),
        (3 * thr, thr),
        (700, 333),
        (1500, 1500),
        // just past a power of two: cyclic transform plus the wrapped top (`FpM::conv_long`)
        (129, 129),
        (160, 160),
        (130, 128),
        (257, 200),
    ];
    // the integer p - 1 as a residue (largest operand: the CRT bound is tight there)
    let mut pm1 = f.p;
    pm1[0] -= 1;
    for (la, lb) in lens {
        for kind in 0..3 {
            let a: Vec<[u64; N]> = (0..la)
                .map(|i| match kind {
                    0 => f.random(&mut rng),
                    1 => pm1,
                    _ => {
                        if i % 3 == 0 {
                            f.zero()
                        } else {
                            f.random(&mut rng)
                        }
                    }
                })
                .collect();
            let b: Vec<[u64; N]> = (0..lb)
                .map(|_| if kind == 1 { pm1 } else { f.random(&mut rng) })
                .collect();
            NTT_N.with(|k| k.set(Some(usize::MAX)));
            let want = poly::mul_raw(f, &a, &b);
            let a2 = a.clone(); // a distinct slice: the general product, not the squaring path
            let want_sq = poly::mul_raw(f, &a, &a2);
            NTT_N.with(|k| k.set(Some(1)));
            let got = poly::mul_raw(f, &a, &b);
            let got_sq = poly::mul_raw(f, &a, &a);
            NTT_N.with(|k| k.set(None));
            assert_eq!(got, want, "N={N} lens ({la},{lb}) kind {kind}");
            assert_eq!(got_sq, want_sq, "N={N} square len {la} kind {kind}");
        }
    }
}

#[test]
fn ntt_equals_karatsuba_one_limb() {
    check(&FpM::<1>::from_dec("2305843009213693951"), 1);
    check(&FpM::<1>::from_dec("18446744073709551557"), 2); // largest 64-bit prime
}

#[test]
fn ntt_equals_karatsuba_two_limbs() {
    check(
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
        3,
    );
    // a 128-bit prime with the top bit set (2^128 - 159)
    check(
        &FpM::<2>::from_dec("340282366920938463463374607431768211297"),
        4,
    );
}

#[test]
fn ntt_equals_karatsuba_four_limbs() {
    check(
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
        5,
    );
}

#[test]
fn ntt_equals_karatsuba_eight_limbs() {
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    check(&FpM::<8>::from_dec(csidh), 6);
}
