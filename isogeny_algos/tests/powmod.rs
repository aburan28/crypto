//! `poly::powmod_big` / `poly::powmod` (sliding windows for non-linear bases, table or Barrett
//! reducer) equal plain square-and-multiply with long division, for linear and non-linear bases,
//! short and long moduli (both reducers, forced either way; Barrett with its prepared and
//! cyclic products), dense and sparse exponents, with the IFMA transforms on and off.
use isogeny_algos::bigint::Big;
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly::{self, Poly, POWMOD_TABLE_N};

fn naive<F: Field>(f: &F, a: &Poly<F>, e: &Big, m: &Poly<F>) -> Poly<F> {
    let mut r = poly::constant(f, f.one());
    let b = poly::divrem(f, a, m).1;
    for i in (0..e.bits()).rev() {
        r = poly::divrem(f, &poly::mul(f, &r, &r), m).1;
        if e.bit(i) {
            r = poly::divrem(f, &poly::mul(f, &r, &b), m).1;
        }
    }
    r
}

fn check<const N: usize>(f: &FpM<N>, seed: u64) {
    let mut rng = Rng::new(seed);
    let q = f.q();
    let exps = [
        q.clone(),
        q.sub_small(1).shr(1), // (q - 1) / 2
        Big::from_u64(1),
        Big::from_u64(2),
        Big::from_u64(0b1000_0000_0001), // sparse
        Big::from_u64(u64::MAX),         // all ones
    ];
    for deg in [3usize, 20, 63, 64, 100, 170] {
        let mut m: Poly<FpM<N>> = (0..deg).map(|_| f.random(&mut rng)).collect();
        m.push(f.random(&mut rng)); // not monic
        for base_len in [2usize, 3, deg + 5] {
            let a: Poly<FpM<N>> = (0..base_len).map(|_| f.random(&mut rng)).collect();
            for e in &exps {
                let want = naive(f, &a, e, &m);
                for force in [None, Some(usize::MAX), Some(0)] {
                    POWMOD_TABLE_N.with(|c| c.set(force));
                    let got = poly::powmod_big(f, &a, e, &m);
                    POWMOD_TABLE_N.with(|c| c.set(None));
                    assert_eq!(
                        got, want,
                        "N={N} deg {deg} base_len {base_len} force {force:?}"
                    );
                }
                if e.bits() <= 64 {
                    let e64 = e.limbs(1)[0] as u128;
                    assert_eq!(poly::powmod(f, &a, e64, &m), want, "N={N} powmod deg {deg}");
                }
            }
        }
    }
}

#[test]
fn powmod_matches_naive() {
    for simd in [true, false] {
        isogeny_algos::ntt::SIMD.with(|s| s.set(simd));
        check_all();
    }
    isogeny_algos::ntt::SIMD.with(|s| s.set(true));
}

#[test]
fn barrett_matches_division_eight_limbs() {
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    let f = FpM::<8>::from_dec(csidh);
    let mut rng = Rng::new(8);
    for simd in [true, false] {
        isogeny_algos::ntt::SIMD.with(|s| s.set(simd));
        for deg in [30usize, 64, 90] {
            let mut m: Poly<FpM<8>> = (0..deg).map(|_| f.random(&mut rng)).collect();
            m.push(f.one());
            let pm = poly::PolyModulus::new(&f, &m);
            for _ in 0..3 {
                let a: Poly<FpM<8>> = (0..2 * deg - 1).map(|_| f.random(&mut rng)).collect();
                assert_eq!(
                    pm.rem(&f, &a),
                    poly::divrem(&f, &a, &m).1,
                    "deg {deg} simd {simd}"
                );
            }
        }
    }
    isogeny_algos::ntt::SIMD.with(|s| s.set(true));
}

fn check_all() {
    check(&FpM::<1>::from_dec("2305843009213693951"), 1);
    check(
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
        2,
    );
    check(
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
        3,
    );
}
