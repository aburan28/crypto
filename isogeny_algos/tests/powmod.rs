//! `poly::powmod_big` / `poly::powmod` (sliding windows for non-linear bases, table or Barrett
//! reducer) equal plain square-and-multiply with long division, for linear and non-linear bases,
//! short and long moduli (both reducers, forced either way), and dense and sparse exponents.
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
