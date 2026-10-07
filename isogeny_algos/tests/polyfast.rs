use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly::{self, PolyModulus};
use isogeny_algos::series;

fn rand_poly<F: Field>(f: &F, n: usize, rng: &mut Rng) -> Vec<F::E> {
    let mut v: Vec<F::E> = (0..n).map(|_| f.random(rng)).collect();
    if n > 0 {
        v[n - 1] = f.one();
    }
    v
}

fn check<F: Field>(f: &F, rng: &mut Rng) {
    for &(n, m) in &[
        (1usize, 1usize),
        (31, 33),
        (32, 32),
        (33, 100),
        (64, 64),
        (100, 37),
        (257, 255),
        (300, 7),
        (513, 300),
    ] {
        let (a, b) = (rand_poly(f, n, rng), rand_poly(f, m, rng));
        assert_eq!(
            poly::mul(f, &a, &b),
            poly::mul_schoolbook(f, &a, &b),
            "karatsuba {n}x{m}"
        );
    }
    for &n in &[2usize, 17, 40, 100, 257] {
        let m = rand_poly(f, n + 1, rng);
        let pm = PolyModulus::new(f, &m);
        for &l in &[0usize, 5, n, n + 1, 2 * n - 1, 2 * n, 2 * n + 1, 2 * n + 3] {
            let a = rand_poly(f, l, rng);
            assert_eq!(
                pm.rem(f, &a),
                poly::rem(f, &a, &poly::monic(f, &m)),
                "rem deg {n} len {l}"
            );
        }
    }
    // poly::rem dispatch boundary (regression: lengths 2n and 2n+1 recursed forever)
    for &n in &[63usize, 64, 100] {
        let b = rand_poly(f, n + 1, rng);
        for l in [2 * n - 1, 2 * n, 2 * n + 1, 2 * n + 2] {
            let a = rand_poly(f, l, rng);
            assert_eq!(poly::rem(f, &a, &b), poly::divrem(f, &a, &b).1);
        }
    }
    for &n in &[10usize, 64, 65, 200, 513] {
        let mut a: Vec<F::E> = (0..n).map(|_| f.random(rng)).collect();
        a[0] = f.one();
        assert_eq!(
            series::inv(f, &a, n),
            series::inv_quadratic(f, &a, n),
            "series inv {n}"
        );
    }
}

#[test]
fn fast_poly_and_series_match_reference() {
    let mut rng = Rng::new(101);
    check(&Zp::new(next_prime(1 << 61)), &mut rng);
    check(&FpM::<1>::from_u64_modulus(next_prime(1 << 40)), &mut rng);
    check(
        &FpM::<4>::from_dec(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        ),
        &mut rng,
    );
    check(&Zp2::new(next_prime(1 << 30)), &mut rng);
}
