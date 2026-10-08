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

/// Brent–Kung composition equals the naive sum h_i xi^i mod g, and iterating it from xi = Y^q
/// gives Y^(q^k) exactly as repeated exponentiation does (Zp at 61 bits and FpM at 127 bits,
/// moduli of degree 5..90, both reducer regimes of PolyModulus).
#[test]
fn composer_matches_naive_and_frobenius_powers() {
    fn check<F: Field>(f: &F, seed: u64) {
        let mut rng = Rng::new(seed);
        for &n in &[5usize, 12, 31, 64, 90] {
            let mut g: Vec<F::E> = (0..n).map(|_| f.random(&mut rng)).collect();
            g.push(f.one());
            let xi: Vec<F::E> = (0..n).map(|_| f.random(&mut rng)).collect();
            let c = poly::Composer::new(f, &xi, &g);
            for _ in 0..3 {
                let h: Vec<F::E> = (0..n).map(|_| f.random(&mut rng)).collect();
                let mut naive: Vec<F::E> = vec![];
                for &hi in h.iter().rev() {
                    naive = poly::add(f, &poly::mulmod(f, &naive, &xi, &g), &poly::constant(f, hi));
                }
                assert_eq!(c.compose(f, &h), poly::rem(f, &naive, &g), "n = {n}");
            }
            let y = poly::x_poly(f);
            let yq = poly::powmod_big(f, &y, &f.q(), &g);
            let frob = poly::Composer::new(f, &yq, &g);
            let (mut by_pow, mut by_comp) = (yq.clone(), yq.clone());
            for _ in 0..4 {
                by_pow = poly::powmod_big(f, &by_pow, &f.q(), &g);
                by_comp = frob.compose(f, &by_comp);
                assert_eq!(by_comp, by_pow, "n = {n}");
            }
        }
    }
    check(&Zp::new((1u64 << 61) - 1), 71);
    check(&FpM::<2>::from_dec("170141183460469231731687303715884105727"), 72);
}

/// Distinct-degree factorisation with composition Frobenius steps returns exactly what the
/// exponentiation-per-degree version returns (random squarefree polynomials of degree 8..96,
/// 61-bit and 127-bit fields).
#[test]
fn ddf_matches_exponentiation_ddf() {
    fn naive<F: Field>(f: &F, p: &[F::E]) -> Vec<(usize, Vec<F::E>)> {
        let mut out = vec![];
        let mut fp = poly::monic(f, &p.to_vec());
        let mut h = poly::x_poly(f);
        let mut k = 1usize;
        while poly::deg(f, &fp) >= 2 * k as isize {
            h = poly::powmod_big(f, &h, &f.q(), &fp);
            let g = poly::gcd(f, &fp, &poly::sub(f, &h, &poly::x_poly(f)));
            if poly::deg(f, &g) > 0 {
                fp = poly::monic(f, &poly::divrem(f, &fp, &g).0);
                h = poly::rem(f, &h, &fp);
                out.push((k, g));
            }
            k += 1;
        }
        if poly::deg(f, &fp) > 0 {
            out.push((poly::deg(f, &fp) as usize, fp));
        }
        out
    }
    fn check<F: Field>(f: &F, seed: u64) {
        let mut rng = Rng::new(seed);
        for &n in &[8usize, 24, 60, 96] {
            for _ in 0..2 {
                let p = rand_poly(f, n + 1, &mut rng);
                if poly::deg(f, &poly::gcd(f, &p, &poly::derivative(f, &p))) > 0 {
                    continue;
                }
                assert_eq!(poly::ddf(f, &p), naive(f, &p), "n = {n}");
            }
        }
    }
    check(&Zp::new((1u64 << 61) - 1), 81);
    check(&FpM::<2>::from_dec("170141183460469231731687303715884105727"), 82);
}
