use isogeny_algos::bigint::{big_mod, Big};
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::{kohel::kohel, velu::velu_cyclic, xonly::velu_xonly};
use isogeny_algos::testdata::*;

fn check_against_big<const N: usize>(f: &FpM<N>, rng: &mut Rng) {
    let p = f.modulus().clone();
    for _ in 0..200 {
        let (a, b) = (f.random(rng), f.random(rng));
        let (ab, bb) = (f.to_big(&a), f.to_big(&b));
        assert!(ab < p && bb < p);
        assert_eq!(f.to_big(&f.mul(a, b)), big_mod(&ab.mul(&bb), &p), "mul");
        assert_eq!(f.to_big(&f.add(a, b)), big_mod(&ab.add(&bb), &p), "add");
        assert_eq!(f.to_big(&f.sub(f.add(a, b), b)), ab, "sub");
        assert_eq!(f.add(a, f.neg(a)), f.zero());
        if !f.is_zero(a) {
            assert_eq!(f.mul(a, f.inv(a)), f.one(), "inv");
        }
        assert_eq!(f.to_big(&f.from_big(&ab)), ab);
    }
    // Fermat and square roots
    let a = f.random(rng);
    assert_eq!(
        f.pow_big(a, &p.sub_small(1)),
        f.one(),
        "Fermat: p must be prime"
    );
    let s = f.sq(a);
    let r = f.sqrt(s).expect("square has a root");
    assert_eq!(f.sq(r), s);
    assert_eq!(
        f.to_big(&f.from_u64(12345)),
        big_mod(&Big::from_u64(12345), &p)
    );
}

#[test]
fn fpm_matches_bigint_arithmetic() {
    let mut rng = Rng::new(91);
    check_against_big(&FpM::<1>::from_u64_modulus(next_prime(1 << 30)), &mut rng);
    check_against_big(&FpM::<1>::from_u64_modulus(next_prime(1 << 62)), &mut rng);
    check_against_big(&FpM::<1>::from_u64_modulus(18446744073709551557), &mut rng); // 2^64 - 59, generic path
    check_against_big(
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
        &mut rng,
    ); // 2^127 - 1
       // P-256 and secp256k1 primes (full 256-bit)
    check_against_big(
        &FpM::<4>::from_dec(
            "115792089210356248762697446949407573530086143415290314195533631308867097853951",
        ),
        &mut rng,
    );
    check_against_big(
        &FpM::<4>::from_dec(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        ),
        &mut rng,
    );
    // CSIDH-512 prime: 4 * (3 * 5 * ... * 373) * 587 - 1
    let mut prod = Big::from_u64(4);
    let mut l = 3u64;
    let mut count = 0;
    while count < 73 {
        if is_prime(l) {
            prod = prod.mul_small(l);
            count += 1;
        }
        l += 2;
    }
    let p512 = prod.mul_small(587).sub_small(1);
    assert_eq!(p512.bits(), 511);
    check_against_big(&FpM::<8>::new(&p512), &mut rng);
}

#[test]
fn fpm1_agrees_with_zp_on_isogenies() {
    let mut rng = Rng::new(92);
    let p = next_prime(1u64 << 40);
    let zp = Zp::new(p);
    let fm = FpM::<1>::from_u64_modulus(p);
    let lift = |x: u64| fm.from_u64(x);
    let down = |x: [u64; 1]| fm.to_canonical(&x)[0];
    for &ell in &[3u64, 5, 7, 13, 31] {
        let (e, pt) = curve_with_point(&zp, ell, &mut rng);
        let (h, _) = kernel_poly_from_point(&zp, &e, &pt, ell);
        let k0 = kohel(&zp, &e, &h, ell);
        let em = Curve::new(lift(e.a), lift(e.b));
        let hm: Vec<_> = h.iter().map(|&c| lift(c)).collect();
        let k1 = kohel(&fm, &em, &hm, ell);
        assert_eq!((down(k1.cod.a), down(k1.cod.b)), (k0.cod.a, k0.cod.b));
        let (x, y) = if let Pt::Aff(x, y) = pt {
            (x, y)
        } else {
            unreachable!()
        };
        let v1 = velu_cyclic(&fm, &em, &Pt::Aff(lift(x), lift(y)), ell);
        let x1 = velu_xonly(&fm, &em, lift(x), ell).unwrap();
        assert_eq!(v1.cod, k1.cod);
        assert_eq!(x1.cod, k1.cod);
    }
}
