use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;

fn inv_checks<const N: usize>(f: &FpM<N>, rng: &mut Rng) {
    let p = f.modulus().clone();
    // random, tiny, and near-p inputs
    let mut xs: Vec<[u64; N]> = (0..2000).map(|_| f.random(rng)).collect();
    for k in 1..50u64 {
        xs.push(f.from_u64(k));
        xs.push(f.neg(f.from_u64(k)));
    }
    xs.push(f.from_big(&p.sub_small(1).shr(1)));
    for x in xs {
        let i = f.inv(x);
        assert_eq!(f.mul(x, i), f.one(), "x = {:?}", f.to_big(&x).to_dec());
        // agrees with Fermat
        assert_eq!(i, f.pow_big(x, &p.sub_small(2)));
    }
}

#[test]
fn bingcd_inverse_matches_fermat() {
    let mut rng = Rng::new(91);
    inv_checks(
        &FpM::<4>::from_dec(
            "115792089210356248762697446949407573530086143415290314195533631308867097853951",
        ),
        &mut rng,
    );
    inv_checks(
        &FpM::<4>::from_dec(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        ),
        &mut rng,
    );
    inv_checks(&isogeny_algos::path::csidh::Csidh::<FpM<8>>::csidh512().fp, &mut rng);
    inv_checks(&FpM::<2>::from_dec("340282366920938463463374607431768211297"), &mut rng);
    // 64-bit prime with the top bit set (N = 1 but not the fast path)
    inv_checks(&FpM::<1>::from_dec("18446744073709551557"), &mut rng);
}

#[test]
fn miller_rabin_known_values() {
    use isogeny_algos::bigint::{is_probable_prime, Big};
    let primes = [
        "115792089210356248762697446949407573530086143415290314195533631308867097853951",
        "115792089237316195423570985008687907853269984665640564039457584007908834671663",
        "340282366920938463463374607431768211297",
        "18446744073709551557",
    ];
    for p in primes {
        assert!(is_probable_prime(&Big::from_dec(p)), "{p}");
    }
    let csidh = isogeny_algos::path::csidh::Csidh::<FpM<8>>::csidh512();
    assert!(is_probable_prime(csidh.fp.modulus()));
    // composites: products of two large primes, a strong pseudoprime to small bases scaled up
    let p = Big::from_dec("18446744073709551557");
    let q = Big::from_dec("340282366920938463463374607431768211297");
    assert!(!is_probable_prime(&p.mul(&q)));
    assert!(!is_probable_prime(&q.mul(&q)));
    // 3825123056546413051: strong pseudoprime to bases 2..23 (below 2^63, so the u64 path decides)
    assert!(!is_probable_prime(&Big::from_dec("3825123056546413051")));
    assert!(!is_probable_prime(&csidh.fp.modulus().mul(&Big::from_dec("18446744073709551557"))));
}
