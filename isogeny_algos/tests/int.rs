use isogeny_algos::bigint::Big;
use isogeny_algos::field::Rng;
use isogeny_algos::int::Int;

fn rand_big(rng: &mut Rng, limbs: usize) -> Big {
    let v: Vec<u64> = (0..limbs)
        .map(|_| rng.next() >> (rng.next() % 64))
        .collect();
    Big::from_limbs(&v)
}

#[test]
fn knuth_division() {
    let mut rng = Rng::new(800);
    for _ in 0..3000 {
        let (la, lb) = (1 + rng.below(8) as usize, 1 + rng.below(5) as usize);
        let a = rand_big(&mut rng, la);
        let d = rand_big(&mut rng, lb);
        if d.is_zero() {
            continue;
        }
        let (q, r) = a.divrem(&d);
        assert!(r < d);
        assert_eq!(q.mul(&d).add(&r), a);
    }
    // adversarial: divisor with top limb 2^63, dividend near multiples
    for _ in 0..500 {
        let d = Big::from_limbs(&[rng.next(), rng.next(), 1u64 << 63]);
        let q0 = rand_big(&mut rng, 3);
        let a = q0.mul(&d).add(&d.sub_small(1));
        let (q, r) = a.divrem(&d);
        assert_eq!(q, q0);
        assert_eq!(r, d.sub_small(1));
    }
    // against u128
    for _ in 0..2000 {
        let a = ((rng.next() as u128) << 64) | rng.next() as u128;
        let d = (((rng.next() >> rng.below(64)) as u128) << 64) | rng.next() as u128;
        if d == 0 {
            continue;
        }
        let (q, r) = Big::from_u128(a).divrem(&Big::from_u128(d));
        assert_eq!((q.to_u128().unwrap(), r.to_u128().unwrap()), (a / d, a % d));
    }
}

#[test]
fn signed_arithmetic() {
    let mut rng = Rng::new(801);
    for _ in 0..3000 {
        let a =
            (rng.next() as i64 >> rng.below(63)) as i128 * if rng.next() & 1 == 1 { -1 } else { 1 };
        let b =
            (rng.next() as i64 >> rng.below(63)) as i128 * if rng.next() & 1 == 1 { -1 } else { 1 };
        let (x, y) = (Int::from(a), Int::from(b));
        assert_eq!((&x + &y).to_i128(), Some(a + b));
        assert_eq!((&x - &y).to_i128(), Some(a - b));
        assert_eq!((&x * &y).to_i128(), Some(a * b));
        assert_eq!(x.cmp(&y), a.cmp(&b));
        if b != 0 {
            let q = a / b;
            let floor = if a % b != 0 && ((a < 0) != (b < 0)) {
                q - 1
            } else {
                q
            };
            assert_eq!(x.div_floor(&y).to_i128(), Some(floor));
            assert_eq!(x.modulo(&y).to_i128(), Some(a.rem_euclid(b)));
        }
        let (g, s, t) = Int::xgcd(&x, &y);
        assert_eq!(&(&x * &s) + &(&y * &t), g);
        if !g.is_zero() {
            assert!(x.modulo(&g).is_zero() && y.modulo(&g).is_zero());
        }
    }
    let big = Int::from(3i64).pow(200);
    assert_eq!((&big * &big).isqrt(), big);
    assert!((&big * &big).is_square());
    assert!(!(&(&big * &big) + &Int::one()).is_square());
    // modular square root
    let m = Int::from(1_000_000_007i64);
    for v in [2i64, 3, 5, 10, 12345] {
        if let Some(r) = Int::sqrt_mod_prime(&Int::from(v), &m) {
            assert_eq!((&r * &r).modulo(&m), Int::from(v));
        }
    }
    let p = Int::from(2i64).pow(127) - Int::one(); // Mersenne prime
    assert!(p.is_probable_prime());
    let r = Int::sqrt_mod_prime(&Int::from(9i64), &p).unwrap();
    assert_eq!((&r * &r).modulo(&p), Int::from(9i64));
}

#[test]
fn miller_rabin_matches_u128_reference() {
    use isogeny_algos::bigint::is_probable_prime;
    fn mulmod(a: u128, b: u128, m: u128) -> u128 {
        // double-and-add to avoid overflow
        let (mut r, mut a, mut b) = (0u128, a % m, b);
        while b > 0 {
            if b & 1 == 1 {
                r = (r + a) % m;
            }
            a = (a << 1) % m;
            b >>= 1;
        }
        r
    }
    fn mr(n: u128) -> bool {
        if n < 2 {
            return false;
        }
        for p in [2u128, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            if n.is_multiple_of(p) {
                return n == p;
            }
        }
        let (mut d, mut s) = (n - 1, 0);
        while d % 2 == 0 {
            d /= 2;
            s += 1;
        }
        'o: for a in [2u128, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            let mut x = 1u128;
            let (mut b, mut e) = (a, d);
            while e > 0 {
                if e & 1 == 1 {
                    x = mulmod(x, b, n);
                }
                b = mulmod(b, b, n);
                e >>= 1;
            }
            if x == 1 || x == n - 1 {
                continue;
            }
            for _ in 1..s {
                x = mulmod(x, x, n);
                if x == n - 1 {
                    continue 'o;
                }
            }
            return false;
        }
        true
    }
    let mut rng = Rng::new(803);
    let mut primes = 0;
    for _ in 0..3000 {
        let bits = 64 + rng.below(62);
        let n = (((rng.next() as u128) << 64 | rng.next() as u128) >> (128 - bits)) | 1;
        let r = mr(n);
        primes += r as usize;
        assert_eq!(is_probable_prime(&Big::from_u128(n)), r, "n = {n}");
    }
    assert!(primes > 10);
}
