//! Field arithmetic against num-bigint, for several limb counts.
use eigen_seta::fp::{Fp, FpConfig};
use eigen_seta::fp2::Fp2;
use num_bigint::BigUint;
use num_traits::{One, Zero};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

fn to_big<const N: usize>(l: &[u64; N]) -> BigUint {
    let mut b = BigUint::zero();
    for i in (0..N).rev() {
        b = (b << 64u32) + BigUint::from(l[i]);
    }
    b
}

fn from_big<const N: usize>(b: &BigUint) -> [u64; N] {
    let mut l = [0u64; N];
    let d = b.to_u64_digits();
    for (i, x) in d.iter().enumerate() {
        l[i] = *x;
    }
    l
}

macro_rules! field_test {
    ($name:ident, $cfg:ident, $n:expr, $hex:expr) => {
        #[derive(Clone, Copy, Debug)]
        struct $cfg;
        impl FpConfig<$n> for $cfg {
            const MODULUS: [u64; $n] = eigen_seta::hexlimbs!($hex, $n);
            const NAME: &'static str = stringify!($cfg);
        }
        #[test]
        fn $name() {
            type F = Fp<$cfg, $n>;
            type F2 = Fp2<$cfg, $n>;
            let p = to_big(&<$cfg as FpConfig<$n>>::MODULUS);
            assert_eq!(&p % 4u32, BigUint::from(3u32));
            let mut rng = Rng(42 + $n);
            let rand_big = |rng: &mut Rng| {
                let mut l = [0u64; $n];
                for i in 0..$n {
                    l[i] = rng.next();
                }
                to_big(&l) % &p
            };
            for _ in 0..200 {
                let a = rand_big(&mut rng);
                let b = rand_big(&mut rng);
                let fa = F::from_canonical_limbs(from_big(&a));
                let fb = F::from_canonical_limbs(from_big(&b));
                assert_eq!(to_big(&fa.to_canonical_limbs()), a);
                assert_eq!(to_big(&(fa + fb).to_canonical_limbs()), (&a + &b) % &p);
                assert_eq!(to_big(&(fa - fb).to_canonical_limbs()), (&a + &p - &b) % &p);
                assert_eq!(to_big(&(fa * fb).to_canonical_limbs()), (&a * &b) % &p);
                assert_eq!(to_big(&(-fa).to_canonical_limbs()), (&p - &a) % &p);
                if !a.is_zero() {
                    assert_eq!(fa * fa.invert(), F::one());
                    assert_eq!(fa.invert(), fa.invert_fermat(), "safegcd vs Fermat");
                }
                let sq = fa.square();
                let (r, ok) = sq.sqrt();
                assert!(bool::from(ok));
                assert_eq!(r.square(), sq);
                assert!(bool::from(sq.is_square()));
                // bytes round trip
                let (back, ok) = F::from_bytes_be(&fa.to_bytes_be());
                assert!(bool::from(ok) && back == fa);
                // Fp2
                let c = rand_big(&mut rng);
                let d = rand_big(&mut rng);
                let x = F2::new(fa, fb);
                let y = F2::new(F::from_canonical_limbs(from_big(&c)), F::from_canonical_limbs(from_big(&d)));
                let xy = x * y;
                let re = (&a * &c + &p * &p - (&b * &d) % &p) % &p;
                let im = (&a * &d + &b * &c) % &p;
                assert_eq!(to_big(&xy.c0.to_canonical_limbs()), re);
                assert_eq!(to_big(&xy.c1.to_canonical_limbs()), im);
                assert_eq!(x.square(), x * x);
                if !(a.is_zero() && b.is_zero()) {
                    assert_eq!(x * x.invert(), F2::one());
                }
                let x2 = x.square();
                let (r2, ok2) = x2.sqrt();
                assert!(bool::from(ok2), "sqrt of a square failed");
                assert_eq!(r2.square(), x2);
                // a random element is a square half the time; check consistency
                let (r3, ok3) = x.sqrt();
                assert_eq!(bool::from(ok3), bool::from(x.is_square()));
                if bool::from(ok3) {
                    assert_eq!(r3.square(), x);
                }
            }
            assert_eq!(F::zero().invert(), F::zero());
            assert_eq!(F2::zero().invert(), F2::zero());
            let (z, okz) = F2::zero().sqrt();
            assert!(bool::from(okz) && z == F2::zero());
            let _ = BigUint::one();
        }
    };
}

// Moduli = 3 mod 4: 2^61 - 1 (1 limb), 2^127 - 1 (2 limbs), 2^521 - 1 (9
// limbs), and the L1 and L5 primes from src/params_generated.rs (21 and 42
// limbs; these exercise the Montgomery-friendly reduction and safegcd).
field_test!(f1, P61, 1, "1fffffffffffffff");
field_test!(f2, P127, 2, "7fffffffffffffffffffffffffffffff");
field_test!(f3, P521ish, 9, "1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
field_test!(fl1, PL1, 21, "1e90dd6696c7cb6ab11dd31b2fc8b8c51a85a0269c88a697887798072b7f6ec1a1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
field_test!(fl5, PL5, 42, "223e047a381431b6bd283938b1136820f011153efd3d97b1741eb046d717d8f9ab92388d45960d699e8a5bc8d98a394926f07f3ec962f520fbe604caa2bbe4afa1dfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
