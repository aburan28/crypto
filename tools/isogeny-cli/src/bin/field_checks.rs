#[path = "../verification_field.rs"]
mod optimized;
#[path = "../../../../src/cryptanalysis/isogeny_walk/field.rs"]
mod original;
use num_bigint::BigUint;
use num_traits::{One, Zero};
fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    fn check<const N: usize>(p: BigUint) {
        let f = optimized::PrimeField::<N>::new(&p).unwrap();
        let reference = (p.bits() <= 256).then(|| original::Field::new(&p).unwrap());
        let mut state = 0xabcde123u64;
        for i in 0..300 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let a = (BigUint::from(state) << (i % p.bits() as usize)) % &p;
            let b = (&p - BigUint::one() + BigUint::from(state)) % &p;
            let aa = f.from_big(&a);
            let bb = f.from_big(&b);
            assert_eq!(f.to_big(&f.mul(&aa, &bb)), (&a * &b) % &p);
            assert_eq!(f.to_big(&f.add(&aa, &bb)), (&a + &b) % &p);
            if !a.is_zero() {
                let inverse = f.inv(&aa).unwrap();
                assert_eq!(f.mul(&aa, &inverse), f.one());
                assert_eq!(f.to_big(&inverse), a.modpow(&(&p - 2u8), &p));
                if let Some(r) = &reference {
                    assert_eq!(
                        r.to_big(&r.inv(&r.from_big(&a)).unwrap()),
                        f.to_big(&inverse)
                    );
                }
            }
        }
        assert!(f.inv(&f.zero()).is_none());
        let last = f.from_big(&(&p - 1u8));
        assert_eq!(f.to_big(&f.inv(&last).unwrap()), &p - 1u8);
    }
    #[test]
    fn small_and_full_width_fields() {
        check::<1>(BigUint::from(1009u32));
        check::<1>(BigUint::from(18446744073709551557u64));
        check::<2>((BigUint::one() << 127) - 1u8);
        check::<3>((BigUint::one() << 192) - (BigUint::one() << 64) - 1u8);
        check::<4>((BigUint::one() << 256) - (BigUint::one() << 32) - 977u32);
    }
    #[test]
    fn wider_standard_fields() {
        check::<6>(
            (BigUint::one() << 384) - (BigUint::one() << 128) - (BigUint::one() << 96)
                + (BigUint::one() << 32)
                - 1u8,
        );
        check::<8>((BigUint::one() << 512) - 569u16);
        check::<9>((BigUint::one() << 521) - 1u8);
    }
}
