mod common;
use common::*;
use eigen_seta::field::Field;
use eigen_seta::kani::{div_add, div_dbl, div_on_curve, Div};
use rand_core::{OsRng, RngCore};

fn rand_f(rng: &mut OsRng) -> F {
    F::from_counter(rng.next_u64()) * F::from_counter(rng.next_u64())
}

fn point_on(h: &[F; 7], rng: &mut OsRng) -> (F, F) {
    loop {
        let x = rand_f(rng);
        let mut y2 = F::zero();
        for i in (0..7).rev() {
            y2 = y2 * x + h[i];
        }
        let (y, ok) = y2.sqrt();
        if bool::from(ok) {
            return (x, y);
        }
    }
}

fn divisor(h: &[F; 7], rng: &mut OsRng) -> Div<F> {
    let (x1, y1) = point_on(h, rng);
    let (x2, y2) = point_on(h, rng);
    // u = (x - x1)(x - x2); v interpolates (x1,y1), (x2,y2)
    let u1 = -(x1 + x2);
    let u0 = x1 * x2;
    let v1 = (y2 - y1) * (x2 - x1).invert();
    let v0 = y1 - v1 * x1;
    Div { u1, u0, v1, v0 }
}

fn deq(a: &Div<F>, b: &Div<F>) -> bool {
    a.u1 == b.u1 && a.u0 == b.u0 && a.v1 == b.v1 && a.v0 == b.v0
}

#[test]
fn jacobian_laws() {
    let mut rng = OsRng;
    let mut h = [F::zero(); 7];
    for c in h.iter_mut() {
        *c = rand_f(&mut rng);
    }
    for _ in 0..20 {
        let a = divisor(&h, &mut rng);
        let b = divisor(&h, &mut rng);
        assert!(div_on_curve(&h, &a) && div_on_curve(&h, &b));
        let s = div_add(&h, &a, &b);
        assert!(div_on_curve(&h, &s), "sum off curve");
        let da = div_dbl(&h, &a);
        assert!(div_on_curve(&h, &da), "double off curve");
        let db = div_dbl(&h, &b);
        let lhs = div_dbl(&h, &s);
        let rhs = div_add(&h, &da, &db);
        assert!(deq(&lhs, &rhs), "2(a+b) != 2a + 2b");
        // (a + b) + a == 2a + b
        let l2 = div_add(&h, &s, &a);
        let r2 = div_add(&h, &da, &b);
        assert!(deq(&l2, &r2), "(a+b)+a != 2a+b");
    }
}

/// Regression: Tonelli-Shanks must not accept z = 0 as a non-residue (it
/// looped forever).  Exercise many small q = 1 mod 8 with many RNG streams.
#[test]
fn sqrt_mod_prime_small_q_terminates() {
    use eigen_seta::ntheory::sqrt_mod_prime;
    use num_bigint::BigUint;
    struct Ctr(u64);
    impl rand_core::RngCore for Ctr {
        fn next_u32(&mut self) -> u32 {
            self.next_u64() as u32
        }
        fn next_u64(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
        fn fill_bytes(&mut self, d: &mut [u8]) {
            for c in d.chunks_mut(8) {
                let v = self.next_u64().to_le_bytes();
                c.copy_from_slice(&v[..c.len()]);
            }
        }
        fn try_fill_bytes(&mut self, d: &mut [u8]) -> Result<(), rand_core::Error> {
            self.fill_bytes(d);
            Ok(())
        }
    }
    for q in [17u32, 41, 73, 89, 97, 113, 137, 193, 233, 241, 257] {
        let qb = BigUint::from(q);
        for seed in 0..200u64 {
            let mut rng = Ctr(seed);
            let a = BigUint::from(q - 1); // -1 is a square mod q = 1 mod 4
            let r = sqrt_mod_prime(&a, &qb, &mut rng).expect("square");
            assert_eq!((&r * &r) % &qb, a);
        }
    }
}
