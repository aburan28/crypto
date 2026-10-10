//! Native finite extensions, validated by the unchanged polynomial arithmetic.
use isogeny_algos::{
    bigint::Big,
    field::{Field, Rng},
    poly,
};

#[derive(Clone)]
pub struct Extension<F: Field, const D: usize> {
    pub base: F,
    pub modulus: Vec<F::E>,
}
impl<F: Field, const D: usize> Extension<F, D> {
    pub fn new(base: F, modulus: Vec<F::E>) -> Result<Self, String> {
        if D < 2 || modulus.len() != D + 1 || modulus[D] != base.one() {
            return Err("invalid monic extension modulus".into());
        }
        let x = vec![base.zero(), base.one()];
        let mut xp = x.clone();
        let p = base.q();
        let divisors: Vec<_> = (2..=D)
            .filter(|r| D % r == 0 && (2..*r).all(|k| r % k != 0))
            .collect();
        for k in 1..=D {
            xp = poly::powmod_big(&base, &xp, &p, &modulus);
            if divisors.iter().any(|r| k == D / r)
                && poly::gcd(&base, &poly::sub(&base, &xp, &x), &modulus).len() != 1
            {
                return Err("reducible extension modulus".into());
            }
        }
        if xp != x {
            return Err("extension fails Frobenius identity".into());
        }
        Ok(Self { base, modulus })
    }
    pub fn search(base: F, seed: u64) -> Self {
        let mut rng = Rng::new(seed);
        for _ in 0..1024 {
            let mut m: Vec<_> = (0..D).map(|_| base.random(&mut rng)).collect();
            m.push(base.one());
            if let Ok(f) = Self::new(base.clone(), m) {
                return f;
            }
        }
        panic!("no irreducible modulus within the declared search bound")
    }
    pub fn embed(&self, a: F::E) -> [F::E; D] {
        let mut out = [self.base.zero(); D];
        out[0] = a;
        out
    }
    fn element(&self, v: &[F::E]) -> [F::E; D] {
        let mut out = [self.base.zero(); D];
        out[..v.len()].copy_from_slice(v);
        out
    }
}
impl<F: Field, const D: usize> Field for Extension<F, D> {
    type E = [F::E; D];
    fn zero(&self) -> Self::E {
        [self.base.zero(); D]
    }
    fn one(&self) -> Self::E {
        self.embed(self.base.one())
    }
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        std::array::from_fn(|i| self.base.add(a[i], b[i]))
    }
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        std::array::from_fn(|i| self.base.sub(a[i], b[i]))
    }
    fn neg(&self, a: Self::E) -> Self::E {
        std::array::from_fn(|i| self.base.neg(a[i]))
    }
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        let mut c = vec![self.base.zero(); 2 * D - 1];
        for i in 0..D {
            for j in 0..D {
                c[i + j] = self.base.add(c[i + j], self.base.mul(a[i], b[j]));
            }
        }
        for k in (D..c.len()).rev() {
            let top = c[k];
            for j in 0..D {
                c[k - D + j] = self
                    .base
                    .sub(c[k - D + j], self.base.mul(top, self.modulus[j]));
            }
        }
        self.element(&c[..D])
    }
    fn inv(&self, a: Self::E) -> Self::E {
        assert!(!self.is_zero(a), "zero extension inverse");
        let mut r0 = self.modulus.clone();
        let mut r1 = a.to_vec();
        poly::trim(&self.base, &mut r1);
        let mut t0 = vec![];
        let mut t1 = vec![self.base.one()];
        while !r1.is_empty() {
            let (q, r) = poly::divrem(&self.base, &r0, &r1);
            let t = poly::sub(&self.base, &t0, &poly::mul(&self.base, &q, &t1));
            r0 = r1;
            r1 = r;
            t0 = t1;
            t1 = t;
        }
        assert_eq!(r0.len(), 1);
        let t = poly::scale(&self.base, &t0, self.base.inv(r0[0]));
        self.element(&poly::rem(&self.base, &t, &self.modulus))
    }
    fn from_u64(&self, n: u64) -> Self::E {
        self.embed(self.base.from_u64(n))
    }
    fn char(&self) -> u64 {
        self.base.char()
    }
    fn characteristic_exceeds(&self, n: u64) -> bool {
        self.base.characteristic_exceeds(n)
    }
    fn q(&self) -> Big {
        let p = self.base.q();
        (0..D).fold(Big::from_u64(1), |v, _| v.mul(&p))
    }
    fn size(&self) -> u128 {
        self.q().to_u128().unwrap_or(u128::MAX)
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        std::array::from_fn(|_| self.base.random(rng))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use isogeny_algos::{field::Zp, fpm::FpM};
    #[test]
    fn rejects_a_reducible_modulus() {
        let f = Zp::new(101);
        assert!(Extension::<_, 3>::new(f, vec![100, 0, 0, 1]).is_err());
    }
    fn exercise<F: Field, const D: usize>(f: Extension<F, D>) {
        let mut rng = Rng::new(991);
        for _ in 0..24 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            let c = f.random(&mut rng);
            assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
            assert_eq!(f.pow_big(a, &f.q()), a);
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a)), f.one());
            }
        }
    }
    #[test]
    fn small_extension_arithmetic() {
        exercise(Extension::<_, 3>::search(Zp::new(103), 11));
        exercise(Extension::<_, 5>::search(Zp::new(101), 12));
    }
    #[test]
    fn large_prime_extension_arithmetic() {
        let base = FpM::<3>::from_dec("6277101735386680763835789423207666416083908700390324961279");
        exercise(Extension::<_, 3>::search(base, 19));
        let base = FpM::<4>::from_dec(
            "26959946667150639794667015087019630673557916260026308143510066298881",
        );
        exercise(Extension::<_, 5>::search(base, 23));
    }

    #[test]
    fn extension_point_count_matches_the_known_base_trace() {
        use isogeny_algos::curve::{rhs, Curve};
        let f = Extension::<_, 3>::search(Zp::new(7), 31);
        let curve = Curve::new(f.from_u64(1), f.from_u64(1));
        // Exhaustive base-field count gives #E(F_7)=5, hence t=3.
        // The third trace is 3^3-3*7*3=-36, so #E(F_343)=380.
        let mut count = 1;
        for a in 0..7 {
            for b in 0..7 {
                for c in 0..7 {
                    let value = rhs(&f, &curve, [a, b, c]);
                    count += if f.is_zero(value) {
                        1
                    } else if f.sqrt(value).is_some() {
                        2
                    } else {
                        0
                    };
                }
            }
        }
        assert_eq!(count, 380);
    }
}
