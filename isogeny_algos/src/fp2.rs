//! F_{p^2} = F_p[i]/(i^2 + 1) over any prime field `F` with p = 3 mod 4 (e.g. `FpM<N>`), so that
//! supersingular-curve work (SIDH-type chains, Kani chains) runs at cryptographic sizes.
//! Products by Karatsuba (3 base multiplications), squares with 2, inversion through the norm.
use crate::bigint::Big;
use crate::field::{Field, Rng};

#[derive(Clone, Debug)]
pub struct Fp2<F: Field> {
    pub base: F,
}

impl<F: Field> Fp2<F> {
    pub fn new(base: F) -> Self {
        let q = base.q();
        assert!(q.bit(0) && q.bit(1), "Fp2: needs p = 3 mod 4 (i^2 = -1)");
        Fp2 { base }
    }
    pub fn embed(&self, a: F::E) -> (F::E, F::E) {
        (a, self.base.zero())
    }
    pub fn conj(&self, a: (F::E, F::E)) -> (F::E, F::E) {
        (a.0, self.base.neg(a.1))
    }
}

impl<F: Field> Field for Fp2<F> {
    type E = (F::E, F::E);
    fn zero(&self) -> Self::E {
        (self.base.zero(), self.base.zero())
    }
    fn one(&self) -> Self::E {
        (self.base.one(), self.base.zero())
    }
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.base.add(a.0, b.0), self.base.add(a.1, b.1))
    }
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.base.sub(a.0, b.0), self.base.sub(a.1, b.1))
    }
    fn neg(&self, a: Self::E) -> Self::E {
        (self.base.neg(a.0), self.base.neg(a.1))
    }
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        let f = &self.base;
        let t0 = f.mul(a.0, b.0);
        let t1 = f.mul(a.1, b.1);
        let t2 = f.mul(f.add(a.0, a.1), f.add(b.0, b.1));
        (f.sub(t0, t1), f.sub(f.sub(t2, t0), t1))
    }
    fn sq(&self, a: Self::E) -> Self::E {
        let f = &self.base;
        let t = f.mul(a.0, a.1);
        (f.mul(f.add(a.0, a.1), f.sub(a.0, a.1)), f.add(t, t))
    }
    fn inv(&self, a: Self::E) -> Self::E {
        let f = &self.base;
        let n = f.add(f.sq(a.0), f.sq(a.1));
        let ni = f.inv(n);
        (f.mul(a.0, ni), f.neg(f.mul(a.1, ni)))
    }
    fn from_u64(&self, n: u64) -> Self::E {
        (self.base.from_u64(n), self.base.zero())
    }
    fn char(&self) -> u64 {
        self.base.char()
    }
    fn size(&self) -> u128 {
        let p = self.base.size();
        p.checked_mul(p)
            .expect("Fp2::size: p^2 exceeds 128 bits; use q()")
    }
    fn q(&self) -> Big {
        let p = self.base.q();
        p.mul(&p)
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        (self.base.random(rng), self.base.random(rng))
    }
}
