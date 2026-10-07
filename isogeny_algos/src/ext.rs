//! F_{p^4} = F_{p^2}[s]/(s^2 - nu) over `field::Zp2` (nu a non-square in F_{p^2}), p < 2^31.
//! Used where full (p^2 - 1)-torsion is needed: a supersingular curve over F_{p^2} with
//! Frobenius -p (such as y^2 = x^3 + x for p = 3 mod 4) has E(F_{p^4}) = (Z/(p^2-1))^2.
use crate::bigint::Big;
use crate::field::{Field, Rng, Zp2};

type E2 = (u64, u64);

#[derive(Clone, Debug)]
pub struct Fp4 {
    pub f2: Zp2,
    pub nu: E2,
    /// s^p = frob_s * s
    pub frob_s: E2,
}

impl Fp4 {
    pub fn new(p: u64) -> Fp4 {
        assert!(p < (1u64 << 31), "Fp4: p^4 must fit in u128");
        let f2 = Zp2::new(p);
        let q2 = p as u128 * p as u128;
        let mut rng = Rng::new(0xF4);
        let nu = loop {
            let a = f2.random(&mut rng);
            if a != f2.zero() && f2.pow(a, (q2 - 1) / 2) != f2.one() {
                break a;
            }
        };
        let frob_s = f2.pow(nu, ((p - 1) / 2) as u128);
        Fp4 { f2, nu, frob_s }
    }
    pub fn p(&self) -> u64 {
        self.f2.base.p
    }
    pub fn embed(&self, a: E2) -> (E2, E2) {
        (a, (0, 0))
    }
    /// The element of F_{p^2} if b = 0.
    pub fn project(&self, a: (E2, E2)) -> Option<E2> {
        if a.1 == (0, 0) {
            Some(a.0)
        } else {
            None
        }
    }
    fn conj2(&self, a: E2) -> E2 {
        (a.0, (self.p() - a.1) % self.p())
    }
    /// x -> x^p.
    pub fn frob(&self, x: (E2, E2)) -> (E2, E2) {
        let f = &self.f2;
        (self.conj2(x.0), f.mul(self.conj2(x.1), self.frob_s))
    }
}

impl Field for Fp4 {
    type E = (E2, E2);
    fn zero(&self) -> Self::E {
        ((0, 0), (0, 0))
    }
    fn one(&self) -> Self::E {
        ((1, 0), (0, 0))
    }
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.f2.add(a.0, b.0), self.f2.add(a.1, b.1))
    }
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.f2.sub(a.0, b.0), self.f2.sub(a.1, b.1))
    }
    fn neg(&self, a: Self::E) -> Self::E {
        (self.f2.neg(a.0), self.f2.neg(a.1))
    }
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        let f = &self.f2;
        let ac = f.mul(a.0, b.0);
        let bd = f.mul(a.1, b.1);
        // (a + b)(c + d) - ac - bd = ad + bc
        let cross = f.sub(f.sub(f.mul(f.add(a.0, a.1), f.add(b.0, b.1)), ac), bd);
        (f.add(ac, f.mul(self.nu, bd)), cross)
    }
    fn inv(&self, a: Self::E) -> Self::E {
        let f = &self.f2;
        let n = f.sub(f.mul(a.0, a.0), f.mul(self.nu, f.mul(a.1, a.1)));
        let ni = f.inv(n);
        (f.mul(a.0, ni), f.neg(f.mul(a.1, ni)))
    }
    fn from_u64(&self, n: u64) -> Self::E {
        (self.f2.from_u64(n), (0, 0))
    }
    fn char(&self) -> u64 {
        self.p()
    }
    fn size(&self) -> u128 {
        let p = self.p() as u128;
        p * p * p * p
    }
    fn q(&self) -> Big {
        Big::from_u128(self.size())
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        (self.f2.random(rng), self.f2.random(rng))
    }
}
