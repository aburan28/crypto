//! `F_p` and `F_{p²} = F_p(i)`, `i² = −1`, for primes `p ≡ 3 (mod 4)` below `2⁶³`.
//!
//! Toy-parameter arithmetic: a `u64` residue with `u128` intermediates.
//! Nothing here is constant-time; see the crate README for the security
//! posture.  `p ≡ 3 (mod 4)` makes `−1` a non-residue, so `x² + 1` is
//! irreducible and `F_{p²} = F_p[i]/(i² + 1)`, and it gives the closed-form
//! square root `a^{(p+1)/4}` in `F_p`.

use std::fmt;

/// The prime field modulus with its derived constants.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Prime {
    pub p: u64,
}

impl Prime {
    /// Construct, checking `p ≡ 3 (mod 4)`, `p < 2⁶³`, and primality
    /// (deterministic Miller–Rabin for 64-bit inputs).
    pub fn new(p: u64) -> Result<Self, String> {
        if p >= 1u64 << 63 {
            return Err("p must be below 2^63".into());
        }
        if p % 4 != 3 {
            return Err("p must be 3 mod 4".into());
        }
        if !is_prime_u64(p) {
            return Err("p is not prime".into());
        }
        Ok(Self { p })
    }

    pub fn elem(&self, v: u64) -> Fp {
        Fp {
            v: v % self.p,
            p: self.p,
        }
    }

    pub fn zero(&self) -> Fp {
        self.elem(0)
    }

    pub fn one(&self) -> Fp {
        self.elem(1)
    }

    pub fn fp2(&self, a: u64, b: u64) -> Fp2 {
        Fp2 {
            a: self.elem(a),
            b: self.elem(b),
        }
    }

    pub fn fp2_zero(&self) -> Fp2 {
        self.fp2(0, 0)
    }

    pub fn fp2_one(&self) -> Fp2 {
        self.fp2(1, 0)
    }

    /// The element `i` with `i² = −1`.
    pub fn fp2_i(&self) -> Fp2 {
        self.fp2(0, 1)
    }
}

/// Deterministic Miller–Rabin, correct for all `n < 2⁶⁴`.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'witness: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow_mod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

#[inline]
pub fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

pub fn pow_mod(mut base: u64, mut e: u64, m: u64) -> u64 {
    let mut acc = 1u64 % m;
    base %= m;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod(acc, base, m);
        }
        base = mul_mod(base, base, m);
        e >>= 1;
    }
    acc
}

/// An element of `F_p`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fp {
    pub v: u64,
    pub p: u64,
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.v)
    }
}

impl Fp {
    #[inline]
    pub fn add(self, o: Fp) -> Fp {
        let s = self.v as u128 + o.v as u128;
        let p = self.p as u128;
        Fp {
            v: (if s >= p { s - p } else { s }) as u64,
            p: self.p,
        }
    }

    #[inline]
    pub fn sub(self, o: Fp) -> Fp {
        Fp {
            v: if self.v >= o.v {
                self.v - o.v
            } else {
                self.v + self.p - o.v
            },
            p: self.p,
        }
    }

    #[inline]
    pub fn neg(self) -> Fp {
        Fp {
            v: if self.v == 0 { 0 } else { self.p - self.v },
            p: self.p,
        }
    }

    #[inline]
    pub fn mul(self, o: Fp) -> Fp {
        Fp {
            v: mul_mod(self.v, o.v, self.p),
            p: self.p,
        }
    }

    pub fn square(self) -> Fp {
        self.mul(self)
    }

    pub fn pow(self, e: u64) -> Fp {
        Fp {
            v: pow_mod(self.v, e, self.p),
            p: self.p,
        }
    }

    /// Multiplicative inverse by Fermat; `None` for zero.
    pub fn inv(self) -> Option<Fp> {
        if self.v == 0 {
            None
        } else {
            Some(self.pow(self.p - 2))
        }
    }

    pub fn is_zero(self) -> bool {
        self.v == 0
    }

    /// Euler criterion: is `self` a square in `F_p`?  Zero counts as a square.
    pub fn is_square(self) -> bool {
        self.v == 0 || self.pow((self.p - 1) / 2).v == 1
    }

    /// A square root, `None` if `self` is a non-residue.  `p ≡ 3 (mod 4)`.
    pub fn sqrt(self) -> Option<Fp> {
        let r = self.pow((self.p + 1) / 4);
        if r.square() == self {
            Some(r)
        } else {
            None
        }
    }

    pub fn from_i64(v: i64, p: u64) -> Fp {
        let m = v.rem_euclid(p as i64) as u64;
        Fp { v: m, p }
    }
}

/// An element `a + b·i` of `F_{p²}`.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fp2 {
    pub a: Fp,
    pub b: Fp,
}

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}+{}i", self.a.v, self.b.v)
    }
}

impl Fp2 {
    pub fn from_fp(a: Fp) -> Fp2 {
        Fp2 {
            a,
            b: Fp { v: 0, p: a.p },
        }
    }

    #[inline]
    pub fn add(self, o: Fp2) -> Fp2 {
        Fp2 {
            a: self.a.add(o.a),
            b: self.b.add(o.b),
        }
    }

    #[inline]
    pub fn sub(self, o: Fp2) -> Fp2 {
        Fp2 {
            a: self.a.sub(o.a),
            b: self.b.sub(o.b),
        }
    }

    pub fn neg(self) -> Fp2 {
        Fp2 {
            a: self.a.neg(),
            b: self.b.neg(),
        }
    }

    /// `(a + bi)(c + di) = (ac − bd) + (ad + bc) i`.
    #[inline]
    pub fn mul(self, o: Fp2) -> Fp2 {
        let ac = self.a.mul(o.a);
        let bd = self.b.mul(o.b);
        let ad = self.a.mul(o.b);
        let bc = self.b.mul(o.a);
        Fp2 {
            a: ac.sub(bd),
            b: ad.add(bc),
        }
    }

    pub fn square(self) -> Fp2 {
        self.mul(self)
    }

    pub fn mul_fp(self, c: Fp) -> Fp2 {
        Fp2 {
            a: self.a.mul(c),
            b: self.b.mul(c),
        }
    }

    pub fn mul_small(self, c: u64) -> Fp2 {
        let c = Fp {
            v: c % self.a.p,
            p: self.a.p,
        };
        self.mul_fp(c)
    }

    pub fn is_zero(self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    /// Conjugate `a − bi`, which is also the Frobenius `x ↦ x^p`.
    pub fn conj(self) -> Fp2 {
        Fp2 {
            a: self.a,
            b: self.b.neg(),
        }
    }

    /// The norm `a² + b²` down to `F_p`.
    pub fn norm(self) -> Fp {
        self.a.square().add(self.b.square())
    }

    pub fn inv(self) -> Option<Fp2> {
        let n = self.norm().inv()?;
        Some(self.conj().mul_fp(n))
    }

    pub fn pow(self, mut e: u128) -> Fp2 {
        let mut acc = Fp2::from_fp(Fp { v: 1, p: self.a.p });
        let mut base = self;
        while e > 0 {
            if e & 1 == 1 {
                acc = acc.mul(base);
            }
            base = base.square();
            e >>= 1;
        }
        acc
    }

    /// Frobenius `x ↦ x^p`; equals conjugation on `F_{p²}`.
    pub fn frobenius(self) -> Fp2 {
        self.conj()
    }

    /// Is `self` a square in `F_{p²}`?  `x` is a square iff `N(x) = x^{p+1}`
    /// is a square in `F_p`.
    pub fn is_square(self) -> bool {
        self.norm().is_square()
    }

    /// A square root in `F_{p²}`, `p ≡ 3 (mod 4)`.
    ///
    /// For `x = a + bi`: `N = a² + b²` must be a square in `F_p`; with
    /// `n = √N`, one of `(a + n)/2` or `(a − n)/2` is a square `α²` in
    /// `F_p` (when it is nonzero), and then `√x = α + (b/(2α)) i`.  The
    /// case `b = 0` with `a` a non-residue gives `√(−a)·i`.
    pub fn sqrt(self) -> Option<Fp2> {
        let p = self.a.p;
        if self.is_zero() {
            return Some(self);
        }
        if self.b.is_zero() {
            if let Some(r) = self.a.sqrt() {
                return Some(Fp2::from_fp(r));
            }
            let r = self.a.neg().sqrt()?;
            return Some(Fp2 {
                a: Fp { v: 0, p },
                b: r,
            });
        }
        let n = self.norm().sqrt()?;
        let half = Fp { v: 2, p }.inv().expect("2 invertible");
        for s in [n, n.neg()] {
            let t = self.a.add(s).mul(half);
            if t.is_zero() {
                continue;
            }
            if let Some(alpha) = t.sqrt() {
                let beta = self.b.mul(half).mul(alpha.inv().expect("alpha nonzero"));
                let r = Fp2 { a: alpha, b: beta };
                if r.square() == self {
                    return Some(r);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 810_647_932_926_689_279; // 2^54 * 3^2 * 5 - 1

    #[test]
    fn prime_constants() {
        let pr = Prime::new(P).unwrap();
        assert_eq!(pr.p % 4, 3);
        assert!(Prime::new(P + 2).is_err());
        assert!(Prime::new(1_000_003).is_ok()); // 1000003 = 3 mod 4
    }

    #[test]
    fn fp_field_laws() {
        let pr = Prime::new(P).unwrap();
        let a = pr.elem(123_456_789_123);
        let b = pr.elem(987_654_321_987);
        assert_eq!(a.add(b), b.add(a));
        assert_eq!(a.mul(b), b.mul(a));
        assert_eq!(a.mul(a.inv().unwrap()), pr.one());
        assert_eq!(a.sub(a), pr.zero());
        let s = a.square();
        assert!(s.is_square());
        let r = s.sqrt().unwrap();
        assert!(r == a || r == a.neg());
        assert!(!pr.elem(P - 1).is_square()); // -1 is a non-residue
    }

    #[test]
    fn fp2_field_laws_and_sqrt() {
        let pr = Prime::new(P).unwrap();
        let i = pr.fp2_i();
        assert_eq!(i.square(), pr.fp2_one().neg());
        let x = pr.fp2(55_555, 777_777_777);
        let y = pr.fp2(1234, 5678);
        assert_eq!(x.mul(y), y.mul(x));
        assert_eq!(x.mul(x.inv().unwrap()), pr.fp2_one());
        assert_eq!(x.pow(P as u128), x.frobenius());
        assert_eq!(x.pow((P as u128) * (P as u128)), x);
        let s = x.square();
        let r = s.sqrt().unwrap();
        assert!(r == x || r == x.neg());
        // every element of F_p is a square in F_{p^2}
        let m1 = Fp2::from_fp(pr.elem(P - 1));
        assert_eq!(m1.sqrt().unwrap().square(), m1);
        // a random-ish non-square
        let mut found = false;
        for t in 1..50u64 {
            let z = pr.fp2(t, t * 7 + 1);
            if !z.is_square() {
                assert!(z.sqrt().is_none());
                found = true;
                break;
            }
        }
        assert!(found);
    }
}
