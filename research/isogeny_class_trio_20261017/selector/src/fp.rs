//! GF(p) for p < 2^64 and the short Weierstrass curve y^2 = x^3 + a x + b.

use crate::bsgs::Group;
use crate::nt::powmod;

#[derive(Clone, Debug)]
pub struct Fp {
    pub p: u64,
}

impl Fp {
    #[inline]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let (s, o) = a.overflowing_add(b);
        if o || s >= self.p {
            s.wrapping_sub(self.p)
        } else {
            s
        }
    }
    #[inline]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a.wrapping_sub(b).wrapping_add(self.p)
        }
    }
    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    pub fn inv(&self, a: u64) -> u64 {
        assert!(a % self.p != 0, "inverse of zero");
        let (mut r0, mut r1) = (self.p as i128, a as i128);
        let (mut s0, mut s1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        debug_assert_eq!(r0, 1);
        s0.rem_euclid(self.p as i128) as u64
    }
    pub fn pow(&self, a: u64, e: u128) -> u64 {
        powmod(a as u128, e, self.p as u128) as u64
    }
    pub fn is_square(&self, a: u64) -> bool {
        a == 0 || self.pow(a, (self.p as u128 - 1) / 2) == 1
    }
    /// Tonelli-Shanks.
    pub fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        if !self.is_square(a) {
            return None;
        }
        let p = self.p as u128;
        let (mut q, mut s) = (p - 1, 0);
        while q % 2 == 0 {
            q /= 2;
            s += 1;
        }
        let mut z = 2u64;
        while self.is_square(z) {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(z, q);
        let mut t = self.pow(a, q);
        let mut r = self.pow(a, (q + 1) / 2);
        while t != 1 {
            let mut i = 0;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.mul(b, b);
            }
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pt {
    Inf,
    A(u64, u64),
}

#[derive(Clone, Debug)]
pub struct PCurve {
    pub k: Fp,
    pub a: u64,
    pub b: u64,
}

impl PCurve {
    pub fn on_curve(&self, p: Pt) -> bool {
        match p {
            Pt::Inf => true,
            Pt::A(x, y) => {
                let k = &self.k;
                let rhs = k.add(k.mul(k.add(k.mul(x, x), self.a), x), self.b);
                k.mul(y, y) == rhs
            }
        }
    }
    pub fn lift_x(&self, x: u64) -> Option<Pt> {
        let k = &self.k;
        let rhs = k.add(k.mul(k.add(k.mul(x, x), self.a), x), self.b);
        k.sqrt(rhs).map(|y| Pt::A(x, y))
    }
    pub fn dbl(&self, p: Pt) -> Pt {
        match p {
            Pt::Inf => Pt::Inf,
            Pt::A(x, y) => {
                if y == 0 {
                    return Pt::Inf;
                }
                let k = &self.k;
                let num = k.add(k.mul(3, k.mul(x, x)), self.a);
                let l = k.mul(num, k.inv(k.add(y, y)));
                let x3 = k.sub(k.mul(l, l), k.add(x, x));
                let y3 = k.sub(k.mul(l, k.sub(x, x3)), y);
                Pt::A(x3, y3)
            }
        }
    }
    pub fn add_pts(&self, p: Pt, q: Pt) -> Pt {
        match (p, q) {
            (Pt::Inf, _) => q,
            (_, Pt::Inf) => p,
            (Pt::A(x1, y1), Pt::A(x2, y2)) => {
                if x1 == x2 {
                    if y1 == y2 {
                        return self.dbl(p);
                    }
                    return Pt::Inf;
                }
                let k = &self.k;
                let l = k.mul(k.sub(y2, y1), k.inv(k.sub(x2, x1)));
                let x3 = k.sub(k.sub(k.mul(l, l), x1), x2);
                let y3 = k.sub(k.mul(l, k.sub(x1, x3)), y1);
                Pt::A(x3, y3)
            }
        }
    }
}

impl Group for PCurve {
    type P = Pt;
    fn inf(&self) -> Pt {
        Pt::Inf
    }
    fn is_inf(&self, p: Pt) -> bool {
        p == Pt::Inf
    }
    fn add(&self, p: Pt, q: Pt) -> Pt {
        self.add_pts(p, q)
    }
    fn neg(&self, p: Pt) -> Pt {
        match p {
            Pt::Inf => Pt::Inf,
            Pt::A(x, y) => Pt::A(x, self.k.sub(0, y)),
        }
    }
    fn mul(&self, p: Pt, n: u128) -> Pt {
        let mut r = Pt::Inf;
        if n == 0 {
            return r;
        }
        for i in (0..128 - n.leading_zeros()).rev() {
            r = self.dbl(r);
            if (n >> i) & 1 == 1 {
                r = self.add_pts(r, p);
            }
        }
        r
    }
    fn key(&self, p: Pt) -> u64 {
        match p {
            Pt::Inf => u64::MAX,
            Pt::A(x, _) => x,
        }
    }
}

impl Group for crate::f2m::BinCurve {
    type P = crate::f2m::BPt;
    fn inf(&self) -> Self::P {
        crate::f2m::BPt::Inf
    }
    fn is_inf(&self, p: Self::P) -> bool {
        p == crate::f2m::BPt::Inf
    }
    fn add(&self, p: Self::P, q: Self::P) -> Self::P {
        crate::f2m::BinCurve::add(self, p, q)
    }
    fn neg(&self, p: Self::P) -> Self::P {
        crate::f2m::BinCurve::neg(self, p)
    }
    fn mul(&self, p: Self::P, n: u128) -> Self::P {
        crate::f2m::BinCurve::mul(self, p, n)
    }
    fn key(&self, p: Self::P) -> u64 {
        match p {
            crate::f2m::BPt::Inf => u64::MAX,
            crate::f2m::BPt::A(x, _) => x,
        }
    }
}
