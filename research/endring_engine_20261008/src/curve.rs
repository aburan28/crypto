//! Elliptic curves `y² = x³ + a·x + b` over `F_{p²}`, affine points with
//! full `(x, y)`, scalar multiplication, torsion bases and smooth-order
//! two-dimensional discrete logarithms.
//!
//! The curves this crate walks are supersingular over `F_{p²}` with
//! `E(F_{p²}) ≅ (Z/(p+1))²` (every curve in the isogeny class of
//! `E₀: y² = x³ + x` reached by an `F_{p²}`-rational isogeny of degree
//! dividing `p + 1` has this group structure), so for `N | p + 1` the full
//! `N`-torsion is rational and a basis of `E[N] ≅ (Z/N)²` can be sampled.

use crate::fp::{Fp2, Prime};
use rand::Rng;

/// A short Weierstrass curve over `F_{p²}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve {
    pub prime: Prime,
    pub a: Fp2,
    pub b: Fp2,
}

/// An affine point or the point at infinity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Point {
    Infinity,
    Affine { x: Fp2, y: Fp2 },
}

impl Point {
    pub fn is_infinity(&self) -> bool {
        matches!(self, Point::Infinity)
    }

    pub fn affine(x: Fp2, y: Fp2) -> Point {
        Point::Affine { x, y }
    }
}

impl Curve {
    pub fn new(prime: Prime, a: Fp2, b: Fp2) -> Result<Curve, String> {
        let c = Curve { prime, a, b };
        if c.discriminant().is_zero() {
            return Err("singular curve".into());
        }
        Ok(c)
    }

    /// `E₀: y² = x³ + x`, `j = 1728`, supersingular for `p ≡ 3 (mod 4)`.
    pub fn e0(prime: Prime) -> Curve {
        Curve {
            prime,
            a: prime.fp2_one(),
            b: prime.fp2_zero(),
        }
    }

    /// `−16(4a³ + 27b²)`, zero exactly for a singular model.
    pub fn discriminant(&self) -> Fp2 {
        let a3 = self.a.square().mul(self.a).mul_small(4);
        let b2 = self.b.square().mul_small(27);
        a3.add(b2).mul_small(16).neg()
    }

    /// `j = 1728 · 4a³ / (4a³ + 27b²)`.
    pub fn j_invariant(&self) -> Fp2 {
        let a3 = self.a.square().mul(self.a).mul_small(4);
        let den = a3.add(self.b.square().mul_small(27));
        a3.mul_small(1728).mul(den.inv().expect("nonsingular"))
    }

    pub fn is_on_curve(&self, pt: &Point) -> bool {
        match pt {
            Point::Infinity => true,
            Point::Affine { x, y } => {
                let rhs = x.square().mul(*x).add(self.a.mul(*x)).add(self.b);
                y.square() == rhs
            }
        }
    }

    pub fn neg(&self, pt: &Point) -> Point {
        match pt {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => Point::Affine { x: *x, y: y.neg() },
        }
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        let (x1, y1) = match p {
            Point::Infinity => return *q,
            Point::Affine { x, y } => (*x, *y),
        };
        let (x2, y2) = match q {
            Point::Infinity => return *p,
            Point::Affine { x, y } => (*x, *y),
        };
        let lambda = if x1 == x2 {
            if y1 != y2 || y1.is_zero() {
                return Point::Infinity;
            }
            // tangent: (3x² + a) / (2y)
            let num = x1.square().mul_small(3).add(self.a);
            let den = y1.mul_small(2);
            num.mul(den.inv().expect("y nonzero"))
        } else {
            let num = y2.sub(y1);
            let den = x2.sub(x1);
            num.mul(den.inv().expect("x1 != x2"))
        };
        let x3 = lambda.square().sub(x1).sub(x2);
        let y3 = lambda.mul(x1.sub(x3)).sub(y1);
        Point::Affine { x: x3, y: y3 }
    }

    pub fn double(&self, p: &Point) -> Point {
        self.add(p, p)
    }

    pub fn sub(&self, p: &Point, q: &Point) -> Point {
        self.add(p, &self.neg(q))
    }

    /// `[k]P` by double-and-add; `k` is a `u128` so that `(p+1)²` fits.
    pub fn mul(&self, p: &Point, mut k: u128) -> Point {
        let mut acc = Point::Infinity;
        let mut base = *p;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            base = self.double(&base);
            k >>= 1;
        }
        acc
    }

    /// `[k]P` for a signed scalar.
    pub fn mul_signed(&self, p: &Point, k: i128) -> Point {
        if k < 0 {
            self.neg(&self.mul(p, k.unsigned_abs()))
        } else {
            self.mul(p, k as u128)
        }
    }

    /// The Frobenius `π: (x, y) ↦ (x^p, y^p)`.  An endomorphism of every
    /// curve defined over `F_p`, and in particular of `E₀`.
    pub fn frobenius(&self, p: &Point) -> Point {
        match p {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => Point::Affine {
                x: x.frobenius(),
                y: y.frobenius(),
            },
        }
    }

    /// A uniformly random affine point.
    pub fn random_point<R: Rng>(&self, rng: &mut R) -> Point {
        loop {
            let x = self.prime.fp2(rng.gen::<u64>(), rng.gen::<u64>());
            let rhs = x.square().mul(x).add(self.a.mul(x)).add(self.b);
            if let Some(y) = rhs.sqrt() {
                let y = if rng.gen::<bool>() { y } else { y.neg() };
                return Point::Affine { x, y };
            }
        }
    }

    /// The order of a point whose order divides the smooth `n` with the
    /// given prime factorisation (as `(ℓ, e)` pairs).
    pub fn order_dividing(&self, p: &Point, n: u128, factors: &[(u128, u32)]) -> u128 {
        let mut order = n;
        for &(ell, e) in factors {
            for _ in 0..e {
                let candidate = order / ell;
                if self.mul(p, candidate).is_infinity() {
                    order = candidate;
                } else {
                    break;
                }
            }
        }
        order
    }

    /// A random point of exact order `n`, where `n | p + 1` and the curve
    /// has `E(F_{p²}) ≅ (Z/(p+1))²`.
    pub fn random_point_of_order<R: Rng>(
        &self,
        rng: &mut R,
        n: u128,
        factors: &[(u128, u32)],
    ) -> Point {
        let cofactor = (self.prime.p as u128 + 1) / n;
        loop {
            let r = self.random_point(rng);
            let q = self.mul(&r, cofactor);
            if self.order_dividing(&q, n, factors) == n {
                return q;
            }
        }
    }

    /// Are `P` and `Q`, both of order `n`, independent, i.e. do they
    /// generate `E[n]`?  For each prime `ℓ | n` the images in `E[ℓ]` must
    /// not be proportional; the check enumerates the `ℓ` multiples, so it
    /// is for small `ℓ`.
    pub fn independent(&self, p: &Point, q: &Point, n: u128, factors: &[(u128, u32)]) -> bool {
        for &(ell, _) in factors {
            let m = n / ell;
            let pl = self.mul(p, m);
            let ql = self.mul(q, m);
            if pl.is_infinity() || ql.is_infinity() {
                return false;
            }
            let mut acc = Point::Infinity;
            for _ in 0..ell {
                if acc == ql {
                    return false;
                }
                acc = self.add(&acc, &pl);
            }
        }
        true
    }

    /// A basis `(P, Q)` of `E[n]`, `n | p + 1`.
    pub fn torsion_basis<R: Rng>(
        &self,
        rng: &mut R,
        n: u128,
        factors: &[(u128, u32)],
    ) -> (Point, Point) {
        let p = self.random_point_of_order(rng, n, factors);
        loop {
            let q = self.random_point_of_order(rng, n, factors);
            if self.independent(&p, &q, n, factors) {
                return (p, q);
            }
        }
    }

    /// Two-dimensional discrete logarithm in `E[ℓ^e]`: with `(P, Q)` a
    /// basis, write `R = [a]P + [b]Q` and return `(a, b)` modulo `ℓ^e`.
    /// Digit by digit (Pohlig–Hellman), each digit by exhausting the `ℓ²`
    /// combinations, so it is for small `ℓ`.
    pub fn dlog2_prime_power(
        &self,
        p: &Point,
        q: &Point,
        r: &Point,
        ell: u128,
        e: u32,
    ) -> Option<(u128, u128)> {
        let n = ell.pow(e);
        // Images in E[ℓ].
        let p1 = self.mul(p, n / ell);
        let q1 = self.mul(q, n / ell);
        let mut a = 0u128;
        let mut b = 0u128;
        let mut scale = 1u128;
        for k in 0..e {
            // Residual R − [a]P − [b]Q has order dividing ℓ^{e−k}; project to E[ℓ].
            let residual = self.sub(r, &self.add(&self.mul(p, a), &self.mul(q, b)));
            let target = self.mul(&residual, n / ell.pow(k + 1));
            let mut found = None;
            'search: for da in 0..ell {
                let pa = self.mul(&p1, da);
                for db in 0..ell {
                    let cand = self.add(&pa, &self.mul(&q1, db));
                    if cand == target {
                        found = Some((da, db));
                        break 'search;
                    }
                }
            }
            let (da, db) = found?;
            a += da * scale;
            b += db * scale;
            scale *= ell;
        }
        Some((a % n, b % n))
    }

    /// Two-dimensional discrete logarithm in `E[n]` for smooth `n`, by
    /// combining the prime-power solutions with the Chinese remainder
    /// theorem.  `(P, Q)` must be a basis of `E[n]`.
    pub fn dlog2(
        &self,
        p: &Point,
        q: &Point,
        r: &Point,
        n: u128,
        factors: &[(u128, u32)],
    ) -> Option<(u128, u128)> {
        let mut a = 0u128;
        let mut b = 0u128;
        let mut modulus = 1u128;
        for &(ell, e) in factors {
            let m = ell.pow(e);
            let cof = n / m;
            let pm = self.mul(p, cof);
            let qm = self.mul(q, cof);
            let rm = self.mul(r, cof);
            let (am, bm) = self.dlog2_prime_power(&pm, &qm, &rm, ell, e)?;
            // Solve [cof·a_m... ]: R_m = cof·R = [a]cof·P + [b]cof·Q, so a ≡ a_m, b ≡ b_m mod m.
            a = crt(a, modulus, am, m);
            b = crt(b, modulus, bm, m);
            modulus *= m;
        }
        Some((a, b))
    }
}

/// Combine `x ≡ a (mod m)` and `x ≡ b (mod n)` with `gcd(m, n) = 1`.
pub fn crt(a: u128, m: u128, b: u128, n: u128) -> u128 {
    if m == 1 {
        return b % n;
    }
    // x = a + m·t, m·t ≡ b − a (mod n)
    let diff = (b % n + n - a % n) % n;
    let inv = mod_inverse(m % n, n).expect("coprime moduli");
    let t = mulmod_u128(diff, inv, n);
    a + m * t
}

/// Modular inverse by the extended Euclid; `None` if not invertible.
pub fn mod_inverse(a: u128, m: u128) -> Option<u128> {
    let (mut old_r, mut r) = (a as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    if old_r != 1 {
        return None;
    }
    Some(old_s.rem_euclid(m as i128) as u128)
}

/// `a·b mod m` without overflow for `m < 2⁶⁴`; for larger moduli it uses
/// the schoolbook double-and-add on `u128`.
pub fn mulmod_u128(a: u128, b: u128, m: u128) -> u128 {
    if m <= u64::MAX as u128 {
        return ((a % m) * (b % m)) % m;
    }
    let mut acc = 0u128;
    let mut base = a % m;
    let mut e = b % m;
    while e > 0 {
        if e & 1 == 1 {
            acc = (acc + base) % m;
        }
        base = (base << 1) % m;
        e >>= 1;
    }
    acc
}

/// The prime factorisation of a smooth `n` by trial division up to `bound`.
pub fn factor_smooth(mut n: u128, bound: u128) -> Option<Vec<(u128, u32)>> {
    let mut out = Vec::new();
    let mut q = 2u128;
    while q <= bound && q * q <= n {
        if n % q == 0 {
            let mut e = 0;
            while n % q == 0 {
                n /= q;
                e += 1;
            }
            out.push((q, e));
        }
        q += if q == 2 { 1 } else { 2 };
    }
    if n > 1 {
        if n > bound {
            return None;
        }
        out.push((n, 1));
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    const P: u64 = 810_647_932_926_689_279; // p + 1 = 2^54 · 3^2 · 5

    fn setup() -> (Curve, rand_chacha::ChaCha8Rng, Vec<(u128, u32)>) {
        let prime = Prime::new(P).unwrap();
        let e0 = Curve::e0(prime);
        let rng = rand_chacha::ChaCha8Rng::seed_from_u64(1);
        let factors = factor_smooth(P as u128 + 1, 100).unwrap();
        assert_eq!(factors, vec![(2, 54), (3, 2), (5, 1)]);
        (e0, rng, factors)
    }

    #[test]
    fn e0_group_structure() {
        let (e0, mut rng, _) = setup();
        assert_eq!(e0.j_invariant(), e0.prime.fp2(1728, 0));
        for _ in 0..8 {
            let r = e0.random_point(&mut rng);
            assert!(e0.is_on_curve(&r));
            // (p+1)·R = O since E(F_{p²}) ≅ (Z/(p+1))²
            assert!(e0.mul(&r, P as u128 + 1).is_infinity());
            // group law: 2R = R + R, (R + S) − S = R
            let s = e0.random_point(&mut rng);
            assert_eq!(e0.double(&r), e0.add(&r, &r));
            assert_eq!(e0.sub(&e0.add(&r, &s), &s), r);
            assert!(e0.is_on_curve(&e0.add(&r, &s)));
            // π² = 1 on F_{p²}-points, and π commutes with the group law
            assert_eq!(e0.frobenius(&e0.frobenius(&r)), r);
            assert_eq!(
                e0.frobenius(&e0.add(&r, &s)),
                e0.add(&e0.frobenius(&r), &e0.frobenius(&s))
            );
        }
    }

    #[test]
    fn torsion_basis_and_dlog() {
        let (e0, mut rng, _) = setup();
        let cases: Vec<(u128, Vec<(u128, u32)>)> = vec![
            (
                2u128.pow(20) * 9 * 5,
                factor_smooth(2u128.pow(20) * 45, 100).unwrap(),
            ),
            (3u128 * 5, vec![(3, 1), (5, 1)]),
            (2u128.pow(30), vec![(2, 30)]),
        ];
        for (n, factors) in cases {
            let factors = &factors;
            let (p, q) = e0.torsion_basis(&mut rng, n, factors);
            assert!(e0.mul(&p, n).is_infinity() && e0.mul(&q, n).is_infinity());
            assert_eq!(e0.order_dividing(&p, n, factors), n);
            for _ in 0..4 {
                let a = rng.gen::<u128>() % n;
                let b = rng.gen::<u128>() % n;
                let r = e0.add(&e0.mul(&p, a), &e0.mul(&q, b));
                let (a2, b2) = e0.dlog2(&p, &q, &r, n, factors).unwrap();
                assert_eq!((a2, b2), (a, b), "n = {n}");
            }
        }
    }

    #[test]
    fn crt_and_inverse() {
        assert_eq!(crt(2, 3, 3, 5), 8);
        assert_eq!(mod_inverse(3, 7), Some(5));
        assert_eq!(mod_inverse(4, 8), None);
        assert_eq!(
            mulmod_u128(u128::from(u64::MAX), 3, u128::from(u64::MAX) + 7),
            (u128::from(u64::MAX) * 3) % (u128::from(u64::MAX) + 7)
        );
    }
}
