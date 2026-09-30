//! The Koblitz curve `K_0 : y² + xy = x³ + 1` over `F_{2^n}`.

use crate::gf2n::Field;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Point {
    Inf,
    Aff(u64, u64),
}

pub struct Curve<'a> {
    pub f: &'a Field,
    pub order: u64,
    pub odd_order: u64,
}

impl<'a> Curve<'a> {
    pub fn new(f: &'a Field) -> Curve<'a> {
        let n = f.n;
        // Koblitz recurrence for a = 0: V_0 = 2, V_1 = -1, V_k = -V_{k-1} - 2 V_{k-2};
        // #E = 2^n + 1 - V_n.
        let (mut v0, mut v1) = (2i128, -1i128);
        for _ in 1..n {
            let v2 = -v1 - 2 * v0;
            v0 = v1;
            v1 = v2;
        }
        let expected = (1i128 << n) + 1 - v1;
        // O, (0,1), and two points per x ≠ 0 with Tr(x + x^{-2}) = 0.
        // The scan is the check that the recurrence matches an enumeration.
        // Frozen cells stop at n = 19, so n ≤ 20 keeps that check. Past it
        // the scan is 2^n field inversions and is not a check that can be run;
        // the order is the recurrence alone.
        let count = if n <= 20 {
            let mut count = 2u64;
            for x in 1..(1u64 << n) {
                let c = x ^ f.sqr(f.inv(x));
                if f.trace(c) == 0 {
                    count += 2;
                }
            }
            assert_eq!(
                count as i128, expected,
                "point count must match the Koblitz recurrence"
            );
            count
        } else {
            assert!(expected > 0 && expected <= u64::MAX as i128);
            expected as u64
        };
        assert_eq!(count % 4, 0);
        Curve {
            f,
            order: count,
            odd_order: count / 4,
        }
    }

    /// The points with abscissa `x`, if any.
    pub fn lift(&self, x: u64) -> Option<[Point; 2]> {
        let f = self.f;
        if x == 0 {
            return Some([Point::Aff(0, 1), Point::Aff(0, 1)]);
        }
        let c = x ^ f.sqr(f.inv(x));
        if f.trace(c) != 0 {
            return None;
        }
        let z = f.half_trace(c);
        debug_assert_eq!(f.sqr(z) ^ z, c);
        let y = f.mul(x, z);
        Some([Point::Aff(x, y), Point::Aff(x, y ^ x)])
    }

    pub fn on_curve(&self, p: Point) -> bool {
        match p {
            Point::Inf => true,
            Point::Aff(x, y) => {
                let f = self.f;
                f.sqr(y) ^ f.mul(x, y) == f.mul(f.sqr(x), x) ^ 1
            }
        }
    }

    pub fn neg(&self, p: Point) -> Point {
        match p {
            Point::Inf => Point::Inf,
            Point::Aff(x, y) => Point::Aff(x, x ^ y),
        }
    }

    pub fn add(&self, p: Point, q: Point) -> Point {
        let f = self.f;
        match (p, q) {
            (Point::Inf, _) => q,
            (_, Point::Inf) => p,
            (Point::Aff(x1, y1), Point::Aff(x2, y2)) => {
                if x1 == x2 {
                    if y1 != y2 || x1 == 0 {
                        return Point::Inf;
                    }
                    // doubling
                    let lam = x1 ^ f.mul(y1, f.inv(x1));
                    let x3 = f.sqr(lam) ^ lam;
                    let y3 = f.sqr(x1) ^ f.mul(lam ^ 1, x3);
                    Point::Aff(x3, y3)
                } else {
                    let lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
                    let x3 = f.sqr(lam) ^ lam ^ x1 ^ x2;
                    let y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1;
                    Point::Aff(x3, y3)
                }
            }
        }
    }

    pub fn sub(&self, p: Point, q: Point) -> Point {
        self.add(p, self.neg(q))
    }

    pub fn mul(&self, mut k: u64, p: Point) -> Point {
        let mut r = Point::Inf;
        let mut base = p;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(r, base);
            }
            base = self.add(base, base);
            k >>= 1;
        }
        r
    }

    pub fn x(&self, p: Point) -> Option<u64> {
        match p {
            Point::Inf => None,
            Point::Aff(x, _) => Some(x),
        }
    }
}
