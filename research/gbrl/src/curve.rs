// Short-Weierstrass curve arithmetic over Fp.
//
// Pulled out of the demos so the spurious-root filter, ECDLP scaffolding, and
// future index-calculus relation collector can share one implementation.

use crate::field::{is_qr, Fp, P};
use crate::semaev::Curve;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Point {
    Infinity,
    Affine(Fp, Fp),
}

pub fn on_curve(curve: &Curve, x: Fp, y: Fp) -> bool {
    y * y == x * x * x + curve.a * x + curve.b
}

pub fn neg(p: Point) -> Point {
    match p {
        Point::Infinity => Point::Infinity,
        Point::Affine(x, y) => Point::Affine(x, -y),
    }
}

pub fn add(curve: &Curve, p: Point, q: Point) -> Point {
    match (p, q) {
        (Point::Infinity, _) => q,
        (_, Point::Infinity) => p,
        (Point::Affine(x1, y1), Point::Affine(x2, y2)) => {
            if x1 == x2 {
                if y1 + y2 == Fp::zero() {
                    return Point::Infinity;
                }
                let two = Fp::new(2);
                let three = Fp::new(3);
                let m = (three * x1 * x1 + curve.a) * (two * y1).inv().unwrap();
                let x3 = m * m - x1 - x1;
                let y3 = m * (x1 - x3) - y1;
                Point::Affine(x3, y3)
            } else {
                let m = (y2 - y1) * (x2 - x1).inv().unwrap();
                let x3 = m * m - x1 - x2;
                let y3 = m * (x1 - x3) - y1;
                Point::Affine(x3, y3)
            }
        }
    }
}

// Scalar multiplication via double-and-add. `n` is interpreted as non-negative.
pub fn scalar_mul(curve: &Curve, p: Point, n: u64) -> Point {
    let mut acc = Point::Infinity;
    let mut base = p;
    let mut k = n;
    while k > 0 {
        if k & 1 == 1 {
            acc = add(curve, acc, base);
        }
        base = add(curve, base, base);
        k >>= 1;
    }
    acc
}

// Find one affine point at x via Tonelli-Shanks sqrt. None if the RHS is
// a non-residue (no curve point at this x). O(log p), feasible at P ~ 10^6+.
pub fn point_at_x(curve: &Curve, x: Fp) -> Option<Point> {
    let rhs = x * x * x + curve.a * x + curve.b;
    rhs.sqrt().map(|y| Point::Affine(x, y))
}

// Compute #E(F_p) without actually finding y for each x. Iterates all x,
// counts 1 + (# y with y^2 = rhs) using Euler. Adds 1 for the point at
// infinity. O(p log p), fine for p ≲ 10^7.
pub fn group_order(curve: &Curve) -> u64 {
    let mut count: u64 = 1; // point at infinity
    for xi in 0..P {
        let x = Fp(xi);
        let rhs = x * x * x + curve.a * x + curve.b;
        if rhs.is_zero() {
            count += 1; // single y = 0
        } else if is_qr(rhs) {
            count += 2;
        }
    }
    count
}

// Spurious-root filter for Semaev S_3 decompositions.
//
// A common zero (x1, x2) of {S_3(x_R, x_1, x_2), F(x_1), F(x_2)} corresponds
// to a *real* decomposition R = P_1 + P_2 iff some sign choice of y_1, y_2
// gives (x_1, y_1) + (x_2, y_2) ∈ {R, -R}. Returns the concrete (P_1, P_2)
// if so, else None. (R, -R both accepted because S_3 forgets the sign of y_R.)
pub fn filter_decomp(curve: &Curve, x1: Fp, x2: Fp, r: Point) -> Option<(Point, Point)> {
    let p1 = point_at_x(curve, x1)?;
    let p2 = point_at_x(curve, x2)?;
    let nr = neg(r);
    for &s1 in &[p1, neg(p1)] {
        for &s2 in &[p2, neg(p2)] {
            let sum = add(curve, s1, s2);
            if sum == r || sum == nr {
                return Some((s1, s2));
            }
        }
    }
    None
}

// Signed variant: given two canonical factor-base points fb1, fb2 and target R,
// returns (s1, s2) ∈ {-1, +1}² such that s1·fb1 + s2·fb2 = R, if any choice
// works. Returns None if no sign combination yields R. Caller uses (s1, s2) as
// the linear-algebra coefficients for the relation row.
pub fn filter_decomp_signed(curve: &Curve, fb1: Point, fb2: Point, r: Point) -> Option<(i64, i64)> {
    for &s1 in &[1i64, -1] {
        let p1 = if s1 == 1 { fb1 } else { neg(fb1) };
        for &s2 in &[1i64, -1] {
            let p2 = if s2 == 1 { fb2 } else { neg(fb2) };
            if add(curve, p1, p2) == r {
                return Some((s1, s2));
            }
        }
    }
    None
}
