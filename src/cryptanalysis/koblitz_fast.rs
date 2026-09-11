//! Single-word point arithmetic for the hot loops of the Koblitz index
//! calculus: the pair-sum table, the `|F|` subtractions of an `m = 3`
//! decomposition, the probe `[a]G + [b]Q` of every trial, and the
//! signed-Frobenius ρ walk the pipeline is measured against.
//!
//! The general [`BinaryPoint`] carries `Vec<u64>`-backed coordinates and
//! allocates on every field operation, which for the `n ≤ 62` fields the
//! pipeline actually runs in costs two orders of magnitude more than the
//! arithmetic.  Every operation here is a few carry-less multiplies on a
//! [`Gf2`] and never allocates; a slice of additions sharing one summand
//! is done with a single field inversion by Montgomery's trick.
//!
//! The formulas are exactly those of [`crate::binary_ecc::curve`]
//! (`y² + xy = x³ + ax² + b`), and the tests below check every operation
//! against that reference on random points, so the two representations
//! are interchangeable: [`FastCurve::lift`] and [`FastCurve::lower`]
//! convert.  Results handed back to the pipeline are still re-verified
//! in the general arithmetic where the code already did so.

use num_bigint::BigUint;

use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement};
use crate::cryptanalysis::semaev_decomp::Gf2;

/// A point with single-word coordinates, or `O`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FastPoint {
    pub x: u64,
    pub y: u64,
    pub infinity: bool,
}

impl FastPoint {
    /// The point at infinity.
    pub const INFINITY: Self = Self {
        x: 0,
        y: 0,
        infinity: true,
    };

    /// An affine point.
    pub fn affine(x: u64, y: u64) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }

    /// Packed identity, identical to
    /// [`crate::cryptanalysis::koblitz_index_calculus::pack_point`] on the
    /// lowered point: `0` for `O`, otherwise `2(x + 1) + s` with the sign
    /// bit `s` separating `P` from `−P`.
    #[inline]
    pub fn pack(self) -> u64 {
        if self.infinity {
            0
        } else {
            let sign = u64::from(self.y > (self.x ^ self.y));
            ((self.x + 1) << 1) | sign
        }
    }
}

/// `y² + xy = x³ + ax² + b` over `F_{2^n}`, `n ≤ 62`, in one word per
/// coordinate.
#[derive(Clone, Debug)]
pub struct FastCurve {
    pub field: Gf2,
    pub n: u32,
    pub a: u64,
    pub b: u64,
}

impl FastCurve {
    /// Largest extension degree handled: the packed identity needs
    /// `2(x + 1)` to fit in a word.
    pub const MAX_DEGREE: u32 = 62;

    /// `None` when the field is too wide.
    pub fn new(curve: &BinaryCurve) -> Option<Self> {
        if curve.m > Self::MAX_DEGREE || curve.irreducible.degree != curve.m {
            return None;
        }
        let field = Gf2::new(&curve.irreducible);
        let a = field.from_element(&curve.a);
        let b = field.from_element(&curve.b);
        Some(Self {
            field,
            n: curve.m,
            a,
            b,
        })
    }

    /// Project a general point down.
    pub fn lift(&self, p: &BinaryPoint) -> FastPoint {
        match p {
            BinaryPoint::Infinity => FastPoint::INFINITY,
            BinaryPoint::Affine { x, y } => {
                FastPoint::affine(self.field.from_element(x), self.field.from_element(y))
            }
        }
    }

    /// Lift back to the general representation.
    pub fn lower(&self, p: FastPoint) -> BinaryPoint {
        if p.infinity {
            BinaryPoint::Infinity
        } else {
            BinaryPoint::Affine {
                x: self.element(p.x),
                y: self.element(p.y),
            }
        }
    }

    fn element(&self, v: u64) -> F2mElement {
        F2mElement::from_biguint(&BigUint::from(v), self.n)
    }

    /// Whether `p` satisfies the curve equation.
    pub fn is_on_curve(&self, p: FastPoint) -> bool {
        if p.infinity {
            return true;
        }
        let f = &self.field;
        let lhs = f.sqr(p.y) ^ f.mul(p.x, p.y);
        let x2 = f.sqr(p.x);
        let rhs = f.mul(x2, p.x) ^ f.mul(self.a, x2) ^ self.b;
        lhs == rhs
    }

    /// `−P = (x, x + y)`.
    #[inline]
    pub fn neg(&self, p: FastPoint) -> FastPoint {
        if p.infinity {
            p
        } else {
            FastPoint::affine(p.x, p.x ^ p.y)
        }
    }

    /// `[2]P`.
    pub fn double(&self, p: FastPoint) -> FastPoint {
        if p.infinity || p.x == 0 {
            return FastPoint::INFINITY;
        }
        let f = &self.field;
        let lambda = p.x ^ f.mul(p.y, f.inv(p.x));
        self.double_with_lambda(p, lambda)
    }

    #[inline]
    fn double_with_lambda(&self, p: FastPoint, lambda: u64) -> FastPoint {
        let f = &self.field;
        let x3 = f.sqr(lambda) ^ lambda ^ self.a;
        let y3 = f.sqr(p.x) ^ f.mul(lambda ^ 1, x3);
        FastPoint::affine(x3, y3)
    }

    #[inline]
    fn add_with_lambda(&self, p: FastPoint, q: FastPoint, lambda: u64) -> FastPoint {
        let f = &self.field;
        let x3 = f.sqr(lambda) ^ lambda ^ p.x ^ q.x ^ self.a;
        let y3 = f.mul(lambda, p.x ^ x3) ^ x3 ^ p.y;
        FastPoint::affine(x3, y3)
    }

    /// `P + Q`.
    pub fn add(&self, p: FastPoint, q: FastPoint) -> FastPoint {
        if p.infinity {
            return q;
        }
        if q.infinity {
            return p;
        }
        if p.x == q.x {
            if p.y ^ q.y == p.x {
                return FastPoint::INFINITY;
            }
            return self.double(p);
        }
        let f = &self.field;
        let lambda = f.mul(p.y ^ q.y, f.inv(p.x ^ q.x));
        self.add_with_lambda(p, q, lambda)
    }

    /// `P + Q_j` for every `j`, appended to `out`, with **one** field
    /// inversion for the whole slice (Montgomery's trick over the
    /// denominators `x_P + x_{Q_j}`).  Degenerate pairs (an `O`
    /// operand, `Q_j = ±P`) are handled individually.
    pub fn add_many(
        &self,
        p: FastPoint,
        qs: &[FastPoint],
        out: &mut Vec<FastPoint>,
        scratch: &mut BatchScratch,
    ) {
        let f = &self.field;
        scratch.dens.clear();
        scratch.dens.extend(qs.iter().map(|q| {
            if p.infinity || q.infinity || q.x == p.x {
                0
            } else {
                p.x ^ q.x
            }
        }));
        f.batch_inv(&mut scratch.dens, &mut scratch.acc);
        let mut double_p: Option<FastPoint> = None;
        out.reserve(qs.len());
        for (q, &inv) in qs.iter().zip(&scratch.dens) {
            if inv == 0 {
                let sum = if p.infinity {
                    *q
                } else if q.infinity {
                    p
                } else if p.y ^ q.y == p.x {
                    FastPoint::INFINITY
                } else {
                    *double_p.get_or_insert_with(|| self.double(p))
                };
                out.push(sum);
                continue;
            }
            let lambda = f.mul(p.y ^ q.y, inv);
            out.push(self.add_with_lambda(p, *q, lambda));
        }
    }

    /// `P_i + Q_i` for every `i`, appended to `out`, with **one** field
    /// inversion for the whole slice.  This is what lets many
    /// independent walks share the cost of the one inversion a point
    /// addition needs: with `w` walks stepping together the inversion
    /// costs `3 + 1/w` multiplications each instead of the `2n` a
    /// Fermat inversion takes.
    pub fn add_pairwise(
        &self,
        ps: &[FastPoint],
        qs: &[FastPoint],
        out: &mut Vec<FastPoint>,
        scratch: &mut BatchScratch,
    ) {
        assert_eq!(ps.len(), qs.len(), "pairwise addition needs equal slices");
        let f = &self.field;
        scratch.dens.clear();
        scratch.dens.extend(ps.iter().zip(qs).map(|(p, q)| {
            if p.infinity || q.infinity || p.x == q.x {
                0
            } else {
                p.x ^ q.x
            }
        }));
        f.batch_inv(&mut scratch.dens, &mut scratch.acc);
        out.reserve(ps.len());
        for ((p, q), &inv) in ps.iter().zip(qs).zip(&scratch.dens) {
            if inv == 0 {
                out.push(self.add(*p, *q));
                continue;
            }
            let lambda = f.mul(p.y ^ q.y, inv);
            out.push(self.add_with_lambda(*p, *q, lambda));
        }
    }

    /// `[k]P` by double-and-add.
    pub fn mul_u64(&self, p: FastPoint, k: u64) -> FastPoint {
        if k == 0 || p.infinity {
            return FastPoint::INFINITY;
        }
        let mut result = FastPoint::INFINITY;
        for i in (0..64 - k.leading_zeros()).rev() {
            result = self.double(result);
            if (k >> i) & 1 == 1 {
                result = self.add(result, p);
            }
        }
        result
    }

    /// `[k]P` for a scalar of any size.
    pub fn mul(&self, p: FastPoint, k: &BigUint) -> FastPoint {
        let bits = k.bits();
        if bits == 0 || p.infinity {
            return FastPoint::INFINITY;
        }
        let mut result = FastPoint::INFINITY;
        for i in (0..bits).rev() {
            result = self.double(result);
            if k.bit(i) {
                result = self.add(result, p);
            }
        }
        result
    }

    /// The `2^k`-power Frobenius `(x, y) ↦ (x^{2^k}, y^{2^k})`.
    #[inline]
    pub fn frobenius_k(&self, p: FastPoint, k: u32) -> FastPoint {
        if p.infinity {
            return p;
        }
        FastPoint::affine(self.field.sqr_k(p.x, k), self.field.sqr_k(p.y, k))
    }
}

/// Reusable buffers for [`FastCurve::add_many`].
#[derive(Clone, Debug, Default)]
pub struct BatchScratch {
    dens: Vec<u64>,
    acc: Vec<u64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::curve::{point_add, point_double, scalar_mul};
    use crate::cryptanalysis::koblitz_index_calculus::{pack_point, KoblitzCurve};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn random_points(kc: &KoblitzCurve, count: usize, seed: u64) -> Vec<BinaryPoint> {
        let mut rng = StdRng::seed_from_u64(seed);
        let r = kc.subgroup_order.to_u64_digits()[0];
        let mut points: Vec<BinaryPoint> = (0..count)
            .map(|_| kc.mul(kc.generator(), &BigUint::from(rng.gen_range(1..r))))
            .collect();
        // A few points outside the prime-order subgroup too.
        for _ in 0..4 {
            let k = BigUint::from(rng.gen_range(1..(1u64 << 20)));
            let x = crate::binary_ecc::F2mElement::from_biguint(&k, kc.n);
            if let Some(p) = crate::cryptanalysis::koblitz_index_calculus::points_with_x(
                &kc.curve, &x,
            )
            .into_iter()
            .next()
            {
                points.push(p);
            }
        }
        points.push(BinaryPoint::Infinity);
        points
    }

    #[test]
    fn matches_the_general_arithmetic_on_random_points() {
        for (a, n) in [(0u8, 9u32), (1, 11), (1, 15), (1, 19), (0, 31)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fc = FastCurve::new(&kc.curve).unwrap();
            let points = random_points(&kc, 24, 7 + u64::from(n));
            for p in &points {
                let fp = fc.lift(p);
                assert_eq!(&fc.lower(fp), p, "lift/lower round trip at n = {n}");
                assert!(fc.is_on_curve(fp));
                assert_eq!(fp.pack(), pack_point(p), "pack at n = {n}");
                assert_eq!(fc.lower(fc.neg(fp)), crate::binary_ecc::curve::point_neg(p));
                assert_eq!(fc.lower(fc.double(fp)), point_double(&kc.curve, p));
                assert_eq!(fc.lower(fc.frobenius_k(fp, 1)), kc.frobenius(p));
                for q in &points {
                    let fq = fc.lift(q);
                    assert_eq!(
                        fc.lower(fc.add(fp, fq)),
                        point_add(&kc.curve, p, q),
                        "add at n = {n}"
                    );
                }
                let k = BigUint::from(0x9e37_79b9u64 + u64::from(n));
                assert_eq!(fc.lower(fc.mul(fp, &k)), scalar_mul(&kc.curve, p, &k));
                assert_eq!(fc.mul(fp, &k), fc.mul_u64(fp, 0x9e37_79b9u64 + u64::from(n)));
            }
        }
    }

    #[test]
    fn batched_additions_equal_individual_ones() {
        let kc = KoblitzCurve::new(0, 19).unwrap();
        let fc = FastCurve::new(&kc.curve).unwrap();
        let points: Vec<FastPoint> = random_points(&kc, 40, 3)
            .iter()
            .map(|p| fc.lift(p))
            .collect();
        let mut scratch = BatchScratch::default();
        let mut out = Vec::new();
        for &p in &points {
            // Include P itself, −P and O among the addends.
            let mut qs = points.clone();
            qs.push(p);
            qs.push(fc.neg(p));
            qs.push(FastPoint::INFINITY);
            out.clear();
            fc.add_many(p, &qs, &mut out, &mut scratch);
            assert_eq!(out.len(), qs.len());
            for (q, sum) in qs.iter().zip(&out) {
                assert_eq!(*sum, fc.add(p, *q));
            }
        }
    }

    #[test]
    fn pairwise_additions_equal_individual_ones() {
        let kc = KoblitzCurve::new(1, 19).unwrap();
        let fc = FastCurve::new(&kc.curve).unwrap();
        let points: Vec<FastPoint> = random_points(&kc, 24, 11)
            .iter()
            .map(|p| fc.lift(p))
            .collect();
        // Degenerate pairs too: P + P, P + (−P), P + O.
        let mut ps: Vec<FastPoint> = Vec::new();
        let mut qs: Vec<FastPoint> = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            ps.push(p);
            qs.push(points[(i * 7 + 3) % points.len()]);
            ps.push(p);
            qs.push(p);
            ps.push(p);
            qs.push(fc.neg(p));
            ps.push(p);
            qs.push(FastPoint::INFINITY);
        }
        let mut out = Vec::new();
        let mut scratch = BatchScratch::default();
        fc.add_pairwise(&ps, &qs, &mut out, &mut scratch);
        assert_eq!(out.len(), ps.len());
        for ((p, q), sum) in ps.iter().zip(&qs).zip(&out) {
            assert_eq!(*sum, fc.add(*p, *q));
        }
    }

    #[test]
    fn refuses_fields_wider_than_a_word() {
        let curve = BinaryCurve::sect113r1();
        assert!(FastCurve::new(&curve).is_none());
    }
}
