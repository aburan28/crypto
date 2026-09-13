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

/// **Canonical names for Frobenius orbits of abscissae**, in a normal
/// basis where the Frobenius is a rotation.
///
/// The folded pair table needs one key per orbit of `x` under `x ↦ x²`,
/// and the obvious one — the least element of the orbit — costs `n − 1`
/// squarings, a dependency chain of spreads and reduction-table lookups.
/// In a normal basis `{β, β², β⁴, …}` the same map is a one-bit cyclic
/// rotation of the coordinate word, because `(Σ c_k β^{2^k})² =
/// Σ c_k β^{2^{k+1}}`.  So the orbit of `x` is the set of rotations of
/// its coordinate word, and the least rotation names it.
///
/// The change of basis is one `F_2`-linear map, applied as eight
/// byte-table lookups; the minimum is `n` rotations of a word.  Nothing
/// is squared and nothing is reduced.
///
/// The key it produces is *not* the least element of the orbit in the
/// polynomial basis — it is a different function of the point.  What
/// matters is only that it is constant on orbits and distinct across
/// them, which a bijective linear map followed by a rotation-invariant
/// minimum is.
#[derive(Clone, Debug)]
pub struct FrobeniusCanon {
    n: u32,
    mask: u64,
    /// `tables[i][b]`: the normal coordinates of the field element whose
    /// `i`-th byte is `b` and whose other bytes are zero.
    tables: Vec<[u64; 256]>,
}

impl FrobeniusCanon {
    /// Build the basis change, or `None` if no normal element turns up —
    /// which the normal basis theorem says will not happen, but the
    /// search is randomised and bounded rather than trusted.
    pub fn new(field: &Gf2, n: u32) -> Option<Self> {
        if n == 0 || n > 63 {
            return None;
        }
        let mask = (1u64 << n) - 1;
        // A normal element: one whose Frobenius orbit is a basis.
        let mut inverse = None;
        let mut candidate = 2u64;
        for _ in 0..4096 {
            candidate = candidate.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ (candidate >> 29);
            let gamma = candidate & mask;
            if gamma == 0 {
                continue;
            }
            // Columns of the map "normal coordinates -> field element".
            let mut column = gamma;
            let columns: Vec<u64> = (0..n)
                .map(|_| {
                    let c = column;
                    column = field.sqr(column);
                    c
                })
                .collect();
            if let Some(inv) = invert_f2(&columns, n) {
                inverse = Some(inv);
                break;
            }
        }
        let inverse = inverse?;
        // Column `j` of the inverse, so a set bit of `x` contributes one
        // XOR rather than one parity.
        let by_bit: Vec<u64> = (0..n)
            .map(|j| {
                (0..n).fold(0u64, |acc, i| acc | (((inverse[i as usize] >> j) & 1) << i))
            })
            .collect();
        let bytes = ((n + 7) / 8) as usize;
        let tables = (0..bytes)
            .map(|bi| {
                let mut table = [0u64; 256];
                for (b, slot) in table.iter_mut().enumerate() {
                    let mut acc = 0u64;
                    for t in 0..8 {
                        let j = bi * 8 + t;
                        if (b >> t) & 1 == 1 && j < n as usize {
                            acc ^= by_bit[j];
                        }
                    }
                    *slot = acc;
                }
                table
            })
            .collect();
        Some(Self { n, mask, tables })
    }

    /// The normal coordinates of a field element.
    #[inline]
    pub fn coords(&self, x: u64) -> u64 {
        let mut acc = 0u64;
        for (i, table) in self.tables.iter().enumerate() {
            acc ^= table[((x >> (8 * i)) & 0xff) as usize];
        }
        acc
    }

    /// The least rotation of the normal coordinates: a name for the
    /// Frobenius orbit of `x`, equal for every element of it and for no
    /// element outside it.
    #[inline]
    pub fn canon(&self, x: u64) -> u64 {
        let c = self.coords(x);
        let n = self.n;
        let mut best = c;
        let mut v = c;
        for _ in 1..n {
            v = ((v << 1) | (v >> (n - 1))) & self.mask;
            if v < best {
                best = v;
            }
        }
        best
    }

    /// The extension degree this was built for.
    pub fn degree(&self) -> u32 {
        self.n
    }
}

/// Invert an `n × n` matrix over `F_2` given as its columns, returning
/// the rows of the inverse; `None` when the columns are dependent.
fn invert_f2(columns: &[u64], n: u32) -> Option<Vec<u64>> {
    // Row `i` carries, in bit `k`, the `i`-th bit of column `k`.
    let mut a: Vec<u64> = (0..n)
        .map(|i| {
            (0..n).fold(0u64, |acc, k| acc | (((columns[k as usize] >> i) & 1) << k))
        })
        .collect();
    let mut inv: Vec<u64> = (0..n).map(|i| 1u64 << i).collect();
    for c in 0..n as usize {
        let pivot = (c..n as usize).find(|&r| (a[r] >> c) & 1 == 1)?;
        a.swap(c, pivot);
        inv.swap(c, pivot);
        for r in 0..n as usize {
            if r != c && (a[r] >> c) & 1 == 1 {
                a[r] ^= a[c];
                inv[r] ^= inv[c];
            }
        }
    }
    Some(inv)
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
            if let Some(p) =
                crate::cryptanalysis::koblitz_index_calculus::points_with_x(&kc.curve, &x)
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
                assert_eq!(
                    fc.mul(fp, &k),
                    fc.mul_u64(fp, 0x9e37_79b9u64 + u64::from(n))
                );
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
