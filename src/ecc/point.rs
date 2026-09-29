//! Elliptic curve point arithmetic over short Weierstrass curves: y² = x³ + ax + b (mod p).
//!
//! # Point addition
//! For P ≠ Q and P ≠ -Q:
//!   λ = (y₂ - y₁) / (x₂ - x₁)
//!   x₃ = λ² - x₁ - x₂
//!   y₃ = λ(x₁ - x₃) - y₁
//!
//! # Point doubling (P = Q)
//!   λ = (3x² + a) / (2y)
//!   x₃ = λ² - 2x
//!   y₃ = λ(x - x₃) - y
//!
//! # Scalar multiplication
//! Uses the double-and-add algorithm: scan bits of k from LSB to MSB,
//! doubling the accumulator each step and adding P when the bit is 1.
//! This is O(log k) point operations.

use super::field::{inv_mod_u64, single_word, FieldElement};
use num_bigint::BigUint;

/// A point on a Weierstrass curve, either the identity (point at infinity)
/// or an affine coordinate pair.
#[derive(Clone, Debug, PartialEq)]
pub enum Point {
    Infinity,
    Affine { x: FieldElement, y: FieldElement },
}

impl Point {
    /// Add two curve points. `a` is the curve coefficient from y²=x³+ax+b.
    pub fn add(&self, other: &Point, a: &FieldElement) -> Point {
        match (self, other) {
            (Point::Infinity, p) | (p, Point::Infinity) => p.clone(),
            (Point::Affine { x: x1, y: y1 }, Point::Affine { x: x2, y: y2 }) => {
                if x1 == x2 {
                    if y1 == y2 {
                        // Same point → double
                        return self.double(a);
                    } else {
                        // P + (-P) = ∞
                        return Point::Infinity;
                    }
                }
                // λ = (y2 - y1) / (x2 - x1)
                let lambda = y2.sub(y1).mul(&x2.sub(x1).inv().unwrap());
                let x3 = lambda.mul(&lambda).sub(x1).sub(x2);
                let y3 = lambda.mul(&x1.sub(&x3)).sub(y1);
                Point::Affine { x: x3, y: y3 }
            }
        }
    }

    /// Double a point: 2P.
    pub fn double(&self, a: &FieldElement) -> Point {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => {
                if y.is_zero() {
                    return Point::Infinity;
                }
                let p = x.modulus.clone();
                let three = FieldElement::new(BigUint::from(3u32), p.clone());
                let two = FieldElement::new(BigUint::from(2u32), p);

                // λ = (3x² + a) / (2y)
                let numerator = three.mul(&x.mul(x)).add(a);
                let denominator = two.mul(y);
                let lambda = numerator.mul(&denominator.inv().unwrap());

                let x3 = lambda.mul(&lambda).sub(&two.mul(x));
                let y3 = lambda.mul(&x.sub(&x3)).sub(y);
                Point::Affine { x: x3, y: y3 }
            }
        }
    }

    /// [`add`](Self::add) in **variable time**, for public points only.
    ///
    /// The same affine formulas and the same identity, doubling and
    /// `P + (−P)` cases, so over a prime `p` the result is exactly
    /// `self.add(other, a)`.  (On a modulus that is not prime, where
    /// neither computes a group law, the point can differ as
    /// [`FieldElement::inv_vartime`] says, but neither call panics where
    /// the other does not.)  What differs is how the field arithmetic
    /// runs, and its running time now depends on the coordinates:
    ///
    /// - the slope's inversion is [`FieldElement::inv_vartime`], a
    ///   Euclidean inversion, instead of the constant-time Fermat ladder
    ///   that was most of an affine addition's cost;
    /// - when `p` fits a machine word the whole formula runs on words
    ///   (the private `Word` arithmetic below) and only the result goes
    ///   back to [`BigUint`]; otherwise subtractions use
    ///   [`FieldElement::sub_vartime`].
    ///
    /// For the cryptanalysis walks, whose points are public.  Protocol
    /// code keeps [`add`](Self::add).
    pub fn add_vartime(&self, other: &Point, a: &FieldElement) -> Point {
        match (self, other) {
            (Point::Infinity, p) | (p, Point::Infinity) => p.clone(),
            (Point::Affine { x: x1, y: y1 }, Point::Affine { x: x2, y: y2 }) => {
                let m = &x1.modulus;
                if let Some(w) = Word::of(m) {
                    if let Some([x1, y1, x2, y2, a]) = w.get([x1, y1, x2, y2, a]) {
                        return w.point(w.affine_add((x1, y1), (x2, y2), a), m);
                    }
                }
                if x1 == x2 {
                    if y1 == y2 {
                        return self.double_vartime(a);
                    } else {
                        return Point::Infinity;
                    }
                }
                let lambda = y2
                    .sub_vartime(y1)
                    .mul(&x2.sub_vartime(x1).inv_vartime().unwrap());
                let x3 = lambda.mul(&lambda).sub_vartime(x1).sub_vartime(x2);
                let y3 = lambda.mul(&x1.sub_vartime(&x3)).sub_vartime(y1);
                Point::Affine { x: x3, y: y3 }
            }
        }
    }

    /// [`double`](Self::double) in **variable time**, for public points
    /// only: exactly `self.double(a)`, computed as in
    /// [`add_vartime`](Self::add_vartime).
    pub fn double_vartime(&self, a: &FieldElement) -> Point {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => {
                let m = &x.modulus;
                if let Some(w) = Word::of(m) {
                    if let Some([x, y, a]) = w.get([x, y, a]) {
                        return w.point(w.affine_double((x, y), a), m);
                    }
                }
                if y.is_zero() {
                    return Point::Infinity;
                }
                let p = x.modulus.clone();
                let three = FieldElement::new(BigUint::from(3u32), p.clone());
                let two = FieldElement::new(BigUint::from(2u32), p);
                let numerator = three.mul(&x.mul(x)).add(a);
                let denominator = two.mul(y);
                let lambda = numerator.mul(&denominator.inv_vartime().unwrap());
                let x3 = lambda.mul(&lambda).sub_vartime(&two.mul(x));
                let y3 = lambda.mul(&x.sub_vartime(&x3)).sub_vartime(y);
                Point::Affine { x: x3, y: y3 }
            }
        }
    }

    /// [`scalar_mul`](Self::scalar_mul) in **variable time**, for public
    /// points and scalars only: over a prime `p`, exactly
    /// `self.scalar_mul(k, a)` (on other moduli, as for
    /// [`add_vartime`](Self::add_vartime)).
    ///
    /// When `p` fits a word the same Jacobian ladder runs on word
    /// arithmetic and inverts once at the end with the word Euclid of
    /// [`FieldElement::inv_vartime`]; otherwise this is `scalar_mul`
    /// itself, whose one constant-time inversion is a small part of a
    /// multi-word ladder.  For the set-up multiplications of the
    /// cryptanalysis walks and searches (walker starts, subgroup
    /// projections, verification of a recovered log).
    pub fn scalar_mul_vartime(&self, k: &BigUint, a: &FieldElement) -> Point {
        if let Point::Affine { x, y } = self {
            let m = &x.modulus;
            if let Some(w) = Word::of(m) {
                if let Some([x, y, a]) = w.get([x, y, a]) {
                    let mut acc: Option<WordJacobian> = None;
                    for i in (0..k.bits()).rev() {
                        if let Some(j) = acc {
                            acc = w.jacobian_double(j, a);
                        }
                        if k.bit(i) {
                            acc = match acc {
                                None => Some(WordJacobian { x, y, z: 1 }),
                                Some(j) => w.jacobian_add_affine(j, (x, y), a),
                            };
                        }
                    }
                    return w.point(acc.map(|j| w.jacobian_to_affine(j)), m);
                }
            }
        }
        self.scalar_mul(k, a)
    }

    /// The least `k < bound` with `k·self = target`, trying
    /// `k = 0, 1, 2, …` in turn — a linear discrete-log search, as in a
    /// Pohlig–Hellman digit — or `None` if there is none below `bound`.
    /// Variable time; for public points only.
    ///
    /// The running multiple `k·self` is kept in Jacobian coordinates and
    /// compared with the affine `target = (x, y)` by cross-multiplying,
    /// `X = x·Z²` and `Y = y·Z³`, so a step costs a mixed addition and
    /// four multiplications and never an inversion.  The multiples, the
    /// identity among them, and so the `k` returned are exactly those of
    /// stepping `current = current.add(self, a)` and testing
    /// `current == *target`.  On word arithmetic when `p` fits a word, on
    /// the [`BigUint`] Jacobian points of [`scalar_mul`](Self::scalar_mul)
    /// otherwise.
    pub fn linear_dlog_vartime(&self, target: &Point, bound: u64, a: &FieldElement) -> Option<u64> {
        let (gx, gy) = match self {
            // Every multiple of O is O.
            Point::Infinity => return (*target == Point::Infinity && bound > 0).then_some(0),
            Point::Affine { x, y } => (x, y),
        };
        let m = &gx.modulus;
        // A base with an unreduced coordinate (no constructor makes one) is
        // itself the multiple `1·self` compared, unreduced, so it is
        // stepped affinely as the search this stands for does.
        if gx.value >= *m || gy.value >= *m {
            let mut current = Point::Infinity;
            for k in 0..bound {
                if current == *target {
                    return Some(k);
                }
                current = current.add_vartime(self, a);
            }
            return None;
        }
        // Every other multiple is reduced, so an affine target with an
        // unreduced coordinate equals none of them.
        let t = match target {
            Point::Infinity => None,
            Point::Affine { x, y } if x.value < *m && y.value < *m => Some((x, y)),
            Point::Affine { .. } => return None,
        };
        if let Some(w) = Word::of(m) {
            let tw = match t {
                None => Some(None),
                Some((x, y)) => w.get([x, y]).map(|[x, y]| Some((x, y))),
            };
            if let (Some([gx, gy, a]), Some(t)) = (w.get([gx, gy, a]), tw) {
                let mut acc: Option<WordJacobian> = None;
                for k in 0..bound {
                    let hit = match (acc, t) {
                        (None, None) => true,
                        (Some(j), Some(t)) => w.jacobian_is(j, t),
                        _ => false,
                    };
                    if hit {
                        return Some(k);
                    }
                    acc = match acc {
                        None => Some(WordJacobian { x: gx, y: gy, z: 1 }),
                        Some(j) => w.jacobian_add_affine(j, (gx, gy), a),
                    };
                }
                return None;
            }
        }
        let mut acc: Option<Jacobian> = None;
        for k in 0..bound {
            let hit = match (&acc, t) {
                (None, None) => true,
                (Some(j), Some((x, y))) => j.is(x, y),
                _ => false,
            };
            if hit {
                return Some(k);
            }
            acc = match &acc {
                None => Some(Jacobian::from_affine(gx, gy)),
                Some(j) => j.add_affine(gx, gy, a),
            };
        }
        None
    }

    /// Scalar multiplication kP using the double-and-add algorithm.
    ///
    /// Security note: this implementation is NOT constant-time and is
    /// vulnerable to timing side-channels.  Use [`Point::scalar_mul_ct`]
    /// for secret scalars.  This variable-time variant is appropriate
    /// only for public-input scalar multiplications, e.g. the verifier's
    /// `u1·G + u2·Q` in ECDSA verification.
    ///
    /// The ladder runs in Jacobian coordinates `(X : Y : Z)`, affine
    /// `(X/Z², Y/Z³)`, and converts back with a single inversion at the
    /// end: an affine step costs an inversion — a full modular
    /// exponentiation here — where a Jacobian step costs a dozen
    /// multiplications.  The result is the same affine point the
    /// double-and-add over [`Point::add`] and [`Point::double`] returns,
    /// including its identity and `2P = O` cases, because the two agree
    /// wherever the affine formulas are defined.
    pub fn scalar_mul(&self, k: &BigUint, a: &FieldElement) -> Point {
        let (x, y) = match self {
            Point::Infinity => return Point::Infinity,
            Point::Affine { x, y } => (x, y),
        };
        let mut acc: Option<Jacobian> = None;
        for i in (0..k.bits()).rev() {
            if let Some(j) = &acc {
                acc = j.double(a);
            }
            if k.bit(i) {
                acc = match &acc {
                    None => Some(Jacobian::from_affine(x, y)),
                    Some(j) => j.add_affine(x, y, a),
                };
            }
        }
        match acc {
            None => Point::Infinity,
            Some(j) => j.to_affine(),
        }
    }

    /// Scalar multiplication kP using the Montgomery ladder, processing
    /// `scalar_bits` bits MSB-first.  Every iteration performs exactly
    /// one point addition and one point doubling, regardless of the bit
    /// value, so the *operation count* no longer depends on the
    /// scalar's bit pattern.  Iteration count is fixed at `scalar_bits`,
    /// so the runtime no longer scales with the position of the
    /// most-significant set bit.
    ///
    /// Security caveats (read carefully):
    ///
    /// - The branch that re-binds `(r0, r1)` after each step is a
    ///   plain Rust `if`.  The two branches assign the same set of
    ///   precomputed values in different orders, but the compiler is
    ///   free to emit a conditional jump.  This still closes the
    ///   coarse-grained "count operations to count set bits" leak that
    ///   double-and-add has.
    /// - The underlying field arithmetic uses `num-bigint`, which is
    ///   *not* constant-time at the limb level.  An attacker with
    ///   fine-grained timing access can still extract information from
    ///   individual `FieldElement::mul` / `inv` operations.  Fixing
    ///   that requires replacing `num-bigint` with a fixed-width,
    ///   constant-time bignum (out of scope for this change).
    ///
    /// In other words: this is a real improvement over double-and-add,
    /// but it does not make ECC operations side-channel safe.  See
    /// `SECURITY.md` for the full picture.
    pub fn scalar_mul_ct(&self, k: &BigUint, a: &FieldElement, scalar_bits: usize) -> Point {
        let mut r0 = Point::Infinity;
        let mut r1 = self.clone();

        for i in (0..scalar_bits).rev() {
            let bit = k.bit(i as u64);

            // Always compute all three potential next states so the
            // operation count is independent of `bit`.
            let sum = r0.add(&r1, a);
            let r0_dbl = r0.double(a);
            let r1_dbl = r1.double(a);

            // Re-bind (r0, r1) according to the current bit.  Both
            // arms assign exactly two values; the work is the same.
            if bit {
                r0 = sum;
                r1 = r1_dbl;
            } else {
                r0 = r0_dbl;
                r1 = sum;
            }
        }
        r0
    }

    /// Negate a point: -P = (x, -y mod p).
    pub fn neg(&self) -> Point {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => Point::Affine {
                x: x.clone(),
                y: y.neg(),
            },
        }
    }

    /// [`neg`](Self::neg) in **variable time**, for public points only:
    /// exactly `self.neg()`, through [`FieldElement::neg_vartime`].
    pub fn neg_vartime(&self) -> Point {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Affine { x, y } => Point::Affine {
                x: x.clone(),
                y: y.neg_vartime(),
            },
        }
    }

    /// Return the x-coordinate as a `BigUint`, or `None` for the point at infinity.
    pub fn x_coord(&self) -> Option<&BigUint> {
        match self {
            Point::Affine { x, .. } => Some(&x.value),
            Point::Infinity => None,
        }
    }
}

/// A finite point in Jacobian coordinates, `(X : Y : Z)` with `Z ≠ 0`
/// standing for the affine `(X/Z², Y/Z³)`; the identity is `None` where
/// these are used.  [`Point::scalar_mul`] and, above one word,
/// [`Point::linear_dlog_vartime`] use it.
struct Jacobian {
    x: FieldElement,
    y: FieldElement,
    z: FieldElement,
}

impl Jacobian {
    fn from_affine(x: &FieldElement, y: &FieldElement) -> Self {
        Jacobian {
            x: x.clone(),
            y: y.clone(),
            z: FieldElement::one(x.modulus.clone()),
        }
    }

    /// `2P` (dbl-2007-bl, general `a`).  `Y = 0` is exactly the affine
    /// `y = 0`, whose double is the identity.
    fn double(&self, a: &FieldElement) -> Option<Self> {
        if self.y.is_zero() {
            return None;
        }
        let xx = self.x.mul(&self.x);
        let yy = self.y.mul(&self.y);
        let yyyy = yy.mul(&yy);
        let zz = self.z.mul(&self.z);
        let t = self.x.add(&yy);
        let s = t.mul(&t).sub(&xx).sub(&yyyy);
        let s = s.add(&s);
        let m = xx.add(&xx).add(&xx).add(&a.mul(&zz.mul(&zz)));
        let x3 = m.mul(&m).sub(&s).sub(&s);
        let yyyy8 = {
            let d = yyyy.add(&yyyy);
            let q = d.add(&d);
            q.add(&q)
        };
        let y3 = m.mul(&s.sub(&x3)).sub(&yyyy8);
        let yz = self.y.add(&self.z);
        let z3 = yz.mul(&yz).sub(&yy).sub(&zz);
        Some(Jacobian {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    /// `P + (x₂, y₂)` for an affine `(x₂, y₂)` (madd-2007-bl).  Equal
    /// abscissae (`H = 0`) mean `Q = ±P`: double when the ordinates agree
    /// too, the identity otherwise — the cases [`Point::add`] branches on.
    fn add_affine(&self, x2: &FieldElement, y2: &FieldElement, a: &FieldElement) -> Option<Self> {
        let z1z1 = self.z.mul(&self.z);
        let u2 = x2.mul(&z1z1);
        let s2 = y2.mul(&self.z).mul(&z1z1);
        let h = u2.sub(&self.x);
        let r = s2.sub(&self.y);
        if h.is_zero() {
            return if r.is_zero() { self.double(a) } else { None };
        }
        let r = r.add(&r);
        let hh = h.mul(&h);
        let i = hh.add(&hh);
        let i = i.add(&i);
        let j = h.mul(&i);
        let v = self.x.mul(&i);
        let x3 = r.mul(&r).sub(&j).sub(&v).sub(&v);
        let y1j = self.y.mul(&j);
        let y3 = r.mul(&v.sub(&x3)).sub(&y1j).sub(&y1j);
        let zh = self.z.add(&h);
        let z3 = zh.mul(&zh).sub(&z1z1).sub(&hh);
        Some(Jacobian {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    /// Is this the affine `(x, y)`?  `X = x·Z²` and `Y = y·Z³`, compared
    /// without inverting `Z`; `x` and `y` must be reduced.
    fn is(&self, x: &FieldElement, y: &FieldElement) -> bool {
        let zz = self.z.mul(&self.z);
        self.x == x.mul(&zz) && self.y == y.mul(&zz.mul(&self.z))
    }

    fn to_affine(&self) -> Point {
        let zi = self.z.inv().expect("a Jacobian point has Z != 0");
        let zi2 = zi.mul(&zi);
        Point::Affine {
            x: self.x.mul(&zi2),
            y: self.y.mul(&zi2.mul(&zi)),
        }
    }
}

/// `F_p` on bare machine words, for a prime `p < 2⁶⁴`: the arithmetic of
/// the `*_vartime` operations when the modulus fits one word.
///
/// A [`FieldElement`] operation goes through `num-bigint`'s general
/// multi-limb code, reduces by a division even for an addition, and
/// clones the modulus into its result, so on a 32-bit curve an affine
/// addition spent far more in that machinery than in the arithmetic.
/// Here an element is its canonical value in a `u64`, a multiplication
/// is one 128-bit product and one remainder, an addition or subtraction
/// a compare, and an inversion [`inv_mod_u64`]; only results are turned
/// back into [`BigUint`].  Every value is exact mod `p`, and an inversion
/// is [`FieldElement::inv_vartime`]'s (see [`Word::inv`]), so the
/// formulas below take the same branches and produce the same points as
/// their [`FieldElement`] counterparts.  Variable time, like its callers.
#[derive(Clone, Copy)]
struct Word {
    p: u64,
}

impl Word {
    /// The word arithmetic mod `modulus`, if `modulus` fits one word.
    fn of(modulus: &BigUint) -> Option<Self> {
        single_word(modulus).map(|p| Word { p })
    }

    /// The elements' values as words, or `None` if one is not a reduced
    /// residue (never true of an element the arithmetic produced; the
    /// caller then keeps to the general path).
    fn get<const N: usize>(self, elems: [&FieldElement; N]) -> Option<[u64; N]> {
        let mut out = [0u64; N];
        for (o, e) in out.iter_mut().zip(elems) {
            *o = single_word(&e.value).filter(|&v| v < self.p)?;
        }
        Some(out)
    }

    /// An affine result (`None` for the identity) as a [`Point`] mod `m`.
    fn point(self, pt: Option<(u64, u64)>, m: &BigUint) -> Point {
        match pt {
            None => Point::Infinity,
            Some((x, y)) => Point::Affine {
                x: FieldElement {
                    value: BigUint::from(x),
                    modulus: m.clone(),
                },
                y: FieldElement {
                    value: BigUint::from(y),
                    modulus: m.clone(),
                },
            },
        }
    }

    fn add(self, a: u64, b: u64) -> u64 {
        let (s, carry) = a.overflowing_add(b);
        if carry || s >= self.p {
            s.wrapping_sub(self.p)
        } else {
            s
        }
    }

    fn sub(self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a.wrapping_sub(b).wrapping_add(self.p)
        }
    }

    fn mul(self, a: u64, b: u64) -> u64 {
        (u128::from(a) * u128::from(b) % u128::from(self.p)) as u64
    }

    /// `a⁻¹`, as [`FieldElement::inv_vartime`] computes it: Euclid, and
    /// [`non_unit_inv`](Self::non_unit_inv) where Euclid finds no inverse.
    fn inv(self, a: u64) -> u64 {
        match inv_mod_u64(a, self.p) {
            Some(i) => i,
            None => self.non_unit_inv(a),
        }
    }

    /// [`FieldElement::inv`]'s value `a^(p−2)` for an `a` with no inverse,
    /// which for nonzero `a` means a non-unit modulo a `p` that is not
    /// prime.  Zero has no value to take and panics, as `inv`'s `None`
    /// does at the `unwrap` in [`Point::add`] and [`Point::double`]; their
    /// callers rule it out the same way these do.  Out of line and cold:
    /// over a prime no walk reaches it, and it would otherwise weigh on
    /// the inlining of every formula that inverts.
    #[cold]
    #[inline(never)]
    fn non_unit_inv(self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero");
        self.pow(a, self.p - 2)
    }

    /// `a^e`, square and multiply; only for
    /// [`non_unit_inv`](Self::non_unit_inv).
    fn pow(self, a: u64, mut e: u64) -> u64 {
        let (mut base, mut acc) = (a, 1 % self.p);
        while e != 0 {
            if e & 1 == 1 {
                acc = self.mul(acc, base);
            }
            base = self.mul(base, base);
            e >>= 1;
        }
        acc
    }
}

/// Affine points in words, `None` for the identity: [`Point::add`] and
/// [`Point::double`] branch for branch.
impl Word {
    fn affine_add(self, (x1, y1): (u64, u64), (x2, y2): (u64, u64), a: u64) -> Option<(u64, u64)> {
        if x1 == x2 {
            return if y1 == y2 {
                self.affine_double((x1, y1), a)
            } else {
                None
            };
        }
        let lambda = self.mul(self.sub(y2, y1), self.inv(self.sub(x2, x1)));
        Some(self.chord(lambda, x1, x2, y1))
    }

    fn affine_double(self, (x, y): (u64, u64), a: u64) -> Option<(u64, u64)> {
        if y == 0 {
            return None;
        }
        let xx = self.mul(x, x);
        let numerator = self.add(self.add(self.add(xx, xx), xx), a);
        let lambda = self.mul(numerator, self.inv(self.add(y, y)));
        Some(self.chord(lambda, x, x, y))
    }

    /// The tail both formulas share: `x₃ = λ² − x₁ − x₂`,
    /// `y₃ = λ(x₁ − x₃) − y₁`.
    fn chord(self, lambda: u64, x1: u64, x2: u64, y1: u64) -> (u64, u64) {
        let x3 = self.sub(self.sub(self.mul(lambda, lambda), x1), x2);
        let y3 = self.sub(self.mul(lambda, self.sub(x1, x3)), y1);
        (x3, y3)
    }
}

/// A finite point `(X : Y : Z)` in Jacobian coordinates on [`Word`]s.
#[derive(Clone, Copy)]
struct WordJacobian {
    x: u64,
    y: u64,
    z: u64,
}

/// [`Jacobian`]'s formulas and branches on words, `None` for the
/// identity.
impl Word {
    /// dbl-2007-bl, as [`Jacobian::double`].
    fn jacobian_double(self, j: WordJacobian, a: u64) -> Option<WordJacobian> {
        if j.y == 0 {
            return None;
        }
        let xx = self.mul(j.x, j.x);
        let yy = self.mul(j.y, j.y);
        let yyyy = self.mul(yy, yy);
        let zz = self.mul(j.z, j.z);
        let t = self.add(j.x, yy);
        let s = self.sub(self.sub(self.mul(t, t), xx), yyyy);
        let s = self.add(s, s);
        let m = self.add(
            self.add(self.add(xx, xx), xx),
            self.mul(a, self.mul(zz, zz)),
        );
        let x3 = self.sub(self.sub(self.mul(m, m), s), s);
        let yyyy8 = {
            let d = self.add(yyyy, yyyy);
            let q = self.add(d, d);
            self.add(q, q)
        };
        let y3 = self.sub(self.mul(m, self.sub(s, x3)), yyyy8);
        let yz = self.add(j.y, j.z);
        let z3 = self.sub(self.sub(self.mul(yz, yz), yy), zz);
        Some(WordJacobian {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    /// madd-2007-bl, as [`Jacobian::add_affine`].
    fn jacobian_add_affine(
        self,
        j: WordJacobian,
        (x2, y2): (u64, u64),
        a: u64,
    ) -> Option<WordJacobian> {
        let z1z1 = self.mul(j.z, j.z);
        let u2 = self.mul(x2, z1z1);
        let s2 = self.mul(self.mul(y2, j.z), z1z1);
        let h = self.sub(u2, j.x);
        let r = self.sub(s2, j.y);
        if h == 0 {
            return if r == 0 {
                self.jacobian_double(j, a)
            } else {
                None
            };
        }
        let r = self.add(r, r);
        let hh = self.mul(h, h);
        let i = self.add(hh, hh);
        let i = self.add(i, i);
        let jj = self.mul(h, i);
        let v = self.mul(j.x, i);
        let x3 = self.sub(self.sub(self.sub(self.mul(r, r), jj), v), v);
        let y1j = self.mul(j.y, jj);
        let y3 = self.sub(self.sub(self.mul(r, self.sub(v, x3)), y1j), y1j);
        let zh = self.add(j.z, h);
        let z3 = self.sub(self.sub(self.mul(zh, zh), z1z1), hh);
        Some(WordJacobian {
            x: x3,
            y: y3,
            z: z3,
        })
    }

    /// As [`Jacobian::is`].
    fn jacobian_is(self, j: WordJacobian, (x, y): (u64, u64)) -> bool {
        let zz = self.mul(j.z, j.z);
        j.x == self.mul(x, zz) && j.y == self.mul(y, self.mul(zz, j.z))
    }

    fn jacobian_to_affine(self, j: WordJacobian) -> (u64, u64) {
        let zi = self.inv(j.z);
        let zi2 = self.mul(zi, zi);
        (self.mul(j.x, zi2), self.mul(j.y, self.mul(zi2, zi)))
    }
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    /// The affine double-and-add that `scalar_mul` replaced.
    fn scalar_mul_affine(p: &Point, k: &BigUint, a: &FieldElement) -> Point {
        let mut result = Point::Infinity;
        let mut addend = p.clone();
        let mut k = k.clone();
        let one = BigUint::one();
        while !k.is_zero() {
            if &k & &one == one {
                result = result.add(&addend, a);
            }
            addend = addend.double(a);
            k >>= 1;
        }
        result
    }

    /// Jacobian `scalar_mul` returns exactly the affine ladder's point,
    /// over every scalar up to past the group order of a small curve
    /// (identity, `2P = O` and `Q = ±P` cases included), and on P-256.
    #[test]
    fn jacobian_scalar_mul_matches_affine_double_and_add() {
        // y² = x³ + 2x + 3 over F_97: every point, every k < 2·#E + 3.
        let p = BigUint::from(97u32);
        let fe = |v: u32| FieldElement::new(BigUint::from(v), p.clone());
        let a = fe(2);
        let mut points = vec![Point::Infinity];
        for x in 0..97u32 {
            for y in 0..97u32 {
                if (y * y) % 97 == (x * x * x + 2 * x + 3) % 97 {
                    points.push(Point::Affine { x: fe(x), y: fe(y) });
                }
            }
        }
        let order = points.len() as u32;
        for pt in &points {
            for k in 0..(2 * order + 3) {
                let k = BigUint::from(k);
                assert_eq!(
                    pt.scalar_mul(&k, &a),
                    scalar_mul_affine(pt, &k, &a),
                    "{pt:?} k={k}"
                );
            }
        }
        let curve = crate::ecc::curve::CurveParams::p256();
        let g = curve.generator();
        let a = curve.a_fe();
        let mut k = BigUint::from(0x1234_5678_9abc_def1u64);
        for _ in 0..8 {
            k = (&k * &k + 12345u32) % &curve.n;
            assert_eq!(g.scalar_mul(&k, &a), scalar_mul_affine(&g, &k, &a));
        }
        let n_minus_1 = &curve.n - 1u32;
        assert_eq!(g.scalar_mul(&n_minus_1, &a), g.neg());
        assert_eq!(g.scalar_mul(&curve.n, &a), Point::Infinity);
    }
    use super::*;

    /// Every point of `y² = x³ + 2x + 3` over `F_97`, identity first.
    fn f97_points() -> (Vec<Point>, FieldElement) {
        let p = BigUint::from(97u32);
        let fe = |v: u32| FieldElement::new(BigUint::from(v), p.clone());
        let mut points = vec![Point::Infinity];
        for x in 0..97u32 {
            for y in 0..97u32 {
                if (y * y) % 97 == (x * x * x + 2 * x + 3) % 97 {
                    points.push(Point::Affine { x: fe(x), y: fe(y) });
                }
            }
        }
        (points, fe(2))
    }

    /// `(G, a)` on a curve over the largest 64-bit prime, `2⁶⁴ − 59`,
    /// with `b` chosen to put `G` on it: coordinates near `2⁶⁴` exercise
    /// the carries of the word arithmetic.
    fn word_edge_curve() -> (Point, FieldElement) {
        let p = BigUint::from(u64::MAX - 58);
        let fe = |v: u64| FieldElement::new(BigUint::from(v), p.clone());
        let g = Point::Affine {
            x: fe(u64::MAX - 1000),
            y: fe(u64::MAX - 77),
        };
        (g, fe(u64::MAX - 60))
    }

    /// The multiples `G, 2G, …, nG` by repeated [`Point::add`].
    fn multiples(g: &Point, a: &FieldElement, n: usize) -> Vec<Point> {
        let mut out = vec![g.clone()];
        while out.len() < n {
            let next = out.last().unwrap().add(g, a);
            out.push(next);
        }
        out
    }

    /// The vartime operations return exactly what the constant-time ones
    /// do: on every pair of points of a small curve (the word path, with
    /// its identity, doubling and `P + (−P)` cases), near the top of a
    /// 64-bit word, and on P-256 (the multi-word path).
    #[test]
    fn vartime_add_double_match_reference() {
        let (points, a) = f97_points();
        for p in &points {
            assert_eq!(p.double_vartime(&a), p.double(&a), "2·{p:?}");
            for q in &points {
                assert_eq!(p.add_vartime(q, &a), p.add(q, &a), "{p:?} + {q:?}");
            }
        }
        for p in &points {
            assert_eq!(p.neg_vartime(), p.neg());
        }
        let (g, a) = word_edge_curve();
        let ms = multiples(&g, &a, 40);
        for p in &ms {
            assert_eq!(p.double_vartime(&a), p.double(&a));
            for q in ms.iter().step_by(3) {
                assert_eq!(p.add_vartime(q, &a), p.add(q, &a));
                assert_eq!(p.add_vartime(&q.neg(), &a), p.add(&q.neg(), &a));
            }
        }
        let curve = crate::ecc::curve::CurveParams::p256();
        let a = curve.a_fe();
        let ms = multiples(&curve.generator(), &a, 6);
        for p in &ms {
            assert_eq!(p.double_vartime(&a), p.double(&a));
            for q in &ms {
                assert_eq!(p.add_vartime(q, &a), p.add(q, &a));
                assert_eq!(p.add_vartime(&q.neg(), &a), p.add(&q.neg(), &a));
            }
        }
    }

    /// `scalar_mul_vartime` is `scalar_mul`: every point and scalar up to
    /// past the group order of a small curve, scalars of every length
    /// near the top of a word, and the multi-word fallback on P-256.
    #[test]
    fn scalar_mul_vartime_matches_scalar_mul() {
        let (points, a) = f97_points();
        let order = points.len() as u32;
        for pt in &points {
            for k in 0..(2 * order + 3) {
                let k = BigUint::from(k);
                assert_eq!(pt.scalar_mul_vartime(&k, &a), pt.scalar_mul(&k, &a));
            }
        }
        let (g, a) = word_edge_curve();
        let mut k = BigUint::from(1u32);
        for _ in 0..70 {
            assert_eq!(g.scalar_mul_vartime(&k, &a), g.scalar_mul(&k, &a), "k={k}");
            k = &k * 3u32 + 1u32;
        }
        let curve = crate::ecc::curve::CurveParams::p256();
        let a = curve.a_fe();
        let k = BigUint::from(0xfeed_beef_dead_f00du64).pow(3);
        let g = curve.generator();
        assert_eq!(g.scalar_mul_vartime(&k, &a), g.scalar_mul(&k, &a));
        assert_eq!(g.scalar_mul_vartime(&curve.n, &a), Point::Infinity);
    }

    /// The linear search the Pohlig–Hellman digits used before, by affine
    /// additions.
    fn linear_dlog_affine(g: &Point, t: &Point, bound: u64, a: &FieldElement) -> Option<u64> {
        let mut current = Point::Infinity;
        for k in 0..bound {
            if current == *t {
                return Some(k);
            }
            current = current.add(g, a);
        }
        None
    }

    /// `linear_dlog_vartime` finds the `k` the affine search finds, or
    /// fails where it fails: every base/target pair of a small curve with
    /// bounds short of, at and past the order (a base of small order
    /// wraps through the identity), an unreduced target, and the
    /// multi-word path on P-256.
    #[test]
    fn linear_dlog_vartime_matches_affine_search() {
        let (points, a) = f97_points();
        let order = points.len() as u64;
        for g in &points {
            for t in &points {
                for bound in [0, 1, 3, order / 2, order + 2] {
                    assert_eq!(
                        g.linear_dlog_vartime(t, bound, &a),
                        linear_dlog_affine(g, t, bound, &a),
                        "{g:?} {t:?} {bound}"
                    );
                }
            }
        }
        if let Point::Affine { x, y } = &points[5] {
            let unreduced = Point::Affine {
                x: FieldElement {
                    value: &x.value + 97u32,
                    modulus: x.modulus.clone(),
                },
                y: y.clone(),
            };
            let g = &points[5];
            assert_eq!(g.linear_dlog_vartime(&unreduced, order + 2, &a), None);
            assert_eq!(linear_dlog_affine(g, &unreduced, order + 2, &a), None);
        }
        let (g, a) = word_edge_curve();
        let t = g.scalar_mul(&BigUint::from(29u32), &a);
        assert_eq!(g.linear_dlog_vartime(&t, 64, &a), Some(29));
        assert_eq!(g.linear_dlog_vartime(&t, 29, &a), None);
        let curve = crate::ecc::curve::CurveParams::p256();
        let a = curve.a_fe();
        let g = curve.generator();
        for k in [0u64, 1, 2, 11] {
            let t = g.scalar_mul(&BigUint::from(k), &a);
            assert_eq!(g.linear_dlog_vartime(&t, 12, &a), Some(k));
            assert_eq!(g.linear_dlog_vartime(&t, k, &a), None);
        }
    }

    /// Small curve for testing: y² = x³ + 2x + 3 mod 97
    /// Generator G = (3, 6), order 5.
    fn field(v: u64) -> FieldElement {
        FieldElement::new(BigUint::from(v), BigUint::from(97u64))
    }

    fn g() -> Point {
        Point::Affine {
            x: field(3),
            y: field(6),
        }
    }

    fn a() -> FieldElement {
        field(2)
    }

    #[test]
    fn double_then_add() {
        // 2G + G = 3G; verify the result is still a finite point
        let g = g();
        let g2 = g.double(&a());
        let g3 = g2.add(&g, &a());
        assert!(matches!(g3, Point::Affine { .. }));
    }

    #[test]
    fn scalar_mul_zero() {
        let g = g();
        let result = g.scalar_mul(&BigUint::zero(), &a());
        assert_eq!(result, Point::Infinity);
    }

    #[test]
    fn scalar_mul_one() {
        let g = g();
        let result = g.scalar_mul(&BigUint::one(), &a());
        assert_eq!(result, g);
    }

    #[test]
    fn ladder_matches_double_and_add() {
        let g = g();
        let a = a();
        // Try every scalar in [0, 5] (curve order) and verify ladder == double-and-add.
        for k in 0u64..=5 {
            let k_bu = BigUint::from(k);
            let v = g.scalar_mul(&k_bu, &a);
            let ct = g.scalar_mul_ct(&k_bu, &a, 4);
            assert_eq!(v, ct, "mismatch at k={}", k);
        }
    }

    #[test]
    fn ladder_matches_double_and_add_p256() {
        use super::super::curve::CurveParams;
        let curve = CurveParams::p256();
        let g = curve.generator();
        let a = curve.a_fe();
        // A handful of arbitrary scalars: the ladder must agree with double-and-add.
        for k_hex in [
            "1",
            "2",
            "deadbeef",
            "ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc63254f",
        ] {
            let k = BigUint::parse_bytes(k_hex.as_bytes(), 16).unwrap();
            let v = g.scalar_mul(&k, &a);
            let ct = g.scalar_mul_ct(&k, &a, 256);
            assert_eq!(v, ct, "mismatch at k={}", k_hex);
        }
    }

    #[test]
    fn neg_point() {
        let g = g();
        let neg_g = g.neg();
        // P + (-P) = ∞
        let sum = g.add(&neg_g, &a());
        assert_eq!(sum, Point::Infinity);
    }
}

/// Random differential tests of the `*_vartime` operations against the
/// operations they stand in for, over primes on both sides of every
/// width the word path distinguishes (8, 16, 31, 32 … 64 bits, and just
/// above one word, where the `BigUint` fallback takes over).  The
/// affine and Jacobian formulas are rational expressions in the
/// coordinates, so the pairs need not lie on a common curve; the
/// special cases (`x₁ = x₂`, `y₂ = ±y₁`, `y = 0`, values at `0`, `1`
/// and `p − 1`) are drawn deliberately.  Results are compared value
/// *and* modulus.
#[cfg(test)]
mod vartime_differential {
    use super::*;
    use num_bigint::RandBigInt;
    use num_traits::Zero;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    const PRIMES: [&str; 13] = [
        "251",
        "65521",
        "2147483647",
        "4294967291",
        "8589934583",
        "281474976710597",
        "4611686018427387847",
        "9223372036854775783",
        "9223372036854775837",
        "18446744073709551557",
        "18446744073709551629",
        "36893488147419103183",
        "340282366920938463463374607431768211297",
    ];

    fn primes() -> Vec<BigUint> {
        PRIMES
            .iter()
            .map(|s| BigUint::parse_bytes(s.as_bytes(), 10).unwrap())
            .collect()
    }

    fn same(a: &Point, b: &Point) -> bool {
        match (a, b) {
            (Point::Infinity, Point::Infinity) => true,
            (Point::Affine { x: x1, y: y1 }, Point::Affine { x: x2, y: y2 }) => {
                x1.value == x2.value
                    && y1.value == y2.value
                    && x1.modulus == x2.modulus
                    && y1.modulus == y2.modulus
            }
            _ => false,
        }
    }

    /// A coordinate: uniform, or one of the edge values.
    fn coord(rng: &mut StdRng, p: &BigUint) -> BigUint {
        match rng.gen_range(0..8) {
            0 => BigUint::zero(),
            1 => BigUint::from(1u32),
            2 => p - 1u32,
            3 => p - 2u32,
            _ => rng.gen_biguint_below(p),
        }
    }

    fn fe(v: BigUint, p: &BigUint) -> FieldElement {
        FieldElement::new(v, p.clone())
    }

    fn random_point(rng: &mut StdRng, p: &BigUint) -> Point {
        if rng.gen_range(0..16) == 0 {
            return Point::Infinity;
        }
        let y = if rng.gen_range(0..8) == 0 {
            BigUint::zero()
        } else {
            coord(rng, p)
        };
        Point::Affine {
            x: fe(coord(rng, p), p),
            y: fe(y, p),
        }
    }

    /// A second operand related to the first as the branches need it:
    /// the same point, its negative, the same `x` with an unrelated `y`,
    /// or unrelated.
    fn partner(rng: &mut StdRng, p: &BigUint, q: &Point) -> Point {
        match (rng.gen_range(0..6), q) {
            (0, _) => q.clone(),
            (1, _) => q.neg(),
            (2, Point::Affine { x, .. }) => Point::Affine {
                x: x.clone(),
                y: fe(coord(rng, p), p),
            },
            _ => random_point(rng, p),
        }
    }

    #[test]
    fn field_vartime_matches_constant_time() {
        let mut rng = StdRng::seed_from_u64(0x7e57_0001);
        for p in primes() {
            for _ in 0..400 {
                let a = fe(coord(&mut rng, &p), &p);
                let b = fe(coord(&mut rng, &p), &p);
                let (s, t) = (a.sub_vartime(&b), a.sub(&b));
                assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus), "{a} - {b}");
                let (s, t) = (a.neg_vartime(), a.neg());
                assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus), "-{a}");
                match (a.inv_vartime(), a.inv()) {
                    (None, None) => {}
                    (Some(s), Some(t)) => {
                        assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus), "1/{a}")
                    }
                    (s, t) => panic!("1/{a} mod {p}: {s:?} vs {t:?}"),
                }
            }
        }
    }

    #[test]
    fn point_vartime_matches_constant_time() {
        let mut rng = StdRng::seed_from_u64(0x7e57_0002);
        for p in primes() {
            for _ in 0..300 {
                let a = fe(coord(&mut rng, &p), &p);
                let p1 = random_point(&mut rng, &p);
                let p2 = partner(&mut rng, &p, &p1);
                assert!(
                    same(&p1.add_vartime(&p2, &a), &p1.add(&p2, &a)),
                    "{p1:?}+{p2:?}"
                );
                assert!(same(&p1.double_vartime(&a), &p1.double(&a)), "2·{p1:?}");
                assert!(same(&p1.neg_vartime(), &p1.neg()), "-{p1:?}");
                let k = match rng.gen_range(0..4) {
                    0 => BigUint::from(rng.gen_range(0..4u32)),
                    1 => BigUint::from(rng.gen::<u64>()),
                    _ => {
                        let bits = rng.gen_range(1..200);
                        rng.gen_biguint(bits)
                    }
                };
                assert!(
                    same(&p1.scalar_mul_vartime(&k, &a), &p1.scalar_mul(&k, &a)),
                    "{k}·{p1:?}"
                );
            }
        }
    }

    /// A base or target with an unreduced coordinate (the fields are
    /// public; no constructor makes one) is found where the affine search
    /// finds it: an unreduced base is its own first multiple, compared
    /// unreduced, and every later multiple is reduced.
    #[test]
    fn linear_dlog_vartime_unreduced_points_match_reference_search() {
        // (0, 10) on y² = x³ + 2x + 3 over F_97 has order 50, so no
        // multiple inside the bounds is ±g again (which would send an
        // unreduced chord to `x₂ − x₁ ≡ 0`, where both searches panic).
        let p = BigUint::from(97u32);
        let a = fe(BigUint::from(2u32), &p);
        let g = Point::Affine {
            x: fe(BigUint::zero(), &p),
            y: fe(BigUint::from(10u32), &p),
        };
        let lift = |pt: &Point| match pt {
            Point::Affine { x, y } => Point::Affine {
                x: FieldElement {
                    value: &x.value + &p,
                    modulus: p.clone(),
                },
                y: y.clone(),
            },
            Point::Infinity => Point::Infinity,
        };
        let g2 = g.add(&g, &a);
        let g3 = g2.add(&g, &a);
        for base in [g.clone(), lift(&g)] {
            for target in [
                Point::Infinity,
                g.clone(),
                lift(&g),
                g2.clone(),
                lift(&g2),
                g3.clone(),
            ] {
                for bound in [0, 1, 2, 3, 5] {
                    assert_eq!(
                        base.linear_dlog_vartime(&target, bound, &a),
                        linear_dlog_reference(&base, &target, bound, &a),
                        "base={base:?} target={target:?} bound={bound}"
                    );
                }
            }
        }
    }

    /// On a modulus that is not prime, outside the fields this type is
    /// for, a non-unit has no Euclidean inverse; the vartime operations
    /// then take `inv`'s value `a^(p−2)`, as the constant-time ones do,
    /// rather than panicking.  On one word (1009 · 2003) and on several
    /// ((2³¹ − 1)(2⁸⁹ − 1)): an inversion, a chord whose `x₂ − x₁` is a
    /// non-unit, a tangent whose `2y` is, and a ladder whose final `Z` is
    /// (1019 is the order of `(5, 7)` mod 1009 on `y² = x³ + 2x + b`).
    #[test]
    fn composite_modulus_non_units_match_constant_time() {
        let q1 = BigUint::from(2_147_483_647u32);
        let q2 = (BigUint::from(1u32) << 89) - 1u32;
        let cases = [
            (BigUint::from(1009u32 * 2003), BigUint::from(1009u32)),
            (&q1 * &q2, q1),
        ];
        for (n, q) in cases {
            let elem = |v: BigUint| fe(v, &n);
            let a = elem(BigUint::from(2u32));
            let z = elem(&q * 5u32);
            let (s, t) = (z.inv_vartime().unwrap(), z.inv().unwrap());
            assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus));
            let p1 = Point::Affine {
                x: elem(BigUint::from(7u32)),
                y: elem(BigUint::from(11u32)),
            };
            let p2 = Point::Affine {
                x: elem(&q * 3u32 + 7u32),
                y: elem(BigUint::from(13u32)),
            };
            assert!(same(&p1.add_vartime(&p2, &a), &p1.add(&p2, &a)), "{n}");
            let p3 = Point::Affine {
                x: elem(BigUint::from(5u32)),
                y: elem(&q * 2u32),
            };
            assert!(same(&p3.double_vartime(&a), &p3.double(&a)), "{n}");
            let g = Point::Affine {
                x: elem(BigUint::from(5u32)),
                y: elem(BigUint::from(7u32)),
            };
            let k = BigUint::from(1019u32);
            assert!(
                same(&g.scalar_mul_vartime(&k, &a), &g.scalar_mul(&k, &a)),
                "{n}"
            );
        }
    }

    /// The search `linear_dlog_vartime` replaced.
    fn linear_dlog_reference(g: &Point, t: &Point, bound: u64, a: &FieldElement) -> Option<u64> {
        let mut current = Point::Infinity;
        for k in 0..bound {
            if current == *t {
                return Some(k);
            }
            current = current.add(g, a);
        }
        None
    }

    /// Random curves through a random base, targets among its multiples
    /// (so small-order bases wrap through the identity on the small
    /// primes) and off them.
    #[test]
    fn linear_dlog_vartime_matches_reference_search() {
        let mut rng = StdRng::seed_from_u64(0x7e57_0003);
        for p in primes() {
            let small = p.bits() <= 16;
            let rounds = if small { 300 } else { 40 };
            for _ in 0..rounds {
                let a = fe(coord(&mut rng, &p), &p);
                let g = random_point(&mut rng, &p);
                let bound: u64 = if small {
                    rng.gen_range(0..600)
                } else {
                    rng.gen_range(0..60)
                };
                let t = match rng.gen_range(0..4) {
                    0 => random_point(&mut rng, &p),
                    1 => Point::Infinity,
                    _ => g.scalar_mul(&BigUint::from(rng.gen_range(0..=bound + 2)), &a),
                };
                assert_eq!(
                    g.linear_dlog_vartime(&t, bound, &a),
                    linear_dlog_reference(&g, &t, bound, &a),
                    "p={p} g={g:?} t={t:?} bound={bound}"
                );
            }
        }
    }
}
