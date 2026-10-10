//! The group interface the NIST ECDLP solvers run on, and its two
//! realisations: prime-field short Weierstrass curves (`P-xxx`) and binary
//! curves `y² + xy = x³ + ax² + b` over `F_{2^m}` (`B-xxx`, `K-xxx`).
//!
//! Everything a collision search needs is here: the group law, a stable
//! 64-bit key per point (for partitioning walks and for distinguished-point
//! tables), the negation map, and — on Koblitz curves — the Frobenius
//! endomorphism `τ(x, y) = (x², y²)` together with its eigenvalue `λ` on the
//! prime-order subgroup, so a walk can fold the `2m` points `±τ^i(P)` into
//! one representative and the coefficient bookkeeping stays exact.

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::binary_ecc::curve::{self as bcurve, BinaryCurve, BinaryPoint};
use crate::binary_ecc::f2m::F2mElement;
use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;

/// SplitMix64 finaliser: a cheap, well-mixed 64-bit hash step.
#[inline]
pub fn mix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Fold a little-endian word slice into one 64-bit key.
#[inline]
fn fold_words(seed: u64, words: impl IntoIterator<Item = u64>) -> u64 {
    let mut h = seed;
    for w in words {
        h = mix64(h ^ w);
    }
    h
}

/// A group of prime order `n` with the extra structure the solvers use.
///
/// Implementors must be `Sync`: the parallel rho shares one group across
/// worker threads.
pub trait EcdlpGroup: Sync {
    /// A group element.  Cloned on every walk step, so keep it compact.
    type Elt: Clone + PartialEq + Send + Sync + std::fmt::Debug;

    /// Display name of the curve.
    fn name(&self) -> &str;
    /// The (prime) order of the generator.
    fn order(&self) -> &BigUint;
    /// The generator `G`.
    fn generator(&self) -> Self::Elt;
    /// The point at infinity.
    fn identity(&self) -> Self::Elt;
    fn is_identity(&self, p: &Self::Elt) -> bool;
    fn add(&self, p: &Self::Elt, q: &Self::Elt) -> Self::Elt;
    fn double(&self, p: &Self::Elt) -> Self::Elt;
    fn neg(&self, p: &Self::Elt) -> Self::Elt;
    /// `[k]P` by double-and-add (variable time; every scalar here is public).
    fn mul(&self, p: &Self::Elt, k: &BigUint) -> Self::Elt {
        let mut acc = self.identity();
        if k.is_zero() || self.is_identity(p) {
            return acc;
        }
        let bits = k.bits();
        for i in (0..bits).rev() {
            acc = self.double(&acc);
            if k.bit(i) {
                acc = self.add(&acc, p);
            }
        }
        acc
    }
    /// A 64-bit key of the point, equal for equal points.  Used to choose
    /// the walk's jump and to decide distinguishedness, so it must depend on
    /// both coordinates.
    fn key(&self, p: &Self::Elt) -> u64;
    /// The negation-map representative of `{P, −P}` and whether it is `−P`.
    fn canonical_sign(&self, p: &Self::Elt) -> (Self::Elt, bool);
    /// Size of the cyclic automorphism group folded by
    /// [`canonical_orbit`](Self::canonical_orbit) beyond negation: `m` on a
    /// Koblitz curve over `F_{2^m}` (Frobenius), `1` otherwise.
    fn automorphism_order(&self) -> u32 {
        1
    }
    /// `λ^i mod n` for `i < automorphism_order()`, where the automorphism
    /// acts as `[λ]` on the subgroup.  `[1]` when there is none.
    fn orbit_multipliers(&self) -> Vec<BigUint> {
        vec![BigUint::one()]
    }
    /// The representative of the orbit `{±σ^i(P)}` together with `(i, neg)`
    /// such that the representative equals `[(−1)^neg · λ^i] P`.
    fn canonical_orbit(&self, p: &Self::Elt) -> (Self::Elt, u32, bool) {
        let (r, neg) = self.canonical_sign(p);
        (r, 0, neg)
    }
    /// Encode a point for display (hex coordinates) — `None` for infinity.
    fn coords_hex(&self, p: &Self::Elt) -> Option<(String, String)>;
}

// ── Prime field curves ──────────────────────────────────────────────────────

/// A prime-field curve `y² = x³ + ax + b` over `F_p` with a generator of
/// prime order `n`.  Thin wrapper over [`CurveParams`] with the curve
/// coefficient `a` pre-lifted into the field.
#[derive(Clone)]
pub struct PrimeGroup {
    pub curve: CurveParams,
    a: FieldElement,
    half_p: BigUint,
}

impl PrimeGroup {
    pub fn new(curve: CurveParams) -> Self {
        let a = curve.a_fe();
        let half_p = &curve.p >> 1u32;
        Self { curve, a, half_p }
    }

    /// The field coefficient `a` as a [`FieldElement`].
    pub fn a_fe(&self) -> &FieldElement {
        &self.a
    }

    /// Build the affine point `(x, y)` without checking it is on the curve.
    pub fn point(&self, x: &BigUint, y: &BigUint) -> Point {
        Point::Affine {
            x: self.curve.fe(x.clone()),
            y: self.curve.fe(y.clone()),
        }
    }

    /// Whether `p` is on the curve (the identity counts).
    pub fn is_on_curve(&self, p: &Point) -> bool {
        self.curve.is_on_curve(p)
    }
}

impl EcdlpGroup for PrimeGroup {
    type Elt = Point;

    fn name(&self) -> &str {
        self.curve.name
    }
    fn order(&self) -> &BigUint {
        &self.curve.n
    }
    fn generator(&self) -> Point {
        self.curve.generator()
    }
    fn identity(&self) -> Point {
        Point::Infinity
    }
    fn is_identity(&self, p: &Point) -> bool {
        matches!(p, Point::Infinity)
    }
    fn add(&self, p: &Point, q: &Point) -> Point {
        p.add_vartime(q, &self.a)
    }
    fn double(&self, p: &Point) -> Point {
        p.double_vartime(&self.a)
    }
    fn neg(&self, p: &Point) -> Point {
        p.neg_vartime()
    }
    fn mul(&self, p: &Point, k: &BigUint) -> Point {
        p.scalar_mul_vartime(k, &self.a)
    }
    fn key(&self, p: &Point) -> u64 {
        match p {
            Point::Infinity => 0x5EED_0000_0000_0001,
            Point::Affine { x, y } => {
                let h = fold_words(0x1234_5678_9ABC_DEF0, x.value.iter_u64_digits());
                fold_words(h, y.value.iter_u64_digits())
            }
        }
    }
    fn canonical_sign(&self, p: &Point) -> (Point, bool) {
        match p {
            Point::Infinity => (Point::Infinity, false),
            Point::Affine { y, .. } => {
                if y.value > self.half_p {
                    (p.neg_vartime(), true)
                } else {
                    (p.clone(), false)
                }
            }
        }
    }
    fn coords_hex(&self, p: &Point) -> Option<(String, String)> {
        match p {
            Point::Infinity => None,
            Point::Affine { x, y } => Some((format!("{:x}", x.value), format!("{:x}", y.value))),
        }
    }
}

// ── Binary curves (random and Koblitz) ──────────────────────────────────────

/// Frobenius data for a Koblitz curve `E_a : y² + xy = x³ + ax² + 1` over
/// `F_{2^m}`: `τ(x, y) = (x², y²)` satisfies `τ² − μτ + 2 = 0` with
/// `μ = (−1)^{1−a}`, and acts on the prime-order subgroup as `[λ]` where
/// `λ` is the matching root of `λ² − μλ + 2 ≡ 0 (mod n)`.
#[derive(Clone, Debug)]
pub struct Frobenius {
    pub m: u32,
    pub mu: i8,
    pub lambda: BigUint,
    /// `λ^i mod n` for `0 ≤ i < m`.
    pub powers: Vec<BigUint>,
}

/// A binary curve with its generator's prime-order subgroup; carries
/// Frobenius data when the curve is Koblitz (`b = 1`, `a ∈ {0, 1}`).
#[derive(Clone)]
pub struct BinaryGroup {
    pub curve: BinaryCurve,
    pub name: String,
    pub frobenius: Option<Frobenius>,
}

impl BinaryGroup {
    /// Wrap a binary curve, detecting Koblitz structure from `(a, b)`.
    ///
    /// Returns `Err` only when the curve is Koblitz but no root of
    /// `λ² − μλ + 2` mod `n` matches `τ(G)` — which would mean `n` is not
    /// the generator's order or the curve constants are inconsistent.
    pub fn new(name: &str, curve: BinaryCurve) -> Result<Self, String> {
        let m = curve.m;
        let is_koblitz =
            curve.b == F2mElement::one(m) && (curve.a.is_zero() || curve.a == F2mElement::one(m));
        let frobenius = if is_koblitz {
            Some(frobenius_data(&curve)?)
        } else {
            None
        };
        Ok(Self {
            curve,
            name: name.to_string(),
            frobenius,
        })
    }

    /// Wrap a binary curve *without* Frobenius folding even if it is Koblitz
    /// (the control arm for measuring the `√m` gain).
    pub fn without_frobenius(name: &str, curve: BinaryCurve) -> Self {
        Self {
            curve,
            name: name.to_string(),
            frobenius: None,
        }
    }

    /// `τ(P) = (x², y²)`.
    pub fn frobenius_map(&self, p: &BinaryPoint) -> BinaryPoint {
        match p {
            BinaryPoint::Infinity => BinaryPoint::Infinity,
            BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
                x: x.square(&self.curve.irreducible),
                y: y.square(&self.curve.irreducible),
            },
        }
    }

    /// Build the affine point `(x, y)` from integer-encoded coordinates.
    pub fn point(&self, x: &BigUint, y: &BigUint) -> BinaryPoint {
        BinaryPoint::Affine {
            x: F2mElement::from_biguint(x, self.curve.m),
            y: F2mElement::from_biguint(y, self.curve.m),
        }
    }

    pub fn is_on_curve(&self, p: &BinaryPoint) -> bool {
        self.curve.is_on_curve(p)
    }

    /// Solve `y² + xy = x³ + ax² + b` for `y` given `x` (odd `m` only, via
    /// the half-trace), returning one of the two points or `None`.
    pub fn lift_x(&self, x: &F2mElement) -> Option<BinaryPoint> {
        let irr = &self.curve.irreducible;
        let m = self.curve.m;
        if x.is_zero() {
            // y² = b  ⇒  y = b^(2^(m−1)).
            let y = self.curve.b.square_k_times(m - 1, irr);
            return Some(BinaryPoint::Affine { x: x.clone(), y });
        }
        if m.is_multiple_of(2) {
            return None;
        }
        // Substitute y = x·z: z² + z = x + a + b/x².
        let x_inv = x.flt_inverse(irr)?;
        let c = x
            .add(&self.curve.a)
            .add(&self.curve.b.mul(&x_inv.square(irr), irr));
        // Half-trace H(c) = Σ_{i=0}^{(m−1)/2} c^(2^(2i)) solves z² + z = c
        // when Tr(c) = 0.
        let mut z = c.clone();
        let mut term = c.clone();
        for _ in 0..(m - 1) / 2 {
            term = term.square_k_times(2, irr);
            z = z.add(&term);
        }
        // Verify (Tr(c) might be 1).
        if z.square(irr).add(&z) != c {
            return None;
        }
        let y = x.mul(&z, irr);
        Some(BinaryPoint::Affine { x: x.clone(), y })
    }
}

/// Compute the Frobenius eigenvalue on the generator's subgroup.
fn frobenius_data(curve: &BinaryCurve) -> Result<Frobenius, String> {
    let m = curve.m;
    let n = &curve.order;
    let mu: i8 = if curve.a.is_zero() { -1 } else { 1 };
    // λ = (μ ± √(μ² − 8)) / 2 = (μ ± √(−7)) / 2  (mod n).
    let minus7 = (n - BigUint::from(7u32) % n) % n;
    let root = sqrt_mod_p(&minus7, n).ok_or_else(|| {
        format!(
            "−7 is not a square mod n = {n}; the subgroup order of a Koblitz curve must split in Z[τ]"
        )
    })?;
    let inv2 = BigUint::from(2u32).modpow(&(n - BigUint::from(2u32)), n);
    let mu_mod = if mu == 1 {
        BigUint::one()
    } else {
        n - BigUint::one()
    };
    let cand1 = (&mu_mod + &root) % n * &inv2 % n;
    let cand2 = ((&mu_mod + n - &root) % n) * &inv2 % n;
    let g = &curve.generator;
    let tau_g = match g {
        BinaryPoint::Infinity => return Err("generator is the identity".into()),
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.square(&curve.irreducible),
            y: y.square(&curve.irreducible),
        },
    };
    let lambda = [cand1, cand2]
        .into_iter()
        .find(|l| bcurve::scalar_mul(curve, g, l) == tau_g)
        .ok_or_else(|| "neither root of λ² − μλ + 2 matches τ(G)".to_string())?;
    let mut powers = Vec::with_capacity(m as usize);
    let mut acc = BigUint::one();
    for _ in 0..m {
        powers.push(acc.clone());
        acc = acc * &lambda % n;
    }
    Ok(Frobenius {
        m,
        mu,
        lambda,
        powers,
    })
}

fn f2m_key(seed: u64, e: &F2mElement) -> u64 {
    fold_words(seed, e.raw_bits().iter().copied())
}

/// `a < b` as integers (bit-vectors compared from the top word).
fn f2m_less(a: &F2mElement, b: &F2mElement) -> bool {
    let (wa, wb) = (a.raw_bits(), b.raw_bits());
    let len = wa.len().max(wb.len());
    for i in (0..len).rev() {
        let (x, y) = (
            wa.get(i).copied().unwrap_or(0),
            wb.get(i).copied().unwrap_or(0),
        );
        if x != y {
            return x < y;
        }
    }
    false
}

impl EcdlpGroup for BinaryGroup {
    type Elt = BinaryPoint;

    fn name(&self) -> &str {
        &self.name
    }
    fn order(&self) -> &BigUint {
        &self.curve.order
    }
    fn generator(&self) -> BinaryPoint {
        self.curve.generator.clone()
    }
    fn identity(&self) -> BinaryPoint {
        BinaryPoint::Infinity
    }
    fn is_identity(&self, p: &BinaryPoint) -> bool {
        matches!(p, BinaryPoint::Infinity)
    }
    fn add(&self, p: &BinaryPoint, q: &BinaryPoint) -> BinaryPoint {
        bcurve::point_add(&self.curve, p, q)
    }
    fn double(&self, p: &BinaryPoint) -> BinaryPoint {
        bcurve::point_double(&self.curve, p)
    }
    fn neg(&self, p: &BinaryPoint) -> BinaryPoint {
        bcurve::point_neg(p)
    }
    fn mul(&self, p: &BinaryPoint, k: &BigUint) -> BinaryPoint {
        bcurve::scalar_mul(&self.curve, p, k)
    }
    fn key(&self, p: &BinaryPoint) -> u64 {
        match p {
            BinaryPoint::Infinity => 0x5EED_0000_0000_0002,
            BinaryPoint::Affine { x, y } => {
                let h = f2m_key(0x0F1E_2D3C_4B5A_6978, x);
                f2m_key(h, y)
            }
        }
    }
    fn canonical_sign(&self, p: &BinaryPoint) -> (BinaryPoint, bool) {
        match p {
            BinaryPoint::Infinity => (BinaryPoint::Infinity, false),
            BinaryPoint::Affine { x, y } => {
                // −P = (x, x + y); pick the smaller y as integer.
                let y_neg = y.add(x);
                if f2m_less(&y_neg, y) {
                    (
                        BinaryPoint::Affine {
                            x: x.clone(),
                            y: y_neg,
                        },
                        true,
                    )
                } else {
                    (p.clone(), false)
                }
            }
        }
    }
    fn automorphism_order(&self) -> u32 {
        self.frobenius.as_ref().map_or(1, |f| f.m)
    }
    fn orbit_multipliers(&self) -> Vec<BigUint> {
        self.frobenius
            .as_ref()
            .map_or_else(|| vec![BigUint::one()], |f| f.powers.clone())
    }
    fn canonical_orbit(&self, p: &BinaryPoint) -> (BinaryPoint, u32, bool) {
        let Some(fr) = &self.frobenius else {
            let (r, neg) = self.canonical_sign(p);
            return (r, 0, neg);
        };
        let BinaryPoint::Affine { x, y } = p else {
            return (BinaryPoint::Infinity, 0, false);
        };
        let irr = &self.curve.irreducible;
        // Walk the Frobenius orbit of x; keep the minimal x and its index.
        let mut best_x = x.clone();
        let mut best_i = 0u32;
        let mut cur = x.clone();
        for i in 1..fr.m {
            cur = cur.square(irr);
            if f2m_less(&cur, &best_x) {
                best_x = cur.clone();
                best_i = i;
            }
        }
        let y_i = y.square_k_times(best_i, irr);
        let rep = BinaryPoint::Affine { x: best_x, y: y_i };
        let (rep, neg) = self.canonical_sign(&rep);
        (rep, best_i, neg)
    }
    fn coords_hex(&self, p: &BinaryPoint) -> Option<(String, String)> {
        match p {
            BinaryPoint::Infinity => None,
            BinaryPoint::Affine { x, y } => Some((
                format!("{:x}", x.to_biguint()),
                format!("{:x}", y.to_biguint()),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecdlp_nist::toy;

    #[test]
    fn koblitz_frobenius_eigenvalue_matches_on_random_multiples() {
        let g = toy::koblitz(17, 1).expect("toy K curve");
        let fr = g.frobenius.as_ref().expect("Koblitz");
        assert_eq!(fr.m, 17);
        let p = g.mul(&g.generator(), &BigUint::from(12345u32));
        assert_eq!(g.frobenius_map(&p), g.mul(&p, &fr.lambda));
        // λ^m ≡ 1 on the subgroup since τ^m = identity on E(F_{2^m}).
        let lm = fr.lambda.modpow(&BigUint::from(fr.m), g.order());
        assert_eq!(lm, BigUint::one());
    }

    #[test]
    fn canonical_orbit_is_a_class_function_with_exact_multiplier() {
        let g = toy::koblitz(19, 0).expect("toy K curve");
        let p = g.mul(&g.generator(), &BigUint::from(777u32));
        let (rep, i, neg) = g.canonical_orbit(&p);
        let mults = g.orbit_multipliers();
        let mut s = mults[i as usize].clone();
        if neg {
            s = g.order() - s;
        }
        assert_eq!(rep, g.mul(&p, &s), "rep = [±λ^i]P");
        // Every element of the orbit maps to the same representative.
        let mut q = p.clone();
        for _ in 0..g.automorphism_order() {
            q = g.frobenius_map(&q);
            assert_eq!(g.canonical_orbit(&q).0, rep);
            assert_eq!(g.canonical_orbit(&g.neg(&q)).0, rep);
        }
    }

    #[test]
    fn prime_negation_map_is_a_class_function() {
        let g = toy::prime_small();
        let p = g.mul(&g.generator(), &BigUint::from(4242u32));
        let (r1, s1) = g.canonical_sign(&p);
        let (r2, s2) = g.canonical_sign(&g.neg(&p));
        assert_eq!(r1, r2);
        assert_ne!(s1, s2);
    }

    #[test]
    fn binary_lift_x_lands_on_curve() {
        let g = toy::koblitz(23, 1).expect("toy");
        let gen = g.generator();
        let BinaryPoint::Affine { x, .. } = &gen else {
            panic!()
        };
        let p = g.lift_x(x).expect("x of a point lifts");
        assert!(g.is_on_curve(&p));
        assert!(p == gen || p == g.neg(&gen));
    }
}
