//! Curve-generic, variable-time Jacobian arithmetic on a short Weierstrass
//! curve `y² = x³ + a·x + b` over [`GenericField`], and the
//! scalar-multiplication arms of the native GLV-chain benchmark
//! (`research/endo_native_bench_20261016`): width-`w` NAF, and GLV-2 with an
//! isogeny-chain endomorphism.
//!
//! **This is a measurement harness, not a production implementation.**
//! Every routine branches on point and scalar values; never use it with a
//! secret scalar.
//!
//! It is the curve-generic form of [`crate::ecc::cryptopro_b_point`] (whose
//! code paths it leaves untouched) with the formulas the autoresearcher's
//! chain-sweep model counts for each curve (`harness/endosweep/curvesweep.py`,
//! `arithmetic`):
//!
//! * doubling — `a = 0`: EFD `dbl-2009-l` (2M + 5S); `a = −3`: `dbl-2001-b`
//!   (3M + 5S); `|a| ≤ 16`: `dbl-2007-bl` with `a·ZZ²` by additions
//!   (1M + 8S); any other `a`: `dbl-2007-bl` with one multiplication by `a`
//!   (2M + 8S).  A curve with `a ≠ −3` that is isomorphic to an `a = −3`
//!   model is refused rather than silently counted with the wrong formula;
//! * addition — `add-2007-bl` (11M + 5S) and `madd-2007-bl` (7M + 4S) with an
//!   affine second operand, for every curve (they do not use `a`);
//! * inversion — Fermat, `bits(p−2)` squarings and `popcount(p−2)`
//!   multiplications.
//!
//! The chain evaluators are those of the CryptoPro-B module, generalised to
//! the curve's own field and coefficients, plus degree-2 steps (Vélu, for
//! curves where 2 splits or ramifies, e.g. Bandersnatch's `√−2`):
//!
//! * [`ChainEval::Optimised`] — Jacobian steps; an odd `ℓ`-step costs
//!   `3ℓ − 4` M on the affine input and `11s + 2ℓ − 2` M + 1 S later
//!   (`s = (ℓ−1)/2`); a 2-step costs 4M + 1S on the affine input and
//!   7M + 3S later (`X' = (X·d + v·Z⁴)·d`, `Y' = Y·(d² − v·Z⁴)·d`,
//!   `Z' = Z·d`, `d = X − x₀·Z²`); the isomorphism back to `E` is `Z·u⁻¹`.
//! * [`ChainEval::Generic`] — the evaluator of crypto#1408: projective `P²`
//!   steps, `11s + 2ℓ + 3` M + 1 S per odd step, then `u²X, u³Y` and the
//!   projective-to-Jacobian conversion (4M + 1S).  A 2-step in projective
//!   coordinates is `X' = (X·d + v·Z²)·d`, `Y' = Y·(d² − v·Z²)`,
//!   `Z' = Z·d²` with `d = X − x₀·Z`: 3M + 1S on the affine input, 6M + 2S
//!   later.

use super::endo_generic_field::{Fe, GenericField};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};

/// Which doubling formula the curve's `a` admits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DoublingKind {
    /// `a = 0`: dbl-2009-l, 2M + 5S.
    AZero,
    /// `a = −3`: dbl-2001-b, 3M + 5S.
    AMinus3,
    /// `|a| ≤ 16`, `a ∉ {0, −3}`: dbl-2007-bl with `a·ZZ²` by additions,
    /// 1M + 8S.
    SmallA(i64),
    /// Any other `a`: dbl-2007-bl with one multiplication by `a`, 2M + 8S.
    General,
}

impl DoublingKind {
    /// `(M, S)` of one doubling.
    pub fn cost(&self) -> (u64, u64) {
        match self {
            DoublingKind::AZero => (2, 5),
            DoublingKind::AMinus3 => (3, 5),
            DoublingKind::SmallA(_) => (1, 8),
            DoublingKind::General => (2, 8),
        }
    }

    pub fn describe(&self) -> String {
        match self {
            DoublingKind::AZero => "a = 0: dbl-2009-l 2M + 5S".into(),
            DoublingKind::AMinus3 => "a = -3: dbl-2001-b 3M + 5S".into(),
            DoublingKind::SmallA(a) => format!("small a = {a}: dbl-2007-bl 1M + 8S"),
            DoublingKind::General => "general a: dbl-2007-bl 1M + 8S + 1M for a".into(),
        }
    }
}

/// `x` is a non-zero fourth power mod the prime `p`.
fn is_fourth_power(x: &BigUint, p: &BigUint) -> bool {
    let x = x % p;
    if x.is_zero() {
        return false;
    }
    let pm1 = p - 1u8;
    let g = if (&pm1 % 4u8).is_zero() { 4u8 } else { 2u8 };
    x.modpow(&(&pm1 / g), p).is_one()
}

/// The doubling formula for `y² = x³ + a·x + b` over `F_p`, with the same
/// case analysis as the model (`curvesweep.arithmetic`).  `Err` for a curve
/// whose `a ≠ −3` is isomorphic to an `a = −3` model: the model counts it
/// with dbl-2001-b on that model, which this harness does not build.
pub fn doubling_kind(p: &BigUint, a: &BigUint) -> Result<DoublingKind, String> {
    let a = a % p;
    if a.is_zero() {
        return Ok(DoublingKind::AZero);
    }
    if a == p - 3u8 {
        return Ok(DoublingKind::AMinus3);
    }
    let a_inv = a.modpow(&(p - 2u8), p);
    let t = ((p - 3u8) * a_inv) % p;
    if is_fourth_power(&t, p) {
        return Err(
            "a is isomorphic to an a = -3 model (-3/a is a fourth power); not supported".into(),
        );
    }
    let neg = p - &a;
    let small = BigUint::from(16u8);
    if a <= small {
        return Ok(DoublingKind::SmallA(
            i64::try_from(a.to_u64_digits().first().copied().unwrap_or(0)).unwrap(),
        ));
    }
    if neg <= small {
        return Ok(DoublingKind::SmallA(
            -i64::try_from(neg.to_u64_digits().first().copied().unwrap_or(0)).unwrap(),
        ));
    }
    Ok(DoublingKind::General)
}

/// An affine point; never the identity.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Affine<const L: usize> {
    pub x: Fe<L>,
    pub y: Fe<L>,
}

/// A Jacobian point `(X : Y : Z)`, `x = X/Z²`, `y = Y/Z³`; the identity is
/// any point with `Z = 0`.
#[derive(Copy, Clone, Debug)]
pub struct Jacobian<const L: usize> {
    pub x: Fe<L>,
    pub y: Fe<L>,
    pub z: Fe<L>,
}

impl<const L: usize> Jacobian<L> {
    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }
}

/// The curve: its field, coefficients and doubling formula.
#[derive(Clone, Debug)]
pub struct Curve<const L: usize> {
    pub f: GenericField<L>,
    pub a: Fe<L>,
    pub b: Fe<L>,
    pub a_big: BigUint,
    pub b_big: BigUint,
    pub dbl: DoublingKind,
}

impl<const L: usize> Curve<L> {
    pub fn new(p: &BigUint, a: &BigUint, b: &BigUint) -> Result<Self, String> {
        let f = GenericField::<L>::new(p)?;
        if a >= p || b >= p {
            return Err("curve coefficients must be canonical (below p)".into());
        }
        let disc = (BigUint::from(4u8) * a.modpow(&BigUint::from(3u8), p)
            + BigUint::from(27u8) * b * b)
            % p;
        if disc.is_zero() {
            return Err("the curve is singular".into());
        }
        let dbl = doubling_kind(p, a)?;
        Ok(Curve {
            a: f.from_biguint(a),
            b: f.from_biguint(b),
            a_big: a.clone(),
            b_big: b.clone(),
            f,
            dbl,
        })
    }

    pub fn identity(&self) -> Jacobian<L> {
        Jacobian {
            x: self.f.one(),
            y: self.f.one(),
            z: self.f.zero(),
        }
    }

    pub fn affine(&self, x: &BigUint, y: &BigUint) -> Result<Affine<L>, String> {
        Ok(Affine {
            x: self.f.from_canonical(x)?,
            y: self.f.from_canonical(y)?,
        })
    }

    pub fn affine_to_biguint(&self, p: &Affine<L>) -> (BigUint, BigUint) {
        (self.f.to_biguint(&p.x), self.f.to_biguint(&p.y))
    }

    /// `y² = x³ + a·x + b` for the given coefficients.
    pub fn on_curve_with(&self, p: &Affine<L>, a: &Fe<L>, b: &Fe<L>) -> bool {
        let f = &self.f;
        let lhs = f.sqr(&p.y);
        let rhs = f.add(&f.add(&f.mul(&f.sqr(&p.x), &p.x), &f.mul(a, &p.x)), b);
        lhs == rhs
    }

    pub fn is_on_curve(&self, p: &Affine<L>) -> bool {
        self.on_curve_with(p, &self.a, &self.b)
    }

    pub fn neg_affine(&self, p: &Affine<L>) -> Affine<L> {
        Affine {
            x: p.x,
            y: self.f.neg(&p.y),
        }
    }

    pub fn neg(&self, p: &Jacobian<L>) -> Jacobian<L> {
        Jacobian {
            x: p.x,
            y: self.f.neg(&p.y),
            z: p.z,
        }
    }

    pub fn from_affine(&self, p: &Affine<L>) -> Jacobian<L> {
        Jacobian {
            x: p.x,
            y: p.y,
            z: self.f.one(),
        }
    }

    /// One inversion; `None` for the identity.
    pub fn to_affine(&self, p: &Jacobian<L>) -> Option<Affine<L>> {
        if p.is_identity() {
            return None;
        }
        let f = &self.f;
        let zi = f.inv(&p.z);
        let zi2 = f.sqr(&zi);
        let zi3 = f.mul(&zi2, &zi);
        Some(Affine {
            x: f.mul(&p.x, &zi2),
            y: f.mul(&p.y, &zi3),
        })
    }

    /// Equality as curve points by cross-multiplication.
    pub fn eq_point(&self, p: &Jacobian<L>, q: &Jacobian<L>) -> bool {
        match (p.is_identity(), q.is_identity()) {
            (true, true) => true,
            (true, false) | (false, true) => false,
            (false, false) => {
                let f = &self.f;
                let z1z1 = f.sqr(&p.z);
                let z2z2 = f.sqr(&q.z);
                let u1 = f.mul(&p.x, &z2z2);
                let u2 = f.mul(&q.x, &z1z1);
                let s1 = f.mul(&f.mul(&p.y, &q.z), &z2z2);
                let s2 = f.mul(&f.mul(&q.y, &p.z), &z1z1);
                u1 == u2 && s1 == s2
            }
        }
    }

    /// Doubling with the curve's formula ([`DoublingKind`]).  `Y = 0` (order
    /// two) gives `Z3 = 0`; the harness only doubles points of odd prime
    /// order and their multiples.
    pub fn double(&self, p: &Jacobian<L>) -> Jacobian<L> {
        if p.is_identity() {
            return *p;
        }
        let f = &self.f;
        match self.dbl {
            DoublingKind::AMinus3 => {
                // dbl-2001-b, exactly the operation sequence of
                // CryptoProBJacobian::double
                let delta = f.sqr(&p.z);
                let gamma = f.sqr(&p.y);
                let beta = f.mul(&p.x, &gamma);
                let t = f.mul(&f.sub(&p.x, &delta), &f.add(&p.x, &delta));
                let alpha = f.add(&f.add(&t, &t), &t);
                let beta2 = f.add(&beta, &beta);
                let beta4 = f.add(&beta2, &beta2);
                let beta8 = f.add(&beta4, &beta4);
                let x3 = f.sub(&f.sqr(&alpha), &beta8);
                let z3 = f.sub(&f.sub(&f.sqr(&f.add(&p.y, &p.z)), &gamma), &delta);
                let gamma2 = f.sqr(&gamma);
                let g2 = f.add(&gamma2, &gamma2);
                let g4 = f.add(&g2, &g2);
                let g8 = f.add(&g4, &g4);
                let y3 = f.sub(&f.mul(&alpha, &f.sub(&beta4, &x3)), &g8);
                Jacobian {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }
            DoublingKind::AZero => {
                // dbl-2009-l
                let a = f.sqr(&p.x);
                let b = f.sqr(&p.y);
                let c = f.sqr(&b);
                let t = f.sub(&f.sub(&f.sqr(&f.add(&p.x, &b)), &a), &c);
                let d = f.add(&t, &t);
                let e = f.add(&f.add(&a, &a), &a);
                let ff = f.sqr(&e);
                let x3 = f.sub(&ff, &f.add(&d, &d));
                let c2 = f.add(&c, &c);
                let c4 = f.add(&c2, &c2);
                let c8 = f.add(&c4, &c4);
                let y3 = f.sub(&f.mul(&e, &f.sub(&d, &x3)), &c8);
                let yz = f.mul(&p.y, &p.z);
                let z3 = f.add(&yz, &yz);
                Jacobian {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }
            DoublingKind::SmallA(_) | DoublingKind::General => {
                // dbl-2007-bl
                let xx = f.sqr(&p.x);
                let yy = f.sqr(&p.y);
                let yyyy = f.sqr(&yy);
                let zz = f.sqr(&p.z);
                let t = f.sub(&f.sub(&f.sqr(&f.add(&p.x, &yy)), &xx), &yyyy);
                let s = f.add(&t, &t);
                let zz2 = f.sqr(&zz);
                let azz2 = match self.dbl {
                    DoublingKind::SmallA(c) => f.mul_small(&zz2, c),
                    _ => f.mul(&self.a, &zz2),
                };
                let m = f.add(&f.add(&f.add(&xx, &xx), &xx), &azz2);
                let x3 = f.sub(&f.sqr(&m), &f.add(&s, &s));
                let y2 = f.add(&yyyy, &yyyy);
                let y4 = f.add(&y2, &y2);
                let y8 = f.add(&y4, &y4);
                let y3 = f.sub(&f.mul(&m, &f.sub(&s, &x3)), &y8);
                let z3 = f.sub(&f.sub(&f.sqr(&f.add(&p.y, &p.z)), &yy), &zz);
                Jacobian {
                    x: x3,
                    y: y3,
                    z: z3,
                }
            }
        }
    }

    /// EFD `add-2007-bl`, 11M + 5S; identity and `P = ±Q` by branches.
    pub fn add(&self, p: &Jacobian<L>, q: &Jacobian<L>) -> Jacobian<L> {
        if p.is_identity() {
            return *q;
        }
        if q.is_identity() {
            return *p;
        }
        let f = &self.f;
        let z1z1 = f.sqr(&p.z);
        let z2z2 = f.sqr(&q.z);
        let u1 = f.mul(&p.x, &z2z2);
        let u2 = f.mul(&q.x, &z1z1);
        let s1 = f.mul(&f.mul(&p.y, &q.z), &z2z2);
        let s2 = f.mul(&f.mul(&q.y, &p.z), &z1z1);
        let h = f.sub(&u2, &u1);
        let r = f.sub(&s2, &s1);
        let r = f.add(&r, &r);
        if h.is_zero() {
            return if r.is_zero() {
                self.double(p)
            } else {
                self.identity()
            };
        }
        let h2 = f.add(&h, &h);
        let i = f.sqr(&h2);
        let j = f.mul(&h, &i);
        let v = f.mul(&u1, &i);
        let x3 = f.sub(&f.sub(&f.sub(&f.sqr(&r), &j), &v), &v);
        let s1j = f.mul(&s1, &j);
        let y3 = f.sub(&f.sub(&f.mul(&r, &f.sub(&v, &x3)), &s1j), &s1j);
        let z3 = f.mul(&f.sub(&f.sub(&f.sqr(&f.add(&p.z, &q.z)), &z1z1), &z2z2), &h);
        Jacobian {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// EFD `madd-2007-bl`, 7M + 4S, affine second operand.
    pub fn add_affine(&self, p: &Jacobian<L>, q: &Affine<L>) -> Jacobian<L> {
        if p.is_identity() {
            return self.from_affine(q);
        }
        let f = &self.f;
        let z1z1 = f.sqr(&p.z);
        let u2 = f.mul(&q.x, &z1z1);
        let s2 = f.mul(&f.mul(&q.y, &p.z), &z1z1);
        let h = f.sub(&u2, &p.x);
        let r = f.sub(&s2, &p.y);
        let r = f.add(&r, &r);
        if h.is_zero() {
            return if r.is_zero() {
                self.double(p)
            } else {
                self.identity()
            };
        }
        let hh = f.sqr(&h);
        let i = f.add(&hh, &hh);
        let i = f.add(&i, &i);
        let j = f.mul(&h, &i);
        let v = f.mul(&p.x, &i);
        let x3 = f.sub(&f.sub(&f.sub(&f.sqr(&r), &j), &v), &v);
        let y1j = f.mul(&p.y, &j);
        let y3 = f.sub(&f.sub(&f.mul(&r, &f.sub(&v, &x3)), &y1j), &y1j);
        let z3 = f.sub(&f.sub(&f.sqr(&f.add(&p.z, &h)), &z1z1), &hh);
        Jacobian {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Affine forms with one inversion (Montgomery's trick), the operation
    /// sequence of `cryptopro_b_point::batch_to_affine`.  Panics on the
    /// identity.
    pub fn batch_to_affine(&self, points: &[Jacobian<L>]) -> Vec<Affine<L>> {
        let n = points.len();
        if n == 0 {
            return Vec::new();
        }
        let f = &self.f;
        let mut prefix = Vec::with_capacity(n);
        let mut acc = f.one();
        for p in points {
            assert!(
                !p.is_identity(),
                "batch_to_affine: the identity has no affine form"
            );
            acc = f.mul(&acc, &p.z);
            prefix.push(acc);
        }
        let mut inv = f.inv(&acc);
        let mut out = vec![
            Affine {
                x: f.zero(),
                y: f.zero()
            };
            n
        ];
        for i in (0..n).rev() {
            let zi = if i == 0 {
                inv
            } else {
                f.mul(&inv, &prefix[i - 1])
            };
            if i > 0 {
                inv = f.mul(&inv, &points[i].z);
            }
            let zi2 = f.sqr(&zi);
            let zi3 = f.mul(&zi2, &zi);
            out[i] = Affine {
                x: f.mul(&points[i].x, &zi2),
                y: f.mul(&points[i].y, &zi3),
            };
        }
        out
    }

    /// `P, 3P, …, (2^(w−1) − 1)·P`, Jacobian: one doubling and
    /// `2^(w−2) − 1` full additions.
    fn odd_multiples(&self, p: &Jacobian<L>, w: u32) -> Vec<Jacobian<L>> {
        let count = 1usize << (w - 2);
        let mut table = Vec::with_capacity(count);
        table.push(*p);
        if count > 1 {
            let p2 = self.double(p);
            for i in 1..count {
                let next = self.add(&table[i - 1], &p2);
                table.push(next);
            }
        }
        table
    }

    #[inline]
    fn entry_affine(&self, table: &[Affine<L>], digit: i8, negate: bool) -> Affine<L> {
        let e = table[(digit.unsigned_abs() >> 1) as usize];
        if (digit < 0) != negate {
            self.neg_affine(&e)
        } else {
            e
        }
    }

    #[inline]
    fn entry_jac(&self, table: &[Jacobian<L>], digit: i8, negate: bool) -> Jacobian<L> {
        let e = table[(digit.unsigned_abs() >> 1) as usize];
        if (digit < 0) != negate {
            self.neg(&e)
        } else {
            e
        }
    }
}

// ── width-w NAF ────────────────────────────────────────────────────────────

/// Width-`w` NAF of `k`, least-significant digit first: the recoding of
/// `cryptopro_b_point::wnaf_digits` for scalars of any length.  `2 ≤ w ≤ 7`.
pub fn wnaf_digits(k: &BigUint, w: u32) -> Vec<i8> {
    assert!((2..=7).contains(&w), "wNAF width must be in 2..=7");
    let mut limbs = k.to_u64_digits();
    limbs.push(0); // `k + (2^w − m)` can carry past the top limb
    let modulus = 1u64 << w;
    let half = 1u64 << (w - 1);
    let mut digits = Vec::with_capacity(64 * limbs.len() + 1);
    while limbs.iter().any(|&l| l != 0) {
        let digit = if limbs[0] & 1 == 1 {
            let m = limbs[0] & (modulus - 1);
            if m >= half {
                add_small(&mut limbs, modulus - m);
                (m as i64 - modulus as i64) as i8
            } else {
                sub_small(&mut limbs, m);
                m as i8
            }
        } else {
            0
        };
        digits.push(digit);
        let n = limbs.len();
        for i in 0..n - 1 {
            limbs[i] = (limbs[i] >> 1) | (limbs[i + 1] << 63);
        }
        limbs[n - 1] >>= 1;
    }
    digits
}

fn add_small(k: &mut [u64], v: u64) {
    let mut carry = v;
    for limb in k.iter_mut() {
        let (s, c) = limb.overflowing_add(carry);
        *limb = s;
        carry = u64::from(c);
        if carry == 0 {
            break;
        }
    }
}

fn sub_small(k: &mut [u64], v: u64) {
    let mut borrow = v;
    for limb in k.iter_mut() {
        let (d, b) = limb.overflowing_sub(borrow);
        *limb = d;
        borrow = u64::from(b);
        if borrow == 0 {
            break;
        }
    }
}

/// How the odd-multiple tables are kept.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TableMode {
    /// One batched inversion, then mixed additions.
    Affine,
    /// No inversion, full additions.
    Jacobian,
}

/// **Baseline arm**: `k·P` by width-`w` NAF with the table kept as `table`
/// says; the result is Jacobian (no final inversion).
pub fn scalar_mul_wnaf<const L: usize>(
    c: &Curve<L>,
    p: &Affine<L>,
    k: &BigUint,
    w: u32,
    table: TableMode,
) -> Jacobian<L> {
    let digits = wnaf_digits(k, w);
    let jt = c.odd_multiples(&c.from_affine(p), w);
    let mut q = c.identity();
    match table {
        TableMode::Affine => {
            let t = c.batch_to_affine(&jt);
            for &d in digits.iter().rev() {
                if !q.is_identity() {
                    q = c.double(&q);
                }
                if d != 0 {
                    q = c.add_affine(&q, &c.entry_affine(&t, d, false));
                }
            }
        }
        TableMode::Jacobian => {
            for &d in digits.iter().rev() {
                if !q.is_identity() {
                    q = c.double(&q);
                }
                if d != 0 {
                    q = c.add(&q, &c.entry_jac(&jt, d, false));
                }
            }
        }
    }
    q
}

// ── the chain endomorphism ─────────────────────────────────────────────────

/// Which chain evaluator computes `φ(P)` (see the module docs).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainEval {
    Generic,
    Optimised,
}

/// One prepared step.
#[derive(Clone, Debug)]
enum Step<const L: usize> {
    /// An odd `ℓ`-step in Kohel form, `x' = N/ψ²`, `y' = y·M/ψ³`, monic
    /// polynomials low degree first.
    Odd {
        psi: Vec<Fe<L>>,
        n: Vec<Fe<L>>,
        m: Vec<Fe<L>>,
    },
    /// A 2-step with kernel `(x₀, 0)`: `x' = x + v/(x − x₀)`,
    /// `y' = y·(1 − v/(x − x₀)²)`, `v = 3x₀² + a`.
    Two { x0: Fe<L>, v: Fe<L> },
}

/// One step's constants as integers (canonical, below `p`), as the loader
/// hands them over: `psi`, `n`, `m` are the Kohel polynomials, low degree
/// first.  For `ell = 2`, `psi = x − x₀` and the step is Vélu's 2-isogeny
/// (see [`PreparedChain::from_steps`] for the accepted forms of `n`, `m`).
#[derive(Clone, Debug)]
pub struct StepInput {
    pub ell: u32,
    pub domain_a: BigUint,
    pub domain_b: BigUint,
    pub codomain_a: BigUint,
    pub codomain_b: BigUint,
    pub psi: Vec<BigUint>,
    pub n: Vec<BigUint>,
    pub m: Vec<BigUint>,
    /// For a 2-step, the exported `(x₀, v)`, checked against `ψ` and `a`.
    pub velu2: Option<(BigUint, BigUint)>,
}

/// A validated chain with its constants in Montgomery form.
#[derive(Clone, Debug)]
pub struct PreparedChain<const L: usize> {
    pub name: String,
    pub degree: u64,
    ells: Vec<u32>,
    steps: Vec<Step<L>>,
    u2: Fe<L>,
    u3: Fe<L>,
    u_inv: Fe<L>,
}

fn poly_fe<const L: usize>(
    f: &GenericField<L>,
    coeffs: &[BigUint],
    want_degree: usize,
    what: &str,
) -> Result<Vec<Fe<L>>, String> {
    if coeffs.len() != want_degree + 1 {
        return Err(format!(
            "polynomial {what} has degree {}, expected {want_degree}",
            coeffs.len() as i64 - 1
        ));
    }
    let v = coeffs
        .iter()
        .map(|c| f.from_canonical(c))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("{what}: {e}"))?;
    if v.last() != Some(&f.one()) {
        return Err(format!("polynomial {what} is not monic"));
    }
    Ok(v)
}

/// `q·(x − x₀)` mod `p`, coefficients low degree first.
fn times_x_minus(q: &[BigUint], x0: &BigUint, p: &BigUint) -> Vec<BigUint> {
    let mut out = vec![BigUint::zero(); q.len() + 1];
    for (i, c) in q.iter().enumerate() {
        out[i + 1] = (&out[i + 1] + c) % p;
        out[i] = sub_mod_big(&out[i], &(c * x0), p);
    }
    out
}

fn modp(v: BigUint, p: &BigUint) -> BigUint {
    v % p
}

fn sub_mod_big(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    ((a % p) + p - (b % p)) % p
}

impl<const L: usize> PreparedChain<L> {
    /// Validate and prepare: steps compose from `E` (the curve's `a`, `b`)
    /// back to a curve that `iso_u` maps onto `E` (`a = u⁴·a_last`,
    /// `b = u⁶·b_last`); `degree` is the product of the steps; odd steps
    /// have monic `ψ`, `N`, `M` of degrees `(ℓ−1)/2`, `ℓ`, `3(ℓ−1)/2`; a
    /// 2-step has `ψ = x − x₀` with `(x₀, 0)` a 2-torsion point of its
    /// domain, Vélu's codomain `a − 5v`, `b − 7x₀v`, and `N`, `M` (when
    /// given) equal to `x² − x₀x + v` and `x² − 2x₀x + x₀² − v` (the forms
    /// `x' = N/ψ`, `y' = y·M/ψ²`).
    pub fn from_steps(
        curve: &Curve<L>,
        name: &str,
        degree: u64,
        steps_in: &[StepInput],
        iso_u: &BigUint,
    ) -> Result<Self, String> {
        let f = &curve.f;
        let p = f.modulus().clone();
        if steps_in.is_empty() {
            return Err("a chain has at least one step".into());
        }
        let mut expect_a = curve.a_big.clone();
        let mut expect_b = curve.b_big.clone();
        let mut product = 1u64;
        let mut steps = Vec::with_capacity(steps_in.len());
        let mut ells = Vec::with_capacity(steps_in.len());
        for (idx, st) in steps_in.iter().enumerate() {
            let ell = st.ell;
            if st.domain_a != expect_a || st.domain_b != expect_b {
                return Err(format!("step {idx}: domain does not compose"));
            }
            for v in [&st.codomain_a, &st.codomain_b] {
                if v >= &p {
                    return Err(format!("step {idx}: codomain coefficient not below p"));
                }
            }
            let l = ell as usize;
            if ell == 2 {
                if st.psi.len() != 2 || !st.psi[1].is_one() || st.psi[0] >= p {
                    return Err(format!("step {idx}: a 2-step needs psi = x - x0"));
                }
                let x0 = (&p - &st.psi[0]) % &p;
                let x0_2 = modp(&x0 * &x0, &p);
                let fx0 = modp(&x0_2 * &x0 + &st.domain_a * &x0 + &st.domain_b, &p);
                if !fx0.is_zero() {
                    return Err(format!("step {idx}: (x0, 0) is not a 2-torsion point"));
                }
                let v = modp(BigUint::from(3u8) * &x0_2 + &st.domain_a, &p);
                let ca = sub_mod_big(&st.domain_a, &(BigUint::from(5u8) * &v), &p);
                let cb = sub_mod_big(&st.domain_b, &(BigUint::from(7u8) * &x0 * &v), &p);
                if ca != st.codomain_a || cb != st.codomain_b {
                    return Err(format!(
                        "step {idx}: codomain is not Velu's (a - 5v, b - 7 x0 v)"
                    ));
                }
                // N and M in either Kohel form: x' = N/ψ, y' = y·M/ψ² (degree
                // 2) or x' = N/ψ², y' = y·M/ψ³ (degree 3, the form the
                // chain sweep exports for every ℓ).
                let n2 = vec![v.clone(), (&p - &x0) % &p, BigUint::one()];
                let m2 = vec![
                    sub_mod_big(&x0_2, &v, &p),
                    sub_mod_big(&BigUint::zero(), &(BigUint::from(2u8) * &x0), &p),
                    BigUint::one(),
                ];
                let (n3, m3) = (times_x_minus(&n2, &x0, &p), times_x_minus(&m2, &x0, &p));
                if !st.n.is_empty() && st.n != n2 && st.n != n3 {
                    return Err(format!("step {idx}: N is not Velu's x-map numerator"));
                }
                if !st.m.is_empty() && st.m != m2 && st.m != m3 {
                    return Err(format!("step {idx}: M is not Velu's y-map numerator"));
                }
                if let Some((vx0, vv)) = &st.velu2 {
                    if vx0 != &x0 || vv != &v {
                        return Err(format!("step {idx}: velu2 (x0, v) disagrees with psi"));
                    }
                }
                steps.push(Step::Two {
                    x0: f.from_biguint(&x0),
                    v: f.from_biguint(&v),
                });
            } else {
                if ell < 3 || ell % 2 == 0 {
                    return Err(format!("step {idx}: ell must be 2 or an odd prime"));
                }
                steps.push(Step::Odd {
                    psi: poly_fe(f, &st.psi, (l - 1) / 2, "psi")?,
                    n: poly_fe(f, &st.n, l, "N")?,
                    m: poly_fe(f, &st.m, (3 * l - 3) / 2, "M")?,
                });
            }
            ells.push(ell);
            expect_a = st.codomain_a.clone();
            expect_b = st.codomain_b.clone();
            product = product
                .checked_mul(u64::from(ell))
                .ok_or("degree overflow")?;
        }
        if product != degree {
            return Err("chain degree is not the product of the steps".into());
        }
        let u = f.from_canonical(iso_u).map_err(|e| format!("iso_u: {e}"))?;
        if u.is_zero() {
            return Err("iso_u must be a unit".into());
        }
        let u2 = f.sqr(&u);
        let u3 = f.mul(&u2, &u);
        let u4 = f.sqr(&u2);
        let u6 = f.mul(&u4, &u2);
        if f.mul(&u4, &f.from_biguint(&expect_a)) != curve.a
            || f.mul(&u6, &f.from_biguint(&expect_b)) != curve.b
        {
            return Err("iso_u does not map the last codomain back to E".into());
        }
        Ok(PreparedChain {
            name: name.to_string(),
            degree,
            ells,
            steps,
            u2,
            u3,
            u_inv: f.inv(&u),
        })
    }

    pub fn step_degrees(&self) -> Vec<u32> {
        self.ells.clone()
    }

    /// `φ(P)` by the chosen evaluator, as a Jacobian point on `E`.
    pub fn apply_with(&self, c: &Curve<L>, p: &Affine<L>, eval: ChainEval) -> Jacobian<L> {
        match eval {
            ChainEval::Generic => self.apply_generic(c, p),
            ChainEval::Optimised => self.apply_optimised(c, p),
        }
    }

    fn apply_generic(&self, c: &Curve<L>, p: &Affine<L>) -> Jacobian<L> {
        let f = &c.f;
        let (mut x, mut y, mut z) = (p.x, p.y, f.one());
        for (i, step) in self.steps.iter().enumerate() {
            (x, y, z) = match step {
                Step::Odd { psi, n, m } => step_projective(f, psi, n, m, &x, &y, &z),
                Step::Two { x0, v } => two_projective(f, x0, v, &x, &y, &z, i == 0),
            };
        }
        x = f.mul(&x, &self.u2);
        y = f.mul(&y, &self.u3);
        let z2 = f.sqr(&z);
        Jacobian {
            x: f.mul(&x, &z),
            y: f.mul(&y, &z2),
            z,
        }
    }

    fn apply_optimised(&self, c: &Curve<L>, p: &Affine<L>) -> Jacobian<L> {
        let f = &c.f;
        let (first, rest) = self.steps.split_first().expect("non-empty chain");
        let (mut x, mut y, mut z) = match first {
            Step::Odd { psi, n, m } => {
                let ps = horner_monic(f, psi, &p.x);
                let nu = horner_monic(f, n, &p.x);
                let mu = horner_monic(f, m, &p.x);
                (nu, f.mul(&p.y, &mu), ps)
            }
            Step::Two { x0, v } => two_jacobian_affine(f, x0, v, &p.x, &p.y),
        };
        for step in rest {
            (x, y, z) = match step {
                Step::Odd { psi, n, m } => step_jacobian(f, psi, n, m, &x, &y, &z),
                Step::Two { x0, v } => two_jacobian(f, x0, v, &x, &y, &z),
            };
        }
        Jacobian {
            x,
            y,
            z: f.mul(&z, &self.u_inv),
        }
    }
}

#[inline]
fn horner_homogeneous<const L: usize>(
    f: &GenericField<L>,
    coeffs: &[Fe<L>],
    x: &Fe<L>,
    zpow: &[Fe<L>],
) -> Fe<L> {
    let d = coeffs.len() - 1;
    let mut acc = coeffs[d];
    for i in (0..d).rev() {
        acc = f.add(&f.mul(&acc, x), &f.mul(&coeffs[i], &zpow[d - i]));
    }
    acc
}

#[inline]
fn horner_monic<const L: usize>(f: &GenericField<L>, coeffs: &[Fe<L>], x: &Fe<L>) -> Fe<L> {
    let d = coeffs.len() - 1;
    let mut acc = f.add(x, &coeffs[d - 1]);
    for i in (0..d - 1).rev() {
        acc = f.add(&f.mul(&acc, x), &coeffs[i]);
    }
    acc
}

#[inline]
fn horner_monic_homogeneous<const L: usize>(
    f: &GenericField<L>,
    coeffs: &[Fe<L>],
    x: &Fe<L>,
    wpow: &[Fe<L>],
) -> Fe<L> {
    let d = coeffs.len() - 1;
    let mut acc = f.add(x, &f.mul(&coeffs[d - 1], &wpow[1]));
    for i in (0..d - 1).rev() {
        acc = f.add(&f.mul(&acc, x), &f.mul(&coeffs[i], &wpow[d - i]));
    }
    acc
}

fn step_jacobian<const L: usize>(
    f: &GenericField<L>,
    psi: &[Fe<L>],
    n: &[Fe<L>],
    m: &[Fe<L>],
    x: &Fe<L>,
    y: &Fe<L>,
    z: &Fe<L>,
) -> (Fe<L>, Fe<L>, Fe<L>) {
    let max_deg = m.len() - 1;
    let w = f.sqr(z);
    let mut wpow = Vec::with_capacity(max_deg + 1);
    wpow.push(f.one());
    wpow.push(w);
    for i in 2..=max_deg {
        let next = f.mul(&wpow[i - 1], &w);
        wpow.push(next);
    }
    let ps = horner_monic_homogeneous(f, psi, x, &wpow);
    let nu = horner_monic_homogeneous(f, n, x, &wpow);
    let mu = horner_monic_homogeneous(f, m, x, &wpow);
    (nu, f.mul(y, &mu), f.mul(z, &ps))
}

fn step_projective<const L: usize>(
    f: &GenericField<L>,
    psi: &[Fe<L>],
    n: &[Fe<L>],
    m: &[Fe<L>],
    x: &Fe<L>,
    y: &Fe<L>,
    z: &Fe<L>,
) -> (Fe<L>, Fe<L>, Fe<L>) {
    let max_deg = m.len() - 1;
    let mut zpow = Vec::with_capacity(max_deg + 1);
    zpow.push(f.one());
    zpow.push(*z);
    for i in 2..=max_deg {
        let next = f.mul(&zpow[i - 1], z);
        zpow.push(next);
    }
    let ps = horner_homogeneous(f, psi, x, &zpow);
    let nu = horner_homogeneous(f, n, x, &zpow);
    let mu = horner_homogeneous(f, m, x, &zpow);
    let ps3 = f.mul(&f.sqr(&ps), &ps);
    (f.mul(&nu, &ps), f.mul(y, &mu), f.mul(z, &ps3))
}

/// A 2-step on an affine input, Jacobian output: `d = x − x₀`,
/// `X' = (x·d + v)·d`, `Y' = y·(d² − v)·d`, `Z' = d`; 4M + 1S.
fn two_jacobian_affine<const L: usize>(
    f: &GenericField<L>,
    x0: &Fe<L>,
    v: &Fe<L>,
    x: &Fe<L>,
    y: &Fe<L>,
) -> (Fe<L>, Fe<L>, Fe<L>) {
    let d = f.sub(x, x0);
    let xn = f.mul(&f.add(&f.mul(x, &d), v), &d);
    let dd = f.sqr(&d);
    let yn = f.mul(&f.mul(y, &f.sub(&dd, v)), &d);
    (xn, yn, d)
}

/// A 2-step on a Jacobian input: `d = X − x₀·Z²`, `X' = (X·d + v·Z⁴)·d`,
/// `Y' = Y·(d² − v·Z⁴)·d`, `Z' = Z·d`; 7M + 3S.
fn two_jacobian<const L: usize>(
    f: &GenericField<L>,
    x0: &Fe<L>,
    v: &Fe<L>,
    x: &Fe<L>,
    y: &Fe<L>,
    z: &Fe<L>,
) -> (Fe<L>, Fe<L>, Fe<L>) {
    let z2 = f.sqr(z);
    let z4 = f.sqr(&z2);
    let d = f.sub(x, &f.mul(x0, &z2));
    let vz4 = f.mul(v, &z4);
    let xn = f.mul(&f.add(&f.mul(x, &d), &vz4), &d);
    let dd = f.sqr(&d);
    let yn = f.mul(&f.mul(y, &f.sub(&dd, &vz4)), &d);
    (xn, yn, f.mul(z, &d))
}

/// A 2-step in projective coordinates (`x = X/Z`): `d = X − x₀·Z`,
/// `X' = (X·d + v·Z²)·d`, `Y' = Y·(d² − v·Z²)`, `Z' = Z·d²`; 6M + 2S, and
/// 3M + 1S when the input is the affine point (`Z = 1`, `first`).
fn two_projective<const L: usize>(
    f: &GenericField<L>,
    x0: &Fe<L>,
    v: &Fe<L>,
    x: &Fe<L>,
    y: &Fe<L>,
    z: &Fe<L>,
    first: bool,
) -> (Fe<L>, Fe<L>, Fe<L>) {
    if first {
        let d = f.sub(x, x0);
        let xn = f.mul(&f.add(&f.mul(x, &d), v), &d);
        let dd = f.sqr(&d);
        let yn = f.mul(y, &f.sub(&dd, v));
        return (xn, yn, dd);
    }
    let d = f.sub(x, &f.mul(x0, z));
    let vz2 = f.mul(v, &f.sqr(z));
    let xn = f.mul(&f.add(&f.mul(x, &d), &vz2), &d);
    let dd = f.sqr(&d);
    let yn = f.mul(y, &f.sub(&dd, &vz2));
    (xn, yn, f.mul(z, &dd))
}

// ── GLV ────────────────────────────────────────────────────────────────────

/// The input-independent part of a GLV arm: the chain, `λ`, `n`, the reduced
/// basis and its determinant (built once, outside every timed arm).
#[derive(Clone, Debug)]
pub struct GlvContext<const L: usize> {
    pub chain: PreparedChain<L>,
    pub lambda: BigUint,
    pub order: BigUint,
    basis: [[BigInt; 2]; 2],
    det: BigInt,
    pub babai_bound_bits: u32,
}

fn round_div(num: &BigInt, den: &BigInt) -> BigInt {
    let (num, den) = if den.sign() == Sign::Minus {
        (-num, -den)
    } else {
        (num.clone(), den.clone())
    };
    (num * 2u8 + &den).div_floor(&(den * 2u8))
}

impl<const L: usize> GlvContext<L> {
    /// Checks `λ < n`, a non-singular basis with `|det| = n`, and both rows in
    /// `{(v1, v2) : v1 + v2·λ ≡ 0 (mod n)}`.
    pub fn from_parts(
        chain: PreparedChain<L>,
        lambda: &BigUint,
        order: &BigUint,
        basis: [[BigInt; 2]; 2],
        babai_bound_bits: u32,
    ) -> Result<Self, String> {
        if lambda >= order {
            return Err("lambda must be reduced below n".into());
        }
        let det = &basis[0][0] * &basis[1][1] - &basis[0][1] * &basis[1][0];
        let n_int = BigInt::from(order.clone());
        if det.magnitude() != order {
            return Err("the GLV basis does not have determinant ±n".into());
        }
        let lam = BigInt::from(lambda.clone());
        for row in &basis {
            if !(&row[0] + &row[1] * &lam).mod_floor(&n_int).is_zero() {
                return Err("a GLV basis row is not in the lattice of lambda".into());
            }
        }
        Ok(GlvContext {
            chain,
            lambda: lambda.clone(),
            order: order.clone(),
            basis,
            det,
            babai_bound_bits,
        })
    }

    /// Babai rounding, the decomposition of `cryptopro_b_point::GlvContext`.
    pub fn decompose(&self, k: &BigUint) -> (BigInt, BigInt) {
        let k = BigInt::from(k.clone());
        let r1 = round_div(&(&k * &self.basis[1][1]), &self.det);
        let r2 = round_div(&(-&k * &self.basis[0][1]), &self.det);
        let k1 = &k - &r1 * &self.basis[0][0] - &r2 * &self.basis[1][0];
        let k2 = -(&r1 * &self.basis[0][1]) - &r2 * &self.basis[1][1];
        (k1, k2)
    }

    pub fn decomposition_is_valid(&self, k: &BigUint, k1: &BigInt, k2: &BigInt) -> bool {
        let n = BigInt::from(self.order.clone());
        let lam = BigInt::from(self.lambda.clone());
        let lhs = (k1 + k2 * lam).mod_floor(&n);
        let rhs = BigInt::from(k.clone()).mod_floor(&n);
        let bound = u64::from(self.babai_bound_bits);
        lhs == rhs && k1.bits() <= bound && k2.bits() <= bound
    }
}

/// `k1·P + k2·φ(P)`: interleaved two-scalar width-`w` NAF with shared
/// doublings, the tables kept as `table` says (both converted in one batched
/// inversion when affine).
pub fn glv_combine<const L: usize>(
    c: &Curve<L>,
    p: &Affine<L>,
    phi_p: &Jacobian<L>,
    k1: &BigInt,
    k2: &BigInt,
    w: u32,
    table: TableMode,
) -> Jacobian<L> {
    let (neg1, neg2) = (k1.sign() == Sign::Minus, k2.sign() == Sign::Minus);
    let d1 = wnaf_digits(k1.magnitude(), w);
    let d2 = wnaf_digits(k2.magnitude(), w);
    let len = d1.len().max(d2.len());
    let mut q = c.identity();
    let t1 = c.odd_multiples(&c.from_affine(p), w);
    let t2 = c.odd_multiples(phi_p, w);
    match table {
        TableMode::Affine => {
            let count = t1.len();
            let mut both = t1;
            both.extend(t2);
            let t = c.batch_to_affine(&both);
            let (a1, a2) = t.split_at(count);
            for i in (0..len).rev() {
                if !q.is_identity() {
                    q = c.double(&q);
                }
                let a = d1.get(i).copied().unwrap_or(0);
                if a != 0 {
                    q = c.add_affine(&q, &c.entry_affine(a1, a, neg1));
                }
                let b = d2.get(i).copied().unwrap_or(0);
                if b != 0 {
                    q = c.add_affine(&q, &c.entry_affine(a2, b, neg2));
                }
            }
        }
        TableMode::Jacobian => {
            for i in (0..len).rev() {
                if !q.is_identity() {
                    q = c.double(&q);
                }
                let a = d1.get(i).copied().unwrap_or(0);
                if a != 0 {
                    q = c.add(&q, &c.entry_jac(&t1, a, neg1));
                }
                let b = d2.get(i).copied().unwrap_or(0);
                if b != 0 {
                    q = c.add(&q, &c.entry_jac(&t2, b, neg2));
                }
            }
        }
    }
    q
}

/// **GLV arm**: decomposition, `φ(P)` by `eval`, then [`glv_combine`].
pub fn scalar_mul_glv<const L: usize>(
    c: &Curve<L>,
    ctx: &GlvContext<L>,
    p: &Affine<L>,
    k: &BigUint,
    w: u32,
    table: TableMode,
    eval: ChainEval,
) -> Jacobian<L> {
    let (k1, k2) = ctx.decompose(k);
    let phi_p = ctx.chain.apply_with(c, p, eval);
    glv_combine(c, p, &phi_p, &k1, &k2, w, table)
}

/// `k·P` by the crate's textbook `BigUint` arithmetic
/// ([`crate::ecc::point::Point::scalar_mul_vartime`]), as affine integers;
/// `None` for the identity.  The independent reference of every check.
pub fn textbook_mul<const L: usize>(
    c: &Curve<L>,
    p: &Affine<L>,
    k: &BigUint,
) -> Option<(BigUint, BigUint)> {
    use crate::ecc::field::FieldElement;
    use crate::ecc::point::Point;
    let m = c.f.modulus().clone();
    let (x, y) = c.affine_to_biguint(p);
    let tp = Point::Affine {
        x: FieldElement::new(x, m.clone()),
        y: FieldElement::new(y, m.clone()),
    };
    match tp.scalar_mul_vartime(k, &FieldElement::new(c.a_big.clone(), m)) {
        Point::Infinity => None,
        Point::Affine { x, y } => Some((x.value, y.value)),
    }
}

/// The Jacobian point's affine integers (one inversion); `None` for the
/// identity.
pub fn jacobian_to_biguint<const L: usize>(
    c: &Curve<L>,
    p: &Jacobian<L>,
) -> Option<(BigUint, BigUint)> {
    c.to_affine(p).map(|a| c.affine_to_biguint(&a))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::cryptopro_b_point as cpb;
    use crate::ecc::CurveParams;
    use num_bigint::RandBigInt;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn cpb_curve() -> (Curve<4>, CurveParams) {
        let cp = CurveParams::gost_cryptopro_b();
        (Curve::<4>::new(&cp.p, &cp.a, &cp.b).unwrap(), cp)
    }

    /// A curve with general `a` over the P-384 prime (not a standard curve;
    /// only its group law is used), and a point on it found by search.
    fn p384_general_curve() -> (Curve<6>, Affine<6>) {
        let p = BigUint::parse_bytes(
            b"fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff",
            16,
        )
        .unwrap();
        let b = BigUint::from(7u8);
        // the first a from a fixed start that is neither small nor
        // isomorphic to an a = -3 model
        let mut a = BigUint::parse_bytes(b"123456789abcdef0123456789abcdef", 16).unwrap();
        let c = loop {
            if let Ok(c) = Curve::<6>::new(&p, &a, &b) {
                if c.dbl == DoublingKind::General {
                    break c;
                }
            }
            a += 1u8;
        };
        // p ≡ 3 (mod 4): sqrt by (p+1)/4
        let e = (&p + 1u8) / 4u8;
        let mut x = BigUint::from(2u8);
        loop {
            let rhs = (&x * &x * &x + &a * &x + &b) % &p;
            let y = rhs.modpow(&e, &p);
            if (&y * &y) % &p == rhs {
                return (c.clone(), c.affine(&x, &y).unwrap());
            }
            x += 1u8;
        }
    }

    fn jac_is(c: &Curve<4>, j: &Jacobian<4>, want: Option<(BigUint, BigUint)>) -> bool {
        jacobian_to_biguint(c, j) == want
    }

    #[test]
    fn doubling_kinds_follow_the_model() {
        let cp = CurveParams::gost_cryptopro_b();
        assert_eq!(doubling_kind(&cp.p, &cp.a).unwrap(), DoublingKind::AMinus3);
        assert_eq!(
            doubling_kind(&cp.p, &BigUint::zero()).unwrap(),
            DoublingKind::AZero
        );
        let s = CurveParams::secp256k1();
        assert_eq!(doubling_kind(&s.p, &s.a).unwrap(), DoublingKind::AZero);
        // a = 1 over a p ≡ 1 (mod 4) prime where -3 is not a fourth power
        // is "small"; the answer is checked against the fourth-power test.
        let p = &cp.p;
        for a in [1u32, 2, 5, 7] {
            let a = BigUint::from(a);
            let t = ((p - 3u8) * a.modpow(&(p - 2u8), p)) % p;
            match doubling_kind(p, &a) {
                Ok(DoublingKind::SmallA(v)) => {
                    assert!(!is_fourth_power(&t, p));
                    assert_eq!(BigUint::from(v as u64), a);
                }
                Err(_) => assert!(is_fourth_power(&t, p)),
                other => panic!("{other:?}"),
            }
        }
    }

    #[test]
    fn group_law_matches_textbook_on_cryptopro_b_and_a_general_a_curve() {
        let (c, cp) = cpb_curve();
        let g = c.affine(&cp.gx, &cp.gy).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..8 {
            let k = rng.gen_biguint_below(&cp.n);
            for w in [3u32, 5, 7] {
                for t in [TableMode::Affine, TableMode::Jacobian] {
                    let got = scalar_mul_wnaf(&c, &g, &k, w, t);
                    assert!(jac_is(&c, &got, textbook_mul(&c, &g, &k)));
                }
            }
            // same result as the specialised CryptoPro-B code
            let g_cpb = cpb::CryptoProBAffine::from_textbook(&cp.generator()).unwrap();
            let want = cpb::scalar_mul_wnaf(&g_cpb, &k, 5).to_affine().unwrap();
            let got = c
                .to_affine(&scalar_mul_wnaf(&c, &g, &k, 5, TableMode::Affine))
                .unwrap();
            assert_eq!(c.f.to_biguint(&got.x), want.x.to_biguint());
            assert_eq!(c.f.to_biguint(&got.y), want.y.to_biguint());
        }
        let (c6, q) = p384_general_curve();
        assert!(c6.is_on_curve(&q));
        let jq = c6.from_affine(&q);
        let d = c6.double(&jq);
        let t = textbook_mul(&c6, &q, &BigUint::from(2u8));
        assert_eq!(jacobian_to_biguint(&c6, &d), t);
        for k in [3u32, 5, 77, 1001, 65537] {
            let k = BigUint::from(k);
            for w in [3u32, 4, 6] {
                for tm in [TableMode::Affine, TableMode::Jacobian] {
                    let got = scalar_mul_wnaf(&c6, &q, &k, w, tm);
                    assert_eq!(jacobian_to_biguint(&c6, &got), textbook_mul(&c6, &q, &k));
                }
            }
        }
        // add and madd agree with Z ≠ 1
        let d2 = c6.double(&d);
        assert!(c6.eq_point(&c6.add(&d2, &jq), &c6.add_affine(&d2, &q)));
        assert!(c6.add(&jq, &c6.neg(&jq)).is_identity());
        assert!(c6.eq_point(&c6.add(&jq, &jq), &d));
    }

    #[test]
    fn small_a_and_a_zero_doublings_match_textbook() {
        // a = 0: secp256k1; small a: a = 5 over P-384's prime, if it is not
        // isomorphic to an a = -3 model there.
        let s = CurveParams::secp256k1();
        let c = Curve::<4>::new(&s.p, &s.a, &s.b).unwrap();
        let g = c.affine(&s.gx, &s.gy).unwrap();
        for k in [2u32, 3, 1234567, 0xfffffff] {
            let k = BigUint::from(k);
            let got = scalar_mul_wnaf(&c, &g, &k, 4, TableMode::Jacobian);
            assert!(jac_is(&c, &got, textbook_mul(&c, &g, &k)));
        }
        let p = BigUint::parse_bytes(
            b"fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff",
            16,
        )
        .unwrap();
        for a in [5u32, 7, 11, 13] {
            let a = BigUint::from(a);
            let b = BigUint::from(3u8);
            let Ok(c) = Curve::<6>::new(&p, &a, &b) else {
                continue;
            };
            let DoublingKind::SmallA(_) = c.dbl else {
                continue;
            };
            let e = (&p + 1u8) / 4u8;
            let mut x = BigUint::from(1u8);
            let q = loop {
                let rhs = (&x * &x * &x + &a * &x + &b) % &p;
                let y = rhs.modpow(&e, &p);
                if (&y * &y) % &p == rhs {
                    break c.affine(&x, &y).unwrap();
                }
                x += 1u8;
            };
            for k in [2u32, 9, 31337] {
                let k = BigUint::from(k);
                let got = scalar_mul_wnaf(&c, &q, &k, 5, TableMode::Affine);
                assert_eq!(jacobian_to_biguint(&c, &got), textbook_mul(&c, &q, &k));
            }
            return;
        }
        panic!("no small-a curve found");
    }

    #[test]
    fn wnaf_matches_the_cryptopro_b_recoding_and_reconstructs() {
        let mut rng = StdRng::seed_from_u64(11);
        for _ in 0..50 {
            let k = rng.gen_biguint(256);
            for w in 2..=7 {
                let d = wnaf_digits(&k, w);
                let want = cpb::wnaf_digits(&crate::ct_bignum::U256::from_biguint(&k), w);
                assert_eq!(d, want);
            }
        }
        for _ in 0..20 {
            let k = rng.gen_biguint(782);
            let d = wnaf_digits(&k, 5);
            let mut acc = BigInt::zero();
            for &x in d.iter().rev() {
                acc = acc * 2 + BigInt::from(x);
            }
            assert_eq!(acc, BigInt::from(k));
        }
    }

    #[test]
    fn two_isogeny_evaluators_agree_with_velu() {
        // Over P-384's prime, find a curve with a rational 2-torsion point:
        // y² = (x − x0)(x² + x0 x + c) = x³ + (c − x0²) x − x0 c.
        let p = BigUint::parse_bytes(
            b"fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff",
            16,
        )
        .unwrap();
        let x0 = BigUint::from(1234567u32);
        let cc = BigUint::from(98765u32);
        let a = sub_mod_big(&cc, &(&x0 * &x0), &p);
        let b = sub_mod_big(&BigUint::zero(), &(&x0 * &cc), &p);
        let c = Curve::<6>::new(&p, &a, &b).unwrap();
        let v = (BigUint::from(3u8) * &x0 * &x0 + &a) % &p;
        let ca = sub_mod_big(&a, &(BigUint::from(5u8) * &v), &p);
        let cb = sub_mod_big(&b, &(BigUint::from(7u8) * &x0 * &v), &p);
        let st = StepInput {
            ell: 2,
            domain_a: a.clone(),
            domain_b: b.clone(),
            codomain_a: ca.clone(),
            codomain_b: cb.clone(),
            psi: vec![(&p - &x0) % &p, BigUint::one()],
            n: vec![],
            m: vec![],
            velu2: None,
        };
        // a "chain" into the codomain, closed with u = 1 only if the
        // codomain is E itself, which it is not: test the step maps
        // directly through two_* instead of PreparedChain.
        let f = &c.f;
        let (x0f, vf) = (f.from_biguint(&x0), f.from_biguint(&v));
        let e = (&p + 1u8) / 4u8;
        let mut x = BigUint::from(5u8);
        let mut done = 0;
        while done < 5 {
            let rhs = (&x * &x * &x + &a * &x + &b) % &p;
            let y = rhs.modpow(&e, &p);
            if (&y * &y) % &p == rhs && x != x0 {
                let pt = c.affine(&x, &y).unwrap();
                // Velu, affine, with BigUint
                let dinv = sub_mod_big(&x, &x0, &p).modpow(&(&p - 2u8), &p);
                let xp = (&x + &v * &dinv) % &p;
                let yp = (&y * sub_mod_big(&BigUint::one(), &(&v * &dinv * &dinv), &p)) % &p;
                let check = |(xx, yy, zz): (Fe<6>, Fe<6>, Fe<6>), proj: bool| {
                    let zi = f.inv(&zz);
                    let (xa, ya) = if proj {
                        (f.mul(&xx, &zi), f.mul(&yy, &zi))
                    } else {
                        let zi2 = f.sqr(&zi);
                        (f.mul(&xx, &zi2), f.mul(&yy, &f.mul(&zi2, &zi)))
                    };
                    assert_eq!(f.to_biguint(&xa), xp);
                    assert_eq!(f.to_biguint(&ya), yp);
                };
                check(two_jacobian_affine(f, &x0f, &vf, &pt.x, &pt.y), false);
                check(
                    two_projective(f, &x0f, &vf, &pt.x, &pt.y, &f.one(), true),
                    true,
                );
                // a non-trivial Z: (λ²x, λ³y, λ) Jacobian and (λx, λy, λ) projective
                let lam = f.from_biguint(&BigUint::from(777u32));
                let l2 = f.sqr(&lam);
                let l3 = f.mul(&l2, &lam);
                check(
                    two_jacobian(f, &x0f, &vf, &f.mul(&pt.x, &l2), &f.mul(&pt.y, &l3), &lam),
                    false,
                );
                check(
                    two_projective(
                        f,
                        &x0f,
                        &vf,
                        &f.mul(&pt.x, &lam),
                        &f.mul(&pt.y, &lam),
                        &lam,
                        false,
                    ),
                    true,
                );
                // the image lies on Velu's codomain
                let img = Affine {
                    x: f.from_biguint(&xp),
                    y: f.from_biguint(&yp),
                };
                assert!(c.on_curve_with(&img, &f.from_biguint(&ca), &f.from_biguint(&cb)));
                done += 1;
            }
            x += 1u8;
        }
        // the loader accepts the step's own Velu N and M and refuses others
        let mut good = st.clone();
        good.n = vec![v.clone(), (&p - &x0) % &p, BigUint::one()];
        good.m = vec![
            sub_mod_big(&((&x0 * &x0) % &p), &v, &p),
            sub_mod_big(&BigUint::zero(), &(BigUint::from(2u8) * &x0), &p),
            BigUint::one(),
        ];
        let err = PreparedChain::<6>::from_steps(&c, "t", 2, &[good.clone()], &BigUint::one())
            .unwrap_err();
        assert!(err.contains("iso_u"), "{err}");
        let mut bad = good;
        bad.n[0] = (&bad.n[0] + 1u8) % &p;
        let err = PreparedChain::<6>::from_steps(&c, "t", 2, &[bad], &BigUint::one()).unwrap_err();
        assert!(err.contains("N is not"), "{err}");
        // the degree-3 forms N·ψ, M·ψ are accepted too
        let mut g3 = st.clone();
        g3.n = times_x_minus(&[v.clone(), (&p - &x0) % &p, BigUint::one()], &x0, &p);
        g3.velu2 = Some((x0.clone(), v.clone()));
        let err = PreparedChain::<6>::from_steps(&c, "t", 2, &[g3], &BigUint::one()).unwrap_err();
        assert!(err.contains("iso_u"), "{err}");
    }

    #[cfg(feature = "endo-chain-opcount")]
    #[test]
    fn formula_costs_are_the_declared_ones() {
        use crate::ecc::endo_generic_field::opcount;
        let (c, cp) = cpb_curve();
        let g = c.from_affine(&c.affine(&cp.gx, &cp.gy).unwrap());
        let g2 = c.double(&g);
        let ga = c.to_affine(&g2).unwrap();
        let g4 = c.double(&g2);
        let count = |f: &dyn Fn()| {
            opcount::reset();
            f();
            let o = opcount::read();
            (o.mul, o.sqr)
        };
        assert_eq!(
            count(&|| {
                c.double(&g2);
            }),
            (3, 5)
        );
        assert_eq!(
            count(&|| {
                c.add(&g4, &g2);
            }),
            (11, 5)
        );
        assert_eq!(
            count(&|| {
                c.add_affine(&g4, &ga);
            }),
            (7, 4)
        );
        let (c6, q) = p384_general_curve();
        let j = c6.double(&c6.from_affine(&q));
        assert_eq!(
            count(&|| {
                c6.double(&j);
            }),
            (2, 8)
        );
        let s = CurveParams::secp256k1();
        let ck = Curve::<4>::new(&s.p, &s.a, &s.b).unwrap();
        let jk = ck.double(&ck.from_affine(&ck.affine(&s.gx, &s.gy).unwrap()));
        assert_eq!(
            count(&|| {
                ck.double(&jk);
            }),
            (2, 5)
        );
    }
}
