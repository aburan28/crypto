//! The wide auxiliary-curve search of the PKM isogeny tower: a curve
//! `E'/F_p` with a rational point of order `2^t`, found by 2-Sylow
//! halving instead of point counting.
//!
//! # Setting
//!
//! `p ≡ 3 (mod 4)` (NIST P-192: `v₂(p − 1) = 1`, `v₂(p + 1) = 64`). Every
//! curve drawn has full rational 2-torsion,
//! `E : y² = (x − e₁)(x − e₂)(x − e₃)`. Because `μ₄ ⊄ F_p`, `E[4] ⊄
//! E(F_p)` (Weil pairing), so the 2-Sylow subgroup is
//! `E(F_p)[2^∞] ≅ Z/2 × Z/2^b` with `b ≥ 1`, and `E` has a point of order
//! `2^t` exactly when `b ≥ t`. Also `v₂(#E) = b + 1`.
//!
//! # Halving
//!
//! A point `P = (x, y)` is in `2E(F_p)` iff `x − e₁, x − e₂, x − e₃` are
//! all squares. With `rᵢ = √(x − eᵢ)`, the four halves of `±P` have
//! `x_H = x + r₁r₂ + r₁r₃ + r₂r₃` over the sign choices of the `rᵢ`
//! (up to a global sign), `x_H − eᵢ = (rᵢ + rⱼ)(rᵢ + r_k)` and
//! `y_H = ±(r₁ + r₂)(r₁ + r₃)(r₂ + r₃)`. Negating `r_n` alone gives
//! `H + T_n`, `T_n = (e_n, 0)`.
//!
//! The descent map `δ(P) = (χ(x − e₁), χ(x − e₂), χ(x − e₃))` (`χ` the
//! Legendre symbol, with the product of the other two in place of a zero)
//! is an isomorphism `E(F_p)/2E(F_p) → {±1}²`. When `b ≥ 2` exactly one
//! 2-torsion point `T_h` is in `2E`; the other two share `d = δ(T_n)`,
//! which has exactly one component `d_k = +1`. Some half of a halvable `P`
//! is itself halvable iff `δ(H) ∈ {1, d}`, iff `χ(x_H − e_k) = +1`: one
//! Legendre symbol per level. Then `H` (if `δ(H) = 1`) or `H + T_n`
//! (otherwise) is the next point of the chain, whose square roots cost two
//! exponentiations: `i₁ = (x − e₁)^((p−3)/4) = 1/r₁`,
//! `u = ((x − e₁)(x − e₂))^((p−3)/4) = 1/(r₁r₂)`, `r₂ = (x − e₂)·u·r₁`,
//! `r₃ = y·u`.
//!
//! # Families
//!
//! - [`Family::Z2xZ2`]: `e₁ = 0, e₂ = α, e₃ = β` with `α, β` uniform.
//! - [`Family::Z2xZ8`]: the genus-0 modular curve `X₁(2, 8)`. With
//!   `A = m² − 1`, `B = 2m`, `C = m² + 1` (so `A² + B² = C²`),
//!   `Y² = X(X + B⁴)(X + A⁴)` has the point `H₄ = (A²B², A²B²C²)` of order
//!   4 with `r = (AB, BC, AC)`, and so `b ≥ 3` for every `m`.
//!
//! [`Arm::Baseline`] is the direct implementation (three square roots per
//! level through [`PrimeField::sqrt`], the four halves tested one by one
//! with [`PrimeField::legendre`]); [`Arm::Fast`] is the one above, with
//! [`HalvingField`]'s addition-chain exponentiation and Jacobi symbol. Both
//! return the same `b` on every draw (tested).

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::fp_backend::{jacobi_limbs, Fp64, FpP192, PrimeField};

// ── Field hooks ───────────────────────────────────────────────────────

/// The two primitives the fast arm specialises: `a^((p−3)/4)` and the
/// Legendre symbol. The defaults are the generic ones.
pub trait HalvingField: PrimeField {
    /// `a^((p − 3)/4)` (for a non-zero square `a`, `±1/√a`). `e` holds
    /// `(p − 3)/4` as little-endian limbs.
    fn pow_isqrt(&self, a: Self::Elem, e: &[u64]) -> Self::Elem {
        self.pow(a, e)
    }
    /// The Legendre symbol, by the fastest available method.
    fn chi(&self, a: Self::Elem) -> i32 {
        self.legendre(a)
    }
}

impl HalvingField for Fp64 {}

impl HalvingField for FpP192 {
    fn pow_isqrt(&self, a: [u64; 3], _e: &[u64]) -> [u64; 3] {
        p192_pow_isqrt(self, a)
    }
    fn chi(&self, a: [u64; 3]) -> i32 {
        jacobi_p192_fast(&a)
    }
}

/// `a^((p−3)/4)` for P-192 by an addition chain: `(p − 3)/4 = 2^190 −
/// 2^62 − 1 = 2^63·(2^127 − 1) + (2^62 − 1)`, 189 squarings and 9
/// multiplications.
pub fn p192_pow_isqrt(f: &FpP192, a: [u64; 3]) -> [u64; 3] {
    let sqn = |mut x: [u64; 3], n: u32| {
        for _ in 0..n {
            x = f.sqr(x);
        }
        x
    };
    // x_k = a^(2^k − 1).
    let x1 = a;
    let x2 = f.mul(f.sqr(x1), x1);
    let x3 = f.mul(f.sqr(x2), x1);
    let x6 = f.mul(sqn(x3, 3), x3);
    let x12 = f.mul(sqn(x6, 6), x6);
    let x24 = f.mul(sqn(x12, 12), x12);
    let x48 = f.mul(sqn(x24, 24), x24);
    let x60 = f.mul(sqn(x48, 12), x12);
    let x62 = f.mul(sqn(x60, 2), x2);
    let x124 = f.mul(sqn(x62, 62), x62);
    let x127 = f.mul(sqn(x124, 3), x3);
    f.mul(sqn(x127, 63), x62)
}

/// The P-192 prime as limbs.
const P192_LIMBS: [u64; 3] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFE,
    0xFFFF_FFFF_FFFF_FFFF,
];

/// 62 "posdivsteps" on the low 64 bits of `f` (odd) and `g`, tracking the
/// Jacobi symbol `(g/f)` in bit 0 of `jac` (Wuille's non-negative variant
/// of the Bernstein–Yang divstep, as in libsecp256k1's `jacobi64`).
/// Returns `(u, v, q, r)` with `2^62·f' = u·f + v·g`, `2^62·g' = q·f + r·g`,
/// all entries non-negative. Every decision reads only the low bits, which
/// are exact.
#[inline]
fn posdivsteps_62(eta: &mut i64, mut f: u64, mut g: u64, jac: &mut u32) -> [u64; 4] {
    let (mut u, mut v, mut q, mut r) = (1u64, 0u64, 0u64, 1u64);
    let mut i: i64 = 62;
    loop {
        // Count zeros only up to i (sentinel bit).
        let zeros = (g | (u64::MAX << i)).trailing_zeros();
        g >>= zeros;
        u <<= zeros;
        v <<= zeros;
        *eta -= i64::from(zeros);
        i -= i64::from(zeros);
        // Dividing g by 2 multiplies (g/f) by (2/f) = −1 iff f ≡ 3, 5 (mod 8).
        *jac ^= zeros & ((f >> 1) ^ (f >> 2)) as u32;
        if i == 0 {
            break;
        }
        // g is odd here. If eta < 0, negate it and swap (f, g), (u, q),
        // (v, r) — branch-free, since the branch is unpredictable.
        // Reciprocity: both odd and positive; flip iff both ≡ 3 (mod 4).
        let mask = (*eta >> 63) as u64;
        *jac ^= (((f & g) >> 1) & mask) as u32;
        *eta = (*eta ^ mask as i64) - mask as i64;
        let x = (f ^ g) & mask;
        f ^= x;
        g ^= x;
        let x = (u ^ q) & mask;
        u ^= x;
        q ^= x;
        let x = (v ^ r) & mask;
        v ^= x;
        r ^= x;
        // Cancel the bottom min(eta + 1, i, 4) bits of g with a multiple of f.
        let limit = (*eta + 1).min(i).min(4) as u32;
        let m = (1u64 << limit) - 1;
        // −f⁻¹ mod 16 from a table indexed by (f mod 16) >> 1.
        const NEG_INV16: [u64; 8] = [15, 5, 3, 9, 7, 13, 11, 1];
        let w = g.wrapping_mul(NEG_INV16[((f >> 1) & 7) as usize]) & m;
        g = g.wrapping_add(f.wrapping_mul(w));
        q = q.wrapping_add(u.wrapping_mul(w));
        r = r.wrapping_add(v.wrapping_mul(w));
        debug_assert_eq!(g & m, 0);
        if cfg!(test) {
            assert!(
                u.max(v).max(q).max(r) < 1 << 63,
                "posdivstep matrix overflow"
            );
        }
    }
    [u, v, q, r]
}

/// `(c0·x + c1·y) >> 62` for 3-limb `x, y` and `c0, c1 < 2^63`; the caller
/// guarantees the result fits 3 limbs and the low 62 bits are zero.
#[inline(always)]
fn lincomb_shr62(x: &[u64; 3], c0: u64, y: &[u64; 3], c1: u64) -> [u64; 3] {
    let mut out = [0u64; 4];
    let mut carry: u128 = 0;
    for k in 0..3 {
        let t = u128::from(x[k]) * u128::from(c0) + u128::from(y[k]) * u128::from(c1) + carry;
        out[k] = t as u64;
        carry = t >> 64;
    }
    out[3] = carry as u64;
    debug_assert!(carry >> 64 == 0);
    debug_assert_eq!(out[0] & ((1u64 << 62) - 1), 0);
    [
        (out[0] >> 62) | (out[1] << 2),
        (out[1] >> 62) | (out[2] << 2),
        (out[2] >> 62) | (out[3] << 2),
    ]
}

/// The Jacobi symbol `(a / p)` for the P-192 prime and canonical `a`, by
/// posdivsteps (62 per round, about 9 rounds). Falls back to the binary
/// algorithm if 40 rounds do not finish (never observed; counted in tests).
pub fn jacobi_p192_fast(a: &[u64; 3]) -> i32 {
    if a.iter().all(|&w| w == 0) {
        return 0;
    }
    let mut f = P192_LIMBS;
    let mut g = *a;
    let mut eta: i64 = -1;
    let mut jac: u32 = 0;
    for _ in 0..40 {
        let [u, v, q, r] = posdivsteps_62(&mut eta, f[0], g[0], &mut jac);
        let nf = lincomb_shr62(&f, u, &g, v);
        let ng = lincomb_shr62(&f, q, &g, r);
        f = nf;
        g = ng;
        if f == [1, 0, 0] {
            return 1 - 2 * (jac & 1) as i32;
        }
        if g == [0, 0, 0] {
            // gcd(a, p) = f > 1 cannot happen for prime p and 0 < a < p.
            return 0;
        }
    }
    #[cfg(test)]
    JACOBI_FALLBACKS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    jacobi_limbs(a, &P192_LIMBS)
}

#[cfg(test)]
static JACOBI_FALLBACKS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

// ── Draws ─────────────────────────────────────────────────────────────

/// SplitMix64, the seed → draw derivation (frozen in the study protocol).
#[derive(Clone, Copy, Debug)]
pub struct SplitMix(pub u64);

impl SplitMix {
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// A field element from a seed: three SplitMix64 words, little-endian, as
/// a 192-bit integer reduced mod `p` (the bias is below `2^-60` at P-192).
pub fn element_from_rng<F: PrimeField>(f: &F, rng: &mut SplitMix) -> F::Elem {
    let w = [rng.next_u64(), rng.next_u64(), rng.next_u64()];
    let mut x = BigUint::zero();
    for &wi in w.iter().rev() {
        x = (x << 64u32) + BigUint::from(wi);
    }
    f.from_biguint(&x)
}

/// The domain-separated generator for draw `seed` of `family`.
pub fn draw_rng(family: Family, seed: u64) -> SplitMix {
    let tag = match family {
        Family::Z2xZ2 => 0x5A32_5A32_0000_0001u64,
        Family::Z2xZ8 => 0x5A32_5A38_0000_0001u64,
    };
    let mut s = SplitMix(seed ^ tag);
    s.next_u64();
    s
}

/// The curve families of the search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Family {
    Z2xZ2,
    Z2xZ8,
}

impl Family {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "z2xz2" => Some(Family::Z2xZ2),
            "z2xz8" => Some(Family::Z2xZ8),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Family::Z2xZ2 => "z2xz2",
            Family::Z2xZ8 => "z2xz8",
        }
    }
    /// The 2-depth every non-singular member has.
    pub fn forced_depth(self) -> u32 {
        match self {
            Family::Z2xZ2 => 1,
            Family::Z2xZ8 => 3,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Arm {
    /// Three square roots per level, halves tested one by one.
    Baseline,
    /// One Legendre symbol per level test, two exponentiations per halving.
    FastChi,
    /// The test is the next level's exponentiation (no Legendre symbol
    /// below the family's start).
    Fast,
}

/// A curve `y² = (x − e₁)(x − e₂)(x − e₃)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve3<E> {
    pub e: [E; 3],
}

/// A point of order `2^depth` on the drawn curve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeepPoint<E> {
    pub x: E,
    pub y: E,
    pub depth: u32,
}

/// The outcome of one draw.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DrawResult<E> {
    pub curve: Curve3<E>,
    /// `b` (the 2-Sylow is `Z/2 × Z/2^b`); `0` for a singular draw.
    pub depth: u32,
    /// The point of order exactly `2^record_t`, when `depth ≥ record_t`.
    pub point: Option<DeepPoint<E>>,
}

/// Operation counts of one search (all draws).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HalvingCost {
    pub draws: u64,
    pub singular: u64,
    pub legendre: u64,
    /// Exponentiations: `(p−3)/4` powers (fast) or square roots (baseline).
    pub exps: u64,
    /// Successful halvings (levels descended past the family's start).
    pub halvings: u64,
}

/// The search context: the field, its exponent and the arm.
pub struct Halver<F: HalvingField> {
    pub f: F,
    pub arm: Arm,
    isqrt_exp: Vec<u64>,
    pub cost: HalvingCost,
}

/// Per-curve constants of the fast test.
struct Plan {
    /// The index of the 2-torsion point in `2E`.
    h: usize,
    /// The component `k` with `d_k = +1`.
    k: usize,
    /// Another component, `≠ k`: `δ(H) = 1` iff it is `+1` (given the
    /// `k` test passed).
    other: usize,
    /// A 2-torsion index not in `2E` (the sign to flip for `H + T_n`).
    n: usize,
}

impl<F: HalvingField> Halver<F> {
    pub fn new(f: F, arm: Arm) -> Self {
        let p = f.modulus().clone();
        assert!(
            (&p % 4u32) == BigUint::from(3u32),
            "the halving search needs p ≡ 3 (mod 4)"
        );
        let e = (&p - 3u32) >> 2u32;
        Halver {
            f,
            arm,
            isqrt_exp: e.to_u64_digits(),
            cost: HalvingCost::default(),
        }
    }

    #[inline]
    fn chi(&mut self, a: F::Elem) -> i32 {
        self.cost.legendre += 1;
        match self.arm {
            Arm::Fast | Arm::FastChi => self.f.chi(a),
            Arm::Baseline => self.f.legendre(a),
        }
    }

    #[inline]
    fn isqrt(&mut self, a: F::Elem) -> F::Elem {
        self.cost.exps += 1;
        self.f.pow_isqrt(a, &self.isqrt_exp)
    }

    #[inline]
    fn sqrt(&mut self, a: F::Elem) -> Option<F::Elem> {
        self.cost.exps += 1;
        self.f.sqrt(a)
    }

    /// The curve of draw `seed`, and for `Z2xZ8` the order-4 point with its
    /// square roots. `None` when singular.
    #[allow(clippy::type_complexity)]
    pub fn draw_curve(
        &self,
        family: Family,
        seed: u64,
    ) -> Option<(Curve3<F::Elem>, Option<(F::Elem, F::Elem, [F::Elem; 3])>)> {
        let f = &self.f;
        let mut rng = draw_rng(family, seed);
        match family {
            Family::Z2xZ2 => {
                let al = element_from_rng(f, &mut rng);
                let be = element_from_rng(f, &mut rng);
                let z = f.zero();
                if al == z || be == z || al == be {
                    return None;
                }
                Some((Curve3 { e: [z, al, be] }, None))
            }
            Family::Z2xZ8 => {
                let m = element_from_rng(f, &mut rng);
                let m2 = f.sqr(m);
                let a = f.sub(m2, f.one());
                let b = f.double(m);
                let c = f.add(m2, f.one());
                let (a2, b2) = (f.sqr(a), f.sqr(b));
                let (a4, b4) = (f.sqr(a2), f.sqr(b2));
                // Singular iff A·B = 0 or A⁴ = B⁴ (C ≠ 0 since −1 is a non-residue).
                if f.is_zero(a) || f.is_zero(b) || a4 == b4 {
                    return None;
                }
                let curve = Curve3 {
                    e: [f.zero(), f.neg(b4), f.neg(a4)],
                };
                let x4 = f.mul(a2, b2);
                let r = [f.mul(a, b), f.mul(b, c), f.mul(a, c)];
                let y4 = f.mul(x4, f.sqr(c));
                Some((curve, Some((x4, y4, r))))
            }
        }
    }

    /// `b` for draw `seed`, and a point of order `2^record_t` if `b ≥
    /// record_t` (`record_t = 0`: none).
    pub fn run_draw(&mut self, family: Family, seed: u64, record_t: u32) -> DrawResult<F::Elem> {
        self.cost.draws += 1;
        let Some((curve, start)) = self.draw_curve(family, seed) else {
            self.cost.singular += 1;
            let z = self.f.zero();
            return DrawResult {
                curve: Curve3 { e: [z; 3] },
                depth: 0,
                point: None,
            };
        };
        let (depth, point) = match self.arm {
            Arm::Fast | Arm::FastChi => self.descend_fast(&curve, start, record_t),
            Arm::Baseline => self.descend_baseline(&curve, start, record_t),
        };
        DrawResult {
            curve,
            depth,
            point,
        }
    }

    /// `x_H` and `y_H` for the signed roots `r`, from `x`.
    #[inline]
    fn half(&self, x: F::Elem, r: &[F::Elem; 3]) -> (F::Elem, F::Elem) {
        let f = &self.f;
        let s01 = f.add(r[0], r[1]);
        let s02 = f.add(r[0], r[2]);
        let s12 = f.add(r[1], r[2]);
        // x + r0r1 + r0r2 + r1r2 = x + r0(r1 + r2) + r1r2.
        let xh = f.add(f.add(x, f.mul(r[0], s12)), f.mul(r[1], r[2]));
        let yh = f.mul(f.mul(s01, s02), s12);
        (xh, yh)
    }

    /// `x_H − e_i = (r_i + r_j)(r_i + r_k)`.
    #[inline]
    fn half_minus_e(&self, r: &[F::Elem; 3], i: usize) -> F::Elem {
        let f = &self.f;
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        f.mul(f.add(r[i], r[j]), f.add(r[i], r[k]))
    }

    /// The plan, from the three difference symbols; `None` if `b = 1`.
    fn plan(&mut self, c: &Curve3<F::Elem>) -> Option<Plan> {
        let f = &self.f;
        // s[i][j] = χ(e_i − e_j); χ(−1) = −1.
        let d01 = f.sub(c.e[0], c.e[1]);
        let d02 = f.sub(c.e[0], c.e[2]);
        let d12 = f.sub(c.e[1], c.e[2]);
        let c01 = self.chi(d01);
        let c02 = self.chi(d02);
        let c12 = self.chi(d12);
        let s = |i: usize, j: usize| -> i32 {
            match (i, j) {
                (0, 1) => c01,
                (1, 0) => -c01,
                (0, 2) => c02,
                (2, 0) => -c02,
                (1, 2) => c12,
                (2, 1) => -c12,
                _ => unreachable!(),
            }
        };
        let h = (0..3).find(|&i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            s(i, j) == 1 && s(i, k) == 1
        })?;
        Some(self.plan_from(h, &s))
    }

    fn plan_from(&self, h: usize, s: &dyn Fn(usize, usize) -> i32) -> Plan {
        let n = (h + 1) % 3;
        // d = δ(T_n): components χ(e_n − e_i) for i ≠ n, product at n.
        let mut d = [0i32; 3];
        for (i, di) in d.iter_mut().enumerate() {
            if i != n {
                *di = s(n, i);
            }
        }
        let (a, b) = ((n + 1) % 3, (n + 2) % 3);
        d[n] = d[a] * d[b];
        let k = (0..3).find(|&i| d[i] == 1).expect("d has one +1 component");
        let other = (k + 1) % 3;
        Plan { h, k, other, n }
    }

    /// Fast descent. Returns `(b, point of order 2^record_t)`.
    #[allow(clippy::type_complexity)]
    fn descend_fast(
        &mut self,
        c: &Curve3<F::Elem>,
        start: Option<(F::Elem, F::Elem, [F::Elem; 3])>,
        record_t: u32,
    ) -> (u32, Option<DeepPoint<F::Elem>>) {
        let f = self.f.clone();
        // Level j: a halvable point (x, y) of order 2^j with roots r.
        let (plan, mut j, mut x, mut y, mut r) = match start {
            Some((x4, y4, r4)) => {
                // Z2xZ8: T₁ = (0,0) ∈ 2E; e₂ − e₁ = −B⁴ and e₃ − e₁ = −A⁴
                // are non-residues, so only χ(e₂ − e₃) is needed.
                let d12 = f.sub(c.e[1], c.e[2]);
                let c12 = self.chi(d12);
                let s = move |i: usize, j: usize| -> i32 {
                    match (i, j) {
                        (0, _) => 1,
                        (_, 0) => -1,
                        (1, 2) => c12,
                        (2, 1) => -c12,
                        _ => unreachable!(),
                    }
                };
                (self.plan_from(0, &s), 2u32, x4, y4, r4)
            }
            None => {
                let Some(plan) = self.plan(c) else {
                    return (1, two_torsion_record(c, f.zero(), record_t));
                };
                // P = T_h (order 2): r_h = 0, r_i = √(e_h − e_i).
                let h = plan.h;
                let mut r = [f.zero(); 3];
                for i in 0..3 {
                    if i != h {
                        let a = f.sub(c.e[h], c.e[i]);
                        let is = self.isqrt(a);
                        r[i] = f.mul(a, is);
                        debug_assert_eq!(f.sqr(r[i]), a);
                    }
                }
                (plan, 1u32, c.e[h], f.zero(), r)
            }
        };
        let mut rec = None;
        if record_t == j {
            rec = Some(DeepPoint { x, y, depth: j });
        }
        let one = f.one();
        let (k, o, m) = (plan.k, plan.other, (plan.k + 2) % 3);
        loop {
            if self.arm == Arm::FastChi {
                // v1: one Legendre symbol decides whether some half is
                // halvable, a second which one; then two exponentiations.
                let hk = self.half_minus_e(&r, k);
                if self.chi(hk) != 1 {
                    return (j + 1, rec.or_else(|| self.record_half(x, &r, j, record_t)));
                }
                let ho = self.half_minus_e(&r, o);
                if self.chi(ho) != 1 {
                    r[plan.n] = f.neg(r[plan.n]);
                }
                let (xh, yh) = self.half(x, &r);
                x = xh;
                y = yh;
                j += 1;
                self.cost.halvings += 1;
                if record_t == j {
                    rec = Some(DeepPoint { x, y, depth: j });
                }
                let a1 = f.sub(x, c.e[0]);
                let a2 = f.sub(x, c.e[1]);
                let i1 = self.isqrt(a1);
                let r0 = f.mul(a1, i1);
                let u = self.isqrt(f.mul(a1, a2));
                let r1 = f.mul(f.mul(a2, u), r0);
                let r2 = f.mul(y, u);
                debug_assert_eq!(f.sqr(r0), a1);
                debug_assert_eq!(f.sqr(r1), a2);
                debug_assert_eq!(f.sqr(r2), f.sub(x, c.e[2]));
                r = [r0, r1, r2];
                continue;
            }
            // v2: the exponentiation the next level needs is the test.
            // i_k = (x_H − e_k)^((p−3)/4); x_H − e_k is a square (some
            // half is halvable) iff i_k²·(x_H − e_k) = 1.
            let mut hk = self.half_minus_e(&r, k);
            let mut ik = self.isqrt(hk);
            if f.mul(f.sqr(ik), hk) != one {
                return (j + 1, rec.or_else(|| self.record_half(x, &r, j, record_t)));
            }
            // u = ((x_H − e_k)(x_H − e_o))^((p−3)/4) decides H or H + T_n.
            let mut ho = self.half_minus_e(&r, o);
            let mut u = self.isqrt(f.mul(hk, ho));
            if f.mul(f.sqr(u), f.mul(hk, ho)) != one {
                r[plan.n] = f.neg(r[plan.n]);
                hk = self.half_minus_e(&r, k);
                ik = self.isqrt(hk);
                ho = self.half_minus_e(&r, o);
                u = self.isqrt(f.mul(hk, ho));
            }
            let (xh, yh) = self.half(x, &r);
            x = xh;
            y = yh;
            j += 1;
            self.cost.halvings += 1;
            if record_t == j {
                rec = Some(DeepPoint { x, y, depth: j });
            }
            let rk = f.mul(hk, ik);
            let ro = f.mul(f.mul(ho, u), rk);
            let rm = f.mul(y, u);
            debug_assert_eq!(f.sqr(rk), f.sub(x, c.e[k]));
            debug_assert_eq!(f.sqr(ro), f.sub(x, c.e[o]));
            debug_assert_eq!(f.sqr(rm), f.sub(x, c.e[m]));
            r[k] = rk;
            r[o] = ro;
            r[m] = rm;
        }
    }

    /// When the descent stops at level `j` (halves exist, none halvable),
    /// a half has order `2^(j+1)`; record it if that is the target.
    fn record_half(
        &self,
        x: F::Elem,
        r: &[F::Elem; 3],
        j: u32,
        record_t: u32,
    ) -> Option<DeepPoint<F::Elem>> {
        if record_t != j + 1 {
            return None;
        }
        let (xh, yh) = self.half(x, r);
        Some(DeepPoint {
            x: xh,
            y: yh,
            depth: j + 1,
        })
    }

    /// Baseline descent: three square roots per level, the four halves
    /// tested in turn with two Legendre symbols each.
    #[allow(clippy::type_complexity)]
    fn descend_baseline(
        &mut self,
        c: &Curve3<F::Elem>,
        start: Option<(F::Elem, F::Elem, [F::Elem; 3])>,
        record_t: u32,
    ) -> (u32, Option<DeepPoint<F::Elem>>) {
        let f = self.f.clone();
        let (mut j, mut x, mut y) = match start {
            Some((x4, y4, _)) => (2u32, x4, y4),
            None => {
                // Find the halvable 2-torsion point, testing each in turn.
                let mut found = None;
                for h in 0..3 {
                    let (a, b) = ((h + 1) % 3, (h + 2) % 3);
                    if self.chi(f.sub(c.e[h], c.e[a])) == 1 && self.chi(f.sub(c.e[h], c.e[b])) == 1
                    {
                        found = Some(h);
                        break;
                    }
                }
                let Some(h) = found else {
                    return (1, two_torsion_record(c, f.zero(), record_t));
                };
                (1u32, c.e[h], f.zero())
            }
        };
        let mut rec = None;
        if record_t == j {
            rec = Some(DeepPoint { x, y, depth: j });
        }
        loop {
            // P = (x, y) is halvable: its roots.
            let mut r = [f.zero(); 3];
            for (i, ri) in r.iter_mut().enumerate() {
                *ri = self
                    .sqrt(f.sub(x, c.e[i]))
                    .expect("a halvable point has square x − e_i");
            }
            // The four halves: flip none, r₀, r₁ or r₂.
            let mut next = None;
            let mut first = None;
            for flip in 0..4 {
                let mut rs = r;
                if flip > 0 {
                    rs[flip - 1] = f.neg(rs[flip - 1]);
                }
                let (xh, yh) = self.half(x, &rs);
                if first.is_none() {
                    first = Some((xh, yh));
                }
                if self.chi(f.sub(xh, c.e[0])) == 1 && self.chi(f.sub(xh, c.e[1])) == 1 {
                    next = Some((xh, yh));
                    break;
                }
            }
            match next {
                None => {
                    if rec.is_none() && record_t == j + 1 {
                        let (xh, yh) = first.unwrap();
                        rec = Some(DeepPoint {
                            x: xh,
                            y: yh,
                            depth: j + 1,
                        });
                    }
                    return (j + 1, rec);
                }
                Some((xh, yh)) => {
                    x = xh;
                    y = yh;
                    j += 1;
                    self.cost.halvings += 1;
                    if record_t == j {
                        rec = Some(DeepPoint { x, y, depth: j });
                    }
                }
            }
        }
    }
}

/// When `b = 1`, any 2-torsion point has the maximal order 2.
fn two_torsion_record<E: Copy>(c: &Curve3<E>, zero: E, record_t: u32) -> Option<DeepPoint<E>> {
    (record_t == 1).then_some(DeepPoint {
        x: c.e[0],
        y: zero,
        depth: 1,
    })
}

// ── Short Weierstrass form and BigUint verification ───────────────────

/// `(a, b)` of `y² = x³ + ax + b` isomorphic to `y² = (x − e₁)(x − e₂)(x −
/// e₃)` by `x ↦ x − s/3` (`s = e₁ + e₂ + e₃`), and that shift.
pub fn short_weierstrass(p: &BigUint, e: [&BigUint; 3]) -> (BigUint, BigUint, BigUint) {
    let m = |a: &BigUint, b: &BigUint| (a * b) % p;
    let inv = |a: &BigUint| a.modpow(&(p - 2u32), p);
    let s1 = (e[0] + e[1] + e[2]) % p;
    let s2 = (m(e[0], e[1]) + m(e[0], e[2]) + m(e[1], e[2])) % p;
    let s3 = m(&m(e[0], e[1]), e[2]);
    let i3 = inv(&BigUint::from(3u32));
    let i27 = inv(&BigUint::from(27u32));
    let s1sq = m(&s1, &s1);
    // a = s2 − s1²/3; b = −2s1³/27 + s1s2/3 − s3.
    let a = (&s2 + p - m(&s1sq, &i3)) % p;
    let t1 = m(&m(&BigUint::from(2u32), &m(&s1sq, &s1)), &i27);
    let t2 = m(&m(&s1, &s2), &i3);
    let b = ((&t2 + p - t1) % p + p - s3) % p;
    let shift = m(&s1, &i3);
    (a, b, shift)
}

/// Affine doubling over `BigUint` on `y² = x³ + ax + b` (`None` = O).
pub fn big_double(
    p: &BigUint,
    a: &BigUint,
    pt: &Option<(BigUint, BigUint)>,
) -> Option<(BigUint, BigUint)> {
    let (x, y) = pt.as_ref()?;
    if y.is_zero() {
        return None;
    }
    let inv = |v: &BigUint| v.modpow(&(p - 2u32), p);
    let num = (BigUint::from(3u32) * x * x + a) % p;
    let l = (num * inv(&((BigUint::from(2u32) * y) % p))) % p;
    let x3 = ((&l * &l) % p + 2u32 * p - (2u32 * x) % p) % p;
    let y3 = ((&l * ((x + p - &x3) % p)) % p + p - y) % p;
    Some((x3, y3))
}

/// An artifact for `backlog_check_aux2t`: `(a, b, Qx, Qy)` on the short
/// Weierstrass model, after checking over `BigUint` that `Q` is on it,
/// `[2^(t−1)]Q ≠ O` and `[2^t]Q = O`.
pub fn certify(
    p: &BigUint,
    e: [&BigUint; 3],
    qx: &BigUint,
    qy: &BigUint,
    t: u32,
) -> Result<(BigUint, BigUint, BigUint, BigUint), String> {
    let (a, b, shift) = short_weierstrass(p, e);
    let x = (qx + p - &shift) % p;
    let y = qy % p;
    let disc = (BigUint::from(4u32) * &a * &a * &a + BigUint::from(27u32) * &b * &b) % p;
    if disc.is_zero() {
        return Err("singular".into());
    }
    if (&y * &y) % p != (&x * &x * &x + &a * &x + &b) % p {
        return Err("Q is not on the curve".into());
    }
    let mut q = Some((x.clone(), y.clone()));
    for i in 1..t {
        q = big_double(p, &a, &q);
        if q.is_none() {
            return Err(format!("[2^{i}]Q = O"));
        }
    }
    if big_double(p, &a, &q).is_some() {
        return Err(format!("[2^{t}]Q != O"));
    }
    Ok((a, b, x, y))
}

/// The P-192 prime.
pub fn p192() -> BigUint {
    (BigUint::one() << 192u32) - (BigUint::one() << 64u32) - BigUint::one()
}

// ── The tower on a found curve (generic short Weierstrass) ────────────

/// An affine point on `y² = x³ + ax + b`; `None` is O.
pub type Aff<E> = Option<(E, E)>;

/// Affine group law over a [`PrimeField`].
#[derive(Clone, Debug)]
pub struct Sw<F: PrimeField> {
    pub f: F,
    pub a: F::Elem,
    pub b: F::Elem,
}

impl<F: PrimeField> Sw<F> {
    pub fn on_curve(&self, pt: &Aff<F::Elem>) -> bool {
        match pt {
            None => true,
            Some((x, y)) => {
                let f = &self.f;
                let rhs = f.add(f.mul(f.add(f.sqr(*x), self.a), *x), self.b);
                f.sqr(*y) == rhs
            }
        }
    }
    pub fn add(&self, p1: &Aff<F::Elem>, p2: &Aff<F::Elem>) -> Aff<F::Elem> {
        let f = &self.f;
        let (Some((x1, y1)), Some((x2, y2))) = (p1, p2) else {
            return if p1.is_none() { *p2 } else { *p1 };
        };
        let l = if x1 == x2 {
            if f.is_zero(f.add(*y1, *y2)) {
                return None;
            }
            let num = f.add(f.mul(f.from_u64(3), f.sqr(*x1)), self.a);
            f.mul(num, f.inv(f.double(*y1))?)
        } else {
            f.mul(f.sub(*y2, *y1), f.inv(f.sub(*x2, *x1))?)
        };
        let x3 = f.sub(f.sub(f.sqr(l), *x1), *x2);
        let y3 = f.sub(f.mul(l, f.sub(*x1, x3)), *y1);
        Some((x3, y3))
    }
    pub fn double(&self, p: &Aff<F::Elem>) -> Aff<F::Elem> {
        self.add(p, p)
    }
    pub fn rhs(&self, x: F::Elem) -> F::Elem {
        let f = &self.f;
        f.add(f.mul(f.add(f.sqr(x), self.a), x), self.b)
    }
}

/// One Vélu 2-isogeny step: kernel abscissa `ξ`, `τ = 3ξ² + a`, and the
/// domain's coefficients.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VeluStepF<E> {
    pub xi: E,
    pub tau: E,
    pub a: E,
    pub b: E,
}

/// `V = x(R' + ⟨Q⟩)`, the Vélu chain and `c = x(Φ(R'))`.
#[derive(Clone, Debug)]
pub struct TowerF<E> {
    pub t: u32,
    pub v: Vec<E>,
    pub steps: Vec<VeluStepF<E>>,
    pub c: E,
    pub r: (E, E),
}

/// The composed x-map `L = L_t ∘ … ∘ L_1` at `x` (`None` at a pole).
pub fn tower_eval<F: PrimeField>(
    f: &F,
    steps: &[VeluStepF<F::Elem>],
    x: F::Elem,
) -> Option<F::Elem> {
    let mut y = x;
    for s in steps {
        let d = f.sub(y, s.xi);
        y = f.add(y, f.mul(s.tau, f.inv(d)?));
    }
    Some(y)
}

/// Build the tower on `E'` with `Q` of order `2^t`; `R'` is drawn from
/// `seed` with `2R' ∉ ⟨Q⟩`. Fails loudly on any broken invariant.
pub fn build_tower<F: PrimeField>(
    e: &Sw<F>,
    q: (F::Elem, F::Elem),
    t: u32,
    seed: u64,
) -> Result<TowerF<F::Elem>, String> {
    let f = &e.f;
    let qa = Some(q);
    if !e.on_curve(&qa) {
        return Err("Q is not on E'".into());
    }
    let n = 1usize << t;
    let mut multiples: Vec<Aff<F::Elem>> = Vec::with_capacity(n);
    multiples.push(None);
    let mut acc = qa;
    for _ in 1..n {
        if acc.is_none() {
            return Err("Q has order below 2^t".into());
        }
        multiples.push(acc);
        acc = e.add(&acc, &qa);
    }
    if acc.is_some() {
        return Err("Q does not have order 2^t".into());
    }
    let keys: std::collections::HashSet<Option<F::Elem>> =
        multiples.iter().map(|m| m.map(|(x, _)| x)).collect();
    let mut rng = SplitMix(seed ^ 0x7E57_0000_0000_0001);
    let r0 = loop {
        let x = element_from_rng(f, &mut rng);
        let Some(y) = f.sqrt(e.rhs(x)) else { continue };
        let r = Some((x, y));
        let r2 = e.double(&r);
        if !keys.contains(&r2.map(|(x, _)| x)) {
            break (x, y);
        }
    };
    let ra = Some(r0);
    let mut v = Vec::with_capacity(n);
    for m in &multiples {
        let s = e.add(&ra, m);
        let (x, _) = s.ok_or("R' + kQ = O")?;
        v.push(x);
    }
    let (mut aj, mut bj) = (e.a, e.b);
    let (mut qj, mut rj) = (q, r0);
    let mut steps = Vec::with_capacity(t as usize);
    let three = f.from_u64(3);
    let five = f.from_u64(5);
    let seven = f.from_u64(7);
    for j in 0..t {
        let cj = Sw {
            f: f.clone(),
            a: aj,
            b: bj,
        };
        let mut kj = Some(qj);
        for _ in 0..(t - 1 - j) {
            kj = cj.double(&kj);
        }
        let (xi, yk) = kj.ok_or(format!("step {j}: kernel point is O"))?;
        if !f.is_zero(yk) {
            return Err(format!("step {j}: the kernel point is not 2-torsion"));
        }
        let tau = f.add(f.mul(three, f.sqr(xi)), aj);
        steps.push(VeluStepF {
            xi,
            tau,
            a: aj,
            b: bj,
        });
        let velu = |pt: (F::Elem, F::Elem)| -> Result<(F::Elem, F::Elem), String> {
            let d = f.inv(f.sub(pt.0, xi)).ok_or("Vélu pole")?;
            let x = f.add(pt.0, f.mul(tau, d));
            let y = f.mul(pt.1, f.sub(f.one(), f.mul(tau, f.sqr(d))));
            Ok((x, y))
        };
        // Q_j has order 2^(t−j) > 2 for j < t − 1, so it is not the kernel;
        // at j = t − 1 the image of Q is O and is not needed.
        if j + 1 < t {
            qj = velu(qj)?;
        }
        rj = velu(rj)?;
        aj = f.sub(aj, f.mul(five, tau));
        bj = f.sub(bj, f.mul(seven, f.mul(xi, tau)));
        let next = Sw {
            f: f.clone(),
            a: aj,
            b: bj,
        };
        if (j + 1 < t && !next.on_curve(&Some(qj))) || !next.on_curve(&Some(rj)) {
            return Err(format!("step {j}: a Vélu image is off the codomain"));
        }
    }
    Ok(TowerF {
        t,
        v,
        steps,
        c: rj.0,
        r: r0,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn is_prime_u64(n: u64) -> bool {
        n >= 2
            && (2..)
                .take_while(|d| d * d <= n)
                .all(|d| !n.is_multiple_of(d))
    }

    /// The exponent `b` of the 2-Sylow `Z/2 × Z/2^b` by exhaustion: the
    /// largest 2-power dividing a point order, over all points, by
    /// repeated addition in naive u64 arithmetic.
    fn exhaustive_b(p: u64, e: [u64; 3]) -> (u32, u64) {
        let pw = |mut b: u64, mut x: u64| {
            let mut r = 1u64;
            b %= p;
            while x > 0 {
                if x & 1 == 1 {
                    r = r * b % p;
                }
                b = b * b % p;
                x >>= 1;
            }
            r
        };
        let inv = |v: u64| pw(v, p - 2);
        let rhs = |x: u64| (x + p - e[0]) % p * ((x + p - e[1]) % p) % p * ((x + p - e[2]) % p) % p;
        // Coefficients of x³ + c2 x² + c1 x + c0.
        let c2 = (3 * p - e[0] - e[1] - e[2]) % p;
        let c1 = (e[0] * e[1] % p + e[0] * e[2] % p + e[1] * e[2] % p) % p;
        let add = |a: Option<(u64, u64)>, b: Option<(u64, u64)>| -> Option<(u64, u64)> {
            let (Some((x1, y1)), Some((x2, y2))) = (a, b) else {
                return a.or(b);
            };
            let l = if x1 == x2 {
                if (y1 + y2) % p == 0 {
                    return None;
                }
                (3 * x1 % p * x1 % p + 2 * c2 % p * x1 % p + c1) % p * inv(2 * y1 % p) % p
            } else {
                (y2 + p - y1) % p * inv((x2 + p - x1) % p) % p
            };
            let x3 = (l * l % p + 3 * p - c2 - x1 - x2) % p;
            let y3 = (l * ((x1 + p - x3) % p) % p + p - y1) % p;
            Some((x3, y3))
        };
        let mut pts = vec![];
        for x in 0..p {
            let r = rhs(x);
            for y in 0..p {
                if y * y % p == r {
                    pts.push((x, y));
                }
            }
        }
        let n = pts.len() as u64 + 1;
        let mut best = 0u32;
        for &pt in &pts {
            let mut acc = Some(pt);
            let mut k = 1u64;
            while acc.is_some() {
                acc = add(acc, Some(pt));
                k += 1;
            }
            best = best.max(k.trailing_zeros());
        }
        (best, n)
    }

    #[test]
    fn halving_depth_matches_exhaustive_group_structure() {
        // Primes ≡ 3 mod 4 with v₂(p + 1) small and large.
        let primes: Vec<u64> = (100u64..360)
            .filter(|&p| p % 4 == 3 && is_prime_u64(p))
            .collect();
        let mut checked = 0;
        let mut depths = std::collections::BTreeMap::new();
        for &p in &primes {
            let f = Fp64::new(p);
            let mut fast = Halver::new(f.clone(), Arm::Fast);
            let mut base = Halver::new(f.clone(), Arm::Baseline);
            let mut v1 = Halver::new(f.clone(), Arm::FastChi);
            let mut rng = SplitMix(p);
            for _ in 0..12 {
                let al = 1 + rng.next_u64() % (p - 1);
                let be = 1 + rng.next_u64() % (p - 1);
                if al == be {
                    continue;
                }
                let e = [0, al, be];
                let (b_true, n) = exhaustive_b(p, e);
                assert_eq!(n.trailing_zeros(), b_true + 1, "v₂(#E) = b + 1");
                let curve = Curve3 {
                    e: [f.from_u64(0), f.from_u64(al), f.from_u64(be)],
                };
                for b_rec in [b_true, 0] {
                    let (bf, pf) = fast.descend_fast(&curve, None, b_rec);
                    let (bb, pb) = base.descend_baseline(&curve, None, b_rec);
                    let (b1, p1) = v1.descend_fast(&curve, None, b_rec);
                    assert_eq!(bf, b_true, "fast p={p} e={e:?}");
                    assert_eq!(bb, b_true, "baseline p={p} e={e:?}");
                    assert_eq!(b1, b_true, "fast-chi p={p} e={e:?}");
                    if b_rec > 0 {
                        for pt in [pf.unwrap(), pb.unwrap(), p1.unwrap()] {
                            assert_eq!(pt.depth, b_true);
                            check_order_u64(p, e, pt, b_true);
                        }
                    }
                }
                *depths.entry(b_true).or_insert(0u32) += 1;
                checked += 1;
            }
        }
        assert!(checked > 250, "{checked}");
        // The test population reaches several depths.
        assert!(depths.keys().max().copied().unwrap_or(0) >= 5, "{depths:?}");
    }

    /// The recorded point has order exactly 2^t (u64 arithmetic).
    fn check_order_u64(p: u64, e: [u64; 3], pt: DeepPoint<u64>, t: u32) {
        let f = Fp64::new(p);
        let (pb, eb) = (BigUint::from(p), e.map(BigUint::from));
        let x = BigUint::from(f.to_u64(pt.x).unwrap());
        let y = BigUint::from(f.to_u64(pt.y).unwrap());
        certify(&pb, [&eb[0], &eb[1], &eb[2]], &x, &y, t).unwrap();
    }

    #[test]
    fn z2xz8_family_forces_depth_three_and_agrees_with_exhaustion() {
        let mut n = 0;
        for p in [103u64, 199, 263, 307, 383, 479, 503] {
            assert!(is_prime_u64(p) && p % 4 == 3);
            let f = Fp64::new(p);
            let mut fast = Halver::new(f.clone(), Arm::Fast);
            let mut base = Halver::new(f.clone(), Arm::Baseline);
            let mut v1 = Halver::new(f.clone(), Arm::FastChi);
            for seed in 0..25 {
                let Some((c, start)) = fast.draw_curve(Family::Z2xZ8, seed) else {
                    continue;
                };
                let e = c.e.map(|v| f.to_u64(v).unwrap());
                let (b_true, _) = exhaustive_b(p, e);
                assert!(b_true >= 3, "the family forces b ≥ 3");
                let (x4, y4, r4) = start.unwrap();
                // H₄ is on the curve with the stated roots.
                let rhs = f.mul(
                    f.mul(f.sub(x4, c.e[0]), f.sub(x4, c.e[1])),
                    f.sub(x4, c.e[2]),
                );
                assert_eq!(f.sqr(y4), rhs);
                for i in 0..3 {
                    assert_eq!(f.sqr(r4[i]), f.sub(x4, c.e[i]));
                }
                let df = fast.run_draw(Family::Z2xZ8, seed, b_true);
                let db = base.run_draw(Family::Z2xZ8, seed, b_true);
                let d1 = v1.run_draw(Family::Z2xZ8, seed, b_true);
                assert_eq!(df.depth, b_true);
                assert_eq!(db.depth, b_true);
                assert_eq!(d1.depth, b_true);
                check_order_u64(p, e, d1.point.unwrap(), b_true);
                check_order_u64(p, e, df.point.unwrap(), b_true);
                check_order_u64(p, e, db.point.unwrap(), b_true);
                n += 1;
            }
        }
        assert!(n > 100);
    }

    #[test]
    fn fast_jacobi_matches_binary_and_euler_at_p192() {
        let f = FpP192::new();
        let p = p192();
        let mut rng = StdRng::seed_from_u64(5);
        let mut vals: Vec<BigUint> = vec![
            BigUint::zero(),
            BigUint::one(),
            BigUint::from(2u32),
            BigUint::from(3u32),
            &p - 1u32,
            &p - 2u32,
            &p >> 1u32,
            (&p >> 1u32) + 1u32,
        ];
        for k in [62u32, 63, 64, 65, 126, 127, 128, 129, 190, 191] {
            vals.push(BigUint::one() << k);
            vals.push((BigUint::one() << k) - 1u32);
            vals.push((BigUint::one() << k) + 1u32);
        }
        for _ in 0..20_000 {
            vals.push(rng.gen_biguint_below(&p));
        }
        for v in &vals {
            let v = v % &p;
            let a = f.from_biguint(&v);
            let fast = jacobi_p192_fast(&a);
            assert_eq!(fast, f.legendre(a), "{v}");
            if !v.is_zero() {
                let euler = v.modpow(&((&p - 1u32) >> 1u32), &p);
                assert_eq!(fast == 1, euler.is_one(), "{v}");
            }
        }
        assert_eq!(
            JACOBI_FALLBACKS.load(std::sync::atomic::Ordering::Relaxed),
            0
        );
    }

    #[test]
    fn isqrt_chain_matches_biguint_at_p192() {
        let f = FpP192::new();
        let p = p192();
        let e = (&p - 3u32) >> 2u32;
        let mut rng = StdRng::seed_from_u64(6);
        let mut vals = vec![BigUint::one(), BigUint::from(2u32), &p - 1u32, &p - 2u32];
        for _ in 0..500 {
            vals.push(rng.gen_biguint_below(&p));
        }
        for v in vals {
            let got = f.to_biguint(p192_pow_isqrt(&f, f.from_biguint(&v)));
            assert_eq!(got, v.modpow(&e, &p));
        }
    }

    /// Positive control at full size: y² = x³ − x is supersingular over
    /// the P-192 prime (p ≡ 3 mod 4), #E = p + 1, v₂(p + 1) = 64, so the
    /// halving must find b = 63 exactly, in both arms.
    #[test]
    fn supersingular_control_has_depth_63_at_p192() {
        let f = FpP192::new();
        let curve = Curve3 {
            e: [f.zero(), f.one(), f.neg(f.one())],
        };
        for arm in [Arm::Fast, Arm::FastChi, Arm::Baseline] {
            let mut h = Halver::new(f.clone(), arm);
            let (b, pt) = match arm {
                Arm::Fast | Arm::FastChi => h.descend_fast(&curve, None, 63),
                Arm::Baseline => h.descend_baseline(&curve, None, 63),
            };
            assert_eq!(b, 63, "{arm:?}");
            let pt = pt.unwrap();
            let p = p192();
            let e = curve.e.map(|v| f.to_biguint(v));
            certify(
                &p,
                [&e[0], &e[1], &e[2]],
                &f.to_biguint(pt.x),
                &f.to_biguint(pt.y),
                63,
            )
            .unwrap();
        }
    }

    #[test]
    fn arms_agree_on_p192_draws_and_found_points_certify() {
        let f = FpP192::new();
        let p = p192();
        for fam in [Family::Z2xZ2, Family::Z2xZ8] {
            let mut fast = Halver::new(f.clone(), Arm::Fast);
            let mut base = Halver::new(f.clone(), Arm::Baseline);
            let mut v1 = Halver::new(f.clone(), Arm::FastChi);
            let mut deep = 0;
            for seed in 0..3000u64 {
                let a = fast.run_draw(fam, seed, 6);
                let b = base.run_draw(fam, seed, 6);
                let c1 = v1.run_draw(fam, seed, 6);
                assert_eq!(a.depth, b.depth, "{fam:?} seed {seed}");
                assert_eq!(a.depth, c1.depth, "{fam:?} seed {seed}");
                if a.depth >= 6 {
                    deep += 1;
                    let e = a.curve.e.map(|v| f.to_biguint(v));
                    for pt in [a.point.unwrap(), b.point.unwrap(), c1.point.unwrap()] {
                        certify(
                            &p,
                            [&e[0], &e[1], &e[2]],
                            &f.to_biguint(pt.x),
                            &f.to_biguint(pt.y),
                            6,
                        )
                        .unwrap();
                    }
                }
            }
            assert!(deep > 10, "{fam:?}: {deep}");
            assert!(fast.cost.exps < base.cost.exps);
        }
    }

    #[test]
    fn tower_on_a_found_p192_curve_has_two_to_the_t_roots() {
        let f = FpP192::new();
        let p = p192();
        let t = 8u32;
        let mut h = Halver::new(f.clone(), Arm::Fast);
        let (curve, pt) = (0..)
            .find_map(|s| {
                let d = h.run_draw(Family::Z2xZ8, s, t);
                d.point.map(|pt| (d.curve, pt))
            })
            .unwrap();
        let e = curve.e.map(|v| f.to_biguint(v));
        let (a, b, qx, qy) = certify(
            &p,
            [&e[0], &e[1], &e[2]],
            &f.to_biguint(pt.x),
            &f.to_biguint(pt.y),
            t,
        )
        .unwrap();
        let sw = Sw {
            f: f.clone(),
            a: f.from_biguint(&a),
            b: f.from_biguint(&b),
        };
        let tower = build_tower(&sw, (f.from_biguint(&qx), f.from_biguint(&qy)), t, 1).unwrap();
        let mut v = tower.v.clone();
        v.sort();
        v.dedup();
        assert_eq!(v.len(), 1 << t, "|V| = 2^t");
        for &x in &tower.v {
            assert_eq!(tower_eval(&f, &tower.steps, x), Some(tower.c));
        }
        // Independently over BigUint: every x in V maps to c.
        let (pp, cc) = (p.clone(), f.to_biguint(tower.c));
        for &x in tower.v.iter().take(64) {
            let mut y = f.to_biguint(x);
            for s in &tower.steps {
                let d = (&y + &pp - f.to_biguint(s.xi)) % &pp;
                let di = d.modpow(&(&pp - 2u32), &pp);
                y = (&y + f.to_biguint(s.tau) * di) % &pp;
            }
            assert_eq!(y, cc);
        }
        // And random abscissae outside V (almost surely) do not.
        let mut rng = SplitMix(9);
        for _ in 0..64 {
            let x = element_from_rng(&f, &mut rng);
            if !tower.v.contains(&x) {
                assert_ne!(tower_eval(&f, &tower.steps, x), Some(tower.c));
            }
        }
    }

    #[test]
    fn small_prime_tower_has_exactly_the_roots_v() {
        // Exhaustive version of the u64 test: L(x) = c has exactly V.
        let p = 10_007u64;
        assert!(p % 4 == 3);
        let f = Fp64::new(p);
        let mut h = Halver::new(f.clone(), Arm::Fast);
        let t = 4;
        let (curve, pt) = (0..)
            .find_map(|s| {
                let d = h.run_draw(Family::Z2xZ8, s, t);
                d.point.map(|pt| (d.curve, pt))
            })
            .unwrap();
        let pb = BigUint::from(p);
        let e = curve.e.map(|v| f.to_biguint(v));
        let (a, b, qx, qy) = certify(
            &pb,
            [&e[0], &e[1], &e[2]],
            &f.to_biguint(pt.x),
            &f.to_biguint(pt.y),
            t,
        )
        .unwrap();
        let sw = Sw {
            f: f.clone(),
            a: f.from_biguint(&a),
            b: f.from_biguint(&b),
        };
        let tower = build_tower(&sw, (f.from_biguint(&qx), f.from_biguint(&qy)), t, 3).unwrap();
        let roots = (0..p)
            .filter(|&x| tower_eval(&f, &tower.steps, f.from_u64(x)) == Some(tower.c))
            .count();
        assert_eq!(roots, 1 << t);
    }
}
