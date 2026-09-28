//! # GLS curves over `F_{p²}` and their `ψ`-invariant line base (type B).
//!
//! Galbraith–Lin–Scott take a curve `E / F_p` and its quadratic twist
//! `E' / F_{p²}`, `E' : y² = x³ + a u² x + b u³` for a non-square
//! `u ∈ F_{p²}`, and conjugate the `p`-power Frobenius through the twist
//! isomorphism `τ_u : (x, y) ↦ (u x, u^{3/2} y)`:
//!
//! ```text
//!     ψ = τ_u ∘ π_p ∘ τ_u⁻¹ :  (x, y) ↦ (u^{1−p} x^p,  u^{3(1−p)/2} y^p),
//!     ψ² = −1 on E'(F_{p²}),   #E'(F_{p²}) = (p − 1)² + t².
//! ```
//!
//! `ψ` is a Frobenius-type endomorphism of degree `p` whose eigenvalue
//! has order **4** modulo `r` (`λ² ≡ −1`), so it folds a factor base
//! exactly as the Koblitz Frobenius does, by `⟨−1, ψ⟩` of order 4.  The
//! base that is invariant by construction is a **line**: writing
//! `F_{p²} = F_p(s)`, `s² = ν`, the abscissa map `x ↦ u^{1−p} x^p` sends
//! `u·s·t ↦ −u·s·t` for `t ∈ F_p` (because `s^p = −s`), so
//!
//! ```text
//!     F = { P ∈ E'(F_{p²}) : x(P) ∈ u·s·F_p }
//! ```
//!
//! is `ψ`-stable with `4` signed points per column.  (The other
//! `ψ`-stable line, `u·F_p`, carries no points: `x = u t` gives
//! `y² = u³(t³ + at + b)`, a non-square times a square.)  Both facts are
//! checked at run time, not assumed: [`gls_line_base`] folds under
//! [`Closure::Strict`], which refuses a base that is not invariant, and
//! [`generate_gls_instance`] verifies `ψ` on random points.
//!
//! This is the same construction as the Koblitz fold in a different
//! field: there the `π`-stable sets are `F_2`-subspaces of `F_{2^n}`
//! that are `π`-stable, here the `ψ`-stable sets are the `F_p`-lines of
//! `F_{p²}` that `x ↦ u^{1−p}x^p` preserves.  A `4`-dimensional GLV+GLS
//! curve (`j = 1728` twisted over `F_{p²}`, FourQ-like) would fold by
//! `⟨ι, ψ⟩` of order 8 through the same [`fold_by_endomorphisms`] call
//! with both generators; that is the next family the plan lists.
//!
//! The group is a [`CountedGroup`], so the framework's oracles, relation
//! loop and reports run on it unchanged.  `p < 2^31` so that a packed
//! `F_{p²}` element fits the framework's `u64` keys.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::glv_invariant_base::{
    addm, factor_u64, fold_by_endomorphisms, mulm, subm, verify_endomorphism, Closure,
    Endomorphism, FoldReport, Negation,
};
use crate::cryptanalysis::ic_boundary::{CountedGroup, FactorBase, GroupOps};
use crate::cryptanalysis::residual_walk::{inv_mod, is_prime_u64, sqrt_mod};

/// `a + b·s` with `s² = ν`, as `[a, b]`.
pub type Fp2El = [u64; 2];

/// `F_{p²} = F_p[s] / (s² − ν)`, `ν` the smallest non-residue.
#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct Fp2 {
    pub p: u64,
    pub nu: u64,
}

fn jacobi(mut a: u64, mut n: u64) -> i32 {
    let mut result = 1i32;
    a %= n;
    while a != 0 {
        while a.is_multiple_of(2) {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a %= n;
    }
    if n == 1 {
        result
    } else {
        0
    }
}

impl Fp2 {
    /// The field for an odd prime `p < 2^31`.
    pub fn new(p: u64) -> Result<Self, String> {
        if !(3..(1 << 31)).contains(&p) || !is_prime_u64(p) {
            return Err(format!("p = {p} must be an odd prime below 2^31"));
        }
        let nu = (2..p)
            .find(|&n| jacobi(n, p) == -1)
            .ok_or("no quadratic non-residue")?;
        Ok(Self { p, nu })
    }
    pub fn zero(&self) -> Fp2El {
        [0, 0]
    }
    pub fn one(&self) -> Fp2El {
        [1, 0]
    }
    /// The generator `s` of `F_{p²}` over `F_p`.
    pub fn s(&self) -> Fp2El {
        [0, 1]
    }
    pub fn from_fp(&self, a: u64) -> Fp2El {
        [a % self.p, 0]
    }
    pub fn is_zero(&self, a: Fp2El) -> bool {
        a == [0, 0]
    }
    pub fn add(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        [addm(a[0], b[0], self.p), addm(a[1], b[1], self.p)]
    }
    pub fn sub(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        [subm(a[0], b[0], self.p), subm(a[1], b[1], self.p)]
    }
    pub fn neg(&self, a: Fp2El) -> Fp2El {
        [subm(0, a[0], self.p), subm(0, a[1], self.p)]
    }
    pub fn mul(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        let p = self.p;
        let a0b0 = mulm(a[0], b[0], p);
        let a1b1 = mulm(a[1], b[1], p);
        let cross = addm(mulm(a[0], b[1], p), mulm(a[1], b[0], p), p);
        [addm(a0b0, mulm(self.nu, a1b1, p), p), cross]
    }
    pub fn scale(&self, a: Fp2El, k: u64) -> Fp2El {
        [mulm(a[0], k, self.p), mulm(a[1], k, self.p)]
    }
    pub fn sqr(&self, a: Fp2El) -> Fp2El {
        self.mul(a, a)
    }
    /// `a₀² − ν a₁²`, the norm to `F_p`.
    pub fn norm(&self, a: Fp2El) -> u64 {
        subm(
            mulm(a[0], a[0], self.p),
            mulm(self.nu, mulm(a[1], a[1], self.p), self.p),
            self.p,
        )
    }
    /// `a^p = a₀ − a₁ s`: the conjugate.
    pub fn frob(&self, a: Fp2El) -> Fp2El {
        [a[0], subm(0, a[1], self.p)]
    }
    pub fn inv(&self, a: Fp2El) -> Option<Fp2El> {
        let n = self.norm(a);
        if n == 0 {
            return None;
        }
        let ni = inv_mod(n, self.p);
        Some(self.scale(self.frob(a), ni))
    }
    pub fn pow(&self, mut a: Fp2El, mut e: u64) -> Fp2El {
        let mut acc = self.one();
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(acc, a);
            }
            a = self.sqr(a);
            e >>= 1;
        }
        acc
    }
    /// Whether `a` is a square in `F_{p²}`: `a^{(p²−1)/2} = 1`.
    pub fn is_square(&self, a: Fp2El) -> bool {
        if self.is_zero(a) {
            return true;
        }
        self.pow(a, (self.p * self.p - 1) / 2) == self.one()
    }
    /// A square root in `F_{p²}`, or `None` for a non-square.
    pub fn sqrt(&self, a: Fp2El) -> Option<Fp2El> {
        let p = self.p;
        if self.is_zero(a) {
            return Some(self.zero());
        }
        if a[1] == 0 {
            return match sqrt_mod(a[0], p) {
                Some(r) => Some([r, 0]),
                None => sqrt_mod(mulm(a[0], inv_mod(self.nu, p), p), p).map(|r| [0, r]),
            };
        }
        // (x₀ + x₁ s)² = a: x₀² − ν x₁² = ±√N and 2 x₀ x₁ = a₁.
        let n = sqrt_mod(self.norm(a), p)?;
        let half = inv_mod(2, p);
        for sign in [n, subm(0, n, p)] {
            let x0_sq = mulm(addm(a[0], sign, p), half, p);
            if let Some(x0) = sqrt_mod(x0_sq, p) {
                if x0 == 0 {
                    continue;
                }
                let x1 = mulm(a[1], inv_mod(mulm(2, x0, p), p), p);
                let root = [x0, x1];
                if self.sqr(root) == a {
                    return Some(root);
                }
            }
        }
        None
    }
    /// `a₀ + a₁ p`: an injective packing into a word.
    pub fn pack(&self, a: Fp2El) -> u64 {
        a[0] + a[1] * self.p
    }
}

/// A point of `E'(F_{p²})`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Fp2Point {
    pub x: Fp2El,
    pub y: Fp2El,
    pub infinity: bool,
}

impl Fp2Point {
    pub const INFINITY: Self = Self {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    };
    pub fn affine(x: Fp2El, y: Fp2El) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }
}

/// `y² = x³ + Ax + B` over `F_{p²}`, every operation charged.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Fp2Curve {
    pub f: Fp2,
    pub a: Fp2El,
    pub b: Fp2El,
}

impl Fp2Curve {
    /// `x³ + Ax + B`.
    pub fn rhs(&self, x: Fp2El) -> Fp2El {
        let f = &self.f;
        f.add(f.add(f.mul(f.sqr(x), x), f.mul(self.a, x)), self.b)
    }
    pub fn is_on_curve(&self, pt: Fp2Point) -> bool {
        pt.infinity || self.f.sqr(pt.y) == self.rhs(pt.x)
    }
    /// The two points above `x`, if any (one if `y = 0`).
    pub fn lift_x(&self, x: Fp2El) -> Vec<Fp2Point> {
        match self.f.sqrt(self.rhs(x)) {
            None => Vec::new(),
            Some(y) if self.f.is_zero(y) => vec![Fp2Point::affine(x, y)],
            Some(y) => vec![Fp2Point::affine(x, y), Fp2Point::affine(x, self.f.neg(y))],
        }
    }
    fn add_raw(&self, a: Fp2Point, b: Fp2Point) -> Fp2Point {
        if a.infinity {
            return b;
        }
        if b.infinity {
            return a;
        }
        let f = &self.f;
        if a.x == b.x {
            if f.is_zero(f.add(a.y, b.y)) {
                return Fp2Point::INFINITY;
            }
            return self.double_raw(a);
        }
        let lambda = f.mul(f.sub(b.y, a.y), f.inv(f.sub(b.x, a.x)).expect("distinct x"));
        let x3 = f.sub(f.sub(f.sqr(lambda), a.x), b.x);
        let y3 = f.sub(f.mul(lambda, f.sub(a.x, x3)), a.y);
        Fp2Point::affine(x3, y3)
    }
    fn double_raw(&self, a: Fp2Point) -> Fp2Point {
        if a.infinity || self.f.is_zero(a.y) {
            return Fp2Point::INFINITY;
        }
        let f = &self.f;
        let num = f.add(f.scale(f.sqr(a.x), 3), self.a);
        let lambda = f.mul(num, f.inv(f.scale(a.y, 2)).expect("y ≠ 0"));
        let x3 = f.sub(f.sqr(lambda), f.scale(a.x, 2));
        let y3 = f.sub(f.mul(lambda, f.sub(a.x, x3)), a.y);
        Fp2Point::affine(x3, y3)
    }
}

impl CountedGroup for Fp2Curve {
    type Elt = Fp2Point;
    fn identity(&self) -> Fp2Point {
        Fp2Point::INFINITY
    }
    fn is_identity(&self, p: &Fp2Point) -> bool {
        p.infinity
    }
    fn add(&self, ops: &mut GroupOps, p: Fp2Point, q: Fp2Point) -> Fp2Point {
        ops.adds += 1;
        self.add_raw(p, q)
    }
    fn double(&self, ops: &mut GroupOps, p: Fp2Point) -> Fp2Point {
        ops.doubles += 1;
        self.double_raw(p)
    }
    fn neg(&self, p: Fp2Point) -> Fp2Point {
        if p.infinity {
            p
        } else {
            Fp2Point::affine(p.x, self.f.neg(p.y))
        }
    }
    fn key(&self, p: &Fp2Point) -> u64 {
        if p.infinity {
            return 0;
        }
        let sign = u64::from(self.f.pack(p.y) > self.f.pack(self.f.neg(p.y)));
        ((self.f.pack(p.x) + 1) << 1) | sign
    }
}

/// `ψ(x, y) = (cx · x^p, cy · y^p)`.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct GlsEndomorphism {
    pub cx: Fp2El,
    pub cy: Fp2El,
    pub eigenvalue: u64,
    pub p: u64,
}

impl Endomorphism<Fp2Curve> for GlsEndomorphism {
    fn name(&self) -> String {
        "gls-psi".into()
    }
    fn degree(&self) -> u64 {
        self.p
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &Fp2Curve, pt: Fp2Point) -> Fp2Point {
        if pt.infinity {
            return pt;
        }
        let f = &g.f;
        Fp2Point::affine(f.mul(self.cx, f.frob(pt.x)), f.mul(self.cy, f.frob(pt.y)))
    }
}

/// A GLS instance: the twist, its certified order, the prime-order
/// subgroup, and the verified `ψ`.
#[derive(Clone, Debug, Serialize)]
pub struct GlsInstance {
    pub name: String,
    pub p: u64,
    /// The base curve `y² = x³ + ax + b` over `F_p` and its trace.
    pub base_a: u64,
    pub base_b: u64,
    pub trace: i64,
    pub u: Fp2El,
    pub curve: Fp2Curve,
    /// `(p − 1)² + t²`.
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Fp2Point,
    pub psi: GlsEndomorphism,
    /// The `ψ`-stable line's direction, `u·s`.
    pub line: Fp2El,
}

fn random_point(curve: &Fp2Curve, rng: &mut StdRng) -> Fp2Point {
    loop {
        let x = [rng.gen_range(0..curve.f.p), rng.gen_range(0..curve.f.p)];
        let pts = curve.lift_x(x);
        if let Some(&pt) = pts.first() {
            if !curve.f.is_zero(pt.y) {
                return pt;
            }
        }
    }
}

/// **A GLS instance with `p` of about `p_bits` bits**, deterministic in
/// `seed`.  The base curve's order is counted in `O(p)`, the twist's
/// order `(p − 1)² + t²` is certified as `[N]P = O` on random points,
/// `r` is its largest prime factor with cofactor at most
/// `max_cofactor`, and `ψ` is verified on random points before it is
/// returned.  `p_bits ≤ 24`.
pub fn generate_gls_instance(
    p_bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<GlsInstance, String> {
    if !(6..=24).contains(&p_bits) {
        return Err(format!("p_bits = {p_bits} outside 6..=24"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x474C_535F_5053_4900 ^ (p_bits as u64));
    for _ in 0..20_000u32 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            if is_prime_u64(c) {
                break c;
            }
        };
        let f = Fp2::new(p)?;
        let a = rng.gen_range(1..p);
        let b = rng.gen_range(1..p);
        let disc = addm(
            mulm(4, mulm(mulm(a, a, p), a, p), p),
            mulm(27, mulm(b, b, p), p),
            p,
        );
        if disc == 0 {
            continue;
        }
        // #E(F_p) by the Legendre symbol at every abscissa.
        let mut count = 1u64;
        for x in 0..p {
            let rhs = addm(addm(mulm(mulm(x, x, p), x, p), mulm(a, x, p), p), b, p);
            count += (1 + jacobi(rhs, p)) as u64;
        }
        let t = p as i64 + 1 - count as i64;
        let order = ((p - 1) as u128 * (p - 1) as u128 + (t as i128 * t as i128) as u128) as u64;
        let factors = factor_u64(order);
        let Some(&(r, _)) = factors.last() else {
            continue;
        };
        let h = order / r;
        if h > max_cofactor || r < 64 || r % 4 != 1 {
            continue;
        }
        // A non-square u ∈ F_{p²} and the twist.
        let u = loop {
            let c = [rng.gen_range(1..p), rng.gen_range(1..p)];
            if !f.is_square(c) {
                break c;
            }
        };
        let u2 = f.sqr(u);
        let u3 = f.mul(u2, u);
        let curve = Fp2Curve {
            f,
            a: f.mul(f.from_fp(a), u2),
            b: f.mul(f.from_fp(b), u3),
        };
        let mut ops = GroupOps::default();
        let pts: Vec<Fp2Point> = (0..3).map(|_| random_point(&curve, &mut rng)).collect();
        if !pts
            .iter()
            .all(|&pt| curve.mul(&mut ops, pt, order).infinity)
        {
            return Err(format!(
                "the twist's order (p − 1)² + t² = {order} is wrong at p = {p}: [N]P ≠ O"
            ));
        }
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        // ψ: cx = u^{1−p} = u / u^p, cy = √(cx³).
        let cx = f.mul(u, f.inv(f.frob(u)).expect("u ≠ 0"));
        let Some(cy0) = f.sqrt(f.mul(f.sqr(cx), cx)) else {
            return Err("cx³ is not a square in F_{p²}; the GLS derivation is wrong".into());
        };
        let s = sqrt_mod(r - 1, r).ok_or("r ≢ 1 (mod 4)")?;
        let mut psi = None;
        'outer: for cy in [cy0, f.neg(cy0)] {
            for lambda in [s, r - s] {
                let cand = GlsEndomorphism {
                    cx,
                    cy,
                    eigenvalue: lambda,
                    p,
                };
                if cand.apply(&curve, generator) == curve.mul(&mut ops, generator, lambda) {
                    psi = Some(cand);
                    break 'outer;
                }
            }
        }
        let Some(psi) = psi else {
            return Err(format!("no sign of ψ acts as √−1 on G at p = {p}"));
        };
        verify_endomorphism(&curve, generator, r, &psi, 20, seed)?;
        return Ok(GlsInstance {
            name: format!("gls-p{p_bits}bit-p{p}"),
            p,
            base_a: a,
            base_b: b,
            trace: t,
            u,
            curve,
            group_order: order,
            r,
            cofactor: h,
            generator,
            psi,
            line: f.mul(u, f.s()),
        });
    }
    Err(format!(
        "no GLS instance found at p_bits = {p_bits} with cofactor ≤ {max_cofactor}"
    ))
}

/// **The `ψ`-stable line base** `{P : x(P) ∈ u·s·F_p^*}`, folded by
/// `⟨−1, ψ⟩` (`fold = true`, four signed points per column) or by
/// negation alone (the control on the same points).  The set is
/// checked to be invariant, not assumed.
pub fn gls_line_base(
    inst: &GlsInstance,
    fold: bool,
) -> Result<(FactorBase<Fp2Point>, FoldReport), String> {
    let start = std::time::Instant::now();
    let curve = &inst.curve;
    let f = &curve.f;
    let mut seed = Vec::with_capacity(inst.p as usize);
    let mut sqrts = 0u64;
    for t in 1..inst.p {
        let x = f.scale(inst.line, t);
        sqrts += 1;
        seed.extend(curve.lift_x(x));
    }
    let neg = Negation { r: inst.r };
    let gens: Vec<&dyn Endomorphism<Fp2Curve>> = if fold {
        vec![&neg, &inst.psi]
    } else {
        vec![&neg]
    };
    let description = format!(
        "the ψ-stable line x ∈ u·s·F_p ({} abscissae), folded by {}",
        inst.p - 1,
        if fold { "⟨−1, ψ⟩" } else { "negation" }
    );
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        seed,
        &gens,
        Closure::Strict,
        8,
        |p| curve.key(p),
        |p| f.pack(p.x),
        description,
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("abscissae_scanned", inst.p - 1);
    fb.cost.count("sqrt_solves", sqrts);
    Ok((fb, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::eigenvalue_order;
    use crate::cryptanalysis::ic_boundary::{collect_and_solve, RestartPool, TargetSource};
    use crate::cryptanalysis::ic_framework::plugins::SubtractOracle;
    use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx};

    #[test]
    fn fp2_arithmetic_is_a_field() {
        let f = Fp2::new(1009).unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..200 {
            let a = [rng.gen_range(0..1009), rng.gen_range(0..1009)];
            let b = [rng.gen_range(0..1009), rng.gen_range(0..1009)];
            assert_eq!(f.mul(a, b), f.mul(b, a));
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a).unwrap()), f.one());
            }
            // Frobenius is a field automorphism of order 2.
            assert_eq!(f.frob(f.mul(a, b)), f.mul(f.frob(a), f.frob(b)));
            assert_eq!(f.frob(f.frob(a)), a);
            assert_eq!(f.pow(a, 1009), f.frob(a), "a^p is the conjugate");
            // Every square has a root, and non-squares have none.
            let sq = f.sqr(a);
            let root = f.sqrt(sq).expect("a square has a root");
            assert!(root == a || root == f.neg(a));
            assert_eq!(f.is_square(a), f.sqrt(a).is_some());
        }
        assert_eq!(f.pow(f.s(), 1009), f.neg(f.s()), "s^p = −s");
    }

    #[test]
    fn a_gls_instance_has_a_verified_psi_of_order_four() {
        let inst = generate_gls_instance(10, 1, 16).unwrap();
        assert_eq!(
            inst.group_order as i128,
            (inst.p as i128 - 1).pow(2) + (inst.trace as i128).pow(2)
        );
        let check =
            verify_endomorphism(&inst.curve, inst.generator, inst.r, &inst.psi, 30, 1).unwrap();
        assert_eq!(check.eigenvalue_order, 4);
        assert_eq!(eigenvalue_order(inst.psi.eigenvalue, inst.r), 4);
        // ψ² = −1 on random points of the whole group, not only ⟨G⟩.
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..20 {
            let pt = random_point(&inst.curve, &mut rng);
            let psi2 = inst.psi.apply(&inst.curve, inst.psi.apply(&inst.curve, pt));
            assert_eq!(psi2, inst.curve.neg(pt), "ψ² = −1");
        }
    }

    /// The line `u·s·F_p` is `ψ`-stable and `u·F_p` carries no points:
    /// the two facts the module header derives.
    #[test]
    fn the_line_is_psi_stable_and_the_other_line_is_empty() {
        let inst = generate_gls_instance(9, 2, 16).unwrap();
        let f = &inst.curve.f;
        for t in 1..inst.p {
            let x = f.scale(inst.line, t);
            let image = f.mul(inst.psi.cx, f.frob(x));
            assert_eq!(image, f.neg(x), "ψ_x(u s t) = −u s t");
        }
        let mut on_other_line = 0;
        for t in 1..inst.p {
            on_other_line += inst.curve.lift_x(f.scale(inst.u, t)).len();
        }
        assert!(
            on_other_line <= 3,
            "u·F_p carries only 2-torsion: {on_other_line}"
        );
    }

    #[test]
    fn the_line_base_folds_four_to_one_and_recovers_the_logarithm() {
        let inst = generate_gls_instance(9, 3, 16).unwrap();
        let (folded, rep) = gls_line_base(&inst, true).unwrap();
        let (control, crep) = gls_line_base(&inst, false).unwrap();
        assert_eq!(folded.points.len(), control.points.len());
        assert!(
            (rep.points_per_orbit - 4.0).abs() < 1e-9,
            "{}",
            rep.points_per_orbit
        );
        assert!((crep.points_per_orbit - 2.0).abs() < 1e-9);
        assert_eq!(folded.columns * 2, control.columns);
        assert!(
            folded.points.len() > inst.p as usize / 2,
            "about p signed points"
        );

        let planted = 4321 % inst.r;
        for fb in [&folded, &control] {
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &inst.curve,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(2),
            };
            let mut oracle = SubtractOracle;
            let out = collect_and_solve(
                &inst.curve,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                fb,
                5,
                2_000_000,
                TargetSource::Walk,
                RestartPool::Lazy,
                |ops, ctr, pt| oracle.decompose(&ctx, fb, ops, ctr, pt),
            );
            assert_eq!(out.recovered, Some(planted), "{}", fb.description);
            assert!(out.verified);
        }
    }
}
