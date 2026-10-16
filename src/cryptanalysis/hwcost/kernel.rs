//! The Pollard rho iteration the model is calibrated on: an r-adding walk
//! of `B` independent affine walks on `y² = x³ + 7` over a prime field,
//! advanced together so that one Montgomery batch inversion serves the
//! whole batch.
//!
//! One step of one walk is `P ← P + R_j`, `j` taken from the low bits of
//! `x(P)`:
//!
//! ```text
//! dx = x_R − x_P                      (1 sub)          per walk
//! batch-invert dx                     (3 mul, 1/B inv) per walk, two of
//!                                     the three on the prefix-product
//!                                     chain (serially dependent)
//! λ  = (y_R − y_P) · dx⁻¹             (1 sub, 1 mul)
//! x₃ = λ² − x_P − x_R                 (1 sqr, 2 sub)
//! y₃ = λ · (x_P − x₃) − y_P           (2 sub, 1 mul)
//! ```
//!
//! so `5M + 1S + 6A + I/B` per iteration, of which `2M` (and the inversion)
//! form one dependent chain through the whole batch and run at latency
//! however many walks are in flight.  [`count_ops`] confirms the counts by
//! running the same code through a counting field instead of asserting
//! them; [`Arith::mul_dep`] marks the chained multiplications in the code.
//! The kernel is generic over [`Arith`]; [`Plain`] wraps any
//! [`PrimeField`] at no cost and [`Counted`] tallies every call.

use std::cell::Cell;

use num_bigint::RandBigInt;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::cryptanalysis::fp_backend::PrimeField;

/// Number of precomputed points `R_j` of the r-adding walk.
pub const TABLE: usize = 32;

/// The least significant 64 bits of an element's representation, used for
/// the partition index and the distinguished-point test.  Any
/// deterministic function of the point serves a walk; this one is free.
pub trait LowWord {
    fn low_word(&self) -> u64;
}

impl LowWord for u32 {
    #[inline(always)]
    fn low_word(&self) -> u64 {
        u64::from(*self)
    }
}

impl LowWord for u64 {
    #[inline(always)]
    fn low_word(&self) -> u64 {
        *self
    }
}

impl<const N: usize> LowWord for [u64; N] {
    #[inline(always)]
    fn low_word(&self) -> u64 {
        self[0]
    }
}

/// The field operations the walk performs.
pub trait Arith {
    type E: Copy + PartialEq;
    fn one(&self) -> Self::E;
    fn is_zero(&self, a: Self::E) -> bool;
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E;
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E;
    /// A multiplication on a serial dependency chain (the batch
    /// inversion's running product): the same operation as [`Arith::mul`],
    /// tallied apart because no amount of independent work hides its
    /// latency.
    #[inline(always)]
    fn mul_dep(&self, a: Self::E, b: Self::E) -> Self::E {
        self.mul(a, b)
    }
    fn sqr(&self, a: Self::E) -> Self::E;
    fn inv(&self, a: Self::E) -> Option<Self::E>;
    fn low(&self, a: Self::E) -> u64;
}

/// A [`PrimeField`] used directly.
pub struct Plain<'a, F: PrimeField>(pub &'a F);

impl<F: PrimeField> Arith for Plain<'_, F>
where
    F::Elem: LowWord,
{
    type E = F::Elem;
    #[inline(always)]
    fn one(&self) -> F::Elem {
        self.0.one()
    }
    #[inline(always)]
    fn is_zero(&self, a: F::Elem) -> bool {
        self.0.is_zero(a)
    }
    #[inline(always)]
    fn sub(&self, a: F::Elem, b: F::Elem) -> F::Elem {
        self.0.sub(a, b)
    }
    #[inline(always)]
    fn mul(&self, a: F::Elem, b: F::Elem) -> F::Elem {
        self.0.mul(a, b)
    }
    #[inline(always)]
    fn sqr(&self, a: F::Elem) -> F::Elem {
        self.0.sqr(a)
    }
    #[inline(always)]
    fn inv(&self, a: F::Elem) -> Option<F::Elem> {
        self.0.inv(a)
    }
    #[inline(always)]
    fn low(&self, a: F::Elem) -> u64 {
        a.low_word()
    }
}

/// Exact tally of the field operations of a run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpTally {
    /// Multiplications with independent work around them.
    pub mul: u64,
    /// Multiplications on the batch's serial chain ([`Arith::mul_dep`]).
    pub mul_dep: u64,
    pub sqr: u64,
    /// Additions and subtractions.
    pub add: u64,
    pub inv: u64,
}

/// A [`PrimeField`] whose every [`Arith`] call is counted.
pub struct Counted<'a, F: PrimeField> {
    pub f: &'a F,
    pub tally: Cell<OpTally>,
}

impl<'a, F: PrimeField> Counted<'a, F> {
    pub fn new(f: &'a F) -> Self {
        Counted {
            f,
            tally: Cell::new(OpTally::default()),
        }
    }
    fn bump(&self, k: impl FnOnce(&mut OpTally)) {
        let mut t = self.tally.get();
        k(&mut t);
        self.tally.set(t);
    }
}

impl<F: PrimeField> Arith for Counted<'_, F>
where
    F::Elem: LowWord,
{
    type E = F::Elem;
    fn one(&self) -> F::Elem {
        self.f.one()
    }
    fn is_zero(&self, a: F::Elem) -> bool {
        self.f.is_zero(a)
    }
    fn sub(&self, a: F::Elem, b: F::Elem) -> F::Elem {
        self.bump(|t| t.add += 1);
        self.f.sub(a, b)
    }
    fn mul(&self, a: F::Elem, b: F::Elem) -> F::Elem {
        self.bump(|t| t.mul += 1);
        self.f.mul(a, b)
    }
    fn mul_dep(&self, a: F::Elem, b: F::Elem) -> F::Elem {
        self.bump(|t| t.mul_dep += 1);
        self.f.mul(a, b)
    }
    fn sqr(&self, a: F::Elem) -> F::Elem {
        self.bump(|t| t.sqr += 1);
        self.f.sqr(a)
    }
    fn inv(&self, a: F::Elem) -> Option<F::Elem> {
        self.bump(|t| t.inv += 1);
        self.f.inv(a)
    }
    fn low(&self, a: F::Elem) -> u64 {
        a.low_word()
    }
}

/// `B` walks and their scratch, laid out as separate arrays.
pub struct Walks<E> {
    pub xs: Vec<E>,
    pub ys: Vec<E>,
    dx: Vec<E>,
    scratch: Vec<E>,
    /// Steps whose `dx` was zero (`x_P = x_R`); the step is then wrong, so
    /// a run that reports any is rejected.  Probability about `B/p` per
    /// step.
    pub degenerate: u64,
    /// Distinguished points seen (low `dp_bits` bits of `x` zero).
    pub distinguished: u64,
}

impl<E: Copy> Walks<E> {
    /// Walks starting at the points `(xs[i], ys[i])`.
    pub fn new<F: PrimeField<Elem = E>>(f: &F, xs: Vec<E>, ys: Vec<E>) -> Self {
        assert_eq!(xs.len(), ys.len());
        let b = xs.len();
        Walks {
            xs,
            ys,
            dx: vec![f.zero(); b],
            scratch: vec![f.zero(); b],
            degenerate: 0,
            distinguished: 0,
        }
    }
}

/// The walk's precomputed points `R_0 … R_{TABLE−1}`.
pub struct Table<E> {
    pub xs: [E; TABLE],
    pub ys: [E; TABLE],
}

/// Element bytes in memory for a field element type.
pub fn elem_bytes<E>() -> usize {
    std::mem::size_of::<E>()
}

/// Bytes the kernel reads plus bytes it writes per iteration, counted from
/// the loops of [`step`]: reads `x` (partition), `dx`, `dx`+`scratch`
/// (backward pass), `x, y, dx` (update): 7 elements; writes `dx`,
/// `scratch`, `dx⁻¹`, `x₃, y₃`: 5 elements.  Table reads stay in L1 and
/// are not counted.
pub const ELEMS_MOVED_PER_ITERATION: usize = 12;

/// Elements resident per walk: `x`, `y`, `dx`, `scratch`.
pub const ELEMS_RESIDENT_PER_WALK: usize = 4;

#[inline(always)]
fn part<A: Arith>(a: &A, x: A::E) -> usize {
    ((a.low(x) >> 1) as usize) & (TABLE - 1)
}

/// One iteration of every walk.
#[inline(always)]
pub fn step<A: Arith>(a: &A, w: &mut Walks<A::E>, t: &Table<A::E>, dp_mask: u64) {
    let b = w.xs.len();
    for i in 0..b {
        let j = part(a, w.xs[i]);
        w.dx[i] = a.sub(t.xs[j], w.xs[i]);
    }
    // Montgomery's trick: 3 multiplications per element, one inversion.
    let mut acc = a.one();
    for i in 0..b {
        w.scratch[i] = acc;
        if a.is_zero(w.dx[i]) {
            w.degenerate += 1;
        } else {
            acc = a.mul_dep(acc, w.dx[i]);
        }
    }
    let mut inv = match a.inv(acc) {
        Some(v) => v,
        None => {
            w.degenerate += 1;
            return;
        }
    };
    for i in (0..b).rev() {
        let d = w.dx[i];
        if a.is_zero(d) {
            continue;
        }
        w.dx[i] = a.mul(inv, w.scratch[i]);
        inv = a.mul_dep(inv, d);
    }
    for i in 0..b {
        let x = w.xs[i];
        let y = w.ys[i];
        let j = part(a, x);
        let (xr, yr) = (t.xs[j], t.ys[j]);
        let l = a.mul(a.sub(yr, y), w.dx[i]);
        let x3 = a.sub(a.sub(a.sqr(l), x), xr);
        let y3 = a.sub(a.mul(l, a.sub(x, x3)), y);
        w.xs[i] = x3;
        w.ys[i] = y3;
        if a.low(x3) & dp_mask == 0 {
            w.distinguished += 1;
        }
    }
}

/// The same step with one inversion per walk and no batching: the
/// reference [`step`] is checked against.
pub fn step_unbatched<F: PrimeField>(
    f: &F,
    xs: &mut [F::Elem],
    ys: &mut [F::Elem],
    t: &Table<F::Elem>,
) where
    F::Elem: LowWord,
{
    for i in 0..xs.len() {
        let (x, y) = (xs[i], ys[i]);
        let j = ((x.low_word() >> 1) as usize) & (TABLE - 1);
        let (xr, yr) = (t.xs[j], t.ys[j]);
        let l = f.mul(f.sub(yr, y), f.inv(f.sub(xr, x)).expect("x_P ≠ x_R"));
        let x3 = f.sub(f.sub(f.sqr(l), x), xr);
        ys[i] = f.sub(f.mul(l, f.sub(x, x3)), y);
        xs[i] = x3;
    }
}

/// `y² = x³ + 7`.
pub fn on_curve<F: PrimeField>(f: &F, x: F::Elem, y: F::Elem) -> bool {
    let rhs = f.add(f.mul(f.sqr(x), x), f.from_u64(7));
    f.sqr(y) == rhs
}

/// A uniformly random point of `y² = x³ + 7`.
pub fn random_point<F: PrimeField>(f: &F, rng: &mut StdRng) -> (F::Elem, F::Elem) {
    loop {
        let x = f.from_biguint(&rng.gen_biguint_below(f.modulus()));
        let rhs = f.add(f.mul(f.sqr(x), x), f.from_u64(7));
        if let Some(y) = f.sqrt(rhs) {
            return (x, y);
        }
    }
}

/// `b` random walks and a random table, from `seed`.
pub fn setup<F: PrimeField>(f: &F, b: usize, seed: u64) -> (Walks<F::Elem>, Table<F::Elem>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut txs = [f.zero(); TABLE];
    let mut tys = [f.zero(); TABLE];
    for j in 0..TABLE {
        let (x, y) = random_point(f, &mut rng);
        txs[j] = x;
        tys[j] = y;
    }
    let mut xs = Vec::with_capacity(b);
    let mut ys = Vec::with_capacity(b);
    for _ in 0..b {
        let (x, y) = random_point(f, &mut rng);
        xs.push(x);
        ys.push(y);
    }
    (Walks::new(f, xs, ys), Table { xs: txs, ys: tys })
}

/// The field operations per iteration of [`step`] at batch `b`, counted
/// over `steps` steps of a real walk.
pub fn count_ops<F: PrimeField>(f: &F, b: usize, steps: usize, seed: u64) -> (OpTally, f64)
where
    F::Elem: LowWord,
{
    let (mut w, t) = setup(f, b, seed);
    let c = Counted::new(f);
    for _ in 0..steps {
        step(&c, &mut w, &t, 0);
    }
    assert_eq!(w.degenerate, 0, "degenerate step while counting");
    (c.tally.get(), (b * steps) as f64)
}

/// Run `steps` batched steps and the unbatched reference side by side and
/// confirm they agree and every point stays on the curve.
pub fn check<F: PrimeField>(f: &F, b: usize, steps: usize, seed: u64) -> Result<(), String>
where
    F::Elem: LowWord,
{
    let (mut w, t) = setup(f, b, seed);
    let mut rx = w.xs.clone();
    let mut ry = w.ys.clone();
    let p = Plain(f);
    for s in 0..steps {
        step(&p, &mut w, &t, 0);
        step_unbatched(f, &mut rx, &mut ry, &t);
        if w.degenerate != 0 {
            return Err(format!("{}: degenerate step at step {s}", f.name()));
        }
        if w.xs != rx || w.ys != ry {
            return Err(format!(
                "{}: batched and unbatched walks differ at step {s}",
                f.name()
            ));
        }
        for i in 0..b {
            if !on_curve(f, w.xs[i], w.ys[i]) {
                return Err(format!("{}: walk {i} left the curve at step {s}", f.name()));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::fp_backend::{Fp64, FpWide};
    use num_bigint::BigUint;

    #[test]
    fn batched_walk_matches_reference_and_stays_on_curve() {
        let f64_ = Fp64::new((1u64 << 63) - 25);
        check(&f64_, 17, 5, 1).unwrap();
        let p127 = (BigUint::from(1u32) << 127u32) - 1u32;
        check(&FpWide::<2>::new(&p127), 9, 4, 2).unwrap();
    }

    #[test]
    fn counted_ops_are_5m_1s_6a_one_inversion_per_batch_two_m_chained() {
        let f = Fp64::new((1u64 << 63) - 25);
        let (t, iters) = count_ops(&f, 64, 3, 7);
        assert_eq!(iters, 192.0);
        assert_eq!(t.mul, 3 * 192);
        assert_eq!(t.mul_dep, 2 * 192);
        assert_eq!(t.sqr, 192);
        assert_eq!(t.add, 6 * 192);
        assert_eq!(t.inv, 3);
    }
}
