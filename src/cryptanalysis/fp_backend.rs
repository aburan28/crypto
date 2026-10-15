//! # Variable-time prime-field backends for research engines
//!
//! One trait, [`PrimeField`], and four implementations, for code that has
//! to run the same algorithm (F4 row reduction, curve arithmetic, Vélu
//! chains) at a 31-bit prime, a 63-bit prime and a NIST prime without
//! rewriting it:
//!
//! | backend | primes | element | representation |
//! |:--|:--|:--|:--|
//! | [`Fp32`] | odd `p < 2³²` | `u32` | canonical, Barrett `⌊2⁶⁴/p⌋` |
//! | [`Fp64`] | odd `p < 2⁶³` | `u64` | Montgomery, `R = 2⁶⁴` |
//! | [`FpWide<N>`] | odd `p < 2^(64N)`, `N ≤ 8` | `[u64; N]` | Montgomery (CIOS), `R = 2^(64N)` |
//! | [`FpWide<N, true>`](FpWide) | odd `p`, top limb `< 2⁶³ − 1` | `[u64; N]` | Montgomery, carry-free CIOS |
//! | [`FpP192`] | NIST P-192 only | `[u64; 3]` | canonical, Solinas reduction |
//!
//! Aliases: [`Fp192`] = `FpWide<3>` (P-192), [`Fp256`] = `FpWide<4>`
//! (P-256, any prime below `2²⁵⁶`), [`Fp224`] = `FpWide<4, true>` (P-224,
//! the 200-bit test prime). [`FpWide`] is the default and the correctness
//! reference for wide primes; [`FpP192`] and the `SPARE` variant exist
//! because they measured faster (`research/fp_backend_20261015/README.md`).
//! NIST P-224 has no Solinas specialisation (its `2¹²⁸`-free form does not
//! align with 64-bit words); use [`Fp224`].
//!
//! **Nothing here is constant-time.** Inversion, square roots, the
//! Legendre symbol, the conditional subtractions and the lazy
//! accumulators branch on data. Use [`crate::ct_bignum`] for secrets.
//!
//! ## Elements
//!
//! A field `F: PrimeField` owns its constants; elements are plain `Copy`
//! values (`F::Elem`) passed by value to `&self` methods. Every element
//! is stored fully reduced (`< p`) in a unique representation, so `==`
//! and `Hash` on `F::Elem` are field equality and hashing. Montgomery
//! backends store `a·R mod p`; never read the raw words as the integer —
//! use [`PrimeField::to_biguint`] / [`PrimeField::to_u64`] and
//! [`PrimeField::from_u64`] / [`PrimeField::from_biguint`].
//!
//! ## Arithmetic
//!
//! `add`, `sub`, `neg`, `double`, `mul`, `sqr`, `inv` (returns `None` on
//! zero), `batch_inv` (Montgomery's trick, zeros are left as zero),
//! `pow` / `pow_u64` / `pow_biguint` (4-bit fixed window), `legendre`
//! (`1`, `−1`, or `0`), `sqrt` (`None` on a non-residue), `is_square`.
//!
//! - **Inverse:** Pornin's optimised binary GCD ([`inv_bingcd`]: 31
//!   divsteps at a time on 64-bit approximations, then one update of the
//!   full values) on every backend, the `u64` ones included.
//! - **Legendre:** a binary Jacobi symbol (shifts and subtractions, no
//!   exponentiation) on every backend.
//! - **Square root:** `p ≡ 3 mod 4`: `a^((p+1)/4)`. Otherwise
//!   Tonelli–Shanks from a non-residue found at construction; when
//!   `s = v₂(p − 1) ≥ 16` (P-224 has `s = 96`) the 2-Sylow logarithm is
//!   read 8 or 12 bits at a time from a table of `2^w` entries
//!   ([`SqrtChunks`]).
//!
//! ## The lazy accumulator (`F::Wide`)
//!
//! Row operations and dot products should not reduce after every
//! product. `F::Wide` is an unreduced accumulator that represents one
//! field element:
//!
//! ```text
//! let mut w = f.wide_zero();             // or f.wide_from(x)   (a load)
//! f.wide_mul_add(&mut w, a, b);           // w += a·b, unreduced
//! ... repeated ...
//! let x = f.wide_reduce(&w);              // one reduction at the end
//! ```
//!
//! and the row kernels apply it to whole slices:
//!
//! ```text
//! f.row_load(&row, &mut acc);             // acc[i] = row[i]
//! f.row_sub(&mut acc, c, &pivot);         // acc[i] -= c·pivot[i]  (lazy)
//! f.row_axpy(&mut acc, c, &pivot);        // acc[i] += c·pivot[i]  (lazy)
//! f.row_reduce(&acc, &mut out);           // out[i] = acc[i] mod p
//! f.row_reduce_scaled(&acc, s, &mut out); // out[i] = s·acc[i]  (F4 makes rows monic)
//! f.row_scatter_axpy(&mut acc, c, &cols, &vals); // acc[cols[k]] += c·vals[k]  (sparse)
//! ```
//!
//! **Contract.** After a load (or `wide_zero`) at most
//! [`PrimeField::lazy_capacity`] products may be added before
//! `wide_reduce`; one `row_axpy`/`row_sub` call adds one product to each
//! entry. Debug builds check the bound inside every backend (an `assert`
//! on the invariant that the bound protects). [`LazyRow`] counts the
//! products for you and renormalises when the bound is reached, so code
//! built on it never needs to know the bound. The capacities, derived:
//!
//! - **`Fp32`, `p < 2³¹`.** `Wide = u64`. A product is `< p² < 2⁶²`;
//!   the entry is kept `< 2⁶³` by subtracting `M`, the largest multiple
//!   of `p` not above `2⁶³`, whenever the sum reaches `2⁶³`
//!   (`s < 2⁶³ + 2⁶²`, so `s − M < 2⁶² + p < 2⁶³`). The invariant holds
//!   after every product: **unbounded** (`u64::MAX`). This is the
//!   `lazy_axpy` kernel of `f4_fp_tower`, which vectorises (AVX2/AVX-512
//!   dispatch at construction).
//! - **`Fp32`, `2³¹ ≤ p < 2³²`.** A product no longer leaves room in a
//!   `u64`. Rows use Shoup multiplication by the fixed `c`
//!   (`c' = ⌊c·2³²/p⌋`, `c·y − ⌊c'·y/2³²⌋·p ∈ [0, 2p)`, three 32×32→64
//!   multiplies, result `< 2p`), single products unreduced; either way the
//!   product is added to the `u64` entry and a carry out of the word
//!   (worth `2⁶⁴`) is folded back in as `2⁶⁴ mod p`. Every `u64` is a
//!   valid entry, read through Barrett reduction: **unbounded**.
//! - **`Fp64`.** `Wide = u128` holding `W` with `W·R⁻¹` the value. A
//!   product of two elements is `< p²`; the invariant `W < p·2⁶⁴` is
//!   kept by subtracting `p·2⁶⁴` (a multiple of `p`) when it is reached
//!   (`W + p² < 2p·2⁶⁴ < 2¹²⁸` because `p < 2⁶³`), and a single
//!   Montgomery reduction of `W < p·2⁶⁴` is fully reduced after one
//!   conditional subtraction: **unbounded**.
//! - **`FpWide<N>`.** `Wide` is `2N + 1` words with no renormalisation.
//!   A product is `< p² < p·R`, a load is `< p·R`; the reduction is a
//!   Montgomery reduction over `N + 1` words, valid for `W < p·R·2⁶⁴`.
//!   So a load plus `k` products needs `(k + 1)·p·R ≤ p·R·2⁶⁴`:
//!   **`k ≤ 2⁶⁴ − 2`** (`lazy_capacity() = u64::MAX − 1`). Debug builds
//!   assert `W < p·R·2⁶⁴` in `wide_reduce` and that the top word never
//!   carries out. The reduction then multiplies by the constant
//!   `2⁶⁴·R mod p` to undo the extra word (one Montgomery product per
//!   reduced entry; `row_reduce_scaled` folds it into the scale).
//! - **`FpP192`.** `Wide` is 7 words (a 384-bit product sum plus a
//!   carry word), folded by `2¹⁹² ≡ 2⁶⁴ + 1` at the end: **`k ≤ 2⁶⁴ − 2`**,
//!   asserted the same way.
//!
//! ## How the engines should use it
//!
//! **F4 (lane F).** Make the engine generic over `F: PrimeField`, store
//! coefficients as `F::Elem`, and replace the eliminations:
//!
//! - dense block `acc −= v·bred[j]`: `Vec<F::Wide>` per row,
//!   `row_load` / `row_sub` / `row_reduce`; one `row_sub` per pivot row.
//! - sparse scatter (`Acc::sub_tail`): keep `vals: Vec<F::Wide>`, call
//!   `row_scatter_axpy(&mut vals, nf, &row.cols[1..], &row.vals[1..])`
//!   and set the bitset bits in a separate pass (per-entry `wide_mul_add`
//!   is correct too, but measured 1.6× slower for `Fp32` above `2³¹`
//!   because the regime test stays inside the loop); read an entry with `wide_reduce` before testing it for zero;
//!   `drain_from` becomes `wide_reduce` + `mul` by the new lead's inverse.
//! - for `Fp32` with `p < 2³¹` this is exactly today's arithmetic
//!   (measured in the report); for `p ≥ 2³¹`, `Fp64` and wider fields the
//!   capacities above are never reached by an F4 row, so no counting is
//!   needed — [`LazyRow`] is there if in doubt.
//!
//! **Curves, square roots and Vélu (lane G).** Use [`FpP192`] for NIST
//! P-192 (`p ≡ 3 mod 4`, so `sqrt` is one exponentiation; `Fp192` gives
//! the same answers more slowly), [`Fp224`] for P-224 and [`Fp256`] for
//! P-256. Batch inversions with `batch_inv`; one `inv` costs on the order
//! of 100–250 multiplications at 192–256 bits (measured in the report), so
//! keep points projective where an inversion per step would dominate.
//! Elements implement `Hash`, so j-invariants can key a `HashMap`.
//!
//! The benchmark driver is `examples/fp_backend_bench.rs`.

use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive};

// ── The trait ─────────────────────────────────────────────────────────

/// A prime field `F_p` with its constants; see the module documentation.
///
/// Conversions take `&self` because the field owns the constants they need.
#[allow(clippy::wrong_self_convention)]
pub trait PrimeField: Clone + Debug + Send + Sync {
    /// A field element in this backend's representation (always `< p`).
    type Elem: Copy + PartialEq + Eq + Hash + Debug + Send + Sync;
    /// An unreduced accumulator holding one field element.
    type Wide: Copy + Debug + Send + Sync;

    /// Short backend name, for reports.
    fn name(&self) -> &'static str;
    /// The prime `p`.
    fn modulus(&self) -> &BigUint;
    /// `⌈log₂ p⌉`.
    fn bits(&self) -> u32 {
        self.modulus().bits() as u32
    }

    fn zero(&self) -> Self::Elem;
    fn one(&self) -> Self::Elem;
    /// `x mod p`.
    fn from_u64(&self, x: u64) -> Self::Elem;
    /// `x mod p`.
    fn from_biguint(&self, x: &BigUint) -> Self::Elem;
    /// The canonical integer in `[0, p)`.
    fn to_biguint(&self, a: Self::Elem) -> BigUint;
    /// The canonical integer, when it fits a `u64`.
    fn to_u64(&self, a: Self::Elem) -> Option<u64> {
        self.to_biguint(a).to_u64()
    }
    /// `x mod p` for a signed `x`.
    fn from_i64(&self, x: i64) -> Self::Elem {
        let a = self.from_u64(x.unsigned_abs());
        if x < 0 {
            self.neg(a)
        } else {
            a
        }
    }
    fn is_zero(&self, a: Self::Elem) -> bool {
        a == self.zero()
    }
    fn is_one(&self, a: Self::Elem) -> bool {
        a == self.one()
    }

    fn add(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;
    fn sub(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;
    fn neg(&self, a: Self::Elem) -> Self::Elem;
    fn double(&self, a: Self::Elem) -> Self::Elem {
        self.add(a, a)
    }
    fn mul(&self, a: Self::Elem, b: Self::Elem) -> Self::Elem;
    #[inline(always)]
    fn sqr(&self, a: Self::Elem) -> Self::Elem {
        self.mul(a, a)
    }
    /// `a⁻¹`, or `None` for zero.
    fn inv(&self, a: Self::Elem) -> Option<Self::Elem>;

    /// Invert every non-zero entry of `xs` in place (zeros stay zero) with
    /// one inversion and `3(k − 1)` multiplications.
    fn batch_inv(&self, xs: &mut [Self::Elem], scratch: &mut Vec<Self::Elem>) {
        scratch.clear();
        scratch.reserve(xs.len());
        let mut acc = self.one();
        for &x in xs.iter() {
            scratch.push(acc);
            if !self.is_zero(x) {
                acc = self.mul(acc, x);
            }
        }
        let mut inv = match self.inv(acc) {
            Some(v) => v,
            None => return,
        };
        for i in (0..xs.len()).rev() {
            let x = xs[i];
            if self.is_zero(x) {
                continue;
            }
            xs[i] = self.mul(inv, scratch[i]);
            inv = self.mul(inv, x);
        }
    }

    /// `a^e` for an exponent given as little-endian 64-bit limbs
    /// (fixed 4-bit window).
    fn pow(&self, a: Self::Elem, e: &[u64]) -> Self::Elem {
        pow_window(self, a, e)
    }
    fn pow_u64(&self, a: Self::Elem, e: u64) -> Self::Elem {
        let mut acc = self.one();
        let mut base = a;
        let mut e = e;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(acc, base);
            }
            e >>= 1;
            if e > 0 {
                base = self.sqr(base);
            }
        }
        acc
    }
    fn pow_biguint(&self, a: Self::Elem, e: &BigUint) -> Self::Elem {
        self.pow(a, &e.to_u64_digits())
    }

    /// Precomputed square-root / Legendre constants.
    fn sqrt_pre(&self) -> &SqrtPre<Self::Elem>;
    /// The Legendre symbol `(a / p)`: `1`, `−1`, or `0` for zero.
    fn legendre(&self, a: Self::Elem) -> i32 {
        if self.is_zero(a) {
            return 0;
        }
        let r = self.pow(a, &self.sqrt_pre().p_minus1_half);
        if self.is_one(r) {
            1
        } else {
            -1
        }
    }
    fn is_square(&self, a: Self::Elem) -> bool {
        self.legendre(a) >= 0
    }
    /// A square root of `a` (either one), or `None` for a non-residue.
    fn sqrt(&self, a: Self::Elem) -> Option<Self::Elem> {
        sqrt_generic(self, a)
    }

    // ── lazy accumulation ────────────────────────────────────────────

    /// How many products may be added to a loaded accumulator before
    /// `wide_reduce` (module documentation, "The lazy accumulator").
    fn lazy_capacity(&self) -> u64;
    fn wide_zero(&self) -> Self::Wide;
    /// An accumulator holding `a` (counts as the load, not a product).
    fn wide_from(&self, a: Self::Elem) -> Self::Wide;
    /// `w += a·b` without reduction.
    fn wide_mul_add(&self, w: &mut Self::Wide, a: Self::Elem, b: Self::Elem);
    /// The element `w` holds.
    fn wide_reduce(&self, w: &Self::Wide) -> Self::Elem;

    /// `acc[i] = row[i]` (a load of every entry).
    fn row_load(&self, row: &[Self::Elem], acc: &mut [Self::Wide]) {
        for (w, &x) in acc.iter_mut().zip(row) {
            *w = self.wide_from(x);
        }
    }
    /// `acc[i] += c·row[i]` for `i < min(len)`: one product per entry.
    fn row_axpy(&self, acc: &mut [Self::Wide], c: Self::Elem, row: &[Self::Elem]) {
        for (w, &y) in acc.iter_mut().zip(row) {
            self.wide_mul_add(w, c, y);
        }
    }
    /// Sparse `acc[cols[k]] += c·vals[k]` (one product per touched entry;
    /// `cols` need not be sorted but must index `acc`).
    fn row_scatter_axpy(
        &self,
        acc: &mut [Self::Wide],
        c: Self::Elem,
        cols: &[u32],
        vals: &[Self::Elem],
    ) {
        for (&col, &y) in cols.iter().zip(vals) {
            self.wide_mul_add(&mut acc[col as usize], c, y);
        }
    }
    /// `acc[i] −= c·row[i]`: one product per entry.
    fn row_sub(&self, acc: &mut [Self::Wide], c: Self::Elem, row: &[Self::Elem]) {
        self.row_axpy(acc, self.neg(c), row)
    }
    /// `out[i] = acc[i]` reduced.
    fn row_reduce(&self, acc: &[Self::Wide], out: &mut [Self::Elem]) {
        for (o, w) in out.iter_mut().zip(acc) {
            *o = self.wide_reduce(w);
        }
    }
    /// `out[i] = s · acc[i]` reduced.
    fn row_reduce_scaled(&self, acc: &[Self::Wide], s: Self::Elem, out: &mut [Self::Elem]) {
        for (o, w) in out.iter_mut().zip(acc) {
            *o = self.mul(self.wide_reduce(w), s);
        }
    }
    /// Eager in-place `dst[i] += c·row[i]`, reduced every step (the
    /// comparison arm for the lazy kernels; also the right call when a row
    /// receives a single update).
    fn row_axpy_eager(&self, dst: &mut [Self::Elem], c: Self::Elem, row: &[Self::Elem]) {
        for (x, &y) in dst.iter_mut().zip(row) {
            *x = self.add(*x, self.mul(c, y));
        }
    }
    /// `Σ a[i]·b[i]` with one reduction per `lazy_capacity` products.
    fn dot(&self, a: &[Self::Elem], b: &[Self::Elem]) -> Self::Elem {
        let cap = self.lazy_capacity();
        let mut res = self.zero();
        let mut w = self.wide_zero();
        let mut used = 0u64;
        for (&x, &y) in a.iter().zip(b) {
            if used == cap {
                res = self.add(res, self.wide_reduce(&w));
                w = self.wide_zero();
                used = 0;
            }
            self.wide_mul_add(&mut w, x, y);
            used += 1;
        }
        self.add(res, self.wide_reduce(&w))
    }
}

/// A dense lazy row that counts the products added since its last
/// reduction and renormalises before the backend's capacity is reached.
#[derive(Clone, Debug)]
pub struct LazyRow<F: PrimeField> {
    pub acc: Vec<F::Wide>,
    used: u64,
}

impl<F: PrimeField> LazyRow<F> {
    pub fn zeros(f: &F, len: usize) -> Self {
        LazyRow {
            acc: vec![f.wide_zero(); len],
            used: 0,
        }
    }
    pub fn from_row(f: &F, row: &[F::Elem]) -> Self {
        let mut acc = vec![f.wide_zero(); row.len()];
        f.row_load(row, &mut acc);
        LazyRow { acc, used: 0 }
    }
    fn make_room(&mut self, f: &F) {
        if self.used >= f.lazy_capacity() {
            for w in self.acc.iter_mut() {
                *w = f.wide_from(f.wide_reduce(w));
            }
            self.used = 0;
        }
        self.used += 1;
        debug_assert!(self.used <= f.lazy_capacity());
    }
    /// `acc[off + i] += c·row[i]`.
    pub fn axpy(&mut self, f: &F, off: usize, c: F::Elem, row: &[F::Elem]) {
        self.make_room(f);
        f.row_axpy(&mut self.acc[off..], c, row);
    }
    /// `acc[off + i] −= c·row[i]`.
    pub fn sub(&mut self, f: &F, off: usize, c: F::Elem, row: &[F::Elem]) {
        self.make_room(f);
        f.row_sub(&mut self.acc[off..], c, row);
    }
    pub fn finish(&self, f: &F, out: &mut Vec<F::Elem>) {
        out.clear();
        out.resize(self.acc.len(), f.zero());
        f.row_reduce(&self.acc, out);
    }
}

// ── Shared generic algorithms ─────────────────────────────────────────

/// Square-root and Legendre constants: `p − 1 = 2^s·q` with `q` odd.
#[derive(Clone, Debug)]
pub struct SqrtPre<E> {
    pub s: u32,
    /// `(p − 1)/2`, little-endian limbs (the Euler exponent).
    pub p_minus1_half: Vec<u64>,
    /// `(p + 1)/4` when `p ≡ 3 mod 4`.
    pub p_plus1_quarter: Option<Vec<u64>>,
    /// `(q − 1)/2`.
    pub q_minus1_half: Vec<u64>,
    /// The smallest non-residue `z ≥ 2`.
    pub non_residue: u64,
    /// `z^q`, a generator of the 2-Sylow subgroup.
    pub z_q: E,
    /// Windowed discrete-log tables for Tonelli–Shanks when `s ≥ 16`.
    pub chunks: Option<SqrtChunks<E>>,
}

/// Tables for the windowed Tonelli–Shanks: the logarithm of `t = a^q` in
/// the 2-Sylow subgroup `⟨g⟩`, `g = z^q`, is read `w` bits at a time from
/// a table of `h = g^(2^(s − w))`, about `s²/(2w)` squarings instead of
/// about `s²/4`.
#[derive(Clone, Debug)]
pub struct SqrtChunks<E> {
    pub w: u32,
    /// `h^d ↦ d` for `0 ≤ d < 2^w`.
    pub table: HashMap<E, u32>,
    /// `half[0] = g⁻¹`, `half[i] = g^(−2^(w·i − 1))` for `i ≥ 1`: the
    /// square root of window `i`'s correction.
    pub half: Vec<E>,
}

impl<E: Copy> SqrtPre<E> {
    fn placeholder(zero: E) -> Self {
        SqrtPre {
            s: 0,
            p_minus1_half: Vec::new(),
            p_plus1_quarter: None,
            q_minus1_half: Vec::new(),
            non_residue: 0,
            z_q: zero,
            chunks: None,
        }
    }
}

/// Fill a backend's [`SqrtPre`] (uses only `pow`, `from_u64`, `is_one`).
fn compute_sqrt_pre<F: PrimeField>(f: &F) -> SqrtPre<F::Elem> {
    let p = f.modulus().clone();
    let one = BigUint::one();
    let pm1 = &p - &one;
    let s = pm1.trailing_zeros().unwrap_or(0) as u32;
    let q = &pm1 >> s;
    let p_minus1_half: Vec<u64> = (&pm1 >> 1u32).to_u64_digits();
    let p_plus1_quarter = if (&p % 4u32) == BigUint::from(3u32) {
        Some(((&p + &one) >> 2u32).to_u64_digits())
    } else {
        None
    };
    let q_minus1_half = ((&q - &one) >> 1u32).to_u64_digits();
    // Smallest non-residue by Euler's criterion.
    let mut z = 2u64;
    loop {
        let r = f.pow(f.from_u64(z), &p_minus1_half);
        if !f.is_one(r) && !f.is_zero(r) {
            break;
        }
        z += 1;
    }
    let z_q = f.pow(f.from_u64(z), &q.to_u64_digits());
    let chunks = (s >= 16).then(|| {
        let w = if s >= 48 { 12 } else { 8 };
        let g_inv = f.inv(z_q).expect("z^q is a unit");
        let mut h = z_q;
        for _ in 0..s - w {
            h = f.sqr(h);
        }
        let mut table = HashMap::with_capacity(1 << w);
        let mut cur = f.one();
        for d in 0..1u32 << w {
            table.insert(cur, d);
            cur = f.mul(cur, h);
        }
        let mut half = vec![g_inv];
        let mut gi = g_inv;
        for i in 1..s.div_ceil(w) {
            // g^(−2^(w·i − 1)) from g^(−2^(w·(i−1))).
            let shift = if i == 1 { w - 1 } else { w };
            for _ in 0..shift {
                gi = f.sqr(gi);
            }
            half.push(gi);
        }
        SqrtChunks { w, table, half }
    });
    SqrtPre {
        s,
        p_minus1_half,
        p_plus1_quarter,
        q_minus1_half,
        non_residue: z,
        z_q,
        chunks,
    }
}

fn pow_window<F: PrimeField>(f: &F, a: F::Elem, e: &[u64]) -> F::Elem {
    let mut top = e.len();
    while top > 0 && e[top - 1] == 0 {
        top -= 1;
    }
    if top == 0 {
        return f.one();
    }
    let nbits = 64 * (top - 1) + (64 - e[top - 1].leading_zeros() as usize);
    if nbits <= 16 {
        return f.pow_u64(a, e[0]);
    }
    let mut table = [a; 16];
    table[0] = f.one();
    for i in 2..16 {
        table[i] = f.mul(table[i - 1], a);
    }
    let window = |pos: usize| -> usize {
        // 4-bit window whose lowest bit is `pos` (pos multiple of 4).
        ((e[pos / 64] >> (pos % 64)) & 0xf) as usize
    };
    let nwin = nbits.div_ceil(4);
    let mut acc = table[window(4 * (nwin - 1))];
    for w in (0..nwin - 1).rev() {
        acc = f.sqr(acc);
        acc = f.sqr(acc);
        acc = f.sqr(acc);
        acc = f.sqr(acc);
        let d = window(4 * w);
        if d != 0 {
            acc = f.mul(acc, table[d]);
        }
    }
    acc
}

/// Tonelli–Shanks, with the `p ≡ 3 mod 4` shortcut.
fn sqrt_generic<F: PrimeField>(f: &F, a: F::Elem) -> Option<F::Elem> {
    if f.is_zero(a) {
        return Some(a);
    }
    let pre = f.sqrt_pre();
    if let Some(e) = &pre.p_plus1_quarter {
        let r = f.pow(a, e);
        return if f.sqr(r) == a { Some(r) } else { None };
    }
    // x = a^((q−1)/2); t = a^q = a·x²; r = a^((q+1)/2) = a·x.
    let x = f.pow(a, &pre.q_minus1_half);
    let r = f.mul(a, x);
    let t = f.mul(r, x);
    if let Some(ch) = &pre.chunks {
        // r² = a·t with t = g^e (e even iff a is a square), so
        // √a = r·g^(−e/2). Read e w bits at a time, lowest first; the
        // window's correction is hp² with hp its share of g^(−e/2).
        let s = pre.s;
        let mut cur = t;
        let mut acc = r;
        for (i, &h) in ch.half.iter().enumerate() {
            let lo = ch.w * i as u32;
            let wi = ch.w.min(s - lo);
            let mut y = cur;
            for _ in 0..s - lo - wi {
                y = f.sqr(y);
            }
            let d = ch.table[&y] >> (ch.w - wi);
            let exp = if i == 0 {
                if d & 1 == 1 {
                    return None;
                }
                d >> 1
            } else {
                d
            };
            if exp != 0 {
                let hp = f.pow_u64(h, u64::from(exp));
                cur = f.mul(cur, f.sqr(hp));
                acc = f.mul(acc, hp);
            }
        }
        debug_assert!(f.is_one(cur));
        return Some(acc);
    }
    ts_classic(f, r, t)
}

/// Tonelli–Shanks without the window tables (comparison arm for the
/// benchmark; [`PrimeField::sqrt`] is the one to use).
pub fn sqrt_tonelli_shanks_plain<F: PrimeField>(f: &F, a: F::Elem) -> Option<F::Elem> {
    if f.is_zero(a) {
        return Some(a);
    }
    let pre = f.sqrt_pre();
    let x = f.pow(a, &pre.q_minus1_half);
    let r = f.mul(a, x);
    let t = f.mul(r, x);
    ts_classic(f, r, t)
}

/// The classic Tonelli–Shanks loop from `r = a^((q+1)/2)`, `t = a^q`.
fn ts_classic<F: PrimeField>(f: &F, mut r: F::Elem, mut t: F::Elem) -> Option<F::Elem> {
    let pre = f.sqrt_pre();
    let mut c = pre.z_q;
    let mut m = pre.s;
    loop {
        if f.is_one(t) {
            return Some(r);
        }
        let mut i = 0u32;
        let mut t2 = t;
        while !f.is_one(t2) {
            t2 = f.sqr(t2);
            i += 1;
            if i == m {
                return None;
            }
        }
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = f.sqr(b);
        }
        m = i;
        c = f.sqr(b);
        t = f.mul(t, c);
        r = f.mul(r, b);
    }
}

/// The Jacobi symbol `(a / n)` for odd `n`, binary algorithm on `u64`
/// (one remainder, then shifts and subtractions only).
pub fn jacobi_u64(a: u64, n: u64) -> i32 {
    debug_assert!(n & 1 == 1);
    let mut a = a % n;
    let mut n = n;
    let mut t = 1i32;
    while a != 0 {
        let z = a.trailing_zeros();
        a >>= z;
        if z & 1 == 1 && (n & 7 == 3 || n & 7 == 5) {
            t = -t;
        }
        if a < n {
            std::mem::swap(&mut a, &mut n);
            if a & 3 == 3 && n & 3 == 3 {
                t = -t;
            }
        }
        a -= n;
    }
    if n == 1 {
        t
    } else {
        0
    }
}

/// Kaliski's almost Montgomery inverse on one word: for odd `p < 2⁶³` and
/// `0 < a < p`, returns `(x, k)` with `x = a⁻¹·2^k mod p` and
/// `k ≤ 2·⌈log₂ p⌉ ≤ 126`; runs of zero bits are shifted out at once.
/// Invariant `p = u·s + v·r`, so `r, s ≤ p` until the last step and
/// `r < 2p < 2⁶⁴` at the end.
#[inline]
fn almost_inv_u64(a: u64, p: u64) -> (u64, u32) {
    debug_assert!(a != 0 && a < p && p & 1 == 1 && p < 1 << 63);
    let (mut u, mut v, mut r, mut s, mut k) = (p, a, 0u64, 1u64, 0u32);
    while v > 0 {
        if u & 1 == 0 {
            let t = u.trailing_zeros();
            u >>= t;
            s <<= t;
            k += t;
        } else if v & 1 == 0 {
            let t = v.trailing_zeros();
            v >>= t;
            r <<= t;
            k += t;
        } else if u > v {
            u = (u - v) >> 1;
            r += s;
            s <<= 1;
            k += 1;
        } else {
            v = (v - u) >> 1;
            s += r;
            r <<= 1;
            k += 1;
        }
    }
    debug_assert_eq!(u, 1, "not invertible");
    if r >= p {
        r -= p;
    }
    (if r == 0 { 0 } else { p - r }, k)
}

/// `table[k] = c·2^(−k) mod p` for `k ≤ 128` (halving mod `p`).
fn halving_table(c: u64, p: u64) -> Vec<u64> {
    let mut t = Vec::with_capacity(129);
    let mut x = c % p;
    for _ in 0..=128 {
        t.push(x);
        x = if x & 1 == 0 {
            x >> 1
        } else {
            ((u128::from(x) + u128::from(p)) >> 1) as u64
        };
    }
    t
}

// ── Fp32: p < 2^32, canonical u32, Barrett ────────────────────────────

/// `F_p` for an odd prime `p < 2³²`: canonical `u32` elements, Barrett
/// reduction with `⌊2⁶⁴/p⌋` (the arithmetic of `f4_fp_tower`).
#[derive(Clone, Debug)]
pub struct Fp32 {
    p: u64,
    /// `⌊2⁶⁴/p⌋`.
    m: u64,
    /// The largest multiple of `p` not above `2⁶³`.
    lazy_m: u64,
    /// `2⁶⁴ mod p` (the carry fold for `p ≥ 2³¹`).
    fold: u64,
    /// `2^(−k) mod p`, the correction of the almost inverse.
    inv_pow2: Vec<u64>,
    /// `p < 2³¹`: the lazy product needs no reduction.
    small: bool,
    lanes: Lanes,
    modulus: BigUint,
    sqrt: SqrtPre<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
enum Lanes {
    Scalar,
    Avx2,
    Avx512,
}

impl Lanes {
    fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512f") {
                return Lanes::Avx512;
            }
            if std::is_x86_feature_detected!("avx2") {
                return Lanes::Avx2;
            }
        }
        Lanes::Scalar
    }
}

/// `acc += a·b` with the `2⁶³` invariant (`p < 2³¹`); vectorises.
#[inline(always)]
fn fp32_lazy_axpy(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    let n = acc.len().min(b.len());
    for (x, &y) in acc[..n].iter_mut().zip(&b[..n]) {
        let s = *x + u64::from(a) * u64::from(y);
        *x = s - ((s >> 63).wrapping_neg() & m);
    }
}

/// `acc += a·b` through Shoup's product (`p < 2³²`), folding a carry out
/// of the word back in as `2⁶⁴ mod p`; vectorises.
#[inline(always)]
fn fp32_shoup_axpy(acc: &mut [u64], a: u32, a_sh: u64, b: &[u32], p: u64, fold: u64) {
    let n = acc.len().min(b.len());
    let a = u64::from(a);
    for (x, &y) in acc[..n].iter_mut().zip(&b[..n]) {
        let y = u64::from(y);
        let q = (a_sh * y) >> 32;
        let r = (a * y).wrapping_sub(q * p);
        let s = x.wrapping_add(r);
        *x = s.wrapping_add(u64::from(s < r).wrapping_neg() & fold);
    }
}

#[inline(always)]
fn barrett_reduce(x: u64, p: u64, m: u64) -> u64 {
    let q = ((u128::from(x) * u128::from(m)) >> 64) as u64;
    let r = x - q * p;
    if r >= p {
        r - p
    } else {
        r
    }
}

#[inline(always)]
fn fp32_reduce_rows(acc: &[u64], out: &mut [u32], p: u64, m: u64) {
    for (o, &w) in out.iter_mut().zip(acc) {
        *o = barrett_reduce(w, p, m) as u32;
    }
}

#[cfg(target_arch = "x86_64")]
mod fp32_simd {
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn lazy_avx2(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
        super::fp32_lazy_axpy(acc, a, b, m)
    }
    #[target_feature(enable = "avx512f")]
    pub(super) unsafe fn lazy_avx512(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
        super::fp32_lazy_axpy(acc, a, b, m)
    }
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn shoup_avx2(
        acc: &mut [u64],
        a: u32,
        a_sh: u64,
        b: &[u32],
        p: u64,
        fold: u64,
    ) {
        super::fp32_shoup_axpy(acc, a, a_sh, b, p, fold)
    }
    #[target_feature(enable = "avx512f")]
    pub(super) unsafe fn shoup_avx512(
        acc: &mut [u64],
        a: u32,
        a_sh: u64,
        b: &[u32],
        p: u64,
        fold: u64,
    ) {
        super::fp32_shoup_axpy(acc, a, a_sh, b, p, fold)
    }
}

impl Fp32 {
    /// The field for an odd prime `3 ≤ p < 2³²` (primality is the
    /// caller's promise).
    pub fn new(p: u64) -> Self {
        assert!(
            p & 1 == 1 && (3..1u64 << 32).contains(&p),
            "Fp32 needs an odd prime below 2^32"
        );
        Self::with_lanes(p, Lanes::detect())
    }

    /// Scalar kernels only (for testing the vector dispatch).
    pub fn new_scalar(p: u64) -> Self {
        let mut f = Self::new(p);
        f.lanes = Lanes::Scalar;
        f
    }

    fn with_lanes(p: u64, lanes: Lanes) -> Self {
        let mut f = Fp32 {
            p,
            m: ((1u128 << 64) / u128::from(p)) as u64,
            lazy_m: ((1u64 << 63) / p) * p,
            fold: ((1u128 << 64) % u128::from(p)) as u64,
            inv_pow2: halving_table(1, p),
            small: p < 1 << 31,
            lanes,
            modulus: BigUint::from(p),
            sqrt: SqrtPre::placeholder(0),
        };
        f.sqrt = compute_sqrt_pre(&f);
        f
    }

    pub fn p(&self) -> u64 {
        self.p
    }

    /// `x mod p` for any `x < 2⁶⁴`.
    #[inline(always)]
    pub fn reduce(&self, x: u64) -> u32 {
        barrett_reduce(x, self.p, self.m) as u32
    }

    #[inline(always)]
    fn renorm(&self, s: u64) -> u64 {
        s - ((s >> 63).wrapping_neg() & self.lazy_m)
    }
}

impl PrimeField for Fp32 {
    type Elem = u32;
    type Wide = u64;

    fn name(&self) -> &'static str {
        "Fp32"
    }
    fn modulus(&self) -> &BigUint {
        &self.modulus
    }
    #[inline(always)]
    fn zero(&self) -> u32 {
        0
    }
    #[inline(always)]
    fn one(&self) -> u32 {
        1
    }
    #[inline(always)]
    fn from_u64(&self, x: u64) -> u32 {
        self.reduce(x)
    }
    fn from_biguint(&self, x: &BigUint) -> u32 {
        (x % self.p).to_u64().unwrap() as u32
    }
    fn to_biguint(&self, a: u32) -> BigUint {
        BigUint::from(a)
    }
    #[inline(always)]
    fn to_u64(&self, a: u32) -> Option<u64> {
        Some(u64::from(a))
    }
    #[inline(always)]
    fn add(&self, a: u32, b: u32) -> u32 {
        let s = u64::from(a) + u64::from(b);
        (if s >= self.p { s - self.p } else { s }) as u32
    }
    #[inline(always)]
    fn sub(&self, a: u32, b: u32) -> u32 {
        if a >= b {
            a - b
        } else {
            (u64::from(a) + self.p - u64::from(b)) as u32
        }
    }
    #[inline(always)]
    fn neg(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (self.p - u64::from(a)) as u32
        }
    }
    #[inline(always)]
    fn mul(&self, a: u32, b: u32) -> u32 {
        self.reduce(u64::from(a) * u64::from(b))
    }
    fn inv(&self, a: u32) -> Option<u32> {
        if a == 0 {
            None
        } else {
            let (x, k) = almost_inv_u64(u64::from(a), self.p);
            Some(self.reduce(x * self.inv_pow2[k as usize]))
        }
    }
    fn sqrt_pre(&self) -> &SqrtPre<u32> {
        &self.sqrt
    }
    fn legendre(&self, a: u32) -> i32 {
        jacobi_u64(u64::from(a), self.p)
    }

    fn lazy_capacity(&self) -> u64 {
        u64::MAX
    }
    #[inline(always)]
    fn wide_zero(&self) -> u64 {
        0
    }
    #[inline(always)]
    fn wide_from(&self, a: u32) -> u64 {
        u64::from(a)
    }
    #[inline(always)]
    fn wide_mul_add(&self, w: &mut u64, a: u32, b: u32) {
        let prod = u64::from(a) * u64::from(b);
        if self.small {
            debug_assert!(*w < 1 << 63, "Fp32 accumulator invariant");
            *w = self.renorm(*w + prod);
        } else {
            // Any u64 is valid; a carry out is 2⁶⁴ ≡ fold. After a carry
            // s < prod < 2⁶⁴ − 2³³, so adding fold cannot carry again.
            let s = w.wrapping_add(prod);
            *w = s.wrapping_add(u64::from(s < prod).wrapping_neg() & self.fold);
        }
    }
    #[inline(always)]
    fn wide_reduce(&self, w: &u64) -> u32 {
        self.reduce(*w)
    }
    #[inline]
    fn row_axpy(&self, acc: &mut [u64], c: u32, row: &[u32]) {
        if self.small {
            debug_assert!(acc.iter().all(|&w| w < 1 << 63));
            match self.lanes {
                #[cfg(target_arch = "x86_64")]
                // SAFETY: `Lanes::detect` found AVX-512F on this CPU.
                Lanes::Avx512 => unsafe { fp32_simd::lazy_avx512(acc, c, row, self.lazy_m) },
                #[cfg(target_arch = "x86_64")]
                // SAFETY: `Lanes::detect` found AVX2 on this CPU.
                Lanes::Avx2 => unsafe { fp32_simd::lazy_avx2(acc, c, row, self.lazy_m) },
                _ => fp32_lazy_axpy(acc, c, row, self.lazy_m),
            }
        } else {
            let c_sh = (u64::from(c) << 32) / self.p;
            match self.lanes {
                #[cfg(target_arch = "x86_64")]
                // SAFETY: as above.
                Lanes::Avx512 => unsafe {
                    fp32_simd::shoup_avx512(acc, c, c_sh, row, self.p, self.fold)
                },
                #[cfg(target_arch = "x86_64")]
                // SAFETY: as above.
                Lanes::Avx2 => unsafe {
                    fp32_simd::shoup_avx2(acc, c, c_sh, row, self.p, self.fold)
                },
                _ => fp32_shoup_axpy(acc, c, c_sh, row, self.p, self.fold),
            }
        }
    }
    #[inline]
    fn row_scatter_axpy(&self, acc: &mut [u64], c: u32, cols: &[u32], vals: &[u32]) {
        let c = u64::from(c);
        if self.small {
            let m = self.lazy_m;
            for (&col, &y) in cols.iter().zip(vals) {
                let w = &mut acc[col as usize];
                let s = *w + c * u64::from(y);
                *w = s - ((s >> 63).wrapping_neg() & m);
            }
        } else {
            let fold = self.fold;
            for (&col, &y) in cols.iter().zip(vals) {
                let w = &mut acc[col as usize];
                let prod = c * u64::from(y);
                let (s, o) = w.overflowing_add(prod);
                *w = if o { s + fold } else { s };
            }
        }
    }
    #[inline]
    fn row_reduce(&self, acc: &[u64], out: &mut [u32]) {
        fp32_reduce_rows(acc, out, self.p, self.m)
    }
    #[inline]
    fn row_axpy_eager(&self, dst: &mut [u32], c: u32, row: &[u32]) {
        // x + c·y ≤ (p − 1) + (p − 1)² < 2⁶⁴: one Barrett per entry.
        let c = u64::from(c);
        for (x, &y) in dst.iter_mut().zip(row) {
            *x = self.reduce(u64::from(*x) + c * u64::from(y));
        }
    }
}

// ── Fp64: p < 2^63, Montgomery u64 ────────────────────────────────────

/// `F_p` for an odd prime `p < 2⁶³`, Montgomery form with `R = 2⁶⁴`.
/// Extends [`crate::cryptanalysis::prime_fast::Fp64`] (`p < 2⁶²`) by one
/// bit and adds the lazy accumulator.
#[derive(Clone, Debug)]
pub struct Fp64 {
    p: u64,
    /// `−p⁻¹ mod 2⁶⁴`.
    n_inv: u64,
    /// `R mod p`.
    one: u64,
    /// `R² mod p`.
    r2: u64,
    /// `R³ mod p`.
    r3: u64,
    /// `p·2⁶⁴`, the accumulator bound.
    p_hi: u128,
    /// `R³·2^(−k) mod p`: one Montgomery product turns the almost inverse
    /// of `a·R` into `a⁻¹·R`.
    inv_fix: Vec<u64>,
    modulus: BigUint,
    sqrt: SqrtPre<u64>,
}

#[inline(always)]
fn neg_inv_u64(p: u64) -> u64 {
    let mut inv: u64 = 1;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
    }
    inv.wrapping_neg()
}

impl Fp64 {
    pub fn new(p: u64) -> Self {
        assert!(
            p & 1 == 1 && (3..1 << 63).contains(&p),
            "Fp64 needs an odd prime below 2^63"
        );
        let r = ((1u128 << 64) % u128::from(p)) as u64;
        let r2 = ((u128::from(r) * u128::from(r)) % u128::from(p)) as u64;
        let mut f = Fp64 {
            p,
            n_inv: neg_inv_u64(p),
            one: r,
            r2,
            r3: 0,
            p_hi: u128::from(p) << 64,
            inv_fix: Vec::new(),
            modulus: BigUint::from(p),
            sqrt: SqrtPre::placeholder(0),
        };
        f.r3 = f.mul(r2, r2);
        // r3 is the integer R³ mod p (mul divides r2·r2 = R⁴ by R).
        f.inv_fix = halving_table(f.r3, p);
        f.sqrt = compute_sqrt_pre(&f);
        f
    }

    pub fn p(&self) -> u64 {
        self.p
    }

    /// Montgomery reduction of `t < p·2⁶⁴`: `t·R⁻¹ mod p`, fully reduced.
    #[inline(always)]
    pub fn redc(&self, t: u128) -> u64 {
        let lo = t as u64;
        let m = lo.wrapping_mul(self.n_inv);
        let mp = u128::from(m) * u128::from(self.p);
        // lo + (mp as u64) ≡ 0 mod 2⁶⁴, carrying exactly when lo ≠ 0.
        let u = ((t >> 64) as u64) + ((mp >> 64) as u64) + u64::from(lo != 0);
        if u >= self.p {
            u - self.p
        } else {
            u
        }
    }

    /// Montgomery form of a canonical `a < p`.
    #[inline(always)]
    pub fn to_mont(&self, a: u64) -> u64 {
        self.redc(u128::from(a) * u128::from(self.r2))
    }

    /// Canonical value of a Montgomery element.
    #[inline(always)]
    pub fn from_mont(&self, a: u64) -> u64 {
        self.redc(u128::from(a))
    }
}

impl PrimeField for Fp64 {
    type Elem = u64;
    type Wide = u128;

    fn name(&self) -> &'static str {
        "Fp64"
    }
    fn modulus(&self) -> &BigUint {
        &self.modulus
    }
    #[inline(always)]
    fn zero(&self) -> u64 {
        0
    }
    #[inline(always)]
    fn one(&self) -> u64 {
        self.one
    }
    #[inline(always)]
    fn from_u64(&self, x: u64) -> u64 {
        // x·R² mod p: x·(R² mod p) < 2⁶⁴·p, inside redc's range.
        self.to_mont(x)
    }
    fn from_biguint(&self, x: &BigUint) -> u64 {
        self.to_mont((x % self.p).to_u64().unwrap())
    }
    fn to_biguint(&self, a: u64) -> BigUint {
        BigUint::from(self.from_mont(a))
    }
    #[inline(always)]
    fn to_u64(&self, a: u64) -> Option<u64> {
        Some(self.from_mont(a))
    }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    #[inline(always)]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.redc(u128::from(a) * u128::from(b))
    }
    fn inv(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return None;
        }
        // (a·R)⁻¹ = a⁻¹·R⁻¹; times R³ through one Montgomery product.
        // (a·R)⁻¹·2^k, times R³·2^(−k) through one Montgomery product.
        let (x, k) = almost_inv_u64(a, self.p);
        Some(self.mul(x, self.inv_fix[k as usize]))
    }
    fn sqrt_pre(&self) -> &SqrtPre<u64> {
        &self.sqrt
    }
    fn legendre(&self, a: u64) -> i32 {
        jacobi_u64(self.from_mont(a), self.p)
    }

    fn lazy_capacity(&self) -> u64 {
        u64::MAX
    }
    #[inline(always)]
    fn wide_zero(&self) -> u128 {
        0
    }
    #[inline(always)]
    fn wide_from(&self, a: u64) -> u128 {
        u128::from(a) << 64
    }
    #[inline(always)]
    fn wide_mul_add(&self, w: &mut u128, a: u64, b: u64) {
        debug_assert!(*w < self.p_hi, "Fp64 accumulator invariant");
        let s = *w + u128::from(a) * u128::from(b);
        *w = s - (self.p_hi & u128::from(s >= self.p_hi).wrapping_neg());
    }
    #[inline(always)]
    fn wide_reduce(&self, w: &u128) -> u64 {
        debug_assert!(*w < self.p_hi);
        self.redc(*w)
    }
    #[inline]
    fn row_axpy(&self, acc: &mut [u128], c: u64, row: &[u64]) {
        let ph = self.p_hi;
        let c = u128::from(c);
        for (w, &y) in acc.iter_mut().zip(row) {
            let s = *w + c * u128::from(y);
            *w = s - (ph & u128::from(s >= ph).wrapping_neg());
        }
    }
    #[inline]
    fn row_axpy_eager(&self, dst: &mut [u64], c: u64, row: &[u64]) {
        // Shoup: with c̄ = from_mont(c), c̄·y ≡ (c·y)·R⁻¹, i.e. the
        // Montgomery product, and c' = ⌊c̄·2⁶⁴/p⌋ gives it in [0, 2p).
        let cb = self.from_mont(c);
        let c_sh = ((u128::from(cb) << 64) / u128::from(self.p)) as u64;
        let p = self.p;
        for (x, &y) in dst.iter_mut().zip(row) {
            let q = ((u128::from(c_sh) * u128::from(y)) >> 64) as u64;
            let r = cb.wrapping_mul(y).wrapping_sub(q.wrapping_mul(p));
            let r = r - (p & u64::from(r >= p).wrapping_neg());
            let s = *x + r;
            *x = s - (p & u64::from(s >= p).wrapping_neg());
        }
    }
}

// ── Limb helpers ──────────────────────────────────────────────────────

#[inline(always)]
fn mac(acc: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = u128::from(acc) + u128::from(a) * u128::from(b) + u128::from(carry);
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = u128::from(a) + u128::from(b) + u128::from(carry);
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let t = u128::from(a)
        .wrapping_sub(u128::from(b))
        .wrapping_sub(u128::from(borrow));
    (t as u64, (t >> 127) as u64)
}

#[inline(always)]
fn limbs_lt<const N: usize>(a: &[u64; N], b: &[u64; N]) -> bool {
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return a[i] < b[i];
        }
    }
    false
}

#[inline(always)]
fn limbs_sub<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], u64) {
    let mut r = [0u64; N];
    let mut br = 0;
    for i in 0..N {
        let (d, b2) = sbb(a[i], b[i], br);
        r[i] = d;
        br = b2;
    }
    (r, br)
}

#[inline(always)]
fn limbs_add<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], u64) {
    let mut r = [0u64; N];
    let mut c = 0;
    for i in 0..N {
        let (s, c2) = adc(a[i], b[i], c);
        r[i] = s;
        c = c2;
    }
    (r, c)
}

#[inline(always)]
fn limbs_is_zero<const N: usize>(a: &[u64; N]) -> bool {
    a.iter().all(|&x| x == 0)
}

/// `t − p` when `extra ≠ 0` or `t ≥ p`, else `t`, without a branch
/// (the final subtraction is a coin flip on random data).
#[inline(always)]
fn limbs_csub<const N: usize>(t: &[u64; N], extra: u64, p: &[u64; N]) -> [u64; N] {
    let (d, borrow) = limbs_sub(t, p);
    // Keep t exactly when the subtraction borrowed and nothing overflowed.
    let keep = (borrow & (extra ^ 1)).wrapping_neg();
    let mut r = [0u64; N];
    for i in 0..N {
        r[i] = (t[i] & keep) | (d[i] & !keep);
    }
    r
}

/// `(a + b) mod p` for `a, b < p`, branch-free.
#[inline(always)]
fn limbs_add_mod<const N: usize>(a: &[u64; N], b: &[u64; N], p: &[u64; N]) -> [u64; N] {
    let (s, c) = limbs_add(a, b);
    limbs_csub(&s, c, p)
}

/// `(a − b) mod p` for `a, b < p`, branch-free.
#[inline(always)]
fn limbs_sub_mod<const N: usize>(a: &[u64; N], b: &[u64; N], p: &[u64; N]) -> [u64; N] {
    let (d, br) = limbs_sub(a, b);
    let mask = br.wrapping_neg();
    let mut pm = [0u64; N];
    for i in 0..N {
        pm[i] = p[i] & mask;
    }
    limbs_add(&d, &pm).0
}

fn limbs_from_biguint<const N: usize>(x: &BigUint) -> [u64; N] {
    let d = x.to_u64_digits();
    assert!(d.len() <= N, "value does not fit {N} limbs");
    let mut r = [0u64; N];
    r[..d.len()].copy_from_slice(&d);
    r
}

fn limbs_to_biguint<const N: usize>(a: &[u64; N]) -> BigUint {
    let mut bytes = Vec::with_capacity(8 * N);
    for w in a {
        bytes.extend_from_slice(&w.to_le_bytes());
    }
    BigUint::from_bytes_le(&bytes)
}

/// Shift right by one, with `top` shifted into the top bit.
#[inline(always)]
fn limbs_shr1<const N: usize>(a: &mut [u64; N], top: u64) {
    for i in 0..N - 1 {
        a[i] = (a[i] >> 1) | (a[i + 1] << 63);
    }
    a[N - 1] = (a[N - 1] >> 1) | (top << 63);
}

/// Variable-time binary extended Euclid: `a⁻¹ mod p` for canonical
/// `0 < a < p`, odd `p` (Guide to ECC, Algorithm 2.22).
fn inv_binary<const N: usize>(a: &[u64; N], p: &[u64; N]) -> [u64; N] {
    let mut u = *a;
    let mut v = *p;
    let mut x1 = [0u64; N];
    x1[0] = 1;
    let mut x2 = [0u64; N];
    let mut one = [0u64; N];
    one[0] = 1;
    let halve = |x: &mut [u64; N]| {
        if x[0] & 1 == 0 {
            limbs_shr1(x, 0);
        } else {
            let (s, c) = limbs_add(x, p);
            *x = s;
            limbs_shr1(x, c);
        }
    };
    while u != one && v != one {
        while u[0] & 1 == 0 {
            limbs_shr1(&mut u, 0);
            halve(&mut x1);
        }
        while v[0] & 1 == 0 {
            limbs_shr1(&mut v, 0);
            halve(&mut x2);
        }
        if !limbs_lt(&u, &v) {
            u = limbs_sub(&u, &v).0;
            let (d, b) = limbs_sub(&x1, &x2);
            x1 = if b != 0 { limbs_add(&d, p).0 } else { d };
        } else {
            v = limbs_sub(&v, &u).0;
            let (d, b) = limbs_sub(&x2, &x1);
            x2 = if b != 0 { limbs_add(&d, p).0 } else { d };
        }
    }
    if u == one {
        x1
    } else {
        x2
    }
}

/// `f·a + g·b` for non-negative `a, b < 2^(64N)` and `|f|, |g| ≤ 2⁶²`:
/// the low `N` limbs and the signed top limb (two's complement).
#[inline(always)]
fn lincomb<const N: usize>(a: &[u64; N], f: i64, b: &[u64; N], g: i64) -> ([u64; N], i64) {
    let mut r = [0u64; N];
    let mut carry: i128 = 0;
    for i in 0..N {
        let t = i128::from(a[i]) * i128::from(f) + i128::from(b[i]) * i128::from(g) + carry;
        r[i] = t as u64;
        carry = t >> 64;
    }
    (r, carry as i64)
}

/// Two's-complement negation of `(lo, top)`.
#[inline(always)]
fn neg_signed<const N: usize>(lo: &mut [u64; N], top: &mut i64) {
    let mut c = 1u64;
    for x in lo.iter_mut() {
        let (s, c2) = adc(!*x, 0, c);
        *x = s;
        c = c2;
    }
    *top = (!*top).wrapping_add(c as i64);
}

/// `(lo, top) >> k` for `0 < k < 64`, arithmetic, keeping `N` limbs and
/// the signed top.
#[inline(always)]
fn shr_signed<const N: usize>(lo: &mut [u64; N], top: &mut i64, k: u32) {
    for i in 0..N - 1 {
        lo[i] = (lo[i] >> k) | (lo[i + 1] << (64 - k));
    }
    lo[N - 1] = (lo[N - 1] >> k) | ((*top as u64) << (64 - k));
    *top >>= k;
}

/// Bit length of `a`.
#[inline(always)]
fn limbs_bits<const N: usize>(a: &[u64; N]) -> u32 {
    for i in (0..N).rev() {
        if a[i] != 0 {
            return 64 * i as u32 + 64 - a[i].leading_zeros();
        }
    }
    0
}

/// Bits `[pos, pos + 33)` of `a` as a `u64` (`pos + 33 ≤ 64N`).
#[inline(always)]
fn limbs_bits_at<const N: usize>(a: &[u64; N], pos: u32) -> u64 {
    let (w, s) = ((pos / 64) as usize, pos % 64);
    let mut v = a[w] >> s;
    if s > 0 && w + 1 < N {
        v |= a[w + 1] << (64 - s);
    }
    v & ((1u64 << 33) - 1)
}

/// `(f·u + g·v) / 2³¹ mod p`, for `u, v < p`: Montgomery division by
/// `2³¹` (`n_inv = −p⁻¹ mod 2⁶⁴`).
#[inline(always)]
fn lincomb_div31<const N: usize>(
    u: &[u64; N],
    f: i64,
    v: &[u64; N],
    g: i64,
    p: &[u64; N],
    n_inv: u64,
) -> [u64; N] {
    let (mut lo, mut top) = lincomb(u, f, v, g);
    // k = −x·p⁻¹ mod 2³¹ makes x + k·p divisible by 2³¹.
    let k = lo[0].wrapping_mul(n_inv) & 0x7fff_ffff;
    let mut c = 0u64;
    for i in 0..N {
        let (s, c2) = mac(lo[i], k, p[i], c);
        lo[i] = s;
        c = c2;
    }
    top = top.wrapping_add(c as i64);
    shr_signed(&mut lo, &mut top, 31);
    // |result| < 4p: fix the sign, then subtract p while ≥ p.
    while top < 0 {
        let (s, c) = limbs_add(&lo, p);
        lo = s;
        top += c as i64;
    }
    while top > 0 || !limbs_lt(&lo, p) {
        let (d, b) = limbs_sub(&lo, p);
        lo = d;
        top -= b as i64;
    }
    lo
}

/// `a⁻¹ mod p` for canonical `0 < a < p`, odd `p`: Pornin's optimised
/// binary GCD (IACR ePrint 2020/972, Algorithm 2) in its variable-time
/// form — 31 divsteps at a time on 64-bit approximations (low 31 bits
/// exact, top 33 bits), then one update of the full values. Falls back to
/// [`inv_binary`] if the iteration bound is ever exceeded.
pub fn inv_bingcd<const N: usize>(y: &[u64; N], p: &[u64; N], n_inv: u64) -> [u64; N] {
    let mut a = *y;
    let mut b = *p;
    let mut u = [0u64; N];
    u[0] = 1;
    let mut v = [0u64; N];
    // Invariants: a ≡ u·y, b ≡ v·y (mod p); b odd.
    let max_outer = (2 * 64 * N as u32).div_ceil(31) + 8;
    for _ in 0..max_outer {
        if limbs_is_zero(&a) {
            return v;
        }
        let n = limbs_bits(&a).max(limbs_bits(&b));
        let (mut xa, mut xb) = if n <= 64 {
            (a[0], b[0])
        } else {
            (
                (a[0] & 0x7fff_ffff) | (limbs_bits_at(&a, n - 33) << 31),
                (b[0] & 0x7fff_ffff) | (limbs_bits_at(&b, n - 33) << 31),
            )
        };
        let (mut f0, mut g0, mut f1, mut g1) = (1i64, 0i64, 0i64, 1i64);
        for _ in 0..31 {
            if xa & 1 == 1 {
                if xa < xb {
                    std::mem::swap(&mut xa, &mut xb);
                    std::mem::swap(&mut f0, &mut f1);
                    std::mem::swap(&mut g0, &mut g1);
                }
                xa -= xb;
                f0 -= f1;
                g0 -= g1;
            }
            xa >>= 1;
            f1 <<= 1;
            g1 <<= 1;
        }
        let (mut na, mut ta) = lincomb(&a, f0, &b, g0);
        let (mut nb, mut tb) = lincomb(&a, f1, &b, g1);
        if ta < 0 {
            neg_signed(&mut na, &mut ta);
            f0 = -f0;
            g0 = -g0;
        }
        if tb < 0 {
            neg_signed(&mut nb, &mut tb);
            f1 = -f1;
            g1 = -g1;
        }
        shr_signed(&mut na, &mut ta, 31);
        shr_signed(&mut nb, &mut tb, 31);
        debug_assert!(ta == 0 && tb == 0);
        let nu = lincomb_div31(&u, f0, &v, g0, p, n_inv);
        let nv = lincomb_div31(&u, f1, &v, g1, p, n_inv);
        a = na;
        b = nb;
        u = nu;
        v = nv;
    }
    if limbs_is_zero(&a) {
        return v;
    }
    #[cfg(test)]
    BINGCD_FALLBACKS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    inv_binary(y, p)
}

#[cfg(test)]
static BINGCD_FALLBACKS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// The Jacobi symbol `(a / n)` for odd `n` on limbs (binary algorithm).
pub fn jacobi_limbs<const N: usize>(a: &[u64; N], n: &[u64; N]) -> i32 {
    let mut a = *a;
    let mut n = *n;
    let mut t = 1i32;
    loop {
        if limbs_is_zero(&a) {
            let mut one = [0u64; N];
            one[0] = 1;
            return if n == one { t } else { 0 };
        }
        // Strip twos: (2/n) = −1 iff n ≡ 3, 5 mod 8.
        let mut z = 0u32;
        while a[0] == 0 {
            a.copy_within(1.., 0);
            a[N - 1] = 0;
            z += 64;
        }
        let tz = a[0].trailing_zeros();
        if tz > 0 {
            for i in 0..N - 1 {
                a[i] = (a[i] >> tz) | (a[i + 1] << (64 - tz));
            }
            a[N - 1] >>= tz;
        }
        z += tz;
        if z & 1 == 1 && (n[0] & 7 == 3 || n[0] & 7 == 5) {
            t = -t;
        }
        // Both odd: make a ≥ n by reciprocity, then a ← a − n.
        if limbs_lt(&a, &n) {
            std::mem::swap(&mut a, &mut n);
            if a[0] & 3 == 3 && n[0] & 3 == 3 {
                t = -t;
            }
        }
        a = limbs_sub(&a, &n).0;
    }
}

/// The largest supported limb count of [`FpWide`].
pub const MAX_LIMBS: usize = 8;

// ── FpWide<N>: Montgomery over N limbs ────────────────────────────────

/// `F_p` for an odd prime `p < 2^(64N)`, `1 ≤ N ≤ 8`: variable-time
/// Montgomery arithmetic and binary-GCD inversion.
///
/// `SPARE = false` (the default) is the general CIOS product with a
/// branch-free final subtraction, for any `p < 2^(64N)`. `SPARE = true`
/// needs a spare top bit (`p[N−1] < 2⁶³ − 1`, asserted): the carry-free
/// CIOS of gnark-crypto with a branch for the final subtraction, which is
/// then rarely taken when `p` is well below `R` (P-224 in four limbs). The
/// choice is a type parameter because a run-time switch inside the
/// product measured 15–40 % slower.
#[derive(Clone, Debug)]
pub struct FpWide<const N: usize, const SPARE: bool = false> {
    p: [u64; N],
    n_inv: u64,
    one: [u64; N],
    r2: [u64; N],
    r3: [u64; N],
    /// `2⁶⁴·R mod p`: undoes the extra word of the lazy reduction.
    k64: [u64; N],
    modulus: BigUint,
    sqrt: SqrtPre<[u64; N]>,
}

/// `FpWide<3>`: primes up to 192 bits (NIST P-192).
pub type Fp192 = FpWide<3>;
/// `FpWide<4>`: primes up to 256 bits (NIST P-256 and any other).
pub type Fp256 = FpWide<4>;
/// `FpWide<4, true>`: primes below `2²⁵⁵` with a spare top bit, the
/// faster choice for NIST P-224 and the 200-bit test prime.
pub type Fp224 = FpWide<4, true>;

/// The lazy accumulator of [`FpWide<N>`]: `lo + hi·R + top·R²`.
#[derive(Clone, Copy, Debug)]
pub struct WideAcc<const N: usize> {
    pub lo: [u64; N],
    pub hi: [u64; N],
    pub top: u64,
}

impl<const N: usize, const SPARE: bool> FpWide<N, SPARE> {
    /// The field for an odd prime `p < 2^(64N)` (primality is the
    /// caller's promise; `p > 3`).
    pub fn new(p: &BigUint) -> Self {
        assert!((1..=MAX_LIMBS).contains(&N));
        assert!(
            p.bit(0) && p > &BigUint::from(3u32),
            "need an odd prime > 3"
        );
        assert!(p.bits() as usize <= 64 * N, "p does not fit {N} limbs");
        let pl: [u64; N] = limbs_from_biguint(p);
        assert!(
            !SPARE || pl[N - 1] < (u64::MAX >> 1) - 1,
            "FpWide<N, true> needs p[N-1] < 2^63 - 1"
        );
        let r = BigUint::one() << (64 * N);
        let one = limbs_from_biguint(&(&r % p));
        let r2 = limbs_from_biguint(&((&r * &r) % p));
        let r3 = limbs_from_biguint(&((&r * &r * &r) % p));
        let k64 = limbs_from_biguint(&(((&r) << 64u32) % p));
        let mut f = FpWide {
            p: pl,
            n_inv: neg_inv_u64(pl[0]),
            one,
            r2,
            r3,
            k64,
            modulus: p.clone(),
            sqrt: SqrtPre::placeholder([0; N]),
        };
        f.sqrt = compute_sqrt_pre(&f);
        f
    }

    /// From a big-endian hexadecimal string.
    pub fn from_hex(hex: &str) -> Self {
        Self::new(&BigUint::parse_bytes(hex.as_bytes(), 16).expect("hex prime"))
    }

    pub fn p_limbs(&self) -> &[u64; N] {
        &self.p
    }

    /// The general CIOS product (`p` may fill the top limb).
    #[inline(always)]
    fn mul_cios(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        let p = &self.p;
        let mut t = [0u64; N];
        let mut t_n = 0u64;
        for &bi in b.iter() {
            let mut c = 0u64;
            for j in 0..N {
                let (lo, hi) = mac(t[j], a[j], bi, c);
                t[j] = lo;
                c = hi;
            }
            let (tn, t_n1) = adc(t_n, c, 0);
            let m = t[0].wrapping_mul(self.n_inv);
            let (_, mut c) = mac(t[0], m, p[0], 0);
            for j in 1..N {
                let (lo, hi) = mac(t[j], m, p[j], c);
                t[j - 1] = lo;
                c = hi;
            }
            let (lo, hi) = adc(tn, c, 0);
            t[N - 1] = lo;
            t_n = t_n1 + hi;
        }
        limbs_csub(&t, t_n, p)
    }

    /// CIOS without the overflow word, for `p[N−1] < 2⁶³ − 1`
    /// (the "no-carry" variant of gnark-crypto).
    #[inline(always)]
    fn mul_nc(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        let p = &self.p;
        let mut t = [0u64; N];
        for &bi in b.iter() {
            let (lo, mut ca) = mac(t[0], a[0], bi, 0);
            t[0] = lo;
            let m = t[0].wrapping_mul(self.n_inv);
            let (_, mut cc) = mac(t[0], m, p[0], 0);
            for j in 1..N {
                let (lo, hi) = mac(t[j], a[j], bi, ca);
                ca = hi;
                let (lo2, hi2) = mac(lo, m, p[j], cc);
                cc = hi2;
                t[j - 1] = lo2;
            }
            t[N - 1] = cc + ca;
        }
        if limbs_lt(&t, p) {
            t
        } else {
            limbs_sub(&t, p).0
        }
    }

    #[inline(always)]
    fn mont_mul(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        if SPARE {
            self.mul_nc(a, b)
        } else {
            self.mul_cios(a, b)
        }
    }

    /// Montgomery reduction of the `2N + 1`-word value `w` over `N + 1`
    /// words: `w·R⁻¹·2⁻⁶⁴ mod p` for `w < p·R·2⁶⁴`, fully reduced.
    #[inline(always)]
    fn redc_wide(&self, w: &WideAcc<N>) -> [u64; N] {
        let p = &self.p;
        let mut t = [0u64; 2 * MAX_LIMBS + 2];
        t[..N].copy_from_slice(&w.lo);
        t[N..2 * N].copy_from_slice(&w.hi);
        t[2 * N] = w.top;
        let mut ex = 0u64;
        for i in 0..=N {
            let m = t[i].wrapping_mul(self.n_inv);
            let mut c = 0u64;
            for j in 0..N {
                let (lo, hi) = mac(t[i + j], m, p[j], c);
                t[i + j] = lo;
                c = hi;
            }
            let (s, c2) = adc(t[i + N], c, ex);
            t[i + N] = s;
            ex = c2;
        }
        let mut r = [0u64; N];
        r.copy_from_slice(&t[N + 1..2 * N + 1]);
        if SPARE {
            if ex == 0 && limbs_lt(&r, p) {
                r
            } else {
                limbs_sub(&r, p).0
            }
        } else {
            limbs_csub(&r, ex, p)
        }
    }

    /// `w < p·R·2⁶⁴`: words `N + 1 ..= 2N` compared with `p`.
    fn wide_in_bound(&self, w: &WideAcc<N>) -> bool {
        let mut hiw = [0u64; N];
        hiw[..N - 1].copy_from_slice(&w.hi[1..]);
        hiw[N - 1] = w.top;
        limbs_lt(&hiw, &self.p)
    }

    /// `a⁻¹ mod p` for a canonical `0 < a < p`.
    pub fn inv_canonical(&self, a: &[u64; N]) -> [u64; N] {
        inv_bingcd(a, &self.p, self.n_inv)
    }

    /// Inverse by the plain binary extended Euclid (comparison arm).
    pub fn inv_plain_binary(&self, a: [u64; N]) -> Option<[u64; N]> {
        if limbs_is_zero(&a) {
            return None;
        }
        Some(self.mont_mul(&inv_binary(&a, &self.p), &self.r3))
    }

    /// The Legendre symbol by Euler's criterion (comparison arm).
    pub fn legendre_euler(&self, a: [u64; N]) -> i32 {
        if limbs_is_zero(&a) {
            return 0;
        }
        let r = self.pow(a, &self.sqrt.p_minus1_half);
        if r == self.one {
            1
        } else {
            -1
        }
    }

    /// Inverse by Fermat's little theorem (comparison arm).
    pub fn inv_fermat(&self, a: [u64; N]) -> Option<[u64; N]> {
        if limbs_is_zero(&a) {
            return None;
        }
        let e = (&self.modulus - BigUint::from(2u32)).to_u64_digits();
        Some(self.pow(a, &e))
    }
}

impl<const N: usize, const SPARE: bool> PrimeField for FpWide<N, SPARE> {
    type Elem = [u64; N];
    type Wide = WideAcc<N>;

    fn name(&self) -> &'static str {
        match (N, SPARE) {
            (3, false) => "FpWide<3>",
            (4, false) => "FpWide<4>",
            (3, true) => "FpWide<3+spare>",
            (4, true) => "FpWide<4+spare>",
            _ => "FpWide<N>",
        }
    }
    fn modulus(&self) -> &BigUint {
        &self.modulus
    }
    #[inline(always)]
    fn zero(&self) -> [u64; N] {
        [0; N]
    }
    #[inline(always)]
    fn one(&self) -> [u64; N] {
        self.one
    }
    fn from_u64(&self, x: u64) -> [u64; N] {
        let mut a = [0u64; N];
        a[0] = x;
        if !limbs_lt(&a, &self.p) {
            a[0] = x % self.p[0];
        }
        self.mont_mul(&a, &self.r2)
    }
    fn from_biguint(&self, x: &BigUint) -> [u64; N] {
        let a = limbs_from_biguint(&(x % &self.modulus));
        self.mont_mul(&a, &self.r2)
    }
    fn to_biguint(&self, a: [u64; N]) -> BigUint {
        let mut one = [0u64; N];
        one[0] = 1;
        limbs_to_biguint(&self.mont_mul(&a, &one))
    }
    #[inline(always)]
    fn is_zero(&self, a: [u64; N]) -> bool {
        limbs_is_zero(&a)
    }
    #[inline(always)]
    fn add(&self, a: [u64; N], b: [u64; N]) -> [u64; N] {
        limbs_add_mod(&a, &b, &self.p)
    }
    #[inline(always)]
    fn sub(&self, a: [u64; N], b: [u64; N]) -> [u64; N] {
        limbs_sub_mod(&a, &b, &self.p)
    }
    #[inline(always)]
    fn neg(&self, a: [u64; N]) -> [u64; N] {
        if limbs_is_zero(&a) {
            a
        } else {
            limbs_sub(&self.p, &a).0
        }
    }
    #[inline(always)]
    fn mul(&self, a: [u64; N], b: [u64; N]) -> [u64; N] {
        self.mont_mul(&a, &b)
    }
    #[inline(always)]
    fn sqr(&self, a: [u64; N]) -> [u64; N] {
        self.mont_mul(&a, &a)
    }
    fn inv(&self, a: [u64; N]) -> Option<[u64; N]> {
        if limbs_is_zero(&a) {
            return None;
        }
        // The inverse of the Montgomery value a·R is a⁻¹·R⁻¹; R³ fixes it.
        let c = inv_bingcd(&a, &self.p, self.n_inv);
        Some(self.mont_mul(&c, &self.r3))
    }
    fn legendre(&self, a: [u64; N]) -> i32 {
        let mut one = [0u64; N];
        one[0] = 1;
        jacobi_limbs(&self.mont_mul(&a, &one), &self.p)
    }
    fn sqrt_pre(&self) -> &SqrtPre<[u64; N]> {
        &self.sqrt
    }

    fn lazy_capacity(&self) -> u64 {
        u64::MAX - 1
    }
    #[inline(always)]
    fn wide_zero(&self) -> WideAcc<N> {
        WideAcc {
            lo: [0; N],
            hi: [0; N],
            top: 0,
        }
    }
    #[inline(always)]
    fn wide_from(&self, a: [u64; N]) -> WideAcc<N> {
        WideAcc {
            lo: [0; N],
            hi: a,
            top: 0,
        }
    }
    #[inline(always)]
    fn wide_mul_add(&self, w: &mut WideAcc<N>, a: [u64; N], b: [u64; N]) {
        // Full product into (pl, ph), operand scanning.
        let mut pl = [0u64; N];
        let mut ph = [0u64; N];
        for i in 0..N {
            let mut c = 0u64;
            for j in 0..N {
                let k = i + j;
                let cur = if k < N { pl[k] } else { ph[k - N] };
                let (lo, hi) = mac(cur, a[j], b[i], c);
                if k < N {
                    pl[k] = lo;
                } else {
                    ph[k - N] = lo;
                }
                c = hi;
            }
            ph[i] = c;
        }
        let mut c = 0u64;
        for k in 0..N {
            let (s, c2) = adc(w.lo[k], pl[k], c);
            w.lo[k] = s;
            c = c2;
        }
        for k in 0..N {
            let (s, c2) = adc(w.hi[k], ph[k], c);
            w.hi[k] = s;
            c = c2;
        }
        debug_assert!(w.top < u64::MAX || c == 0, "FpWide accumulator overflow");
        w.top += c;
    }
    #[inline(always)]
    fn wide_reduce(&self, w: &WideAcc<N>) -> [u64; N] {
        debug_assert!(self.wide_in_bound(w), "FpWide accumulator past capacity");
        let r = self.redc_wide(w);
        self.mont_mul(&r, &self.k64)
    }
    #[inline]
    fn row_reduce_scaled(&self, acc: &[WideAcc<N>], s: [u64; N], out: &mut [[u64; N]]) {
        let sk = self.mont_mul(&s, &self.k64);
        for (o, w) in out.iter_mut().zip(acc) {
            debug_assert!(self.wide_in_bound(w));
            *o = self.mont_mul(&self.redc_wide(w), &sk);
        }
    }
}

// ── FpP192: NIST P-192 with Solinas reduction ─────────────────────────

/// NIST P-192, `p = 2¹⁹² − 2⁶⁴ − 1`, canonical limbs with the Solinas
/// reduction `2¹⁹² ≡ 2⁶⁴ + 1`. Same interface as `FpWide<3>`.
#[derive(Clone, Debug)]
pub struct FpP192 {
    modulus: BigUint,
    sqrt: SqrtPre<[u64; 3]>,
}

const P192: [u64; 3] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFE,
    0xFFFF_FFFF_FFFF_FFFF,
];

/// The 7-word lazy accumulator of [`FpP192`] (`w[6]` counts overflows).
#[derive(Clone, Copy, Debug)]
pub struct P192Acc(pub [u64; 7]);

impl Default for FpP192 {
    fn default() -> Self {
        Self::new()
    }
}

impl FpP192 {
    pub fn new() -> Self {
        let mut f = FpP192 {
            modulus: limbs_to_biguint(&P192),
            sqrt: SqrtPre::placeholder([0; 3]),
        };
        f.sqrt = compute_sqrt_pre(&f);
        f
    }

    /// Reduce `x = Σ x_i 2^(64i)`, six words plus a small seventh, mod p.
    #[inline(always)]
    fn reduce7(x: &[u64; 7]) -> [u64; 3] {
        // 2^192 ≡ 2^64 + 1, 2^256 ≡ 2^128 + 2^64, 2^320 ≡ 2^128 + 2^64 + 1,
        // 2^384 ≡ 2^128 + 2·2^64 + 1 (= (2^64 + 1)²).
        let [x0, x1, x2, x3, x4, x5, x6] = *x;
        // Word sums per position, as u128 accumulators.
        // pos0: x0 + x3 + x5 + x6
        // pos1: x1 + x3 + x4 + x5 + 2·x6
        // pos2: x2 + x4 + x5 + x6
        let s0 = u128::from(x0) + u128::from(x3) + u128::from(x5) + u128::from(x6);
        let s1 = u128::from(x1)
            + u128::from(x3)
            + u128::from(x4)
            + u128::from(x5)
            + 2 * u128::from(x6)
            + (s0 >> 64);
        let s2 = u128::from(x2) + u128::from(x4) + u128::from(x5) + u128::from(x6) + (s1 >> 64);
        // Carry c = s2 >> 64 < 6: c·2^192 ≡ c·(2^64 + 1). Two folds:
        // the first leaves a carry of at most one, the second (adding
        // 2^64 + 1 to a value below 2^67) none.
        // Both branches are rarely taken on random data (p is within
        // 2^65 of 2^192), so they predict well.
        let c = (s2 >> 64) as u64;
        let (a0, k0) = adc(s0 as u64, c, 0);
        let (a1, k1) = adc(s1 as u64, c, k0);
        let (a2, c) = adc(s2 as u64, 0, k1);
        let mut r = [a0, a1, a2];
        if c != 0 {
            let (a0, k0) = adc(r[0], 1, 0);
            let (a1, k1) = adc(r[1], 1, k0);
            let (a2, c2) = adc(r[2], 0, k1);
            debug_assert_eq!(c2, 0);
            r = [a0, a1, a2];
        }
        if limbs_lt(&r, &P192) {
            r
        } else {
            limbs_sub(&r, &P192).0
        }
    }

    #[inline(always)]
    fn mul_wide(a: &[u64; 3], b: &[u64; 3]) -> [u64; 6] {
        let mut t = [0u64; 6];
        for i in 0..3 {
            let mut c = 0u64;
            for j in 0..3 {
                let (lo, hi) = mac(t[i + j], a[j], b[i], c);
                t[i + j] = lo;
                c = hi;
            }
            t[i + 3] = c;
        }
        t
    }
}

impl PrimeField for FpP192 {
    type Elem = [u64; 3];
    type Wide = P192Acc;

    fn name(&self) -> &'static str {
        "FpP192"
    }
    fn modulus(&self) -> &BigUint {
        &self.modulus
    }
    #[inline(always)]
    fn zero(&self) -> [u64; 3] {
        [0; 3]
    }
    #[inline(always)]
    fn one(&self) -> [u64; 3] {
        [1, 0, 0]
    }
    fn from_u64(&self, x: u64) -> [u64; 3] {
        [x, 0, 0]
    }
    fn from_biguint(&self, x: &BigUint) -> [u64; 3] {
        limbs_from_biguint(&(x % &self.modulus))
    }
    fn to_biguint(&self, a: [u64; 3]) -> BigUint {
        limbs_to_biguint(&a)
    }
    fn to_u64(&self, a: [u64; 3]) -> Option<u64> {
        if a[1] == 0 && a[2] == 0 {
            Some(a[0])
        } else {
            None
        }
    }
    #[inline(always)]
    fn is_zero(&self, a: [u64; 3]) -> bool {
        limbs_is_zero(&a)
    }
    #[inline(always)]
    fn add(&self, a: [u64; 3], b: [u64; 3]) -> [u64; 3] {
        limbs_add_mod(&a, &b, &P192)
    }
    #[inline(always)]
    fn sub(&self, a: [u64; 3], b: [u64; 3]) -> [u64; 3] {
        limbs_sub_mod(&a, &b, &P192)
    }
    #[inline(always)]
    fn neg(&self, a: [u64; 3]) -> [u64; 3] {
        if limbs_is_zero(&a) {
            a
        } else {
            limbs_sub(&P192, &a).0
        }
    }
    #[inline(always)]
    fn mul(&self, a: [u64; 3], b: [u64; 3]) -> [u64; 3] {
        let t = Self::mul_wide(&a, &b);
        Self::reduce7(&[t[0], t[1], t[2], t[3], t[4], t[5], 0])
    }
    fn inv(&self, a: [u64; 3]) -> Option<[u64; 3]> {
        if limbs_is_zero(&a) {
            return None;
        }
        Some(inv_bingcd(&a, &P192, 1))
    }
    fn sqrt_pre(&self) -> &SqrtPre<[u64; 3]> {
        &self.sqrt
    }
    fn legendre(&self, a: [u64; 3]) -> i32 {
        jacobi_limbs(&a, &P192)
    }

    fn lazy_capacity(&self) -> u64 {
        u64::MAX - 1
    }
    #[inline(always)]
    fn wide_zero(&self) -> P192Acc {
        P192Acc([0; 7])
    }
    #[inline(always)]
    fn wide_from(&self, a: [u64; 3]) -> P192Acc {
        P192Acc([a[0], a[1], a[2], 0, 0, 0, 0])
    }
    #[inline(always)]
    fn wide_mul_add(&self, w: &mut P192Acc, a: [u64; 3], b: [u64; 3]) {
        let t = Self::mul_wide(&a, &b);
        let mut c = 0u64;
        for (wk, &tk) in w.0[..6].iter_mut().zip(&t) {
            let (s, c2) = adc(*wk, tk, c);
            *wk = s;
            c = c2;
        }
        debug_assert!(
            w.0[6] < u64::MAX - 1 || c == 0,
            "FpP192 accumulator overflow"
        );
        w.0[6] += c;
    }
    #[inline(always)]
    fn wide_reduce(&self, w: &P192Acc) -> [u64; 3] {
        // reduce7 handles x6 up to 2^64 − 1: all sums stay below 2^67.
        Self::reduce7(&w.0)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use num_traits::Zero;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn big(hex: &str) -> BigUint {
        BigUint::parse_bytes(hex.as_bytes(), 16).unwrap()
    }

    pub(super) fn p192() -> BigUint {
        big("fffffffffffffffffffffffffffffffeffffffffffffffff")
    }
    pub(super) fn p224() -> BigUint {
        big("ffffffffffffffffffffffffffffffff000000000000000000000001")
    }
    pub(super) fn p256() -> BigUint {
        big("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff")
    }
    /// A 200-bit prime: the first prime above 2^199 + 0x1234567.
    pub(super) fn p200() -> BigUint {
        prime_above(199, 0x1234567)
    }

    /// The first prime above `2^bits + off`.
    pub(super) fn prime_above(bits: u32, off: u64) -> BigUint {
        let mut c = (BigUint::one() << bits) + BigUint::from(off);
        if !c.bit(0) {
            c += 1u32;
        }
        while !is_probable_prime(&c) {
            c += 2u32;
        }
        c
    }

    fn is_probable_prime(n: &BigUint) -> bool {
        let one = BigUint::one();
        let nm1 = n - &one;
        let s = nm1.trailing_zeros().unwrap();
        let d = &nm1 >> s;
        for a in [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            let mut x = BigUint::from(a).modpow(&d, n);
            if x == one || x == nm1 {
                continue;
            }
            let mut ok = false;
            for _ in 1..s {
                x = (&x * &x) % n;
                if x == nm1 {
                    ok = true;
                    break;
                }
            }
            if !ok {
                return false;
            }
        }
        true
    }

    /// Values every op is checked on: 0, 1, 2, p − 1, p − 2, (p ± 1)/2,
    /// limb boundaries 2^(64k) ± 1 below p, then random.
    fn samples(p: &BigUint, rng: &mut StdRng, n_random: usize) -> Vec<BigUint> {
        let one = BigUint::one();
        let mut v = vec![
            BigUint::zero(),
            one.clone(),
            BigUint::from(2u32),
            p - &one,
            p - BigUint::from(2u32),
            p >> 1u32,
            (p >> 1u32) + &one,
        ];
        let mut k = 32u32;
        while u64::from(k) < p.bits() {
            let b = &one << k;
            v.push(&b - &one);
            v.push(b.clone());
            v.push(&b + &one);
            k += 32;
        }
        for _ in 0..n_random {
            v.push(rng.gen_biguint_below(p));
        }
        v.into_iter().filter(|x| x < p).collect()
    }

    fn check_field<F: PrimeField>(f: &F, seed: u64, n_random: usize) {
        let p = f.modulus().clone();
        let mut rng = StdRng::seed_from_u64(seed);
        let xs = samples(&p, &mut rng, n_random);
        let el: Vec<F::Elem> = xs.iter().map(|x| f.from_biguint(x)).collect();
        for (x, &e) in xs.iter().zip(&el) {
            assert_eq!(&f.to_biguint(e), x, "{} round trip", f.name());
            if let Some(v) = x.to_u64() {
                assert_eq!(f.from_u64(v), e);
                assert_eq!(f.to_u64(e), Some(v));
            }
        }
        assert_eq!(f.to_biguint(f.zero()), BigUint::zero());
        assert_eq!(f.to_biguint(f.one()), BigUint::one());
        assert_eq!(
            f.from_u64(u64::MAX),
            f.from_biguint(&BigUint::from(u64::MAX))
        );
        assert_eq!(f.to_biguint(f.from_i64(-5)), (&p * 5u32 - 5u32) % &p);
        let red = |x: BigUint| x % &p;
        for (i, (x, &a)) in xs.iter().zip(&el).enumerate() {
            assert_eq!(f.to_biguint(f.neg(a)), red(&p - x));
            assert_eq!(f.to_biguint(f.sqr(a)), red(x * x));
            assert_eq!(f.to_biguint(f.double(a)), red(x * 2u32));
            if !x.is_zero() {
                let ia = f.inv(a).unwrap();
                assert_eq!(f.to_biguint(f.mul(a, ia)), BigUint::one(), "inv");
            } else {
                assert!(f.inv(a).is_none());
            }
            // A window of partners, not all pairs.
            for (y, &b) in xs.iter().zip(&el).skip(i).take(12) {
                assert_eq!(f.to_biguint(f.add(a, b)), red(x + y), "add");
                assert_eq!(f.to_biguint(f.sub(a, b)), red(x + &p - y), "sub");
                assert_eq!(f.to_biguint(f.mul(a, b)), red(x * y), "mul");
            }
        }
        // pow against modpow.
        for (x, &a) in xs.iter().zip(&el).take(20) {
            let e: BigUint = rng.gen_biguint(f.bits() as u64 + 7);
            assert_eq!(f.to_biguint(f.pow_biguint(a, &e)), x.modpow(&e, &p), "pow");
            let e64 = rng.gen::<u64>();
            assert_eq!(
                f.to_biguint(f.pow_u64(a, e64)),
                x.modpow(&BigUint::from(e64), &p)
            );
        }
        // Legendre and sqrt.
        let half = (&p - 1u32) >> 1u32;
        for (x, &a) in xs.iter().zip(&el) {
            let euler = x.modpow(&half, &p);
            let want = if x.is_zero() {
                0
            } else if euler.is_one() {
                1
            } else {
                -1
            };
            assert_eq!(f.legendre(a), want, "legendre");
            match f.sqrt(a) {
                Some(r) => {
                    assert!(want >= 0, "sqrt of a non-residue");
                    assert_eq!(f.sqr(r), a, "sqrt");
                }
                None => assert_eq!(want, -1, "sqrt missed a residue"),
            }
            // Squares always have a root.
            let s = f.sqr(a);
            let r = f.sqrt(s).expect("square");
            assert_eq!(f.sqr(r), s);
        }
        let z = f.from_u64(f.sqrt_pre().non_residue);
        assert_eq!(f.legendre(z), -1);
        assert!(f.sqrt(z).is_none());
        // Batch inverse, zeros left in place.
        let mut v = el.clone();
        let mut scratch = Vec::new();
        f.batch_inv(&mut v, &mut scratch);
        for (&a, &b) in el.iter().zip(&v) {
            if f.is_zero(a) {
                assert!(f.is_zero(b));
            } else {
                assert_eq!(Some(b), f.inv(a));
            }
        }
        check_lazy(f, &mut rng);
    }

    fn check_lazy<F: PrimeField>(f: &F, rng: &mut StdRng) {
        let p = f.modulus().clone();
        let n = 300;
        let rand_el = |rng: &mut StdRng| f.from_biguint(&rng.gen_biguint_below(&p));
        let pm1 = f.from_biguint(&(&p - 1u32));
        // Worst case: every entry p − 1.
        let worst = vec![pm1; n];
        let a: Vec<F::Elem> = (0..n).map(|_| rand_el(rng)).collect();
        let b: Vec<F::Elem> = (0..n).map(|_| rand_el(rng)).collect();
        for (u, v) in [(&a, &b), (&worst, &worst), (&a, &worst)] {
            let naive = u
                .iter()
                .zip(v.iter())
                .fold(f.zero(), |s, (&x, &y)| f.add(s, f.mul(x, y)));
            assert_eq!(f.dot(u, v), naive, "{} dot", f.name());
        }
        // Row ops: load, many axpy/sub with p − 1 and random scalars, reduce.
        let start: Vec<F::Elem> = (0..n).map(|_| rand_el(rng)).collect();
        let mut acc = vec![f.wide_zero(); n];
        f.row_load(&start, &mut acc);
        let mut eager = start.clone();
        for k in 0..200 {
            let c = if k % 3 == 0 { pm1 } else { rand_el(rng) };
            let row: Vec<F::Elem> = if k % 2 == 0 {
                worst.clone()
            } else {
                (0..n).map(|_| rand_el(rng)).collect()
            };
            if k % 5 == 0 {
                f.row_sub(&mut acc, c, &row);
                let nc = f.neg(c);
                for (x, &y) in eager.iter_mut().zip(&row) {
                    *x = f.add(*x, f.mul(nc, y));
                }
            } else {
                f.row_axpy(&mut acc, c, &row);
                let mut chk = eager.clone();
                f.row_axpy_eager(&mut chk, c, &row);
                for (x, &y) in eager.iter_mut().zip(&row) {
                    *x = f.add(*x, f.mul(c, y));
                }
                assert_eq!(chk, eager, "{} eager axpy", f.name());
            }
        }
        let mut out = vec![f.zero(); n];
        f.row_reduce(&acc, &mut out);
        assert_eq!(out, eager, "{} lazy rows", f.name());
        let s = rand_el(rng);
        f.row_reduce_scaled(&acc, s, &mut out);
        let want: Vec<F::Elem> = eager.iter().map(|&x| f.mul(x, s)).collect();
        assert_eq!(out, want, "{} scaled reduce", f.name());
        // Sparse scatter rows against the eager reference.
        let mut sacc = vec![f.wide_zero(); n];
        f.row_load(&start, &mut sacc);
        let mut se = start.clone();
        for _ in 0..300 {
            let cols: Vec<u32> = (0..20).map(|_| rng.gen_range(0..n as u32)).collect();
            let vals: Vec<F::Elem> = (0..20).map(|_| rand_el(rng)).collect();
            let c = if rng.gen::<bool>() { pm1 } else { rand_el(rng) };
            f.row_scatter_axpy(&mut sacc, c, &cols, &vals);
            for (&col, &y) in cols.iter().zip(&vals) {
                se[col as usize] = f.add(se[col as usize], f.mul(c, y));
            }
        }
        let mut so = vec![f.zero(); n];
        f.row_reduce(&sacc, &mut so);
        assert_eq!(so, se, "{} scatter", f.name());
        // Scalar accumulator, scattered like the sparse F4 kernel.
        let mut w = f.wide_from(start[0]);
        let mut e = start[0];
        for _ in 0..1000 {
            let (x, y) = (rand_el(rng), rand_el(rng));
            f.wide_mul_add(&mut w, x, y);
            e = f.add(e, f.mul(x, y));
        }
        assert_eq!(f.wide_reduce(&w), e);
        // LazyRow with offsets.
        let mut lr = LazyRow::from_row(f, &start);
        let mut e2 = start.clone();
        for k in 0..50 {
            let off = k % 7;
            let c = rand_el(rng);
            let row: Vec<F::Elem> = (0..n - off).map(|_| rand_el(rng)).collect();
            lr.sub(f, off, c, &row);
            for (x, &y) in e2[off..].iter_mut().zip(&row) {
                *x = f.sub(*x, f.mul(c, y));
            }
        }
        let mut o = Vec::new();
        lr.finish(f, &mut o);
        assert_eq!(o, e2);
    }

    const SMALL: [u64; 5] = [3, 5, 7, 13, 17];
    // 31-bit and 32-bit primes, including 2^31 − 1 (≡ 3 mod 4), 2013265921
    // (= 15·2^27 + 1, v2 = 27) and the largest prime below 2^32.
    const P32: [u64; 5] = [
        2_147_483_647,
        2_013_265_921,
        1_000_000_007,
        2_147_483_659,
        4_294_967_291,
    ];
    // 61/62/63-bit: 2^61 − 1, 2^62 − 57, 2^63 − 25, and two primes
    // c·2^40 + 1 (62 and 63 bits, v2(p − 1) = 40, long Tonelli–Shanks).
    const P64: [u64; 5] = [
        (1 << 61) - 1,
        (1 << 62) - 57,
        (1 << 63) - 25,
        0x2000_4600_0000_0001,
        0x4000_1200_0000_0001,
    ];

    #[test]
    fn fp32_matches_biguint() {
        for &p in SMALL.iter().chain(&P32) {
            check_field(&Fp32::new(p), p, 120);
            check_field(&Fp32::new_scalar(p), p ^ 1, 40);
        }
    }

    #[test]
    fn fp64_matches_biguint() {
        for &p in SMALL.iter().chain(&P32).chain(&P64) {
            assert!(crate::cryptanalysis::residual_walk::is_prime_u64(p), "{p}");
            check_field(&Fp64::new(p), p, 120);
        }
    }

    #[test]
    fn fpwide_matches_biguint() {
        // P-192 fills three limbs; a 190-bit prime takes the carry-free path.
        for p in [p192(), prime_above(189, 0xabcdef)] {
            check_field(&FpWide::<3>::new(&p), 3, 150);
        }
        check_field(&FpWide::<3, true>::new(&prime_above(189, 0xabcdef)), 3, 150);
        for p in [p200(), p224()] {
            check_field(&Fp224::new(&p), 4, 150);
        }
        for p in [p192(), p200(), p224(), p256()] {
            check_field(&FpWide::<4>::new(&p), 4, 150);
        }
        // Small and one-limb primes through the wide code.
        for p in [
            5u64,
            2_147_483_647,
            (1 << 61) - 1,
            18_446_744_073_709_551_557,
        ] {
            check_field(&FpWide::<1>::new(&BigUint::from(p)), p, 60);
            check_field(&FpWide::<2>::new(&BigUint::from(p)), p, 60);
        }
        // A 127-bit and a full 128-bit prime in two limbs.
        let m127 = (BigUint::one() << 127u32) - 1u32;
        check_field(&FpWide::<2>::new(&m127), 2, 100);
        let p128 = (BigUint::one() << 128u32) - 159u32;
        assert!(is_probable_prime(&p128));
        check_field(&FpWide::<2>::new(&p128), 2, 100);
    }

    #[test]
    fn p192_solinas_matches_biguint() {
        let f = FpP192::new();
        assert_eq!(f.modulus(), &p192());
        check_field(&f, 192, 300);
        let g = Fp192::new(&p192());
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..2000 {
            let x = rng.gen_biguint_below(&p192());
            let y = rng.gen_biguint_below(&p192());
            let a = f.mul(f.from_biguint(&x), f.from_biguint(&y));
            let b = g.mul(g.from_biguint(&x), g.from_biguint(&y));
            assert_eq!(f.to_biguint(a), g.to_biguint(b));
        }
    }

    #[test]
    fn variant_paths_agree() {
        // The carry-free CIOS is used for P-224 in 4 limbs; the general one
        // for P-192 and P-256. Force both on one prime.
        let mut rng = StdRng::seed_from_u64(77);
        for p in [p224(), p200()] {
            let f = Fp224::new(&p);
            for _ in 0..2000 {
                let a = f.from_biguint(&rng.gen_biguint_below(&p));
                let b = f.from_biguint(&rng.gen_biguint_below(&p));
                assert_eq!(f.mul_nc(&a, &b), f.mul_cios(&a, &b));
            }
            let a = f.from_biguint(&rng.gen_biguint_below(&p));
            assert_eq!(f.inv(a), f.inv_fermat(a));
        }
        // Big-p Fp32 Shoup rows against scalar dispatch.
        let p = 4_294_967_291u64;
        let (fv, fs) = (Fp32::new(p), Fp32::new_scalar(p));
        let row: Vec<u32> = (0..1000).map(|_| rng.gen_range(0..p as u32)).collect();
        let mut a1 = vec![0u64; 1000];
        let mut a2 = vec![0u64; 1000];
        for k in 0..100 {
            let c = (rng.gen::<u32>() as u64 % p) as u32;
            fv.row_axpy(&mut a1, c, &row[k..]);
            fs.row_axpy(&mut a2, c, &row[k..]);
        }
        let (mut o1, mut o2) = (vec![0; 1000], vec![0; 1000]);
        fv.row_reduce(&a1, &mut o1);
        fs.row_reduce(&a2, &mut o2);
        assert_eq!(o1, o2);
    }

    #[test]
    fn two_adic_primes_take_tonelli_shanks() {
        // P-224: v2(p − 1) = 96, the worst case for Tonelli–Shanks.
        let f = Fp256::new(&p224());
        assert_eq!(f.sqrt_pre().s, 96);
        let mut rng = StdRng::seed_from_u64(5);
        let mut seen_nr = 0;
        for _ in 0..40 {
            let a = f.from_biguint(&rng.gen_biguint_below(&p224()));
            match f.sqrt(a) {
                Some(r) => assert_eq!(f.sqr(r), a),
                None => {
                    assert_eq!(f.legendre(a), -1);
                    seen_nr += 1;
                }
            }
        }
        assert!(seen_nr > 5);
        for _ in 0..40 {
            let a = f.sqr(f.from_biguint(&rng.gen_biguint_below(&p224())));
            let r1 = f.sqrt(a).unwrap();
            let r2 = sqrt_tonelli_shanks_plain(&f, a).unwrap();
            assert!(r1 == r2 || r1 == f.neg(r2));
        }
        // P-192: p ≡ 3 mod 4 fast path.
        assert!(Fp192::new(&p192()).sqrt_pre().p_plus1_quarter.is_some());
        assert_eq!(Fp192::new(&p192()).sqrt_pre().s, 1);
    }

    #[test]
    fn bingcd_and_jacobi_match_reference() {
        let mut rng = StdRng::seed_from_u64(31);
        fn run<const N: usize>(p: &BigUint, rng: &mut StdRng) {
            let f = FpWide::<N>::new(p);
            let pl = *f.p_limbs();
            let one = BigUint::one();
            let mut xs: Vec<BigUint> = vec![one.clone(), BigUint::from(2u32), p - &one, p - 2u32];
            for k in 1..(64 * N as u32) {
                let v = &one << k;
                if &v < p {
                    xs.push(v.clone());
                    xs.push(&v - &one);
                }
            }
            for _ in 0..3000 {
                xs.push(rng.gen_biguint_below(p));
            }
            let pm2 = p - 2u32;
            let half = (p - &one) >> 1u32;
            for x in xs.iter().filter(|x| !x.is_zero()) {
                let a: [u64; N] = limbs_from_biguint(x);
                let inv = limbs_to_biguint(&inv_bingcd(&a, &pl, neg_inv_u64(pl[0])));
                assert_eq!(inv, x.modpow(&pm2, p), "bingcd {x:x} mod {p:x}");
                let e = x.modpow(&half, p);
                let want = if e.is_one() { 1 } else { -1 };
                assert_eq!(jacobi_limbs(&a, &pl), want, "jacobi {x:x}");
            }
        }
        run::<3>(&p192(), &mut rng);
        run::<3>(&prime_above(189, 0xabcdef), &mut rng);
        run::<4>(&p224(), &mut rng);
        run::<4>(&p256(), &mut rng);
        run::<4>(&p200(), &mut rng);
        run::<2>(&((BigUint::one() << 127u32) - 1u32), &mut rng);
        run::<1>(&BigUint::from((1u64 << 61) - 1), &mut rng);
        run::<8>(&prime_above(511, 0x1f), &mut rng);
        assert_eq!(
            BINGCD_FALLBACKS.load(std::sync::atomic::Ordering::Relaxed),
            0,
            "the binary-GCD iteration bound was reached"
        );
    }

    #[test]
    fn jacobi_matches_euler() {
        let mut rng = StdRng::seed_from_u64(1);
        for &p in P32.iter().chain(&P64) {
            for _ in 0..200 {
                let a = rng.gen::<u64>() % p;
                let e = BigUint::from(a).modpow(&BigUint::from((p - 1) / 2), &BigUint::from(p));
                let want = if a == 0 {
                    0
                } else if e.is_one() {
                    1
                } else {
                    -1
                };
                assert_eq!(jacobi_u64(a, p), want);
            }
        }
    }
}
