//! # Relation filtering and block Wiedemann for the Koblitz logarithm precompute.
//!
//! The factor-base logarithm database
//! ([`super::koblitz_index_calculus::solve_factor_base_logs`]) is the
//! solution of one linear system over `Z/rZ`:
//!
//! ```text
//!     Σ_o c_{i,o} · x_o  ≡  h·a_i   (mod r)         one row per relation,
//! ```
//!
//! with `x_o` the unknown logarithm of projected column `o`.  Each row
//! has at most `m` nonzero entries (`m` summands per decomposition), so
//! for a factor base with thousands of columns the matrix is extremely
//! sparse and dense Gaussian elimination — `O(rows · cols²)` big-integer
//! operations, repeated on every attempt — is the wrong tool.  This
//! module gives the ECDLP pipeline the two stages every number-field
//! sieve driver puts between relation collection and the logarithm
//! solve:
//!
//! 1. **Filtering** ([`filter_relations`]): duplicate rows are dropped;
//!    a column that occurs in exactly one row (a *singleton*) is
//!    eliminated together with that row and recovered afterwards by
//!    back-substitution; surplus rows beyond a target excess are
//!    removed, choosing rows whose removal cascades through the
//!    weight-2 columns (the clique rule); and light columns are merged
//!    away by structured Gaussian elimination under a fill-in bound.
//!    Every elimination is recorded so the eliminated unknowns are
//!    reconstructed exactly from the core solution.
//! 2. **Block Wiedemann** ([`block_wiedemann_kernel`]): the reduced core
//!    is made square by folding its excess rows into random earlier
//!    rows, homogenised to `M (x, 1)ᵀ = 0`, and a kernel vector is found
//!    from the Krylov sequence `X Mⁱ Y` (`m × n` blocks) through a
//!    shifted minimal approximant basis — the matrix Berlekamp–Massey
//!    step of Coppersmith's algorithm.  Only sparse matrix-times-block
//!    products touch the matrix, and they run in parallel over rows.
//!
//! The solution is checked against every original row before it is
//! returned, and the caller certifies each logarithm in the group
//! (`[x_o]G == R_o`), so a wrong answer cannot escape.  The dense solver
//! in [`super::ec_index_calculus::gaussian_eliminate_mod_n`] remains the
//! reference the tests cross-check against.
//!
//! ## Word widths
//!
//! Everything is generic over a modular word type ([`ModWord`]).  The
//! `u64` instance is the original path for subgroup orders below `2⁶³`
//! (`u128` products, Shoup multipliers, lazily reduced accumulators) and
//! is unchanged.  [`Wide<L>`] (`L = 2, 3, 4` limbs: [`U128`], [`U192`],
//! [`U256`]) covers odd moduli below `2^(64L−1)`, which includes the
//! 129-bit ECC2K-130 subgroup order; it keeps canonical residues in the
//! public row and vector types and uses Montgomery products internally:
//! a matrix coefficient is stored as `a·R mod m`, so a row of the sparse
//! product is one exact multi-limb sum of products reduced once by REDC,
//! and a fixed multiplier (a filtering pivot factor, an approximant-basis
//! update) is one Montgomery product per element.  [`SparseRows`] picks
//! the width for a modulus given as a `BigUint`, so the index-calculus
//! driver never names a width.
//!
//! ## References
//!
//! - D. Wiedemann, *Solving sparse linear equations over finite fields*,
//!   IEEE Trans. Inf. Theory 32 (1986).
//! - D. Coppersmith, *Solving homogeneous linear equations over GF(2)
//!   via block Wiedemann algorithm*, Math. Comp. 62 (1994).
//! - E. Kaltofen, *Analysis of Coppersmith's block Wiedemann algorithm
//!   for the parallel solution of sparse linear systems*, Math. Comp. 64
//!   (1995).
//! - B. Beckermann, G. Labahn, *A uniform approach for the fast
//!   computation of matrix-type Padé approximants*, SIAM J. Matrix Anal.
//!   Appl. 15 (1994); P. Giorgi, C.-P. Jeannerod, G. Villard, *On the
//!   complexity of polynomial matrix computations*, ISSAC 2003 (the
//!   iterative order basis used for the matrix Berlekamp–Massey step).
//! - C. Bouvier, *The filtering step of discrete logarithm and integer
//!   factorization algorithms*, 2013 (singleton, clique and merge
//!   rules as run by CADO-NFS).
//! - P. Montgomery, *Modular multiplication without trial division*,
//!   Math. Comp. 44 (1985) (the wide word type's products).

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;
use std::io::Write;
use std::process::Command;

// ── Modular arithmetic over Z/rZ with r < 2^63 ─────────────────────

#[inline]
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

#[inline]
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    let s = a + b;
    if s >= m {
        s - m
    } else {
        s
    }
}

#[inline]
fn submod(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + (m - b)
    }
}

#[inline]
fn negmod(a: u64, m: u64) -> u64 {
    if a == 0 {
        0
    } else {
        m - a
    }
}

/// Modular inverse by the extended Euclidean algorithm; `None` when
/// `a` and `m` are not coprime (impossible for prime `m` and `a ≠ 0`).
fn invmod(a: u64, m: u64) -> Option<u64> {
    // The remainders never go negative, so they and the quotients stay
    // in `u64`: a hardware division per step instead of the `i128`
    // library call (`__divti3`), which filtering pays once per pivot.
    let (mut old_r, mut r) = (a, m);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q as i128 * s);
    }
    if old_r != 1 {
        return None;
    }
    Some(old_s.rem_euclid(m as i128) as u64)
}

fn to_u64_mod(v: &BigUint, m: u64) -> u64 {
    let reduced = v % BigUint::from(m);
    reduced.to_u64_digits().first().copied().unwrap_or(0)
}

/// The largest modulus the `u64` arithmetic supports: sums of two
/// residues must not overflow.
pub const U64_MODULUS_BITS: u64 = 63;

/// The largest modulus any word width here supports (`U256`, one bit of
/// headroom for the lazily reduced sums).
pub const MAX_MODULUS_BITS: u64 = 4 * 64 - 1;

/// The word width a modulus is handled in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum WordWidth {
    U64,
    U128,
    U192,
    U256,
}

impl WordWidth {
    /// The narrowest width that handles `r`, or `None` when no width
    /// does (`r < 2`, `r` past [`MAX_MODULUS_BITS`], or an even `r`
    /// above the `u64` range, which Montgomery reduction excludes).
    pub fn for_modulus(r: &BigUint) -> Option<Self> {
        if *r < BigUint::from(2u32) {
            return None;
        }
        let bits = r.bits();
        if bits <= U64_MODULUS_BITS {
            return Some(Self::U64);
        }
        if r.is_even() {
            return None;
        }
        match bits {
            0..=127 => Some(Self::U128),
            128..=191 => Some(Self::U192),
            192..=MAX_MODULUS_BITS => Some(Self::U256),
            _ => None,
        }
    }

    pub fn bits(self) -> u32 {
        match self {
            Self::U64 => 64,
            Self::U128 => 128,
            Self::U192 => 192,
            Self::U256 => 256,
        }
    }
}

/// Whether a subgroup order can be handled by this module.
pub fn modulus_supported(r: &BigUint) -> bool {
    WordWidth::for_modulus(r).is_some()
}

// ── The modular word abstraction ───────────────────────────────────

/// A residue type over `Z/mZ` with the operations the filtering and
/// Krylov code need.  Values are canonical residues (`0 ≤ v < m`)
/// wherever they are visible outside this module; the *lazy form*
/// ([`ModWord::lazy_form`]) is the representation a matrix coefficient
/// is stored in for the sparse products and never leaves a
/// [`CsrMatrix`].
pub trait ModWord: Copy + Eq + Hash + Debug + Default + Send + Sync + 'static {
    /// The modulus with whatever constants the arithmetic precomputes.
    type Modulus: Copy + Debug + Send + Sync + 'static;
    /// A multiplier prepared for many products by the same residue.
    type Fixed: Copy + Send + Sync;
    /// An exact accumulator for sums of products.
    type Acc: Copy;
    /// The context of the lazily reduced sums ([`LazyLanes`]), including
    /// how many products an accumulator absorbs between reductions.
    type Lazy: Copy + Send + Sync;

    /// The width this type is.
    const WIDTH: WordWidth;

    /// The modulus for `r`, or `None` when this width cannot hold it.
    fn modulus_from_biguint(r: &BigUint) -> Option<Self::Modulus>;
    fn modulus_to_biguint(m: &Self::Modulus) -> BigUint;

    fn from_biguint(v: &BigUint, m: &Self::Modulus) -> Self;
    fn to_biguint(self) -> BigUint;
    /// `v mod m`.
    fn from_u64(v: u64, m: &Self::Modulus) -> Self;
    /// The canonical residue of any stored value.
    fn reduce(self, m: &Self::Modulus) -> Self;

    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(self) -> bool;
    /// `m − 1`.
    fn max_residue(m: &Self::Modulus) -> Self;
    /// The larger of two stored values, as integers.
    fn larger(a: Self, b: Self) -> Self;

    fn add(self, b: Self, m: &Self::Modulus) -> Self;
    fn sub(self, b: Self, m: &Self::Modulus) -> Self;
    fn neg(self, m: &Self::Modulus) -> Self;
    fn mul(self, b: Self, m: &Self::Modulus) -> Self;
    /// `None` when `self` and `m` are not coprime.
    fn inv(self, m: &Self::Modulus) -> Option<Self>;

    /// Uniform in `[0, m)`.
    fn random<R: Rng + ?Sized>(rng: &mut R, m: &Self::Modulus) -> Self;
    /// Uniform in `[1, m)`.
    fn random_nonzero<R: Rng + ?Sized>(rng: &mut R, m: &Self::Modulus) -> Self;

    /// Prepare `w` (a residue) for repeated products.
    fn fixed(w: Self, m: &Self::Modulus) -> Self::Fixed;
    /// `w · t`.
    fn fixed_mul(f: &Self::Fixed, t: Self, m: &Self::Modulus) -> Self;
    /// `a − w · t` for a residue `a`.
    fn fixed_sub_from(f: &Self::Fixed, a: Self, t: Self, m: &Self::Modulus) -> Self;

    /// The lazy context for sums whose coefficients are at most `a_max`
    /// (in lazy form, products of the [`ModWord::lazy_form`] of a
    /// coefficient with a residue); `None` when the width has no lazy
    /// path for `m` (the `u64` path past `2⁶³`), in which case
    /// coefficients stay canonical and products are reduced one by one.
    fn lazy(m: &Self::Modulus, a_max: Self) -> Option<Self::Lazy>;
    /// Products an accumulator absorbs between reductions (`≥ 1`).
    fn lazy_terms(l: &Self::Lazy) -> usize;
    /// The stored form of a coefficient for the lazy products.
    fn lazy_form(a: Self, m: &Self::Modulus) -> Self;
    /// The canonical residue a lazy-form coefficient stands for.
    fn from_lazy_form(a: Self, m: &Self::Modulus) -> Self;
    /// `a · x` for a lazy-form `a` and a residue `x`.
    fn lazy_mul(a: Self, x: Self, m: &Self::Modulus) -> Self;

    fn acc_zero() -> Self::Acc;
    /// `acc += a · x`, exactly.
    fn acc_mac(acc: &mut Self::Acc, a: Self, x: Self);
    /// The residue of a sum of lazy-form coefficients times residues.
    fn acc_reduce(acc: Self::Acc, l: &Self::Lazy) -> Self;
    /// The residue of a sum of products of residues (both canonical).
    fn acc_reduce_plain(acc: Self::Acc, l: &Self::Lazy) -> Self;
    /// The accumulator standing for residue `r` in a sum of lazy-form
    /// products, so that [`ModWord::acc_reduce`] of it plus further
    /// products is `r` plus theirs.  It counts as one residue (`< m`)
    /// towards the accumulator bound, like a fresh product's worth.
    fn acc_from_residue(r: Self, l: &Self::Lazy) -> Self::Acc;
    /// The same for a sum of plain products ([`ModWord::acc_reduce_plain`]).
    fn acc_from_residue_plain(r: Self, l: &Self::Lazy) -> Self::Acc;
}

// ── Division-free products by a fixed modulus ──────────────────────
//
// `mulmod` divides a `u128` by the modulus (`__umodti3`, tens of cycles
// on x86-64), and that division was most of the block Wiedemann time:
// once per matrix entry and lane in the block products, and once per
// element update in the approximant basis.  The helpers below return
// the same canonical residue from multiplications alone, so every
// vector, generator and report is unchanged.  They need `m < 2⁶³`
// (`U64_MODULUS_BITS`); each caller keeps the plain `mulmod` path for a
// modulus outside that range.

/// `x mod m` for `x < 2m`: `x − m` wraps past `x` exactly when `x < m`,
/// so the smaller of the two is the residue.  The comparison is a coin
/// flip on the data; written as `min` it compiles to a conditional move,
/// where the `if`s of `submod` became branches in the approximant basis
/// whose mispredictions cost more than the division these helpers save.
#[inline(always)]
fn fold_once(x: u64, m: u64) -> u64 {
    x.min(x.wrapping_sub(m))
}

/// `w·t mod m` for a fixed `w < m < 2⁶³` and any `t`, by Shoup's method
/// with `ws = ⌊w·2⁶⁴/m⌋`.  Writing `ws = (w·2⁶⁴ − e)/m` with `0 ≤ e < m`,
/// `ws·t/2⁶⁴ = w·t/m − e·t/(m·2⁶⁴)` falls short of `w·t/m` by less than
/// one, so `q = ⌊ws·t/2⁶⁴⌋` is `⌊w·t/m⌋` or one less and `w·t − q·m`
/// lies in `[0, 2m)`, which fits a `u64` because `m < 2⁶³`; one
/// conditional subtraction makes it canonical.
#[inline]
fn mulmod_shoup(w: u64, ws: u64, t: u64, m: u64) -> u64 {
    let q = ((ws as u128 * t as u128) >> 64) as u64;
    fold_once(w.wrapping_mul(t).wrapping_sub(q.wrapping_mul(m)), m)
}

/// Multiplication by one fixed residue: Shoup's method when `m < 2⁶³`,
/// the `u128` division otherwise.  Its set-up costs one division, so it
/// pays wherever the same factor scales more than a couple of entries.
#[derive(Clone, Copy, Debug)]
pub struct FixedMul {
    w: u64,
    /// `⌊w·2⁶⁴/m⌋`, or `None` when the modulus is too wide for it.
    ws: Option<u64>,
}

impl FixedMul {
    #[inline]
    fn new(w: u64, m: u64) -> Self {
        debug_assert!(w < m, "a fixed multiplier must be reduced");
        let ws = (m >> U64_MODULUS_BITS == 0).then(|| (((w as u128) << 64) / m as u128) as u64);
        Self { w, ws }
    }

    #[inline]
    fn mul(self, t: u64, m: u64) -> u64 {
        match self.ws {
            Some(ws) => mulmod_shoup(self.w, ws, t, m),
            None => mulmod(self.w, t, m),
        }
    }

    /// `a − w·t mod m` for a residue `a`: `submod(a, mulmod(w, t, m), m)`
    /// with the subtraction as a conditional move too
    /// (`a − p + m ∈ [1, 2m)` for residues `a`, `p`).
    #[inline]
    fn sub_from(self, a: u64, t: u64, m: u64) -> u64 {
        match self.ws {
            Some(ws) => fold_once(
                a.wrapping_sub(mulmod_shoup(self.w, ws, t, m))
                    .wrapping_add(m),
                m,
            ),
            None => submod(a, mulmod(self.w, t, m), m),
        }
    }
}

/// Reduction of a `u128` modulo a fixed `m` with `2 ≤ m < 2⁶³`:
/// `v = hi·2⁶⁴ + lo ≡ hi·(2⁶⁴ mod m) + lo·1`, both products by Shoup's
/// method, so five multiplications replace the `u128` division.
#[derive(Clone, Copy, Debug)]
pub struct Reducer {
    m: u64,
    /// `2⁶⁴ mod m` and its Shoup multiplier.
    r64: u64,
    r64_shoup: u64,
    /// `⌊2⁶⁴/m⌋`, the Shoup multiplier of `1`.
    one_shoup: u64,
}

impl Reducer {
    /// `None` outside `2 ≤ m < 2⁶³`, where the Shoup bound fails (or,
    /// for `m = 1`, the multiplier of `1` does not fit); the callers keep
    /// their per-entry `mulmod` path for such moduli.
    fn new(m: u64) -> Option<Self> {
        if m < 2 || m >> U64_MODULUS_BITS != 0 {
            return None;
        }
        let shoup = |w: u64| (((w as u128) << 64) / m as u128) as u64;
        let r64 = ((1u128 << 64) % m as u128) as u64;
        Some(Self {
            m,
            r64,
            r64_shoup: shoup(r64),
            one_shoup: shoup(1),
        })
    }

    #[inline]
    fn reduce(&self, v: u128) -> u64 {
        let m = self.m;
        let high = mulmod_shoup(self.r64, self.r64_shoup, (v >> 64) as u64, m);
        let low = mulmod_shoup(1, self.one_shoup, v as u64, m);
        fold_once(high + low, m)
    }

    /// How many products `a·t` with `a ≤ a_max` and any `t < 2⁶⁴` a
    /// `u128` accumulator holding a residue (`≤ m − 1`) absorbs before
    /// it must be reduced again: the largest `K` with
    /// `(m − 1) + K·a_max·(2⁶⁴ − 1) < 2¹²⁸`.  At least 1 for any `a_max`,
    /// since `(2⁶⁴ − 1)² + 2⁶³ < 2¹²⁸`; at least 2 for a reduced
    /// `a_max < 2⁶³`; and over 800,000 for the 45-bit `k0n53` order, far
    /// above any row weight, so there the reduction happens once a lane.
    fn lazy_terms(&self, a_max: u64) -> usize {
        let a_max = a_max.max(1) as u128;
        let room = u128::MAX - (self.m - 1) as u128;
        (room / (a_max * u64::MAX as u128)).min(usize::MAX as u128) as usize
    }
}

impl ModWord for u64 {
    type Modulus = u64;
    type Fixed = FixedMul;
    type Acc = u128;
    type Lazy = (Reducer, usize);

    const WIDTH: WordWidth = WordWidth::U64;

    fn modulus_from_biguint(r: &BigUint) -> Option<u64> {
        (r.bits() <= U64_MODULUS_BITS && *r >= BigUint::from(2u32))
            .then(|| r.to_u64_digits().first().copied().unwrap_or(0))
    }
    fn modulus_to_biguint(m: &u64) -> BigUint {
        BigUint::from(*m)
    }
    fn from_biguint(v: &BigUint, m: &u64) -> u64 {
        to_u64_mod(v, *m)
    }
    fn to_biguint(self) -> BigUint {
        BigUint::from(self)
    }
    #[inline]
    fn from_u64(v: u64, m: &u64) -> u64 {
        v % *m
    }
    #[inline]
    fn reduce(self, m: &u64) -> u64 {
        self % *m
    }
    #[inline]
    fn zero() -> u64 {
        0
    }
    #[inline]
    fn one() -> u64 {
        1
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
    #[inline]
    fn max_residue(m: &u64) -> u64 {
        m.saturating_sub(1)
    }
    #[inline]
    fn larger(a: u64, b: u64) -> u64 {
        a.max(b)
    }
    #[inline]
    fn add(self, b: u64, m: &u64) -> u64 {
        addmod(self, b, *m)
    }
    #[inline]
    fn sub(self, b: u64, m: &u64) -> u64 {
        submod(self, b, *m)
    }
    #[inline]
    fn neg(self, m: &u64) -> u64 {
        negmod(self, *m)
    }
    #[inline]
    fn mul(self, b: u64, m: &u64) -> u64 {
        mulmod(self, b, *m)
    }
    fn inv(self, m: &u64) -> Option<u64> {
        invmod(self, *m)
    }
    #[inline]
    fn random<R: Rng + ?Sized>(rng: &mut R, m: &u64) -> u64 {
        rng.gen_range(0..*m)
    }
    #[inline]
    fn random_nonzero<R: Rng + ?Sized>(rng: &mut R, m: &u64) -> u64 {
        rng.gen_range(1..*m)
    }
    #[inline]
    fn fixed(w: u64, m: &u64) -> FixedMul {
        FixedMul::new(w, *m)
    }
    #[inline]
    fn fixed_mul(f: &FixedMul, t: u64, m: &u64) -> u64 {
        f.mul(t, *m)
    }
    #[inline]
    fn fixed_sub_from(f: &FixedMul, a: u64, t: u64, m: &u64) -> u64 {
        f.sub_from(a, t, *m)
    }
    fn lazy(m: &u64, a_max: u64) -> Option<(Reducer, usize)> {
        Reducer::new(*m).map(|r| (r, r.lazy_terms(a_max)))
    }
    #[inline(always)]
    fn lazy_terms(l: &(Reducer, usize)) -> usize {
        l.1
    }
    #[inline(always)]
    fn lazy_form(a: u64, _m: &u64) -> u64 {
        a
    }
    #[inline(always)]
    fn from_lazy_form(a: u64, _m: &u64) -> u64 {
        a
    }
    #[inline]
    fn lazy_mul(a: u64, x: u64, m: &u64) -> u64 {
        mulmod(a, x, *m)
    }
    #[inline(always)]
    fn acc_zero() -> u128 {
        0
    }
    #[inline(always)]
    fn acc_mac(acc: &mut u128, a: u64, x: u64) {
        *acc += a as u128 * x as u128;
    }
    #[inline(always)]
    fn acc_reduce(acc: u128, l: &(Reducer, usize)) -> u64 {
        l.0.reduce(acc)
    }
    #[inline(always)]
    fn acc_reduce_plain(acc: u128, l: &(Reducer, usize)) -> u64 {
        l.0.reduce(acc)
    }
    #[inline(always)]
    fn acc_from_residue(r: u64, _l: &(Reducer, usize)) -> u128 {
        r as u128
    }
    #[inline(always)]
    fn acc_from_residue_plain(r: u64, _l: &(Reducer, usize)) -> u128 {
        r as u128
    }
}

// ── Wide residues: `L` limbs, Montgomery products ──────────────────

/// The most limbs a [`Wide`] may have (the exact accumulators are sized
/// for twice this).
pub const MAX_LIMBS: usize = 4;

/// A canonical residue of `L` little-endian 64-bit limbs.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Wide<const L: usize>(pub [u64; L]);

impl<const L: usize> Default for Wide<L> {
    fn default() -> Self {
        Self([0; L])
    }
}

/// Residues modulo an odd `m < 2^(64L−1)`, for `L` up to four limbs.
pub type U128 = Wide<2>;
pub type U192 = Wide<3>;
pub type U256 = Wide<4>;

impl<const L: usize> Debug for Wide<L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_biguint())
    }
}

impl<const L: usize> Wide<L> {
    const ZERO: Self = Self([0; L]);

    fn from_biguint_unchecked(v: &BigUint) -> Self {
        let mut limbs = [0u64; L];
        for (slot, digit) in limbs.iter_mut().zip(v.iter_u64_digits()) {
            *slot = digit;
        }
        Self(limbs)
    }

    fn to_biguint(self) -> BigUint {
        let mut digits: Vec<u32> = Vec::with_capacity(2 * L);
        for limb in self.0 {
            digits.push(limb as u32);
            digits.push((limb >> 32) as u32);
        }
        BigUint::new(digits)
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self.0.iter().all(|&l| l == 0)
    }

    /// `self < other` as integers.
    #[inline(always)]
    fn lt(&self, other: &Self) -> bool {
        for i in (0..L).rev() {
            if self.0[i] != other.0[i] {
                return self.0[i] < other.0[i];
            }
        }
        false
    }

    /// `self + b` (no overflow for two residues of an `m < 2^(64L−1)`).
    #[inline(always)]
    fn add_raw(self, b: Self) -> Self {
        let mut out = [0u64; L];
        let mut carry = 0u64;
        for i in 0..L {
            let (s, c1) = self.0[i].overflowing_add(b.0[i]);
            let (s, c2) = s.overflowing_add(carry);
            out[i] = s;
            carry = u64::from(c1) + u64::from(c2);
        }
        Self(out)
    }

    /// `self − b` for `b ≤ self`.
    #[inline(always)]
    fn sub_raw(self, b: Self) -> Self {
        let mut out = [0u64; L];
        let mut borrow = 0u64;
        for i in 0..L {
            let (d, b1) = self.0[i].overflowing_sub(b.0[i]);
            let (d, b2) = d.overflowing_sub(borrow);
            out[i] = d;
            borrow = u64::from(b1) + u64::from(b2);
        }
        Self(out)
    }

    /// `self − m` if `self ≥ m`, else `self`.
    #[inline(always)]
    fn fold(self, m: &Self) -> Self {
        if self.lt(m) {
            self
        } else {
            self.sub_raw(*m)
        }
    }
}

/// An odd modulus for [`Wide<L>`] with its Montgomery constants
/// (`R = 2^(64L)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WideModulus<const L: usize> {
    m: Wide<L>,
    /// `−m⁻¹ mod 2⁶⁴`.
    neg_inv: u64,
    /// `R² mod m`: the Montgomery form of `R`, so `REDC(a · r2) = a·R`.
    r2: Wide<L>,
    /// `R mod m`, the Montgomery form of `1`.
    r1: Wide<L>,
    bits: u64,
    /// `⌊(R − 1) / m⌋`: products of a residue by a lazy-form coefficient
    /// (both `< m`) an accumulator takes while staying below `m·R`,
    /// which REDC needs.  At least 2 because `m < R/2`.
    lazy_terms: usize,
}

impl<const L: usize> WideModulus<L> {
    /// `None` unless `m` is odd, at least 3, and below `2^(64L−1)`.
    pub fn new(m: &BigUint) -> Option<Self> {
        assert!(
            (1..=MAX_LIMBS).contains(&L),
            "Wide<{L}>: between one and {MAX_LIMBS} limbs"
        );
        if m.is_even() || *m < BigUint::from(3u32) || m.bits() > 64 * L as u64 - 1 {
            return None;
        }
        let r = BigUint::one() << (64 * L);
        let m0 = m.iter_u64_digits().next().unwrap_or(1);
        // Newton iteration for m⁻¹ mod 2⁶⁴ (each step doubles the bits).
        let mut inv = m0;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(m0.wrapping_mul(inv)));
        }
        debug_assert_eq!(m0.wrapping_mul(inv), 1);
        let r1 = &r % m;
        let r2 = (&r1 * &r1) % m;
        let lazy_terms = ((&r - BigUint::one()) / m)
            .to_u64_digits()
            .first()
            .copied()
            .map_or(usize::MAX, |k| k.min(usize::MAX as u64) as usize);
        Some(Self {
            m: Wide::from_biguint_unchecked(m),
            neg_inv: inv.wrapping_neg(),
            r2: Wide::from_biguint_unchecked(&r2),
            r1: Wide::from_biguint_unchecked(&r1),
            bits: m.bits(),
            lazy_terms,
        })
    }

    pub fn to_biguint(&self) -> BigUint {
        self.m.to_biguint()
    }

    pub fn bits(&self) -> u64 {
        self.bits
    }

    /// `t · R⁻¹ mod m` for `t < m·R` held in `2L` limbs (the top limbs
    /// of `t` beyond `2L` are zero), canonical.
    #[inline(always)]
    fn redc(&self, mut t: [u64; 2 * MAX_LIMBS]) -> Wide<L> {
        let m = &self.m.0;
        for i in 0..L {
            let u = t[i].wrapping_mul(self.neg_inv);
            let mut carry = 0u128;
            for j in 0..L {
                let s = t[i + j] as u128 + u as u128 * m[j] as u128 + carry;
                t[i + j] = s as u64;
                carry = s >> 64;
            }
            let mut k = i + L;
            while carry != 0 {
                // `t + u·m·2^(64i) < 2mR < R²` for `m < R/2`, so the
                // carry dies inside the `2L` limbs.
                let s = t[k] as u128 + carry;
                t[k] = s as u64;
                carry = s >> 64;
                k += 1;
            }
        }
        debug_assert!(
            t[2 * L..].iter().all(|&v| v == 0),
            "REDC input exceeded m·R"
        );
        let mut out = [0u64; L];
        out.copy_from_slice(&t[L..2 * L]);
        Wide(out).fold(&self.m)
    }

    /// `a · b · R⁻¹ mod m` for `a < R`, `b < m` (either order).
    #[inline(always)]
    fn mont_mul(&self, a: &Wide<L>, b: &Wide<L>) -> Wide<L> {
        let mut t = [0u64; 2 * MAX_LIMBS];
        Wide::<L>::mac_into(&mut t, a, b);
        self.redc(t)
    }
}

impl<const L: usize> Wide<L> {
    /// `acc += a · b`, exactly, in `2L` limbs (the caller bounds the sum
    /// below `R²`).
    #[inline(always)]
    fn mac_into(acc: &mut [u64; 2 * MAX_LIMBS], a: &Self, b: &Self) {
        for i in 0..L {
            let ai = a.0[i] as u128;
            if ai == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..L {
                let s = acc[i + j] as u128 + ai * b.0[j] as u128 + carry;
                acc[i + j] = s as u64;
                carry = s >> 64;
            }
            let mut k = i + L;
            while carry != 0 {
                let s = acc[k] as u128 + carry;
                acc[k] = s as u64;
                carry = s >> 64;
                k += 1;
            }
        }
    }
}

impl<const L: usize> ModWord for Wide<L> {
    type Modulus = WideModulus<L>;
    type Fixed = Wide<L>;
    type Acc = [u64; 2 * MAX_LIMBS];
    type Lazy = (WideModulus<L>, usize);

    const WIDTH: WordWidth = match L {
        2 => WordWidth::U128,
        3 => WordWidth::U192,
        _ => WordWidth::U256,
    };

    fn modulus_from_biguint(r: &BigUint) -> Option<WideModulus<L>> {
        WideModulus::new(r)
    }
    fn modulus_to_biguint(m: &WideModulus<L>) -> BigUint {
        m.to_biguint()
    }
    fn from_biguint(v: &BigUint, m: &WideModulus<L>) -> Self {
        Self::from_biguint_unchecked(&(v % m.to_biguint()))
    }
    fn to_biguint(self) -> BigUint {
        Wide::to_biguint(self)
    }
    fn from_u64(v: u64, m: &WideModulus<L>) -> Self {
        let mut limbs = [0u64; L];
        limbs[0] = v;
        Self(limbs).reduce(m)
    }
    fn reduce(self, m: &WideModulus<L>) -> Self {
        if self.lt(&m.m) {
            self
        } else {
            Self::from_biguint(&self.to_biguint(), m)
        }
    }
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }
    #[inline(always)]
    fn one() -> Self {
        let mut limbs = [0u64; L];
        limbs[0] = 1;
        Self(limbs)
    }
    #[inline(always)]
    fn is_zero(self) -> bool {
        Wide::is_zero(self)
    }
    fn max_residue(m: &WideModulus<L>) -> Self {
        m.m.sub_raw(Self::one())
    }
    fn larger(a: Self, b: Self) -> Self {
        if a.lt(&b) {
            b
        } else {
            a
        }
    }
    #[inline(always)]
    fn add(self, b: Self, m: &WideModulus<L>) -> Self {
        self.add_raw(b).fold(&m.m)
    }
    #[inline(always)]
    fn sub(self, b: Self, m: &WideModulus<L>) -> Self {
        if self.lt(&b) {
            self.add_raw(m.m).sub_raw(b)
        } else {
            self.sub_raw(b)
        }
    }
    #[inline(always)]
    fn neg(self, m: &WideModulus<L>) -> Self {
        if Wide::is_zero(self) {
            self
        } else {
            m.m.sub_raw(self)
        }
    }
    #[inline]
    fn mul(self, b: Self, m: &WideModulus<L>) -> Self {
        // (a·R) · b · R⁻¹ = a·b.
        m.mont_mul(&m.mont_mul(&self, &m.r2), &b)
    }
    fn inv(self, m: &WideModulus<L>) -> Option<Self> {
        self.to_biguint()
            .modinv(&m.to_biguint())
            .map(|v| Self::from_biguint_unchecked(&v))
    }
    fn random<R: Rng + ?Sized>(rng: &mut R, m: &WideModulus<L>) -> Self {
        // Rejection sampling on `bits` random bits: fewer than two draws
        // on average, uniform exactly.
        let top = ((m.bits - 1) / 64) as usize;
        let mask = u64::MAX >> ((64 - m.bits % 64) % 64);
        loop {
            let mut limbs = [0u64; L];
            for (i, limb) in limbs.iter_mut().enumerate().take(top + 1) {
                *limb = rng.gen();
                if i == top {
                    *limb &= mask;
                }
            }
            let v = Self(limbs);
            if v.lt(&m.m) {
                return v;
            }
        }
    }
    fn random_nonzero<R: Rng + ?Sized>(rng: &mut R, m: &WideModulus<L>) -> Self {
        loop {
            let v = Self::random(rng, m);
            if !Wide::is_zero(v) {
                return v;
            }
        }
    }
    #[inline]
    fn fixed(w: Self, m: &WideModulus<L>) -> Self {
        m.mont_mul(&w, &m.r2)
    }
    #[inline(always)]
    fn fixed_mul(f: &Self, t: Self, m: &WideModulus<L>) -> Self {
        m.mont_mul(f, &t)
    }
    #[inline(always)]
    fn fixed_sub_from(f: &Self, a: Self, t: Self, m: &WideModulus<L>) -> Self {
        a.sub(m.mont_mul(f, &t), m)
    }
    fn lazy(m: &WideModulus<L>, _a_max: Self) -> Option<(WideModulus<L>, usize)> {
        // Coefficients are canonicalised by `lazy_form`, so the bound
        // depends on the modulus alone.
        Some((*m, m.lazy_terms))
    }
    #[inline(always)]
    fn lazy_terms(l: &(WideModulus<L>, usize)) -> usize {
        l.1
    }
    #[inline]
    fn lazy_form(a: Self, m: &WideModulus<L>) -> Self {
        m.mont_mul(&a, &m.r2)
    }
    #[inline]
    fn from_lazy_form(a: Self, m: &WideModulus<L>) -> Self {
        m.mont_mul(&a, &Self::one())
    }
    #[inline(always)]
    fn lazy_mul(a: Self, x: Self, m: &WideModulus<L>) -> Self {
        m.mont_mul(&a, &x)
    }
    #[inline(always)]
    fn acc_zero() -> [u64; 2 * MAX_LIMBS] {
        [0; 2 * MAX_LIMBS]
    }
    #[inline(always)]
    fn acc_mac(acc: &mut [u64; 2 * MAX_LIMBS], a: Self, x: Self) {
        Self::mac_into(acc, &a, &x);
    }
    #[inline(always)]
    fn acc_reduce(acc: [u64; 2 * MAX_LIMBS], l: &(WideModulus<L>, usize)) -> Self {
        l.0.redc(acc)
    }
    #[inline(always)]
    fn acc_reduce_plain(acc: [u64; 2 * MAX_LIMBS], l: &(WideModulus<L>, usize)) -> Self {
        // Σ a·x · R⁻¹, then times R² · R⁻¹.
        let m = &l.0;
        m.mont_mul(&m.redc(acc), &m.r2)
    }
    #[inline(always)]
    fn acc_from_residue(r: Self, l: &(WideModulus<L>, usize)) -> [u64; 2 * MAX_LIMBS] {
        // A lazy sum `T` stands for `T·R⁻¹`, so `r` re-enters as
        // `r·R mod m`, which is below `m` like any single product term.
        let mut acc = [0u64; 2 * MAX_LIMBS];
        acc[..L].copy_from_slice(&l.0.mont_mul(&r, &l.0.r2).0);
        acc
    }
    #[inline(always)]
    fn acc_from_residue_plain(r: Self, _l: &(WideModulus<L>, usize)) -> [u64; 2 * MAX_LIMBS] {
        let mut acc = [0u64; 2 * MAX_LIMBS];
        acc[..L].copy_from_slice(&r.0);
        acc
    }
}

/// `W` lanes of a sum of products `a·t`, each summed exactly in the
/// word's accumulator and reduced once at the end, instead of one
/// division (or Montgomery reduction) per product.  A lane is folded
/// back to a residue whenever `terms` products
/// ([`ModWord::lazy_terms`]) have accumulated, which only a modulus near
/// the top of its width ever needs; the slice methods split their input
/// at those points so that the inner loop is only multiply-and-add.
/// The sum is exact and the reduction canonical, so the result is the
/// per-product `add(mul(..))` chain's to the bit.
///
/// Coefficients (`a`) are in lazy form and the data (`src`, `x`)
/// residues, unless the lanes were opened with [`LazyLanes::plain`],
/// where both are residues.
struct LazyLanes<'r, W: ModWord, const N: usize> {
    acc: [W::Acc; N],
    /// Products the lanes can still absorb before a fold.
    room: usize,
    terms: usize,
    plain: bool,
    lazy: &'r W::Lazy,
}

impl<'r, W: ModWord, const N: usize> LazyLanes<'r, W, N> {
    #[inline(always)]
    fn new(lazy: &'r W::Lazy) -> Self {
        let terms = W::lazy_terms(lazy);
        Self {
            acc: [W::acc_zero(); N],
            room: terms,
            terms,
            plain: false,
            lazy,
        }
    }

    /// Lanes for sums of products of two residues.
    #[inline(always)]
    fn plain(lazy: &'r W::Lazy) -> Self {
        Self {
            plain: true,
            ..Self::new(lazy)
        }
    }

    #[inline(always)]
    fn mac(&mut self, a: W, src: &[W]) {
        for (acc, &s) in self.acc.iter_mut().zip(&src[..N]) {
            W::acc_mac(acc, a, s);
        }
    }

    #[inline(always)]
    fn fold(&mut self) {
        let lazy = self.lazy;
        if self.plain {
            self.acc
                .iter_mut()
                .for_each(|v| *v = W::acc_from_residue_plain(W::acc_reduce_plain(*v, lazy), lazy));
        } else {
            self.acc
                .iter_mut()
                .for_each(|v| *v = W::acc_from_residue(W::acc_reduce(*v, lazy), lazy));
        }
        self.room = self.terms;
    }

    /// Lane `l` gains `a · src[l]`.
    #[inline(always)]
    fn add(&mut self, a: W, src: &[W]) {
        if self.room == 0 {
            self.fold();
        }
        self.room -= 1;
        self.mac(a, src);
    }

    /// Lane `l` gains `Σ_k vals[k] · x[cols[k]·n + l]`.
    #[inline(always)]
    fn add_sparse(&mut self, cols: &[u32], vals: &[W], x: &[W], n: usize) {
        let (mut cols, mut vals) = (cols, vals);
        loop {
            let take = vals.len().min(self.room);
            for (&c, &a) in cols[..take].iter().zip(&vals[..take]) {
                self.mac(a, &x[c as usize * n..]);
            }
            self.room -= take;
            (cols, vals) = (&cols[take..], &vals[take..]);
            if vals.is_empty() {
                break;
            }
            self.fold();
        }
    }

    /// Lane `l` gains `Σ_k coeffs[k] · x[k·stride + l]`.
    #[inline(always)]
    fn add_strided(&mut self, coeffs: &[W], x: &[W], stride: usize) {
        let (mut coeffs, mut x) = (coeffs, x);
        loop {
            let take = coeffs.len().min(self.room);
            for (k, &a) in coeffs[..take].iter().enumerate() {
                self.mac(a, &x[k * stride..]);
            }
            self.room -= take;
            coeffs = &coeffs[take..];
            if coeffs.is_empty() {
                break;
            }
            x = &x[take * stride..];
            self.fold();
        }
    }

    #[inline(always)]
    fn finish(self, out: &mut [W]) {
        if self.plain {
            for (o, &a) in out.iter_mut().zip(&self.acc) {
                *o = W::acc_reduce_plain(a, self.lazy);
            }
        } else {
            for (o, &a) in out.iter_mut().zip(&self.acc) {
                *o = W::acc_reduce(a, self.lazy);
            }
        }
    }
}

// ── Sparse rows ────────────────────────────────────────────────────

/// One relation row `Σ c · x_col ≡ rhs (mod r)`: `(column, coefficient)`
/// pairs sorted by column, coefficients nonzero and reduced.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SparseRow<W: ModWord = u64> {
    pub entries: Vec<(u32, W)>,
    pub rhs: W,
}

impl<W: ModWord> SparseRow<W> {
    /// Build a row from unsorted, possibly repeated, unreduced entries.
    pub fn new(mut entries: Vec<(u32, W)>, rhs: W, modulus: W::Modulus) -> Self {
        let m = &modulus;
        entries.sort_by_key(|&(c, _)| c);
        let mut merged: Vec<(u32, W)> = Vec::with_capacity(entries.len());
        for (c, v) in entries {
            let v = v.reduce(m);
            match merged.last_mut() {
                Some((lc, lv)) if *lc == c => *lv = lv.add(v, m),
                _ => merged.push((c, v)),
            }
        }
        merged.retain(|&(_, v)| !v.is_zero());
        Self {
            entries: merged,
            rhs: rhs.reduce(m),
        }
    }

    /// A row from the dense big-integer form the relation collector
    /// produces.
    pub fn from_dense(row: &[BigUint], rhs: &BigUint, modulus: W::Modulus) -> Self {
        let m = &modulus;
        let entries = row
            .iter()
            .enumerate()
            .filter_map(|(c, v)| {
                let v = W::from_biguint(v, m);
                (!v.is_zero()).then_some((c as u32, v))
            })
            .collect();
        Self {
            entries,
            rhs: W::from_biguint(rhs, m),
        }
    }

    /// Number of nonzero coefficients.
    pub fn weight(&self) -> usize {
        self.entries.len()
    }

    fn coefficient(&self, col: u32) -> Option<W> {
        self.entries
            .binary_search_by_key(&col, |&(c, _)| c)
            .ok()
            .map(|i| self.entries[i].1)
    }

    /// `Σ c · x_col mod r` for a full assignment `x`.
    fn evaluate(&self, x: &[W], modulus: W::Modulus) -> W {
        let m = &modulus;
        self.entries.iter().fold(W::zero(), |acc, &(c, v)| {
            acc.add(v.mul(x[c as usize], m), m)
        })
    }

    /// `self − factor · other`, entries merged by column (`factor < r`).
    fn axpy(&self, factor: W, other: &SparseRow<W>, modulus: W::Modulus) -> SparseRow<W> {
        let m = &modulus;
        let factor = W::fixed(factor, m);
        let mut entries = Vec::with_capacity(self.entries.len() + other.entries.len());
        let (mut i, mut j) = (0, 0);
        while i < self.entries.len() || j < other.entries.len() {
            let take_self = j >= other.entries.len()
                || (i < self.entries.len() && self.entries[i].0 < other.entries[j].0);
            let take_other = i >= self.entries.len()
                || (j < other.entries.len() && other.entries[j].0 < self.entries[i].0);
            if take_self {
                entries.push(self.entries[i]);
                i += 1;
            } else if take_other {
                let (c, v) = other.entries[j];
                entries.push((c, W::fixed_mul(&factor, v, m).neg(m)));
                j += 1;
            } else {
                let (c, a) = self.entries[i];
                let b = other.entries[j].1;
                let v = W::fixed_sub_from(&factor, a, b, m);
                if !v.is_zero() {
                    entries.push((c, v));
                }
                i += 1;
                j += 1;
            }
        }
        SparseRow {
            entries,
            rhs: W::fixed_sub_from(&factor, self.rhs, other.rhs, m),
        }
    }
}

// ── Filtering ──────────────────────────────────────────────────────

/// Controls for [`filter_relations`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FilterOptions {
    /// Drop rows identical to an earlier one.
    pub remove_duplicates: bool,
    /// Eliminate singleton columns with their rows (recovered by
    /// back-substitution).
    pub remove_singletons: bool,
    /// Remove surplus rows until `rows − columns` is at most this.
    /// `usize::MAX` keeps every row.
    pub target_excess: usize,
    /// Merge away columns of weight at most this by structured Gaussian
    /// elimination (`1` restricts the stage to singletons, `0` disables
    /// it).
    pub merge_max_weight: usize,
    /// A merge is skipped when it would leave a row heavier than this.
    pub max_row_weight: usize,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            remove_duplicates: true,
            remove_singletons: true,
            target_excess: 32,
            merge_max_weight: 8,
            max_row_weight: 32,
        }
    }
}

/// What filtering did to the system.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct FilterReport {
    pub rows_in: usize,
    pub columns_in: usize,
    pub nonzeros_in: usize,
    pub duplicates_removed: usize,
    pub singletons_removed: usize,
    pub excess_rows_removed: usize,
    pub merged_columns: usize,
    /// Rows that became `0 ≡ 0` during merging (linear combinations of
    /// others) and were dropped.
    pub dependent_rows_dropped: usize,
    /// Columns left in no row (never determined by the core).
    pub uncovered_columns: usize,
    pub rows_out: usize,
    pub columns_out: usize,
    pub nonzeros_out: usize,
}

/// The reduced system plus everything needed to undo the reduction.
#[derive(Clone, Debug)]
pub struct FilteredSystem<W: ModWord = u64> {
    modulus: W::Modulus,
    columns: usize,
    /// Rows of the core, over the original column numbering.
    pub core_rows: Vec<SparseRow<W>>,
    /// Original ids of the columns that remain in the core, ascending.
    pub core_columns: Vec<u32>,
    /// `(pivot row, pivot column)` in elimination order.
    eliminated: Vec<(SparseRow<W>, u32)>,
    pub report: FilterReport,
}

/// The filter found a contradiction, which only a wrong relation can
/// produce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Inconsistent;

struct Workspace<W: ModWord> {
    modulus: W::Modulus,
    rows: Vec<Option<SparseRow<W>>>,
    /// Active rows containing each column (unordered).
    col_rows: Vec<Vec<usize>>,
    col_eliminated: Vec<bool>,
    active_rows: usize,
    eliminated: Vec<(SparseRow<W>, u32)>,
    report: FilterReport,
}

impl<W: ModWord> Workspace<W> {
    fn column_weight(&self, c: u32) -> usize {
        self.col_rows[c as usize].len()
    }

    fn active_columns(&self) -> usize {
        (0..self.col_rows.len())
            .filter(|&c| !self.col_eliminated[c] && !self.col_rows[c].is_empty())
            .count()
    }

    fn detach(&mut self, row_id: usize) -> Option<SparseRow<W>> {
        let row = self.rows[row_id].take()?;
        for &(c, _) in &row.entries {
            let list = &mut self.col_rows[c as usize];
            if let Some(pos) = list.iter().position(|&r| r == row_id) {
                list.swap_remove(pos);
            }
        }
        self.active_rows -= 1;
        Some(row)
    }

    fn attach(&mut self, row_id: usize, row: SparseRow<W>) {
        for &(c, _) in &row.entries {
            self.col_rows[c as usize].push(row_id);
        }
        self.rows[row_id] = Some(row);
        self.active_rows += 1;
    }

    /// Eliminate column `c` (weight ≥ 1) with its lightest row as pivot.
    /// Returns the columns whose weight changed, or `Err` on a
    /// contradiction.
    fn eliminate(
        &mut self,
        c: u32,
        max_row_weight: usize,
    ) -> Result<Option<Vec<u32>>, Inconsistent> {
        let m = self.modulus;
        let members = self.col_rows[c as usize].clone();
        if members.is_empty() {
            return Ok(None);
        }
        let pivot_id = members
            .iter()
            .copied()
            .min_by_key(|&r| {
                (
                    self.rows[r].as_ref().map_or(usize::MAX, SparseRow::weight),
                    r,
                )
            })
            .expect("nonempty");
        let pivot_weight = self.rows[pivot_id].as_ref().expect("active").weight();
        // Fill-in bound: every other row grows by at most pivot_weight − 2.
        for &r in &members {
            if r == pivot_id {
                continue;
            }
            let w = self.rows[r].as_ref().expect("active").weight();
            if w + pivot_weight - 2 > max_row_weight {
                return Ok(None);
            }
        }
        let pivot = self.detach(pivot_id).expect("active");
        let pc = pivot.coefficient(c).expect("member");
        let inv = pc.inv(&m).ok_or(Inconsistent)?;
        let mut touched: Vec<u32> = pivot.entries.iter().map(|&(c, _)| c).collect();
        for &r in &members {
            if r == pivot_id {
                continue;
            }
            let row = self.detach(r).expect("active");
            let factor = row.coefficient(c).expect("member").mul(inv, &m);
            let reduced = row.axpy(factor, &pivot, m);
            debug_assert!(reduced.coefficient(c).is_none());
            if reduced.entries.is_empty() {
                if !reduced.rhs.is_zero() {
                    return Err(Inconsistent);
                }
                self.report.dependent_rows_dropped += 1;
                continue;
            }
            touched.extend(reduced.entries.iter().map(|&(c, _)| c));
            self.attach(r, reduced);
        }
        self.col_eliminated[c as usize] = true;
        self.eliminated.push((pivot, c));
        touched.sort_unstable();
        touched.dedup();
        touched.retain(|&t| t != c);
        Ok(Some(touched))
    }

    /// Eliminate every singleton column, cascading.
    fn cascade_singletons(&mut self, seeds: Vec<u32>) -> Result<(), Inconsistent> {
        let mut stack = seeds;
        while let Some(c) = stack.pop() {
            if self.col_eliminated[c as usize] || self.column_weight(c) != 1 {
                continue;
            }
            if let Some(touched) = self.eliminate(c, usize::MAX)? {
                self.report.singletons_removed += 1;
                stack.extend(touched.into_iter().filter(|&t| self.column_weight(t) == 1));
            }
        }
        Ok(())
    }

    fn remove_row(&mut self, row_id: usize) -> Result<(), Inconsistent> {
        let row = self.detach(row_id).expect("active");
        self.report.excess_rows_removed += 1;
        let seeds = row.entries.iter().map(|&(c, _)| c).collect();
        self.cascade_singletons(seeds)
    }
}

fn find_root(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Reduce a relation system before the sparse solve.  `columns` is the
/// number of unknowns; every row's columns must be below it.
pub fn filter_relations<W: ModWord>(
    rows: Vec<SparseRow<W>>,
    columns: usize,
    modulus: W::Modulus,
    opts: &FilterOptions,
) -> Result<FilteredSystem<W>, Inconsistent> {
    let mut report = FilterReport {
        rows_in: rows.len(),
        columns_in: columns,
        nonzeros_in: rows.iter().map(SparseRow::weight).sum(),
        ..FilterReport::default()
    };
    // Duplicates and trivial rows.
    let mut kept: Vec<SparseRow<W>> = Vec::with_capacity(rows.len());
    let mut seen: HashSet<SparseRow<W>> = HashSet::new();
    for row in rows {
        if row.entries.is_empty() {
            if !row.rhs.is_zero() {
                return Err(Inconsistent);
            }
            report.dependent_rows_dropped += 1;
            continue;
        }
        if opts.remove_duplicates {
            if seen.contains(&row) {
                report.duplicates_removed += 1;
                continue;
            }
            seen.insert(row.clone());
        }
        kept.push(row);
    }
    drop(seen);

    let mut ws = Workspace {
        modulus,
        rows: Vec::with_capacity(kept.len()),
        col_rows: vec![Vec::new(); columns],
        col_eliminated: vec![false; columns],
        active_rows: 0,
        eliminated: Vec::new(),
        report,
    };
    for (id, row) in kept.into_iter().enumerate() {
        ws.rows.push(None);
        ws.attach(id, row);
    }

    // Singletons.
    if opts.remove_singletons {
        let seeds: Vec<u32> = (0..columns as u32)
            .filter(|&c| ws.column_weight(c) == 1)
            .collect();
        ws.cascade_singletons(seeds)?;
    }

    // Excess: remove rows, preferring those whose removal cascades
    // through the largest connected component of weight-2 columns
    // (the clique rule), heaviest row first.
    if opts.target_excess != usize::MAX && opts.remove_singletons {
        loop {
            let active_cols = ws.active_columns();
            if ws.active_rows <= active_cols.saturating_add(opts.target_excess) {
                break;
            }
            let to_remove = ws.active_rows - active_cols - opts.target_excess;
            let n = ws.rows.len();
            let mut parent: Vec<usize> = (0..n).collect();
            for c in 0..columns {
                if ws.col_eliminated[c] || ws.col_rows[c].len() != 2 {
                    continue;
                }
                let (a, b) = (ws.col_rows[c][0], ws.col_rows[c][1]);
                let (ra, rb) = (find_root(&mut parent, a), find_root(&mut parent, b));
                if ra != rb {
                    parent[ra] = rb;
                }
            }
            let mut size: HashMap<usize, usize> = HashMap::new();
            let mut heaviest: HashMap<usize, (usize, usize)> = HashMap::new();
            for id in 0..n {
                let Some(row) = ws.rows[id].as_ref() else {
                    continue;
                };
                let root = find_root(&mut parent, id);
                *size.entry(root).or_insert(0) += 1;
                let entry = heaviest.entry(root).or_insert((0, id));
                if row.weight() > entry.0 {
                    *entry = (row.weight(), id);
                }
            }
            let mut components: Vec<(usize, usize, usize)> = size
                .iter()
                .map(|(&root, &s)| (s, heaviest[&root].0, heaviest[&root].1))
                .collect();
            components.sort_unstable_by(|a, b| b.cmp(a));
            let mut removed = 0;
            for &(_, _, row_id) in components.iter().take(to_remove) {
                if ws.rows[row_id].is_some() {
                    ws.remove_row(row_id)?;
                    removed += 1;
                }
            }
            if removed == 0 {
                break;
            }
        }
    }

    // Merge light columns, lightest first, under the fill-in bound.
    if opts.merge_max_weight >= 2 {
        let mut heap: BinaryHeap<Reverse<(usize, u32)>> = BinaryHeap::new();
        for c in 0..columns {
            let w = ws.col_rows[c].len();
            if !ws.col_eliminated[c] && (1..=opts.merge_max_weight).contains(&w) {
                heap.push(Reverse((w, c as u32)));
            }
        }
        while let Some(Reverse((w, c))) = heap.pop() {
            if ws.col_eliminated[c as usize] || ws.column_weight(c) != w {
                continue; // stale entry
            }
            let Some(touched) = ws.eliminate(c, opts.max_row_weight)? else {
                continue; // fill-in bound; retried if its weight changes
            };
            if w == 1 {
                ws.report.singletons_removed += 1;
            } else {
                ws.report.merged_columns += 1;
            }
            for t in touched {
                let tw = ws.column_weight(t);
                if !ws.col_eliminated[t as usize] && (1..=opts.merge_max_weight).contains(&tw) {
                    heap.push(Reverse((tw, t)));
                }
            }
        }
    }

    let core_rows: Vec<SparseRow<W>> = ws.rows.iter().flatten().cloned().collect();
    let core_columns: Vec<u32> = (0..columns as u32)
        .filter(|&c| !ws.col_eliminated[c as usize] && !ws.col_rows[c as usize].is_empty())
        .collect();
    let mut report = ws.report;
    report.uncovered_columns = (0..columns)
        .filter(|&c| !ws.col_eliminated[c] && ws.col_rows[c].is_empty())
        .count();
    report.rows_out = core_rows.len();
    report.columns_out = core_columns.len();
    report.nonzeros_out = core_rows.iter().map(SparseRow::weight).sum();
    Ok(FilteredSystem {
        modulus,
        columns,
        core_rows,
        core_columns,
        eliminated: ws.eliminated,
        report,
    })
}

impl<W: ModWord> FilteredSystem<W> {
    /// Number of unknowns of the original system.
    pub fn columns(&self) -> usize {
        self.columns
    }

    /// Recover the full solution from the core one (indexed like
    /// [`Self::core_columns`]) by back-substituting the elimination
    /// stack, then propagating through `original_rows` while some row
    /// has exactly one unknown, then solving whatever remains densely.
    /// `None` when the original system does not determine every column
    /// (the caller collects more relations).
    pub fn reconstruct(
        &self,
        core_solution: &[W],
        original_rows: &[SparseRow<W>],
    ) -> Option<Vec<W>> {
        let m = &self.modulus;
        assert_eq!(
            core_solution.len(),
            self.core_columns.len(),
            "core solution width"
        );
        let mut value: Vec<Option<W>> = vec![None; self.columns];
        for (&c, &v) in self.core_columns.iter().zip(core_solution) {
            value[c as usize] = Some(v);
        }
        let solve_row = |row: &SparseRow<W>, col: u32, value: &[Option<W>]| -> Option<W> {
            let mut acc = row.rhs;
            let mut pivot = None;
            for &(c, v) in &row.entries {
                if c == col {
                    pivot = Some(v);
                } else {
                    acc = acc.sub(v.mul(value[c as usize]?, m), m);
                }
            }
            Some(acc.mul(pivot?.inv(m)?, m))
        };
        // Stack in reverse: a pivot row's other columns are core columns
        // or were eliminated later.
        let mut deferred = Vec::new();
        for (row, col) in self.eliminated.iter().rev() {
            match solve_row(row, *col, &value) {
                Some(v) => value[*col as usize] = Some(v),
                None => deferred.push((row, *col)),
            }
        }
        // Propagation over the original rows (and deferred pivots) while
        // some row has exactly one unknown.
        let mut progress = true;
        while progress && value.iter().any(Option::is_none) {
            progress = false;
            for (row, col) in &deferred {
                if value[*col as usize].is_none() {
                    if let Some(v) = solve_row(row, *col, &value) {
                        value[*col as usize] = Some(v);
                        progress = true;
                    }
                }
            }
            for row in original_rows {
                let mut unknown = None;
                let mut count = 0;
                for &(c, _) in &row.entries {
                    if value[c as usize].is_none() {
                        unknown = Some(c);
                        count += 1;
                        if count > 1 {
                            break;
                        }
                    }
                }
                if count == 1 {
                    let c = unknown.expect("one unknown");
                    if let Some(v) = solve_row(row, c, &value) {
                        value[c as usize] = Some(v);
                        progress = true;
                    }
                }
            }
        }
        // Residual dense solve over the columns still unknown.
        let unknown: Vec<u32> = (0..self.columns as u32)
            .filter(|&c| value[c as usize].is_none())
            .collect();
        if !unknown.is_empty() {
            let index: HashMap<u32, usize> =
                unknown.iter().enumerate().map(|(i, &c)| (c, i)).collect();
            let mut rows = Vec::new();
            let mut rhs = Vec::new();
            for row in original_rows {
                if !row.entries.iter().any(|&(c, _)| index.contains_key(&c)) {
                    continue;
                }
                let mut dense = vec![W::zero(); unknown.len()];
                let mut b = row.rhs;
                for &(c, v) in &row.entries {
                    match index.get(&c) {
                        Some(&i) => dense[i] = v,
                        None => b = b.sub(v.mul(value[c as usize].expect("known"), m), m),
                    }
                }
                rows.push(dense);
                rhs.push(b);
            }
            let solved = dense_solve_unique(rows, rhs, self.modulus)?;
            for (&c, v) in unknown.iter().zip(solved) {
                value[c as usize] = Some(v);
            }
        }
        value.into_iter().collect()
    }
}

/// Gauss–Jordan over `Z/rZ` in the word type; `Some` only for a unique
/// solution.
pub fn dense_solve_unique<W: ModWord>(
    mut rows: Vec<Vec<W>>,
    mut rhs: Vec<W>,
    modulus: W::Modulus,
) -> Option<Vec<W>> {
    let m = &modulus;
    let cols = rows.first().map_or(0, Vec::len);
    if cols == 0 {
        return rhs.iter().all(|b| b.is_zero()).then_some(Vec::new());
    }
    let mut pivot_row = 0;
    let mut pivots = Vec::with_capacity(cols);
    for col in 0..cols {
        let Some(p) = (pivot_row..rows.len()).find(|&r| !rows[r][col].is_zero()) else {
            return None; // free column: not unique
        };
        rows.swap(pivot_row, p);
        rhs.swap(pivot_row, p);
        let inv = rows[pivot_row][col].inv(m)?;
        for v in rows[pivot_row].iter_mut() {
            *v = v.mul(inv, m);
        }
        rhs[pivot_row] = rhs[pivot_row].mul(inv, m);
        let pivot = rows[pivot_row].clone();
        let pb = rhs[pivot_row];
        for r in 0..rows.len() {
            if r == pivot_row || rows[r][col].is_zero() {
                continue;
            }
            let f = rows[r][col];
            for k in col..cols {
                if !pivot[k].is_zero() {
                    rows[r][k] = rows[r][k].sub(f.mul(pivot[k], m), m);
                }
            }
            rhs[r] = rhs[r].sub(f.mul(pb, m), m);
        }
        pivots.push(pivot_row);
        pivot_row += 1;
    }
    // Remaining rows must read 0 ≡ 0.
    if (pivot_row..rows.len()).any(|r| !rhs[r].is_zero()) {
        return None;
    }
    Some(pivots.iter().map(|&r| rhs[r]).collect())
}

// ── Sparse matrices and block operators ────────────────────────────

/// Compressed-sparse-row matrix over `Z/rZ`.
#[derive(Clone, Debug)]
pub struct CsrMatrix<W: ModWord = u64> {
    pub n_rows: usize,
    pub n_cols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
    /// Coefficients, in lazy form ([`ModWord::lazy_form`]) when `lazy`
    /// is set and canonical otherwise.
    vals: Vec<W>,
    modulus: W::Modulus,
    /// The context of the lazily reduced row sums ([`ModWord::lazy`]
    /// for the largest stored coefficient); `None` when the width has
    /// no lazy path for the modulus (`u64` outside `[2, 2⁶³)`).
    lazy: Option<W::Lazy>,
    spmv: SpmvOptions,
}

/// Sparse matrix-times-block execution backend.  `Worker` uses the stable
/// `SPMV1` contract implemented by `gpu/spmv`; it may be a CUDA executable or
/// a distributed launcher.  Worker output is always cross-checked on the CPU.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SpmvBackend {
    #[default]
    Auto,
    Serial,
    Rayon,
    Sharded,
    Worker,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SpmvOptions {
    pub backend: SpmvBackend,
    /// Deterministic row shards used by `sharded`; zero selects the Rayon
    /// thread count.  A distributed launcher can use the same row ranges.
    pub shards: usize,
}

impl Default for SpmvOptions {
    fn default() -> Self {
        Self {
            backend: SpmvBackend::Auto,
            shards: 0,
        }
    }
}

/// Rows below this many are multiplied serially.
const PARALLEL_ROWS: usize = 2048;

impl<W: ModWord> CsrMatrix<W> {
    /// From rows over columns `0..n_cols` (right-hand sides ignored).
    pub fn from_rows(rows: &[SparseRow<W>], n_cols: usize, modulus: W::Modulus) -> Self {
        let m = &modulus;
        let mut row_ptr = Vec::with_capacity(rows.len() + 1);
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        row_ptr.push(0);
        for row in rows {
            for &(c, v) in &row.entries {
                debug_assert!((c as usize) < n_cols);
                col_idx.push(c);
                vals.push(v);
            }
            row_ptr.push(col_idx.len());
        }
        // Rows are built reduced, but `entries` is public: bound the
        // products by the largest coefficient actually stored, and by
        // `m − 1` for the homogenising column `Homogenised` appends.
        // (Only the `u64` width needs the bound; a wide coefficient is
        // canonicalised by its lazy form.)
        let a_max = vals
            .iter()
            .fold(W::max_residue(m), |acc: W, &v| W::larger(acc, v));
        let lazy = W::lazy(m, a_max);
        if lazy.is_some() {
            vals.iter_mut().for_each(|v| *v = W::lazy_form(*v, m));
        }
        Self {
            n_rows: rows.len(),
            n_cols,
            row_ptr,
            col_idx,
            vals,
            modulus,
            lazy,
            spmv: SpmvOptions::default(),
        }
    }

    pub fn with_spmv(mut self, options: SpmvOptions) -> Self {
        self.spmv = options;
        self
    }

    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.vals.len()
    }

    /// The modulus.
    pub fn modulus(&self) -> W::Modulus {
        self.modulus
    }

    /// Stored coefficient `k` as a canonical residue.
    fn coefficient(&self, k: usize) -> W {
        if self.lazy.is_some() {
            W::from_lazy_form(self.vals[k], &self.modulus)
        } else {
            self.vals[k]
        }
    }

    /// Stored coefficient `k` times a residue.
    #[inline]
    fn coefficient_times(&self, k: usize, x: W) -> W {
        if self.lazy.is_some() {
            W::lazy_mul(self.vals[k], x, &self.modulus)
        } else {
            self.vals[k].mul(x, &self.modulus)
        }
    }

    #[inline]
    fn row_into(&self, i: usize, x: &[W], n: usize, out: &mut [W]) {
        match &self.lazy {
            Some(lazy) => self.row_lazy(i, x, n, W::zero(), lazy, out),
            None => self.row_into_per_entry(i, x, n, out),
        }
    }

    /// Row `i` of `A · x` over `n` lanes, plus `extra · x[n_cols·n + j]` in
    /// lane `j` when `extra ≠ 0` (the homogenising column, in lazy
    /// form), with one reduction per lane ([`LazyLanes`]).  Lanes go
    /// four at a time (the default block width), then singly, so the
    /// accumulators stay in registers.
    #[inline]
    fn row_lazy(&self, i: usize, x: &[W], n: usize, extra: W, lazy: &W::Lazy, out: &mut [W]) {
        let span = self.row_ptr[i]..self.row_ptr[i + 1];
        let row = (&self.col_idx[span.clone()], &self.vals[span]);
        let extra = (self.n_cols * n, extra);
        let mut j = 0;
        while j + 4 <= n {
            let acc = LazyLanes::<W, 4>::new(lazy);
            csr_lanes(acc, row, &x[j..], n, extra, &mut out[j..j + 4]);
            j += 4;
        }
        while j < n {
            let acc = LazyLanes::<W, 1>::new(lazy);
            csr_lanes(acc, row, &x[j..], n, extra, &mut out[j..j + 1]);
            j += 1;
        }
    }

    /// One reduced product per entry and lane: the path for a modulus
    /// the width has no lazy context for, and the reference the lazy
    /// path is tested against.
    fn row_into_per_entry(&self, i: usize, x: &[W], n: usize, out: &mut [W]) {
        let m = &self.modulus;
        out.iter_mut().for_each(|v| *v = W::zero());
        for k in self.row_ptr[i]..self.row_ptr[i + 1] {
            let c = self.col_idx[k] as usize;
            let src = &x[c * n..(c + 1) * n];
            for (o, &s) in out.iter_mut().zip(src) {
                *o = o.add(self.coefficient_times(k, s), m);
            }
        }
    }

    /// `row(i, out)` for every `n`-lane row of `y`, in parallel over rows
    /// for a large matrix.
    fn for_each_row(&self, y: &mut [W], n: usize, row: impl Fn(usize, &mut [W]) + Sync) {
        if self.n_rows >= PARALLEL_ROWS {
            y.par_chunks_mut(n)
                .enumerate()
                .for_each(|(i, out)| row(i, out));
        } else {
            for (i, out) in y.chunks_mut(n).enumerate() {
                row(i, out);
            }
        }
    }

    fn mul_block_serial(&self, x: &[W], n: usize, y: &mut [W]) {
        for (i, out) in y.chunks_mut(n).enumerate() {
            self.row_into(i, x, n, out);
        }
    }

    fn mul_block_rayon(&self, x: &[W], n: usize, y: &mut [W]) {
        y.par_chunks_mut(n)
            .enumerate()
            .for_each(|(i, out)| self.row_into(i, x, n, out));
    }

    fn mul_block_sharded(&self, x: &[W], n: usize, y: &mut [W]) {
        let shards = if self.spmv.shards == 0 {
            rayon::current_num_threads()
        } else {
            self.spmv.shards
        }
        .max(1)
        .min(self.n_rows.max(1));
        let rows_per_shard = self.n_rows.div_ceil(shards);
        y.par_chunks_mut(rows_per_shard * n)
            .enumerate()
            .for_each(|(shard, output)| {
                let first = shard * rows_per_shard;
                for (offset, row) in output.chunks_mut(n).enumerate() {
                    self.row_into(first + offset, x, n, row);
                }
            });
    }

    fn worker_contract(&self, x: &[W], lanes: usize) -> String {
        let mut contract = format!(
            "SPMV1 {} {} {} {} {}\n",
            self.n_rows,
            self.n_cols,
            lanes,
            W::modulus_to_biguint(&self.modulus),
            self.vals.len()
        );
        let append = |target: &mut String, values: &mut dyn Iterator<Item = String>| {
            for (index, value) in values.enumerate() {
                if index != 0 {
                    target.push(' ');
                }
                target.push_str(&value);
            }
            target.push('\n');
        };
        append(
            &mut contract,
            &mut self.row_ptr.iter().map(|v| v.to_string()),
        );
        append(
            &mut contract,
            &mut self.col_idx.iter().map(|v| v.to_string()),
        );
        append(
            &mut contract,
            &mut (0..self.vals.len()).map(|k| self.coefficient(k).to_biguint().to_string()),
        );
        append(
            &mut contract,
            &mut x.iter().map(|v| v.to_biguint().to_string()),
        );
        contract
    }

    fn mul_block_worker(&self, x: &[W], lanes: usize, y: &mut [W]) -> Result<(), String> {
        let worker = std::env::var_os("IC_SPMV_WORKER").ok_or("IC_SPMV_WORKER is not set")?;
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_nanos())
            .unwrap_or(0);
        let mut input = std::env::temp_dir();
        input.push(format!("ic-spmv-{}-{nonce}.txt", std::process::id()));
        let mut input_file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&input)
            .map_err(|error| format!("create SPMV1 input: {error}"))?;
        let write_result = input_file.write_all(self.worker_contract(x, lanes).as_bytes());
        drop(input_file);
        if let Err(error) = write_result {
            let _ = std::fs::remove_file(&input);
            return Err(format!("write SPMV1 input: {error}"));
        }
        let result = Command::new(worker).arg("--in").arg(&input).output();
        let _ = std::fs::remove_file(&input);
        let output = result.map_err(|error| format!("run SPMV1 worker: {error}"))?;
        if !output.status.success() {
            return Err(format!("SPMV1 worker exited with {}", output.status));
        }
        let text = std::str::from_utf8(&output.stdout)
            .map_err(|error| format!("SPMV1 output is not UTF-8: {error}"))?;
        let mut tokens = text.split_whitespace();
        if tokens.next() != Some("SPMV1") {
            return Err("SPMV1 worker returned the wrong header".into());
        }
        let rows: usize = tokens
            .next()
            .ok_or("SPMV1 output is missing rows")?
            .parse()
            .map_err(|_| "SPMV1 rows are invalid")?;
        let returned_lanes: usize = tokens
            .next()
            .ok_or("SPMV1 output is missing lanes")?
            .parse()
            .map_err(|_| "SPMV1 lanes are invalid")?;
        if rows != self.n_rows || returned_lanes != lanes {
            return Err("SPMV1 output dimensions do not match the request".into());
        }
        for value in y.iter_mut() {
            let parsed: BigUint = tokens
                .next()
                .ok_or("SPMV1 output is truncated")?
                .parse()
                .map_err(|_| "SPMV1 output value is invalid")?;
            *value = W::from_biguint(&parsed, &self.modulus);
        }
        if tokens.next().is_some() {
            return Err("SPMV1 output has trailing values".into());
        }
        // An accelerator proposes a product; the portable implementation is
        // the authority until the worker has a separately audited certificate.
        let mut reference = vec![W::zero(); y.len()];
        self.mul_block_serial(x, lanes, &mut reference);
        if reference != y {
            return Err("SPMV1 worker disagrees with the CPU reference".into());
        }
        Ok(())
    }

    /// `y = A · x` for `n` column vectors stored row-major
    /// (`x[c * n + j]`), in parallel over rows.  The wide widths need
    /// `x` reduced; the `u64` width accepts any values.
    pub fn mul_block(&self, x: &[W], n: usize, y: &mut [W]) {
        debug_assert_eq!(x.len(), self.n_cols * n);
        debug_assert_eq!(y.len(), self.n_rows * n);
        match self.spmv.backend {
            SpmvBackend::Auto => self.for_each_row(y, n, |i, out| self.row_into(i, x, n, out)),
            SpmvBackend::Serial => self.mul_block_serial(x, n, y),
            SpmvBackend::Rayon => self.mul_block_rayon(x, n, y),
            SpmvBackend::Sharded => self.mul_block_sharded(x, n, y),
            SpmvBackend::Worker => {
                if let Err(error) = self.mul_block_worker(x, n, y) {
                    static WARNING: std::sync::Once = std::sync::Once::new();
                    WARNING.call_once(|| {
                        eprintln!(
                            "warning: external SPMV1 worker failed; using serial CPU: {error}"
                        );
                    });
                    self.mul_block_serial(x, n, y);
                }
            }
        }
    }

    /// `y = A^T · x` for `n` lanes.  This is kept beside the forward
    /// product so block Lanczos can form the verified normal operator
    /// `A^T D A` without materialising it.
    pub fn mul_block_transpose(&self, x: &[W], n: usize, y: &mut [W]) {
        debug_assert_eq!(x.len(), self.n_rows * n);
        debug_assert_eq!(y.len(), self.n_cols * n);
        y.fill(W::zero());
        let modulus = &self.modulus;
        for row in 0..self.n_rows {
            let source = &x[row * n..(row + 1) * n];
            for index in self.row_ptr[row]..self.row_ptr[row + 1] {
                let column = self.col_idx[index] as usize;
                let target = &mut y[column * n..(column + 1) * n];
                for (out, &value) in target.iter_mut().zip(source) {
                    *out = out.add(self.coefficient_times(index, value), modulus);
                }
            }
        }
    }
}

/// `N` lanes of one CSR row `(columns, coefficients)` times `x`, where
/// `x` starts at the first of those lanes and holds `n` lanes per column,
/// plus `extra.1` times the lanes at offset `extra.0` when nonzero.
#[inline(always)]
fn csr_lanes<W: ModWord, const N: usize>(
    mut acc: LazyLanes<'_, W, N>,
    (cols, vals): (&[u32], &[W]),
    x: &[W],
    n: usize,
    extra: (usize, W),
    out: &mut [W],
) {
    acc.add_sparse(cols, vals, x, n);
    if !extra.1.is_zero() {
        acc.add(extra.1, &x[extra.0..]);
    }
    acc.finish(out);
}

/// A square linear operator on block vectors, as block Wiedemann sees it.
pub trait BlockOperator<W: ModWord>: Sync {
    /// Dimension `N` of the square operator.
    fn dim(&self) -> usize;
    /// `y = M · x` on `n` column vectors stored row-major.
    fn apply(&self, x: &[W], n: usize, y: &mut [W]);
    /// The modulus.
    fn modulus(&self) -> W::Modulus;
}

impl<W: ModWord> BlockOperator<W> for CsrMatrix<W> {
    fn dim(&self) -> usize {
        assert_eq!(self.n_rows, self.n_cols, "square operator");
        self.n_rows
    }
    fn apply(&self, x: &[W], n: usize, y: &mut [W]) {
        self.mul_block(x, n, y)
    }
    fn modulus(&self) -> W::Modulus {
        self.modulus
    }
}

/// `M = [[A, −b], [0, 0]]` of dimension `N + 1`: a kernel vector `(x, t)`
/// with `t ≠ 0` gives `A (x / t) = b`.  Never materialised.
pub struct Homogenised<'a, W: ModWord = u64> {
    pub a: &'a CsrMatrix<W>,
    pub b: &'a [W],
    /// `−b`, reduced, in the lazy form the fused row product multiplies
    /// by (canonical when `a` has no lazy path).
    neg_b: Vec<W>,
}

impl<'a, W: ModWord> Homogenised<'a, W> {
    /// `b` may hold unreduced values; they are reduced here.
    pub fn new(a: &'a CsrMatrix<W>, b: &'a [W]) -> Self {
        let m = &a.modulus;
        let neg_b = b
            .iter()
            .map(|&bi| {
                let v = bi.reduce(m).neg(m);
                if a.lazy.is_some() {
                    W::lazy_form(v, m)
                } else {
                    v
                }
            })
            .collect();
        Self { a, b, neg_b }
    }
}

impl<W: ModWord> BlockOperator<W> for Homogenised<'_, W> {
    fn dim(&self) -> usize {
        self.a.dim() + 1
    }
    fn apply(&self, x: &[W], n: usize, y: &mut [W]) {
        let big_n = self.a.dim();
        let m = &self.a.modulus;
        let (ya, yt) = y.split_at_mut(big_n * n);
        yt.iter_mut().for_each(|v| *v = W::zero());
        if let Some(lazy) = &self.a.lazy {
            // `−b_i · t_j` is one more product in row `i`'s sum, against
            // the `t` lanes that follow `A`'s `N·n` in `x`: one reduction
            // per lane for the whole row of `M`.
            self.a.for_each_row(ya, n, |i, out| {
                self.a.row_lazy(i, x, n, self.neg_b[i], lazy, out);
            });
            return;
        }
        self.a.mul_block(&x[..big_n * n], n, ya);
        let t = &x[big_n * n..];
        for (i, out) in ya.chunks_mut(n).enumerate() {
            let nbi = self.neg_b[i];
            if nbi.is_zero() {
                continue;
            }
            for (o, &tj) in out.iter_mut().zip(t) {
                *o = o.add(nbi.mul(tj, m), m);
            }
        }
    }
    fn modulus(&self) -> W::Modulus {
        self.a.modulus
    }
}

// ── Block Wiedemann ────────────────────────────────────────────────

/// Block sizes and safety margin for [`block_wiedemann_kernel`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BlockWiedemannOptions {
    /// Rows of the left projection `X` (`m`).
    pub block_m: usize,
    /// Columns of the right block `Y` (`n`).
    pub block_n: usize,
    /// Sequence terms beyond `⌈N/m⌉ + ⌈N/n⌉`.
    pub margin: usize,
}

impl Default for BlockWiedemannOptions {
    fn default() -> Self {
        Self {
            block_m: 4,
            block_n: 4,
            margin: 8,
        }
    }
}

/// What one block Wiedemann run did.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct BlockWiedemannReport {
    pub dimension: usize,
    pub block_m: usize,
    pub block_n: usize,
    /// Krylov terms computed (`X Mⁱ Y`, `i < sequence_length`).
    pub sequence_length: usize,
    /// Nominal degree of the generator column that produced the kernel
    /// vector.
    pub generator_degree: usize,
    /// Extra applications of `M` needed after evaluating that column.
    pub extra_steps: usize,
    /// Matrix–block products performed in total.
    pub products: usize,
}

/// One column of the approximant basis: `f` (coefficients over the `n`
/// right vectors), the residual `[a | I] · [f; g] mod λ^L` (column-major
/// `m`-vectors per order), and the shifted degree bound `δ`.
struct BasisColumn<W: ModWord> {
    f: Vec<Vec<W>>,
    residual: Vec<W>,
    delta: usize,
}

/// Shifted minimal approximant basis of order `L` for `[a(λ) | I_m]`
/// with shift `(0ⁿ, 1ᵐ)`: every column `[f; g]` satisfies
/// `a·f + g ≡ 0 (mod λ^L)` with `deg f ≤ δ`, `deg g ≤ δ − 1`, hence the
/// block recurrence `Σ_{k≤δ} S_{s−k} f_k = 0` for `δ ≤ s < L`.
fn minimal_approximant_basis<W: ModWord>(
    seq: &[Vec<W>],
    m: usize,
    n: usize,
    modulus: W::Modulus,
) -> Vec<BasisColumn<W>> {
    let md = &modulus;
    let big_l = seq.len();
    let mut cols: Vec<BasisColumn<W>> = Vec::with_capacity(n + m);
    for j in 0..n {
        let mut e = vec![W::zero(); n];
        e[j] = W::one();
        let mut residual = vec![W::zero(); big_l * m];
        for (t, s) in seq.iter().enumerate() {
            for i in 0..m {
                residual[t * m + i] = s[i * n + j];
            }
        }
        cols.push(BasisColumn {
            f: vec![e],
            residual,
            delta: 0,
        });
    }
    for i in 0..m {
        let mut residual = vec![W::zero(); big_l * m];
        residual[i] = W::one();
        cols.push(BasisColumn {
            f: vec![vec![W::zero(); n]],
            residual,
            delta: 1,
        });
    }
    let mut order: Vec<usize> = (0..n + m).collect();
    let mut pivots: Vec<(usize, usize)> = Vec::with_capacity(m);
    for t in 0..big_l {
        order.sort_by_key(|&j| (cols[j].delta, j));
        pivots.clear();
        for &j in &order {
            for &(p, row) in &pivots {
                let c = cols[j].residual[t * m + row];
                if c.is_zero() {
                    continue;
                }
                // column j -= c · column p  (δ_p ≤ δ_j by the ordering)
                let (pj, pp) = if j < p {
                    let (lo, hi) = cols.split_at_mut(p);
                    (&mut lo[j], &hi[0])
                } else {
                    let (lo, hi) = cols.split_at_mut(j);
                    (&mut hi[0], &lo[p])
                };
                let width = pp.f.len();
                if pj.f.len() < width {
                    pj.f.resize(width, vec![W::zero(); n]);
                }
                // `c` scales a whole column: one fixed-multiplier set-up,
                // then a few multiplications per element instead of a
                // division (and no test for a zero entry, which the
                // product leaves be).
                let c = W::fixed(c, md);
                for (fj, fp) in pj.f.iter_mut().zip(&pp.f) {
                    for (a, &b) in fj.iter_mut().zip(fp) {
                        *a = W::fixed_sub_from(&c, *a, b, md);
                    }
                }
                let span = t * m..big_l * m;
                for (a, &b) in pj.residual[span.clone()].iter_mut().zip(&pp.residual[span]) {
                    *a = W::fixed_sub_from(&c, *a, b, md);
                }
            }
            if let Some(row) = (0..m).find(|&i| !cols[j].residual[t * m + i].is_zero()) {
                // Normalise the pivot so later reductions are one product.
                let inv = cols[j].residual[t * m + row]
                    .inv(md)
                    .expect("prime modulus");
                if inv != W::one() {
                    let inv = W::fixed(inv, md);
                    let col = &mut cols[j];
                    for fk in col.f.iter_mut() {
                        for a in fk.iter_mut() {
                            *a = W::fixed_mul(&inv, *a, md);
                        }
                    }
                    for a in &mut col.residual[t * m..big_l * m] {
                        *a = W::fixed_mul(&inv, *a, md);
                    }
                }
                pivots.push((j, row));
            }
        }
        for &(p, _) in &pivots {
            let col = &mut cols[p];
            col.f.insert(0, vec![W::zero(); n]);
            col.residual
                .copy_within(t * m..(big_l - 1) * m, (t + 1) * m);
            col.residual[t * m..(t + 1) * m]
                .iter_mut()
                .for_each(|v| *v = W::zero());
            col.delta += 1;
        }
    }
    cols
}

/// A nonzero kernel vector of a singular square operator by block
/// Wiedemann, or `None` if none of the generator columns produced one
/// (retry with another seed, or the operator is nonsingular).
pub fn block_wiedemann_kernel<W: ModWord>(
    op: &impl BlockOperator<W>,
    opts: &BlockWiedemannOptions,
    seed: u64,
) -> Option<(Vec<W>, BlockWiedemannReport)> {
    let big_n = op.dim();
    let modulus = op.modulus();
    let md = &modulus;
    let m = opts.block_m.max(1);
    let n = opts.block_n.max(1);
    if big_n == 0 {
        return None;
    }
    let big_l = big_n.div_ceil(m) + big_n.div_ceil(n) + opts.margin.max(2);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x424c_4f43_4b57_4945);
    // X, Z and the basis are reduced, so their products are bounded by
    // `m − 1` times whatever the operator returns.
    let lazy = W::lazy(md, W::max_residue(md));
    // X is only ever a coefficient in the Krylov projection, so it is
    // kept in lazy form (the identity for `u64`).
    let x: Vec<W> = (0..m * big_n)
        .map(|_| {
            let v = W::random(&mut rng, md);
            if lazy.is_some() {
                W::lazy_form(v, md)
            } else {
                v
            }
        })
        .collect();
    // Y = M Z for random Z (Coppersmith): a generator of the Y-sequence
    // evaluated on Z lands in the kernel of M instead of vanishing.
    let mut z: Vec<W> = (0..big_n * n).map(|_| W::random(&mut rng, md)).collect();
    let mut y = vec![W::zero(); big_n * n];
    op.apply(&z, n, &mut y);
    // From here on Z is a coefficient block too (`z_times`).
    if lazy.is_some() {
        z.iter_mut().for_each(|v| *v = W::lazy_form(*v, md));
    }
    let mut report = BlockWiedemannReport {
        dimension: big_n,
        block_m: m,
        block_n: n,
        sequence_length: big_l,
        products: 1,
        ..BlockWiedemannReport::default()
    };

    // Krylov sequence S_t = X Mᵗ Y = X Mᵗ⁺¹ Z.
    let mut seq: Vec<Vec<W>> = Vec::with_capacity(big_l);
    let mut v = y.clone();
    let mut w = vec![W::zero(); big_n * n];
    for t in 0..big_l {
        let mut s = vec![W::zero(); m * n];
        for i in 0..m {
            let xi = &x[i * big_n..(i + 1) * big_n];
            let out = &mut s[i * n..(i + 1) * n];
            if let Some(lazy) = &lazy {
                // Row i of X against the lanes of v, like a dense matrix
                // row: lanes four at a time, one reduction each.
                let mut j = 0;
                while j + 4 <= n {
                    let mut acc = LazyLanes::<W, 4>::new(lazy);
                    acc.add_strided(xi, &v[j..], n);
                    acc.finish(&mut out[j..j + 4]);
                    j += 4;
                }
                while j < n {
                    let mut acc = LazyLanes::<W, 1>::new(lazy);
                    acc.add_strided(xi, &v[j..], n);
                    acc.finish(&mut out[j..j + 1]);
                    j += 1;
                }
                continue;
            }
            for (j, o) in out.iter_mut().enumerate() {
                let mut acc = W::zero();
                for (r, &xr) in xi.iter().enumerate() {
                    acc = acc.add(xr.mul(v[r * n + j], md), md);
                }
                *o = acc;
            }
        }
        seq.push(s);
        if t + 1 < big_l {
            op.apply(&v, n, &mut w);
            report.products += 1;
            std::mem::swap(&mut v, &mut w);
        }
    }

    let basis = minimal_approximant_basis(&seq, m, n, modulus);
    let mut order: Vec<usize> = (0..basis.len()).collect();
    order.sort_by_key(|&j| basis[j].delta);
    drop(y);
    let mut buf = vec![W::zero(); big_n];
    let z_times = |f: &[W], out: &mut [W]| {
        for (r, o) in out.iter_mut().enumerate() {
            let zr = &z[r * n..(r + 1) * n];
            if let Some(lazy) = &lazy {
                // Z is the lazy-form side here, f the residues.
                let mut acc = LazyLanes::<W, 1>::new(lazy);
                for (k, &zk) in zr.iter().enumerate() {
                    acc.add(zk, &f[k..k + 1]);
                }
                acc.finish(std::slice::from_mut(o));
                continue;
            }
            let mut acc = W::zero();
            for (j, &fj) in f.iter().enumerate() {
                if !fj.is_zero() {
                    acc = acc.add(zr[j].mul(fj, md), md);
                }
            }
            *o = acc;
        }
    };
    for &j in &order {
        let col = &basis[j];
        let Some(e) = col.f.iter().rposition(|fk| fk.iter().any(|a| !a.is_zero())) else {
            continue;
        };
        // Horner from the top power: u = Σ_{k≤e} M^{e−k} Z f_k, so f_0
        // multiplies Mᵉ and f_e multiplies the identity.  The recurrence
        // gives M^{δ−e+1} u = 0, so some Mⁱu below that is in the kernel.
        let mut u = vec![W::zero(); big_n];
        z_times(&col.f[0], &mut u);
        for fk in &col.f[1..=e] {
            op.apply(&u, 1, &mut buf);
            report.products += 1;
            z_times(fk, &mut u);
            for (a, &b) in u.iter_mut().zip(&buf) {
                *a = a.add(b, md);
            }
        }
        if u.iter().all(|a| a.is_zero()) {
            continue;
        }
        // Push through M until it dies; the last nonzero image is in
        // the kernel.
        let limit = col.delta.saturating_sub(e) + 4;
        for step in 0..=limit {
            op.apply(&u, 1, &mut buf);
            report.products += 1;
            if buf.iter().all(|a| a.is_zero()) {
                report.generator_degree = col.delta;
                report.extra_steps = step;
                return Some((u, report));
            }
            std::mem::swap(&mut u, &mut buf);
        }
    }
    None
}

// ── Block Lanczos ─────────────────────────────────────────────────

/// Controls for the finite-field block-Lanczos core solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BlockLanczosOptions {
    /// Simultaneous Krylov directions.  The first right-hand side is the
    /// ECDLP system; the remaining lanes are deterministic random probes that
    /// reduce scalar Lanczos breakdowns.
    pub block_size: usize,
    /// Iterations beyond `ceil(N / block_size)` before declaring breakdown.
    pub margin: usize,
}

impl Default for BlockLanczosOptions {
    fn default() -> Self {
        Self {
            block_size: 4,
            margin: 8,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct BlockLanczosReport {
    pub dimension: usize,
    pub block_size: usize,
    pub iterations: usize,
    /// Applications of the original sparse matrix or its transpose.
    pub sparse_products: usize,
}

/// The symmetric, nonsingular (when `A` is nonsingular) congruence
/// `S = A^T D A`, where `D` is a random nonzero diagonal.  Working with this
/// operator is sound over the odd subgroup-order fields used by ECDLP: a
/// returned vector is accepted only after checking `A x = b` itself.
struct NormalOperator<'a, W: ModWord> {
    a: &'a CsrMatrix<W>,
    diagonal: Vec<W::Fixed>,
}

impl<W: ModWord> NormalOperator<'_, W> {
    fn apply(&self, x: &[W], lanes: usize, y: &mut [W], scratch: &mut [W]) {
        self.a.mul_block(x, lanes, scratch);
        let modulus = &self.a.modulus;
        for (row, block) in scratch.chunks_mut(lanes).enumerate() {
            let scale = &self.diagonal[row];
            for value in block {
                *value = W::fixed_mul(scale, *value, modulus);
            }
        }
        self.a.mul_block_transpose(scratch, lanes, y);
    }

    fn rhs(&self, b: &[W], lanes: usize, out: &mut [W]) {
        let modulus = &self.a.modulus;
        let mut scaled = b.to_vec();
        for (row, block) in scaled.chunks_mut(lanes).enumerate() {
            let scale = &self.diagonal[row];
            for value in block {
                *value = W::fixed_mul(scale, *value, modulus);
            }
        }
        self.a.mul_block_transpose(&scaled, lanes, out);
    }
}

fn small_inverse<W: ModWord>(matrix: &[W], n: usize, modulus: &W::Modulus) -> Option<Vec<W>> {
    let mut augmented = vec![W::zero(); n * 2 * n];
    for row in 0..n {
        for column in 0..n {
            augmented[row * 2 * n + column] = matrix[row * n + column];
        }
        augmented[row * 2 * n + n + row] = W::one();
    }
    for column in 0..n {
        let pivot = (column..n).find(|&row| !augmented[row * 2 * n + column].is_zero())?;
        if pivot != column {
            for j in 0..2 * n {
                augmented.swap(column * 2 * n + j, pivot * 2 * n + j);
            }
        }
        let inverse = augmented[column * 2 * n + column].inv(modulus)?;
        for j in 0..2 * n {
            let slot = &mut augmented[column * 2 * n + j];
            *slot = slot.mul(inverse, modulus);
        }
        for row in 0..n {
            if row == column {
                continue;
            }
            let factor = augmented[row * 2 * n + column];
            if factor.is_zero() {
                continue;
            }
            for j in 0..2 * n {
                let product = factor.mul(augmented[column * 2 * n + j], modulus);
                let slot = &mut augmented[row * 2 * n + j];
                *slot = slot.sub(product, modulus);
            }
        }
    }
    let mut inverse = vec![W::zero(); n * n];
    for row in 0..n {
        inverse[row * n..(row + 1) * n]
            .copy_from_slice(&augmented[row * 2 * n + n..row * 2 * n + 2 * n]);
    }
    Some(inverse)
}

/// `left^T * right` for two row-major `rows × lanes` blocks of residues.
/// The `rows` products of each entry are summed lazily when the width
/// has a lazy context for the modulus.
fn block_gram<W: ModWord>(
    left: &[W],
    right: &[W],
    rows: usize,
    lanes: usize,
    modulus: &W::Modulus,
) -> Vec<W> {
    let mut out = vec![W::zero(); lanes * lanes];
    if let Some(lazy) = W::lazy(modulus, W::max_residue(modulus)) {
        for i in 0..lanes {
            for j in 0..lanes {
                let mut acc = LazyLanes::<W, 1>::plain(&lazy);
                for row in 0..rows {
                    acc.add(left[row * lanes + i], &right[row * lanes + j..]);
                }
                acc.finish(&mut out[i * lanes + j..i * lanes + j + 1]);
            }
        }
        return out;
    }
    for row in 0..rows {
        for i in 0..lanes {
            let a = left[row * lanes + i];
            for j in 0..lanes {
                out[i * lanes + j] =
                    out[i * lanes + j].add(a.mul(right[row * lanes + j], modulus), modulus);
            }
        }
    }
    out
}

fn small_product<W: ModWord>(a: &[W], b: &[W], n: usize, modulus: &W::Modulus) -> Vec<W> {
    let mut out = vec![W::zero(); n * n];
    for i in 0..n {
        for k in 0..n {
            if a[i * n + k].is_zero() {
                continue;
            }
            for j in 0..n {
                out[i * n + j] =
                    out[i * n + j].add(a[i * n + k].mul(b[k * n + j], modulus), modulus);
            }
        }
    }
    out
}

/// `target += sign * block * small`, all row-major. `sign` is `1` or `-1`.
/// `small` is `lanes × lanes`; its entries are prepared as fixed
/// multipliers once, so the `rows · lanes²` products are each one
/// multiplication by a fixed residue.
fn block_axpy_small<W: ModWord>(
    target: &mut [W],
    block: &[W],
    small: &[W],
    rows: usize,
    lanes: usize,
    subtract: bool,
    modulus: &W::Modulus,
) {
    let fixed: Vec<W::Fixed> = small.iter().map(|&s| W::fixed(s, modulus)).collect();
    for row in 0..rows {
        for j in 0..lanes {
            let mut value = W::zero();
            for k in 0..lanes {
                value = value.add(
                    W::fixed_mul(&fixed[k * lanes + j], block[row * lanes + k], modulus),
                    modulus,
                );
            }
            let slot = &mut target[row * lanes + j];
            *slot = if subtract {
                slot.sub(value, modulus)
            } else {
                slot.add(value, modulus)
            };
        }
    }
}

/// Solve a nonsingular square sparse system with a block-Lanczos/conjugate-
/// direction iteration on `A^T D A`.  Breakdown is reported as `None`; callers
/// retry with a fresh diagonal.  The first lane is always the requested `b`,
/// and the returned solution is verified against `A x = b`.
pub fn block_lanczos_solve<W: ModWord>(
    a: &CsrMatrix<W>,
    b: &[W],
    options: &BlockLanczosOptions,
    seed: u64,
) -> Option<(Vec<W>, BlockLanczosReport)> {
    if a.n_rows == 0 || a.n_rows != a.n_cols || b.len() != a.n_rows {
        return None;
    }
    let dimension = a.n_rows;
    let modulus = &a.modulus;
    // A final short block needs Montgomery-style mask deflation.  Until that
    // specialised mask path is useful at the repository's matrix sizes, pick
    // the largest requested width dividing N; this preserves exact block
    // recurrences and falls back to scalar Lanczos only for prime dimensions.
    let requested_lanes = options.block_size.max(1).min(dimension);
    let lanes = (1..=requested_lanes)
        .rev()
        .find(|width| dimension.is_multiple_of(*width))
        .unwrap_or(1);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x424c_4f43_4b4c_414e);
    let normal = NormalOperator {
        a,
        diagonal: (0..dimension)
            .map(|_| W::fixed(W::random_nonzero(&mut rng, modulus), modulus))
            .collect(),
    };
    let mut right = vec![W::zero(); dimension * lanes];
    for row in 0..dimension {
        right[row * lanes] = b[row].reduce(modulus);
        for lane in 1..lanes {
            right[row * lanes + lane] = W::random(&mut rng, modulus);
        }
    }
    let mut residual = vec![W::zero(); dimension * lanes];
    normal.rhs(&right, lanes, &mut residual);
    let mut direction = residual.clone();
    let mut solution = vec![W::zero(); dimension * lanes];
    let mut image = vec![W::zero(); dimension * lanes];
    let mut scratch = vec![W::zero(); dimension * lanes];
    let mut report = BlockLanczosReport {
        dimension,
        block_size: lanes,
        sparse_products: 1,
        ..BlockLanczosReport::default()
    };
    let limit = dimension.div_ceil(lanes) + options.margin.max(2);
    for iteration in 0..limit {
        normal.apply(&direction, lanes, &mut image, &mut scratch);
        report.sparse_products += 2;
        let h = block_gram(&direction, &image, dimension, lanes, modulus);
        let h_inverse = small_inverse(&h, lanes, modulus)?;
        let projected = block_gram(&direction, &residual, dimension, lanes, modulus);
        let alpha = small_product(&h_inverse, &projected, lanes, modulus);
        block_axpy_small(
            &mut solution,
            &direction,
            &alpha,
            dimension,
            lanes,
            false,
            modulus,
        );
        block_axpy_small(
            &mut residual,
            &image,
            &alpha,
            dimension,
            lanes,
            true,
            modulus,
        );
        report.iterations = iteration + 1;
        // The first lane is the requested ECDLP right-hand side.  Random
        // companion lanes may lose rank in the final partial block when
        // `dimension` is not divisible by `lanes`; do not turn their expected
        // terminal deflation into a breakdown once lane zero has converged.
        if residual.chunks(lanes).all(|block| block[0].is_zero()) {
            let answer: Vec<W> = solution.chunks(lanes).map(|block| block[0]).collect();
            let mut check = vec![W::zero(); dimension];
            a.mul_block(&answer, 1, &mut check);
            if check.iter().zip(b).all(|(&x, &y)| x == y.reduce(modulus)) {
                return Some((answer, report));
            }
            return None;
        }
        // P_new = R_new - P H^{-1} (A P)^T R_new. This makes the new
        // block A-conjugate to the previous direction block.
        let cross = block_gram(&image, &residual, dimension, lanes, modulus);
        let beta = small_product(&h_inverse, &cross, lanes, modulus);
        let mut next_direction = residual.clone();
        block_axpy_small(
            &mut next_direction,
            &direction,
            &beta,
            dimension,
            lanes,
            true,
            modulus,
        );
        direction = next_direction;
    }
    None
}

// ── The sparse solve: filter → fold → block Wiedemann → reconstruct ──

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SparseCoreSolver {
    #[default]
    BlockWiedemann,
    BlockLanczos,
}

fn is_default<T: Default + PartialEq>(value: &T) -> bool {
    value == &T::default()
}

/// Controls for [`solve_sparse_system`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SparseSolveOptions {
    pub filter: FilterOptions,
    /// Omitted from JSON at its historical block-Wiedemann default so frozen
    /// job configurations retain their exact schema and digest.
    #[serde(skip_serializing_if = "is_default")]
    pub solver: SparseCoreSolver,
    #[serde(skip_serializing_if = "is_default")]
    pub spmv: SpmvOptions,
    pub wiedemann: BlockWiedemannOptions,
    #[serde(skip_serializing_if = "is_default")]
    pub lanczos: BlockLanczosOptions,
    /// Base rows each excess row is folded into.
    pub fold: usize,
    /// Independent `(fold, X, Y)` draws before giving up.
    pub attempts: usize,
    pub seed: u64,
}

impl Default for SparseSolveOptions {
    fn default() -> Self {
        Self {
            filter: FilterOptions::default(),
            solver: SparseCoreSolver::default(),
            spmv: SpmvOptions::default(),
            wiedemann: BlockWiedemannOptions::default(),
            lanczos: BlockLanczosOptions::default(),
            fold: 3,
            attempts: 3,
            seed: 0x5350_4152_5345_4c41,
        }
    }
}

/// What [`solve_sparse_system`] did.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct SparseSolveReport {
    pub filter: FilterReport,
    #[serde(skip_serializing_if = "is_default")]
    pub solver: SparseCoreSolver,
    #[serde(skip_serializing_if = "is_default")]
    pub spmv: SpmvOptions,
    /// Rows and columns of the square system handed to the Krylov solver.
    pub core_dimension: usize,
    pub core_nonzeros: usize,
    pub wiedemann: Option<BlockWiedemannReport>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lanczos: Option<BlockLanczosReport>,
    /// Draws of `(fold, X, Y)` used.
    pub attempts: usize,
    /// Columns recovered by back-substitution, propagation and the
    /// residual dense solve after the core.
    pub reconstructed_columns: usize,
}

/// Outcome of [`solve_sparse_system`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SparseSolveOutcome<T = u64> {
    /// A solution over every column, verified against every row.  It is
    /// the unique one whenever the rows have full column rank; the
    /// solver does not certify rank itself (a consistent rank-deficient
    /// system whose dependency survives filtering yields one of its
    /// solutions), which is why the logarithm precompute gates the
    /// solve on an independent rank tracker and certifies every
    /// logarithm in the group.
    Solved(Vec<T>),
    /// The rows do not (yet) determine every column, or the solve did
    /// not converge within the attempts.
    Undetermined,
    /// A contradiction: some relation is wrong.
    Inconsistent,
}

/// Fold rows `C..R` of a full-column-rank system into rows `< C` at
/// random with random coefficients, giving a square system with the same
/// solution.
fn fold_square<W: ModWord>(
    rows: &[SparseRow<W>],
    n_cols: usize,
    modulus: W::Modulus,
    fold: usize,
    rng: &mut StdRng,
) -> Vec<SparseRow<W>> {
    let m = &modulus;
    let mut square: Vec<SparseRow<W>> = rows[..n_cols].to_vec();
    for extra in &rows[n_cols..] {
        for _ in 0..fold.max(1) {
            let i = rng.gen_range(0..n_cols);
            let c = W::random_nonzero(rng, m);
            square[i] = square[i].axpy(c.neg(m), extra, modulus);
        }
    }
    square
}

/// Solve `Σ c_{i,o} x_o ≡ rhs_i (mod r)` for every `x_o` by filtering,
/// block Wiedemann on the reduced core, and exact reconstruction.
/// `columns` is the number of unknowns.  The solution is verified against
/// every input row before it is returned.
pub fn solve_sparse_system<W: ModWord>(
    rows: &[SparseRow<W>],
    columns: usize,
    modulus: W::Modulus,
    opts: &SparseSolveOptions,
) -> (SparseSolveOutcome<W>, SparseSolveReport) {
    let m = &modulus;
    let mut report = SparseSolveReport {
        solver: opts.solver,
        spmv: opts.spmv,
        ..SparseSolveReport::default()
    };
    // Cheap necessary condition: every column occurs somewhere.
    let mut covered = vec![false; columns];
    for row in rows {
        for &(c, _) in &row.entries {
            covered[c as usize] = true;
        }
    }
    if covered.iter().any(|&c| !c) || rows.len() < columns {
        report.filter.rows_in = rows.len();
        report.filter.columns_in = columns;
        report.filter.uncovered_columns = covered.iter().filter(|&&c| !c).count();
        return (SparseSolveOutcome::Undetermined, report);
    }
    let filtered = match filter_relations(rows.to_vec(), columns, modulus, &opts.filter) {
        Ok(f) => f,
        Err(Inconsistent) => return (SparseSolveOutcome::Inconsistent, report),
    };
    report.filter = filtered.report.clone();
    let core_cols = filtered.core_columns.len();
    if filtered.core_rows.len() < core_cols {
        return (SparseSolveOutcome::Undetermined, report);
    }
    report.core_dimension = core_cols;

    // Solve the core.
    let mut core_solution: Option<Vec<W>> = None;
    if core_cols == 0 {
        core_solution = Some(Vec::new());
    } else {
        // Renumber core columns densely.
        let index: HashMap<u32, u32> = filtered
            .core_columns
            .iter()
            .enumerate()
            .map(|(i, &c)| (c, i as u32))
            .collect();
        let core: Vec<SparseRow<W>> = filtered
            .core_rows
            .iter()
            .map(|row| SparseRow {
                entries: row.entries.iter().map(|&(c, v)| (index[&c], v)).collect(),
                rhs: row.rhs,
            })
            .collect();
        report.core_nonzeros = core.iter().map(SparseRow::weight).sum();
        let mut rng = StdRng::seed_from_u64(opts.seed ^ 0x464f_4c44_5351_5552);
        for attempt in 0..opts.attempts.max(1) {
            report.attempts = attempt + 1;
            let square = fold_square(&core, core_cols, modulus, opts.fold, &mut rng);
            let a = CsrMatrix::from_rows(&square, core_cols, modulus).with_spmv(opts.spmv);
            let b: Vec<W> = square.iter().map(|r| r.rhs).collect();
            let x = match opts.solver {
                SparseCoreSolver::BlockWiedemann => {
                    let op = Homogenised::new(&a, &b);
                    let Some((z, wr)) = block_wiedemann_kernel(&op, &opts.wiedemann, rng.gen())
                    else {
                        continue;
                    };
                    report.wiedemann = Some(wr);
                    let t = z[core_cols];
                    let Some(inv) = (!t.is_zero()).then(|| t.inv(m)).flatten() else {
                        continue;
                    };
                    z[..core_cols].iter().map(|&v| v.mul(inv, m)).collect()
                }
                SparseCoreSolver::BlockLanczos => {
                    let Some((x, lr)) = block_lanczos_solve(&a, &b, &opts.lanczos, rng.gen())
                    else {
                        continue;
                    };
                    report.lanczos = Some(lr);
                    x
                }
            };
            if core.iter().all(|row| row.evaluate(&x, modulus) == row.rhs) {
                core_solution = Some(x);
                break;
            }
        }
    }
    let Some(core_solution) = core_solution else {
        return (SparseSolveOutcome::Undetermined, report);
    };
    let Some(full) = filtered.reconstruct(&core_solution, rows) else {
        return (SparseSolveOutcome::Undetermined, report);
    };
    report.reconstructed_columns = columns - core_cols;
    if rows
        .iter()
        .all(|row| row.evaluate(&full, modulus) == row.rhs)
    {
        (SparseSolveOutcome::Solved(full), report)
    } else {
        (SparseSolveOutcome::Inconsistent, report)
    }
}

// ── Rows over whichever width a modulus needs ──────────────────────

/// Relation rows in the word width [`WordWidth::for_modulus`] picks for
/// a modulus given as a `BigUint`, with the operations the logarithm
/// precompute needs, so that caller never names a width.  Values cross
/// the boundary as `BigUint`s.
#[derive(Clone, Debug)]
pub enum SparseRows {
    U64 {
        modulus: u64,
        rows: Vec<SparseRow<u64>>,
    },
    U128 {
        modulus: WideModulus<2>,
        rows: Vec<SparseRow<U128>>,
    },
    U192 {
        modulus: WideModulus<3>,
        rows: Vec<SparseRow<U192>>,
    },
    U256 {
        modulus: WideModulus<4>,
        rows: Vec<SparseRow<U256>>,
    },
}

/// Runs `$body` with `$rows`/`$modulus` bound to the active variant's
/// fields, generically.
macro_rules! with_rows {
    ($self:expr, |$modulus:ident, $rows:ident| $body:expr) => {
        match $self {
            SparseRows::U64 {
                modulus: $modulus,
                rows: $rows,
            } => $body,
            SparseRows::U128 {
                modulus: $modulus,
                rows: $rows,
            } => $body,
            SparseRows::U192 {
                modulus: $modulus,
                rows: $rows,
            } => $body,
            SparseRows::U256 {
                modulus: $modulus,
                rows: $rows,
            } => $body,
        }
    };
}

impl SparseRows {
    /// An empty row set for `r`, or `None` when no width supports it
    /// ([`modulus_supported`]).
    pub fn for_modulus(r: &BigUint) -> Option<Self> {
        Some(match WordWidth::for_modulus(r)? {
            WordWidth::U64 => Self::U64 {
                modulus: <u64 as ModWord>::modulus_from_biguint(r)?,
                rows: Vec::new(),
            },
            WordWidth::U128 => Self::U128 {
                modulus: WideModulus::new(r)?,
                rows: Vec::new(),
            },
            WordWidth::U192 => Self::U192 {
                modulus: WideModulus::new(r)?,
                rows: Vec::new(),
            },
            WordWidth::U256 => Self::U256 {
                modulus: WideModulus::new(r)?,
                rows: Vec::new(),
            },
        })
    }

    pub fn width(&self) -> WordWidth {
        match self {
            Self::U64 { .. } => WordWidth::U64,
            Self::U128 { .. } => WordWidth::U128,
            Self::U192 { .. } => WordWidth::U192,
            Self::U256 { .. } => WordWidth::U256,
        }
    }

    pub fn modulus(&self) -> BigUint {
        match self {
            Self::U64 { modulus, .. } => BigUint::from(*modulus),
            Self::U128 { modulus, .. } => modulus.to_biguint(),
            Self::U192 { modulus, .. } => modulus.to_biguint(),
            Self::U256 { modulus, .. } => modulus.to_biguint(),
        }
    }

    pub fn len(&self) -> usize {
        with_rows!(self, |_m, rows| rows.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn clear(&mut self) {
        with_rows!(self, |_m, rows| rows.clear())
    }

    /// Append a row given densely over the columns.
    pub fn push_dense(&mut self, row: &[BigUint], rhs: &BigUint) {
        with_rows!(self, |modulus, rows| rows
            .push(SparseRow::from_dense(row, rhs, *modulus)))
    }

    /// Append a row from `(column, coefficient)` pairs (unsorted,
    /// unreduced, repeats allowed).
    pub fn push_entries(&mut self, entries: &[(u32, BigUint)], rhs: &BigUint) {
        with_rows!(self, |modulus, rows| {
            let m = &*modulus;
            let entries = entries
                .iter()
                .map(|(c, v)| (*c, ModWord::from_biguint(v, m)))
                .collect();
            rows.push(SparseRow::new(entries, ModWord::from_biguint(rhs, m), *m))
        })
    }

    /// Row `i` as `(column, coefficient)` pairs and its right-hand side,
    /// as big integers.
    pub fn row_biguint(&self, i: usize) -> (Vec<(u32, BigUint)>, BigUint) {
        with_rows!(self, |_m, rows| {
            let row = &rows[i];
            (
                row.entries
                    .iter()
                    .map(|&(c, v)| (c, v.to_biguint()))
                    .collect(),
                row.rhs.to_biguint(),
            )
        })
    }

    /// The last row densely over `columns` columns, as big integers.
    pub fn last_dense(&self, columns: usize) -> Option<Vec<BigUint>> {
        with_rows!(self, |_m, rows| rows.last().map(|row| {
            let mut dense = vec![BigUint::zero(); columns];
            for &(c, v) in &row.entries {
                dense[c as usize] = v.to_biguint();
            }
            dense
        }))
    }

    /// How many rows each of `columns` columns occurs in.
    pub fn column_weights(&self, columns: usize) -> Vec<usize> {
        let mut weights = vec![0usize; columns];
        with_rows!(self, |_m, rows| {
            for row in rows.iter() {
                for &(c, _) in &row.entries {
                    weights[c as usize] += 1;
                }
            }
        });
        weights
    }

    /// [`solve_sparse_system`] in the active width, the solution as big
    /// integers.
    pub fn solve(
        &self,
        columns: usize,
        opts: &SparseSolveOptions,
    ) -> (SparseSolveOutcome<BigUint>, SparseSolveReport) {
        with_rows!(self, |modulus, rows| {
            let (outcome, report) = solve_sparse_system(rows, columns, *modulus, opts);
            let outcome = match outcome {
                SparseSolveOutcome::Solved(x) => {
                    SparseSolveOutcome::Solved(x.into_iter().map(ModWord::to_biguint).collect())
                }
                SparseSolveOutcome::Undetermined => SparseSolveOutcome::Undetermined,
                SparseSolveOutcome::Inconsistent => SparseSolveOutcome::Inconsistent,
            };
            (outcome, report)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n;

    const P: u64 = 1_439_393; // r for K_1 / 2^31
    const Q: u64 = (1u64 << 61) - 1; // Mersenne prime, near the width limit

    /// Rank of the coefficient matrix by plain forward elimination in
    /// `u128`, independent of the solver under test.
    fn reference_rank(rows: &[SparseRow], columns: usize, modulus: u64) -> usize {
        let p = modulus as u128;
        let mut m: Vec<Vec<u128>> = rows
            .iter()
            .map(|r| {
                let mut d = vec![0u128; columns];
                for &(c, v) in &r.entries {
                    d[c as usize] = v as u128;
                }
                d
            })
            .collect();
        let mut rank = 0;
        for col in 0..columns {
            let Some(piv) = (rank..m.len()).find(|&r| m[r][col] != 0) else {
                continue;
            };
            m.swap(rank, piv);
            let inv = invmod(m[rank][col] as u64, modulus).unwrap() as u128;
            for r in rank + 1..m.len() {
                if m[r][col] == 0 {
                    continue;
                }
                let f = m[r][col] * inv % p;
                for k in col..columns {
                    m[r][k] = (m[r][k] + p - f * m[rank][k] % p) % p;
                }
            }
            rank += 1;
        }
        rank
    }

    /// The unique solution by the big-integer reference solver, or
    /// `None` when the system is inconsistent or has free columns.
    fn dense_reference(rows: &[SparseRow], columns: usize, modulus: u64) -> Option<Vec<u64>> {
        if reference_rank(rows, columns, modulus) < columns {
            return None;
        }
        let mut m: Vec<Vec<BigUint>> = rows
            .iter()
            .map(|r| {
                let mut d = vec![BigUint::from(0u32); columns];
                for &(c, v) in &r.entries {
                    d[c as usize] = BigUint::from(v);
                }
                d
            })
            .collect();
        let mut b: Vec<BigUint> = rows.iter().map(|r| BigUint::from(r.rhs)).collect();
        let sol = gaussian_eliminate_mod_n(&mut m, &mut b, &BigUint::from(modulus))?;
        let sol: Vec<u64> = sol.iter().map(|v| to_u64_mod(v, modulus)).collect();
        // The reference is rank-checked; re-check every row anyway.
        rows.iter()
            .all(|r| r.evaluate(&sol, modulus) == r.rhs)
            .then_some(sol)
    }

    /// Random relation-shaped system (row weight ≤ `weight`) with a
    /// planted solution.
    fn planted_system(
        rng: &mut StdRng,
        rows: usize,
        columns: usize,
        weight: usize,
        modulus: u64,
    ) -> (Vec<SparseRow>, Vec<u64>) {
        let x: Vec<u64> = (0..columns).map(|_| rng.gen_range(0..modulus)).collect();
        let mut out = Vec::with_capacity(rows);
        for _ in 0..rows {
            let w = rng.gen_range(1..=weight);
            let entries: Vec<(u32, u64)> = (0..w)
                .map(|_| (rng.gen_range(0..columns) as u32, rng.gen_range(1..modulus)))
                .collect();
            let mut row = SparseRow::new(entries, 0, modulus);
            row.rhs = row.evaluate(&x, modulus);
            out.push(row);
        }
        (out, x)
    }

    #[test]
    fn modular_helpers_agree_with_big_integers() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..200 {
            let a = rng.gen_range(0..Q);
            let b = rng.gen_range(0..Q);
            let big = |v: u64| BigUint::from(v);
            assert_eq!(big(mulmod(a, b, Q)), (big(a) * big(b)) % big(Q));
            assert_eq!(big(addmod(a, b, Q)), (big(a) + big(b)) % big(Q));
            assert_eq!(submod(a, b, Q), addmod(a, negmod(b, Q), Q));
            if a != 0 {
                assert_eq!(mulmod(a, invmod(a, Q).unwrap(), Q), 1);
            }
        }
        assert_eq!(invmod(0, P), None);
        assert_eq!(invmod(6, 9), None);
    }

    /// The largest prime below `2⁶³`: the widest modulus the module
    /// admits, where a lazy accumulator absorbs only two products.
    const Q63: u64 = (1u64 << 63) - 25;

    #[test]
    fn division_free_reduction_matches_the_division() {
        let mut rng = StdRng::seed_from_u64(12);
        let mut moduli = vec![2, 3, P, Q, Q63, (1 << 63) - 1, 21_044_858_204_113];
        moduli.extend(
            (2..=63).map(|bits| rng.gen_range(1u64 << (bits - 1)..=u64::MAX >> (64 - bits))),
        );
        for &m in &moduli {
            let red = Reducer::new(m).expect("in range");
            let mut values: Vec<u128> = vec![
                0,
                1,
                m as u128 - 1,
                m as u128,
                (m as u128 - 1) * (m as u128 - 1),
                u64::MAX as u128,
                1u128 << 64,
                u128::MAX,
                u128::MAX - 1,
            ];
            values.extend((0..200).map(|_| rng.gen::<u128>() >> rng.gen_range(0..128)));
            for &v in &values {
                assert_eq!(red.reduce(v) as u128, v % m as u128, "m = {m}, v = {v}");
            }
            for _ in 0..200 {
                let w = rng.gen_range(0..m);
                let t = if rng.gen_bool(0.1) {
                    u64::MAX
                } else {
                    rng.gen()
                };
                let fixed = FixedMul::new(w, m);
                assert!(fixed.ws.is_some());
                assert_eq!(
                    fixed.mul(t, m),
                    mulmod(w, t, m),
                    "m = {m}, w = {w}, t = {t}"
                );
                let a = rng.gen_range(0..m);
                assert_eq!(fixed.sub_from(a, t, m), submod(a, mulmod(w, t, m), m));
            }
            // The accumulation bound, checked against the exact inequality.
            for a_max in [0, 1, m - 1, u64::MAX] {
                let k = red.lazy_terms(a_max) as u128;
                let per = a_max.max(1) as u128 * u64::MAX as u128;
                let fits = |k: u128| {
                    k.checked_mul(per)
                        .and_then(|s| s.checked_add(m as u128 - 1))
                        .is_some()
                };
                assert!(k >= 1 && fits(k), "m = {m}, a_max = {a_max}");
                if k < usize::MAX as u128 {
                    assert!(!fits(k + 1), "m = {m}, a_max = {a_max}: not the largest");
                }
            }
        }
        assert_eq!(Reducer::new(Q63).unwrap().lazy_terms(Q63 - 1), 2);
        assert_eq!(
            Reducer::new(21_044_858_204_113)
                .unwrap()
                .lazy_terms(21_044_858_204_112),
            876_543
        );
        for m in [0, 1, 1 << 63, u64::MAX] {
            assert!(Reducer::new(m).is_none(), "m = {m}");
        }
        // Past 2⁶³ a fixed multiplier falls back to the division.
        let wide = u64::MAX - 58; // prime
        for _ in 0..200 {
            let (w, t) = (rng.gen_range(0..wide), rng.gen());
            let fixed = FixedMul::new(w, wide);
            assert!(fixed.ws.is_none());
            assert_eq!(fixed.mul(t, wide), mulmod(w, t, wide));
            let a = rng.gen_range(0..wide);
            assert_eq!(
                fixed.sub_from(a, t, wide),
                submod(a, mulmod(w, t, wide), wide)
            );
        }
    }

    /// Random CSR rows up to `max_weight` entries, coefficients below
    /// `coeff_bound` (unreduced when it exceeds the modulus; `entries` is
    /// public, so the products must hold for any stored value).
    fn random_csr(
        rng: &mut StdRng,
        rows: usize,
        cols: usize,
        max_weight: usize,
        coeff_bound: u64,
        modulus: u64,
    ) -> CsrMatrix {
        let rows: Vec<SparseRow> = (0..rows)
            .map(|_| {
                let w = rng.gen_range(0..=max_weight);
                let mut entries: Vec<(u32, u64)> = (0..w)
                    .map(|_| (rng.gen_range(0..cols) as u32, rng.gen_range(1..coeff_bound)))
                    .collect();
                entries.sort_unstable_by_key(|&(c, _)| c);
                entries.dedup_by_key(|&mut (c, _)| c);
                SparseRow { entries, rhs: 0 }
            })
            .collect();
        CsrMatrix::from_rows(&rows, cols, modulus)
    }

    #[test]
    fn lazy_block_products_match_per_entry_products() {
        let mut rng = StdRng::seed_from_u64(13);
        let mut chunked = false;
        for &modulus in &[P, 21_044_858_204_113, Q, Q63] {
            for &coeff_bound in &[modulus, u64::MAX] {
                let big_n = rng.gen_range(20..60);
                // Rows up to 40 entries, far past the two products a
                // 63-bit accumulator holds between reductions.
                let a = random_csr(&mut rng, big_n, big_n, 40, coeff_bound, modulus);
                let (_, terms) = a.lazy.expect("modulus in range");
                chunked |= terms < 40;
                for n in 1..=9 {
                    let x: Vec<u64> = (0..big_n * n)
                        .map(|_| match rng.gen_range(0..3) {
                            0 => rng.gen_range(0..modulus),
                            1 => rng.gen(),
                            _ => u64::MAX - rng.gen_range(0..2),
                        })
                        .collect();
                    let mut y = vec![0u64; big_n * n];
                    a.mul_block(&x, n, &mut y);
                    for (i, out) in y.chunks(n).enumerate() {
                        let mut expect = vec![0u64; n];
                        a.row_into_per_entry(i, &x, n, &mut expect);
                        assert_eq!(out, &expect[..], "modulus {modulus}, n {n}, row {i}");
                    }
                    // The homogenised operator, fused, against A·x
                    // followed by the per-entry subtraction of b ⊗ t.
                    let b: Vec<u64> = (0..big_n)
                        .map(|_| match rng.gen_range(0..4) {
                            0 => 0,
                            1 => modulus,
                            2 => rng.gen(),
                            _ => rng.gen_range(0..modulus),
                        })
                        .collect();
                    let t: Vec<u64> = (0..n).map(|_| rng.gen()).collect();
                    let xt: Vec<u64> = x.iter().chain(&t).copied().collect();
                    let op = Homogenised::new(&a, &b);
                    let mut z = vec![1u64; (big_n + 1) * n];
                    op.apply(&xt, n, &mut z);
                    for (i, out) in z[..big_n * n].chunks(n).enumerate() {
                        for (j, &o) in out.iter().enumerate() {
                            let ax = y[i * n + j];
                            let expect = submod(ax, mulmod(b[i], t[j], modulus), modulus);
                            assert_eq!(o, expect, "modulus {modulus}, n {n}, row {i}");
                        }
                    }
                    assert!(z[big_n * n..].iter().all(|&v| v == 0));
                }
            }
        }
        assert!(chunked, "some row outgrew the accumulator");
        // Past 2⁶³ the per-entry path runs; small values keep its sums
        // inside a u64, so it can be checked against plain u128 sums.
        let wide = (1u64 << 63) + 29;
        let a = random_csr(&mut rng, 30, 30, 8, 1 << 40, wide);
        assert!(a.lazy.is_none());
        let n = 3;
        let x: Vec<u64> = (0..30 * n).map(|_| rng.gen_range(0..1 << 20)).collect();
        let mut y = vec![0u64; 30 * n];
        a.mul_block(&x, n, &mut y);
        for (i, out) in y.chunks(n).enumerate() {
            for (j, &o) in out.iter().enumerate() {
                let sum: u128 = (a.row_ptr[i]..a.row_ptr[i + 1])
                    .map(|k| a.vals[k] as u128 * x[a.col_idx[k] as usize * n + j] as u128)
                    .sum();
                assert_eq!(o as u128, sum % wide as u128);
            }
        }
    }

    #[test]
    fn lazy_strided_sums_match_per_product_sums() {
        // The Krylov projection's and `z_times`' shape: reduced
        // coefficients against arbitrary lanes at a stride, from empty to
        // long enough that at 63 bits the sum, and the slice with it, is
        // split at every second product.
        let mut rng = StdRng::seed_from_u64(15);
        for &modulus in &[2, 3, P, 21_044_858_204_113, Q, Q63] {
            let red = Reducer::new(modulus).expect("in range");
            let lazy = (red, red.lazy_terms(modulus - 1));
            for len in 0..24 {
                let coeffs: Vec<u64> = (0..len).map(|_| rng.gen_range(0..modulus)).collect();
                for stride in [1, 4, 5] {
                    let x: Vec<u64> = (0..len * stride + 4)
                        .map(|_| {
                            if rng.gen_bool(0.3) {
                                u64::MAX
                            } else {
                                rng.gen()
                            }
                        })
                        .collect();
                    let expect: Vec<u64> = (0..4)
                        .map(|l| {
                            coeffs.iter().enumerate().fold(0, |acc, (k, &a)| {
                                addmod(acc, mulmod(a, x[k * stride + l], modulus), modulus)
                            })
                        })
                        .collect();
                    let mut four = [0u64; 4];
                    let mut acc = LazyLanes::<u64, 4>::new(&lazy);
                    acc.add_strided(&coeffs, &x, stride);
                    acc.finish(&mut four);
                    assert_eq!(four[..], expect[..], "m = {modulus}, len {len}");
                    let mut one = 0u64;
                    let mut acc = LazyLanes::<u64, 1>::new(&lazy);
                    acc.add_strided(&coeffs, &x, stride);
                    acc.finish(std::slice::from_mut(&mut one));
                    assert_eq!(one, expect[0], "m = {modulus}, len {len}");
                }
            }
        }
    }

    #[test]
    fn block_wiedemann_solves_systems_modulo_a_63_bit_prime() {
        // Two products per accumulator: every row, the homogenising
        // column and the Krylov projection reduce part-way.
        let mut rng = StdRng::seed_from_u64(14);
        for trial in 0..6 {
            let big_n = rng.gen_range(10..60);
            let (rows, x) = planted_system(&mut rng, big_n, big_n, 12, Q63);
            let a = CsrMatrix::from_rows(&rows, big_n, Q63);
            let b: Vec<u64> = rows.iter().map(|r| r.rhs).collect();
            let op = Homogenised::new(&a, &b);
            let nonsingular = dense_reference(&rows, big_n, Q63).is_some();
            let Some((z, _)) =
                block_wiedemann_kernel(&op, &BlockWiedemannOptions::default(), rng.gen())
            else {
                assert!(!nonsingular, "trial {trial}");
                continue;
            };
            let mut img = vec![0u64; big_n + 1];
            op.apply(&z, 1, &mut img);
            assert!(
                img.iter().all(|&v| v == 0),
                "trial {trial}: not a kernel vector"
            );
            if nonsingular {
                let inv = invmod(z[big_n], Q63).expect("t ≠ 0");
                let sol: Vec<u64> = z[..big_n].iter().map(|&v| mulmod(v, inv, Q63)).collect();
                assert_eq!(sol, x, "trial {trial}");
            }
        }
    }

    #[test]
    fn sparse_row_construction_merges_and_reduces() {
        let row = SparseRow::new(vec![(3, P + 2), (1, 5), (3, P - 2), (0, P)], 2 * P + 7, P);
        assert_eq!(row.entries, vec![(1, 5)]);
        assert_eq!(row.rhs, 7);
        let dense = vec![
            BigUint::from(0u32),
            BigUint::from(P + 4),
            BigUint::from(3u32),
        ];
        let from_dense = SparseRow::from_dense(&dense, &BigUint::from(P), P);
        assert_eq!(from_dense.entries, vec![(1, 4), (2, 3)]);
        assert_eq!(from_dense.rhs, 0);
        let a = SparseRow::new(vec![(0, 2), (2, 5)], 9, P);
        let b = SparseRow::new(vec![(0, 1), (1, 7), (2, 5)], 4, P);
        let c = a.axpy(1, &b, P);
        assert_eq!(c.entries, vec![(0, 1), (1, P - 7)]);
        assert_eq!(c.rhs, 5);
    }

    #[test]
    fn dense_u64_solver_matches_the_big_integer_reference() {
        let mut rng = StdRng::seed_from_u64(2);
        for trial in 0..40 {
            let n = rng.gen_range(1..12);
            let extra = rng.gen_range(0..3);
            let (rows, x) = planted_system(&mut rng, n + extra, n, n.min(4), P);
            let dense: Vec<Vec<u64>> = rows
                .iter()
                .map(|r| {
                    let mut d = vec![0u64; n];
                    for &(c, v) in &r.entries {
                        d[c as usize] = v;
                    }
                    d
                })
                .collect();
            let rhs: Vec<u64> = rows.iter().map(|r| r.rhs).collect();
            let ours = dense_solve_unique(dense, rhs, P);
            let reference = dense_reference(&rows, n, P);
            assert_eq!(ours.is_some(), reference.is_some(), "trial {trial}");
            if let Some(sol) = ours {
                assert_eq!(sol, x, "trial {trial}");
            }
        }
    }

    #[test]
    fn filtering_preserves_the_solution_and_reconstructs_every_column() {
        let mut rng = StdRng::seed_from_u64(3);
        let mut merged_any = false;
        let mut purged_any = false;
        for trial in 0..60 {
            let columns = rng.gen_range(5..80);
            let rows_n = columns + rng.gen_range(0..columns);
            let weight = rng.gen_range(2..=4);
            let (rows, x) = planted_system(&mut rng, rows_n, columns, weight, P);
            let opts = FilterOptions {
                target_excess: rng.gen_range(0..8),
                merge_max_weight: rng.gen_range(0..10),
                max_row_weight: rng.gen_range(4..40),
                ..FilterOptions::default()
            };
            let filtered = filter_relations(rows.clone(), columns, P, &opts).expect("consistent");
            let rep = &filtered.report;
            merged_any |= rep.merged_columns > 0;
            purged_any |= rep.excess_rows_removed > 0;
            assert_eq!(rep.rows_in, rows_n);
            assert_eq!(rep.rows_out, filtered.core_rows.len());
            assert_eq!(rep.columns_out, filtered.core_columns.len());
            // The planted solution satisfies every core row (they are
            // combinations of the originals).
            for row in &filtered.core_rows {
                assert_eq!(
                    row.evaluate(&x, P),
                    row.rhs,
                    "trial {trial}: core row broken"
                );
                for &(c, _) in &row.entries {
                    assert!(filtered.core_columns.contains(&c));
                }
            }
            // Reconstruct from the planted core values: exact whenever the
            // original system determines every column.
            let core: Vec<u64> = filtered
                .core_columns
                .iter()
                .map(|&c| x[c as usize])
                .collect();
            let full = filtered.reconstruct(&core, &rows);
            let determined = dense_reference(&rows, columns, P).is_some();
            assert_eq!(full.is_some(), determined, "trial {trial}: determinacy");
            if let Some(full) = full {
                assert_eq!(full, x, "trial {trial}");
            }
        }
        assert!(
            merged_any && purged_any,
            "the random trials exercised merge and purge"
        );
    }

    #[test]
    fn filtering_detects_duplicates_singletons_and_contradictions() {
        let r = |e: Vec<(u32, u64)>, rhs| SparseRow::new(e, rhs, P);
        let rows = vec![
            r(vec![(0, 1), (1, 1)], 3),
            r(vec![(0, 1), (1, 1)], 3), // duplicate
            r(vec![(1, 2), (2, 1)], 4),
            r(vec![(2, 1), (3, 5)], 6), // column 3 is a singleton
            r(vec![(0, 1), (2, 1)], 2),
        ];
        let f = filter_relations(
            rows.clone(),
            4,
            P,
            &FilterOptions {
                merge_max_weight: 0,
                ..FilterOptions::default()
            },
        )
        .unwrap();
        assert_eq!(f.report.duplicates_removed, 1);
        assert_eq!(f.report.singletons_removed, 1);
        assert_eq!(f.core_columns, vec![0, 1, 2]);
        assert_eq!(f.core_rows.len(), 3);
        // A contradictory pair.
        let bad = vec![r(vec![(0, 1)], 1), r(vec![(0, 1)], 2)];
        assert!(filter_relations(bad, 1, P, &FilterOptions::default()).is_err());
        let empty = vec![r(vec![], 5)];
        assert!(filter_relations(empty, 1, P, &FilterOptions::default()).is_err());
    }

    #[test]
    fn approximant_basis_reproduces_scalar_berlekamp_massey() {
        // A scalar linear recurrence of order 5: the basis must contain
        // a column of degree 5 whose f is the connection polynomial.
        let modulus = P;
        let coeffs = [3u64, 7, 11, 13, 17];
        let mut s: Vec<u64> = vec![1, 2, 3, 4, 5];
        for t in 5..40 {
            let mut v = 0;
            for (k, &c) in coeffs.iter().enumerate() {
                v = addmod(v, mulmod(c, s[t - 1 - k], modulus), modulus);
            }
            s.push(v);
        }
        let seq: Vec<Vec<u64>> = s.iter().map(|&v| vec![v]).collect();
        let basis = minimal_approximant_basis(&seq, 1, 1, modulus);
        let best = basis.iter().min_by_key(|c| c.delta).unwrap();
        assert_eq!(best.delta, 5);
        // f(λ) = f_0 + f_1 λ + … with Σ_k s_{t−k} f_k = 0 for t ≥ 5.
        for t in 5..40 {
            let mut acc = 0;
            for (k, fk) in best.f.iter().enumerate() {
                acc = addmod(acc, mulmod(s[t - k], fk[0], modulus), modulus);
            }
            assert_eq!(acc, 0, "recurrence at t = {t}");
        }
        // And it is the true connection polynomial up to scaling.
        let scale = invmod(best.f[0][0], modulus).unwrap();
        for (k, &c) in coeffs.iter().enumerate() {
            assert_eq!(mulmod(best.f[k + 1][0], scale, modulus), negmod(c, modulus));
        }
    }

    #[test]
    fn approximant_basis_satisfies_the_block_recurrence() {
        let mut rng = StdRng::seed_from_u64(4);
        for &(m, n, big_n) in &[
            (1usize, 1usize, 9usize),
            (2, 2, 13),
            (3, 2, 20),
            (2, 4, 17),
            (4, 4, 40),
        ] {
            // A random dense square matrix; the sequence X Mᵗ Y.
            let mat: Vec<u64> = (0..big_n * big_n).map(|_| rng.gen_range(0..P)).collect();
            let x: Vec<u64> = (0..m * big_n).map(|_| rng.gen_range(0..P)).collect();
            let y: Vec<u64> = (0..big_n * n).map(|_| rng.gen_range(0..P)).collect();
            let big_l = big_n.div_ceil(m) + big_n.div_ceil(n) + 6;
            let mut v = y.clone();
            let mut seq = Vec::new();
            for _ in 0..big_l {
                let mut s = vec![0u64; m * n];
                for i in 0..m {
                    for j in 0..n {
                        for r in 0..big_n {
                            s[i * n + j] =
                                addmod(s[i * n + j], mulmod(x[i * big_n + r], v[r * n + j], P), P);
                        }
                    }
                }
                seq.push(s);
                let mut w = vec![0u64; big_n * n];
                for r in 0..big_n {
                    for j in 0..n {
                        for c in 0..big_n {
                            w[r * n + j] = addmod(
                                w[r * n + j],
                                mulmod(mat[r * big_n + c], v[c * n + j], P),
                                P,
                            );
                        }
                    }
                }
                v = w;
            }
            let basis = minimal_approximant_basis(&seq, m, n, P);
            assert_eq!(basis.len(), m + n);
            let total: usize = basis.iter().map(|c| c.delta).sum();
            assert!(total <= m + m * big_l, "degree budget");
            for col in &basis {
                for t in col.delta..big_l {
                    let mut acc = vec![0u64; m];
                    for (k, fk) in col.f.iter().enumerate() {
                        if k > t {
                            break;
                        }
                        for i in 0..m {
                            for j in 0..n {
                                acc[i] = addmod(acc[i], mulmod(seq[t - k][i * n + j], fk[j], P), P);
                            }
                        }
                    }
                    assert!(
                        acc.iter().all(|&a| a == 0),
                        "m={m} n={n}: recurrence fails at t={t}"
                    );
                }
            }
            // The n lowest-degree columns carry the Kalman degree N.
            let mut deltas: Vec<usize> = basis.iter().map(|c| c.delta).collect();
            deltas.sort_unstable();
            let low: usize = deltas[..n].iter().sum();
            assert!(
                low >= big_n && low <= big_n + n,
                "m={m} n={n}: low degrees {deltas:?}"
            );
        }
    }

    #[test]
    fn block_wiedemann_solves_random_sparse_systems() {
        let mut rng = StdRng::seed_from_u64(5);
        for &(m, n) in &[(1usize, 1usize), (2, 2), (4, 2), (2, 4), (4, 4), (8, 8)] {
            for &modulus in &[P, Q] {
                for trial in 0..4 {
                    let big_n = rng.gen_range(1..70);
                    let (rows, x) = planted_system(&mut rng, big_n, big_n, 4, modulus);
                    let a = CsrMatrix::from_rows(&rows, big_n, modulus);
                    let b: Vec<u64> = rows.iter().map(|r| r.rhs).collect();
                    let nonsingular = dense_reference(&rows, big_n, modulus).is_some();
                    let op = Homogenised::new(&a, &b);
                    let opts = BlockWiedemannOptions {
                        block_m: m,
                        block_n: n,
                        margin: 8,
                    };
                    let found = block_wiedemann_kernel(&op, &opts, rng.gen());
                    let Some((z, report)) = found else {
                        // Only a singular A can defeat the kernel search.
                        assert!(
                            !nonsingular,
                            "m={m} n={n} trial {trial}: no kernel vector for a nonsingular system"
                        );
                        continue;
                    };
                    assert!(z.iter().any(|&v| v != 0));
                    // z is in the kernel of the homogenised operator.
                    let mut img = vec![0u64; big_n + 1];
                    op.apply(&z, 1, &mut img);
                    assert!(
                        img.iter().all(|&v| v == 0),
                        "m={m} n={n}: not a kernel vector"
                    );
                    if nonsingular {
                        let t = z[big_n];
                        assert_ne!(t, 0, "kernel vector must have t ≠ 0 for nonsingular A");
                        let inv = invmod(t, modulus).unwrap();
                        let sol: Vec<u64> = z[..big_n]
                            .iter()
                            .map(|&v| mulmod(v, inv, modulus))
                            .collect();
                        assert_eq!(sol, x, "m={m} n={n} trial {trial} dim {big_n}");
                    }
                    assert_eq!(report.dimension, big_n + 1);
                    assert!(
                        report.sequence_length >= (big_n + 1).div_ceil(m) + (big_n + 1).div_ceil(n)
                    );
                }
            }
        }
    }

    #[test]
    fn block_wiedemann_finds_kernel_vectors_of_singular_matrices() {
        let mut rng = StdRng::seed_from_u64(6);
        for trial in 0..10 {
            let big_n = rng.gen_range(3..40);
            // Rank-deficient: the last row is a combination of two others.
            let (mut rows, _) = planted_system(&mut rng, big_n, big_n, 3, P);
            let (i, j) = (rng.gen_range(0..big_n - 1), rng.gen_range(0..big_n - 1));
            let combo = rows[i].axpy(rng.gen_range(1..P), &rows[j], P);
            let last = rows.len() - 1;
            rows[last] = combo;
            let a = CsrMatrix::from_rows(&rows, big_n, P);
            let opts = BlockWiedemannOptions::default();
            let mut found = false;
            for attempt in 0..3 {
                if let Some((z, _)) = block_wiedemann_kernel(&a, &opts, trial * 10 + attempt) {
                    assert!(z.iter().any(|&v| v != 0));
                    let mut img = vec![0u64; big_n];
                    a.apply(&z, 1, &mut img);
                    assert!(
                        img.iter().all(|&v| v == 0),
                        "trial {trial}: not in the kernel"
                    );
                    found = true;
                    break;
                }
            }
            assert!(
                found,
                "trial {trial}: singular matrix of dimension {big_n} yielded no kernel vector"
            );
        }
    }

    #[test]
    fn parallel_and_serial_block_products_agree() {
        let mut rng = StdRng::seed_from_u64(7);
        let big_n = PARALLEL_ROWS + 17;
        let (rows, _) = planted_system(&mut rng, big_n, big_n, 4, P);
        let a = CsrMatrix::from_rows(&rows, big_n, P);
        let n = 3;
        let x: Vec<u64> = (0..big_n * n).map(|_| rng.gen_range(0..P)).collect();
        let mut y = vec![0u64; big_n * n];
        a.mul_block(&x, n, &mut y);
        for (i, out) in y.chunks(n).enumerate() {
            let mut expect = vec![0u64; n];
            a.row_into(i, &x, n, &mut expect);
            assert_eq!(out, &expect[..], "row {i}");
        }
        assert_eq!(a.nnz(), rows.iter().map(SparseRow::weight).sum::<usize>());
    }

    #[test]
    fn sparse_solve_matches_dense_on_relation_shaped_systems() {
        let mut rng = StdRng::seed_from_u64(8);
        let mut solved = 0;
        let mut undetermined = 0;
        for trial in 0..50 {
            let columns = rng.gen_range(1..120);
            let rows_n = columns + rng.gen_range(0..columns + 4);
            let weight = rng.gen_range(2..=4);
            let modulus = if trial % 2 == 0 { P } else { Q };
            let (rows, x) = planted_system(&mut rng, rows_n, columns, weight, modulus);
            let opts = SparseSolveOptions {
                filter: FilterOptions {
                    target_excess: rng.gen_range(0..40),
                    merge_max_weight: rng.gen_range(0..10),
                    ..FilterOptions::default()
                },
                wiedemann: BlockWiedemannOptions {
                    block_m: rng.gen_range(1..5),
                    block_n: rng.gen_range(1..5),
                    margin: 8,
                },
                seed: trial,
                ..SparseSolveOptions::default()
            };
            let (outcome, report) = solve_sparse_system(&rows, columns, modulus, &opts);
            let reference = dense_reference(&rows, columns, modulus);
            match outcome {
                SparseSolveOutcome::Solved(sol) => {
                    assert_eq!(sol, x, "trial {trial}");
                    assert!(
                        reference.is_some(),
                        "trial {trial}: sparse solved an underdetermined system?"
                    );
                    assert_eq!(report.filter.rows_in, rows_n);
                    assert_eq!(
                        report.core_dimension + report.reconstructed_columns,
                        columns
                    );
                    solved += 1;
                }
                SparseSolveOutcome::Undetermined => {
                    // The sparse path may fail to converge on a determined
                    // system only through the random fold; with three
                    // attempts that should be rare, so demand agreement.
                    assert!(
                        reference.is_none(),
                        "trial {trial}: dense determined it but sparse gave up ({report:?})"
                    );
                    undetermined += 1;
                }
                SparseSolveOutcome::Inconsistent => {
                    panic!("trial {trial}: consistent system reported inconsistent")
                }
            }
        }
        assert!(
            solved >= 10 && undetermined >= 1,
            "solved {solved}, undetermined {undetermined}"
        );
    }

    #[test]
    fn sparse_solve_flags_a_wrong_relation() {
        let mut rng = StdRng::seed_from_u64(9);
        let (mut rows, _) = planted_system(&mut rng, 40, 20, 3, P);
        // Make sure it is determined before corrupting.
        let (before, _) = solve_sparse_system(&rows, 20, P, &SparseSolveOptions::default());
        assert!(matches!(before, SparseSolveOutcome::Solved(_)));
        rows[7].rhs = addmod(rows[7].rhs, 1, P);
        let (after, _) = solve_sparse_system(&rows, 20, P, &SparseSolveOptions::default());
        assert!(
            matches!(
                after,
                SparseSolveOutcome::Inconsistent | SparseSolveOutcome::Undetermined
            ),
            "a corrupted relation must not yield a verified solution: {after:?}"
        );
    }

    #[test]
    fn sparse_solve_handles_a_larger_system() {
        // Rows of weight exactly 3 with a sixfold excess, so every column
        // is covered many times over and the purge has work to do.
        let mut rng = StdRng::seed_from_u64(10);
        let columns = 3000usize;
        let x: Vec<u64> = (0..columns).map(|_| rng.gen_range(0..P)).collect();
        let rows: Vec<SparseRow> = (0..6 * columns)
            .map(|_| {
                let entries: Vec<(u32, u64)> = (0..3)
                    .map(|_| (rng.gen_range(0..columns) as u32, rng.gen_range(1..P)))
                    .collect();
                let mut row = SparseRow::new(entries, 0, P);
                row.rhs = row.evaluate(&x, P);
                row
            })
            .collect();
        let begin = std::time::Instant::now();
        let (outcome, report) =
            solve_sparse_system(&rows, columns, P, &SparseSolveOptions::default());
        println!(
            "3000-column sparse solve in {:.3}s: {report:?}",
            begin.elapsed().as_secs_f64()
        );
        match outcome {
            SparseSolveOutcome::Solved(sol) => assert_eq!(sol, x),
            other => panic!("expected a solution: {other:?} ({report:?})"),
        }
        assert!(
            report.filter.excess_rows_removed > 0,
            "purge ran: {report:?}"
        );
        assert!(
            report.core_dimension < columns,
            "filtering shrank the core: {report:?}"
        );
        assert_eq!(
            report.core_dimension + report.reconstructed_columns,
            columns
        );
        assert!(report.wiedemann.is_some());
    }

    #[test]
    fn block_lanczos_solves_and_verifies_a_sparse_prime_field_system() {
        let columns = 12usize;
        let planted: Vec<u64> = (0..columns).map(|i| (17 * i as u64 + 9) % P).collect();
        // A cyclic, nonsymmetric three-diagonal matrix.  Lanczos sees only
        // A^T D A; the final assertion checks the original equations.
        let mut rows = Vec::new();
        for i in 0..columns {
            let entries = vec![
                (i as u32, 3),
                (((i + 1) % columns) as u32, 5),
                (((i + 3) % columns) as u32, 7),
            ];
            let mut row = SparseRow::new(entries, 0, P);
            row.rhs = row.evaluate(&planted, P);
            rows.push(row);
        }
        assert!(dense_reference(&rows, columns, P).is_some());
        let matrix = CsrMatrix::from_rows(&rows, columns, P);
        let rhs: Vec<u64> = rows.iter().map(|row| row.rhs).collect();
        let options = BlockLanczosOptions {
            block_size: 3,
            margin: 12,
        };
        let mut found = None;
        for seed in 0..64 {
            if let Some(result) = block_lanczos_solve(&matrix, &rhs, &options, seed) {
                found = Some(result);
                break;
            }
        }
        let (solution, report) =
            found.expect("a randomized block-Lanczos draw must not break down");
        assert_eq!(solution, planted);
        assert!(report.iterations > 0);
        assert!(rows.iter().all(|row| row.evaluate(&solution, P) == row.rhs));
    }

    #[test]
    fn new_sparse_defaults_preserve_the_frozen_json_shape() {
        let defaults = SparseSolveOptions::default();
        let value = serde_json::to_value(defaults).unwrap();
        assert!(value.get("solver").is_none());
        assert!(value.get("spmv").is_none());
        assert!(value.get("lanczos").is_none());
        assert_eq!(
            serde_json::from_value::<SparseSolveOptions>(value).unwrap(),
            defaults
        );

        let selected = SparseSolveOptions {
            solver: SparseCoreSolver::BlockLanczos,
            spmv: SpmvOptions {
                backend: SpmvBackend::Sharded,
                shards: 8,
            },
            lanczos: BlockLanczosOptions {
                block_size: 8,
                ..BlockLanczosOptions::default()
            },
            ..defaults
        };
        let selected_json = serde_json::to_value(selected).unwrap();
        assert_eq!(selected_json["solver"], "block-lanczos");
        assert_eq!(selected_json["spmv"]["backend"], "sharded");
        assert_eq!(selected_json["lanczos"]["block_size"], 8);
        assert_eq!(
            serde_json::from_value::<SparseSolveOptions>(selected_json).unwrap(),
            selected
        );

        let report = serde_json::to_value(SparseSolveReport::default()).unwrap();
        assert!(report.get("solver").is_none());
        assert!(report.get("spmv").is_none());
        assert!(report.get("lanczos").is_none());
    }

    #[test]
    fn serial_rayon_and_sharded_spmv_are_identical() {
        let rows = vec![
            SparseRow::new(vec![(0, 3), (2, 5)], 0, P),
            SparseRow::new(vec![(1, 7), (3, 11)], 0, P),
            SparseRow::new(vec![(0, 13), (3, 17)], 0, P),
            SparseRow::new(vec![(1, 19), (2, 23)], 0, P),
        ];
        let x: Vec<u64> = (0..4 * 3).map(|i| (29 * i as u64 + 1) % P).collect();
        let multiply = |backend| {
            let matrix =
                CsrMatrix::from_rows(&rows, 4, P).with_spmv(SpmvOptions { backend, shards: 3 });
            let mut y = vec![0u64; 4 * 3];
            matrix.mul_block(&x, 3, &mut y);
            y
        };
        let serial = multiply(SpmvBackend::Serial);
        assert_eq!(multiply(SpmvBackend::Rayon), serial);
        assert_eq!(multiply(SpmvBackend::Sharded), serial);
    }

    #[test]
    fn configured_spmv_worker_matches_the_cpu_reference() {
        if std::env::var_os("IC_SPMV_WORKER").is_none() {
            return;
        }
        let rows = vec![
            SparseRow::new(vec![(0, 91), (2, 17)], 0, P),
            SparseRow::new(vec![(1, 31), (2, 43)], 0, P),
            SparseRow::new(vec![(0, 59), (1, 61)], 0, P),
        ];
        let x: Vec<u64> = (0..3 * 4).map(|i| (67 * i as u64 + 5) % P).collect();
        let multiply = |backend| {
            let matrix =
                CsrMatrix::from_rows(&rows, 3, P).with_spmv(SpmvOptions { backend, shards: 0 });
            let mut y = vec![0u64; 3 * 4];
            matrix.mul_block(&x, 4, &mut y);
            y
        };
        assert_eq!(multiply(SpmvBackend::Worker), multiply(SpmvBackend::Serial));
    }

    #[test]
    fn sparse_pipeline_selects_block_lanczos() {
        let columns = 12usize;
        let planted: Vec<u64> = (0..columns).map(|i| (31 * i as u64 + 4) % P).collect();
        let mut rows = Vec::new();
        for i in 0..columns {
            let entries = vec![
                (i as u32, 3),
                (((i + 1) % columns) as u32, 5),
                (((i + 3) % columns) as u32, 7),
            ];
            let mut row = SparseRow::new(entries, 0, P);
            row.rhs = row.evaluate(&planted, P);
            rows.push(row);
        }
        assert!(dense_reference(&rows, columns, P).is_some());
        let options = SparseSolveOptions {
            solver: SparseCoreSolver::BlockLanczos,
            lanczos: BlockLanczosOptions {
                block_size: 3,
                margin: 8,
            },
            attempts: 64,
            ..SparseSolveOptions::default()
        };
        let (outcome, report) = solve_sparse_system(&rows, columns, P, &options);
        assert_eq!(outcome, SparseSolveOutcome::Solved(planted));
        assert_eq!(report.solver, SparseCoreSolver::BlockLanczos);
        assert_eq!(report.spmv, options.spmv);
        // Filtering is allowed to solve the entire system before a Krylov
        // backend is needed; otherwise the selected backend must be Lanczos.
        assert!(report.core_dimension == 0 || report.lanczos.is_some());
        assert!(report.wiedemann.is_none());
    }

    /// Timing on a relation-shaped system at a chosen size, against the
    /// dense reference where that is affordable.  Run with
    /// `KOBLITZ_SPARSE_LA_COLUMNS=20000 cargo test --release --lib
    /// sparse_solve_benchmark -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn sparse_solve_benchmark() {
        let columns: usize = std::env::var("KOBLITZ_SPARSE_LA_COLUMNS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(20_000);
        let weight: usize = std::env::var("KOBLITZ_SPARSE_LA_WEIGHT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(3);
        let excess = columns / 10 + 64;
        let mut rng = StdRng::seed_from_u64(11);
        let x: Vec<u64> = (0..columns).map(|_| rng.gen_range(0..P)).collect();
        // Cover every column at least once, then random rows of the weight.
        let mut rows: Vec<SparseRow> = Vec::with_capacity(columns + excess);
        for c in 0..columns as u32 {
            let mut entries = vec![(c, rng.gen_range(1..P))];
            for _ in 1..weight {
                entries.push((rng.gen_range(0..columns) as u32, rng.gen_range(1..P)));
            }
            let mut row = SparseRow::new(entries, 0, P);
            row.rhs = row.evaluate(&x, P);
            rows.push(row);
        }
        for _ in 0..excess {
            let entries: Vec<(u32, u64)> = (0..weight)
                .map(|_| (rng.gen_range(0..columns) as u32, rng.gen_range(1..P)))
                .collect();
            let mut row = SparseRow::new(entries, 0, P);
            row.rhs = row.evaluate(&x, P);
            rows.push(row);
        }
        for block in [1usize, 2, 4, 8] {
            let opts = SparseSolveOptions {
                wiedemann: BlockWiedemannOptions {
                    block_m: block,
                    block_n: block,
                    margin: 8,
                },
                ..SparseSolveOptions::default()
            };
            let begin = std::time::Instant::now();
            let (outcome, report) = solve_sparse_system(&rows, columns, P, &opts);
            let secs = begin.elapsed().as_secs_f64();
            let ok = matches!(&outcome, SparseSolveOutcome::Solved(sol) if *sol == x);
            println!(
                "columns {columns} rows {} weight {weight} block {block}×{block}: {} in {secs:.3}s; core {} (nnz {}), attempts {}, {:?}",
                rows.len(),
                if ok { "solved" } else { "FAILED" },
                report.core_dimension,
                report.core_nonzeros,
                report.attempts,
                report.wiedemann
            );
            assert!(ok, "{outcome:?}");
        }
        if columns <= 3000 {
            let begin = std::time::Instant::now();
            let reference = dense_reference(&rows, columns, P);
            println!(
                "dense big-integer reference: {:.3}s",
                begin.elapsed().as_secs_f64()
            );
            assert_eq!(reference.as_ref(), Some(&x));
        }
    }
}

/// The wide word widths: arithmetic against `num-bigint`, and the whole
/// pipeline (filtering, block Wiedemann, block Lanczos, dense fallback)
/// modulo primes of 61, 127, 130 (the ECC2K-130 subgroup order), 190
/// and 250 bits, on relation-shaped systems with planted solutions.
#[cfg(test)]
mod wide_tests {
    use super::*;
    use num_traits::{One, Zero};
    use std::str::FromStr;

    /// `2⁶¹ − 1`: the `u64` path's regression fixture, also run in `U128`.
    const P61: &str = "2305843009213693951";
    /// `2¹²⁷ − 1`.
    const P127: &str = "170141183460469231731687303715884105727";
    /// The ECC2K-130 subgroup order (prime, 130 bits).
    const R130: &str = "680564733841876926932320129493409985129";
    /// `2¹⁹⁰ − 11`.
    const P190: &str = "1569275433846670190958947355801916604025588861116008628213";
    /// `2¹⁹¹ − 19`: the top of `U192`, where an accumulator absorbs two products.
    const P191: &str = "3138550867693340381917894711603833208051177722232017256429";
    /// `2²⁵⁰ − 207`.
    const P250: &str =
        "1809251394333065553493296640760748560207343510400633813116524750123642650417";

    fn big(s: &str) -> BigUint {
        BigUint::from_str(s).unwrap()
    }

    fn modulus<W: ModWord>(s: &str) -> W::Modulus {
        W::modulus_from_biguint(&big(s)).expect("modulus fits the width")
    }

    /// `rows` relation-shaped rows of exactly `weight` random columns
    /// (merged if they collide), the first `columns` rows each covering
    /// their own column, with a planted solution.
    fn planted<W: ModWord>(
        rng: &mut StdRng,
        rows: usize,
        columns: usize,
        weight: usize,
        m: &W::Modulus,
    ) -> (Vec<SparseRow<W>>, Vec<W>) {
        let x: Vec<W> = (0..columns).map(|_| W::random(rng, m)).collect();
        let mut out = Vec::with_capacity(rows);
        for i in 0..rows {
            let mut entries: Vec<(u32, W)> = Vec::with_capacity(weight);
            if i < columns {
                entries.push((i as u32, W::random_nonzero(rng, m)));
            }
            while entries.len() < weight {
                entries.push((rng.gen_range(0..columns) as u32, W::random_nonzero(rng, m)));
            }
            let mut row = SparseRow::new(entries, W::zero(), *m);
            row.rhs = row.evaluate(&x, *m);
            out.push(row);
        }
        (out, x)
    }

    fn check_arithmetic<W: ModWord>(s: &str, rng: &mut StdRng) {
        let m = modulus::<W>(s);
        let mb = W::modulus_to_biguint(&m);
        assert_eq!(mb, big(s));
        let lazy = W::lazy(&m, W::max_residue(&m)).expect("wide widths always have a lazy path");
        let terms = W::lazy_terms(&lazy);
        assert!(terms >= 2, "{s}: lazy terms {terms}");
        assert_eq!(W::max_residue(&m).to_biguint(), &mb - BigUint::one());
        for _ in 0..300 {
            let a = W::random(rng, &m);
            let b = W::random(rng, &m);
            let (ab, bb) = (a.to_biguint(), b.to_biguint());
            assert!(ab < mb && bb < mb);
            assert_eq!(W::from_biguint(&ab, &m), a);
            assert_eq!(a.add(b, &m).to_biguint(), (&ab + &bb) % &mb);
            assert_eq!(a.sub(b, &m).to_biguint(), (&ab + &mb - &bb) % &mb);
            assert_eq!(a.neg(&m).to_biguint(), (&mb - &ab) % &mb);
            assert_eq!(a.mul(b, &m).to_biguint(), (&ab * &bb) % &mb, "{s}: mul");
            let f = W::fixed(a, &m);
            assert_eq!(W::fixed_mul(&f, b, &m).to_biguint(), (&ab * &bb) % &mb);
            assert_eq!(
                W::fixed_sub_from(&f, b, a, &m).to_biguint(),
                (&bb + &mb - (&ab * &ab) % &mb) % &mb
            );
            let al = W::lazy_form(a, &m);
            assert_eq!(W::from_lazy_form(al, &m), a);
            assert_eq!(W::lazy_mul(al, b, &m).to_biguint(), (&ab * &bb) % &mb);
            if !a.is_zero() {
                let inv = a.inv(&m).expect("prime modulus");
                assert_eq!(a.mul(inv, &m), W::one(), "{s}: inverse");
            }
            // Unreduced input to `reduce` and `from_u64`.
            let u: u64 = rng.gen();
            assert_eq!(W::from_u64(u, &m).to_biguint(), BigUint::from(u) % &mb);
        }
        assert!(W::zero().inv(&m).is_none());
        // Exact accumulation of the longest sums the module forms (a
        // row of up to 40 entries, well past `terms` for the top-of-width
        // primes), both coefficient forms, against big-integer sums.
        for len in [0usize, 1, 2, 3, 7, 40] {
            let coeffs: Vec<W> = (0..len).map(|_| W::random(rng, &m)).collect();
            let data: Vec<W> = (0..len).map(|_| W::random(rng, &m)).collect();
            let expect = coeffs
                .iter()
                .zip(&data)
                .fold(BigUint::zero(), |acc, (c, d)| {
                    acc + c.to_biguint() * d.to_biguint()
                })
                % &mb;
            let lazy_coeffs: Vec<W> = coeffs.iter().map(|&c| W::lazy_form(c, &m)).collect();
            let mut acc = LazyLanes::<W, 1>::new(&lazy);
            for (c, d) in lazy_coeffs.iter().zip(&data) {
                acc.add(*c, std::slice::from_ref(d));
            }
            let mut out = [W::zero()];
            acc.finish(&mut out);
            assert_eq!(out[0].to_biguint(), expect, "{s}: lazy sum of {len}");
            let mut acc = LazyLanes::<W, 1>::plain(&lazy);
            for (c, d) in coeffs.iter().zip(&data) {
                acc.add(*c, std::slice::from_ref(d));
            }
            acc.finish(&mut out);
            assert_eq!(out[0].to_biguint(), expect, "{s}: plain sum of {len}");
            // Four lanes at a stride, the Krylov projection's shape.
            let x: Vec<W> = (0..len * 4).map(|_| W::random(rng, &m)).collect();
            let mut acc = LazyLanes::<W, 4>::new(&lazy);
            acc.add_strided(&lazy_coeffs, &x, 4);
            let mut four = [W::zero(); 4];
            acc.finish(&mut four);
            for (l, o) in four.iter().enumerate() {
                let expect = coeffs
                    .iter()
                    .enumerate()
                    .fold(BigUint::zero(), |acc, (k, c)| {
                        acc + c.to_biguint() * x[k * 4 + l].to_biguint()
                    })
                    % &mb;
                assert_eq!(o.to_biguint(), expect, "{s}: strided lane {l} of {len}");
            }
        }
    }

    #[test]
    fn wide_arithmetic_agrees_with_big_integers() {
        let mut rng = StdRng::seed_from_u64(0x57_49_44_45);
        check_arithmetic::<U128>(P61, &mut rng);
        check_arithmetic::<U128>(P127, &mut rng);
        check_arithmetic::<U192>(R130, &mut rng);
        check_arithmetic::<U192>(P190, &mut rng);
        check_arithmetic::<U192>(P191, &mut rng);
        check_arithmetic::<U256>(P250, &mut rng);
        check_arithmetic::<U256>("1439393", &mut rng);
        // Top-of-width primes leave an accumulator two products.
        let m = modulus::<U192>(P191);
        assert_eq!(U192::lazy_terms(&U192::lazy(&m, U192::zero()).unwrap()), 2);
        let m = modulus::<U128>(P127);
        assert_eq!(U128::lazy_terms(&U128::lazy(&m, U128::zero()).unwrap()), 2);
        // The ECC2K-130 order leaves it over 2⁶² of them.
        let m = modulus::<U192>(R130);
        assert!(U192::lazy_terms(&U192::lazy(&m, U192::zero()).unwrap()) > 1 << 61);
    }

    #[test]
    fn word_width_selection() {
        let two = BigUint::from(2u32);
        assert_eq!(WordWidth::for_modulus(&two), Some(WordWidth::U64));
        assert_eq!(WordWidth::for_modulus(&BigUint::one()), None);
        assert_eq!(WordWidth::for_modulus(&big(P61)), Some(WordWidth::U64));
        assert_eq!(
            WordWidth::for_modulus(&((BigUint::one() << 63u32) - BigUint::from(25u32))),
            Some(WordWidth::U64)
        );
        // Past 2⁶³ the wide widths need an odd modulus.
        assert_eq!(
            WordWidth::for_modulus(&((BigUint::one() << 63u32) + BigUint::from(29u32))),
            Some(WordWidth::U128)
        );
        assert_eq!(WordWidth::for_modulus(&(BigUint::one() << 64u32)), None);
        assert_eq!(WordWidth::for_modulus(&big(P127)), Some(WordWidth::U128));
        assert_eq!(WordWidth::for_modulus(&big(R130)), Some(WordWidth::U192));
        assert_eq!(big(R130).bits(), 130);
        assert_eq!(WordWidth::for_modulus(&big(P190)), Some(WordWidth::U192));
        assert_eq!(WordWidth::for_modulus(&big(P191)), Some(WordWidth::U192));
        assert_eq!(WordWidth::for_modulus(&big(P250)), Some(WordWidth::U256));
        assert_eq!(
            WordWidth::for_modulus(&((BigUint::one() << 255u32) - BigUint::from(19u32))),
            Some(WordWidth::U256)
        );
        assert_eq!(
            WordWidth::for_modulus(&((BigUint::one() << 255u32) + BigUint::from(95u32))),
            None
        );
        assert!(modulus_supported(&big(R130)));
        assert!(!modulus_supported(&(BigUint::one() << 256u32)));
        // A width refuses what it cannot hold.
        assert!(WideModulus::<2>::new(&big(R130)).is_none());
        assert!(WideModulus::<3>::new(&(BigUint::one() << 191u32)).is_none());
        assert!(WideModulus::<3>::new(&(big(P190) + BigUint::one())).is_none());
    }

    /// Filtering + block Wiedemann, filtering + block Lanczos, and the
    /// dense solver all recover the planted solution; the verified
    /// outcome is `Solved` only then.
    fn check_pipeline<W: ModWord>(s: &str, columns: usize, seed: u64) -> SparseSolveReport {
        let m = modulus::<W>(s);
        let mut rng = StdRng::seed_from_u64(seed);
        let rows_n = columns + columns / 10 + 64;
        let (rows, x) = planted::<W>(&mut rng, rows_n, columns, 5, &m);
        let mut last = None;
        for solver in [
            SparseCoreSolver::BlockWiedemann,
            SparseCoreSolver::BlockLanczos,
        ] {
            let opts = SparseSolveOptions {
                solver,
                seed,
                ..SparseSolveOptions::default()
            };
            let (outcome, report) = solve_sparse_system(&rows, columns, m, &opts);
            match outcome {
                SparseSolveOutcome::Solved(sol) => assert_eq!(sol, x, "{s} {solver:?}"),
                other => panic!("{s} {solver:?}: {other:?} ({report:?})"),
            }
            assert_eq!(
                report.core_dimension + report.reconstructed_columns,
                columns
            );
            match solver {
                SparseCoreSolver::BlockWiedemann => assert!(report.wiedemann.is_some()),
                SparseCoreSolver::BlockLanczos => assert!(report.lanczos.is_some()),
            }
            last = Some(report);
        }
        // The dense solver on a small instance.
        let small = 120.min(columns);
        let (rows_s, x_s) = planted::<W>(&mut rng, small + 8, small, 5, &m);
        let dense: Vec<Vec<W>> = rows_s
            .iter()
            .map(|r| {
                let mut d = vec![W::zero(); small];
                for &(c, v) in &r.entries {
                    d[c as usize] = v;
                }
                d
            })
            .collect();
        let rhs: Vec<W> = rows_s.iter().map(|r| r.rhs).collect();
        assert_eq!(
            dense_solve_unique(dense, rhs, m),
            Some(x_s.clone()),
            "{s} dense"
        );
        // The pipeline on the small instance agrees with the dense solve.
        let (outcome, _) = solve_sparse_system(&rows_s, small, m, &SparseSolveOptions::default());
        assert_eq!(
            outcome,
            SparseSolveOutcome::Solved(x_s),
            "{s} small pipeline"
        );
        // Singular: a column no row covers, and a rank-deficient system.
        let (mut rows_u, _) = planted::<W>(&mut rng, small + 8, small, 5, &m);
        rows_u.retain(|r| r.entries.iter().all(|&(c, _)| c != 3));
        let (outcome, report) =
            solve_sparse_system(&rows_u, small, m, &SparseSolveOptions::default());
        assert_eq!(outcome, SparseSolveOutcome::Undetermined, "{s} uncovered");
        assert_eq!(report.filter.uncovered_columns, 1);
        // Rank-deficient, consistent, every column covered: a square
        // system whose row 2 is a combination of rows 0 and 1.  The dense
        // solver reports it (no unique solution).  The sparse pipeline
        // reports `Undetermined` when filtering exposes the dependency
        // (the row reduces to `0 ≡ 0`) and otherwise returns one of the
        // system's solutions, verified against every row, as the `u64`
        // path always has: `Solved` certifies the rows, and uniqueness
        // is the caller's rank gate and group check.  It is never
        // `Inconsistent`, and a duplicate row is always caught.
        let (mut rows_d, _) = planted::<W>(&mut rng, small, small, 5, &m);
        let combo = rows_d[0].axpy(W::random_nonzero(&mut rng, &m), &rows_d[1], m);
        rows_d[2] = combo; // row 2 now depends on rows 0 and 1
        let dense_d: Vec<Vec<W>> = rows_d
            .iter()
            .map(|r| {
                let mut d = vec![W::zero(); small];
                for &(c, v) in &r.entries {
                    d[c as usize] = v;
                }
                d
            })
            .collect();
        let rhs_d: Vec<W> = rows_d.iter().map(|r| r.rhs).collect();
        assert_eq!(
            dense_solve_unique(dense_d, rhs_d, m),
            None,
            "{s} dense rank-deficient"
        );
        for solver in [
            SparseCoreSolver::BlockWiedemann,
            SparseCoreSolver::BlockLanczos,
        ] {
            let opts = SparseSolveOptions {
                solver,
                ..SparseSolveOptions::default()
            };
            let (outcome, report) = solve_sparse_system(&rows_d, small, m, &opts);
            match outcome {
                SparseSolveOutcome::Undetermined => {}
                SparseSolveOutcome::Solved(sol) => {
                    assert!(
                        rows_d.iter().all(|r| r.evaluate(&sol, m) == r.rhs),
                        "{s} {solver:?} rank-deficient: unverified solution"
                    );
                    assert_eq!(report.filter.dependent_rows_dropped, 0, "{s} {solver:?}");
                }
                SparseSolveOutcome::Inconsistent => {
                    panic!("{s} {solver:?}: consistent rank-deficient system reported inconsistent")
                }
            }
        }
        rows_d[2] = rows_d[1].clone(); // an exact duplicate: removed, rows < columns
        let (outcome, report) =
            solve_sparse_system(&rows_d, small, m, &SparseSolveOptions::default());
        assert_eq!(
            outcome,
            SparseSolveOutcome::Undetermined,
            "{s} duplicate row"
        );
        assert!(
            report.filter.duplicates_removed == 1 || report.filter.uncovered_columns == 1,
            "{s} duplicate row: {report:?}"
        );
        // Inconsistent: one corrupted right-hand side is never "solved".
        let (mut rows_i, _) = planted::<W>(&mut rng, small + 8, small, 5, &m);
        rows_i[5].rhs = rows_i[5].rhs.add(W::one(), &m);
        let (outcome, _) = solve_sparse_system(&rows_i, small, m, &SparseSolveOptions::default());
        assert!(
            matches!(
                outcome,
                SparseSolveOutcome::Inconsistent | SparseSolveOutcome::Undetermined
            ),
            "{s} corrupted: {outcome:?}"
        );
        last.expect("two solvers ran")
    }

    /// The `u64` path through the same checks as the wide widths, on the
    /// 61-bit regression prime.
    #[test]
    fn the_u64_pipeline_passes_the_wide_checks_modulo_a_61_bit_prime() {
        check_pipeline::<u64>(P61, 400, 4);
        check_pipeline::<u64>(P61, 1500, 2);
    }

    #[test]
    fn wide_pipeline_solves_planted_systems_modulo_a_127_bit_prime() {
        check_pipeline::<U128>(P127, 1500, 1);
    }

    #[test]
    fn wide_pipeline_solves_planted_systems_modulo_the_ecc2k130_order() {
        check_pipeline::<U192>(R130, 1500, 2);
    }

    #[test]
    fn wide_pipeline_solves_planted_systems_modulo_190_and_191_bit_primes() {
        check_pipeline::<U192>(P190, 1000, 3);
        check_pipeline::<U192>(P191, 400, 4);
    }

    #[test]
    fn wide_pipeline_solves_planted_systems_modulo_a_250_bit_prime() {
        check_pipeline::<U256>(P250, 800, 5);
    }

    /// The same 61-bit system through the `u64` path and through `U128`:
    /// identical filtering and identical (unique) solutions.  The Krylov
    /// draws differ between widths, so the Krylov reports may.
    #[test]
    fn the_u64_and_wide_paths_agree_modulo_a_61_bit_prime() {
        let p = big(P61);
        let m64 = <u64 as ModWord>::modulus_from_biguint(&p).unwrap();
        let m128 = modulus::<U128>(P61);
        let mut rng = StdRng::seed_from_u64(61);
        let columns = 2000;
        let (rows64, x64) = planted::<u64>(&mut rng, columns + columns / 10 + 64, columns, 5, &m64);
        let rows128: Vec<SparseRow<U128>> = rows64
            .iter()
            .map(|r| SparseRow {
                entries: r
                    .entries
                    .iter()
                    .map(|&(c, v)| (c, U128::from_u64(v, &m128)))
                    .collect(),
                rhs: U128::from_u64(r.rhs, &m128),
            })
            .collect();
        for solver in [
            SparseCoreSolver::BlockWiedemann,
            SparseCoreSolver::BlockLanczos,
        ] {
            let opts = SparseSolveOptions {
                solver,
                ..SparseSolveOptions::default()
            };
            let (o64, r64) = solve_sparse_system(&rows64, columns, m64, &opts);
            let (o128, r128) = solve_sparse_system(&rows128, columns, m128, &opts);
            assert_eq!(
                o64,
                SparseSolveOutcome::Solved(x64.clone()),
                "{solver:?} u64"
            );
            let SparseSolveOutcome::Solved(s128) = o128 else {
                panic!("{solver:?} U128: {o128:?}");
            };
            let s128: Vec<u64> = s128
                .iter()
                .map(|v| v.to_biguint().to_u64_digits().first().copied().unwrap_or(0))
                .collect();
            assert_eq!(s128, x64, "{solver:?} U128");
            assert_eq!(r64.filter, r128.filter, "{solver:?}: filtering differs");
            assert_eq!(r64.core_dimension, r128.core_dimension);
            assert_eq!(r64.core_nonzeros, r128.core_nonzeros);
        }
        // The CSR products agree lane for lane.
        let a64 = CsrMatrix::from_rows(&rows64, columns, m64);
        let a128 = CsrMatrix::from_rows(&rows128, columns, m128);
        let n = 4;
        let v64: Vec<u64> = (0..columns * n).map(|_| rng.gen_range(0..m64)).collect();
        let v128: Vec<U128> = v64.iter().map(|&v| U128::from_u64(v, &m128)).collect();
        let mut y64 = vec![0u64; rows64.len() * n];
        let mut y128 = vec![U128::zero(); rows64.len() * n];
        a64.mul_block(&v64, n, &mut y64);
        a128.mul_block(&v128, n, &mut y128);
        assert!(y64
            .iter()
            .zip(&y128)
            .all(|(&a, b)| U128::from_u64(a, &m128) == *b));
        let mut t64 = vec![0u64; columns * n];
        let mut t128 = vec![U128::zero(); columns * n];
        a64.mul_block_transpose(&y64, n, &mut t64);
        a128.mul_block_transpose(&y128, n, &mut t128);
        assert!(t64
            .iter()
            .zip(&t128)
            .all(|(&a, b)| U128::from_u64(a, &m128) == *b));
    }

    #[test]
    fn sparse_rows_pick_the_width_and_solve_as_big_integers() {
        for (s, width) in [
            (P61, WordWidth::U64),
            (P127, WordWidth::U128),
            (R130, WordWidth::U192),
            (P250, WordWidth::U256),
        ] {
            let p = big(s);
            let mut rows = SparseRows::for_modulus(&p).expect("supported");
            assert_eq!(rows.width(), width);
            assert_eq!(rows.modulus(), p);
            assert!(rows.is_empty());
            let columns = 300;
            let mut rng = StdRng::seed_from_u64(7);
            let x: Vec<BigUint> = (0..columns)
                .map(|_| BigUint::from(rng.gen::<u128>()) % &p)
                .collect();
            for i in 0..columns + 40 {
                let mut entries: Vec<(u32, BigUint)> = vec![];
                if i < columns {
                    entries.push((i as u32, BigUint::from(rng.gen_range(1..1u64 << 40))));
                }
                while entries.len() < 5 {
                    entries.push((
                        rng.gen_range(0..columns) as u32,
                        BigUint::from(rng.gen_range(1..1u64 << 40)),
                    ));
                }
                let rhs = entries
                    .iter()
                    .fold(BigUint::zero(), |acc, (c, v)| acc + v * &x[*c as usize])
                    % &p;
                if i % 2 == 0 {
                    rows.push_entries(&entries, &rhs);
                } else {
                    let mut dense = vec![BigUint::zero(); columns];
                    for (c, v) in &entries {
                        dense[*c as usize] += v;
                    }
                    rows.push_dense(&dense, &rhs);
                }
            }
            assert_eq!(rows.len(), columns + 40);
            let weights = rows.column_weights(columns);
            assert!(weights.iter().all(|&w| w >= 1));
            let (entries, rhs) = rows.row_biguint(0);
            assert!(entries.iter().any(|(c, _)| *c == 0));
            assert!(rhs < p);
            assert_eq!(rows.last_dense(columns).unwrap().len(), columns);
            let (outcome, report) = rows.solve(columns, &SparseSolveOptions::default());
            assert_eq!(outcome, SparseSolveOutcome::Solved(x), "{s}: {report:?}");
            rows.clear();
            assert!(rows.is_empty());
        }
        assert!(SparseRows::for_modulus(&(BigUint::one() << 64u32)).is_none());
    }

    /// Timing of the sparse product and of full solves per width.  Set
    /// `KOBLITZ_SPARSE_LA_COLUMNS` (default 10000),
    /// `KOBLITZ_SPARSE_LA_WIDTH` (`u64`, `u128`, `u192`, `u256`; the
    /// modulus is the fixture prime of that width), and
    /// `KOBLITZ_SPARSE_LA_SOLVERS` (comma list of `wiedemann`, `lanczos`,
    /// `none`); run with `--ignored --nocapture`.  Prints one
    /// `BENCH_JSON {...}` line per measurement.
    #[test]
    #[ignore]
    fn wide_sparse_la_benchmark() {
        let columns: usize = std::env::var("KOBLITZ_SPARSE_LA_COLUMNS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(10_000);
        let width = std::env::var("KOBLITZ_SPARSE_LA_WIDTH").unwrap_or_else(|_| "u64".into());
        let solvers = std::env::var("KOBLITZ_SPARSE_LA_SOLVERS")
            .unwrap_or_else(|_| "wiedemann,lanczos".into());
        match width.as_str() {
            "u64" => bench::<u64>(P61, "u64", columns, &solvers),
            "u128" => bench::<U128>(P127, "u128", columns, &solvers),
            "u192" => bench::<U192>(R130, "u192", columns, &solvers),
            "u256" => bench::<U256>(P250, "u256", columns, &solvers),
            other => panic!("unknown width {other}"),
        }
    }

    /// CPU seconds this process has used over all its threads: on a
    /// loaded shared host wall time overstates the work, so each
    /// measurement reports both.
    fn cpu_seconds() -> f64 {
        #[cfg(unix)]
        // SAFETY: clock_gettime writes into the zeroed timespec.
        unsafe {
            let mut ts: libc::timespec = std::mem::zeroed();
            libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut ts);
            ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
        }
        #[cfg(not(unix))]
        f64::NAN
    }

    /// The host's one-minute load average (`NaN` where unavailable).
    fn load_average() -> f64 {
        #[cfg(unix)]
        // SAFETY: getloadavg writes at most one f64 into the array.
        unsafe {
            let mut avg = [0f64; 1];
            if libc::getloadavg(avg.as_mut_ptr(), 1) == 1 {
                return avg[0];
            }
        }
        f64::NAN
    }

    fn bench<W: ModWord>(s: &str, width: &str, columns: usize, solvers: &str) {
        let m = modulus::<W>(s);
        let bits = W::modulus_to_biguint(&m).bits();
        let threads = rayon::current_num_threads();
        let mut rng = StdRng::seed_from_u64(0xbe4c);
        let rows_n = columns + columns / 10 + 64;
        let (rows, x) = planted::<W>(&mut rng, rows_n, columns, 5, &m);
        let a = CsrMatrix::from_rows(&rows, columns, m);
        let nnz = a.nnz();
        // SpMV: block of 4 lanes (the default Krylov block) and 1 lane,
        // serial and with the default (parallel over rows) backend.
        for lanes in [4usize, 1] {
            let v: Vec<W> = (0..columns * lanes)
                .map(|_| W::random(&mut rng, &m))
                .collect();
            let mut y = vec![W::zero(); rows_n * lanes];
            for backend in [SpmvBackend::Serial, SpmvBackend::Auto] {
                let a = a.clone().with_spmv(SpmvOptions { backend, shards: 0 });
                a.mul_block(&v, lanes, &mut y);
                let reps = (2_000_000_000usize / (nnz * lanes).max(1)).clamp(3, 200);
                let load = load_average();
                let cpu_begin = cpu_seconds();
                let begin = std::time::Instant::now();
                for _ in 0..reps {
                    a.mul_block(&v, lanes, &mut y);
                }
                let secs = begin.elapsed().as_secs_f64() / reps as f64;
                let cpu = (cpu_seconds() - cpu_begin) / reps as f64;
                // Verified against the per-entry reference.
                let mut reference = vec![W::zero(); rows_n * lanes];
                for i in 0..rows_n {
                    a.row_into_per_entry(i, &v, lanes, &mut reference[i * lanes..(i + 1) * lanes]);
                }
                assert_eq!(y, reference, "SpMV mismatch");
                println!(
                    "BENCH_JSON {{\"kind\":\"spmv\",\"width\":\"{width}\",\"modulus_bits\":{bits},\"columns\":{columns},\"rows\":{rows_n},\"nnz\":{nnz},\"lanes\":{lanes},\"backend\":\"{backend:?}\",\"threads\":{threads},\"reps\":{reps},\"seconds_per_product\":{secs:.9},\"ns_per_nnz_lane\":{:.3},\"cpu_seconds_per_product\":{cpu:.9},\"cpu_ns_per_nnz_lane\":{:.3},\"load_average\":{load:.1},\"verified\":true}}",
                    secs * 1e9 / (nnz * lanes) as f64,
                    cpu * 1e9 / (nnz * lanes) as f64
                );
            }
        }
        for solver in solvers.split(',') {
            let solver = match solver.trim() {
                "wiedemann" => SparseCoreSolver::BlockWiedemann,
                "lanczos" => SparseCoreSolver::BlockLanczos,
                _ => continue,
            };
            let opts = SparseSolveOptions {
                solver,
                ..SparseSolveOptions::default()
            };
            let load = load_average();
            let cpu_begin = cpu_seconds();
            let begin = std::time::Instant::now();
            let (outcome, report) = solve_sparse_system(&rows, columns, m, &opts);
            let secs = begin.elapsed().as_secs_f64();
            let cpu = cpu_seconds() - cpu_begin;
            let verified = matches!(&outcome, SparseSolveOutcome::Solved(sol) if *sol == x);
            let core_products = report
                .wiedemann
                .as_ref()
                .map(|w| w.products * w.block_n)
                .or_else(|| {
                    report
                        .lanczos
                        .as_ref()
                        .map(|l| l.sparse_products * l.block_size)
                })
                .unwrap_or(0);
            println!(
                "BENCH_JSON {{\"kind\":\"solve\",\"width\":\"{width}\",\"modulus_bits\":{bits},\"columns\":{columns},\"rows\":{rows_n},\"nnz\":{nnz},\"solver\":\"{solver:?}\",\"threads\":{threads},\"seconds\":{secs:.3},\"cpu_seconds\":{cpu:.3},\"load_average\":{load:.1},\"verified\":{verified},\"core_dimension\":{},\"core_nonzeros\":{},\"attempts\":{},\"core_lane_products\":{core_products},\"report\":{}}}",
                report.core_dimension,
                report.core_nonzeros,
                report.attempts,
                serde_json::to_string(&report).unwrap()
            );
            assert!(verified, "{outcome:?}");
        }
    }
}
