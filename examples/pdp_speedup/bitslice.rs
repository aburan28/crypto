//! Bit-sliced GF(2^n) arithmetic over generic lanes (PS-2 gate counting and
//! PS-3 execution).  An element is `n` planes; plane `i` holds bit `i` of
//! every lane's element.  The same code runs on `u64` (64 lanes), wider
//! vectors, a counting lane ([`Cnt`]) that tallies gates, and a recording
//! lane ([`Tape`]) that emits the straight-line program for hashing and for
//! the GPU kernel generator.

// The recording lane, the lane accessors and `from_planes_u64` are the interface
// PS-2 (program export and hash) and PS-3 (hit extraction) build on; PS-1 only
// counts gates and runs 64 lanes.
#![allow(dead_code)]

use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use std::cell::{Cell, RefCell};

pub const MAXN: usize = 40;

/// One bit-plane vector.  `xand(acc, a, b) = acc ^ (a & b)` is the fused
/// operation that a ternary-logic instruction (`vpternlogd`, `lop3`) does
/// in one issue.
pub trait Lane: Copy {
    const LANES: usize;
    /// 64-bit words per vector (`LANES / 64`, or 1 for the counting and recording lanes).
    const WORDS: usize;
    /// Load a vector from `WORDS` words, word `w` holding lanes `64w … 64w + 63`.
    fn from_words(w: &[u64]) -> Self;
    /// Store a vector as `WORDS` words.
    fn to_words(self, w: &mut [u64]);
    /// The field product; vector lanes override it with a copy compiled for their
    /// target feature so that the intrinsics inline.
    #[inline(always)]
    fn mul<const N: usize>(f: &BsField<N>, a: &E<Self>, b: &E<Self>) -> E<Self>
    where
        Self: Sized,
    {
        generic_mul_noinline::<Self, N>(f, a, b)
    }
    fn zero() -> Self;
    fn ones() -> Self;
    fn xor(self, o: Self) -> Self;
    fn and(self, o: Self) -> Self;
    fn or(self, o: Self) -> Self;
    fn not(self) -> Self;
    fn xand(acc: Self, a: Self, b: Self) -> Self;
    fn is_zero(self) -> bool;
}

impl Lane for u64 {
    const LANES: usize = 64;
    const WORDS: usize = 1;
    #[inline(always)]
    fn from_words(w: &[u64]) -> Self {
        w[0]
    }
    #[inline(always)]
    fn to_words(self, w: &mut [u64]) {
        w[0] = self;
    }
    #[inline(always)]
    fn zero() -> Self {
        0
    }
    #[inline(always)]
    fn ones() -> Self {
        !0
    }
    #[inline(always)]
    fn xor(self, o: Self) -> Self {
        self ^ o
    }
    #[inline(always)]
    fn and(self, o: Self) -> Self {
        self & o
    }
    #[inline(always)]
    fn or(self, o: Self) -> Self {
        self | o
    }
    #[inline(always)]
    fn not(self) -> Self {
        !self
    }
    #[inline(always)]
    fn xand(acc: Self, a: Self, b: Self) -> Self {
        acc ^ (a & b)
    }
    #[inline(always)]
    fn is_zero(self) -> bool {
        self == 0
    }
    #[inline(always)]
    fn mul<const N: usize>(f: &BsField<N>, a: &E<Self>, b: &E<Self>) -> E<Self> {
        generic_mul_noinline::<u64, N>(f, a, b)
    }
}

// ── gate counting lane ──────────────────────────────────────────────────

/// Gate tallies.  `fused` counts every operation as one issue (an XOR-AND
/// is one ternary-logic issue); `unfused` counts an XOR-AND as two.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Gates {
    pub xor: u64,
    pub and: u64,
    pub or: u64,
    pub not: u64,
    pub xand: u64,
}

impl Gates {
    pub fn fused(&self) -> u64 {
        self.xor + self.and + self.or + self.not + self.xand
    }
    pub fn unfused(&self) -> u64 {
        self.xor + self.and + self.or + self.not + 2 * self.xand
    }
}

thread_local! {
    static GATES: Cell<Gates> = const { Cell::new(Gates { xor: 0, and: 0, or: 0, not: 0, xand: 0 }) };
}

pub fn gates_reset() {
    GATES.with(|g| g.set(Gates::default()));
}

pub fn gates_get() -> Gates {
    GATES.with(|g| g.get())
}

#[inline(always)]
fn bump(f: impl FnOnce(&mut Gates)) {
    GATES.with(|g| {
        let mut v = g.get();
        f(&mut v);
        g.set(v);
    });
}

/// A lane that does no arithmetic and counts the gates asked of it.  Its
/// `is_zero` is `false`, so a data-dependent branch is counted on the
/// not-taken side only: use it on branch-free code.
#[derive(Clone, Copy)]
pub struct Cnt;

impl Lane for Cnt {
    const LANES: usize = 1;
    const WORDS: usize = 1;
    fn from_words(_: &[u64]) -> Self {
        Cnt
    }
    fn to_words(self, w: &mut [u64]) {
        w[0] = 0;
    }
    fn zero() -> Self {
        Cnt
    }
    fn ones() -> Self {
        Cnt
    }
    fn xor(self, _: Self) -> Self {
        bump(|g| g.xor += 1);
        Cnt
    }
    fn and(self, _: Self) -> Self {
        bump(|g| g.and += 1);
        Cnt
    }
    fn or(self, _: Self) -> Self {
        bump(|g| g.or += 1);
        Cnt
    }
    fn not(self) -> Self {
        bump(|g| g.not += 1);
        Cnt
    }
    fn xand(_: Self, _: Self, _: Self) -> Self {
        bump(|g| g.xand += 1);
        Cnt
    }
    fn is_zero(self) -> bool {
        false
    }
}

// ── program-recording lane ──────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TapeOp {
    /// Program input number `k` (a plane supplied by the caller).
    Input(u32),
    Zero,
    Ones,
    Xor(u32, u32),
    And(u32, u32),
    Or(u32, u32),
    Not(u32),
    Xand(u32, u32, u32),
}

thread_local! {
    static TAPE: RefCell<Vec<TapeOp>> = const { RefCell::new(Vec::new()) };
}

pub fn tape_reset() {
    TAPE.with(|t| t.borrow_mut().clear());
}

pub fn tape_take() -> Vec<TapeOp> {
    TAPE.with(|t| std::mem::take(&mut *t.borrow_mut()))
}

/// A wire of the recorded straight-line program.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Tape(pub u32);

fn emit(op: TapeOp) -> Tape {
    TAPE.with(|t| {
        let mut v = t.borrow_mut();
        v.push(op);
        Tape((v.len() - 1) as u32)
    })
}

impl Tape {
    /// A fresh program input.
    pub fn input(k: u32) -> Self {
        emit(TapeOp::Input(k))
    }
}

impl Lane for Tape {
    const LANES: usize = 1;
    const WORDS: usize = 1;
    fn from_words(_: &[u64]) -> Self {
        Tape(0)
    }
    fn to_words(self, w: &mut [u64]) {
        w[0] = 0;
    }
    fn zero() -> Self {
        emit(TapeOp::Zero)
    }
    fn ones() -> Self {
        emit(TapeOp::Ones)
    }
    fn xor(self, o: Self) -> Self {
        emit(TapeOp::Xor(self.0, o.0))
    }
    fn and(self, o: Self) -> Self {
        emit(TapeOp::And(self.0, o.0))
    }
    fn or(self, o: Self) -> Self {
        emit(TapeOp::Or(self.0, o.0))
    }
    fn not(self) -> Self {
        emit(TapeOp::Not(self.0))
    }
    fn xand(acc: Self, a: Self, b: Self) -> Self {
        emit(TapeOp::Xand(acc.0, a.0, b.0))
    }
    fn is_zero(self) -> bool {
        false
    }
}

// ── the field ───────────────────────────────────────────────────────────

pub type E<L> = [L; MAXN];

pub fn zero_elt<L: Lane>() -> E<L> {
    [L::zero(); MAXN]
}

/// `GF(2^N)` in the polynomial basis with `x^N = Σ_{t ∈ tail} x^t`.
pub struct BsField<const N: usize> {
    /// Low terms of the modulus (exponents below `N`), at most four.
    pub tail: [usize; 4],
    pub ntail: usize,
    gf: Gf2,
    /// `sqr_k` maps: `(k, rows)` with `rows[j]` the input planes XORed into output plane `j`.
    sq: Vec<(u32, Vec<Vec<u8>>)>,
}

impl<const N: usize> BsField<N> {
    pub fn new(gf: &Gf2) -> Self {
        assert_eq!(gf.n as usize, N, "field degree differs from the type's");
        assert!(N <= MAXN);
        let mut tail = [0usize; 4];
        let mut ntail = 0;
        for t in 0..N {
            if (gf.irr >> t) & 1 == 1 {
                assert!(ntail < 4, "modulus has more than four low terms");
                tail[ntail] = t;
                ntail += 1;
            }
        }
        let mut f = Self {
            tail,
            ntail,
            gf: gf.clone(),
            sq: Vec::new(),
        };
        // Every `k` an inversion chain and the pipeline use, built once.
        let mut ks = vec![1u32];
        let e = (N - 1) as u32;
        let mut len = 1u32;
        for bit in (0..(31 - e.leading_zeros())).rev() {
            ks.push(len);
            len *= 2;
            if (e >> bit) & 1 == 1 {
                len += 1;
            }
        }
        ks.sort_unstable();
        ks.dedup();
        for k in ks {
            let mut rows: Vec<Vec<u8>> = vec![Vec::new(); N];
            for i in 0..N {
                let col = gf.sqr_k(1u64 << i, k);
                for (j, row) in rows.iter_mut().enumerate() {
                    if (col >> j) & 1 == 1 {
                        row.push(i as u8);
                    }
                }
            }
            f.sq.push((k, rows));
        }
        f
    }

    pub fn scalar(&self) -> &Gf2 {
        &self.gf
    }

    #[inline(always)]
    pub fn add<L: Lane>(&self, a: &E<L>, b: &E<L>) -> E<L> {
        let mut o = zero_elt::<L>();
        for i in 0..N {
            o[i] = a[i].xor(b[i]);
        }
        o
    }

    /// `a + 1`: flips plane 0 (one gate).
    #[inline(always)]
    pub fn add_one<L: Lane>(&self, a: &E<L>) -> E<L> {
        let mut o = *a;
        o[0] = a[0].not();
        o
    }

    /// Product: `n²` fused XOR-AND gates, then the reduction (see [`generic_mul`]).
    #[inline(always)]
    pub fn mul<L: Lane>(&self, a: &E<L>, b: &E<L>) -> E<L> {
        L::mul::<N>(self, a, b)
    }

    /// `a^(2^k)`, a linear map: `popcount − 1` XORs per output plane.
    #[inline(always)]
    pub fn sqr_k<L: Lane>(&self, a: &E<L>, k: u32) -> E<L> {
        let rows = &self
            .sq
            .iter()
            .find(|(kk, _)| *kk == k)
            .unwrap_or_else(|| panic!("sqr_k map for k = {k} was not prepared"))
            .1;
        let mut o = zero_elt::<L>();
        for j in 0..N {
            let ins = &rows[j];
            let mut acc = match ins.first() {
                Some(&i) => a[i as usize],
                None => L::zero(),
            };
            for &i in ins.iter().skip(1) {
                acc = acc.xor(a[i as usize]);
            }
            o[j] = acc;
        }
        o
    }

    #[inline(always)]
    pub fn sqr<L: Lane>(&self, a: &E<L>) -> E<L> {
        self.sqr_k(a, 1)
    }

    /// Inverse by the Itoh–Tsujii chain, as `Gf2::inv_inline` (zero maps to zero).
    #[inline(always)]
    pub fn inv<L: Lane>(&self, a: &E<L>) -> E<L> {
        let e = (N - 1) as u32;
        let mut beta = *a;
        let mut len = 1u32;
        for bit in (0..(31 - e.leading_zeros())).rev() {
            let s = self.sqr_k(&beta, len);
            beta = self.mul(&s, &beta);
            len *= 2;
            if (e >> bit) & 1 == 1 {
                let s = self.sqr(&beta);
                beta = self.mul(&s, a);
                len += 1;
            }
        }
        debug_assert_eq!(len, e);
        self.sqr(&beta)
    }

    /// Lanes whose element is zero (a mask lane); `N − 1` OR gates and a NOT.
    #[inline(always)]
    pub fn zero_mask<L: Lane>(&self, a: &E<L>) -> L {
        let mut acc = a[0];
        for i in 1..N {
            acc = acc.or(a[i]);
        }
        acc.not()
    }

    /// Lanes whose element lies in `V` (all planes from `l` up are zero).
    #[inline(always)]
    pub fn in_subspace_mask<L: Lane>(&self, a: &E<L>, l: usize) -> L {
        let mut acc = a[l];
        for i in l + 1..N {
            acc = acc.or(a[i]);
        }
        acc.not()
    }
}

/// [`generic_mul`] as one function per `(lane, N)`, so that the portable, counting and
/// recording lanes do not inline a 1,089-operation product at every call site.
#[inline(never)]
pub fn generic_mul_noinline<L: Lane, const N: usize>(f: &BsField<N>, a: &E<L>, b: &E<L>) -> E<L> {
    generic_mul::<L, N>(f, a, b)
}

/// The schoolbook product and its reduction, the default [`Lane::mul`]: `n²`
/// operations (`2n − 1` plain ANDs, the rest fused XOR-ANDs) in the
/// output-stationary order `i` outermost, then the fold of the high half.
#[inline(always)]
pub fn generic_mul<L: Lane, const N: usize>(f: &BsField<N>, a: &E<L>, b: &E<L>) -> E<L> {
    let mut c = [L::zero(); 2 * MAXN];
    // Column `k` is first touched at `i = 0` (k < N) or at `j = N − 1` (k ≥ N).
    for i in 0..N {
        for j in 0..N {
            c[i + j] = if i == 0 || j == N - 1 {
                a[i].and(b[j])
            } else {
                L::xand(c[i + j], a[i], b[j])
            };
        }
    }
    for k in (N..2 * N - 1).rev() {
        let hi = c[k];
        for t in 0..f.ntail {
            let p = k - N + f.tail[t];
            c[p] = c[p].xor(hi);
        }
    }
    let mut o = zero_elt::<L>();
    o[..N].copy_from_slice(&c[..N]);
    o
}

/// `c[i + j] ^= a_i · b_j` for `i ∈ is`, `j ∈ js` (edge strips of the tiled product).
#[inline(always)]
fn accumulate_rect<L: Lane>(
    c: &mut [L; 2 * MAXN],
    a: &[L],
    b: &[L],
    is: std::ops::Range<usize>,
    js: std::ops::Range<usize>,
) {
    for i in is {
        let ai = a[i];
        for j in js.clone() {
            c[i + j] = L::xand(c[i + j], ai, b[j]);
        }
    }
}

// The `T × T` tiles with named accumulators, written out gate by gate (an array of accumulators is
// not promoted to registers before the loops unroll, and LLVM then issues a load and a store per
// gate).  `c` is the product array, `o = ib + jb` the tile's first column; row `i` adds
// `a_i · b_j` into column `o + i + j`.
macro_rules! tile8_gates {
    ($L:ty, $a:expr, $b:expr, $c:ident, $ib:expr, $jb:expr) => {{
        let o = $ib + $jb;
        let mut c0 = $c[o + 0];
        let mut c1 = $c[o + 1];
        let mut c2 = $c[o + 2];
        let mut c3 = $c[o + 3];
        let mut c4 = $c[o + 4];
        let mut c5 = $c[o + 5];
        let mut c6 = $c[o + 6];
        let mut c7 = $c[o + 7];
        let mut c8 = $c[o + 8];
        let mut c9 = $c[o + 9];
        let mut c10 = $c[o + 10];
        let mut c11 = $c[o + 11];
        let mut c12 = $c[o + 12];
        let mut c13 = $c[o + 13];
        let mut c14 = $c[o + 14];
        let a0 = $a[$ib + 0];
        let a1 = $a[$ib + 1];
        let a2 = $a[$ib + 2];
        let a3 = $a[$ib + 3];
        let a4 = $a[$ib + 4];
        let a5 = $a[$ib + 5];
        let a6 = $a[$ib + 6];
        let a7 = $a[$ib + 7];
        let b0 = $b[$jb + 0];
        let b1 = $b[$jb + 1];
        let b2 = $b[$jb + 2];
        let b3 = $b[$jb + 3];
        let b4 = $b[$jb + 4];
        let b5 = $b[$jb + 5];
        let b6 = $b[$jb + 6];
        let b7 = $b[$jb + 7];
        c0 = <$L as Lane>::xand(c0, a0, b0);
        c1 = <$L as Lane>::xand(c1, a0, b1);
        c2 = <$L as Lane>::xand(c2, a0, b2);
        c3 = <$L as Lane>::xand(c3, a0, b3);
        c4 = <$L as Lane>::xand(c4, a0, b4);
        c5 = <$L as Lane>::xand(c5, a0, b5);
        c6 = <$L as Lane>::xand(c6, a0, b6);
        c7 = <$L as Lane>::xand(c7, a0, b7);
        c1 = <$L as Lane>::xand(c1, a1, b0);
        c2 = <$L as Lane>::xand(c2, a1, b1);
        c3 = <$L as Lane>::xand(c3, a1, b2);
        c4 = <$L as Lane>::xand(c4, a1, b3);
        c5 = <$L as Lane>::xand(c5, a1, b4);
        c6 = <$L as Lane>::xand(c6, a1, b5);
        c7 = <$L as Lane>::xand(c7, a1, b6);
        c8 = <$L as Lane>::xand(c8, a1, b7);
        c2 = <$L as Lane>::xand(c2, a2, b0);
        c3 = <$L as Lane>::xand(c3, a2, b1);
        c4 = <$L as Lane>::xand(c4, a2, b2);
        c5 = <$L as Lane>::xand(c5, a2, b3);
        c6 = <$L as Lane>::xand(c6, a2, b4);
        c7 = <$L as Lane>::xand(c7, a2, b5);
        c8 = <$L as Lane>::xand(c8, a2, b6);
        c9 = <$L as Lane>::xand(c9, a2, b7);
        c3 = <$L as Lane>::xand(c3, a3, b0);
        c4 = <$L as Lane>::xand(c4, a3, b1);
        c5 = <$L as Lane>::xand(c5, a3, b2);
        c6 = <$L as Lane>::xand(c6, a3, b3);
        c7 = <$L as Lane>::xand(c7, a3, b4);
        c8 = <$L as Lane>::xand(c8, a3, b5);
        c9 = <$L as Lane>::xand(c9, a3, b6);
        c10 = <$L as Lane>::xand(c10, a3, b7);
        c4 = <$L as Lane>::xand(c4, a4, b0);
        c5 = <$L as Lane>::xand(c5, a4, b1);
        c6 = <$L as Lane>::xand(c6, a4, b2);
        c7 = <$L as Lane>::xand(c7, a4, b3);
        c8 = <$L as Lane>::xand(c8, a4, b4);
        c9 = <$L as Lane>::xand(c9, a4, b5);
        c10 = <$L as Lane>::xand(c10, a4, b6);
        c11 = <$L as Lane>::xand(c11, a4, b7);
        c5 = <$L as Lane>::xand(c5, a5, b0);
        c6 = <$L as Lane>::xand(c6, a5, b1);
        c7 = <$L as Lane>::xand(c7, a5, b2);
        c8 = <$L as Lane>::xand(c8, a5, b3);
        c9 = <$L as Lane>::xand(c9, a5, b4);
        c10 = <$L as Lane>::xand(c10, a5, b5);
        c11 = <$L as Lane>::xand(c11, a5, b6);
        c12 = <$L as Lane>::xand(c12, a5, b7);
        c6 = <$L as Lane>::xand(c6, a6, b0);
        c7 = <$L as Lane>::xand(c7, a6, b1);
        c8 = <$L as Lane>::xand(c8, a6, b2);
        c9 = <$L as Lane>::xand(c9, a6, b3);
        c10 = <$L as Lane>::xand(c10, a6, b4);
        c11 = <$L as Lane>::xand(c11, a6, b5);
        c12 = <$L as Lane>::xand(c12, a6, b6);
        c13 = <$L as Lane>::xand(c13, a6, b7);
        c7 = <$L as Lane>::xand(c7, a7, b0);
        c8 = <$L as Lane>::xand(c8, a7, b1);
        c9 = <$L as Lane>::xand(c9, a7, b2);
        c10 = <$L as Lane>::xand(c10, a7, b3);
        c11 = <$L as Lane>::xand(c11, a7, b4);
        c12 = <$L as Lane>::xand(c12, a7, b5);
        c13 = <$L as Lane>::xand(c13, a7, b6);
        c14 = <$L as Lane>::xand(c14, a7, b7);
        $c[o + 0] = c0;
        $c[o + 1] = c1;
        $c[o + 2] = c2;
        $c[o + 3] = c3;
        $c[o + 4] = c4;
        $c[o + 5] = c5;
        $c[o + 6] = c6;
        $c[o + 7] = c7;
        $c[o + 8] = c8;
        $c[o + 9] = c9;
        $c[o + 10] = c10;
        $c[o + 11] = c11;
        $c[o + 12] = c12;
        $c[o + 13] = c13;
        $c[o + 14] = c14;
    }};
}

macro_rules! tile4_gates {
    ($L:ty, $a:expr, $b:expr, $c:ident, $ib:expr, $jb:expr) => {{
        let o = $ib + $jb;
        let mut c0 = $c[o + 0];
        let mut c1 = $c[o + 1];
        let mut c2 = $c[o + 2];
        let mut c3 = $c[o + 3];
        let mut c4 = $c[o + 4];
        let mut c5 = $c[o + 5];
        let mut c6 = $c[o + 6];
        let a0 = $a[$ib + 0];
        let a1 = $a[$ib + 1];
        let a2 = $a[$ib + 2];
        let a3 = $a[$ib + 3];
        let b0 = $b[$jb + 0];
        let b1 = $b[$jb + 1];
        let b2 = $b[$jb + 2];
        let b3 = $b[$jb + 3];
        c0 = <$L as Lane>::xand(c0, a0, b0);
        c1 = <$L as Lane>::xand(c1, a0, b1);
        c2 = <$L as Lane>::xand(c2, a0, b2);
        c3 = <$L as Lane>::xand(c3, a0, b3);
        c1 = <$L as Lane>::xand(c1, a1, b0);
        c2 = <$L as Lane>::xand(c2, a1, b1);
        c3 = <$L as Lane>::xand(c3, a1, b2);
        c4 = <$L as Lane>::xand(c4, a1, b3);
        c2 = <$L as Lane>::xand(c2, a2, b0);
        c3 = <$L as Lane>::xand(c3, a2, b1);
        c4 = <$L as Lane>::xand(c4, a2, b2);
        c5 = <$L as Lane>::xand(c5, a2, b3);
        c3 = <$L as Lane>::xand(c3, a3, b0);
        c4 = <$L as Lane>::xand(c4, a3, b1);
        c5 = <$L as Lane>::xand(c5, a3, b2);
        c6 = <$L as Lane>::xand(c6, a3, b3);
        $c[o + 0] = c0;
        $c[o + 1] = c1;
        $c[o + 2] = c2;
        $c[o + 3] = c3;
        $c[o + 4] = c4;
        $c[o + 5] = c5;
        $c[o + 6] = c6;
    }};
}

/// The product in `T × T` register tiles, `T ∈ {4, 8}`.  The `2T − 1` partial columns of a tile
/// live in registers while its `T²` fused gates run; `T` loads of `a`, `T` of `b` and `2T − 1`
/// loads and stores of `c` serve `T²` gates.  The edge strips (`N mod T`) run untiled.  The
/// gates are those of [`generic_mul`] (`n²` fused, then the fold) in a different order; the
/// recorded program keeps the output-stationary order.
#[inline(always)]
pub fn tiled_mul<L: Lane, const N: usize, const T: usize>(
    f: &BsField<N>,
    a: &E<L>,
    b: &E<L>,
) -> E<L> {
    assert!(T == 4 || T == 8);
    let mut c = [L::zero(); 2 * MAXN];
    let full = N - N % T;
    let mut ib = 0;
    while ib < full {
        let mut jb = 0;
        while jb < full {
            if T == 8 {
                tile8_gates!(L, a, b, c, ib, jb);
            } else {
                tile4_gates!(L, a, b, c, ib, jb);
            }
            jb += T;
        }
        accumulate_rect(&mut c, a, b, ib..ib + T, full..N);
        ib += T;
    }
    accumulate_rect(&mut c, a, b, full..N, 0..N);
    for k in (N..2 * N - 1).rev() {
        let hi = c[k];
        for t in 0..f.ntail {
            let p = k - N + f.tail[t];
            c[p] = c[p].xor(hi);
        }
    }
    let mut o = zero_elt::<L>();
    o[..N].copy_from_slice(&c[..N]);
    o
}

/// [`tiled_mul`] as one function per `(lane, N, T)`.
#[inline(never)]
pub fn tiled_mul_noinline<L: Lane, const N: usize, const T: usize>(
    f: &BsField<N>,
    a: &E<L>,
    b: &E<L>,
) -> E<L> {
    tiled_mul::<L, N, T>(f, a, b)
}

/// [`colwise_mul`] as one function per `(lane, N)`.
#[inline(never)]
pub fn colwise_mul_noinline<L: Lane, const N: usize>(f: &BsField<N>, a: &E<L>, b: &E<L>) -> E<L> {
    colwise_mul::<L, N>(f, a, b)
}

/// The same product scheduled column by column: each output column `k` is one accumulator chain
/// held in a register (`a₀b_k`, then `acc ^= a_i b_{k−i}`), written once; no zero-filled
/// scratch and one store per column instead of one per gate.  The operations and their count
/// (`n²` fused, then the fold) are those of [`generic_mul`], in a different order; the
/// recorded program (whose hash PS-2 froze) keeps the output-stationary order.
#[inline(always)]
pub fn colwise_mul<L: Lane, const N: usize>(f: &BsField<N>, a: &E<L>, b: &E<L>) -> E<L> {
    use std::mem::MaybeUninit;
    let mut c = [MaybeUninit::<L>::uninit(); 2 * MAXN];
    for k in 0..2 * N - 1 {
        let lo = k.saturating_sub(N - 1);
        let hi = if k < N { k } else { N - 1 };
        let mut acc = a[lo].and(b[k - lo]);
        for i in lo + 1..=hi {
            acc = L::xand(acc, a[i], b[k - i]);
        }
        c[k] = MaybeUninit::new(acc);
    }
    for k in (N..2 * N - 1).rev() {
        // SAFETY: every `c[k]`, `k < 2N − 1`, was written in the loop above.
        let hi = unsafe { c[k].assume_init() };
        for t in 0..f.ntail {
            let p = k - N + f.tail[t];
            // SAFETY: as above (`p < 2N − 1`).
            c[p] = MaybeUninit::new(unsafe { c[p].assume_init() }.xor(hi));
        }
    }
    let mut o = zero_elt::<L>();
    for i in 0..N {
        // SAFETY: as above.
        o[i] = unsafe { c[i].assume_init() };
    }
    o
}

// ── transposition for 64 lanes ──────────────────────────────────────────

/// Words (one per lane) to planes.
pub fn to_planes_u64(words: &[u64; 64], n: usize) -> E<u64> {
    let mut o = zero_elt::<u64>();
    for (j, &w) in words.iter().enumerate() {
        for (i, plane) in o.iter_mut().enumerate().take(n) {
            *plane |= ((w >> i) & 1) << j;
        }
    }
    o
}

/// Planes to words (one per lane).
pub fn from_planes_u64(e: &E<u64>, n: usize) -> [u64; 64] {
    let mut w = [0u64; 64];
    for (j, word) in w.iter_mut().enumerate() {
        for (i, plane) in e.iter().enumerate().take(n) {
            *word |= ((plane >> j) & 1) << i;
        }
    }
    w
}

// ── exporting, hashing and interpreting a recorded program ──────────────

/// Keep only the operations an output depends on, renumbering wires.  Returns
/// the live program and the outputs in it.
pub fn eliminate_dead(ops: &[TapeOp], outputs: &[Tape]) -> (Vec<TapeOp>, Vec<u32>) {
    let mut live = vec![false; ops.len()];
    let mut stack: Vec<u32> = outputs.iter().map(|t| t.0).collect();
    while let Some(w) = stack.pop() {
        if live[w as usize] {
            continue;
        }
        live[w as usize] = true;
        match ops[w as usize] {
            TapeOp::Input(_) | TapeOp::Zero | TapeOp::Ones => {}
            TapeOp::Not(a) => stack.push(a),
            TapeOp::Xor(a, b) | TapeOp::And(a, b) | TapeOp::Or(a, b) => {
                stack.push(a);
                stack.push(b);
            }
            TapeOp::Xand(a, b, c) => {
                stack.push(a);
                stack.push(b);
                stack.push(c);
            }
        }
    }
    let mut renum = vec![u32::MAX; ops.len()];
    let mut out = Vec::new();
    for (i, op) in ops.iter().enumerate() {
        if !live[i] {
            continue;
        }
        let r = |w: u32| renum[w as usize];
        let n = match *op {
            TapeOp::Input(k) => TapeOp::Input(k),
            TapeOp::Zero => TapeOp::Zero,
            TapeOp::Ones => TapeOp::Ones,
            TapeOp::Not(a) => TapeOp::Not(r(a)),
            TapeOp::Xor(a, b) => TapeOp::Xor(r(a), r(b)),
            TapeOp::And(a, b) => TapeOp::And(r(a), r(b)),
            TapeOp::Or(a, b) => TapeOp::Or(r(a), r(b)),
            TapeOp::Xand(a, b, c) => TapeOp::Xand(r(a), r(b), r(c)),
        };
        renum[i] = out.len() as u32;
        out.push(n);
    }
    let outs = outputs.iter().map(|t| renum[t.0 as usize]).collect();
    (out, outs)
}

/// Canonical bytes of a program: per operation a tag byte and its operands
/// as `u32` little-endian, then the output wires.
pub fn program_bytes(ops: &[TapeOp], outputs: &[u32]) -> Vec<u8> {
    let mut b = Vec::with_capacity(ops.len() * 13);
    for op in ops {
        match *op {
            TapeOp::Input(k) => {
                b.push(0);
                b.extend_from_slice(&k.to_le_bytes());
            }
            TapeOp::Zero => b.push(1),
            TapeOp::Ones => b.push(2),
            TapeOp::Not(a) => {
                b.push(3);
                b.extend_from_slice(&a.to_le_bytes());
            }
            TapeOp::Xor(x, y) => {
                b.push(4);
                b.extend_from_slice(&x.to_le_bytes());
                b.extend_from_slice(&y.to_le_bytes());
            }
            TapeOp::And(x, y) => {
                b.push(5);
                b.extend_from_slice(&x.to_le_bytes());
                b.extend_from_slice(&y.to_le_bytes());
            }
            TapeOp::Or(x, y) => {
                b.push(6);
                b.extend_from_slice(&x.to_le_bytes());
                b.extend_from_slice(&y.to_le_bytes());
            }
            TapeOp::Xand(x, y, z) => {
                b.push(7);
                b.extend_from_slice(&x.to_le_bytes());
                b.extend_from_slice(&y.to_le_bytes());
                b.extend_from_slice(&z.to_le_bytes());
            }
        }
    }
    b.push(255);
    for o in outputs {
        b.extend_from_slice(&o.to_le_bytes());
    }
    b
}

/// Run a recorded program on `u64` words (64 lanes): inputs by number.
pub fn eval_program(ops: &[TapeOp], outputs: &[u32], inputs: &[u64]) -> Vec<u64> {
    let mut v: Vec<u64> = Vec::with_capacity(ops.len());
    for op in ops {
        let x = match *op {
            TapeOp::Input(k) => inputs[k as usize],
            TapeOp::Zero => 0,
            TapeOp::Ones => !0,
            TapeOp::Not(a) => !v[a as usize],
            TapeOp::Xor(a, b) => v[a as usize] ^ v[b as usize],
            TapeOp::And(a, b) => v[a as usize] & v[b as usize],
            TapeOp::Or(a, b) => v[a as usize] | v[b as usize],
            TapeOp::Xand(a, b, c) => v[a as usize] ^ (v[b as usize] & v[c as usize]),
        };
        v.push(x);
    }
    outputs.iter().map(|&o| v[o as usize]).collect()
}
