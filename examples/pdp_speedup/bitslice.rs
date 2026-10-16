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
use std::collections::HashMap;

pub const MAXN: usize = 40;

/// One bit-plane vector.  `xand(acc, a, b) = acc ^ (a & b)` is the fused
/// operation that a ternary-logic instruction (`vpternlogd`, `lop3`) does
/// in one issue.
pub trait Lane: Copy {
    const LANES: usize;
    fn zero() -> Self;
    fn ones() -> Self;
    fn xor(self, o: Self) -> Self;
    fn and(self, o: Self) -> Self;
    fn or(self, o: Self) -> Self;
    fn not(self) -> Self;
    fn xand(acc: Self, a: Self, b: Self) -> Self;
    fn is_zero(self) -> bool;
    /// Lane `j`'s bit.
    fn bit(self, j: usize) -> bool;
    fn set_bit(self, j: usize) -> Self;
}

impl Lane for u64 {
    const LANES: usize = 64;
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
    fn bit(self, j: usize) -> bool {
        (self >> j) & 1 == 1
    }
    #[inline(always)]
    fn set_bit(self, j: usize) -> Self {
        self | (1u64 << j)
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
    fn bit(self, _: usize) -> bool {
        false
    }
    fn set_bit(self, _: usize) -> Self {
        Cnt
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
    fn bit(self, _: usize) -> bool {
        false
    }
    fn set_bit(self, _: usize) -> Self {
        self
    }
}

// ── the field ───────────────────────────────────────────────────────────

pub type E<L> = [L; MAXN];

pub fn zero_elt<L: Lane>() -> E<L> {
    [L::zero(); MAXN]
}

/// `GF(2^n)` in the polynomial basis with `x^n = Σ_{t ∈ tail} x^t`.
pub struct BsField {
    pub n: usize,
    pub tail: Vec<usize>,
    gf: Gf2,
    /// `sqr_k` matrices: row `j` is the mask of input planes XORed into
    /// output plane `j`.
    sq: RefCell<HashMap<u32, Vec<u64>>>,
}

impl BsField {
    pub fn new(gf: &Gf2) -> Self {
        let n = gf.n as usize;
        assert!(n <= MAXN);
        let tail: Vec<usize> = (0..n).filter(|&t| (gf.irr >> t) & 1 == 1).collect();
        Self {
            n,
            tail,
            gf: gf.clone(),
            sq: RefCell::new(HashMap::new()),
        }
    }

    pub fn scalar(&self) -> &Gf2 {
        &self.gf
    }

    #[inline(always)]
    pub fn add<L: Lane>(&self, a: &E<L>, b: &E<L>) -> E<L> {
        let mut o = zero_elt::<L>();
        for i in 0..self.n {
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

    /// Schoolbook product, `n²` fused XOR-AND gates, then the reduction.
    pub fn mul<L: Lane>(&self, a: &E<L>, b: &E<L>) -> E<L> {
        let n = self.n;
        let mut c = [L::zero(); 2 * MAXN];
        // The first term of each column is a plain AND, the rest are fused
        // XOR-ANDs: column `k` is first touched at `i = 0` (k < n) or at
        // `j = n − 1` (k ≥ n), looping `i` outermost.
        for i in 0..n {
            for j in 0..n {
                c[i + j] = if i == 0 || j == n - 1 {
                    a[i].and(b[j])
                } else {
                    L::xand(c[i + j], a[i], b[j])
                };
            }
        }
        self.reduce_wide(&mut c)
    }

    #[inline(always)]
    fn reduce_wide<L: Lane>(&self, c: &mut [L; 2 * MAXN]) -> E<L> {
        let n = self.n;
        for k in (n..2 * n - 1).rev() {
            let hi = c[k];
            for &t in &self.tail {
                c[k - n + t] = c[k - n + t].xor(hi);
            }
        }
        let mut o = zero_elt::<L>();
        o[..n].copy_from_slice(&c[..n]);
        o
    }

    fn sq_masks(&self, k: u32) -> std::cell::Ref<'_, Vec<u64>> {
        if !self.sq.borrow().contains_key(&k) {
            let n = self.n;
            let mut rows = vec![0u64; n];
            for i in 0..n {
                let col = self.gf.sqr_k(1u64 << i, k);
                for (j, row) in rows.iter_mut().enumerate() {
                    *row |= ((col >> j) & 1) << i;
                }
            }
            self.sq.borrow_mut().insert(k, rows);
        }
        std::cell::Ref::map(self.sq.borrow(), |m| &m[&k])
    }

    /// `a^(2^k)`, a linear map: `popcount − 1` XORs per output plane.
    pub fn sqr_k<L: Lane>(&self, a: &E<L>, k: u32) -> E<L> {
        let masks = self.sq_masks(k);
        let mut o = zero_elt::<L>();
        for j in 0..self.n {
            let mut m = masks[j];
            let mut acc: Option<L> = None;
            while m != 0 {
                let i = m.trailing_zeros() as usize;
                m &= m - 1;
                acc = Some(match acc {
                    None => a[i],
                    Some(v) => v.xor(a[i]),
                });
            }
            o[j] = acc.unwrap_or_else(L::zero);
        }
        o
    }

    pub fn sqr<L: Lane>(&self, a: &E<L>) -> E<L> {
        self.sqr_k(a, 1)
    }

    /// Inverse by the Itoh–Tsujii chain, as `Gf2::inv_inline` (zero maps to zero).
    pub fn inv<L: Lane>(&self, a: &E<L>) -> E<L> {
        let e = (self.n - 1) as u32;
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

    /// Lanes whose element is zero (a mask lane); `n − 1` OR gates and a NOT.
    pub fn zero_mask<L: Lane>(&self, a: &E<L>) -> L {
        let mut acc = a[0];
        for i in 1..self.n {
            acc = acc.or(a[i]);
        }
        acc.not()
    }

    /// Lanes whose element lies in `V` (all planes from `l` up are zero).
    pub fn in_subspace_mask<L: Lane>(&self, a: &E<L>, l: usize) -> L {
        let mut acc = a[l];
        for i in l + 1..self.n {
            acc = acc.or(a[i]);
        }
        acc.not()
    }
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
