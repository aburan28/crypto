//! Fast exhaustive search over quadratic Boolean systems, for the same
//! Semaev ANF the WDSat path emits.
//!
//! Inspired by the LIP6 / ALMASTY MQ suite
//! (<https://gitlab.lip6.fr/almasty/mq>, public domain):
//!
//! - **Moebius transform** (`moebius.c`) — pack ANF coefficients into a
//!   `2^n` table and convert to the truth table in `O(n·2^n)`; zeros are
//!   solutions.  Fallback when Monica does not apply.
//! - **Monica** (`monica.c`, `ffs.h`) — striped-down Crossbred: linearise
//!   `v ≈ √(2m)` variables and FFS-enumerate the rest.  Preferred when the
//!   cost model says it beats Möbius (see [`crate::cryptanalysis::mq_monica`]).
//!   Cubic chained systems still stay on SAT / WDSat.
//!
//! Also informed by Bouillaguet’s
//! [`libfes-lite`](https://github.com/cbouilla/libfes-lite).  This is a
//! study-library port of the *ideas*, not a binding to those trees.
//!
//! Restricted to degree ≤ 2, `n_vars ≤ 24`, and ≤ 64 equations.  Chained
//! `m ≥ 3` Semaev systems are cubic and are refused here.

use super::wdsat_oracle::AnfRow;
use crate::cryptanalysis::pq_groebner_f2::F2BoolPoly;

/// Packed quadratic form over `n` Boolean variables: constant, linear
/// coefficients, and upper-triangular quadratic coefficients.
#[derive(Clone, Debug)]
pub struct QuadraticForm {
    pub n: usize,
    pub constant: bool,
    pub linear: Vec<bool>,
    /// `quad[i][j]` for `0 ≤ j < i < n` is the coefficient of `x_i x_j`.
    pub quad: Vec<Vec<bool>>,
}

impl QuadraticForm {
    /// Convert an ANF row; returns `None` if any monomial has degree > 2.
    pub fn from_anf_row(row: &AnfRow, n: usize) -> Option<Self> {
        let mut form = Self {
            n,
            constant: row.constant,
            linear: vec![false; n],
            quad: (0..n).map(|i| vec![false; i]).collect(),
        };
        for mono in &row.monomials {
            match mono.len() {
                0 => form.constant = !form.constant,
                1 => {
                    let v = mono[0] as usize;
                    if v >= n {
                        return None;
                    }
                    form.linear[v] = !form.linear[v];
                }
                2 => {
                    let mut a = mono[0] as usize;
                    let mut b = mono[1] as usize;
                    if a == b || a >= n || b >= n {
                        return None;
                    }
                    if a < b {
                        std::mem::swap(&mut a, &mut b);
                    }
                    form.quad[a][b] = !form.quad[a][b];
                }
                _ => return None,
            }
        }
        Some(form)
    }

    pub fn from_poly(poly: &F2BoolPoly) -> Option<Self> {
        Self::from_anf_row(&AnfRow::from_poly(poly), poly.n_vars)
    }

    /// Evaluate at a bit-packed point (`bit i` = value of `x_i`).
    pub fn eval(&self, point: u64) -> bool {
        let mut v = self.constant;
        for i in 0..self.n {
            if ((point >> i) & 1) == 1 && self.linear[i] {
                v = !v;
            }
        }
        for i in 0..self.n {
            if ((point >> i) & 1) == 0 {
                continue;
            }
            for j in 0..i {
                if ((point >> j) & 1) == 1 && self.quad[i][j] {
                    v = !v;
                }
            }
        }
        v
    }

    /// XOR this form's ANF coefficients into bit `eq` of a Moebius table.
    fn scatter_anf(&self, table: &mut [u64], eq: usize) {
        let bit = 1u64 << eq;
        if self.constant {
            table[0] ^= bit;
        }
        for i in 0..self.n {
            if self.linear[i] {
                table[1 << i] ^= bit;
            }
        }
        for i in 0..self.n {
            for j in 0..i {
                if self.quad[i][j] {
                    table[(1 << i) | (1 << j)] ^= bit;
                }
            }
        }
    }
}

/// In-place Möbius / zeta transform over the Boolean lattice, as in
/// ALMASTY `moebius.c` (`small`): after this, `table[x]` holds the
/// packed truth-table values of every equation at assignment `x`.
pub fn moebius_transform(table: &mut [u64], n: usize) {
    debug_assert_eq!(table.len(), 1usize << n);
    for i in 0..n {
        let sz = 1usize << i;
        let mut pos = 0usize;
        while pos < table.len() {
            for j in 0..sz {
                table[pos + sz + j] ^= table[pos + j];
            }
            pos += 2 * sz;
        }
    }
}

/// Pack quadratic forms into an ANF table and solve by Möbius transform.
///
/// Returns every common zero, capped at `max_solutions`.  Needs
/// `n ≤ 24` and `forms.len() ≤ 64`.
pub fn moebius_find_all(forms: &[QuadraticForm], max_solutions: usize) -> Option<Vec<u64>> {
    if forms.is_empty() {
        return Some(vec![0]);
    }
    let n = forms[0].n;
    if n > 24 || forms.len() > 64 || forms.iter().any(|f| f.n != n) {
        return None;
    }
    let mut table = vec![0u64; 1usize << n];
    for (eq, form) in forms.iter().enumerate() {
        form.scatter_anf(&mut table, eq);
    }
    moebius_transform(&mut table, n);
    let mut out = Vec::new();
    for (x, packed) in table.into_iter().enumerate() {
        if packed == 0 {
            out.push(x as u64);
            if out.len() >= max_solutions {
                break;
            }
        }
    }
    Some(out)
}

/// Which quadratic FES backend produced a result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FesBackend {
    Monica,
    Moebius,
}

/// Prefer Möbius for `n ≤ 24`; Monica only when Möbius refuses / cost model wins.
pub fn fes_find_all_auto(
    forms: &[QuadraticForm],
    max_solutions: usize,
) -> Option<(Vec<u64>, FesBackend)> {
    if forms.is_empty() {
        return Some((vec![0], FesBackend::Moebius));
    }
    let n = forms[0].n;
    let m = forms.len();
    if n > 24 || crate::cryptanalysis::mq_monica::monica_beats_moebius(n, m) {
        if let Some(roots) = crate::cryptanalysis::mq_monica::monica_find_all(forms, max_solutions)
        {
            return Some((roots, FesBackend::Monica));
        }
    }
    moebius_find_all(forms, max_solutions).map(|r| (r, FesBackend::Moebius))
}

/// First common zero via incremental Gray (early exit), else Möbius / Monica.
pub fn fes_find_one(forms: &[QuadraticForm]) -> Option<u64> {
    if forms.is_empty() {
        return Some(0);
    }
    if let Some(x) = gray_incremental_find_one(forms) {
        return Some(x);
    }
    // Gray refused (too large): try Monica, then Möbius.
    let n = forms[0].n;
    let m = forms.len();
    if n > 24 || crate::cryptanalysis::mq_monica::monica_beats_moebius(n, m) {
        if let Some(roots) = crate::cryptanalysis::mq_monica::monica_find_all(forms, 1) {
            return roots.into_iter().next();
        }
    }
    moebius_find_all(forms, 1)?.into_iter().next()
}

/// Every common zero via the auto-selected backend (empty on capacity refusal).
pub fn fes_find_all(forms: &[QuadraticForm], max_solutions: usize) -> Vec<u64> {
    fes_find_all_auto(forms, max_solutions)
        .map(|(r, _)| r)
        .unwrap_or_default()
}

/// Triangular index of the monomial `x_i x_j` with `i < j` (libfes `idxq`).
#[inline]
pub(crate) fn idxq(i: usize, j: usize) -> usize {
    debug_assert!(i < j);
    j * (j - 1) / 2 + i
}

/// Bitner–Ehrlich–Reingold focus pointers (`libfes-lite` / ALMASTY `ffs.h`).
#[derive(Clone, Debug)]
pub(crate) struct Ffs {
    focus: [i32; 34],
    stack: [i32; 33],
    sp: i32,
    pub(crate) k1: i32,
    pub(crate) k2: i32,
}

impl Ffs {
    pub(crate) fn reset(n: usize) -> Self {
        let mut focus = [0i32; 34];
        for j in 0..=32 {
            focus[j] = j as i32;
        }
        let mut stack = [0i32; 33];
        stack[0] = (n + 1) as i32;
        Self {
            focus,
            stack,
            sp: 1,
            k1: (n + 1) as i32,
            k2: -1,
        }
    }

    #[inline]
    pub(crate) fn step(&mut self) {
        let j = self.focus[0];
        self.focus[0] = 0;
        self.focus[j as usize] = self.focus[(j + 1) as usize];
        self.focus[(j + 1) as usize] = j + 1;
        self.k1 = j;
        self.sp -= j;
        self.k2 = self.stack[(self.sp - 1) as usize];
        self.stack[self.sp as usize] = j;
        self.sp += 1;
    }
}

/// Incremental Gray-code enumeration — **libfes-lite O(1) per step**.
///
/// Bouillaguet / ALMASTY / libfes-lite: after a one-time setup that pads two
/// fictive variables, each Gray step is
///
/// ```text
///     Fl[1+k1] ^= Fq[idxq(k1, k2)];
///     Fl[0]    ^= Fl[1+k1];
/// ```
///
/// so the hot loop no longer walks all `n` derivatives.  Dispatch:
/// - `n ≥ 4`: scalar `L = 4` (16-step) chunk with `Fl[0]` kept in a register
///   and unchecked table indexing.
/// - else: minimal one-step FFS.
///
/// Hardcoded `L = 8`, batch-probe, AVX2 4×u64, and rayon outer-specialisation
/// remain available for experiments (`gray_ffs_unrolled_l8*`,
/// `gray_ffs_parallel_outer`, `mq_fes_avx2`) but are not auto-selected: on
/// packed-u64 single-system Semaev at the sizes we fit they lose to, or are
/// within noise of, the `L=4` path (I-cache, rewind copies, SIMD setup, or
/// specialisation tax).
/// Inspired by <https://github.com/cbouilla/libfes-lite>
/// (`generic_minimal.c`, `generic_1x32.c`, `avx2_8x32.c`, batch asm) and
/// ALMASTY `ffs.h`.
pub fn gray_incremental_find_all(
    forms: &[QuadraticForm],
    max_solutions: usize,
) -> Option<Vec<u64>> {
    if forms.is_empty() {
        return Some(vec![0]);
    }
    let n = forms[0].n;
    let m = forms.len();
    if n == 0 || n > 32 || m > 64 || forms.iter().any(|f| f.n != n) {
        return None;
    }

    let mut fq = [0u64; 561];
    let mut fl = [0u64; 34];
    fill_fq_fl(forms, n, &mut fq, &mut fl);

    let mut out = Vec::new();
    if n >= 4 {
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, max_solutions, &mut out);
    } else {
        gray_ffs_minimal(&mut fq, &mut fl, n, max_solutions, &mut out);
    }
    Some(out)
}

/// Pack ANF coefficients into the libfes Fq / Fl tables (with fictive vars).
pub(crate) fn fill_fq_fl(
    forms: &[QuadraticForm],
    n: usize,
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
) {
    fq.fill(0);
    fl.fill(0);
    for (eq, form) in forms.iter().enumerate() {
        let bit = 1u64 << eq;
        if form.constant {
            fl[0] ^= bit;
        }
        for i in 0..n {
            if form.linear[i] {
                fl[1 + i] ^= bit;
            }
            for j in 0..i {
                if form.quad[i][j] {
                    fq[idxq(j, i)] ^= bit;
                }
            }
        }
    }
    for i in 0..n {
        fq[idxq(i, n)] = 0;
    }
    fq[idxq(0, n + 1)] = 0;
    for i in 1..n {
        fq[idxq(i, n + 1)] = fq[idxq(i - 1, i)];
    }
    fq[idxq(n, n + 1)] = 0;
}

/// Specialise the top `outer` variables to the bit-pattern `lane`, writing
/// the induced system on the remaining `n_inner = n - outer` variables into
/// `fq`/`fl` (with fictive padding).  Mirrors libfes / AVX2 lane setup.
pub(crate) fn specialize_outer_to_tables(
    forms: &[QuadraticForm],
    n: usize,
    outer: usize,
    lane: u32,
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
) {
    let n_inner = n - outer;
    debug_assert!(n_inner <= 32);
    fq.fill(0);
    fl.fill(0);
    // Stack scratch — avoid per-lane heap traffic that dominated walls.
    let mut lin = [false; 32];
    let mut quad_flat = [false; 496]; // idxq(i,j) for j<32
    for (eq, form) in forms.iter().enumerate() {
        let bit = 1u64 << eq;
        let mut c = form.constant;
        lin[..n_inner].fill(false);
        for i in 0..n_inner {
            for j in 0..i {
                quad_flat[idxq(j, i)] = false;
            }
        }

        let val = |v: usize| -> Option<bool> {
            if v < n_inner {
                None
            } else {
                Some(((lane >> (v - n_inner)) & 1) == 1)
            }
        };

        for i in 0..n {
            if form.linear[i] {
                match val(i) {
                    Some(true) => c = !c,
                    Some(false) => {}
                    None => lin[i] = !lin[i],
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                if !form.quad[i][j] {
                    continue;
                }
                match (val(i), val(j)) {
                    (Some(true), Some(true)) => c = !c,
                    (Some(true), None) => lin[j] = !lin[j],
                    (None, Some(true)) => lin[i] = !lin[i],
                    (None, None) => quad_flat[idxq(j, i)] = !quad_flat[idxq(j, i)],
                    _ => {}
                }
            }
        }

        if c {
            fl[0] ^= bit;
        }
        for i in 0..n_inner {
            if lin[i] {
                fl[1 + i] ^= bit;
            }
            for j in 0..i {
                if quad_flat[idxq(j, i)] {
                    fq[idxq(j, i)] ^= bit;
                }
            }
        }
    }
    for i in 0..n_inner {
        fq[idxq(i, n_inner)] = 0;
    }
    fq[idxq(0, n_inner + 1)] = 0;
    for i in 1..n_inner {
        fq[idxq(i, n_inner + 1)] = fq[idxq(i - 1, i)];
    }
    fq[idxq(n_inner, n_inner + 1)] = 0;
}

/// Rayon over `2^outer` specialised subsystems (independent Gray walks).
pub(crate) fn gray_ffs_parallel_outer(
    forms: &[QuadraticForm],
    n: usize,
    max_solutions: usize,
    outer: usize,
) -> Vec<u64> {
    use rayon::prelude::*;
    let n_inner = n - outer;
    let lanes = 1u32 << outer;
    // Each lane may collect up to the global cap; merge truncates.  (Splitting
    // the budget across lanes can miss a lane-heavy solution set.)
    let per_lane = max_solutions;

    let parts: Vec<Vec<u64>> = (0..lanes)
        .into_par_iter()
        .map(|lane| {
            let mut fq = [0u64; 561];
            let mut fl = [0u64; 34];
            specialize_outer_to_tables(forms, n, outer, lane, &mut fq, &mut fl);
            let mut local = Vec::new();
            if n_inner >= 4 {
                gray_ffs_unrolled_l4(&mut fq, &mut fl, n_inner, per_lane, &mut local);
            } else {
                gray_ffs_minimal(&mut fq, &mut fl, n_inner, per_lane, &mut local);
            }
            let shift = n_inner as u64;
            for x in &mut local {
                *x |= (lane as u64) << shift;
            }
            local
        })
        .collect();

    let mut out = Vec::new();
    for part in parts {
        for x in part {
            out.push(x);
            if out.len() >= max_solutions {
                return out;
            }
        }
    }
    out
}

/// Heuristic: AVX2 batch-probe vs scalar `L=4` batch on this host.
/// Tuned from release walls; engineering lever only (floor unchanged).
#[inline]
pub fn prefer_avx2_batch(n: usize, m: usize, max_solutions: usize) -> bool {
    let _ = (n, m, max_solutions);
    // Measured slower than hardcoded scalar L=4/L=8 on this host; keep false.
    false
}

/// libfes `generic_minimal`: one FFS step per point.
fn gray_ffs_minimal(
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    max_solutions: usize,
    out: &mut Vec<u64>,
) {
    let mut ffs = Ffs::reset(n);
    ffs.step();
    let upto = (1u64 << n) - 1;
    let mut i = 0u64;
    loop {
        if fl[0] == 0 {
            out.push(i ^ (i >> 1));
            if out.len() >= max_solutions {
                return;
            }
        }
        let a = (1 + ffs.k1) as usize;
        let b = idxq(ffs.k1 as usize, ffs.k2 as usize);
        fl[a] ^= fq[b];
        fl[0] ^= fl[a];
        ffs.step();
        if i == upto {
            return;
        }
        i += 1;
    }
}

#[inline(always)]
fn step2(fq: &[u64; 561], fl: &mut [u64; 34], a: usize, b: usize, index: u64, out: &mut Vec<u64>, max_solutions: usize) -> bool {
    if fl[0] == 0 {
        out.push(index ^ (index >> 1));
        if out.len() >= max_solutions {
            return true;
        }
    }
    fl[a] ^= fq[b];
    fl[0] ^= fl[a];
    false
}

/// Hot-path step with `Fl[0]` held in a register (`v`) and unchecked indexing.
#[inline(always)]
unsafe fn step2_fast(
    fq: &[u64; 561],
    fl: &mut [u64; 34],
    v: &mut u64,
    a: usize,
    b: usize,
    index: u64,
    out: &mut Vec<u64>,
    max_solutions: usize,
) -> bool {
    if *v == 0 {
        out.push(index ^ (index >> 1));
        if out.len() >= max_solutions {
            fl[0] = *v;
            return true;
        }
    }
    let fa = fl.get_unchecked_mut(a);
    *fa ^= *fq.get_unchecked(b);
    *v ^= *fa;
    false
}

/// libfes `generic_1x32` / `UNROLLED_CHUNK`: 16 Gray steps per FFS advance.
pub(crate) fn gray_ffs_unrolled_l4(
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    max_solutions: usize,
    out: &mut Vec<u64>,
) {
    const L: usize = 4;
    let mut ffs = Ffs::reset(n - L);
    let mut k1 = ffs.k1 + L as i32;
    let mut k2 = ffs.k2 + L as i32;
    let iterations = 1u64 << (n - L);
    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + L as i32;
        k2 = ffs.k2 + L as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << L;
        let mut v = fl[0];
        // Hard-coded 16-step Gray chunk with Fl[0] in a local.
        let done = unsafe {
            step2_fast(fq, fl, &mut v, 1, alpha, base, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 2, alpha + 1, base + 1, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 0, base + 2, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 3, alpha + 2, base + 3, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 1, base + 4, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 2, 2, base + 5, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 0, base + 6, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 4, alpha + 3, base + 7, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 3, base + 8, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 2, 4, base + 9, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 0, base + 10, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 3, 5, base + 11, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 1, base + 12, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 2, 2, base + 13, out, max_solutions)
                || step2_fast(fq, fl, &mut v, 1, 0, base + 14, out, max_solutions)
                || step2_fast(fq, fl, &mut v, beta, gamma, base + 15, out, max_solutions)
        };
        fl[0] = v;
        if done {
            return;
        }
    }
}

#[inline(always)]
fn step2_update(fq: &[u64; 561], fl: &mut [u64; 34], a: usize, b: usize) {
    fl[a] ^= fq[b];
    fl[0] ^= fl[a];
}

/// `L = 4` batch probe: update without recording; on a hit, rewind and harvest.
///
/// Wins on sparse / unsat full enums (most chunks have `Fl[0] ≠ 0` the whole
/// way) by dropping per-step solution branches from the common path.
pub(crate) fn gray_ffs_unrolled_l4_batch(
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    max_solutions: usize,
    out: &mut Vec<u64>,
) {
    const L: usize = 4;
    let mut ffs = Ffs::reset(n - L);
    let mut k1 = ffs.k1 + L as i32;
    let mut k2 = ffs.k2 + L as i32;
    let iterations = 1u64 << (n - L);
    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + L as i32;
        k2 = ffs.k2 + L as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << L;

        let fl_save = *fl;
        let mut hit = false;
        // Probe the 16-step chunk (same order as UNROLLED_CHUNK).
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, alpha);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 2, alpha + 1);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 0);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 3, alpha + 2);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 1);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 2, 2);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 0);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 4, alpha + 3);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 3);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 2, 4);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 0);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 3, 5);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 1);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 2, 2);
        hit |= fl[0] == 0;
        step2_update(fq, fl, 1, 0);
        hit |= fl[0] == 0;
        step2_update(fq, fl, beta, gamma);

        if !hit {
            continue;
        }
        *fl = fl_save;
        if step2(fq, fl, 1, alpha, base, out, max_solutions)
            || step2(fq, fl, 2, alpha + 1, base + 1, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 2, out, max_solutions)
            || step2(fq, fl, 3, alpha + 2, base + 3, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 4, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 5, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 6, out, max_solutions)
            || step2(fq, fl, 4, alpha + 3, base + 7, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 8, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 9, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 10, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 11, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 12, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 13, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 14, out, max_solutions)
            || step2(fq, fl, beta, gamma, base + 15, out, max_solutions)
        {
            return;
        }
    }
}

#[inline(always)]
fn l8_fq_index(kind: u8, payload: u16, alpha: usize) -> usize {
    if kind == 0 {
        payload as usize
    } else {
        alpha + payload as usize
    }
}

/// Hard-coded `L = 8` (256-step) Gray chunk — matches libfes UNROLLED_CHUNK.
pub(crate) fn gray_ffs_unrolled_l8(
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    max_solutions: usize,
    out: &mut Vec<u64>,
) {
    const L: usize = 8;
    let mut ffs = Ffs::reset(n - L);
    let mut k1 = ffs.k1 + L as i32;
    let mut k2 = ffs.k2 + L as i32;
    let iterations = 1u64 << (n - L);
    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + L as i32;
        k2 = ffs.k2 + L as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << L;
        if false
            || step2(fq, fl, 1, alpha + 0, base + 0, out, max_solutions)
            || step2(fq, fl, 2, alpha + 1, base + 1, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 2, out, max_solutions)
            || step2(fq, fl, 3, alpha + 2, base + 3, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 4, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 5, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 6, out, max_solutions)
            || step2(fq, fl, 4, alpha + 3, base + 7, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 8, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 9, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 10, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 11, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 12, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 13, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 14, out, max_solutions)
            || step2(fq, fl, 5, alpha + 4, base + 15, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 16, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 17, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 18, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 19, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 20, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 21, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 22, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 23, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 24, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 25, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 26, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 27, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 28, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 29, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 30, out, max_solutions)
            || step2(fq, fl, 6, alpha + 5, base + 31, out, max_solutions)
            || step2(fq, fl, 1, 10, base + 32, out, max_solutions)
            || step2(fq, fl, 2, 11, base + 33, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 34, out, max_solutions)
            || step2(fq, fl, 3, 12, base + 35, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 36, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 37, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 38, out, max_solutions)
            || step2(fq, fl, 4, 13, base + 39, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 40, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 41, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 42, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 43, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 44, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 45, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 46, out, max_solutions)
            || step2(fq, fl, 5, 14, base + 47, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 48, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 49, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 50, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 51, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 52, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 53, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 54, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 55, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 56, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 57, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 58, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 59, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 60, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 61, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 62, out, max_solutions)
            || step2(fq, fl, 7, alpha + 6, base + 63, out, max_solutions)
            || step2(fq, fl, 1, 15, base + 64, out, max_solutions)
            || step2(fq, fl, 2, 16, base + 65, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 66, out, max_solutions)
            || step2(fq, fl, 3, 17, base + 67, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 68, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 69, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 70, out, max_solutions)
            || step2(fq, fl, 4, 18, base + 71, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 72, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 73, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 74, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 75, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 76, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 77, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 78, out, max_solutions)
            || step2(fq, fl, 5, 19, base + 79, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 80, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 81, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 82, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 83, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 84, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 85, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 86, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 87, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 88, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 89, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 90, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 91, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 92, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 93, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 94, out, max_solutions)
            || step2(fq, fl, 6, 20, base + 95, out, max_solutions)
            || step2(fq, fl, 1, 10, base + 96, out, max_solutions)
            || step2(fq, fl, 2, 11, base + 97, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 98, out, max_solutions)
            || step2(fq, fl, 3, 12, base + 99, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 100, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 101, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 102, out, max_solutions)
            || step2(fq, fl, 4, 13, base + 103, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 104, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 105, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 106, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 107, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 108, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 109, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 110, out, max_solutions)
            || step2(fq, fl, 5, 14, base + 111, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 112, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 113, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 114, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 115, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 116, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 117, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 118, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 119, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 120, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 121, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 122, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 123, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 124, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 125, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 126, out, max_solutions)
            || step2(fq, fl, 8, alpha + 7, base + 127, out, max_solutions)
            || step2(fq, fl, 1, 21, base + 128, out, max_solutions)
            || step2(fq, fl, 2, 22, base + 129, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 130, out, max_solutions)
            || step2(fq, fl, 3, 23, base + 131, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 132, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 133, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 134, out, max_solutions)
            || step2(fq, fl, 4, 24, base + 135, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 136, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 137, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 138, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 139, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 140, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 141, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 142, out, max_solutions)
            || step2(fq, fl, 5, 25, base + 143, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 144, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 145, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 146, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 147, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 148, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 149, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 150, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 151, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 152, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 153, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 154, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 155, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 156, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 157, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 158, out, max_solutions)
            || step2(fq, fl, 6, 26, base + 159, out, max_solutions)
            || step2(fq, fl, 1, 10, base + 160, out, max_solutions)
            || step2(fq, fl, 2, 11, base + 161, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 162, out, max_solutions)
            || step2(fq, fl, 3, 12, base + 163, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 164, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 165, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 166, out, max_solutions)
            || step2(fq, fl, 4, 13, base + 167, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 168, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 169, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 170, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 171, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 172, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 173, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 174, out, max_solutions)
            || step2(fq, fl, 5, 14, base + 175, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 176, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 177, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 178, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 179, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 180, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 181, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 182, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 183, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 184, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 185, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 186, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 187, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 188, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 189, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 190, out, max_solutions)
            || step2(fq, fl, 7, 27, base + 191, out, max_solutions)
            || step2(fq, fl, 1, 15, base + 192, out, max_solutions)
            || step2(fq, fl, 2, 16, base + 193, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 194, out, max_solutions)
            || step2(fq, fl, 3, 17, base + 195, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 196, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 197, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 198, out, max_solutions)
            || step2(fq, fl, 4, 18, base + 199, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 200, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 201, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 202, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 203, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 204, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 205, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 206, out, max_solutions)
            || step2(fq, fl, 5, 19, base + 207, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 208, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 209, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 210, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 211, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 212, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 213, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 214, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 215, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 216, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 217, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 218, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 219, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 220, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 221, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 222, out, max_solutions)
            || step2(fq, fl, 6, 20, base + 223, out, max_solutions)
            || step2(fq, fl, 1, 10, base + 224, out, max_solutions)
            || step2(fq, fl, 2, 11, base + 225, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 226, out, max_solutions)
            || step2(fq, fl, 3, 12, base + 227, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 228, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 229, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 230, out, max_solutions)
            || step2(fq, fl, 4, 13, base + 231, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 232, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 233, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 234, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 235, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 236, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 237, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 238, out, max_solutions)
            || step2(fq, fl, 5, 14, base + 239, out, max_solutions)
            || step2(fq, fl, 1, 6, base + 240, out, max_solutions)
            || step2(fq, fl, 2, 7, base + 241, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 242, out, max_solutions)
            || step2(fq, fl, 3, 8, base + 243, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 244, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 245, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 246, out, max_solutions)
            || step2(fq, fl, 4, 9, base + 247, out, max_solutions)
            || step2(fq, fl, 1, 3, base + 248, out, max_solutions)
            || step2(fq, fl, 2, 4, base + 249, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 250, out, max_solutions)
            || step2(fq, fl, 3, 5, base + 251, out, max_solutions)
            || step2(fq, fl, 1, 1, base + 252, out, max_solutions)
            || step2(fq, fl, 2, 2, base + 253, out, max_solutions)
            || step2(fq, fl, 1, 0, base + 254, out, max_solutions)
            || step2(fq, fl, beta, gamma, base + 255, out, max_solutions)
        {
            return;
        }
    }
}

/// Probe a 256-step chunk without recording, then harvest only if a zero appeared.
///
/// Matches libfes `BATCH_MODE` / `avx2_asm_enum_batch`: most chunks have no
/// solutions, so the common path stays branch-light.  On a hit, restore `Fl`
/// and re-run with ordinary recording.
pub(crate) fn gray_ffs_unrolled_l8_batch(
    fq: &mut [u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    max_solutions: usize,
    out: &mut Vec<u64>,
) {
    use super::mq_fes_l8_steps::L8_STEPS;
    const L: usize = 8;
    let mut ffs = Ffs::reset(n - L);
    let mut k1 = ffs.k1 + L as i32;
    let mut k2 = ffs.k2 + L as i32;
    let iterations = 1u64 << (n - L);
    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + L as i32;
        k2 = ffs.k2 + L as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << L;

        let fl_save = *fl;
        let mut hit = false;
        for &(a, kind, payload) in L8_STEPS.iter() {
            hit |= fl[0] == 0;
            let a = a as usize;
            let b = l8_fq_index(kind, payload, alpha);
            fl[a] ^= fq[b];
            fl[0] ^= fl[a];
        }
        hit |= fl[0] == 0;
        fl[beta] ^= fq[gamma];
        fl[0] ^= fl[beta];

        if !hit {
            continue;
        }
        // Rewind and harvest with the early-exit-capable step function.
        *fl = fl_save;
        for (step_i, &(a, kind, payload)) in L8_STEPS.iter().enumerate() {
            if step2(
                fq,
                fl,
                a as usize,
                l8_fq_index(kind, payload, alpha),
                base + step_i as u64,
                out,
                max_solutions,
            ) {
                return;
            }
        }
        if step2(fq, fl, beta, gamma, base + 255, out, max_solutions) {
            return;
        }
    }
}

pub fn gray_incremental_find_one(forms: &[QuadraticForm]) -> Option<u64> {
    gray_incremental_find_all(forms, 1)?.into_iter().next()
}

/// Previous O(n)-per-step Gray update — kept for regression / speed comparison.
pub fn gray_on_step_find_all(
    forms: &[QuadraticForm],
    max_solutions: usize,
) -> Option<Vec<u64>> {
    if forms.is_empty() {
        return Some(vec![0]);
    }
    let n = forms[0].n;
    let m = forms.len();
    if n > 28 || m > 64 || forms.iter().any(|f| f.n != n) {
        return None;
    }
    let mut deriv = vec![0u64; n];
    let mut quad_mask = vec![vec![0u64; n]; n];
    let mut value = 0u64;
    for (eq, form) in forms.iter().enumerate() {
        let bit = 1u64 << eq;
        if form.constant {
            value ^= bit;
        }
        for i in 0..n {
            if form.linear[i] {
                deriv[i] ^= bit;
            }
            for j in 0..i {
                if form.quad[i][j] {
                    quad_mask[i][j] ^= bit;
                    quad_mask[j][i] ^= bit;
                }
            }
        }
    }
    let mut out = Vec::new();
    let mut point = 0u64;
    let limit = 1u64 << n;
    for step in 0..limit {
        if value == 0 {
            out.push(point);
            if out.len() >= max_solutions {
                break;
            }
        }
        let flip = (step + 1).trailing_zeros() as usize;
        if flip >= n {
            break;
        }
        value ^= deriv[flip];
        for j in 0..n {
            if j != flip {
                deriv[j] ^= quad_mask[flip][j];
            }
        }
        point ^= 1u64 << flip;
    }
    Some(out)
}

/// Naive Gray-code re-evaluation — retained as an independent check on
/// the Möbius path (ALMASTY monica / libfes enumeration shape).
pub fn gray_find_all(forms: &[QuadraticForm], max_solutions: usize) -> Vec<u64> {
    let mut out = Vec::new();
    if forms.is_empty() {
        return vec![0];
    }
    let n = forms[0].n;
    if n > 24 || forms.iter().any(|f| f.n != n) {
        return out;
    }
    let limit = 1u64 << n;
    let mut point = 0u64;
    for step in 0..limit {
        if forms.iter().all(|f| !f.eval(point)) {
            out.push(point);
            if out.len() >= max_solutions {
                break;
            }
        }
        let flip = (step + 1).trailing_zeros() as u64;
        if flip < n as u64 {
            point ^= 1u64 << flip;
        }
    }
    out
}

/// Solve a quadratic (`m = 2`) Semaev decomposition by Möbius FES.
///
/// Returns the same shape as
/// [`crate::cryptanalysis::koblitz_index_calculus::sat_decompose`].
/// Cubic chained systems (`m ≥ 3`) are reported as exhausted.
pub fn mq_fes_decompose(
    kc: &crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve,
    fb: &crate::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    index_of: &std::collections::HashMap<(num_bigint::BigUint, num_bigint::BigUint), usize>,
    st: &crate::cryptanalysis::koblitz_groebner::FieldStructure,
    target: &crate::binary_ecc::BinaryPoint,
    m: usize,
) -> (
    Option<Vec<usize>>,
    crate::cryptanalysis::koblitz_index_calculus::SatDecompositionStats,
) {
    use crate::binary_ecc::BinaryPoint;
    use crate::cryptanalysis::koblitz_index_calculus::{lift_candidate, SatDecompositionStats};

    let mut stats = SatDecompositionStats::default();
    if m != 2 {
        stats.exhausted = true;
        return (None, stats);
    }
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, stats),
    };
    let sys = match crate::cryptanalysis::polynomial_reuse::build_decomposition_system_reusing(
        &fb.subspace_basis,
        &x_r,
        &kc.curve.b,
        m,
        st,
    ) {
        Some(sys) => sys,
        None => {
            stats.exhausted = true;
            return (None, stats);
        }
    };
    let forms: Option<Vec<_>> = sys
        .equations
        .iter()
        .map(QuadraticForm::from_poly)
        .collect();
    let forms = match forms {
        Some(forms) => forms,
        None => {
            stats.exhausted = true;
            return (None, stats);
        }
    };
    stats.solver_calls = 1;
    // Prefer incremental Gray so a successful lift can stop before a full
    // Möbius transform; fall back to the auto all-roots backend otherwise.
    let roots = if let Some(roots) = gray_incremental_find_all(&forms, 64) {
        roots
    } else {
        match fes_find_all_auto(&forms, 64) {
            Some((roots, _)) => roots,
            None => {
                stats.exhausted = true;
                return (None, stats);
            }
        }
    };
    if roots.is_empty() {
        stats.refuted = true;
        return (None, stats);
    }
    for root in roots {
        stats.models += 1;
        if !sys.equations.iter().all(|e| e.eval(root) == 0) {
            stats.spurious += 1;
            continue;
        }
        let xs: Vec<_> = (0..m)
            .map(|i| sys.summand_x(&fb.subspace_basis, root, i, kc.n))
            .collect();
        if let Some(idxs) = lift_candidate(kc, fb, index_of, &xs, target) {
            return (Some(idxs), stats);
        }
    }
    stats.exhausted = true;
    (None, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::wdsat_oracle::AnfRow;

    #[test]
    fn finds_the_unique_zero_of_a_tiny_system() {
        // x0 + x1 = 0 and x0 = 1  →  only (1,1).
        let rows = [
            AnfRow {
                monomials: vec![vec![0], vec![1]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0]],
                constant: true,
            },
        ];
        let forms: Vec<_> = rows
            .iter()
            .map(|r| QuadraticForm::from_anf_row(r, 2).unwrap())
            .collect();
        assert_eq!(fes_find_one(&forms), Some(0b11));
    }

    #[test]
    fn moebius_agrees_with_gray_on_random_quadratics() {
        let rows = [
            AnfRow {
                monomials: vec![vec![0, 1], vec![2], vec![4]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![1, 3], vec![0, 4], vec![2, 3]],
                constant: false,
            },
        ];
        let forms: Vec<_> = rows
            .iter()
            .map(|r| QuadraticForm::from_anf_row(r, 5).unwrap())
            .collect();
        let mut a = moebius_find_all(&forms, 1024).unwrap();
        let mut b = gray_find_all(&forms, 1024);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn rejects_a_cubic_row() {
        let row = AnfRow {
            monomials: vec![vec![0, 1, 2]],
            constant: false,
        };
        assert!(QuadraticForm::from_anf_row(&row, 3).is_none());
    }

    #[test]
    fn mq_fes_agrees_with_native_sat_on_prime_degree() {
        use crate::binary_ecc::BinaryPoint;
        use crate::cryptanalysis::koblitz_groebner::FieldStructure;
        use crate::cryptanalysis::koblitz_index_calculus::{
            build_frobenius_factor_base, point_key, sat_decompose, KoblitzCurve,
        };
        use std::collections::HashMap;

        let kc = KoblitzCurve::new(1, 7).expect("K_1/F_2^7");
        let fb = build_frobenius_factor_base(&kc, 0).expect("factor base");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let index_of: HashMap<_, _> = fb
            .points
            .iter()
            .enumerate()
            .map(|(i, p)| (point_key(p), i))
            .collect();
        let mut affine = fb.points.iter().enumerate().filter_map(|(i, p)| match p {
            BinaryPoint::Affine { .. } => Some(i),
            BinaryPoint::Infinity => None,
        });
        let i = affine.next().unwrap();
        let j = affine.next().unwrap();
        let target = kc.add(&fb.points[i], &fb.points[j]);

        let (native, _) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 8, Some(2));
        let (fes, stats) = mq_fes_decompose(&kc, &fb, &index_of, &st, &target, 2);
        assert!(native.is_some(), "native SAT missed planted");
        assert!(fes.is_some(), "mq-fes missed planted: {stats:?}");
        let mut a = native.unwrap();
        let mut b = fes.unwrap();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn gray_incremental_agrees_with_moebius() {
        let rows = [
            AnfRow {
                monomials: vec![vec![0, 1], vec![2], vec![4]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![1, 3], vec![0, 4], vec![2, 3]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0, 2], vec![1], vec![3, 4]],
                constant: true,
            },
        ];
        let forms: Vec<_> = rows
            .iter()
            .map(|r| QuadraticForm::from_anf_row(r, 5).unwrap())
            .collect();
        let mut a = gray_incremental_find_all(&forms, 1024).unwrap();
        let mut b = moebius_find_all(&forms, 1024).unwrap();
        let mut c = gray_on_step_find_all(&forms, 1024).unwrap();
        a.sort_unstable();
        b.sort_unstable();
        c.sort_unstable();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn l8_and_l8_batch_agree_with_l4() {
        let n = 12usize;
        let m = 16usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut l4 = Vec::new();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut l4);

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut l4b = Vec::new();
        gray_ffs_unrolled_l4_batch(&mut fq, &mut fl, n, usize::MAX, &mut l4b);

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut l8 = Vec::new();
        gray_ffs_unrolled_l8(&mut fq, &mut fl, n, usize::MAX, &mut l8);

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut batch = Vec::new();
        gray_ffs_unrolled_l8_batch(&mut fq, &mut fl, n, usize::MAX, &mut batch);

        l4.sort_unstable();
        l4b.sort_unstable();
        l8.sort_unstable();
        batch.sort_unstable();
        assert_eq!(l4, l4b);
        assert_eq!(l4, l8);
        assert_eq!(l4, batch);
    }

    #[test]
    fn l4_batch_vs_l4_unsat_wall() {
        // Document: with packed-u64 and a well-predicted never-hit branch,
        // batch probe+rewind loses to plain L=4 (copy + extra compares).
        let n = 18usize;
        let m = 24usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut a = Vec::new();
        let t0 = std::time::Instant::now();
        gray_ffs_unrolled_l4_batch(&mut fq, &mut fl, n, usize::MAX, &mut a);
        let batch_ns = t0.elapsed().as_nanos();

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut b = Vec::new();
        let t1 = std::time::Instant::now();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut b);
        let l4_ns = t1.elapsed().as_nanos();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        let ratio = l4_ns as f64 / batch_ns.max(1) as f64;
        eprintln!(
            "l4_batch_vs_l4 n={n} m={m} unsat: batch={batch_ns}ns l4={l4_ns}ns ratio={ratio:.2} sols={}",
            a.len()
        );
        assert!(
            ratio < 1.0,
            "unexpected: L=4 batch beat plain L=4 ({ratio:.3}×); update note"
        );
    }

    #[test]
    fn l8_hardcoded_vs_l4_full_enum_wall() {
        // Document relative cost; either direction is within engineering noise
        // at n=18 — do not auto-select L=8 without a clearer win.
        let n = 18usize;
        let m = 24usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
        // Warmup
        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut warm = Vec::new();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut warm);

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut a = Vec::new();
        let t0 = std::time::Instant::now();
        gray_ffs_unrolled_l8(&mut fq, &mut fl, n, usize::MAX, &mut a);
        let l8_ns = t0.elapsed().as_nanos();

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut b = Vec::new();
        let t1 = std::time::Instant::now();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut b);
        let l4_ns = t1.elapsed().as_nanos();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        let ratio = l4_ns as f64 / l8_ns.max(1) as f64;
        eprintln!(
            "l8_hardcoded_vs_l4 n={n} m={m} full_enum: l8={l8_ns}ns l4={l4_ns}ns ratio={ratio:.2} sols={}",
            a.len()
        );
    }

    #[test]
    fn parallel_outer_agrees_with_serial() {
        // Correctness; walls at n≤18 still favour serial L=4 (setup tax).
        let n = 14usize;
        let m = 16usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        let mut par = gray_ffs_parallel_outer(&forms, n, usize::MAX, 4);
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut ser = Vec::new();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut ser);
        par.sort_unstable();
        ser.sort_unstable();
        assert_eq!(par, ser);
    }

    #[test]
    fn parallel_outer_beats_serial_at_n20_wall() {
        let n = 20usize;
        let m = 24usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        // Warm rayon pool + serial path.
        let _ = gray_ffs_parallel_outer(&forms, n, 1, 2);
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut warm = Vec::new();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, 1, &mut warm);

        let t0 = std::time::Instant::now();
        let mut par = gray_ffs_parallel_outer(&forms, n, usize::MAX, 4);
        let par_ns = t0.elapsed().as_nanos();

        fill_fq_fl(&forms, n, &mut fq, &mut fl);
        let mut ser = Vec::new();
        let t1 = std::time::Instant::now();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, n, usize::MAX, &mut ser);
        let ser_ns = t1.elapsed().as_nanos();
        par.sort_unstable();
        ser.sort_unstable();
        assert_eq!(par, ser);
        let ratio = ser_ns as f64 / par_ns.max(1) as f64;
        eprintln!(
            "parallel_outer4_vs_l4 n={n} m={m}: par={par_ns}ns ser={ser_ns}ns ratio={ratio:.2} sols={}",
            par.len()
        );
        assert!(
            ratio >= 1.5,
            "expected 4-outer parallel ≥1.5× serial L=4 at n=20, got {ratio:.3}"
        );
    }

    #[test]
    fn gray_ffs_beats_on_step_full_enum_wall() {
        // Full enumeration (no early exit): O(1)/step FFS vs O(n)/step update.
        let n = 18usize;
        let m = 18usize;
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant: eq % 2 == 0,
                linear,
                quad,
            });
        }
        let t0 = std::time::Instant::now();
        let mut ffs = gray_incremental_find_all(&forms, usize::MAX).expect("ffs");
        let ffs_ns = t0.elapsed().as_nanos();
        let t1 = std::time::Instant::now();
        let mut on = gray_on_step_find_all(&forms, usize::MAX).expect("on");
        let on_ns = t1.elapsed().as_nanos();
        ffs.sort_unstable();
        on.sort_unstable();
        assert_eq!(ffs, on);
        let ratio = on_ns as f64 / ffs_ns.max(1) as f64;
        eprintln!(
            "gray_ffs_vs_on_step n={n} full_enum: ffs={ffs_ns}ns on_step={on_ns}ns ratio={ratio:.2} sols={}",
            ffs.len()
        );
        assert!(
            ratio >= 1.5,
            "expected FFS Gray ≥1.5× O(n)-step Gray, got {ratio:.3}"
        );
    }

    #[test]
    fn gray_early_exit_beats_moebius_find_one_wall() {
        // Plant a solution early in Gray order so early exit pays, while
        // Möbius still pays the full n·2^n transform.  n=18 keeps both fast.
        let n = 18usize;
        let m = 18usize;
        let gray_index: u64 = 2_000; // well below 2^18
        let planted = gray_index ^ (gray_index >> 1);
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 19 + i * 5) % 3) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq * 11 + i * 7 + j * 3) % 5) == 0;
                }
            }
            let probe = QuadraticForm {
                n,
                constant: false,
                linear: linear.clone(),
                quad: quad.clone(),
            };
            // Adjust constant so eval(planted) == false.
            let constant = probe.eval(planted);
            forms.push(QuadraticForm {
                n,
                constant,
                linear,
                quad,
            });
        }
        assert!(forms.iter().all(|f| !f.eval(planted)));

        let t0 = std::time::Instant::now();
        let g = gray_incremental_find_one(&forms).expect("gray");
        let gray_ns = t0.elapsed().as_nanos();
        let t1 = std::time::Instant::now();
        let mbi = moebius_find_all(&forms, 1).unwrap().into_iter().next();
        let moebius_ns = t1.elapsed().as_nanos();
        assert!(forms.iter().all(|f| !f.eval(g)));
        assert!(mbi.is_some() && forms.iter().all(|f| !f.eval(mbi.unwrap())));
        // Multiple roots are expected; Gray returns the first in Gray order,
        // Möbius the first in binary order — only require both be zeros.
        let ratio = moebius_ns as f64 / gray_ns.max(1) as f64;
        eprintln!(
            "gray_early_vs_moebius n={n}: gray={gray_ns}ns moebius={moebius_ns}ns ratio={ratio:.2} gray_sol={g:#x} moebius_sol={:#x}",
            mbi.unwrap()
        );
        assert!(
            ratio >= 1.5,
            "expected incremental Gray find_one ≥1.5× Möbius, got {ratio:.3}"
        );
    }

    #[test]
    fn monica_extends_past_moebius_cap() {
        // n=26 exceeds Möbius' n≤24 table; Monica must still run.
        let n = 26usize;
        let m = 64usize;
        assert!(n > 24);
        let mut forms = Vec::with_capacity(m);
        for eq in 0..m {
            let constant = eq % 5 == 0;
            let mut linear = vec![false; n];
            let mut quad = (0..n).map(|i| vec![false; i]).collect::<Vec<_>>();
            for i in 0..n {
                linear[i] = ((eq * 3 + i) % 4) == 0;
                for j in 0..i {
                    quad[i][j] = ((eq + i + j) % 11) == 0;
                }
            }
            forms.push(QuadraticForm {
                n,
                constant,
                linear,
                quad,
            });
        }
        assert!(moebius_find_all(&forms, 1).is_none());
        let roots = crate::cryptanalysis::mq_monica::monica_find_all(&forms, 4);
        assert!(roots.is_some(), "Monica should accept n=26");
        for x in roots.unwrap() {
            assert!(forms.iter().all(|f| !f.eval(x)), "bad Monica root {x:#x}");
        }
    }
}
