//! PS-3: the bit-sliced group-law kernel on `u64` (64 lanes), AVX2 (256) and
//! AVX-512F (512) vectors, with runtime dispatch and `std::thread` rows.
//!
//! Layout: the rational base points `(x, y)` (sorted by `x`) are laid across the
//! lanes in blocks of `LANES`; a task is a (row `i₁`, block) cell with
//! `block ≥ i₁/LANES`.  The tasks of a chunk of rows are taken in order, `BATCH`
//! at a time, each batch sharing one bit-sliced inversion per chain
//! (Montgomery).  A hit or exceptional lane is re-derived exactly with
//! [`thirds`]: the cold path, which also confirms every hit.

use crate::bitslice::*;
use crate::bspairs::*;
use crate::group::*;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

pub const BATCH: usize = 8;
pub const ROW_CHUNK: usize = 8;

pub struct BsKernel<L: Lane, const N: usize> {
    pub f: BsField<N>,
    gf: Gf2,
    pub l: usize,
    pub base: Vec<(u64, u64)>,
    xb: Vec<E<L>>,
    yb: Vec<E<L>>,
}

/// The planes of a block of up to `LANES` values (lane `j` = value `j`).
fn planes_of<L: Lane, const N: usize>(vals: &[u64]) -> E<L> {
    let mut e = zero_elt::<L>();
    let mut words = vec![0u64; L::WORDS];
    for (i, plane) in e.iter_mut().enumerate().take(N) {
        words.iter_mut().for_each(|w| *w = 0);
        for (j, &v) in vals.iter().enumerate() {
            words[j / 64] |= ((v >> i) & 1) << (j % 64);
        }
        *plane = L::from_words(&words);
    }
    e
}

impl<L: Lane, const N: usize> BsKernel<L, N> {
    pub fn new(gf: &Gf2, solver: &ArtinSchreier, l: u32) -> Self {
        let base = base_points(gf, solver, l);
        let lanes = L::WORDS * 64;
        let blocks = base.len().div_ceil(lanes);
        let mut xb = Vec::with_capacity(blocks);
        let mut yb = Vec::with_capacity(blocks);
        for b in 0..blocks {
            let chunk = &base[b * lanes..((b + 1) * lanes).min(base.len())];
            let xs: Vec<u64> = chunk.iter().map(|p| p.0).collect();
            let ys: Vec<u64> = chunk.iter().map(|p| p.1).collect();
            xb.push(planes_of::<L, N>(&xs));
            yb.push(planes_of::<L, N>(&ys));
        }
        Self {
            f: BsField::<N>::new(gf),
            gf: gf.clone(),
            l: l as usize,
            base,
            xb,
            yb,
        }
    }

    pub fn pairs_total(&self) -> u64 {
        let m = self.base.len() as u64;
        m * (m + 1) / 2
    }

    /// Lanes `i₁ ≤ idx < m` of block `blk` as a vector.
    #[inline(always)]
    fn valid_lanes(&self, i1: usize, blk: usize) -> L {
        let lo = blk * L::WORDS * 64;
        let mut words = vec![0u64; L::WORDS];
        for (w, word) in words.iter_mut().enumerate() {
            for b in 0..64 {
                let idx = lo + w * 64 + b;
                if idx >= i1 && idx < self.base.len() {
                    *word |= 1u64 << b;
                }
            }
        }
        L::from_words(&words)
    }

    /// Process the rows `rows` for target `r`; returns whether a root was found.
    /// Unless `all`, stops at the first confirmed one (and tells the other
    /// threads through `stop`).  `pairs` counts the pairs visited.
    #[inline(always)]
    pub fn run_rows(
        &self,
        r: (u64, u64),
        rows: std::ops::Range<usize>,
        all: bool,
        out: &mut Vec<[u64; 3]>,
        stop: &AtomicBool,
        pairs: &mut u64,
    ) -> bool {
        let lanes = L::WORDS * 64;
        let blocks = self.xb.len();
        let rp = Pt::Aff(r.0, r.1);
        let nr = neg(rp);
        let span = 1u64 << self.l;
        // Row points and planes for the chunk (exact, one inversion per row point).
        let mut row_store: Vec<RowPlanes<L>> = Vec::with_capacity(ROW_CHUNK);
        let mut tasks: Vec<(usize, usize, usize)> = Vec::new(); // (slot, i1, blk)
        for i1 in rows {
            let p1 = Pt::Aff(self.base[i1].0, self.base[i1].1);
            let (Pt::Aff(xu, yu), Pt::Aff(xv, yv)) = (add(&self.gf, p1, rp), add(&self.gf, p1, nr))
            else {
                continue; // P₁ = ±R: x_R ∈ V, excluded by the caller
            };
            let slot = row_store.len();
            row_store.push(RowPlanes::<L>::new(N, xu, yu, xv, yv));
            for blk in (i1 / lanes)..blocks {
                tasks.push((slot, i1, blk));
            }
        }
        let mut found = false;
        let mut words = vec![0u64; L::WORDS];
        for group in tasks.chunks(BATCH) {
            if stop.load(Ordering::Relaxed) {
                break;
            }
            let rs: Vec<&RowPlanes<L>> = group.iter().map(|t| &row_store[t.0]).collect();
            let xs: Vec<&E<L>> = group.iter().map(|t| &self.xb[t.2]).collect();
            let ys: Vec<&E<L>> = group.iter().map(|t| &self.yb[t.2]).collect();
            let res = batch::<L, N>(&self.f, &rs, &xs, &ys, self.l);
            for (k, &(hit, exc)) in res.iter().enumerate() {
                let (_, i1, blk) = group[k];
                let valid = self.valid_lanes(i1, blk);
                valid.to_words(&mut words);
                *pairs += words.iter().map(|w| w.count_ones() as u64).sum::<u64>();
                let hit_v = hit.and(valid).and(exc.not());
                let exc_v = exc.and(valid);
                for (mask, is_hit) in [(hit_v, true), (exc_v, false)] {
                    if mask.is_zero() {
                        continue;
                    }
                    mask.to_words(&mut words);
                    for (w, &word) in words.clone().iter().enumerate() {
                        let mut bits = word;
                        while bits != 0 {
                            let b = bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let i2 = blk * lanes + w * 64 + b;
                            let mut any = false;
                            for t in thirds(&self.gf, self.base[i1], self.base[i2], r)
                                .into_iter()
                                .flatten()
                            {
                                if t < span {
                                    any = true;
                                    let mut tr = [self.base[i1].0, self.base[i2].0, t];
                                    tr.sort_unstable();
                                    out.push(tr);
                                }
                            }
                            // A hit lane the exact path does not confirm is a kernel bug.
                            assert!(!is_hit || any, "bit-sliced hit not confirmed exactly");
                            found |= any;
                        }
                    }
                }
            }
            if found && !all {
                stop.store(true, Ordering::Relaxed);
                break;
            }
        }
        found
    }

    /// All rows with `threads` workers (rows in chunks of `ROW_CHUNK`, dynamic).
    pub fn run(
        &self,
        r: (u64, u64),
        all: bool,
        threads: usize,
        out: &mut Vec<[u64; 3]>,
    ) -> (bool, u64)
    where
        L: Dispatch,
    {
        let m = self.base.len();
        let stop = AtomicBool::new(false);
        let next = AtomicUsize::new(0);
        let worker = |out: &mut Vec<[u64; 3]>| -> (bool, u64) {
            let (mut found, mut pairs) = (false, 0u64);
            loop {
                let s = next.fetch_add(ROW_CHUNK, Ordering::Relaxed);
                if s >= m || stop.load(Ordering::Relaxed) {
                    break;
                }
                found |= L::rows::<N>(
                    self,
                    r,
                    s..(s + ROW_CHUNK).min(m),
                    all,
                    out,
                    &stop,
                    &mut pairs,
                );
            }
            (found, pairs)
        };
        if threads <= 1 {
            return worker(out);
        }
        let results: Vec<(bool, u64, Vec<[u64; 3]>)> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..threads)
                .map(|_| {
                    sc.spawn(|| {
                        let mut o = Vec::new();
                        let (f, p) = worker(&mut o);
                        (f, p, o)
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let (mut found, mut pairs) = (false, 0u64);
        for (f, p, o) in results {
            found |= f;
            pairs += p;
            out.extend(o);
        }
        (found, pairs)
    }
}

/// Dependent chain of `iters` products `a ← a·b` on lane-varying operands; returns a checksum.
#[inline(always)]
fn micro_body<L: Lane, const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
    let mut s = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let mut words = vec![0u64; L::WORDS];
    let mut a = zero_elt::<L>();
    let mut b = zero_elt::<L>();
    for i in 0..N {
        for w in words.iter_mut() {
            *w = next();
        }
        a[i] = L::from_words(&words);
        for w in words.iter_mut() {
            *w = next();
        }
        b[i] = L::from_words(&words);
    }
    for _ in 0..iters {
        a = match mode {
            1 => f.sqr(&a),
            2 => generic_mul::<L, N>(f, &a, &b),
            3 => colwise_mul::<L, N>(f, &a, &b),
            4 => tiled_mul::<L, N, 4>(f, &a, &b),
            5 => tiled_mul::<L, N, 8>(f, &a, &b),
            _ => f.mul(&a, &b),
        };
    }
    if (7..=9).contains(&mode) {
        let row = RowPlanes::<L>::new(
            N,
            0x1234_5678 & ((1u64 << N) - 1),
            0x0F0F_0F0F,
            0x2222_2222,
            0x5555_5555 & ((1u64 << N) - 1),
        );
        let xs: Vec<E<L>> = (0..BATCH)
            .map(|_| {
                let mut e = zero_elt::<L>();
                for i in 0..N {
                    for w in words.iter_mut() {
                        *w = next();
                    }
                    e[i] = L::from_words(&words);
                }
                e
            })
            .collect();
        let xr: Vec<&E<L>> = xs.iter().collect();
        let mut acc = 0u64;
        let sa: Vec<StageA<L>> = xs.iter().map(|x| stage_a(f, &row, x)).collect();
        let dus: Vec<&E<L>> = sa.iter().map(|s| &s.du).collect();
        for _ in 0..iters {
            match mode {
                7 => {
                    let r = montgomery(f, &dus);
                    r[3][0].to_words(&mut words);
                }
                8 => {
                    for k in 0..BATCH {
                        let h = stage_c(
                            f,
                            &row,
                            xr[k],
                            xr[k],
                            &sa[k],
                            xr[k],
                            xr[(k + 1) % BATCH],
                            N / 3,
                        );
                        h.to_words(&mut words);
                    }
                }
                _ => {
                    for k in 0..BATCH {
                        let s = stage_a(f, &row, xr[k]);
                        s.exceptional.to_words(&mut words);
                    }
                }
            }
            acc ^= words[0];
        }
        return acc;
    }
    if mode == 6 {
        // A batch of `BATCH` tasks against one row (all lanes lifted-looking random values).
        let row = RowPlanes::<L>::new(
            N,
            0x1234_5678 & ((1u64 << N) - 1),
            0x0F0F_0F0F,
            0x2222_2222,
            0x5555_5555 & ((1u64 << N) - 1),
        );
        let xs: Vec<E<L>> = (0..BATCH)
            .map(|_| {
                let mut e = zero_elt::<L>();
                for i in 0..N {
                    for w in words.iter_mut() {
                        *w = next();
                    }
                    e[i] = L::from_words(&words);
                }
                e
            })
            .collect();
        let ys = xs.clone();
        let xr: Vec<&E<L>> = xs.iter().collect();
        let yr: Vec<&E<L>> = ys.iter().collect();
        let rr: Vec<&RowPlanes<L>> = (0..BATCH).map(|_| &row).collect();
        let mut acc = 0u64;
        for _ in 0..iters {
            let res = batch::<L, N>(f, &rr, &xr, &yr, N / 3);
            for (h, e) in res {
                h.to_words(&mut words);
                acc ^= words[0];
                e.to_words(&mut words);
                acc ^= words[0];
            }
        }
        return acc;
    }
    a[0].to_words(&mut words);
    words.iter().fold(0, |x, &w| x ^ w)
}

/// Runs the row loop inside a function compiled for the lane's target feature.
pub trait Dispatch: Lane + Send + Sync {
    /// The micro-benchmark body (`bs-micro`), inside the lane's target feature.
    fn micro<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64
    where
        Self: Sized;
    #[allow(clippy::too_many_arguments)]
    fn rows<const N: usize>(
        k: &BsKernel<Self, N>,
        r: (u64, u64),
        rows: std::ops::Range<usize>,
        all: bool,
        out: &mut Vec<[u64; 3]>,
        stop: &AtomicBool,
        pairs: &mut u64,
    ) -> bool
    where
        Self: Sized;
}

impl Dispatch for u64 {
    #[inline(always)]
    fn micro<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
        micro_body::<u64, N>(f, iters, mode)
    }
    #[inline(always)]
    fn rows<const N: usize>(
        k: &BsKernel<u64, N>,
        r: (u64, u64),
        rows: std::ops::Range<usize>,
        all: bool,
        out: &mut Vec<[u64; 3]>,
        stop: &AtomicBool,
        pairs: &mut u64,
    ) -> bool {
        k.run_rows(r, rows, all, out, stop, pairs)
    }
}

#[cfg(target_arch = "x86_64")]
mod x86_dispatch {
    use super::*;
    use crate::lanes_x86::{W256, W512};

    #[target_feature(enable = "avx2")]
    unsafe fn rows_avx2<const N: usize>(
        k: &BsKernel<W256, N>,
        r: (u64, u64),
        rows: std::ops::Range<usize>,
        all: bool,
        out: &mut Vec<[u64; 3]>,
        stop: &AtomicBool,
        pairs: &mut u64,
    ) -> bool {
        k.run_rows(r, rows, all, out, stop, pairs)
    }

    #[target_feature(enable = "avx512f")]
    unsafe fn rows_avx512<const N: usize>(
        k: &BsKernel<W512, N>,
        r: (u64, u64),
        rows: std::ops::Range<usize>,
        all: bool,
        out: &mut Vec<[u64; 3]>,
        stop: &AtomicBool,
        pairs: &mut u64,
    ) -> bool {
        k.run_rows(r, rows, all, out, stop, pairs)
    }

    #[target_feature(enable = "avx2")]
    unsafe fn micro_avx2<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
        micro_body::<W256, N>(f, iters, mode)
    }

    #[target_feature(enable = "avx512f")]
    unsafe fn micro_avx512<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
        micro_body::<W512, N>(f, iters, mode)
    }

    impl Dispatch for W256 {
        fn micro<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
            assert!(std::arch::is_x86_feature_detected!("avx2"));
            // SAFETY: avx2 was just detected.
            unsafe { micro_avx2::<N>(f, iters, mode) }
        }
        fn rows<const N: usize>(
            k: &BsKernel<W256, N>,
            r: (u64, u64),
            rows: std::ops::Range<usize>,
            all: bool,
            out: &mut Vec<[u64; 3]>,
            stop: &AtomicBool,
            pairs: &mut u64,
        ) -> bool {
            assert!(std::arch::is_x86_feature_detected!("avx2"));
            // SAFETY: avx2 was just detected.
            unsafe { rows_avx2::<N>(k, r, rows, all, out, stop, pairs) }
        }
    }

    impl Dispatch for W512 {
        fn micro<const N: usize>(f: &BsField<N>, iters: usize, mode: u8) -> u64 {
            assert!(std::arch::is_x86_feature_detected!("avx512f"));
            // SAFETY: avx512f was just detected.
            unsafe { micro_avx512::<N>(f, iters, mode) }
        }
        fn rows<const N: usize>(
            k: &BsKernel<W512, N>,
            r: (u64, u64),
            rows: std::ops::Range<usize>,
            all: bool,
            out: &mut Vec<[u64; 3]>,
            stop: &AtomicBool,
            pairs: &mut u64,
        ) -> bool {
            assert!(std::arch::is_x86_feature_detected!("avx512f"));
            // SAFETY: avx512f was just detected.
            unsafe { rows_avx512::<N>(k, r, rows, all, out, stop, pairs) }
        }
    }
}
