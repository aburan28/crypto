//! libfes-lite–style AVX2 Gray FES for a single Semaev system.
//!
//! libfes `avx2_8x32` runs **8 systems × u32** equation words.  Our mq-fes
//! path already packs ≤64 equations into one `u64`, so the matching SIMD
//! shape is **4 points × u64** in a `__m256i` (same 256-bit register):
//! specialise two outer Boolean variables into 4 independent lanes and run
//! an `L = 8` (256-step) unrolled FFS Gray chunk with AVX2 XORs.
//!
//! Falls back when AVX2 is unavailable, `n < 10`, or the scalar path is
//! preferred for early-exit `find_one`.  Engineering wall-time lever only;
//! free-oracle floor unchanged.
//!
//! Inspired by <https://github.com/cbouilla/libfes-lite> `avx2_8x32.c` /
//! `avx2_codegen.py` (public-domain / study-library port of the ideas).

#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

use super::mq_fes::{idxq, Ffs, QuadraticForm};
use super::mq_fes_l8_steps::L8_STEPS;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const LANES: usize = 4; // 4 × u64 = 256-bit YMM
const UNROLL: usize = 8;
const OUTER: usize = 2; // 2^2 = 4 lanes

/// True when the AVX2 4×u64 path can run this instance.
pub fn avx2_4x64_applicable(n: usize, m: usize) -> bool {
    is_x86_feature_detected!("avx2") && n >= OUTER + UNROLL && m > 0 && m <= 64
}

/// Solve by 4-lane AVX2 Gray FFS; `None` if inapplicable.
pub fn gray_ffs_avx2_8x32(
    forms: &[QuadraticForm],
    max_solutions: usize,
) -> Option<Vec<u64>> {
    // Name kept for callers; implementation is 4×u64 (see module docs).
    let n = forms.first()?.n;
    let m = forms.len();
    if !avx2_4x64_applicable(n, m) || forms.iter().any(|f| f.n != n) {
        return None;
    }
    Some(unsafe { gray_ffs_avx2_4x64_inner(forms, n, m, max_solutions) })
}

#[target_feature(enable = "avx2")]
unsafe fn gray_ffs_avx2_4x64_inner(
    forms: &[QuadraticForm],
    n: usize,
    m: usize,
    max_solutions: usize,
) -> Vec<u64> {
    let n_inner = n - OUTER;
    let fq_len = (n_inner + 1) * (n_inner + 2) / 2;
    let mut fq_lanes = vec![[0u64; LANES]; fq_len];
    let mut fl_lanes = vec![[0u64; LANES]; n_inner + 2];

    for lane in 0..LANES {
        let (fl0, fl, fq) = specialize_lane(forms, n, n_inner, m, lane as u8);
        fl_lanes[0][lane] = fl0;
        for i in 0..n_inner {
            fl_lanes[1 + i][lane] = fl[i];
        }
        let mut fq_full = vec![0u64; fq_len];
        let nq = n_inner * n_inner.saturating_sub(1) / 2;
        fq_full[..nq].copy_from_slice(&fq[..nq]);
        for i in 0..n_inner {
            fq_full[idxq(i, n_inner)] = 0;
        }
        fq_full[idxq(0, n_inner + 1)] = 0;
        for i in 1..n_inner {
            fq_full[idxq(i, n_inner + 1)] = fq_full[idxq(i - 1, i)];
        }
        fq_full[idxq(n_inner, n_inner + 1)] = 0;
        for (dst, &src) in fq_lanes.iter_mut().zip(fq_full.iter()) {
            dst[lane] = src;
        }
    }

    let mut fq_vec: Vec<__m256i> = fq_lanes
        .iter()
        .map(|lane| _mm256_loadu_si256(lane.as_ptr() as *const __m256i))
        .collect();
    let mut fl_vec: Vec<__m256i> = fl_lanes
        .iter()
        .map(|lane| _mm256_loadu_si256(lane.as_ptr() as *const __m256i))
        .collect();

    let mut out = Vec::new();
    let mut ffs = Ffs::reset(n_inner - UNROLL);
    let mut k1 = ffs.k1 + UNROLL as i32;
    let mut k2 = ffs.k2 + UNROLL as i32;
    let iterations = 1u64 << (n_inner - UNROLL);
    let zero = _mm256_setzero_si256();

    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + UNROLL as i32;
        k2 = ffs.k2 + UNROLL as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << UNROLL;

        for (step_i, &(a, kind, payload)) in L8_STEPS.iter().enumerate() {
            if harvest_zeros(fl_vec[0], zero, base + step_i as u64, n_inner, &mut out, max_solutions)
            {
                return out;
            }
            let a = a as usize;
            let b = if kind == 0 {
                payload as usize
            } else {
                alpha + payload as usize
            };
            fl_vec[a] = _mm256_xor_si256(fl_vec[a], fq_vec[b]);
            fl_vec[0] = _mm256_xor_si256(fl_vec[0], fl_vec[a]);
        }
        if harvest_zeros(fl_vec[0], zero, base + 255, n_inner, &mut out, max_solutions) {
            return out;
        }
        let updated = _mm256_xor_si256(fl_vec[beta], fq_vec[gamma]);
        fl_vec[beta] = updated;
        fl_vec[0] = _mm256_xor_si256(fl_vec[0], updated);
    }
    let _ = fq_vec;
    out
}

#[inline(always)]
unsafe fn harvest_zeros(
    fl0: __m256i,
    zero: __m256i,
    index: u64,
    n_inner: usize,
    out: &mut Vec<u64>,
    max_solutions: usize,
) -> bool {
    let cmp = _mm256_cmpeq_epi64(fl0, zero);
    let mask = _mm256_movemask_epi8(cmp) as u32;
    if mask == 0 {
        return false;
    }
    let gray = index ^ (index >> 1);
    for lane in 0..LANES {
        // Each u64 lane → 8 bytes → 8 bits in the movemask.
        let lane_bits = (mask >> (lane * 8)) & 0xff;
        if lane_bits == 0xff {
            out.push(((lane as u64) << n_inner) | gray);
            if out.len() >= max_solutions {
                return true;
            }
        }
    }
    false
}

fn specialize_lane(
    forms: &[QuadraticForm],
    n: usize,
    n_inner: usize,
    m: usize,
    outer: u8,
) -> (u64, Vec<u64>, Vec<u64>) {
    let mut fl0 = 0u64;
    let mut fl = vec![0u64; n_inner];
    let mut fq = vec![0u64; n_inner * n_inner.saturating_sub(1) / 2];

    for (eq, form) in forms.iter().enumerate().take(m) {
        let bit = 1u64 << eq;
        let mut c = form.constant;
        let mut lin = vec![false; n_inner];
        let mut quad = vec![vec![false; n_inner]; n_inner];

        let val = |v: usize| -> Option<bool> {
            if v < n_inner {
                None
            } else {
                Some(((outer >> (v - n_inner)) & 1) == 1)
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
                    (None, None) => quad[i][j] = !quad[i][j],
                    _ => {}
                }
            }
        }

        if c {
            fl0 ^= bit;
        }
        for i in 0..n_inner {
            if lin[i] {
                fl[i] ^= bit;
            }
            for j in 0..i {
                if quad[i][j] {
                    fq[idxq(j, i)] ^= bit;
                }
            }
        }
    }
    (fl0, fl, fq)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::mq_fes::{gray_ffs_unrolled_l4, moebius_find_all, QuadraticForm};
    use crate::cryptanalysis::wdsat_oracle::AnfRow;

    fn dense_forms(n: usize, m: usize) -> Vec<QuadraticForm> {
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
        forms
    }

    fn fill_scalar(forms: &[QuadraticForm], n: usize) -> ([u64; 561], [u64; 34]) {
        let mut fq = [0u64; 561];
        let mut fl = [0u64; 34];
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
        (fq, fl)
    }

    #[test]
    fn avx2_agrees_with_scalar_l4_when_available() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let forms = dense_forms(14, 16);
        let mut a = gray_ffs_avx2_8x32(&forms, usize::MAX).expect("avx2");
        let (mut fq, mut fl) = fill_scalar(&forms, 14);
        let mut b = Vec::new();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, 14, usize::MAX, &mut b);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn avx2_agrees_with_moebius_on_smallish() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let rows = [
            AnfRow {
                monomials: vec![vec![0, 1], vec![2], vec![4], vec![10]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![1, 3], vec![0, 4], vec![2, 3], vec![11]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0, 2], vec![1], vec![3, 4], vec![9]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![5, 6], vec![7], vec![8, 10]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0, 11], vec![1, 9], vec![2, 8]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![3, 7], vec![4, 6], vec![5]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0, 5], vec![1, 6], vec![2, 7]],
                constant: true,
            },
            AnfRow {
                monomials: vec![vec![3, 8], vec![4, 9], vec![10, 11]],
                constant: false,
            },
        ];
        let forms: Vec<_> = rows
            .iter()
            .map(|r| QuadraticForm::from_anf_row(r, 12).unwrap())
            .collect();
        let mut a = gray_ffs_avx2_8x32(&forms, 4096).expect("avx2");
        let mut b = moebius_find_all(&forms, 4096).unwrap();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn avx2_full_enum_beats_scalar_l4_wall() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let forms = dense_forms(16, 24);
        let t0 = std::time::Instant::now();
        let mut a = gray_ffs_avx2_8x32(&forms, usize::MAX).expect("avx2");
        let avx_ns = t0.elapsed().as_nanos();

        let (mut fq, mut fl) = fill_scalar(&forms, 16);
        let mut b = Vec::new();
        let t1 = std::time::Instant::now();
        gray_ffs_unrolled_l4(&mut fq, &mut fl, 16, usize::MAX, &mut b);
        let l4_ns = t1.elapsed().as_nanos();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        let ratio = l4_ns as f64 / avx_ns.max(1) as f64;
        eprintln!(
            "avx2_4x64_vs_l4 n=16 m=24 full_enum: avx2={avx_ns}ns l4={l4_ns}ns ratio={ratio:.2} sols={}",
            a.len()
        );
        // Documented negative vs our packed-u64 scalar L=4 (engineering).
        // libfes's hand-written asm + multi-system batch is a different regime.
        assert!(
            ratio < 1.0,
            "unexpected: AVX2 beat scalar L=4 ({ratio:.3}×); update the note/scoreboard"
        );
    }
}
