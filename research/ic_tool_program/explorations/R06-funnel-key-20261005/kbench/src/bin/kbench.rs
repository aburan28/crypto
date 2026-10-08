//! Per-key cost of three AVX-512 kernels on 16 normal-coordinate words:
//! the chained least rotation, the funnel least rotation, and a
//! rotation-invariant popcount hash (popcount and K autocorrelations).
use std::arch::x86_64::*;
use std::time::Instant;

#[target_feature(enable = "avx512f")]
unsafe fn chained(c: &[u64; 16], n: u32, mask: u64) -> [u64; 16] {
    let m = _mm512_set1_epi64(mask as i64);
    let sl = _mm_cvtsi64_si128(1);
    let sr = _mm_cvtsi64_si128(i64::from(n - 1));
    let mut va = _mm512_loadu_si512(c.as_ptr().cast());
    let mut vb = _mm512_loadu_si512(c.as_ptr().add(8).cast());
    let (mut ba, mut bb) = (va, vb);
    for _ in 1..n {
        va = _mm512_ternarylogic_epi64::<0xA8>(_mm512_sll_epi64(va, sl), _mm512_srl_epi64(va, sr), m);
        vb = _mm512_ternarylogic_epi64::<0xA8>(_mm512_sll_epi64(vb, sl), _mm512_srl_epi64(vb, sr), m);
        ba = _mm512_min_epu64(ba, va);
        bb = _mm512_min_epu64(bb, vb);
    }
    let mut out = [0u64; 16];
    _mm512_storeu_si512(out.as_mut_ptr().cast(), ba);
    _mm512_storeu_si512(out.as_mut_ptr().add(8).cast(), bb);
    out
}

#[inline(always)]
unsafe fn doubled(v: __m512i, n: u32) -> (__m512i, __m512i) {
    let count = |s: u32| _mm_cvtsi64_si128(i64::from(s));
    if n > 32 {
        (_mm512_or_si512(_mm512_sll_epi64(v, count(64 - n)), _mm512_srl_epi64(v, count(2 * n - 64))),
         _mm512_sll_epi64(v, count(128 - 2 * n)))
    } else {
        (_mm512_sll_epi64(_mm512_or_si512(_mm512_sll_epi64(v, count(n)), v), count(64 - 2 * n)),
         _mm512_setzero_si512())
    }
}

#[target_feature(enable = "avx512f,avx512vbmi2")]
unsafe fn funnel(c: &[u64; 16], n: u32) -> [u64; 16] {
    let count = |s: u32| _mm_cvtsi64_si128(i64::from(s));
    let (hi_a, lo_a) = doubled(_mm512_loadu_si512(c.as_ptr().cast()), n);
    let (hi_b, lo_b) = doubled(_mm512_loadu_si512(c.as_ptr().add(8).cast()), n);
    let (mut best_a, mut best_b) = (hi_a, hi_b);
    macro_rules! rotations {
        ($($t:literal)*) => {$(
            if $t < n {
                best_a = _mm512_min_epu64(best_a, _mm512_shldi_epi64::<$t>(hi_a, lo_a));
                best_b = _mm512_min_epu64(best_b, _mm512_shldi_epi64::<$t>(hi_b, lo_b));
            }
        )*};
    }
    rotations!(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
        31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
        60 61 62);
    let mut out = [0u64; 16];
    _mm512_storeu_si512(out.as_mut_ptr().cast(), _mm512_srl_epi64(best_a, count(64 - n)));
    _mm512_storeu_si512(out.as_mut_ptr().add(8).cast(), _mm512_srl_epi64(best_b, count(64 - n)));
    out
}

/// popcount(c) and popcount(c & rotl^k c) for k = 1..=K, folded by
/// rotate-6-and-xor.  Rotations by funnel shifts of the doubled word,
/// masked to the top n bits by AND with the top-aligned c.
#[target_feature(enable = "avx512f,avx512vbmi2,avx512vpopcntdq")]
unsafe fn invariant<const K: u32>(c: &[u64; 16], n: u32) -> [u64; 16] {
    let top = _mm512_set1_epi64((u64::MAX << (64 - n)) as i64);
    let (hi_a, lo_a) = doubled(_mm512_loadu_si512(c.as_ptr().cast()), n);
    let (hi_b, lo_b) = doubled(_mm512_loadu_si512(c.as_ptr().add(8).cast()), n);
    let ca = _mm512_and_si512(hi_a, top);
    let cb = _mm512_and_si512(hi_b, top);
    let mut acc_a = _mm512_popcnt_epi64(ca);
    let mut acc_b = _mm512_popcnt_epi64(cb);
    macro_rules! ks {
        ($($t:literal)*) => {$(
            if $t <= K {
                acc_a = _mm512_xor_si512(_mm512_rol_epi64::<6>(acc_a),
                    _mm512_popcnt_epi64(_mm512_and_si512(_mm512_shldi_epi64::<$t>(hi_a, lo_a), ca)));
                acc_b = _mm512_xor_si512(_mm512_rol_epi64::<6>(acc_b),
                    _mm512_popcnt_epi64(_mm512_and_si512(_mm512_shldi_epi64::<$t>(hi_b, lo_b), cb)));
            }
        )*};
    }
    ks!(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16);
    let mut out = [0u64; 16];
    _mm512_storeu_si512(out.as_mut_ptr().cast(), acc_a);
    _mm512_storeu_si512(out.as_mut_ptr().add(8).cast(), acc_b);
    out
}

fn main() {
    let ns: Vec<u32> = std::env::args().skip(1).map(|a| a.parse().unwrap()).collect();
    let mut s = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    const KEYS: usize = 1 << 20;
    for n in ns {
        let mask = (1u64 << n) - 1;
        let words: Vec<u64> = (0..KEYS).map(|_| next() & mask).collect();
        let mut results = Vec::new();
        for kernel in ["chained", "funnel", "inv10", "inv12", "inv14"] {
            let mut best = f64::MAX;
            let mut check = 0u64;
            for _ in 0..7 {
                let t = Instant::now();
                let mut x = 0u64;
                for chunk in words.chunks_exact(16) {
                    let c: [u64; 16] = chunk.try_into().unwrap();
                    let out = unsafe {
                        match kernel {
                            "chained" => chained(&c, n, mask),
                            "funnel" => funnel(&c, n),
                            "inv10" => invariant::<10>(&c, n),
                            "inv12" => invariant::<12>(&c, n),
                            _ => invariant::<14>(&c, n),
                        }
                    };
                    for v in out { x = x.wrapping_add(v); }
                }
                let ns_key = t.elapsed().as_nanos() as f64 / KEYS as f64;
                best = best.min(ns_key);
                check = x;
            }
            results.push(format!("{kernel} {best:.2} ns/key (sum {check:#x})"));
        }
        println!("n={n}: {}", results.join(", "));
    }
}
