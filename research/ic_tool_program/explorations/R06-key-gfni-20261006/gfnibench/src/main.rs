//! The normal-basis change of `FrobeniusCanon::coords`, by byte tables (as
//! main does) and by GFNI (`vgf2p8affineqb` on a byte-transposed block of
//! eight keys).  Same random n x n matrix for both; every output compared.
use std::arch::x86_64::*;
use std::time::Instant;

struct Tables {
    bytes: usize,
    t: Vec<[u64; 256]>,
}

fn splitmix(s: &mut u64) -> u64 {
    *s = s.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *s;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

impl Tables {
    /// Column `j` of the map is `by_bit[j]`, an n-bit word.
    fn new(by_bit: &[u64], n: usize) -> Self {
        let bytes = n.div_ceil(8);
        let t = (0..bytes)
            .map(|bi| {
                let mut table = [0u64; 256];
                for (b, slot) in table.iter_mut().enumerate() {
                    let mut acc = 0;
                    for c in 0..8 {
                        let j = bi * 8 + c;
                        if (b >> c) & 1 == 1 && j < n {
                            acc ^= by_bit[j];
                        }
                    }
                    *slot = acc;
                }
                table
            })
            .collect();
        Tables { bytes, t }
    }
    #[inline]
    fn coords(&self, x: u64) -> u64 {
        let mut acc = 0;
        for (i, table) in self.t.iter().enumerate() {
            acc ^= table[((x >> (8 * i)) & 0xff) as usize];
        }
        acc
    }
    /// The GFNI matrices: `m[i]`'s lane `k` is the 8x8 block from input
    /// byte `i` to output byte `k`, row `r` (output bit) in byte `7 - r`.
    fn gfni(&self) -> [[u64; 8]; 8] {
        let mut m = [[0u64; 8]; 8];
        for i in 0..self.bytes {
            for k in 0..8 {
                let mut q = 0u64;
                for r in 0..8 {
                    let mut row = 0u64;
                    for c in 0..8 {
                        row |= ((self.t[i][1 << c] >> (8 * k + r)) & 1) << c;
                    }
                    q |= row << (8 * (7 - r));
                }
                m[i][k] = q;
            }
        }
        m
    }
}

/// Byte `8i + k` of the result is byte `8k + i` of the input.
fn transpose_index() -> [u8; 64] {
    let mut idx = [0u8; 64];
    for i in 0..8 {
        for k in 0..8 {
            idx[8 * i + k] = (8 * k + i) as u8;
        }
    }
    idx
}

#[target_feature(enable = "avx512f,avx512bw,avx512vbmi,gfni")]
unsafe fn coords8_gfni(x: __m512i, m: &[__m512i; 8], bytes: usize, tr: __m512i) -> __m512i {
    let t = _mm512_permutexvar_epi8(tr, x);
    let mut acc = _mm512_setzero_si512();
    for (i, mi) in m.iter().enumerate().take(bytes) {
        let bi = _mm512_permutexvar_epi64(_mm512_set1_epi64(i as i64), t);
        acc = _mm512_xor_si512(acc, _mm512_gf2p8affine_epi64_epi8::<0>(bi, *mi));
    }
    _mm512_permutexvar_epi8(tr, acc)
}

#[target_feature(enable = "avx512f,avx512bw,avx512vbmi,gfni")]
unsafe fn coords_gfni_all(xs: &[u64], out: &mut [u64], mq: &[[u64; 8]; 8], bytes: usize) {
    let tr = _mm512_loadu_si512(transpose_index().as_ptr().cast());
    let m: [__m512i; 8] = std::array::from_fn(|i| _mm512_loadu_si512(mq[i].as_ptr().cast()));
    for (c, o) in xs.chunks_exact(8).zip(out.chunks_exact_mut(8)) {
        let x = _mm512_loadu_si512(c.as_ptr().cast());
        _mm512_storeu_si512(o.as_mut_ptr().cast(), coords8_gfni(x, &m, bytes, tr));
    }
}

fn main() {
    assert!(is_x86_feature_detected!("gfni") && is_x86_feature_detected!("avx512vbmi"));
    let len = 1 << 20;
    for n in [41usize, 53, 59, 61, 63] {
        let mut s = 0x6f6e_6900 + n as u64;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let by_bit: Vec<u64> = (0..n).map(|_| splitmix(&mut s) & mask).collect();
        let tables = Tables::new(&by_bit, n);
        let mq = tables.gfni();
        let xs: Vec<u64> = (0..len).map(|_| splitmix(&mut s) & mask).collect();
        let mut a = vec![0u64; len];
        let mut b = vec![0u64; len];
        let (mut best_t, mut best_g) = (f64::MAX, f64::MAX);
        for _ in 0..7 {
            let t0 = Instant::now();
            for (o, &x) in a.iter_mut().zip(&xs) {
                *o = tables.coords(x);
            }
            best_t = best_t.min(t0.elapsed().as_nanos() as f64 / len as f64);
            std::hint::black_box(&a);
            let t0 = Instant::now();
            unsafe { coords_gfni_all(&xs, &mut b, &mq, tables.bytes) };
            best_g = best_g.min(t0.elapsed().as_nanos() as f64 / len as f64);
            std::hint::black_box(&b);
        }
        let same = a == b;
        let sum = a.iter().fold(0u64, |h, &v| h.rotate_left(5) ^ v);
        println!("{{\"n\": {n}, \"tables_ns\": {best_t:.3}, \"gfni_ns\": {best_g:.3}, \"identical\": {same}, \"checksum\": \"{sum:016x}\"}}");
        assert!(same, "n = {n}: GFNI coords differ from the tables'");
    }
}
