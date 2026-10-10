//! Four Keccak-f[1600] permutations at once, one per 64-bit lane of a 256-bit
//! vector, and the few sponge operations ML-KEM needs on top of them.
//!
//! ML-KEM's hashing is dominated by many short, independent sponges: one
//! SHAKE128 per entry of `Â` and one SHAKE256 per noise polynomial. They share
//! no state, so running four side by side in SIMD lanes costs roughly what one
//! costs in scalar code, a trick that goes back to the Kyber AVX2 submission.
//!
//! Three backends, chosen at run time so the binary stays portable:
//!
//! * **AVX-512VL** (`vprolq`, `vpternlogq`): native 64-bit rotates, and the
//!   χ step `a ^ (!b & c)` and θ's five-way parity are single instructions.
//! * **AVX2**: rotates as shift-shift-or, χ as and-not then xor.
//! * **portable**: four calls to the scalar permutation. Callers ask
//!   [`accelerated`] first and keep their one-at-a-time code when it is false,
//!   so hosts without either extension (Arm64, older x86) pay nothing for
//!   this module's existence.
//!
//! The state is lane-interleaved, `s[word][instance]`, so absorbing and
//! squeezing are plain scalar byte operations on one column and only the
//! permutation touches vectors. The tests check every backend against the
//! scalar permutation on random states.

use super::keccak::keccak_f1600;

/// Round constants, as in the scalar permutation.
#[cfg(target_arch = "x86_64")]
const RC: [u64; 24] = [
    0x0000_0000_0000_0001,
    0x0000_0000_0000_8082,
    0x8000_0000_0000_808a,
    0x8000_0000_8000_8000,
    0x0000_0000_0000_808b,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8009,
    0x0000_0000_0000_008a,
    0x0000_0000_0000_0088,
    0x0000_0000_8000_8009,
    0x0000_0000_8000_000a,
    0x0000_0000_8000_808b,
    0x8000_0000_0000_008b,
    0x8000_0000_0000_8089,
    0x8000_0000_0000_8003,
    0x8000_0000_0000_8002,
    0x8000_0000_0000_0080,
    0x0000_0000_0000_800a,
    0x8000_0000_8000_000a,
    0x8000_0000_8000_8081,
    0x8000_0000_0000_8080,
    0x0000_0000_8000_0001,
    0x8000_0000_8000_8008,
];

/// Four Keccak states, interleaved: `0[w][i]` is word `w` of instance `i`.
#[derive(Clone, Copy)]
#[repr(C, align(32))]
pub struct State4(pub [[u64; 4]; 25]);

impl State4 {
    pub const fn zero() -> Self {
        State4([[0; 4]; 25])
    }
}

/// Which permutation backend this process uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    Avx512,
    Avx2,
    Portable,
}

#[cfg(test)]
thread_local! {
    /// Test-only: force a backend on this thread so the ML-KEM differential
    /// tests can run every code path the host supports, not just the best.
    pub static FORCE: std::cell::Cell<Option<Backend>> = const { std::cell::Cell::new(None) };
}

/// The fastest backend the running CPU supports.
#[inline]
pub fn backend() -> Backend {
    #[cfg(test)]
    if let Some(b) = FORCE.with(|f| f.get()) {
        return b;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f")
            && std::arch::is_x86_feature_detected!("avx512vl")
        {
            return Backend::Avx512;
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            return Backend::Avx2;
        }
    }
    Backend::Portable
}

/// Whether four-way hashing is faster than one-at-a-time on this CPU.
#[inline]
pub fn accelerated() -> bool {
    backend() != Backend::Portable
}

/// Keccak-f[1600] on all four states, with the fastest available backend.
#[inline]
pub fn keccak_f1600_x4(s: &mut State4) {
    keccak_f1600_x4_with(s, backend())
}

/// Keccak-f[1600] on all four states with a named backend. Used directly only
/// by the tests, which run every backend the host supports.
pub fn keccak_f1600_x4_with(s: &mut State4, b: Backend) {
    match b {
        #[cfg(target_arch = "x86_64")]
        Backend::Avx512 => {
            assert!(
                std::arch::is_x86_feature_detected!("avx512f")
                    && std::arch::is_x86_feature_detected!("avx512vl")
            );
            // SAFETY: the features this function is compiled for were just
            // checked on the running CPU.
            unsafe { x86::permute_avx512(s) }
        }
        #[cfg(target_arch = "x86_64")]
        Backend::Avx2 => {
            assert!(std::arch::is_x86_feature_detected!("avx2"));
            // SAFETY: as above.
            unsafe { x86::permute_avx2(s) }
        }
        _ => permute_portable(s),
    }
}

/// Four scalar permutations, through a transpose. The reference the vector
/// backends are tested against, and the fallback.
fn permute_portable(s: &mut State4) {
    for i in 0..4 {
        let mut a = [0u64; 25];
        for w in 0..25 {
            a[w] = s.0[w][i];
        }
        keccak_f1600(&mut a);
        for w in 0..25 {
            s.0[w][i] = a[w];
        }
    }
}

/// The permutation body, written once over a small set of lane operations so
/// that the AVX2 and AVX-512 backends cannot drift apart. `$rol!(x, n)` is a
/// constant rotate left, `$xor5!` five-way xor, `$chi!(a, b, c)` is
/// `a ^ (!b & c)`.
#[cfg(target_arch = "x86_64")]
macro_rules! keccak_rounds {
    ($ld:ident, $st:ident, $xor:path, $bcast:path,
     $rol:ident, $xor5:ident, $chi:ident) => {{
        // `$ld!(w)` loads word `w` of every instance into one vector and
        // `$st!(w, v)` stores it back; the caller is an `unsafe fn` compiled
        // for the features the intrinsics need.
        let mut a00 = $ld!(0);
        let mut a01 = $ld!(1);
        let mut a02 = $ld!(2);
        let mut a03 = $ld!(3);
        let mut a04 = $ld!(4);
        let mut a05 = $ld!(5);
        let mut a06 = $ld!(6);
        let mut a07 = $ld!(7);
        let mut a08 = $ld!(8);
        let mut a09 = $ld!(9);
        let mut a10 = $ld!(10);
        let mut a11 = $ld!(11);
        let mut a12 = $ld!(12);
        let mut a13 = $ld!(13);
        let mut a14 = $ld!(14);
        let mut a15 = $ld!(15);
        let mut a16 = $ld!(16);
        let mut a17 = $ld!(17);
        let mut a18 = $ld!(18);
        let mut a19 = $ld!(19);
        let mut a20 = $ld!(20);
        let mut a21 = $ld!(21);
        let mut a22 = $ld!(22);
        let mut a23 = $ld!(23);
        let mut a24 = $ld!(24);

        for &rc in RC.iter() {
            // θ
            let c0 = $xor5!(a00, a05, a10, a15, a20);
            let c1 = $xor5!(a01, a06, a11, a16, a21);
            let c2 = $xor5!(a02, a07, a12, a17, a22);
            let c3 = $xor5!(a03, a08, a13, a18, a23);
            let c4 = $xor5!(a04, a09, a14, a19, a24);
            let d0 = $xor(c4, $rol!(c1, 1));
            let d1 = $xor(c0, $rol!(c2, 1));
            let d2 = $xor(c1, $rol!(c3, 1));
            let d3 = $xor(c2, $rol!(c4, 1));
            let d4 = $xor(c3, $rol!(c0, 1));

            // θ application fused into ρ and π: b[π(i)] = rot(a[i] ^ d, ρ(i)).
            let b00 = $xor(a00, d0);
            let b01 = $rol!($xor(a06, d1), 44);
            let b02 = $rol!($xor(a12, d2), 43);
            let b03 = $rol!($xor(a18, d3), 21);
            let b04 = $rol!($xor(a24, d4), 14);
            let b05 = $rol!($xor(a03, d3), 28);
            let b06 = $rol!($xor(a09, d4), 20);
            let b07 = $rol!($xor(a10, d0), 3);
            let b08 = $rol!($xor(a16, d1), 45);
            let b09 = $rol!($xor(a22, d2), 61);
            let b10 = $rol!($xor(a01, d1), 1);
            let b11 = $rol!($xor(a07, d2), 6);
            let b12 = $rol!($xor(a13, d3), 25);
            let b13 = $rol!($xor(a19, d4), 8);
            let b14 = $rol!($xor(a20, d0), 18);
            let b15 = $rol!($xor(a04, d4), 27);
            let b16 = $rol!($xor(a05, d0), 36);
            let b17 = $rol!($xor(a11, d1), 10);
            let b18 = $rol!($xor(a17, d2), 15);
            let b19 = $rol!($xor(a23, d3), 56);
            let b20 = $rol!($xor(a02, d2), 62);
            let b21 = $rol!($xor(a08, d3), 55);
            let b22 = $rol!($xor(a14, d4), 39);
            let b23 = $rol!($xor(a15, d0), 41);
            let b24 = $rol!($xor(a21, d1), 2);

            // χ, with ι folded into the first lane.
            a00 = $xor($chi!(b00, b01, b02), $bcast(rc as i64));
            a01 = $chi!(b01, b02, b03);
            a02 = $chi!(b02, b03, b04);
            a03 = $chi!(b03, b04, b00);
            a04 = $chi!(b04, b00, b01);
            a05 = $chi!(b05, b06, b07);
            a06 = $chi!(b06, b07, b08);
            a07 = $chi!(b07, b08, b09);
            a08 = $chi!(b08, b09, b05);
            a09 = $chi!(b09, b05, b06);
            a10 = $chi!(b10, b11, b12);
            a11 = $chi!(b11, b12, b13);
            a12 = $chi!(b12, b13, b14);
            a13 = $chi!(b13, b14, b10);
            a14 = $chi!(b14, b10, b11);
            a15 = $chi!(b15, b16, b17);
            a16 = $chi!(b16, b17, b18);
            a17 = $chi!(b17, b18, b19);
            a18 = $chi!(b18, b19, b15);
            a19 = $chi!(b19, b15, b16);
            a20 = $chi!(b20, b21, b22);
            a21 = $chi!(b21, b22, b23);
            a22 = $chi!(b22, b23, b24);
            a23 = $chi!(b23, b24, b20);
            a24 = $chi!(b24, b20, b21);
        }

        $st!(0, a00);
        $st!(1, a01);
        $st!(2, a02);
        $st!(3, a03);
        $st!(4, a04);
        $st!(5, a05);
        $st!(6, a06);
        $st!(7, a07);
        $st!(8, a08);
        $st!(9, a09);
        $st!(10, a10);
        $st!(11, a11);
        $st!(12, a12);
        $st!(13, a13);
        $st!(14, a14);
        $st!(15, a15);
        $st!(16, a16);
        $st!(17, a17);
        $st!(18, a18);
        $st!(19, a19);
        $st!(20, a20);
        $st!(21, a21);
        $st!(22, a22);
        $st!(23, a23);
        $st!(24, a24);
    }};
}

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::{State4, RC};
    use core::arch::x86_64::*;

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn permute_avx2(s: &mut State4) {
        macro_rules! rol {
            ($x:expr, $n:literal) => {{
                let x = $x;
                _mm256_or_si256(
                    _mm256_slli_epi64::<$n>(x),
                    _mm256_srli_epi64::<{ 64 - $n }>(x),
                )
            }};
        }
        macro_rules! xor5 {
            ($a:expr, $b:expr, $c:expr, $d:expr, $e:expr) => {
                _mm256_xor_si256(
                    _mm256_xor_si256(_mm256_xor_si256($a, $b), _mm256_xor_si256($c, $d)),
                    $e,
                )
            };
        }
        macro_rules! chi {
            ($a:expr, $b:expr, $c:expr) => {
                _mm256_xor_si256($a, _mm256_andnot_si256($b, $c))
            };
        }
        // State4 is 25 contiguous rows of 32 bytes, aligned to 32, so each
        // row is one aligned vector.
        let p = s.0.as_mut_ptr() as *mut __m256i;
        macro_rules! ld {
            ($w:literal) => {
                _mm256_load_si256(p.add($w) as *const __m256i)
            };
        }
        macro_rules! st {
            ($w:literal, $v:expr) => {
                _mm256_store_si256(p.add($w), $v)
            };
        }
        keccak_rounds!(ld, st, _mm256_xor_si256, _mm256_set1_epi64x, rol, xor5, chi);
    }

    macro_rules! rol512 {
        ($x:expr, $n:literal) => {
            _mm256_rol_epi64::<$n>($x)
        };
    }
    // 0x96 is three-way xor; 0xd2 is a ^ (!b & c).
    macro_rules! xor5_512 {
        ($a:expr, $b:expr, $c:expr, $d:expr, $e:expr) => {
            _mm256_ternarylogic_epi64::<0x96>(_mm256_ternarylogic_epi64::<0x96>($a, $b, $c), $d, $e)
        };
    }
    macro_rules! chi512 {
        ($a:expr, $b:expr, $c:expr) => {
            _mm256_ternarylogic_epi64::<0xd2>($a, $b, $c)
        };
    }

    #[target_feature(enable = "avx2,avx512f,avx512vl")]
    pub(super) unsafe fn permute_avx512(s: &mut State4) {
        // State4 is 25 contiguous rows of 32 bytes, aligned to 32, so each
        // row is one aligned vector.
        let p = s.0.as_mut_ptr() as *mut __m256i;
        macro_rules! ld {
            ($w:literal) => {
                _mm256_load_si256(p.add($w) as *const __m256i)
            };
        }
        macro_rules! st {
            ($w:literal, $v:expr) => {
                _mm256_store_si256(p.add($w), $v)
            };
        }
        keccak_rounds!(
            ld,
            st,
            _mm256_xor_si256,
            _mm256_set1_epi64x,
            rol512,
            xor5_512,
            chi512
        );
    }

    /// One state in lane 0; the other lanes compute on whatever the
    /// broadcast loads put there and are never stored.
    #[target_feature(enable = "avx2,avx512f,avx512vl")]
    pub(super) unsafe fn permute1_avx512(a: &mut [u64; 25]) {
        let p = a.as_mut_ptr();
        macro_rules! ld {
            ($w:literal) => {
                _mm256_set1_epi64x(*p.add($w) as i64)
            };
        }
        macro_rules! st {
            ($w:literal, $v:expr) => {
                *p.add($w) = _mm256_extract_epi64::<0>($v) as u64
            };
        }
        keccak_rounds!(
            ld,
            st,
            _mm256_xor_si256,
            _mm256_set1_epi64x,
            rol512,
            xor5_512,
            chi512
        );
    }
}

// ── sponge helpers ───────────────────────────────────────────────────────────

/// XOR one byte into instance `i` at byte offset `off`.
#[inline(always)]
fn xor_byte(s: &mut State4, i: usize, off: usize, b: u8) {
    s.0[off / 8][i] ^= (b as u64) << (8 * (off % 8));
}

/// Absorb one short message per instance (each shorter than the rate), pad
/// with domain byte `suffix`, and permute: the state is then ready to
/// squeeze its first block.
pub fn absorb_short_x4<const RATE: usize>(inputs: [&[u8]; 4], suffix: u8) -> State4 {
    let mut s = State4::zero();
    for (i, inp) in inputs.iter().enumerate() {
        assert!(inp.len() < RATE);
        let whole = inp.len() / 8;
        for w in 0..whole {
            let mut b = [0u8; 8];
            b.copy_from_slice(&inp[8 * w..8 * w + 8]);
            s.0[w][i] ^= u64::from_le_bytes(b);
        }
        for (off, &b) in inp.iter().enumerate().skip(8 * whole) {
            xor_byte(&mut s, i, off, b);
        }
        xor_byte(&mut s, i, inp.len(), suffix);
        xor_byte(&mut s, i, RATE - 1, 0x80);
    }
    keccak_f1600_x4(&mut s);
    s
}

/// Copy the current block of instance `i` (the first `RATE` bytes of its
/// state) into `out`.
#[inline]
pub fn extract_block<const RATE: usize>(s: &State4, i: usize, out: &mut [u8]) {
    debug_assert!(out.len() >= RATE && RATE.is_multiple_of(8));
    for w in 0..RATE / 8 {
        out[8 * w..8 * w + 8].copy_from_slice(&s.0[w][i].to_le_bytes());
    }
}

/// Four SHAKE256 outputs of equal length from four short inputs.
///
/// `out[i]` receives `SHAKE256(inputs[i])[..out[i].len()]`; all outputs have
/// to be the same length, which is the shape ML-KEM's PRF calls take.
pub fn shake256_x4(out: [&mut [u8]; 4], inputs: [&[u8]; 4]) {
    const R: usize = super::keccak::SHAKE256_RATE;
    let len = out[0].len();
    assert!(out.iter().all(|o| o.len() == len));
    let mut s = absorb_short_x4::<R>(inputs, 0x1f);
    let mut block = [0u8; R];
    let mut done = 0;
    let [o0, o1, o2, o3] = out;
    let outs: [&mut [u8]; 4] = [o0, o1, o2, o3];
    let mut outs = outs;
    loop {
        let take = (len - done).min(R);
        for (i, o) in outs.iter_mut().enumerate() {
            extract_block::<R>(&s, i, &mut block);
            o[done..done + take].copy_from_slice(&block[..take]);
        }
        done += take;
        if done == len {
            break;
        }
        keccak_f1600_x4(&mut s);
    }
}

/// One Keccak-f[1600] through the four-way AVX-512 permutation, using lane 0.
///
/// With `vprolq` and `vpternlogq` four lanes cost about 610 cycles on the
/// development host against about 880 for one scalar permutation, so even a
/// single sponge is faster in one lane of four. Returns false, leaving the
/// state alone, when AVX-512VL is not available.
#[inline]
pub fn keccak_f1600_lane0(a: &mut [u64; 25]) -> bool {
    if backend() != Backend::Avx512 {
        return false;
    }
    #[cfg(target_arch = "x86_64")]
    {
        // SAFETY: backend() returned Avx512 only after detecting the
        // features permute1_avx512 is compiled for.
        unsafe { x86::permute1_avx512(a) };
        true
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = a;
        false
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::pqc::fast::keccak::{shake256_into, SHAKE256_RATE};

    pub(crate) fn backends() -> Vec<Backend> {
        #[cfg_attr(not(target_arch = "x86_64"), allow(unused_mut))]
        let mut v = vec![Backend::Portable];
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx2") {
                v.push(Backend::Avx2);
            }
            if std::arch::is_x86_feature_detected!("avx512f")
                && std::arch::is_x86_feature_detected!("avx512vl")
            {
                v.push(Backend::Avx512);
            }
        }
        v
    }

    /// Every backend the host has, against the scalar permutation applied to
    /// each instance separately, from unstructured states and through several
    /// chained permutations so a wrong rotation cannot cancel out.
    #[test]
    fn every_backend_matches_four_scalar_permutations() {
        let mut x = 0x0123_4567_89ab_cdefu64;
        for b in backends() {
            for _ in 0..32 {
                let mut s = State4::zero();
                for w in 0..25 {
                    for i in 0..4 {
                        x ^= x << 13;
                        x ^= x >> 7;
                        x ^= x << 17;
                        s.0[w][i] = x;
                    }
                }
                let mut r = s;
                for _ in 0..3 {
                    keccak_f1600_x4_with(&mut s, b);
                    for i in 0..4 {
                        let mut a = [0u64; 25];
                        for w in 0..25 {
                            a[w] = r.0[w][i];
                        }
                        keccak_f1600(&mut a);
                        for w in 0..25 {
                            r.0[w][i] = a[w];
                        }
                    }
                }
                assert_eq!(s.0, r.0, "{b:?}");
            }
        }
    }

    /// The four-way SHAKE256 is four SHAKE256s, at the lengths ML-KEM uses
    /// and across a block boundary, with inputs of different lengths so the
    /// padding position differs between lanes.
    #[test]
    fn shake256_x4_is_four_shake256s() {
        let msgs: [Vec<u8>; 4] = [
            (0..33).map(|i| i as u8).collect(),
            (0..33).map(|i| (i * 7 + 1) as u8).collect(),
            (0..64).map(|i| (i * 3) as u8).collect(),
            vec![],
        ];
        for len in [1usize, 32, 128, 135, 136, 137, 192, 300] {
            let mut o = [
                vec![0u8; len],
                vec![0u8; len],
                vec![0u8; len],
                vec![0u8; len],
            ];
            {
                let [a, b, c, d] = &mut o;
                shake256_x4([a, b, c, d], [&msgs[0], &msgs[1], &msgs[2], &msgs[3]]);
            }
            for i in 0..4 {
                let mut want = vec![0u8; len];
                shake256_into(&mut want, &msgs[i]);
                assert_eq!(o[i], want, "lane {i}, length {len}");
            }
        }
        let _ = SHAKE256_RATE;
    }

    #[test]
    fn lane0_permutation_matches_scalar() {
        let mut a = [0u64; 25];
        let mut x = 0xdead_beef_0bad_f00du64;
        for lane in a.iter_mut() {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            *lane = x;
        }
        let mut b = a;
        for _ in 0..16 {
            if !keccak_f1600_lane0(&mut a) {
                return; // no AVX-512 on this host: nothing to compare
            }
            keccak_f1600(&mut b);
            assert_eq!(a, b);
        }
    }
}
