//! # Challenge 55 — MD4 Collisions (Wang's attack)
//!
//! Xiaoyun Wang & Hongbo Yu's 2005 attack on MD4 finds a collision
//! in expected `2^8` MD4 evaluations — eight orders of magnitude
//! below the `2^64` birthday bound.
//!
//! The attack is a *differential* attack: fix a precise XOR
//! difference between two messages `M` and `M'`, then force the
//! internal state to follow a specific path so the difference
//! cancels out by the end of round 3.
//!
//! ## Differential (Wang 2005)
//!
//! ```text
//!     m'[1]  = m[1]  + 2^31
//!     m'[2]  = m[2]  + 2^31 - 2^28
//!     m'[12] = m[12] - 2^16
//! ```
//!
//! ## Sufficient conditions
//!
//! For the differential to propagate to a collision, the
//! intermediate state at every step must satisfy a list of bit-level
//! conditions (set/clear/equal-another-bit).  Wang's table lists
//! ~120 conditions across 3 rounds — too many to bake into a fast
//! attack, but the first-round ones (~50) can all be enforced
//! deterministically by *modifying the corresponding message word*
//! after each step.  Remaining (round-2/3) conditions hold
//! probabilistically, giving a success rate ≈ `2^-25` per random
//! retry once round-1 is satisfied.
//!
//! This implementation enforces the round-1 conditions via the
//! "single-step modification" recipe and then loops until the
//! differential rounds out to a real collision.  We do **not**
//! implement the round-2 multi-step modifications described in
//! Wang's §4.3 — those are necessary to hit the headline `2^8`
//! complexity.  Without them, success rate per random trial is
//! dominated by ~12 unenforced round-2 conditions, so the
//! collision search may take minutes (≈ 2^20–2^28 trials).
//! The `#[ignore]`d test exercises the full search.
//!
//! ## Reference
//!
//! - Wang, Lai, Feng, Chen, Yu — *Cryptanalysis of the Hash
//!   Functions MD4 and RIPEMD*, EUROCRYPT 2005.
//! - Naito, Sasaki, Kunihiro, Ohta — *Improved Collision Search
//!   for MD4*, IWSEC 2006.

use crate::cryptopals::Report;
use crate::hash::md4::md4_compress;

#[inline] fn f(x: u32, y: u32, z: u32) -> u32 { (x & y) | (!x & z) }

const IV: [u32; 4] = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476];
const SHIFTS_R1: [u32; 16] = [3, 7, 11, 19, 3, 7, 11, 19, 3, 7, 11, 19, 3, 7, 11, 19];

/// Set bit `b` of `x` to 0.
#[inline] fn clr(x: u32, b: u32) -> u32 { x & !(1 << b) }
/// Set bit `b` of `x` to 1.
#[inline] fn set(x: u32, b: u32) -> u32 { x | (1 << b) }
/// Copy bit `b` from `src` into `x`.
#[inline]
fn copy_bit(x: u32, src: u32, b: u32) -> u32 {
    (x & !(1 << b)) | (src & (1 << b))
}

/// Apply Wang's round-1 sufficient conditions to register `q_new`
/// (the *output* of step i), using already-computed registers
/// (a,b,c,d) at this step.  Bits are 0-indexed.
///
/// Returns the masked `q_new`.  The caller then solves for the
/// message word that would have produced this value.
fn apply_conditions(step: usize, q: u32, a: u32, b: u32, c: u32, d: u32) -> u32 {
    // a,b,c,d are this step's input registers (i.e. the previous
    // (a,b,c,d) tuple at the start of step `step`).
    //
    // We use Wang's table (Table 6 of the EUROCRYPT 2005 paper),
    // 0-indexed: bit_k in Wang = bit_(k-1) for us.
    let mut q = q;
    match step {
        0 => {
            // a1: a1[6] = b0[6]
            q = copy_bit(q, b, 6);
        }
        1 => {
            // d1: d1[6]=0, d1[7]=a1[7], d1[10]=a1[10]
            // After step 0, a-slot in step 1 is d0; b-slot is a1.
            // So `b` here is a1.
            q = clr(q, 6);
            q = copy_bit(q, b, 7);
            q = copy_bit(q, b, 10);
        }
        2 => {
            // c1: bit6=1, bit7=1, bit10=0, bit25=d1[25]
            // b slot here = d1.
            q = set(q, 6);
            q = set(q, 7);
            q = clr(q, 10);
            q = copy_bit(q, b, 25);
        }
        3 => {
            // b1: bit6=1, bit7=0, bit10=0, bit25=0
            q = set(q, 6);
            q = clr(q, 7);
            q = clr(q, 10);
            q = clr(q, 25);
        }
        4 => {
            // a2: bit7=1, bit10=1, bit25=0, bit13=b1[13]
            // b slot here = b1.
            q = set(q, 7);
            q = set(q, 10);
            q = clr(q, 25);
            q = copy_bit(q, b, 13);
        }
        5 => {
            // d2: bit13=0, bit18=a2[18], bit19=a2[19], bit20=a2[20],
            //     bit21=a2[21], bit25=1
            // b slot = a2.
            q = clr(q, 13);
            q = copy_bit(q, b, 18);
            q = copy_bit(q, b, 19);
            q = copy_bit(q, b, 20);
            q = copy_bit(q, b, 21);
            q = set(q, 25);
        }
        6 => {
            // c2: bit12=d2[12], bit13=0, bit14=d2[14], bit18=0,
            //     bit19=0, bit20=1, bit21=0
            // b slot = d2.
            q = copy_bit(q, b, 12);
            q = clr(q, 13);
            q = copy_bit(q, b, 14);
            q = clr(q, 18);
            q = clr(q, 19);
            q = set(q, 20);
            q = clr(q, 21);
        }
        7 => {
            // b2: bit12=1, bit13=1, bit14=0, bit16=c2[16], bit18=0,
            //     bit19=0, bit20=0, bit21=0
            // b slot = c2.
            q = set(q, 12);
            q = set(q, 13);
            q = clr(q, 14);
            q = copy_bit(q, b, 16);
            q = clr(q, 18);
            q = clr(q, 19);
            q = clr(q, 20);
            q = clr(q, 21);
        }
        8 => {
            // a3: bit12=1, bit13=1, bit14=1, bit16=0, bit18=0,
            //     bit19=0, bit20=0, bit22=b2[22], bit21=1, bit25=b2[25]
            // b slot = b2.
            q = set(q, 12);
            q = set(q, 13);
            q = set(q, 14);
            q = clr(q, 16);
            q = clr(q, 18);
            q = clr(q, 19);
            q = clr(q, 20);
            q = copy_bit(q, b, 22);
            q = set(q, 21);
            q = copy_bit(q, b, 25);
        }
        9 => {
            // d3: bit16=0, bit19=a3[19], bit22=0, bit12=a3[12],
            //     bit13=a3[13], bit14=a3[14], bit20=0, bit21=1,
            //     bit25=1, bit29=a3[29]
            // b slot = a3.
            q = clr(q, 16);
            q = copy_bit(q, b, 19);
            q = clr(q, 22);
            q = copy_bit(q, b, 12);
            q = copy_bit(q, b, 13);
            q = copy_bit(q, b, 14);
            q = clr(q, 20);
            q = set(q, 21);
            q = set(q, 25);
            q = copy_bit(q, b, 29);
        }
        10 => {
            // c3: bit16=1, bit19=0, bit20=0, bit21=0, bit22=0,
            //     bit25=1, bit29=1, bit31=d3[31]
            // b slot = d3.
            q = set(q, 16);
            q = clr(q, 19);
            q = clr(q, 20);
            q = clr(q, 21);
            q = clr(q, 22);
            q = set(q, 25);
            q = set(q, 29);
            q = copy_bit(q, b, 31);
        }
        11 => {
            // b3: bit19=0, bit20=1, bit21=1, bit22=c3[22], bit25=1,
            //     bit29=0, bit31=0
            // b slot = c3.
            q = clr(q, 19);
            q = set(q, 20);
            q = set(q, 21);
            q = copy_bit(q, b, 22);
            q = set(q, 25);
            q = clr(q, 29);
            q = clr(q, 31);
        }
        12 => {
            // a4: bit22=0, bit25=0, bit26=b3[26], bit28=b3[28],
            //     bit29=1, bit31=0
            // b slot = b3.
            q = clr(q, 22);
            q = clr(q, 25);
            q = copy_bit(q, b, 26);
            q = copy_bit(q, b, 28);
            q = set(q, 29);
            q = clr(q, 31);
        }
        13 => {
            // d4: bit22=0, bit25=0, bit26=1, bit28=1, bit29=0, bit31=1
            q = clr(q, 22);
            q = clr(q, 25);
            q = set(q, 26);
            q = set(q, 28);
            q = clr(q, 29);
            q = set(q, 31);
        }
        14 => {
            // c4: bit18=d4[18], bit22=1, bit25=1, bit26=0,
            //     bit28=0, bit29=0  (Wang's table 6, 1-indexed bits
            //     19, 23, 26, 27, 29, 30).
            // b slot = d4.
            q = copy_bit(q, b, 18);
            q = set(q, 22);
            q = set(q, 25);
            q = clr(q, 26);
            q = clr(q, 28);
            q = clr(q, 29);
        }
        15 => {
            // b4: bit18=0, bit25=c4[25], bit26=c4[26], bit28=c4[28],
            //     bit29=c4[29]  (Wang's table 6 — five equal-to-c4
            //     conditions; bit 19 zero).
            // b slot = c4.
            q = clr(q, 18);
            q = copy_bit(q, b, 25);
            q = copy_bit(q, b, 26);
            q = copy_bit(q, b, 28);
            q = copy_bit(q, b, 29);
        }
        _ => {}
    }
    let _ = (a, c, d);
    q
}

/// Run one MD4 step (round 1) and return `(q_new, m_used)`.
/// Always solves for the message word that makes `q_new` satisfy
/// the conditions.
fn forward_step(
    step: usize,
    a: u32, b: u32, c: u32, d: u32,
    m: u32,
) -> (u32, u32) {
    // First, compute what the natural next register would be.
    let s = SHIFTS_R1[step];
    let raw = a.wrapping_add(f(b, c, d)).wrapping_add(m).rotate_left(s);
    // Apply conditions.
    let target = apply_conditions(step, raw, a, b, c, d);
    // Solve for m: m = (target >>> s) - a - F(b,c,d).
    let new_m = target
        .rotate_right(s)
        .wrapping_sub(a)
        .wrapping_sub(f(b, c, d));
    (target, new_m)
}

/// Apply Wang's round-1 single-step message modification to a
/// 16-word random message.  Returns the modified message.
pub fn round1_modify(msg: &[u32; 16]) -> [u32; 16] {
    let mut out = *msg;
    let (mut a, mut b, mut c, mut d) = (IV[0], IV[1], IV[2], IV[3]);
    for step in 0..16 {
        let (q, new_m) = forward_step(step, a, b, c, d, out[step]);
        out[step] = new_m;
        // MD4 register rotation: (a,b,c,d) ← (d, q, b, c)
        let next_a = d;
        let next_d = c;
        let next_c = b;
        let next_b = q;
        a = next_a;
        b = next_b;
        c = next_c;
        d = next_d;
    }
    out
}

/// Differential: produce `M'` from `M` per Wang's table.
pub fn apply_differential(msg: &[u32; 16]) -> [u32; 16] {
    let mut out = *msg;
    out[1] = out[1].wrapping_add(1u32 << 31);
    out[2] = out[2].wrapping_add((1u32 << 31).wrapping_sub(1u32 << 28));
    out[12] = out[12].wrapping_sub(1u32 << 16);
    out
}

/// Hash a 16-word message as one MD4 compression-function call
/// (no padding, no length-encoding — we operate directly on the
/// block as Wang does).
fn compress_words(m: &[u32; 16]) -> [u32; 4] {
    let mut block = [0u8; 64];
    for (i, w) in m.iter().enumerate() {
        block[4 * i..4 * i + 4].copy_from_slice(&w.to_le_bytes());
    }
    let mut state = IV;
    md4_compress(&mut state, &block);
    state
}

/// Try `n` times to find a Wang collision.  Returns `Some((M, M'))`
/// on success or `None` if budget runs out.
///
/// Each iteration: pick fresh random message, apply round-1
/// modification, derive M', test for collision.  Success rate per
/// iteration depends on how many round-2/3 conditions the random
/// message happens to satisfy — typically a few in 2^25, but
/// occasionally much luckier.
pub fn find_collision(max_tries: u64) -> Option<([u32; 16], [u32; 16])> {
    // Use the same RNG as the rest of the repo to avoid pulling new
    // crates.
    use rand::{Rng, SeedableRng};
    use rand::rngs::StdRng;
    let mut rng = StdRng::from_entropy();

    for _ in 0..max_tries {
        let mut m = [0u32; 16];
        for w in m.iter_mut() {
            *w = rng.gen();
        }
        let m = round1_modify(&m);
        let m_prime = apply_differential(&m);
        if compress_words(&m) == compress_words(&m_prime) && m != m_prime {
            return Some((m, m_prime));
        }
    }
    None
}

pub fn run() -> Report {
    let mut r = Report::new(55, "MD4 Collisions (Wang)");
    r.line("Searching for Wang collision (single-block, MD4 differential)...");
    r.line("(Round-1 single-step modifications only; success per trial is");
    r.line(" probabilistic.  See module docs.)");
    // Generous budget — round-1-only Wang search is much slower
    // than the headline 2^8 figure.  ~30 s on modern hardware.
    let result = find_collision(1 << 26);
    match result {
        Some((m, m_prime)) => {
            let m_bytes = words_to_bytes(&m);
            let mp_bytes = words_to_bytes(&m_prime);
            r.line(format!("M  = {}", hex::encode(&m_bytes)));
            r.line(format!("M' = {}", hex::encode(&mp_bytes)));
            r.line(format!(
                "MD4(M)  = {}",
                hex::encode_upper(state_to_bytes(&compress_words(&m)))
            ));
            r.line(format!(
                "MD4(M') = {}",
                hex::encode_upper(state_to_bytes(&compress_words(&m_prime)))
            ));
            r.line(format!(
                "differences M⊕M' = {:08x} (word 1), {:08x} (word 2), {:08x} (word 12)",
                m[1] ^ m_prime[1],
                m[2] ^ m_prime[2],
                m[12] ^ m_prime[12]
            ));
            r.succeed()
        }
        None => {
            r.line("(no collision in budget — try again; success rate ≈ 2^-25)");
            r
        }
    }
}

fn words_to_bytes(w: &[u32; 16]) -> [u8; 64] {
    let mut out = [0u8; 64];
    for (i, x) in w.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&x.to_le_bytes());
    }
    out
}

fn state_to_bytes(s: &[u32; 4]) -> [u8; 16] {
    let mut out = [0u8; 16];
    for (i, x) in s.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&x.to_le_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round1_modification_satisfies_conditions() {
        // After applying round1_modify, the chain of intermediate
        // registers must match Wang's table on every step.
        let msg = [0xdeadbeefu32; 16];
        let m = round1_modify(&msg);
        let (mut a, mut b, mut c, mut d) = (IV[0], IV[1], IV[2], IV[3]);
        for step in 0..16 {
            let s = SHIFTS_R1[step];
            let q = a
                .wrapping_add(f(b, c, d))
                .wrapping_add(m[step])
                .rotate_left(s);
            let expected = apply_conditions(step, q, a, b, c, d);
            assert_eq!(q, expected, "step {step}: condition not met");
            // Rotate registers like MD4 does.
            let nb = q;
            let na = d;
            let nd = c;
            let nc = b;
            a = na;
            b = nb;
            c = nc;
            d = nd;
        }
    }

    #[test]
    fn differential_words_match() {
        // The differential is purely additive on words 1, 2, 12.
        let m = [1u32; 16];
        let mp = apply_differential(&m);
        for i in 0..16 {
            if i == 1 {
                assert_eq!(mp[1], m[1].wrapping_add(1 << 31));
            } else if i == 2 {
                assert_eq!(
                    mp[2],
                    m[2].wrapping_add((1u32 << 31).wrapping_sub(1u32 << 28))
                );
            } else if i == 12 {
                assert_eq!(mp[12], m[12].wrapping_sub(1u32 << 16));
            } else {
                assert_eq!(mp[i], m[i]);
            }
        }
    }

    // This test runs the full attack and is allowed to be slow.
    // Round-1-only modifications give success rate ≈ 2^-25, so
    // we mark this `#[ignore]` for the default test run and let
    // CI / curious users invoke it explicitly with `--ignored`.
    #[test]
    #[ignore]
    fn finds_real_collision() {
        let r = find_collision(1 << 28).expect("Wang attack must find collision");
        let (m, mp) = r;
        assert_ne!(m, mp);
        assert_eq!(compress_words(&m), compress_words(&mp));
    }
}
