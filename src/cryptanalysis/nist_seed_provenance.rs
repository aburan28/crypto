//! **NIST seed provenance** — the published X9.62 / FIPS 186-2 curve seeds
//! regenerated from the generator that produced them.
//!
//! On 2026-10-15 Anthropic's Frontier Red Team showed (and F. Valsorda's
//! bounty write-up of 2026-10-16 confirmed) that the seeds were not hashes
//! of lost English phrases: each prime-curve seed is one 20-byte block of
//! `(int)(drand48() * 256)` output after `srand48(t)`, with `t` a UNIX
//! time between May 1997 and April 1998, and four of the five binary-curve
//! seeds are five `lrand48()` words after `srand48(t / 2)` on 1998-04-17.
//! This module re-executes that claim natively:
//!
//! 1. [`Drand48`] — the POSIX 48-bit LCG (`x ← 0x5DEECE66D·x + 0xB mod 2^48`),
//!    seeded the way `srand48` does (`x = t << 16 | 0x330E`), with the
//!    inverse step that makes the walk back to the seeding state possible.
//! 2. [`prime_seed`] / [`binary_seed`] — the `n`-th candidate seed of a run.
//! 3. [`PRIME_SEEDS`] / [`BINARY_SEEDS`] — the published seeds with their
//!    run timestamp and candidate index; the tests regenerate every one and
//!    re-check the X9.62 relation `r(seed)·b² ≡ −27 (mod p)` for the prime
//!    curves (`a = −3`).
//!
//! **Status: verified provenance, not an attack.**  Nothing here changes the
//! ECDLP on any curve.  What it settles is the seed-grinding question
//! [`super::p256_backdoor_map`] left open: the freedom a generator had was
//! the choice of start time (about 2^26 seconds in 1997–98) and of which
//! accepted candidate to publish, roughly 30 bits per curve.  Still open:
//! why about half of the valid candidates in each run were skipped, the
//! B-163 seed, and the base points of P-224, P-256, P-384, P-521 and the
//! binary curves.  Report and census: `research/nist_seed_provenance_20261016/`.

use num_bigint::BigUint;

use super::p256_backdoor_map::x962_candidate_r;

// ── The generator ─────────────────────────────────────────────────────

const LCG_A: u64 = 0x5_DEEC_E66D;
const LCG_C: u64 = 0xB;
const LCG_MASK: u64 = (1 << 48) - 1;

/// POSIX `drand48` / `lrand48` state (48 bits).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Drand48 {
    pub state: u64,
}

impl Drand48 {
    /// `srand48(seedval)`: the high 32 bits of the state are the low 32 bits
    /// of `seedval`, the low 16 bits are `0x330E`.
    pub fn srand48(seedval: u32) -> Self {
        Drand48 {
            state: (u64::from(seedval) << 16) | 0x330E,
        }
    }

    /// Advance one step and return the new state.
    pub fn step(&mut self) -> u64 {
        self.state = LCG_A.wrapping_mul(self.state).wrapping_add(LCG_C) & LCG_MASK;
        self.state
    }

    /// Undo one step (the multiplier is odd, so the LCG is a bijection).
    pub fn step_back(&mut self) {
        self.state = self.state.wrapping_sub(LCG_C).wrapping_mul(lcg_a_inverse()) & LCG_MASK;
    }

    /// `(int)(drand48() * 256)`: the top 8 bits of the new state.
    pub fn byte(&mut self) -> u8 {
        (self.step() >> 40) as u8
    }

    /// `lrand48()`: the top 31 bits of the new state.
    pub fn lrand48(&mut self) -> u32 {
        (self.step() >> 17) as u32
    }

    /// Is this a state `srand48` could have produced?
    pub fn is_srand48_form(&self) -> bool {
        self.state & 0xFFFF == 0x330E
    }
}

/// `A^-1 mod 2^48` by Newton iteration (each round doubles the correct bits).
fn lcg_a_inverse() -> u64 {
    let mut inv: u64 = 1;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(LCG_A.wrapping_mul(inv)));
    }
    inv & LCG_MASK
}

/// The `n`-th (1-based) prime-curve candidate of the run `srand48(t)`:
/// bytes `20(n−1) .. 20n` of `(int)(drand48() * 256)`.
pub fn prime_seed(t: u32, n: u32) -> [u8; 20] {
    let mut g = Drand48::srand48(t);
    for _ in 0..20 * (n - 1) {
        g.step();
    }
    let mut out = [0u8; 20];
    for b in out.iter_mut() {
        *b = g.byte();
    }
    out
}

/// The `n`-th (1-based) binary-curve candidate of the run `srand48(seedval)`:
/// five `lrand48()` words, big-endian.
pub fn binary_seed(seedval: u32, n: u32) -> [u8; 20] {
    let mut g = Drand48::srand48(seedval);
    for _ in 0..5 * (n - 1) {
        g.step();
    }
    let mut out = [0u8; 20];
    for chunk in out.chunks_mut(4) {
        chunk.copy_from_slice(&g.lrand48().to_be_bytes());
    }
    out
}

// ── The published seeds ───────────────────────────────────────────────

/// A published prime-curve seed with the run that produced it.
#[derive(Clone, Copy, Debug)]
pub struct PrimeSeed {
    pub name: &'static str,
    /// The `srand48` argument: a UNIX time.
    pub srand48: u32,
    /// 1-based candidate index within the run.
    pub candidate: u32,
    pub seed_hex: &'static str,
    pub p_hex: &'static str,
    pub b_hex: &'static str,
}

/// A published binary-curve seed with the run that produced it
/// (`srand48` argument = UNIX time / 2, if the halved-clock reading is right).
#[derive(Clone, Copy, Debug)]
pub struct BinarySeed {
    pub name: &'static str,
    pub srand48: u32,
    pub candidate: u32,
    pub seed_hex: &'static str,
}

pub const PRIME_SEEDS: [PrimeSeed; 10] = [
    PrimeSeed { name: "prime239v1", srand48: 863711343, candidate: 2185, seed_hex: "e43bb460f0b80cc0c0b075798e948060f8321b7d", p_hex: "7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff", b_hex: "6b016c3bdcf18941d0d654921475ca71a9db2fb27d1d37796185c2942c0a" },
    PrimeSeed { name: "prime239v2", srand48: 863711343, candidate: 2320, seed_hex: "e8b4011604095303ca3b8099982be09fcb9ae616", p_hex: "7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff", b_hex: "617fab6832576cbbfed50d99f0249c3fee58b94ba0038c7ae84c8c832f2c" },
    PrimeSeed { name: "prime239v3", srand48: 863711343, candidate: 3048, seed_hex: "7d7374168ffe3471b60a857686a19475d3bfa2ff", p_hex: "7fffffffffffffffffffffff7fffffffffff8000000000007fffffffffff", b_hex: "255705fa2a306654b1f4cb03d6a750a30c250102d4988717d9ba15ab6d3e" },
    PrimeSeed { name: "P-192", srand48: 863966734, candidate: 139, seed_hex: "3045ae6fc8422f64ed579528d38120eae12196d5", p_hex: "fffffffffffffffffffffffffffffffeffffffffffffffff", b_hex: "64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1" },
    PrimeSeed { name: "prime192v3", srand48: 863966734, candidate: 5228, seed_hex: "c469684435deb378c4b65ca9591e2a5763059a2e", p_hex: "fffffffffffffffffffffffffffffffeffffffffffffffff", b_hex: "22123dc2395a05caa7423daeccc94760a7d462256bd56916" },
    PrimeSeed { name: "prime192v2", srand48: 863966734, candidate: 6013, seed_hex: "31a92ee2029fd10d901b113e990710f0d21ac6b6", p_hex: "fffffffffffffffffffffffffffffffeffffffffffffffff", b_hex: "cc22d6dfb95c6b25e49c0d6364a4e5980c393aa21668d953" },
    PrimeSeed { name: "P-521", srand48: 871002724, candidate: 1136, seed_hex: "d09e8800291cb85396cc6717393284aaa0da64ba", p_hex: "1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", b_hex: "51953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00" },
    PrimeSeed { name: "P-256", srand48: 871696354, candidate: 1842, seed_hex: "c49d360886e704936a6678e1139d26b7819f7e90", p_hex: "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", b_hex: "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b" },
    PrimeSeed { name: "P-384", srand48: 872952608, candidate: 2688, seed_hex: "a335926aa319a27a1d00896a6773a4827acdac73", p_hex: "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff", b_hex: "b3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef" },
    PrimeSeed { name: "P-224", srand48: 891607927, candidate: 7026, seed_hex: "bd71344799d5c7fcdc45b59fa3b9ab8f6a948bc5", p_hex: "ffffffffffffffffffffffffffffffff000000000000000000000001", b_hex: "b4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4" },
];
pub const BINARY_SEEDS: [BinarySeed; 4] = [
    BinarySeed {
        name: "B-233",
        srand48: 446413329,
        candidate: 155,
        seed_hex: "74d59ff07f6b413d0ea14b344b20a2db049b50c3",
    },
    BinarySeed {
        name: "B-283",
        srand48: 446413432,
        candidate: 32,
        seed_hex: "77e2b07370eb0f832a6dd5b62dfc88cd06bb84be",
    },
    BinarySeed {
        name: "B-409",
        srand48: 446413509,
        candidate: 60,
        seed_hex: "4099b5a457f9d69f79213d094c4bcd4d4262210b",
    },
    BinarySeed {
        name: "B-571",
        srand48: 446413564,
        candidate: 73,
        seed_hex: "2aa058f73a0e33ab486b0f610410c53a7f132310",
    },
];

/// Parse a 40-hex-digit seed.
pub fn seed_bytes(hex: &str) -> [u8; 20] {
    let mut out = [0u8; 20];
    for (i, b) in out.iter_mut().enumerate() {
        *b = u8::from_str_radix(&hex[2 * i..2 * i + 2], 16).expect("seed hex");
    }
    out
}

/// `r(seed)·b² + 27 ≡ 0 (mod p)`: the X9.62 A.3.3.1 relation with `a = −3`.
pub fn x962_relation_holds(s: &PrimeSeed) -> bool {
    let p = BigUint::parse_bytes(s.p_hex.as_bytes(), 16).expect("p hex");
    let b = BigUint::parse_bytes(s.b_hex.as_bytes(), 16).expect("b hex");
    let (_, _, r) = x962_candidate_r(&seed_bytes(s.seed_hex), p.bits() as u32);
    (r * &b % &p * &b + BigUint::from(27u8)) % &p == BigUint::from(0u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 1997-01-01T00:00:00Z and 1999-01-01T00:00:00Z.
    const Y1997: u32 = 852_076_800;
    const Y1999: u32 = 915_148_800;

    #[test]
    fn every_prime_seed_is_a_drand48_block() {
        for s in PRIME_SEEDS.iter() {
            assert_eq!(
                prime_seed(s.srand48, s.candidate),
                seed_bytes(s.seed_hex),
                "{}",
                s.name
            );
            assert!(
                (Y1997..Y1999).contains(&s.srand48),
                "{} start time outside 1997-98",
                s.name
            );
        }
    }

    #[test]
    fn every_prime_seed_satisfies_x962_with_published_b() {
        for s in PRIME_SEEDS.iter() {
            assert!(x962_relation_holds(s), "{}", s.name);
        }
    }

    #[test]
    fn four_binary_seeds_are_lrand48_words() {
        for s in BINARY_SEEDS.iter() {
            assert_eq!(
                binary_seed(s.srand48, s.candidate),
                seed_bytes(s.seed_hex),
                "{}",
                s.name
            );
            // Doubled, the srand48 arguments fall on 1998-04-17 (15:24-15:32 UTC).
            let t = 2 * s.srand48;
            assert!((892_771_200..892_857_600).contains(&t), "{}", s.name);
        }
    }

    #[test]
    fn walking_back_from_a_seed_block_reaches_its_srand48_state() {
        // P-256: the state that emitted the seed's first byte, stepped back
        // 20·(1842−1)+1 times, is srand48(871696354) and of srand48 form.
        let s = PRIME_SEEDS.iter().find(|s| s.name == "P-256").unwrap();
        let mut g = Drand48::srand48(s.srand48);
        let steps = 20 * (s.candidate - 1) + 1;
        for _ in 0..steps {
            g.step();
        }
        assert_eq!((g.state >> 40) as u8, seed_bytes(s.seed_hex)[0]);
        for _ in 0..steps {
            g.step_back();
        }
        assert!(g.is_srand48_form());
        assert_eq!(g, Drand48::srand48(s.srand48));
    }

    #[test]
    fn a_published_seed_is_not_an_earlier_block_of_its_run() {
        // Control: the match is at the stated index, not anywhere nearby.
        for s in PRIME_SEEDS.iter() {
            assert_ne!(
                prime_seed(s.srand48, s.candidate - 1),
                seed_bytes(s.seed_hex),
                "{}",
                s.name
            );
            assert_ne!(
                prime_seed(s.srand48 + 1, s.candidate),
                seed_bytes(s.seed_hex),
                "{}",
                s.name
            );
        }
    }
}
