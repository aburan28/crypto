//! The `secpkk1` ECDLP challenge ladder of Dallaire-Demers, Doyle and Foo
//! (arXiv 2508.14011, "Brace for impact: ECDLP challenges for quantum
//! cryptanalysis", Section 2, Algorithm 3), rebuilt in Rust for
//! `k ≤ 62` bits so the repository's generic solvers can be run on it.
//!
//! For a bit length `k` the instance is:
//!
//! - **Prime.**  Walk the `k`-bit primes `p < 2^k` downwards and take the
//!   first for which `E: y² = x³ + 7` over `F_p` has prime order `n ≠ p`.
//!   For `p ≡ 2 (mod 3)` the curve is supersingular with `p + 1` points
//!   (even), so only `p ≡ 1 (mod 3)` can qualify; there `E` has CM by
//!   `Z[ω]` and its order is one of the six `p + 1 − t` with
//!   `4p = t² + 3v²`, found here by Cornacchia and settled by checking
//!   which candidate annihilates points of `E` (the paper uses SEA; for a
//!   `j = 0` curve the CM count is exact).
//! - **Points.**  For the labels `"Quantum"` (`P`) and `"Challenge"`
//!   (`Q`): `h = SHA256(label) mod p` (digest read big-endian),
//!   `x = h + δ` for the least `δ ≥ 0` with `x³ + 7` a square (or zero),
//!   and `y` the root in `[0, (p−1)/2]`.
//! - **Task.**  Recover `κ` with `P = [κ]Q`.
//!
//! The published cards at 16, 24, 32 and 48 bits are checked in the tests,
//! so the derivation is pinned to the paper's own values.

use num_bigint::BigUint;

use crate::hash::sha256::sha256;

/// `a·b mod m`.
#[inline]
fn mulm(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// `a^e mod m`.
fn powm(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, a, m);
        }
        a = mulm(a, a, m);
        e >>= 1;
    }
    r
}

/// Deterministic Miller–Rabin for `u64`.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    'witness: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// Square root mod an odd prime by Tonelli–Shanks; `None` for a
/// non-residue.
pub fn sqrt_mod(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if powm(a, (p - 1) / 2, p) != 1 {
        return None;
    }
    if p % 4 == 3 {
        return Some(powm(a, (p + 1) / 4, p));
    }
    let mut q = p - 1;
    let mut s = 0u32;
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let mut z = 2u64;
    while powm(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = powm(z, q, p);
    let mut t = powm(a, q, p);
    let mut r = powm(a, q.div_ceil(2), p);
    while t != 1 {
        let mut i = 0u32;
        let mut tt = t;
        while tt != 1 {
            tt = mulm(tt, tt, p);
            i += 1;
        }
        let b = powm(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mulm(b, b, p);
        t = mulm(t, c, p);
        r = mulm(r, b, p);
    }
    Some(r)
}

fn isqrt(v: u64) -> u64 {
    let mut s = (v as f64).sqrt() as u64;
    while s > 0 && s.saturating_mul(s) > v {
        s -= 1;
    }
    while (s + 1).saturating_mul(s + 1) <= v {
        s += 1;
    }
    s
}

/// Cornacchia: `(a, b)` with `a² + 3b² = p`, for a prime `p ≡ 1 (mod 3)`.
pub fn cornacchia3(p: u64) -> Option<(u64, u64)> {
    let mut r0 = sqrt_mod(p - 3, p)?;
    if r0 <= p / 2 {
        r0 = p - r0;
    }
    let (mut a, mut b) = (p, r0);
    let lim = isqrt(p);
    while b > lim {
        let t = a % b;
        a = b;
        b = t;
    }
    let rest = p - b * b;
    if !rest.is_multiple_of(3) {
        return None;
    }
    let c = isqrt(rest / 3);
    (c * c == rest / 3).then_some((b, c))
}

/// An affine point of `y² = x³ + 7` over `F_p` (canonical integers), or
/// `None` for `O`.
pub type Pt = Option<(u64, u64)>;

/// `P + Q` on `y² = x³ + 7` over `F_p`.
pub fn add(p: u64, a: Pt, b: Pt) -> Pt {
    let (x1, y1) = match a {
        None => return b,
        Some(v) => v,
    };
    let (x2, y2) = match b {
        None => return a,
        Some(v) => v,
    };
    let lam = if x1 == x2 {
        if (y1 + y2) % p == 0 {
            return None;
        }
        // 3x² / 2y
        let num = mulm(3, mulm(x1, x1, p), p);
        mulm(num, powm(mulm(2, y1, p), p - 2, p), p)
    } else {
        let num = (y2 + p - y1) % p;
        let den = (x2 + p - x1) % p;
        mulm(num, powm(den, p - 2, p), p)
    };
    let x3 = (mulm(lam, lam, p) + 2 * p - x1 - x2) % p;
    let y3 = (mulm(lam, (x1 + p - x3) % p, p) + p - y1) % p;
    Some((x3, y3))
}

/// `[k]P`.
pub fn mul(p: u64, pt: Pt, mut k: u64) -> Pt {
    let mut acc = None;
    let mut base = pt;
    while k > 0 {
        if k & 1 == 1 {
            acc = add(p, acc, base);
        }
        base = add(p, base, base);
        k >>= 1;
    }
    acc
}

/// `[k]P` for a 128-bit `k` (the order candidates can exceed `u64` only
/// in principle; `p + 1 + 2√p` fits for `p < 2^62`).
fn mul_u128(p: u64, pt: Pt, k: u128) -> Pt {
    let mut acc = None;
    let mut base = pt;
    let mut k = k;
    while k > 0 {
        if k & 1 == 1 {
            acc = add(p, acc, base);
        }
        base = add(p, base, base);
        k >>= 1;
    }
    acc
}

/// The six CM order candidates `p + 1 − t` of `y² = x³ + b` over `F_p`,
/// `p ≡ 1 (mod 3)`.
pub fn j0_order_candidates(p: u64) -> Option<Vec<u128>> {
    let (a, b) = cornacchia3(p)?;
    let (a, b) = (a as i128, b as i128);
    let ts = [
        2 * a,
        -2 * a,
        a + 3 * b,
        -(a + 3 * b),
        a - 3 * b,
        -(a - 3 * b),
    ];
    Some(ts.iter().map(|t| (p as i128 + 1 - t) as u128).collect())
}

/// The order of `y² = x³ + 7` over `F_p`: the unique CM candidate that
/// annihilates a handful of deterministic points.  `None` when `p` is not
/// `1 mod 3` or the candidates cannot be separated.
pub fn order_y2_x3_7(p: u64) -> Option<u128> {
    if p % 3 != 1 {
        return None;
    }
    let cands = j0_order_candidates(p)?;
    let mut alive: Vec<u128> = cands.clone();
    alive.sort();
    alive.dedup();
    let mut x = 1u64;
    let mut tested = 0;
    while alive.len() > 1 || tested < 3 {
        let rhs = (mulm(mulm(x, x, p), x, p) + 7) % p;
        if let Some(y) = sqrt_mod(rhs, p) {
            let pt = Some((x, y));
            alive.retain(|&n| mul_u128(p, pt, n).is_none());
            tested += 1;
            if alive.is_empty() || tested > 40 {
                return None;
            }
        }
        x += 1;
    }
    alive.first().copied()
}

/// One rung of the ladder.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChallengeInstance {
    pub bits: u32,
    pub p: u64,
    /// `|E(F_p)|`, prime.
    pub n: u64,
    /// The `"Quantum"` point (the target).
    pub p_pt: (u64, u64),
    /// The `"Challenge"` point (the generator).
    pub q_pt: (u64, u64),
    /// `δ` used for `P` and `Q`.
    pub delta_p: u64,
    pub delta_q: u64,
    /// Primes below `2^k` rejected before `p` was accepted.
    pub primes_rejected: u64,
    /// `ord_n(p)` if it is at most 64, else `None` (reported as `> 64`).
    pub embedding_degree_le_64: Option<u32>,
}

impl ChallengeInstance {
    /// Compressed SEC1 of a point, `x` padded to `⌈k/8⌉` bytes.
    pub fn sec1(&self, pt: (u64, u64)) -> String {
        let len = (self.bits as usize).div_ceil(8);
        let pre = if pt.1.is_multiple_of(2) { "02" } else { "03" };
        format!("{pre}{:0width$X}", pt.0, width = 2 * len)
    }
}

/// The NUMS point for `label` on `y² = x³ + 7` over `F_p`.
pub fn nums_point(p: u64, label: &str) -> ((u64, u64), u64) {
    let h = BigUint::from_bytes_be(&sha256(label.as_bytes())) % BigUint::from(p);
    let h = h.iter_u64_digits().next().unwrap_or(0);
    let mut delta = 0u64;
    loop {
        let x = ((h as u128 + delta as u128) % p as u128) as u64;
        let rho = (mulm(mulm(x, x, p), x, p) + 7) % p;
        if let Some(y) = sqrt_mod(rho, p) {
            let y = y.min((p - y) % p);
            return ((x, y), delta);
        }
        delta += 1;
    }
}

/// Which way the prime search runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrimeSearch {
    /// Algorithm 3 as printed: `k`-bit primes in descending order from `2^k`.
    Descending,
    /// Diagnostic: ascending from `2^{k−1}`.  The paper's published 8- and
    /// 16-bit cards (`p = 163`, `32803`) are not the first descending hits
    /// (`211`, `65419`); this variant tests whether they are ascending ones.
    Ascending,
}

/// Algorithm 3 of arXiv 2508.14011 at `bits ≤ 62`.
pub fn derive(bits: u32) -> Option<ChallengeInstance> {
    derive_with(bits, PrimeSearch::Descending)
}

/// [`derive`] with the prime search direction chosen.
pub fn derive_with(bits: u32, dir: PrimeSearch) -> Option<ChallengeInstance> {
    assert!((6..=62).contains(&bits), "u64 ladder covers 6..=62 bits");
    let lo = 1u64 << (bits - 1);
    let hi = (1u64 << bits) - 1;
    let mut p = if dir == PrimeSearch::Descending {
        hi
    } else {
        lo + 1
    };
    let mut rejected = 0u64;
    while p >= lo && p <= hi {
        if is_prime_u64(p) {
            if let Some(n) = order_y2_x3_7(p) {
                if n <= u64::MAX as u128 && is_prime_u64(n as u64) && n as u64 != p {
                    let n = n as u64;
                    let (p_pt, dp) = nums_point(p, "Quantum");
                    let (q_pt, dq) = nums_point(p, "Challenge");
                    let pm = p % n;
                    let mut acc = pm;
                    let mut r = None;
                    for t in 1..=64u32 {
                        if acc == 1 {
                            r = Some(t);
                            break;
                        }
                        acc = mulm(acc, pm, n);
                    }
                    return Some(ChallengeInstance {
                        bits,
                        p,
                        n,
                        p_pt,
                        q_pt,
                        delta_p: dp,
                        delta_q: dq,
                        primes_rejected: rejected,
                        embedding_degree_le_64: r,
                    });
                }
            }
            rejected += 1;
        }
        if dir == PrimeSearch::Descending {
            p -= 2;
        } else {
            p += 2;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The paper's Appendix A cards (p, n, P, Q).  The 6-, 24-, 32- and 48-bit cards are the first descending hits, as
    /// Algorithm 3 is printed; the 8-, 12- and 16-bit cards are the first
    /// *ascending* hits from `2^{k−1}` (the descending search gives
    /// different primes, e.g. `211` and `65419`), so the published ladder
    /// mixes two search directions.
    #[test]
    fn matches_published_cards() {
        for (bits, p, n, pp, qq) in [
            (8u32, 163u64, 139u64, "0205", "0261"),
            (12, 2089, 2143, "0207C0", "03053E"),
            (16, 32803, 32497, "0260A3", "036D88"),
        ] {
            let c = derive_with(bits, PrimeSearch::Ascending).expect("instance");
            assert_eq!((c.p, c.n), (p, n), "{bits}-bit p, n (ascending)");
            assert_eq!(c.sec1(c.p_pt), pp, "{bits}-bit P");
            assert_eq!(c.sec1(c.q_pt), qq, "{bits}-bit Q");
            let d = derive(bits).unwrap();
            assert_ne!(d.p, p, "{bits}: the descending search differs");
        }
        let cards: [(u32, u64, u64, &str, &str); 4] = [
            (6, 43, 31, "0307", "0320"),
            (24, 16_777_213, 16_770_451, "0276BA79", "0385303B"),
            (32, 4_294_966_177, 4_294_835_173, "021E0EF3EB", "02EF01BFF4"),
            (
                48,
                281_474_976_707_983,
                281_474_961_723_409,
                "03769EC7F35857",
                "02210644F01A41",
            ),
        ];
        for (bits, p, n, pp, qq) in cards {
            let c = derive(bits).expect("instance");
            assert_eq!((c.p, c.n), (p, n), "{bits}-bit p, n");
            assert_eq!(c.sec1(c.p_pt), pp, "{bits}-bit P");
            assert_eq!(c.sec1(c.q_pt), qq, "{bits}-bit Q");
            assert!(mul(p, Some(c.q_pt), n).is_none(), "{bits}: nQ = O");
        }
    }

    #[test]
    fn cornacchia_and_sqrt() {
        for p in [7u64, 13, 31, 4_294_966_177, 281_474_976_707_983] {
            if p % 3 == 1 {
                let (a, b) = cornacchia3(p).unwrap();
                assert_eq!(a * a + 3 * b * b, p);
            }
        }
        let p = 1_000_000_007u64; // p ≡ 3 mod 4 is not this one: 1e9+7 ≡ 3
        for a in 1..200u64 {
            if let Some(r) = sqrt_mod(a, p) {
                assert_eq!(mulm(r, r, p), a);
            }
        }
        let p = 998_244_353u64; // 2^23 · 119 + 1, exercises Tonelli–Shanks
        for a in 1..200u64 {
            if let Some(r) = sqrt_mod(a, p) {
                assert_eq!(mulm(r, r, p), a);
            }
        }
    }
}
