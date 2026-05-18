//! # Challenge 57 — Diffie-Hellman Subgroup-Confinement Attack
//!
//! Classic small-subgroup attack on multiplicative-group DH.  Alice
//! and Bob have agreed on `(p, g, q)` where `q | p − 1` and `g` has
//! order `q`.  The catch: `j = (p−1)/q` has lots of *small* prime
//! factors.
//!
//! Eve plays the role of Bob's peer and, for each small prime
//! `r | j`:
//!
//! 1. Builds a point `h` of order `r` by computing
//!    `h = rand(1, p)^((p−1)/r) mod p` (re-rolling if it lands on 1).
//! 2. Sends `h` as her "public key".  Bob computes
//!    `K = h^x mod p` where `x` is *Bob's* secret.
//! 3. Because `h` has order `r`, `K` has at most `r` possible
//!    values.  Eve brute-forces them against the MAC Bob returned,
//!    recovering `x mod r`.
//!
//! After collecting enough congruences `x ≡ b_i (mod r_i)` such
//! that `∏ r_i > q`, Eve applies CRT to lift them into `x mod ∏ r_i`,
//! which (since the product exceeds `q`) is the full `x`.
//!
//! ## Defence
//!
//! Validate that every received public key has order `q`.  Or
//! switch to a group whose order is itself prime (e.g. ECC over a
//! prime-order curve).

use crate::cryptopals::set8_util::{
    biguint_to_bytes_be, crt_combine, hmac_sha256, parse_big, small_factors,
};
use crate::cryptopals::Report;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::SeedableRng;

/// Cryptopals' published parameters.
fn params() -> (BigUint, BigUint, BigUint, BigUint) {
    let p = parse_big(
        "7199773997391911030609999317773941274322764333428698921736339643928346453700085358802973900485592910475480089726140708102474957429903531369589969318716771",
    );
    let g = parse_big(
        "4565356397095740655436854503483826832136106141639563487732438195343690437606117828318042418238184896212352329118608100083187535033402010599512641674644143",
    );
    let q = parse_big("236234353446506858198510045061214171961");
    let j = parse_big(
        "30477252323177606811760882179058908038824640750610513771646768011063128035873508507547741559514324673960576895059570",
    );
    (p, g, q, j)
}

/// Bob's role: holds secret `x ∈ [1, q)` and, when handed a peer
/// public key `h`, computes the shared secret `K = h^x mod p` and
/// returns `(msg, HMAC-SHA256(K, msg))`.
///
/// In a *correct* implementation Bob would also verify that `h^q = 1`
/// before using `h`.  The whole point of this challenge is that Bob
/// *doesn't*.
pub struct Bob {
    pub(crate) p: BigUint,
    pub(crate) x: BigUint,
}

impl Bob {
    pub fn new(p: BigUint, x: BigUint) -> Self {
        Self { p, x }
    }

    /// Sign a fixed message using a shared secret derived from `h^x`.
    /// Mirrors the cryptopals protocol: same `m`, MAC under
    /// `HMAC-SHA256` keyed by the shared secret's big-endian bytes.
    pub fn oracle(&self, h: &BigUint) -> ([u8; 32], &'static [u8]) {
        let k = h.modpow(&self.x, &self.p);
        let msg: &[u8] = b"crazy flamboyant for the rap enjoyment";
        let key = biguint_to_bytes_be(&k, 192);
        (hmac_sha256(&key, msg), msg)
    }

    /// Test-only accessor: in the real protocol Bob would never
    /// reveal `x`.  The Challenge 58 demo uses this to compute
    /// Bob's public key for the kangaroo step (in a real attack
    /// that key would be observed on the wire).
    pub fn x_clone(&self) -> BigUint {
        self.x.clone()
    }
}

/// Find an element of order exactly `r` in `Zp*`, where `r | (p-1)`.
///
/// Standard recipe: pick a random `g_rand ∈ [2, p-1]`, set
/// `h := g_rand^((p-1)/r) mod p`, and re-roll if `h == 1`.  Almost
/// every randomisation succeeds on first try because the number of
/// `r`-th non-residues mod `p` vastly outnumbers the residues.
pub fn element_of_order(p: &BigUint, r: &BigUint, rng_seed: u64) -> BigUint {
    use rand::rngs::StdRng;
    let mut rng = StdRng::seed_from_u64(rng_seed);
    let exp = (p - BigUint::one()) / r;
    loop {
        let g_rand = rng.gen_biguint_range(&BigUint::from(2u32), &(p - BigUint::one()));
        let h = g_rand.modpow(&exp, p);
        if !h.is_one() {
            return h;
        }
    }
}

/// Brute-force a small-subgroup MAC: try every `b ∈ [0, r)` and
/// check whether `HMAC-SHA256(h^b, msg) == tag`.  Returns `b` when
/// it lands.
pub fn recover_residue(
    p: &BigUint,
    h: &BigUint,
    r: &BigUint,
    msg: &[u8],
    tag: &[u8; 32],
) -> Option<BigUint> {
    let mut b = BigUint::zero();
    // Iterate by multiplication so each step is one modmul.
    let mut k = BigUint::one();
    while &b < r {
        let key = biguint_to_bytes_be(&k, 192);
        let got = hmac_sha256(&key, msg);
        if &got == tag {
            return Some(b);
        }
        k = (&k * h) % p;
        b += 1u32;
    }
    None
}

/// Run the full attack against `bob`, returning the recovered key.
pub fn attack(bob: &Bob, p: &BigUint, q: &BigUint, j: &BigUint) -> Option<BigUint> {
    // 1. Find distinct small prime factors of `j` up to 2^16.
    let factors = small_factors(j, 1 << 16);
    // 2. For each factor r, build a point of order r and read off
    //    `x mod r` from Bob.  Stop once ∏ r > q.
    let mut residues: Vec<(BigUint, BigUint)> = Vec::new();
    let mut product = BigUint::one();
    let mut seed: u64 = 1;
    for r in &factors {
        if &product > q {
            break;
        }
        // Keep r small enough for brute force.
        if *r > BigUint::from(1u64 << 20) {
            continue;
        }
        let h = element_of_order(p, r, seed);
        seed = seed.wrapping_add(1);
        let (tag, msg) = bob.oracle(&h);
        let b = recover_residue(p, &h, r, msg, &tag)?;
        residues.push((b, r.clone()));
        product *= r;
    }
    if &product <= q {
        return None;
    }
    let (a, _) = crt_combine(&residues);
    Some(a % q)
}

pub fn run() -> Report {
    let mut r = Report::new(57, "Diffie-Hellman Subgroup-Confinement Attack");
    let (p, g, q, j) = params();
    // Pick Bob's secret deterministically so the demo is reproducible.
    let secret_x = parse_big("123456789012345678901234567890123456789");
    let bob = Bob::new(p.clone(), secret_x.clone());
    let _ = g;
    r.line(format!("Bob's secret x  : {}", secret_x));
    r.line(format!("q (order of g)  : {}", q));
    match attack(&bob, &p, &q, &j) {
        Some(x) => {
            r.line(format!("Recovered x     : {}", x));
            r.line(format!("Match           : {}", x == secret_x));
            assert_eq!(x, secret_x);
            r.succeed()
        }
        None => {
            r.line("Attack failed — not enough small-subgroup coverage.");
            r
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_actually_has_target_order() {
        let (p, _g, q, j) = params();
        let factors = small_factors(&j, 1 << 16);
        let r = factors[3].clone(); // a non-trivial small factor
        let h = element_of_order(&p, &r, 42);
        // h^r = 1 (mod p), and h ≠ 1.
        assert_eq!(h.modpow(&r, &p), BigUint::one());
        assert_ne!(h, BigUint::one());
        let _ = q;
    }

    #[test]
    fn attack_recovers_full_key() {
        let (p, _g, q, j) = params();
        let secret_x = parse_big("123456789012345678901234567890123456789");
        let bob = Bob::new(p.clone(), secret_x.clone());
        let x = attack(&bob, &p, &q, &j).unwrap();
        assert_eq!(x, secret_x);
    }
}
