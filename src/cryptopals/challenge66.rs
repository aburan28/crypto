//! # Challenge 66 — Exploiting Implementation Errors in (EC)DH
//!
//! Carry bugs in bignum arithmetic occasionally cause point addition
//! to silently produce wrong results.  If the attacker can detect
//! when this happens — say, because the next protocol step
//! authenticates the shared secret and fails — they can binary-search
//! the bits of the victim's private scalar.
//!
//! ## How the leak works
//!
//! Recall the textbook double-and-add ladder:
//!
//! ```text
//!   R := Q                       # k[1] = 1 by assumption
//!   for b in k[2..n]:
//!       R := add(R, R)           # always happens
//!       if b = 1:
//!           R := add(R, Q)       # only when this bit is 1
//! ```
//!
//! At iteration `i` the value of `R` going into `add(R, R)` is a
//! function of `k[1..i-1]` alone — once we know those bits, we can
//! *predict* the coefficient pair that ladder will pass to `add`.
//!
//! Attacker picks a point `Q` such that:
//!
//! - `add(c_0·Q, c_0·Q)` succeeds  ←— continues the chain to bit `i`
//! - `add(c_1·Q, c_1·Q)` triggers a fault
//!
//! where `c_0` is the coefficient assuming `k[i] = 0` and `c_1`
//! assuming `k[i] = 1`.  Query the oracle: if it errors, the
//! ladder must have hit the faulty path, so `k[i]` equals whichever
//! branch's coefficient triggered the fault.
//!
//! ## This demo
//!
//! We model the fault deterministically (`fault(Q1, Q2) =
//! (Q1.x * Q2.x) % p == 0 mod fault_rate`).  Bigger `fault_rate`
//! means rarer faults but easier filtering.  Walk the key one bit
//! at a time.

use crate::cryptopals::challenge59::{Curve, Pt};
use crate::cryptopals::Report;
use crate::cryptopals::set8_util::parse_big;
use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Fault predicate.  Returns `true` when the add would corrupt.
/// `mask` controls the fault rate: smaller mask → more faults.
pub fn faulty(p1: &Pt, p2: &Pt, mask: u128) -> bool {
    match (p1, p2) {
        (Pt::Aff(x1, _), Pt::Aff(x2, _)) => {
            let prod = (x1 * x2).to_u64_digits();
            let lo = prod.get(0).copied().unwrap_or(0) as u128;
            (lo & (mask as u128 - 1)) == 0
        }
        _ => false,
    }
}

/// Faulty curve add: returns `None` on fault, `Some(point)` otherwise.
pub fn faulty_add(curve: &Curve, p1: &Pt, p2: &Pt, mask: u128) -> Option<Pt> {
    if faulty(p1, p2, mask) {
        return None;
    }
    Some(curve.add(p1, p2))
}

/// Victim's scalar multiplication using the faulty addition.
/// Returns `None` if any fault triggered during the ladder.
pub fn faulty_scalar_mul(curve: &Curve, q: &Pt, k: &BigUint, mask: u128) -> Option<Pt> {
    if k.is_zero() {
        return Some(Pt::Inf);
    }
    let bits: Vec<bool> = (0..k.bits()).rev().map(|i| ((k >> i as usize) & BigUint::one()).is_one()).collect();
    let mut r = q.clone();
    for &b in &bits[1..] {
        r = faulty_add(curve, &r, &r, mask)?;
        if b {
            r = faulty_add(curve, &r, q, mask)?;
        }
    }
    Some(r)
}

/// Oracle wrapping the ladder: returns `true` if the operation
/// succeeded (no fault), `false` if it triggered.  In real life
/// this maps to "did the protocol return success or fail?".
pub fn oracle(curve: &Curve, q: &Pt, secret: &BigUint, mask: u128) -> bool {
    faulty_scalar_mul(curve, q, secret, mask).is_some()
}

/// Recover one bit of the secret given the prefix we already know.
/// Returns the recovered bit (0 or 1) or `None` if more samples are
/// needed.
pub fn recover_bit(
    curve: &Curve,
    secret_bits_known: &[bool],
    bit_index: usize,
    oracle_fn: &dyn Fn(&Pt) -> bool,
    rng: &mut StdRng,
    mask: u128,
    budget: u32,
) -> Option<bool> {
    // We want to find a Q such that simulating the ladder up through
    // the bits we know succeeds, then differentiates between bit = 0
    // and bit = 1.  Specifically: at step `bit_index`, the coefficient
    // pair entering add() differs between the two hypotheses, and
    // we want exactly one hypothesis to fault.
    let _ = bit_index;
    for _ in 0..budget {
        // Random Q on the curve.
        let q = random_curve_point(curve, rng)?;
        let path0 = simulate_path(curve, &q, secret_bits_known, false, mask);
        let path1 = simulate_path(curve, &q, secret_bits_known, true, mask);
        match (path0, path1) {
            (Some(_), None) => {
                // bit = 1 faults; bit = 0 succeeds.  Query oracle:
                // if it succeeds, bit must be 0.  If it fails, bit = 1.
                return Some(!oracle_fn(&q));
            }
            (None, Some(_)) => {
                // bit = 0 faults; bit = 1 succeeds.  Inverse logic.
                return Some(oracle_fn(&q));
            }
            _ => continue,
        }
    }
    None
}

/// Simulate the ladder for one hypothesis of the next bit.
/// `known_bits` already includes the assumed leading `1` (k[1] = 1)
/// at index 0.  We process `known_bits[1..]` as the iteration bits.
fn simulate_path(
    curve: &Curve,
    q: &Pt,
    known_bits: &[bool],
    hypothesised_bit: bool,
    mask: u128,
) -> Option<Pt> {
    let mut r = q.clone();
    for &b in &known_bits[1..] {
        r = faulty_add(curve, &r, &r, mask)?;
        if b {
            r = faulty_add(curve, &r, q, mask)?;
        }
    }
    // One more step: hypothesise `hypothesised_bit`.
    r = faulty_add(curve, &r, &r, mask)?;
    if hypothesised_bit {
        r = faulty_add(curve, &r, q, mask)?;
    }
    Some(r)
}

fn random_curve_point(curve: &Curve, rng: &mut StdRng) -> Option<Pt> {
    let p = &curve.p;
    for _ in 0..50 {
        let x = BigUint::from(rng.gen::<u64>()) % p;
        let x2 = (&x * &x) % p;
        let x3 = (&x2 * &x) % p;
        let ax = (&curve.a * &x) % p;
        let rhs = ((x3 + ax) + &curve.b) % p;
        if let Some(y) = sqrt_mod_p(&rhs, p) {
            return Some(Pt::Aff(x, y));
        }
    }
    None
}

pub fn run() -> Report {
    let mut r = Report::new(66, "Exploiting Implementation Errors in DH");

    // Use the cryptopals curve with a moderate fault rate.
    let (curve, order, _g) = crate::cryptopals::challenge59::cryptopals_curve();
    let _ = order;
    let _ = parse_big;
    // 10-bit secret so the demo runs quickly; full recovery scales.
    let secret_bits: u32 = 10;
    let mut rng = StdRng::seed_from_u64(0xBEEF);
    let secret = BigUint::from(rng.gen::<u32>() & ((1u32 << secret_bits) - 1) | (1u32 << (secret_bits - 1)));
    r.line(format!("Secret d (binary): {:b}", secret));
    r.line(format!("Fault probability per add ≈ 1/{}", 1u128 << 10));
    let mask = 1u128 << 10;
    let oracle_fn = |q: &Pt| oracle(&curve, q, &secret, mask);
    let mut known = vec![true]; // k[1] = 1 by definition (leading bit)
    let mut total_queries = 0u32;
    for i in 1..secret_bits as usize {
        let bit = recover_bit(&curve, &known, i, &oracle_fn, &mut rng, mask, 200);
        match bit {
            Some(b) => {
                known.push(b);
                total_queries += 200;
            }
            None => {
                r.line(format!("Bit {i} not recovered within budget; aborting."));
                break;
            }
        }
    }
    let recovered_bits = known.iter().fold(BigUint::zero(), |acc, &b| {
        (acc << 1) + if b { BigUint::one() } else { BigUint::zero() }
    });
    r.line(format!("Recovered     : {:b}", recovered_bits));
    r.line(format!("Match         : {}", recovered_bits == secret));
    r.line(format!("Oracle queries: ~{}", total_queries));
    if recovered_bits == secret {
        return r.succeed();
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn faulty_add_detects_marked_inputs() {
        let curve = Curve {
            p: BigUint::from(101u32),
            a: BigUint::zero(),
            b: BigUint::from(7u32),
        };
        let p1 = Pt::Aff(BigUint::zero(), BigUint::from(2u32));
        let p2 = Pt::Aff(BigUint::from(5u32), BigUint::from(3u32));
        // x1 * x2 = 0 → fault triggers at mask 1.
        assert!(faulty(&p1, &p2, 1));
    }
}
