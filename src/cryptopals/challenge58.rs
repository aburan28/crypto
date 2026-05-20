//! # Challenge 58 — Pollard's Kangaroo (Lambda) Method
//!
//! Generic discrete-log algorithm whose work factor is `O(√(b - a))`
//! where `[a, b]` is a known *contiguous* range containing the
//! log.  Constant space — no big lookup table.
//!
//! ## How it works
//!
//! Pick a pseudo-random jump function `f: G → ℤ`.  Then:
//!
//! - **Tame kangaroo**: start at `y_T = g^b`, accumulate
//!   `x_T = Σ f(y_T)` over `N` jumps.  End at `y_T = g^(b + x_T)`.
//! - **Wild kangaroo**: start at `y_W = y` (the target).  Take jumps
//!   `y_W ← y_W · g^f(y_W)` while accumulating `x_W = Σ f(y_W)`.
//!   On each step check if `y_W == y_T`; if so, the wild has
//!   joined the tame's path and `index(y) = b + x_T − x_W`.
//!
//! The mean jump size derived from `f(y) = 2^(y mod k)` is
//! `(2^k − 1) / k`, and one tunes `N` proportional to that.
//!
//! Application to subgroup confinement: after Challenge 57 we know
//! `x ≡ n (mod r)`.  Rewriting `y = g^x = g^n · (g^r)^m` lets us
//! reduce to a kangaroo search for `m` in `[0, (q-1)/r)`.

use crate::cryptopals::challenge57::Bob;
use crate::cryptopals::set8_util::{
    biguint_to_bytes_be, crt_combine, hmac_sha256, parse_big, small_factors,
};
use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

/// Parameters specific to Challenge 58 (the "less accommodating" group).
fn params() -> (BigUint, BigUint, BigUint, BigUint) {
    let p = parse_big(
        "11470374874925275658116663507232161402086650258453896274534991676898999262641581519101074740642369848233294239851519212341844337347119899874391456329785623",
    );
    let q = parse_big("335062023296420808191071248367701059461");
    let j = parse_big(
        "34233586850807404623475048381328686211071196701374230492615844865929237417097514638999377942356150481334217896204702",
    );
    let g = parse_big(
        "622952335333961296978159266084741085889881358738459939978290179936063635566740258555167783009058567397963466103140082647486611657350811560630587013183357",
    );
    (p, g, q, j)
}

/// Index of the jump: `i(y) = y mod k`, in `[0, k)`.
fn jump_idx(y: &BigUint, k: u32) -> usize {
    let r = y % BigUint::from(k);
    let limbs = r.to_u64_digits();
    if limbs.is_empty() {
        0
    } else {
        limbs[0] as usize
    }
}

/// Pollard's kangaroo for discrete logs in `Zp*`.
///
/// `g` is the base, `y = g^x mod p` is the target, and `x` is
/// promised to lie in `[a, b)`.  Returns `Some(x)` on success and
/// `None` if the tame–wild trails fail to meet within budget.
///
/// The jump function is `f(y) = 2^(y mod k)` (so jumps come from the
/// fixed set `{1, 2, 4, …, 2^(k-1)}`).  We precompute `g^(2^i) mod p`
/// for `i ∈ [0, k)` so each iteration costs one mod-mul, not a
/// full modpow.
pub fn kangaroo(
    g: &BigUint,
    y: &BigUint,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    k: u32,
) -> Option<BigUint> {
    // Precompute g^(2^i) for i = 0..k.
    let mut g_pow: Vec<BigUint> = Vec::with_capacity(k as usize);
    {
        let mut cur = g.clone() % p;
        for _ in 0..k {
            g_pow.push(cur.clone());
            cur = (&cur * &cur) % p;
        }
    }
    // mean(f) = (2^k - 1)/k ≈ 2^k/k.  Set N := 4 · mean.
    let mean_jump = ((1u64 << k) - 1) / (k as u64).max(1);
    let n_iters = 4u64 * mean_jump;

    // Tame: starts at g^b.
    let mut x_t = BigUint::zero();
    let mut y_t = g.modpow(b, p);
    for _ in 0..n_iters {
        let idx = jump_idx(&y_t, k);
        let f = 1u64 << idx;
        x_t += BigUint::from(f);
        y_t = (&y_t * &g_pow[idx]) % p;
    }

    // Wild: starts at y, runs until x_w would overshoot.
    let mut x_w = BigUint::zero();
    let mut y_w = y.clone();
    let range = b - a;
    let max_w = &range + &x_t;
    while x_w < max_w {
        let idx = jump_idx(&y_w, k);
        let f = 1u64 << idx;
        x_w += BigUint::from(f);
        y_w = (&y_w * &g_pow[idx]) % p;
        if y_w == y_t {
            return Some(b + &x_t - x_w);
        }
    }
    None
}

/// Apply the subgroup-confinement attack to recover `x mod r_total`,
/// then close the residual gap with the kangaroo algorithm.
pub fn attack(bob: &Bob, p: &BigUint, g: &BigUint, q: &BigUint, j: &BigUint) -> Option<BigUint> {
    use crate::cryptopals::challenge57::{element_of_order, recover_residue};
    // Step 1: harvest residues mod small primes of j.
    let factors = small_factors(j, 1 << 16);
    let mut residues: Vec<(BigUint, BigUint)> = Vec::new();
    let mut product = BigUint::one();
    let mut seed: u64 = 7;
    for r in &factors {
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
    if residues.is_empty() {
        return None;
    }
    let (n, r_total) = crt_combine(&residues);
    // Step 2: kangaroo over the residual range.  We seek `m` with
    //   y' = (g')^m, where g' = g^r_total, y' = y · g^(-n).
    // y = g^x is Bob's public; we don't have it directly, but we
    // can ask Bob to MAC with `g` itself — the public is
    // `g^x mod p`.  For test purposes the caller will pass the
    // already-known public; here we compute it via an oracle call.
    let pub_y = g.modpow(&bob_secret(bob, p, g), p);
    let n_inv_exp = q - (&n % q); // -n mod q
    let g_pow_neg_n = g.modpow(&n_inv_exp, p);
    let y_prime = (&pub_y * &g_pow_neg_n) % p;
    let g_prime = g.modpow(&r_total, p);
    let upper = (q - BigUint::one()) / &r_total + BigUint::one();
    let m = kangaroo(&g_prime, &y_prime, &BigUint::zero(), &upper, p, 22)?;
    let x = n + r_total * m;
    Some(x % q)
}

/// Convenience getter — Bob's secret leaks here only because we
/// (the demo) own Bob.  In a real attack the public key comes from
/// the network.
fn bob_secret(bob: &Bob, _p: &BigUint, _g: &BigUint) -> BigUint {
    bob.x_clone()
}

/// One-shot sample lookup that demonstrates the algorithm on the
/// cryptopals example `y` (index in `[0, 2^20]`).
pub fn run() -> Report {
    let mut r = Report::new(58, "Pollard's Kangaroo (Lambda) Method");
    let (p, g, q, j) = params();

    // ── Small example: y with index in [0, 2^20]. ──
    let y_small =
        parse_big("7760073848032689505395005705677365876654629189298052775754597607446617558600394076764814236081991643094239886772481052254010323780165093955236429914607119");
    let small = kangaroo(
        &g,
        &y_small,
        &BigUint::zero(),
        &BigUint::from(1u64 << 20),
        &p,
        12,
    );
    r.line(format!(
        "kangaroo over [0, 2^20): {}",
        match &small {
            Some(x) => format!("found x = {}", x),
            None => "not found".to_string(),
        }
    ));
    // Verify the recovered index.
    if let Some(x) = &small {
        let lhs = g.modpow(x, &p);
        r.line(format!("g^x mod p == y       : {}", lhs == y_small));
        assert_eq!(lhs, y_small);
    }

    // ── End-to-end: subgroup confinement + kangaroo on bigger range. ──
    // Bob's secret is mod q ≈ 2^128; we recover ~80 bits of it from
    // the small-subgroup attack and run the kangaroo on the residue.
    let secret_x = parse_big("214581739081370940921376148093017584711");
    let bob = Bob::new(p.clone(), secret_x.clone());
    r.line("");
    r.line("End-to-end: subgroup confinement → kangaroo on residual");
    r.line(format!("Bob's secret x : {}", secret_x));
    match attack(&bob, &p, &g, &q, &j) {
        Some(x) => {
            r.line(format!("Recovered  x   : {}", x));
            r.line(format!("Match           : {}", x == secret_x));
            if x == secret_x {
                return r.succeed();
            }
        }
        None => r.line("Attack did not converge (try a bigger kangaroo budget)."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kangaroo_small_range_finds_log() {
        let (p, g, _q, _j) = params();
        let y = parse_big("7760073848032689505395005705677365876654629189298052775754597607446617558600394076764814236081991643094239886772481052254010323780165093955236429914607119");
        let x = kangaroo(&g, &y, &BigUint::zero(), &BigUint::from(1u64 << 20), &p, 12)
            .expect("kangaroo must find x in [0, 2^20]");
        assert_eq!(g.modpow(&x, &p), y);
    }

    /// Larger range; allowed to be slow.
    #[test]
    #[ignore]
    fn kangaroo_40_bit_range() {
        let (p, g, _q, _j) = params();
        let y = parse_big("9388897478013399550694114614498790691034187453089355259602614074132918843899833277397448144245883225611726912025846772975325932794909655215329941809013733");
        let x = kangaroo(&g, &y, &BigUint::zero(), &BigUint::from(1u64 << 40), &p, 22)
            .expect("kangaroo must find x in [0, 2^40]");
        assert_eq!(g.modpow(&x, &p), y);
    }
}
