//! Integer factorisation and modular square roots for the Cornacchia step of the fixed-degree
//! algorithm (Benčina–Kutas–Merz–Petit–Stopar–Weitkämper 2023).
//!
//! * [`factor`] / [`factor_bounded`]: trial division by the primes below 1000, then (for what is
//!   left) `field::factor_u64` when it fits in 64 bits, otherwise Pollard rho with Brent's cycle
//!   detection and batched gcds over [`Int`], with Miller–Rabin (`Int::is_probable_prime`) as the
//!   primality test. `factor_bounded` gives up (returns `None`) once the total number of rho
//!   iterations exceeds the budget, which is how callers abandon hard-to-factor targets.
//!   Primality is probabilistic (Miller–Rabin), so a "prime" factor above 2^64 is a probable prime.
//! * [`sqrt_mod_prime_power`]: *all* square roots of `a` modulo `p^k`, including `p = 2` and
//!   `p | a` (non-unit roots), by Tonelli–Shanks + Hensel lifting on the unit part.
//! * [`sqrt_mod_composite`]: all square roots modulo `n` from the factorisation of `n`, by CRT.
//!   Restriction: the number of roots is returned in full; it is `prod_i r_i` where `r_i` is the
//!   number of roots modulo `p_i^{k_i}`, which is at most 4 for unit `a` but can be as large as
//!   `p^{floor(k/2)}` when `p^k | a`. Callers that pass such `a` get all of them.
use crate::field::factor_u64;
use crate::int::Int;

/// Prime factorisation as `(prime, exponent)` pairs, primes increasing.
pub type Factors = Vec<(Int, u32)>;

fn small_primes() -> Vec<u64> {
    let n = 1000usize;
    let mut sieve = vec![true; n + 1];
    sieve[0] = false;
    sieve[1] = false;
    let mut i = 2;
    while i * i <= n {
        if sieve[i] {
            let mut j = i * i;
            while j <= n {
                sieve[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    (2..=n).filter(|&k| sieve[k]).map(|k| k as u64).collect()
}

/// Number of Pollard-rho iterations allowed by [`factor`] (effectively unbounded for the sizes
/// used here: a 2^40 prime factor needs about 2^20 iterations).
pub const DEFAULT_RHO_BUDGET: u64 = 1 << 26;

/// Full factorisation of `|n|` (n != 0). Panics if the default rho budget is exhausted.
pub fn factor(n: &Int) -> Factors {
    factor_bounded(n, DEFAULT_RHO_BUDGET).expect("factor: rho budget exhausted")
}

/// Factorisation of `|n|` (n != 0), or `None` once more than `max_rho_iters` Pollard-rho
/// iterations (summed over all splits) have been spent. `|n| = 1` gives the empty list.
pub fn factor_bounded(n: &Int, max_rho_iters: u64) -> Option<Factors> {
    assert!(!n.is_zero(), "factor of zero");
    let mut m = n.abs();
    let mut primes: Vec<Int> = vec![];
    for q in small_primes() {
        if m.bits() <= 64 {
            break;
        }
        while m.mod_u64(q) == 0 {
            primes.push(Int::from(q));
            m = m.div_floor(&Int::from(q));
        }
    }
    let mut budget = max_rho_iters;
    if m.bits() <= 64 {
        let v = m.to_i128().unwrap() as u64;
        if v > 1 {
            for (q, k) in factor_u64(v) {
                for _ in 0..k {
                    primes.push(Int::from(q));
                }
            }
        }
    } else {
        split_into(&m, &mut primes, &mut budget)?;
    }
    primes.sort();
    let mut out: Factors = vec![];
    for q in primes {
        match out.last_mut() {
            Some((r, k)) if *r == q => *k += 1,
            _ => out.push((q, 1)),
        }
    }
    Some(out)
}

fn split_into(n: &Int, out: &mut Vec<Int>, budget: &mut u64) -> Option<()> {
    if *n == Int::one() {
        return Some(());
    }
    if n.bits() <= 64 {
        for (q, k) in factor_u64(n.to_i128().unwrap() as u64) {
            for _ in 0..k {
                out.push(Int::from(q));
            }
        }
        return Some(());
    }
    if n.is_probable_prime() {
        out.push(n.clone());
        return Some(());
    }
    if n.is_square() {
        let r = n.isqrt();
        split_into(&r, out, budget)?;
        return split_into(&r, out, budget);
    }
    let mut c = 1u64;
    loop {
        if let Some(g) = brent(n, &Int::from(c), budget) {
            split_into(&g, out, budget)?;
            return split_into(&n.div_floor(&g), out, budget);
        }
        if *budget == 0 {
            return None;
        }
        c += 1;
    }
}

/// One run of Brent's variant of Pollard rho with f(y) = y^2 + c: a non-trivial factor, or `None`
/// (failure with this c, or budget exhausted -- the caller distinguishes by `*budget == 0`).
fn brent(n: &Int, c: &Int, budget: &mut u64) -> Option<Int> {
    let f = |y: &Int| (&(y * y) + c).modulo(n);
    let one = Int::one();
    let m = 64u64;
    let mut y = Int::from(2i64);
    let mut x = y.clone();
    let mut ys = y.clone();
    let mut g = one.clone();
    let mut q = one.clone();
    let mut r = 1u64;
    while g == one {
        x = y.clone();
        for _ in 0..r {
            y = f(&y);
        }
        let mut k = 0u64;
        while k < r && g == one {
            ys = y.clone();
            let lim = m.min(r - k);
            if *budget < lim {
                *budget = 0;
                return None;
            }
            *budget -= lim;
            for _ in 0..lim {
                y = f(&y);
                q = (&q * &(&x - &y).abs()).modulo(n);
            }
            g = q.gcd(n);
            k += m;
        }
        r *= 2;
    }
    if g == *n {
        loop {
            if *budget == 0 {
                return None;
            }
            *budget -= 1;
            ys = f(&ys);
            g = (&x - &ys).abs().gcd(n);
            if g != one {
                break;
            }
        }
    }
    if g == *n {
        None
    } else {
        Some(g)
    }
}

/// Number of distinct prime factors of a factorisation.
pub fn omega_of(f: &Factors) -> usize {
    f.len()
}

/// Number of distinct prime factors of `n` (n != 0), with the default budget.
pub fn omega(n: &Int) -> usize {
    factor(n).len()
}

/// The product of a factorisation.
pub fn product(f: &[(Int, u32)]) -> Int {
    let mut r = Int::one();
    for (p, k) in f {
        r = &r * &p.pow(*k);
    }
    r
}

/// All x in [0, 2^m) with x^2 = a (mod 2^m), for odd a.
fn sqrt_unit_mod_2pow(a: &Int, m: u32) -> Vec<Int> {
    let modulus = Int::one().shl(m as usize);
    let a = a.modulo(&modulus);
    match m {
        0 => return vec![Int::zero()],
        1 => return vec![Int::one()],
        2 => {
            return if a.mod_u64(4) == 1 {
                vec![Int::one(), Int::from(3i64)]
            } else {
                vec![]
            }
        }
        _ => {}
    }
    if a.mod_u64(8) != 1 {
        return vec![];
    }
    // invariant: x^2 = a mod 2^j
    let mut x = Int::one();
    for j in 3..m {
        let mj1 = Int::one().shl(j as usize + 1);
        if !(&(&x * &x) - &a).modulo(&mj1).is_zero() {
            x = &x + &Int::one().shl(j as usize - 1);
        }
    }
    let half = Int::one().shl(m as usize - 1);
    let mut v = vec![
        x.modulo(&modulus),
        (-&x).modulo(&modulus),
        (&x + &half).modulo(&modulus),
        (&(-&x) + &half).modulo(&modulus),
    ];
    v.sort();
    v.dedup();
    v
}

/// All x in [0, p^k) with x^2 = a (mod p^k), for a prime to the odd prime p.
fn sqrt_unit_mod_odd_prime_power(a: &Int, p: &Int, k: u32) -> Vec<Int> {
    if k == 0 {
        return vec![Int::zero()];
    }
    let Some(mut r) = Int::sqrt_mod_prime(a, p) else {
        return vec![];
    };
    let mut pk = p.clone();
    let two = Int::from(2i64);
    for _ in 1..k {
        pk = &pk * p;
        // Newton step r <- r - (r^2 - a) / (2r) mod p^{i+1}
        let inv = (&two * &r).inv_mod(&pk).expect("unit root");
        r = (&r - &(&(&(&r * &r) - a) * &inv)).modulo(&pk);
    }
    let mut v = vec![r.modulo(&pk), (-&r).modulo(&pk)];
    v.sort();
    v.dedup();
    v
}

/// All square roots of `a` modulo `p^k` (p prime, k >= 1), in [0, p^k), sorted. Handles p = 2
/// and p | a: if a = p^v a' with v < k, roots exist iff v is even and a' is a square mod
/// p^{k-v}; they are p^{v/2} y + s p^{k - v/2} for the unit roots y mod p^{k-v} and
/// s in [0, p^{v/2}). If p^k | a the roots are the multiples of p^{ceil(k/2)}.
pub fn sqrt_mod_prime_power(a: &Int, p: &Int, k: u32) -> Vec<Int> {
    assert!(k >= 1);
    let pk = p.pow(k);
    let a = a.modulo(&pk);
    if a.is_zero() {
        let step = p.pow(k.div_ceil(2));
        let cnt = p.pow(k / 2);
        let mut out = vec![];
        let mut t = Int::zero();
        while t < cnt {
            out.push(&t * &step);
            t = &t + &Int::one();
        }
        return out;
    }
    let mut v = 0u32;
    let mut a1 = a.clone();
    while a1.modulo(p).is_zero() {
        a1 = a1.div_floor(p);
        v += 1;
    }
    if v % 2 == 1 {
        return vec![];
    }
    let rem = k - v;
    let unit = if *p == Int::from(2i64) {
        sqrt_unit_mod_2pow(&a1, rem)
    } else {
        sqrt_unit_mod_odd_prime_power(&a1, p, rem)
    };
    let base_mul = p.pow(v / 2);
    let step = p.pow(k - v / 2);
    let mut out = vec![];
    for y in unit {
        let base = &base_mul * &y;
        let mut s = Int::zero();
        while s < base_mul {
            out.push((&base + &(&s * &step)).modulo(&pk));
            s = &s + &Int::one();
        }
    }
    out.sort();
    out.dedup();
    out
}

/// All square roots of `a` modulo `n = prod p_i^{k_i}` (given by `factors`, n >= 1), in [0, n),
/// sorted, combined by CRT from [`sqrt_mod_prime_power`]. n = 1 gives `[0]`.
pub fn sqrt_mod_composite(a: &Int, n: &Int, factors: &[(Int, u32)]) -> Vec<Int> {
    debug_assert_eq!(product(factors), *n);
    let mut res = vec![Int::zero()];
    let mut m = Int::one();
    for (p, k) in factors {
        let pk = p.pow(*k);
        let rs = sqrt_mod_prime_power(a, p, *k);
        if rs.is_empty() {
            return vec![];
        }
        let minv = m.inv_mod(&pk).expect("coprime moduli");
        let mut next = Vec::with_capacity(res.len() * rs.len());
        for r1 in &res {
            for r2 in &rs {
                let t = (&(r2 - r1) * &minv).modulo(&pk);
                next.push(r1 + &(&m * &t));
            }
        }
        m = &m * &pk;
        res = next;
    }
    let _ = n;
    res.sort();
    res
}
