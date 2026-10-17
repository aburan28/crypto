//! Integer arithmetic for the selector: u128 modular arithmetic, deterministic
//! Miller-Rabin, Pollard-Brent factorisation, Kronecker symbols and the
//! decomposition of a Frobenius discriminant t^2 - 4q = f^2 d_K.

pub fn mulmod(a: u128, b: u128, n: u128) -> u128 {
    if n <= u64::MAX as u128 {
        return (a % n) * (b % n) % n;
    }
    debug_assert!(n < (1u128 << 127));
    let (mut a, mut b, mut r) = (a % n, b % n, 0u128);
    while b > 0 {
        if b & 1 == 1 {
            r += a;
            if r >= n {
                r -= n;
            }
        }
        a <<= 1;
        if a >= n {
            a -= n;
        }
        b >>= 1;
    }
    r
}

pub fn powmod(mut a: u128, mut e: u128, n: u128) -> u128 {
    let mut r = 1 % n;
    a %= n;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, n);
        }
        a = mulmod(a, a, n);
        e >>= 1;
    }
    r
}

const SMALL_PRIMES: [u128; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];

/// Deterministic for n < 3.3e24 (bases 2..41); every input here is < 2^72.
pub fn is_prime(n: u128) -> bool {
    if n < 2 {
        return false;
    }
    for p in SMALL_PRIMES {
        if n == p {
            return true;
        }
        if n % p == 0 {
            return false;
        }
    }
    assert!(n < 3_317_044_064_679_887_385_961_981u128, "MR bases not proved for n");
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'base: for a in SMALL_PRIMES {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn brent(n: u128, c: u128) -> Option<u128> {
    let f = |x: u128| (mulmod(x, x, n) + c) % n;
    let (mut y, mut r, mut q, m) = (2u128, 1u64, 1u128, 128u64);
    let (mut g, mut x, mut ys) = (1u128, 0u128, 0u128);
    while g == 1 {
        x = y;
        for _ in 0..r {
            y = f(y);
        }
        let mut k = 0;
        while k < r && g == 1 {
            ys = y;
            for _ in 0..m.min(r - k) {
                y = f(y);
                q = mulmod(q, x.abs_diff(y), n);
            }
            g = gcd(q, n);
            k += m;
        }
        r *= 2;
        if r > 1 << 40 {
            return None;
        }
    }
    if g == n {
        loop {
            ys = f(ys);
            g = gcd(x.abs_diff(ys), n);
            if g > 1 {
                break;
            }
        }
    }
    if g == n {
        None
    } else {
        Some(g)
    }
}

fn split(n: u128, out: &mut Vec<u128>) {
    if n == 1 {
        return;
    }
    if is_prime(n) {
        out.push(n);
        return;
    }
    let r = isqrt(n);
    if r * r == n {
        split(r, out);
        split(r, out);
        return;
    }
    for c in 1..200u128 {
        if let Some(d) = brent(n, c) {
            split(d, out);
            split(n / d, out);
            return;
        }
    }
    panic!("factorisation failed for {n}");
}

/// Prime factorisation as sorted (prime, exponent) pairs.
pub fn factor(mut n: u128) -> Vec<(u128, u32)> {
    assert!(n > 0);
    let mut ps = Vec::new();
    let mut p = 2u128;
    while p < 10_000 && p * p <= n {
        while n % p == 0 {
            ps.push(p);
            n /= p;
        }
        p += if p == 2 { 1 } else { 2 };
    }
    split(n, &mut ps);
    ps.sort_unstable();
    let mut out: Vec<(u128, u32)> = Vec::new();
    for p in ps {
        match out.last_mut() {
            Some((q, e)) if *q == p => *e += 1,
            _ => out.push((p, 1)),
        }
    }
    out
}

pub fn isqrt(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u128;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

/// Kronecker symbol (d / n) for n > 0.
pub fn kronecker(d: i128, n: u128) -> i32 {
    let mut n = n;
    let mut res = 1i32;
    let mut a = d;
    while n % 2 == 0 {
        n /= 2;
        let am8 = a.rem_euclid(8);
        if am8 % 2 == 0 {
            return 0;
        }
        if am8 == 3 || am8 == 5 {
            res = -res;
        }
    }
    // Jacobi (a / n), n odd
    let mut a = a.rem_euclid(n as i128) as u128;
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                res = -res;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            res = -res;
        }
        a %= n;
    }
    if n == 1 {
        res
    } else {
        0
    }
}

/// Decomposition t^2 - 4q = f^2 * d_K with d_K fundamental.
#[derive(Clone, Debug)]
pub struct FrobDisc {
    pub t: i128,
    pub delta: i128,
    pub d_k: i128,
    pub f: u128,
    pub f_factors: Vec<(u128, u32)>,
}

pub fn frob_disc(t: i128, q: u128) -> FrobDisc {
    let delta = t * t - 4 * q as i128;
    assert!(delta < 0, "not ordinary-imaginary");
    let fac = factor(delta.unsigned_abs());
    let mut f: u128 = 1;
    let mut core: u128 = 1;
    for &(p, e) in &fac {
        f *= p.pow(e / 2);
        if e % 2 == 1 {
            core *= p;
        }
    }
    let mut d0 = -(core as i128);
    if d0.rem_euclid(4) != 1 {
        assert!(f % 2 == 0, "discriminant bookkeeping");
        f /= 2;
        d0 *= 4;
    }
    assert_eq!(d0 * (f * f) as i128, delta);
    let f_factors = if f > 1 { factor(f) } else { vec![] };
    FrobDisc { t, delta, d_k: d0, f, f_factors }
}

/// h(c^2 d_K) / h(d_K) for an imaginary quadratic order of conductor c.
pub fn order_class_ratio(d_k: i128, c_factors: &[(u128, u32)]) -> (u128, u128) {
    // h(O_c) = h(O_K) * c * prod_{l|c} (1 - (d_K/l)/l) / [O_K^* : O_c^*]
    let mut num: u128 = 1;
    let mut den: u128 = 1;
    let mut c: u128 = 1;
    for &(l, e) in c_factors {
        c *= l.pow(e);
        let chi = kronecker(d_k, l) as i128;
        num *= (l as i128 - chi) as u128;
        den *= l;
    }
    let mut num = c * num;
    if c > 1 {
        if d_k == -3 {
            den *= 3;
        } else if d_k == -4 {
            den *= 2;
        }
    }
    let g = gcd(num, den);
    num /= g;
    den /= g;
    (num, den)
}

/// Every divisor c of f (as factorisations), in increasing order of c.
pub fn divisors(f_factors: &[(u128, u32)]) -> Vec<(u128, Vec<(u128, u32)>)> {
    let mut out: Vec<(u128, Vec<(u128, u32)>)> = vec![(1, vec![])];
    for &(p, e) in f_factors {
        let mut next = Vec::new();
        for (d, fac) in &out {
            let mut pk = 1;
            for k in 0..=e {
                let mut fk = fac.clone();
                if k > 0 {
                    fk.push((p, k));
                }
                next.push((d * pk, fk));
                pk *= p;
            }
        }
        out = next;
    }
    out.sort_by_key(|x| x.0);
    out
}

/// Euler-product estimate of h(D), D < -4 fundamental (bound B).
pub fn class_number_estimate(d: i128, bound: u64) -> f64 {
    let mut l = 1.0f64;
    let sieve = primes_upto(bound);
    for p in sieve {
        let chi = kronecker(d, p as u128) as f64;
        l *= 1.0 / (1.0 - chi / p as f64);
    }
    (d.unsigned_abs() as f64).sqrt() / std::f64::consts::PI * l
}

pub fn primes_upto(n: u64) -> Vec<u64> {
    let n = n as usize;
    let mut is = vec![true; n + 1];
    is[0] = false;
    if n >= 1 {
        is[1] = false;
    }
    let mut i = 2;
    while i * i <= n {
        if is[i] {
            let mut j = i * i;
            while j <= n {
                is[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    (0..=n).filter(|&i| is[i]).map(|i| i as u64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_and_primes() {
        assert!(is_prime(2u128.pow(61) - 1));
        assert!(!is_prime((2u128.pow(31) - 1) * 4_294_967_291));
        let n = (2u128.pow(31) - 1) * 4_294_967_291 * 9;
        assert_eq!(factor(n), vec![(3, 2), (2u128.pow(31) - 1, 1), (4_294_967_291, 1)]);
        assert_eq!(factor(1_000_000_007 * 998_244_353), vec![(998_244_353, 1), (1_000_000_007, 1)]);
    }

    #[test]
    fn kronecker_values() {
        assert_eq!(kronecker(-7, 2), 1);
        assert_eq!(kronecker(-3, 2), -1);
        assert_eq!(kronecker(-4, 3), -1);
        assert_eq!(kronecker(-7, 11), 1);
        assert_eq!(kronecker(-7, 7), 0);
    }

    #[test]
    fn disc_decomposition() {
        // t = 1, q = 2: 1 - 8 = -7
        let d = frob_disc(1, 2);
        assert_eq!((d.d_k, d.f), (-7, 1));
        // -7 * 9 * 4 = -252 = t^2 - 4q with t=2? 4 - 4q = -252 -> q = 64
        let d = frob_disc(2, 64);
        assert_eq!((d.d_k, d.f), (-7, 6));
        // h(-7 * 9) = 3 - (-7/3) = 3 - (-1) = 4 -> wait chi(-7,3) = (2/3) = -1, so h = 4
        assert_eq!(order_class_ratio(-7, &[(3, 1)]), (4, 1));
        // h(-3 * 49) = 7*(1 - (-3/7)/7)/3 = (7 - 1)/3 = 2
        assert_eq!(order_class_ratio(-3, &[(7, 1)]), (2, 1));
        assert_eq!(divisors(&[(3, 2), (5, 1)]).iter().map(|d| d.0).collect::<Vec<_>>(), vec![1, 3, 5, 9, 15, 45]);
    }
}
