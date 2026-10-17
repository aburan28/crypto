//! Exact class number h(D) of a negative fundamental discriminant by counting
//! reduced primitive forms (a, b, c): |b| <= a <= c, b >= 0 when |b| = a or
//! a = c. For fundamental D every form of discriminant D is primitive.
//! Work and memory are O(sqrt(|D|)); intended for |D| < 2^56.

use crate::nt::{isqrt, kronecker};

fn sqrt_mod_prime(d: u64, l: u64) -> Option<u64> {
    // Tonelli-Shanks with u128 products
    let d = d % l;
    if d == 0 {
        return Some(0);
    }
    let pw = |mut a: u64, mut e: u64| {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = (r as u128 * a as u128 % l as u128) as u64;
            }
            a = (a as u128 * a as u128 % l as u128) as u64;
            e >>= 1;
        }
        r
    };
    if l == 2 {
        return Some(d);
    }
    if pw(d, (l - 1) / 2) != 1 {
        return None;
    }
    let (mut q, mut s) = (l - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while pw(z, (l - 1) / 2) == 1 {
        z += 1;
    }
    let mul = |a: u64, b: u64| (a as u128 * b as u128 % l as u128) as u64;
    let (mut m, mut c, mut t, mut r) = (s, pw(z, q), pw(d, q), pw(d, (q + 1) / 2));
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mul(tt, tt);
            i += 1;
        }
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = mul(b, b);
        }
        m = i;
        c = mul(b, b);
        t = mul(t, c);
        r = mul(r, b);
    }
    Some(r)
}

fn crt(rs1: &[u128], m1: u128, rs2: &[u128], m2: u128) -> Vec<u128> {
    // m1, m2 coprime
    let inv = {
        let (mut r0, mut r1) = (m2 as i128, (m1 % m2) as i128);
        let (mut s0, mut s1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        s0.rem_euclid(m2 as i128) as u128
    };
    let mut out = Vec::with_capacity(rs1.len() * rs2.len());
    for &a in rs1 {
        for &b in rs2 {
            let diff = (b + m2 - a % m2) % m2;
            let k = diff * inv % m2;
            out.push(a + m1 * k);
        }
    }
    out
}

/// Exact h(D) for a fundamental discriminant D < 0.
pub fn class_number(d: i128) -> u64 {
    assert!(d < 0);
    let nd = d.unsigned_abs();
    let amax = isqrt(nd / 3) as usize;
    // smallest prime factor sieve and sqrt(D) mod each odd prime (u32::MAX = none)
    let mut spf = vec![0u32; amax + 1];
    for i in 2..=amax {
        if spf[i] == 0 {
            let mut j = i;
            while j <= amax {
                if spf[j] == 0 {
                    spf[j] = i as u32;
                }
                j += i;
            }
        }
    }
    let dmod = |l: u64| -> u64 { d.rem_euclid(l as i128) as u64 };
    let mut sq = vec![u32::MAX; amax + 1];
    for l in 3..=amax {
        if spf[l] as usize == l {
            if let Some(r) = sqrt_mod_prime(dmod(l as u64), l as u64) {
                sq[l] = r as u32;
            }
        }
    }
    // 2-adic part: S[e] = residues b mod 2^(e+1) with b^2 = D mod 2^(e+2)
    let mut s2: Vec<Vec<u128>> = Vec::new();
    {
        let mut cur: Vec<u128> = (0..2u128)
            .filter(|&b| ((b * b) as i128 - d).rem_euclid(4) == 0)
            .collect();
        s2.push(cur.clone());
        for e in 1..70u32 {
            let modb = 1u128 << (e + 1);
            let modq = 1i128 << (e + 2);
            let mut nxt = Vec::new();
            for &b in &cur {
                for cand in [b, b + (1u128 << e)] {
                    let c = cand % modb;
                    if ((c * c) as i128 - d).rem_euclid(modq) == 0 && !nxt.contains(&c) {
                        nxt.push(c);
                    }
                }
            }
            s2.push(nxt.clone());
            cur = nxt;
            if (1u128 << e) > amax as u128 {
                break;
            }
        }
    }
    let mut h: u64 = 0;
    'a: for a in 1..=amax {
        // factor a
        let mut x = a;
        let mut e2 = 0usize;
        while x % 2 == 0 {
            x /= 2;
            e2 += 1;
        }
        let mut roots: Vec<u128> = s2[e2].clone();
        let mut modulus: u128 = 1u128 << (e2 + 1);
        if roots.is_empty() {
            continue;
        }
        while x > 1 {
            let l = spf[x] as usize;
            let mut k = 0;
            while x % l == 0 {
                x /= l;
                k += 1;
            }
            let l64 = l as u64;
            let lk = (l as u128).pow(k);
            let rl: Vec<u128> = if dmod(l64) == 0 {
                if k == 1 {
                    vec![0]
                } else {
                    // l^2 | D impossible for odd l and fundamental D; check anyway
                    (0..lk).filter(|&b| ((b * b) as i128 - d).rem_euclid(lk as i128) == 0).collect()
                }
            } else {
                if sq[l] == u32::MAX {
                    continue 'a;
                }
                let r = sq[l] as u128;
                let mut sols = vec![r, l as u128 - r];
                // Hensel lift to l^k
                let mut mk = l as u128;
                for _ in 1..k {
                    let mk2 = mk * l as u128;
                    sols = sols
                        .iter()
                        .map(|&s| {
                            // s' = s - (s^2 - D) / (2s) mod mk2
                            let f = ((s * s) as i128 - d).rem_euclid(mk2 as i128) as u128;
                            let two_s = (2 * s) % mk2;
                            let inv = {
                                let (mut r0, mut r1) = (mk2 as i128, two_s as i128);
                                let (mut s0, mut s1) = (0i128, 1i128);
                                while r1 != 0 {
                                    let q = r0 / r1;
                                    (r0, r1) = (r1, r0 - q * r1);
                                    (s0, s1) = (s1, s0 - q * s1);
                                }
                                s0.rem_euclid(mk2 as i128) as u128
                            };
                            (s + mk2 - (f % mk2) * inv % mk2) % mk2
                        })
                        .collect();
                    mk = mk2;
                }
                sols.sort_unstable();
                sols.dedup();
                sols
            };
            roots = crt(&roots, modulus, &rl, lk);
            modulus *= lk;
        }
        let a128 = a as u128;
        debug_assert_eq!(modulus, 2 * a128);
        for r in roots {
            let r = r % (2 * a128);
            let b: i128 = if r <= a128 { r as i128 } else { r as i128 - 2 * a as i128 };
            let num = b * b - d;
            debug_assert_eq!(num.rem_euclid(4 * a as i128), 0);
            let c = (num / (4 * a as i128)) as u128;
            if c < a128 {
                continue;
            }
            if (c == a128 || b.unsigned_abs() == a128) && b < 0 {
                continue;
            }
            h += 1;
        }
    }
    h
}

/// Brute-force reference (O(|D|)), for tests.
pub fn class_number_brute(d: i128) -> u64 {
    let nd = d.unsigned_abs();
    let amax = isqrt(nd / 3);
    let mut h = 0;
    for a in 1..=amax as i128 {
        for b in -a + 1..=a {
            let num = b * b - d;
            if num % (4 * a) != 0 {
                continue;
            }
            let c = num / (4 * a);
            if c < a || ((c == a) && b < 0) {
                continue;
            }
            // primitivity
            let g = gcd3(a, b.abs(), c);
            if g == 1 {
                h += 1;
            }
        }
    }
    h
}

fn gcd3(a: i128, b: i128, c: i128) -> i128 {
    fn g(mut x: i128, mut y: i128) -> i128 {
        while y != 0 {
            let t = x % y;
            x = y;
            y = t;
        }
        x.abs()
    }
    g(g(a, b), c)
}

/// True when d is a fundamental discriminant.
pub fn is_fundamental(d: i128) -> bool {
    let m4 = d.rem_euclid(4);
    let sqfree = |n: u128| crate::nt::factor(n).iter().all(|&(_, e)| e == 1);
    if m4 == 1 {
        return sqfree(d.unsigned_abs());
    }
    if m4 == 0 {
        let e = d / 4;
        let em = e.rem_euclid(4);
        return (em == 2 || em == 3) && sqfree(e.unsigned_abs());
    }
    let _ = kronecker;
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values() {
        for (d, h) in [(-3, 1), (-4, 1), (-7, 1), (-8, 1), (-23, 3), (-47, 5), (-71, 7), (-199, 9), (-20, 2), (-56, 4), (-163, 1)] {
            assert_eq!(class_number(d), h, "D={d}");
        }
    }

    #[test]
    fn matches_brute_force() {
        let mut checked = 0;
        for n in 3..4000i128 {
            let d = -n;
            if is_fundamental(d) {
                assert_eq!(class_number(d), class_number_brute(d), "D={d}");
                checked += 1;
            }
        }
        assert!(checked > 1000);
        for d in [-5923i128, -10007, -99991, -65539] {
            if is_fundamental(d) {
                assert_eq!(class_number(d), class_number_brute(d), "D={d}");
            }
        }
    }
}
