//! Dense univariate polynomials over [`Field`]: the arithmetic the root
//! finder, the kernel construction and the verifier share.
//!
//! A polynomial is a `Vec<Fe>` of coefficients, lowest degree first, with
//! no trailing zeros (the zero polynomial is empty).

use num_bigint::BigUint;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use super::field::{Fe, Field};

pub type Poly = Vec<Fe>;

pub fn trim(f: &Field, mut a: Poly) -> Poly {
    while a.last().is_some_and(|c| f.is_zero(c)) {
        a.pop();
    }
    a
}

/// Degree, with `None` for the zero polynomial.
pub fn degree(a: &Poly) -> Option<usize> {
    a.len().checked_sub(1)
}

pub fn constant(f: &Field, c: Fe) -> Poly {
    trim(f, vec![c])
}

/// `X`.
pub fn x(f: &Field) -> Poly {
    vec![f.zero(), f.one()]
}

pub fn add(f: &Field, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let z = f.zero();
    let out = (0..n)
        .map(|i| f.add(a.get(i).unwrap_or(&z), b.get(i).unwrap_or(&z)))
        .collect();
    trim(f, out)
}

pub fn sub(f: &Field, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let z = f.zero();
    let out = (0..n)
        .map(|i| f.sub(a.get(i).unwrap_or(&z), b.get(i).unwrap_or(&z)))
        .collect();
    trim(f, out)
}

pub fn scale(f: &Field, a: &Poly, c: &Fe) -> Poly {
    trim(f, a.iter().map(|v| f.mul(v, c)).collect())
}

pub fn mul(f: &Field, a: &Poly, b: &Poly) -> Poly {
    trim(f, fast_product(f, a, b))
}

fn fast_product(f: &Field, a: &[Fe], b: &[Fe]) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    if b.len() >= 48 {
        if a.len() >= 2 * b.len() {
            let mut out = vec![f.zero(); a.len() + b.len() - 1];
            for (chunk_index, chunk) in a.chunks(b.len()).enumerate() {
                let p = fast_product(f, chunk, b);
                let offset = chunk_index * b.len();
                for (i, c) in p.iter().enumerate() {
                    out[offset + i] = f.add(&out[offset + i], c);
                }
            }
            return out;
        }
        let middle = a.len().div_ceil(2);
        let (a0, a1) = a.split_at(middle);
        let (b0, b1) = b.split_at(middle.min(b.len()));
        let z0 = fast_product(f, a0, b0);
        let z2 = fast_product(f, a1, b1);
        let sum_a: Poly = (0..middle)
            .map(|i| f.add(&a0[i], a1.get(i).unwrap_or(&f.zero())))
            .collect();
        let sum_b: Poly = (0..middle)
            .map(|i| {
                f.add(
                    b0.get(i).unwrap_or(&f.zero()),
                    b1.get(i).unwrap_or(&f.zero()),
                )
            })
            .collect();
        let mut z1 = fast_product(f, &sum_a, &sum_b);
        for (i, c) in z0.iter().enumerate() {
            z1[i] = f.sub(&z1[i], c);
        }
        for (i, c) in z2.iter().enumerate() {
            z1[i] = f.sub(&z1[i], c);
        }
        let mut out = vec![f.zero(); a.len() + b.len() - 1];
        for (offset, p) in [(0, z0), (middle, z1), (2 * middle, z2)] {
            for (i, c) in p.iter().enumerate() {
                if offset + i < out.len() {
                    out[offset + i] = f.add(&out[offset + i], c);
                }
            }
        }
        return out;
    }
    let mut out = vec![f.zero(); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        if f.is_zero(ai) {
            continue;
        }
        for (j, bj) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(ai, bj));
        }
    }
    out
}

/// `(q, r)` with `a = q·b + r`, `deg r < deg b`.  Panics on `b = 0`.
pub fn divrem(f: &Field, a: &Poly, b: &Poly) -> (Poly, Poly) {
    let db = degree(b).expect("division by the zero polynomial");
    let inv_lead = f.inv(&b[db]).expect("trimmed polynomial has a unit lead");
    let mut r = a.clone();
    if r.len() <= db {
        return (Vec::new(), r);
    }
    let mut q = vec![f.zero(); r.len() - db];
    for k in (0..q.len()).rev() {
        let c = f.mul(&r[k + db], &inv_lead);
        q[k] = c;
        if f.is_zero(&c) {
            continue;
        }
        for (j, bj) in b.iter().enumerate() {
            r[k + j] = f.sub(&r[k + j], &f.mul(&c, bj));
        }
    }
    r.truncate(db);
    (trim(f, q), trim(f, r))
}

pub fn rem(f: &Field, a: &Poly, b: &Poly) -> Poly {
    divrem(f, a, b).1
}

pub fn mulmod(f: &Field, a: &Poly, b: &Poly, m: &Poly) -> Poly {
    rem(f, &mul(f, a, b), m)
}

pub fn monic(f: &Field, a: &Poly) -> Poly {
    match a.last() {
        None => Vec::new(),
        Some(lead) => scale(f, a, &f.inv(lead).expect("nonzero lead")),
    }
}

/// The monic gcd.
pub fn gcd(f: &Field, a: &Poly, b: &Poly) -> Poly {
    let (mut a, mut b) = (a.clone(), b.clone());
    while !b.is_empty() {
        let r = rem(f, &a, &b);
        a = b;
        b = r;
    }
    monic(f, &a)
}

/// `a^{-1} mod m`, or `None` when `gcd(a, m) ≠ 1`.
pub fn invmod(f: &Field, a: &Poly, m: &Poly) -> Option<Poly> {
    let (mut r0, mut r1) = (m.clone(), rem(f, a, m));
    let (mut s0, mut s1): (Poly, Poly) = (Vec::new(), vec![f.one()]);
    while !r1.is_empty() {
        let (q, r) = divrem(f, &r0, &r1);
        let s = sub(f, &s0, &mul(f, &q, &s1));
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s);
    }
    if degree(&r0) != Some(0) {
        return None;
    }
    let c = f.inv(&r0[0])?;
    Some(rem(f, &scale(f, &s0, &c), m))
}

pub fn powmod(f: &Field, base: &Poly, e: &BigUint, m: &Poly) -> Poly {
    let base = rem(f, base, m);
    let mut acc = rem(f, &vec![f.one()], m);
    for i in (0..e.bits()).rev() {
        acc = mulmod(f, &acc, &acc, m);
        if e.bit(i) {
            acc = mulmod(f, &acc, &base, m);
        }
    }
    acc
}

pub fn derivative(f: &Field, a: &Poly) -> Poly {
    trim(
        f,
        a.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| f.mul(c, &f.from_u64(i as u64)))
            .collect(),
    )
}

pub fn eval(f: &Field, a: &Poly, x: &Fe) -> Fe {
    a.iter()
        .rev()
        .fold(f.zero(), |acc, c| f.add(&f.mul(&acc, x), c))
}

/// Evaluate `a` at an element `x` of `F_p[X]/(m)`.
pub fn eval_mod(f: &Field, a: &Poly, x: &Poly, m: &Poly) -> Poly {
    a.iter().rev().fold(Vec::new(), |acc, c| {
        add(f, &mulmod(f, &acc, x, m), &constant(f, *c))
    })
}

/// The distinct roots of `a` in `F_p`, sorted by canonical integer value.
/// `seed` drives the Cantor–Zassenhaus splits; the result does not depend
/// on it.
pub fn roots(f: &Field, a: &Poly, seed: u64) -> Vec<Fe> {
    let a = monic(f, a);
    if degree(&a).unwrap_or(0) == 0 {
        return Vec::new();
    }
    let xp = powmod(f, &x(f), f.modulus(), &a);
    let g = gcd(f, &a, &sub(f, &xp, &x(f)));
    let mut out = Vec::new();
    let mut rng = SmallRng::seed_from_u64(seed);
    split_linear(f, g, &mut rng, &mut out);
    out.sort_by_key(|r| f.to_big(r));
    out
}

/// Split a squarefree product of distinct linear factors.
fn split_linear(f: &Field, g: Poly, rng: &mut SmallRng, out: &mut Vec<Fe>) {
    match degree(&g) {
        None | Some(0) => {}
        Some(1) => out.push(f.neg(&f.mul(&g[0], &f.inv(&g[1]).unwrap()))),
        Some(_) => {
            let half = (f.modulus() - 1u8) >> 1;
            loop {
                let d = f.from_u64(rng.gen());
                let t = vec![d, f.one()];
                let h = sub(f, &powmod(f, &t, &half, &g), &vec![f.one()]);
                let s = gcd(f, &g, &h);
                let ds = degree(&s).unwrap_or(0);
                if ds > 0 && ds < degree(&g).unwrap() {
                    let other = divrem(f, &g, &s).0;
                    split_linear(f, s, rng, out);
                    split_linear(f, monic(f, &other), rng, out);
                    return;
                }
            }
        }
    }
}

/// Power sums `Σ r_i^k`, `k = 1..=n`, of the roots of the monic `h`, by
/// Newton's identities.
pub fn power_sums(f: &Field, h: &Poly, n: usize) -> Vec<Fe> {
    let d = degree(h).unwrap_or(0);
    // e_k with h = Σ (-1)^k e_k X^{d-k}.
    let e: Vec<Fe> = (0..=d)
        .map(|k| {
            let c = h[d - k];
            if k % 2 == 1 {
                f.neg(&c)
            } else {
                c
            }
        })
        .collect();
    let mut p = vec![f.zero(); n + 1];
    for k in 1..=n {
        // p_k = Σ_{i=1}^{k-1} (-1)^{i-1} e_i p_{k-i} + (-1)^{k-1} k e_k.
        let mut acc = f.zero();
        for i in 1..k {
            if i > d {
                break;
            }
            let t = f.mul(&e[i], &p[k - i]);
            acc = if i % 2 == 1 {
                f.add(&acc, &t)
            } else {
                f.sub(&acc, &t)
            };
        }
        if k <= d {
            let t = f.mul(&e[k], &f.from_u64(k as u64));
            acc = if k % 2 == 1 {
                f.add(&acc, &t)
            } else {
                f.sub(&acc, &t)
            };
        }
        p[k] = acc;
    }
    p.remove(0);
    p
}

/// The monic polynomial of degree `d` whose root power sums are
/// `s[0..d]` (`s[k-1] = Σ r^k`).  Needs `p > d`.
pub fn from_power_sums(f: &Field, s: &[Fe], d: usize) -> Poly {
    let mut e = vec![f.one()];
    for k in 1..=d {
        let mut acc = f.zero();
        for i in 1..=k {
            let t = f.mul(&e[k - i], &s[i - 1]);
            acc = if i % 2 == 1 {
                f.add(&acc, &t)
            } else {
                f.sub(&acc, &t)
            };
        }
        e.push(f.div(&acc, &f.from_u64(k as u64)).expect("p > d"));
    }
    let mut h = vec![f.zero(); d + 1];
    for (k, ek) in e.iter().enumerate() {
        h[d - k] = if k % 2 == 1 { f.neg(ek) } else { *ek };
    }
    h
}
