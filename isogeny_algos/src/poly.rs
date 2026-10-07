//! Dense univariate polynomials over a generic field, schoolbook arithmetic.
//! Coefficients low-to-high; the zero polynomial is the empty vec.
use crate::field::{Field, Rng};

pub type Poly<F> = Vec<<F as Field>::E>;

pub fn trim<F: Field>(f: &F, p: &mut Poly<F>) {
    while let Some(&l) = p.last() {
        if f.is_zero(l) {
            p.pop();
        } else {
            break;
        }
    }
}
pub fn deg<F: Field>(_f: &F, p: &Poly<F>) -> isize {
    p.len() as isize - 1
}
pub fn constant<F: Field>(f: &F, c: F::E) -> Poly<F> {
    let mut v = vec![c];
    trim(f, &mut v);
    v
}
pub fn x_poly<F: Field>(f: &F) -> Poly<F> {
    vec![f.zero(), f.one()]
}
pub fn add<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    let n = a.len().max(b.len());
    let mut r = Vec::with_capacity(n);
    for i in 0..n {
        let x = if i < a.len() { a[i] } else { f.zero() };
        let y = if i < b.len() { b[i] } else { f.zero() };
        r.push(f.add(x, y));
    }
    trim(f, &mut r);
    r
}
pub fn sub<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    let n = a.len().max(b.len());
    let mut r = Vec::with_capacity(n);
    for i in 0..n {
        let x = if i < a.len() { a[i] } else { f.zero() };
        let y = if i < b.len() { b[i] } else { f.zero() };
        r.push(f.sub(x, y));
    }
    trim(f, &mut r);
    r
}
pub fn scale<F: Field>(f: &F, a: &Poly<F>, c: F::E) -> Poly<F> {
    let mut r: Poly<F> = a.iter().map(|&x| f.mul(x, c)).collect();
    trim(f, &mut r);
    r
}
pub fn mul<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![f.zero(); a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if f.is_zero(x) {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = f.add(r[i + j], f.mul(x, y));
        }
    }
    trim(f, &mut r);
    r
}
pub fn derivative<F: Field>(f: &F, a: &Poly<F>) -> Poly<F> {
    let mut r: Poly<F> = (1..a.len()).map(|i| f.mul(a[i], f.from_u64(i as u64))).collect();
    trim(f, &mut r);
    r
}
pub fn eval<F: Field>(f: &F, a: &Poly<F>, x: F::E) -> F::E {
    let mut r = f.zero();
    for &c in a.iter().rev() {
        r = f.add(f.mul(r, x), c);
    }
    r
}
pub fn monic<F: Field>(f: &F, a: &Poly<F>) -> Poly<F> {
    if a.is_empty() {
        return vec![];
    }
    let li = f.inv(*a.last().unwrap());
    scale(f, a, li)
}
pub fn divrem<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> (Poly<F>, Poly<F>) {
    assert!(!b.is_empty(), "division by zero polynomial");
    if a.len() < b.len() {
        return (vec![], a.clone());
    }
    let li = f.inv(*b.last().unwrap());
    let mut r = a.clone();
    let mut q = vec![f.zero(); a.len() - b.len() + 1];
    for i in (0..q.len()).rev() {
        let c = f.mul(r[i + b.len() - 1], li);
        q[i] = c;
        if !f.is_zero(c) {
            for j in 0..b.len() {
                r[i + j] = f.sub(r[i + j], f.mul(c, b[j]));
            }
        }
    }
    r.truncate(b.len() - 1);
    trim(f, &mut r);
    trim(f, &mut q);
    (q, r)
}
pub fn rem<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    divrem(f, a, b).1
}
pub fn gcd<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    let (mut a, mut b) = (a.clone(), b.clone());
    while !b.is_empty() {
        let r = rem(f, &a, &b);
        a = b;
        b = r;
    }
    monic(f, &a)
}
pub fn mulmod<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>, m: &Poly<F>) -> Poly<F> {
    rem(f, &mul(f, a, b), m)
}
pub fn powmod<F: Field>(f: &F, a: &Poly<F>, mut e: u128, m: &Poly<F>) -> Poly<F> {
    let mut r = constant(f, f.one());
    let mut b = rem(f, a, m);
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(f, &r, &b, m);
        }
        e >>= 1;
        if e > 0 {
            b = mulmod(f, &b, &b, m);
        }
    }
    r
}
/// Product of (x - r) over the given roots.
pub fn from_roots<F: Field>(f: &F, roots: &[F::E]) -> Poly<F> {
    let mut p = vec![f.one()];
    for &r in roots {
        let mut n = vec![f.zero(); p.len() + 1];
        for (i, &c) in p.iter().enumerate() {
            n[i + 1] = f.add(n[i + 1], c);
            n[i] = f.sub(n[i], f.mul(c, r));
        }
        p = n;
    }
    p
}

fn split_roots<F: Field>(f: &F, g: &Poly<F>, rng: &mut Rng, out: &mut Vec<F::E>) {
    let d = deg(f, g);
    if d <= 0 {
        return;
    }
    if d == 1 {
        out.push(f.neg(f.div(g[0], g[1])));
        return;
    }
    let half = (f.size() - 1) / 2;
    loop {
        let r = f.random(rng);
        let base = vec![r, f.one()];
        let w = sub(f, &powmod(f, &base, half, g), &constant(f, f.one()));
        let dd = gcd(f, g, &w);
        let dd_deg = deg(f, &dd);
        if dd_deg > 0 && dd_deg < d {
            let (q, _) = divrem(f, g, &dd);
            split_roots(f, &dd, rng, out);
            split_roots(f, &monic(f, &q), rng, out);
            return;
        }
    }
}

/// Distinct roots in the field (Cantor-Zassenhaus).
pub fn roots<F: Field>(f: &F, p: &Poly<F>, rng: &mut Rng) -> Vec<F::E> {
    let p = monic(f, p);
    if deg(f, &p) <= 0 {
        return vec![];
    }
    let xq = powmod(f, &x_poly(f), f.size(), &p);
    let g = gcd(f, &p, &sub(f, &xq, &x_poly(f)));
    let mut out = vec![];
    split_roots(f, &g, rng, &mut out);
    out
}

/// Distinct-degree factorisation of a squarefree polynomial.
pub fn ddf<F: Field>(f: &F, p: &Poly<F>) -> Vec<(usize, Poly<F>)> {
    let mut out = vec![];
    let mut fp = monic(f, p);
    let mut h = x_poly(f);
    let mut k = 1usize;
    while deg(f, &fp) >= 2 * k as isize {
        h = powmod(f, &h, f.size(), &fp);
        let g = gcd(f, &fp, &sub(f, &h, &x_poly(f)));
        if deg(f, &g) > 0 {
            fp = monic(f, &divrem(f, &fp, &g).0);
            h = rem(f, &h, &fp);
            out.push((k, g));
        }
        k += 1;
    }
    if deg(f, &fp) > 0 {
        out.push((deg(f, &fp) as usize, fp));
    }
    out
}

/// Equal-degree factorisation: split a product of degree-k irreducibles.
pub fn edf<F: Field>(f: &F, g: &Poly<F>, k: usize, rng: &mut Rng, out: &mut Vec<Poly<F>>) {
    let d = deg(f, g);
    if d as usize == k {
        out.push(monic(f, g));
        return;
    }
    let half = (f.size() - 1) / 2;
    loop {
        let a: Poly<F> = {
            let mut a: Poly<F> = (0..d as usize).map(|_| f.random(rng)).collect();
            trim(f, &mut a);
            a
        };
        if a.is_empty() {
            continue;
        }
        // a^{(q^k-1)/2} = (a^{1+q+..+q^{k-1}})^{(q-1)/2}
        let mut t = rem(f, &a, g);
        let mut s = t.clone();
        for _ in 1..k {
            t = powmod(f, &t, f.size(), g);
            s = mulmod(f, &s, &t, g);
        }
        let b = powmod(f, &s, half, g);
        let dd = gcd(f, g, &sub(f, &b, &constant(f, f.one())));
        let dd_deg = deg(f, &dd);
        if dd_deg > 0 && dd_deg < d {
            let (q, _) = divrem(f, g, &dd);
            edf(f, &dd, k, rng, out);
            edf(f, &monic(f, &q), k, rng, out);
            return;
        }
    }
}

/// Full factorisation of a squarefree polynomial into monic irreducibles.
pub fn factor_squarefree<F: Field>(f: &F, p: &Poly<F>, rng: &mut Rng) -> Vec<Poly<F>> {
    let mut out = vec![];
    for (k, g) in ddf(f, p) {
        edf(f, &g, k, rng, &mut out);
    }
    out
}

/// Inverse of `a` modulo `m` (extended Euclid); None if not coprime.
pub fn invmod<F: Field>(f: &F, a: &Poly<F>, m: &Poly<F>) -> Option<Poly<F>> {
    let (mut r0, mut r1) = (m.clone(), rem(f, a, m));
    let (mut t0, mut t1): (Poly<F>, Poly<F>) = (vec![], vec![f.one()]);
    while !r1.is_empty() {
        let (q, r) = divrem(f, &r0, &r1);
        let t2 = sub(f, &t0, &mul(f, &q, &t1));
        r0 = r1;
        r1 = r;
        t0 = t1;
        t1 = t2;
    }
    if r0.len() != 1 {
        return None;
    }
    Some(rem(f, &scale(f, &t0, f.inv(r0[0])), m))
}
