//! Dense univariate polynomials over a generic field, schoolbook arithmetic.
//! Coefficients low-to-high; the zero polynomial is the empty vec.
use crate::bigint::Big;
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
/// Karatsuba threshold (schoolbook below this many coefficients).
pub const KARATSUBA_THRESHOLD: usize = 32;

fn school<F: Field>(f: &F, a: &[F::E], b: &[F::E]) -> Vec<F::E> {
    f.conv(a, b)
}

/// Untrimmed product of two non-empty coefficient slices, length a.len() + b.len() - 1.
pub fn mul_raw<F: Field>(f: &F, a: &[F::E], b: &[F::E]) -> Vec<F::E> {
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    if b.len() < KARATSUBA_THRESHOLD {
        return school(f, a, b);
    }
    if a.len() >= 2 * b.len() {
        // unbalanced: chunk the longer operand
        let mut r = vec![f.zero(); a.len() + b.len() - 1];
        for (ci, chunk) in a.chunks(b.len()).enumerate() {
            let pr = mul_raw(f, chunk, b);
            let off = ci * b.len();
            for (k, &v) in pr.iter().enumerate() {
                r[off + k] = f.add(r[off + k], v);
            }
        }
        return r;
    }
    let h = a.len().div_ceil(2);
    let (a0, a1) = a.split_at(h.min(a.len()));
    let (b0, b1) = b.split_at(h.min(b.len()));
    let z0 = mul_raw(f, a0, b0);
    let mut r = vec![f.zero(); a.len() + b.len() - 1];
    for (k, &v) in z0.iter().enumerate() {
        r[k] = f.add(r[k], v);
    }
    if b1.is_empty() {
        // b fits in the low half: a1 * b
        let z = mul_raw(f, a1, b);
        for (k, &v) in z.iter().enumerate() {
            r[h + k] = f.add(r[h + k], v);
        }
        return r;
    }
    let z2 = mul_raw(f, a1, b1);
    let sa: Vec<F::E> = (0..h)
        .map(|i| f.add(a0[i], if i < a1.len() { a1[i] } else { f.zero() }))
        .collect();
    let sb: Vec<F::E> = (0..h)
        .map(|i| {
            f.add(
                if i < b0.len() { b0[i] } else { f.zero() },
                if i < b1.len() { b1[i] } else { f.zero() },
            )
        })
        .collect();
    let mut z1 = mul_raw(f, &sa, &sb);
    for (k, &v) in z0.iter().enumerate() {
        z1[k] = f.sub(z1[k], v);
    }
    for (k, &v) in z2.iter().enumerate() {
        z1[k] = f.sub(z1[k], v);
    }
    for (k, &v) in z1.iter().enumerate() {
        if h + k < r.len() {
            r[h + k] = f.add(r[h + k], v);
        }
    }
    for (k, &v) in z2.iter().enumerate() {
        r[2 * h + k] = f.add(r[2 * h + k], v);
    }
    r
}

pub fn mul<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = mul_raw(f, a, b);
    trim(f, &mut r);
    r
}

/// Schoolbook product (kept for reference / benchmarks).
pub fn mul_schoolbook<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = school(f, a, b);
    trim(f, &mut r);
    r
}

/// Reducer for a fixed modulus: remainder by two multiplications with a precomputed
/// inverse of the reversed modulus (Newton/Barrett style division).
pub struct PolyModulus<F: Field> {
    pub m: Poly<F>,
    n: usize,
    rinv: Vec<F::E>,
}

impl<F: Field> PolyModulus<F> {
    pub fn new(f: &F, m: &Poly<F>) -> Self {
        let m = monic(f, m);
        let n = m.len() - 1;
        let rev: Vec<F::E> = m.iter().rev().copied().collect();
        let rinv = if n >= 2 {
            crate::series::inv(f, &rev, n - 1)
        } else {
            vec![f.one()]
        };
        PolyModulus { m, n, rinv }
    }
    pub fn rem(&self, f: &F, a: &Poly<F>) -> Poly<F> {
        let n = self.n;
        if a.len() <= n {
            let mut r = a.clone();
            trim(f, &mut r);
            return r;
        }
        if a.len() > 2 * n - 1 || n < 2 {
            // outside the reducer's range: plain division (never back through `rem`, which may
            // dispatch here again)
            return divrem(f, a, &self.m).1;
        }
        let k = a.len() - n; // number of quotient coefficients (<= n-1)
        let rev_a: Vec<F::E> = a.iter().rev().take(k).copied().collect();
        let mut qr = mul_raw(f, &rev_a, &self.rinv[..k]);
        qr.truncate(k);
        let q: Vec<F::E> = qr.into_iter().rev().collect();
        let qm = mul_raw(f, &q, &self.m);
        let mut r: Vec<F::E> = (0..n).map(|i| f.sub(a[i], qm[i])).collect();
        trim(f, &mut r);
        r
    }
    pub fn mulmod(&self, f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        self.rem(f, &mul(f, a, b))
    }
    pub fn powmod_big(&self, f: &F, a: &Poly<F>, e: &Big) -> Poly<F> {
        let mut r = constant(f, f.one());
        let b = self.rem(f, &rem(f, a, &self.m));
        for i in (0..e.bits()).rev() {
            r = self.mulmod(f, &r, &r);
            if e.bit(i) {
                r = self.mulmod(f, &r, &b);
            }
        }
        r
    }
}

pub fn derivative<F: Field>(f: &F, a: &Poly<F>) -> Poly<F> {
    let mut r: Poly<F> = (1..a.len())
        .map(|i| f.mul(a[i], f.from_u64(i as u64)))
        .collect();
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
/// (a(x), a'(x)) in one Horner pass.
pub fn eval_with_derivative<F: Field>(f: &F, a: &Poly<F>, x: F::E) -> (F::E, F::E) {
    let mut v = f.zero();
    let mut d = f.zero();
    for &c in a.iter().rev() {
        d = f.add(f.mul(d, x), v);
        v = f.add(f.mul(v, x), c);
    }
    (v, d)
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
    // large quotients: Newton/Barrett division (two multiplications); small: schoolbook
    // reducer range: deg b = n, a.len() <= 2n - 1
    if b.len() >= 64 && a.len() >= b.len() + 32 && a.len() + 3 <= 2 * b.len() {
        return PolyModulus::new(f, b).rem(f, a);
    }
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
    if m.len() > 64 {
        return PolyModulus::new(f, m).powmod_big(f, a, &Big::from_u128(e));
    }
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
pub fn powmod_big<F: Field>(f: &F, a: &Poly<F>, e: &Big, m: &Poly<F>) -> Poly<F> {
    if m.len() > 64 {
        return PolyModulus::new(f, m).powmod_big(f, a, e);
    }
    let mut r = constant(f, f.one());
    let b = rem(f, a, m);
    for i in (0..e.bits()).rev() {
        r = mulmod(f, &r, &r, m);
        if e.bit(i) {
            r = mulmod(f, &r, &b, m);
        }
    }
    r
}

/// Subproduct tree of the linear factors (x - r_i): level 0 = leaves, last level = product.
pub fn subproduct_tree<F: Field>(f: &F, roots: &[F::E]) -> Vec<Vec<Poly<F>>> {
    let mut levels: Vec<Vec<Poly<F>>> =
        vec![roots.iter().map(|&r| vec![f.neg(r), f.one()]).collect()];
    while levels.last().unwrap().len() > 1 {
        let prev = levels.last().unwrap();
        let next: Vec<Poly<F>> = prev
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    mul(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
        levels.push(next);
    }
    levels
}

/// Values of `a` at all points (remainder tree; Horner for few points).
pub fn multipoint_eval<F: Field>(f: &F, a: &Poly<F>, points: &[F::E]) -> Vec<F::E> {
    if points.len() <= 8 || a.len() <= 8 {
        return points.iter().map(|&x| eval(f, a, x)).collect();
    }
    let tree = subproduct_tree(f, points);
    multipoint_eval_tree(f, a, &tree)
}

pub fn multipoint_eval_tree<F: Field>(f: &F, a: &Poly<F>, tree: &[Vec<Poly<F>>]) -> Vec<F::E> {
    let top = tree.len() - 1;
    let mut rems: Vec<Poly<F>> = vec![rem(f, a, &tree[top][0])];
    for lvl in (0..top).rev() {
        let mut next = Vec::with_capacity(tree[lvl].len());
        for (i, r) in rems.iter().enumerate() {
            for c in 0..2 {
                let idx = 2 * i + c;
                if idx < tree[lvl].len() {
                    next.push(rem(f, r, &tree[lvl][idx]));
                }
            }
        }
        rems = next;
    }
    rems.iter()
        .map(|r| if r.is_empty() { f.zero() } else { r[0] })
        .collect()
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
    let half = f.q().sub_small(1).shr(1);
    loop {
        let r = f.random(rng);
        let base = vec![r, f.one()];
        let w = sub(f, &powmod_big(f, &base, &half, g), &constant(f, f.one()));
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
    let xq = powmod_big(f, &x_poly(f), &f.q(), &p);
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
        h = powmod_big(f, &h, &f.q(), &fp);
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
    let half = f.q().sub_small(1).shr(1);
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
            t = powmod_big(f, &t, &f.q(), g);
            s = mulmod(f, &s, &t, g);
        }
        let b = powmod_big(f, &s, &half, g);
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
