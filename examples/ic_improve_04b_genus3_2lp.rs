//! Genus-3 index calculus with double large primes (IC-2LP) on the hyperelliptic GHS cubic cover,
//! with the factor-base and linear-algebra arms of
//! `research/ic_improve_20261014/04b-genus3-2lp-rebalance/PROTOCOL.md`.
//!
//! Dependency-free: `rustc --edition 2021 -O examples/ic_improve_04b_genus3_2lp.rs`
//! or `cargo run --release --example ic_improve_04b_genus3_2lp -- <cmd> ...`.
//!
//! Commands
//!   curves --k K --count C --seed S           build C accepted GHS genus-3 curves over q = prevprime(2^K), JSON lines
//!   solve  --curves FILE --idx I --target T --kappa X --la bfs,sge [--excess E]
//!                                     one relation phase, then each listed linear algebra, JSON line
//!   rho    --curves FILE --idx I --target T   Pollard rho reference on the same instance, JSON line
//!
//! Counted unit: multiplications in F_q (squarings and the multiplications inside Fermat inversions
//! included) plus multiplications mod ell (Montgomery products, Fermat inversions included, and one
//! per scalar-by-entry product in large-prime elimination and merging).
use std::cell::Cell;
use std::collections::HashMap;
use std::time::Instant;

// The base seed is a required argument of `curves` and is stored in every curve record.

// ------------------------------------------------------------------ rng and hashing
fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}
fn mix(a: u64, b: u64) -> u64 {
    splitmix(a ^ splitmix(b.wrapping_add(0x632BE59BD9B4E019)))
}
struct Rng(u64);
impl Rng {
    fn new(s: u64) -> Rng {
        Rng(splitmix(s))
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        splitmix(self.0)
    }
    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

// ------------------------------------------------------------------ integers
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn isqrt(n: u128) -> u128 {
    if n == 0 {
        return 0;
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

// ------------------------------------------------------------------ F_q, q < 2^21 prime, counted
struct Fq {
    p: u64,
    bm: u64,
    cnt: Cell<u64>,
}
impl Fq {
    fn new(p: u64) -> Fq {
        Fq {
            p,
            bm: ((1u128 << 64) / p as u128) as u64,
            cnt: Cell::new(0),
        }
    }
    #[inline]
    fn red(&self, x: u64) -> u64 {
        let qt = ((x as u128 * self.bm as u128) >> 64) as u64;
        let r = x - qt * self.p;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.cnt.set(self.cnt.get() + 1);
        self.red(a * b)
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    #[inline]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn pow(&self, a: u64, mut e: u64) -> u64 {
        let mut r = 1;
        let mut b = a;
        let mut first = true;
        while e > 0 {
            if e & 1 == 1 {
                if first {
                    r = b;
                    first = false;
                } else {
                    r = self.mul(r, b);
                }
            }
            e >>= 1;
            if e > 0 {
                b = self.mul(b, b);
            }
        }
        r
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "F_q inverse of zero");
        self.pow(a, self.p - 2)
    }
    fn legendre(&self, a: u64) -> i32 {
        if a == 0 {
            0
        } else if self.pow(a, (self.p - 1) / 2) == 1 {
            1
        } else {
            -1
        }
    }
    fn int(&self, c: i64) -> u64 {
        c.rem_euclid(self.p as i64) as u64
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        // Tonelli-Shanks, counted
        if a == 0 {
            return Some(0);
        }
        if self.legendre(a) != 1 {
            return None;
        }
        let p = self.p;
        let (mut qq, mut s) = (p - 1, 0u32);
        while qq % 2 == 0 {
            qq /= 2;
            s += 1;
        }
        let mut z = 2;
        while self.legendre(z) != -1 {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(z, qq);
        let mut t = self.pow(a, qq);
        let mut r = self.pow(a, qq.div_ceil(2));
        while t != 1 {
            let mut i = 0;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.mul(b, b);
            }
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
}

// ------------------------------------------------------------------ polynomials over F_q, fixed capacity
const PN: usize = 20;
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Po {
    c: [u64; PN],
    d: i32,
}
impl Po {
    fn zero() -> Po {
        Po { c: [0; PN], d: -1 }
    }
    fn cst(a: u64) -> Po {
        let mut p = Po::zero();
        if a != 0 {
            p.c[0] = a;
            p.d = 0;
        }
        p
    }
    fn x() -> Po {
        let mut p = Po::zero();
        p.c[1] = 1;
        p.d = 1;
        p
    }
    fn from(v: &[u64]) -> Po {
        let mut p = Po::zero();
        for (i, &a) in v.iter().enumerate() {
            p.c[i] = a;
        }
        p.d = v.len() as i32 - 1;
        p.trim();
        p
    }
    fn trim(&mut self) {
        while self.d >= 0 && self.c[self.d as usize] == 0 {
            self.d -= 1;
        }
    }
    fn lead(&self) -> u64 {
        if self.d < 0 {
            0
        } else {
            self.c[self.d as usize]
        }
    }
    fn is_zero(&self) -> bool {
        self.d < 0
    }
}
fn padd(f: &Fq, a: &Po, b: &Po) -> Po {
    let d = a.d.max(b.d);
    let mut r = Po::zero();
    if d >= 0 {
        for i in 0..=d as usize {
            r.c[i] = f.add(a.c[i], b.c[i]);
        }
    }
    r.d = d;
    r.trim();
    r
}
fn psub(f: &Fq, a: &Po, b: &Po) -> Po {
    let d = a.d.max(b.d);
    let mut r = Po::zero();
    if d >= 0 {
        for i in 0..=d as usize {
            r.c[i] = f.sub(a.c[i], b.c[i]);
        }
    }
    r.d = d;
    r.trim();
    r
}
fn pneg(f: &Fq, a: &Po) -> Po {
    let mut r = *a;
    if a.d >= 0 {
        for i in 0..=a.d as usize {
            r.c[i] = f.neg(a.c[i]);
        }
    }
    r
}
fn pmul(f: &Fq, a: &Po, b: &Po) -> Po {
    if a.d < 0 || b.d < 0 {
        return Po::zero();
    }
    let d = (a.d + b.d) as usize;
    assert!(d < PN, "poly capacity");
    let mut r = Po::zero();
    for i in 0..=a.d as usize {
        if a.c[i] == 0 {
            continue;
        }
        for j in 0..=b.d as usize {
            r.c[i + j] = f.add(r.c[i + j], f.mul(a.c[i], b.c[j]));
        }
    }
    r.d = d as i32;
    r.trim();
    r
}
fn psqr(f: &Fq, a: &Po) -> Po {
    if a.d < 0 {
        return Po::zero();
    }
    let n = a.d as usize;
    assert!(2 * n < PN);
    let mut r = Po::zero();
    for i in 0..=n {
        r.c[2 * i] = f.add(r.c[2 * i], f.mul(a.c[i], a.c[i]));
        for j in (i + 1)..=n {
            let t = f.mul(a.c[i], a.c[j]);
            r.c[i + j] = f.add(r.c[i + j], f.add(t, t));
        }
    }
    r.d = 2 * n as i32;
    r.trim();
    r
}
fn pscale(f: &Fq, a: &Po, s: u64) -> Po {
    if s == 1 {
        return *a;
    }
    let mut r = *a;
    if a.d >= 0 {
        for i in 0..=a.d as usize {
            r.c[i] = f.mul(a.c[i], s);
        }
    }
    r.trim();
    r
}
fn pdivrem(f: &Fq, a: &Po, b: &Po) -> (Po, Po) {
    assert!(!b.is_zero(), "poly division by zero");
    let mut r = *a;
    let mut qt = Po::zero();
    if a.d < b.d {
        return (qt, r);
    }
    let lb = b.lead();
    let inv = if lb == 1 { 1 } else { f.inv(lb) };
    let bd = b.d as usize;
    let mut i = a.d as usize;
    loop {
        let top = r.c[i];
        if top != 0 {
            let co = if inv == 1 { top } else { f.mul(top, inv) };
            qt.c[i - bd] = co;
            for j in 0..bd {
                r.c[i - bd + j] = f.sub(r.c[i - bd + j], f.mul(co, b.c[j]));
            }
            r.c[i] = 0;
        }
        if i == bd {
            break;
        }
        i -= 1;
    }
    qt.d = a.d - b.d;
    qt.trim();
    r.d = b.d - 1;
    r.trim();
    (qt, r)
}
fn prem(f: &Fq, a: &Po, b: &Po) -> Po {
    if a.d < b.d {
        *a
    } else {
        pdivrem(f, a, b).1
    }
}
fn pexact(f: &Fq, a: &Po, b: &Po) -> Po {
    let (qt, r) = pdivrem(f, a, b);
    assert!(r.is_zero(), "inexact poly division");
    qt
}
fn pmonic(f: &Fq, a: &Po) -> Po {
    let l = a.lead();
    if l == 0 || l == 1 {
        *a
    } else {
        pscale(f, a, f.inv(l))
    }
}
/// returns (g monic, s, t) with s a + t b = g
fn pxgcd(f: &Fq, a: &Po, b: &Po) -> (Po, Po, Po) {
    let (mut r0, mut r1) = (*a, *b);
    let (mut s0, mut s1) = (Po::cst(1), Po::zero());
    let (mut t0, mut t1) = (Po::zero(), Po::cst(1));
    while !r1.is_zero() {
        let (qt, r) = pdivrem(f, &r0, &r1);
        let s2 = psub(f, &s0, &pmul(f, &qt, &s1));
        let t2 = psub(f, &t0, &pmul(f, &qt, &t1));
        r0 = r1;
        r1 = r;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    let l = r0.lead();
    if l == 0 {
        return (r0, s0, t0);
    }
    let li = if l == 1 { 1 } else { f.inv(l) };
    (pscale(f, &r0, li), pscale(f, &s0, li), pscale(f, &t0, li))
}
/// returns (g monic, s) with s a = g mod b
fn pxgcd_s(f: &Fq, a: &Po, b: &Po) -> (Po, Po) {
    let (mut r0, mut r1) = (*a, *b);
    let (mut s0, mut s1) = (Po::cst(1), Po::zero());
    while !r1.is_zero() {
        let (qt, r) = pdivrem(f, &r0, &r1);
        let s2 = psub(f, &s0, &pmul(f, &qt, &s1));
        r0 = r1;
        r1 = r;
        s0 = s1;
        s1 = s2;
    }
    let l = r0.lead();
    let li = if l == 1 { 1 } else { f.inv(l) };
    (pscale(f, &r0, li), pscale(f, &s0, li))
}
fn pgcd(f: &Fq, a: &Po, b: &Po) -> Po {
    let (mut x, mut y) = (*a, *b);
    while !y.is_zero() {
        let r = prem(f, &x, &y);
        x = y;
        y = r;
    }
    pmonic(f, &x)
}
fn peval(f: &Fq, a: &Po, x: u64) -> u64 {
    let mut r = 0;
    if a.d < 0 {
        return 0;
    }
    for i in (0..=a.d as usize).rev() {
        r = f.add(f.mul(r, x), a.c[i]);
    }
    r
}
fn pderiv(f: &Fq, a: &Po) -> Po {
    let mut r = Po::zero();
    if a.d >= 1 {
        for i in 1..=a.d as usize {
            r.c[i - 1] = f.mul(a.c[i], f.int(i as i64));
        }
        r.d = a.d - 1;
        r.trim();
    }
    r
}
/// base^e mod m
fn ppowmod(f: &Fq, base: &Po, mut e: u64, m: &Po) -> Po {
    let mut r = Po::cst(1);
    let mut b = prem(f, base, m);
    let mut first = true;
    while e > 0 {
        if e & 1 == 1 {
            r = if first {
                b
            } else {
                prem(f, &pmul(f, &r, &b), m)
            };
            first = false;
        }
        e >>= 1;
        if e > 0 {
            b = prem(f, &psqr(f, &b), m);
        }
    }
    r
}
/// x^e mod m (left-to-right; multiplication by x is a shift, then reduction)
fn pxpow(f: &Fq, e: u64, m: &Po) -> Po {
    let bits = 64 - e.leading_zeros();
    let mut r = prem(f, &Po::x(), m);
    for i in (0..bits - 1).rev() {
        r = prem(f, &psqr(f, &r), m);
        if (e >> i) & 1 == 1 {
            let mut s = Po::zero();
            if r.d >= 0 {
                for j in 0..=r.d as usize {
                    s.c[j + 1] = r.c[j];
                }
                s.d = r.d + 1;
            }
            r = prem(f, &s, m);
        }
    }
    r
}
/// roots of a squarefree polynomial that splits into distinct linear factors over F_q
fn split_roots(f: &Fq, u: &Po, rng: &mut Rng, out: &mut Vec<u64>) {
    let u = pmonic(f, u);
    match u.d {
        d if d <= 0 => {}
        1 => out.push(f.neg(u.c[0])),
        2 => {
            let (b, c) = (u.c[1], u.c[0]);
            let disc = f.sub(f.mul(b, b), f.mul(4 % f.p, c));
            let s = f.sqrt(disc).expect("split quadratic");
            let i2 = f.p.div_ceil(2);
            out.push(f.mul(f.sub(s, b), i2));
            out.push(f.mul(f.sub(f.neg(s), b), i2));
        }
        _ => loop {
            let dl = rng.below(f.p);
            let base = Po::from(&[dl, 1]);
            let g = pgcd(
                f,
                &u,
                &psub(f, &ppowmod(f, &base, (f.p - 1) / 2, &u), &Po::cst(1)),
            );
            if g.d >= 1 && g.d < u.d {
                let h = pexact(f, &u, &g);
                split_roots(f, &g, rng, out);
                split_roots(f, &h, rng, out);
                break;
            }
        },
    }
}

// ------------------------------------------------------------------ Jacobian of y^2 = f(x), deg f = 7 monic (Cantor)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Dv {
    u: Po,
    v: Po,
}
fn jzero() -> Dv {
    Dv {
        u: Po::cst(1),
        v: Po::zero(),
    }
}
fn jadd(f: &Fq, fc: &Po, a: &Dv, b: &Dv) -> Dv {
    if a.u.d == 0 {
        return *b;
    }
    if b.u.d == 0 {
        return *a;
    }
    let (u1, v1, u2, v2) = (&a.u, &a.v, &b.u, &b.v);
    let (d1, e1) = pxgcd_s(f, u1, u2);
    let (mut u, mut v);
    if d1.d == 0 {
        let k = prem(f, &pmul(f, &psub(f, v2, v1), &e1), u2);
        u = pmul(f, u1, u2);
        v = padd(f, v1, &pmul(f, &k, u1));
    } else {
        let e2 = pexact(f, &psub(f, &d1, &pmul(f, &e1, u1)), u2);
        let (d, c1, c2) = pxgcd(f, &d1, &padd(f, v1, v2));
        let (s1, s2, s3) = (pmul(f, &c1, &e1), pmul(f, &c1, &e2), c2);
        u = pexact(f, &pmul(f, u1, u2), &psqr(f, &d));
        let num = padd(
            f,
            &padd(
                f,
                &pmul(f, &pmul(f, &s1, u1), v2),
                &pmul(f, &pmul(f, &s2, u2), v1),
            ),
            &pmul(f, &s3, &padd(f, &pmul(f, v1, v2), fc)),
        );
        v = prem(f, &pexact(f, &num, &d), &u);
    }
    v = prem(f, &v, &u);
    while u.d > 3 {
        u = pmonic(f, &pexact(f, &psub(f, fc, &psqr(f, &v)), &u));
        v = prem(f, &pneg(f, &v), &u);
    }
    Dv { u, v }
}
fn jmul(f: &Fq, fc: &Po, mut n: u64, a: &Dv) -> Dv {
    let mut r = jzero();
    let mut b = *a;
    while n > 0 {
        if n & 1 == 1 {
            r = jadd(f, fc, &r, &b);
        }
        n >>= 1;
        if n > 0 {
            b = jadd(f, fc, &b, &b);
        }
    }
    r
}
fn jvalid(f: &Fq, fc: &Po, a: &Dv) -> bool {
    a.u.lead() == 1
        && a.u.d <= 3
        && a.v.d < a.u.d.max(0) + if a.u.d == 0 { 1 } else { 0 }
        && prem(f, &psub(f, &psqr(f, &a.v), fc), &a.u).is_zero()
}
fn jhash(a: &Dv) -> u64 {
    let mut h = 0x1234567u64;
    for i in 0..4 {
        h = mix(h, a.u.c[i]);
        h = mix(h, a.v.c[i]);
    }
    h
}
/// random divisor: sum of three random affine points
fn jrandom(f: &Fq, fc: &Po, rng: &mut Rng) -> Dv {
    let mut r = jzero();
    for _ in 0..3 {
        loop {
            let x = rng.below(f.p);
            let y2 = peval(f, fc, x);
            if let Some(y) = f.sqrt(y2) {
                let y = if rng.next() & 1 == 1 { f.neg(y) } else { y };
                r = jadd(
                    f,
                    fc,
                    &r,
                    &Dv {
                        u: Po::from(&[f.neg(x), 1]),
                        v: Po::cst(y),
                    },
                );
                break;
            }
        }
    }
    r
}

// ------------------------------------------------------------------ F_{q^3} = F_q[t]/(t^3 - s1 t - s0), construction only
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct K3([u64; 3]);
struct Fk<'a> {
    f: &'a Fq,
    s0: u64,
    s1: u64,
}
impl<'a> Fk<'a> {
    fn add(&self, a: K3, b: K3) -> K3 {
        K3([
            self.f.add(a.0[0], b.0[0]),
            self.f.add(a.0[1], b.0[1]),
            self.f.add(a.0[2], b.0[2]),
        ])
    }
    fn sub(&self, a: K3, b: K3) -> K3 {
        K3([
            self.f.sub(a.0[0], b.0[0]),
            self.f.sub(a.0[1], b.0[1]),
            self.f.sub(a.0[2], b.0[2]),
        ])
    }
    fn mul(&self, a: K3, b: K3) -> K3 {
        let f = self.f;
        let mut c = [0u64; 5];
        for i in 0..3 {
            for j in 0..3 {
                c[i + j] = f.add(c[i + j], f.mul(a.0[i], b.0[j]));
            }
        }
        // t^4 = s1 t^2 + s0 t ; t^3 = s1 t + s0
        let (c4, c3) = (c[4], c[3]);
        c[2] = f.add(c[2], f.mul(c4, self.s1));
        c[1] = f.add(c[1], f.mul(c4, self.s0));
        c[1] = f.add(c[1], f.mul(c3, self.s1));
        c[0] = f.add(c[0], f.mul(c3, self.s0));
        K3([c[0], c[1], c[2]])
    }
    fn pow(&self, a: K3, mut e: u64) -> K3 {
        let mut r = K3([1, 0, 0]);
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            e >>= 1;
            if e > 0 {
                b = self.mul(b, b);
            }
        }
        r
    }
    fn frob(&self, a: K3) -> K3 {
        self.pow(a, self.f.p)
    }
    fn inv(&self, a: K3) -> K3 {
        let q = self.f.p;
        self.pow(a, q * q * q - 2)
    }
    fn inq(&self, a: K3) -> bool {
        a.0[1] == 0 && a.0[2] == 0
    }
    fn tr(&self, a: K3) -> u64 {
        let s = self.add(self.add(a, self.frob(a)), self.frob(self.frob(a)));
        assert!(self.inq(s), "trace not in F_q");
        s.0[0]
    }
}
type KP = Vec<K3>; // low -> high, trimmed
fn ktrim(mut a: KP) -> KP {
    while let Some(l) = a.last() {
        if *l == K3([0, 0, 0]) {
            a.pop();
        } else {
            break;
        }
    }
    a
}
fn kmul(k: &Fk, a: &KP, b: &KP) -> KP {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![K3([0, 0, 0]); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            r[i + j] = k.add(r[i + j], k.mul(a[i], b[j]));
        }
    }
    ktrim(r)
}
fn ksub(k: &Fk, a: &KP, b: &KP) -> KP {
    let n = a.len().max(b.len());
    ktrim(
        (0..n)
            .map(|i| {
                k.sub(
                    *a.get(i).unwrap_or(&K3([0, 0, 0])),
                    *b.get(i).unwrap_or(&K3([0, 0, 0])),
                )
            })
            .collect(),
    )
}
fn kscale(k: &Fk, a: &KP, s: K3) -> KP {
    ktrim(a.iter().map(|&c| k.mul(c, s)).collect())
}
fn kdivrem(k: &Fk, a: &KP, b: &KP) -> (KP, KP) {
    let mut r = a.clone();
    if r.len() < b.len() {
        return (vec![], r);
    }
    let li = k.inv(*b.last().unwrap());
    let mut qt = vec![K3([0, 0, 0]); r.len() - b.len() + 1];
    for i in (0..qt.len()).rev() {
        let co = k.mul(r[i + b.len() - 1], li);
        qt[i] = co;
        for j in 0..b.len() {
            r[i + j] = k.sub(r[i + j], k.mul(co, b[j]));
        }
    }
    (ktrim(qt), ktrim(r))
}
fn kgcd(k: &Fk, a: &KP, b: &KP) -> KP {
    let (mut x, mut y) = (a.clone(), b.clone());
    while !y.is_empty() {
        let r = kdivrem(k, &x, &y).1;
        x = y;
        y = r;
    }
    let li = k.inv(*x.last().unwrap());
    kscale(k, &x, li)
}

// ------------------------------------------------------------------ the GHS m = 3, r = 5 hyperelliptic cover (T-2463 construction)
struct Curve {
    k: u32,
    q: u64,
    idx: u32,
    f: Po,
    n: u64,
    ell: u64,
    h: u64,
    seed: u64,
}

/// One candidate. Returns the monic degree-7 model and construction facts, or why it was rejected.
fn ghs_candidate(f: &Fq, rng: &mut Rng) -> Result<(Po, String), String> {
    let q = f.p;
    // K = F_q[t]/(t^3 - s1 t - s0) irreducible: no root in F_q
    let (s0, s1) = loop {
        let s0 = 1 + rng.below(q - 1);
        let s1 = rng.below(q);
        let m = Po::from(&[f.neg(s0), f.neg(s1), 0, 1]);
        let xq = pxpow(f, q, &m);
        if pgcd(f, &m, &psub(f, &xq, &Po::x())).d == 0 {
            break (s0, s1);
        }
    };
    let k = Fk { f, s0, s1 };
    let (a, b) = loop {
        let a = rng.below(q);
        let b = rng.below(q);
        if a != b {
            break (a, b);
        }
    };
    let p = loop {
        let p = K3([rng.below(q), rng.below(q), rng.below(q)]);
        if !k.inq(p) {
            break p;
        }
    };
    let r0 = p;
    let r1 = k.frob(r0);
    let r2 = k.frob(r1);
    let c0 = k.sub(r2, r0);
    // Gram matrix of the F_q quadratic form Tr(c0 xi^2), xi = Z0 + Z1 t + Z2 t^2
    let tp = [
        K3([1, 0, 0]),
        K3([0, 1, 0]),
        K3([0, 0, 1]),
        k.mul(K3([0, 0, 1]), K3([0, 1, 0])),
        k.mul(K3([0, 0, 1]), K3([0, 0, 1])),
    ];
    let mut g = [[0u64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            g[i][j] = k.tr(k.mul(c0, tp[i + j]));
        }
    }
    let bil = |u: &[u64; 3], v: &[u64; 3]| -> u64 {
        let mut s = 0;
        for i in 0..3 {
            for j in 0..3 {
                s = f.add(s, f.mul(f.mul(u[i], g[i][j]), v[j]));
            }
        }
        s
    };
    // a rational point of the conic
    let ps = loop {
        let u = [rng.below(q), rng.below(q), rng.below(q)];
        let v = [rng.below(q), rng.below(q), rng.below(q)];
        let (qu, quv, qv) = (bil(&u, &u), bil(&u, &v), bil(&v, &v));
        if qv == 0 {
            if v != [0, 0, 0] {
                break v;
            }
            continue;
        }
        let disc = f.sub(f.mul(quv, quv), f.mul(qu, qv));
        if let Some(sq) = f.sqrt(disc) {
            let s = f.mul(f.sub(sq, quv), f.inv(qv));
            let w = [
                f.add(u[0], f.mul(s, v[0])),
                f.add(u[1], f.mul(s, v[1])),
                f.add(u[2], f.mul(s, v[2])),
            ];
            if w != [0, 0, 0] {
                break w;
            }
        }
    };
    let l = [
        bil(&ps, &[1, 0, 0]),
        bil(&ps, &[0, 1, 0]),
        bil(&ps, &[0, 0, 1]),
    ];
    let kk = (0..3).find(|&i| l[i] != 0).ok_or("degenerate conic")?;
    let mut v1 = [0u64; 3];
    v1[kk] = 1;
    let v0 = loop {
        let r = [rng.below(q), rng.below(q), rng.below(q)];
        let lr = f.add(
            f.add(f.mul(l[0], r[0]), f.mul(l[1], r[1])),
            f.mul(l[2], r[2]),
        );
        let mut w = r;
        w[kk] = f.sub(w[kk], f.mul(lr, f.inv(l[kk])));
        let par = (0..3).all(|i| {
            let j = (i + 1) % 3;
            f.mul(w[i], ps[j]) == f.mul(w[j], ps[i])
        });
        if w != [0, 0, 0] && !par {
            break w;
        }
    };
    let (qv0, bv01, qv1, b1) = (bil(&v0, &v0), bil(&v0, &v1), bil(&v1, &v1), bil(&ps, &v1));
    let two = 2 % q;
    let mut zp = [[0u64; 3]; 3]; // zp[j][deg]
    for j in 0..3 {
        zp[j] = [
            f.sub(f.mul(qv1, ps[j]), f.mul(f.mul(two, b1), v1[j])),
            f.sub(f.mul(f.mul(two, bv01), ps[j]), f.mul(f.mul(two, b1), v0[j])),
            f.mul(qv0, ps[j]),
        ];
    }
    for _ in 0..4 {
        let z = rng.below(q);
        let pt = [
            peval(f, &Po::from(&zp[0]), z),
            peval(f, &Po::from(&zp[1]), z),
            peval(f, &Po::from(&zp[2]), z),
        ];
        assert!(bil(&pt, &pt) == 0, "conic parametrization");
    }
    // xi(z) in K[z]: coefficient of z^d is sum_j zp[j][d] t^j
    let xi: KP = ktrim((0..3).map(|d| K3([zp[0][d], zp[1][d], zp[2][d]])).collect());
    let s2 = |a: K3| k.frob(k.frob(a));
    let xi2: KP = xi.iter().map(|&c| s2(c)).collect();
    let (a0, a2) = (kmul(&k, &xi, &xi), kmul(&k, &xi2, &xi2));
    let num = ksub(&k, &kscale(&k, &a0, r0), &kscale(&k, &a2, r1));
    let den = ksub(&k, &a0, &a2);
    if den.is_empty() {
        return Err("x denominator zero".into());
    }
    let gg = kgcd(&k, &num, &den);
    let (mut xn, mut xd) = (kdivrem(&k, &num, &gg).0, kdivrem(&k, &den, &gg).0);
    let li = k.inv(*xd.last().unwrap());
    xn = kscale(&k, &xn, li);
    xd = kscale(&k, &xd, li);
    if !xn.iter().chain(xd.iter()).all(|c| k.inq(*c)) {
        return Err("P1: x(z) not over F_q".into());
    }
    let xnp = Po::from(&xn.iter().map(|c| c.0[0]).collect::<Vec<_>>());
    let xdp = Po::from(&xd.iter().map(|c| c.0[0]).collect::<Vec<_>>());
    let ha = psub(f, &xnp, &pscale(f, &xdp, a));
    let hb = psub(f, &xnp, &pscale(f, &xdp, b));
    let h = pmul(f, &ha, &hb);
    if h.d < 7 || h.d > 8 {
        return Err(format!("deg h = {}", h.d));
    }
    if pgcd(f, &h, &pderiv(f, &h)).d != 0 {
        return Err("h not squarefree".into());
    }
    let lin = |p: &Po| -> Po {
        if p.d < 1 {
            return Po::cst(1);
        }
        let m = pmonic(f, p);
        pgcd(f, &m, &psub(f, &pxpow(f, q, &m), &Po::x()))
    };
    let (la, lb) = (lin(&ha), lin(&hb));
    let mut roots = vec![];
    let mut rr = Rng::new(rng.next());
    split_roots(f, &la, &mut rr, &mut roots);
    split_roots(f, &lb, &mut rr, &mut roots);
    let (na, nb) = (la.d.max(0), lb.d.max(0));
    let f1 = if h.d == 7 {
        h
    } else {
        if roots.is_empty() {
            return Err("no rational Weierstrass point".into());
        }
        let e = roots[0]; // h(z + e), then reverse: z^8 h(e + 1/z)
        let mut sh = Po::zero();
        let mut pw = Po::cst(1);
        let lin_e = Po::from(&[e, 1]);
        for i in 0..=h.d as usize {
            sh = padd(f, &sh, &pscale(f, &pw, h.c[i]));
            pw = pmul(f, &pw, &lin_e);
        }
        assert!(sh.c[0] == 0, "shift root");
        let mut rv = Po::zero();
        for i in 0..=8usize {
            rv.c[8 - i] = sh.c[i];
        }
        rv.d = 8;
        rv.trim();
        assert!(rv.d == 7);
        rv
    };
    let c = f1.lead();
    let ci = f.inv(c);
    let mut fm = Po::zero();
    let mut cp = f.mul(f.mul(ci, ci), f.mul(ci, ci));
    cp = f.mul(cp, cp); // c^{-8}
    for i in 0..=7usize {
        fm.c[i] = f.mul(f1.c[i], cp);
        cp = f.mul(cp, c);
    }
    fm.d = 7;
    fm.trim();
    assert!(fm.lead() == 1 && fm.d == 7);
    if pgcd(f, &fm, &pderiv(f, &fm)).d != 0 {
        return Err("model not squarefree".into());
    }
    let info = format!("\"K_mod\":[{},{}],\"a\":{},\"b\":{},\"p\":[{},{},{}],\"deg_h\":{},\"rat_weier_over_a\":{},\"rat_weier_over_b\":{}", s0, s1, a, b, p.0[0], p.0[1], p.0[2], h.d, na, nb);
    Ok((fm, info))
}

/// #Jac(C)(F_q) = q^3 + 1 - t with |t| <= 2 q^{3/2}; baby-step giant-step on random divisors, then 4 confirmations.
fn jac_order(f: &Fq, fc: &Po, rng: &mut Rng) -> Option<u64> {
    let q = f.p as u128;
    let center = q * q * q + 1;
    let bnd = isqrt(4 * q * q * q) + 1;
    let lo = (center - bnd) as u64;
    let w = (2 * bnd) as u64;
    let m = isqrt(w as u128) as u64 + 1;
    let mut cands: Option<Vec<u64>> = None;
    for _ in 0..3 {
        let d = jrandom(f, fc, rng);
        let mut tab: HashMap<u64, Vec<u32>> = HashMap::with_capacity(m as usize * 2);
        let mut cur = jzero();
        for j in 0..m {
            tab.entry(jhash(&cur)).or_default().push(j as u32);
            cur = jadd(f, fc, &cur, &d);
        }
        let step = cur;
        let mut g = jmul(f, fc, lo, &d);
        let mut found = vec![];
        for i in 0..=(w / m + 1) {
            if let Some(js) = tab.get(&jhash(&g)) {
                for &j in js {
                    let mm = lo as i128 + (i * m) as i128 - j as i128;
                    if mm >= lo as i128
                        && mm <= (lo + w) as i128
                        && jmul(f, fc, mm as u64, &d).u.d == 0
                    {
                        found.push(mm as u64);
                    }
                }
            }
            g = jadd(f, fc, &g, &step);
        }
        found.sort();
        found.dedup();
        cands = Some(match cands {
            None => found,
            Some(c) => c.into_iter().filter(|x| found.contains(x)).collect(),
        });
        if cands.as_ref().unwrap().len() <= 1 {
            break;
        }
    }
    let c = cands?;
    if c.len() != 1 {
        return None;
    }
    let n = c[0];
    for _ in 0..4 {
        let d = jrandom(f, fc, rng);
        if jmul(f, fc, n, &d).u.d != 0 {
            return None;
        }
    }
    Some(n)
}

fn prev_prime(n: u64) -> u64 {
    let mut p = n - 1;
    while !is_prime(p) {
        p -= 1;
    }
    p
}

fn cmd_curves(kbits: u32, count: u32, seed: u64) {
    let q = prev_prime(1u64 << kbits);
    let f = Fq::new(q);
    let mut rng = Rng::new(mix(seed, mix(0xC0FFEE, kbits as u64)));
    let (mut tried, mut acc) = (0u64, 0u32);
    let mut rej: HashMap<String, u64> = HashMap::new();
    let t0 = Instant::now();
    while acc < count {
        tried += 1;
        let (fc, info) = match ghs_candidate(&f, &mut rng) {
            Ok(x) => x,
            Err(e) => {
                *rej.entry(e.split(':').next().unwrap().to_string())
                    .or_default() += 1;
                continue;
            }
        };
        let n = match jac_order(&f, &fc, &mut rng) {
            Some(n) => n,
            None => {
                *rej.entry("order not unique".into()).or_default() += 1;
                continue;
            }
        };
        let (h, ell) = match [4u64, 8].iter().find(|&&h| n % h == 0 && is_prime(n / h)) {
            Some(&h) => (h, n / h),
            None => {
                *rej.entry("cofactor not 4 or 8".into()).or_default() += 1;
                continue;
            }
        };
        // GHS shape check: N1 = #C(F_q) = q + 1 (a1 = 0 in L(T) = 1 - t T^3 + q^3 T^6)
        let mut n1 = 1u64;
        for x in 0..q {
            n1 += (1 + f.legendre(peval(&f, &fc, x))) as u64;
        }
        let t = (q as i128 * q as i128 * q as i128 + 1 - n as i128) as i64;
        assert!(
            n1 == q + 1,
            "N1 != q+1 on an accepted curve: construction error"
        );
        acc += 1;
        let fs: Vec<String> = (0..8).map(|i| fc.c[i].to_string()).collect();
        println!("{{\"tool\":\"ic04b\",\"cmd\":\"curves\",\"seed\":{},\"k\":{},\"q\":{},\"idx\":{},\"f\":[{}],\"N\":{},\"t\":{},\"h\":{},\"ell\":{},\"N1\":{},{},\"candidates_tried\":{},\"wall_s\":{:.3}}}",
            seed, kbits, q, acc, fs.join(","), n, t, h, ell, n1, info, tried, t0.elapsed().as_secs_f64());
    }
    let mut rv: Vec<_> = rej.into_iter().collect();
    rv.sort();
    eprintln!(
        "k={} q={} accepted={} tried={} rejected={:?}",
        kbits, q, acc, tried, rv
    );
}

// ------------------------------------------------------------------ minimal JSON field access for our own records
fn jnum(line: &str, key: &str) -> Option<i128> {
    let k = format!("\"{}\":", key);
    let i = line.find(&k)? + k.len();
    let s: String = line[i..]
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-')
        .collect();
    s.parse().ok()
}
fn jarr(line: &str, key: &str) -> Option<Vec<u64>> {
    let k = format!("\"{}\":[", key);
    let i = line.find(&k)? + k.len();
    let j = line[i..].find(']')? + i;
    Some(
        line[i..j]
            .split(',')
            .map(|s| s.trim().parse().unwrap())
            .collect(),
    )
}
fn load_curve(path: &str, idx: u32) -> Curve {
    let txt = std::fs::read_to_string(path).expect("curves file");
    for line in txt.lines() {
        if jnum(line, "idx") == Some(idx as i128) && line.contains("\"cmd\":\"curves\"") {
            let q = jnum(line, "q").unwrap() as u64;
            let f = Fq::new(q);
            let fc = Po::from(&jarr(line, "f").unwrap());
            let (n, ell, h) = (
                jnum(line, "N").unwrap() as u64,
                jnum(line, "ell").unwrap() as u64,
                jnum(line, "h").unwrap() as u64,
            );
            assert!(
                fc.d == 7 && fc.lead() == 1 && pgcd(&f, &fc, &pderiv(&f, &fc)).d == 0,
                "curve model"
            );
            assert!(n == h * ell && is_prime(ell), "order factorisation");
            let seed = jnum(line, "seed").unwrap() as u64;
            let mut rng = Rng::new(mix(seed, 77));
            for _ in 0..2 {
                assert!(
                    jmul(&f, &fc, n, &jrandom(&f, &fc, &mut rng)).u.d == 0,
                    "group order check"
                );
            }
            return Curve {
                k: jnum(line, "k").unwrap() as u32,
                q,
                idx,
                f: fc,
                n,
                ell,
                h,
                seed,
            };
        }
    }
    panic!("curve {} not in {}", idx, path);
}

// ------------------------------------------------------------------ F_ell, Montgomery, counted
struct Ml {
    n: u64,
    ninv: u64,
    r2: u64,
    cnt: Cell<u64>,
    small: Vec<u64>,
}
impl Ml {
    fn new(n: u64) -> Ml {
        assert!(n % 2 == 1 && n < (1 << 62));
        let mut x = n;
        for _ in 0..6 {
            x = x.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(x)));
        }
        let r1 = ((1u128 << 64) % n as u128) as u64;
        let r2 = mulmod(r1, r1, n);
        let mut m = Ml {
            n,
            ninv: x.wrapping_neg(),
            r2,
            cnt: Cell::new(0),
            small: vec![],
        };
        m.small = (-16i64..=16)
            .map(|c| m.to_m(c.rem_euclid(n as i64) as u64))
            .collect();
        m
    }
    #[inline]
    fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.ninv);
        let u = ((t + m as u128 * self.n as u128) >> 64) as u64;
        if u >= self.n {
            u - self.n
        } else {
            u
        }
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.cnt.set(self.cnt.get() + 1);
        self.redc(a as u128 * b as u128)
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.n {
            s - self.n
        } else {
            s
        }
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.n - b
        }
    }
    #[inline]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.n - a
        }
    }
    fn to_m(&self, x: u64) -> u64 {
        self.mul(x % self.n, self.r2)
    }
    fn decode_m(&self, x: u64) -> u64 {
        self.redc(x as u128)
    }
    fn sm(&self, c: i64) -> u64 {
        if (-16..=16).contains(&c) {
            self.small[(c + 16) as usize]
        } else {
            self.to_m(c.rem_euclid(self.n as i64) as u64)
        }
    }
    fn pow(&self, a: u64, mut e: u64) -> u64 {
        let mut r = self.small[17];
        let mut b = a;
        let mut first = true;
        while e > 0 {
            if e & 1 == 1 {
                if first {
                    r = b;
                    first = false;
                } else {
                    r = self.mul(r, b);
                }
            }
            e >>= 1;
            if e > 0 {
                b = self.mul(b, b);
            }
        }
        r
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "F_ell inverse of zero");
        self.pow(a, self.n - 2)
    }
}

// ------------------------------------------------------------------ relation phase
#[derive(Clone, Copy)]
struct Rel {
    fb: [(u32, i8); 3],
    nfb: u8,
    lp: [(u32, i8); 2],
    nlp: u8,
    a: u64,
    b: u64,
}
struct RelSet {
    fulls: Vec<Rel>,
    edges: Vec<Rel>,
    fsize: usize,
    comp0_v: u64,
    comp0_e: u64,
    uf: Vec<u32>,
}

struct Inst<'a> {
    f: &'a Fq,
    fc: Po,
    ell: u64,
    p: Dv,
    qd: Dv,
    planted: u64,
}

fn uf_find(par: &mut [u32], mut x: u32) -> u32 {
    let mut r = x;
    while par[r as usize] != r {
        r = par[r as usize];
    }
    while par[x as usize] != r {
        let n = par[x as usize];
        par[x as usize] = r;
        x = n;
    }
    r
}

struct RelStats {
    trials: u64,
    split: u64,
    full: u64,
    p1: u64,
    p2: u64,
    disc: u64,
    fq_setup: u64,
    fq_walk: u64,
    fq_decomp: u64,
    wall: f64,
    target: u64,
}

fn relation_phase(inst: &Inst, kappa: f64, walk_seed: u64) -> (RelSet, RelStats) {
    let f = inst.f;
    let fc = &inst.fc;
    let q = f.p;
    let t0 = Instant::now();
    let c0 = f.cnt.get();
    let fsize = (kappa * (q as f64).powf(2.0 / 3.0)).ceil() as usize;
    let mut fbcol = vec![u32::MAX; q as usize];
    let mut nf = 0usize;
    let mut x = 0u64;
    while nf < fsize {
        assert!(x < q, "factor base larger than available places");
        if f.legendre(peval(f, fc, x)) == 1 {
            fbcol[x as usize] = nf as u32;
            nf += 1;
        }
        x += 1;
    }
    let ml = |a: u64, b: u64| mulmod(a, b, inst.ell);
    let mut rng = Rng::new(walk_seed);
    const R: usize = 32;
    let (a0, b0) = (rng.below(inst.ell), rng.below(inst.ell));
    let mut cur = jadd(f, fc, &jmul(f, fc, a0, &inst.p), &jmul(f, fc, b0, &inst.qd));
    let (mut ca, mut cb) = (a0, b0);
    let mut mults = vec![];
    for _ in 0..R {
        let (al, be) = (rng.below(inst.ell), rng.below(inst.ell));
        mults.push((
            jadd(f, fc, &jmul(f, fc, al, &inst.p), &jmul(f, fc, be, &inst.qd)),
            al,
            be,
        ));
    }
    let _ = ml;
    let fq_setup = f.cnt.get() - c0;
    let target = (fsize + 10) as u64;
    let mut par: Vec<u32> = (0..=q as u32).collect();
    let mut ecnt = vec![0u32; q as usize + 1];
    let mut ncnt = vec![1u32; q as usize + 1];
    ncnt[0] = 0;
    let mut st = RelStats {
        trials: 0,
        split: 0,
        full: 0,
        p1: 0,
        p2: 0,
        disc: 0,
        fq_setup,
        fq_walk: 0,
        fq_decomp: 0,
        wall: 0.0,
        target,
    };
    let mut fulls = vec![];
    let mut edges = vec![];
    let mut roots = vec![];
    let mut rrng = Rng::new(mix(walk_seed, 99));
    let (mut fq_walk, mut fq_dec) = (0u64, 0u64);
    let mut j: u64 = 0;
    loop {
        let r0 = uf_find(&mut par, 0) as usize;
        let have = st.full + (ecnt[r0] as u64).saturating_sub(ncnt[r0] as u64);
        if have >= target {
            break;
        }
        // decompose the current element
        st.trials += 1;
        let cd = f.cnt.get();
        let u = cur.u;
        let mut ok = true;
        if u.d >= 2 {
            let xq = pxpow(f, q, &u);
            ok = xq == prem(f, &Po::x(), &u);
        }
        if ok && u.d >= 1 {
            st.split += 1;
            roots.clear();
            split_roots(f, &u, &mut rrng, &mut roots);
            let mut r = Rel {
                fb: [(0, 0); 3],
                nfb: 0,
                lp: [(0, 0); 2],
                nlp: 0,
                a: ca,
                b: cb,
            };
            let mut nlp = 0;
            for &al in roots.iter() {
                let yv = peval(f, &cur.v, al);
                if yv == 0 {
                    continue;
                }
                let s: i8 = if yv <= q - yv { 1 } else { -1 };
                let col = fbcol[al as usize];
                if col != u32::MAX {
                    r.fb[r.nfb as usize] = (col, s);
                    r.nfb += 1;
                } else {
                    if nlp < 2 {
                        r.lp[nlp] = (al as u32 + 1, s);
                    }
                    nlp += 1;
                }
            }
            fq_dec += f.cnt.get() - cd;
            if nlp == 0 {
                st.full += 1;
                fulls.push(r);
            } else if nlp <= 2 {
                r.nlp = nlp as u8;
                if nlp == 1 {
                    st.p1 += 1;
                } else {
                    st.p2 += 1;
                }
                let (va, vb) = if nlp == 1 {
                    (0u32, r.lp[0].0)
                } else {
                    (r.lp[0].0, r.lp[1].0)
                };
                let (ra, rb) = (uf_find(&mut par, va), uf_find(&mut par, vb));
                if ra != rb {
                    let (big, small) = if ncnt[ra as usize] + ecnt[ra as usize]
                        >= ncnt[rb as usize] + ecnt[rb as usize]
                    {
                        (ra, rb)
                    } else {
                        (rb, ra)
                    };
                    par[small as usize] = big;
                    ecnt[big as usize] += ecnt[small as usize] + 1;
                    ncnt[big as usize] += ncnt[small as usize];
                } else {
                    ecnt[ra as usize] += 1;
                }
                edges.push(r);
            } else {
                st.disc += 1;
            }
        } else {
            fq_dec += f.cnt.get() - cd;
        }
        // walk step: index from the step counter, not from the element (no rho cycling)
        let cw = f.cnt.get();
        let i = (mix(walk_seed, j) % R as u64) as usize;
        j += 1;
        cur = jadd(f, fc, &cur, &mults[i].0);
        ca = (ca + mults[i].1) % inst.ell;
        cb = (cb + mults[i].2) % inst.ell;
        fq_walk += f.cnt.get() - cw;
    }
    let r0 = uf_find(&mut par, 0) as usize;
    st.fq_walk = fq_walk;
    st.fq_decomp = fq_dec;
    st.wall = t0.elapsed().as_secs_f64();
    let (cv, ce) = (ncnt[r0] as u64, ecnt[r0] as u64);
    (
        RelSet {
            fulls,
            edges,
            fsize,
            comp0_v: cv,
            comp0_e: ce,
            uf: par,
        },
        st,
    )
}

// ------------------------------------------------------------------ sparse matrix and Lanczos on the normal equations
struct Mat {
    rows: Vec<Vec<(u32, u64)>>,
    a: Vec<u64>,
    b: Vec<u64>,
    ncols: usize,
} // values and (a, b) in Montgomery form

struct LaOut {
    method: String,
    rows: usize,
    cols: usize,
    nnz: usize,
    elim_mults: u64,
    lanczos_mults: u64,
    iters: u64,
    attempts: u32,
    log: Option<u64>,
    wall: f64,
    extra: String,
}

fn compact_cols(m: &mut Mat) {
    let mut map = HashMap::new();
    for r in m.rows.iter_mut() {
        for e in r.iter_mut() {
            let n = map.len() as u32;
            let c = *map.entry(e.0).or_insert(n);
            e.0 = c;
        }
        r.sort_unstable_by_key(|e| e.0);
    }
    m.ncols = map.len();
}

/// returns (log, iterations, attempts); all products counted in ml.cnt
fn lanczos_log(m: &Mat, ml: &Ml, seed: u64) -> (Option<u64>, u64, u32) {
    let n = m.rows.len();
    let nc = m.ncols;
    let mut iters = 0u64;
    let apply = |v: &[u64], out: &mut Vec<u64>, tmp: &mut Vec<u64>| {
        for t in tmp.iter_mut() {
            *t = 0;
        }
        for (j, r) in m.rows.iter().enumerate() {
            let vj = v[j];
            if vj == 0 {
                ml.cnt.set(ml.cnt.get() + r.len() as u64);
                continue;
            }
            for &(c, val) in r {
                tmp[c as usize] = ml.add(tmp[c as usize], ml.mul(val, vj));
            }
        }
        for (j, r) in m.rows.iter().enumerate() {
            let mut s = 0;
            for &(c, val) in r {
                s = ml.add(s, ml.mul(val, tmp[c as usize]));
            }
            out[j] = s;
        }
    };
    let dot = |x: &[u64], y: &[u64]| -> u64 {
        let mut s = 0;
        for i in 0..x.len() {
            s = ml.add(s, ml.mul(x[i], y[i]));
        }
        s
    };
    let mut rng = Rng::new(seed);
    for attempt in 1..=5u32 {
        let y: Vec<u64> = (0..n).map(|_| rng.below(ml.n)).collect();
        let mut tmp = vec![0u64; nc];
        let mut bvec = vec![0u64; n];
        apply(&y, &mut bvec, &mut tmp);
        let mut x = vec![0u64; n];
        let mut w = bvec.clone();
        let mut wprev = vec![0u64; n];
        let mut vprev = vec![0u64; n];
        let mut tprev_inv = 0u64;
        let mut v = vec![0u64; n];
        let mut ok = false;
        for it in 0..(n + 10) {
            if w.iter().all(|&c| c == 0) {
                ok = true;
                break;
            }
            iters += 1;
            apply(&w, &mut v, &mut tmp);
            let t = dot(&w, &v);
            if t == 0 {
                break;
            }
            let tinv = ml.inv(t);
            let cx = ml.mul(dot(&w, &bvec), tinv);
            for i in 0..n {
                x[i] = ml.add(x[i], ml.mul(cx, w[i]));
            }
            let c1 = ml.mul(dot(&v, &v), tinv);
            let c2 = if it > 0 {
                ml.mul(dot(&v, &vprev), tprev_inv)
            } else {
                0
            };
            let mut wn = vec![0u64; n];
            for i in 0..n {
                let mut s = ml.sub(v[i], ml.mul(c1, w[i]));
                if it > 0 {
                    s = ml.sub(s, ml.mul(c2, wprev[i]));
                }
                wn[i] = s;
            }
            if it == 0 { /* no c2 products performed */ }
            wprev = std::mem::replace(&mut w, wn);
            vprev = v.clone();
            tprev_inv = tinv;
        }
        if !ok {
            continue;
        }
        let z: Vec<u64> = (0..n).map(|i| ml.sub(x[i], y[i])).collect();
        if z.iter().all(|&c| c == 0) {
            continue;
        }
        // check M^T z = 0
        for t in tmp.iter_mut() {
            *t = 0;
        }
        for (j, r) in m.rows.iter().enumerate() {
            for &(c, val) in r {
                tmp[c as usize] = ml.add(tmp[c as usize], ml.mul(val, z[j]));
            }
        }
        if tmp.iter().any(|&c| c != 0) {
            continue;
        }
        let sa = dot(&z, &m.a);
        let sb = dot(&z, &m.b);
        if sb == 0 {
            continue;
        }
        let lg = ml.decode_m(ml.mul(ml.neg(sa), ml.inv(sb)));
        return (Some(lg), iters, attempt);
    }
    (None, iters, 5)
}

/// Baseline large-prime elimination: BFS spanning tree of the component of vertex 0, one relation per non-tree edge.
fn la_bfs(rs: &RelSet, ml: &Ml, seed: u64) -> LaOut {
    let t0 = Instant::now();
    let c0 = ml.cnt.get();
    let mut uf = rs.uf.clone();
    let r0 = uf_find(&mut uf, 0);
    let ein: Vec<usize> = (0..rs.edges.len())
        .filter(|&i| {
            let e = &rs.edges[i];
            uf_find(&mut uf, e.lp[0].0) == r0
        })
        .collect();
    // local vertex ids
    let mut loc: HashMap<u32, u32> = HashMap::new();
    loc.insert(0, 0);
    let ends = |e: &Rel| -> (u32, u32) {
        if e.nlp == 1 {
            (0, e.lp[0].0)
        } else {
            (e.lp[0].0, e.lp[1].0)
        }
    };
    for &i in &ein {
        let (x, y) = ends(&rs.edges[i]);
        for v in [x, y] {
            let n = loc.len() as u32;
            loc.entry(v).or_insert(n);
        }
    }
    let nv = loc.len();
    let mut adj: Vec<Vec<(u32, u32)>> = vec![vec![]; nv];
    for (k, &i) in ein.iter().enumerate() {
        let (x, y) = ends(&rs.edges[i]);
        let (lx, ly) = (loc[&x], loc[&y]);
        adj[lx as usize].push((ly, k as u32));
        if lx != ly {
            adj[ly as usize].push((lx, k as u32));
        }
    }
    let mut depth = vec![u32::MAX; nv];
    let mut pedge = vec![u32::MAX; nv];
    let mut intree = vec![false; ein.len()];
    depth[0] = 0;
    let mut qu = std::collections::VecDeque::new();
    qu.push_back(0u32);
    while let Some(x) = qu.pop_front() {
        for &(y, k) in &adj[x as usize] {
            if depth[y as usize] == u32::MAX {
                depth[y as usize] = depth[x as usize] + 1;
                pedge[y as usize] = k;
                intree[k as usize] = true;
                qu.push_back(y);
            }
        }
    }
    let coef_at = |e: &Rel, v: u32| -> i64 {
        let mut s = 0;
        for t in 0..e.nlp as usize {
            if e.lp[t].0 == v {
                s += e.lp[t].1 as i64;
            }
        }
        s
    };
    let mut rows: Vec<Vec<(u32, u64)>> = vec![];
    let mut av = vec![];
    let mut bv = vec![];
    let mut acc: HashMap<u32, i64> = HashMap::new();
    let mut prods = 0u64;
    let mut wsum = 0u64;
    for f in &rs.fulls {
        let mut r: Vec<(u32, u64)> = (0..f.nfb as usize)
            .map(|t| (f.fb[t].0, ml.sm(f.fb[t].1 as i64)))
            .collect();
        r.sort_unstable_by_key(|e| e.0);
        rows.push(r);
        av.push(ml.to_m(f.a));
        bv.push(ml.to_m(f.b));
    }
    let mut maxw = 0usize;
    for (k, &i) in ein.iter().enumerate() {
        if intree[k] {
            continue;
        }
        acc.clear();
        let e = &rs.edges[i];
        let mut ra = ml.to_m(e.a);
        let mut rb = ml.to_m(e.b);
        for t in 0..e.nfb as usize {
            *acc.entry(e.fb[t].0).or_default() += e.fb[t].1 as i64;
        }
        // residual on LP vertices: (global vertex id, coefficient)
        let mut res: Vec<(u32, i64)> = vec![];
        for t in 0..e.nlp as usize {
            res.push((e.lp[t].0, e.lp[t].1 as i64));
        }
        let mut steps = 0u32;
        loop {
            res.retain(|&(v, c)| v != 0 && c != 0);
            // merge duplicate vertices
            res.sort_unstable_by_key(|x| x.0);
            let mut m2: Vec<(u32, i64)> = vec![];
            for &(v, c) in &res {
                if let Some(l) = m2.last_mut() {
                    if l.0 == v {
                        l.1 += c;
                        continue;
                    }
                }
                m2.push((v, c));
            }
            res = m2.into_iter().filter(|x| x.1 != 0).collect();
            if res.is_empty() {
                break;
            }
            let (pi, _) = res
                .iter()
                .enumerate()
                .max_by_key(|(_, &(v, _))| depth[loc[&v] as usize])
                .unwrap();
            let (w, cw) = res[pi];
            let lw = loc[&w] as usize;
            let te = &rs.edges[ein[pedge[lw] as usize]];
            let cte = coef_at(te, w);
            assert!(cte == 1 || cte == -1, "tree edge coefficient");
            let s = -cw * cte; // add s * te
            for t in 0..te.nfb as usize {
                *acc.entry(te.fb[t].0).or_default() += s * te.fb[t].1 as i64;
            }
            for t in 0..te.nlp as usize {
                res.push((te.lp[t].0, s * te.lp[t].1 as i64));
            }
            let sm = ml.sm(s);
            ra = ml.add(ra, ml.mul(sm, ml.to_m(te.a)));
            rb = ml.add(rb, ml.mul(sm, ml.to_m(te.b)));
            prods += te.nfb as u64 + te.nlp as u64; // scalar x entry products (the two (a, b) products are counted by ml)
            steps += 1;
            assert!(steps < 1_000_000, "cycle walk runaway");
        }
        let mut r: Vec<(u32, u64)> = acc
            .iter()
            .filter(|(_, &c)| c != 0)
            .map(|(&col, &c)| (col, ml.sm(c)))
            .collect();
        r.sort_unstable_by_key(|e| e.0);
        wsum += r.len() as u64;
        maxw = maxw.max(r.len());
        rows.push(r);
        av.push(ra);
        bv.push(rb);
    }
    ml.cnt.set(ml.cnt.get() + prods);
    let mut m = Mat {
        rows,
        a: av,
        b: bv,
        ncols: 0,
    };
    compact_cols(&mut m);
    let elim = ml.cnt.get() - c0;
    let c1 = ml.cnt.get();
    let nnz: usize = m.rows.iter().map(|r| r.len()).sum();
    let (lg, iters, att) = lanczos_log(&m, ml, seed);
    let ncyc = m.rows.len() - rs.fulls.len();
    LaOut {
        method: "bfs_lanczos".into(),
        rows: m.rows.len(),
        cols: m.ncols,
        nnz,
        elim_mults: elim,
        lanczos_mults: ml.cnt.get() - c1,
        iters,
        attempts: att,
        log: lg,
        wall: t0.elapsed().as_secs_f64(),
        extra: format!(
            "\"cycles\":{},\"cycle_mean_weight\":{:.3},\"cycle_max_weight\":{}",
            ncyc,
            if ncyc > 0 {
                wsum as f64 / ncyc as f64
            } else {
                0.0
            },
            maxw
        ),
    }
}

/// Structured Gaussian elimination: singleton removal, excess pruning, Markowitz merges while the Lanczos cost estimate falls.
fn la_sge(rs: &RelSet, ml: &Ml, seed: u64, keep_excess: usize) -> LaOut {
    let t0 = Instant::now();
    let c0 = ml.cnt.get();
    let fsz = rs.fsize as u32;
    let mut lpcol: HashMap<u32, u32> = HashMap::new();
    let mut rows: Vec<Vec<(u32, u64)>> = vec![];
    let mut av = vec![];
    let mut bv = vec![];
    for r in rs.fulls.iter().chain(rs.edges.iter()) {
        let mut e: Vec<(u32, u64)> = (0..r.nfb as usize)
            .map(|t| (r.fb[t].0, ml.sm(r.fb[t].1 as i64)))
            .collect();
        for t in 0..r.nlp as usize {
            let n = fsz + lpcol.len() as u32;
            let c = *lpcol.entry(r.lp[t].0).or_insert(n);
            e.push((c, ml.sm(r.lp[t].1 as i64)));
        }
        e.sort_unstable_by_key(|x| x.0);
        rows.push(e);
        av.push(ml.to_m(r.a));
        bv.push(ml.to_m(r.b));
    }
    let ncol = fsz as usize + lpcol.len();
    let nrow0 = rows.len();
    let mut alive = vec![true; rows.len()];
    let mut cw = vec![0u32; ncol];
    let mut crow: Vec<Vec<u32>> = vec![vec![]; ncol];
    for (i, r) in rows.iter().enumerate() {
        for &(c, _) in r {
            cw[c as usize] += 1;
            crow[c as usize].push(i as u32);
        }
    }
    let mut nr = rows.len() as i64;
    let mut nc = cw.iter().filter(|&&w| w > 0).count() as i64;
    let mut nnz: i64 = rows.iter().map(|r| r.len() as i64).sum();
    let col_rows = |crow: &mut Vec<Vec<u32>>,
                    rows: &Vec<Vec<(u32, u64)>>,
                    alive: &Vec<bool>,
                    c: u32|
     -> Vec<u32> {
        let mut v: Vec<u32> = crow[c as usize]
            .iter()
            .cloned()
            .filter(|&i| {
                alive[i as usize] && rows[i as usize].binary_search_by_key(&c, |e| e.0).is_ok()
            })
            .collect();
        v.sort_unstable();
        v.dedup();
        crow[c as usize] = v.clone();
        v
    };
    // delete a row, return columns whose weight changed
    macro_rules! kill {
        ($i:expr) => {{
            let i = $i as usize;
            alive[i] = false;
            nr -= 1;
            nnz -= rows[i].len() as i64;
            let mut ch = vec![];
            for &(c, _) in &rows[i] {
                cw[c as usize] -= 1;
                if cw[c as usize] == 0 {
                    nc -= 1;
                }
                ch.push(c);
            }
            rows[i].clear();
            ch
        }};
    }
    // 1. singletons
    let mut stack: Vec<u32> = (0..ncol as u32).filter(|&c| cw[c as usize] == 1).collect();
    let mut singles = 0u64;
    while let Some(c) = stack.pop() {
        if cw[c as usize] != 1 {
            continue;
        }
        let rs1 = col_rows(&mut crow, &rows, &alive, c);
        if rs1.len() != 1 {
            continue;
        }
        singles += 1;
        for c2 in kill!(rs1[0]) {
            if cw[c2 as usize] == 1 {
                stack.push(c2);
            }
        }
    }
    // 2. excess pruning: drop heaviest rows while rows - cols > keep_excess
    let mut order: Vec<usize> = (0..rows.len()).filter(|&i| alive[i]).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(rows[i].len()));
    let mut pruned = 0u64;
    let mut oi = 0;
    while nr - nc > keep_excess as i64 && oi < order.len() {
        let i = order[oi];
        oi += 1;
        if !alive[i] {
            continue;
        }
        pruned += 1;
        let mut st: Vec<u32> = kill!(i)
            .into_iter()
            .filter(|&c| cw[c as usize] == 1)
            .collect();
        while let Some(c) = st.pop() {
            if cw[c as usize] != 1 {
                continue;
            }
            let rs1 = col_rows(&mut crow, &rows, &alive, c);
            if rs1.len() != 1 {
                continue;
            }
            singles += 1;
            for c2 in kill!(rs1[0]) {
                if cw[c2 as usize] == 1 {
                    st.push(c2);
                }
            }
        }
    }
    // 3. merges. Cost estimate of Lanczos: cols * (2 nnz + 7 rows)
    let est =
        |nr: i64, nc: i64, nnz: i64| -> f64 { nc as f64 * (2.0 * nnz as f64 + 7.0 * nr as f64) };
    let mut inv_cache: HashMap<u64, u64> = HashMap::new();
    let mut heap = std::collections::BinaryHeap::new(); // (Reverse(key), col)
    let key_of = |w: u32, lp: usize| -> i64 {
        if w == 1 {
            i64::MIN
        } else {
            (w as i64 - 1) * (lp as i64 - 2) - lp as i64
        }
    };
    for c in 0..ncol as u32 {
        if cw[c as usize] >= 1 {
            let rr = col_rows(&mut crow, &rows, &alive, c);
            if rr.is_empty() {
                continue;
            }
            let lp = rr.iter().map(|&i| rows[i as usize].len()).min().unwrap();
            heap.push((std::cmp::Reverse(key_of(cw[c as usize], lp)), c));
        }
    }
    let (mut merges, mut merge_mults) = (0u64, 0u64);
    let cm0 = ml.cnt.get();
    while let Some((std::cmp::Reverse(k), c)) = heap.pop() {
        let w = cw[c as usize];
        if w == 0 {
            continue;
        }
        let rr = col_rows(&mut crow, &rows, &alive, c);
        if rr.is_empty() {
            continue;
        }
        let piv = *rr.iter().min_by_key(|&&i| rows[i as usize].len()).unwrap();
        let lp = rows[piv as usize].len();
        let kn = key_of(w, lp);
        if kn != k {
            heap.push((std::cmp::Reverse(kn), c));
            continue;
        }
        if w == 1 {
            singles += 1;
            for c2 in kill!(piv) {
                if cw[c2 as usize] >= 1 {
                    let r2 = col_rows(&mut crow, &rows, &alive, c2);
                    if let Some(l) = r2.iter().map(|&i| rows[i as usize].len()).min() {
                        heap.push((std::cmp::Reverse(key_of(cw[c2 as usize], l)), c2));
                    }
                }
            }
            continue;
        }
        let dnnz = (w as i64 - 1) * (lp as i64 - 2) - lp as i64;
        if est(nr - 1, nc - 1, nnz + dnnz) >= est(nr, nc, nnz) {
            break;
        }
        // perform the merge
        let prow = rows[piv as usize].clone();
        let pv = prow[prow.binary_search_by_key(&c, |e| e.0).unwrap()].1;
        let pinv = match inv_cache.get(&pv) {
            Some(&x) => x,
            None => {
                let x = ml.inv(pv);
                inv_cache.insert(pv, x);
                x
            }
        };
        let (pa, pb) = (av[piv as usize], bv[piv as usize]);
        let mut touched: Vec<u32> = vec![];
        for &i in rr.iter().filter(|&&i| i != piv) {
            let i = i as usize;
            let rv = rows[i][rows[i].binary_search_by_key(&c, |e| e.0).unwrap()].1;
            let fct = ml.mul(rv, pinv);
            let old = std::mem::take(&mut rows[i]);
            let mut nw: Vec<(u32, u64)> = Vec::with_capacity(old.len() + prow.len());
            let (mut x, mut y) = (0, 0);
            while x < old.len() || y < prow.len() {
                if y >= prow.len() || (x < old.len() && old[x].0 < prow[y].0) {
                    nw.push(old[x]);
                    x += 1;
                } else if x >= old.len() || prow[y].0 < old[x].0 {
                    let col = prow[y].0;
                    let val = ml.neg(ml.mul(fct, prow[y].1));
                    if col != c && val != 0 {
                        nw.push((col, val));
                        cw[col as usize] += 1;
                        if cw[col as usize] == 1 {
                            nc += 1;
                        }
                        crow[col as usize].push(i as u32);
                        touched.push(col);
                    }
                    y += 1;
                } else {
                    let col = old[x].0;
                    let val = ml.sub(old[x].1, ml.mul(fct, prow[y].1));
                    if val != 0 {
                        nw.push((col, val));
                    } else {
                        cw[col as usize] -= 1;
                        if cw[col as usize] == 0 {
                            nc -= 1;
                        }
                        touched.push(col);
                    }
                    x += 1;
                    y += 1;
                }
            }
            nnz += nw.len() as i64 - old.len() as i64;
            rows[i] = nw;
            av[i] = ml.sub(av[i], ml.mul(fct, pa));
            bv[i] = ml.sub(bv[i], ml.mul(fct, pb));
        }
        for c2 in kill!(piv) {
            touched.push(c2);
        }
        merges += 1;
        touched.sort_unstable();
        touched.dedup();
        for c2 in touched {
            if cw[c2 as usize] >= 1 {
                let r2 = col_rows(&mut crow, &rows, &alive, c2);
                if let Some(l) = r2.iter().map(|&i| rows[i as usize].len()).min() {
                    heap.push((std::cmp::Reverse(key_of(cw[c2 as usize], l)), c2));
                }
            }
        }
    }
    merge_mults += ml.cnt.get() - cm0;
    let mut m = Mat {
        rows: vec![],
        a: vec![],
        b: vec![],
        ncols: 0,
    };
    for i in 0..rows.len() {
        if alive[i] {
            m.rows.push(std::mem::take(&mut rows[i]));
            m.a.push(av[i]);
            m.b.push(bv[i]);
        }
    }
    compact_cols(&mut m);
    let elim = ml.cnt.get() - c0;
    let c1 = ml.cnt.get();
    let nnzf: usize = m.rows.iter().map(|r| r.len()).sum();
    let (lg, iters, att) = lanczos_log(&m, ml, seed);
    LaOut { method: "sge_lanczos".into(), rows: m.rows.len(), cols: m.ncols, nnz: nnzf, elim_mults: elim, lanczos_mults: ml.cnt.get() - c1, iters, attempts: att, log: lg, wall: t0.elapsed().as_secs_f64(),
        extra: format!("\"rows_in\":{},\"cols_in\":{},\"singletons\":{},\"pruned\":{},\"merges\":{},\"merge_mults\":{},\"keep_excess\":{}", nrow0, ncol, singles, pruned, merges, merge_mults, keep_excess) }
}

fn make_instance<'a>(f: &'a Fq, c: &Curve, target: u32) -> (Inst<'a>, u64) {
    let ts = mix(mix(c.seed, c.k as u64), mix(c.idx as u64, target as u64));
    let mut rng = Rng::new(mix(ts, 1));
    let p = loop {
        let d = jrandom(f, &c.f, &mut rng);
        let p = jmul(f, &c.f, c.h, &d);
        if p.u.d > 0 {
            break p;
        }
    };
    assert!(jmul(f, &c.f, c.ell, &p).u.d == 0, "P not of order ell");
    let k = 1 + rng.below(c.ell - 1);
    let qd = jmul(f, &c.f, k, &p);
    assert!(jvalid(f, &c.f, &p) && jvalid(f, &c.f, &qd));
    (
        Inst {
            f,
            fc: c.f,
            ell: c.ell,
            p,
            qd,
            planted: k,
        },
        ts,
    )
}

fn rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find(|l| l.starts_with("VmHWM")).and_then(|l| {
                l.split_whitespace()
                    .nth(1)
                    .and_then(|x| x.parse::<f64>().ok())
            })
        })
        .map(|kb| kb / 1024.0)
        .unwrap_or(-1.0)
}

fn cmd_solve(path: &str, idx: u32, target: u32, kappa: f64, las: &[String], keep_excess: usize) {
    let c = load_curve(path, idx);
    let f = Fq::new(c.q);
    let t0 = Instant::now();
    let (inst, ts) = make_instance(&f, &c, target);
    f.cnt.set(0);
    let walk_seed = mix(ts, 2);
    let (rs, st) = relation_phase(&inst, kappa, walk_seed);
    let rel_mults = st.fq_setup + st.fq_walk + st.fq_decomp;
    let mut la_json = vec![];
    for la in las {
        let ml = Ml::new(c.ell);
        ml.cnt.set(0);
        let lseed = mix(ts, 3);
        let out = match la.as_str() {
            "bfs" => la_bfs(&rs, &ml, lseed),
            "sge" => la_sge(&rs, &ml, lseed, keep_excess),
            _ => panic!("unknown la {}", la),
        };
        let cv = f.cnt.get();
        let (ok, okp) = match out.log {
            Some(x) => (jmul(&f, &c.f, x, &inst.p) == inst.qd, x == inst.planted),
            None => (false, false),
        };
        let ver = f.cnt.get() - cv;
        f.cnt.set(cv);
        let total = rel_mults + out.elim_mults + out.lanczos_mults + ver;
        la_json.push(format!("{{\"method\":\"{}\",\"rows\":{},\"cols\":{},\"nnz\":{},\"w_row\":{:.4},\"w_col\":{:.4},\"elim_mults\":{},\"lanczos_mults\":{},\"lanczos_iters\":{},\"attempts\":{},\"verify_fq_mults\":{},\"log\":{},\"log_verified\":{},\"log_equals_planted\":{},\"total_mults\":{},\"S\":{:.6},\"wall_s\":{:.3},{}}}",
            out.method, out.rows, out.cols, out.nnz, out.nnz as f64 / out.rows.max(1) as f64, out.nnz as f64 / out.cols.max(1) as f64, out.elim_mults, out.lanczos_mults, out.iters, out.attempts, ver,
            out.log.map(|x| x.to_string()).unwrap_or("null".into()), ok, okp, total, total as f64 / (c.ell as f64).sqrt(), out.wall, out.extra));
    }
    println!("{{\"tool\":\"ic04b\",\"cmd\":\"solve\",\"seed\":{},\"k\":{},\"q\":{},\"curve\":{},\"target\":{},\"kappa\":{},\"F\":{},\"ell\":{},\"N\":{},\"planted\":{},\"rel\":{{\"trials\":{},\"split\":{},\"full\":{},\"p1\":{},\"p2\":{},\"disc3\":{},\"edges\":{},\"comp0_lp_vertices\":{},\"comp0_edges\":{},\"target_relations\":{},\"fq_mults_setup\":{},\"fq_mults_walk\":{},\"fq_mults_decomp\":{},\"rel_mults\":{},\"wall_s\":{:.3}}},\"la\":[{}],\"rss_mb\":{:.1},\"wall_s\":{:.3}}}",
        c.seed, c.k, c.q, idx, target, kappa, rs.fsize, c.ell, c.n, inst.planted, st.trials, st.split, st.full, st.p1, st.p2, st.disc, rs.edges.len(), rs.comp0_v, rs.comp0_e, st.target,
        st.fq_setup, st.fq_walk, st.fq_decomp, rel_mults, st.wall, la_json.join(","), rss_mb(), t0.elapsed().as_secs_f64());
}

/// Pollard rho reference: r-adding walk (r = 20) with Brent cycle detection, same instance, same counted unit.
fn cmd_rho(path: &str, idx: u32, target: u32) {
    let c = load_curve(path, idx);
    let f = Fq::new(c.q);
    let t0 = Instant::now();
    let (inst, ts) = make_instance(&f, &c, target);
    let fc = &c.f;
    let ell = c.ell;
    f.cnt.set(0);
    let mut rng = Rng::new(mix(ts, 4));
    const R: usize = 20;
    let mut mults = vec![];
    for _ in 0..R {
        let (al, be) = (rng.below(ell), rng.below(ell));
        mults.push((
            jadd(
                &f,
                fc,
                &jmul(&f, fc, al, &inst.p),
                &jmul(&f, fc, be, &inst.qd),
            ),
            al,
            be,
        ));
    }
    let mut steps = 0u64;
    let mut restarts = 0u32;
    let res = 'outer: loop {
        let (a0, b0) = (rng.below(ell), rng.below(ell));
        let mut x = jadd(
            &f,
            fc,
            &jmul(&f, fc, a0, &inst.p),
            &jmul(&f, fc, b0, &inst.qd),
        );
        let (mut xa, mut xb) = (a0, b0);
        let (mut s, mut sa, mut sb) = (x, xa, xb);
        let mut pw = 1u64;
        let mut lam = 0u64;
        loop {
            let i = (jhash(&x) % R as u64) as usize;
            x = jadd(&f, fc, &x, &mults[i].0);
            xa = (xa + mults[i].1) % ell;
            xb = (xb + mults[i].2) % ell;
            steps += 1;
            lam += 1;
            if x == s {
                if xb != sb {
                    let num = (sa + ell - xa) % ell;
                    let den = (xb + ell - sb) % ell;
                    break 'outer Some(mulmod(num, powmod(den, ell - 2, ell), ell));
                }
                restarts += 1;
                break;
            }
            if lam == pw {
                s = x;
                sa = xa;
                sb = xb;
                pw *= 2;
                lam = 0;
            }
        }
    };
    let total = f.cnt.get();
    let cv = f.cnt.get();
    let ok = res
        .map(|x| jmul(&f, fc, x, &inst.p) == inst.qd)
        .unwrap_or(false);
    let ver = f.cnt.get() - cv;
    println!("{{\"tool\":\"ic04b\",\"cmd\":\"rho\",\"seed\":{},\"k\":{},\"q\":{},\"curve\":{},\"target\":{},\"ell\":{},\"steps\":{},\"restarts\":{},\"fq_mults\":{},\"verify_fq_mults\":{},\"log\":{},\"log_verified\":{},\"log_equals_planted\":{},\"S\":{:.6},\"steps_per_sqrt_ell\":{:.4},\"wall_s\":{:.3}}}",
        c.seed, c.k, c.q, idx, target, ell, steps, restarts, total, ver, res.map(|x| x.to_string()).unwrap_or("null".into()), ok, res == Some(inst.planted), (total + ver) as f64 / (ell as f64).sqrt(), steps as f64 / (ell as f64).sqrt(), t0.elapsed().as_secs_f64());
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let get =
        |k: &str| -> Option<String> { a.iter().position(|x| x == k).map(|i| a[i + 1].clone()) };
    match a.get(1).map(|s| s.as_str()) {
        Some("curves") => cmd_curves(
            get("--k").unwrap().parse().unwrap(),
            get("--count").unwrap_or("7".into()).parse().unwrap(),
            get("--seed").expect("--seed required").parse().unwrap(),
        ),
        Some("solve") => {
            let las: Vec<String> = get("--la")
                .unwrap_or("bfs,sge".into())
                .split(',')
                .map(|s| s.to_string())
                .collect();
            cmd_solve(
                &get("--curves").unwrap(),
                get("--idx").unwrap().parse().unwrap(),
                get("--target").unwrap().parse().unwrap(),
                get("--kappa").unwrap().parse().unwrap(),
                &las,
                get("--excess").unwrap_or("20".into()).parse().unwrap(),
            )
        }
        Some("rho") => cmd_rho(
            &get("--curves").unwrap(),
            get("--idx").unwrap().parse().unwrap(),
            get("--target").unwrap().parse().unwrap(),
        ),
        _ => {
            eprintln!("usage: curves --k K --count C --seed S | solve --curves F --idx I --target T --kappa X [--la bfs,sge] [--excess E] | rho --curves F --idx I --target T");
            std::process::exit(2);
        }
    }
}
