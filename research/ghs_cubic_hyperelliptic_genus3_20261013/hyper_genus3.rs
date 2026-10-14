// Hyperelliptic genus-3 GHS cover of type m = 3, r = 5 (PROTOCOL.md).
// Builds C : w^2 = h(z) over F_q from E0 : y0^2 = Q(x)(x-r0)(x-r1) over K = F_{q^3}, checks descent,
// Frobenius and the transfer phi = (1 + Fr + Fr^2) o pi0^* with Cantor arithmetic.
// No index calculus, no relation collection, no discrete-log solver.
// Dependency-free: rustc --edition 2021 -O hyper_genus3.rs ; ./hyper_genus3 <q> <census jsonl>
// One JSON object per line on stdout. Single-threaded, deterministic (seed 20261013).
use std::collections::{HashMap, HashSet};

const SEED: u64 = 20261013;
struct Rng(u64);
impl Rng {
    fn new(tag: u64) -> Self { Rng(SEED ^ tag.wrapping_mul(0x9E3779B97F4A7C15)) }
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: u64) -> u64 { self.next() % n }
}

// ---------------------------------------------------------------- finite fields, log representation
// 0 is zero, v >= 1 is g^(v-1). Index = coefficients (base q, little endian) in the polynomial basis.
struct F { q: u32, size: u32, m: u32, zech: Vec<u32>, idx_of: Vec<u32>, el_of: Vec<u32>, modulus: Vec<u32> }
fn digits(q: u32, n: usize, mut i: u32) -> Vec<u32> { let mut v = vec![0; n]; for k in 0..n { v[k] = i % q; i /= q; } v }
fn index(q: u32, v: &[u32]) -> u32 { let mut i = 0; for k in (0..v.len()).rev() { i = i * q + v[k]; } i }
impl F {
    fn finish(q: u32, n: u32, exp_idx: Vec<u32>, modulus: Vec<u32>) -> F {
        let size = q.pow(n); let m = size - 1;
        let mut el_of = vec![0u32; size as usize]; let mut idx_of = vec![0u32; size as usize];
        for k in 0..m as usize { el_of[exp_idx[k] as usize] = k as u32 + 1; idx_of[k + 1] = exp_idx[k]; }
        let mut zech = vec![0u32; m as usize];
        for k in 0..m as usize { let i = exp_idx[k]; let d0 = i % q; let j = i - d0 + (d0 + 1) % q; zech[k] = el_of[j as usize]; }
        F { q, size, m, zech, idx_of, el_of, modulus }
    }
    // L = F_{q^n}: first monic x^n + sum c_i x^i (c in index order, c0 != 0) with x primitive
    fn build_l(q: u32, n: u32) -> F {
        let size = q.pow(n); let m = size - 1; let nn = n as usize;
        for cand in 0..size { let c = digits(q, nn, cand); if c[0] == 0 { continue; }
            let mut exp_idx = vec![0u32; m as usize]; let mut cur = vec![0u32; nn]; cur[0] = 1; let mut ok = true;
            for k in 0..m as usize { let ci = index(q, &cur); if k > 0 && ci == 1 { ok = false; break; } exp_idx[k] = ci;
                let top = cur[nn - 1]; for i in (1..nn).rev() { cur[i] = (cur[i - 1] + (q - c[i]) * top) % q; } cur[0] = ((q - c[0]) * top) % q; }
            if !ok || index(q, &cur) != 1 { continue; }
            return F::finish(q, n, exp_idx, c);
        }
        panic!("no primitive polynomial");
    }
    // verbatim field search of research/ghs_cubic_cover_census_20261012/cover_census.rs (same representation)
    fn build_census(q: u32, n: u32) -> F {
        let size = q.pow(n); let m = size - 1; let nn = n as usize;
        for cand in 0..size { let md = digits(q, nn, cand);
            for s in 0..q {
                let mut exp_idx = vec![0u32; m as usize]; let mut cur = vec![0u32; nn]; cur[0] = 1; let mut ok = true;
                for k in 0..m as usize { let ci = index(q, &cur); if ci == 0 || (k > 0 && ci == 1) { ok = false; break; } exp_idx[k] = ci;
                    let mut nx = vec![0u32; nn]; let top = cur[nn - 1];
                    for i in 0..nn { let sh = if i > 0 { cur[i - 1] } else { 0 }; nx[i] = (sh + (q - md[i]) * top + s * cur[i]) % q; } cur = nx; }
                if !ok || index(q, &cur) != 1 { continue; }
                return F::finish(q, n, exp_idx, md);
            } }
        panic!("no field");
    }
    #[inline] fn mul(&self, a: u32, b: u32) -> u32 { if a == 0 || b == 0 { 0 } else { let s = a + b - 2; (if s >= self.m { s - self.m } else { s }) + 1 } }
    #[inline] fn inv(&self, a: u32) -> u32 { assert!(a != 0, "inverse of zero"); (self.m - (a - 1)) % self.m + 1 }
    #[inline] fn div(&self, a: u32, b: u32) -> u32 { self.mul(a, self.inv(b)) }
    #[inline] fn neg(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1 + self.m / 2) % self.m + 1 } }
    #[inline] fn add(&self, a: u32, b: u32) -> u32 { if a == 0 { return b; } if b == 0 { return a; }
        let dlt = if b >= a { b - a } else { b + self.m - a }; let z = self.zech[dlt as usize]; if z == 0 { 0 } else { let s = a + z - 2; (if s >= self.m { s - self.m } else { s }) + 1 } }
    #[inline] fn sub(&self, a: u32, b: u32) -> u32 { self.add(a, self.neg(b)) }
    #[inline] fn frob(&self, a: u32) -> u32 { if a == 0 { 0 } else { (((a - 1) as u64 * self.q as u64) % self.m as u64) as u32 + 1 } }
    fn pow(&self, a: u32, e: u64) -> u32 { if e == 0 { return 1; } if a == 0 { 0 } else { (((a - 1) as u64 * (e % self.m as u64)) % self.m as u64) as u32 + 1 } }
    fn int(&self, c: i64) -> u32 { self.el_of[c.rem_euclid(self.q as i64) as usize] }
    fn to_int(&self, v: u32) -> u32 { self.idx_of[v as usize] } // prime-field elements only
    fn step(&self, k: u32) -> u32 { self.m / (self.q.pow(k) - 1) }
    fn in_sub(&self, v: u32, k: u32) -> bool { v == 0 || (v - 1) % self.step(k) == 0 }
    fn chi(&self, v: u32, k: u32) -> i32 { if v == 0 { return 0; } let s = self.step(k); assert!((v - 1) % s == 0, "chi outside subfield"); if ((v - 1) / s) % 2 == 0 { 1 } else { -1 } }
    fn sqrt(&self, v: u32) -> Option<u32> { if v == 0 { Some(0) } else if (v - 1) % 2 == 0 { Some((v - 1) / 2 + 1) } else { None } }
    fn elems(&self, k: u32) -> Vec<u32> { let s = self.step(k); let mut v = vec![0u32]; for t in 0..(self.q.pow(k) - 1) { v.push(t * s + 1); } v }
}

// ---------------------------------------------------------------- polynomials over L, low -> high, trimmed
type P = Vec<u32>;
fn tr(mut a: P) -> P { while let Some(&0) = a.last() { a.pop(); } a }
fn dg(a: &P) -> i32 { a.len() as i32 - 1 }
fn padd(f: &F, a: &P, b: &P) -> P { let n = a.len().max(b.len()); tr((0..n).map(|i| f.add(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0))).collect()) }
fn pneg(f: &F, a: &P) -> P { a.iter().map(|&c| f.neg(c)).collect() }
fn psub(f: &F, a: &P, b: &P) -> P { padd(f, a, &pneg(f, b)) }
fn pmul(f: &F, a: &P, b: &P) -> P { if a.is_empty() || b.is_empty() { return vec![]; } let mut r = vec![0u32; a.len() + b.len() - 1];
    for i in 0..a.len() { if a[i] == 0 { continue; } for j in 0..b.len() { if b[j] != 0 { r[i + j] = f.add(r[i + j], f.mul(a[i], b[j])); } } } tr(r) }
fn psc(f: &F, a: &P, c: u32) -> P { tr(a.iter().map(|&x| f.mul(x, c)).collect()) }
fn pdivrem(f: &F, a: &P, b: &P) -> (P, P) { assert!(!b.is_empty(), "division by zero polynomial"); let mut r = a.clone(); if r.len() < b.len() { return (vec![], r); }
    let mut qv = vec![0u32; r.len() - b.len() + 1]; let li = f.inv(*b.last().unwrap());
    while r.len() >= b.len() { let k = r.len() - b.len(); let c = f.mul(*r.last().unwrap(), li); qv[k] = c;
        for j in 0..b.len() { r[k + j] = f.sub(r[k + j], f.mul(c, b[j])); } r = tr(r); if r.is_empty() { break; } }
    (tr(qv), r) }
fn pmod(f: &F, a: &P, b: &P) -> P { pdivrem(f, a, b).1 }
fn pexact(f: &F, a: &P, b: &P) -> P { let (qq, r) = pdivrem(f, a, b); assert!(r.is_empty(), "inexact division"); qq }
fn pmonic(f: &F, a: &P) -> P { if a.is_empty() { return vec![]; } psc(f, a, f.inv(*a.last().unwrap())) }
fn pgcd(f: &F, a: &P, b: &P) -> P { let (mut x, mut y) = (a.clone(), b.clone()); while !y.is_empty() { let r = pmod(f, &x, &y); x = y; y = r; } pmonic(f, &x) }
fn pxgcd(f: &F, a: &P, b: &P) -> (P, P, P) { // s a + t b = g, g monic
    let (mut r0, mut r1) = (a.clone(), b.clone()); let (mut s0, mut s1) = (vec![1u32], vec![]); let (mut t0, mut t1) = (vec![], vec![1u32]);
    while !r1.is_empty() { let (qq, r) = pdivrem(f, &r0, &r1); r0 = r1; r1 = r;
        let s2 = psub(f, &s0, &pmul(f, &qq, &s1)); s0 = s1; s1 = s2; let t2 = psub(f, &t0, &pmul(f, &qq, &t1)); t0 = t1; t1 = t2; }
    let li = f.inv(*r0.last().unwrap()); (psc(f, &r0, li), psc(f, &s0, li), psc(f, &t0, li)) }
fn pev(f: &F, a: &P, x: u32) -> u32 { let mut r = 0; for &c in a.iter().rev() { r = f.add(f.mul(r, x), c); } r }
fn pder(f: &F, a: &P) -> P { tr((1..a.len()).map(|i| f.mul(f.int(i as i64), a[i])).collect()) }
fn p_in(f: &F, a: &P, k: u32) -> bool { a.iter().all(|&c| f.in_sub(c, k)) }
fn plin(f: &F, r: u32) -> P { vec![f.neg(r), 1] }
fn pints(f: &F, a: &P) -> String { format!("[{}]", a.iter().map(|&c| f.to_int(c).to_string()).collect::<Vec<_>>().join(",")) }

// ---------------------------------------------------------------- Cantor arithmetic, w^2 = h(z), deg h = 2g+1
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Dv { u: P, v: P }
fn dzero() -> Dv { Dv { u: vec![1], v: vec![] } }
fn dneg(f: &F, d: &Dv) -> Dv { Dv { u: d.u.clone(), v: pneg(f, &d.v) } }
fn dpt(f: &F, z: u32, w: u32) -> Dv { Dv { u: plin(f, z), v: tr(vec![w]) } }
fn dred(f: &F, h: &P, g: i32, u: P, v: P) -> Dv {
    let mut u = pmonic(f, &u); let mut v = pmod(f, &v, &u);
    while dg(&u) > g { let num = psub(f, h, &pmul(f, &v, &v)); u = pmonic(f, &pexact(f, &num, &u)); v = pmod(f, &pneg(f, &v), &u); }
    if dg(&u) <= 0 { return dzero(); } Dv { u, v } }
fn dadd(f: &F, h: &P, g: i32, a: &Dv, b: &Dv) -> Dv {
    let (d0, e1, e2) = pxgcd(f, &a.u, &b.u); let vs = padd(f, &a.v, &b.v);
    let (d, c1, c2) = pxgcd(f, &d0, &vs);
    let (s1, s2, s3) = (pmul(f, &c1, &e1), pmul(f, &c1, &e2), c2);
    let u = pexact(f, &pmul(f, &a.u, &b.u), &pmul(f, &d, &d));
    let t = padd(f, &padd(f, &pmul(f, &pmul(f, &s1, &a.u), &b.v), &pmul(f, &pmul(f, &s2, &b.u), &a.v)), &pmul(f, &s3, &padd(f, &pmul(f, &a.v, &b.v), h)));
    let v = pexact(f, &t, &d); dred(f, h, g, u, v) }
fn dmul(f: &F, h: &P, g: i32, mut n: u64, a: &Dv) -> Dv { let mut r = dzero(); let mut b = a.clone(); while n > 0 { if n & 1 == 1 { r = dadd(f, h, g, &r, &b); } b = dadd(f, h, g, &b, &b); n >>= 1; } r }
fn dv_in(f: &F, d: &Dv) -> bool { p_in(f, &d.u, 1) && p_in(f, &d.v, 1) }
fn dv_ok(f: &F, h: &P, d: &Dv) -> bool { pmod(f, &psub(f, &pmul(f, &d.v, &d.v), h), &d.u).is_empty() && dg(&d.v) < dg(&d.u).max(0) || d == &dzero() }

// ---------------------------------------------------------------- elliptic curve Y^2 = S^3 + a2 S^2 + a4 S + a6
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum WP { O, A(u32, u32) }
struct Wc { a2: u32, a4: u32, a6: u32 }
fn wrhs(f: &F, w: &Wc, s: u32) -> u32 { let s2 = f.mul(s, s); f.add(f.add(f.add(f.mul(s2, s), f.mul(w.a2, s2)), f.mul(w.a4, s)), w.a6) }
fn wneg(f: &F, p: WP) -> WP { match p { WP::O => WP::O, WP::A(x, y) => WP::A(x, f.neg(y)) } }
fn wadd(f: &F, w: &Wc, p: WP, q: WP) -> WP {
    let (x1, y1, x2, y2) = match (p, q) { (WP::O, _) => return q, (_, WP::O) => return p, (WP::A(a, b), WP::A(c, d)) => (a, b, c, d) };
    let lam = if x1 == x2 { if f.add(y1, y2) == 0 { return WP::O; }
        f.div(f.add(f.add(f.mul(f.int(3), f.mul(x1, x1)), f.mul(f.int(2), f.mul(w.a2, x1))), w.a4), f.mul(f.int(2), y1)) } else { f.div(f.sub(y2, y1), f.sub(x2, x1)) };
    let x3 = f.sub(f.sub(f.sub(f.mul(lam, lam), w.a2), x1), x2); let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1); WP::A(x3, y3) }
fn wmul(f: &F, w: &Wc, mut n: u64, p: WP) -> WP { let mut r = WP::O; let mut b = p; while n > 0 { if n & 1 == 1 { r = wadd(f, w, r, b); } b = wadd(f, w, b, b); n >>= 1; } r }

// ---------------------------------------------------------------- 3x3 linear algebra over L
type M3 = [[u32; 3]; 3];
fn mv(f: &F, m: &M3, v: &[u32; 3]) -> [u32; 3] { let mut r = [0u32; 3]; for i in 0..3 { for j in 0..3 { r[i] = f.add(r[i], f.mul(m[i][j], v[j])); } } r }
fn inv3(f: &F, a: &M3) -> Option<M3> {
    let mut c = [[0u32; 3]; 3];
    for i in 0..3 { for j in 0..3 { c[i][j] = f.sub(f.mul(a[(j + 1) % 3][(i + 1) % 3], a[(j + 2) % 3][(i + 2) % 3]), f.mul(a[(j + 1) % 3][(i + 2) % 3], a[(j + 2) % 3][(i + 1) % 3])); } }
    let det = f.add(f.add(f.mul(a[0][0], c[0][0]), f.mul(a[0][1], c[1][0])), f.mul(a[0][2], c[2][0])); if det == 0 { return None; }
    let di = f.inv(det); for i in 0..3 { for j in 0..3 { c[i][j] = f.mul(c[i][j], di); } } Some(c) }
fn bil(f: &F, g: &M3, a: &[u32; 3], b: &[u32; 3]) -> u32 { let gb = mv(f, g, b); f.add(f.add(f.mul(a[0], gb[0]), f.mul(a[1], gb[1])), f.mul(a[2], gb[2])) }
fn parallel(f: &F, a: &[u32; 3], b: &[u32; 3]) -> bool { (0..3).all(|i| { let j = (i + 1) % 3; f.mul(a[i], b[j]) == f.mul(a[j], b[i]) }) }

// ---------------------------------------------------------------- the cover
struct Inst { split: bool, c0: u32, c1: u32, ab: Option<(u32, u32)>, r: [u32; 3], qp: P, fi: [P; 3], np: P,
    mm: M3, g: M3, ps: [u32; 3], v0: [u32; 3], v1: [u32; 3], binv: M3, xp: [P; 3], xn: P, xd: P, h: P, gen: i32, pts_conic: u32, pts_over_a: u32 }
impl Inst {
    fn special(&self, f: &F, x: u32) -> bool { pev(f, &self.qp, x) == 0 || self.r.contains(&x) }
    fn fwd(&self, f: &F, pt: [u32; 4]) -> Option<(u32, u32)> {
        let [x, y0, y1, y2] = pt; let qx = pev(f, &self.qp, x); if qx == 0 || self.r.contains(&x) || y0 == 0 || y1 == 0 || y2 == 0 { return None; }
        let u01 = f.div(f.mul(y0, y1), f.mul(qx, f.sub(x, self.r[1]))); let u12 = f.div(f.mul(y1, y2), f.mul(qx, f.sub(x, self.r[2]))); let u20 = f.div(f.mul(y2, y0), f.mul(qx, f.sub(x, self.r[0])));
        let xx = [f.mul(u12, u20), f.mul(u20, u01), f.mul(u01, u12)]; let zz = mv(f, &self.mm, &xx); let co = mv(f, &self.binv, &zz);
        if co[2] == 0 { return None; } let z = f.div(co[1], co[2]);
        let w = f.div(f.mul(f.mul(f.mul(y0, y1), y2), pev(f, &self.xd, z)), f.mul(qx, pev(f, &self.np, x))); Some((z, w)) }
    // addendum: a generic point lying over z = infinity (the conic base point P*): returns s = lim w/z^4
    fn fwd_inf(&self, f: &F, pt: [u32; 4]) -> Option<u32> {
        let [x, y0, y1, y2] = pt; let qx = pev(f, &self.qp, x); if qx == 0 || self.r.contains(&x) || y0 == 0 || y1 == 0 || y2 == 0 || dg(&self.xd) != 4 { return None; }
        let co = mv(f, &self.binv, &mv(f, &self.mm, &self.xtrue(f, pt))); if co[2] != 0 || co[1] != 0 { return None; }
        Some(f.div(f.mul(f.mul(f.mul(y0, y1), y2), self.xd[4]), f.mul(qx, pev(f, &self.np, x)))) }
    fn xtrue(&self, f: &F, pt: [u32; 4]) -> [u32; 3] { let [x, y0, y1, y2] = pt; let qx = pev(f, &self.qp, x);
        let u01 = f.div(f.mul(y0, y1), f.mul(qx, f.sub(x, self.r[1]))); let u12 = f.div(f.mul(y1, y2), f.mul(qx, f.sub(x, self.r[2]))); let u20 = f.div(f.mul(y2, y0), f.mul(qx, f.sub(x, self.r[0])));
        [f.mul(u12, u20), f.mul(u20, u01), f.mul(u01, u12)] }
    // (z, w) -> (x, y0, y1, y2): the unique sign pattern matching w and the conic point
    fn bwd(&self, f: &F, z: u32, w: u32) -> Result<[u32; 4], &'static str> {
        let xdz = pev(f, &self.xd, z); if xdz == 0 { return Err("x=inf"); } let x = f.div(pev(f, &self.xn, z), xdz); if self.special(f, x) { return Err("special x"); }
        let mut s = [0u32; 3]; for i in 0..3 { s[i] = match f.sqrt(pev(f, &self.fi[i], x)) { Some(v) => v, None => return Err("y not in L") }; }
        let target = f.div(f.mul(f.mul(w, pev(f, &self.qp, x)), pev(f, &self.np, x)), xdz);
        let xz = [pev(f, &self.xp[0], z), pev(f, &self.xp[1], z), pev(f, &self.xp[2], z)]; let mut found = vec![];
        for e in 0..8u32 { let y: Vec<u32> = (0..3).map(|i| if e >> i & 1 == 1 { f.neg(s[i]) } else { s[i] }).collect();
            if f.mul(f.mul(y[0], y[1]), y[2]) != target { continue; } let pt = [x, y[0], y[1], y[2]];
            if parallel(f, &self.xtrue(f, pt), &xz) { found.push(pt); } }
        if found.len() == 1 { Ok(found[0]) } else { Err("sign match not unique") } }
}
fn frc(f: &F, p: [u32; 4]) -> [u32; 4] { [f.frob(p[0]), f.frob(p[3]), f.frob(p[1]), f.frob(p[2])] }

fn build(f: &F, split: bool, c0: u32, c1: u32, ab: Option<(u32, u32)>, p: u32) -> Result<Inst, String> {
    let q = f.q; let r = [p, f.frob(p), f.frob(f.frob(p))];
    let qp: P = vec![c0, c1, 1];
    let fi = [pmul(f, &qp, &pmul(f, &plin(f, r[0]), &plin(f, r[1]))), pmul(f, &qp, &pmul(f, &plin(f, r[1]), &plin(f, r[2]))), pmul(f, &qp, &pmul(f, &plin(f, r[2]), &plin(f, r[0])))];
    let np = pmul(f, &plin(f, r[0]), &pmul(f, &plin(f, r[1]), &plin(f, r[2])));
    if !p_in(f, &np, 1) { return Err("N not over F_q".into()); }
    let gam = (q * q * q + 1) + 1; let th = [1u32, gam, f.mul(gam, gam)];
    let mut mm = [[0u32; 3]; 3]; for j in 0..3 { let mut t = th[j]; for i in 0..3 { mm[j][i] = t; t = f.frob(t); } }
    let minv = inv3(f, &mm).ok_or("M singular")?;
    let c = [f.sub(r[2], r[0]), f.sub(r[0], r[1]), f.sub(r[1], r[2])];
    let mut g = [[0u32; 3]; 3]; for j in 0..3 { for k in 0..3 { for i in 0..3 { g[j][k] = f.add(g[j][k], f.mul(f.mul(minv[i][j], c[i]), minv[i][k])); } } }
    if !(0..9).all(|t| f.in_sub(g[t / 3][t % 3], 1)) { return Err("P1 fail: conic not over F_q".into()); }
    let xof = |zv: &[u32; 3]| -> Option<u32> { let xx = mv(f, &minv, zv); let (a0, a2) = (f.mul(xx[0], xx[0]), f.mul(xx[2], xx[2])); let d = f.sub(a0, a2);
        if d == 0 { None } else { Some(f.div(f.sub(f.mul(r[0], a0), f.mul(r[1], a2)), d)) } };
    // F_q points of the conic, normalised, in enumeration order
    let mut cpts: Vec<[u32; 3]> = vec![]; for t in 0..(q * q + q + 1) { let v = if t < q * q { [1, t % q, t / q] } else if t < q * q + q { [0, 1, t - q * q] } else { [0, 0, 1] };
        let v = [f.int(v[0] as i64), f.int(v[1] as i64), f.int(v[2] as i64)]; if bil(f, &g, &v, &v) == 0 { cpts.push(v); } }
    let over_a: Vec<[u32; 3]> = match ab { Some((a, _)) => cpts.iter().filter(|z| xof(z) == Some(a)).cloned().collect(), None => vec![] };
    let ps = if split { if over_a.len() != 1 { return Err(format!("P1 fail: {} F_q points over a", over_a.len())); } over_a[0] } else { cpts[0] };
    let l = mv(f, &g, &ps);
    let dot = |a: &[u32; 3], b: &[u32; 3]| f.add(f.add(f.mul(a[0], b[0]), f.mul(a[1], b[1])), f.mul(a[2], b[2]));
    let mut v0 = [0u32; 3]; for t in 1..q * q * q { let d = digits(q, 3, t); let v = [f.int(d[0] as i64), f.int(d[1] as i64), f.int(d[2] as i64)]; if dot(&l, &v) == 0 && !parallel(f, &v, &ps) { v0 = v; break; } }
    let mut v1 = [0u32; 3]; for k in 0..3 { let mut v = [0u32; 3]; v[k] = 1; if dot(&l, &v) != 0 { v1 = v; break; } }
    let (qv0, bv01, qv1, b1) = (bil(f, &g, &v0, &v0), bil(f, &g, &v0, &v1), bil(f, &g, &v1, &v1), bil(f, &g, &ps, &v1));
    let two = f.int(2); let mut zp: [P; 3] = [vec![], vec![], vec![]];
    for j in 0..3 { zp[j] = tr(vec![f.sub(f.mul(qv1, ps[j]), f.mul(f.mul(two, b1), v1[j])), f.sub(f.mul(f.mul(two, bv01), ps[j]), f.mul(f.mul(two, b1), v0[j])), f.mul(qv0, ps[j])]); }
    let mut xp: [P; 3] = [vec![], vec![], vec![]]; for i in 0..3 { for j in 0..3 { xp[i] = padd(f, &xp[i], &psc(f, &zp[j], minv[i][j])); } }
    let (s0, s2) = (pmul(f, &xp[0], &xp[0]), pmul(f, &xp[2], &xp[2]));
    let num = psub(f, &psc(f, &s0, r[0]), &psc(f, &s2, r[1])); let den = psub(f, &s0, &s2);
    let gg = pgcd(f, &num, &den); let (mut xn, mut xd) = (pexact(f, &num, &gg), pexact(f, &den, &gg));
    let li = f.inv(*xd.last().unwrap()); xn = psc(f, &xn, li); xd = psc(f, &xd, li);
    let h = padd(f, &padd(f, &pmul(f, &xn, &xn), &psc(f, &pmul(f, &xn, &xd), c1)), &psc(f, &pmul(f, &xd, &xd), c0));
    let bm = [[ps[0], v0[0], v1[0]], [ps[1], v0[1], v1[1]], [ps[2], v0[2], v1[2]]];
    let binv = inv3(f, &bm).ok_or("basis singular")?;
    let pts_conic = cpts.len() as u32; let pts_over_a = over_a.len() as u32;
    Ok(Inst { split, c0, c1, ab, r, qp, fi, np, mm, g, ps, v0, v1, binv, xp, xn, xd, h, gen: 3, pts_conic, pts_over_a })
}

// ---------------------------------------------------------------- real model (addendum 2026-10-14): y^2 = f, f monic deg 8, genus 3
// class (u, v, n) = div(u, v) + n inf+ + (3 - deg u - n) inf- - Dinf, Dinf = 2 inf+ + inf-; reduced iff deg u <= 3, 0 <= n <= 3 - deg u
struct Rm { f: P, vp: P }
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Rd { u: P, v: P, n: i32 }
fn rzero() -> Rd { Rd { u: vec![1], v: vec![], n: 2 } }
// one step with v adapted to +V+ (plus) or -V+ (minus). Returns u', v' and the pole order e of y - va at the
// other point at infinity (inf- for plus, inf+ for minus): plus gives n' = n + deg u - e, minus gives m' = m + deg u - e.
fn rstep(fl: &F, rm: &Rm, u: &P, v: &P, plus: bool) -> (P, P, i32) {
    let va = if plus { psub(fl, &rm.vp, &pmod(fl, &psub(fl, &rm.vp, v), u)) } else { padd(fl, &pneg(fl, &rm.vp), &pmod(fl, &padd(fl, &rm.vp, v), u)) };
    let up = pmonic(fl, &pexact(fl, &psub(fl, &rm.f, &pmul(fl, &va, &va)), u)); let vn = pmod(fl, &pneg(fl, &va), &up);
    let e = if plus { dg(&padd(fl, &rm.vp, &va)) } else { dg(&psub(fl, &rm.vp, &va)) }; (up, vn, e) }
fn rred(fl: &F, rm: &Rm, u: P, v: P, n: i32) -> Rd {
    let mut u = pmonic(fl, &u); let mut v = pmod(fl, &v, &u); let mut n = n;
    while dg(&u) > 3 { let du = dg(&u); let (a, b, e) = rstep(fl, rm, &u, &v, true); u = a; v = b; n += du - e; }
    let mut guard = 0;
    loop { guard += 1; assert!(guard < 100, "real-model reduction does not terminate");
        let m = 3 - dg(&u) - n;
        if n < 0 { let du = dg(&u); let (a, b, e) = rstep(fl, rm, &u, &v, false); u = a; v = b; let m2 = m + du - e; n = 3 - dg(&u) - m2; }
        else if m < 0 { let du = dg(&u); let (a, b, e) = rstep(fl, rm, &u, &v, true); u = a; v = b; n += du - e; }
        else { break; } }
    if dg(&u) <= 0 { u = vec![1]; v = vec![]; }
    Rd { u, v, n } }
fn radd(fl: &F, rm: &Rm, a: &Rd, b: &Rd) -> Rd {
    let (d0, e1, e2) = pxgcd(fl, &a.u, &b.u); let vs = padd(fl, &a.v, &b.v); let (d, c1, c2) = pxgcd(fl, &d0, &vs);
    let (s1, s2, s3) = (pmul(fl, &c1, &e1), pmul(fl, &c1, &e2), c2);
    let u = pexact(fl, &pmul(fl, &a.u, &b.u), &pmul(fl, &d, &d));
    let t = padd(fl, &padd(fl, &pmul(fl, &pmul(fl, &s1, &a.u), &b.v), &pmul(fl, &pmul(fl, &s2, &b.u), &a.v)), &pmul(fl, &s3, &padd(fl, &pmul(fl, &a.v, &b.v), &rm.f)));
    let v = pexact(fl, &t, &d); rred(fl, rm, u, v, a.n + b.n + dg(&d) - 2) }
fn rneg(fl: &F, rm: &Rm, a: &Rd) -> Rd { rred(fl, rm, a.u.clone(), pneg(fl, &a.v), 4 - dg(&a.u) - a.n) }
fn rmul(fl: &F, rm: &Rm, mut n: u64, a: &Rd) -> Rd { let mut r = rzero(); let mut b = a.clone(); while n > 0 { if n & 1 == 1 { r = radd(fl, rm, &r, &b); } b = radd(fl, rm, &b, &b); n >>= 1; } r }
fn rvalid(fl: &F, rm: &Rm, d: &Rd) -> bool { d.u.last() == Some(&1) && dg(&d.u) <= 3 && d.n >= 0 && d.n <= 3 - dg(&d.u) && dg(&d.v) < dg(&d.u)
    && pmod(fl, &psub(fl, &pmul(fl, &d.v, &d.v), &rm.f), &d.u).is_empty() && p_in(fl, &d.u, 1) && p_in(fl, &d.v, 1) }
// roots over L of u in F_q[z], deg <= 3, with multiplicities
fn roots3(f: &F, u: &P, kel: &[u32]) -> Vec<(u32, u64)> {
    let mut rem = u.clone(); let mut roots = vec![];
    for c in 0..f.q { let cv = f.int(c as i64); let mut mlt = 0; while dg(&rem) >= 1 && pev(f, &rem, cv) == 0 { rem = pexact(f, &rem, &plin(f, cv)); mlt += 1; } if mlt > 0 { roots.push((cv, mlt)); } }
    if dg(&rem) == 2 { let (c0r, c1r) = (rem[0], rem[1]); let disc = f.sub(f.mul(c1r, c1r), f.mul(f.int(4), c0r)); let sd = f.sqrt(disc).unwrap(); let i2 = f.inv(f.int(2));
        roots.push((f.mul(f.sub(sd, c1r), i2), 1)); roots.push((f.mul(f.sub(f.neg(sd), c1r), i2), 1)); }
    else if dg(&rem) == 3 { let rho = *kel.iter().find(|&&t| pev(f, &rem, t) == 0).unwrap(); roots.push((rho, 1)); roots.push((f.frob(rho), 1)); roots.push((f.frob(f.frob(rho)), 1)); }
    roots }

// ---------------------------------------------------------------- helpers for output
struct Iso { cen2l: Vec<u32>, l2cen: HashMap<u32, u32> }
fn kc(iso: &Iso, q: u32, v: u32) -> String { let i = iso.l2cen[&v]; let d = digits(q, 3, i); format!("[{},{},{}]", d[0], d[1], d[2]) }
fn is_prime(n: u64) -> bool { if n < 2 { return false; } let mut d = 2; while d * d <= n { if n % d == 0 { return false; } d += 1; } true }
fn cof(n: u64) -> Option<(u64, u64)> { for h in [4u64, 8] { if n % h == 0 && is_prime(n / h) && n / h >= 5 { return Some((h, n / h)); } } None }
fn jnum(line: &str, key: &str) -> Option<i64> { let k = format!("\"{}\":", key); let i = line.find(&k)? + k.len(); let s: String = line[i..].chars().take_while(|c| c.is_ascii_digit() || *c == '-').collect(); s.parse().ok() }

// census lookup of the K-isomorphism class of y^2 = x^3 + A x + B (census field, its index basis)
fn census_lookup(cf: &F, a_idx: u32, b_idx: u32, recs: &[(u32, u32, i64, i64, bool, String)]) -> Vec<(u32, u32, i64, i64, bool, String)> {
    let (a, b) = (cf.el_of[a_idx as usize], cf.el_of[b_idx as usize]); let mut cls = HashSet::new();
    for u in 1..=cf.m { cls.insert((cf.idx_of[cf.mul(a, cf.pow(u, 4)) as usize], cf.idx_of[cf.mul(b, cf.pow(u, 6)) as usize])); }
    recs.iter().filter(|r| cls.contains(&(r.0, r.1))).cloned().collect() }

fn main() {
    let args: Vec<String> = std::env::args().collect(); let q: u32 = args[1].parse().unwrap(); let census_path = &args[2];
    let t0 = std::time::Instant::now();
    let f = F::build_l(q, 6); let cf = F::build_census(q, 3); let q3 = q * q * q;
    // isomorphism census K -> K inside L: t |-> beta, beta the first root (log order) of t^3 + md2 t^2 + md1 t + md0
    let md = &cf.modulus; let kel = f.elems(3);
    let beta = *kel.iter().find(|&&b| b != 0 && { let b2 = f.mul(b, b); f.add(f.add(f.add(f.mul(b2, b), f.mul(f.int(md[2] as i64), b2)), f.mul(f.int(md[1] as i64), b)), f.int(md[0] as i64)) == 0 }).expect("no root");
    let mut cen2l = vec![0u32; q3 as usize]; let mut l2cen = HashMap::new();
    for i in 0..q3 { let d = digits(q, 3, i); let v = f.add(f.add(f.int(d[0] as i64), f.mul(f.int(d[1] as i64), beta)), f.mul(f.int(d[2] as i64), f.mul(beta, beta))); cen2l[i as usize] = v; l2cen.insert(v, i); }
    assert_eq!(l2cen.len(), q3 as usize);
    // the map is a ring isomorphism: check on all products of a sample
    let mut iso_ok = true; for i in (0..q3).step_by(7) { for j in (0..q3).step_by(11) { let p = cf.idx_of[cf.mul(cf.el_of[i as usize], cf.el_of[j as usize]) as usize];
        let s = cf.idx_of[cf.add(cf.el_of[i as usize], cf.el_of[j as usize]) as usize];
        if f.mul(cen2l[i as usize], cen2l[j as usize]) != cen2l[p as usize] || f.add(cen2l[i as usize], cen2l[j as usize]) != cen2l[s as usize] { iso_ok = false; } } }
    let iso = Iso { cen2l, l2cen };
    println!("{{\"rec\":\"field\",\"q\":{},\"L_modulus_low\":[{}],\"census_K_modulus_low\":[{}],\"iso_beta_log\":{},\"iso_ring_check\":{}}}", q,
        f.modulus.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","), md.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","), beta - 1, iso_ok);
    // census records
    let mut recs = vec![]; if let Ok(txt) = std::fs::read_to_string(census_path) { for line in txt.lines() {
        if !line.contains("\"rec\":\"curve\"") || jnum(line, "q") != Some(q as i64) { continue; }
        let role = if line.contains("\"role\":\"class\"") { "class" } else { "sample" };
        recs.push((jnum(line, "A_idx").unwrap() as u32, jnum(line, "B_idx").unwrap() as u32, jnum(line, "order").unwrap(), jnum(line, "min_genus").unwrap(),
            line.contains("{\"m\":3,\"r\":5,\"hyperelliptic\":true"), role.to_string())); } }
    // Q lists
    let mut qs: Vec<(bool, u32, u32, Option<(u32, u32)>, String)> = vec![];
    for a in 0..q { for b in a + 1..q { qs.push((true, f.int((a * b) as i64), f.int(-((a + b) as i64)), Some((f.int(a as i64), f.int(b as i64))), format!("[{},{}]", a, b))); } }
    for c1 in 0..q { for c0 in 0..q { let disc = f.sub(f.mul(f.int(c1 as i64), f.int(c1 as i64)), f.mul(f.int(4), f.int(c0 as i64))); if f.chi(disc, 1) == -1 { qs.push((false, f.int(c0 as i64), f.int(c1 as i64), None, format!("[{},{}]", c1, c0))); } } }
    let chiq: Vec<Vec<i8>> = qs.iter().map(|(_, c0, c1, _, _)| kel.iter().map(|&x| f.chi(f.add(f.add(f.mul(x, x), f.mul(*c1, x)), *c0), 3) as i8).collect()).collect();
    let quota = |split: bool| -> usize { match (q, split) { (5, true) => 3, (7, true) => 4, (11, true) => 3, (13, true) => 2, _ => 2 } };
    let mut sel: Vec<(usize, u32, u64, u64, u64, bool)> = vec![]; let mut conj_first: Vec<(usize, u32, u64)> = vec![]; let mut nsel = [0usize; 2]; let mut tally: [HashMap<u64, u64>; 2] = [HashMap::new(), HashMap::new()]; let mut ncand = [0u64; 2];
    let mut cand_lines = vec![];
    for i in 0..q3 { let v = iso.cen2l[i as usize]; if f.in_sub(v, 1) { continue; }
        let (v1, v2) = (f.frob(v), f.frob(f.frob(v))); if iso.l2cen[&v1] < i || iso.l2cen[&v2] < i { continue; }
        let chip: Vec<i8> = kel.iter().map(|&x| f.chi(f.mul(f.sub(x, v), f.sub(x, v1)), 3) as i8).collect();
        for (k, qd) in qs.iter().enumerate() { let arm = if qd.0 { 0 } else { 1 };
            let s: i64 = chiq[k].iter().zip(chip.iter()).map(|(&a, &b)| (a * b) as i64).sum(); let n = (q3 as i64 + 2 + s) as u64;
            ncand[arm] += 1; *tally[arm].entry(n % 4).or_insert(0) += 1; if arm == 1 && conj_first.len() < 2 { conj_first.push((k, i, n)); }
            if nsel[arm] < quota(qd.0) { let c = cof(n); cand_lines.push(format!("{{\"rec\":\"candidate\",\"q\":{},\"arm\":\"{}\",\"p_census_idx\":{},\"Q\":{},\"order\":{},\"selected\":{}}}", q, if qd.0 { "split" } else { "conjugate" }, i, qd.4, n, c.is_some()));
                if let Some((hh, ll)) = c { sel.push((k, i, n, hh, ll, false)); nsel[arm] += 1; } } } }
    // amendment 2026-10-13 (PROTOCOL.md): an arm with no candidate meeting the cofactor rule takes its first 2 candidates, P1-P3 and P5 only
    if nsel[1] == 0 { for &(k, i, n) in &conj_first { let mut lp = 1; let mut t = n; let mut d = 2; while d * d <= t { while t % d == 0 { lp = d; t /= d; } d += 1; } if t > 1 { lp = t; } sel.push((k, i, n, n / lp, lp, true)); } }
    for l in &cand_lines { println!("{}", l); }
    for arm in 0..2 { let t = &tally[arm]; println!("{{\"rec\":\"census_mod4\",\"q\":{},\"arm\":\"{}\",\"candidates\":{},\"order_0mod4\":{},\"order_2mod4\":{},\"order_odd\":{}}}", q, if arm == 0 { "split" } else { "conjugate" }, ncand[arm],
        t.get(&0).unwrap_or(&0), t.get(&2).unwrap_or(&0), t.get(&1).unwrap_or(&0) + t.get(&3).unwrap_or(&0)); }
    let mut all_pass = true;
    for (inst_no, &(k, pidx, n, hh, ll, amended)) in sel.iter().enumerate() {
        let (split, c0, c1, ab, ref qs_s) = qs[k]; let p = iso.cen2l[pidx as usize];
        let mut rng = Rng::new((q as u64) << 32 | inst_no as u64);
        let mut out = format!("{{\"rec\":\"instance\",\"q\":{},\"arm\":\"{}\",\"Q\":{},\"Q_coeffs_c0_c1\":[{},{}],\"p_census\":{},\"r\":[{},{},{}],\"order_E0\":{},\"cofactor\":{},\"ell\":{},\"t_E\":{},\"selected_by_amendment\":{}",
            q, if split { "split" } else { "conjugate" }, qs_s, f.to_int(c0), f.to_int(c1), kc(&iso, q, p), kc(&iso, q, p), kc(&iso, q, f.frob(p)), kc(&iso, q, f.frob(f.frob(p))), n, hh, ll, q3 as i64 + 1 - n as i64, amended);
        let inst = match build(&f, split, c0, c1, ab, p) { Ok(x) => x, Err(e) => { println!("{},\"build_error\":\"{}\"}}", out, e); all_pass = false; continue; } };
        let mut pass: Vec<(&str, bool)> = vec![];
        // ---- P1 / P2
        let p1_x = p_in(&f, &inst.xn, 1) && p_in(&f, &inst.xd, 1); let p1_h = p_in(&f, &inst.h, 1);
        let degx = dg(&inst.xn).max(dg(&inst.xd)); let hd = dg(&inst.h); let hder = pder(&f, &inst.h);
        let sqf = !hder.is_empty() && dg(&pgcd(&f, &inst.h, &hder)) == 0;
        let gstr: Vec<String> = (0..9).map(|t| f.to_int(inst.g[t / 3][t % 3]).to_string()).collect();
        let v3 = |v: &[u32; 3]| format!("[{},{},{}]", f.to_int(v[0]), f.to_int(v[1]), f.to_int(v[2]));
        // random points of C over L
        let (mut tested, mut ok_curve, mut ok_x, mut ok_frob, mut ok_rt) = (0, 0, 0, 0, 0);
        while tested < 200 { let x = rng.below(f.m as u64) as u32 + 1; if inst.special(&f, x) { continue; }
            let ys: Vec<Option<u32>> = (0..3).map(|i| f.sqrt(pev(&f, &inst.fi[i], x))).collect(); if ys.iter().any(|y| y.is_none()) { continue; }
            let mut pt = [x, ys[0].unwrap(), ys[1].unwrap(), ys[2].unwrap()]; for i in 1..4 { if rng.below(2) == 1 { pt[i] = f.neg(pt[i]); } }
            let Some((z, w)) = inst.fwd(&f, pt) else { continue }; tested += 1;
            if f.mul(w, w) == pev(&f, &inst.h, z) { ok_curve += 1; }
            if pev(&f, &inst.xd, z) != 0 && f.div(pev(&f, &inst.xn, z), pev(&f, &inst.xd, z)) == x { ok_x += 1; }
            if inst.fwd(&f, frc(&f, pt)) == Some((f.frob(z), f.frob(w))) { ok_frob += 1; }
            if inst.bwd(&f, z, w) == Ok(pt) { ok_rt += 1; } }
        let p1 = p1_x && p1_h && ok_curve == 200 && ok_x == 200 && ok_frob == 200 && ok_rt == 200 && inst.pts_conic == q + 1 && (!split || inst.pts_over_a == 1);
        let p2 = degx == 4 && sqf && hd == if split { 7 } else { 8 };
        pass.push(("P1", p1)); pass.push(("P2", p2));
        out += &format!(",\"conic_G_Fq\":[{}],\"conic_Fq_points\":{},\"Fq_points_of_D_over_a\":{},\"P_star\":{},\"V0\":{},\"V1\":{},\"x_num\":{},\"x_den\":{},\"h\":{},\"deg_h\":{},\"deg_x\":{},\"h_squarefree\":{},\
\"P1_coeffs_in_Fq\":{},\"P1_random_points\":{},\"P1_on_curve\":{},\"P1_x_consistent\":{},\"P1_frobenius_commutes\":{},\"P1_roundtrip\":{},\"P1\":{},\"P2\":{}",
            gstr.join(","), inst.pts_conic, inst.pts_over_a, v3(&inst.ps), v3(&inst.v0), v3(&inst.v1), pints(&f, &inst.xn), pints(&f, &inst.xd), pints(&f, &inst.h), hd, degx, sqf, p1_x && p1_h, tested, ok_curve, ok_x, ok_frob, ok_rt, p1, p2);
        // ---- P3: point counts
        let mut nh = [0i64; 3]; let mut nk = [0i64; 3];
        for k in 1..=3u32 { let el = f.elems(k); let mut c: i64 = el.iter().map(|&z| 1 + f.chi(pev(&f, &inst.h, z), k) as i64).sum();
            c += if hd == 7 { 1 } else { 1 + f.chi(*inst.h.last().unwrap(), k) as i64 }; nh[k as usize - 1] = c; }
        nk[0] = f.elems(1).iter().map(|&x| 1 + f.chi(pev(&f, &inst.fi[0], x), 3) as i64).sum::<i64>() + 2;
        nk[1] = f.elems(2).iter().map(|&x| 1 + f.chi(pev(&f, &inst.fi[0], x), 6) as i64).sum::<i64>() + 2;
        nk[2] = 8 + kel.iter().map(|&x| { let r = &inst.r;
            if let Some(i) = r.iter().position(|&ri| ri == x) { (1 + f.chi(pev(&f, &inst.fi[(i + 1) % 3], x), 3) as i64) * (1 + f.chi(f.mul(f.sub(x, r[(i + 1) % 3]), f.sub(x, r[(i + 2) % 3])), 3) as i64) }
            else if pev(&f, &inst.qp, x) == 0 { (1 + f.chi(f.mul(f.sub(x, r[0]), f.sub(x, r[2])), 3) as i64) * (1 + f.chi(f.mul(f.sub(x, r[1]), f.sub(x, r[0])), 3) as i64) }
            else { (0..3).map(|i| 1 + f.chi(pev(&f, &inst.fi[i], x), 3) as i64).product::<i64>() } }).sum::<i64>();
        let qi = q as i64; let te = q3 as i64 + 1 - n as i64; let pred = [qi + 1, qi * qi + 1, qi * qi * qi + 1 - 3 * te];
        let pk: Vec<i64> = (0..3).map(|k| qi.pow(k as u32 + 1) + 1 - nh[k]).collect();
        let e1 = pk[0]; let e2x2 = e1 * pk[0] - pk[1]; let e2 = e2x2 / 2; let e3x3 = e2 * pk[0] - e1 * pk[1] + pk[2]; let e3 = e3x3 / 3;
        let (a1, a2, a3) = (-e1, e2, -e3); let lpoly_exact = e2x2 % 2 == 0 && e3x3 % 3 == 0;
        let l1 = 1 + a1 + a2 + a3 + qi * a2 + qi * qi * a1 + qi * qi * qi;
        let mut menum: i64 = -1;
        if split { menum = 0; for d in 0..=3u32 { for uc in 0..q.pow(d) { let mut u: P = digits(q, d as usize, uc).iter().map(|&c| f.int(c as i64)).collect(); u.push(1);
            let hm = pmod(&f, &inst.h, &u);
            for vc in 0..q.pow(d) { let v = tr(digits(q, d as usize, vc).iter().map(|&c| f.int(c as i64)).collect()); if pmod(&f, &pmul(&f, &v, &v), &u) == hm { menum += 1; } } } } }
        let p3 = nh.iter().zip(pred.iter()).all(|(a, b)| a == b) && nk == nh && lpoly_exact && a1 == 0 && a2 == 0 && a3 == -te && l1 == n as i64 && (!split || menum == n as i64);
        pass.push(("P3", p3));
        out += &format!(",\"N_hyperelliptic\":[{},{},{}],\"N_K_model\":[{},{},{}],\"N_predicted\":[{},{},{}],\"L_coeffs_a1_a2_a3\":[{},{},{}],\"L_of_1\":{},\"Jac_Fq_by_Mumford_enumeration\":{},\"P3\":{}",
            nh[0], nh[1], nh[2], nk[0], nk[1], nk[2], pred[0], pred[1], pred[2], a1, a2, a3, l1, if menum >= 0 { menum.to_string() } else { "null".into() }, p3);
        // ---- P5: census cross-check
        // Weierstrass model at r0: F0(r0 + t) = e1 t + e2 t^2 + e3 t^3 + t^4
        let mut tay = inst.fi[0].clone(); let mut ec = [0u32; 5]; for j in 0..5 { let (qq, rr) = pdivrem(&f, &tay, &plin(&f, inst.r[0])); ec[j] = *rr.first().unwrap_or(&0); tay = qq; }
        let (te1, te2, te3) = (ec[1], ec[2], ec[3]); let wc = Wc { a2: te2, a4: f.mul(te1, te3), a6: f.mul(te1, te1) };
        let three = f.int(3); let aa = f.sub(wc.a4, f.div(f.mul(wc.a2, wc.a2), three));
        let bb = f.add(f.sub(f.div(f.mul(f.int(2), f.mul(wc.a2, f.mul(wc.a2, wc.a2))), f.int(27)), f.div(f.mul(wc.a2, wc.a4), three)), wc.a6);
        let mut p5_lines = vec![]; let mut p5 = true; let census_scope = q == 5 || q == 7;
        for j in 0..3 { let (ja, jb) = { let mut x = (aa, bb); for _ in 0..j { x = (f.frob(x.0), f.frob(x.1)); } x };
            let hits = census_lookup(&cf, iso.l2cen[&ja], iso.l2cen[&jb], &recs);
            if census_scope && hits.is_empty() { p5 = false; }
            for hcr in &hits { if hcr.2 != n as i64 || hcr.3 != 3 || !hcr.4 { p5 = false; } }
            p5_lines.push(format!("{{\"conjugate\":{},\"A\":{},\"B\":{},\"records\":[{}]}}", j, kc(&iso, q, ja), kc(&iso, q, jb),
                hits.iter().map(|h| format!("{{\"role\":\"{}\",\"A_idx\":{},\"B_idx\":{},\"order\":{},\"min_genus\":{},\"has_m3_r5_hyperelliptic\":{}}}", h.5, h.0, h.1, h.2, h.3, h.4)).collect::<Vec<_>>().join(","))); }
        let wpts: Vec<WP> = { let mut v = vec![WP::O]; for &s in &kel { let rh = wrhs(&f, &wc, s); if rh == 0 { v.push(WP::A(s, 0)); } else if f.chi(rh, 3) == 1 { let y = f.sqrt(rh).unwrap(); v.push(WP::A(s, y)); v.push(WP::A(s, f.neg(y))); } } v };
        let w_order_ok = wpts.len() as u64 == n; if !w_order_ok { p5 = false; }
        pass.push(("P5a", p5));
        out += &format!(",\"weierstrass_at_r0\":{{\"a2\":{},\"a4\":{},\"a6\":{}}},\"short_weierstrass_census_basis\":{{\"A\":{},\"B\":{}}},\"weierstrass_order_matches\":{},\"census_scope_exhaustive\":{},\"census_lookup\":[{}],\"P5a\":{}",
            kc(&iso, q, wc.a2), kc(&iso, q, wc.a4), kc(&iso, q, wc.a6), kc(&iso, q, aa), kc(&iso, q, bb), w_order_ok, census_scope, p5_lines.join(","), p5);
        // ---- P4: transfer (split arm)
        if split {
            let (a, b) = ab.unwrap(); let r = inst.r; let h = &inst.h; let gg = inst.gen;
            let to_e = |pp: WP| -> Option<(u32, u32)> { match pp { WP::O => None, WP::A(s, y) => if s == 0 { None } else { Some((f.add(r[0], f.div(te1, s)), f.div(f.mul(y, te1), f.mul(s, s)))) } } };
            let to_w = |x: u32, y0: u32| -> WP { let d = f.sub(x, r[0]); WP::A(f.div(te1, d), f.div(f.mul(te1, y0), f.mul(d, d))) };
            let elig = |pp: WP| -> bool { match to_e(pp) { None => false, Some((x, _)) => !inst.special(&f, x) } };
            let phi = |pp: WP| -> Dv { let (x, y0) = to_e(pp).unwrap(); let s1 = f.sqrt(pev(&f, &inst.fi[1], x)).unwrap(); let s2 = f.sqrt(pev(&f, &inst.fi[2], x)).unwrap();
                let mut d = dzero(); for e in 0..4 { let mut pt = [x, y0, if e & 1 == 1 { f.neg(s1) } else { s1 }, if e & 2 == 2 { f.neg(s2) } else { s2 }];
                    for _ in 0..3 { let (z, w) = inst.fwd(&f, pt).expect("generic point maps"); d = dadd(&f, h, gg, &d, &dpt(&f, z, w)); pt = frc(&f, pt); } } d };
            let mut val: HashMap<WP, Dv> = HashMap::new(); let mut fq_ok = 0; let mut mumford_ok = 0; let el: Vec<WP> = wpts.iter().cloned().filter(|&pp| elig(pp)).collect();
            for &pp in &el { let d = phi(pp); if dv_in(&f, &d) { fq_ok += 1; } if dv_ok(&f, h, &d) { mumford_ok += 1; } val.insert(pp, d); }
            let exc: Vec<WP> = wpts.iter().cloned().filter(|&pp| !elig(pp)).collect();
            // direct values for exceptional points
            let ua = pmonic(&f, &psub(&f, &inst.xn, &psc(&f, &inst.xd, a))); let ub = pmonic(&f, &psub(&f, &inst.xn, &psc(&f, &inst.xd, b)));
            let ua_ok = dg(&ua) == 3 && pmod(&f, h, &ua).is_empty(); let ub_ok = dg(&ub) == 4 && pmod(&f, h, &ub).is_empty();
            let direct = |pp: WP| -> Option<Dv> { match to_e(pp) { None => if pp == WP::O { Some(dzero()) } else { None },
                Some((x, _)) => if x == a { Some(Dv { u: ua.clone(), v: vec![] }) } else if x == b { Some(dred(&f, h, gg, ub.clone(), vec![])) } else if x == r[1] { Some(dzero()) } else { None } } };
            let mut exc_lines = vec![]; let mut exc_ok = true;
            for &pp in &exc { let mut vals = vec![]; let mut tries = 0; while vals.len() < 3 && tries < 10000 { tries += 1; let rr = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, pp, rr);
                    if elig(s) { vals.push(dadd(&f, h, gg, &val[&s], &dneg(&f, &val[&rr]))); } }
                let agree = vals.len() == 3 && vals.iter().all(|v| v == &vals[0]); let dir = direct(pp); let dir_agree = dir.as_ref().map(|d| d == &vals[0]);
                if !agree || dir_agree == Some(false) { exc_ok = false; }
                let kind = match to_e(pp) { None => if pp == WP::O { "O=(r0,0)".to_string() } else { "x=inf".to_string() }, Some((x, _)) => if x == a { "(a,0)".into() } else if x == b { "(b,0)".into() } else if x == r[1] { "(r1,0)".into() } else if x == r[2] { "x=r2".into() } else { "other".into() } };
                exc_lines.push(format!("{{\"kind\":\"{}\",\"R_route_values_agree\":{},\"direct_value_agrees\":{},\"phi_is_zero\":{}}}", kind, agree, match dir_agree { Some(t) => t.to_string(), None => "null".into() }, vals[0] == dzero()));
                val.insert(pp, vals[0].clone()); }
            // homomorphism
            let (mut hom_t, mut hom_ok, mut dbl_t, mut dbl_ok) = (0, 0, 0, 0);
            while hom_t < 400 { let p1 = el[rng.below(el.len() as u64) as usize]; let p2 = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, p1, p2); if !elig(s) { continue; }
                hom_t += 1; if dadd(&f, h, gg, &val[&p1], &val[&p2]) == val[&s] { hom_ok += 1; } }
            let mut guard = 0; while dbl_t < 100 && guard < 100000 { guard += 1; let p1 = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, p1, p1); if !elig(s) { continue; }
                dbl_t += 1; if dadd(&f, h, gg, &val[&p1], &val[&p1]) == val[&s] { dbl_ok += 1; } }
            // generator of E0(K)[ell]
            let mut gpt = WP::O; while gpt == WP::O { gpt = wmul(&f, &wc, n / ll, wpts[rng.below(wpts.len() as u64) as usize]); }
            let gv = val[&gpt].clone(); let g_nonzero = gv != dzero(); let g_ell = dmul(&f, h, gg, ll, &gv) == dzero(); let g_elig = elig(gpt);
            // kernel and image over all of E0(K)
            let ker: Vec<WP> = wpts.iter().cloned().filter(|pp| val[pp] == dzero()).collect();
            let ker4 = ker.iter().all(|&pp| wmul(&f, &wc, 4, pp) == WP::O); let t1 = to_w(r[1], 0); let t1_in = ker.contains(&t1);
            let img: HashSet<Dv> = wpts.iter().map(|pp| val[pp].clone()).collect(); let img_ok = img.len() as u64 * ker.len() as u64 == n;
            let ker_orders: Vec<u64> = ker.iter().map(|&pp| (1..=4u64).find(|&k| wmul(&f, &wc, k, pp) == WP::O).unwrap()).collect();
            let all_vals_fq = wpts.iter().all(|pp| dv_in(&f, &val[pp]));
            // pi0_* phi = [4]
            let ta = to_w(a, 0); let (mut pf_ok, mut pf_bad, mut pf_skip) = (0, 0, 0); let (mut pfx_ok, mut pfx_bad) = (0, 0); let mut pf_why: HashMap<&str, u32> = HashMap::new();
            for &pp in &wpts { let is_el = elig(pp); let d = &val[&pp]; let mut rem = d.u.clone(); let mut roots: Vec<(u32, u64)> = vec![];
                for c in 0..q { let cv = f.int(c as i64); let mut mlt = 0; while dg(&rem) >= 1 && pev(&f, &rem, cv) == 0 { rem = pexact(&f, &rem, &plin(&f, cv)); mlt += 1; } if mlt > 0 { roots.push((cv, mlt)); } }
                if dg(&rem) == 2 { let (c0r, c1r) = (rem[0], rem[1]); let disc = f.sub(f.mul(c1r, c1r), f.mul(f.int(4), c0r)); let sd = f.sqrt(disc).unwrap(); let i2 = f.inv(f.int(2));
                    roots.push((f.mul(f.sub(sd, c1r), i2), 1)); roots.push((f.mul(f.sub(f.neg(sd), c1r), i2), 1)); }
                else if dg(&rem) == 3 { let rho = *kel.iter().find(|&&t| pev(&f, &rem, t) == 0).unwrap(); roots.push((rho, 1)); roots.push((f.frob(rho), 1)); roots.push((f.frob(f.frob(rho)), 1)); }
                let mut acc = WP::O; let mut skip = false; let mut acc_x = WP::O; let mut bad_x = false;
                for &(z, mlt) in &roots { let wv = pev(&f, &d.v, z); match inst.bwd(&f, z, wv) { Ok(pt) => { let t = wmul(&f, &wc, mlt, to_w(pt[0], pt[1])); acc = wadd(&f, &wc, acc, t); acc_x = wadd(&f, &wc, acc_x, t); }
                    Err(e) => { skip = true; if is_el { *pf_why.entry(e).or_insert(0) += 1; }
                        // extension beyond the protocol: images of points in the special fibres (derived formulas, see README)
                        let xdz = pev(&f, &inst.xd, z); let az = [pev(&f, &inst.xp[0], z), pev(&f, &inst.xp[1], z), pev(&f, &inst.xp[2], z)];
                        let img = if xdz == 0 { let prod = f.div(wv, pev(&f, &inst.xn, z)); let (r1v, r2v) = (f.div(az[0], az[1]), f.div(az[2], az[1]));
                                let v0 = f.div(prod, f.mul(r1v, r2v)); if f.mul(v0, v0) != 1 { None } else { Some(WP::A(0, f.mul(te1, v0))) } }
                            else { let x = f.div(pev(&f, &inst.xn, z), xdz);
                                if x == a || x == b || x == r[1] { Some(to_w(x, 0)) } else if x == r[0] { Some(WP::O) }
                                else if x == r[2] { let u12 = f.div(f.mul(f.sub(r[2], r[1]), az[2]), az[0]); let y0 = f.div(f.mul(wv, f.mul(f.sub(r[2], r[0]), f.sub(r[2], r[1]))), f.mul(u12, xdz));
                                    if f.mul(y0, y0) != pev(&f, &inst.fi[0], x) { None } else { Some(to_w(x, y0)) } } else { None } };
                        match img { Some(t) => { acc_x = wadd(&f, &wc, acc_x, wmul(&f, &wc, mlt, t)); } None => { bad_x = true; } } } } }
                let dpart = wmul(&f, &wc, dg(&d.u).max(0) as u64, ta); let four = wmul(&f, &wc, 4, pp);
                if bad_x || wadd(&f, &wc, acc_x, dpart) != four { pfx_bad += 1; } else { pfx_ok += 1; }
                if !is_el { continue; }
                if skip { pf_skip += 1; continue; }
                acc = wadd(&f, &wc, acc, dpart);
                if acc == four { pf_ok += 1; } else { pf_bad += 1; } }
            let p4 = fq_ok == el.len() && mumford_ok == el.len() && all_vals_fq && exc_ok && ua_ok && ub_ok && hom_ok == hom_t && dbl_ok == dbl_t && dbl_t == 100 && g_nonzero && g_ell && ker4 && t1_in && img_ok && pf_bad == 0 && pfx_bad == 0;
            pass.push(("P4", p4));
            let gpt_s = match gpt { WP::A(s, y) => format!("[{},{}]", kc(&iso, q, s), kc(&iso, q, y)), WP::O => "null".into() };
            out += &format!(",\"P4_detail\":{{\"E0_points\":{},\"eligible\":{},\"exceptional\":{},\"exceptional_detail\":[{}],\"phi_Fq_coeffs\":{},\"phi_valid_mumford\":{},\"u_a_deg3_divides_h\":{},\"u_b_deg4_divides_h\":{},\
\"hom_pairs\":{},\"hom_pass\":{},\"doublings\":{},\"doubling_pass\":{},\"G_weierstrass_SY\":{},\"G_eligible\":{},\"phi_G_nonzero\":{},\"ell_phi_G_zero\":{},\"kernel_size\":{},\"kernel_orders\":{:?},\"kernel_in_E0_4\":{},\"T1_in_kernel\":{},\"image_size\":{},\"image_times_kernel_eq_order\":{},\
\"pushforward_eq_4P\":{},\"pushforward_fail\":{},\"pushforward_not_evaluable\":{},\"pushforward_skip_reasons\":{:?},\"ext_pushforward_all_E0K_eq_4P\":{},\"ext_pushforward_fail\":{},\"phi_G\":{{\"u\":{},\"v\":{}}}}},\"P4\":{}",
                wpts.len(), el.len(), exc.len(), exc_lines.join(","), fq_ok, mumford_ok, ua_ok, ub_ok, hom_t, hom_ok, dbl_t, dbl_ok, gpt_s, g_elig, g_nonzero, g_ell, ker.len(), ker_orders, ker4, t1_in, img.len(), img_ok,
                pf_ok, pf_bad, pf_skip, { let mut v: Vec<_> = pf_why.iter().map(|(k, c)| format!("{}:{}", k, c)).collect(); v.sort(); v }, pfx_ok, pfx_bad, pints(&f, &gv.u), pints(&f, &gv.v), p4);
        } else {
            // ---- addendum 2026-10-14: transfer on the conjugate arm with real-model (balanced) arithmetic
            let r = inst.r; let h = &inst.h;
            // F_q-rational base point (z0, w0) of the model change
            let mut z0w0 = None; for c in 0..q { let z = f.int(c as i64); let hz = pev(&f, h, z); let xdz = pev(&f, &inst.xd, z);
                if hz != 0 && f.chi(hz, 1) == 1 && xdz != 0 && !inst.special(&f, f.div(pev(&f, &inst.xn, z), xdz)) { z0w0 = Some((z, f.sqrt(hz).unwrap())); break; } }
            let (z0, w0) = z0w0.expect("no base point");
            // f(z') = z'^8 h(z0 + 1/z') / h(z0)
            let mut hs = vec![0u32; 9]; { let mut t = h.clone(); for k in 0..9 { let (qq, rr) = pdivrem(&f, &t, &plin(&f, z0)); hs[k] = *rr.first().unwrap_or(&0); t = qq; } }
            let fr: P = tr((0..9).map(|k| f.div(hs[8 - k], hs[0])).collect());
            let mut vp = vec![0u32; 5]; vp[4] = 1; for k in (0..4).rev() { let sq = pmul(&f, &vp, &vp); vp[k] = f.div(f.sub(*fr.get(4 + k).unwrap_or(&0), *sq.get(4 + k).unwrap_or(&0)), f.int(2)); }
            let vp = tr(vp); let rm = Rm { f: fr.clone(), vp: vp.clone() };
            let fr_ok = dg(&fr) == 8 && fr[8] == 1 && p_in(&f, &fr, 1) && dg(&pgcd(&f, &fr, &pder(&f, &fr))) == 0 && dg(&psub(&f, &fr, &pmul(&f, &vp, &vp))) <= 3;
            let jp = Rd { u: vec![1], v: vec![], n: 3 }; // [inf+ - inf-]
            let ptc = |z: u32, w: u32| -> Rd { if z == z0 { if w == w0 { rzero() } else { Rd { u: vec![1], v: vec![], n: 1 } } }
                else { let d = f.sub(z, z0); let zp = f.inv(d); let wp = f.div(w, f.mul(f.pow(d, 4), w0)); Rd { u: plin(&f, zp), v: tr(vec![wp]), n: 1 } } };
            // V: enumeration of reduced triples, associativity, order, inverse
            let mut all: Vec<Rd> = vec![]; for d in 0..=3u32 { for uc in 0..q.pow(d) { let mut u: P = digits(q, d as usize, uc).iter().map(|&c| f.int(c as i64)).collect(); u.push(1);
                let fm = pmod(&f, &fr, &u);
                for vc in 0..q.pow(d) { let v = tr(digits(q, d as usize, vc).iter().map(|&c| f.int(c as i64)).collect()); if pmod(&f, &pmul(&f, &v, &v), &u) == fm { for nn in 0..=(3 - d as i32) { all.push(Rd { u: u.clone(), v: v.clone(), n: nn }); } } } } }
            let v1 = all.len() as u64 == n;
            let pick = |rng: &mut Rng| all[rng.below(all.len() as u64) as usize].clone();
            let (mut v2_ok, mut v3_ok, mut v4_ok) = (0, 0, 0);
            for _ in 0..300 { let (a1, b1, c1) = (pick(&mut rng), pick(&mut rng), pick(&mut rng)); if radd(&f, &rm, &radd(&f, &rm, &a1, &b1), &c1) == radd(&f, &rm, &a1, &radd(&f, &rm, &b1, &c1)) { v2_ok += 1; } }
            for _ in 0..100 { let d = pick(&mut rng); if rmul(&f, &rm, n, &d) == rzero() { v3_ok += 1; }
                if radd(&f, &rm, &d, &rneg(&f, &rm, &d)) == rzero() && radd(&f, &rm, &d, &rzero()) == d { v4_ok += 1; } }
            let vpass = fr_ok && v1 && v2_ok == 300 && v3_ok == 100 && v4_ok == 100;
            pass.push(("V", vpass));
            // E0 <-> W
            let to_e = |pp: WP| -> Option<(u32, u32)> { match pp { WP::O => None, WP::A(s, y) => if s == 0 { None } else { Some((f.add(r[0], f.div(te1, s)), f.div(f.mul(y, te1), f.mul(s, s)))) } } };
            let to_w = |x: u32, y0: u32| -> WP { let d = f.sub(x, r[0]); WP::A(f.div(te1, d), f.div(f.mul(te1, y0), f.mul(d, d))) };
            let elig = |pp: WP| -> bool { match to_e(pp) { None => false, Some((x, _)) => !inst.special(&f, x) } };
            let phi = |pp: WP| -> Rd { let (x, y0) = to_e(pp).unwrap(); let s1 = f.sqrt(pev(&f, &inst.fi[1], x)).unwrap(); let s2 = f.sqrt(pev(&f, &inst.fi[2], x)).unwrap();
                let mut d = rzero(); for e in 0..4 { let mut pt = [x, y0, if e & 1 == 1 { f.neg(s1) } else { s1 }, if e & 2 == 2 { f.neg(s2) } else { s2 }];
                    for _ in 0..3 { let c = match inst.fwd(&f, pt) { Some((z, w)) => ptc(z, w),
                            // point over z = infinity of the original model (fibre of the conic base point P*): z' = 0, w' = s/w0, s = lim w/z^4
                            None => { let s = inst.fwd_inf(&f, pt).expect("generic point maps"); Rd { u: vec![0, 1], v: tr(vec![f.div(s, w0)]), n: 1 } } };
                        d = radd(&f, &rm, &d, &c); pt = frc(&f, pt); } }
                radd(&f, &rm, &d, &rmul(&f, &rm, 6, &jp)) };
            let mut val: HashMap<WP, Rd> = HashMap::new(); let el: Vec<WP> = wpts.iter().cloned().filter(|&pp| elig(pp)).collect(); let mut valid_ok = 0;
            for &pp in &el { let d = phi(pp); if rvalid(&f, &rm, &d) { valid_ok += 1; } val.insert(pp, d); }
            let exc: Vec<WP> = wpts.iter().cloned().filter(|&pp| !elig(pp)).collect();
            let mut exc_lines = vec![]; let mut exc_ok = true;
            for &pp in &exc { let mut vals = vec![]; let mut tries = 0; while vals.len() < 3 && tries < 10000 { tries += 1; let rr = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, pp, rr);
                    if elig(s) { vals.push(radd(&f, &rm, &val[&s], &rneg(&f, &rm, &val[&rr]))); } }
                let agree = vals.len() == 3 && vals.iter().all(|v| v == &vals[0]);
                let dir = match to_e(pp) { None => if pp == WP::O { Some(rzero()) } else { None }, Some((x, _)) => if x == r[1] { Some(rzero()) } else { None } };
                let dir_agree = dir.as_ref().map(|d| d == &vals[0]); if !agree || dir_agree == Some(false) || !rvalid(&f, &rm, &vals[0]) { exc_ok = false; }
                let kind = match to_e(pp) { None => if pp == WP::O { "O=(r0,0)".to_string() } else { "x=inf".to_string() }, Some((x, _)) => if x == r[1] { "(r1,0)".into() } else if x == r[2] { "x=r2".into() } else { "other".into() } };
                exc_lines.push(format!("{{\"kind\":\"{}\",\"R_route_values_agree\":{},\"direct_value_agrees\":{},\"phi_is_zero\":{}}}", kind, agree, match dir_agree { Some(t) => t.to_string(), None => "null".into() }, vals[0] == rzero()));
                val.insert(pp, vals[0].clone()); }
            let (mut hom_t, mut hom_ok, mut dbl_t, mut dbl_ok) = (0, 0, 0, 0);
            while hom_t < 400 { let p1 = el[rng.below(el.len() as u64) as usize]; let p2 = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, p1, p2); if !elig(s) { continue; }
                hom_t += 1; if radd(&f, &rm, &val[&p1], &val[&p2]) == val[&s] { hom_ok += 1; } }
            let mut guard = 0; while dbl_t < 100 && guard < 100000 { guard += 1; let p1 = el[rng.below(el.len() as u64) as usize]; let s = wadd(&f, &wc, p1, p1); if !elig(s) { continue; }
                dbl_t += 1; if radd(&f, &rm, &val[&p1], &val[&p1]) == val[&s] { dbl_ok += 1; } }
            let mut gpt = WP::O; while gpt == WP::O { gpt = wmul(&f, &wc, n / ll, wpts[rng.below(wpts.len() as u64) as usize]); }
            let gv = val[&gpt].clone(); let g_nonzero = gv != rzero(); let g_ell = rmul(&f, &rm, ll, &gv) == rzero();
            let ker: Vec<WP> = wpts.iter().cloned().filter(|pp| val[pp] == rzero()).collect();
            let ker4 = ker.iter().all(|&pp| wmul(&f, &wc, 4, pp) == WP::O); let t1_in = ker.contains(&to_w(r[1], 0));
            let img: HashSet<Rd> = wpts.iter().map(|pp| val[pp].clone()).collect(); let img_ok = img.len() as u64 * ker.len() as u64 == n;
            let ker_orders: Vec<u64> = ker.iter().map(|&pp| (1..=n).find(|&k| wmul(&f, &wc, k, pp) == WP::O).unwrap()).collect();
            let all_valid = wpts.iter().all(|pp| rvalid(&f, &rm, &val[pp]));
            // pi0_* : images of points of C given homogeneous data (X at z, xn, xd, w), special fibres included
            let pimg = |az: [u32; 3], xnum: u32, xden: u32, wv: u32| -> Option<WP> {
                if xden == 0 { let prod = f.div(wv, xnum); let (r1v, r2v) = (f.div(az[0], az[1]), f.div(az[2], az[1])); let v0 = f.div(prod, f.mul(r1v, r2v));
                    return if f.mul(v0, v0) != 1 { None } else { Some(WP::A(0, f.mul(te1, v0))) }; }
                let x = f.div(xnum, xden);
                if pev(&f, &inst.qp, x) == 0 || x == r[1] { return Some(to_w(x, 0)); } if x == r[0] { return Some(WP::O); }
                if x == r[2] { let u12 = f.div(f.mul(f.sub(r[2], r[1]), az[2]), az[0]); let y0 = f.div(f.mul(wv, f.mul(f.sub(r[2], r[0]), f.sub(r[2], r[1]))), f.mul(u12, xden));
                    return if f.mul(y0, y0) != pev(&f, &inst.fi[0], x) { None } else { Some(to_w(x, y0)) }; }
                let mut s = [0u32; 3]; for i in 0..3 { s[i] = f.sqrt(pev(&f, &inst.fi[i], x))?; }
                let target = f.div(f.mul(f.mul(wv, pev(&f, &inst.qp, x)), pev(&f, &inst.np, x)), xden); let mut found = vec![];
                for e in 0..8u32 { let y: Vec<u32> = (0..3).map(|i| if e >> i & 1 == 1 { f.neg(s[i]) } else { s[i] }).collect();
                    if f.mul(f.mul(y[0], y[1]), y[2]) != target { continue; } let pt = [x, y[0], y[1], y[2]]; if parallel(&f, &inst.xtrue(&f, pt), &az) { found.push(pt); } }
                if found.len() == 1 { Some(to_w(found[0][0], found[0][1])) } else { None } };
            let img_fin = |z: u32, w: u32| -> Option<WP> { pimg([pev(&f, &inst.xp[0], z), pev(&f, &inst.xp[1], z), pev(&f, &inst.xp[2], z)], pev(&f, &inst.xn, z), pev(&f, &inst.xd, z), w) };
            let top = |p: &P, k: usize| *p.get(k).unwrap_or(&0);
            let img_inf = |s: u32| -> Option<WP> { pimg([top(&inst.xp[0], 2), top(&inst.xp[1], 2), top(&inst.xp[2], 2)], top(&inst.xn, 4), top(&inst.xd, 4), s) };
            let (ip_plus, ip_minus) = (img_fin(z0, w0), img_fin(z0, f.neg(w0)));
            let (mut pf_ok, mut pf_bad, mut pf_ne) = (0, 0, 0);
            for &pp in &wpts { let d = &val[&pp]; let mut acc = WP::O; let mut bad = ip_plus.is_none() || ip_minus.is_none();
                for (zp, mlt) in roots3(&f, &d.u, &kel) { let wp = pev(&f, &d.v, zp);
                    let im = if zp == 0 { img_inf(f.mul(wp, w0)) } else { let zi = f.inv(zp); img_fin(f.add(z0, zi), f.div(f.mul(wp, w0), f.pow(zp, 4))) };
                    match im { Some(t) => { acc = wadd(&f, &wc, acc, wmul(&f, &wc, mlt, t)); } None => { bad = true; } } }
                if bad { pf_ne += 1; continue; }
                let (nn, mm) = (d.n as i64 - 2, 3 - dg(&d.u) as i64 - d.n as i64 - 1);
                let sc = |k: i64, t: WP| -> WP { if k >= 0 { wmul(&f, &wc, k as u64, t) } else { wneg(&f, wmul(&f, &wc, (-k) as u64, t)) } };
                acc = wadd(&f, &wc, acc, wadd(&f, &wc, sc(nn, ip_plus.unwrap()), sc(mm, ip_minus.unwrap())));
                if acc == wmul(&f, &wc, 4, pp) { pf_ok += 1; } else { pf_bad += 1; } }
            let p4c = vpass && valid_ok == el.len() && all_valid && exc_ok && hom_ok == hom_t && dbl_ok == dbl_t && dbl_t == 100 && g_nonzero && g_ell && ker4 && t1_in && img_ok && pf_bad == 0 && pf_ne == 0;
            pass.push(("P4c", p4c));
            let gpt_s = match gpt { WP::A(s, y) => format!("[{},{}]", kc(&iso, q, s), kc(&iso, q, y)), WP::O => "null".into() };
            out += &format!(",\"real_model\":{{\"z0\":{},\"w0\":{},\"f\":{},\"V_plus\":{},\"f_monic_deg8_squarefree\":{},\"reduced_triples_enumerated\":{},\"V1_enum_eq_order\":{},\"V2_assoc_pass\":{},\"V3_order_kills_pass\":{},\"V4_inverse_identity_pass\":{},\"V\":{}}},\
\"P4c_detail\":{{\"E0_points\":{},\"eligible\":{},\"exceptional\":{},\"exceptional_detail\":[{}],\"phi_valid_reduced_Fq\":{},\"all_values_valid\":{},\"hom_pairs\":{},\"hom_pass\":{},\"doublings\":{},\"doubling_pass\":{},\"G_weierstrass_SY\":{},\"G_eligible\":{},\"phi_G_nonzero\":{},\"ell_phi_G_zero\":{},\
\"kernel_size\":{},\"kernel_orders\":{:?},\"kernel_in_E0_4\":{},\"T1_in_kernel\":{},\"image_size\":{},\"image_times_kernel_eq_order\":{},\"pushforward_all_E0K_eq_4P\":{},\"pushforward_fail\":{},\"pushforward_not_evaluable\":{},\"phi_G\":{{\"u\":{},\"v\":{},\"n\":{}}}}},\"V\":{},\"P4c\":{}",
                f.to_int(z0), f.to_int(w0), pints(&f, &fr), pints(&f, &vp), fr_ok, all.len(), v1, v2_ok, v3_ok, v4_ok, vpass,
                wpts.len(), el.len(), exc.len(), exc_lines.join(","), valid_ok, all_valid, hom_t, hom_ok, dbl_t, dbl_ok, gpt_s, elig(gpt), g_nonzero, g_ell,
                ker.len(), ker_orders, ker4, t1_in, img.len(), img_ok, pf_ok, pf_bad, pf_ne, pints(&f, &gv.u), pints(&f, &gv.v), gv.n, vpass, p4c);
        }
        let ip = pass.iter().all(|x| x.1); all_pass &= ip;
        out += &format!(",\"all_pass\":{}}}", ip); println!("{}", out);
    }
    let p5b = tally.iter().all(|t| t.get(&2).is_none() && t.get(&1).is_none() && t.get(&3).is_none());
    println!("{{\"rec\":\"totals\",\"q\":{},\"instances\":{},\"all_instances_pass\":{},\"P5b_all_candidates_0mod4\":{},\"wall_s\":{:.1}}}", q, sel.len(), all_pass, p5b, t0.elapsed().as_secs_f64());
}
