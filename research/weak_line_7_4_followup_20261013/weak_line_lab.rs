// Follow-up lab: why the (7,4) weak curve of weak_curve_lab.rs E3 has two weak lines (PROTOCOL.md).
// Everything up to "follow-up code (new)" is a verbatim fork of weak_curve_lab.rs (field, polynomials, Cantor, census).
// rustc --edition 2021 -O weak_line_lab.rs ; ./weak_line_lab <f1|f2|f3|f4|f5>
#![allow(dead_code)]
use std::collections::HashSet;

const SEED: u64 = 20261012;
const CAP: usize = 12_000_000;

struct Rng(u64);
impl Rng {
    fn new(tag: u64) -> Self { Rng(SEED ^ tag.wrapping_mul(0x9E3779B97F4A7C15)) }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u32) -> u32 { (self.next() % n as u64) as u32 }
}

// ---------------------------------------------------------------- field GF(q^n)
// Element encoding: 0 is zero, e >= 1 is alpha^(e-1). "Index" encoding: base-q digits.
struct Fld { q: u32, n: u32, size: u32, m: u32, d: u32, zech: Vec<u32>, idx_of: Vec<u32>, el_of: Vec<u32>, modulus: Vec<u32>, gen_shift: u32 }
impl Fld {
    fn build(q: u32, n: u32) -> Fld {
        let size = q.pow(n); let m = size - 1; let nn = n as usize;
        let digits = |mut i: u32| -> Vec<u32> { let mut v = vec![0; nn]; for k in 0..nn { v[k] = i % q; i /= q; } v };
        let index = |v: &[u32]| -> u32 { let mut i = 0; for k in (0..nn).rev() { i = i * q + v[k]; } i };
        for cand in 0..size { let md = digits(cand);
            for s in 0..q {
                if n == 1 && (s == 0 || cand != 0) { continue; }
                let mut exp_idx = vec![0u32; m as usize]; let mut cur = vec![0u32; nn]; cur[0] = 1; let mut ok = true;
                for k in 0..m as usize {
                    let ci = index(&cur);
                    if ci == 0 || (k > 0 && ci == 1) { ok = false; break; }
                    exp_idx[k] = ci;
                    let mut nx = vec![0u32; nn];
                    if n == 1 { nx[0] = cur[0] * s % q; } else {
                        let top = cur[nn - 1];
                        for i in 0..nn { let sh = if i > 0 { cur[i - 1] } else { 0 }; nx[i] = (sh + (q - md[i]) * top + s * cur[i]) % q; }
                    }
                    cur = nx;
                }
                if !ok || index(&cur) != 1 { continue; }
                let mut el_of = vec![0u32; size as usize]; let mut idx_of = vec![0u32; size as usize];
                for k in 0..m as usize { el_of[exp_idx[k] as usize] = k as u32 + 1; idx_of[k + 1] = exp_idx[k]; }
                let mut zech = vec![0u32; m as usize];
                for k in 0..m as usize { let mut dv = digits(exp_idx[k]); dv[0] = (dv[0] + 1) % q; zech[k] = el_of[index(&dv) as usize]; }
                return Fld { q, n, size, m, d: m / (q - 1), zech, idx_of, el_of, modulus: md, gen_shift: s };
            }
        }
        panic!("no field");
    }
    #[inline] fn mul(&self, a: u32, b: u32) -> u32 { if a == 0 || b == 0 { 0 } else { (a - 1 + b - 1) % self.m + 1 } }
    #[inline] fn inv(&self, a: u32) -> u32 { assert!(a != 0); (self.m - (a - 1)) % self.m + 1 }
    #[inline] fn neg(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1 + self.m / 2) % self.m + 1 } }
    #[inline] fn add(&self, a: u32, b: u32) -> u32 {
        if a == 0 { return b; } if b == 0 { return a; }
        let z = self.zech[((b + self.m - a) % self.m) as usize];
        if z == 0 { 0 } else { (a - 1 + z - 1) % self.m + 1 }
    }
    #[inline] fn sub(&self, a: u32, b: u32) -> u32 { self.add(a, self.neg(b)) }
    fn pow(&self, a: u32, e: u64) -> u32 { if a == 0 { return if e == 0 { 1 } else { 0 }; } (((a - 1) as u64 * (e % self.m as u64)) % self.m as u64) as u32 + 1 }
    #[inline] fn is_sq(&self, a: u32) -> bool { a == 0 || (a - 1) % 2 == 0 }
    #[inline] fn sqrt(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1) / 2 + 1 } }
    #[inline] fn in_fq(&self, a: u32) -> bool { a == 0 || (a - 1) % self.d == 0 }
    fn in_sub(&self, a: u32, k: u32) -> bool { a == 0 || (a - 1) % (self.m / (self.q.pow(k) - 1)) == 0 }
    fn int(&self, c: u32) -> u32 { self.el_of[(c % self.q) as usize] }
    fn to_int(&self, a: u32) -> u32 { assert!(self.in_fq(a)); self.idx_of[a as usize] }
    fn rand(&self, r: &mut Rng) -> u32 { r.below(self.size) }
    fn rand_nz(&self, r: &mut Rng) -> u32 { 1 + r.below(self.m) }
    fn two_inv(&self) -> u32 { self.inv(self.int(2)) }
}

// ---------------------------------------------------------------- polynomials (low -> high, trimmed)
type Poly = Vec<u32>;
fn trim(mut a: Poly) -> Poly { while a.last() == Some(&0) { a.pop(); } a }
fn padd(f: &Fld, a: &Poly, b: &Poly) -> Poly { let n = a.len().max(b.len()); let mut r = vec![0; n]; for i in 0..n { r[i] = f.add(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)); } trim(r) }
fn pneg(f: &Fld, a: &Poly) -> Poly { a.iter().map(|&x| f.neg(x)).collect() }
fn psub(f: &Fld, a: &Poly, b: &Poly) -> Poly { padd(f, a, &pneg(f, b)) }
fn pmul(f: &Fld, a: &Poly, b: &Poly) -> Poly { if a.is_empty() || b.is_empty() { return vec![]; } let mut r = vec![0; a.len() + b.len() - 1]; for i in 0..a.len() { if a[i] == 0 { continue; } for j in 0..b.len() { r[i + j] = f.add(r[i + j], f.mul(a[i], b[j])); } } trim(r) }
fn pscale(f: &Fld, a: &Poly, c: u32) -> Poly { trim(a.iter().map(|&x| f.mul(x, c)).collect()) }
fn pdivrem(f: &Fld, a: &Poly, b: &Poly) -> (Poly, Poly) {
    assert!(!b.is_empty()); let mut r = a.clone(); if r.len() < b.len() { return (vec![], r); }
    let mut qv = vec![0; r.len() - b.len() + 1]; let li = f.inv(*b.last().unwrap());
    while r.len() >= b.len() { let k = r.len() - b.len(); let c = f.mul(*r.last().unwrap(), li); qv[k] = c;
        for j in 0..b.len() { r[k + j] = f.sub(r[k + j], f.mul(c, b[j])); } r = trim(r); if r.is_empty() { break; } }
    (trim(qv), r)
}
fn pmod(f: &Fld, a: &Poly, b: &Poly) -> Poly { pdivrem(f, a, b).1 }
fn pxgcd(f: &Fld, a: &Poly, b: &Poly) -> (Poly, Poly, Poly) {
    let (mut r0, mut s0, mut t0) = (a.clone(), vec![1u32], vec![]); let (mut r1, mut s1, mut t1) = (b.clone(), vec![], vec![1u32]);
    while !r1.is_empty() { let (qq, r) = pdivrem(f, &r0, &r1);
        let s = psub(f, &s0, &pmul(f, &qq, &s1)); let t = psub(f, &t0, &pmul(f, &qq, &t1));
        r0 = r1; s0 = s1; t0 = t1; r1 = r; s1 = s; t1 = t; }
    let c = f.inv(*r0.last().unwrap()); (pscale(f, &r0, c), pscale(f, &s0, c), pscale(f, &t0, c))
}
fn peval(f: &Fld, a: &Poly, x: u32) -> u32 { let mut r = 0; for &c in a.iter().rev() { r = f.add(f.mul(r, x), c); } r }
fn pderiv(f: &Fld, a: &Poly) -> Poly { let mut r = vec![]; for i in 1..a.len() { r.push(f.mul(a[i], f.int(i as u32 % f.q))); } trim(r) }
fn squarefree(f: &Fld, a: &Poly) -> bool { let d = pderiv(f, a); !d.is_empty() && pxgcd(f, a, &d).0.len() == 1 }
// f(beta*t + b)
fn paffine(f: &Fld, a: &Poly, beta: u32, b: u32) -> Poly { let lin = vec![b, beta]; let mut r: Poly = vec![]; for &c in a.iter().rev() { r = pmul(f, &r, &lin); r = padd(f, &r, &trim(vec![c])); } r }

// ---------------------------------------------------------------- Jacobian, Mumford + Cantor
#[derive(Clone, PartialEq)]
struct Div { u: Poly, v: Poly }
fn dzero() -> Div { Div { u: vec![1], v: vec![] } }
fn key(d: &Div) -> u128 { let g = |p: &Poly, i: usize| *p.get(i).unwrap_or(&0) as u128; let deg = d.u.len() as u128 - 1;
    let (u0, u1) = if deg == 2 { (g(&d.u, 0), g(&d.u, 1)) } else if deg == 1 { (g(&d.u, 0), 0) } else { (0, 0) };
    deg | u0 << 2 | u1 << 20 | g(&d.v, 0) << 38 | g(&d.v, 1) << 56 }
fn unkey(k: u128) -> Div { let m = (1u128 << 18) - 1; let deg = (k & 3) as usize; let u0 = ((k >> 2) & m) as u32; let u1 = ((k >> 20) & m) as u32;
    let u = match deg { 0 => vec![1], 1 => vec![u0, 1], _ => vec![u0, u1, 1] }; Div { u, v: trim(vec![((k >> 38) & m) as u32, ((k >> 56) & m) as u32]) } }
fn dadd(f: &Fld, cf: &Poly, a: &Div, b: &Div) -> Div {
    let (d1, e1, e2) = pxgcd(f, &a.u, &b.u);
    let (d, c1, c2) = if d1.len() == 1 { (vec![1u32], vec![1u32], vec![]) } else { pxgcd(f, &d1, &padd(f, &a.v, &b.v)) };
    let (s1, s2, s3) = (pmul(f, &c1, &e1), pmul(f, &c1, &e2), c2);
    let mut u = pdivrem(f, &pmul(f, &a.u, &b.u), &pmul(f, &d, &d)).0;
    let num = padd(f, &padd(f, &pmul(f, &pmul(f, &s1, &a.u), &b.v), &pmul(f, &pmul(f, &s2, &b.u), &a.v)), &pmul(f, &s3, &padd(f, &pmul(f, &a.v, &b.v), cf)));
    let mut v = pmod(f, &pdivrem(f, &num, &d).0, &u);
    while u.len() > 3 { let (u2, rem) = pdivrem(f, &psub(f, cf, &pmul(f, &v, &v)), &u); debug_assert!(rem.is_empty());
        u = pscale(f, &u2, f.inv(*u2.last().unwrap())); v = pmod(f, &pneg(f, &v), &u); }
    Div { u, v }
}
fn dmul(f: &Fld, cf: &Poly, a: &Div, mut k: u128) -> Div { let mut r = dzero(); let mut b = a.clone(); while k > 0 { if k & 1 == 1 { r = dadd(f, cf, &r, &b); } b = dadd(f, cf, &b, &b); k >>= 1; } r }
fn dvalid(f: &Fld, cf: &Poly, a: &Div) -> bool { pmod(f, &psub(f, &pmul(f, &a.v, &a.v), cf), &a.u).is_empty() }
fn point_div(f: &Fld, cf: &Poly, x: u32) -> Option<Div> { let fx = peval(f, cf, x); if !f.is_sq(fx) { return None; } Some(Div { u: vec![f.neg(x), 1], v: trim(vec![f.sqrt(fx)]) }) }
fn rand_div(f: &Fld, cf: &Poly, r: &mut Rng) -> Div { let mut acc = dzero(); let mut got = 0; while got < 2 { if let Some(p) = point_div(f, cf, f.rand(r)) { acc = dadd(f, cf, &acc, &p); got += 1; } } acc }

// Exact subgroup generated by `gens`, by coset extension. Returns (keys, lower_bound, capped).
fn gen_subgroup(f: &Fld, cf: &Poly, gens: &[Div], full: Option<u128>) -> (Vec<u128>, HashSet<u128>, u128, bool) {
    let mut elems = vec![key(&dzero())]; let mut set: HashSet<u128> = elems.iter().cloned().collect();
    for g in gens {
        if set.contains(&key(g)) { continue; }
        let mut mult = vec![key(g)]; let mut cur = g.clone();
        loop { cur = dadd(f, cf, &cur, g); let k = key(&cur); if set.contains(&k) { break; } mult.push(k); if mult.len() > CAP { break; } }
        let newsize = elems.len() as u128 * (mult.len() as u128 + 1);
        if newsize > CAP as u128 { return (elems, set, newsize, true); }
        let old = elems.clone();
        for mk in &mult { let md = unkey(*mk); for h in &old { let e = key(&dadd(f, cf, &unkey(*h), &md)); if set.insert(e) { elems.push(e); } } }
        if Some(elems.len() as u128) == full { break; }
    }
    let n = elems.len() as u128; (elems, set, n, false)
}

// ---------------------------------------------------------------- point counts and group orders
fn chi(f: &Fld, a: u32) -> i64 { if a == 0 { 0 } else if f.is_sq(a) { 1 } else { -1 } }
fn count_k(f: &Fld, cf: &Poly) -> i128 { let mut c = 1i128; for x in 0..f.size { c += 1 + chi(f, peval(f, cf, x)) as i128; } c }
fn count_k2(f: &Fld, cf: &Poly) -> i128 { // over K[t]/(t^2 - alpha)
    let ns = 2u32; let mut c = 1i128;
    for xa in 0..f.size { for xb in 0..f.size { let (mut ra, mut rb) = (0u32, 0u32);
        for &co in cf.iter().rev() { let na = f.add(f.add(f.mul(ra, xa), f.mul(ns, f.mul(rb, xb))), co); let nb = f.add(f.mul(ra, xb), f.mul(rb, xa)); ra = na; rb = nb; }
        if ra == 0 && rb == 0 { c += 1; } else { let nm = f.sub(f.mul(ra, ra), f.mul(ns, f.mul(rb, rb))); c += if f.is_sq(nm) { 2 } else { 0 }; } } }
    c
}
fn jac_order_counts(f: &Fld, cf: &Poly) -> u128 { let n1 = count_k(f, cf); let n2 = count_k2(f, cf); ((n1 * n1 + n2) / 2 - f.size as i128) as u128 }
type M4 = [[i128; 4]; 4];
fn m4mul(a: &M4, b: &M4) -> M4 { let mut r = [[0i128; 4]; 4]; for i in 0..4 { for j in 0..4 { for k in 0..4 { r[i][j] += a[i][k] * b[k][j]; } } } r }
fn det(mut a: Vec<Vec<i128>>) -> i128 { // Bareiss
    let n = a.len(); let mut sign = 1i128; let mut prev = 1i128;
    for k in 0..n - 1 { if a[k][k] == 0 { let mut sw = None; for i in k + 1..n { if a[i][k] != 0 { sw = Some(i); break; } } match sw { None => return 0, Some(i) => { a.swap(i, k); sign = -sign; } } }
        for i in k + 1..n { for j in k + 1..n { a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev; } } prev = a[k][k]; }
    sign * a[n - 1][n - 1]
}
// Frobenius data of a genus-2 curve y^2 = psi(x) over F_q (psi given by integer coefficients, degree 5, any leading coeff).
struct Frob { a1: i128, a2: i128, q: i128 }
impl Frob {
    fn of(q: u32, psi: &[u32]) -> Frob {
        let f1 = Fld::build(q, 1); let f2 = Fld::build(q, 2);
        let p1: Poly = psi.iter().map(|&c| f1.int(c)).collect(); let p2: Poly = psi.iter().map(|&c| f2.int(c)).collect();
        // odd degree model: one point at infinity
        let n1 = count_k(&f1, &p1); let n2 = count_k(&f2, &p2); let qq = q as i128;
        let a1 = n1 - qq - 1; let a2 = (n2 - qq * qq - 1 + a1 * a1) / 2; Frob { a1, a2, q: qq }
    }
    fn jac(&self, k: u32) -> u128 { // det(1 - M^k), M companion of T^4 + a1 T^3 + a2 T^2 + q a1 T + q^2
        let c = [self.q * self.q, self.q * self.a1, self.a2, self.a1];
        let mut m: M4 = [[0; 4]; 4]; for i in 0..3 { m[i + 1][i] = 1; } for i in 0..4 { m[i][3] = -c[i]; }
        let mut p: M4 = [[0; 4]; 4]; for i in 0..4 { p[i][i] = 1; } for _ in 0..k { p = m4mul(&p, &m); }
        let v: Vec<Vec<i128>> = (0..4).map(|i| (0..4).map(|j| (if i == j { 1 } else { 0 }) - p[i][j]).collect()).collect(); det(v) as u128
    }
}
fn factor(mut n: u128) -> Vec<(u128, u32)> { let mut r = vec![]; let mut p = 2u128; while p * p <= n { let mut e = 0; while n % p == 0 { n /= p; e += 1; } if e > 0 { r.push((p, e)); } p += if p == 2 { 1 } else { 2 }; } if n > 1 { r.push((n, 1)); } r }
fn maxprime(n: u128) -> u128 { factor(n).last().map(|x| x.0).unwrap_or(1) }
fn fstr(n: u128) -> String { factor(n).iter().map(|(p, e)| if *e > 1 { format!("{}^{}", p, e) } else { format!("{}", p) }).collect::<Vec<_>>().join("*") }
fn gcd(a: u128, b: u128) -> u128 { if b == 0 { a } else { gcd(b, a % b) } }
fn bound(q: u32) -> f64 { (1.0 + (q as f64).sqrt()).powi(4) }

// ---------------------------------------------------------------- curves and lines
fn rand_monic_quintic(f: &Fld, r: &mut Rng, over_fq: bool) -> Poly { loop { let mut p = vec![0u32; 6]; p[5] = 1;
    for i in 0..5 { p[i] = if over_fq { f.int(1 + r.below(f.q - 1)) } else { f.rand_nz(r) }; } if squarefree(f, &p) { return p; } } }
// A line: x = beta*t + b, t in F_q.  Weak data: c with c*f(beta t + b) in F_q[t].
struct Weak { rational: bool, square_ok: bool, c: u32, psi: Vec<u32> }
fn weak_data(f: &Fld, cf: &Poly, beta: u32, b: u32) -> Weak {
    let g = paffine(f, cf, beta, b); let c0 = f.inv(*g.last().unwrap());
    if !g.iter().all(|&x| f.in_fq(f.mul(x, c0))) { return Weak { rational: false, square_ok: false, c: 0, psi: vec![] }; }
    let mut c = c0; let mut ok = f.is_sq(c0);
    if !ok && f.n % 2 == 1 { for e in 2..f.q { let ee = f.int(e); if !f.is_sq(ee) { c = f.mul(c0, ee); ok = true; break; } } }
    Weak { rational: true, square_ok: ok, c, psi: g.iter().map(|&x| f.to_int(f.mul(x, c))).collect() }
}
fn line_xs(f: &Fld, beta: u32, b: u32) -> Vec<u32> { (0..f.q).map(|t| f.add(f.mul(beta, f.int(t)), b)).collect() }
fn fb_from_xs(f: &Fld, cf: &Poly, xs: &[u32]) -> Vec<Div> { xs.iter().filter_map(|&x| point_div(f, cf, x)).collect() }
// paper condition (16): exists a with j_i + i a == 5 a (mod d) for i = 0..4, j_i = log_alpha f_i
fn paper_condition(f: &Fld, cf: &Poly) -> Option<u32> { if cf.iter().any(|&c| c == 0) { return None; }
    'a: for a in 0..f.d { for i in 0..5u64 { let lhs = ((cf[i as usize] - 1) as u64 + i * a as u64) % f.d as u64; if lhs != (5 * a as u64) % f.d as u64 { continue 'a; } } return Some(a); } None }
// all affine weak lines by brute force over (beta class, b); returns unique (j, b_canon, square_ok)
fn census(f: &Fld, cf: &Poly, shortcut: bool) -> Vec<(u32, u32, bool)> {
    let mut out: HashSet<(u32, u32, bool)> = HashSet::new();
    let inv5 = if f.q != 5 { f.inv(f.int(5)) } else { 0 };
    for j in 0..f.d { let beta = j + 1;
        let bs: Vec<u32> = if shortcut { vec![f.neg(f.mul(cf[4], inv5))] } else { (0..f.size).collect() };
        for b in bs { let mut arr = cf.clone(); let mut ok = true;
            for k in 0..6usize { for i in (k..5).rev() { arr[i] = f.add(arr[i], f.mul(b, arr[i + 1])); }
                // coefficient k of f(beta t + b) scaled by beta^-5 is arr[k]*beta^(k-5)
                let e = f.mul(arr[k], f.pow(f.inv(beta), 5 - k as u64)); if !f.in_fq(e) { ok = false; break; } }
            if ok { let w = weak_data(f, cf, beta, b); assert!(w.rational);
                let canon = (0..f.q).map(|s| f.idx_of[f.add(b, f.mul(beta, f.int(s))) as usize]).min().unwrap(); out.insert((j, canon, w.square_ok)); } } }
    let mut v: Vec<_> = out.into_iter().collect(); v.sort(); v
}

// relation candidates (c, lambda), lambda of exact degree m, lambda^2 - c*S in F_q[t]; S of degree 5 (m=2) or 6 (m=3)
fn candidates(f: &Fld, s: &Poly) -> u64 {
    let deg = s.len() - 1; let ti = f.two_inv(); let fq: Vec<u32> = (0..f.q).map(|e| f.int(e)).collect(); let mut count = 0u64;
    let cs: Vec<u32> = if deg == 5 { let c0 = f.inv(s[5]); fq[1..].iter().map(|&e| f.mul(c0, e)).filter(|&c| f.is_sq(c)).collect() } else { (0..f.m / 2).map(|k| 2 * k + 1).collect() };
    for &c in &cs { let cs_: Vec<u32> = s.iter().map(|&x| f.mul(c, x)).collect(); let m = deg / 2;
        for &e in &fq { let top2 = f.add(cs_[2 * m], e); if top2 == 0 || !f.is_sq(top2) { continue; } let r = f.sqrt(top2);
            for lm in [r, f.neg(r)] { let h = f.mul(ti, f.inv(lm));
                if m == 2 { for &e1 in &fq { let l1 = f.mul(f.add(e1, cs_[3]), h);
                    for &e0 in &fq { let l0 = f.mul(f.sub(f.add(e0, cs_[2]), f.mul(l1, l1)), h);
                        let h1 = f.sub(f.mul(f.int(2), f.mul(l1, l0)), cs_[1]); let h0 = f.sub(f.mul(l0, l0), cs_[0]);
                        if f.in_fq(h1) && f.in_fq(h0) { count += 1; } } }
                } else { for &e2 in &fq { let l2 = f.mul(f.add(e2, cs_[5]), h);
                    for &e1 in &fq { let l1 = f.mul(f.sub(f.add(e1, cs_[4]), f.mul(l2, l2)), h);
                        for &e0 in &fq { let l0 = f.mul(f.sub(f.add(e0, cs_[3]), f.mul(f.int(2), f.mul(l2, l1))), h);
                            let h2 = f.sub(f.add(f.mul(l1, l1), f.mul(f.int(2), f.mul(l2, l0))), cs_[2]);
                            let h1 = f.sub(f.mul(f.int(2), f.mul(l1, l0)), cs_[1]); let h0 = f.sub(f.mul(l0, l0), cs_[0]);
                            if f.in_fq(h2) && f.in_fq(h1) && f.in_fq(h0) { count += 1; } } } } }
            } } }
    count
}

fn selftest(f: &Fld, cf: &Poly, n: u128, r: &mut Rng) -> bool { for _ in 0..20 { let a = rand_div(f, cf, r); let b = rand_div(f, cf, r); let c = rand_div(f, cf, r);
    if !dvalid(f, cf, &a) || dmul(f, cf, &a, n) != dzero() { return false; }
    if dadd(f, cf, &dadd(f, cf, &a, &b), &c) != dadd(f, cf, &a, &dadd(f, cf, &b, &c)) { return false; } } true }
fn polystr(p: &[u32]) -> String { format!("[{}]", p.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")) }
fn idxpoly(f: &Fld, p: &Poly) -> String { polystr(&p.iter().map(|&x| f.idx_of[x as usize]).collect::<Vec<_>>()) }
fn field_json(f: &Fld) -> String { format!("\"q\":{},\"n\":{},\"modulus_low\":{},\"gen\":\"x+{}\"", f.q, f.n, polystr(&f.modulus), f.gen_shift) }

// ================================================================ follow-up code (new)
// Affine maps x -> z x + c over K, as (z, c); composition (a o b)(x) = a(b(x)).
type Aff = (u32, u32);
fn acomp(f: &Fld, a: Aff, b: Aff) -> Aff { (f.mul(a.0, b.0), f.add(f.mul(a.0, b.1), a.1)) }
fn ainv(f: &Fld, a: Aff) -> Aff { let zi = f.inv(a.0); (zi, f.neg(f.mul(zi, a.1))) }
fn afrob(f: &Fld, a: Aff) -> Aff { (f.pow(a.0, f.q as u64), f.pow(a.1, f.q as u64)) }
fn aid() -> Aff { (1, 0) }
// N(nu) = phi^{n-1}(nu) o ... o phi(nu) o nu
fn anorm(f: &Fld, nu: Aff) -> Aff { let mut acc = nu; let mut cur = nu; for _ in 1..f.n { cur = afrob(f, cur); acc = acomp(f, cur, acc); } acc }
fn aorder(f: &Fld, a: Aff) -> u32 { let mut c = a; let mut k = 1; while c != aid() && k < 100 { c = acomp(f, a, c); k += 1; } k }
fn defk(f: &Fld, a: Aff) -> u32 { for k in 1..=f.n { if f.n % k == 0 && f.in_sub(a.0, k) && f.in_sub(a.1, k) { return k; } } f.n }
fn inv5(f: &Fld) -> u32 { f.inv(f.int(5)) }
// all (z, c) with p(z x + c) = z^5 qq(x), p, qq monic quintics; c forced by the x^4 coefficient. zs: candidate scalars.
fn affine_equiv(f: &Fld, p: &Poly, qq: &Poly, zs: &[u32]) -> Vec<Aff> { assert!(f.q != 5); let i5 = inv5(f); let mut out = vec![];
    for &z in zs { let c = f.mul(f.sub(f.mul(z, qq[4]), p[4]), i5);
        let lhs = paffine(f, p, z, c); let rhs = pscale(f, qq, f.pow(z, 5)); if lhs == rhs { out.push((z, c)); } }
    out }
fn all_scalars(f: &Fld) -> Vec<u32> { (1..=f.m).collect() }
// scalars that can occur in a non-trivial Aut of a quintic: roots of unity of order 2, 4, 5 (plus 1)
fn root_scalars(f: &Fld) -> Vec<u32> { let mut v: Vec<u32> = vec![1]; for r in [2u32, 4, 5] { if f.m % r == 0 { for j in 1..r { let z = (f.m / r) * j + 1; if !v.contains(&z) { v.push(z); } } } } v }
fn aut(f: &Fld, g: &Poly, full: bool) -> Vec<Aff> { affine_equiv(f, g, g, &if full { all_scalars(f) } else { root_scalars(f) }) }
fn lemma_a(f: &Fld, auts: &[Aff]) -> usize { auts.iter().filter(|&&nu| anorm(f, nu) == aid()).count() }
fn make_weak(f: &Fld, g: &Poly, beta: u32) -> Poly { (0..6).map(|k| f.mul(g[k], f.pow(beta, 5 - k as u64))).collect() }
fn rand_beta(f: &Fld, r: &mut Rng) -> u32 { loop { let b = f.rand_nz(r); if !f.in_fq(b) { return b; } } }
fn gints(f: &Fld, g: &Poly) -> String { polystr(&g.iter().map(|&x| f.to_int(x)).collect::<Vec<_>>()) }
fn affstr(f: &Fld, a: Aff) -> String { format!("[{},{}]", f.idx_of[a.0 as usize], f.idx_of[a.1 as usize]) }
fn autjson(f: &Fld, auts: &[Aff]) -> String { format!("[{}]", auts.iter().map(|&a| format!("{{\"nu_idx\":{},\"order\":{},\"def_k\":{},\"norm_is_id\":{}}}", affstr(f, a), aorder(f, a), defk(f, a), anorm(f, a) == aid())).collect::<Vec<_>>().join(",")) }
// monic g(s + x0) centred: coefficients of g(x0 + s), x0 = -g4/5
fn centred(f: &Fld, g: &Poly) -> Poly { let x0 = f.neg(f.mul(g[4], inv5(f))); paffine(f, g, 1, x0) }
fn is_odd(f: &Fld, g: &Poly) -> bool { let c = centred(f, g); c.len() == 6 && (0..6).step_by(2).all(|i| c[i] == 0) }

// ---------------------------------------------------------------- F1: reproduce the E3 weak-curve record at (7,4)
fn e3_weak_curve(q: u32, n: u32) -> (Fld, Poly, Poly, u32) { let f = Fld::build(q, n); let mut r = Rng::new(300 + (q * 10 + n) as u64);
    let g = rand_monic_quintic(&f, &mut r, true); let beta = loop { let b = f.rand_nz(&mut r); if !f.in_fq(b) { break b; } };
    let weak = make_weak(&f, &g, beta); (f, g, weak, beta) }
fn f1() { let (f, _g, cf, _b) = e3_weak_curve(7, 4); let lines = census(&f, &cf, false);
    println!("{{\"exp\":\"E3\",{},\"curve\":\"{}\",\"f_idx\":{},\"affine_lines_total\":{},\"weak_lines\":{},\"weak_lines_square_ok\":{},\"weak_scaling_lines\":{}}}",
        field_json(&f), "weak", idxpoly(&f, &cf), f.d as u64 * (f.size / f.q) as u64, lines.len(), lines.iter().filter(|l| l.2).count(), lines.iter().filter(|l| l.1 == 0).count()); }

// ---------------------------------------------------------------- F2: the two lines
fn f2() { let (f, g, cf, beta) = e3_weak_curve(7, 4); let q = f.q;
    let lines = census(&f, &cf, false); assert!(lines.len() == 2);
    let ncount = jac_order_counts(&f, &cf);
    let gc = centred(&f, &g);
    println!("{{\"exp\":\"F2\",\"part\":\"curve\",{},\"g\":{},\"g_centred_ints\":{},\"g_odd_about_Fq_point\":{},\"x0_idx\":{},\"beta_idx\":{},\"f_idx\":{},\"J_K_counts\":{},\"lines\":{}}}",
        field_json(&f), gints(&f, &g), gints(&f, &gc), is_odd(&f, &g), f.idx_of[f.neg(f.mul(g[4], inv5(&f))) as usize], f.idx_of[beta as usize], idxpoly(&f, &cf), ncount, lines.len());
    let mut data: Vec<(u32, u32, Weak)> = vec![];
    for (i, ln) in lines.iter().enumerate() { let (be, b) = (ln.0 + 1, f.el_of[ln.1 as usize]); let w = weak_data(&f, &cf, be, b);
        let monic: Poly = w.psi.iter().map(|&x| f.int(x)).collect(); let fr = Frob::of(q, &w.psi);
        let xs = line_xs(&f, be, b); let scaling = xs.contains(&0);
        let (j1, j2, jn) = (fr.jac(1), fr.jac(2), fr.jac(f.n));
        // the monic psi on this line is psi(t) = f(be t + b)/be^5 up to the F_q scalar in w.c; check
        let direct = pscale(&f, &paffine(&f, &cf, be, b), f.inv(f.pow(be, 5)));
        println!("{{\"exp\":\"F2\",\"part\":\"line\",\"i\":{},\"beta_idx\":{},\"b_idx\":{},\"contains_0_scaling\":{},\"c_idx\":{},\"square_ok\":{},\"psi_monic_equals_f_over_beta5\":{},\"psi\":{},\"a1\":{},\"a2\":{},\"J_Fq\":{},\"J_Fq2\":{},\"J_Fq4_lpoly\":{},\"J_K_equals_counts\":{},\"J_Fq_factored\":\"{}\",\"J_Fq2_factored\":\"{}\"}}",
            i, f.idx_of[be as usize], f.idx_of[b as usize], scaling, f.idx_of[w.c as usize], w.square_ok, direct == monic, polystr(&w.psi), fr.a1, fr.a2, j1, j2, jn, jn == ncount, fstr(j1), fstr(j2));
        data.push((be, b, w)); }
    let (p1, p2): (Poly, Poly) = (data[0].2.psi.iter().map(|&x| f.int(x)).collect(), data[1].2.psi.iter().map(|&x| f.int(x)).collect());
    // K-affine map mu = sigma1^{-1} sigma2 : t -> (be2 t + b2 - b1)/be1, and nu = mu^phi mu^{-1}
    let (be1, b1, be2, b2) = (data[0].0, data[0].1, data[1].0, data[1].1);
    let mu: Aff = (f.mul(be2, f.inv(be1)), f.mul(f.sub(b2, b1), f.inv(be1)));
    let nu = acomp(&f, afrob(&f, mu), ainv(&f, mu)); let aut1 = aut(&f, &p1, true);
    println!("{{\"exp\":\"F2\",\"part\":\"mu\",\"mu_idx\":{},\"mu_scalar_in_Fq2\":{},\"mu_scalar_in_Fq\":{},\"mu_c_in_Fq2\":{},\"nu_idx\":{},\"nu_order\":{},\"nu_in_Aut_psi1\":{},\"nu_def_k\":{},\"psi1_o_mu_equals_scal_psi2\":{}}}",
        affstr(&f, mu), f.in_sub(mu.0, 2), f.in_fq(mu.0), f.in_sub(mu.1, 2), affstr(&f, nu), aorder(&f, nu), aut1.contains(&nu), defk(&f, nu),
        paffine(&f, &p1, mu.0, mu.1) == pscale(&f, &p2, f.pow(mu.0, 5)));
    // equivalence psi2 ~ psi1 over F_q, F_{q^2}, K ; psi1 ~ g
    let eq_all = affine_equiv(&f, &p1, &p2, &all_scalars(&f));
    let eq_k = |k: u32| eq_all.iter().filter(|a| f.in_sub(a.0, k) && f.in_sub(a.1, k)).count();
    let (gaut_full, gaut_fast) = (aut(&f, &g, true), aut(&f, &g, false));
    let g_vs_p1 = affine_equiv(&f, &g, &p1, &all_scalars(&f)); let g_vs_p2 = affine_equiv(&f, &g, &p2, &all_scalars(&f));
    println!("{{\"exp\":\"F2\",\"part\":\"equiv\",\"psi1_to_psi2_maps_over_K\":{},\"over_Fq2\":{},\"over_Fq\":{},\"psi1_to_psi2_maps\":{},\"g_to_psi1_over_Fq\":{},\"g_to_psi2_over_Fq\":{},\"g_to_psi1_over_K\":{},\"g_to_psi2_over_K\":{}}}",
        eq_all.len(), eq_k(2), eq_k(1), format!("[{}]", eq_all.iter().map(|&a| format!("{{\"map_idx\":{},\"def_k\":{}}}", affstr(&f, a), defk(&f, a))).collect::<Vec<_>>().join(",")),
        g_vs_p1.iter().filter(|a| f.in_fq(a.0) && f.in_fq(a.1)).count(), g_vs_p2.iter().filter(|a| f.in_fq(a.0) && f.in_fq(a.1)).count(), g_vs_p1.len(), g_vs_p2.len());
    println!("{{\"exp\":\"F2\",\"part\":\"aut_g\",\"aut_full\":{},\"aut_fast_agrees\":{},\"size\":{},\"lemma_A_count\":{},\"census_lines\":{}}}",
        autjson(&f, &gaut_full), { let mut a = gaut_full.clone(); a.sort(); let mut b = gaut_fast.clone(); b.sort(); a == b }, gaut_full.len(), lemma_a(&f, &gaut_full), lines.len());
    let aut2 = aut(&f, &p2, true);
    println!("{{\"exp\":\"F2\",\"part\":\"aut_psi\",\"aut_psi1_size\":{},\"aut_psi2_size\":{},\"aut_psi1\":{},\"aut_psi2\":{}}}", aut1.len(), aut2.len(), autjson(&f, &aut1), autjson(&f, &aut2));
}

// ---------------------------------------------------------------- F3: control census
fn f3() { let fields: [(u32, u32, usize); 13] = [(3, 2, 300), (7, 2, 300), (11, 2, 300), (13, 2, 300), (3, 4, 300), (7, 4, 2000), (11, 4, 300), (13, 4, 300), (3, 6, 300), (3, 3, 300), (7, 3, 300), (11, 3, 300), (3, 5, 300)];
    for (q, n, ncurves) in fields { let f = Fld::build(q, n); let mut r = Rng::new(7000 + (q * 10 + n) as u64);
        let (mut mism, mut xcheck, mut xcheck_bad, mut autx_bad) = (0usize, 0usize, 0usize, 0usize);
        let mut hist: std::collections::BTreeMap<(usize, usize), usize> = Default::default(); let mut odd_cnt = 0usize; let mut odd_even_extra = 0usize;
        for i in 0..ncurves {
            let g = rand_monic_quintic(&f, &mut r, true); let beta = rand_beta(&f, &mut r); let cf = make_weak(&f, &g, beta);
            let small = f.size <= 2401;
            let census_n = if small { census(&f, &cf, false).len() } else { census(&f, &cf, true).len() };
            let full = f.size <= 30000 && (small || i < 5);
            let auts = aut(&f, &g, full);
            if full { let mut a = auts.clone(); a.sort(); let mut b = aut(&f, &g, false); b.sort(); if a != b { autx_bad += 1; } }
            let la = lemma_a(&f, &auts);
            if !small && i < 5 { xcheck += 1; if census(&f, &cf, false).len() != census_n { xcheck_bad += 1; } }
            if small { xcheck += 1; if census(&f, &cf, true).len() != census_n { xcheck_bad += 1; } }
            if la != census_n { mism += 1; }
            *hist.entry((census_n, auts.len())).or_insert(0) += 1;
            if is_odd(&f, &g) { odd_cnt += 1; if census_n > 1 { odd_even_extra += 1; } }
            if census_n > 1 || auts.len() > 1 { println!("{{\"exp\":\"F3\",\"part\":\"curve\",\"q\":{},\"n\":{},\"i\":{},\"g\":{},\"census_lines\":{},\"lemmaA\":{},\"aut_size\":{},\"odd_about_Fq_point\":{},\"aut\":{}}}", q, n, i, gints(&f, &g), census_n, la, auts.len(), is_odd(&f, &g), autjson(&f, &auts)); }
        }
        println!("{{\"exp\":\"F3\",\"part\":\"field\",\"q\":{},\"n\":{},\"curves\":{},\"mismatch_lemmaA_vs_census\":{},\"shortcut_vs_exhaustive_checked\":{},\"shortcut_vs_exhaustive_bad\":{},\"aut_fast_vs_full_bad\":{},\"odd_g\":{},\"odd_g_with_extra_line\":{},\"hist_census_lines_autsize\":{:?}}}",
            q, n, ncurves, mism, xcheck, xcheck_bad, autx_bad, odd_cnt, odd_even_extra, hist.iter().map(|(k, v)| format!("({},{}):{}", k.0, k.1, v)).collect::<Vec<_>>());
    } }

// ---------------------------------------------------------------- F4: population by g (Lemma A with root-of-unity scalars)
fn f4() { for (q, n) in [(3u32, 2u32), (3, 3), (3, 4), (3, 6), (7, 2), (7, 3), (7, 4), (11, 2), (11, 3), (11, 4), (13, 2), (13, 3), (13, 4)] { let f = Fld::build(q, n);
        let mut total = 0u64; let mut lines_hist: std::collections::BTreeMap<usize, u64> = Default::default(); let mut kinds: std::collections::BTreeMap<String, u64> = Default::default();
        let mut idx = vec![1u32; 5]; let nz = q - 1;
        loop { let mut g: Poly = idx.iter().map(|&e| f.int(e)).collect(); g.push(1);
            if squarefree(&f, &g) { total += 1; let auts = aut(&f, &g, false); let la = lemma_a(&f, &auts); *lines_hist.entry(la).or_insert(0) += 1;
                if auts.len() > 1 { let mut key: Vec<String> = auts.iter().filter(|&&a| a != aid()).map(|&a| format!("o{}k{}{}", aorder(&f, a), defk(&f, a), if anorm(&f, a) == aid() { "N1" } else { "" })).collect(); key.sort();
                    *kinds.entry(format!("{}|{}", if is_odd(&f, &g) { "odd" } else { "-" }, key.join(","))).or_insert(0) += 1; } }
            let mut k = 0; while k < 5 { idx[k] += 1; if idx[k] <= nz { break; } idx[k] = 1; k += 1; } if k == 5 { break; } }
        println!("{{\"exp\":\"F4\",\"q\":{},\"n\":{},\"g_total\":{},\"hist_lines\":{:?},\"nontrivial_aut_kinds\":{:?}}}", q, n, total,
            lines_hist.iter().map(|(k, v)| format!("{}:{}", k, v)).collect::<Vec<_>>(), kinds.iter().map(|(k, v)| format!("{}:{}", k, v)).collect::<Vec<_>>());
    } }

// ---------------------------------------------------------------- F5: positive controls (constructed symmetric g)
fn f5() {
    for (q, n, arm) in [(3u32, 2u32, "odd"), (7, 2, "odd"), (3, 4, "odd"), (7, 4, "odd"), (3, 6, "odd"), (3, 3, "odd"), (7, 3, "odd"), (3, 5, "odd"), (11, 2, "x5+c"), (11, 3, "x5+c"), (11, 4, "x5+c"), (11, 5, "x5+c-supp")] {
        let f = Fld::build(q, n); let mut r = Rng::new(8000 + (q * 10 + n) as u64); let (mut ok, mut tot) = (0, 0); let mut hist: std::collections::BTreeMap<(usize, usize), u32> = Default::default();
        for _ in 0..20 {
            let g: Poly = if arm == "odd" { loop { let x0 = f.int(r.below(q)); let (a, b) = (f.int(1 + r.below(q - 1)), f.int(1 + r.below(q - 1)));
                let s = vec![0, b, 0, a, 0, 1]; let g = paffine(&f, &s, 1, f.neg(x0)); if squarefree(&f, &g) { break g; } } }
              else { let c = f.int(1 + r.below(q - 1)); vec![c, 0, 0, 0, 0, 1] };
            let cf = make_weak(&f, &g, rand_beta(&f, &mut r)); let small = f.size <= 2401;
            let cn = census(&f, &cf, !small).len(); let auts = aut(&f, &g, small); let la = lemma_a(&f, &auts);
            tot += 1; if la == cn { ok += 1; } *hist.entry((cn, auts.len())).or_insert(0) += 1; }
        println!("{{\"exp\":\"F5\",\"q\":{},\"n\":{},\"arm\":\"{}\",\"curves\":{},\"lemmaA_equals_census\":{},\"hist_census_lines_autsize\":{:?},\"census_kind\":\"{}\"}}", q, n, arm, tot, ok, hist.iter().map(|(k, v)| format!("({},{}):{}", k.0, k.1, v)).collect::<Vec<_>>(), if f.size <= 2401 { "exhaustive" } else { "shortcut" });
    } }

fn main() { let arg = std::env::args().nth(1).unwrap_or("all".into());
    let run = |n: &str| arg == "all" || arg == n;
    if run("f1") { f1(); } if run("f2") { f2(); } if run("f3") { f3(); } if run("f4") { f4(); } if run("f5") { f5(); } }
