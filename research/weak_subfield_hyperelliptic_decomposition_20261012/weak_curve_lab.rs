// Toy-field laboratory for the Bhati-Singh "weak curve" decomposition (PROTOCOL.md).
// Genus 2, C/K : y^2 = f(x), f monic quintic, K = GF(q^n), q an odd prime.
// Dependency-free: rustc -O weak_curve_lab.rs ; ./weak_curve_lab <e1|e2|e3|e3b|e4|e4x|e5|all>
// Every result is one JSON object per line on stdout. e4x (a, b >= 1) was added after E4 was seen; it is supplementary.
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

// ---------------------------------------------------------------- E1: collapse + control
fn e1() { for (q, n) in [(7u32, 3u32), (11, 3), (13, 3), (5, 4), (7, 4)] { let f = Fld::build(q, n); let mut r = Rng::new(100 + (q * 10 + n) as u64);
    let sub = if n % 2 == 0 { 2 } else { 1 };
    for i in 0..3 { // weak curves
        let g = rand_monic_quintic(&f, &mut r, true);
        let beta = loop { let b = f.rand_nz(&mut r); if f.in_fq(b) { continue; } if n % 2 == 0 && !f.is_sq(b) { continue; } break b; };
        let cf: Poly = (0..6).map(|k| f.mul(g[k], f.pow(beta, 5 - k as u64))).collect();
        let w = weak_data(&f, &cf, beta, 0); assert!(w.rational && w.square_ok);
        let fr = Frob::of(q, &w.psi); let (j1, j2, jn) = (fr.jac(1), fr.jac(2), fr.jac(n));
        let ncount = jac_order_counts(&f, &cf);
        let fb = fb_from_xs(&f, &cf, &line_xs(&f, beta, 0));
        let (elems, _, hsz, capped) = gen_subgroup(&f, &cf, &fb, Some(ncount));
        // coordinates of every element of H in the C0 model: u(beta t) monic and mu*v(beta t)
        let mu = f.sqrt(w.c); let mut coords_ok = true;
        for k in &elems { let dv = unkey(*k); let ut = paffine(&f, &dv.u, beta, 0); let ut = pscale(&f, &ut, f.inv(*ut.last().unwrap())); let vt = pscale(&f, &paffine(&f, &dv.v, beta, 0), mu);
            if !ut.iter().chain(vt.iter()).all(|&x| f.in_sub(x, sub)) { coords_ok = false; } }
        let jsub = if sub == 2 { j2 } else { j1 };
        println!("{{\"exp\":\"E1\",\"kind\":\"weak\",{},\"i\":{},\"g\":{},\"beta_idx\":{},\"f_idx\":{},\"paper_a\":{},\"psi\":{},\"J_K_counts\":{},\"J_K_lpoly\":{},\"orders_agree\":{},\"J_C0_Fq\":{},\"J_C0_Fq2\":{},\"fb\":{},\"H\":{},\"capped\":{},\"H_divides_sub\":{},\"coords_in_subfield\":{},\"index\":{},\"maxprime_J\":{},\"maxprime_H\":{},\"bound\":{:.1},\"selftest\":{}}}",
            field_json(&f), i, polystr(&g.iter().map(|&x| f.to_int(x)).collect::<Vec<_>>()), f.idx_of[beta as usize], idxpoly(&f, &cf),
            paper_condition(&f, &cf).map(|a| a.to_string()).unwrap_or("null".into()), polystr(&w.psi), ncount, jn, ncount == jn, j1, j2, fb.len(), hsz, capped,
            jsub % hsz == 0, coords_ok, ncount / hsz, maxprime(ncount), maxprime(hsz), bound(q), selftest(&f, &cf, ncount, &mut r));
    }
    for i in 0..3 { // random controls
        let cf = loop { let c = rand_monic_quintic(&f, &mut r, false); if !weak_data(&f, &c, 1, 0).rational { break c; } };
        let ncount = jac_order_counts(&f, &cf); let fb = fb_from_xs(&f, &cf, &line_xs(&f, 1, 0));
        let (_, _, hsz, capped) = gen_subgroup(&f, &cf, &fb, Some(ncount));
        println!("{{\"exp\":\"E1\",\"kind\":\"random\",{},\"i\":{},\"f_idx\":{},\"paper_a\":{},\"J_K_counts\":{},\"fb\":{},\"H\":{},\"capped\":{},\"index\":{},\"H_divides_J\":{},\"maxprime_J\":{},\"maxprime_H\":{},\"bound\":{:.1},\"selftest\":{}}}",
            field_json(&f), i, idxpoly(&f, &cf), paper_condition(&f, &cf).map(|a| a.to_string()).unwrap_or("null".into()), ncount, fb.len(), hsz, capped,
            if capped { 0 } else { ncount / hsz }, ncount % hsz == 0, maxprime(ncount), maxprime(hsz), bound(q), selftest(&f, &cf, ncount, &mut r));
    } } }

// ---------------------------------------------------------------- E2: factor-base lines on one weak curve
fn e2() { for (q, n) in [(7u32, 3u32), (11, 3), (13, 3), (7, 4)] { let f = Fld::build(q, n);
    let mut r = Rng::new(100 + (q * 10 + n) as u64); // same stream as E1: first weak curve
    let g = rand_monic_quintic(&f, &mut r, true);
    let beta = loop { let b = f.rand_nz(&mut r); if f.in_fq(b) { continue; } if n % 2 == 0 && !f.is_sq(b) { continue; } break b; };
    let cf: Poly = (0..6).map(|k| f.mul(g[k], f.pow(beta, 5 - k as u64))).collect();
    let ncount = jac_order_counts(&f, &cf); let mut r = Rng::new(200 + (q * 10 + n) as u64);
    let mut lines: Vec<(&str, u32, u32)> = vec![("weak", beta, 0)];
    for _ in 0..3 { let b2 = loop { let x = f.rand_nz(&mut r); if !f.in_fq(f.mul(x, f.inv(beta))) { break x; } }; lines.push(("scaling", b2, 0)); }
    for _ in 0..3 { let b = loop { let x = f.rand_nz(&mut r); if !f.in_fq(f.mul(x, f.inv(beta))) { break x; } }; lines.push(("shift", beta, b)); }
    for _ in 0..3 { lines.push(("affine", f.rand_nz(&mut r), f.rand(&mut r))); }
    for (kind, be, b) in lines { let w = weak_data(&f, &cf, be, b); let model = paffine(&f, &cf, be, b);
        let fb = fb_from_xs(&f, &cf, &line_xs(&f, be, b)); let (_, _, hsz, capped) = gen_subgroup(&f, &cf, &fb, Some(ncount));
        println!("{{\"exp\":\"E2\",{},\"line\":\"{}\",\"model_deg\":5,\"beta_idx\":{},\"b_idx\":{},\"weak_rational\":{},\"candidates\":{},\"fb\":{},\"H\":{},\"capped\":{},\"J_K\":{},\"index\":{},\"maxprime_H\":{},\"maxprime_J\":{},\"bound\":{:.1}}}",
            field_json(&f), kind, f.idx_of[be as usize], f.idx_of[b as usize], w.rational, candidates(&f, &model), fb.len(), hsz, capped, ncount, if capped { 0 } else { ncount / hsz }, maxprime(hsz), maxprime(ncount), bound(q));
    }
    for _ in 0..3 { // Moebius lines x = (a t + b)/(c t + d), t in P^1(F_q)
        let (a, b, c, d) = loop { let (a, b, c, d) = (f.rand(&mut r), f.rand(&mut r), f.rand_nz(&mut r), f.rand(&mut r)); if f.sub(f.mul(a, d), f.mul(b, c)) != 0 && peval(&f, &cf, f.mul(a, f.inv(c))) != 0 { break (a, b, c, d); } };
        let mut xs = vec![f.mul(a, f.inv(c))];
        for t in 0..q { let tt = f.int(t); let den = f.add(f.mul(c, tt), d); if den != 0 { xs.push(f.mul(f.add(f.mul(a, tt), b), f.inv(den))); } }
        let (num, den) = (vec![b, a], trim(vec![d, c])); let mut s: Poly = vec![]; let mut dp: Vec<Poly> = vec![vec![1]]; for _ in 0..6 { let l = dp.last().unwrap().clone(); dp.push(pmul(&f, &l, &den)); }
        let mut np: Poly = vec![1]; for i in 0..6 { s = padd(&f, &s, &pscale(&f, &pmul(&f, &np, &dp[6 - i]), cf[i])); np = pmul(&f, &np, &num); }
        let fb = fb_from_xs(&f, &cf, &xs); let (_, _, hsz, capped) = gen_subgroup(&f, &cf, &fb, Some(ncount));
        println!("{{\"exp\":\"E2\",{},\"line\":\"moebius\",\"model_deg\":{},\"abcd_idx\":{},\"weak_rational\":false,\"candidates\":{},\"fb\":{},\"H\":{},\"capped\":{},\"J_K\":{},\"index\":{},\"maxprime_H\":{},\"maxprime_J\":{},\"bound\":{:.1}}}",
            field_json(&f), s.len() - 1, polystr(&[a, b, c, d].iter().map(|&x| f.idx_of[x as usize]).collect::<Vec<_>>()), candidates(&f, &s), fb.len(), hsz, capped, ncount, if capped { 0 } else { ncount / hsz }, maxprime(hsz), maxprime(ncount), bound(q));
    } } }

// ---------------------------------------------------------------- E3: census of weak affine lines
fn special(_f: &Fld, kind: &str, c: u32) -> Poly { match kind { "x^5+c" => vec![c, 0, 0, 0, 0, 1], _ => vec![0, c, 0, 0, 0, 1] } }
fn e3() { for (q, n) in [(7u32, 3u32), (11, 3), (13, 3), (5, 4), (7, 4), (3, 4)] { let f = Fld::build(q, n); let mut r = Rng::new(300 + (q * 10 + n) as u64);
    let g = rand_monic_quintic(&f, &mut r, true); let beta = loop { let b = f.rand_nz(&mut r); if !f.in_fq(b) { break b; } };
    let weak: Poly = (0..6).map(|k| f.mul(g[k], f.pow(beta, 5 - k as u64))).collect();
    let rnd = rand_monic_quintic(&f, &mut r, false);
    let mut curves = vec![("weak", weak), ("random", rnd)];
    for kind in ["x^5+c", "x^5+cx"] { if kind == "x^5+c" && q == 5 { continue; } let c = loop { let c = f.rand_nz(&mut r); if squarefree(&f, &special(&f, kind, c)) { break c; } }; curves.push((kind, special(&f, kind, c))); }
    for (kind, cf) in curves { let lines = census(&f, &cf, false);
        println!("{{\"exp\":\"E3\",{},\"curve\":\"{}\",\"f_idx\":{},\"affine_lines_total\":{},\"weak_lines\":{},\"weak_lines_square_ok\":{},\"weak_scaling_lines\":{}}}",
            field_json(&f), kind, idxpoly(&f, &cf), f.d as u64 * (f.size / q) as u64, lines.len(), lines.iter().filter(|l| l.2).count(), lines.iter().filter(|l| l.1 == 0).count());
    } } }

// ---------------------------------------------------------------- E3b: union of all weak lines on automorphism-rich curves
fn e3b() { for (kind, q, n) in [("x^5+c", 3u32, 4u32), ("x^5+c", 7, 4), ("x^5+c", 11, 5), ("x^5+cx", 5, 4), ("x^5+cx", 7, 4), ("x^5+cx", 13, 4)] {
    let f = Fld::build(q, n); let mut r = Rng::new(400 + (q * 10 + n) as u64 + if kind == "x^5+c" { 0 } else { 7 });
    let shortcut = f.size > 3000; let mut tries = 0;
    let (cf, lines) = loop { tries += 1; let c = f.rand_nz(&mut r); let cf = special(&f, kind, c); if !squarefree(&f, &cf) { continue; }
        let l: Vec<_> = census(&f, &cf, shortcut).into_iter().filter(|l| l.2).collect(); if !l.is_empty() { break (cf, l); } if tries > 200 { break (cf, l); } };
    if lines.is_empty() { println!("{{\"exp\":\"E3b\",{},\"curve\":\"{}\",\"tries\":{},\"weak_lines_square_ok\":0}}", field_json(&f), kind, tries); continue; }
    let b_of = |l: &(u32, u32, bool)| f.el_of[l.1 as usize];
    let w0 = weak_data(&f, &cf, lines[0].0 + 1, b_of(&lines[0])); let fr = Frob::of(q, &w0.psi); let jn = fr.jac(n);
    let mut all: Vec<Div> = vec![]; let mut per = vec![]; let mut l = 1u128; let mut prod = 1u128; let mut forms: HashSet<(i128, i128)> = HashSet::new();
    for ln in &lines { let (be, b) = (ln.0 + 1, b_of(ln)); let fb = fb_from_xs(&f, &cf, &line_xs(&f, be, b)); let (_, _, h, _) = gen_subgroup(&f, &cf, &fb, None);
        let w = weak_data(&f, &cf, be, b); let fw = Frob::of(q, &w.psi); forms.insert((fw.a1, fw.a2));
        per.push(h); l = l / gcd(l, h) * h; prod = prod.saturating_mul(h); all.extend(fb); }
    let (_, _, hu, capped) = gen_subgroup(&f, &cf, &all, Some(jn));
    let coprime = (0..per.len()).all(|i| (0..i).all(|j| gcd(per[i], per[j]) == 1));
    println!("{{\"exp\":\"E3b\",{},\"curve\":\"{}\",\"f_idx\":{},\"tries\":{},\"census\":\"{}\",\"weak_lines_square_ok\":{},\"distinct_frobenius_forms\":{},\"per_line_H\":{:?},\"lcm_per_line\":{},\"product_per_line\":{},\"pairwise_coprime\":{},\"H_union\":{},\"capped\":{},\"J_K_lpoly\":{},\"J_K_factored\":\"{}\",\"H_union_factored\":\"{}\",\"maxprime_J\":{},\"maxprime_H_union\":{},\"maxprime_per_line\":{},\"bound\":{:.1}}}",
        field_json(&f), kind, idxpoly(&f, &cf), tries, if shortcut { "b=-f4/5 per scaling class" } else { "exhaustive affine" }, lines.len(), forms.len(), per, l, prod, coprime, hu, capped, jn, fstr(jn),
        if capped { "capped".into() } else { fstr(hu) }, maxprime(jn), if capped { 0 } else { maxprime(hu) }, per.iter().map(|&h| maxprime(h)).max().unwrap(), bound(q));
} }

// ---------------------------------------------------------------- E4: (2,2)-split weak curve and the pullback E1(K) -> J_C(K)
type Pt = Option<(u32, u32)>;
fn ec_add(f: &Fld, a: u32, b: u32, p: &Pt, q: &Pt) -> Pt { match (p, q) { (None, _) => *q, (_, None) => *p, (Some((x1, y1)), Some((x2, y2))) => {
    let lam = if x1 == x2 { if f.add(*y1, *y2) == 0 { return None; } f.mul(f.add(f.add(f.mul(f.int(3), f.mul(*x1, *x1)), f.mul(f.mul(f.int(2), a), *x1)), b), f.inv(f.mul(f.int(2), *y1))) } else { f.mul(f.sub(*y2, *y1), f.inv(f.sub(*x2, *x1))) };
    let x3 = f.sub(f.sub(f.sub(f.mul(lam, lam), a), *x1), *x2); Some((x3, f.sub(f.mul(lam, f.sub(*x1, x3)), *y1))) } } }
fn ec_mul(f: &Fld, a: u32, b: u32, p: &Pt, mut k: u128) -> Pt { let mut r = None; let mut t = *p; while k > 0 { if k & 1 == 1 { r = ec_add(f, a, b, &r, &t); } t = ec_add(f, a, b, &t, &t); k >>= 1; } r }
fn e4(lo: u32, tag: &str) { for (q, n) in [(7u32, 3u32), (11, 3), (13, 3)] { let f = Fld::build(q, n); let mut r = Rng::new(500 + (q * 10 + n) as u64);
    let mut found = None;
    'srch: for a in lo..q { for b in lo..q { for rr in 1..q { let r2 = rr * rr % q; let c = (q - (r2 * r2 % q * r2 + a * (r2 * r2 % q) + b * r2) % q) % q; if c == 0 { continue; }
        let s: Poly = trim(vec![f.int(c), 0, f.int(b), 0, f.int(a), 0, 1]); if squarefree(&f, &s) { found = Some((a, b, c, rr, s)); break 'srch; } } } }
    let (a, b, c, rr, s) = found.unwrap(); let (ea, eb, ec, er) = (f.int(a), f.int(b), f.int(c), f.int(rr));
    // quintic in t = 1/(x-r): t^6 S(r + 1/t) = sum_k S_k (r t + 1)^k t^(6-k)
    let mut q5: Poly = vec![]; let lin = vec![1u32, er]; let mut lp: Poly = vec![1];
    for k in 0..7 { let mut term = vec![0u32; 6 - k]; term.extend(pscale(&f, &lp, s[k])); q5 = padd(&f, &q5, &trim(term)); lp = pmul(&f, &lp, &lin); }
    assert!(q5.len() == 6); let ld = q5[5]; // = S'(r)
    let psi0: Poly = (0..6).map(|k| f.mul(q5[k], f.pow(ld, 4)) ).enumerate().map(|(k, x)| f.mul(x, f.pow(f.inv(ld), k as u64))).collect(); assert!(psi0[5] == 1);
    let gamma = loop { let g = f.rand_nz(&mut r); if !f.in_fq(f.mul(g, g)) { break g; } }; let beta = f.mul(gamma, gamma);
    let cf: Poly = (0..6).map(|k| f.mul(psi0[k], f.pow(beta, 5 - k as u64))).collect();
    let w = weak_data(&f, &cf, beta, 0); assert!(w.rational && w.square_ok);
    let fr = Frob::of(q, &w.psi); let jn = fr.jac(n); let ncount = jac_order_counts(&f, &cf);
    // E1 : v^2 = u^3 + a u^2 + b u + c ; E2 : v^2 = c u^3 + b u^2 + a u + 1
    let e1p: Poly = vec![ec, eb, ea, 1]; let e2p: Poly = vec![1, ea, eb, ec]; let (n_e1, n_e2) = (count_k(&f, &e1p) as u128, count_k(&f, &e2p) as u128);
    let phi = |p: &Pt| -> Option<Div> { match p { None => Some(dzero()), Some((u0, v0)) => { let dd = f.sub(*u0, f.mul(er, er)); if dd == 0 { return None; } let wv = f.inv(dd);
        // u_T = T^2 - 2 r w L T - w L^2 ; v_T = (v0/L) T^3 mod u_T ; then x = beta T, y = gamma^5 Y'
        let ut: Poly = vec![f.neg(f.mul(wv, f.mul(ld, ld))), f.neg(f.mul(f.int(2), f.mul(er, f.mul(wv, ld)))), 1];
        let vt = pmod(&f, &vec![0, 0, 0, f.mul(*v0, f.inv(ld))], &ut);
        let bi = f.inv(beta); let ux: Poly = (0..3).map(|k| f.mul(ut[k], f.pow(beta, 2 - k as u64))).collect();
        let g5 = f.pow(gamma, 5); let vx: Poly = trim((0..vt.len()).map(|k| f.mul(f.mul(vt[k], g5), f.pow(bi, k as u64))).collect());
        let dv = Div { u: ux, v: vx }; assert!(dvalid(&f, &cf, &dv)); Some(dadd(&f, &cf, &dv, &dzero())) } } };
    let rand_pt = |r: &mut Rng| -> Pt { loop { let u = f.rand(r); let rhs = peval(&f, &e1p, u); if rhs != 0 && f.is_sq(rhs) { return Some((u, f.sqrt(rhs))); } } };
    let (mut hom_ok, mut hom_tested) = (0, 0);
    for _ in 0..200 { let (p, qq) = (rand_pt(&mut r), rand_pt(&mut r)); let s = ec_add(&f, ea, eb, &p, &qq);
        if let (Some(dp), Some(dq), Some(ds)) = (phi(&p), phi(&qq), phi(&s)) { hom_tested += 1; if dadd(&f, &cf, &dp, &dq) == ds { hom_ok += 1; } } }
    let ell = maxprime(n_e1); let cof = n_e1 / ell;
    let (pl, dl) = loop { let p = ec_mul(&f, ea, eb, &rand_pt(&mut r), cof); if p.is_none() { continue; } if let Some(dv) = phi(&p) { break (p, dv); } };
    let fb = fb_from_xs(&f, &cf, &line_xs(&f, beta, 0)); let (_, set, hsz, _) = gen_subgroup(&f, &cf, &fb, Some(ncount));
    let mut hits = 0u64; let mut cur = dl.clone(); for _ in 1..ell { if set.contains(&key(&cur)) { hits += 1; } cur = dadd(&f, &cf, &cur, &dl); }
    // reachable isogeny classes: traces of pi^n for |t| <= 2 sqrt q, against all |t| <= 2 sqrt Q
    let mut tr: HashSet<i128> = HashSet::new(); let bq = (2.0 * (q as f64).sqrt()).floor() as i128;
    for t in -bq..=bq { let (mut s0, mut s1) = (2i128, t); for _ in 1..n { let s2 = t * s1 - q as i128 * s0; s0 = s1; s1 = s2; } tr.insert(s1); }
    let all_tr = 2 * (2.0 * (f.size as f64).sqrt()).floor() as i128 + 1;
    println!("{{\"exp\":\"{}\",{},\"abc_r\":{},\"S\":\"x^6+{}x^4+{}x^2+{}\",\"psi0\":{},\"beta_idx\":{},\"f_idx\":{},\"E1_K\":{},\"E2_K\":{},\"J_K_counts\":{},\"J_K_lpoly\":{},\"J_equals_E1E2\":{},\"hom_tested\":{},\"hom_ok\":{},\"ell\":{},\"E1_K_factored\":\"{}\",\"point_order_ell\":{},\"image_nonzero\":{},\"image_order_ell\":{},\"fb\":{},\"H\":{},\"H_factored\":\"{}\",\"ell_divides_H\":{},\"nonzero_multiples_in_H\":{},\"J_C0_Fq\":{},\"J_C0_Fq2\":{},\"reachable_traces\":{},\"all_traces\":{},\"E1_point_idx\":{}}}",
        tag, field_json(&f), polystr(&[a, b, c, rr]), a, b, c, polystr(&psi0.iter().map(|&x| f.to_int(x)).collect::<Vec<_>>()), f.idx_of[beta as usize], idxpoly(&f, &cf), n_e1, n_e2, ncount, jn,
        ncount == n_e1 * n_e2, hom_tested, hom_ok, ell, fstr(n_e1), ec_mul(&f, ea, eb, &pl, ell).is_none(), dl != dzero(), dmul(&f, &cf, &dl, ell) == dzero(), fb.len(), hsz, fstr(hsz), hsz % ell == 0, hits,
        fr.jac(1), fr.jac(2), tr.len(), all_tr, pl.map(|(x, y)| format!("[{},{}]", f.idx_of[x as usize], f.idx_of[y as usize])).unwrap());
} }

// ---------------------------------------------------------------- E5: the paper's printed examples
struct Big { q: u64, md: Vec<u64> } // K = F_q[w]/(w^n + md), elements low -> high
impl Big {
    fn mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> { let n = self.md.len(); let q = self.q as u128; let mut t = vec![0u128; 2 * n - 1];
        for i in 0..n { for j in 0..n { t[i + j] = (t[i + j] + a[i] as u128 * b[j] as u128) % q; } }
        for k in (n..2 * n - 1).rev() { let c = t[k]; t[k] = 0; for i in 0..n { t[k - n + i] = (t[k - n + i] + c * (q - self.md[i] as u128)) % q; } }
        t[..n].iter().map(|&x| x as u64).collect() }
    fn add(&self, a: &[u64], b: &[u64]) -> Vec<u64> { a.iter().zip(b).map(|(x, y)| (x + y) % self.q).collect() }
    fn scal(&self, a: &[u64], s: u64) -> Vec<u64> { a.iter().map(|&x| (x as u128 * s as u128 % self.q as u128) as u64).collect() }
    fn in_fq(&self, a: &[u64]) -> bool { a[1..].iter().all(|&x| x == 0) }
}
fn hi(v: &[u64]) -> Vec<u64> { v.iter().rev().cloned().collect() } // printed order is high -> low
fn e5_one(name: &str, k: &Big, ft: &[Vec<u64>], mu: &[u64], psi: &[u64], lam0: u64, pts: &[(u64, Vec<u64>)]) {
    let q = k.q; let mu2 = k.mul(mu, mu); let n = k.md.len();
    let psi_ok = (0..6).all(|i| { let p = k.mul(&mu2, &ft[i]); k.in_fq(&p) && p[0] == psi[i] });
    let (mut on, mut infq, mut eq, mut root) = (0, 0, 0, 0);
    for (x, y) in pts { let mut fx = vec![0u64; n]; for i in (0..6).rev() { fx = k.add(&k.scal(&fx, *x), &ft[i]); }
        if k.mul(y, y) == fx { on += 1; }
        let yp = k.mul(mu, y); let lam = (x + lam0) % q; if k.in_fq(&yp) { infq += 1; if yp[0] == lam || yp[0] == (q - lam) % q { eq += 1; } }
        let mut px = 0u128; for i in (0..6).rev() { px = (px * *x as u128 + psi[i] as u128) % q as u128; } if (lam as u128 * lam as u128) % q as u128 == px { root += 1; } }
    println!("{{\"exp\":\"E5\",\"example\":\"{}\",\"q\":{},\"n\":{},\"mu2_ftilde_equals_psi_in_Fq\":{},\"points\":{},\"points_on_Ctilde\":{},\"mu_y_in_Fq\":{},\"mu_y_equals_pm_lambda\":{},\"lambda2_equals_psi_at_x\":{}}}", name, q, n, psi_ok, pts.len(), on, infq, eq, root);
}
fn e5() {
    let k = Big { q: 1354523531, md: vec![1354523482, 7, 104, 38, 1354523496, 1354523513, 1] };
    let top = hi(&[780022053, 760928178, 170989435, 662222802, 702422979, 601651399, 636757065]);
    let ft = vec![hi(&[1250291512, 541668418, 648336894, 93981211, 1209366817, 123242898, 285752852]), hi(&[986949060, 603254005, 469173590, 369803482, 1083057270, 1021615666, 1036052397]),
        hi(&[615249757, 538160793, 134609138, 287039557, 47910319, 159765934, 734765645]), hi(&[789992515, 271076474, 641593893, 60273850, 1118544582, 1168055942, 523388883]), top.clone(), top];
    let mu = hi(&[1096004451, 72204545, 1031970086, 850013668, 667154540, 50674261, 1212050995]);
    let psi = [952059020u64, 506022088, 288749419, 163649504, 140224146, 140224146];
    let pts = vec![(807542330u64, hi(&[716687265, 843761534, 828265313, 159145829, 686345735, 700861992, 509324809])), (620447586, hi(&[813283073, 318977187, 267032198, 35254994, 844429980, 369850187, 214843875])),
        (582043617, hi(&[134707935, 32178915, 41003762, 373212771, 69947485, 533130081, 1094103612])), (548956658, hi(&[314837163, 1345447603, 380080922, 807574817, 985485077, 893517804, 1139064555])),
        (150056870, hi(&[85834289, 675835336, 951227639, 362195701, 972755872, 328519770, 389619817]))];
    e5_one("1", &k, &ft, &mu, &psi, 1133867401, &pts);
    let k = Big { q: 13545235267, md: vec![1, 1, 1, 1] };
    let top = hi(&[8311487827, 668903640, 746585135, 10112627890]);
    let ft = vec![hi(&[8555796245, 3896712384, 2849733199, 7361159113]), hi(&[1613278559, 12653224429, 9436869140, 9250569911]), hi(&[7967576419, 7786992992, 10961177735, 857524719]),
        hi(&[11887251939, 4448133116, 3493295058, 10674514901]), top.clone(), top];
    let mu = hi(&[8480343835, 10345332833, 2496233379, 11776077527]);
    let psi = [9900534212u64, 163035135, 12633136871, 533270342, 10076073713, 10076073713];
    let pts = vec![(13539352170u64, hi(&[9808490675, 11848496349, 10783891627, 11360816798])), (7708342780, hi(&[8476725686, 12354179980, 5399924307, 10535120119])),
        (3683205295, hi(&[2470303628, 6406203858, 9996767999, 350144386])), (1428292391, hi(&[12661012320, 9191091394, 3111319514, 771379433])), (731277897, hi(&[1467582925, 1476361776, 6522982528, 1852301214]))];
    e5_one("2", &k, &ft, &mu, &psi, 5631025385, &pts);
}

fn main() { let arg = std::env::args().nth(1).unwrap_or("all".into());
    // `all` covers the frozen experiments only; the supplementary E4x arm runs only when named ("e4x"),
    // so `all && e4x` records it once (run 20261013T0055Z ran it twice: see PROTOCOL.md amendment 1).
    let run = |n: &str| (arg == "all" && n != "e4x") || arg == n;
    if run("e5") { e5(); } if run("e1") { e1(); } if run("e2") { e2(); } if run("e3") { e3(); } if run("e3b") { e3b(); } if run("e4") { e4(0, "E4"); } if run("e4x") { e4(1, "E4x"); }
}
