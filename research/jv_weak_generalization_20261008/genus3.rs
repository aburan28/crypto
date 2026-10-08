//! Which isogeny classes of E/F_{q^3} have a genus-3 Jacobian isogenous to
//! their Weil restriction, and which of those are Joux–Vitse weak.
//! Standalone: `rustc -O genus3.rs -o genus3 && ./genus3 13 results 200000 200000`.
//! q prime, q ≡ 1 (mod 12).  See PROTOCOL.md.
use std::collections::{HashMap, HashSet};
use std::io::Write;

// ---------------------------------------------------------------- F_q
#[derive(Clone, Copy)] struct Fq { q: u64 }
impl Fq {
    fn add(&self, a: u64, b: u64) -> u64 { let s = a + b; if s >= self.q { s - self.q } else { s } }
    fn sub(&self, a: u64, b: u64) -> u64 { if a >= b { a - b } else { a + self.q - b } }
    fn mul(&self, a: u64, b: u64) -> u64 { (a * b) % self.q }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 { let mut r = 1; while e > 0 { if e & 1 == 1 { r = self.mul(r, a); } a = self.mul(a, a); e >>= 1; } r }
    fn inv(&self, a: u64) -> u64 { self.pow(a, self.q - 2) }
    fn chi(&self, a: u64) -> i64 { if a == 0 { 0 } else if self.pow(a, (self.q - 1) / 2) == 1 { 1 } else { -1 } }
}
// ---------------------------------------------------------------- F_{q^k} = F_q[s]/(s^k - c), k = 2 (c non-square) or 3 (c non-cube)
#[derive(Clone)] struct Ext { f: Fq, k: usize, c: u64, frob_s: Vec<Vec<u64>> /* s^(q*i) as vectors, i<k */, order_m1_u: u64, order_m1_s: u32, nonres: Vec<u64> }
type El = Vec<u64>;
impl Ext {
    fn new(q: u64, k: usize, c: u64) -> Ext {
        let f = Fq { q }; let mut e = Ext { f, k, c, frob_s: vec![], order_m1_u: 0, order_m1_s: 0, nonres: vec![] };
        // s^q = s * (s^k)^((q-1)/k) ... compute by repeated squaring of the element s
        let mut fr = vec![vec![0u64; k]; k]; fr[0][0] = 1;
        if k > 1 { let s: El = { let mut v = vec![0; k]; v[1] = 1; v }; let sq = e.pow_el(&s, q); let mut cur = vec![0u64; k]; cur[0] = 1; for i in 1..k { cur = e.mul(&cur, &sq); fr[i] = cur.clone(); } }
        e.frob_s = fr;
        let order = (q as u128).pow(k as u32) as u64 - 1; let mut u = order; let mut sv = 0; while u % 2 == 0 { u /= 2; sv += 1; }
        e.order_m1_u = u; e.order_m1_s = sv;
        let mut t = 2u64; loop { let x = e.unpack(t); if !e.is_zero(&x) && !e.is_square(&x) { e.nonres = x; break; } t += 1; }
        e
    }
    fn zero(&self) -> El { vec![0; self.k] }
    fn one(&self) -> El { let mut v = vec![0; self.k]; v[0] = 1; v }
    fn from_u(&self, a: u64) -> El { let mut v = vec![0; self.k]; v[0] = a % self.f.q; v }
    fn is_zero(&self, a: &El) -> bool { a.iter().all(|&x| x == 0) }
    fn in_base(&self, a: &El) -> bool { a[1..].iter().all(|&x| x == 0) }
    fn add(&self, a: &El, b: &El) -> El { (0..self.k).map(|i| self.f.add(a[i], b[i])).collect() }
    fn sub(&self, a: &El, b: &El) -> El { (0..self.k).map(|i| self.f.sub(a[i], b[i])).collect() }
    fn neg(&self, a: &El) -> El { (0..self.k).map(|i| self.f.sub(0, a[i])).collect() }
    fn scal(&self, s: u64, a: &El) -> El { (0..self.k).map(|i| self.f.mul(s % self.f.q, a[i])).collect() }
    fn mul(&self, a: &El, b: &El) -> El {
        let k = self.k; let mut t = vec![0u64; 2 * k - 1];
        for i in 0..k { if a[i] == 0 { continue; } for j in 0..k { t[i + j] = self.f.add(t[i + j], self.f.mul(a[i], b[j])); } }
        for d in (k..2 * k - 1).rev() { let v = t[d]; if v != 0 { t[d - k] = self.f.add(t[d - k], self.f.mul(v, self.c)); t[d] = 0; } }
        t.truncate(k); t
    }
    fn sq(&self, a: &El) -> El { self.mul(a, a) }
    fn pow_el(&self, a: &El, mut e: u64) -> El { let mut r = self.one(); let mut b = a.clone(); while e > 0 { if e & 1 == 1 { r = self.mul(&r, &b); } b = self.sq(&b); e >>= 1; } r }
    fn frob(&self, a: &El) -> El { let mut r = self.zero(); for i in 0..self.k { if a[i] != 0 { r = self.add(&r, &self.scal(a[i], &self.frob_s[i])); } } r }
    /// N_{F_{q^k}/F_q}
    fn norm(&self, a: &El) -> u64 { let mut n = a.clone(); let mut cur = a.clone(); for _ in 1..self.k { cur = self.frob(&cur); n = self.mul(&n, &cur); } debug_assert!(self.in_base(&n)); n[0] }
    fn inv(&self, a: &El) -> El { let mut num = self.one(); let mut cur = a.clone(); for _ in 1..self.k { cur = self.frob(&cur); num = self.mul(&num, &cur); } let n = self.mul(a, &num)[0]; self.scal(self.f.inv(n), &num) }
    fn is_square(&self, a: &El) -> bool { if self.is_zero(a) { return true; } let n = self.norm(a); if self.k % 2 == 1 { self.f.chi(n) == 1 } else { self.f.chi(n) >= 0 && self.pow_el(a, ((self.f.q as u128).pow(self.k as u32) as u64 - 1) / 2) == self.one() } }
    fn chi(&self, a: &El) -> i64 { if self.is_zero(a) { 0 } else if self.is_square(a) { 1 } else { -1 } }
    fn sqrt(&self, x: &El) -> El {
        if self.is_zero(x) { return x.clone(); }
        let (u, s) = (self.order_m1_u, self.order_m1_s); let mut m = s; let mut c = self.pow_el(&self.nonres, u); let mut t = self.pow_el(x, u); let mut r = self.pow_el(x, (u + 1) / 2);
        loop { if t == self.one() { return r; } let mut i = 0u32; let mut tt = t.clone(); while tt != self.one() { tt = self.sq(&tt); i += 1; if i == m { panic!("sqrt non-square"); } } let mut b = c.clone(); for _ in 0..(m - i - 1) { b = self.sq(&b); } m = i; c = self.sq(&b); t = self.mul(&t, &c); r = self.mul(&r, &b); }
    }
    fn pack(&self, a: &El) -> u64 { let mut v = 0u64; for i in (0..self.k).rev() { v = v * self.f.q + a[i]; } v }
    fn unpack(&self, mut v: u64) -> El { let mut a = vec![0; self.k]; for i in 0..self.k { a[i] = v % self.f.q; v /= self.f.q; } a }
    fn size(&self) -> u64 { (self.f.q as u128).pow(self.k as u32) as u64 }
}
struct Rng(u64);
impl Rng { fn next(&mut self) -> u64 { let mut x = self.0; x ^= x << 13; x ^= x >> 7; x ^= x << 17; self.0 = x; x } }

// ---------------------------------------------------------------- elliptic curves y^2 = x^3 + a2 x^2 + a4 x + a6 over an Ext
struct Curve { a2: El, a4: El, a6: El }
type Pt = Option<(El, El)>;
fn ec_add(e: &Ext, c: &Curve, p: &Pt, q: &Pt) -> Pt {
    let (x1, y1) = match p { None => return q.clone(), Some(v) => v }; let (x2, y2) = match q { None => return p.clone(), Some(v) => v };
    let m = if x1 == x2 { if y1 != y2 || e.is_zero(y1) { return None; } let num = e.add(&e.add(&e.scal(3, &e.sq(x1)), &e.scal(2, &e.mul(&c.a2, x1))), &c.a4); e.mul(&num, &e.inv(&e.scal(2, y1))) } else { e.mul(&e.sub(y2, y1), &e.inv(&e.sub(x2, x1))) };
    let x3 = e.sub(&e.sub(&e.sub(&e.sq(&m), &c.a2), x1), x2); let y3 = e.sub(&e.mul(&m, &e.sub(x1, &x3)), y1); Some((x3, y3))
}
fn ec_mul(e: &Ext, c: &Curve, p: &Pt, mut k: u64) -> Pt { let mut r: Pt = None; let mut b = p.clone(); while k > 0 { if k & 1 == 1 { r = ec_add(e, c, &r, &b); } b = ec_add(e, c, &b, &b); k >>= 1; } r }
fn rhs(e: &Ext, c: &Curve, x: &El) -> El { e.add(&e.add(&e.add(&e.mul(x, &e.sq(x)), &e.mul(&c.a2, &e.sq(x))), &e.mul(&c.a4, x)), &c.a6) }
fn random_point(e: &Ext, c: &Curve, rng: &mut Rng) -> Pt { loop { let x = e.unpack(rng.next() % e.size()); let r = rhs(e, c, &x); if e.is_zero(&r) { continue; } if e.is_square(&r) { return Some((x, e.sqrt(&r))); } } }
fn order_candidates(e: &Ext, c: &Curve, p: &Pt) -> Vec<u64> {
    let n = e.size(); let n0 = n + 1; let half = 2 * ((n as f64).sqrt() as u64 + 1); let lo = n0 - half; let w = 2 * half + 1; let m = ((w as f64).sqrt() as u64) + 2;
    let mut baby: HashMap<u64, u64> = HashMap::new(); let mut cur: Pt = None; let mut small: Option<u64> = None;
    for b in 0..m { if b > 0 && cur.is_none() { small = Some(b); break; } if let Some((x, _)) = &cur { baby.entry(e.pack(x)).or_insert(b); } cur = ec_add(e, c, &cur, p); }
    if let Some(o) = small { let first = ((lo + o - 1) / o) * o; return (0..).map(|k| first + k * o).take_while(|&v| v < lo + w).collect(); }
    let step = ec_mul(e, c, p, m); let mut r = ec_mul(e, c, p, lo); if let Some((x, y)) = &r { r = Some((x.clone(), e.neg(y))); }
    let mut cands = Vec::new(); let amax = w / m + 1;
    for a in 0..=amax { match &r { None => cands.push(a * m), Some((x, _)) => if let Some(&b) = baby.get(&e.pack(x)) { cands.push(a * m + b); if a * m >= b { cands.push(a * m - b); } } } let ns = match &step { Some((x, y)) => Some((x.clone(), e.neg(y))), None => None }; r = ec_add(e, c, &r, &ns); }
    let mut out: Vec<u64> = cands.into_iter().filter(|&k| k < w).map(|k| lo + k).filter(|&v| ec_mul(e, c, p, v).is_none()).collect(); out.sort(); out.dedup(); out
}
fn curve_order(e: &Ext, c: &Curve, rng: &mut Rng) -> u64 {
    let mut inter: Option<Vec<u64>> = None;
    for _ in 0..8 { let pt = random_point(e, c, rng); let cs = order_candidates(e, c, &pt); inter = Some(match inter { None => cs, Some(v) => v.into_iter().filter(|n| cs.contains(n)).collect() }); let v = inter.as_ref().unwrap(); if v.len() == 1 { return v[0]; } if v.is_empty() { panic!("no candidate"); } }
    let v = inter.unwrap(); let d = v.windows(2).map(|w| w[1] - w[0]).min().unwrap(); let ok: Vec<u64> = v.iter().cloned().filter(|&n| n % d == 0 && d % (n / d) == 0).collect(); if ok.len() == 1 { ok[0] } else { panic!("ambiguous {:?}", v) }
}

// ---------------------------------------------------------------- polynomials over F_q (for squarefree checks) and root counting over Ext
fn fq_poly_rem(f: &Fq, a: &[u64], m: &[u64]) -> Vec<u64> { let mut r = a.to_vec(); let dm = m.len() - 1; let il = f.inv(m[dm]); while r.len() > dm { let dr = r.len() - 1; let c = f.mul(r[dr], il); for i in 0..=dm { r[dr - dm + i] = f.sub(r[dr - dm + i], f.mul(c, m[i])); } while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); } if dr == 0 { break; } } while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); } r }
fn fq_squarefree(f: &Fq, h: &[u64]) -> bool { let d: Vec<u64> = (1..h.len()).map(|i| f.mul(i as u64 % f.q, h[i])).collect(); let mut a = h.to_vec(); let mut b = d; while b.len() > 1 && *b.last().unwrap() == 0 { b.pop(); } if b.len() == 1 && b[0] == 0 { return false; } loop { let r = fq_poly_rem(f, &a, &b); a = b; b = r; if b.len() == 1 && b[0] == 0 { break; } } a.len() == 1 }
/// number of roots in Ext of a polynomial with Ext coefficients (degree ≤ 8): x^N - x mod p, then gcd degree.
fn ext_root_count(e: &Ext, poly: &[El]) -> usize {
    let mut p: Vec<El> = poly.to_vec(); while p.len() > 1 && e.is_zero(p.last().unwrap()) { p.pop(); }
    if p.len() <= 1 { return 0; }
    let dm = p.len() - 1; let il = e.inv(&p[dm]); let m: Vec<El> = p.iter().map(|c| e.mul(c, &il)).collect();
    let rem = |a: &Vec<El>| -> Vec<El> { let mut r = a.clone(); while r.len() > dm { let dr = r.len() - 1; let c = r[dr].clone(); if !e.is_zero(&c) { for i in 0..=dm { r[dr - dm + i] = e.sub(&r[dr - dm + i], &e.mul(&c, &m[i])); } } r.pop(); } r };
    let pmul = |a: &Vec<El>, b: &Vec<El>| -> Vec<El> { let mut r = vec![e.zero(); a.len() + b.len() - 1]; for (i, x) in a.iter().enumerate() { if e.is_zero(x) { continue; } for (j, y) in b.iter().enumerate() { r[i + j] = e.add(&r[i + j], &e.mul(x, y)); } } rem(&r) };
    // x^N mod m
    let n = e.size(); let mut res: Vec<El> = vec![e.one()]; let mut base: Vec<El> = vec![e.zero(), e.one()]; base = rem(&{ let mut b = base; b.resize(dm.max(2), e.zero()); b });
    let mut ex = n; while ex > 0 { if ex & 1 == 1 { res = pmul(&res, &base); } base = pmul(&base, &base); ex >>= 1; }
    // g = x^N - x ; gcd(m, g)
    let mut g = res; g.resize(dm.max(2), e.zero()); g[1] = e.sub(&g[1], &e.one()); while g.len() > 1 && e.is_zero(g.last().unwrap()) { g.pop(); }
    let (mut a, mut b) = (m.clone(), g);
    loop { if b.len() == 1 && e.is_zero(&b[0]) { break; } let db = b.len() - 1; if db == 0 { a = b; b = vec![e.zero()]; continue; } let ilb = e.inv(&b[db]); let mut r = a.clone(); while r.len() > db { let dr = r.len() - 1; let c = e.mul(&r[dr], &ilb); if !e.is_zero(&c) { for i in 0..=db { r[dr - db + i] = e.sub(&r[dr - db + i], &e.mul(&c, &b[i])); } } r.pop(); if r.is_empty() { r.push(e.zero()); } while r.len() > 1 && e.is_zero(r.last().unwrap()) { r.pop(); } } a = b; b = r; }
    a.len() - 1
}

// ---------------------------------------------------------------- genus-3 point counts
/// hyperelliptic y^2 = h(x), coefficients in F_q, over the field e (k = 1 via Ext with k=1 not supported: use count with chi on Ext)
fn hyp_count(e: &Ext, h: &[u64]) -> u64 {
    let n = e.size(); let mut cnt = 0u64; let hc: Vec<El> = h.iter().map(|&c| e.from_u(c)).collect();
    for xv in 0..n { let x = e.unpack(xv); let mut v = e.zero(); for c in hc.iter().rev() { v = e.add(&e.mul(&v, &x), c); } cnt += (1 + e.chi(&v)) as u64; }
    let deg = h.len() - 1; if deg % 2 == 1 { cnt += 1; } else { cnt += (1 + e.chi(&e.from_u(h[deg]))) as u64; }
    cnt
}
/// plane quartic F(x, y, z) = sum c[i][j] x^i y^j z^(4-i-j); count projective points over e by counting roots in y for each x (z = 1), then the line z = 0.
fn quartic_count(e: &Ext, c: &[[u64; 5]; 5]) -> u64 {
    let n = e.size(); let mut cnt = 0u64;
    let cc: Vec<Vec<El>> = (0..5).map(|i| (0..5).map(|j| e.from_u(c[i][j])).collect()).collect();
    for xv in 0..n {
        let x = e.unpack(xv); let mut xp = vec![e.one()]; for _ in 1..5 { let l = xp.last().unwrap().clone(); xp.push(e.mul(&l, &x)); }
        let mut poly: Vec<El> = (0..5).map(|j| { let mut s = e.zero(); for i in 0..(5 - j) { s = e.add(&s, &e.mul(&cc[i][j], &xp[i])); } s }).collect();
        while poly.len() > 1 && e.is_zero(poly.last().unwrap()) { poly.pop(); }
        if poly.len() == 1 { if e.is_zero(&poly[0]) { cnt += n; } continue; }
        cnt += ext_root_count(e, &poly) as u64;
    }
    // z = 0: F(x, y, 0) = sum_{i+j=4} c[i][j] x^i y^j; points [x:y:0]: y = 1: roots in x of sum c[i][4-i] x^i ; plus [1:0:0] if c[4][0] == 0
    let mut pz: Vec<El> = (0..5).map(|i| e.from_u(c[i][4 - i])).collect(); while pz.len() > 1 && e.is_zero(pz.last().unwrap()) { pz.pop(); }
    if pz.len() == 1 { if e.is_zero(&pz[0]) { cnt += n; } } else { cnt += ext_root_count(e, &pz) as u64; }
    if c[4][0] == 0 { cnt += 1; }
    cnt
}
/// Evaluate F and its three partials at a projective point of e, without pow_el.
fn quartic_grad(e: &Ext, c: &[[u64; 5]; 5], x: &El, y: &El, z: &El) -> (El, El, El, El) {
    let mut xp = vec![e.one()]; let mut yp = vec![e.one()]; let mut zp = vec![e.one()];
    for _ in 1..5 { xp.push(e.mul(xp.last().unwrap(), x)); yp.push(e.mul(yp.last().unwrap(), y)); zp.push(e.mul(zp.last().unwrap(), z)); }
    let (mut v, mut dx, mut dy, mut dz) = (e.zero(), e.zero(), e.zero(), e.zero());
    for i in 0..5 { for j in 0..(5 - i) { let k = 4 - i - j; let co = c[i][j]; if co == 0 { continue; }
        v = e.add(&v, &e.scal(co, &e.mul(&e.mul(&xp[i], &yp[j]), &zp[k])));
        if i > 0 { dx = e.add(&dx, &e.scal(co * i as u64, &e.mul(&e.mul(&xp[i - 1], &yp[j]), &zp[k]))); }
        if j > 0 { dy = e.add(&dy, &e.scal(co * j as u64, &e.mul(&e.mul(&xp[i], &yp[j - 1]), &zp[k]))); }
        if k > 0 { dz = e.add(&dz, &e.scal(co * k as u64, &e.mul(&e.mul(&xp[i], &yp[j]), &zp[k - 1]))); } } }
    (v, dx, dy, dz)
}
/// Singular over e?  Affine chart: for each x, the quartic f(y) = F(x, y, 1) and f' share a root only at
/// O(1) values of x (vertical tangents and singular points); brute-force y there.  Then the line z = 0.
fn quartic_singular_over(e: &Ext, c: &[[u64; 5]; 5]) -> bool {
    let n = e.size();
    let cc: Vec<Vec<El>> = (0..5).map(|i| (0..5).map(|j| e.from_u(c[i][j])).collect()).collect();
    let gcd_nontrivial = |a: &Vec<El>, b: &Vec<El>| -> bool {
        let (mut p, mut r) = (a.clone(), b.clone());
        let trim = |v: &mut Vec<El>| { while v.len() > 1 && e.is_zero(v.last().unwrap()) { v.pop(); } };
        trim(&mut p); trim(&mut r);
        if r.len() == 1 && e.is_zero(&r[0]) { return p.len() > 1; }
        loop { if r.len() == 1 && e.is_zero(&r[0]) { break; } let dr = r.len() - 1; if dr == 0 { return false; } let il = e.inv(&r[dr]); let mut rem = p.clone();
            while rem.len() > dr { let d = rem.len() - 1; let q = e.mul(&rem[d], &il); if !e.is_zero(&q) { for i in 0..=dr { rem[d - dr + i] = e.sub(&rem[d - dr + i], &e.mul(&q, &r[i])); } } rem.pop(); if rem.is_empty() { rem.push(e.zero()); } trim(&mut rem); }
            p = r; r = rem; }
        p.len() > 1
    };
    for xv in 0..n {
        let x = e.unpack(xv); let mut xp = vec![e.one()]; for _ in 1..5 { xp.push(e.mul(xp.last().unwrap(), &x)); }
        let fy: Vec<El> = (0..5).map(|j| { let mut s = e.zero(); for i in 0..(5 - j) { s = e.add(&s, &e.mul(&cc[i][j], &xp[i])); } s }).collect();
        let dfy: Vec<El> = (1..5).map(|j| e.scal(j as u64, &fy[j])).collect();
        if !gcd_nontrivial(&fy, &dfy) { continue; }
        for yv in 0..n { let y = e.unpack(yv); let (v, dx, dy, dz) = quartic_grad(e, c, &x, &y, &e.one()); if e.is_zero(&v) && e.is_zero(&dx) && e.is_zero(&dy) && e.is_zero(&dz) { return true; } }
    }
    for xv in 0..n { let x = e.unpack(xv); let (v, dx, dy, dz) = quartic_grad(e, c, &x, &e.one(), &e.zero()); if e.is_zero(&v) && e.is_zero(&dx) && e.is_zero(&dy) && e.is_zero(&dz) { return true; } }
    let (v, dx, dy, dz) = quartic_grad(e, c, &e.one(), &e.zero(), &e.zero()); e.is_zero(&v) && e.is_zero(&dx) && e.is_zero(&dy) && e.is_zero(&dz)
}

fn v2f_of(t: i64, q3: u64) -> (u32, i64) { // (v2 of conductor, D_K mod 8)
    let d = t * t - 4 * q3 as i64; if d == 0 { return (64, 99); }
    let mut n = (-d) as u64; let mut f0 = 1u64; let mut sqf: i64 = -1; let mut dd = 2u64;
    while dd * dd <= n { let mut e = 0; while n % dd == 0 { n /= dd; e += 1; } for _ in 0..(e / 2) { f0 *= dd; } if e % 2 == 1 { sqf *= dd as i64; } dd += if dd == 2 { 1 } else { 2 }; }
    if n > 1 { sqf *= n as i64; }
    let (dk, cond) = if sqf.rem_euclid(4) == 1 { (sqf, f0) } else { (4 * sqf, f0 / 2) };
    (cond.trailing_zeros(), dk.rem_euclid(8))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let q: u64 = args[1].parse().unwrap(); let outdir = args.get(2).cloned().unwrap_or("results".into());
    let n_h: u64 = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(100000); let n_c: u64 = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(100000);
    let seed: u64 = args.get(5).map(|s| s.parse().unwrap()).unwrap_or(20261008);
    assert!(q % 12 == 1, "q ≡ 1 mod 12");
    let t0 = std::time::Instant::now();
    let f = Fq { q };
    let nonsq = (2..q).find(|&x| f.chi(x) == -1).unwrap();
    let noncube = (2..q).find(|&x| f.pow(x, (q - 1) / 3) != 1).unwrap();
    let e2 = Ext::new(q, 2, nonsq); let e3 = Ext::new(q, 3, noncube);
    // tower self-checks
    { let x = e3.unpack(12345 % e3.size()); assert!(e3.mul(&x, &e3.inv(&x)) == e3.one()); assert!(e3.frob(&e3.frob(&e3.frob(&x))) == x); assert!(e3.pow_el(&x, q) == e3.frob(&x)); let s = e3.sqrt(&e3.sq(&x)); assert!(s == x || s == e3.neg(&x));
      let y = e2.unpack(777 % e2.size()); assert!(e2.mul(&y, &e2.inv(&y)) == e2.one()); assert!(e2.pow_el(&y, q) == e2.frob(&y)); }
    let q3 = e3.size();
    eprintln!("q={q}: F_q^3 size {q3}, tower ok ({:.1}s)", t0.elapsed().as_secs_f64());

    // 1. full-2-torsion classes over F_{q^3}: Legendre enumeration, weak flag, traces by components
    let mut lambda_of: Vec<u32> = vec![u32::MAX; q3 as usize]; let mut nodes: Vec<(u32, u32)> = Vec::new(); // (j index, lambda index)
    let one = e3.one();
    let j_of = |lam: &El| -> El { let l2 = e3.sq(lam); let num = e3.add(&e3.sub(&l2, lam), &one); let num3 = e3.mul(&num, &e3.sq(&num)); let den = e3.mul(&l2, &e3.sq(&e3.sub(lam, &one))); e3.mul(&e3.scal(256, &num3), &e3.inv(&den)) };
    for lk in 2..q3 { let lam = e3.unpack(lk); if lam == one { continue; } let jk = e3.pack(&j_of(&lam)) as usize; if lambda_of[jk] == u32::MAX { lambda_of[jk] = lk as u32; nodes.push((jk as u32, lk as u32)); } }
    let nn = nodes.len();
    let mut weak = vec![false; nn]; let mut node_of: HashMap<u32, u32> = HashMap::new(); for (i, (jk, _)) in nodes.iter().enumerate() { node_of.insert(*jk, i as u32); }
    let mut edges: Vec<Vec<u32>> = vec![Vec::new(); nn];
    for (id, (_, lk)) in nodes.iter().enumerate() {
        let lam = e3.unpack(*lk as u64);
        let c1 = lam.clone(); let c2 = e3.sub(&one, &lam); let c3 = e3.mul(&lam, &e3.inv(&e3.sub(&lam, &one)));
        weak[id] = [c1, c2, c3].iter().any(|c| e3.norm(c) == 1);
        let e = [e3.zero(), one.clone(), lam.clone()];
        for i in 0..3 { let (ei, ej, ek) = (&e[i], &e[(i + 1) % 3], &e[(i + 2) % 3]); let a = e3.sub(&e3.sub(&e3.scal(2, ei), ej), ek); let b = e3.mul(&e3.sub(ei, ej), &e3.sub(ei, ek)); if !e3.is_square(&b) { continue; }
            let a2 = e3.neg(&e3.scal(2, &a)); let b2 = e3.sub(&e3.sq(&a), &e3.scal(4, &b));
            // j of y^2 = x(x^2 + a2 x + b2) = 256 (a2^2 - 3 b2)^3 / (b2^2 (a2^2 - 4 b2))
            let a22 = e3.sq(&a2); let num = e3.scal(256, &e3.pow_el(&e3.sub(&a22, &e3.scal(3, &b2)), 3)); let den = e3.mul(&e3.sq(&b2), &e3.sub(&a22, &e3.scal(4, &b2)));
            let jj = e3.pack(&e3.mul(&num, &e3.inv(&den))) as u32; if let Some(&nb) = node_of.get(&jj) { edges[id].push(nb); } }
    }
    let mut parent: Vec<u32> = (0..nn as u32).collect(); fn find(p: &mut Vec<u32>, mut x: u32) -> u32 { while p[x as usize] != x { let g = p[p[x as usize] as usize]; p[x as usize] = g; x = g; } x }
    for id in 0..nn { for &nb in &edges[id] { let (a, b) = (find(&mut parent, id as u32), find(&mut parent, nb)); if a != b { parent[a as usize] = b; } } }
    let mut comp_rep: HashMap<u32, u32> = HashMap::new(); for id in 0..nn as u32 { let r = find(&mut parent, id); comp_rep.entry(r).or_insert(id); }
    let mut rng = Rng(seed | 1); let mut trace_of_root: HashMap<u32, i64> = HashMap::new();
    for (r, &rep) in &comp_rep { let lam = e3.unpack(nodes[rep as usize].1 as u64); let c = Curve { a2: e3.neg(&e3.add(&one, &lam)), a4: lam.clone(), a6: e3.zero() }; let n = curve_order(&e3, &c, &mut rng); trace_of_root.insert(*r, q3 as i64 + 1 - n as i64); }
    let mut w_set: HashSet<i64> = HashSet::new(); let mut t2_set: HashSet<i64> = HashSet::new();
    for id in 0..nn as u32 { let t = trace_of_root[&find(&mut parent, id)].abs(); t2_set.insert(t); if weak[id as usize] { w_set.insert(t); } }
    eprintln!("full-2-torsion j: {nn}, weak {}, classes |T2| = {}, weak classes |W| = {} ({:.1}s)", weak.iter().filter(|&&w| w).count(), t2_set.len(), w_set.len(), t0.elapsed().as_secs_f64());
    // all traces: sample random curves y^2 = x^3 + a x + b over F_{q^3}
    let mut t_all: HashSet<i64> = HashSet::new();
    for _ in 0..4000 { let a = e3.unpack(rng.next() % q3); let b = e3.unpack(rng.next() % q3); let disc = e3.add(&e3.scal(4, &e3.pow_el(&a, 3)), &e3.scal(27, &e3.sq(&b))); if e3.is_zero(&disc) { continue; } let c = Curve { a2: e3.zero(), a4: a, a6: b }; let n = curve_order(&e3, &c, &mut rng); t_all.insert((q3 as i64 + 1 - n as i64).abs()); }
    eprintln!("sampled distinct |t| over all curves: {} ({:.1}s)", t_all.len(), t0.elapsed().as_secs_f64());

    // 2. hyperelliptic genus 3 over F_q
    let mut h_set: HashMap<i64, u64> = HashMap::new(); let mut h_sig = 0u64; let mut h_tried = 0u64; let mut h_a1 = 0u64;
    let e1 = Ext::new(q, 1, 0); // k = 1 field wrapper: s^1 = 0 → elements are base; mul works with k=1
    for _ in 0..n_h {
        let deg = if rng.next() & 1 == 0 { 7 } else { 8 }; let mut h: Vec<u64> = (0..=deg).map(|_| rng.next() % q).collect(); if h[deg] == 0 { h[deg] = 1; }
        if !fq_squarefree(&f, &h) { continue; } h_tried += 1;
        let n1 = hyp_count(&e1, &h); if n1 != q + 1 { continue; } h_a1 += 1;
        let n2 = hyp_count(&e2, &h); if n2 != q * q + 1 { continue; } h_sig += 1;
        let n3 = hyp_count(&e3, &h); let num = q3 as i64 + 1 - n3 as i64; if num % 3 != 0 { eprintln!("hyp: #C(F_q3) not ≡ q^3+1 mod 3: {n3}"); continue; } *h_set.entry((num / 3).abs()).or_default() += 1;
    }
    eprintln!("hyperelliptic: tried {h_tried}, a1=0: {h_a1}, signature: {h_sig}, distinct |t|: {} ({:.1}s)", h_set.len(), t0.elapsed().as_secs_f64());
    // 3. plane quartics over F_q
    let mut qset: HashMap<i64, u64> = HashMap::new(); let mut q_sig = 0u64; let mut q_tried = 0u64; let mut q_a1 = 0u64; let mut q_sing = 0u64;
    for _ in 0..n_c {
        let mut c = [[0u64; 5]; 5]; for i in 0..5 { for j in 0..(5 - i) { c[i][j] = rng.next() % q; } }
        q_tried += 1;
        let n1 = quartic_count(&e1, &c); if n1 != q + 1 { continue; } q_a1 += 1;
        if quartic_singular_over(&e1, &c) { q_sing += 1; continue; }
        let n2 = quartic_count(&e2, &c); if n2 != q * q + 1 { continue; }
        if quartic_singular_over(&e2, &c) { q_sing += 1; continue; }
        q_sig += 1;
        let n3 = quartic_count(&e3, &c); let num = q3 as i64 + 1 - n3 as i64; if num % 3 != 0 { eprintln!("quartic: #C(F_q3) not ≡ q^3+1 mod 3: {n3}"); continue; } *qset.entry((num / 3).abs()).or_default() += 1;
    }
    eprintln!("quartics: tried {q_tried}, a1=0: {q_a1}, singular (over F_q or F_q2): {q_sing}, signature: {q_sig}, distinct |t|: {} ({:.1}s)", qset.len(), t0.elapsed().as_secs_f64());
    // 4. report
    std::fs::create_dir_all(&outdir).unwrap();
    let hk: HashSet<i64> = h_set.keys().cloned().collect(); let qk: HashSet<i64> = qset.keys().cloned().collect();
    let h_minus_w: Vec<i64> = { let mut v: Vec<i64> = hk.difference(&w_set).cloned().collect(); v.sort(); v };
    let q_minus_w: Vec<i64> = { let mut v: Vec<i64> = qk.difference(&w_set).cloned().collect(); v.sort(); v };
    let q_minus_wh: Vec<i64> = { let mut v: Vec<i64> = qk.iter().filter(|t| !w_set.contains(t) && !hk.contains(t)).cloned().collect(); v.sort(); v };
    let q_depth1: Vec<i64> = qk.iter().filter(|&&t| v2f_of(t, q3).0 == 1).cloned().collect();
    let h_depth1: Vec<i64> = hk.iter().filter(|&&t| v2f_of(t, q3).0 == 1).cloned().collect();
    let union_wq: HashSet<i64> = w_set.union(&qk).cloned().collect();
    let mut lines = vec![format!("## genus-3 signature census, q = {q}, |F_q^3| = {q3}\n"),
        format!("| set | size |\n|:--|--:|\n| W (JV-weak traces, exhaustive) | {} |\n| T2 (traces of full-2-torsion classes, exhaustive) | {} |\n| T (distinct traces in a 4,000-curve sample) | {} |\n| H (hyperelliptic genus-3 signature traces, {h_sig} curves) | {} |\n| Q (plane-quartic signature traces, {q_sig} curves) | {} |\n| H ∖ W | {} |\n| Q ∖ W | {} |\n| Q ∖ (W ∪ H) | {} |\n| traces in Q with v2(f) = 1 | {} |\n| traces in H with v2(f) = 1 | {} |\n| (W ∪ Q) ∩ T / T | {:.3} |",
            w_set.len(), t2_set.len(), t_all.len(), hk.len(), qk.len(), h_minus_w.len(), q_minus_w.len(), q_minus_wh.len(), q_depth1.len(), h_depth1.len(), union_wq.intersection(&t_all).count() as f64 / t_all.len().max(1) as f64)];
    lines.push(format!("\nhyperelliptic: tried {h_tried}, a1 = 0: {h_a1}, signature {h_sig}.  quartics: tried {q_tried}, a1 = 0: {q_a1}, singular over F_q or F_q²: {q_sing}, signature {q_sig}."));
    lines.push("\n| trace in Q ∖ W | v2(f) | D_K mod 8 | in T2 | quartics |\n|--:|--:|--:|:--|--:|".into());
    for t in &q_minus_w { let (v, dk) = v2f_of(*t, q3); lines.push(format!("| {t} | {v} | {dk} | {} | {} |", t2_set.contains(t), qset[t])); }
    lines.push("\n| trace in H ∖ W | v2(f) | D_K mod 8 | in T2 | curves |\n|--:|--:|--:|:--|--:|".into());
    for t in &h_minus_w { let (v, dk) = v2f_of(*t, q3); lines.push(format!("| {t} | {v} | {dk} | {} | {} |", t2_set.contains(t), h_set[t])); }
    let text = lines.join("\n");
    std::fs::write(format!("{outdir}/genus3_q{q}.md"), &text).unwrap();
    let mut fo = std::fs::File::create(format!("{outdir}/genus3_q{q}_sets.json")).unwrap();
    let mut wv: Vec<i64> = w_set.iter().cloned().collect(); wv.sort(); let mut hv: Vec<(i64, u64)> = h_set.iter().map(|(a, b)| (*a, *b)).collect(); hv.sort(); let mut qv: Vec<(i64, u64)> = qset.iter().map(|(a, b)| (*a, *b)).collect(); qv.sort(); let mut tv: Vec<i64> = t_all.iter().cloned().collect(); tv.sort();
    let j = |v: &Vec<(i64, u64)>| v.iter().map(|(t, c)| format!("[{t},{c}]")).collect::<Vec<_>>().join(",");
    writeln!(fo, "{{\"q\":{q},\"W\":{:?},\"H\":[{}],\"Q\":[{}],\"T_sample\":{:?},\"seconds\":{:.1}}}", wv, j(&hv), j(&qv), tv, t0.elapsed().as_secs_f64()).unwrap();
    println!("{text}");
}
