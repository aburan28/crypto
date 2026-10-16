//! Exhaustive census of full-2-torsion curves over F_{p^6}: which isogeny
//! classes hold a Joux–Vitse weak curve, and where in the 2-volcano the weak
//! curves sit.  See PROTOCOL.md (frozen before this file was written).
//!
//! Standalone: `rustc -O census.rs -o census && ./census 7 results/`.
//! Modes: `P DIR` (census), `P --selftest`, `tabulate DIR P…`, `patterns P DIR`,
//! `walk P DIR [starts]` (PROTOCOL-2-mechanism-and-walk.md).
//! Field tower: F_{p^2} = F_p[i]/(i^2 - n), F_{p^6} = F_{p^2}[s]/(s^3 - c),
//! sigma = Frobenius of F_{p^2}-degree, s -> omega s.

use std::collections::HashMap;
use std::io::Write;

// ---------------------------------------------------------------- F_p
#[derive(Clone, Copy)]
struct Fp { p: u64 }
impl Fp {
    fn add(&self, a: u64, b: u64) -> u64 { let s = a + b; if s >= self.p { s - self.p } else { s } }
    fn sub(&self, a: u64, b: u64) -> u64 { if a >= b { a - b } else { a + self.p - b } }
    fn mul(&self, a: u64, b: u64) -> u64 { (a * b) % self.p }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 { let mut r = 1; while e > 0 { if e & 1 == 1 { r = self.mul(r, a); } a = self.mul(a, a); e >>= 1; } r }
    fn inv(&self, a: u64) -> u64 { self.pow(a, self.p - 2) }
    fn is_qr(&self, a: u64) -> bool { a != 0 && self.pow(a, (self.p - 1) / 2) == 1 }
}

// ---------------------------------------------------------------- F_{p^2}
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct E2 { a: u64, b: u64 } // a + b i, i^2 = n
#[derive(Clone, Copy)]
struct F2 { f: Fp, n: u64 }
impl F2 {
    fn zero() -> E2 { E2 { a: 0, b: 0 } }
    fn one() -> E2 { E2 { a: 1, b: 0 } }
    fn add(&self, x: E2, y: E2) -> E2 { E2 { a: self.f.add(x.a, y.a), b: self.f.add(x.b, y.b) } }
    fn sub(&self, x: E2, y: E2) -> E2 { E2 { a: self.f.sub(x.a, y.a), b: self.f.sub(x.b, y.b) } }
    fn neg(&self, x: E2) -> E2 { E2 { a: self.f.sub(0, x.a), b: self.f.sub(0, x.b) } }
    fn mul(&self, x: E2, y: E2) -> E2 {
        let ac = self.f.mul(x.a, y.a); let bd = self.f.mul(x.b, y.b);
        let ad = self.f.mul(x.a, y.b); let bc = self.f.mul(x.b, y.a);
        E2 { a: self.f.add(ac, self.f.mul(self.n, bd)), b: self.f.add(ad, bc) }
    }
    fn scal(&self, k: u64, x: E2) -> E2 { E2 { a: self.f.mul(k, x.a), b: self.f.mul(k, x.b) } }
    fn norm(&self, x: E2) -> u64 { self.f.sub(self.f.mul(x.a, x.a), self.f.mul(self.n, self.f.mul(x.b, x.b))) }
    fn inv(&self, x: E2) -> E2 { let ni = self.f.inv(self.norm(x)); E2 { a: self.f.mul(x.a, ni), b: self.f.mul(self.f.sub(0, x.b), ni) } }
    fn pow(&self, mut x: E2, mut e: u64) -> E2 { let mut r = Self::one(); while e > 0 { if e & 1 == 1 { r = self.mul(r, x); } x = self.mul(x, x); e >>= 1; } r }
    fn is_zero(&self, x: E2) -> bool { x.a == 0 && x.b == 0 }
}

// ---------------------------------------------------------------- F_{p^6}
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct E6 { c: [E2; 3] } // c0 + c1 s + c2 s^2, s^3 = cc
struct F6 { f2: F2, cc: E2, omega: E2, omega2: E2, p: u64, order_m1_u: u64, order_m1_s: u32, nonres: E6 }
impl F6 {
    fn zero(&self) -> E6 { E6 { c: [F2::zero(); 3] } }
    fn one(&self) -> E6 { E6 { c: [F2::one(), F2::zero(), F2::zero()] } }
    fn from2(&self, x: E2) -> E6 { E6 { c: [x, F2::zero(), F2::zero()] } }
    fn from_u(&self, k: u64) -> E6 { self.from2(E2 { a: k % self.p, b: 0 }) }
    fn is_zero(&self, x: E6) -> bool { x.c.iter().all(|z| self.f2.is_zero(*z)) }
    fn add(&self, x: E6, y: E6) -> E6 { E6 { c: [self.f2.add(x.c[0], y.c[0]), self.f2.add(x.c[1], y.c[1]), self.f2.add(x.c[2], y.c[2])] } }
    fn sub(&self, x: E6, y: E6) -> E6 { E6 { c: [self.f2.sub(x.c[0], y.c[0]), self.f2.sub(x.c[1], y.c[1]), self.f2.sub(x.c[2], y.c[2])] } }
    fn neg(&self, x: E6) -> E6 { E6 { c: [self.f2.neg(x.c[0]), self.f2.neg(x.c[1]), self.f2.neg(x.c[2])] } }
    fn mul(&self, x: E6, y: E6) -> E6 {
        let g = &self.f2;
        let (x0, x1, x2) = (x.c[0], x.c[1], x.c[2]); let (y0, y1, y2) = (y.c[0], y.c[1], y.c[2]);
        let z0 = g.add(g.mul(x0, y0), g.mul(self.cc, g.add(g.mul(x1, y2), g.mul(x2, y1))));
        let z1 = g.add(g.add(g.mul(x0, y1), g.mul(x1, y0)), g.mul(self.cc, g.mul(x2, y2)));
        let z2 = g.add(g.add(g.mul(x0, y2), g.mul(x1, y1)), g.mul(x2, y0));
        E6 { c: [z0, z1, z2] }
    }
    fn sq(&self, x: E6) -> E6 { self.mul(x, x) }
    fn scal(&self, k: u64, x: E6) -> E6 { E6 { c: [self.f2.scal(k, x.c[0]), self.f2.scal(k, x.c[1]), self.f2.scal(k, x.c[2])] } }
    fn frob(&self, x: E6) -> E6 { E6 { c: [x.c[0], self.f2.mul(self.omega, x.c[1]), self.f2.mul(self.omega2, x.c[2])] } }
    /// N_{F_{p^6}/F_{p^2}}(x) = x * sigma(x) * sigma^2(x).
    fn norm2(&self, x: E6) -> E2 {
        let s1 = self.frob(x); let s2 = self.frob(s1);
        let n = self.mul(x, self.mul(s1, s2));
        debug_assert!(self.f2.is_zero(n.c[1]) && self.f2.is_zero(n.c[2]));
        n.c[0]
    }
    fn inv(&self, x: E6) -> E6 {
        let s1 = self.frob(x); let s2 = self.frob(s1);
        let num = self.mul(s1, s2);
        let n = self.mul(x, num).c[0];
        let ni = self.f2.inv(n);
        E6 { c: [self.f2.mul(num.c[0], ni), self.f2.mul(num.c[1], ni), self.f2.mul(num.c[2], ni)] }
    }
    fn is_square(&self, x: E6) -> bool {
        if self.is_zero(x) { return true; }
        let n2 = self.norm2(x); let n1 = self.f2.norm(n2); self.f2.f.is_qr(n1)
    }
    fn pow(&self, mut x: E6, mut e: u64) -> E6 { let mut r = self.one(); while e > 0 { if e & 1 == 1 { r = self.mul(r, x); } x = self.sq(x); e >>= 1; } r }
    /// Tonelli–Shanks; caller guarantees x is a square.
    fn sqrt(&self, x: E6) -> E6 {
        if self.is_zero(x) { return x; }
        let (u, s) = (self.order_m1_u, self.order_m1_s);
        let mut m = s;
        let mut c = self.pow(self.nonres, u);
        let mut t = self.pow(x, u);
        let mut r = self.pow(x, (u + 1) / 2);
        loop {
            if t == self.one() { return r; }
            let mut i = 0u32; let mut tt = t;
            while tt != self.one() { tt = self.sq(tt); i += 1; if i == m { panic!("sqrt of non-square"); } }
            let mut b = c; for _ in 0..(m - i - 1) { b = self.sq(b); }
            m = i; c = self.sq(b); t = self.mul(t, c); r = self.mul(r, b);
        }
    }
    fn pack(&self, x: E6) -> u64 {
        let p = self.p;
        (((((x.c[0].a * p + x.c[0].b) * p + x.c[1].a) * p + x.c[1].b) * p + x.c[2].a) * p) + x.c[2].b
    }
    fn unpack(&self, mut k: u64) -> E6 {
        let p = self.p;
        let c2b = k % p; k /= p; let c2a = k % p; k /= p; let c1b = k % p; k /= p; let c1a = k % p; k /= p; let c0b = k % p; k /= p; let c0a = k % p;
        E6 { c: [E2 { a: c0a, b: c0b }, E2 { a: c1a, b: c1b }, E2 { a: c2a, b: c2b }] }
    }
}

fn build_field(p: u64) -> F6 {
    let f = Fp { p };
    let n = (2..p).find(|&x| !f.is_qr(x)).expect("qnr");
    let f2 = F2 { f, n };
    let m2 = p * p - 1; assert!(m2 % 3 == 0);
    // a non-cube in F_{p^2}
    let mut cc = E2 { a: 0, b: 0 };
    'outer: for a in 0..p { for b in 0..p { let x = E2 { a, b }; if f2.is_zero(x) { continue; } if f2.pow(x, m2 / 3) != F2::one() { cc = x; break 'outer; } } }
    assert!(!f2.is_zero(cc));
    let omega = f2.pow(cc, m2 / 3); let omega2 = f2.mul(omega, omega);
    assert!(f2.mul(omega2, omega) == F2::one() && omega != F2::one());
    let order = p.pow(6) - 1; let mut u = order; let mut s = 0u32; while u % 2 == 0 { u /= 2; s += 1; }
    let mut fld = F6 { f2, cc, omega, omega2, p, order_m1_u: u, order_m1_s: s, nonres: E6 { c: [F2::zero(); 3] } };
    // a quadratic non-residue of F_{p^6}
    let mut k = 2u64;
    loop { let x = fld.unpack(k); if !fld.is_zero(x) && !fld.is_square(x) { fld.nonres = x; break; } k += 1; }
    // self-checks of the tower
    let x = fld.unpack(12345 % order); let y = fld.unpack(6789 % order);
    assert!(fld.mul(x, fld.inv(x)) == fld.one());
    assert!(fld.mul(fld.add(x, y), fld.sub(x, y)) == fld.sub(fld.sq(x), fld.sq(y)));
    let f3 = fld.frob(fld.frob(fld.frob(x))); assert!(f3 == x);
    assert!(fld.pow(x, p * p) == fld.frob(x), "frobenius mismatch");
    let sqx = fld.sqrt(fld.sq(x)); assert!(sqx == x || sqx == fld.neg(x));
    fld
}

// ---------------------------------------------------------------- curve y^2 = x^3 + A x^2 + B x
#[derive(Clone, Copy)]
struct Curve { a: E6, b: E6 }
type Pt = Option<(E6, E6)>;
fn ec_neg(f: &F6, p: Pt) -> Pt { p.map(|(x, y)| (x, f.neg(y))) }
fn ec_add(f: &F6, c: &Curve, p: Pt, q: Pt) -> Pt {
    let (x1, y1) = match p { None => return q, Some(v) => v };
    let (x2, y2) = match q { None => return p, Some(v) => v };
    let m = if x1 == x2 {
        if y1 != y2 || f.is_zero(y1) { return None; }
        // (3x^2 + 2Ax + B) / (2y)
        let num = f.add(f.add(f.scal(3, f.sq(x1)), f.scal(2, f.mul(c.a, x1))), c.b);
        f.mul(num, f.inv(f.scal(2, y1)))
    } else { f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1))) };
    let x3 = f.sub(f.sub(f.sub(f.sq(m), c.a), x1), x2);
    let y3 = f.sub(f.mul(m, f.sub(x1, x3)), y1);
    Some((x3, y3))
}
fn ec_mul(f: &F6, c: &Curve, p: Pt, mut k: u64) -> Pt {
    let mut r: Pt = None; let mut base = p;
    while k > 0 { if k & 1 == 1 { r = ec_add(f, c, r, base); } base = ec_add(f, c, base, base); k >>= 1; }
    r
}

struct Rng(u64);
impl Rng { fn next(&mut self) -> u64 { let mut x = self.0; x ^= x << 13; x ^= x >> 7; x ^= x << 17; self.0 = x; x } }

fn random_point(f: &F6, c: &Curve, rng: &mut Rng) -> Pt {
    let order = f.p.pow(6);
    loop {
        let x = f.unpack(rng.next() % order);
        let rhs = f.mul(x, f.add(f.add(f.sq(x), f.mul(c.a, x)), c.b));
        if f.is_zero(rhs) { continue; }
        if f.is_square(rhs) { return Some((x, f.sqrt(rhs))); }
    }
}

/// All N in the Hasse interval with [N]P = O (candidate group orders).
fn order_candidates(f: &F6, c: &Curve, p: Pt) -> Vec<u64> {
    let pp = f.p; let q3 = pp.pow(6); let n0 = q3 + 1; let half = 2 * pp.pow(3);
    let lo = n0 - half; let w = 2 * half + 1;
    let m = ((w as f64).sqrt() as u64) + 2;
    let mut baby: HashMap<u64, (u64, u64)> = HashMap::with_capacity(m as usize * 2);
    let mut cur: Pt = None;
    let mut small_order: Option<u64> = None;
    for b in 0..m {
        if b > 0 && cur.is_none() { small_order = Some(b); break; }
        if let Some((x, y)) = cur { baby.entry(f.pack(x)).or_insert((b, f.pack(y))); }
        cur = ec_add(f, c, cur, p);
    }
    if let Some(ord) = small_order {
        // every N in the interval that ord divides
        let first = ((lo + ord - 1) / ord) * ord;
        return (0..).map(|k| first + k * ord).take_while(|&n| n < lo + w).collect();
    }
    let step = ec_mul(f, c, p, m);
    let mut r = ec_neg(f, ec_mul(f, c, p, lo));
    let mut cands = Vec::new();
    let amax = w / m + 1;
    for a in 0..=amax {
        match r {
            None => cands.push(a * m),
            Some((x, _y)) => if let Some(&(b, _yb)) = baby.get(&f.pack(x)) {
                cands.push(a * m + b); if a * m >= b { cands.push(a * m - b); }
            },
        }
        r = ec_add(f, c, r, ec_neg(f, step));
    }
    let mut out: Vec<u64> = cands.into_iter().filter(|&k| k < w).map(|k| lo + k).filter(|&n| ec_mul(f, c, p, n).is_none()).collect();
    out.sort(); out.dedup(); out
}

fn curve_order(f: &F6, c: &Curve, rng: &mut Rng) -> u64 {
    let mut inter: Option<Vec<u64>> = None;
    for _ in 0..8 {
        let pt = random_point(f, c, rng);
        let cs = order_candidates(f, c, pt);
        inter = Some(match inter { None => cs, Some(v) => v.into_iter().filter(|n| cs.contains(n)).collect() });
        let v = inter.as_ref().unwrap();
        if v.len() == 1 { return v[0]; }
        if v.is_empty() { panic!("no order candidate"); }
    }
    // Several candidates remain: every sampled point has order dividing the
    // spacing e, so E(F) = Z/n1 x Z/n2 with n2 | e... and #E = e*d with d | e.
    let v = inter.unwrap();
    let e = v.windows(2).map(|w| w[1] - w[0]).min().unwrap();
    let ok: Vec<u64> = v.iter().cloned().filter(|&n| n % e == 0 && e % (n / e) == 0).collect();
    if ok.len() == 1 { return ok[0]; }
    panic!("ambiguous order: candidates {:?}, spacing {e}, structure-filtered {:?}", v, ok);
}

/// Brute-force #E for the self-test.
fn brute_order(f: &F6, c: &Curve) -> u64 {
    let q3 = f.p.pow(6); let mut n = 1u64; // infinity
    for k in 0..q3 { let x = f.unpack(k); let rhs = f.mul(x, f.add(f.add(f.sq(x), f.mul(c.a, x)), c.b)); if f.is_zero(rhs) { n += 1; } else if f.is_square(rhs) { n += 2; } }
    n
}

/// j-invariant of y^2 = x^3 + a x^2 + b x.
fn j_inv(f: &F6, a: E6, b: E6) -> E6 {
    let a2 = f.sq(a);
    let num = f.scal(256, f.pow(f.sub(a2, f.scal(3, b)), 3));
    let den = f.mul(f.sq(b), f.sub(a2, f.scal(4, b)));
    f.mul(num, f.inv(den))
}

fn factor_u(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new(); let mut d = 2u64;
    while d * d <= n { let mut e = 0; while n % d == 0 { n /= d; e += 1; } if e > 0 { out.push((d, e)); } d += if d == 2 { 1 } else { 2 }; }
    if n > 1 { out.push((n, 1)); } out
}



// ---------------------------------------------------------------- polynomials over F6 (for the 3-division polynomial)
type Poly = Vec<E6>; // little-endian coefficients
fn ptrim(f: &F6, a: &mut Poly) { while a.len() > 1 && f.is_zero(*a.last().unwrap()) { a.pop(); } }
fn pmul(f: &F6, a: &Poly, b: &Poly) -> Poly {
    let mut r = vec![f.zero(); a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() { if f.is_zero(x) { continue; } for (j, &y) in b.iter().enumerate() { r[i + j] = f.add(r[i + j], f.mul(x, y)); } }
    ptrim(f, &mut r); r
}
fn psub(f: &F6, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len()); let mut r = vec![f.zero(); n];
    for i in 0..n { let x = if i < a.len() { a[i] } else { f.zero() }; let y = if i < b.len() { b[i] } else { f.zero() }; r[i] = f.sub(x, y); }
    ptrim(f, &mut r); r
}
fn prem(f: &F6, a: &Poly, m: &Poly) -> Poly {
    let mut r = a.clone(); ptrim(f, &mut r);
    let dm = m.len() - 1; let inv_lead = f.inv(m[dm]);
    while r.len() - 1 >= dm && !(r.len() == 1 && f.is_zero(r[0])) {
        let dr = r.len() - 1; let c = f.mul(r[dr], inv_lead);
        for i in 0..=dm { r[dr - dm + i] = f.sub(r[dr - dm + i], f.mul(c, m[i])); }
        ptrim(f, &mut r); if dr == 0 { break; }
    }
    r
}
fn pmonic(f: &F6, a: &Poly) -> Poly { let il = f.inv(*a.last().unwrap()); a.iter().map(|&c| f.mul(c, il)).collect() }
fn pgcd(f: &F6, a: &Poly, b: &Poly) -> Poly {
    let (mut x, mut y) = (a.clone(), b.clone()); ptrim(f, &mut x); ptrim(f, &mut y);
    while !(y.len() == 1 && f.is_zero(y[0])) { let r = prem(f, &x, &y); x = y; y = r; }
    pmonic(f, &x)
}
fn ppowmod(f: &F6, base: &Poly, mut e: u64, m: &Poly) -> Poly {
    let mut r = vec![f.one()]; let mut b = prem(f, base, m);
    while e > 0 { if e & 1 == 1 { r = prem(f, &pmul(f, &r, &b), m); } b = prem(f, &pmul(f, &b, &b), m); e >>= 1; }
    r
}
/// Roots in F of a polynomial over F (Cantor–Zassenhaus), degree ≤ 30.
fn roots_in_f(f: &F6, poly: &Poly, rng: &mut Rng) -> Vec<E6> {
    let q3 = f.p.pow(6);
    let m = pmonic(f, poly);
    let x = vec![f.zero(), f.one()];
    let xq = ppowmod(f, &x, q3, &m);
    let g = pgcd(f, &m, &psub(f, &xq, &x)); // product of the distinct linear factors
    let mut out = Vec::new();
    fn split(f: &F6, g: &Poly, q3: u64, rng: &mut Rng, out: &mut Vec<E6>) {
        if g.len() <= 1 { return; }
        if g.len() == 2 { out.push(f.neg(g[0])); return; } // monic x + c
        loop {
            let r = f.unpack(rng.next() % q3);
            let base = vec![r, f.one()];
            let pw = ppowmod(f, &base, (q3 - 1) / 2, g);
            let h = pgcd(f, g, &psub(f, &pw, &vec![f.one()]));
            if h.len() > 1 && h.len() < g.len() {
                let other = pmonic(f, &pdiv_exact(f, g, &h));
                split(f, &h, q3, rng, out); split(f, &other, q3, rng, out); return;
            }
        }
    }
    split(f, &g, q3, rng, &mut out);
    out
}
fn pdiv_exact(f: &F6, a: &Poly, b: &Poly) -> Poly {
    let mut r = a.clone(); let db = b.len() - 1; let il = f.inv(b[db]);
    let mut qv = vec![f.zero(); a.len() - db];
    for k in (0..qv.len()).rev() { let c = f.mul(r[k + db], il); qv[k] = c; for i in 0..=db { r[k + i] = f.sub(r[k + i], f.mul(c, b[i])); } }
    qv
}

// ---------------------------------------------------------------- the census as a struct
struct Census {
    p: u64, q3: u64, f: F6,
    lambda_of: Vec<u32>, node_of: Vec<u32>, nodes_j: Vec<u32>,
    weak: Vec<bool>, edges: Vec<[(u32, bool); 3]>, trace: Vec<i64>, height: Vec<u32>,
    v2f_of: HashMap<i64, u32>, class_nodes: HashMap<i64, Vec<u32>>,
}
impl Census {
    fn lambda(&self, id: u32) -> E6 { self.f.unpack(self.lambda_of[self.nodes_j[id as usize] as usize] as u64) }
    fn class_key(&self, id: u32) -> i64 { self.trace[id as usize].abs() }
    /// Legendre node of the curve with 2-torsion abscissae e (full 2-torsion assumed).
    fn node_of_roots(&self, e: [E6; 3]) -> Option<u32> {
        let f = &self.f;
        let lam = f.mul(f.sub(e[2], e[0]), f.inv(f.sub(e[1], e[0])));
        let a = f.neg(f.add(f.one(), lam)); let j = j_inv(f, a, lam);
        let n = self.node_of[f.pack(j) as usize]; if n == u32::MAX { None } else { Some(n) }
    }
    /// Legendre curve of node id: y^2 = x(x-1)(x-lam), a2 = -(1+lam), a4 = lam.
    fn curve(&self, id: u32) -> (Curve, E6) { let lam = self.lambda(id); (Curve { a: self.f.neg(self.f.add(self.f.one(), lam)), b: lam }, lam) }
    /// Rational 3-isogeny codomains (as node ids) of node id.
    fn three_moves(&self, id: u32, rng: &mut Rng) -> (Vec<u32>, usize) {
        let f = &self.f; let (c, lam) = self.curve(id);
        // psi_3 = 3x^4 + 4 a2 x^3 + 6 a4 x^2 - a4^2
        let psi = vec![f.neg(f.sq(c.b)), f.zero(), f.scal(6, c.b), f.scal(4, c.a), f.from_u(3)];
        let roots = roots_in_f(f, &psi, rng);
        let e = [f.zero(), f.one(), lam];
        let mut out = Vec::new(); let mut missing = 0usize;
        for x0 in roots {
            let y0sq = f.mul(x0, f.add(f.add(f.sq(x0), f.mul(c.a, x0)), c.b));
            let mut e2 = [f.zero(); 3];
            for i in 0..3 {
                let d = f.sub(e[i], x0);
                let term = f.mul(f.scal(2, y0sq), f.inv(f.sq(d)));
                e2[i] = f.add(f.sub(f.sub(f.neg(e[i]), f.scal(2, c.a)), f.scal(4, x0)), term);
            }
            match self.node_of_roots(e2) { Some(n) => out.push(n), None => missing += 1 }
        }
        (out, missing)
    }
    fn two_moves(&self, id: u32) -> Vec<u32> {
        self.edges[id as usize].iter().filter(|(_, k)| *k).map(|(jk, _)| self.node_of[*jk as usize]).collect()
    }
    /// P = (x, y) halvable iff x - e_j are all squares.
    fn halvable_x(&self, e: &[E6; 3], x: E6) -> bool { e.iter().all(|&ej| self.f.is_square(self.f.sub(x, ej))) }
    /// abscissae of the halves of the point with abscissa x (both signs of y): up to 4 values.
    fn halves_x(&self, e: &[E6; 3], x: E6) -> Vec<E6> {
        let f = &self.f;
        let s: Vec<E6> = e.iter().map(|&ej| f.sqrt(f.sub(x, ej))).collect();
        let mut out: Vec<E6> = Vec::new();
        for mask in 0..8u32 {
            let t: Vec<E6> = (0..3).map(|j| if (mask >> j) & 1 == 1 { f.neg(s[j]) } else { s[j] }).collect();
            let v = f.add(x, f.add(f.add(f.mul(t[0], t[1]), f.mul(t[0], t[2])), f.mul(t[1], t[2])));
            if !out.contains(&v) { out.push(v); }
        }
        out
    }
    /// largest k ≤ kmax with T_i ∈ 2^k E(F), by breadth-first halving on abscissae.
    fn depth2(&self, e: &[E6; 3], i: usize, kmax: u32) -> u32 {
        let mut cur = vec![e[i]]; let mut k = 0u32;
        while k < kmax {
            let mut next = Vec::new();
            for &x in &cur { if self.halvable_x(e, x) { for h in self.halves_x(e, x) { if !next.contains(&h) { next.push(h); } } } }
            if next.is_empty() { break; }
            k += 1; cur = next; if cur.len() > 64 { cur.truncate(64); }
        }
        k
    }
    /// intrinsic height class: 1 if a 2-isogeny descends to the floor, 2 if some codomain is at intrinsic height 1, else 3 (meaning ≥ 3).
    fn h_intr(&self, id: u32) -> u32 {
        let km = self.two_moves(id); if km.len() < 3 { return 1; }
        if km.iter().any(|&nb| self.two_moves(nb).len() < 3) { return 2; }
        3
    }
    fn pattern(&self, id: u32) -> u32 {
        let lam = self.lambda(id); let e = [self.f.zero(), self.f.one(), lam];
        (0..3).filter(|&i| self.halvable_x(&e, e[i]) ).count() as u32
    }
    /// intrinsic 2-adic height estimate: 1 + min_i depth2(T_i), capped.
    fn a_est(&self, id: u32, kmax: u32) -> u32 {
        let lam = self.lambda(id); let e = [self.f.zero(), self.f.one(), lam];
        1 + (0..3).map(|i| self.depth2(&e, i, kmax)).min().unwrap()
    }
}

fn build_census(p: u64, seed: u64, verbose: bool) -> Census {
    let t0 = std::time::Instant::now();
    let f = build_field(p);
    let q3 = p.pow(6);
    let mut lambda_of: Vec<u32> = vec![u32::MAX; q3 as usize];
    let mut node_of: Vec<u32> = vec![u32::MAX; q3 as usize];
    let mut nodes_j: Vec<u32> = Vec::new();
    let one = f.one();
    for lk in 2..q3 {
        let lam = f.unpack(lk); if lam == one { continue; }
        let a = f.neg(f.add(one, lam)); let j = j_inv(&f, a, lam); let jk = f.pack(j) as usize;
        if lambda_of[jk] == u32::MAX { lambda_of[jk] = lk as u32; node_of[jk] = nodes_j.len() as u32; nodes_j.push(jk as u32); }
    }
    let nn = nodes_j.len();
    if verbose { eprintln!("p={p}: full-2-torsion j: {nn} ({:.1}s)", t0.elapsed().as_secs_f64()); }
    let mut weak = vec![false; nn];
    let mut edges: Vec<[(u32, bool); 3]> = vec![[(u32::MAX, false); 3]; nn];
    for (id, &jk) in nodes_j.iter().enumerate() {
        let lam = f.unpack(lambda_of[jk as usize] as u64);
        let c1 = lam; let c2 = f.sub(one, lam); let c3 = f.mul(lam, f.inv(f.sub(lam, one)));
        weak[id] = [c1, c2, c3].iter().any(|&c| f.norm2(c) == F2::one());
        let e = [f.zero(), one, lam];
        for i in 0..3 {
            let (ei, ej, ek) = (e[i], e[(i + 1) % 3], e[(i + 2) % 3]);
            let a = f.sub(f.sub(f.scal(2, ei), ej), ek); let b = f.mul(f.sub(ei, ej), f.sub(ei, ek));
            let a2 = f.neg(f.scal(2, a)); let b2 = f.sub(f.sq(a), f.scal(4, b));
            edges[id][i] = (f.pack(j_inv(&f, a2, b2)) as u32, f.is_square(b));
        }
    }
    let mut parent: Vec<u32> = (0..nn as u32).collect();
    fn find(par: &mut Vec<u32>, mut x: u32) -> u32 { while par[x as usize] != x { let g = par[par[x as usize] as usize]; par[x as usize] = g; x = g; } x }
    let mut floor_adjacent = vec![false; nn];
    for id in 0..nn { for i in 0..3 { let (jk, keeps) = edges[id][i]; if keeps { let nb = node_of[jk as usize]; assert!(nb != u32::MAX); let (ra, rb) = (find(&mut parent, id as u32), find(&mut parent, nb)); if ra != rb { parent[ra as usize] = rb; } } else { floor_adjacent[id] = true; } } }
    let mut comp_members: HashMap<u32, Vec<u32>> = HashMap::new();
    for id in 0..nn as u32 { let r = find(&mut parent, id); comp_members.entry(r).or_default().push(id); }
    let comps: Vec<(u32, Vec<u32>)> = comp_members.into_iter().collect();
    let nthreads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(14);
    let chunk = (comps.len() + nthreads - 1) / nthreads;
    let mut trace_of_comp: HashMap<u32, i64> = HashMap::new();
    std::thread::scope(|sc| {
        let mut handles = Vec::new();
        for (ti, part) in comps.chunks(chunk.max(1)).enumerate() {
            let fref = &f; let lam_ref = &lambda_of; let nodes_ref = &nodes_j;
            handles.push(sc.spawn(move || {
                let mut out = Vec::with_capacity(part.len());
                let mut rng = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ ((ti as u64 + 1) << 32) | 1);
                for (root, members) in part {
                    let lam = fref.unpack(lam_ref[nodes_ref[members[0] as usize] as usize] as u64);
                    let c = Curve { a: fref.neg(fref.add(fref.one(), lam)), b: lam };
                    let n = curve_order(fref, &c, &mut rng);
                    out.push((*root, (q3 as i64 + 1) - n as i64));
                }
                out
            }));
        }
        for h in handles { for (r, t) in h.join().unwrap() { trace_of_comp.insert(r, t); } }
    });
    let mut trace: Vec<i64> = vec![0; nn];
    for id in 0..nn as u32 { let r = find(&mut parent, id); trace[id as usize] = trace_of_comp[&r]; }
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); nn];
    for id in 0..nn { for i in 0..3 { let (jk, keeps) = edges[id][i]; if keeps { let nb = node_of[jk as usize]; adj[id].push(nb); adj[nb as usize].push(id as u32); } } }
    let mut height: Vec<u32> = vec![u32::MAX; nn];
    let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
    for id in 0..nn { if floor_adjacent[id] { height[id] = 1; queue.push_back(id as u32); } }
    while let Some(v) = queue.pop_front() { let h = height[v as usize]; for &w in &adj[v as usize] { if height[w as usize] == u32::MAX { height[w as usize] = h + 1; queue.push_back(w); } } }
    let mut v2f_of: HashMap<i64, u32> = HashMap::new(); let mut class_nodes: HashMap<i64, Vec<u32>> = HashMap::new();
    for id in 0..nn as u32 {
        let t = trace[id as usize].abs(); class_nodes.entry(t).or_default().push(id);
        v2f_of.entry(t).or_insert_with(|| { let d = t * t - 4 * q3 as i64; if d == 0 { return 64; } let fac = factor_u((-d) as u64); let mut f0 = 1u64; let mut sqf: i64 = -1; for (pr, e) in &fac { for _ in 0..(e / 2) { f0 *= pr; } if e % 2 == 1 { sqf *= *pr as i64; } } let cond = if sqf.rem_euclid(4) == 1 { f0 } else { f0 / 2 }; cond.trailing_zeros() });
    }
    if verbose { eprintln!("p={p}: traces and heights done ({:.1}s)", t0.elapsed().as_secs_f64()); }
    Census { p, q3, f, lambda_of, node_of, nodes_j, weak, edges, trace, height, v2f_of, class_nodes }
}

// ---------------------------------------------------------------- mode: patterns
fn run_patterns(c: &Census, outdir: &str) {
    let p = c.p; let nn = c.nodes_j.len();
    let mut by_hp: HashMap<(u32, u32), (u64, u64)> = HashMap::new();
    let mut one_halv_dir: HashMap<(bool, i64), u64> = HashMap::new(); // (below crater?, direction) for the unique halvable point's isogeny
    let mut pat3_min_h = u32::MAX; let mut pat3 = 0u64;
    let mut aest_eq = 0u64; let mut aest_ne = 0u64; let mut aest_tab: HashMap<(u32, u32), u64> = HashMap::new();
    let mut weak_h1_pat0 = 0u64;
    let mut crater_by_dk8: HashMap<(i64, bool), (u64, u64)> = HashMap::new(); // (D_K mod 8, at crater?) -> (nodes, weak)
    let mut h1_by_depth: HashMap<(bool, u32), (u64, u64)> = HashMap::new();   // (depth==1?, pattern) at height 1 -> (nodes, weak)
    let mut pat3_rule: HashMap<(bool, &'static str), u64> = HashMap::new(); let mut pat3_by_h: HashMap<(u32, &'static str), u64> = HashMap::new();     // pattern-3 nodes: (below crater?, rule outcome)
    for id in 0..nn as u32 {
        let h = c.height[id as usize]; if h == u32::MAX { continue; }
        let lam = c.lambda(id); let e = [c.f.zero(), c.f.one(), lam];
        let halv: Vec<bool> = (0..3).map(|i| c.halvable_x(&e, e[i])).collect();
        let pat = halv.iter().filter(|&&b| b).count() as u32;
        let key = c.class_key(id); let v2f = c.v2f_of[&key];
        let ent = by_hp.entry((h, pat)).or_default(); ent.0 += 1; if c.weak[id as usize] { ent.1 += 1; }
        if pat == 3 { pat3 += 1; pat3_min_h = pat3_min_h.min(h); }
        if c.weak[id as usize] && h == 1 && pat == 0 && v2f >= 2 { weak_h1_pat0 += 1; }
        let dk8 = { let d = key * key - 4 * c.q3 as i64; if d == 0 { 99 } else { let fac = factor_u((-d) as u64); let mut sqf: i64 = -1; for (pr, e) in &fac { if e % 2 == 1 { sqf *= *pr as i64; } } (if sqf.rem_euclid(4) == 1 { sqf } else { 4 * sqf }).rem_euclid(8) } };
        { let ent = crater_by_dk8.entry((dk8, h == v2f)).or_default(); ent.0 += 1; if c.weak[id as usize] { ent.1 += 1; } }
        if h == 1 { let ent = h1_by_depth.entry((v2f == 1, pat)).or_default(); ent.0 += 1; if c.weak[id as usize] { ent.1 += 1; } }
        if pat == 3 {
            let depths: Vec<u32> = (0..3).map(|i| c.depth2(&e, i, 4)).collect();
            let mx = *depths.iter().max().unwrap(); let nmax = depths.iter().filter(|&&d| d == mx).count();
            let outcome: &'static str = if nmax != 1 { "no unique max" } else {
                let i = depths.iter().position(|&d| d == mx).unwrap();
                let (jk, keeps) = c.edges[id as usize][i];
                if !keeps { "max-depth point -> floor" } else { let dh = c.height[c.node_of[jk as usize] as usize] as i64 - h as i64; if dh == 1 { "max-depth point ascends" } else if dh == 0 { "max-depth point level" } else { "max-depth point descends" } }
            };
            let outcome2: &'static str = if outcome != "no unique max" { outcome } else {
                // fallback: the codomain with strictly larger a_est than the other two
                let mut best: Vec<(u32, usize)> = (0..3).filter_map(|i| { let (jk, keeps) = c.edges[id as usize][i]; if keeps { Some((c.a_est(c.node_of[jk as usize], 4), i)) } else { None } }).collect();
                best.sort(); best.reverse();
                if best.len() >= 2 && best[0].0 > best[1].0 { let i = best[0].1; let (jk, _) = c.edges[id as usize][i]; let dh = c.height[c.node_of[jk as usize] as usize] as i64 - h as i64; if dh == 1 { "tie -> a_est fallback ascends" } else if dh == 0 { "tie -> a_est fallback level" } else { "tie -> a_est fallback descends" } } else { "tie -> unresolved" }
            };
            *pat3_rule.entry((h < v2f, outcome2)).or_default() += 1;
            *pat3_by_h.entry((h, outcome2)).or_default() += 1;
        }
        if pat == 1 {
            let i = halv.iter().position(|&b| b).unwrap();
            let (jk, keeps) = c.edges[id as usize][i];
            let dir: i64 = if !keeps { -99 } else { c.height[c.node_of[jk as usize] as usize] as i64 - h as i64 };
            *one_halv_dir.entry((h < v2f, dir)).or_default() += 1;
        }
        let a = c.a_est(id, 3.min(v2f.max(1)));
        *aest_tab.entry((h.min(4), a.min(4))).or_default() += 1;
        if a == h.min(4) { aest_eq += 1; } else { aest_ne += 1; }
    }
    let mut lines = Vec::new();
    lines.push(format!("## patterns, p = {p}\n\n| height | pattern | nodes | weak | weak fraction × q |\n|--:|--:|--:|--:|--:|"));
    let mut keys: Vec<(u32, u32)> = by_hp.keys().cloned().collect(); keys.sort();
    for k in keys { let (n, w) = by_hp[&k]; lines.push(format!("| {} | {} | {n} | {w} | {:.2} |", k.0, k.1, w as f64 / n as f64 * (p * p) as f64)); }
    lines.push(format!("\nM1 — unique halvable point's isogeny, (below crater, direction: +1 up, −1 down, 0 level, −99 floor) → count:"));
    let mut dk: Vec<(bool, i64)> = one_halv_dir.keys().cloned().collect(); dk.sort();
    for k in dk { lines.push(format!("- below_crater={} dir={} : {}", k.0, k.1, one_halv_dir[&k])); }
    lines.push(format!("\nM2 — weak nodes at height 1 with pattern 0 in depth ≥ 2 classes: {weak_h1_pat0}"));
    lines.push("\nCrater test — (D_K mod 8, at crater) → nodes, weak, weak fraction × q:".into());
    let mut ck: Vec<(i64, bool)> = crater_by_dk8.keys().cloned().collect(); ck.sort();
    for k in ck { let (n, w) = crater_by_dk8[&k]; lines.push(format!("- D_K mod 8 = {}, at crater = {} : {n} nodes, {w} weak, {:.2}", k.0, k.1, w as f64 / n as f64 * (p * p) as f64)); }
    lines.push("\nHeight 1 by class depth — (depth = 1, pattern) → nodes, weak:".into());
    let mut hk: Vec<(bool, u32)> = h1_by_depth.keys().cloned().collect(); hk.sort();
    for k in hk { let (n, w) = h1_by_depth[&k]; lines.push(format!("- depth1 = {}, pattern {} : {n} nodes, {w} weak", k.0, k.1)); }
    lines.push("\nPattern-3 ascent rule — (below crater, outcome) → count:".into());
    let mut pk: Vec<(bool, &str)> = pat3_rule.keys().cloned().collect(); pk.sort();
    for k in pk { lines.push(format!("- below_crater = {} : {} : {}", k.0, k.1, pat3_rule[&k])); }
    lines.push("\nPattern-3 ascent rule by height — (height, outcome) → count:".into());
    let mut hk2: Vec<(u32, &str)> = pat3_by_h.keys().cloned().collect(); hk2.sort();
    for k in hk2 { lines.push(format!("- height {} : {} : {}", k.0, k.1, pat3_by_h[&k])); }
    lines.push(format!("M3 — pattern-3 nodes: {pat3}, minimum height {}", if pat3_min_h == u32::MAX { 0 } else { pat3_min_h }));
    lines.push(format!("\nIntrinsic estimate a_est = 1 + min_i depth2(T_i) (cap 3) against height (both capped at 4): equal {aest_eq}, different {aest_ne}"));
    let mut ak: Vec<(u32, u32)> = aest_tab.keys().cloned().collect(); ak.sort();
    lines.push("| height | a_est | nodes |\n|--:|--:|--:|".into());
    for k in ak { lines.push(format!("| {} | {} | {} |", k.0, k.1, aest_tab[&k])); }
    let text = lines.join("\n");
    std::fs::create_dir_all(outdir).unwrap();
    std::fs::write(format!("{outdir}/patterns_p{p}.md"), &text).unwrap();
    println!("{text}");
}

// ---------------------------------------------------------------- mode: walk
#[derive(Default, Clone)]
struct WalkStats { starts: u64, success: u64, refused: u64, steps_success: Vec<u64>, steps_all: Vec<u64>, ascent_fail: u64, ascent_steps_max: u64, missing_codomain: u64, trace_mismatch: u64 }

fn ascend_rule(c: &Census, node: u32) -> Option<u32> {
    let lam = c.lambda(node); let e = [c.f.zero(), c.f.one(), lam];
    let halv: Vec<bool> = (0..3).map(|i| c.halvable_x(&e, e[i])).collect(); let pat = halv.iter().filter(|&&b| b).count();
    let kernel = if pat == 1 { halv.iter().position(|&b| b) } else if pat == 3 {
        let depths: Vec<u32> = (0..3).map(|i| c.depth2(&e, i, 4)).collect(); let mx = *depths.iter().max().unwrap();
        if depths.iter().filter(|&&d| d == mx).count() == 1 { depths.iter().position(|&d| d == mx) } else {
            let mut best: Vec<(u32, usize)> = (0..3).filter_map(|i| { let (jk, keeps) = c.edges[node as usize][i]; if keeps { Some((c.a_est(c.node_of[jk as usize], 4), i)) } else { None } }).collect();
            best.sort(); best.reverse(); if best.len() >= 2 && best[0].0 > best[1].0 { Some(best[0].1) } else { None } }
    } else { None };
    let mut next = kernel.and_then(|i| { let (jk, keeps) = c.edges[node as usize][i]; if keeps { Some(c.node_of[jk as usize]) } else { None } });
    if next.is_none() { let km = c.two_moves(node); if km.len() == 1 { next = Some(km[0]); } }
    next
}

fn walk_once(c: &Census, start: u32, policy: u8, rng: &mut Rng, st: &mut WalkStats) {
    let policy_n = policy == 1; let policy_a = policy == 2; let q = c.p * c.p; let cap = 3 * q; let key = c.class_key(start); let v2f = c.v2f_of[&key];
    st.starts += 1;
    if (policy_n || policy_a) && v2f == 1 { st.refused += 1; st.steps_all.push(0); return; }
    let mut cur = start; let mut steps = 0u64; let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new(); seen.insert(cur);
    let mut since_new = 0u64;
    let check = |c: &Census, n: u32, st: &mut WalkStats| { if c.trace[n as usize].abs() != key { st.trace_mismatch += 1; } };
    macro_rules! visit { ($n:expr) => {{ cur = $n; steps += 1; check(c, cur, st); if seen.insert(cur) { since_new = 0; } else { since_new += 1; } if c.weak[cur as usize] { st.success += 1; st.steps_success.push(steps); st.steps_all.push(steps); return; } }} }
    if c.weak[cur as usize] { st.success += 1; st.steps_success.push(0); st.steps_all.push(0); return; }
    let mut target = if v2f == 64 { 1 } else { 3.min(v2f) };
    let mut asc_total = 0u64;
    loop {
        if policy_n || policy_a {
            // ascend to the target level by the intrinsic rules
            let mut guard = 0;
            while c.h_intr(cur) < target && steps < cap {
                match ascend_rule(c, cur) { Some(nb) => { asc_total += 1; visit!(nb); }, None => { st.ascent_fail += 1; break; } }
                guard += 1; if guard > 6 { st.ascent_fail += 1; break; }
            }
            st.ascent_steps_max = st.ascent_steps_max.max(asc_total);
        }
        // level walk (N) or free walk (R)
        let mut exhausted = false;
        while steps < cap {
            let (three, missing) = c.three_moves(cur, rng); st.missing_codomain += missing as u64;
            let mut moves: Vec<(u32, Option<u32>)> = three.into_iter().map(|n| (n, None)).collect();
            if policy_n {
                let at_crater = c.h_intr(cur) >= 3 && target == v2f || (target < 3 && target == v2f);
                let dist = ascend_rule(c, cur);
                for nb in c.two_moves(cur) {
                    if Some(nb) == dist { if at_crater { moves.push((nb, None)); } continue; }
                    if let Some(back) = ascend_rule(c, nb) { if back != cur { moves.push((back, Some(nb))); } }
                }
            } else { for nb in c.two_moves(cur) { moves.push((nb, None)); } }
            if moves.is_empty() { exhausted = true; break; }
            let (dest, via) = moves[(rng.next() % moves.len() as u64) as usize];
            if let Some(mid) = via { visit!(mid); }
            visit!(dest);
            if since_new > 50 * seen.len() as u64 { exhausted = true; break; }
        }
        if !policy_n || !exhausted || target <= 1 || steps >= cap { break; }
        // N: this level is exhausted; lower the target and descend by a non-distinguished 2-move
        target -= 1; since_new = 0;
        let dist = ascend_rule(c, cur);
        let down: Vec<u32> = c.two_moves(cur).into_iter().filter(|&nb| Some(nb) != dist).collect();
        if down.is_empty() { break; }
        visit!(down[(rng.next() % down.len() as u64) as usize]);
    }
    st.steps_all.push(steps);
}

fn median(v: &mut Vec<u64>) -> f64 { if v.is_empty() { return f64::NAN; } v.sort(); let n = v.len(); if n % 2 == 1 { v[n / 2] as f64 } else { (v[n / 2 - 1] + v[n / 2]) as f64 / 2.0 } }

fn run_walk(c: &Census, outdir: &str, starts_per_class: u64, seed: u64) {
    let p = c.p;
    let mut keys: Vec<i64> = c.class_nodes.keys().cloned().collect(); keys.sort();
    let names = ["R (random)", "N (refuse, ascend, walk level, lower on exhaustion)", "A (refuse, ascend, then walk freely)"];
    let mut on_weak: Vec<WalkStats> = vec![WalkStats::default(); 3]; let mut on_not: Vec<WalkStats> = vec![WalkStats::default(); 3];
    let mut refused_weak = [0u64; 3];
    let mut both: Vec<(Vec<u64>, Vec<u64>)> = vec![(Vec::new(), Vec::new()); 3]; // vs R
    for &t in &keys {
        let members = &c.class_nodes[&t]; let holds_weak = members.iter().any(|&id| c.weak[id as usize]);
        let mut rng = Rng(seed ^ ((t as u64).wrapping_mul(0x9E3779B97F4A7C15)) | 1);
        for _ in 0..starts_per_class {
            let start = members[(rng.next() % members.len() as u64) as usize];
            let base = rng.next() | 1;
            let mut res: Vec<WalkStats> = Vec::new();
            for pol in 0..3u8 { let mut rr = Rng(base.wrapping_add(7919 * pol as u64) | 1); let mut st = WalkStats::default(); walk_once(c, start, pol, &mut rr, &mut st); res.push(st); }
            for pol in 0..3 {
                if holds_weak && res[pol].refused > 0 { refused_weak[pol] += 1; }
                if holds_weak && res[0].success == 1 && res[pol].success == 1 { both[pol].0.push(res[0].steps_all[0]); both[pol].1.push(res[pol].steps_all[0]); }
                let acc = if holds_weak { &mut on_weak[pol] } else { &mut on_not[pol] }; let one = &res[pol];
                acc.starts += one.starts; acc.success += one.success; acc.refused += one.refused; acc.ascent_fail += one.ascent_fail; acc.missing_codomain += one.missing_codomain; acc.trace_mismatch += one.trace_mismatch; acc.ascent_steps_max = acc.ascent_steps_max.max(one.ascent_steps_max); acc.steps_success.extend(one.steps_success.iter()); acc.steps_all.extend(one.steps_all.iter());
            }
        }
    }
    let fmt = |name: &str, s: &mut WalkStats| format!("| {name} | {} | {} ({:.3}) | {} | {:.1} | {:.1} | {} | {} | {} | {} |", s.starts, s.success, s.success as f64 / s.starts.max(1) as f64, s.refused, median(&mut s.steps_success), median(&mut s.steps_all), s.ascent_fail, s.ascent_steps_max, s.missing_codomain, s.trace_mismatch);
    let mut lines = vec![format!("## walk, p = {p}, q = {}, {starts_per_class} starts per class, cap 3q = {}\n", p * p, 3 * p * p),
        "| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |\n|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|".into()];
    for pol in 0..3 { lines.push(fmt(&format!("{} on classes with a weak curve", names[pol]), &mut on_weak[pol])); }
    for pol in 0..3 { lines.push(fmt(&format!("{} on classes with no weak curve", names[pol]), &mut on_not[pol])); }
    lines.push(format!("\nW-1: starts on weak-holding classes refused: N {}, A {}", refused_weak[1], refused_weak[2]));
    for pol in 1..3 { let (mut r, mut x) = both[pol].clone(); let mr = median(&mut r); let mx = median(&mut x); lines.push(format!("W-3 ({} vs R): starts where both succeed {}; median steps R {:.1}, this {:.1}, ratio {:.3}", names[pol], both[pol].0.len(), mr, mx, mx / mr.max(1e-9))); }
    let text = lines.join("\n");
    std::fs::create_dir_all(outdir).unwrap();
    std::fs::write(format!("{outdir}/walk_p{p}.md"), &text).unwrap();
    println!("{text}");
}

fn selftest_walk(c: &Census) {
    // halving formula: for random points P with halvable abscissa, some half Q satisfies 2Q = P
    let f = &c.f; let mut rng = Rng(99); let mut ok = 0; let mut tried = 0; let mut iso_ok = 0; let mut iso_tried = 0;
    for id in (0..c.nodes_j.len() as u32).step_by((c.nodes_j.len() / 40).max(1)) {
        let (cv, lam) = c.curve(id); let e = [f.zero(), f.one(), lam];
        for _ in 0..3 {
            let pt = random_point(f, &cv, &mut rng); let (x, _y) = pt.unwrap();
            if !c.halvable_x(&e, x) { continue; }
            tried += 1;
            let mut found = false;
            for hx in c.halves_x(&e, x) {
                let rhs = f.mul(hx, f.add(f.add(f.sq(hx), f.mul(cv.a, hx)), cv.b));
                if !f.is_square(rhs) { continue; }
                let hy = f.sqrt(rhs); let q2 = ec_add(f, &cv, Some((hx, hy)), Some((hx, hy)));
                if let Some((x2, _)) = q2 { if x2 == x { found = true; break; } }
            }
            if found { ok += 1; }
        }
        // 3-isogeny: codomain in the same class and a 3-isogeny back exists
        let (three, _missing) = c.three_moves(id, &mut rng);
        for nb in three { iso_tried += 1; if c.trace[nb as usize].abs() == c.class_key(id) && c.height[nb as usize] == c.height[id as usize] { iso_ok += 1; } }
    }
    println!("selftest p={}: halving {ok}/{tried}; 3-isogeny codomain same class and height {iso_ok}/{iso_tried}", c.p);
}
/// `census tabulate results/ 7 11 13 ...`: per-size tables from the frozen JSONL
/// (no recomputation): (v2(f), D_K mod 8) rows and pooled height rows.
fn tabulate(dir: &str, ps: &[u64]) {
    for &p in ps {
        let q = (p * p) as f64;
        let text = match std::fs::read_to_string(format!("{dir}/census_p{p}.jsonl")) { Ok(t) => t, Err(_) => { println!("p={p}: no file"); continue; } };
        let mut rows: Vec<(u32, i64, u64, u64, Vec<(u64, u64)>)> = Vec::new(); // v2f, dk8, n, nw, by_height
        for line in text.lines() {
            let get = |k: &str| -> i64 { let i = line.find(&format!("\"{k}\":")).unwrap() + k.len() + 3; let rest = &line[i..]; let e = rest.find(|c: char| c == ',' || c == '}').unwrap(); rest[..e].parse().unwrap() };
            let v2f = get("v2_f") as u32; let dk8 = get("D_K_mod8"); let n = get("n_full2t_j") as u64; let nw = get("n_weak_j") as u64;
            let hi = line.find("\"by_height\":[").unwrap() + 13; let hs = &line[hi..line.len() - 2];
            let mut byh = Vec::new();
            for trip in hs.split("],[") { let t = trip.trim_matches(|c| c == '[' || c == ']'); let v: Vec<u64> = t.split(',').map(|x| x.parse().unwrap()).collect(); byh.push((v[1], v[2])); }
            rows.push((v2f, dk8, n, nw, byh));
        }
        println!("\n=== p = {p}, q = {}, classes = {} ===", p * p, rows.len());
        println!("| v2(f) | D_K mod 8 | classes | weak classes | nodes | weak nodes | weak fraction × q |");
        println!("|--:|--:|--:|--:|--:|--:|--:|");
        let mut keys: Vec<(u32, i64)> = rows.iter().map(|r| (r.0, r.1)).collect(); keys.sort(); keys.dedup();
        for (v2f, dk8) in keys {
            let (mut c, mut wc, mut n, mut nw) = (0u64, 0u64, 0u64, 0u64);
            for r in &rows { if r.0 == v2f && r.1 == dk8 { c += 1; if r.3 > 0 { wc += 1; } n += r.2; nw += r.3; } }
            let label = if v2f == 64 { "D = 0".to_string() } else { v2f.to_string() };
            println!("| {label} | {dk8} | {c} | {wc} | {n} | {nw} | {:.2} |", if n > 0 { nw as f64 / n as f64 * q } else { 0.0 });
        }
        let mut byh: Vec<(u64, u64)> = Vec::new();
        for r in &rows { for (h, (a, w)) in r.4.iter().enumerate() { if byh.len() <= h { byh.resize(h + 1, (0, 0)); } byh[h].0 += a; byh[h].1 += w; } }
        println!("\n| height | nodes | weak nodes | weak fraction × q |");
        println!("|--:|--:|--:|--:|");
        for (h, (a, w)) in byh.iter().enumerate() { if *a > 0 { println!("| {h} | {a} | {w} | {:.2} |", *w as f64 / *a as f64 * q); } }
        let (d1c, d1w): (u64, u64) = rows.iter().filter(|r| r.0 == 1).fold((0, 0), |acc, r| (acc.0 + 1, acc.1 + (r.3 > 0) as u64));
        let (d2c, d2w): (u64, u64) = rows.iter().filter(|r| r.0 >= 2 && r.0 != 64).fold((0, 0), |acc, r| (acc.0 + 1, acc.1 + (r.3 > 0) as u64));
        println!("\nA1: v2(f)=1 classes {d1c}, with a weak j {d1w}.   A2: v2(f)>=2 classes {d2c}, with a weak j {d2w} ({:.3}).", if d2c > 0 { d2w as f64 / d2c as f64 } else { 0.0 });
    }
}


fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(|s| s == "tabulate").unwrap_or(false) {
        let dir = args.get(2).cloned().unwrap_or_else(|| "results".into());
        let ps: Vec<u64> = args[3..].iter().map(|s| s.parse().unwrap()).collect();
        tabulate(&dir, &ps); return;
    }
    if let Some(mode) = args.get(1).filter(|m| *m == "patterns" || *m == "walk" || *m == "walkselftest") {
        let p: u64 = args[2].parse().unwrap(); let outdir = args.get(3).cloned().unwrap_or_else(|| "results".into());
        let c = build_census(p, 20261007, true);
        match mode.as_str() {
            "patterns" => run_patterns(&c, &outdir),
            "walkselftest" => selftest_walk(&c),
            _ => { let starts: u64 = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(20); run_walk(&c, &outdir, starts, 20261008); }
        }
        return;
    }
    let p: u64 = args.get(1).expect("p").parse().unwrap();
    let outdir = args.get(2).cloned().unwrap_or_else(|| "results".into());
    let seed: u64 = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(20261007);
    assert!(p % 2 == 1 && p != 3 && p < 64);
    let t0 = std::time::Instant::now();
    let f = build_field(p);

    if args.get(2).map(|s| s == "--selftest").unwrap_or(false) {
        let mut rng = Rng(seed | 1); let mut bad = 0;
        for i in 0..20 {
            let lk = 2 + rng.next() % (p.pow(6) - 2); let lam = f.unpack(lk); if lam == f.one() { continue; }
            let c = Curve { a: f.neg(f.add(f.one(), lam)), b: lam };
            let nb = brute_order(&f, &c); let nq = curve_order(&f, &c, &mut rng);
            if nb != nq { bad += 1; eprintln!("MISMATCH lambda#{i}: brute {nb} bsgs {nq}"); }
        }
        println!("selftest p={p}: mismatches {bad}/20"); return;
    }
    let q = p * p; let q3 = p.pow(6);
    eprintln!("p={p} q={q} |F_{{p^6}}|={q3}  field built in {:.1}s", t0.elapsed().as_secs_f64());

    // 1. Legendre enumeration: one lambda per j, weak flag, 2-isogeny edges.
    let mut lambda_of: Vec<u32> = vec![u32::MAX; q3 as usize];
    let mut node_of: Vec<u32> = vec![u32::MAX; q3 as usize];
    let mut nodes_j: Vec<u32> = Vec::new();
    let one = f.one();
    for lk in 2..q3 {
        let lam = f.unpack(lk);
        if lam == one { continue; }
        let a = f.neg(f.add(one, lam)); let b = lam;
        let j = j_inv(&f, a, b); let jk = f.pack(j) as usize;
        if lambda_of[jk] == u32::MAX { lambda_of[jk] = lk as u32; node_of[jk] = nodes_j.len() as u32; nodes_j.push(jk as u32); }
    }
    let nn = nodes_j.len();
    eprintln!("full-2-torsion j-invariants: {nn}  ({:.1}s)", t0.elapsed().as_secs_f64());

    // weak flag and edges
    let mut weak = vec![false; nn];
    let mut edges: Vec<[(u32, bool); 3]> = vec![[(u32::MAX, false); 3]; nn]; // (codomain j index, keeps full 2-torsion)
    let mut n_weak = 0usize;
    for (id, &jk) in nodes_j.iter().enumerate() {
        let lam = f.unpack(lambda_of[jk as usize] as u64);
        // cross ratios lam, 1-lam, lam/(lam-1)
        let c1 = lam; let c2 = f.sub(one, lam); let c3 = f.mul(lam, f.inv(f.sub(lam, one)));
        let w = [c1, c2, c3].iter().any(|&c| f.norm2(c) == F2::one());
        weak[id] = w; if w { n_weak += 1; }
        let e = [f.zero(), one, lam];
        for i in 0..3 {
            let (ei, ej, ek) = (e[i], e[(i + 1) % 3], e[(i + 2) % 3]);
            let a = f.sub(f.sub(f.scal(2, ei), ej), ek);
            let b = f.mul(f.sub(ei, ej), f.sub(ei, ek));
            let a2 = f.neg(f.scal(2, a)); let b2 = f.sub(f.sq(a), f.scal(4, b));
            let jj = j_inv(&f, a2, b2);
            edges[id][i] = (f.pack(jj) as u32, f.is_square(b));
        }
    }
    eprintln!("weak j: {n_weak}  fraction {:.5}  (3/q = {:.5})  ({:.1}s)", n_weak as f64 / nn as f64, 3.0 / q as f64, t0.elapsed().as_secs_f64());

    // 2. components under 2-isogenies among full-2-torsion nodes (union-find)
    let mut parent: Vec<u32> = (0..nn as u32).collect();
    fn find(par: &mut Vec<u32>, mut x: u32) -> u32 { while par[x as usize] != x { let g = par[par[x as usize] as usize]; par[x as usize] = g; x = g; } x }
    let mut floor_adjacent = vec![false; nn];
    let mut unresolved_edges = 0usize;
    for id in 0..nn {
        for i in 0..3 {
            let (jk, keeps) = edges[id][i];
            if keeps {
                let nb = node_of[jk as usize];
                if nb == u32::MAX { unresolved_edges += 1; continue; } // should not happen: a full-2-torsion codomain is a Legendre j
                let (ra, rb) = (find(&mut parent, id as u32), find(&mut parent, nb));
                if ra != rb { parent[ra as usize] = rb; }
            } else { floor_adjacent[id] = true; }
        }
    }
    assert!(unresolved_edges == 0, "codomain with full 2-torsion missing from the Legendre set: {unresolved_edges}");
    let mut comp_members: HashMap<u32, Vec<u32>> = HashMap::new();
    for id in 0..nn as u32 { let r = find(&mut parent, id); comp_members.entry(r).or_default().push(id); }
    let ncomp = comp_members.len();
    eprintln!("2-isogeny components: {ncomp}  ({:.1}s)", t0.elapsed().as_secs_f64());

    // 3. traces per component by BSGS (parallel over components), with a second check
    let comps: Vec<(u32, Vec<u32>)> = comp_members.into_iter().collect();
    let nthreads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(14);
    let chunk = (comps.len() + nthreads - 1) / nthreads;
    let mut trace_of_comp: HashMap<u32, i64> = HashMap::new();
    let mut check_fail = 0usize;
    std::thread::scope(|sc| {
        let mut handles = Vec::new();
        for (ti, part) in comps.chunks(chunk.max(1)).enumerate() {
            let fref = &f; let lam_ref = &lambda_of; let nodes_ref = &nodes_j;
            handles.push(sc.spawn(move || {
                let mut out = Vec::with_capacity(part.len()); let mut fails = 0usize;
                let mut rng = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ ((ti as u64 + 1) << 32) | 1);
                for (root, members) in part {
                    let m0 = members[0];
                    let lam = fref.unpack(lam_ref[nodes_ref[m0 as usize] as usize] as u64);
                    let c = Curve { a: fref.neg(fref.add(fref.one(), lam)), b: lam };
                    let n = curve_order(fref, &c, &mut rng);
                    let t = (q3 as i64 + 1) - n as i64;
                    // second check on another member (or the same curve with fresh points)
                    let m1 = members[members.len() / 2];
                    let lam1 = fref.unpack(lam_ref[nodes_ref[m1 as usize] as usize] as u64);
                    let c1 = Curve { a: fref.neg(fref.add(fref.one(), lam1)), b: lam1 };
                    let n1 = curve_order(fref, &c1, &mut rng);
                    if n1 != n && n1 != 2 * (q3 + 1) - n { fails += 1; } // twist pairs share a j-component
                    out.push((*root, t));
                }
                (out, fails)
            }));
        }
        for h in handles { let (out, fails) = h.join().unwrap(); check_fail += fails; for (r, t) in out { trace_of_comp.insert(r, t); } }
    });
    eprintln!("traces done, second-BSGS mismatches: {check_fail}  ({:.1}s)", t0.elapsed().as_secs_f64());
    let mut trace: Vec<i64> = vec![0; nn];
    for id in 0..nn as u32 { let r = find(&mut parent, id); trace[id as usize] = trace_of_comp[&r]; }

    // 4. heights: BFS from floor-adjacent nodes over keep-edges (undirected via symmetric adjacency)
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); nn];
    for id in 0..nn { for i in 0..3 { let (jk, keeps) = edges[id][i]; if keeps { let nb = node_of[jk as usize]; adj[id].push(nb); adj[nb as usize].push(id as u32); } } }
    let mut height: Vec<u32> = vec![u32::MAX; nn];
    let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
    for id in 0..nn { if floor_adjacent[id] { height[id] = 1; queue.push_back(id as u32); } }
    while let Some(v) = queue.pop_front() { let h = height[v as usize]; for &w in &adj[v as usize] { if height[w as usize] == u32::MAX { height[w as usize] = h + 1; queue.push_back(w); } } }
    let no_floor = height.iter().filter(|&&h| h == u32::MAX).count();

    // 5. per-class aggregation by |t|
    #[derive(Default)]
    struct Cls { n: u64, nw: u64, by_h: Vec<(u64, u64)>, comps: u64 }
    let mut classes: HashMap<i64, Cls> = HashMap::new();
    for id in 0..nn {
        let key = trace[id].abs(); let e = classes.entry(key).or_default();
        e.n += 1; if weak[id] { e.nw += 1; }
        let h = if height[id] == u32::MAX { 0 } else { height[id] as usize };
        if e.by_h.len() <= h { e.by_h.resize(h + 1, (0, 0)); }
        e.by_h[h].0 += 1; if weak[id] { e.by_h[h].1 += 1; }
    }
    for (r, members) in &comps { let key = trace_of_comp[r].abs(); classes.get_mut(&key).unwrap().comps += 1; let _ = members; }

    std::fs::create_dir_all(&outdir).unwrap();
    let mut fo = std::fs::File::create(format!("{outdir}/census_p{p}.jsonl")).unwrap();
    let mut keys: Vec<i64> = classes.keys().cloned().collect(); keys.sort();
    let mut n_in_weak_class = 0u64; let mut n_classes_weak = 0usize;
    let mut triples: HashMap<(u32, i64, i64), (usize, usize)> = HashMap::new();
    let mut w4_h1_low = 0usize; let mut w4_top_high = 0usize; let mut w4_eval = 0usize;
    let mut w5_height_ok = 0usize; let mut w5_height_bad = 0usize;
    for &t in &keys {
        let cl = &classes[&t];
        let d: i64 = t * t - 4 * q3 as i64; // negative
        let fac = factor_u((-d) as u64);
        let mut f0: u64 = 1; let mut sqf: i64 = -1;
        for (pr, e) in &fac { for _ in 0..(e / 2) { f0 *= pr; } if e % 2 == 1 { sqf *= *pr as i64; } }
        let (dk, cond) = if sqf.rem_euclid(4) == 1 { (sqf, f0) } else { (4 * sqf, f0 / 2) };
        assert!(cond as i64 * cond as i64 * dk == d, "D decomposition");
        let v2f = cond.trailing_zeros(); let dk8 = dk.rem_euclid(8); let t16 = t.rem_euclid(16);
        let is_weak_class = cl.nw > 0;
        if is_weak_class { n_in_weak_class += cl.n; n_classes_weak += 1; }
        let ent = triples.entry((v2f, dk8, t16)).or_default(); if is_weak_class { ent.0 += 1 } else { ent.1 += 1 }
        let maxh = cl.by_h.len().saturating_sub(1);
        if maxh > 0 { if maxh as u32 == v2f { w5_height_ok += 1 } else { w5_height_bad += 1 } }
        if is_weak_class && cl.n >= 20 && cl.by_h.len() > 2 {
            w4_eval += 1;
            let pooled = cl.nw as f64 / cl.n as f64;
            let (a1, w1) = cl.by_h[1]; let (at, wt) = cl.by_h[maxh];
            if a1 > 0 && (w1 as f64 / a1 as f64) < 0.5 * pooled { w4_h1_low += 1; }
            if at > 0 && (wt as f64 / at as f64) > 1.5 * pooled { w4_top_high += 1; }
        }
        let byh: Vec<String> = cl.by_h.iter().enumerate().map(|(h, (a, w))| format!("[{h},{a},{w}]")).collect();
        writeln!(fo, "{{\"p\":{p},\"t_abs\":{t},\"n_full2t_j\":{},\"n_weak_j\":{},\"components\":{},\"D\":{d},\"D_K\":{dk},\"f\":{cond},\"v2_f\":{v2f},\"D_K_mod8\":{dk8},\"t_mod16\":{t16},\"v2_N_plus\":{},\"v2_N_minus\":{},\"by_height\":[{}]}}",
            cl.n, cl.nw, cl.comps,
            (q3 as i64 + 1 - t).trailing_zeros(), (q3 as i64 + 1 + t).trailing_zeros(), byh.join(",")).unwrap();
    }
    // W3: conflicts
    let mut conflicts: Vec<String> = Vec::new(); let mut n_triples = 0usize;
    for ((v2f, dk8, t16), (nw, nnw)) in &triples { n_triples += 1; if *nw > 0 && *nnw > 0 { conflicts.push(format!("(v2f={v2f},DK8={dk8},t16={t16}):weak={nw},not={nnw}")); } }
    conflicts.sort();
    let frac_in_weak = n_in_weak_class as f64 / nn as f64;
    let summary = format!("{{\"p\":{p},\"q\":{q},\"nodes\":{nn},\"weak_nodes\":{n_weak},\"weak_fraction\":{:.6},\"three_over_q\":{:.6},\"classes\":{},\"weak_classes\":{},\"fraction_nodes_in_weak_class\":{:.4},\"components\":{ncomp},\"second_bsgs_mismatch\":{check_fail},\"nodes_without_floor\":{no_floor},\"triples\":{n_triples},\"triples_with_conflict\":{},\"conflicts\":{:?},\"w4_classes_evaluated\":{w4_eval},\"w4_h1_below_half\":{w4_h1_low},\"w4_top_above_1p5\":{w4_top_high},\"w5_maxheight_eq_v2f\":{w5_height_ok},\"w5_maxheight_ne_v2f\":{w5_height_bad},\"seconds\":{:.1}}}",
        n_weak as f64 / nn as f64, 3.0 / q as f64, keys.len(), n_classes_weak, frac_in_weak, conflicts.len(), conflicts, t0.elapsed().as_secs_f64());
    std::fs::write(format!("{outdir}/summary_p{p}.json"), &summary).unwrap();
    println!("{summary}");
}
