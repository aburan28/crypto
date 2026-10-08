//! Exhaustive census of full-2-torsion curves over F_{p^6}: which isogeny
//! classes hold a Joux–Vitse weak curve, and where in the 2-volcano the weak
//! curves sit.  See PROTOCOL.md (frozen before this file was written).
//!
//! Standalone: `rustc -O census.rs -o census && ./census 7 results/`.
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
