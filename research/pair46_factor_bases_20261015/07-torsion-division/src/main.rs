// 07-torsion-division: torsion / division-polynomial factor bases.
// Standalone, no dependencies. Usage: torsion_division > results.json
use std::collections::HashMap;
use std::fmt::Write as _;
use std::time::Instant;

// ---------------- modular arithmetic ----------------
#[inline]
fn mulm(a: u64, b: u64, p: u64) -> u64 { ((a as u128 * b as u128) % p as u128) as u64 }
#[inline]
fn addm(a: u64, b: u64, p: u64) -> u64 { let s = a as u128 + b as u128; (s % p as u128) as u64 }
#[inline]
fn subm(a: u64, b: u64, p: u64) -> u64 { if a >= b { a - b } else { p - (b - a) } }
fn powm(mut a: u64, mut e: u128, p: u64) -> u64 {
    let mut r = 1u64 % p; a %= p;
    while e > 0 { if e & 1 == 1 { r = mulm(r, a, p); } a = mulm(a, a, p); e >>= 1; }
    r
}
fn invm(a: u64, p: u64) -> u64 { powm(a, (p - 2) as u128, p) }
fn legendre(a: u64, p: u64) -> i32 {
    let a = a % p; if a == 0 { return 0; }
    if powm(a, ((p - 1) / 2) as u128, p) == 1 { 1 } else { -1 }
}
fn sqrtm(a: u64, p: u64) -> Option<u64> {
    let a = a % p; if a == 0 { return Some(0); }
    if legendre(a, p) != 1 { return None; }
    if p % 4 == 3 { return Some(powm(a, ((p + 1) / 4) as u128, p)); }
    // Tonelli-Shanks
    let mut q = p - 1; let mut s = 0; while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while legendre(z, p) != -1 { z += 1; }
    let mut m = s; let mut c = powm(z, q as u128, p); let mut t = powm(a, q as u128, p);
    let mut r = powm(a, ((q + 1) / 2) as u128, p);
    while t != 1 {
        let mut i = 0; let mut tt = t; while tt != 1 { tt = mulm(tt, tt, p); i += 1; }
        let mut b = c; for _ in 0..(m - i - 1) { b = mulm(b, b, p); }
        m = i; c = mulm(b, b, p); t = mulm(t, c, p); r = mulm(r, b, p);
    }
    Some(r)
}
fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % q == 0 { return n == q; } }
    let mut d = n - 1; let mut s = 0; while d % 2 == 0 { d /= 2; s += 1; }
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d as u128, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulm(x, x, n); if x == n - 1 { continue 'w; } }
        return false;
    }
    true
}

// ---------------- rng ----------------
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; self.0 }
    fn below(&mut self, n: u64) -> u64 { ((self.next() as u128 * n as u128) >> 64) as u64 }
    fn new(seed: u64) -> Self { let mut r = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0xD1B54A32D192ED03); for _ in 0..8 { r.next(); } r }
}

// ---------------- polynomials over F_p (coeffs low->high) ----------------
type Poly = Vec<u64>;
fn trim(mut a: Poly) -> Poly { while a.len() > 0 && *a.last().unwrap() == 0 { a.pop(); } a }
fn padd(a: &Poly, b: &Poly, p: u64) -> Poly {
    let n = a.len().max(b.len()); let mut r = vec![0; n];
    for i in 0..n { r[i] = addm(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0), p); }
    trim(r)
}
fn psub(a: &Poly, b: &Poly, p: u64) -> Poly {
    let n = a.len().max(b.len()); let mut r = vec![0; n];
    for i in 0..n { r[i] = subm(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0), p); }
    trim(r)
}
fn pmul(a: &Poly, b: &Poly, p: u64) -> Poly {
    if a.is_empty() || b.is_empty() { return vec![]; }
    let mut acc = vec![0u128; a.len() + b.len() - 1];
    let pp = p as u128; let lim = u128::MAX - pp * pp;
    for (i, &x) in a.iter().enumerate() { if x == 0 { continue; }
        for (j, &y) in b.iter().enumerate() { let t = &mut acc[i + j]; *t += x as u128 * y as u128; if *t > lim { *t %= pp; } } }
    trim(acc.into_iter().map(|v| (v % pp) as u64).collect())
}
fn pscale(a: &Poly, c: u64, p: u64) -> Poly { trim(a.iter().map(|&x| mulm(x, c, p)).collect()) }
fn pdivrem(a: &Poly, b: &Poly, p: u64) -> (Poly, Poly) {
    let a = trim(a.clone()); let b = trim(b.clone());
    if a.len() < b.len() { return (vec![], a); }
    let inv = invm(*b.last().unwrap(), p); let mut r = a.clone(); let db = b.len() - 1;
    let mut q = vec![0; a.len() - db];
    for i in (0..q.len()).rev() {
        let c = mulm(r[i + db], inv, p); q[i] = c;
        if c != 0 { for j in 0..=db { r[i + j] = subm(r[i + j], mulm(c, b[j], p), p); } }
    }
    (trim(q), trim(r))
}
fn pgcd(a: &Poly, b: &Poly, p: u64) -> Poly {
    let mut a = trim(a.clone()); let mut b = trim(b.clone());
    while !b.is_empty() { let (_, r) = pdivrem(&a, &b, p); a = b; b = r; }
    if a.is_empty() { return a; } let inv = invm(*a.last().unwrap(), p); pscale(&a, inv, p)
}
fn peval(a: &Poly, x: u64, p: u64) -> u64 { let mut r = 0; for &c in a.iter().rev() { r = addm(mulm(r, x, p), c, p); } r }
fn deg(a: &Poly) -> i64 { a.len() as i64 - 1 }
// x^e mod g
fn xpow_mod(e: u128, g: &Poly, p: u64) -> Poly {
    let mut r: Poly = vec![1]; let mut base: Poly = pdivrem(&vec![0, 1], g, p).1;
    let mut e = e;
    while e > 0 { if e & 1 == 1 { r = pdivrem(&pmul(&r, &base, p), g, p).1; } base = pdivrem(&pmul(&base, &base, p), g, p).1; e >>= 1; }
    r
}
fn count_fp_roots(g: &Poly, p: u64) -> (usize, Poly) {
    let g = trim(g.clone()); if deg(&g) < 1 { return (0, vec![1]); }
    let xp = xpow_mod(p as u128, &g, p);
    let h = pgcd(&g, &psub(&xp, &vec![0, 1], p), p);
    (deg(&h).max(0) as usize, h)
}

// ---------------- division polynomials ----------------
// psi_n = g_n (n odd) or y*g_n (n even), y^2 = f.
struct Div { p: u64, f: Poly, g: Vec<Poly> }
impl Div {
    fn new(a: u64, b: u64, p: u64, nmax: usize) -> Div {
        let f: Poly = trim(vec![b % p, a % p, 0, 1]);
        let mut g: Vec<Poly> = vec![vec![]; nmax + 3];
        g[0] = vec![]; g[1] = vec![1]; g[2] = vec![2 % p];
        // psi_3 = 3x^4 + 6a x^2 + 12 b x - a^2
        g[3] = trim(vec![subm(0, mulm(a, a, p), p), mulm(12, b, p), mulm(6, a, p), 0, 3]);
        // psi_4 = 4y(x^6 + 5a x^4 + 20 b x^3 - 5a^2 x^2 - 4ab x - 8b^2 - a^3)
        if nmax + 2 >= 4 {
            let c0 = subm(0, addm(mulm(8, mulm(b, b, p), p), mulm(a, mulm(a, a, p), p), p), p);
            g[4] = pscale(&trim(vec![c0, subm(0, mulm(4, mulm(a, b, p), p), p), subm(0, mulm(5, mulm(a, a, p), p), p), mulm(20, b, p), mulm(5, a, p), 0, 1]), 4, p);
        }
        let mut d = Div { p, f, g };
        for n in 5..=nmax + 2 { let v = d.compute(n); d.g[n] = v; }
        d
    }
    // returns psi_n as (poly, y-exponent)
    fn psi(&self, n: usize) -> (Poly, u32) { (self.g[n].clone(), if n % 2 == 0 { 1 } else { 0 }) }
    fn mul(&self, a: &(Poly, u32), b: &(Poly, u32)) -> (Poly, u32) {
        let mut r = pmul(&a.0, &b.0, self.p); let e = a.1 + b.1;
        for _ in 0..(e / 2) { r = pmul(&r, &self.f, self.p); }
        (r, e % 2)
    }
    fn cube(&self, a: &(Poly, u32)) -> (Poly, u32) { let s = self.mul(a, a); self.mul(&s, a) }
    fn compute(&self, n: usize) -> Poly {
        let p = self.p;
        if n % 2 == 1 {
            let k = (n - 1) / 2; // psi_{2k+1} = psi_{k+2} psi_k^3 - psi_{k-1} psi_{k+1}^3
            let t1 = self.mul(&self.psi(k + 2), &self.cube(&self.psi(k)));
            let t2 = self.mul(&self.psi(k - 1), &self.cube(&self.psi(k + 1)));
            assert_eq!(t1.1, t2.1); assert_eq!(t1.1, 0);
            psub(&t1.0, &t2.0, p)
        } else {
            let k = n / 2; // 2y psi_{2k} = psi_k (psi_{k+2} psi_{k-1}^2 - psi_{k-2} psi_{k+1}^2)
            let a1 = self.mul(&self.psi(k + 2), &self.mul(&self.psi(k - 1), &self.psi(k - 1)));
            let a2 = self.mul(&self.psi(k - 2), &self.mul(&self.psi(k + 1), &self.psi(k + 1)));
            assert_eq!(a1.1, a2.1);
            let inner = (psub(&a1.0, &a2.0, p), a1.1);
            let x = self.mul(&self.psi(k), &inner); // = 2 y^2 g = 2 f g, y-exp 0
            assert_eq!(x.1, 0);
            let denom = pscale(&self.f, 2 % p, p);
            let (q, r) = pdivrem(&x.0, &denom, p); assert!(r.is_empty(), "psi_{} division not exact", n);
            q
        }
    }
    // psi_m^2 as a polynomial in x
    fn psi_sq(&self, m: usize) -> Poly { let s = self.mul(&self.psi(m), &self.psi(m)); assert_eq!(s.1, 0); s.0 }
    // phi_m = x psi_m^2 - psi_{m+1} psi_{m-1}
    fn phi(&self, m: usize) -> Poly {
        let t = self.mul(&self.psi(m + 1), &self.psi(m - 1)); assert_eq!(t.1, 0);
        psub(&pmul(&vec![0, 1], &self.psi_sq(m), self.p), &t.0, self.p)
    }
    // psi_m as a polynomial in x for odd m
    fn psi_odd(&self, m: usize) -> Poly { assert!(m % 2 == 1); self.g[m].clone() }
}

// ---------------- curve arithmetic ----------------
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Pt { O, A(u64, u64) }
#[derive(Clone, Copy)]
struct Curve { a: u64, b: u64, p: u64 }
impl Curve {
    fn rhs(&self, x: u64) -> u64 { let p = self.p; addm(addm(mulm(mulm(x, x, p), x, p), mulm(self.a, x, p), p), self.b, p) }
    fn neg(&self, q: Pt) -> Pt { match q { Pt::O => Pt::O, Pt::A(x, y) => Pt::A(x, subm(0, y, self.p)) } }
    fn add(&self, u: Pt, v: Pt) -> Pt {
        let p = self.p;
        match (u, v) {
            (Pt::O, _) => v, (_, Pt::O) => u,
            (Pt::A(x1, y1), Pt::A(x2, y2)) => {
                let l = if x1 == x2 {
                    if addm(y1, y2, p) == 0 { return Pt::O; }
                    mulm(addm(mulm(3, mulm(x1, x1, p), p), self.a, p), invm(mulm(2, y1, p), p), p)
                } else { mulm(subm(y2, y1, p), invm(subm(x2, x1, p), p), p) };
                let x3 = subm(subm(mulm(l, l, p), x1, p), x2, p);
                let y3 = subm(mulm(l, subm(x1, x3, p), p), y1, p);
                Pt::A(x3, y3)
            }
        }
    }
    fn mul(&self, k: u128, q: Pt) -> Pt { let mut r = Pt::O; let mut b = q; let mut k = k; while k > 0 { if k & 1 == 1 { r = self.add(r, b); } b = self.add(b, b); k >>= 1; } r }
    fn count(&self) -> u64 { let p = self.p; let mut s: i64 = 0; for x in 0..p { s += legendre(self.rhs(x), p) as i64; } (p as i64 + 1 + s) as u64 }
    fn random_point(&self, r: &mut Rng) -> Pt { loop { let x = r.below(self.p); if let Some(y) = sqrtm(self.rhs(x), self.p) { if y != 0 || self.rhs(x) == 0 { return Pt::A(x, y); } } } }
    fn singular(&self) -> bool { let p = self.p; addm(mulm(4, mulm(self.a, mulm(self.a, self.a, p), p), p), mulm(27, mulm(self.b, self.b, p), p), p) == 0 }
}

// ---------------- toy setup ----------------
struct Toy { c: Curve, n: u64, gen: Pt, logx: HashMap<u64, u64>, pts: Vec<Pt> }
fn build_toy(c: Curve, n: u64, gen: Pt) -> Toy {
    let mut logx = HashMap::new(); let mut pts = Vec::with_capacity(n as usize);
    let mut q = Pt::O;
    for k in 0..n { pts.push(q); if let Pt::A(x, _) = q { logx.entry(x).or_insert(k); } q = c.add(q, gen); }
    assert_eq!(q, Pt::O, "generator order mismatch");
    Toy { c, n, gen, logx, pts }
}
fn next_prime(mut x: u64) -> u64 { while !is_prime(x) { x += 1; } x }

fn find_prime_order_curve(bits: u32, seed: u64) -> (Curve, u64) {
    let mut r = Rng::new(1000 + bits as u64 * 100 + seed);
    let p = next_prime((1u64 << (bits - 1)) + r.below(1u64 << (bits - 1)));
    loop {
        let c = Curve { a: r.below(p), b: r.below(p), p };
        if c.singular() { continue; }
        let n = c.count(); if is_prime(n) && n != p { return (c, n); }
    }
}
fn find_lq_curve(bits: u32, l: u64, seed: u64) -> (Curve, u64, u64) {
    let mut r = Rng::new(5000 + bits as u64 * 100 + l * 10 + seed);
    let p = next_prime((1u64 << (bits - 1)) + r.below(1u64 << (bits - 1)));
    loop {
        let c = Curve { a: r.below(p), b: r.below(p), p };
        if c.singular() { continue; }
        let n = c.count(); if n % l == 0 && is_prime(n / l) && n / l != l { return (c, n, n / l); }
    }
}

fn coverage_bits(logs: &[u64], n: u64) -> Vec<bool> {
    let mut t = vec![false; n as usize];
    for &a in logs { for &b in logs {
        for (s1, s2) in [(1i64, 1i64), (1, -1), (-1, 1), (-1, -1)] {
            let v = ((s1 * a as i64 + s2 * b as i64).rem_euclid(n as i64)) as usize; t[v] = true;
        }
    } }
    t
}
fn coverage(logs: &[u64], n: u64) -> f64 { let t = coverage_bits(logs, n); t.iter().filter(|&&b| b).count() as f64 / n as f64 }

fn bootstrap_mean(v: &[f64], seed: u64) -> (f64, f64, f64) {
    let n = v.len(); let mean = v.iter().sum::<f64>() / n as f64; let mut r = Rng::new(seed);
    let mut ms: Vec<f64> = (0..2000).map(|_| { let mut s = 0.0; for _ in 0..n { s += v[r.below(n as u64) as usize]; } s / n as f64 }).collect();
    ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (mean, ms[50], ms[1949])
}

fn jf(x: f64) -> String { if x.is_finite() { format!("{:.6}", x) } else { "null".into() } }

fn main() {
    let t0 = Instant::now();
    let mut out = String::new();
    out.push_str("{\n\"item\": \"07-torsion-division\",\n\"protocol_sha256\": \"97f5763067d1a35b1d699b4c5972269dc7dc5408f10d6a03d5232614768f62b6\",\n");

    // ================= Pair 46 exact checks =================
    let p46: u64 = 5947698030651343; let n46: u64 = 5947697996161643; let v: u64 = 34489703;
    let curves46 = [("crater", 5757116383569871u64, 3560173328417442u64), ("floor", 2087370704375809u64, 1402887937802154u64)];
    out.push_str("\"pair46\": [\n");
    for (ci, &(name, a, b)) in curves46.iter().enumerate() {
        let c = Curve { a, b, p: p46 };
        let mut r = Rng::new(46 + ci as u64);
        let mut order_ok = true; for _ in 0..5 { let q = c.random_point(&mut r); if c.mul(n46 as u128, q) != Pt::O { order_ok = false; } }
        let div = Div::new(a, b, p46, 13);
        let mut rows = String::new();
        // m = 2: roots of f
        let (r2, _) = count_fp_roots(&div.f, p46);
        let _ = write!(rows, "{{\"m\": 2, \"deg\": 3, \"fp_roots\": {}, \"on_E\": 0, \"on_twist\": 0}}", r2);
        for m in [3usize, 5, 7, 11, 13] {
            let ts = Instant::now();
            let g = div.psi_odd(m);
            let (nr, h) = count_fp_roots(&g, p46);
            // find the actual roots when few (deg h <= 2) and classify E vs twist
            let mut on_e = 0; let mut on_t = 0; let mut roots = vec![];
            if nr == 1 { let x = subm(0, mulm(h[0], invm(h[1], p46), p46), p46); roots.push(x); }
            if nr == 2 {
                let (c0, c1) = (h[0], h[1]); // x^2 + c1 x + c0
                let disc = subm(mulm(c1, c1, p46), mulm(4, c0, p46), p46);
                if let Some(s) = sqrtm(disc, p46) { let i2 = invm(2, p46);
                    roots.push(mulm(subm(s, c1, p46), i2, p46)); roots.push(mulm(subm(subm(0, s, p46), c1, p46), i2, p46)); }
            }
            for &x in &roots { assert_eq!(peval(&g, x, p46), 0); match legendre(c.rhs(x), p46) { 1 => on_e += 1, -1 => on_t += 1, _ => {} } }
            let _ = write!(rows, ",\n    {{\"m\": {}, \"deg\": {}, \"fp_roots\": {}, \"on_E\": {}, \"on_twist\": {}, \"roots\": {:?}, \"secs\": {:.3}}}", m, deg(&g), nr, on_e, on_t, roots, ts.elapsed().as_secs_f64());
        }
        // twist exponent check
        let mut d = 2u64; while legendre(d, p46) != -1 { d += 1; }
        let ct = Curve { a: mulm(a, mulm(d, d, p46), p46), b: mulm(b, mulm(d, mulm(d, d, p46), p46), p46), p: p46 };
        let nt: u128 = 5 * v as u128 * v as u128;
        assert_eq!(nt, 2 * (p46 as u128 + 1) - n46 as u128);
        let mut r = Rng::new(4600 + ci as u64);
        let mut ord_ok = 0; let mut exp5v = 0; let mut exp_v_in_5 = 0;
        for _ in 0..20 { let q = ct.random_point(&mut r);
            if ct.mul(nt, q) == Pt::O { ord_ok += 1; }
            if ct.mul(5 * v as u128, q) == Pt::O { exp5v += 1; }
            if ct.mul(v as u128, q) == Pt::O { exp_v_in_5 += 1; } }
        let _ = write!(out, "  {{\"curve\": \"{}\", \"a\": {}, \"b\": {}, \"N_check_5pts\": {}, \"psi_roots\": [\n    {}\n  ],\n   \"twist\": {{\"order\": {}, \"nonresidue_d\": {}, \"points\": 20, \"killed_by_order\": {}, \"killed_by_5v\": {}, \"killed_by_v\": {}}}}}{}\n",
            name, a, b, order_ok, rows, nt, d, ord_ok, exp5v, exp_v_in_5, if ci == 0 { "," } else { "" });
    }
    out.push_str("],\n");
    eprintln!("pair46 done {:.1}s", t0.elapsed().as_secs_f64());

    // ================= Toy prime-order curves =================
    let ms_pre = [2usize, 3, 5, 7]; let ms_lvl = [3usize, 5, 7];
    out.push_str("\"toys\": [\n");
    let mut summary: HashMap<(u32, String), Vec<f64>> = HashMap::new();
    let mut first = true;
    for bits in [12u32, 14, 16] {
        for seed in 0..5u64 {
            let ts = Instant::now();
            let (c, n) = find_prime_order_curve(bits, seed);
            let p = c.p;
            let mut r = Rng::new(77 + seed * 31 + bits as u64);
            let gen = c.random_point(&mut r);
            let toy = build_toy(c, n, gen);
            let bsz = ((n as f64).sqrt() / 2.0).ceil() as usize;
            // SX
            let mut sx: Vec<u64> = vec![]; let mut x = 0; while sx.len() < bsz { if toy.logx.contains_key(&x) { sx.push(x); } x += 1; }
            let logs_of = |xs: &[u64]| -> Vec<u64> { xs.iter().map(|x| toy.logx[x]).collect() };
            let c_sx = coverage(&logs_of(&sx), n);
            // RND
            let mut allx: Vec<u64> = toy.logx.keys().cloned().collect(); allx.sort();
            let mut rr = Rng::new(900 + seed + bits as u64 * 7);
            for i in (1..allx.len()).rev() { let j = rr.below(i as u64 + 1) as usize; allx.swap(i, j); }
            let rnd: Vec<u64> = allx[..bsz].to_vec();
            let c_rnd = coverage(&logs_of(&rnd), n);
            // KL
            let kl: Vec<u64> = (1..=bsz as u64).collect();
            let c_kl = coverage(&kl, n);
            // division polynomials
            let div = Div::new(c.a, c.b, p, 8);
            let sxset: std::collections::HashSet<u64> = sx.iter().cloned().collect();
            let mut pre_json = String::new();
            for (mi, &m) in ms_pre.iter().enumerate() {
                let phi = div.phi(m); let ps2 = div.psi_sq(m);
                let mut base = vec![]; let mut off_e_hits = 0; let mut log_mismatch = 0;
                for x in 0..p {
                    let den = peval(&ps2, x, p); if den == 0 { continue; }
                    let xm = mulm(peval(&phi, x, p), invm(den, p), p);
                    if sxset.contains(&xm) {
                        if toy.logx.contains_key(&x) { base.push(x);
                            // cross-check with log table: m*log(x) = +-log(xm)
                            let k = toy.logx[&x]; let km = (k as u128 * m as u128 % n as u128) as u64; let t = toy.logx[&xm];
                            if km != t && km != n - t { log_mismatch += 1; }
                        } else { off_e_hits += 1; }
                    }
                }
                let lg = logs_of(&base); let c_pre = coverage(&lg, n);
                // exact transport check: T(pre) = m^{-1} T(sx)
                let tsx = coverage_bits(&logs_of(&sx), n); let tpre = coverage_bits(&lg, n);
                let mut transport_ok = true; for k in 0..n as usize { if tpre[k] != tsx[(k as u128 * m as u128 % n as u128) as usize] { transport_ok = false; break; } }
                // ideal build + univariate solve, 12-bit only
                let mut ideal_json = String::from("null");
                if bits == 12 {
                    let ti = Instant::now();
                    let mut g: Poly = vec![1];
                    for &s in &sx { g = pmul(&g, &psub(&phi, &pscale(&ps2, s, p), p), p); }
                    let gdeg = deg(&g);
                    let mut roots_e = vec![]; let mut roots_off = 0;
                    for x in 0..p { if peval(&g, x, p) == 0 { if toy.logx.contains_key(&x) { roots_e.push(x); } else { roots_off += 1; } } }
                    let mut bs = base.clone(); bs.sort(); let mut re = roots_e.clone(); re.sort();
                    let roots_match = bs == re;
                    // M4: 20 targets
                    let mut agree = 0; let mut tot_sol = 0u64; let mut rt = Rng::new(31337 + seed);
                    for _ in 0..20 {
                        let k = 1 + rt.below(n - 1); let xr = match toy.pts[k as usize] { Pt::A(x, _) => x, _ => unreachable!() };
                        // univariate route
                        let mut sols = 0u64;
                        for &al in &roots_e {
                            let aa = mulm(subm(al, xr, p), subm(al, xr, p), p);
                            let bb = subm(0, mulm(2, addm(mulm(addm(al, xr, p), addm(mulm(al, xr, p), c.a, p), p), mulm(2, c.b, p), p), p), p);
                            let t1 = subm(mulm(al, xr, p), c.a, p);
                            let cc = subm(mulm(t1, t1, p), mulm(mulm(4, c.b, p), addm(al, xr, p), p), p);
                            let mut cand = vec![];
                            if aa == 0 { if bb != 0 { cand.push(mulm(subm(0, cc, p), invm(bb, p), p)); } }
                            else { let disc = subm(mulm(bb, bb, p), mulm(4, mulm(aa, cc, p), p), p);
                                if let Some(sq) = sqrtm(disc, p) { let i2a = invm(mulm(2, aa, p), p);
                                    let r1 = mulm(subm(sq, bb, p), i2a, p); let r2 = mulm(subm(subm(0, sq, p), bb, p), i2a, p);
                                    cand.push(r1); if r2 != r1 { cand.push(r2); } } }
                            for x2 in cand { if peval(&g, x2, p) == 0 && toy.logx.contains_key(&x2) { sols += 1; } }
                        }
                        // log-table count: ordered pairs (i,j) with +-ki +-kj = +-k
                        let mut lt = 0u64;
                        for &a1 in &lg { for &a2 in &lg { let mut hit = false;
                            for (s1, s2) in [(1i64, 1i64), (1, -1), (-1, 1), (-1, -1)] { let vv = ((s1 * a1 as i64 + s2 * a2 as i64).rem_euclid(n as i64)) as u64; if vv == k || vv == n - k { hit = true; } }
                            if hit { lt += 1; } } }
                        if lt == sols { agree += 1; } tot_sol += sols;
                    }
                    ideal_json = format!("{{\"g_degree\": {}, \"g_roots_on_E\": {}, \"g_roots_off_E\": {}, \"roots_equal_base\": {}, \"bezout_bound\": {}, \"targets\": 20, \"solution_count_agrees\": {}, \"total_solutions\": {}, \"secs\": {:.3}}}",
                        gdeg, roots_e.len(), roots_off, roots_match, 2 * gdeg * gdeg, agree, tot_sol, ti.elapsed().as_secs_f64());
                }
                summary.entry((bits, format!("PRE_{}_minus_SX", m))).or_default().push(c_pre - c_sx);
                let _ = write!(pre_json, "{}{{\"m\": {}, \"size\": {}, \"membership_degree\": {}, \"rho\": {}, \"coverage\": {}, \"coverage_minus_sx\": {}, \"transport_exact\": {}, \"off_E_hits\": {}, \"log_mismatch\": {}, \"ideal\": {}}}",
                    if mi == 0 { "" } else { ",\n      " }, m, base.len(), m * m * sx.len(), jf(base.len() as f64 / (m * m * sx.len()) as f64), jf(c_pre), jf(c_pre - c_sx), transport_ok, off_e_hits, log_mismatch, ideal_json);
            }
            // LVL_m
            let mut lvl_json = String::new();
            for (mi, &m) in ms_lvl.iter().enumerate() {
                let g = div.psi_odd(m); let dg = deg(&g);
                let vals: Vec<(u64, bool)> = (0..p).map(|x| (peval(&g, x, p), toy.logx.contains_key(&x))).collect();
                let mut rl = Rng::new(4242 + seed * 13 + m as u64 + bits as u64);
                let mut sizes = vec![]; let mut covs = vec![];
                for _ in 0..50 { let cst = rl.below(p);
                    let base: Vec<u64> = (0..p).filter(|&x| vals[x as usize].0 == cst && vals[x as usize].1).collect();
                    sizes.push(base.len() as f64); covs.push(if base.is_empty() { 0.0 } else { coverage(&logs_of(&base), n) }); }
                let ms = sizes.iter().sum::<f64>() / 50.0; let mc = covs.iter().sum::<f64>() / 50.0; let mx = sizes.iter().cloned().fold(0.0, f64::max);
                summary.entry((bits, format!("LVL_{}_rho_mean", m))).or_default().push(ms / dg as f64);
                let _ = write!(lvl_json, "{}{{\"m\": {}, \"membership_degree\": {}, \"mean_size\": {}, \"max_size\": {}, \"rho\": {}, \"mean_coverage\": {}}}", if mi == 0 { "" } else { ",\n      " }, m, dg, jf(ms), mx, jf(ms / dg as f64), jf(mc));
            }
            summary.entry((bits, "RND_minus_SX".into())).or_default().push(c_rnd - c_sx);
            summary.entry((bits, "KL_minus_SX".into())).or_default().push(c_kl - c_sx);
            summary.entry((bits, "SX_coverage".into())).or_default().push(c_sx);
            let _ = write!(out, "{}  {{\"bits\": {}, \"seed\": {}, \"p\": {}, \"a\": {}, \"b\": {}, \"N\": {}, \"B\": {}, \"gen\": \"{:?}\",\n    \"SX\": {{\"coverage\": {}, \"rho\": 1.0}}, \"RND\": {{\"coverage\": {}, \"rho\": 1.0}}, \"KL\": {{\"coverage\": {}, \"rho\": 1.0}},\n    \"PRE\": [{}],\n    \"LVL\": [{}], \"secs\": {:.2}}}",
                if first { "" } else { ",\n" }, bits, seed, p, c.a, c.b, n, bsz, toy.gen, jf(c_sx), jf(c_rnd), jf(c_kl), pre_json, lvl_json, ts.elapsed().as_secs_f64());
            first = false;
            eprintln!("toy bits={} seed={} done {:.1}s", bits, seed, t0.elapsed().as_secs_f64());
        }
    }
    out.push_str("\n],\n");

    // ================= TOR_l positive control =================
    out.push_str("\"torsion\": [\n"); let mut firstt = true;
    for bits in [12u32, 14, 16] { for l in [3u64, 5, 7] { for seed in 0..3u64 {
        let (c, n, q) = find_lq_curve(bits, l, seed); let p = c.p;
        let mut r = Rng::new(555 + seed + l + bits as u64);
        let gen = loop { let g = c.random_point(&mut r); if c.mul(l as u128, g) != Pt::O && c.mul(q as u128, g) != Pt::O { break g; } };
        let toy = build_toy(c, n, gen);
        let div = Div::new(c.a, c.b, p, 8);
        let g = div.psi_odd(l as usize);
        let mut on_e = vec![]; let mut on_t = 0;
        for x in 0..p { if peval(&g, x, p) == 0 { if toy.logx.contains_key(&x) { on_e.push(x); } else if legendre(c.rhs(x), p) == -1 { on_t += 1; } } }
        let tors = toy.pts.iter().filter(|&&pt| c.mul(l as u128, pt) == Pt::O).count() as u64; // includes O
        let lg: Vec<u64> = on_e.iter().map(|x| toy.logx[x]).collect();
        let t = coverage_bits(&lg, n);
        let in_q: u64 = (0..n).filter(|&k| k % l == 0 && t[k as usize]).count() as u64;
        let ttot = t.iter().filter(|&&b| b).count();
        let _ = write!(out, "{}  {{\"bits\": {}, \"l\": {}, \"seed\": {}, \"p\": {}, \"N\": {}, \"q\": {}, \"psi_deg\": {}, \"roots_on_E\": {}, \"roots_on_twist\": {}, \"E_l_torsion_incl_O\": {}, \"predicted_roots\": {}, \"targets_covered_total\": {}, \"q_subgroup_covered\": {}, \"q_subgroup_coverage\": {}}}",
            if firstt { "" } else { ",\n" }, bits, l, seed, p, n, q, deg(&g), on_e.len(), on_t, tors, (tors - 1) / 2, ttot, in_q, jf(in_q as f64 / q as f64));
        firstt = false;
    } } }
    out.push_str("\n],\n");
    eprintln!("torsion done {:.1}s", t0.elapsed().as_secs_f64());

    // ================= summary with bootstrap =================
    out.push_str("\"summary\": [\n"); let mut keys: Vec<_> = summary.keys().cloned().collect(); keys.sort();
    for (i, k) in keys.iter().enumerate() { let (m, lo, hi) = bootstrap_mean(&summary[k], 7);
        let _ = write!(out, "  {{\"bits\": {}, \"quantity\": \"{}\", \"n\": {}, \"mean\": {}, \"ci95\": [{}, {}]}}{}\n", k.0, k.1, summary[k].len(), jf(m), jf(lo), jf(hi), if i + 1 < keys.len() { "," } else { "" }); }
    let _ = write!(out, "],\n\"wall_secs\": {:.2}\n}}\n", t0.elapsed().as_secs_f64());
    print!("{}", out);
}
