// 04-family-sweep: factor-base properties of the family p = 5t^2+19t+19
// (CM disc -19, t^2-4p = -19 (t+2)^2) against matched nulls.
// Usage: family_sweep <out_dir> <T_targets> <min_bits> <max_bits> <perms> <boots>
// Everything (generation, instrument, statistics) is computed here; Python only plots.
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::io::Write;
use std::time::Instant;

const SEED: u64 = 46004;
const J19: i64 = -884736; // -96^3

// ---------- arithmetic ----------
#[inline]
fn mulm(a: u64, b: u64, p: u64) -> u64 { if p < (1u64 << 32) { (a * b) % p } else { ((a as u128 * b as u128) % p as u128) as u64 } }
#[inline]
fn addm(a: u64, b: u64, p: u64) -> u64 { let s = a as u128 + b as u128; (s % p as u128) as u64 }
#[inline]
fn subm(a: u64, b: u64, p: u64) -> u64 { if a >= b { a - b } else { p - (b - a) } }
fn powm(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p; a %= p;
    while e > 0 { if e & 1 == 1 { r = mulm(r, a, p); } a = mulm(a, a, p); e >>= 1; }
    r
}
fn invm(a: u64, p: u64) -> u64 {
    if p < (1u64 << 62) {
        let (mut r0, mut r1) = (p as i64, (a % p) as i64);
        let (mut s0, mut s1) = (0i64, 1i64);
        while r1 != 0 { let q = r0 / r1; let t = r0 - q * r1; r0 = r1; r1 = t; let t = s0 - q * s1; s0 = s1; s1 = t; }
        assert!(r0 == 1, "non-invertible");
        return s0.rem_euclid(p as i64) as u64;
    }
    let (mut r0, mut r1) = (p as i128, (a % p) as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 { let q = r0 / r1; let t = r0 - q * r1; r0 = r1; r1 = t; let t = s0 - q * s1; s0 = s1; s1 = t; }
    assert!(r0 == 1, "non-invertible");
    (((s0 % p as i128) + p as i128) % p as i128) as u64
}
fn modi(x: i64, p: u64) -> u64 { let m = (x as i128).rem_euclid(p as i128); m as u64 }
fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % q == 0 { return n == q; } }
    let mut d = n - 1; let mut s = 0; while d % 2 == 0 { d /= 2; s += 1; }
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulm(x, x, n); if x == n - 1 { continue 'w; } }
        return false;
    }
    true
}
fn jacobi(mut a: u64, mut n: u64) -> i32 {
    a %= n; let mut r = 1;
    while a != 0 {
        let tz = a.trailing_zeros(); a >>= tz;
        if tz & 1 == 1 && (n % 8 == 3 || n % 8 == 5) { r = -r; }
        if a % 4 == 3 && n % 4 == 3 { r = -r; }
        std::mem::swap(&mut a, &mut n); a %= n;
    }
    if n == 1 { r } else { 0 }
}
fn sqrtm(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 { return Some(0); }
    if jacobi(a, p) != 1 { return None; }
    if p % 4 == 3 { return Some(powm(a, (p + 1) / 4, p)); }
    let mut q = p - 1; let mut s = 0; while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while jacobi(z, p) != -1 { z += 1; }
    let mut m = s; let mut c = powm(z, q, p); let mut t = powm(a, q, p); let mut r = powm(a, (q + 1) / 2, p);
    while t != 1 {
        let mut i = 0; let mut tt = t; while tt != 1 { tt = mulm(tt, tt, p); i += 1; }
        let mut b = c; for _ in 0..(m - i - 1) { b = mulm(b, b, p); }
        m = i; c = mulm(b, b, p); t = mulm(t, c, p); r = mulm(r, b, p);
    }
    Some(r)
}
fn isqrt(n: u64) -> u64 { let mut x = (n as f64).sqrt() as u64; while x * x > n { x -= 1; } while (x + 1) * (x + 1) <= n { x += 1; } x }

// ---------- rng ----------
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: u64) -> u64 { (self.next() as u128 * n as u128 >> 64) as u64 }
    fn f64(&mut self) -> f64 { (self.next() >> 11) as f64 / (1u64 << 53) as f64 }
}

// ---------- curves ----------
type Pt = Option<(u64, u64)>;
#[derive(Clone, Copy)]
struct Curve { p: u64, a: u64, b: u64 }
impl Curve {
    fn rhs(&self, x: u64) -> u64 { let p = self.p; addm(addm(mulm(mulm(x, x, p), x, p), mulm(self.a, x, p), p), self.b, p) }
    fn add(&self, P: Pt, Q: Pt) -> Pt {
        let p = self.p;
        match (P, Q) {
            (None, _) => Q, (_, None) => P,
            (Some((x1, y1)), Some((x2, y2))) => {
                let l = if x1 == x2 {
                    if addm(y1, y2, p) == 0 { return None; }
                    mulm(addm(mulm(3, mulm(x1, x1, p), p), self.a, p), invm(mulm(2, y1, p), p), p)
                } else { mulm(subm(y2, y1, p), invm(subm(x2, x1, p), p), p) };
                let x3 = subm(subm(mulm(l, l, p), x1, p), x2, p);
                let y3 = subm(mulm(l, subm(x1, x3, p), p), y1, p);
                Some((x3, y3))
            }
        }
    }
    fn neg(&self, P: Pt) -> Pt { P.map(|(x, y)| (x, if y == 0 { 0 } else { self.p - y })) }
    fn mul(&self, mut k: u64, P: Pt) -> Pt { let mut R = None; let mut Q = P; while k > 0 { if k & 1 == 1 { R = self.add(R, Q); } Q = self.add(Q, Q); k >>= 1; } R }
    fn rand_point(&self, rng: &mut Rng) -> (u64, u64) {
        loop { let x = rng.below(self.p); let r = self.rhs(x); if r != 0 { if let Some(y) = sqrtm(r, self.p) { return (x, y); } } }
    }
    fn disc_ok(&self) -> bool { let p = self.p; let d = addm(mulm(4, powm(self.a, 3, p), p), mulm(27, mulm(self.b, self.b, p), p), p); d != 0 }
    fn j(&self) -> u64 { let p = self.p; let a3 = mulm(4, powm(self.a, 3, p), p); let d = addm(a3, mulm(27, mulm(self.b, self.b, p), p), p); mulm(mulm(1728, a3, p), invm(d, p), p) }
    fn rescale(&self, c: u64) -> Curve { let p = self.p; let c2 = mulm(c, c, p); Curve { p, a: mulm(self.a, mulm(c2, c2, p), p), b: mulm(self.b, mulm(c2, mulm(c2, c2, p), p), p) } }
}

fn order_is(E: &Curve, N: u64, rng: &mut Rng) -> bool {
    for _ in 0..3 { let P = Some(E.rand_point(rng)); if E.mul(N, P).is_some() { return false; } }
    true
}

/// j = -96^3 curve over F_p with #E = N (N prime), random model.
fn cm19_curve(p: u64, N: u64, rng: &mut Rng) -> Curve {
    let j = modi(J19, p);
    let k = mulm(j, invm(subm(1728 % p, j, p), p), p);
    let E0 = Curve { p, a: mulm(3, k, p), b: mulm(2, k, p) };
    assert!(E0.j() == j);
    let E = if order_is(&E0, N, rng) { E0 } else {
        let mut d = 2; while jacobi(d, p) != -1 { d += 1; }
        let d2 = mulm(d, d, p);
        let Et = Curve { p, a: mulm(E0.a, d2, p), b: mulm(E0.b, mulm(d2, d, p), p) };
        assert!(order_is(&Et, N, rng), "neither twist has order N");
        Et
    };
    let c = 1 + rng.below(p - 1);
    E.rescale(c)
}

fn rand_prime(bits: u32, rng: &mut Rng) -> u64 { loop { let p = (1u64 << (bits - 1)) | rng.below(1u64 << (bits - 1)) | 1; if is_prime(p) { return p; } } }

/// 4p = a^2 + 19 y^2 with a,y >= 0 (unique when it exists, h(-19)=1)
fn rep19(p: u64) -> Option<(u64, u64)> {
    let mut y = 1u64; while 19 * y * y <= 4 * p { let r = 4 * p - 19 * y * y; let a = isqrt(r); if a * a == r { return Some((a, y)); } y += 1; } None
}
fn in_family(a: u64, y: u64) -> bool { y == a + 2 || (a >= 2 && y == a - 2) || (a < 2 && y == 2 - a) }

fn null_cm19(bits: u32, rng: &mut Rng) -> (Curve, u64, i64) {
    loop {
        let p = rand_prime(bits, rng);
        if jacobi(modi(-19, p), p) != 1 { continue; }
        let (a, y) = match rep19(p) { Some(v) => v, None => continue };
        if in_family(a, y) { continue; }
        let mut cands = vec![]; if is_prime(p + 1 - a) { cands.push(a as i64); } if is_prime(p + 1 + a) { cands.push(-(a as i64)); }
        if cands.is_empty() { continue; }
        let t = cands[rng.below(cands.len() as u64) as usize];
        let N = (p as i64 + 1 - t) as u64;
        return (cm19_curve(p, N, rng), N, t);
    }
}

/// all k in Hasse interval with kP = O (None if small order detected)
fn bsgs_orders(E: &Curve, P: Pt) -> Option<Vec<u64>> {
    let p = E.p; let s = isqrt(p);
    let L = p + 1 - 2 * s - 2; let U = p + 1 + 2 * s + 2; let W = U - L;
    let m = isqrt(W) + 1;
    let mut tab: HashMap<(u64, u64), u64> = HashMap::with_capacity(m as usize + 1);
    let mut R: Pt = None;
    for j in 0..m { let key = R.unwrap_or((u64::MAX, u64::MAX)); if tab.insert(key, j).is_some() { return None; } R = E.add(R, P); }
    let S = E.mul(m, P); let nS = E.neg(S);
    let mut T = E.neg(E.mul(L, P));
    let mut hits = vec![];
    for g in 0..=(W / m + 1) { let key = T.unwrap_or((u64::MAX, u64::MAX)); if let Some(&j) = tab.get(&key) { let k = L + g * m + j; if k <= U { hits.push(k); } } T = E.add(T, nS); }
    Some(hits)
}
fn null_random(bits: u32, rng: &mut Rng) -> (Curve, u64) {
    loop {
        let p = rand_prime(bits, rng);
        let E = Curve { p, a: rng.below(p), b: rng.below(p) };
        if !E.disc_ok() { continue; }
        let P = Some(E.rand_point(rng));
        if let Some(h) = bsgs_orders(&E, P) { if h.len() == 1 && is_prime(h[0]) && order_is(&E, h[0], rng) { return (E, h[0]); } }
    }
}

// ---------- instrument (identical on every arm) ----------
struct Meas { fx_nom: u64, b_nom: u64, fx_mem: u64, dec: u64, rel: u64, targets: u64, lam: f64 }
/// Semaev m=2 decomposition over the small-x factor base {x < b_mem}; metrics normalised with |{x<b_nom}|.
fn measure(E: &Curve, N: u64, dope: f64, T: u64, rng: &mut Rng) -> Meas {
    let p = E.p;
    let b_nom = ((2.0 * N as f64).sqrt()).ceil() as u64;
    let b_mem = ((b_nom as f64) * dope).ceil() as u64;
    let mut member = vec![false; b_mem as usize];
    let mut fx = vec![]; let mut fx_nom = 0u64;
    for x in 0..b_mem { if jacobi(E.rhs(x), p) == 1 { member[x as usize] = true; fx.push(x); if x < b_nom { fx_nom += 1; } } }
    let (mut dec, mut rel) = (0u64, 0u64);
    let mut pend: Vec<(u64, u64, u64, u64, bool)> = Vec::with_capacity(fx.len());
    let mut pre: Vec<u64> = Vec::with_capacity(fx.len());
    for _ in 0..T {
        let (x3, _) = E.rand_point(rng);
        let mut c = 0u64;
        for &x1 in &fx {
            if x1 == x3 { continue; }
            let d13 = subm(x1, x3, p); let A = mulm(d13, d13, p);
            let s13 = addm(x1, x3, p); let p13 = mulm(x1, x3, p);
            let Bc = subm(0, mulm(2, addm(mulm(s13, addm(p13, E.a, p), p), mulm(2, E.b, p), p), p), p);
            let t = subm(p13, E.a, p);
            let C = subm(mulm(t, t, p), mulm(mulm(4, E.b, p), s13, p), p);
            let disc = subm(mulm(Bc, Bc, p), mulm(4, mulm(A, C, p), p), p);
            if disc == 0 { pend.push((x1, subm(0, Bc, p), 0, A, false)); }
            else if jacobi(disc, p) == 1 { pend.push((x1, subm(0, Bc, p), sqrtm(disc, p).unwrap(), A, true)); }
        }
        // batch inversion of 2A (Montgomery trick); identical arithmetic result to per-root inversion
        let n = pend.len();
        if n > 0 {
            pre.clear(); let mut acc = 1u64;
            for e in &pend { pre.push(acc); acc = mulm(acc, mulm(2, e.3, p), p); }
            let mut inv = invm(acc, p);
            for k in (0..n).rev() {
                let e = pend[k]; let i = mulm(inv, pre[k], p); inv = mulm(inv, mulm(2, e.3, p), p);
                let x1 = e.0;
                let mut chk = |x2: u64| { if x2 >= x1 && x2 < b_mem && member[x2 as usize] { c += 1; } };
                if e.4 { chk(mulm(addm(e.1, e.2, p), i, p)); chk(mulm(subm(e.1, e.2, p), i, p)); } else { chk(mulm(e.1, i, p)); }
            }
        }
        pend.clear();
        dec += c; if c > 0 { rel += 1; }
    }
    let lam = 2.0 * (fx_nom as f64).powi(2) / N as f64;
    Meas { fx_nom, b_nom, fx_mem: fx.len() as u64, dec, rel, targets: T, lam }
}

/// exact self-test of the Semaev S3 formula used above
fn selftest(rng: &mut Rng) {
    let (E, _) = null_random(24, rng);
    for _ in 0..200 {
        let P = Some(E.rand_point(rng)); let Q = Some(E.rand_point(rng)); let R = E.add(P, Q);
        if let (Some((x1, _)), Some((x2, _)), Some((x3, _))) = (P, Q, R) {
            let p = E.p; let d = subm(x1, x2, p);
            let s = addm(x1, x2, p); let pr = mulm(x1, x2, p); let t = subm(pr, E.a, p);
            let v = addm(subm(mulm(mulm(d, d, p), mulm(x3, x3, p), p), mulm(mulm(2, addm(mulm(s, addm(pr, E.a, p), p), mulm(2, E.b, p), p), p), x3, p), p), subm(mulm(t, t, p), mulm(mulm(4, E.b, p), s, p), p), p);
            assert!(v == 0, "S3 self-test failed");
        }
    }
}

// ---------- statistics ----------
#[derive(Clone)]
struct Row { arm: &'static str, bits: u32, p: u64, N: u64, t: i64, v_prime: Option<bool>, m: Meas_s }
#[derive(Clone, Copy)]
struct Meas_s { m1: f64, m2: f64, m3: f64, fx_nom: u64, b_nom: u64, fx_mem: u64, dec: u64, rel: u64, targets: u64, lam: f64, a: u64, b: u64 }

fn mean(v: &[f64]) -> f64 { v.iter().sum::<f64>() / v.len() as f64 }
/// stratified permutation test + stratified bootstrap for mean(X)-mean(Y)
fn compare(x: &[(u32, f64)], y: &[(u32, f64)], perms: u32, boots: u32, rng: &mut Rng) -> (f64, f64, f64, f64, f64, f64, f64) {
    let xv: Vec<f64> = x.iter().map(|r| r.1).collect(); let yv: Vec<f64> = y.iter().map(|r| r.1).collect();
    let obs = mean(&xv) - mean(&yv);
    let mut strata: HashMap<u32, (Vec<f64>, usize)> = HashMap::new();
    for r in x { let e = strata.entry(r.0).or_insert((vec![], 0)); e.0.push(r.1); e.1 += 1; }
    for r in y { strata.entry(r.0).or_insert((vec![], 0)).0.push(r.1); }
    let nx = x.len() as f64; let ny = y.len() as f64;
    let mut keys: Vec<u32> = strata.keys().cloned().collect(); keys.sort();
    let mut ge = 0u32;
    for _ in 0..perms {
        let (mut sx, mut sy) = (0.0, 0.0);
        for k in &keys { let (vals, kx) = strata.get_mut(k).unwrap(); let n = vals.len();
            for i in 0..*kx { let j = i + rng.below((n - i) as u64) as usize; vals.swap(i, j); }
            for (i, v) in vals.iter().enumerate() { if i < *kx { sx += v; } else { sy += v; } } }
        if ((sx / nx - sy / ny).abs()) >= obs.abs() - 1e-15 { ge += 1; }
    }
    let pval = (1.0 + ge as f64) / (1.0 + perms as f64);
    // stratified bootstrap
    let mut bx: HashMap<u32, Vec<f64>> = HashMap::new(); let mut by: HashMap<u32, Vec<f64>> = HashMap::new();
    for r in x { bx.entry(r.0).or_default().push(r.1); } for r in y { by.entry(r.0).or_default().push(r.1); }
    let mut diffs = vec![]; let mut ratios = vec![];
    for _ in 0..boots {
        let (mut sx, mut sy) = (0.0, 0.0);
        for v in bx.values() { for _ in 0..v.len() { sx += v[rng.below(v.len() as u64) as usize]; } }
        for v in by.values() { for _ in 0..v.len() { sy += v[rng.below(v.len() as u64) as usize]; } }
        diffs.push(sx / nx - sy / ny); ratios.push((sx / nx) / (sy / ny));
    }
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap()); ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |v: &Vec<f64>, f: f64| v[((v.len() as f64 - 1.0) * f).round() as usize];
    (obs, q(&diffs, 0.025), q(&diffs, 0.975), mean(&xv) / mean(&yv), q(&ratios, 0.025), q(&ratios, 0.975), pval)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = &args[1]; let T: u64 = args[2].parse().unwrap();
    let minb: u32 = args[3].parse().unwrap(); let maxb: u32 = args[4].parse().unwrap();
    let perms: u32 = args[5].parse().unwrap(); let boots: u32 = args[6].parse().unwrap();
    let t0 = Instant::now();
    let mut rng = Rng(SEED);
    selftest(&mut rng);
    let mut log = std::fs::File::create(format!("{}/run.log", out)).unwrap();
    macro_rules! say { ($($a:tt)*) => {{ let s = format!($($a)*); println!("{}", s); writeln!(log, "{}", s).unwrap(); log.flush().unwrap(); }} }
    say!("family_sweep seed={} T={} bits={}..={} perms={} boots={}", SEED, T, minb, maxb, perms, boots);
    say!("S3 self-test passed (200 random sums on a 24-bit curve)");

    // pair-46 sanity: identity and orders (exact checks)
    let p46: u64 = 5947698030651343; let t46: u64 = 34489701; let N46: u64 = 5947697996161643;
    assert!(p46 == 5 * t46 * t46 + 19 * t46 + 19 && N46 == p46 + 1 - t46);
    let crater = Curve { p: p46, a: 5757116383569871, b: 3560173328417442 };
    let floor = Curve { p: p46, a: 2087370704375809, b: 1402887937802154 };
    let jc = crater.j(); let crater_j_ok = jc == modi(J19, p46);
    let co = order_is(&crater, N46, &mut rng); let fo = order_is(&floor, N46, &mut rng);
    say!("pair46: p=5t^2+19t+19 ok; crater j==-96^3: {}; crater #E=N: {}; floor #E=N: {}; floor j==-96^3: {}", crater_j_ok, co, fo, floor.j() == modi(J19, p46));

    // family members
    let mut fam = vec![];
    let mut t = 1u64;
    loop { let p = 5 * t * t + 19 * t + 19; let b = 64 - p.leading_zeros(); if b > maxb { break; }
        if b >= minb && is_prime(p) && is_prime(p + 1 - t) { fam.push((b, p, t)); } t += 1; }
    say!("family members with p, N=p+1-t prime in {}..={} bits: {}", minb, maxb, fam.len());

    let mut rows: Vec<Row> = vec![];
    let mk = |arm: &'static str, bits: u32, E: &Curve, N: u64, t: i64, vp: Option<bool>, m: Meas| Row { arm, bits, p: E.p, N, t, v_prime: vp,
        m: Meas_s { m1: m.dec as f64 / (m.targets as f64 * m.lam), m2: m.rel as f64 / m.targets as f64, m3: m.fx_nom as f64 / (m.b_nom as f64 / 2.0),
            fx_nom: m.fx_nom, b_nom: m.b_nom, fx_mem: m.fx_mem, dec: m.dec, rel: m.rel, targets: m.targets, lam: m.lam, a: E.a, b: E.b } };
    for (i, &(b, p, t)) in fam.iter().enumerate() {
        let N = p + 1 - t;
        let mut r = Rng(SEED ^ (0x1000 + i as u64) * 0x9E37);
        let E = cm19_curve(p, N, &mut r);
        let m = measure(&E, N, 1.0, T, &mut r); rows.push(mk("F", b, &E, N, t as i64, Some(is_prime(t + 2)), m));
        for k in 0..2u64 {
            let mut r1 = Rng(SEED ^ (0x2000 + 2 * i as u64 + k) * 0x9E37);
            let (E1, N1, t1) = null_cm19(b, &mut r1);
            let m = measure(&E1, N1, 1.0, T, &mut r1); rows.push(mk("N1", b, &E1, N1, t1, None, m));
            let mut r2 = Rng(SEED ^ (0x3000 + 2 * i as u64 + k) * 0x9E37);
            let (E2, N2) = null_random(b, &mut r2);
            let m = measure(&E2, N2, 1.0, T, &mut r2); rows.push(mk("N2", b, &E2, N2, p as i64 + 1 - N2 as i64, None, m));
            if k == 0 {
                // positive control: same N2 curve, best-of-2 model by small-x count, factor base doped by 3%
                let mut r3 = Rng(SEED ^ (0x4000 + i as u64) * 0x9E37);
                let alt = E2.rescale(1 + r3.below(E2.p - 1));
                let bn = ((2.0 * N2 as f64).sqrt()).ceil() as u64;
                let cnt = |E: &Curve| (0..bn).filter(|&x| jacobi(E.rhs(x), E.p) == 1).count();
                let Eb = if cnt(&alt) > cnt(&E2) { alt } else { E2 };
                let m = measure(&Eb, N2, 1.03, T, &mut r3); rows.push(mk("PC", b, &Eb, N2, p as i64 + 1 - N2 as i64, None, m));
            }
        }
        if i % 10 == 9 { say!("  done {} / {} family members, elapsed {:.1}s", i + 1, fam.len(), t0.elapsed().as_secs_f64()); }
    }
    let gen_secs = t0.elapsed().as_secs_f64();

    // pair46-scale small-x add-on (M3 only, B = 2^16)
    let B46 = 1u64 << 16;
    let cnt46 = |E: &Curve| (0..B46).filter(|&x| jacobi(E.rhs(x), E.p) == 1).count() as f64 / (B46 as f64 / 2.0);
    let c_val = cnt46(&crater); let f_val = cnt46(&floor);
    let mut r4 = Rng(SEED ^ 0x5000);
    let mut ctrl_models = vec![]; let mut ctrl_rand = vec![];
    for _ in 0..200 { let c = 1 + r4.below(p46 - 1); ctrl_models.push(cnt46(&crater.rescale(c))); }
    for _ in 0..200 { let E = Curve { p: p46, a: r4.below(p46), b: r4.below(p46) }; ctrl_rand.push(cnt46(&E)); }
    let rank = |v: &Vec<f64>, x: f64| v.iter().filter(|&&y| y < x).count();

    // statistics
    let get = |arm: &str, f: &dyn Fn(&Meas_s) -> f64, lo: u32, hi: u32| -> Vec<(u32, f64)> { rows.iter().filter(|r| r.arm == arm && r.bits >= lo && r.bits <= hi).map(|r| (r.bits, f(&r.m))).collect() };
    let metrics: Vec<(&str, Box<dyn Fn(&Meas_s) -> f64>)> = vec![("M1_norm_decomp_yield", Box::new(|m: &Meas_s| m.m1)), ("M2_relation_rate", Box::new(|m: &Meas_s| m.m2)), ("M3_smallx_yield", Box::new(|m: &Meas_s| m.m3))];
    let mut js = String::new();
    write!(js, "{{\n \"experiment\": \"04-family-sweep\", \"seed\": {}, \"targets_per_curve\": {}, \"bits\": [{}, {}], \"perms\": {}, \"boots\": {},\n", SEED, T, minb, maxb, perms, boots).unwrap();
    write!(js, " \"pair46_checks\": {{\"crater_j_is_-96^3\": {}, \"crater_order_N\": {}, \"floor_order_N\": {}, \"floor_j_is_-96^3\": {}}},\n", crater_j_ok, co, fo, floor.j() == modi(J19, p46)).unwrap();
    write!(js, " \"n_family\": {}, \"counts\": {{\"F\": {}, \"N1\": {}, \"N2\": {}, \"PC\": {}}},\n", fam.len(), rows.iter().filter(|r| r.arm == "F").count(), rows.iter().filter(|r| r.arm == "N1").count(), rows.iter().filter(|r| r.arm == "N2").count(), rows.iter().filter(|r| r.arm == "PC").count()).unwrap();
    js.push_str(" \"tests\": [\n");
    let alpha = 0.05 / 6.0;
    let mut first = true;
    let bands = [(minb, maxb, "all"), (20, 24, "20-24"), (25, 27, "25-27"), (28, 30, "28-30")];
    for (lo, hi, bname) in bands.iter() {
        for (mname, f) in metrics.iter() {
            for (xa, ya, role) in [("F", "N1", "primary"), ("F", "N2", "primary"), ("PC", "N2", "positive_control"), ("N1", "N2", "secondary")] {
                let x = get(xa, f.as_ref(), *lo, *hi); let y = get(ya, f.as_ref(), *lo, *hi);
                if x.is_empty() || y.is_empty() { continue; }
                let role = if *bname != "all" { "per_band_descriptive" } else { role };
                let (d, dl, dh, ra, rl, rh, pv) = compare(&x, &y, perms, boots, &mut rng);
                let xm = mean(&x.iter().map(|r| r.1).collect::<Vec<_>>()); let ym = mean(&y.iter().map(|r| r.1).collect::<Vec<_>>());
                let sig = pv < alpha;
                say!("[{}] {:6} {:22} {:>2} vs {:>2}: n={}/{} mean {:.5} vs {:.5} diff {:+.5} [{:+.5},{:+.5}] ratio {:.4} [{:.4},{:.4}] p={:.4} {}", role, bname, mname, xa, ya, x.len(), y.len(), xm, ym, d, dl, dh, ra, rl, rh, pv, if sig { "SIG(bonf 6)" } else { "" });
                if !first { js.push_str(",\n"); } first = false;
                write!(js, "  {{\"role\": \"{}\", \"band\": \"{}\", \"metric\": \"{}\", \"x\": \"{}\", \"y\": \"{}\", \"nx\": {}, \"ny\": {}, \"mean_x\": {:.6}, \"mean_y\": {:.6}, \"diff\": {:.6}, \"diff_ci95\": [{:.6}, {:.6}], \"ratio\": {:.6}, \"ratio_ci95\": [{:.6}, {:.6}], \"p_perm\": {:.6}, \"bonferroni_alpha\": {:.6}, \"significant\": {}}}", role, bname, mname, xa, ya, x.len(), y.len(), xm, ym, d, dl, dh, ra, rl, rh, pv, alpha, sig).unwrap();
            }
        }
    }
    js.push_str("\n ],\n");
    let med = |v: &Vec<f64>| { let mut s = v.clone(); s.sort_by(|a, b| a.partial_cmp(b).unwrap()); (s[s.len() / 2], s[(s.len() as f64 * 0.025) as usize], s[(s.len() as f64 * 0.975) as usize]) };
    let (m1, l1, h1) = med(&ctrl_models); let (m2, l2, h2) = med(&ctrl_rand);
    say!("pair46 M3 (B=2^16): crater {:.5} rank {}/200 among crater models (median {:.5} [{:.5},{:.5}]), rank {}/200 among random (a,b) (median {:.5} [{:.5},{:.5}]); floor {:.5} rank {}/200 random (a,b)", c_val, rank(&ctrl_models, c_val), m1, l1, h1, rank(&ctrl_rand, c_val), m2, l2, h2, f_val, rank(&ctrl_rand, f_val));
    write!(js, " \"pair46_smallx\": {{\"B\": {}, \"crater\": {:.6}, \"floor\": {:.6}, \"crater_rank_models\": {}, \"crater_rank_random\": {}, \"floor_rank_random\": {}, \"n_ctrl\": 200, \"models_median_ci\": [{:.6}, {:.6}, {:.6}], \"random_median_ci\": [{:.6}, {:.6}, {:.6}], \"ctrl_models\": {:?}, \"ctrl_random\": {:?}}},\n", B46, c_val, f_val, rank(&ctrl_models, c_val), rank(&ctrl_rand, c_val), rank(&ctrl_rand, f_val), m1, l1, h1, m2, l2, h2, ctrl_models, ctrl_rand).unwrap();
    write!(js, " \"runtime_s\": {{\"generation_and_instrument\": {:.1}, \"total\": {:.1}}}\n}}\n", gen_secs, t0.elapsed().as_secs_f64()).unwrap();
    std::fs::write(format!("{}/results.json", out), js).unwrap();
    let mut csv = String::from("arm,bits,p,N,t,v_prime,a,b,b_nom,fx_nom,fx_mem,targets,decomps,relations,lambda,M1,M2,M3\n");
    for r in &rows { writeln!(csv, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6}", r.arm, r.bits, r.p, r.N, r.t, r.v_prime.map(|b| b.to_string()).unwrap_or("".into()), r.m.a, r.m.b, r.m.b_nom, r.m.fx_nom, r.m.fx_mem, r.m.targets, r.m.dec, r.m.rel, r.m.lam, r.m.m1, r.m.m2, r.m.m3).unwrap(); }
    std::fs::write(format!("{}/curves.csv", out), csv).unwrap();
    say!("total runtime {:.1}s", t0.elapsed().as_secs_f64());
}
