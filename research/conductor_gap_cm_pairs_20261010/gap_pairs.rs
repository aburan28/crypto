// Prototype: (crater, floor) isogenous curve pairs with a certified prime conductor gap.
// D = -19, h(D)=1, j_D = -884736.  q = (t^2 + 19 l^2)/4 prime, N = q+1-t prime, so t^2-4q = -19 l^2.
// Floor curves found by trace test with ONE scalar multiplication per random curve.
use std::time::Instant;

fn mulmod(a: u64, b: u64, m: u64) -> u64 { ((a as u128 * b as u128) % m as u128) as u64 }
fn addmod(a: u64, b: u64, m: u64) -> u64 { let s = a as u128 + b as u128; (s % m as u128) as u64 }
fn submod(a: u64, b: u64, m: u64) -> u64 { if a >= b { a - b } else { m - (b - a) } }
fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 { let mut r = 1; a %= m; while e > 0 { if e & 1 == 1 { r = mulmod(r, a, m); } a = mulmod(a, a, m); e >>= 1; } r }
fn inv(a: u64, m: u64) -> u64 { powmod(a, m - 2, m) }
fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % p == 0 { return n == p; } }
    let d0 = n - 1; let s = d0.trailing_zeros(); let d = d0 >> s;
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulmod(x, x, n); if x == n - 1 { continue 'w; } }
        return false;
    }
    true
}
fn legendre(a: u64, p: u64) -> i32 { let r = powmod(a, (p - 1) / 2, p); if r == 0 { 0 } else if r == 1 { 1 } else { -1 } }
fn sqrtmod(a: u64, p: u64) -> Option<u64> { // Tonelli-Shanks
    if a == 0 { return Some(0); } if legendre(a, p) != 1 { return None; }
    let mut q = p - 1; let mut s = 0; while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while legendre(z, p) != -1 { z += 1; }
    let mut m = s; let mut c = powmod(z, q, p); let mut t = powmod(a, q, p); let mut r = powmod(a, (q + 1) / 2, p);
    while t != 1 { let mut i = 0; let mut tt = t; while tt != 1 { tt = mulmod(tt, tt, p); i += 1; }
        let b = powmod(c, 1u64 << (m - i - 1), p); m = i; c = mulmod(b, b, p); t = mulmod(t, c, p); r = mulmod(r, b, p); }
    Some(r)
}
// affine point (None = infinity)
type Pt = Option<(u64, u64)>;
fn add(p: &Pt, q: &Pt, a: u64, m: u64) -> Pt {
    match (p, q) { (None, _) => *q, (_, None) => *p, (Some((x1, y1)), Some((x2, y2))) => {
        let lam = if x1 == x2 { if addmod(*y1, *y2, m) == 0 { return None; }
            mulmod(addmod(mulmod(3, mulmod(*x1, *x1, m), m), a, m), inv(mulmod(2, *y1, m), m), m)
        } else { mulmod(submod(*y2, *y1, m), inv(submod(*x2, *x1, m), m), m) };
        let x3 = submod(submod(mulmod(lam, lam, m), *x1, m), *x2, m);
        let y3 = submod(mulmod(lam, submod(*x1, x3, m), m), *y1, m);
        Some((x3, y3)) } }
}
fn mul(p: &Pt, mut k: u64, a: u64, m: u64) -> Pt { let mut r = None; let mut b = *p; while k > 0 { if k & 1 == 1 { r = add(&r, &b, a, m); } b = add(&b, &b, a, m); k >>= 1; } r }
fn jinv(a: u64, b: u64, p: u64) -> u64 { let a3 = mulmod(mulmod(a, a, p), a, p); let num = mulmod(1728, mulmod(4, a3, p), p);
    let den = addmod(mulmod(4, a3, p), mulmod(27, mulmod(b, b, p), p), p); mulmod(num, inv(den, p), p) }
struct Rng(u64);
impl Rng { fn next(&mut self) -> u64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; self.0 } }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let k: u32 = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(17);   // ell ~ 2^k
    let want: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(3);
    let d: u64 = 19; let jd_neg: u64 = 884736; // D = -19, h = 1, j = -96^3; D=-7 forces even t,f; D=-11 forces 3 | qN
    // search ell prime near 2^k, odd t small, q prime, N prime, N not anomalous / low embedding degree
    let mut ell = (1u64 << k) + 1; let (mut t, mut q) = (0u64, 0u64);
    'outer: loop { ell += 2; if !is_prime(ell) { continue; }
        let mut tt = 3u64; while tt < 4000 { let num = tt * tt + d * ell * ell; if num % 4 == 0 { let qq = num / 4;
            if is_prime(qq) { let n = qq + 1 - tt; if is_prime(n) {
                let mut ok = true; let mut acc = 1u64; for _ in 0..20 { acc = mulmod(acc, qq % n, n); if acc == 1 { ok = false; break; } }
                if ok { t = tt; q = qq; break 'outer; } } } } tt += 2; } }
    let n = q + 1 - t; let nt = q + 1 + t;
    println!("q = {q} ({} bits)  t = {t}  ell = {ell}  N = #E = {n} (prime)   t^2-4q = -19*ell^2: {}", 64 - q.leading_zeros(), (t as i128 * t as i128 - 4 * q as i128) == -(19 * ell as i128 * ell as i128));
    // crater curve from j = -3375
    let j = submod(0, jd_neg % q, q); let kk = mulmod(j, inv(submod(1728, j, q), q), q);
    let (mut a0, mut b0) = (mulmod(3, kk, q), mulmod(2, kk, q));
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let rand_pt = |a: u64, b: u64, rng: &mut Rng| -> Pt { loop { let x = rng.next() % q; let rhs = addmod(addmod(mulmod(mulmod(x, x, q), x, q), mulmod(a, x, q), q), b, q); if let Some(y) = sqrtmod(rhs, q) { return Some((x, y)); } } };
    let has_order = |a: u64, b: u64, ord: u64, rng: &mut Rng| -> bool { (0..8).all(|_| mul(&rand_pt(a, b, rng), ord, a, q).is_none()) };
    if !has_order(a0, b0, n, &mut rng) { let mut g = 2; while legendre(g, q) != -1 { g += 1; } a0 = mulmod(a0, mulmod(g, g, q), q); b0 = mulmod(b0, mulmod(mulmod(g, g, q), g, q), q); }
    println!("crater E0: a = {a0} b = {b0} j = {} (= -884736 mod q)  [N]P = O on 8 random P: {}   f_E = 1 (h(-19)=1)", jinv(a0, b0, q), has_order(a0, b0, n, &mut rng));
    // floor curves
    let t0 = Instant::now(); let mut tries = 0u64; let mut found = 0;
    while found < want { tries += 1;
        let a = rng.next() % q; let b = rng.next() % q;
        if addmod(mulmod(4, mulmod(mulmod(a, a, q), a, q), q), mulmod(27, mulmod(b, b, q), q), q) == 0 { continue; }
        let p = rand_pt(a, b, &mut rng); let r = mul(&p, n, a, q);
        let (ord, twist) = if r.is_none() { (n, false) } else if mul(&p, nt, a, q).is_none() { (nt, true) } else { continue; };
        if !has_order(a, b, ord, &mut rng) { continue; }
        if jinv(a, b, q) == j { continue; }
        // twist of a trace -t floor curve: twist it back to trace t so all curves share #E = N
        let (a, b) = if twist { let mut g = 2; while legendre(g, q) != -1 { g += 1; } (mulmod(a, mulmod(g, g, q), q), mulmod(b, mulmod(mulmod(g, g, q), g, q), q)) } else { (a, b) };
        found += 1;
        println!("floor E{found}: a = {a} b = {b} j = {}  #E = N confirmed by 8 random [N]P = O{}   f_E = ell = {ell}   ({tries} samples, {:.1} s)", jinv(a, b, q), if twist { " (via twist)" } else { "" }, t0.elapsed().as_secs_f64());
    }
    println!("samples/hit observed {:.0}, predicted |D|*ell/4 (both twists) = {:.0}", tries as f64 / want as f64, d as f64 * ell as f64 / 4.0);
    println!("Order certification note: [N]P=O for 8 random P is a probabilistic certificate; confirm with ellcard/Sage before use.");
}
