// Arithmetic, RNG and curve helpers copied verbatim from ../04-family-sweep/src/main.rs (lines 13-122, CM-19 builder omitted).

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

