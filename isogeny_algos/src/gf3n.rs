//! Ternary fields GF(3^n) = F_3[t]/(f), elements as a packed base-3 digit vector in a u64
//! (coefficient of t^i is the i-th ternary digit; n <= 40 so 3^n fits in u64 and the packed form
//! fits too). Addition/subtraction digitwise mod 3, multiplication by schoolbook convolution mod
//! the reduction polynomial, inversion by Fermat (a^(3^n - 2)), square roots when they exist.
//! Characteristic 3: the elliptic-curve work on top uses the general Weierstrass model.
use crate::bigint::Big;
use crate::field::{Field, Rng};

/// A digit vector c[0..n] (each 0, 1, 2) packed as sum c[i] 3^i.
fn to_digits(mut v: u64, n: u32) -> Vec<u8> {
    let mut d = vec![0u8; n as usize];
    for x in d.iter_mut() {
        *x = (v % 3) as u8;
        v /= 3;
    }
    d
}
fn from_digits(d: &[u8]) -> u64 {
    let mut v = 0u64;
    for &x in d.iter().rev() {
        v = v * 3 + x as u64;
    }
    v
}

#[derive(Clone, Debug)]
pub struct GF3n {
    pub n: u32,
    /// reduction polynomial f = t^n - red(t); red as digit vector of degree < n
    red: Vec<u8>,
    pow3n: u64,
}

impl GF3n {
    /// GF(3^n) with a low-weight monic irreducible t^n - red found by search (Rabin's test).
    pub fn new(n: u32) -> Self {
        assert!((1..=40).contains(&n));
        let pow3n = 3u64.pow(n);
        if n == 1 {
            return GF3n { n, red: vec![0], pow3n };
        }
        // search reduction polynomials of small weight: t^n = a + b t^k
        let mut best: Option<Vec<u8>> = None;
        'outer: for weight in 1..=3u32 {
            // try t^n - (c0 + c1 t) ... grow degree of the tail
            for packed in 1..pow3n {
                let d = to_digits(packed, n);
                if d.iter().map(|&x| (x != 0) as u32).sum::<u32>() > weight {
                    continue;
                }
                let cand = GF3n { n, red: d.clone(), pow3n };
                if cand.is_irreducible() {
                    best = Some(d);
                    break 'outer;
                }
            }
        }
        let red = best.expect("no irreducible found");
        GF3n { n, red, pow3n }
    }

    fn reduce(&self, mut c: Vec<u8>) -> u64 {
        let n = self.n as usize;
        // c has degree up to 2n-2; reduce top coefficients using t^n = red
        for i in (n..c.len()).rev() {
            let ci = c[i];
            if ci == 0 {
                continue;
            }
            c[i] = 0;
            for (k, &r) in self.red.iter().enumerate() {
                if r == 0 {
                    continue;
                }
                let idx = i - n + k;
                c[idx] = (c[idx] + ci * r) % 3;
            }
        }
        from_digits(&c[..n])
    }

    fn mul_digits(&self, a: u64, b: u64) -> u64 {
        let n = self.n as usize;
        let da = to_digits(a, self.n);
        let db = to_digits(b, self.n);
        let mut prod = vec![0u8; 2 * n];
        for i in 0..n {
            if da[i] == 0 {
                continue;
            }
            for j in 0..n {
                prod[i + j] = (prod[i + j] + da[i] * db[j]) % 3;
            }
        }
        self.reduce(prod)
    }

    fn pow(&self, a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul_digits(r, b);
            }
            b = self.mul_digits(b, b);
            e >>= 1;
        }
        r
    }

    /// Rabin irreducibility test for the monic degree-n polynomial t^n - red over F_3.
    fn is_irreducible(&self) -> bool {
        // x^(3^n) = x (mod f), and for each maximal proper divisor n/q: gcd(x^(3^(n/q)) - x, f) = 1
        let n = self.n;
        // represent polynomials over F_3 as digit vectors (variable-length); work with t = x
        let xq = |e: u32| self.frob_x(e); // x^(3^e) mod f as an element (digit vector length n)
        // x^(3^n) == x
        if xq(n) != self.coeff_x() {
            return false;
        }
        for p in distinct_primes(n) {
            let m = n / p;
            let d = self.sub(xq(m), self.coeff_x());
            if self.poly_gcd_is_unit(d) {
                // gcd(x^(3^m)-x, f) must be 1
            } else {
                return false;
            }
        }
        true
    }

    /// x^(3^e) mod f, computed by e-fold cubing of x.
    fn frob_x(&self, e: u32) -> u64 {
        let mut v = self.coeff_x();
        for _ in 0..e {
            v = self.pow(v, 3);
        }
        v
    }
    fn coeff_x(&self) -> u64 {
        if self.n == 1 {
            0 // x = t reduces to... deg-1 field: x is 0-th? handle n=1 specially
        } else {
            3 // digit vector [0,1] = t
        }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        let da = to_digits(a, self.n);
        let db = to_digits(b, self.n);
        let d: Vec<u8> = (0..self.n as usize).map(|i| (da[i] + 3 - db[i]) % 3).collect();
        from_digits(&d)
    }

    /// Whether gcd(g, f) is a unit (degree 0), where g is given as an element (deg < n) — here g
    /// is x^(3^m) - x reduced mod f, so this asks whether that reduced value is non-zero and f is
    /// coprime to it; for the Rabin test it suffices that the reduced g is a unit in the ring,
    /// i.e. invertible mod f.
    fn poly_gcd_is_unit(&self, g: u64) -> bool {
        if g == 0 {
            return false;
        }
        // invertible mod irreducible-candidate f iff non-zero; but f may be reducible here, so do
        // a real gcd over F_3 between f and g (as true polynomials)
        let mut a = self.full_f();
        let mut b = to_digits(g, self.n);
        trim3(&mut b);
        while !b.is_empty() {
            let r = polyrem3(&a, &b);
            a = b;
            b = r;
        }
        a.len() == 1 // gcd degree 0
    }
    fn full_f(&self) -> Vec<u8> {
        let mut f = vec![0u8; self.n as usize + 1];
        f[self.n as usize] = 1;
        for (k, &r) in self.red.iter().enumerate() {
            f[k] = (3 - r) % 3; // f = t^n - red
        }
        trim3(&mut f);
        f
    }
}

fn trim3(v: &mut Vec<u8>) {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    if v.len() == 1 && v[0] == 0 {
        v.clear();
    }
}

/// remainder of a mod b over F_3 (b monic-izable).
fn polyrem3(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut r = a.to_vec();
    let bl = b.len();
    let binv = inv3(b[bl - 1]);
    while r.len() >= bl && !(r.len() == 1 && r[0] == 0) {
        trim3(&mut r);
        if r.len() < bl {
            break;
        }
        let shift = r.len() - bl;
        let c = (r[r.len() - 1] * binv) % 3;
        if c == 0 {
            r.pop();
            continue;
        }
        for i in 0..bl {
            r[shift + i] = (r[shift + i] + 3 - (c * b[i]) % 3) % 3;
        }
        trim3(&mut r);
    }
    trim3(&mut r);
    r
}
fn inv3(x: u8) -> u8 {
    match x {
        1 => 1,
        2 => 2,
        _ => 0,
    }
}
fn distinct_primes(mut n: u32) -> Vec<u32> {
    let mut ps = vec![];
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            ps.push(d);
            while n % d == 0 {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        ps.push(n);
    }
    ps
}

impl Field for GF3n {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        let da = to_digits(a, self.n);
        let db = to_digits(b, self.n);
        let d: Vec<u8> = (0..self.n as usize).map(|i| (da[i] + db[i]) % 3).collect();
        from_digits(&d)
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        GF3n::sub(self, a, b)
    }
    fn neg(&self, a: u64) -> u64 {
        GF3n::sub(self, 0, a)
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.mul_digits(a, b)
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero in GF(3^n)");
        self.pow(a, self.pow3n - 2)
    }
    fn from_u64(&self, n: u64) -> u64 {
        n % 3
    }
    fn char(&self) -> u64 {
        3
    }
    fn size(&self) -> u128 {
        self.pow3n as u128
    }
    fn q(&self) -> Big {
        Big::from_u64(self.pow3n)
    }
    fn random(&self, rng: &mut Rng) -> u64 {
        rng.below(self.pow3n)
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        // q = 3^n is odd; a is a square iff a^((q-1)/2) = 1, then r = a^((q+1)/4) if q = 3 mod 4,
        // else Tonelli-Shanks. q = 3^n: 3^n mod 4 is 3 if n odd, 1 if n even.
        if self.pow(a, (self.pow3n - 1) / 2) != 1 {
            return None;
        }
        if self.pow3n % 4 == 3 {
            let r = self.pow(a, (self.pow3n + 1) / 4);
            return Some(r);
        }
        // Tonelli-Shanks over F_q
        let qm1 = self.pow3n - 1;
        let mut s = 0u32;
        let mut m = qm1;
        while m % 2 == 0 {
            m /= 2;
            s += 1;
        }
        let mut rng = Rng::new(0x3333);
        let z = loop {
            let z = 1 + rng.below(self.pow3n - 1);
            if self.pow(z, qm1 / 2) != 1 {
                break z;
            }
        };
        let mut c = self.pow(z, m);
        let mut t = self.pow(a, m);
        let mut r = self.pow(a, (m + 1) / 2);
        let mut mm = s;
        while t != 1 {
            let mut i = 0u32;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul_digits(tt, tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(mm - i - 1) {
                b = self.mul_digits(b, b);
            }
            mm = i;
            c = self.mul_digits(b, b);
            t = self.mul_digits(t, c);
            r = self.mul_digits(r, b);
        }
        Some(r)
    }
}
