//! Degree-r extension field F_{p^r} = F_p[z]/(m(z)) for a general prime p and a runtime degree
//! r <= RMAX, with a low-weight monic irreducible m found by Rabin's test. Elements are a packed
//! coefficient array [u64; RMAX] (only the first r entries are used), so the type is `Copy` and
//! the generic `poly`/SEA machinery runs over it unchanged. Used by the Atkin-prime trace
//! recovery, where a curve E/F_p that has no rational l-isogeny gains one over F_{p^r} (r = the
//! Atkin order): Frobenius^r then acts on E[l] as a scalar, and the ordinary Elkies eigenvalue
//! computation over F_{p^r} returns that scalar lambda^r mod l.
use crate::bigint::Big;
use crate::field::{is_prime, Field, Rng, Zp};
use crate::poly;

/// Maximum extension degree (Atkin orders above this are skipped by the caller).
pub const RMAX: usize = 16;
pub type ER = [u64; RMAX];

#[derive(Clone, Debug)]
pub struct FpR {
    pub p: u64,
    pub r: usize,
    /// z^r = sum_i red[i] z^i (i < r).
    red: ER,
}

impl FpR {
    /// F_{p^r} with a low-weight monic irreducible z^r - red(z) found by Rabin's test over F_p.
    pub fn new(p: u64, r: usize) -> FpR {
        assert!(is_prime(p), "FpR: p must be prime");
        assert!((1..=RMAX).contains(&r), "FpR: 1 <= r <= {RMAX}");
        if r == 1 {
            return FpR { p, r, red: [0; RMAX] };
        }
        let zp = Zp::new(p);
        // search tails of increasing weight: z^r = c0 + c1 z + ... (few nonzero terms, small coeffs)
        for weight in 1..=3usize {
            if let Some(red) = search_tail(&zp, r, weight) {
                return FpR { p, r, red };
            }
        }
        // dense fallback (rare): brute small coefficients
        let zp = Zp::new(p);
        let mut rng = Rng::new(0x1F9C ^ p ^ r as u64);
        for _ in 0..100_000 {
            let mut red = [0u64; RMAX];
            for ri in red.iter_mut().take(r) {
                *ri = rng.below(p);
            }
            if is_irreducible(&zp, r, &red) {
                return FpR { p, r, red };
            }
        }
        panic!("FpR: no irreducible of degree {r} found for p={p}");
    }

    fn reduce_vec(&self, mut c: Vec<u64>) -> ER {
        let p = self.p;
        let r = self.r;
        // c has degree up to 2r-2; fold top coefficients with z^r = red
        for i in (r..c.len()).rev() {
            let ci = c[i] % p;
            if ci == 0 {
                continue;
            }
            c[i] = 0;
            for k in 0..r {
                if self.red[k] != 0 {
                    let idx = i - r + k;
                    c[idx] = (c[idx] + (ci as u128 * self.red[k] as u128 % p as u128) as u64) % p;
                }
            }
        }
        let mut out = [0u64; RMAX];
        for (i, slot) in out.iter_mut().take(r).enumerate() {
            *slot = c.get(i).copied().unwrap_or(0) % p;
        }
        out
    }
}

/// Coefficients of the monic reduction polynomial m(z) = z^r - red as a full F_p coeff vector.
fn modulus_poly(zp: &Zp, r: usize, red: &ER) -> Vec<u64> {
    let p = zp.p;
    let mut m = vec![0u64; r + 1];
    m[r] = 1;
    for k in 0..r {
        m[k] = (p - red[k] % p) % p; // -red[k]
    }
    m
}

fn is_irreducible(zp: &Zp, r: usize, red: &ER) -> bool {
    // Rabin: X^(p^r) == X mod m, and for each prime q | r, gcd(X^(p^(r/q)) - X, m) == 1.
    let m = modulus_poly(zp, r, red);
    let p = zp.p;
    let x = poly::x_poly(zp);
    let frob = |e: u32| {
        // X^(p^e) mod m
        let mut v = x.clone();
        for _ in 0..e {
            v = poly::powmod_big(zp, &v, &Big::from_u64(p), &m);
        }
        v
    };
    // X^(p^r) must reduce to X mod m
    let top = poly::sub(zp, &frob(r as u32), &x);
    if !is_zero_poly(zp, &top) {
        return false;
    }
    for q in distinct_primes(r as u64) {
        let e = (r as u64 / q) as u32;
        let d = poly::sub(zp, &frob(e), &x);
        let g = poly::gcd(zp, &m, &d);
        if poly::deg(zp, &g) > 0 {
            return false;
        }
    }
    true
}

fn is_zero_poly(zp: &Zp, a: &[u64]) -> bool {
    a.iter().all(|&c| c % zp.p == 0)
}

fn search_tail(zp: &Zp, r: usize, weight: usize) -> Option<ER> {
    // enumerate reduction tails with exactly `weight` nonzero coefficients, small values
    let p = zp.p;
    let smalls: Vec<u64> = (1..p.min(7)).collect();
    // choose `weight` positions among 0..r and assign small nonzero values
    let mut red = [0u64; RMAX];
    fn rec(zp: &Zp, r: usize, pos: usize, left: usize, start: usize, smalls: &[u64], red: &mut ER) -> Option<ER> {
        if left == 0 {
            if is_irreducible(zp, r, red) {
                return Some(*red);
            }
            return None;
        }
        for position in start..r {
            for &v in smalls {
                red[position] = v;
                if let Some(ok) = rec(zp, r, pos + 1, left - 1, position + 1, smalls, red) {
                    return Some(ok);
                }
                red[position] = 0;
            }
        }
        None
    }
    rec(zp, r, 0, weight, 0, &smalls, &mut red)
}

fn distinct_primes(mut n: u64) -> Vec<u64> {
    let mut ps = vec![];
    let mut d = 2u64;
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

impl Field for FpR {
    type E = ER;
    fn zero(&self) -> ER {
        [0; RMAX]
    }
    fn one(&self) -> ER {
        let mut e = [0; RMAX];
        e[0] = 1 % self.p;
        e
    }
    fn add(&self, a: ER, b: ER) -> ER {
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = (a[i] + b[i]) % self.p;
        }
        c
    }
    fn sub(&self, a: ER, b: ER) -> ER {
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = (a[i] + self.p - b[i] % self.p) % self.p;
        }
        c
    }
    fn neg(&self, a: ER) -> ER {
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = (self.p - a[i] % self.p) % self.p;
        }
        c
    }
    fn mul(&self, a: ER, b: ER) -> ER {
        let r = self.r;
        let p = self.p as u128;
        let mut prod = vec![0u64; 2 * r];
        for i in 0..r {
            if a[i] == 0 {
                continue;
            }
            let ai = a[i] as u128;
            for j in 0..r {
                if b[j] != 0 {
                    prod[i + j] = ((prod[i + j] as u128 + ai * b[j] as u128) % p) as u64;
                }
            }
        }
        self.reduce_vec(prod)
    }
    fn inv(&self, a: ER) -> ER {
        assert!(!self.is_zero(a), "FpR: inverse of zero");
        // extended Euclid in F_p[z] between m and a
        let zp = Zp::new(self.p);
        let m = modulus_poly(&zp, self.r, &self.red);
        let mut av: Vec<u64> = a[..self.r].iter().map(|&c| c % self.p).collect();
        while av.len() > 1 && *av.last().unwrap() == 0 {
            av.pop();
        }
        let inv = poly::invmod(&zp, &av, &m).expect("FpR inverse exists for nonzero element");
        let mut out = [0u64; RMAX];
        for (i, slot) in out.iter_mut().take(self.r).enumerate() {
            *slot = inv.get(i).copied().unwrap_or(0) % self.p;
        }
        out
    }
    fn from_u64(&self, n: u64) -> ER {
        let mut e = [0; RMAX];
        e[0] = n % self.p;
        e
    }
    fn char(&self) -> u64 {
        self.p
    }
    fn size(&self) -> u128 {
        let mut s = 1u128;
        for _ in 0..self.r {
            s = s.saturating_mul(self.p as u128);
        }
        s
    }
    fn q(&self) -> Big {
        let mut s = Big::from_u64(1);
        let bp = Big::from_u64(self.p);
        for _ in 0..self.r {
            s = s.mul(&bp);
        }
        s
    }
    fn random(&self, rng: &mut Rng) -> ER {
        let mut e = [0u64; RMAX];
        for i in 0..self.r {
            e[i] = rng.below(self.p);
        }
        e
    }
}

impl FpR {
    /// Embed a base-field element a in F_p as the constant a in F_{p^r}.
    pub fn embed(&self, a: u64) -> ER {
        self.from_u64(a)
    }
    /// If the element lies in the prime field (coeffs 1..r all zero), return it.
    pub fn project(&self, a: ER) -> Option<u64> {
        if a[1..self.r].iter().all(|&c| c % self.p == 0) {
            Some(a[0] % self.p)
        } else {
            None
        }
    }
}
