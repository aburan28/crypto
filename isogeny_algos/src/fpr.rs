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
    /// the nonzero entries (k, red[k]) of the reduction tail, where z^r = sum_k red[k] z^k
    /// (low weight: 1-3 terms usually)
    nz: Vec<(usize, u64)>,
    /// true when sums of the products formed in one `mul` fit in u128 without reduction
    lazy: bool,
    /// Frobenius matrix: frob[i] = (z^p)^i, so (sum a_i z^i)^p = sum a_i frob[i]
    frob: Vec<ER>,
}

impl FpR {
    /// F_{p^r} with a low-weight monic irreducible z^r - red(z) found by Rabin's test over F_p.
    pub fn new(p: u64, r: usize) -> FpR {
        assert!(is_prime(p), "FpR: p must be prime");
        assert!((1..=RMAX).contains(&r), "FpR: 1 <= r <= {RMAX}");
        if r == 1 {
            return Self::build(p, r, [0; RMAX]);
        }
        let zp = Zp::new(p);
        // search tails of increasing weight: z^r = c0 + c1 z + ... (few nonzero terms, small coeffs)
        for weight in 1..=3usize {
            if let Some(red) = search_tail(&zp, r, weight) {
                return Self::build(p, r, red);
            }
        }
        // dense fallback (rare): brute small coefficients
        let mut rng = Rng::new(0x1F9C ^ p ^ r as u64);
        for _ in 0..100_000 {
            let mut red = [0u64; RMAX];
            for ri in red.iter_mut().take(r) {
                *ri = rng.below(p);
            }
            if is_irreducible(&zp, r, &red) {
                return Self::build(p, r, red);
            }
        }
        panic!("FpR: no irreducible of degree {r} found for p={p}");
    }

    /// F_{p^r} = F_p[z]/(f) for a monic irreducible f of degree r >= 2 (coefficients low to
    /// high), so that z is a root of f: a root of an F_p-polynomial without root finding in the
    /// extension. Irreducibility is checked (Rabin), since the arithmetic is wrong without it.
    pub fn from_modulus(p: u64, f: &[u64]) -> FpR {
        let r = f.len() - 1;
        assert!((2..=RMAX).contains(&r) && f[r] == 1, "FpR::from_modulus: monic, 2 <= degree <= {RMAX}");
        let mut red = [0u64; RMAX];
        for k in 0..r {
            red[k] = (p - f[k] % p) % p; // z^r = -(f_0 + ... + f_(r-1) z^(r-1))
        }
        assert!(is_irreducible(&Zp::new(p), r, &red), "FpR::from_modulus: modulus is reducible");
        Self::build(p, r, red)
    }

    /// The generator z of F_p[z]/(modulus).
    pub fn gen(&self) -> ER {
        let mut z = [0u64; RMAX];
        z[1] = 1;
        z
    }

    fn build(p: u64, r: usize, red: ER) -> FpR {
        let nz: Vec<(usize, u64)> = (0..r).filter(|&k| red[k] % p != 0).map(|k| (k, red[k] % p)).collect();
        // worst case per accumulator entry: r products (2x for squaring) plus (r-1)*|nz| folds,
        // each < p^2; lazy accumulation is safe when that bound stays below 2^127
        let bits = 64 - p.leading_zeros();
        let terms = (2 * r + r.saturating_sub(1) * nz.len()).max(1) as u32;
        let lazy = 2 * bits + (32 - terms.leading_zeros()) <= 127;
        let mut f = FpR { p, r, nz, lazy, frob: vec![] };
        if r > 1 {
            let mut z = [0u64; RMAX];
            z[1] = 1;
            let zp = f.pow_big(z, &Big::from_u64(p)); // z^p
            let mut fr = vec![f.one()];
            for i in 1..r {
                let next = f.mul(fr[i - 1], zp);
                fr.push(next);
            }
            f.frob = fr;
        }
        f
    }

    /// Fold the coefficients of degree >= r with z^r = red and reduce mod p.
    #[inline]
    fn fold(&self, acc: &mut [u128; 2 * RMAX]) -> ER {
        let r = self.r;
        let pp = self.p as u128;
        for i in (r..2 * r - 1).rev() {
            let ci = acc[i] % pp;
            if ci == 0 {
                continue;
            }
            for &(k, rk) in &self.nz {
                let idx = i - r + k;
                if self.lazy {
                    acc[idx] += ci * rk as u128;
                } else {
                    acc[idx] = (acc[idx] % pp + ci * rk as u128 % pp) % pp;
                }
            }
        }
        let mut out = [0u64; RMAX];
        for i in 0..r {
            out[i] = (acc[i] % pp) as u64;
        }
        out
    }

    /// Frobenius a -> a^p (F_p-linear: a matrix-vector product).
    pub fn frobenius(&self, a: ER) -> ER {
        let r = self.r;
        let pp = self.p as u128;
        let mut acc = [0u128; RMAX];
        for i in 0..r {
            if a[i] == 0 {
                continue;
            }
            let ai = a[i] as u128;
            for k in 0..r {
                if self.lazy {
                    acc[k] += ai * self.frob[i][k] as u128;
                } else {
                    acc[k] = (acc[k] + ai * self.frob[i][k] as u128 % pp) % pp;
                }
            }
        }
        let mut out = [0u64; RMAX];
        for k in 0..r {
            out[k] = (acc[k] % pp) as u64;
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
    // Elements are kept reduced (every coefficient < p), so add/sub need one conditional
    // correction instead of a division.
    fn add(&self, a: ER, b: ER) -> ER {
        let p = self.p as u128;
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            let s = a[i] as u128 + b[i] as u128;
            c[i] = if s >= p { (s - p) as u64 } else { s as u64 };
        }
        c
    }
    fn sub(&self, a: ER, b: ER) -> ER {
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = if a[i] >= b[i] { a[i] - b[i] } else { a[i] + (self.p - b[i]) };
        }
        c
    }
    fn neg(&self, a: ER) -> ER {
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = if a[i] == 0 { 0 } else { self.p - a[i] };
        }
        c
    }
    fn mul(&self, a: ER, b: ER) -> ER {
        let r = self.r;
        if r == 1 {
            let mut c = [0u64; RMAX];
            c[0] = (a[0] as u128 * b[0] as u128 % self.p as u128) as u64;
            return c;
        }
        let pp = self.p as u128;
        let mut acc = [0u128; 2 * RMAX];
        for i in 0..r {
            if a[i] == 0 {
                continue;
            }
            let ai = a[i] as u128;
            for j in 0..r {
                if self.lazy {
                    acc[i + j] += ai * b[j] as u128;
                } else {
                    acc[i + j] = (acc[i + j] + ai * b[j] as u128 % pp) % pp;
                }
            }
        }
        self.fold(&mut acc)
    }
    fn sq(&self, a: ER) -> ER {
        let r = self.r;
        if r == 1 || !self.lazy {
            return self.mul(a, a);
        }
        // a_i a_j for i < j once, doubled; squares on the diagonal
        let mut acc = [0u128; 2 * RMAX];
        for i in 0..r {
            if a[i] == 0 {
                continue;
            }
            let ai = a[i] as u128;
            acc[2 * i] += ai * ai;
            let ai2 = 2 * ai;
            for j in (i + 1)..r {
                acc[i + j] += ai2 * a[j] as u128;
            }
        }
        self.fold(&mut acc)
    }
    /// Inverse through the norm: with b = a^p a^(p^2) ... a^(p^(r-1)) (Frobenius powers, each a
    /// matrix-vector product), N(a) = a b lies in F_p, and a^-1 = b N(a)^-1.
    fn inv(&self, a: ER) -> ER {
        assert!(!self.is_zero(a), "FpR: inverse of zero");
        let p = self.p;
        if self.r == 1 {
            let mut c = [0u64; RMAX];
            c[0] = crate::fpm::inv_u64(a[0], p);
            return c;
        }
        let mut f = self.frobenius(a);
        let mut b = f;
        for _ in 2..self.r {
            f = self.frobenius(f);
            b = self.mul(b, f);
        }
        let n = self.mul(a, b);
        debug_assert!(n[1..self.r].iter().all(|&c| c == 0), "norm must lie in F_p");
        let ninv = crate::fpm::inv_u64(n[0], p) as u128;
        let mut c = [0u64; RMAX];
        for i in 0..self.r {
            c[i] = (b[i] as u128 * ninv % p as u128) as u64;
        }
        c
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
