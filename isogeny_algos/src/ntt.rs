//! Exact integer convolution of multi-limb vectors by number-theoretic transforms over several
//! 62-bit primes p_i = c 2^26 + 1, recombined by Garner's mixed-radix CRT.
//!
//! The multi-limb prime fields use it for long products (`FpM::conv` above its threshold): the
//! Montgomery residues are convolved as integers, and each output coefficient is returned as its
//! mixed-radix digits d_0, d_1, .. with value d_0 + d_1 m_0 + d_2 m_0 m_1 + .., which the field
//! folds back with one wide Montgomery reduction. The result equals the lazily accumulated
//! schoolbook sum exactly; only the cost changes, from O(n^1.58) field multiplications (Karatsuba)
//! to O(n log n) word multiplications per prime.
use std::sync::OnceLock;

/// Largest transform: 2^26 points (products up to 2^26 coefficients).
pub const MAX_LOG: u32 = 26;
/// Number of primes available (enough for 1024-bit fields at 2^26 points: 2 * 1024 + 26 bits).
const NPRIMES: usize = 36;

/// One NTT prime p in (2^61, 2^62) with 2^MAX_LOG | p - 1, in Montgomery form (R = 2^64).
#[derive(Clone, Copy, Debug)]
pub struct Prime {
    pub p: u64,
    /// -p^{-1} mod 2^64
    pinv: u64,
    /// R^2 mod p
    r2: u64,
    /// Montgomery form of a primitive 2^MAX_LOG-th root of unity
    root: u64,
}

impl Prime {
    fn new(p: u64) -> Self {
        let mut inv = 1u64; // p^{-1} mod 2^64 by Newton
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let r1 = ((1u128 << 64) % p as u128) as u64;
        let r2 = ((r1 as u128 * r1 as u128) % p as u128) as u64;
        let mut pr = Prime {
            p,
            pinv: inv.wrapping_neg(),
            r2,
            root: 0,
        };
        let c = (p - 1) >> MAX_LOG;
        for x in 2u64.. {
            let w = pr.pow(pr.mont(x), c);
            // order exactly 2^MAX_LOG: w^(2^(MAX_LOG-1)) = -1
            if pr.pow(w, 1u64 << (MAX_LOG - 1)) != pr.mont(1) {
                pr.root = w;
                break;
            }
        }
        pr
    }
    /// a b R^{-1} mod p, for a b < p 2^64 (in particular a < 2^64, b < p).
    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let t = a as u128 * b as u128;
        let m = (t as u64).wrapping_mul(self.pinv);
        let u = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if u >= self.p {
            u - self.p
        } else {
            u
        }
    }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    fn mont(&self, a: u64) -> u64 {
        self.mul(a % self.p, self.r2)
    }
    fn pow(&self, mut b: u64, mut e: u64) -> u64 {
        let mut r = self.mont(1);
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
    /// Montgomery form of a^{-1} for a not divisible by p.
    fn inv_mont(&self, a: u64) -> u64 {
        self.pow(self.mont(a), self.p - 2)
    }
}

/// The NTT primes, largest first.
pub fn primes() -> &'static [Prime] {
    static P: OnceLock<Vec<Prime>> = OnceLock::new();
    P.get_or_init(|| {
        let mut v = vec![];
        let mut c = ((1u64 << 62) - 1) >> MAX_LOG;
        while v.len() < NPRIMES {
            let p = (c << MAX_LOG) + 1;
            if p < (1u64 << 62) && crate::field::is_prime(p) {
                v.push(Prime::new(p));
            }
            c -= 1;
        }
        assert!(v.last().unwrap().p > 1u64 << 61);
        v
    })
}

/// Shoup's product: w x mod p up to one p, i.e. in [0, 2p), for any x < 2^64, given
/// wp = floor(w 2^64 / p) and w < p.
#[inline(always)]
fn mul_shoup(x: u64, w: u64, wp: u64, p: u64) -> u64 {
    let q = ((x as u128 * wp as u128) >> 64) as u64;
    w.wrapping_mul(x).wrapping_sub(q.wrapping_mul(p))
}

/// Twiddle tables of one prime for transforms up to 2^logn points: w^k and its Shoup companion
/// for k < 2^(logn-1), w a primitive 2^logn-th root (plain residues), and the same for w^{-1}.
/// A transform of 2^m points reads every 2^(logn-m)-th entry.
struct Tw {
    logn: u32,
    fw: Vec<u64>,
    fwp: Vec<u64>,
    iw: Vec<u64>,
    iwp: Vec<u64>,
}

impl Tw {
    fn new(pr: &Prime, logn: u32) -> Tw {
        let half = 1usize << (logn - 1);
        let wm = pr.pow(pr.root, 1u64 << (MAX_LOG - logn));
        let plain = |x: u64| pr.mul(x, 1);
        let w = plain(wm);
        let wi = plain(pr.pow(wm, (1u64 << logn) - 1));
        let gen = |g: u64| -> (Vec<u64>, Vec<u64>) {
            let mut v = Vec::with_capacity(half);
            let mut vp = Vec::with_capacity(half);
            let mut x = 1u64;
            for _ in 0..half {
                v.push(x);
                vp.push((((x as u128) << 64) / pr.p as u128) as u64);
                x = ((x as u128 * g as u128) % pr.p as u128) as u64;
            }
            (v, vp)
        };
        let (fw, fwp) = gen(w);
        let (iw, iwp) = gen(wi);
        Tw {
            logn,
            fw,
            fwp,
            iw,
            iwp,
        }
    }
}

thread_local! {
    /// Twiddle tables per prime index, grown to the largest transform used so far.
    static TW: std::cell::RefCell<Vec<Option<std::rc::Rc<Tw>>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

fn tw(i: usize, pr: &Prime, logn: u32) -> std::rc::Rc<Tw> {
    TW.with(|c| {
        let mut v = c.borrow_mut();
        if v.len() <= i {
            v.resize(i + 1, None);
        }
        match &v[i] {
            Some(t) if t.logn >= logn => t.clone(),
            _ => {
                let t = std::rc::Rc::new(Tw::new(pr, logn.max(10)));
                v[i] = Some(t.clone());
                t
            }
        }
    })
}

/// Gentleman–Sande (decimation in frequency), Harvey's lazy form: input and output in [0, 2p),
/// natural order in, bit-reversed order out.
fn dif(a: &mut [u64], t: &Tw, p: u64) {
    let n = a.len();
    let p2 = 2 * p;
    let base = (1usize << t.logn) / n;
    let mut len = n;
    while len >= 2 {
        let half = len / 2;
        let stride = base * (n / len);
        for blk in a.chunks_exact_mut(len) {
            let (lo, hi) = blk.split_at_mut(half);
            for (j, (x, y)) in lo.iter_mut().zip(hi.iter_mut()).enumerate() {
                let (u, v) = (*x, *y);
                let s = u + v;
                *x = if s >= p2 { s - p2 } else { s };
                let k = j * stride;
                *y = mul_shoup(u + p2 - v, t.fw[k], t.fwp[k], p);
            }
        }
        len = half;
    }
}

/// Cooley–Tukey (decimation in time), Harvey's lazy form: input in [0, 2p), output in [0, 4p),
/// bit-reversed order in, natural order out.
fn dit(a: &mut [u64], t: &Tw, p: u64) {
    let n = a.len();
    let p2 = 2 * p;
    let base = (1usize << t.logn) / n;
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let stride = base * (n / len);
        for blk in a.chunks_exact_mut(len) {
            let (lo, hi) = blk.split_at_mut(half);
            for (j, (x, y)) in lo.iter_mut().zip(hi.iter_mut()).enumerate() {
                let u = if *x >= p2 { *x - p2 } else { *x };
                let k = j * stride;
                let v = mul_shoup(*y, t.iw[k], t.iwp[k], p);
                *x = u + v;
                *y = u + p2 - v;
            }
        }
        len *= 2;
    }
}

/// x mod p for an N-limb integer x: sum of limb_k 2^(64 k) mod p, one Montgomery product per
/// limb against R^(k+1) mod p.
#[inline(always)]
fn reduce_limbs<const N: usize>(pr: &Prime, rk: &[u64; N], x: &[u64; N]) -> u64 {
    let mut s = 0u64;
    for k in 0..N {
        s = pr.add(s, pr.mul(x[k], rk[k]));
    }
    s
}

/// Number of primes whose product exceeds every coefficient of the integer convolution of
/// vectors with entries below 2^entry_bits and shorter length `min_len`.
pub fn primes_needed(entry_bits: u32, min_len: usize) -> usize {
    let bits = 2 * entry_bits as usize + (usize::BITS - min_len.leading_zeros()) as usize + 1;
    bits.div_ceil(61)
}

/// Garner constants m_j^{-1} mod m_i (Montgomery form for prime i), at [j * NPRIMES + i].
fn garner() -> &'static [u64] {
    static G: OnceLock<Vec<u64>> = OnceLock::new();
    G.get_or_init(|| {
        let ps = primes();
        let mut c = vec![0u64; NPRIMES * NPRIMES];
        for i in 0..NPRIMES {
            for j in 0..i {
                c[j * NPRIMES + i] = ps[i].inv_mont(ps[j].p % ps[i].p);
            }
        }
        c
    })
}

/// The exact integer convolution of `a` and `b` (entries read as N-limb integers below
/// 2^entry_bits), as mixed-radix digits: `out[k * t + i]` is digit d_i of coefficient k, with
/// value sum_i d_i prod_{j < i} m_j, m_j = `primes()[j].p`, t = `primes_needed(entry_bits,
/// min(len))`. A squaring (`a` and `b` the same slice) transforms once. Returns (digits, t).
pub fn conv_digits<const N: usize>(
    a: &[[u64; N]],
    b: &[[u64; N]],
    entry_bits: u32,
) -> (Vec<u64>, usize) {
    assert!(!a.is_empty() && !b.is_empty());
    let square = std::ptr::eq(a, b);
    let nout = a.len() + b.len() - 1;
    let logn = (usize::BITS - (nout - 1).leading_zeros()).max(1);
    assert!(logn <= MAX_LOG, "NTT product too long");
    let n = 1usize << logn;
    let t = primes_needed(entry_bits, a.len().min(b.len()));
    assert!(t <= NPRIMES, "NTT: too few primes for this field size");
    let ps = &primes()[..t];
    let mut res: Vec<Vec<u64>> = Vec::with_capacity(t);
    let mut fb = if square { vec![] } else { vec![0u64; n] };
    for (pi, pr) in ps.iter().enumerate() {
        let tw = tw(pi, pr, logn);
        let mut rk = [0u64; N];
        let r1 = ((1u128 << 64) % pr.p as u128) as u64;
        let mut x = r1;
        for slot in rk.iter_mut() {
            *slot = x; // R^(k+1) mod p
            x = ((x as u128 * r1 as u128) % pr.p as u128) as u64;
        }
        let mut fa = vec![0u64; n];
        for (d, s) in fa.iter_mut().zip(a) {
            *d = reduce_limbs(pr, &rk, s);
        }
        dif(&mut fa, &tw, pr.p);
        // pointwise: mul(mul(A, B), n^{-1} R^2) = A B / n (A, B < 2p: A B < p 2^64), so the
        // unnormalised inverse transform returns the convolution itself
        let scale = pr.mul(pr.inv_mont(n as u64), pr.r2); // (n^{-1} R) R^2 R^{-1}
        if square {
            for x in fa.iter_mut() {
                *x = pr.mul(pr.mul(*x, *x), scale);
            }
        } else {
            fb.iter_mut().for_each(|x| *x = 0);
            for (d, s) in fb.iter_mut().zip(b) {
                *d = reduce_limbs(pr, &rk, s);
            }
            dif(&mut fb, &tw, pr.p);
            for (x, y) in fa.iter_mut().zip(&fb) {
                *x = pr.mul(pr.mul(*x, *y), scale);
            }
        }
        dit(&mut fa, &tw, pr.p);
        fa.truncate(nout);
        let p2 = 2 * pr.p;
        for x in fa.iter_mut() {
            if *x >= p2 {
                *x -= p2;
            }
            if *x >= pr.p {
                *x -= pr.p;
            }
        }
        res.push(fa);
    }
    // Garner: d_0 = r_0; d_i = (..((r_i - d_0) c_0i - d_1) c_1i .. - d_{i-1}) c_{i-1,i} mod m_i
    let cinv = garner();
    let mut out = vec![0u64; nout * t];
    for k in 0..nout {
        let d = &mut out[k * t..(k + 1) * t];
        for i in 0..t {
            let pr = &ps[i];
            let mut x = res[i][k];
            for j in 0..i {
                let dj = if d[j] >= pr.p { d[j] - pr.p } else { d[j] }; // d_j < 2^62 < 2 m_i
                x = pr.mul(pr.sub(x, dj), cinv[j * NPRIMES + i]);
            }
            d[i] = x;
        }
    }
    (out, t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primes_have_roots() {
        for pr in primes() {
            assert_eq!((pr.p - 1) % (1u64 << MAX_LOG), 0);
            let one = pr.mont(1);
            assert_eq!(pr.pow(pr.root, 1u64 << MAX_LOG), one);
            assert_ne!(pr.pow(pr.root, 1u64 << (MAX_LOG - 1)), one);
        }
    }

    #[test]
    fn single_limb_against_u128() {
        // small entries: the digits recombine to the exact u128 convolution
        let a: Vec<[u64; 1]> = (0..37u64).map(|i| [i * 1_000_003 + 7]).collect();
        let b: Vec<[u64; 1]> = (0..50u64).map(|i| [i * i + 1]).collect();
        let (d, t) = conv_digits(&a, &b, 40);
        let ps = primes();
        for k in 0..a.len() + b.len() - 1 {
            let mut want = 0u128;
            for i in 0..a.len() {
                if k >= i && k - i < b.len() {
                    want += a[i][0] as u128 * b[k - i][0] as u128;
                }
            }
            let mut got = 0u128;
            let mut m = 1u128;
            for i in 0..t {
                got += d[k * t + i] as u128 * m;
                if i + 1 < t {
                    m *= ps[i].p as u128;
                }
            }
            assert_eq!(got, want, "coefficient {k}");
        }
    }
}
