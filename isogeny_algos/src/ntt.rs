//! Exact integer convolution of multi-limb vectors by number-theoretic transforms over several
//! word-size primes p_i = c 2^26 + 1, recombined by the explicit CRT.
//!
//! The multi-limb prime fields use it for long products (`FpM::conv` above its threshold): the
//! Montgomery residues are convolved as integers, and each output coefficient V is returned as
//! its explicit-CRT residues u_i, with V = sum_i u_i (M / m_i) - q M and q read off a
//! floating-point sum, which the field folds back mod p with one wide Montgomery reduction. The
//! result equals the lazily accumulated schoolbook sum exactly; only the cost changes, from
//! O(n^1.58) field multiplications (Karatsuba) to O(n log n) word multiplications per prime.
//!
//! Two prime families: 62-bit primes with portable 64-bit Shoup butterflies, and, on x86-64
//! CPUs with AVX-512 IFMA (detected at run time), 50-bit primes whose transforms run eight
//! lanes at a time on the 52-bit multiply-add instructions. The vector path needs a few more
//! primes (50 instead of 62 bits each) and does the same work per prime in about a quarter of
//! the time. [`SIMD`] switches it off per thread (tests), the environment variable
//! `ISOGENY_NO_SIMD` for the whole process (A/B runs).
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::OnceLock;

/// Largest transform: 2^26 points (products up to 2^26 coefficients).
pub const MAX_LOG: u32 = 26;

thread_local! {
    /// Whether products may use the AVX-512 IFMA transforms when the CPU has them (switched off
    /// by tests and A/B runs to exercise the portable path).
    pub static SIMD: Cell<bool> = const { Cell::new(true) };
}

/// Whether products on this thread go through the AVX-512 IFMA transforms (the CPU has them
/// and [`SIMD`] is on).
pub fn simd_active() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        SIMD.with(|s| s.get()) && ifma::available()
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// One NTT prime p with 2^MAX_LOG | p - 1, in Montgomery form (R = 2^64).
#[derive(Clone, Copy, Debug)]
pub struct Prime {
    pub p: u64,
    /// -p^{-1} mod 2^64
    pinv: u64,
    /// R^2 mod p
    r2: u64,
    /// Montgomery form of a primitive 2^MAX_LOG-th root of unity
    root: u64,
    /// division-free reduction of u128 sums
    red: crate::field::Red128,
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
            red: crate::field::Red128::new(p),
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
    /// R^(k+1) mod p for k < N: reduces limb k of an N-limb integer (`reduce_limbs`).
    fn limb_weights<const N: usize>(&self) -> [u64; N] {
        let mut rk = [0u64; N];
        let r1 = ((1u128 << 64) % self.p as u128) as u64;
        let mut x = r1;
        for slot in rk.iter_mut() {
            *slot = x;
            x = ((x as u128 * r1 as u128) % self.p as u128) as u64;
        }
        rk
    }
}

/// A family of NTT primes p_i = c 2^26 + 1 in (2^bits, 2^(bits+1)), largest first, with the
/// constants of the explicit CRT.
pub struct PrimeSet {
    pub primes: Vec<Prime>,
    /// every prime exceeds 2^bits
    bits: u32,
    /// 1 / p_i as f64, for the quotient
    inv: Vec<f64>,
    /// at [t * len + i] (i < t): the Montgomery form, for prime i, of (M_t / p_i)^{-1} mod p_i,
    /// with M_t = p_0 .. p_{t-1}
    crt: Vec<u64>,
}

impl PrimeSet {
    fn new(bits: u32, count: usize) -> Self {
        let mut primes = vec![];
        let mut c = ((1u64 << (bits + 1)) - 1) >> MAX_LOG;
        while primes.len() < count {
            let p = (c << MAX_LOG) + 1;
            if p < (1u64 << (bits + 1)) && crate::field::is_prime(p) {
                primes.push(Prime::new(p));
            }
            c -= 1;
        }
        assert!(primes.last().unwrap().p > 1u64 << bits);
        let inv = primes.iter().map(|pr| 1.0 / pr.p as f64).collect();
        let mut crt = vec![0u64; (count + 1) * count];
        for t in 1..=count {
            for (i, pr) in primes[..t].iter().enumerate() {
                let mut q = 1u64; // M_t / p_i mod p_i
                for (j, pj) in primes[..t].iter().enumerate() {
                    if j != i {
                        q = ((q as u128 * pj.p as u128) % pr.p as u128) as u64;
                    }
                }
                crt[t * count + i] = pr.inv_mont(q);
            }
        }
        PrimeSet {
            primes,
            bits,
            inv,
            crt,
        }
    }
    /// Number of primes whose product M exceeds twice every coefficient of the integer
    /// convolution when entries are below 2^entry_bits and each coefficient sums at most
    /// `terms` products (the spare bit keeps V < M / 2 for [`PrimeSet::quotient`]).
    pub fn needed(&self, entry_bits: u32, terms: usize) -> usize {
        let bits = 2 * entry_bits as usize + (usize::BITS - terms.leading_zeros()) as usize + 1;
        bits.div_ceil(self.bits as usize)
    }
    /// The quotient q of the explicit CRT for one coefficient with residues `u` (u_i < p_i):
    /// V = sum_i u_i (M / p_i) - q M, M the product of the `u.len()` first primes. Exact
    /// because V < M / 2: the sum of u_i / p_i is q + V / M, in [q, q + 1/2), and its f64 value
    /// is off by less than t 2^-50, so floor(sum + 1/4) = q. Always q < t.
    #[inline(always)]
    pub fn quotient(&self, u: &[u64]) -> u64 {
        let mut s = 0.25f64;
        for (x, i) in u.iter().zip(&self.inv) {
            s += *x as f64 * i;
        }
        s as u64
    }
    /// The pointwise constant n^{-1} (M_t / p_i)^{-1} of prime i, Montgomery form (n a power
    /// of two dividing p - 1, so n^{-1} = p - (p - 1) / n).
    fn scale(&self, i: usize, t: usize, n: usize) -> u64 {
        let pr = &self.primes[i];
        let ninv = pr.mont(pr.p - (pr.p - 1) / n as u64);
        pr.mul(ninv, self.crt[t * self.primes.len() + i])
    }
}

/// The portable family: 62-bit primes (36, enough for 1024-bit fields at 2^26 points).
pub fn scalar_set() -> &'static PrimeSet {
    static S: OnceLock<PrimeSet> = OnceLock::new();
    S.get_or_init(|| PrimeSet::new(61, 36))
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
    /// Twiddle tables per prime, grown to the largest transform used so far.
    static TW: RefCell<HashMap<u64, Rc<Tw>>> = RefCell::new(HashMap::new());
}

fn tw(pr: &Prime, logn: u32) -> Rc<Tw> {
    TW.with(|c| {
        let mut m = c.borrow_mut();
        match m.get(&pr.p) {
            Some(t) if t.logn >= logn => t.clone(),
            _ => {
                let t = Rc::new(Tw::new(pr, logn.max(10)));
                m.insert(pr.p, t.clone());
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

/// x mod p for an N-limb integer x, up to 2p (N <= 8; the lazy transforms take inputs below
/// 2p): sum of limb_k R^(k+1) mod p accumulated in a u128 four limbs at a time (each product
/// below 2^126), one Montgomery reduction (x 2^-64) per group.
#[inline(always)]
fn reduce_limbs<const N: usize>(pr: &Prime, rk: &[u64; N], x: &[u64; N]) -> u64 {
    let mut s = 0u64;
    for (xs, rs) in x.chunks(4).zip(rk.chunks(4)) {
        let mut acc: u128 = 0;
        for (a, r) in xs.iter().zip(rs) {
            acc += *a as u128 * *r as u128;
        }
        s += pr.red.redc(acc); // below p < 2^62 each: at most four groups fit in a u64
    }
    if N > 8 {
        s = pr.red.mod_word(s);
    }
    s
}

/// The transform input: entries reduced mod p (below 2p), folded mod x^len - 1.
fn load<const N: usize>(f: &mut [u64], a: &[[u64; N]], pr: &Prime, rk: &[u64; N]) {
    f.iter_mut().for_each(|x| *x = 0);
    let (n, p2) = (f.len(), 2 * pr.p);
    for (j, s) in a.iter().enumerate() {
        let d = &mut f[j & (n - 1)];
        *d += reduce_limbs(pr, rk, s);
        if *d >= p2 {
            *d -= p2;
        }
    }
}

/// The second operand of a product: the first one again (a squaring, transformed once), a
/// plain vector, or one transformed beforehand ([`prepare`]).
enum Operand<'a, const N: usize> {
    Square,
    Plain(&'a [[u64; N]]),
    Prepared(&'a [Vec<u64>]),
}

/// The prime family for a product with `terms` products per coefficient: the IFMA one when it
/// is in use and has enough primes, otherwise the portable one. Returns (set, t, simd).
fn choose(entry_bits: u32, terms: usize, logn: u32) -> (&'static PrimeSet, usize, bool) {
    #[cfg(target_arch = "x86_64")]
    if logn >= 4 && simd_active() {
        let set = ifma::set();
        let t = set.needed(entry_bits, terms);
        if t <= set.primes.len() {
            return (set, t, true);
        }
    }
    let _ = logn;
    let set = scalar_set();
    let t = set.needed(entry_bits, terms);
    assert!(
        t <= set.primes.len(),
        "NTT: too few primes for this field size"
    );
    (set, t, false)
}

/// Transform length and coefficient bounds of a product of lengths la, lb, mod x^cyc - 1 when
/// `cyc` is below the product length: (n, nout, terms).
fn shape(la: usize, lb: usize, cyc: usize) -> (usize, usize, usize) {
    let full = la + lb - 1;
    let wrapped = cyc < full;
    // a cyclic coefficient sums at most twice as many products
    let terms = la.min(lb) << wrapped as u32;
    let nout = full.min(cyc);
    let logn = (usize::BITS - (nout - 1).leading_zeros()).max(1);
    assert!(logn <= MAX_LOG, "NTT product too long");
    let n = 1usize << logn;
    assert!(!wrapped || n == cyc, "cyclic length must be a power of two");
    (n, nout, terms)
}

/// The exact integer convolution V of `a` and `b` (entries read as N-limb integers below
/// 2^entry_bits) in explicit-CRT form over the returned prime set: `out[k * t + i]` is
/// u_i = V_k (M / p_i)^{-1} mod p_i for coefficient k, so that V_k = sum_i u_i (M / p_i) - q M
/// with q = [`PrimeSet::quotient`] of the row, M = p_0 .. p_{t-1}, p_i = `set.primes[i].p`.
/// The factor (M / p_i)^{-1} rides on the pointwise scaling, so it costs nothing. A squaring
/// (`a` and `b` the same slice) transforms once. Returns (residues, t, set).
///
/// With `cyc` (a power of two) below the product length, the transforms have `cyc` points and
/// the result is the cyclic convolution, mod x^cyc - 1: coefficient k < cyc is V_k + V_(k+cyc)
/// (the caller recovers the wrapped top from a short product).
pub fn conv_crt<const N: usize>(
    a: &[[u64; N]],
    b: &[[u64; N]],
    entry_bits: u32,
    cyc: usize,
) -> (Vec<u64>, usize, &'static PrimeSet) {
    assert!(!a.is_empty() && !b.is_empty());
    let (n, nout, terms) = shape(a.len(), b.len(), cyc);
    let (set, t, simd) = choose(entry_bits, terms, n.trailing_zeros());
    let op = if std::ptr::eq(a, b) {
        Operand::Square
    } else {
        Operand::Plain(b)
    };
    (run(a, op, n, nout, set, t, simd), t, set)
}

/// A fixed operand transformed once, for repeated products with operands of at most
/// `max_len` entries ([`conv_crt_prepared`]): the forward transform per prime.
pub struct Prepared {
    n: usize,
    blen: usize,
    max_len: usize,
    cyc: usize,
    t: usize,
    set: &'static PrimeSet,
    simd: bool,
    tr: Vec<Vec<u64>>,
}

/// Transforms `b` for products with operands of at most `max_len` entries, mod x^cyc - 1 when
/// `cyc` is below the product length (`usize::MAX`: plain products).
pub fn prepare<const N: usize>(
    b: &[[u64; N]],
    entry_bits: u32,
    max_len: usize,
    cyc: usize,
) -> Prepared {
    assert!(!b.is_empty() && max_len > 0);
    let (n, _, terms) = shape(max_len, b.len(), cyc);
    let logn = n.trailing_zeros();
    let (set, t, simd) = choose(entry_bits, terms, logn);
    #[cfg(target_arch = "x86_64")]
    let dg = simd.then(|| ifma::digits52(b));
    let tr = (0..t)
        .map(|pi| {
            let pr = &set.primes[pi];
            let mut f = vec![0u64; n];
            #[cfg(target_arch = "x86_64")]
            if let Some(dg) = &dg {
                // SAFETY: `choose` picks the IFMA family only when the CPU has it
                unsafe { ifma::forward(&mut f, dg, pr) };
                return f;
            }
            load(&mut f, b, pr, &pr.limb_weights::<N>());
            dif(&mut f, &tw(pr, logn), pr.p);
            f
        })
        .collect();
    Prepared {
        n,
        blen: b.len(),
        max_len,
        cyc,
        t,
        set,
        simd,
        tr,
    }
}

/// [`conv_crt`] of `a` and a prepared operand (same output form), or None when `a` is longer
/// than the preparation allows or the IFMA path it was prepared for is off on this thread.
pub fn conv_crt_prepared<const N: usize>(
    a: &[[u64; N]],
    b: &Prepared,
) -> Option<(Vec<u64>, usize, &'static PrimeSet)> {
    if a.is_empty() || a.len() > b.max_len || (b.simd && !simd_active()) {
        return None;
    }
    let (n, nout, _) = shape(a.len(), b.blen, b.cyc);
    if n != b.n {
        return None;
    }
    let out = run(a, Operand::Prepared(&b.tr), n, nout, b.set, b.t, b.simd);
    Some((out, b.t, b.set))
}

/// The product loop over the primes of `set` (portable or IFMA transforms).
fn run<const N: usize>(
    a: &[[u64; N]],
    b: Operand<N>,
    n: usize,
    nout: usize,
    set: &PrimeSet,
    t: usize,
    simd: bool,
) -> Vec<u64> {
    #[cfg(target_arch = "x86_64")]
    if simd {
        // SAFETY: `choose` picks the IFMA family only when the CPU has it
        return unsafe { ifma::conv(a, b, n, nout, set, t) };
    }
    let _ = simd;
    let logn = n.trailing_zeros();
    let mut out = vec![0u64; nout * t];
    let mut fa = vec![0u64; n];
    let mut fb = match b {
        Operand::Plain(_) => vec![0u64; n],
        _ => vec![],
    };
    for (pi, pr) in set.primes[..t].iter().enumerate() {
        let tw = tw(pr, logn);
        let rk = pr.limb_weights::<N>();
        load(&mut fa, a, pr, &rk);
        dif(&mut fa, &tw, pr.p);
        // pointwise: mul(mul(A, B), c) with c = n^{-1} (M/p_i)^{-1} R^2 is A B n^{-1} (M/p_i)^{-1}
        // (A, B < 2p: A B < p 2^64), so the unnormalised inverse transform returns u_i itself
        let scale = pr.mul(set.scale(pi, t, n), pr.r2);
        let fbv: &[u64] = match &b {
            Operand::Square => {
                for x in fa.iter_mut() {
                    *x = pr.mul(pr.mul(*x, *x), scale);
                }
                &[]
            }
            Operand::Plain(b) => {
                load(&mut fb, b, pr, &rk);
                dif(&mut fb, &tw, pr.p);
                &fb
            }
            Operand::Prepared(tr) => &tr[pi],
        };
        if !fbv.is_empty() {
            for (x, y) in fa.iter_mut().zip(fbv) {
                *x = pr.mul(pr.mul(*x, *y), scale);
            }
        }
        dit(&mut fa, &tw, pr.p);
        let p2 = 2 * pr.p;
        for (k, &v) in fa[..nout].iter().enumerate() {
            let mut v = v;
            if v >= p2 {
                v -= p2;
            }
            if v >= pr.p {
                v -= pr.p;
            }
            out[k * t + pi] = v;
        }
    }
    out
}

/// The transforms on AVX-512 IFMA: 50-bit primes, eight 64-bit lanes, products from the 52-bit
/// multiply-adds (vpmadd52luq / vpmadd52huq). Same lazy butterflies as the portable path:
/// Shoup products with companions floor(w 2^52 / p), values below 4p < 2^52.
#[cfg(target_arch = "x86_64")]
mod ifma {
    use super::*;
    use std::arch::x86_64::*;

    const MASK52: u64 = (1 << 52) - 1;

    pub fn available() -> bool {
        static A: OnceLock<bool> = OnceLock::new();
        *A.get_or_init(|| {
            std::env::var_os("ISOGENY_NO_SIMD").is_none()
                && std::arch::is_x86_feature_detected!("avx512f")
                && std::arch::is_x86_feature_detected!("avx512ifma")
        })
    }

    /// 50-bit primes (44: 1024-bit fields at 2^26 points need 2 * 1024 + 28 bits).
    pub fn set() -> &'static PrimeSet {
        static S: OnceLock<PrimeSet> = OnceLock::new();
        S.get_or_init(|| PrimeSet::new(49, 44))
    }

    /// Per-layer twiddles of one prime for transforms up to 2^logn points: the layer pairing
    /// entries h apart (h = 1, 2, 4, ..) uses w_(2h)^j, j < h, stored at [h + j] with its
    /// 52-bit Shoup companion; forward and inverse. `small[l]` holds the twiddles of layer
    /// h = 2^l (l < 3) in the lane order of the in-register butterflies (`LANES`).
    struct TwV {
        logn: u32,
        f: Vec<u64>,
        fp: Vec<u64>,
        i: Vec<u64>,
        ip: Vec<u64>,
        small: [[[u64; 8]; 4]; 3],
    }

    thread_local! {
        static TWV: RefCell<HashMap<u64, Rc<TwV>>> = RefCell::new(HashMap::new());
    }

    /// For the layers h = 1, 2, 4 on sixteen entries (two vectors): the positions of the upper
    /// (`X`) halves of the butterflies, in lane order; the lower halves are at X + h.
    const LANES: [[u64; 8]; 3] = [
        [0, 2, 4, 6, 8, 10, 12, 14],
        [0, 1, 4, 5, 8, 9, 12, 13],
        [0, 1, 2, 3, 8, 9, 10, 11],
    ];

    fn twv(pr: &Prime, logn: u32) -> Rc<TwV> {
        TWV.with(|c| {
            let mut m = c.borrow_mut();
            if let Some(t) = m.get(&pr.p) {
                if t.logn >= logn {
                    return t.clone();
                }
            }
            let logn = logn.max(10);
            let n = 1usize << logn;
            let p = pr.p;
            let plain = |x: u64| pr.mul(x, 1);
            let wm = pr.pow(pr.root, 1u64 << (MAX_LOG - logn));
            let roots = [plain(wm), plain(pr.pow(wm, (1u64 << logn) - 1))];
            let mut tabs: Vec<(Vec<u64>, Vec<u64>)> = vec![];
            for w in roots {
                let (mut v, mut vp) = (vec![0u64; n], vec![0u64; n]);
                let mut h = n / 2;
                let mut wh = w; // primitive 2h-th root
                while h >= 1 {
                    let mut x = 1u64;
                    for j in 0..h {
                        v[h + j] = x;
                        vp[h + j] = (((x as u128) << 52) / p as u128) as u64;
                        x = ((x as u128 * wh as u128) % p as u128) as u64;
                    }
                    wh = ((wh as u128 * wh as u128) % p as u128) as u64;
                    h /= 2;
                }
                tabs.push((v, vp));
            }
            let mut small = [[[0u64; 8]; 4]; 3];
            for (l, s) in small.iter_mut().enumerate() {
                let h = 1usize << l;
                for (k, &x) in LANES[l].iter().enumerate() {
                    let j = h + (x as usize % (2 * h));
                    s[0][k] = tabs[0].0[j];
                    s[1][k] = tabs[0].1[j];
                    s[2][k] = tabs[1].0[j];
                    s[3][k] = tabs[1].1[j];
                }
            }
            let (f, fp) = tabs.remove(0);
            let (i, ip) = tabs.remove(0);
            let t = Rc::new(TwV {
                logn,
                f,
                fp,
                i,
                ip,
                small,
            });
            m.insert(pr.p, t.clone());
            t
        })
    }

    /// Index vectors for `_mm512_permutex2var_epi64` (entries 0..7 from the first operand,
    /// 8..15 from the second): gather the X and Y = X + h halves of layer 2^l from two
    /// consecutive vectors, and scatter them back.
    struct Perm {
        gx: __m512i,
        gy: __m512i,
        s0: __m512i,
        s1: __m512i,
    }

    #[target_feature(enable = "avx512f")]
    unsafe fn perm(l: usize) -> Perm {
        let h = 1u64 << l;
        let x = LANES[l];
        let y = x.map(|v| v + h);
        let mut back = [0u64; 16];
        for k in 0..8 {
            back[x[k] as usize] = k as u64;
            back[y[k] as usize] = 8 + k as u64;
        }
        let ld = |v: &[u64]| _mm512_loadu_si512(v.as_ptr() as *const _);
        Perm {
            gx: ld(&x),
            gy: ld(&y),
            s0: ld(&back[..8]),
            s1: ld(&back[8..]),
        }
    }

    /// w x mod p in [0, 2p) for x < 2^52 (lanewise), wp = floor(w 2^52 / p).
    #[inline(always)]
    unsafe fn shoup(x: __m512i, w: __m512i, wp: __m512i, p: __m512i) -> __m512i {
        let z = _mm512_setzero_si512();
        let q = _mm512_madd52hi_epu64(z, x, wp);
        let r = _mm512_sub_epi64(
            _mm512_madd52lo_epu64(z, x, w),
            _mm512_madd52lo_epu64(z, q, p),
        );
        _mm512_and_si512(r, _mm512_set1_epi64(MASK52 as i64))
    }

    /// x reduced from [0, 2m) to [0, m) (lanewise; m < 2^63).
    #[inline(always)]
    unsafe fn red(x: __m512i, m: __m512i) -> __m512i {
        _mm512_min_epu64(x, _mm512_sub_epi64(x, m))
    }

    /// DIF butterfly: (u + v, (u - v) w), inputs and outputs in [0, 2p).
    #[inline(always)]
    unsafe fn bf_dif(
        u: __m512i,
        v: __m512i,
        w: __m512i,
        wp: __m512i,
        p: __m512i,
        p2: __m512i,
    ) -> (__m512i, __m512i) {
        let x = red(_mm512_add_epi64(u, v), p2);
        let y = shoup(_mm512_sub_epi64(_mm512_add_epi64(u, p2), v), w, wp, p);
        (x, y)
    }

    /// DIT butterfly: (u + v w, u - v w), inputs in [0, 4p) (u) and [0, 4p) (v), outputs in
    /// [0, 4p).
    #[inline(always)]
    unsafe fn bf_dit(
        u: __m512i,
        v: __m512i,
        w: __m512i,
        wp: __m512i,
        p: __m512i,
        p2: __m512i,
    ) -> (__m512i, __m512i) {
        let u = red(u, p2);
        let t = shoup(v, w, wp, p);
        (
            _mm512_add_epi64(u, t),
            _mm512_sub_epi64(_mm512_add_epi64(u, p2), t),
        )
    }

    #[inline(always)]
    unsafe fn ld(p: *const u64) -> __m512i {
        _mm512_loadu_si512(p as *const _)
    }
    #[inline(always)]
    unsafe fn st(p: *mut u64, v: __m512i) {
        _mm512_storeu_si512(p as *mut _, v)
    }

    /// Forward transform (n >= 16 points), as the portable `dif`.
    #[target_feature(enable = "avx512f,avx512ifma")]
    unsafe fn dif(a: &mut [u64], t: &TwV, pm: &[Perm; 3], p: u64) {
        let n = a.len();
        let (vp, vp2) = (_mm512_set1_epi64(p as i64), _mm512_set1_epi64(2 * p as i64));
        let ptr = a.as_mut_ptr();
        let mut h = n / 2;
        while h >= 8 {
            for blk in (0..n).step_by(2 * h) {
                for j in (0..h).step_by(8) {
                    let (pu, pv) = (ptr.add(blk + j), ptr.add(blk + j + h));
                    let w = ld(t.f.as_ptr().add(h + j));
                    let wp = ld(t.fp.as_ptr().add(h + j));
                    let (x, y) = bf_dif(ld(pu), ld(pv), w, wp, vp, vp2);
                    st(pu, x);
                    st(pv, y);
                }
            }
            h /= 2;
        }
        // h = 4, 2, 1 on sixteen entries held in two registers
        let tw: [(__m512i, __m512i); 3] =
            std::array::from_fn(|l| (ld(t.small[l][0].as_ptr()), ld(t.small[l][1].as_ptr())));
        for c in (0..n).step_by(16) {
            let (mut v0, mut v1) = (ld(ptr.add(c)), ld(ptr.add(c + 8)));
            for l in (0..3).rev() {
                let m = &pm[l];
                let u = _mm512_permutex2var_epi64(v0, m.gx, v1);
                let v = _mm512_permutex2var_epi64(v0, m.gy, v1);
                let (x, y) = bf_dif(u, v, tw[l].0, tw[l].1, vp, vp2);
                v0 = _mm512_permutex2var_epi64(x, m.s0, y);
                v1 = _mm512_permutex2var_epi64(x, m.s1, y);
            }
            st(ptr.add(c), v0);
            st(ptr.add(c + 8), v1);
        }
    }

    /// Inverse transform (n >= 16 points), as the portable `dit`: output in [0, 4p).
    #[target_feature(enable = "avx512f,avx512ifma")]
    unsafe fn dit(a: &mut [u64], t: &TwV, pm: &[Perm; 3], p: u64) {
        let n = a.len();
        let (vp, vp2) = (_mm512_set1_epi64(p as i64), _mm512_set1_epi64(2 * p as i64));
        let ptr = a.as_mut_ptr();
        let tw: [(__m512i, __m512i); 3] =
            std::array::from_fn(|l| (ld(t.small[l][2].as_ptr()), ld(t.small[l][3].as_ptr())));
        for c in (0..n).step_by(16) {
            let (mut v0, mut v1) = (ld(ptr.add(c)), ld(ptr.add(c + 8)));
            for (l, m) in pm.iter().enumerate() {
                let u = _mm512_permutex2var_epi64(v0, m.gx, v1);
                let v = _mm512_permutex2var_epi64(v0, m.gy, v1);
                let (x, y) = bf_dit(u, v, tw[l].0, tw[l].1, vp, vp2);
                v0 = _mm512_permutex2var_epi64(x, m.s0, y);
                v1 = _mm512_permutex2var_epi64(x, m.s1, y);
            }
            st(ptr.add(c), v0);
            st(ptr.add(c + 8), v1);
        }
        let mut h = 8;
        while h < n {
            for blk in (0..n).step_by(2 * h) {
                for j in (0..h).step_by(8) {
                    let (pu, pv) = (ptr.add(blk + j), ptr.add(blk + j + h));
                    let w = ld(t.i.as_ptr().add(h + j));
                    let wp = ld(t.ip.as_ptr().add(h + j));
                    let (x, y) = bf_dit(ld(pu), ld(pv), w, wp, vp, vp2);
                    st(pu, x);
                    st(pv, y);
                }
            }
            h *= 2;
        }
    }

    /// a_k <- a_k b_k c mod p in [0, 2p) for a_k, b_k < 2p: a Montgomery product in radix
    /// 2^52 (a b < 4p^2 < p 2^52, result below 2p), then a Shoup product by the constant
    /// c = 2^52 (the Montgomery factor) times the scale.
    #[target_feature(enable = "avx512f,avx512ifma")]
    unsafe fn pointwise(a: &mut [u64], b: Option<&[u64]>, p: u64, pneg: u64, c: u64) {
        let z = _mm512_setzero_si512();
        let (vp, vneg) = (_mm512_set1_epi64(p as i64), _mm512_set1_epi64(pneg as i64));
        let vc = _mm512_set1_epi64(c as i64);
        let vcp = _mm512_set1_epi64(((c as u128) << 52).div_euclid(p as u128) as i64);
        let one = _mm512_set1_epi64(1);
        let pa = a.as_mut_ptr();
        for k in (0..a.len()).step_by(8) {
            let x = ld(pa.add(k));
            let y = match b {
                Some(b) => ld(b.as_ptr().add(k)),
                None => x,
            };
            let lo = _mm512_madd52lo_epu64(z, x, y);
            let hi = _mm512_madd52hi_epu64(z, x, y);
            let m = _mm512_madd52lo_epu64(z, lo, vneg);
            // (lo + m p) / 2^52 = hi + high(m p) + carry, carry = [lo != 0]
            let t = _mm512_madd52hi_epu64(hi, m, vp);
            let t = _mm512_mask_add_epi64(t, _mm512_cmpneq_epu64_mask(lo, z), t, one);
            st(pa.add(k), shoup(t, vc, vcp, vp));
        }
    }

    /// [0, 4p) to [0, p), in place.
    #[target_feature(enable = "avx512f")]
    unsafe fn normalize(a: &mut [u64], p: u64) {
        let (vp, vp2) = (_mm512_set1_epi64(p as i64), _mm512_set1_epi64(2 * p as i64));
        for k in (0..a.len()).step_by(8) {
            let x = ld(a.as_ptr().add(k));
            st(a.as_mut_ptr().add(k), red(red(x, vp2), vp));
        }
    }

    /// Radix-2^52 digits of the entries, digit-major (digit j of entry k at [j * len8 + k],
    /// len8 = the length rounded up to 8, padding entries zero): shared by every prime.
    pub fn digits52<const N: usize>(a: &[[u64; N]]) -> (Vec<u64>, usize, usize) {
        let d = (64 * N).div_ceil(52);
        let len8 = a.len().div_ceil(8) * 8;
        let mut out = vec![0u64; d * len8];
        for (k, x) in a.iter().enumerate() {
            for j in 0..d {
                let (bit, w) = (52 * j, 52 * j / 64);
                let sh = bit % 64;
                let mut v = x[w] >> sh;
                if sh > 12 && w + 1 < N {
                    v |= x[w + 1] << (64 - sh);
                }
                out[j * len8 + k] = v & MASK52;
            }
        }
        (out, d, len8)
    }

    /// The transform input from the digits of `digits52`: entry x (D digits) becomes
    /// x 2^(-52 D) mod p, below 2p, by D digit-serial Montgomery steps t <- (t + d_j + m p) /
    /// 2^52 (m = -(t + d_j) p^{-1} mod 2^52; t stays below p + 2), folded mod x^len - 1.
    #[target_feature(enable = "avx512f,avx512ifma")]
    unsafe fn load(f: &mut [u64], dg: &(Vec<u64>, usize, usize), p: u64, pneg: u64) {
        let (dg, d, len8) = (&dg.0, dg.1, dg.2);
        let n = f.len();
        let z = _mm512_setzero_si512();
        let (vp, vp2) = (_mm512_set1_epi64(p as i64), _mm512_set1_epi64(2 * p as i64));
        let vneg = _mm512_set1_epi64(pneg as i64);
        if len8 < n {
            f[len8..].iter_mut().for_each(|x| *x = 0);
        }
        let fp = f.as_mut_ptr();
        for c in (0..len8).step_by(8) {
            let mut t = z;
            for j in 0..d {
                let x = _mm512_add_epi64(t, ld(dg.as_ptr().add(j * len8 + c)));
                let m = _mm512_madd52lo_epu64(z, x, vneg);
                let s = _mm512_add_epi64(x, _mm512_madd52lo_epu64(z, m, vp));
                t = _mm512_add_epi64(_mm512_madd52hi_epu64(z, m, vp), _mm512_srli_epi64(s, 52));
            }
            if c < n {
                st(fp.add(c), t);
            } else {
                let q = fp.add(c & (n - 1));
                st(q, red(_mm512_add_epi64(ld(q), t), vp2));
            }
        }
    }

    /// The forward transform of a loaded operand (from its `digits52`), as stored by
    /// `ntt::prepare`.
    ///
    /// # Safety
    /// The CPU must support AVX-512F and AVX-512 IFMA.
    pub unsafe fn forward(f: &mut [u64], dg: &(Vec<u64>, usize, usize), pr: &Prime) {
        let tv = twv(pr, f.len().trailing_zeros());
        let pm: [Perm; 3] = std::array::from_fn(|l| perm(l));
        load(f, dg, pr.p, pr.pinv & MASK52);
        dif(f, &tv, &pm, pr.p);
    }

    /// `ntt::run` over the 50-bit primes; n >= 16 a power of two.
    ///
    /// # Safety
    /// The CPU must support AVX-512F and AVX-512 IFMA.
    pub(super) unsafe fn conv<const N: usize>(
        a: &[[u64; N]],
        b: Operand<N>,
        n: usize,
        nout: usize,
        set: &PrimeSet,
        t: usize,
    ) -> Vec<u64> {
        let logn = n.trailing_zeros();
        let pm: [Perm; 3] = std::array::from_fn(|l| perm(l));
        let da = digits52(a);
        let db = match b {
            Operand::Plain(b) => Some(digits52(b)),
            _ => None,
        };
        let mut out = vec![0u64; nout * t];
        let mut fa = vec![0u64; n];
        let mut fb = if db.is_some() { vec![0u64; n] } else { vec![] };
        // each loaded operand carries 2^(-52 D); a prepared one was loaded the same way
        let d = (64 * N).div_ceil(52) as u64;
        for (pi, pr) in set.primes[..t].iter().enumerate() {
            let p = pr.p;
            let tv = twv(pr, logn);
            let pneg = pr.pinv & MASK52;
            load(&mut fa, &da, p, pneg);
            dif(&mut fa, &tv, &pm, p);
            // the plain scale times 2^52 (undoes the radix-2^52 Montgomery division of the
            // pointwise product) times 2^(104 D) (undoes the 2^(-52 D) of both operands)
            let s = pr.mul(set.scale(pi, t, n), 1);
            let f = pr.mul(pr.pow(pr.mont(2), 52 + 104 * d), 1);
            let c = ((s as u128 * f as u128) % p as u128) as u64;
            match &b {
                Operand::Square => pointwise(&mut fa, None, p, pneg, c),
                Operand::Plain(_) => {
                    load(&mut fb, db.as_ref().unwrap(), p, pneg);
                    dif(&mut fb, &tv, &pm, p);
                    pointwise(&mut fa, Some(&fb), p, pneg, c);
                }
                Operand::Prepared(tr) => pointwise(&mut fa, Some(&tr[pi]), p, pneg, c),
            }
            dit(&mut fa, &tv, &pm, p);
            normalize(&mut fa, p);
            for (k, &v) in fa[..nout].iter().enumerate() {
                out[k * t + pi] = v;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bigint::Big;

    #[test]
    fn primes_have_roots() {
        #[cfg_attr(not(target_arch = "x86_64"), allow(unused_mut))]
        let mut sets = vec![scalar_set()];
        #[cfg(target_arch = "x86_64")]
        sets.push(ifma::set());
        for set in sets {
            for pr in &set.primes {
                assert_eq!((pr.p - 1) % (1u64 << MAX_LOG), 0);
                let one = pr.mont(1);
                assert_eq!(pr.pow(pr.root, 1u64 << MAX_LOG), one);
                assert_ne!(pr.pow(pr.root, 1u64 << (MAX_LOG - 1)), one);
            }
        }
    }

    /// V_k = sum_i u_i (M / p_i) - q M, as a Big.
    fn recombine(u: &[u64], set: &PrimeSet) -> Big {
        let ps = &set.primes[..u.len()];
        let m = ps.iter().fold(Big::from_u64(1), |m, pr| m.mul_small(pr.p));
        let mut s = Big::zero();
        for (&x, pr) in u.iter().zip(ps) {
            s = s.add(&m.divrem_small(pr.p).0.mul_small(x));
        }
        s.sub(&m.mul_small(set.quotient(u)))
    }

    /// The exact (cyclic when `cyc` is short) convolution of small vectors, as Bigs.
    fn naive<const N: usize>(a: &[[u64; N]], b: &[[u64; N]], cyc: usize) -> Vec<Big> {
        let full = a.len() + b.len() - 1;
        let mut r = vec![Big::zero(); full.min(cyc)];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                let k = (i + j) % cyc.min(full);
                r[k] = r[k].add(&Big::from_limbs(x).mul(&Big::from_limbs(y)));
            }
        }
        r
    }

    fn check<const N: usize>(a: &[[u64; N]], b: &[[u64; N]], bits: u32, cyc: usize) {
        let want = naive(a, b, cyc);
        #[cfg_attr(not(target_arch = "x86_64"), allow(unused_mut))]
        let mut modes = vec![false];
        #[cfg(target_arch = "x86_64")]
        if ifma::available() {
            modes.push(true);
        }
        for simd in modes {
            SIMD.with(|s| s.set(simd));
            let (u, t, set) = conv_crt(a, b, bits, cyc);
            SIMD.with(|s| s.set(true));
            for (k, w) in want.iter().enumerate() {
                let got = recombine(&u[k * t..(k + 1) * t], set);
                assert_eq!(
                    &got,
                    w,
                    "{}x{} cyc {cyc} simd {simd} coefficient {k}",
                    a.len(),
                    b.len()
                );
            }
        }
    }

    #[test]
    fn single_limb_against_naive() {
        let a: Vec<[u64; 1]> = (0..37u64).map(|i| [i * 1_000_003 + 7]).collect();
        let b: Vec<[u64; 1]> = (0..50u64).map(|i| [i * i + 1]).collect();
        check(&a, &b, 40, usize::MAX);
        // pseudo-random full-size entries, products and squares, plain and cyclic
        let mut x = 12345u64;
        let mut r = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            x >> 3
        };
        for (la, lb) in [(16usize, 16usize), (100, 29), (319, 319), (64, 70)] {
            let a: Vec<[u64; 1]> = (0..la).map(|_| [r()]).collect();
            let b: Vec<[u64; 1]> = (0..lb).map(|_| [r()]).collect();
            check(&a, &b, 61, usize::MAX);
            check(&a, &a, 61, usize::MAX);
            let cyc = (la + lb - 1).next_power_of_two() / 2;
            if la.max(lb) <= cyc + 1 {
                check(&a, &b, 61, cyc);
            }
        }
    }

    #[test]
    fn multi_limb_extremes() {
        // all-ones entries (the largest products) and zeros: quotient near both ends
        let ones = [u64::MAX; 4];
        for (la, lb) in [(1usize, 1usize), (3, 200), (64, 64)] {
            let a = vec![ones; la];
            let mut b = vec![ones; lb];
            b[lb / 2] = [0; 4];
            check(&a, &b, 256, usize::MAX);
        }
        let ones8 = [u64::MAX; 8];
        let a = vec![ones8; 40];
        check(&a, &a, 512, usize::MAX);
        check(&a, &a, 512, 64);
    }
}
