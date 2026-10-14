//! IC improvement round, item 04: genus-3 double-large-prime index calculus after
//! the cubic GHS transfer, ported from
//! `research/ic_next10_20261014/07-cubic-ghs-genus3/ghs3.py` to Rust and run beyond
//! q = 2^15.
//!
//! The port is operation-for-operation: the same `F_q` operation counters
//! (`CNT[0]` multiplications and squarings, `CNT[1]` inversions, `CNT[2]`
//! multiplications mod ell, `CNT[3]` inversions mod ell, `CNT[4]` group
//! operations) and the same random stream (CPython's MT19937 `random.Random`,
//! `randrange` and `sample` are reproduced bit for bit), so `--mode python`
//! regenerates item 7's instances, walks and counted costs exactly.
//!
//! Protocol: `research/ic_improve_20261014/04-genus3-2lp/PROTOCOL.md`.
//!
//! ```sh
//! cargo run --release --example ic_improve_04_genus3_2lp -- \
//!     --bits 12 --curve 0 --mode paired --arms basic,ic2lp,ic1lp,aa --out out.jsonl
//! ```
//!
//! Arms: `basic` (full factor base, incumbent), `ic2lp` (kappa = 2, q^{2/3}
//! factor base, item 7's candidate), `ic1lp` (single large prime, kappa1 q^{5/7}),
//! `ic2lpt` (kappa2 q^{2/3}), `aa` (ic2lp on an independent walk stream), `rho`
//! (negation-map rho on E0(F_{q^3}), reference only).

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::process::{Command, Stdio};
use std::time::Instant;

// ------------------------------------------------------------------ counters
static mut CNT: [u64; 5] = [0; 5];

#[inline(always)]
fn cnt(i: usize, n: u64) {
    unsafe {
        let p = std::ptr::addr_of_mut!(CNT);
        (*p)[i] += n;
    }
}
fn snap() -> [u64; 5] {
    unsafe { *std::ptr::addr_of!(CNT) }
}
fn restore(s: [u64; 5]) {
    unsafe {
        *std::ptr::addr_of_mut!(CNT) = s;
    }
}
fn reset() {
    restore([0; 5]);
}

// ----------------------------------------------------- CPython random.Random
#[derive(Clone)]
struct PyRng {
    mt: [u32; 624],
    idx: usize,
}

impl PyRng {
    fn init_genrand(s: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = s;
        for i in 1..624 {
            mt[i] = 1812433253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        PyRng { mt, idx: 624 }
    }
    fn from_seed(seed: u64) -> Self {
        let mut key: Vec<u32> = vec![seed as u32];
        if seed >> 32 != 0 {
            key.push((seed >> 32) as u32);
        }
        let mut r = Self::init_genrand(19650218);
        let n = 624usize;
        let (mut i, mut j) = (1usize, 0usize);
        let mut k = n.max(key.len());
        while k > 0 {
            let prev = r.mt[i - 1];
            r.mt[i] = (r.mt[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= n {
                r.mt[0] = r.mt[n - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }
        k = n - 1;
        while k > 0 {
            let prev = r.mt[i - 1];
            r.mt[i] = (r.mt[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1566083941)).wrapping_sub(i as u32);
            i += 1;
            if i >= n {
                r.mt[0] = r.mt[n - 1];
                i = 1;
            }
            k -= 1;
        }
        r.mt[0] = 0x80000000;
        r
    }
    fn genrand(&mut self) -> u32 {
        const N: usize = 624;
        const M: usize = 397;
        if self.idx >= N {
            let mag = [0u32, 0x9908b0df];
            let mt = &mut self.mt;
            for kk in 0..N {
                let y = (mt[kk] & 0x80000000) | (mt[(kk + 1) % N] & 0x7fffffff);
                mt[kk] = mt[(kk + M) % N] ^ (y >> 1) ^ mag[(y & 1) as usize];
            }
            self.idx = 0;
        }
        let mut y = self.mt[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }
    fn getrandbits(&mut self, k: u32) -> u64 {
        if k == 0 {
            return 0;
        }
        if k <= 32 {
            return (self.genrand() >> (32 - k)) as u64;
        }
        let mut r: u64 = 0;
        let mut kk = k as i64;
        let mut sh = 0;
        while kk > 0 {
            let mut w = self.genrand();
            if kk < 32 {
                w >>= 32 - kk;
            }
            r |= (w as u64) << sh;
            sh += 32;
            kk -= 32;
        }
        r
    }
    fn randbelow(&mut self, n: u64) -> u64 {
        let k = 64 - n.leading_zeros();
        let mut r = self.getrandbits(k);
        while r >= n {
            r = self.getrandbits(k);
        }
        r
    }
    fn randrange(&mut self, a: u64, b: u64) -> u64 {
        a + self.randbelow(b - a)
    }
    /// `rng.sample(range(n), 2)` for n > 21.
    fn sample2(&mut self, n: u64) -> (u64, u64) {
        let a = self.randbelow(n);
        let mut b = self.randbelow(n);
        while b == a {
            b = self.randbelow(n);
        }
        (a, b)
    }
}

// ------------------------------------------------------------------ integers
fn powmod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r: u64 = 1 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = ((r as u128 * b as u128) % m as u128) as u64;
        }
        b = ((b as u128 * b as u128) % m as u128) as u64;
        e >>= 1;
    }
    r
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &sp in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % sp == 0 {
            return n == sp;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'outer: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = ((x as u128 * x as u128) % n as u128) as u64;
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn next_prime(mut n: u64) -> u64 {
    while !is_prime(n) {
        n += 1;
    }
    n
}

fn sqrt_fp(a: u64, p: u64) -> u64 {
    let a = a % p;
    if a == 0 {
        return 0;
    }
    if p % 4 == 3 {
        return powmod(a, (p + 1) / 4, p);
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while powmod(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, powmod(z, q, p), powmod(a, q, p), powmod(a, (q + 1) / 2, p));
    while t != 1 {
        let (mut i, mut tt) = (0u32, t);
        while tt != 1 {
            tt = tt * tt % p;
            i += 1;
        }
        let b = powmod(c, 1u64 << (m - i - 1), p);
        m = i;
        c = b * b % p;
        t = t * b % p * b % p;
        r = r * b % p;
    }
    r
}

fn bitlen(x: u64) -> u64 {
    64 - x.leading_zeros() as u64
}

// ---------------------------------------------------------------- F_p polys
type P = Vec<u64>;

fn tr(mut a: P) -> P {
    while let Some(&0) = a.last() {
        a.pop();
    }
    a
}
fn padd(a: &P, b: &P, p: u64) -> P {
    let (a, b) = if a.len() < b.len() { (b, a) } else { (a, b) };
    let mut r = a.clone();
    for (i, &c) in b.iter().enumerate() {
        r[i] = (r[i] + c) % p;
    }
    tr(r)
}
fn psub(a: &P, b: &P, p: u64) -> P {
    let n = a.len().max(b.len());
    let mut r = vec![0u64; n];
    r[..a.len()].copy_from_slice(a);
    for (i, &c) in b.iter().enumerate() {
        r[i] = (r[i] + p - c) % p;
    }
    tr(r)
}
fn pneg(a: &P, p: u64) -> P {
    a.iter().map(|&c| (p - c) % p).collect()
}
fn pmul(a: &P, b: &P, p: u64) -> P {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    cnt(0, (a.len() * b.len()) as u64);
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x != 0 {
            for (j, &y) in b.iter().enumerate() {
                r[i + j] = (r[i + j] + x * y) % p;
            }
        }
    }
    tr(r)
}
fn pscal(a: &P, c: u64, p: u64) -> P {
    if c == 1 {
        return a.clone();
    }
    cnt(0, a.len() as u64);
    tr(a.iter().map(|&x| x * c % p).collect())
}
fn pdivmod(a: &P, b: &P, p: u64) -> (P, P) {
    let lb = b.len();
    if a.len() < lb {
        return (vec![], a.clone());
    }
    let mut r = a.clone();
    let lead = b[lb - 1];
    let mut inv = 1;
    if lead != 1 {
        cnt(1, 1);
        inv = powmod(lead, p - 2, p);
    }
    let mut q = vec![0u64; a.len() - lb + 1];
    for i in (0..=a.len() - lb).rev() {
        let mut c = r[i + lb - 1];
        if c != 0 {
            if inv != 1 {
                c = c * inv % p;
                cnt(0, 1);
            }
            q[i] = c;
            cnt(0, (lb - 1) as u64);
            for j in 0..lb - 1 {
                r[i + j] = (r[i + j] + p - c * b[j] % p) % p;
            }
            r[i + lb - 1] = 0;
        }
    }
    r.truncate(lb - 1);
    (tr(q), tr(r))
}
fn pmod(a: &P, b: &P, p: u64) -> P {
    pdivmod(a, b, p).1
}
fn pmonic(a: &P, p: u64) -> P {
    if a.is_empty() || *a.last().unwrap() == 1 {
        return a.clone();
    }
    cnt(1, 1);
    let inv = powmod(*a.last().unwrap(), p - 2, p);
    cnt(0, (a.len() - 1) as u64);
    let mut r: P = a[..a.len() - 1].iter().map(|&c| c * inv % p).collect();
    r.push(1);
    r
}
fn pxgcd(a: &P, b: &P, p: u64) -> (P, P, P) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1): (P, P) = (vec![1], vec![]);
    let (mut t0, mut t1): (P, P) = (vec![], vec![1]);
    while !r1.is_empty() {
        let (q, r) = pdivmod(&r0, &r1, p);
        r0 = std::mem::replace(&mut r1, r);
        let ns = psub(&s0, &pmul(&q, &s1, p), p);
        s0 = std::mem::replace(&mut s1, ns);
        let nt = psub(&t0, &pmul(&q, &t1, p), p);
        t0 = std::mem::replace(&mut t1, nt);
    }
    if r0.is_empty() {
        return (vec![], vec![], vec![]);
    }
    let lead = *r0.last().unwrap();
    if lead != 1 {
        cnt(1, 1);
        let inv = powmod(lead, p - 2, p);
        r0 = pscal(&r0, inv, p);
        s0 = pscal(&s0, inv, p);
        t0 = pscal(&t0, inv, p);
    }
    (r0, s0, t0)
}
fn pev(a: &P, x: u64, p: u64) -> u64 {
    let mut r = 0u64;
    for &c in a.iter().rev() {
        r = (r * x + c) % p;
    }
    cnt(0, a.len().saturating_sub(1) as u64);
    r
}
fn pder(a: &P, p: u64) -> P {
    tr((1..a.len()).map(|i| (i as u64 % p) * a[i] % p).collect())
}

// ---------------------------------------------------------- Cantor over F_p
const G3: usize = 3;
type Div = (P, P);

fn creduce(mut u: P, mut v: P, h: &P, p: u64) -> Div {
    while u.len() - 1 > G3 {
        let (u2, rem) = pdivmod(&psub(h, &pmul(&v, &v, p), p), &u, p);
        assert!(rem.is_empty());
        u = pmonic(&u2, p);
        v = pmod(&pneg(&v, p), &u, p);
    }
    let u = pmonic(&u, p);
    let v = pmod(&v, &u, p);
    (u, v)
}

fn is_one(u: &P) -> bool {
    u.len() == 1 && u[0] == 1
}

fn cadd(d1: &Div, d2: &Div, h: &P, p: u64) -> Div {
    cnt(4, 1);
    let (u1, v1) = d1;
    let (u2, v2) = d2;
    if is_one(u1) {
        return d2.clone();
    }
    if is_one(u2) {
        return d1.clone();
    }
    let (dd1, e1, e2) = pxgcd(u1, u2, p);
    if is_one(&dd1) {
        let t = pmod(&pmul(&e1, &psub(v2, v1, p), p), u2, p);
        let u = pmul(u1, u2, p);
        let v = padd(v1, &pmul(u1, &t, p), p);
        return creduce(u, v, h, p);
    }
    let (d, c1, c2) = pxgcd(&dd1, &padd(v1, v2, p), p);
    let s1 = pmul(&c1, &e1, p);
    let s2 = pmul(&c1, &e2, p);
    let s3 = c2;
    let u = pdivmod(&pmul(u1, u2, p), &pmul(&d, &d, p), p).0;
    let num = padd(
        &padd(&pmul(&pmul(&s1, u1, p), v2, p), &pmul(&pmul(&s2, u2, p), v1, p), p),
        &pmul(&s3, &padd(&pmul(v1, v2, p), h, p), p),
        p,
    );
    let (v, rem) = pdivmod(&num, &d, p);
    assert!(rem.is_empty());
    let v = pmod(&v, &u, p);
    creduce(u, v, h, p)
}

fn cmul(mut n: u64, d: &Div, h: &P, p: u64) -> Div {
    let mut r: Div = (vec![1], vec![]);
    let mut b = d.clone();
    while n > 0 {
        if n & 1 == 1 {
            r = cadd(&r, &b, h, p);
        }
        n >>= 1;
        if n > 0 {
            b = cadd(&b, &b, h, p);
        }
    }
    r
}

// ------------------------------------------------------------ K = F_{p^3}
type KE = [u64; 3];
const KZ: KE = [0, 0, 0];
const KONE: KE = [1, 0, 0];

#[derive(Clone)]
struct KField {
    p: u64,
    c1: u64,
    c0: u64,
    t1: KE,
    t2: KE,
}

impl KField {
    fn new(p: u64, rng: &mut PyRng) -> Self {
        let mut c1s: Vec<u64> = (1..20).collect();
        for _ in 0..200 {
            c1s.push(rng.randrange(1, p));
        }
        for &c1 in &c1s {
            for _ in 0..200 {
                let c0 = rng.randrange(1, p);
                if Self::irreducible(p, c1, c0) {
                    let mut k = KField { p, c1, c0, t1: KZ, t2: KZ };
                    let old = snap();
                    k.t1 = k.pw([0, 1, 0], p);
                    k.t2 = k.mul(k.t1, k.t1);
                    restore(old);
                    return k;
                }
            }
        }
        panic!("no modulus");
    }
    fn irreducible(p: u64, c1: u64, c0: u64) -> bool {
        let f: P = vec![(p - c0) % p, (p - c1) % p, 0, 1];
        let old = snap();
        let mut res: P = vec![1];
        let mut base: P = vec![0, 1];
        let mut e = p;
        while e > 0 {
            if e & 1 == 1 {
                res = pmod(&pmul(&res, &base, p), &f, p);
            }
            base = pmod(&pmul(&base, &base, p), &f, p);
            e >>= 1;
        }
        let g = pxgcd(&f, &psub(&res, &vec![0, 1], p), p).0;
        restore(old);
        is_one(&g)
    }
    #[inline]
    fn red(&self, d: [u64; 5]) -> KE {
        let p = self.p;
        let d3 = d[3] % p;
        let d4 = d[4] % p;
        [
            (d[0] + self.c0 * d3) % p,
            (d[1] % p + self.c1 * d3 + self.c0 * d4) % p,
            (d[2] % p + self.c1 * d4) % p,
        ]
    }
    #[inline]
    fn mul(&self, a: KE, b: KE) -> KE {
        let d = [
            a[0] * b[0],
            a[0] * b[1] + a[1] * b[0],
            a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
            a[1] * b[2] + a[2] * b[1],
            a[2] * b[2],
        ];
        cnt(0, 9 + 2 + if self.c1 == 1 { 0 } else { 2 });
        self.red(d)
    }
    #[inline]
    fn sqr(&self, a: KE) -> KE {
        let d = [
            a[0] * a[0],
            2 * a[0] * a[1],
            2 * a[0] * a[2] + a[1] * a[1],
            2 * a[1] * a[2],
            a[2] * a[2],
        ];
        cnt(0, 6 + 2 + if self.c1 == 1 { 0 } else { 2 });
        self.red(d)
    }
    #[inline]
    fn add(&self, a: KE, b: KE) -> KE {
        let p = self.p;
        [(a[0] + b[0]) % p, (a[1] + b[1]) % p, (a[2] + b[2]) % p]
    }
    #[inline]
    fn sub(&self, a: KE, b: KE) -> KE {
        let p = self.p;
        [(a[0] + p - b[0]) % p, (a[1] + p - b[1]) % p, (a[2] + p - b[2]) % p]
    }
    #[inline]
    fn neg(&self, a: KE) -> KE {
        let p = self.p;
        [(p - a[0]) % p, (p - a[1]) % p, (p - a[2]) % p]
    }
    #[inline]
    fn smul(&self, a: KE, c: u64) -> KE {
        let p = self.p;
        cnt(0, 3);
        let c = c % p;
        [a[0] * c % p, a[1] * c % p, a[2] * c % p]
    }
    #[inline]
    fn frob(&self, a: KE) -> KE {
        let p = self.p;
        let (t1, t2) = (self.t1, self.t2);
        cnt(0, 6);
        [
            (a[0] + a[1] * t1[0] % p + a[2] * t2[0] % p) % p,
            (a[1] * t1[1] % p + a[2] * t2[1] % p) % p,
            (a[1] * t1[2] % p + a[2] * t2[2] % p) % p,
        ]
    }
    fn inv(&self, a: KE) -> KE {
        let p = self.p;
        let b1 = self.frob(a);
        let b2 = self.frob(b1);
        let b = self.mul(b1, b2);
        let n = (a[0] * b[0] % p + self.c0 * ((a[1] * b[2] + a[2] * b[1]) % p)) % p;
        cnt(0, 4);
        cnt(1, 1);
        let ni = powmod(n, p - 2, p);
        self.smul(b, ni)
    }
    fn pw(&self, a: KE, mut e: u64) -> KE {
        let mut r = KONE;
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            e >>= 1;
            if e > 0 {
                b = self.sqr(b);
            }
        }
        r
    }
    fn in_fp(a: KE) -> bool {
        a[1] == 0 && a[2] == 0
    }
    fn emb(&self, c: u64) -> KE {
        [c % self.p, 0, 0]
    }
    fn rand(&self, rng: &mut PyRng) -> KE {
        let p = self.p;
        let a = rng.randbelow(p);
        let b = rng.randbelow(p);
        let c = rng.randbelow(p);
        [a, b, c]
    }
    fn sqrt(&self, a: KE, rng: &mut PyRng) -> Option<KE> {
        let old = snap();
        let r = (|| {
            if a == KZ {
                return Some(a);
            }
            let n = self.p * self.p * self.p;
            if self.pw(a, (n - 1) / 2) != KONE {
                return None;
            }
            let (mut q, mut s) = (n - 1, 0u32);
            while q % 2 == 0 {
                q /= 2;
                s += 1;
            }
            let z = loop {
                let z = self.rand(rng);
                if z != KZ && self.pw(z, (n - 1) / 2) != KONE {
                    break z;
                }
            };
            let (mut m, mut c, mut t, mut r) = (s, self.pw(z, q), self.pw(a, q), self.pw(a, (q + 1) / 2));
            while t != KONE {
                let (mut i, mut tt) = (0u32, t);
                while tt != KONE {
                    tt = self.sqr(tt);
                    i += 1;
                }
                let mut b = c;
                for _ in 0..(m - i - 1) {
                    b = self.sqr(b);
                }
                m = i;
                c = self.sqr(b);
                t = self.mul(t, c);
                r = self.mul(r, b);
            }
            Some(r)
        })();
        restore(old);
        r
    }
}

// generic polys over K, low -> high
type KP = Vec<KE>;

fn ktr(mut a: KP) -> KP {
    while let Some(&x) = a.last() {
        if x == KZ {
            a.pop();
        } else {
            break;
        }
    }
    a
}
fn kpadd(k: &KField, a: &KP, b: &KP) -> KP {
    let n = a.len().max(b.len());
    ktr((0..n)
        .map(|i| k.add(*a.get(i).unwrap_or(&KZ), *b.get(i).unwrap_or(&KZ)))
        .collect())
}
fn kpsub(k: &KField, a: &KP, b: &KP) -> KP {
    let n = a.len().max(b.len());
    ktr((0..n)
        .map(|i| k.sub(*a.get(i).unwrap_or(&KZ), *b.get(i).unwrap_or(&KZ)))
        .collect())
}
fn kpneg(k: &KField, a: &KP) -> KP {
    a.iter().map(|&c| k.neg(c)).collect()
}
fn kpmul(k: &KField, a: &KP, b: &KP) -> KP {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![KZ; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x != KZ {
            for (j, &y) in b.iter().enumerate() {
                if y != KZ {
                    r[i + j] = k.add(r[i + j], k.mul(x, y));
                }
            }
        }
    }
    ktr(r)
}
fn kpscal(k: &KField, a: &KP, c: KE) -> KP {
    ktr(a.iter().map(|&x| k.mul(x, c)).collect())
}
fn kpdivmod(k: &KField, a: &KP, b: &KP) -> (KP, KP) {
    let lb = b.len();
    if a.len() < lb {
        return (vec![], a.clone());
    }
    let mut r = a.clone();
    let inv = k.inv(b[lb - 1]);
    let mut q = vec![KZ; a.len() - lb + 1];
    for i in (0..=a.len() - lb).rev() {
        let mut c = r[i + lb - 1];
        if c != KZ {
            c = k.mul(c, inv);
            q[i] = c;
            for j in 0..lb - 1 {
                r[i + j] = k.sub(r[i + j], k.mul(c, b[j]));
            }
            r[i + lb - 1] = KZ;
        }
    }
    r.truncate(lb - 1);
    (ktr(q), ktr(r))
}
fn kpmod(k: &KField, a: &KP, b: &KP) -> KP {
    kpdivmod(k, a, b).1
}
fn kpmonic(k: &KField, a: &KP) -> KP {
    if a.is_empty() {
        return a.clone();
    }
    kpscal(k, a, k.inv(*a.last().unwrap()))
}
fn kpxgcd(k: &KField, a: &KP, b: &KP) -> (KP, KP, KP) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1): (KP, KP) = (vec![KONE], vec![]);
    let (mut t0, mut t1): (KP, KP) = (vec![], vec![KONE]);
    while !r1.is_empty() {
        let (q, r) = kpdivmod(k, &r0, &r1);
        r0 = std::mem::replace(&mut r1, r);
        let ns = kpsub(k, &s0, &kpmul(k, &q, &s1));
        s0 = std::mem::replace(&mut s1, ns);
        let nt = kpsub(k, &t0, &kpmul(k, &q, &t1));
        t0 = std::mem::replace(&mut t1, nt);
    }
    let inv = k.inv(*r0.last().unwrap());
    (kpscal(k, &r0, inv), kpscal(k, &s0, inv), kpscal(k, &t0, inv))
}
fn kcreduce(k: &KField, mut u: KP, mut v: KP, h: &KP) -> (KP, KP) {
    while u.len() - 1 > G3 {
        let (u2, rem) = kpdivmod(k, &kpsub(k, h, &kpmul(k, &v, &v)), &u);
        assert!(rem.is_empty());
        u = kpmonic(k, &u2);
        v = kpmod(k, &kpneg(k, &v), &u);
    }
    let u = kpmonic(k, &u);
    let v = kpmod(k, &v, &u);
    (u, v)
}
fn kcadd(k: &KField, d1: &(KP, KP), d2: &(KP, KP), h: &KP) -> (KP, KP) {
    cnt(4, 1);
    let (u1, v1) = d1;
    let (u2, v2) = d2;
    if u1.len() == 1 && u1[0] == KONE {
        return d2.clone();
    }
    if u2.len() == 1 && u2[0] == KONE {
        return d1.clone();
    }
    let (dd1, e1, e2) = kpxgcd(k, u1, u2);
    if dd1.len() == 1 && dd1[0] == KONE {
        let t = kpmod(k, &kpmul(k, &e1, &kpsub(k, v2, v1)), u2);
        return kcreduce(k, kpmul(k, u1, u2), kpadd(k, v1, &kpmul(k, u1, &t)), h);
    }
    let (d, c1, c2) = kpxgcd(k, &dd1, &kpadd(k, v1, v2));
    let s1 = kpmul(k, &c1, &e1);
    let s2 = kpmul(k, &c1, &e2);
    let s3 = c2;
    let u = kpdivmod(k, &kpmul(k, u1, u2), &kpmul(k, &d, &d)).0;
    let num = kpadd(
        k,
        &kpadd(k, &kpmul(k, &kpmul(k, &s1, u1), v2), &kpmul(k, &kpmul(k, &s2, u2), v1)),
        &kpmul(k, &s3, &kpadd(k, &kpmul(k, v1, v2), h)),
    );
    let (v, rem) = kpdivmod(k, &num, &d);
    assert!(rem.is_empty());
    let vm = kpmod(k, &v, &u);
    kcreduce(k, u, vm, h)
}

// -------------------------------------------- elliptic curve over K (short W.)
type EP = Option<(KE, KE)>;

#[derive(Clone)]
struct Ecw {
    k: KField,
    a: KE,
    b: KE,
}

impl Ecw {
    fn add(&self, pp: EP, qq: EP) -> EP {
        cnt(4, 1);
        let k = &self.k;
        let (x1, y1) = match pp {
            None => return qq,
            Some(v) => v,
        };
        let (x2, y2) = match qq {
            None => return pp,
            Some(v) => v,
        };
        let lam = if x1 == x2 {
            if k.add(y1, y2) == KZ {
                return None;
            }
            k.mul(k.add(k.smul(k.sqr(x1), 3), self.a), k.inv(k.add(y1, y1)))
        } else {
            k.mul(k.sub(y2, y1), k.inv(k.sub(x2, x1)))
        };
        let x3 = k.sub(k.sub(k.sqr(lam), x1), x2);
        let y3 = k.sub(k.mul(lam, k.sub(x1, x3)), y1);
        Some((x3, y3))
    }
    fn mul(&self, mut n: u64, pt: EP) -> EP {
        let mut r: EP = None;
        let mut b = pt;
        while n > 0 {
            if n & 1 == 1 {
                r = self.add(r, b);
            }
            n >>= 1;
            if n > 0 {
                b = self.add(b, b);
            }
        }
        r
    }
}

// ---------------------------------------------------------- instance build
fn gp_eval(cmd: &str) -> String {
    let mut child = Command::new("gp")
        .args(["-q", "-f", "-s", "400000000"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("gp");
    child.stdin.take().unwrap().write_all(cmd.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    String::from_utf8(out.stdout).unwrap().trim().to_string()
}

fn gp_order_hyperell(h: &P, p: u64) -> u64 {
    let coeffs: Vec<String> = h.iter().rev().map(|c| c.to_string()).collect();
    let cmd = format!(
        "h=Pol([{}])*Mod(1,{});P=hyperellcharpoly(h);print(subst(P,'x,1))\n",
        coeffs.join(","),
        p
    );
    gp_eval(&cmd).parse().unwrap()
}

fn gp_order_e(k: &KField, a: KE, b: KE) -> u64 {
    let ks = |x: KE| format!("({}+{}*t+{}*t^2)", x[0], x[1], x[2]);
    let cmd = format!(
        "t=ffgen(Mod(1,{})*(x^3-{}*x-{}),'t);E=ellinit([{},{}]);print(ellcard(E))\n",
        k.p,
        k.c1,
        k.c0,
        ks(a),
        ks(b)
    );
    gp_eval(&cmd).parse().unwrap()
}

#[derive(Clone)]
struct Instance {
    p: u64,
    k: KField,
    a: u64,
    b: u64,
    r: [KE; 3],
    h: P,
    xn: P,
    xd: P,
    num_r: KP,
    den_r: KP,
    l3: KE,
    shift: KE,
    e: Ecw,
    n: u64,
    ell: u64,
    c: u64,
}

fn build_instance(p: u64, rng: &mut PyRng, cmax: u64) -> Option<Instance> {
    let old = snap();
    let k = KField::new(p, rng);
    let r0 = loop {
        let pr = k.rand(rng);
        if !KField::in_fp(pr) {
            break pr;
        }
    };
    let r1 = k.frob(r0);
    let r2 = k.frob(r1);
    let (a, b) = rng.sample2(p);
    let c0 = k.sub(r2, r0);
    let mut trb = [0u64; 3];
    for (idx, e) in [[1, 0, 0], [0, 1, 0], [0, 0, 1]].iter().enumerate() {
        let mut x: KE = *e;
        let mut s = KZ;
        for _ in 0..3 {
            s = k.add(s, x);
            x = k.frob(x);
        }
        assert!(KField::in_fp(s));
        trb[idx] = s[0];
    }
    let tr_ = |x: KE| (x[0] * trb[0] % p + x[1] * trb[1] % p + x[2] * trb[2] % p) % p;
    let qf = |s: KE| tr_(k.mul(c0, k.sqr(s)));
    let bf = |s1: KE, s2: KE| tr_(k.mul(c0, k.mul(s1, s2)));
    let sstar = loop {
        let s1 = k.rand(rng);
        let s2 = k.rand(rng);
        let (a2, a1, a0) = (qf(s2), 2 * bf(s1, s2) % p, qf(s1));
        if a2 == 0 {
            continue;
        }
        let disc = (a1 * a1 % p + p - 4 * a2 % p * a0 % p) % p;
        if disc == 0 || powmod(disc, (p - 1) / 2, p) != 1 {
            continue;
        }
        let sq = sqrt_fp(disc, p);
        let lam = (p - a1 + sq) % p * powmod(2 * a2 % p, p - 2, p) % p;
        let ss = k.add(s1, k.smul(s2, lam));
        if ss != KZ {
            break ss;
        }
    };
    assert_eq!(qf(sstar), 0);
    let param = |e1: KE, e2: KE| -> KP {
        let (q0, q1, q2) = (qf(e1), 2 * bf(e1, e2) % p, qf(e2));
        let (b0, b1) = (bf(sstar, e1), bf(sstar, e2));
        ktr(vec![
            k.sub(k.smul(e1, 2 * b0 % p), k.smul(sstar, q0)),
            k.sub(k.add(k.smul(e1, 2 * b1 % p), k.smul(e2, 2 * b0 % p)), k.smul(sstar, q1)),
            k.sub(k.smul(e2, 2 * b1 % p), k.smul(sstar, q2)),
        ])
    };
    struct Bf {
        x0: KP,
        x2: KP,
        delta: KP,
        kap: KE,
        xd: P,
        xn: P,
        h: P,
    }
    let build_from = |sz: &KP| -> Option<Bf> {
        let x0 = sz.clone();
        let x1: KP = x0.iter().map(|&c| k.frob(c)).collect();
        let x2: KP = x1.iter().map(|&c| k.frob(c)).collect();
        let x0s = kpmul(&k, &x0, &x0);
        let x2s = kpmul(&k, &x2, &x2);
        let delta = kpsub(&k, &x0s, &x2s);
        let num = kpsub(&k, &kpscal(&k, &x0s, r0), &kpscal(&k, &x2s, r1));
        let kap = k.inv(*delta.last().unwrap());
        let xd = kpscal(&k, &delta, kap);
        let xn = kpscal(&k, &num, kap);
        if !xd.iter().chain(xn.iter()).all(|&c| KField::in_fp(c)) {
            return None;
        }
        let xdp = tr(xd.iter().map(|c| c[0]).collect());
        let xnp = tr(xn.iter().map(|c| c[0]).collect());
        let h = pmul(&psub(&xnp, &pscal(&xdp, a, p), p), &psub(&xnp, &pscal(&xdp, b, p), p), p);
        Some(Bf { x0, x2, delta, kap, xd: xdp, xn: xnp, h })
    };
    let e1 = k.rand(rng);
    let e2 = k.rand(rng);
    let sz = param(e1, e2);
    let bf0 = build_from(&sz).expect("x(z) not F_q-rational");
    let roots: Vec<u64> = (0..p).filter(|&z| pev(&bf0.h, z, p) == 0).collect();
    if roots.is_empty() {
        restore(old);
        return None;
    }
    let mut sw: Option<KE> = None;
    for &z0 in &roots {
        let mut acc = KZ;
        for &c in sz.iter().rev() {
            acc = k.add(if acc != KZ { k.smul(acc, z0) } else { KZ }, c);
        }
        if bf(sstar, acc) != 0 {
            sw = Some(acc);
            break;
        }
    }
    let sw = match sw {
        Some(s) => s,
        None => {
            restore(old);
            return None;
        }
    };
    assert_eq!(qf(sw), 0);
    let mut found: Option<Bf> = None;
    for _ in 0..1000 {
        let e1 = k.rand(rng);
        let sz2 = param(e1, sw);
        if sz2.len() == 3 {
            if let Some(bb) = build_from(&sz2) {
                if bb.h.len() == 8 {
                    found = Some(bb);
                    break;
                }
            }
        }
    }
    let bf1 = match found {
        Some(b) => b,
        None => {
            restore(old);
            return None;
        }
    };
    let h = bf1.h.clone();
    if !is_one(&pxgcd(&h, &pder(&h, p), p).0) {
        restore(old);
        return None;
    }
    let num_r = kpscal(&k, &kpmul(&k, &bf1.delta, &bf1.delta), bf1.kap);
    let den_r = kpscal(&k, &kpmul(&k, &bf1.x0, &bf1.x2), k.sub(r0, r1));
    let rho = [r1, k.emb(a), k.emb(b)];
    let mut poly: KP = vec![KONE];
    for &rr in &rho {
        poly = kpmul(&k, &poly, &vec![KONE, k.sub(r0, rr)]);
    }
    let (l1, l2, l3) = (poly[1], poly[2], poly[3]);
    let i3 = powmod(3, p - 2, p);
    let s = k.smul(l2, i3);
    let ea = k.sub(k.mul(l1, l3), k.mul(k.smul(l2, i3), l2));
    let c1 = k.mul(l1, l3);
    let cc0 = k.sqr(l3);
    let ms = k.neg(s);
    let eb = k.add(k.add(k.add(k.mul(k.sqr(ms), ms), k.mul(l2, k.sqr(ms))), k.mul(c1, ms)), cc0);
    restore(old);
    let n = gp_order_e(&k, ea, eb);
    let mut sel = None;
    let mut cc = 4;
    while cc <= cmax {
        if n % cc == 0 && is_prime(n / cc) && (n / cc) % 2 == 1 && n / cc > cmax {
            sel = Some((n / cc, cc));
            break;
        }
        cc += 4;
    }
    let (ell, c) = sel?;
    let e = Ecw { k: k.clone(), a: ea, b: eb };
    Some(Instance {
        p,
        k,
        a,
        b,
        r: [r0, r1, r2],
        h,
        xn: bf1.xn,
        xd: bf1.xd,
        num_r,
        den_r,
        l3,
        shift: s,
        e,
        n,
        ell,
        c,
    })
}

fn w_to_quartic(inst: &Instance, pt: (KE, KE)) -> Option<(KE, KE)> {
    let k = &inst.k;
    let (x, y) = pt;
    let xp = k.sub(x, inst.shift);
    if xp == KZ {
        return None;
    }
    let xx = k.mul(xp, k.inv(inst.l3));
    let yy = k.mul(y, k.inv(inst.l3));
    let xi = k.inv(xx);
    Some((k.add(inst.r[0], xi), k.mul(yy, k.sqr(xi))))
}

fn random_point(inst: &Instance, rng: &mut PyRng) -> (KE, KE) {
    let k = &inst.k;
    let e = &inst.e;
    let old = snap();
    loop {
        let x = k.rand(rng);
        let rhs = k.add(k.add(k.mul(k.sqr(x), x), k.mul(e.a, x)), e.b);
        if let Some(y) = k.sqrt(rhs, rng) {
            restore(old);
            return (x, y);
        }
    }
}

fn transfer(inst: &Instance, pt: EP) -> Result<Div, &'static str> {
    let k = &inst.k;
    let pt = pt.ok_or("infinity")?;
    let (xp, yp) = w_to_quartic(inst, pt).ok_or("exceptional")?;
    let xn: KP = inst.xn.iter().map(|&c| [c, 0, 0]).collect();
    let xd: KP = inst.xd.iter().map(|&c| [c, 0, 0]).collect();
    let u = kpsub(k, &xn, &kpscal(k, &xd, xp));
    if u.len() != 5 {
        return Err("exceptional");
    }
    let u = kpmonic(k, &u);
    let num = kpmod(k, &inst.num_r, &u);
    let den = kpmod(k, &inst.den_r, &u);
    if den.is_empty() {
        return Err("exceptional");
    }
    let (g, s, _) = kpxgcd(k, &den, &u);
    if !(g.len() == 1 && g[0] == KONE) {
        return Err("exceptional");
    }
    let v = kpmod(k, &kpscal(k, &kpmul(k, &num, &s), yp), &u);
    let hk: KP = inst.h.iter().map(|&c| [c, 0, 0]).collect();
    assert!(kpmod(k, &kpsub(k, &kpmul(k, &v, &v), &hk), &u).is_empty(), "v^2 != h mod U");
    let d = kcreduce(k, u, v, &hk);
    let d1: (KP, KP) = (d.0.iter().map(|&c| k.frob(c)).collect(), d.1.iter().map(|&c| k.frob(c)).collect());
    let d2: (KP, KP) = (d1.0.iter().map(|&c| k.frob(c)).collect(), d1.1.iter().map(|&c| k.frob(c)).collect());
    let s = kcadd(k, &kcadd(k, &d, &d1, &hk), &d2, &hk);
    assert!(s.0.iter().chain(s.1.iter()).all(|&c| KField::in_fp(c)), "phi(P) not F_q-rational");
    Ok((tr(s.0.iter().map(|c| c[0]).collect()), tr(s.1.iter().map(|c| c[0]).collect())))
}

// --------------------------------------------------------------- IC helpers
fn sqrmod3(a: [u64; 3], u: &P, p: u64) -> [u64; 3] {
    let (a0, a1, a2) = (a[0] as i128, a[1] as i128, a[2] as i128);
    let mut d0 = a0 * a0;
    let mut d1 = 2 * a0 * a1;
    let mut d2 = 2 * a0 * a2 + a1 * a1;
    let mut d3 = 2 * a1 * a2;
    let d4 = a2 * a2;
    let (u0, u1, u2) = (u[0] as i128, u[1] as i128, u[2] as i128);
    d3 -= d4 * u2;
    d2 -= d4 * u1;
    d1 -= d4 * u0;
    d2 -= d3 * u2;
    d1 -= d3 * u1;
    d0 -= d3 * u0;
    cnt(0, 12);
    let pp = p as i128;
    [d0.rem_euclid(pp) as u64, d1.rem_euclid(pp) as u64, d2.rem_euclid(pp) as u64]
}

fn xmul3(a: [u64; 3], u: &P, p: u64) -> [u64; 3] {
    cnt(0, 3);
    let pp = p as i128;
    let (a0, a1, a2) = (a[0] as i128, a[1] as i128, a[2] as i128);
    [
        (-a2 * u[0] as i128).rem_euclid(pp) as u64,
        (a0 - a2 * u[1] as i128).rem_euclid(pp) as u64,
        (a1 - a2 * u[2] as i128).rem_euclid(pp) as u64,
    ]
}

fn splits(u: &P, p: u64) -> bool {
    let d = u.len() - 1;
    if d <= 1 {
        return true;
    }
    if d == 3 {
        let mut r = [0u64, 1, 0];
        let bl = bitlen(p);
        for i in (0..bl - 1).rev() {
            r = sqrmod3(r, u, p);
            if (p >> i) & 1 == 1 {
                r = xmul3(r, u, p);
            }
        }
        return r == [0, 1, 0];
    }
    let disc = (u[1] * u[1] % p + p - 4 * u[0] % p) % p;
    cnt(0, 2);
    cnt(0, 3 * bitlen(p) / 2);
    disc != 0 && powmod(disc, (p - 1) / 2, p) == 1
}

fn ppowmod(base: &P, mut e: u64, m: &P, p: u64) -> P {
    let mut r: P = vec![1];
    let mut b = pmod(base, m, p);
    while e > 0 {
        if e & 1 == 1 {
            r = pmod(&pmul(&r, &b, p), m, p);
        }
        e >>= 1;
        if e > 0 {
            b = pmod(&pmul(&b, &b, p), m, p);
        }
    }
    r
}

fn roots_split(u: &P, p: u64, rng: &mut PyRng) -> Vec<u64> {
    let d = u.len() - 1;
    if d == 0 {
        return vec![];
    }
    if d == 1 {
        return vec![(p - u[0]) % p];
    }
    if d == 2 {
        let disc = (u[1] * u[1] % p + p - 4 * u[0] % p) % p;
        let sq = sqrt_fp(disc, p);
        cnt(0, 2 + 3 * bitlen(p) / 2 + 2);
        cnt(1, 1);
        let i2 = (p + 1) / 2;
        return vec![(2 * p - u[1] + sq) % p * i2 % p, (2 * p - u[1] - sq) % p * i2 % p];
    }
    loop {
        let delta = rng.randbelow(p);
        let w = ppowmod(&vec![delta, 1], (p - 1) / 2, u, p);
        let g = pxgcd(u, &psub(&w, &vec![1], p), p).0;
        if !g.is_empty() && g.len() - 1 > 0 && g.len() - 1 < d {
            let g2 = pdivmod(u, &g, p).0;
            let mut r = roots_split(&g, p, rng);
            r.extend(roots_split(&g2, p, rng));
            return r;
        }
    }
}

// ------------------------------------------------------- arithmetic mod ell
#[derive(Clone, Copy)]
struct Mont {
    m: u64,
    minv: u64, // -m^{-1} mod 2^64
    r2: u64,   // 2^128 mod m
}

impl Mont {
    fn new(m: u64) -> Self {
        assert!(m % 2 == 1 && m < (1u64 << 62));
        let mut inv: u64 = 1;
        for _ in 0..7 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(m.wrapping_mul(inv)));
        }
        let r1 = ((1u128 << 64) % m as u128) as u64;
        let r2 = ((r1 as u128 * r1 as u128) % m as u128) as u64;
        Mont { m, minv: inv.wrapping_neg(), r2 }
    }
    #[inline(always)]
    fn redc(&self, t: u128) -> u64 {
        let u = (t as u64).wrapping_mul(self.minv);
        let s = ((t + u as u128 * self.m as u128) >> 64) as u64;
        if s >= self.m {
            s - self.m
        } else {
            s
        }
    }
    #[inline(always)]
    fn to_mont(&self, a: u64) -> u64 {
        self.redc(a as u128 * self.r2 as u128)
    }
    /// a (normal) times bm (Montgomery form of b) = a*b mod m.
    #[inline(always)]
    fn mulm(&self, a: u64, bm: u64) -> u64 {
        self.redc(a as u128 * bm as u128)
    }
    /// full reduction of any u128.
    #[inline(always)]
    fn reduce(&self, t: u128) -> u64 {
        let mut hi = (t >> 64) as u64;
        let lo = t as u64;
        if hi >= self.m {
            hi %= self.m;
        }
        let x = self.redc(((hi as u128) << 64) | lo as u128);
        self.redc(x as u128 * self.r2 as u128)
    }
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

#[derive(Clone, Copy)]
struct Ent {
    col: u32,
    s: i64,  // small signed value, or i64::MIN when big
    bm: u64, // Montgomery form of the value when big
}

fn lanczos(rows: &[Vec<(u32, u64)>], rhs: &[u64], ncols: usize, ell: u64) -> Vec<u64> {
    let mt = Mont::new(ell);
    let lim: u64 = 1 << 40;
    let mat: Vec<Vec<Ent>> = rows
        .iter()
        .map(|r| {
            r.iter()
                .map(|&(j, c)| {
                    if c < lim {
                        Ent { col: j, s: c as i64, bm: 0 }
                    } else if ell - c < lim {
                        Ent { col: j, s: -((ell - c) as i64), bm: 0 }
                    } else {
                        Ent { col: j, s: i64::MIN, bm: mt.to_mont(c) }
                    }
                })
                .collect()
        })
        .collect();
    let nnz: u64 = rows.iter().map(|r| r.len() as u64).sum();
    let ncu = ncols as u64;
    // pos - neg mod ell with one reduction (both below ell * 2^64 by the row-weight bound)
    let big: u128 = (ell as u128) << 64;
    let fin = move |pos: u128, neg: u128| -> u64 {
        if pos < (1u128 << 126) && neg <= big {
            mt.reduce(pos + (big - neg))
        } else {
            let (a, b) = (mt.reduce(pos), mt.reduce(neg));
            if a >= b {
                a - b
            } else {
                a + ell - b
            }
        }
    };
    let sv = |v: &[u64]| -> Vec<u64> {
        cnt(2, nnz);
        mat.iter()
            .map(|r| {
                let (mut pos, mut neg) = (0u128, 0u128);
                for e in r {
                    let x = v[e.col as usize];
                    match e.s {
                        1 => pos += x as u128,
                        -1 => neg += x as u128,
                        i64::MIN => pos += mt.mulm(x, e.bm) as u128,
                        s if s >= 0 => pos += s as u128 * x as u128,
                        s => neg += (-s) as u128 * x as u128,
                    }
                }
                fin(pos, neg)
            })
            .collect()
    };
    let stv = |w: &[u64]| -> Vec<u64> {
        cnt(2, nnz);
        let mut acc = vec![(0u128, 0u128); ncols];
        for (r, &wi) in mat.iter().zip(w.iter()) {
            if wi != 0 {
                for e in r {
                    let a = &mut acc[e.col as usize];
                    match e.s {
                        1 => a.0 += wi as u128,
                        -1 => a.1 += wi as u128,
                        i64::MIN => a.0 += mt.mulm(wi, e.bm) as u128,
                        s if s >= 0 => a.0 += s as u128 * wi as u128,
                        s => a.1 += (-s) as u128 * wi as u128,
                    }
                }
            }
        }
        acc.iter().map(|&(a, b)| fin(a, b)).collect()
    };
    let dot = |a: &[u64], b: &[u64]| -> u64 {
        cnt(2, ncu);
        let mut acc: u128 = 0;
        let mut res: u64 = 0;
        for (i, (&x, &y)) in a.iter().zip(b.iter()).enumerate() {
            acc += x as u128 * y as u128;
            if i & 63 == 63 {
                res = (res + mt.reduce(acc)) % ell;
                acc = 0;
            }
        }
        (res + mt.reduce(acc)) % ell
    };
    let inv = |x: u64| -> u64 {
        cnt(3, 1);
        powmod(x, ell - 2, ell)
    };
    let b = stv(rhs);
    let mut x = vec![0u64; ncols];
    let mut w_prev: Option<(Vec<u64>, Vec<u64>, u64)> = None;
    let mut w = b.clone();
    for _it in 0..ncols + 5 {
        let aw = stv(&sv(&w));
        let waw = dot(&w, &aw);
        if waw == 0 {
            break;
        }
        let iw = inv(waw);
        let coef = mulmod(dot(&w, &b), iw, ell);
        let cm = mt.to_mont(coef);
        for (xi, &wi) in x.iter_mut().zip(w.iter()) {
            let t = *xi + mt.mulm(wi, cm);
            *xi = if t >= ell { t - ell } else { t };
        }
        cnt(2, ncu);
        let awaw = dot(&aw, &aw);
        let alpha = mulmod(awaw, iw, ell);
        let am = mt.to_mont(alpha);
        let mut nw: Vec<u64> = aw
            .iter()
            .zip(w.iter())
            .map(|(&a, &wi)| {
                let t = mt.mulm(wi, am);
                if a >= t {
                    a - t
                } else {
                    a + ell - t
                }
            })
            .collect();
        cnt(2, ncu);
        if let Some((wp, awp, wawp)) = &w_prev {
            let beta = mulmod(dot(&aw, awp), inv(*wawp), ell);
            let bmn = mt.to_mont(beta);
            for (n, &wi) in nw.iter_mut().zip(wp.iter()) {
                let t = mt.mulm(wi, bmn);
                *n = if *n >= t { *n - t } else { *n + ell - t };
            }
            cnt(2, ncu);
        }
        let done = nw.iter().all(|&v| v == 0);
        w_prev = Some((w, aw, waw));
        w = nw;
        if done {
            break;
        }
    }
    x
}

// ---------------------------------------------------- sparse rows mod ell
type Row = Vec<(u32, u64)>; // sorted by column

fn sub_expr(d: &Row, rhs: u64, s: i8, e: &(Row, u64), ell: u64) -> (Row, u64) {
    // d + s*lambda = rhs, lambda = e.0 . x + e.1  ->  d + s e.0 = rhs - s e.1
    let mut out = Vec::with_capacity(d.len() + e.0.len());
    let (mut i, mut j) = (0, 0);
    let sc = |c: u64| if s == 1 { c } else { (ell - c) % ell };
    while i < d.len() || j < e.0.len() {
        if j >= e.0.len() || (i < d.len() && d[i].0 < e.0[j].0) {
            out.push(d[i]);
            i += 1;
        } else if i >= d.len() || e.0[j].0 < d[i].0 {
            let v = sc(e.0[j].1);
            if v != 0 {
                out.push((e.0[j].0, v));
            }
            j += 1;
        } else {
            let v = (d[i].1 + sc(e.0[j].1)) % ell;
            if v != 0 {
                out.push((d[i].0, v));
            }
            i += 1;
            j += 1;
        }
    }
    let r = (rhs + ell - sc(e.1)) % ell;
    (out, r)
}

fn solve_for(d: &Row, rhs: u64, s: i8, ell: u64) -> (Row, u64) {
    // d + s*lambda = rhs  ->  lambda = s*(rhs - d)
    if s == 1 {
        (d.iter().map(|&(j, c)| (j, (ell - c) % ell)).collect(), rhs)
    } else {
        (d.clone(), (ell - rhs) % ell)
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Variant {
    Basic,
    Lp2,
    Lp1,
}

#[derive(Default, Clone, Debug)]
struct IcStats {
    steps: u64,
    cand_split: u64,
    full: u64,
    lp1: u64,
    lp2: u64,
    lp3: u64,
    restarts: u64,
    cycles: u64,
    weier: u64,
}

struct IcOut {
    k: Option<u64>,
    nfb: usize,
    rows: usize,
    lp_vertices: usize,
    checks: u64,
    counters: [u64; 5],
    setup_counters: [u64; 5],
    la_counters: [u64; 5],
    la_ncols: usize,
    la_nnz: u64,
    st: IcStats,
    wall: f64,
}

fn key_of(u: &P, v: &P) -> u128 {
    let mut k: u128 = (u.len() as u128) & 3;
    for i in 0..3 {
        k = (k << 21) | (*u.get(i).unwrap_or(&0) as u128);
    }
    for i in 0..3 {
        k = (k << 21) | (*v.get(i).unwrap_or(&0) as u128);
    }
    k
}

fn diff(a: [u64; 5], b: [u64; 5]) -> [u64; 5] {
    let mut r = [0; 5];
    for i in 0..5 {
        r[i] = a[i] - b[i];
    }
    r
}

fn ic_solve(inst: &Instance, d1: &Div, d2: &Div, rng: &mut PyRng, variant: Variant, kappa: f64) -> IcOut {
    let (p, h, ell) = (inst.p, &inst.h, inst.ell);
    let t0 = Instant::now();
    let c_before = snap();
    let nfb_target: Option<usize> = match variant {
        Variant::Basic => None,
        Variant::Lp2 => Some((kappa * (p as f64).powf(2.0 / 3.0)).ceil() as usize),
        Variant::Lp1 => Some((kappa * (p as f64).powf(5.0 / 7.0)).ceil() as usize),
    };
    let mut fbidx: Vec<i32> = vec![-1; p as usize];
    let mut nfb = 0usize;
    let blcost = 3 * bitlen(p) / 2;
    for x in 0..p {
        let hx = pev(h, x, p);
        cnt(0, blcost);
        if hx != 0 && powmod(hx, (p - 1) / 2, p) == 1 {
            match nfb_target {
                None => {
                    fbidx[x as usize] = nfb as i32;
                    nfb += 1;
                }
                Some(t) if nfb < t => {
                    fbidx[x as usize] = nfb as i32;
                    nfb += 1;
                }
                Some(_) => break,
            }
        }
    }
    let kcol = nfb as u32;
    let use_lp = variant != Variant::Basic;
    let r = 20usize;
    let a_s = rng.randrange(1, ell);
    let b_s = rng.randrange(1, ell);
    let m0 = cmul(a_s, d1, h, p);
    let m1 = cmul(b_s, d1, h, p);
    let mut ms = vec![m0, m1];
    let mut ma = vec![a_s, b_s];
    for _ in 2..r {
        let n = cadd(&ms[ms.len() - 1], &ms[ms.len() - 2], h, p);
        ms.push(n);
        ma.push((ma[ma.len() - 1] + ma[ma.len() - 2]) % ell);
    }
    let mut al = rng.randbelow(ell);
    let mut be = rng.randrange(1, ell);
    let mut rr = cadd(&cmul(al, d1, h, p), &cmul(be, d2, h, p), h, p);
    let c_setup = snap();
    let mut seen: HashSet<u128> = HashSet::new();
    let mut full_rows: Vec<(Row, u64)> = vec![];
    let mut expr: HashMap<u64, (Row, u64)> = HashMap::new();
    let mut rels: Vec<Rel> = vec![];
    let mut pend: HashMap<u64, Vec<(u64, usize)>> = HashMap::new();
    let mut st = IcStats::default();
    let mut next_check = nfb + 10;
    let mut kfound: Option<u64> = None;
    let mut last_new = 0u64;
    let mut checks = 0u64;
    let mut la_counters = [0u64; 5];
    let (mut la_ncols, mut la_nnz) = (0usize, 0u64);

    fn emit(full_rows: &mut Vec<(Row, u64)>, d: Row, rhs: u64) {
        let d: Row = d.into_iter().filter(|&(_, c)| c != 0).collect();
        full_rows.push((d, rhs));
    }

    while kfound.is_none() {
        st.steps += 1;
        let (u, v) = (&rr.0, &rr.1);
        if splits(u, p) {
            let key = key_of(u, v);
            if seen.contains(&key) {
                st.restarts += 1;
                rr = cadd(&cadd(&rr, &ms[0], h, p), d2, h, p);
                al = (al + ma[0]) % ell;
                be = (be + 1) % ell;
                last_new = st.steps;
                continue;
            }
            seen.insert(key);
            last_new = st.steps;
            st.cand_split += 1;
            let rts = roots_split(u, p, rng);
            let mut rowd: Row = vec![];
            let mut lps: Vec<(u64, i8)> = vec![];
            let mut ok = true;
            for &x0 in &rts {
                let y0 = pev(v, x0, p);
                if y0 == 0 {
                    st.weier += 1;
                    continue;
                }
                let s: i8 = if y0 <= p - y0 { 1 } else { -1 };
                let j = fbidx[x0 as usize];
                if j >= 0 {
                    rowd.push((j as u32, if s == 1 { 1 } else { ell - 1 }));
                } else {
                    if !use_lp {
                        ok = false;
                        break;
                    }
                    lps.push((x0, s));
                }
            }
            if ok {
                rowd.sort();
                rowd.push((kcol, (ell - be % ell) % ell));
                let nl = lps.len();
                if nl == 0 {
                    st.full += 1;
                    emit(&mut full_rows, rowd, al);
                } else if nl == 1 {
                    st.lp1 += 1;
                    let (xl, sl) = lps[0];
                    if let Some(e) = expr.get(&xl) {
                        let (d, r2) = sub_expr(&rowd, al, sl, e, ell);
                        emit(&mut full_rows, d, r2);
                        st.cycles += 1;
                    } else {
                        expr.insert(xl, solve_for(&rowd, al, sl, ell));
                        propagate(xl, &mut expr, &mut pend, &rels, &mut full_rows, &mut st, ell);
                    }
                } else if nl == 2 {
                    st.lp2 += 1;
                    if variant == Variant::Lp2 {
                        let (mut xa, mut sa) = lps[0];
                        let (mut xb, mut sb) = lps[1];
                        let ina = expr.contains_key(&xa);
                        let inb = expr.contains_key(&xb);
                        if ina && inb {
                            let (d, r2) = sub_expr(&rowd, al, sa, &expr[&xa], ell);
                            let (d, r2) = sub_expr(&d, r2, sb, &expr[&xb], ell);
                            emit(&mut full_rows, d, r2);
                            st.cycles += 1;
                        } else if ina || inb {
                            if inb {
                                std::mem::swap(&mut xa, &mut xb);
                                std::mem::swap(&mut sa, &mut sb);
                            }
                            let (d, r2) = sub_expr(&rowd, al, sa, &expr[&xa], ell);
                            expr.insert(xb, solve_for(&d, r2, sb, ell));
                            propagate(xb, &mut expr, &mut pend, &rels, &mut full_rows, &mut st, ell);
                        } else {
                            let id = rels.len();
                            rels.push(Rel { row: rowd, rhs: al, lps: [(xa, sa), (xb, sb)] });
                            pend.entry(xa).or_default().push((xb, id));
                            pend.entry(xb).or_default().push((xa, id));
                        }
                    }
                } else {
                    st.lp3 += 1;
                }
            }
        }
        if st.steps - last_new > 60 {
            st.restarts += 1;
            rr = cadd(&cadd(&rr, &ms[0], h, p), d2, h, p);
            al = (al + ma[0]) % ell;
            be = (be + 1) % ell;
            last_new = st.steps;
            continue;
        }
        if full_rows.len() >= next_check {
            let nrows_all = full_rows.len();
            let mut alive = vec![true; nrows_all];
            loop {
                let mut changed = false;
                let mut cntc = vec![0u32; nfb + 1];
                for (i, (rw, _)) in full_rows.iter().enumerate() {
                    if alive[i] {
                        for &(j, _) in rw {
                            cntc[j as usize] += 1;
                        }
                    }
                }
                for (i, (rw, _)) in full_rows.iter().enumerate() {
                    if alive[i] && rw.iter().any(|&(j, _)| cntc[j as usize] == 1 && j != kcol) {
                        alive[i] = false;
                        changed = true;
                    }
                }
                if !changed {
                    break;
                }
            }
            let mut cols: Vec<i64> = vec![-1; nfb + 1];
            let mut ncols = 0usize;
            for (i, (rw, _)) in full_rows.iter().enumerate() {
                if alive[i] {
                    for &(j, _) in rw {
                        if cols[j as usize] < 0 {
                            cols[j as usize] = ncols as i64;
                            ncols += 1;
                        }
                    }
                }
            }
            let nrows = alive.iter().filter(|&&a| a).count();
            if cols[kcol as usize] >= 0 && nrows >= ncols + 5 {
                checks += 1;
                let mut srows = Vec::with_capacity(nrows);
                let mut srhs = Vec::with_capacity(nrows);
                for (i, (rw, rh)) in full_rows.iter().enumerate() {
                    if alive[i] {
                        srows.push(rw.iter().map(|&(j, c)| (cols[j as usize] as u32, c)).collect::<Vec<_>>());
                        srhs.push(*rh);
                    }
                }
                let c_la0 = snap();
                let sol = lanczos(&srows, &srhs, ncols, ell);
                let c_la1 = snap();
                let dl = diff(c_la1, c_la0);
                for i in 0..5 {
                    la_counters[i] += dl[i];
                }
                la_ncols = ncols;
                la_nnz = srows.iter().map(|r| r.len() as u64).sum();
                let kk = sol[cols[kcol as usize] as usize];
                let old = snap();
                let okk = cmul(kk, d1, h, p) == *d2;
                restore(old);
                if okk {
                    kfound = Some(kk);
                    break;
                }
            }
            next_check = (full_rows.len() as f64 * 1.1) as usize + 5;
        }
        if kfound.is_some() {
            break;
        }
        let (u, v) = (&rr.0, &rr.1);
        let u1 = if u.len() > 1 { u[1] } else { 0 };
        let v0 = if !v.is_empty() { v[0] } else { 0 };
        let j = ((u[0] * 31 + u1 * 7 + v0) % r as u64) as usize;
        rr = cadd(&rr, &ms[j], h, p);
        al = (al + ma[j]) % ell;
    }
    let c_end = snap();
    IcOut {
        k: kfound,
        nfb,
        rows: full_rows.len(),
        lp_vertices: expr.len(),
        checks,
        counters: diff(c_end, c_before),
        setup_counters: diff(c_setup, c_before),
        la_counters,
        la_ncols,
        la_nnz,
        st,
        wall: t0.elapsed().as_secs_f64(),
    }
}

fn propagate(
    start: u64,
    expr: &mut HashMap<u64, (Row, u64)>,
    pend: &mut HashMap<u64, Vec<(u64, usize)>>,
    rels: &[Rel],
    full_rows: &mut Vec<(Row, u64)>,
    st: &mut IcStats,
    ell: u64,
) {
    let mut stack = vec![start];
    while let Some(v) = stack.pop() {
        let list = pend.remove(&v).unwrap_or_default();
        for (w, id) in list {
            let rel = &rels[id];
            let mut d = rel.row.clone();
            let mut rr = rel.rhs;
            for &(xl, sl) in rel.lps.iter() {
                if xl == v {
                    let (d2, r2) = sub_expr(&d, rr, sl, &expr[&v], ell);
                    d = d2;
                    rr = r2;
                }
            }
            let sw = rel.lps.iter().find(|t| t.0 == w).unwrap().1;
            if let Some(e) = expr.get(&w) {
                let (d2, r2) = sub_expr(&d, rr, sw, e, ell);
                let d2: Row = d2.into_iter().filter(|&(_, c)| c != 0).collect();
                full_rows.push((d2, r2));
                st.cycles += 1;
            } else {
                expr.insert(w, solve_for(&d, rr, sw, ell));
                stack.push(w);
            }
            if let Some(lst) = pend.get_mut(&w) {
                lst.retain(|t| t.1 != id);
            }
        }
    }
}

struct Rel {
    row: Row,
    rhs: u64,
    lps: [(u64, i8); 2],
}

// ------------------------------------------------------------------- rho
fn canon(k: &KField, pt: (KE, KE)) -> ((KE, KE), bool) {
    let ny = k.neg(pt.1);
    // canonical: lexicographically smaller of (y2, y1, y0) and (-y2, -y1, -y0)
    let key = |y: KE| (y[2], y[1], y[0]);
    if key(ny) < key(pt.1) {
        ((pt.0, ny), true)
    } else {
        (pt, false)
    }
}

struct RhoOut {
    k: Option<u64>,
    steps: u64,
    counters: [u64; 5],
    setup_counters: [u64; 5],
    restarts: u64,
    escapes: u64,
    wall: f64,
}

fn rho_neg(inst: &Instance, g: EP, t: EP, rng: &mut PyRng, r: usize, max_steps: u64) -> RhoOut {
    let e = &inst.e;
    let k = &inst.k;
    let ell = inst.ell;
    let t0 = Instant::now();
    let c0 = snap();
    let mut coeffs = vec![];
    let mut mpts = vec![];
    for _ in 0..r {
        let a = rng.randbelow(ell);
        let b = rng.randrange(1, ell);
        let pnt = e.add(e.mul(a, g), e.mul(b, t));
        coeffs.push((a, b));
        mpts.push(pnt.expect("multiplier at infinity"));
    }
    let c_setup = snap();
    let sq = (ell as f64).sqrt();
    let dpbits = ((sq.log2() as i64) - 6).max(2) as u32;
    let mask: u64 = (1u64 << dpbits) - 1;
    let hx = |x: KE| x[0].wrapping_mul(1000003).wrapping_add(x[1].wrapping_mul(7919)).wrapping_add(x[2]);
    let ordk = |x: KE| (x[2], x[1], x[0]);
    // one walk step on classes {P, -P}: returns the canonical successor and its coefficients
    let step = |cur: (KE, KE), a: u64, b: u64| -> Option<((KE, KE), u64, u64)> {
        let j = ((hx(cur.0) >> 1) % r as u64) as usize;
        let nxt = e.add(Some(cur), Some(mpts[j]))?;
        let (c2, n2) = canon(k, nxt);
        let (mut na, mut nb) = ((a + coeffs[j].0) % ell, (b + coeffs[j].1) % ell);
        if n2 {
            na = (ell - na) % ell;
            nb = (ell - nb) % ell;
        }
        Some((c2, na, nb))
    };
    let mut seen: HashMap<(KE, KE), (u64, u64)> = HashMap::new();
    let mut steps = 0u64;
    let (mut restarts, mut escapes) = (0u64, 0u64);
    let maxwalk = 40u64 << dpbits;
    const L: u64 = 64; // fruitless cycles up to length L are caught within 2L steps
    'outer: loop {
        let a0 = rng.randbelow(ell);
        let b0 = rng.randrange(1, ell);
        let Some(st) = e.add(e.mul(a0, g), e.mul(b0, t)) else { continue };
        let (mut cur, neg) = canon(k, st);
        let (mut a, mut b) = if neg { ((ell - a0) % ell, (ell - b0) % ell) } else { (a0, b0) };
        let mut walk = 0u64;
        let mut saved: Option<(KE, KE)> = None;
        loop {
            if (hx(cur.0) >> 3) & mask == 0 {
                if let Some(&(a2, b2)) = seen.get(&cur) {
                    if (b2 + ell - b) % ell != 0 {
                        let kk = mulmod((a + ell - a2) % ell, powmod((b2 + ell - b) % ell, ell - 2, ell), ell);
                        let c_end = snap();
                        return RhoOut {
                            k: Some(kk),
                            steps,
                            counters: diff(c_end, c0),
                            setup_counters: diff(c_setup, c0),
                            restarts,
                            escapes,
                            wall: t0.elapsed().as_secs_f64(),
                        };
                    }
                    restarts += 1;
                    continue 'outer;
                }
                seen.insert(cur, (a, b));
            }
            if steps >= max_steps {
                break 'outer;
            }
            if walk > maxwalk {
                restarts += 1;
                continue 'outer;
            }
            if saved == Some(cur) {
                // fruitless cycle: walk it once, escape by doubling its smallest point
                escapes += 1;
                let (mut x, mut xa, mut xb) = (cur, a, b);
                let (mut mn, mut ma_, mut mb) = (cur, a, b);
                loop {
                    let Some((y, ya, yb)) = step(x, xa, xb) else { continue 'outer };
                    steps += 1;
                    walk += 1;
                    if y == cur {
                        break;
                    }
                    if ordk(y.0) < ordk(mn.0) {
                        mn = y;
                        ma_ = ya;
                        mb = yb;
                    }
                    x = y;
                    xa = ya;
                    xb = yb;
                }
                let Some(dbl) = e.add(Some(mn), Some(mn)) else { continue 'outer };
                steps += 1;
                walk += 1;
                let (c2, n2) = canon(k, dbl);
                cur = c2;
                a = 2 * ma_ % ell;
                b = 2 * mb % ell;
                if n2 {
                    a = (ell - a) % ell;
                    b = (ell - b) % ell;
                }
                saved = None;
                continue;
            }
            if walk % L == 0 {
                saved = Some(cur);
            }
            let Some((c2, na, nb)) = step(cur, a, b) else { continue 'outer };
            steps += 1;
            walk += 1;
            cur = c2;
            a = na;
            b = nb;
        }
    }
    let c_end = snap();
    RhoOut {
        k: None,
        steps,
        counters: diff(c_end, c0),
        setup_counters: diff(c_setup, c0),
        restarts,
        escapes,
        wall: t0.elapsed().as_secs_f64(),
    }
}

// ------------------------------------------------------------------ driver
fn cost_units(c: &[u64; 5], p: u64, cl: u64) -> u64 {
    c[0] + bitlen(p - 1) * c[1] + cl * c[2]
}

struct Args {
    bits: u32,
    curve: u64,
    arms: Vec<String>,
    mode: String,
    kappa1: f64,
    kappa2: f64,
    out: String,
    rho_targets: usize,
    rho_r: usize,
    gp_check: bool,
}

fn parse_args() -> Args {
    let mut a = Args {
        bits: 8,
        curve: 0,
        arms: vec!["ic2lp".into()],
        mode: "paired".into(),
        kappa1: 1.0,
        kappa2: 2.0,
        out: "out.jsonl".into(),
        rho_targets: 1,
        rho_r: 64,
        gp_check: true,
    };
    let v: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < v.len() {
        let val = v.get(i + 1).cloned().unwrap_or_default();
        match v[i].as_str() {
            "--bits" => a.bits = val.parse().unwrap(),
            "--curve" => a.curve = val.parse().unwrap(),
            "--arms" => a.arms = val.split(',').map(|s| s.to_string()).collect(),
            "--mode" => a.mode = val,
            "--kappa1" => a.kappa1 = val.parse().unwrap(),
            "--kappa2" => a.kappa2 = val.parse().unwrap(),
            "--out" => a.out = val,
            "--rho-targets" => a.rho_targets = val.parse().unwrap(),
            "--rho-r" => a.rho_r = val.parse().unwrap(),
            "--no-gp-check" => {
                a.gp_check = false;
                i += 1;
                continue;
            }
            x => panic!("unknown arg {x}"),
        }
        i += 2;
    }
    a
}

fn ic_json(arm: &str, kappa: f64, res: &IcOut, ktrue: u64, ctr: [u64; 5], p: u64, wall_total: f64) -> serde_json::Value {
    let mut tot = [0u64; 5];
    for i in 0..5 {
        tot[i] = ctr[i] + res.counters[i];
    }
    serde_json::json!({
        "arm": arm, "kappa": kappa, "k": res.k, "k_true": ktrue, "correct": res.k == Some(ktrue),
        "nfb": res.nfb, "rows": res.rows, "lp_vertices": res.lp_vertices, "checks": res.checks,
        "wall": res.wall, "counters": res.counters, "setup_counters": res.setup_counters,
        "la_counters": res.la_counters, "la_ncols": res.la_ncols, "la_nnz": res.la_nnz,
        "steps": res.st.steps, "cand_split": res.st.cand_split, "full": res.st.full,
        "lp1": res.st.lp1, "lp2": res.st.lp2, "lp3": res.st.lp3, "restarts": res.st.restarts,
        "cycles": res.st.cycles, "weier": res.st.weier, "transfer_counters": ctr,
        "total_counters": tot, "cost": cost_units(&tot, p, 1), "cost_cl9": cost_units(&tot, p, 9),
        "wall_total": wall_total,
    })
}

fn main() {
    let args = parse_args();
    let bits = args.bits;
    let p = next_prime((1u64 << bits) + 1);
    let seed = 20261014u64 * 1000 + bits as u64 * 10 + args.curve;
    let mut rng = PyRng::from_seed(seed);
    let t0 = Instant::now();
    let mut tries = 0;
    let inst = loop {
        tries += 1;
        if let Some(i) = build_instance(p, &mut rng, 64) {
            break i;
        }
    };
    let build_wall = t0.elapsed().as_secs_f64();
    let e = inst.e.clone();
    let g: EP = loop {
        let pt = random_point(&inst, &mut rng);
        let old = snap();
        let gg = e.mul(inst.c, Some(pt));
        restore(old);
        if gg.is_some() {
            break gg;
        }
    };
    let old = snap();
    assert!(e.mul(inst.ell, g).is_none());
    let phig = transfer(&inst, g).expect("phi(G)");
    let id: Div = (vec![1], vec![]);
    let mut chk = serde_json::Map::new();
    chk.insert("phiG_nonzero".into(), (phig != id).into());
    chk.insert("ell_phiG_zero".into(), (cmul(inst.ell, &phig, &inst.h, p) == id).into());
    chk.insert("N_phiG_zero".into(), (cmul(inst.n, &phig, &inst.h, p) == id).into());
    let mut hom = 0;
    for _ in 0..3 {
        let p1 = random_point(&inst, &mut rng);
        let p2 = random_point(&inst, &mut rng);
        let p3 = e.add(Some(p1), Some(p2));
        if let (Ok(a), Ok(b), Ok(c)) = (transfer(&inst, p3), transfer(&inst, Some(p1)), transfer(&inst, Some(p2))) {
            if a == cadd(&b, &c, &inst.h, p) {
                hom += 1;
            }
        }
    }
    chk.insert("hom_ok".into(), hom.into());
    if bits <= 11 && args.gp_check {
        chk.insert("jac_order_hyperellcharpoly_equals_N".into(), (gp_order_hyperell(&inst.h, p) == inst.n).into());
    }
    chk.insert("deg_h".into(), (inst.h.len() - 1).into());
    restore(old);
    let mut results = vec![];
    let run_ic = |arm: &str, rng: &mut PyRng, k: u64, t: EP| -> serde_json::Value {
        let (variant, kappa) = match arm {
            "basic" => (Variant::Basic, 0.0),
            "ic2lp" | "aa" => (Variant::Lp2, 2.0),
            "ic2lpt" => (Variant::Lp2, args.kappa2),
            "ic1lp" => (Variant::Lp1, args.kappa1),
            _ => panic!("arm {arm}"),
        };
        reset();
        let tt = Instant::now();
        let dd1 = transfer(&inst, g).unwrap();
        let dd2 = transfer(&inst, t).unwrap();
        let ctr = snap();
        let res = ic_solve(&inst, &dd1, &dd2, rng, variant, kappa);
        let j = ic_json(arm, kappa, &res, k, ctr, p, tt.elapsed().as_secs_f64());
        eprintln!(
            "{{\"bits\":{},\"curve\":{},\"arm\":\"{}\",\"correct\":{},\"cost\":{},\"steps\":{},\"wall\":{:.1}}}",
            bits, args.curve, arm, j["correct"], j["cost"], j["steps"], j["wall_total"].as_f64().unwrap()
        );
        j
    };
    if args.mode == "python" {
        for arm in args.arms.iter().filter(|a| *a == "ic2lp" || *a == "basic") {
            let k = rng.randrange(1, inst.ell);
            let old = snap();
            let t = e.mul(k, g);
            restore(old);
            results.push(run_ic(arm, &mut rng, k, t));
        }
    } else {
        let k = rng.randrange(1, inst.ell);
        let old = snap();
        let t = e.mul(k, g);
        restore(old);
        let base = rng.clone();
        for arm in &args.arms {
            if arm == "rho" {
                continue;
            }
            let mut r2 = base.clone();
            if arm == "aa" {
                r2.genrand();
            }
            results.push(run_ic(arm, &mut r2, k, t));
        }
        if args.arms.iter().any(|a| a == "rho") {
            let mut r2 = base.clone();
            for ti in 0..args.rho_targets {
                // target 0 is the shared k; later targets are fresh
                let kk = if ti == 0 { k } else { r2.randrange(1, inst.ell) };
                let old = snap();
                let tt = e.mul(kk, g);
                restore(old);
                reset();
                let max_steps = 40 * (inst.ell as f64).sqrt() as u64 + 1000;
                let ro = rho_neg(&inst, g, tt, &mut r2, args.rho_r, max_steps);
                let cost = cost_units(&ro.counters, p, 1);
                let walk_cost = cost - cost_units(&ro.setup_counters, p, 1);
                let per_step = walk_cost as f64 / ro.steps.max(1) as f64;
                let model = per_step * (std::f64::consts::PI * inst.ell as f64 / 4.0).sqrt();
                let j = serde_json::json!({
                    "arm": "rho", "target": ti, "k": ro.k, "k_true": kk, "correct": ro.k == Some(kk),
                    "steps": ro.steps, "counters": ro.counters, "setup_counters": ro.setup_counters,
                    "restarts": ro.restarts, "escapes": ro.escapes, "r": args.rho_r,
                    "cost": cost, "cost_cl9": cost, "cost_per_step": per_step,
                    "rho_model_neg": model, "steps_over_sqrt_ell": ro.steps as f64 / (inst.ell as f64).sqrt(),
                    "wall_total": ro.wall,
                });
                eprintln!(
                    "{{\"bits\":{},\"curve\":{},\"arm\":\"rho\",\"correct\":{},\"cost\":{},\"steps\":{},\"wall\":{:.1}}}",
                    bits, args.curve, j["correct"], cost, ro.steps, ro.wall
                );
                results.push(j);
            }
        }
    }
    let rec = serde_json::json!({
        "bits": bits, "q": p, "curve": args.curve, "seed": seed, "mode": args.mode,
        "build_tries": tries, "N": inst.n, "ell": inst.ell, "c": inst.c, "h": inst.h,
        "a": inst.a, "b": inst.b, "r0": inst.r[0], "Kmod": [inst.k.c1, inst.k.c0],
        "build_wall": build_wall, "checks": chk, "results": results,
    });
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(&args.out).unwrap();
    writeln!(f, "{}", rec).unwrap();
}
