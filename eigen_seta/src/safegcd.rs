//! Constant-time modular inversion by Bernstein–Yang "safegcd" divsteps,
//! generalising the structure of libsecp256k1's `modinv64` (signed 62-bit
//! limbs, 62 divsteps per transition matrix) to arbitrary limb counts.
//!
//! The number of divsteps is the proven bound floor((49 B + 57) / 17) for
//! B-bit inputs (Bernstein–Yang, Theorem 11.2), rounded up to a multiple of
//! 62; every loop runs a public number of times and every data-dependent
//! decision is a mask.  inv(0) = 0.

const M62: u64 = u64::MAX >> 2;
pub const MAXL: usize = 48; // up to 48*62 = 2976 bits

#[derive(Clone, Copy)]
struct S62 {
    v: [i64; MAXL],
}

/// Split nonnegative little-endian u64 limbs into L 62-bit limbs.
fn to_s62(a: &[u64], l: usize) -> S62 {
    let mut v = [0i64; MAXL];
    for i in 0..l {
        let bit = 62 * i;
        let w = bit / 64;
        let o = bit % 64;
        let lo = if w < a.len() { a[w] >> o } else { 0 };
        let hi = if o > 2 && w + 1 < a.len() {
            a[w + 1] << (64 - o)
        } else {
            0
        };
        v[i] = ((lo | hi) & M62) as i64;
    }
    S62 { v }
}

/// Pack nonnegative normalised 62-bit limbs back into n u64 limbs.
fn from_s62(s: &S62, l: usize, n: usize) -> Vec<u64> {
    let mut out = vec![0u64; n];
    for i in 0..l {
        let x = s.v[i] as u64 & M62;
        let bit = 62 * i;
        let w = bit / 64;
        let o = bit % 64;
        if w < n {
            out[w] |= x << o;
        }
        if o > 2 && w + 1 < n {
            out[w + 1] |= x >> (64 - o);
        }
    }
    out
}

/// 62 divsteps on the low bits of f and g; returns delta and the transition
/// matrix (scaled by 2^62).
#[inline(always)]
fn divsteps_62(mut delta: i64, f0: u64, g0: u64) -> (i64, i64, i64, i64, i64) {
    let (mut u, mut v, mut q, mut r) = (1i64, 0i64, 0i64, 1i64);
    let (mut f, mut g) = (f0, g0);
    for _ in 0..62 {
        // c1 = all ones iff delta > 0 and g odd
        let dpos = ((delta.wrapping_neg()) >> 63) as u64;
        let c1 = core::hint::black_box(dpos & (g & 1).wrapping_neg());
        // conditional swap with negation: (f, g, u, v, q, r, delta) -> (g, -f, q, r, -u, -v, -delta)
        let x = (f ^ g) & c1;
        f ^= x;
        g ^= x;
        g = (g ^ c1).wrapping_sub(c1);
        let c1s = c1 as i64;
        let xu = (u ^ q) & c1s;
        u ^= xu;
        q ^= xu;
        q = (q ^ c1s).wrapping_sub(c1s);
        let xv = (v ^ r) & c1s;
        v ^= xv;
        r ^= xv;
        r = (r ^ c1s).wrapping_sub(c1s);
        delta = (delta ^ c1s).wrapping_sub(c1s);
        // if g odd: g += f, (q, r) += (u, v)
        let c2 = core::hint::black_box((g & 1).wrapping_neg());
        let c2s = c2 as i64;
        g = g.wrapping_add(f & c2);
        q = q.wrapping_add(u & c2s);
        r = r.wrapping_add(v & c2s);
        // delta += 1, g /= 2, (u, v) *= 2
        delta = delta.wrapping_add(1);
        g = ((g as i64) >> 1) as u64;
        u = u.wrapping_shl(1);
        v = v.wrapping_shl(1);
    }
    (delta, u, v, q, r)
}

fn update_fg(f: &mut S62, g: &mut S62, l: usize, t: (i64, i64, i64, i64)) {
    let (u, v, q, r) = (t.0 as i128, t.1 as i128, t.2 as i128, t.3 as i128);
    let mut cf = u * f.v[0] as i128 + v * g.v[0] as i128;
    let mut cg = q * f.v[0] as i128 + r * g.v[0] as i128;
    cf >>= 62;
    cg >>= 62;
    for i in 1..l {
        cf += u * f.v[i] as i128 + v * g.v[i] as i128;
        cg += q * f.v[i] as i128 + r * g.v[i] as i128;
        f.v[i - 1] = (cf as u64 & M62) as i64;
        g.v[i - 1] = (cg as u64 & M62) as i64;
        cf >>= 62;
        cg >>= 62;
    }
    f.v[l - 1] = cf as i64;
    g.v[l - 1] = cg as i64;
}

fn update_de(d: &mut S62, e: &mut S62, l: usize, t: (i64, i64, i64, i64), m: &S62, minv62: u64) {
    let (u, v, q, r) = t;
    let sd = core::hint::black_box(d.v[l - 1] >> 63);
    let se = core::hint::black_box(e.v[l - 1] >> 63);
    let mut md = (u & sd).wrapping_add(v & se);
    let mut me = (q & sd).wrapping_add(r & se);
    let (ui, vi, qi, ri) = (u as i128, v as i128, q as i128, r as i128);
    let mut cd = ui * d.v[0] as i128 + vi * e.v[0] as i128;
    let mut ce = qi * d.v[0] as i128 + ri * e.v[0] as i128;
    md = md.wrapping_sub((minv62.wrapping_mul(cd as u64).wrapping_add(md as u64) & M62) as i64);
    me = me.wrapping_sub((minv62.wrapping_mul(ce as u64).wrapping_add(me as u64) & M62) as i64);
    cd += m.v[0] as i128 * md as i128;
    ce += m.v[0] as i128 * me as i128;
    cd >>= 62;
    ce >>= 62;
    for i in 1..l {
        cd += ui * d.v[i] as i128 + vi * e.v[i] as i128 + m.v[i] as i128 * md as i128;
        ce += qi * d.v[i] as i128 + ri * e.v[i] as i128 + m.v[i] as i128 * me as i128;
        d.v[i - 1] = (cd as u64 & M62) as i64;
        e.v[i - 1] = (ce as u64 & M62) as i64;
        cd >>= 62;
        ce >>= 62;
    }
    d.v[l - 1] = cd as i64;
    e.v[l - 1] = ce as i64;
}

fn normalize(r: &mut S62, l: usize, sign: i64, m: &S62) {
    // add the modulus if negative, then negate if requested: (-2m, m) -> (-m, m)
    // The masks go through black_box: without it LLVM turned the masked adds
    // into conditionals (skipped when the mask is zero), which the
    // instruction-count probe (examples/ct_probe.rs) detected.
    let cond_add = core::hint::black_box(r.v[l - 1] >> 63);
    for i in 0..l {
        r.v[i] = r.v[i].wrapping_add(m.v[i] & cond_add);
    }
    let cond_neg = core::hint::black_box(sign >> 63);
    for i in 0..l {
        r.v[i] = (r.v[i] ^ cond_neg).wrapping_sub(cond_neg);
    }
    for i in 0..l - 1 {
        r.v[i + 1] = r.v[i + 1].wrapping_add(r.v[i] >> 62);
        r.v[i] &= M62 as i64;
    }
    // add the modulus again if still negative: -> [0, m)
    let cond_add = core::hint::black_box(r.v[l - 1] >> 63);
    for i in 0..l {
        r.v[i] = r.v[i].wrapping_add(m.v[i] & cond_add);
    }
    for i in 0..l - 1 {
        r.v[i + 1] = r.v[i + 1].wrapping_add(r.v[i] >> 62);
        r.v[i] &= M62 as i64;
    }
}

/// a^{-1} mod p for canonical (non-Montgomery) limbs; p odd, `pinv64` =
/// -p^{-1} mod 2^64, `bits` = bit length of p (public).
pub fn modinv(a: &[u64], p: &[u64], pinv64: u64, bits: usize) -> Vec<u64> {
    let n = p.len();
    let l = (64 * n) / 62 + 2;
    assert!(l <= MAXL);
    let m = to_s62(p, l);
    let minv62 = pinv64.wrapping_neg() & M62; // p^{-1} mod 2^62
    let mut f = m;
    let mut g = to_s62(a, l);
    let mut d = S62 { v: [0i64; MAXL] };
    let mut e = S62 { v: [0i64; MAXL] };
    e.v[0] = 1;
    let mut delta = 1i64;
    let total = (49 * bits + 57) / 17;
    let batches = total.div_ceil(62);
    for _ in 0..batches {
        let f0 = (f.v[0] as u64) | ((f.v[1] as u64) << 62);
        let g0 = (g.v[0] as u64) | ((g.v[1] as u64) << 62);
        let (nd, u, v, q, r) = divsteps_62(delta, f0, g0);
        delta = nd;
        update_de(&mut d, &mut e, l, (u, v, q, r), &m, minv62);
        update_fg(&mut f, &mut g, l, (u, v, q, r));
    }
    // f = +-1 now (for invertible a); d * f is the inverse
    normalize(&mut d, l, f.v[l - 1], &m);
    from_s62(&d, l, n)
}
