//! Montgomery curves y^2 = x^3 + A x^2 + x over F_p with x-only arithmetic and the x-only Vélu
//! formula for odd-degree kernels (Costello–Hisil 2017, Renes 2018, Meyer–Reith 2018):
//!   codomain via the twisted-Edwards correspondence (a, d) = (A+2, A-2),
//!   x-map  x' = x * prod_{i=1}^{(l-1)/2} ((x x_i - 1)/(x - x_i))^2.
use crate::curve::Curve;
use crate::field::{Field, Zp};

pub type XZ = (u64, u64);

/// Weierstrass model of E_A: x_W = x + A/3, y^2 = x_W^3 + (1 - A^2/3) x_W + A(2A^2 - 9)/27.
pub fn to_weierstrass(fp: &Zp, a: u64) -> Curve<u64> {
    let c = |v: u64| fp.from_u64(v);
    let a2 = fp.mul(a, a);
    let wa = fp.sub(c(1), fp.div(a2, c(3)));
    let wb = fp.div(fp.mul(a, fp.sub(fp.mul(c(2), a2), c(9))), c(27));
    Curve::new(wa, wb)
}

/// j = 256 (A^2 - 3)^3 / (A^2 - 4)
pub fn j_invariant(fp: &Zp, a: u64) -> u64 {
    let a2 = fp.mul(a, a);
    let n = fp.sub(a2, 3);
    fp.div(fp.mul(256, fp.mul(n, fp.mul(n, n))), fp.sub(a2, 4))
}

/// Constant (A+2)/4 used by the doubling formula.
pub fn a24(fp: &Zp, a: u64) -> u64 {
    fp.div(fp.add(a, 2), 4)
}

pub fn xdbl(fp: &Zp, a24: u64, p: XZ) -> XZ {
    let t0 = fp.mul(fp.add(p.0, p.1), fp.add(p.0, p.1));
    let t1 = fp.mul(fp.sub(p.0, p.1), fp.sub(p.0, p.1));
    let t2 = fp.sub(t0, t1);
    (fp.mul(t0, t1), fp.mul(t2, fp.add(t1, fp.mul(a24, t2))))
}

/// Differential addition: P + Q given P - Q.
pub fn xadd(fp: &Zp, p: XZ, q: XZ, diff: XZ) -> XZ {
    let v0 = fp.mul(fp.add(p.0, p.1), fp.sub(q.0, q.1));
    let v1 = fp.mul(fp.sub(p.0, p.1), fp.add(q.0, q.1));
    let s = fp.add(v0, v1);
    let d = fp.sub(v0, v1);
    (fp.mul(diff.1, fp.mul(s, s)), fp.mul(diff.0, fp.mul(d, d)))
}

/// Montgomery ladder: [k]P for affine x(P) = x.
pub fn ladder(fp: &Zp, a24: u64, x: u64, k: u128) -> XZ {
    let p = (x, 1);
    let (mut r0, mut r1) = ((1u64, 0u64), p);
    for i in (0..128).rev() {
        let bit = (k >> i) & 1;
        if bit == 0 {
            r1 = xadd(fp, r0, r1, p);
            r0 = xdbl(fp, a24, r0);
        } else {
            r0 = xadd(fp, r0, r1, p);
            r1 = xdbl(fp, a24, r1);
        }
    }
    r0
}

/// [1]K, [2]K, ..., [d]K as projective x-coordinates.
pub fn multiples(fp: &Zp, a24: u64, k: XZ, d: usize) -> Vec<XZ> {
    let mut v = vec![k];
    if d >= 2 {
        v.push(xdbl(fp, a24, k));
    }
    for i in 2..d {
        let next = xadd(fp, v[i - 1], k, v[i - 2]);
        v.push(next);
    }
    v
}

/// Codomain coefficient A' of E_A / <K> for K of odd order l = 2d+1, from [1]K..[d]K.
pub fn velu_codomain(fp: &Zp, a: u64, mults: &[XZ], ell: u64) -> u64 {
    let (mut pp, mut pm) = (1u64, 1u64);
    for &(x, z) in mults {
        pp = fp.mul(pp, fp.add(x, z));
        pm = fp.mul(pm, fp.sub(x, z));
    }
    let p8 = |v: u64| {
        let v2 = fp.mul(v, v);
        let v4 = fp.mul(v2, v2);
        fp.mul(v4, v4)
    };
    let ed_a = fp.pow(fp.add(a, 2), ell as u128);
    let ed_d = fp.pow(fp.sub(a, 2), ell as u128);
    // a' = (A+2)^l prod(X_i+Z_i)^8, d' = (A-2)^l prod(X_i-Z_i)^8, A' = 2(a'+d')/(a'-d').
    // (Verified against Weierstrass Velu for A != 0, and against the x-map `isog_x` for the
    // twist/sign: the other pairing gives a curve with the wrong j once A != 0.)
    let a_new = fp.mul(ed_a, p8(pp));
    let d_new = fp.mul(ed_d, p8(pm));
    fp.div(fp.mul(2, fp.add(a_new, d_new)), fp.sub(a_new, d_new))
}

/// x-coordinate of the image of the point with abscissa u (u not in the kernel).
pub fn isog_x(fp: &Zp, mults: &[XZ], u: u64) -> u64 {
    let (mut num, mut den) = (1u64, 1u64);
    for &(x, z) in mults {
        num = fp.mul(num, fp.sub(fp.mul(u, x), z));
        den = fp.mul(den, fp.sub(fp.mul(u, z), x));
    }
    let r = fp.div(num, den);
    fp.mul(u, fp.mul(r, r))
}
