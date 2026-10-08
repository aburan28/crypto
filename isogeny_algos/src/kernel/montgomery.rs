//! Montgomery curves y^2 = x^3 + A x^2 + x with x-only arithmetic and the x-only Vélu formula for
//! odd-degree kernels (Costello–Hisil 2017, Renes 2018, Meyer–Reith 2018), over any prime field
//! (`Zp`, `FpM<N>`):
//!   codomain via the twisted-Edwards correspondence (a, d) = (A+2, A-2),
//!   x-map  x' = x * prod_{i=1}^{(l-1)/2} ((x x_i - 1)/(x - x_i))^2.
use crate::bigint::Big;
use crate::curve::Curve;
use crate::field::Field;

pub type XZ<E> = (E, E);

/// Weierstrass model of E_A: x_W = x + A/3, y^2 = x_W^3 + (1 - A^2/3) x_W + A(2A^2 - 9)/27.
pub fn to_weierstrass<F: Field>(f: &F, a: F::E) -> Curve<F::E> {
    let c = |v: u64| f.from_u64(v);
    let a2 = f.mul(a, a);
    let wa = f.sub(c(1), f.div(a2, c(3)));
    let wb = f.div(f.mul(a, f.sub(f.mul(c(2), a2), c(9))), c(27));
    Curve::new(wa, wb)
}

/// j = 256 (A^2 - 3)^3 / (A^2 - 4)
pub fn j_invariant<F: Field>(f: &F, a: F::E) -> F::E {
    let a2 = f.mul(a, a);
    let n = f.sub(a2, f.from_u64(3));
    f.div(
        f.mul(f.from_u64(256), f.mul(n, f.mul(n, n))),
        f.sub(a2, f.from_u64(4)),
    )
}

/// Constant (A+2)/4 used by the doubling formula.
pub fn a24<F: Field>(f: &F, a: F::E) -> F::E {
    f.div(f.add(a, f.from_u64(2)), f.from_u64(4))
}

pub fn infinity<F: Field>(f: &F) -> XZ<F::E> {
    (f.one(), f.zero())
}

pub fn xdbl<F: Field>(f: &F, a24: F::E, p: XZ<F::E>) -> XZ<F::E> {
    let s = f.add(p.0, p.1);
    let d = f.sub(p.0, p.1);
    let t0 = f.sq(s);
    let t1 = f.sq(d);
    let t2 = f.sub(t0, t1);
    (f.mul(t0, t1), f.mul(t2, f.add(t1, f.mul(a24, t2))))
}

/// Differential addition: P + Q given P - Q.
pub fn xadd<F: Field>(f: &F, p: XZ<F::E>, q: XZ<F::E>, diff: XZ<F::E>) -> XZ<F::E> {
    let v0 = f.mul(f.add(p.0, p.1), f.sub(q.0, q.1));
    let v1 = f.mul(f.sub(p.0, p.1), f.add(q.0, q.1));
    let s = f.add(v0, v1);
    let d = f.sub(v0, v1);
    (f.mul(diff.1, f.sq(s)), f.mul(diff.0, f.sq(d)))
}

/// Montgomery ladder on a projective point: [k]P.
pub fn ladder_xz<F: Field>(f: &F, a24: F::E, p: XZ<F::E>, k: &Big) -> XZ<F::E> {
    if k.is_zero() {
        return infinity(f);
    }
    let (mut r0, mut r1) = (p, xdbl(f, a24, p));
    for i in (0..k.bits() - 1).rev() {
        if !k.bit(i) {
            r1 = xadd(f, r0, r1, p);
            r0 = xdbl(f, a24, r0);
        } else {
            r0 = xadd(f, r0, r1, p);
            r1 = xdbl(f, a24, r1);
        }
    }
    r0
}

/// Montgomery ladder: [k]P for affine x(P) = x.
pub fn ladder<F: Field>(f: &F, a24: F::E, x: F::E, k: u128) -> XZ<F::E> {
    ladder_xz(f, a24, (x, f.one()), &Big::from_u128(k))
}

/// Projective image of (X:Z) under the odd-degree isogeny with kernel multiples `mults`:
/// (X prod(X X_i - Z Z_i)^2 : Z prod(X Z_i - Z X_i)^2).
pub fn isog_xz<F: Field>(f: &F, mults: &[XZ<F::E>], p: XZ<F::E>) -> XZ<F::E> {
    let (mut num, mut den) = (f.one(), f.one());
    for &(xi, zi) in mults {
        num = f.mul(num, f.sub(f.mul(p.0, xi), f.mul(p.1, zi)));
        den = f.mul(den, f.sub(f.mul(p.0, zi), f.mul(p.1, xi)));
    }
    (f.mul(p.0, f.sq(num)), f.mul(p.1, f.sq(den)))
}

/// [1]K, [2]K, ..., [d]K as projective x-coordinates.
pub fn multiples<F: Field>(f: &F, a24: F::E, k: XZ<F::E>, d: usize) -> Vec<XZ<F::E>> {
    let mut v = vec![k];
    if d >= 2 {
        v.push(xdbl(f, a24, k));
    }
    for i in 2..d {
        let next = xadd(f, v[i - 1], k, v[i - 2]);
        v.push(next);
    }
    v
}

/// Codomain coefficient A' of E_A / <K> for K of odd order l = 2d+1, from [1]K..[d]K:
/// a' = (A+2)^l prod(X_i+Z_i)^8, d' = (A-2)^l prod(X_i-Z_i)^8, A' = 2(a'+d')/(a'-d').
/// (Pairing and sign verified against Weierstrass Vélu for A != 0 and against the x-map `isog_x`:
/// the other pairing gives a curve with the wrong j once A != 0, the other sign gives the twist.)
pub fn velu_codomain<F: Field>(f: &F, a: F::E, mults: &[XZ<F::E>], ell: u64) -> F::E {
    let (mut pp, mut pm) = (f.one(), f.one());
    for &(x, z) in mults {
        pp = f.mul(pp, f.add(x, z));
        pm = f.mul(pm, f.sub(x, z));
    }
    let p8 = |v: F::E| f.sq(f.sq(f.sq(v)));
    let two = f.from_u64(2);
    let ed_a = f.pow(f.add(a, two), ell as u128);
    let ed_d = f.pow(f.sub(a, two), ell as u128);
    let a_new = f.mul(ed_a, p8(pp));
    let d_new = f.mul(ed_d, p8(pm));
    f.div(f.mul(two, f.add(a_new, d_new)), f.sub(a_new, d_new))
}

/// x-coordinate of the image of the point with abscissa u (u not in the kernel).
pub fn isog_x<F: Field>(f: &F, mults: &[XZ<F::E>], u: F::E) -> F::E {
    let (mut num, mut den) = (f.one(), f.one());
    for &(x, z) in mults {
        num = f.mul(num, f.sub(f.mul(u, x), z));
        den = f.mul(den, f.sub(f.mul(u, z), x));
    }
    let r = f.div(num, den);
    f.mul(u, f.sq(r))
}

// ---------------------------------------------------------------------------------------------
// Projective coefficients: E is stored as (A24 : C24) = (A + 2C : 4C), so no inversion is needed
// per isogeny (Costello–Hisil 2017, Meyer–Reith 2018). `KernelPre` holds (X_i + Z_i, X_i - Z_i)
// for the kernel multiples, shared by the codomain and by every point pushed through.

pub type Proj24<E> = (E, E);

pub fn proj24<F: Field>(f: &F, a: F::E) -> Proj24<F::E> {
    (f.add(a, f.from_u64(2)), f.from_u64(4))
}

/// Affine A = (4 A24 - 2 C24) / C24 (one inversion).
pub fn affine_a<F: Field>(f: &F, k: Proj24<F::E>) -> F::E {
    let four_a24 = f.add(f.add(k.0, k.0), f.add(k.0, k.0));
    f.div(f.sub(four_a24, f.add(k.1, k.1)), k.1)
}

/// Doubling with projective constant: (C24 t0 t1 : t2 (C24 t1 + A24 t2)), 4M + 2S.
pub fn xdbl_p<F: Field>(f: &F, k: Proj24<F::E>, p: XZ<F::E>) -> XZ<F::E> {
    let t0 = f.sq(f.add(p.0, p.1));
    let t1 = f.sq(f.sub(p.0, p.1));
    let t2 = f.sub(t0, t1);
    let c1 = f.mul(k.1, t1);
    (f.mul(c1, t0), f.mul(t2, f.add(c1, f.mul(k.0, t2))))
}

pub fn ladder_p<F: Field>(f: &F, k: Proj24<F::E>, p: XZ<F::E>, n: &Big) -> XZ<F::E> {
    if n.is_zero() {
        return infinity(f);
    }
    let (mut r0, mut r1) = (p, xdbl_p(f, k, p));
    for i in (0..n.bits() - 1).rev() {
        if !n.bit(i) {
            r1 = xadd(f, r0, r1, p);
            r0 = xdbl_p(f, k, r0);
        } else {
            r0 = xadd(f, r0, r1, p);
            r1 = xdbl_p(f, k, r1);
        }
    }
    r0
}

pub fn multiples_p<F: Field>(f: &F, k: Proj24<F::E>, kp: XZ<F::E>, d: usize) -> Vec<XZ<F::E>> {
    let mut v = vec![kp];
    if d >= 2 {
        v.push(xdbl_p(f, k, kp));
    }
    for i in 2..d {
        let next = xadd(f, v[i - 1], kp, v[i - 2]);
        v.push(next);
    }
    v
}

pub struct KernelPre<E> {
    pub sd: Vec<(E, E)>,
}

pub fn kernel_pre<F: Field>(f: &F, mults: &[XZ<F::E>]) -> KernelPre<F::E> {
    KernelPre {
        sd: mults
            .iter()
            .map(|&(x, z)| (f.add(x, z), f.sub(x, z)))
            .collect(),
    }
}

/// Image of (X:Z): with t0 = (X-Z)(X_i+Z_i), t1 = (X+Z)(X_i-Z_i), t0+t1 = 2(X X_i - Z Z_i) and
/// t0-t1 = 2(X Z_i - Z X_i); 4M per kernel multiple instead of 6M.
pub fn isog_xz_pre<F: Field>(f: &F, pre: &KernelPre<F::E>, p: XZ<F::E>) -> XZ<F::E> {
    let (s, d) = (f.add(p.0, p.1), f.sub(p.0, p.1));
    let (mut num, mut den) = (f.one(), f.one());
    for &(si, di) in &pre.sd {
        let t0 = f.mul(d, si);
        let t1 = f.mul(s, di);
        num = f.mul(num, f.add(t0, t1));
        den = f.mul(den, f.sub(t0, t1));
    }
    (f.mul(p.0, f.sq(num)), f.mul(p.1, f.sq(den)))
}

/// Projective codomain: Edwards (a : d) = (A24 : A24 - C24); a' = a^l prod(X_i+Z_i)^8,
/// d' = d^l prod(X_i-Z_i)^8; (A24' : C24') = (a' : a' - d').
pub fn velu_codomain_p<F: Field>(
    f: &F,
    k: Proj24<F::E>,
    pre: &KernelPre<F::E>,
    ell: u64,
) -> Proj24<F::E> {
    let (mut pp, mut pm) = (f.one(), f.one());
    for &(s, d) in &pre.sd {
        pp = f.mul(pp, s);
        pm = f.mul(pm, d);
    }
    let p8 = |v: F::E| f.sq(f.sq(f.sq(v)));
    let a_new = f.mul(f.pow(k.0, ell as u128), p8(pp));
    let d_new = f.mul(f.pow(f.sub(k.0, k.1), ell as u128), p8(pm));
    (a_new, f.sub(a_new, d_new))
}
