//! Twisted Hessian curves a X^3 + Y^3 + Z^3 = d X Y Z (Bernstein–Kohel–Lange 2015), identity
//! (0 : -1 : 1), -(X : Y : Z) = (X : Z : Y), and their l-isogenies for odd l prime to 3 as the
//! product of translates phi(P) = (prod_Q X_{P+Q} : prod_Q Y_{P+Q} : prod_Q Z_{P+Q}) over the
//! kernel (the shape of Moody–Shumow's Edwards/Huff isogenies; Dang–Moody 2019 for Hessians).
//!
//! Explicit form: with the rotated addition law
//!   P + Q = (X_P^2 Y_Q Z_Q - Y_P Z_P X_Q^2 : Z_P^2 X_Q Y_Q - X_P Y_P Z_Q^2 : Y_P^2 X_Q Z_Q - X_P Z_P Y_Q^2)
//! and Q, -Q paired, phi(P) is a product over half the kernel of quadratic forms in the
//! monomials of P with constants from Q: no point additions.
//!
//! Codomain (derived here, checked by tests): a' = a^l and, from the tangent at the identity
//! (Y/Z + 1 ~ -(d/3) X/Z on the domain, pulled back through phi with the invariant derivation
//! (d y - 3 a x^2)/3 at each kernel point),
//!   d' = (l d - 6 a sum_Q s_Q) / prod_Q s_Q,   s_Q = X_Q^2 / (Y_Q Z_Q),  Q over half the kernel.
use crate::field::Field;

/// a X^3 + Y^3 + Z^3 = d X Y Z
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hessian<E> {
    pub a: E,
    pub d: E,
}

pub type HPt<E> = [E; 3];

impl<E: Copy + PartialEq> Hessian<E> {
    pub fn identity<F: Field<E = E>>(&self, f: &F) -> HPt<E> {
        [f.zero(), f.neg(f.one()), f.one()]
    }
    pub fn on_curve<F: Field<E = E>>(&self, f: &F, p: &HPt<E>) -> bool {
        let [x, y, z] = *p;
        let lhs = f.add(
            f.add(f.mul(self.a, f.mul(x, f.sq(x))), f.mul(y, f.sq(y))),
            f.mul(z, f.sq(z)),
        );
        lhs == f.mul(self.d, f.mul(x, f.mul(y, z)))
            && !(f.is_zero(x) && f.is_zero(y) && f.is_zero(z))
    }
    pub fn neg(&self, p: &HPt<E>) -> HPt<E> {
        [p[0], p[2], p[1]]
    }
    /// Sum by the rotated addition law, falling back to the twisted one (which also doubles).
    pub fn add<F: Field<E = E>>(&self, f: &F, p: &HPt<E>, q: &HPt<E>) -> HPt<E> {
        let [x1, y1, z1] = *p;
        let [x2, y2, z2] = *q;
        let r = [
            f.sub(
                f.mul(f.sq(x1), f.mul(y2, z2)),
                f.mul(f.mul(y1, z1), f.sq(x2)),
            ),
            f.sub(
                f.mul(f.sq(z1), f.mul(x2, y2)),
                f.mul(f.mul(x1, y1), f.sq(z2)),
            ),
            f.sub(
                f.mul(f.sq(y1), f.mul(x2, z2)),
                f.mul(f.mul(x1, z1), f.sq(y2)),
            ),
        ];
        if !r.iter().all(|&v| f.is_zero(v)) {
            return r;
        }
        [
            f.sub(
                f.mul(f.sq(z2), f.mul(x1, z1)),
                f.mul(f.sq(y1), f.mul(x2, y2)),
            ),
            f.sub(
                f.mul(f.sq(y2), f.mul(y1, z1)),
                f.mul(self.a, f.mul(f.sq(x1), f.mul(x2, z2))),
            ),
            f.sub(
                f.mul(self.a, f.mul(f.sq(x2), f.mul(x1, y1))),
                f.mul(f.sq(z1), f.mul(y2, z2)),
            ),
        ]
    }
    pub fn mul<F: Field<E = E>>(&self, f: &F, p: &HPt<E>, mut k: u64) -> HPt<E> {
        let mut r = self.identity(f);
        let mut b = *p;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(f, &r, &b);
            }
            b = self.add(f, &b, &b);
            k >>= 1;
        }
        r
    }
    pub fn eq<F: Field<E = E>>(&self, f: &F, p: &HPt<E>, q: &HPt<E>) -> bool {
        (0..3).all(|i| (0..3).all(|j| f.mul(p[i], q[j]) == f.mul(p[j], q[i])))
    }
    /// j = d^3 (d^3 + 216 a)^3 / (a (d^3 - 27 a)^3) (Hesse pencil, twisted by X -> a^(1/3) X).
    pub fn j<F: Field<E = E>>(&self, f: &F) -> E {
        let d3 = f.mul(self.d, f.sq(self.d));
        let n = f.add(d3, f.mul(f.from_u64(216), self.a));
        let m = f.sub(d3, f.mul(f.from_u64(27), self.a));
        f.div(
            f.mul(d3, f.mul(n, f.sq(n))),
            f.mul(self.a, f.mul(m, f.sq(m))),
        )
    }
}

/// The l-isogeny with kernel <k>, l odd and prime to 3.
#[derive(Clone, Debug)]
pub struct HessianIso<E> {
    pub dom: Hessian<E>,
    pub cod: Hessian<E>,
    /// per kernel pair {Q, -Q}: (Y_Q Z_Q, X_Q^2, X_Q Y_Q, Z_Q^2, X_Q Z_Q, Y_Q^2)
    pub consts: Vec<[E; 6]>,
}

pub fn hessian_isogeny<F: Field>(
    f: &F,
    e: &Hessian<F::E>,
    k: &HPt<F::E>,
    ell: u64,
) -> Option<HessianIso<F::E>> {
    if ell.is_multiple_of(2) || ell.is_multiple_of(3) {
        return None;
    }
    let mut consts = Vec::with_capacity((ell / 2) as usize);
    let mut q = *k;
    let mut den_prod = f.one();
    // s_Q = X^2 / (Y Z): sum and product projectively (common denominator prod Y Z)
    let mut yz_prod = f.one();
    let mut sum_num = f.zero(); // sum_Q X_Q^2 prod_{Q' != Q} Y_Q' Z_Q'
    for i in 0..ell / 2 {
        if i > 0 {
            q = e.add(f, &q, k);
        }
        let [x, y, z] = q;
        let (yz, x2) = (f.mul(y, z), f.sq(x));
        consts.push([yz, x2, f.mul(x, y), f.sq(z), f.mul(x, z), f.sq(y)]);
        sum_num = f.add(f.mul(sum_num, yz), f.mul(x2, yz_prod));
        yz_prod = f.mul(yz_prod, yz);
        den_prod = f.mul(den_prod, x2);
    }
    // d' = (l d - 6 a sum s) / prod s = (l d yz_prod - 6 a sum_num) / den_prod
    let num_sum = f.sub(
        f.mul(f.mul(f.from_u64(ell), e.d), yz_prod),
        f.mul(f.mul(f.from_u64(6), e.a), sum_num),
    );
    if f.is_zero(den_prod) {
        return None;
    }
    let d2 = f.div(num_sum, den_prod);
    let cod = Hessian {
        a: f.pow(e.a, ell as u128),
        d: d2,
    };
    Some(HessianIso {
        dom: *e,
        cod,
        consts,
    })
}

impl<E: Copy + PartialEq> HessianIso<E> {
    /// phi(P) by the explicit product (two factors per kernel pair and coordinate).
    pub fn eval<F: Field<E = E>>(&self, f: &F, p: &HPt<E>) -> HPt<E> {
        let [x, y, z] = *p;
        let (x2, y2, z2) = (f.sq(x), f.sq(y), f.sq(z));
        let (xy, xz, yz) = (f.mul(x, y), f.mul(x, z), f.mul(y, z));
        let (mut px, mut py, mut pz) = (x, y, z);
        for c in &self.consts {
            let [qyz, qx2, qxy, qz2, qxz, qy2] = *c;
            let tx = f.sub(f.mul(x2, qyz), f.mul(yz, qx2));
            px = f.mul(px, f.sq(tx));
            let ty1 = f.sub(f.mul(z2, qxy), f.mul(xy, qz2));
            let ty2 = f.sub(f.mul(z2, qxz), f.mul(xy, qy2));
            py = f.mul(py, f.mul(ty1, ty2));
            let tz1 = f.sub(f.mul(y2, qxz), f.mul(xz, qy2));
            let tz2 = f.sub(f.mul(y2, qxy), f.mul(xz, qz2));
            pz = f.mul(pz, f.mul(tz1, tz2));
        }
        [px, py, pz]
    }
}
