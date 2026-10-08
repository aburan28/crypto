//! √élu on Montgomery curves (Bernstein–De Feo–Leroux–Smith 2020, "Faster computation of
//! isogenies of large prime degree"): the codomain and the x-map of an odd prime-degree isogeny
//! from Õ(√l) operations.
//!
//! With S = {1, 3, ..., l-2} (x([s]P) over S is every kernel abscissa once) and
//! h_S(X) = prod_{s in S} (X - x([s]P)):
//!   codomain  d = ((A-2)/(A+2))^l (h_S(1)/h_S(-1))^8,   A' = 2(1+d)/(1-d);
//!   x-map     x' = a (prod_S (1 - a x_s) / h_S(a))^2.
//! S is split as (I ± J) ∪ K with I = {2b(2i+1) : i < b'}, J = {2j+1 : j < b}, so
//! prod_{i,j} (a - x(P_i+P_j))(a - x(P_i-P_j)) = prod_i E_J(x_i) / Δ with
//! E_J(Z) = prod_j (F0(Z,x_j) a^2 + F1(Z,x_j) a + F2(Z,x_j)) for the Montgomery biquadratics
//!   F0 = (Z - x)^2,  F1 = -2((Z x + 1)(Z + x) + 2 A Z x),  F2 = (Z x - 1)^2,
//! and Δ = prod F0(x_i, x_j) cancels in every ratio above. prod_i E_J(x_i) is a multipoint
//! evaluation of a degree-2b polynomial at the b' roots of h_I (remainder tree), and the
//! "1 - a x" products use the reversed polynomial of E_J at the same points.
//! Points are affine after one batch inversion; the generic polynomial arithmetic (Karatsuba
//! products, Newton remainders) does the rest. Cross-checked against x-only Vélu.
use super::montgomery::{a24, xadd, xdbl, ladder_xz, XZ};
use crate::bigint::Big;
use crate::field::Field;
use crate::poly::{self, Poly};

pub struct SqrtVeluMont<F: Field> {
    pub a: F::E,
    pub ell: u64,
    pub b: usize,
    pub bp: usize,
    pub xs_j: Vec<F::E>,
    pub xs_i: Vec<F::E>,
    pub xs_k: Vec<F::E>,
    tree_i: Vec<Vec<Poly<F>>>,
}

/// (b, b') for degree l: b = floor(sqrt(l-1)/2), b' = floor((l-1)/(4b)); (0, 0) when l is small.
pub fn sizes(ell: u64) -> (usize, usize) {
    let b = ((ell - 1) as f64).sqrt() as usize / 2;
    if b == 0 {
        return (0, 0);
    }
    let bp = (ell as usize - 1) / (4 * b);
    if bp == 0 {
        (0, 0)
    } else {
        (b, bp)
    }
}

impl<F: Field> SqrtVeluMont<F> {
    /// Kernel point set from the generator (projective x on E_A) of order l. None if a needed
    /// multiple is the point at infinity (wrong order).
    pub fn new(f: &F, a: F::E, kp: XZ<F::E>, ell: u64) -> Option<Self> {
        let k24 = a24(f, a);
        let (b, bp) = sizes(ell);
        let two = xdbl(f, k24, kp);
        let mut pts: Vec<XZ<F::E>> = vec![];
        // J: odd multiples 1, 3, ..., 2b-1
        if b > 0 {
            pts.push(kp);
            if b > 1 {
                pts.push(xadd(f, two, kp, kp));
            }
            for j in 2..b {
                let nx = xadd(f, pts[j - 1], two, pts[j - 2]);
                pts.push(nx);
            }
            // I: 2b(2i+1), step 4b
            let p2b = ladder_xz(f, k24, kp, &Big::from_u64(2 * b as u64));
            let p4b = xdbl(f, k24, p2b);
            pts.push(p2b);
            if bp > 1 {
                pts.push(xadd(f, p4b, p2b, p2b));
            }
            for i in 2..bp {
                let nx = xadd(f, pts[b + i - 1], p4b, pts[b + i - 2]);
                pts.push(nx);
            }
        }
        // K: odd k in [4 b b' + 1, l - 2]
        let k0 = (4 * b * bp + 1) as u64;
        let nk = if k0 + 1 <= ell - 1 { ((ell - 2 - k0) / 2 + 1) as usize } else { 0 };
        if nk > 0 {
            pts.push(ladder_xz(f, k24, kp, &Big::from_u64(k0)));
            if nk > 1 {
                pts.push(ladder_xz(f, k24, kp, &Big::from_u64(k0 + 2)));
            }
            let base = b + bp;
            for t in 2..nk {
                let nx = xadd(f, pts[base + t - 1], two, pts[base + t - 2]);
                pts.push(nx);
            }
        }
        let zs: Vec<F::E> = pts.iter().map(|p| p.1).collect();
        if zs.iter().any(|&z| f.is_zero(z)) {
            return None;
        }
        let zi = f.batch_inv(&zs);
        let xs: Vec<F::E> = pts.iter().zip(zi).map(|(p, iz)| f.mul(p.0, iz)).collect();
        let xs_j = xs[..b].to_vec();
        let xs_i = xs[b..b + bp].to_vec();
        let xs_k = xs[b + bp..].to_vec();
        let tree_i = if bp > 0 { poly::subproduct_tree(f, &xs_i) } else { vec![] };
        Some(SqrtVeluMont { a, ell, b, bp, xs_j, xs_i, xs_k, tree_i })
    }

    /// E_J(Z) for the evaluation point alpha (degree 2b).
    fn e_j(&self, f: &F, alpha: F::E) -> Poly<F> {
        let two = f.from_u64(2);
        let al2 = f.sq(alpha);
        let two_a = f.add(self.a, self.a);
        let mut layer: Vec<Poly<F>> = self
            .xs_j
            .iter()
            .map(|&x| {
                let c2 = f.sq(f.sub(alpha, x));
                let c0 = f.sq(f.sub(f.mul(x, alpha), f.one()));
                // Z^1: -2 x a^2 - 2 (x^2 + 1 + 2 A x) a - 2 x
                let t = f.add(f.add(f.sq(x), f.one()), f.mul(two_a, x));
                let c1 = f.neg(f.mul(two, f.add(f.add(f.mul(x, al2), f.mul(t, alpha)), x)));
                vec![c0, c1, c2]
            })
            .collect();
        while layer.len() > 1 {
            layer = layer
                .chunks(2)
                .map(|c| if c.len() == 2 { poly::mul(f, &c[0], &c[1]) } else { c[0].clone() })
                .collect();
        }
        layer.pop().unwrap_or_else(|| vec![f.one()])
    }

    /// prod_i g(x_i) over I.
    fn prod_at_i(&self, f: &F, g: &Poly<F>) -> F::E {
        let vals = if self.bp <= 8 {
            self.xs_i.iter().map(|&x| poly::eval(f, g, x)).collect()
        } else {
            poly::multipoint_eval_tree(f, g, &self.tree_i)
        };
        vals.into_iter().fold(f.one(), |acc, v| f.mul(acc, v))
    }

    /// h_S(alpha) * Δ (Δ = prod F0(x_i, x_j), common to every call).
    fn h_scaled(&self, f: &F, alpha: F::E) -> F::E {
        let mut r = if self.bp > 0 { self.prod_at_i(f, &self.e_j(f, alpha)) } else { f.one() };
        for &x in &self.xs_k {
            r = f.mul(r, f.sub(alpha, x));
        }
        r
    }

    /// Codomain as projective doubling constants (A24' : C24') = (a' : a' - d') with the Edwards
    /// a' = (A+2)^l h_S(-1)^8, d' = (A-2)^l h_S(1)^8 (both times Δ^8): no inversion.
    pub fn codomain_proj(&self, f: &F) -> (F::E, F::E) {
        let one = f.one();
        let h1 = self.h_scaled(f, one);
        let hm1 = self.h_scaled(f, f.neg(one));
        let two = f.from_u64(2);
        let d_new = f.mul(f.pow(f.sub(self.a, two), self.ell as u128), f.sq(f.sq(f.sq(h1))));
        let a_new = f.mul(f.pow(f.add(self.a, two), self.ell as u128), f.sq(f.sq(f.sq(hm1))));
        (a_new, f.sub(a_new, d_new))
    }

    /// Codomain coefficient A' (affine: one inversion).
    pub fn codomain(&self, f: &F) -> F::E {
        super::montgomery::affine_a(f, self.codomain_proj(f))
    }

    /// x-coordinate of the image of the point with abscissa alpha (not in the kernel).
    pub fn eval(&self, f: &F, alpha: F::E) -> F::E {
        let (num, den) = self.eval_nd(f, alpha);
        f.mul(alpha, f.sq(f.div(num, den)))
    }

    /// (prod_S (1 - alpha x_s), h_S(alpha)), both times Δ.
    pub fn eval_nd(&self, f: &F, alpha: F::E) -> (F::E, F::E) {
        let (mut num, mut den) = (f.one(), f.one());
        if self.bp > 0 {
            let e = self.e_j(f, alpha);
            let mut r = e.clone();
            r.reverse();
            den = self.prod_at_i(f, &e);
            num = self.prod_at_i(f, &r);
        }
        for &x in &self.xs_k {
            den = f.mul(den, f.sub(alpha, x));
            num = f.mul(num, f.sub(f.one(), f.mul(alpha, x)));
        }
        (num, den)
    }
}
