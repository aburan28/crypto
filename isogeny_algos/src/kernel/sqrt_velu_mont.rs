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
use super::montgomery::{ladder_p, proj24, xadd, xdbl_p, Proj24, XZ};
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
    /// scaled remainder tree over the x_i (point images: E_J has no symmetry there), built on
    /// the first image
    tree_i: std::cell::OnceCell<poly::EvalTree<F>>,
    /// w_i = x_i + 1/x_i over I (codomain: at alpha = +-1 E_J is palindromic in Z)
    pub ws: Vec<F::E>,
    /// scaled remainder tree over the w_i, built when b' is large enough to use it
    tree_w: std::cell::OnceCell<poly::EvalTree<F>>,
    /// x_i^2 over I (direct point images)
    xs_i2: Vec<F::E>,
    /// Some(true) / Some(false): always / never use the remainder trees (tests, benchmarks);
    /// None: by the size of b' (`TREE_MIN_BP`)
    pub force_tree: Option<bool>,
}

/// From this b' on, the products over I x J are taken by building the product over J as a
/// polynomial and evaluating it with a scaled remainder tree over I (Õ(sqrt l)); below it they
/// are taken pair by pair (b b' multiplications per resultant). Pair by pair was faster for both
/// the codomain and point images at every size measured (l <= 250007, 256 and 511 bits;
/// `results/p13-large-ell.txt`), so the trees are used only on request (`force_tree`).
pub const TREE_MIN_BP: usize = 1 << 20;

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
        Self::new_proj(f, proj24(f, a), kp, ell)
    }

    /// As `new`, from the projective doubling constant (A24 : C24) = (A + 2C : 4C); C24 joins the
    /// batch inversion of the points, so no separate inversion is spent on A.
    pub fn new_proj(f: &F, k: Proj24<F::E>, kp: XZ<F::E>, ell: u64) -> Option<Self> {
        let (b, bp) = sizes(ell);
        let two = xdbl_p(f, k, kp);
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
            let p2b = ladder_p(f, k, kp, &Big::from_u64(2 * b as u64));
            let p4b = xdbl_p(f, k, p2b);
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
        let nk = if k0 < ell - 1 {
            ((ell - 2 - k0) / 2 + 1) as usize
        } else {
            0
        };
        if nk > 0 {
            pts.push(ladder_p(f, k, kp, &Big::from_u64(k0)));
            if nk > 1 {
                pts.push(ladder_p(f, k, kp, &Big::from_u64(k0 + 2)));
            }
            let base = b + bp;
            for t in 2..nk {
                let nx = xadd(f, pts[base + t - 1], two, pts[base + t - 2]);
                pts.push(nx);
            }
        }
        // one batch inversion for every Z and for the X of the I points (1/x_i = Z/X; x = 0 is
        // the 2-torsion point, never a multiple of an odd-order generator)
        let mut zs: Vec<F::E> = pts.iter().map(|p| p.1).collect();
        zs.extend(pts[b..b + bp].iter().map(|p| p.0));
        zs.push(k.1);
        if zs.iter().any(|&z| f.is_zero(z)) {
            return None;
        }
        let zi = f.batch_inv(&zs);
        // A = (4 A24 - 2 C24) / C24
        let four_a24 = f.add(f.add(k.0, k.0), f.add(k.0, k.0));
        let a = f.sub(f.mul(four_a24, zi[zs.len() - 1]), f.from_u64(2));
        let xs: Vec<F::E> = pts.iter().zip(&zi).map(|(p, &iz)| f.mul(p.0, iz)).collect();
        let xs_j = xs[..b].to_vec();
        let xs_i = xs[b..b + bp].to_vec();
        let xs_k = xs[b + bp..].to_vec();
        let np = pts.len();
        let xs_i2: Vec<F::E> = xs_i.iter().map(|&x| f.sq(x)).collect();
        let ws: Vec<F::E> = (0..bp)
            .map(|i| f.add(xs_i[i], f.mul(pts[b + i].1, zi[np + i])))
            .collect();
        Some(SqrtVeluMont {
            a,
            ell,
            b,
            bp,
            xs_j,
            xs_i,
            xs_k,
            tree_i: std::cell::OnceCell::new(),
            xs_i2,
            ws,
            tree_w: std::cell::OnceCell::new(),
            force_tree: None,
        })
    }

    /// E_J(Z) for the evaluation point alpha (degree 2b).
    pub fn e_j(&self, f: &F, alpha: F::E) -> Poly<F> {
        let al2 = f.sq(alpha);
        let layer: Vec<Poly<F>> = self
            .xs_j
            .iter()
            .map(|&x| self.quad(f, alpha, al2, x).to_vec())
            .collect();
        product(f, layer)
    }

    /// [c0, c1, c2] of the factor F0(Z, x) alpha^2 + F1(Z, x) alpha + F2(Z, x) of E_J:
    /// c2 = (alpha - x)^2, c0 = (x alpha - 1)^2, c1 = -2 (x alpha^2 + (x^2 + 1 + 2 A x) alpha + x).
    fn quad(&self, f: &F, alpha: F::E, al2: F::E, x: F::E) -> [F::E; 3] {
        let two_a = f.add(self.a, self.a);
        let c2 = f.sq(f.sub(alpha, x));
        let c0 = f.sq(f.sub(f.mul(x, alpha), f.one()));
        let t = f.add(f.add(f.sq(x), f.one()), f.mul(two_a, x));
        let s = f.add(f.add(f.mul(x, al2), f.mul(t, alpha)), x);
        [c0, f.neg(f.add(s, s)), c2]
    }

    fn use_tree(&self) -> bool {
        self.force_tree.unwrap_or(self.bp >= TREE_MIN_BP)
    }

    /// G(W) with E_J(Z) = Z^b G(Z + 1/Z) for alpha = +-1, where every factor
    /// c2 Z^2 + c1 Z + c0 has c0 = c2 = (alpha - x_j)^2 and so equals Z (c2 W + c1): degree b
    /// instead of 2b, built from linear factors.
    fn g_w(&self, f: &F, alpha: F::E) -> Poly<F> {
        let layer: Vec<Poly<F>> = self
            .xs_j
            .iter()
            .map(|&x| {
                let (c2, c1) = self.c21(f, alpha, x);
                vec![c1, c2]
            })
            .collect();
        product(f, layer)
    }

    /// (c2, c1) of the factor c2 W + c1 for x_j at alpha = +-1: c2 = (alpha - x)^2,
    /// c1 = -2 (x alpha^2 + (x^2 + 1 + 2 A x) alpha + x) with alpha^2 = 1.
    fn c21(&self, f: &F, alpha: F::E, x: F::E) -> (F::E, F::E) {
        let two_a = f.add(self.a, self.a);
        let c2 = f.sq(f.sub(alpha, x));
        let t = f.add(f.add(f.sq(x), f.one()), f.mul(two_a, x));
        let s = f.add(f.add(x, x), f.mul(t, alpha));
        (c2, f.neg(f.add(s, s)))
    }

    /// prod_{i,j} (c2_j w_i + c1_j) = (prod_j c2_j)^b' prod_{i,j} (w_i + c1_j / c2_j): one
    /// multiplication per pair, four independent product chains.
    fn resultant_direct(&self, f: &F, alpha: F::E) -> F::E {
        let (c2, c1): (Vec<F::E>, Vec<F::E>) =
            self.xs_j.iter().map(|&x| self.c21(f, alpha, x)).unzip();
        let c2i = f.batch_inv(&c2);
        let lead = c2.iter().fold(f.one(), |acc, &c| f.mul(acc, c));
        let mut acc = [f.one(); 4];
        for j in 0..self.b {
            let s = f.mul(c1[j], c2i[j]);
            for (i, &w) in self.ws.iter().enumerate() {
                acc[i & 3] = f.mul(acc[i & 3], f.add(w, s));
            }
        }
        let r = f.mul(f.mul(acc[0], acc[1]), f.mul(acc[2], acc[3]));
        f.mul(r, f.pow(lead, self.bp as u128))
    }

    /// prod_i g(x_i) over I.
    fn prod_at_i(&self, f: &F, g: &Poly<F>) -> F::E {
        self.tree_i
            .get_or_init(|| poly::EvalTree::new(f, &self.xs_i))
            .prod(f, g)
    }

    /// h_S(alpha) * Δ * (prod_I x_i)^b for alpha = +-1 (Δ = prod F0(x_i, x_j)): both factors
    /// are common to h_S(1) and h_S(-1) and cancel in the codomain.
    fn h_scaled(&self, f: &F, alpha: F::E) -> F::E {
        let mut r = if self.bp == 0 {
            f.one()
        } else if !self.use_tree() {
            self.resultant_direct(f, alpha)
        } else {
            self.tree_w
                .get_or_init(|| poly::EvalTree::new(f, &self.ws))
                .prod(f, &self.g_w(f, alpha))
        };
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
        let d_new = f.mul(
            f.pow(f.sub(self.a, two), self.ell as u128),
            f.sq(f.sq(f.sq(h1))),
        );
        let a_new = f.mul(
            f.pow(f.add(self.a, two), self.ell as u128),
            f.sq(f.sq(f.sq(hm1))),
        );
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

    /// Image of a projective point (X : Z) not in the kernel, without inversion: the factors of
    /// the pair-by-pair products are homogeneous in (X, Z) of equal degree in numerator and
    /// denominator, so (X num^2 : Z den^2) with alpha = X / Z never formed. Uses the remainder
    /// trees (and one inversion) when b' is large.
    pub fn eval_xz(&self, f: &F, p: XZ<F::E>) -> XZ<F::E> {
        if self.bp > 0 && self.use_tree() {
            let alpha = f.div(p.0, p.1);
            let (num, den) = self.eval_nd(f, alpha);
            return (f.mul(p.0, f.sq(num)), f.mul(p.1, f.sq(den)));
        }
        let (x, z) = p;
        let (x2, z2, xz) = (f.sq(x), f.sq(z), f.mul(x, z));
        let two_a = f.add(self.a, self.a);
        let (mut num, mut den) = (f.one(), f.one());
        for &xj in &self.xs_j {
            // Z^2 times [c0, c1, c2] at alpha = X / Z
            let c2 = f.sq(f.sub(x, f.mul(xj, z)));
            let c0 = f.sq(f.sub(f.mul(xj, x), z));
            let t = f.add(f.add(f.sq(xj), f.one()), f.mul(two_a, xj));
            let s = f.add(f.add(f.mul(xj, x2), f.mul(t, xz)), f.mul(xj, z2));
            let c1 = f.neg(f.add(s, s));
            for (&xi, &xi2) in self.xs_i.iter().zip(&self.xs_i2) {
                let t = f.mul(c1, xi);
                den = f.mul(den, f.add(f.add(f.mul(c2, xi2), t), c0));
                num = f.mul(num, f.add(f.add(f.mul(c0, xi2), t), c2));
            }
        }
        for &xk in &self.xs_k {
            den = f.mul(den, f.sub(x, f.mul(xk, z)));
            num = f.mul(num, f.sub(z, f.mul(x, xk)));
        }
        (f.mul(x, f.sq(num)), f.mul(z, f.sq(den)))
    }

    /// (prod_S (1 - alpha x_s), h_S(alpha)), both times Δ.
    pub fn eval_nd(&self, f: &F, alpha: F::E) -> (F::E, F::E) {
        let (mut num, mut den) = (f.one(), f.one());
        if self.bp > 0 && !self.use_tree() {
            // pair by pair: q_j(x_i) = c2 x_i^2 + c1 x_i + c0 and its reverse c0 x_i^2 + c1 x_i + c2
            // share c1 x_i: five multiplications per pair
            let al2 = f.sq(alpha);
            for &xj in &self.xs_j {
                let [c0, c1, c2] = self.quad(f, alpha, al2, xj);
                for (&x, &x2) in self.xs_i.iter().zip(&self.xs_i2) {
                    let t = f.mul(c1, x);
                    den = f.mul(den, f.add(f.add(f.mul(c2, x2), t), c0));
                    num = f.mul(num, f.add(f.add(f.mul(c0, x2), t), c2));
                }
            }
        } else if self.bp > 0 {
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

/// Product of the polynomials by a balanced tree.
fn product<F: Field>(f: &F, mut layer: Vec<Poly<F>>) -> Poly<F> {
    while layer.len() > 1 {
        layer = layer
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    poly::mul(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    layer.pop().unwrap_or_else(|| vec![f.one()])
}
