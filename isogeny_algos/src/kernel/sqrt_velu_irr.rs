//! √élu with an irrational kernel generator (Bernstein–De Feo–Leroux–Smith 2020, "Faster
//! computation of isogenies of large prime degree", §4.14), on Montgomery curves
//! E_A: y^2 = x^3 + A x^2 + x with A in F_q, plus the Appendix A.5 time-saving techniques that
//! the Montgomery √élu (`sqrt_velu_mont`) did not yet have.
//!
//! # §4.14
//! The kernel G = <P> (prime order l) is Galois-stable but P is not rational: x(P) lies in
//! F_{q^e} for a minimal e >= 1, and the q-power Frobenius acts as pi(P) = [lambda] P with lambda
//! of order e in (Z/lZ)^* / {±1} (the x-line forgets the sign). With L = <lambda>, the index set
//! S = (Z/lZ)^* / {±1} (one class per pair ±s) splits into L-orbits, and
//!     h_S(X) = prod_{s in S} (X - x([s]P)) = prod_{orbits O} g_O(X),
//!     g_O(X) = prod_{a < e} (X - x([s_O]P)^{q^a}) in F_q[X],
//! because x([lambda s]P) = x([s]P)^q. For alpha in F_q, g_O(alpha) = N(alpha - x([s_O]P)) with N
//! the norm F_{q^e} -> F_q, so only one kernel abscissa per orbit is ever computed.
//!
//! Algorithm 2 needs index systems of the form (I, J) = (L I', L J'): then h_I, E_J and D_J are
//! products of F_q-polynomials, one per orbit,
//!     h_I(Z) = prod_{O in I} g_O(Z),   E_J(alpha, Z) = prod_{O in J} prod_{a<e} Q(alpha; Z, x_O^{q^a}),
//! (Q = F0 alpha^2 + F1 alpha + F2 with the Montgomery biquadratics of Example 4.4, all of whose
//! coefficients are in F_q, so the orbit product is Frobenius-invariant), and the resultant
//! Res_Z(h_I, E_J) = prod_{O in I} N(r_O(x_O)) with r_O = E_J mod g_O (a remainder tree over the
//! monic degree-e leaves g_O, `strassen::RemTree`). All polynomial work is in F_q; the extension
//! is used only for the orbit products of degree e and the norms.
//!
//! The paper warns that such L-stable index systems with #I, #J in O~(√l) typically do not
//! exist (e is usually of size l). There is no closed form for them, so [`index_system`] builds
//! one greedily: J = the first n_J orbits, then every orbit whose 2 e #J sum/difference classes
//! are new joins I; K is the rest. The sizes it reaches are reported in [`SqrtVeluIrr::sizes`];
//! for small e it covers most of S (see the tests), for e near l it degenerates to K = S, which
//! is the orbit-norm product [`SqrtVeluIrr::h_orbits`] with (l-1)/(2e) norms.
//!
//! Outputs: the codomain coefficient A' (via the Edwards form, as in §5.2) and the x-map
//! x' = x (prod_S (1 - x x_s) / h_S(x))^2 (§5.1), both in F_q, for points of E(F_q) or of its twist.
//!
//! # A.5
//! See [`fused_codomain_eval`] (merging Q(alpha), Q(1/alpha), Q(1), Q(-1)) and
//! `strassen::RemTree` (precomputed reciprocals of tree nodes, Karatsuba low products). The
//! module tests compare both with the unmerged code and with Vélu.
use super::montgomery::{a24, ladder_xz, xadd, XZ};
use super::sqrt_velu_mont::SqrtVeluMont;
use crate::bigint::Big;
use crate::ext::Fp4;
use crate::field::{Field, Zp, Zp2};
use crate::fp2::Fp2;
use crate::fpr::{FpR, ER};
use crate::poly::{self, Poly};
use crate::strassen::RemTree;
use std::collections::HashSet;

/// A finite extension F_{q^e} of a base field F_q with its q-power Frobenius.
pub trait ExtField: Field {
    type B: Field;
    /// The base field F_q.
    fn base_field(&self) -> Self::B;
    /// F_q -> F_{q^e}.
    fn lift(&self, a: <Self::B as Field>::E) -> Self::E;
    /// a -> a^q.
    fn frob(&self, a: Self::E) -> Self::E;
    /// The element of F_q, if a lies in it.
    fn down(&self, a: Self::E) -> Option<<Self::B as Field>::E>;
}

impl ExtField for Zp2 {
    type B = Zp;
    fn base_field(&self) -> Zp {
        self.base
    }
    fn lift(&self, a: u64) -> (u64, u64) {
        (a, 0)
    }
    /// i^p = -i for i^2 = nr a non-residue.
    fn frob(&self, a: (u64, u64)) -> (u64, u64) {
        (a.0, self.base.neg(a.1))
    }
    fn down(&self, a: (u64, u64)) -> Option<u64> {
        (a.1 == 0).then_some(a.0)
    }
}

impl<F: Field> ExtField for Fp2<F> {
    type B = F;
    fn base_field(&self) -> F {
        self.base.clone()
    }
    fn lift(&self, a: F::E) -> (F::E, F::E) {
        self.embed(a)
    }
    fn frob(&self, a: (F::E, F::E)) -> (F::E, F::E) {
        self.conj(a)
    }
    fn down(&self, a: (F::E, F::E)) -> Option<F::E> {
        self.base.is_zero(a.1).then_some(a.0)
    }
}

impl ExtField for Fp4 {
    type B = Zp;
    fn base_field(&self) -> Zp {
        self.f2.base
    }
    fn lift(&self, a: u64) -> ((u64, u64), (u64, u64)) {
        ((a, 0), (0, 0))
    }
    fn frob(&self, a: ((u64, u64), (u64, u64))) -> ((u64, u64), (u64, u64)) {
        Fp4::frob(self, a)
    }
    fn down(&self, a: ((u64, u64), (u64, u64))) -> Option<u64> {
        (a.1 == (0, 0) && a.0 .1 == 0).then_some(a.0 .0)
    }
}

impl ExtField for FpR {
    type B = Zp;
    fn base_field(&self) -> Zp {
        Zp::new(self.p)
    }
    fn lift(&self, a: u64) -> ER {
        self.embed(a)
    }
    fn frob(&self, a: ER) -> ER {
        self.frobenius(a)
    }
    fn down(&self, a: ER) -> Option<u64> {
        self.project(a)
    }
}

/// The class of s in (Z/lZ) / {±1}, as a representative in [0, (l-1)/2].
pub fn canon(s: i64, ell: u64) -> u64 {
    let r = s.rem_euclid(ell as i64) as u64;
    r.min(ell - r)
}

/// L-orbits of (Z/lZ)^* / {±1} under multiplication by lambda: each orbit as its list of
/// classes c, lambda c, lambda^2 c, ... (first element = smallest class), orbits by first element.
pub fn orbits(ell: u64, lambda: u64) -> Vec<Vec<u64>> {
    let half = (ell - 1) / 2;
    let mut seen = vec![false; half as usize + 1];
    let mut out = vec![];
    for c in 1..=half {
        if seen[c as usize] {
            continue;
        }
        let mut o = vec![];
        let mut x = c;
        while !seen[x as usize] {
            seen[x as usize] = true;
            o.push(x);
            x = canon((x as u128 * lambda as u128 % ell as u128) as i64, ell);
        }
        out.push(o);
    }
    out
}

/// An L-stable index system: one representative class per orbit of I, J and K (the orbits
/// themselves follow from lambda). (I, J) satisfies Definition 4.6 on classes: the classes of
/// i + j and i - j over all (i, j) in I x J are non-zero and pairwise distinct; K is the rest.
#[derive(Clone, Debug)]
pub struct IndexSystem {
    pub e: usize,
    pub i: Vec<u64>,
    pub j: Vec<u64>,
    pub k: Vec<u64>,
}

/// Greedy L-stable index system with J = the first `nj` orbits (see the module documentation).
pub fn index_system_with(ell: u64, lambda: u64, nj: usize) -> IndexSystem {
    let orbs = orbits(ell, lambda);
    let e = orbs[0].len();
    assert!(orbs.iter().all(|o| o.len() == e), "orbits of unequal size");
    let nj = nj.min(orbs.len());
    let jcl: Vec<i64> = orbs[..nj].iter().flatten().map(|&c| c as i64).collect();
    let mut used: HashSet<u64> = HashSet::new();
    let mut i = vec![];
    for o in &orbs {
        let mut t: HashSet<u64> = HashSet::with_capacity(2 * e * jcl.len());
        let mut ok = true;
        'pairs: for &c in o {
            for &jj in &jcl {
                for v in [canon(c as i64 + jj, ell), canon(c as i64 - jj, ell)] {
                    if v == 0 || used.contains(&v) || !t.insert(v) {
                        ok = false;
                        break 'pairs;
                    }
                }
            }
        }
        if ok {
            used.extend(t);
            i.push(o[0]);
        }
    }
    let k = orbs
        .iter()
        .filter(|o| !used.contains(&o[0]))
        .map(|o| o[0])
        .collect();
    if i.is_empty() {
        // no I orbit fits: everything is in K
        return IndexSystem {
            e,
            i: vec![],
            j: vec![],
            k: orbs.iter().map(|o| o[0]).collect(),
        };
    }
    IndexSystem {
        e,
        i,
        j: orbs[..nj].iter().map(|o| o[0]).collect(),
        k,
    }
}

/// The greedy system minimising max(e #I', 2 e #J', #K classes) over a few n_J around
/// sqrt((l-1) / (4 e^2)).
pub fn index_system(ell: u64, lambda: u64) -> IndexSystem {
    let half = (ell - 1) / 2;
    let e = orbits(ell, lambda)[0].len() as u64;
    let t = (((half as f64) / (2.0 * (e * e) as f64)).sqrt() as usize).max(1);
    let mut best: Option<(u64, IndexSystem)> = None;
    for nj in [t / 2, (2 * t) / 3, t, (3 * t) / 2]
        .into_iter()
        .filter(|&v| v > 0)
    {
        let s = index_system_with(ell, lambda, nj);
        let cost = (e * s.i.len() as u64)
            .max(2 * e * s.j.len() as u64)
            .max(e * s.k.len() as u64);
        if best.as_ref().is_none_or(|b| cost < b.0) {
            best = Some((cost, s));
        }
    }
    best.unwrap().1
}

/// The eigenvalue class of Frobenius on <P>: the s in [1, (l-1)/2] with x([s]P) = x(P)^q, by
/// walking the multiples (O(l) extension operations: a setup/test helper, not part of the
/// Õ(√l) evaluation; in applications lambda is known from the Frobenius eigenvalue).
pub fn find_eigenvalue<G: ExtField>(
    g: &G,
    a: <G::B as Field>::E,
    xp: G::E,
    ell: u64,
) -> Option<u64> {
    let e24 = a24(g, g.lift(a));
    let p1: XZ<G::E> = (xp, g.one());
    let target = g.frob(xp);
    let p2 = super::montgomery::xdbl(g, e24, p1);
    let (mut prev, mut cur) = (p1, p1);
    for s in 1..=(ell - 1) / 2 {
        if s == 2 {
            (prev, cur) = (cur, p2);
        } else if s > 2 {
            let nx = xadd(g, cur, p1, prev);
            (prev, cur) = (cur, nx);
        }
        if !g.is_zero(cur.1) && g.mul(target, cur.1) == cur.0 {
            return Some(s);
        }
    }
    None
}

/// √élu for a Galois-stable kernel <P> of odd prime order l with x(P) in an extension G of the
/// base field F_q, A in F_q (§4.14).
pub struct SqrtVeluIrr<G: ExtField> {
    pub g: G,
    pub base: G::B,
    /// Montgomery coefficient (in F_q)
    pub a: <G::B as Field>::E,
    pub ell: u64,
    /// Frobenius eigenvalue class (pi(P) = [±lambda] P)
    pub lambda: u64,
    pub sys: IndexSystem,
    /// x([s]P) for the representative s of each orbit of I, J, K (in G)
    xi: Vec<G::E>,
    xj: Vec<G::E>,
    xk: Vec<G::E>,
    /// x([s]P) for one class per orbit of the whole of S (for the orbit-norm product)
    xs: Vec<G::E>,
    /// remainder tree over the leaves g_O, O in I (its root is h_I), None when I is empty
    tree: Option<RemTree<G::B>>,
}

impl<G: ExtField> SqrtVeluIrr<G> {
    /// Setup from x(P) in G (affine) on E_A with A in F_q, the degree l and the eigenvalue class
    /// lambda; the index system is [`index_system`]. None if P does not have order l or
    /// x([lambda]P) != x(P)^q.
    pub fn new(g: &G, a: <G::B as Field>::E, xp: G::E, ell: u64, lambda: u64) -> Option<Self> {
        Self::with_system(g, a, xp, ell, lambda, index_system(ell, lambda))
    }

    pub fn with_system(
        g: &G,
        a: <G::B as Field>::E,
        xp: G::E,
        ell: u64,
        lambda: u64,
        sys: IndexSystem,
    ) -> Option<Self> {
        let base = g.base_field();
        let e24 = a24(g, g.lift(a));
        let p1: XZ<G::E> = (xp, g.one());
        if !g.is_zero(ladder_xz(g, e24, p1, &Big::from_u64(ell)).1) {
            return None;
        }
        let xmul = |s: u64| -> Option<G::E> {
            let r = ladder_xz(g, e24, p1, &Big::from_u64(s));
            (!g.is_zero(r.1)).then(|| g.div(r.0, r.1))
        };
        if xmul(lambda)? != g.frob(xp) {
            return None;
        }
        let orbs = orbits(ell, lambda);
        let xi: Vec<G::E> = sys.i.iter().map(|&s| xmul(s)).collect::<Option<_>>()?;
        let xj: Vec<G::E> = sys.j.iter().map(|&s| xmul(s)).collect::<Option<_>>()?;
        let xk: Vec<G::E> = sys.k.iter().map(|&s| xmul(s)).collect::<Option<_>>()?;
        let xs: Vec<G::E> = orbs.iter().map(|o| xmul(o[0])).collect::<Option<_>>()?;
        let e = sys.e;
        let tree = if xi.is_empty() {
            None
        } else {
            let leaves: Vec<Poly<G::B>> = xi
                .iter()
                .map(|&x| {
                    let roots = conjugates(g, x, e);
                    let lin: Vec<Poly<G>> =
                        roots.iter().map(|&r| vec![g.neg(r), g.one()]).collect();
                    down_poly(g, &product(g, lin))
                })
                .collect();
            Some(RemTree::new(&base, leaves))
        };
        Some(SqrtVeluIrr {
            g: g.clone(),
            base,
            a,
            ell,
            lambda,
            sys,
            xi,
            xj,
            xk,
            xs,
            tree,
        })
    }

    /// (#I, #J, #K) in classes of S (each orbit has e classes).
    pub fn sizes(&self) -> (usize, usize, usize) {
        let e = self.sys.e;
        (e * self.xi.len(), e * self.xj.len(), e * self.xk.len())
    }

    fn norm(&self, v: G::E) -> <G::B as Field>::E {
        let g = &self.g;
        let n = conjugates(g, v, self.sys.e)
            .into_iter()
            .fold(g.one(), |acc, c| g.mul(acc, c));
        g.down(n).expect("norm not in the base field")
    }

    /// E_J(alpha, Z) in F_q[Z] (length 2 e #J' + 1, untrimmed).
    fn e_j(&self, alpha: <G::B as Field>::E) -> Poly<G::B> {
        let g = &self.g;
        let (al, two_a) = (g.lift(alpha), g.lift(self.base.add(self.a, self.a)));
        let al2 = g.sq(al);
        let layer: Vec<Poly<G::B>> = self
            .xj
            .iter()
            .map(|&x| {
                // c2 = (alpha - x)^2, c0 = (x alpha - 1)^2,
                // c1 = -2 (x alpha^2 + (x^2 + 1 + 2 A x) alpha + x)
                let c2 = g.sq(g.sub(al, x));
                let c0 = g.sq(g.sub(g.mul(x, al), g.one()));
                let t = g.add(g.add(g.sq(x), g.one()), g.mul(two_a, x));
                let s = g.add(g.add(g.mul(x, al2), g.mul(t, al)), x);
                let q: Poly<G> = vec![c0, g.neg(g.add(s, s)), c2];
                let mut prod = q.clone();
                let mut qa = q;
                for _ in 1..self.sys.e {
                    qa = qa.iter().map(|&c| g.frob(c)).collect();
                    prod = poly::mul_raw(g, &prod, &qa);
                }
                down_poly(g, &prod)
            })
            .collect();
        product(&self.base, layer)
    }

    /// prod_{O in I} N(r(x_O)) for every r in `polys` (r reduced modulo the g_O by the tree).
    fn res_i(&self, polys: &[Poly<G::B>]) -> Vec<<G::B as Field>::E> {
        let (g, b) = (&self.g, &self.base);
        let Some(tree) = &self.tree else {
            return vec![b.one(); polys.len()];
        };
        polys
            .iter()
            .map(|p| {
                let mut p = p.clone();
                poly::trim(b, &mut p);
                tree.rems(b, &p)
                    .iter()
                    .zip(&self.xi)
                    .fold(b.one(), |acc, (r, &x)| {
                        let v = r
                            .iter()
                            .rev()
                            .fold(g.zero(), |s, &c| g.add(g.mul(s, x), g.lift(c)));
                        b.mul(acc, self.norm(v))
                    })
            })
            .collect()
    }

    /// (prod_S (1 - alpha x_s), h_S(alpha)), both times Δ_{I,J} (which cancels in every use).
    pub fn eval_nd(&self, alpha: <G::B as Field>::E) -> (<G::B as Field>::E, <G::B as Field>::E) {
        let (g, b) = (&self.g, &self.base);
        let e = self.e_j(alpha);
        let mut r = e.clone();
        r.reverse();
        let v = self.res_i(&[r, e]);
        let (mut num, mut den) = (v[0], v[1]);
        let al = g.lift(alpha);
        for &x in &self.xk {
            den = b.mul(den, self.norm(g.sub(al, x)));
            num = b.mul(num, self.norm(g.sub(g.one(), g.mul(al, x))));
        }
        (num, den)
    }

    /// h_S(alpha) Δ_{I,J} for alpha = ±1 (K included).
    fn h_scaled(&self, alpha: <G::B as Field>::E) -> <G::B as Field>::E {
        let (g, b) = (&self.g, &self.base);
        let mut r = self.res_i(&[self.e_j(alpha)])[0];
        let al = g.lift(alpha);
        for &x in &self.xk {
            r = b.mul(r, self.norm(g.sub(al, x)));
        }
        r
    }

    /// h_S(alpha) = prod over orbits of N(alpha - x_O): (l-1)/(2e) norms, no index system.
    pub fn h_orbits(&self, alpha: <G::B as Field>::E) -> <G::B as Field>::E {
        let g = &self.g;
        let al = g.lift(alpha);
        self.xs.iter().fold(self.base.one(), |acc, &x| {
            self.base.mul(acc, self.norm(g.sub(al, x)))
        })
    }

    /// Codomain coefficient A' in F_q from h_S(±1) (§5.2, Edwards form).
    pub fn codomain(&self) -> <G::B as Field>::E {
        let b = &self.base;
        let one = b.one();
        codomain_from_h(
            b,
            self.a,
            self.ell,
            self.h_scaled(one),
            self.h_scaled(b.neg(one)),
        )
    }

    /// Codomain from the orbit-norm products (reference path, O(l / e) norms).
    pub fn codomain_orbits(&self) -> <G::B as Field>::E {
        let b = &self.base;
        let one = b.one();
        codomain_from_h(
            b,
            self.a,
            self.ell,
            self.h_orbits(one),
            self.h_orbits(b.neg(one)),
        )
    }

    /// x-coordinate of the image of the point with abscissa alpha in F_q (not in the kernel).
    pub fn eval(&self, alpha: <G::B as Field>::E) -> <G::B as Field>::E {
        let b = &self.base;
        let (num, den) = self.eval_nd(alpha);
        b.mul(alpha, b.sq(b.div(num, den)))
    }
}

/// A' = 2 (a' + d') / (a' - d') with a' = (A+2)^l h(-1)^8, d' = (A-2)^l h(1)^8.
fn codomain_from_h<F: Field>(f: &F, a: F::E, ell: u64, h1: F::E, hm1: F::E) -> F::E {
    let two = f.from_u64(2);
    let p8 = |v: F::E| f.sq(f.sq(f.sq(v)));
    let d_new = f.mul(f.pow(f.sub(a, two), ell as u128), p8(h1));
    let a_new = f.mul(f.pow(f.add(a, two), ell as u128), p8(hm1));
    f.div(f.mul(two, f.add(a_new, d_new)), f.sub(a_new, d_new))
}

/// v, v^q, ..., v^{q^{e-1}}.
fn conjugates<G: ExtField>(g: &G, v: G::E, e: usize) -> Vec<G::E> {
    let mut out = Vec::with_capacity(e);
    let mut c = v;
    for _ in 0..e {
        out.push(c);
        c = g.frob(c);
    }
    out
}

/// Coefficients of a polynomial over G that lie in F_q.
fn down_poly<G: ExtField>(g: &G, p: &Poly<G>) -> Poly<G::B> {
    p.iter()
        .map(|&c| {
            g.down(c)
                .expect("orbit product not defined over the base field")
        })
        .collect()
}

/// Product of polynomials by a balanced tree (untrimmed: exact degrees add).
fn product<F: Field>(f: &F, mut layer: Vec<Poly<F>>) -> Poly<F> {
    while layer.len() > 1 {
        layer = layer
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    poly::mul_raw(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    layer.pop().unwrap_or_else(|| vec![f.one()])
}

// ---------------------------------------------------------------------------------------------
// Appendix A.5

/// Codomain and the image of one point from a single pass over I x J (A.5, first bullet: "merge
/// the four computations of Q(alpha), Q(1/alpha), Q(1), Q(-1)"). For each x_j the shared
/// t_j = x_j^2 + 1 + 2 A x_j feeds all four quadratics, and each pair (i, j) updates the four
/// products h(1), h(-1) (in the self-reciprocal W = x_i + 1/x_i form, A.5 second bullet, as in
/// `sqrt_velu_mont`), h(alpha) and prod (1 - alpha x) in one loop. Pair-by-pair products, as
/// `SqrtVeluMont` uses below `TREE_MIN_BP`. Returns (A', x(phi(alpha))).
pub fn fused_codomain_eval<F: Field>(f: &F, sv: &SqrtVeluMont<F>, alpha: F::E) -> (F::E, F::E) {
    let two_a = f.add(sv.a, sv.a);
    let (one, two) = (f.one(), f.from_u64(2));
    let al2 = f.sq(alpha);
    let xs_i2: Vec<F::E> = sv.xs_i.iter().map(|&x| f.sq(x)).collect();
    // per j: coefficients of the four quadratics; for alpha = ±1, c0 = c2 = (x -+ 1)^2 and the
    // factor is Z (c2 W + c1)
    let (mut hp, mut hm, mut num, mut den) = (one, one, one, one);
    for &xj in &sv.xs_j {
        let t = f.add(f.add(f.sq(xj), one), f.mul(two_a, xj));
        // alpha = 1: c1 = -2 (2 x + t), alpha = -1: c1 = -2 (2 x - t)
        let (xx, tt) = (f.add(xj, xj), t);
        let c1p = f.neg(f.mul(two, f.add(xx, tt)));
        let c1m = f.neg(f.mul(two, f.sub(xx, tt)));
        let c2p = f.sq(f.sub(one, xj));
        let c2m = f.sq(f.add(one, xj));
        // general alpha
        let c2 = f.sq(f.sub(alpha, xj));
        let c0 = f.sq(f.sub(f.mul(xj, alpha), one));
        let s = f.add(f.add(f.mul(xj, al2), f.mul(t, alpha)), xj);
        let c1 = f.neg(f.add(s, s));
        for (&w, (&x, &x2)) in sv.ws.iter().zip(sv.xs_i.iter().zip(&xs_i2)) {
            hp = f.mul(hp, f.add(f.mul(c2p, w), c1p));
            hm = f.mul(hm, f.add(f.mul(c2m, w), c1m));
            let u = f.mul(c1, x);
            den = f.mul(den, f.add(f.add(f.mul(c2, x2), u), c0));
            num = f.mul(num, f.add(f.add(f.mul(c0, x2), u), c2));
        }
    }
    for &xk in &sv.xs_k {
        hp = f.mul(hp, f.sub(one, xk));
        hm = f.mul(hm, f.sub(f.neg(one), xk));
        den = f.mul(den, f.sub(alpha, xk));
        num = f.mul(num, f.sub(one, f.mul(alpha, xk)));
    }
    let a_new = codomain_from_h(f, sv.a, sv.ell, hp, hm);
    (a_new, f.mul(alpha, f.sq(f.div(num, den))))
}

/// The Montgomery √élu resultants by remainder trees with precomputed node reciprocals and
/// Karatsuba low products (`strassen::RemTree`; A.5 bullets "precompute the product trees",
/// "precomputing reciprocals of tree nodes", "low products"), instead of `poly::EvalTree`, whose
/// scaled-remainder-tree middle products are schoolbook dot products (Θ(b'^2) per evaluation)
/// and whose root quotient is a quadratic series recurrence.
pub struct RemTreeSqrtVelu<F: Field> {
    /// tree over the x_i (point images)
    tree_x: Option<RemTree<F>>,
    /// tree over w_i = x_i + 1/x_i (codomain, self-reciprocal form)
    tree_w: Option<RemTree<F>>,
}

impl<F: Field> RemTreeSqrtVelu<F> {
    pub fn new(f: &F, sv: &SqrtVeluMont<F>) -> Self {
        if sv.bp == 0 {
            return RemTreeSqrtVelu {
                tree_x: None,
                tree_w: None,
            };
        }
        RemTreeSqrtVelu {
            tree_x: Some(RemTree::from_roots(f, &sv.xs_i)),
            tree_w: Some(RemTree::from_roots(f, &sv.ws)),
        }
    }

    fn prod(f: &F, t: &RemTree<F>, p: &Poly<F>) -> F::E {
        let mut p = p.clone();
        poly::trim(f, &mut p);
        t.eval(f, &p).into_iter().fold(f.one(), |a, v| f.mul(a, v))
    }

    /// (prod_S (1 - alpha x_s), h_S(alpha)) times Δ, as `SqrtVeluMont::eval_nd`.
    pub fn eval_nd(&self, f: &F, sv: &SqrtVeluMont<F>, alpha: F::E) -> (F::E, F::E) {
        let (mut num, mut den) = (f.one(), f.one());
        if let Some(t) = &self.tree_x {
            let e = sv.e_j(f, alpha);
            let mut r = e.clone();
            r.reverse();
            den = Self::prod(f, t, &e);
            num = Self::prod(f, t, &r);
        }
        for &x in &sv.xs_k {
            den = f.mul(den, f.sub(alpha, x));
            num = f.mul(num, f.sub(f.one(), f.mul(alpha, x)));
        }
        (num, den)
    }

    pub fn eval(&self, f: &F, sv: &SqrtVeluMont<F>, alpha: F::E) -> F::E {
        let (num, den) = self.eval_nd(f, sv, alpha);
        f.mul(alpha, f.sq(f.div(num, den)))
    }

    /// h_S(alpha) Δ (prod x_i)^b for alpha = ±1 via G(W), E_J(Z) = Z^b G(Z + 1/Z).
    fn h_scaled(&self, f: &F, sv: &SqrtVeluMont<F>, alpha: F::E) -> F::E {
        let mut r = f.one();
        if let Some(t) = &self.tree_w {
            let two_a = f.add(sv.a, sv.a);
            let layer: Vec<Poly<F>> = sv
                .xs_j
                .iter()
                .map(|&x| {
                    let c2 = f.sq(f.sub(alpha, x));
                    let tt = f.add(f.add(f.sq(x), f.one()), f.mul(two_a, x));
                    let s = f.add(f.add(x, x), f.mul(tt, alpha));
                    vec![f.neg(f.add(s, s)), c2]
                })
                .collect();
            r = Self::prod(f, t, &product(f, layer));
        }
        for &x in &sv.xs_k {
            r = f.mul(r, f.sub(alpha, x));
        }
        r
    }

    /// Codomain coefficient A'.
    pub fn codomain(&self, f: &F, sv: &SqrtVeluMont<F>) -> F::E {
        let one = f.one();
        let h1 = self.h_scaled(f, sv, one);
        let hm1 = self.h_scaled(f, sv, f.neg(one));
        codomain_from_h(f, sv.a, sv.ell, h1, hm1)
    }
}
