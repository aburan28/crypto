//! (2, ..., 2)-isogenies of principally polarised abelian varieties of any dimension g in level-2
//! theta coordinates (2^g coordinates, index bit k = characteristic component k), the
//! generalisation of `theta.rs` (g = 2) used for dimension-4 Kani embeddings (Robert 2022;
//! Dartois–Leroux–Robert–Wesolowski 2023; Dartois 2024).
//!
//! The formulas are those of `theta.rs` with the Hadamard transform over (Z/2)^g:
//!  * f: A -> B with kernel K2 (points acting by sign changes): H(S(x)) = H(theta^B(f x)) * D,
//!    D = H(theta^B(0)); D_{chi + e_j} / D_chi is read off H(S(T''_j)) for points T''_j of order 8
//!    above the kernel (no square roots). When A is a product some D_chi vanish; the missing
//!    dual coordinates of f(x) are recovered from f(x + T'_m) (T'_m = sum of the 4-torsion points
//!    above the kernel over the bits of m), which permutes the dual coordinates by m.
//!  * Doubling: theta(2x) = H(S(H(S(x))) / H(S(null))) / null.
//!  * Starting structure on E_1 x ... x E_g: the product of one-dimensional structures, then a
//!    change of coordinates making the kernel the K2 subgroup: the theta group of level 2 acts
//!    on coordinates by signed permutations (per factor: Z for the "1/2" point, X for the
//!    "tau/2" point); the new basis is the common +1 eigenvector w_0 of the kernel generators'
//!    operators, moved by the operators of a symplectic complement (w_chi = B^chi w_0). The
//!    remaining sign choices (which fix the level-4 structure) are searched until the first
//!    step is consistent (zero pattern of the 4-torsion points, D^2 = H(S(null))).
use crate::curve::{padd, Curve, Pt};
use crate::field::Field;
use crate::theta::{theta1_structure, Theta1};

/// Walsh–Hadamard transform: H(u)_chi = sum_c (-1)^{popcount(chi & c)} u_c.
pub fn wht<F: Field>(f: &F, x: &[F::E]) -> Vec<F::E> {
    let mut v = x.to_vec();
    let n = v.len();
    let mut h = 1;
    while h < n {
        for i in (0..n).step_by(2 * h) {
            for j in i..i + h {
                let (a, b) = (v[j], v[j + h]);
                v[j] = f.add(a, b);
                v[j + h] = f.sub(a, b);
            }
        }
        h *= 2;
    }
    v
}

pub fn sq_all<F: Field>(f: &F, x: &[F::E]) -> Vec<F::E> {
    x.iter().map(|&v| f.sq(v)).collect()
}

pub fn proj_eq_v<F: Field>(f: &F, x: &[F::E], y: &[F::E]) -> bool {
    let Some(k) = (0..x.len()).find(|&i| !f.is_zero(x[i])) else {
        return y.iter().all(|&v| f.is_zero(v));
    };
    if f.is_zero(y[k]) {
        return false;
    }
    (0..x.len()).all(|i| f.mul(x[k], y[i]) == f.mul(y[k], x[i]))
}

/// Projective inverse skipping zero coordinates (whose entry stays zero).
fn proj_inv_skip<F: Field>(f: &F, v: &[F::E]) -> Vec<F::E> {
    let n = v.len();
    let w: Vec<F::E> = v
        .iter()
        .map(|&c| if f.is_zero(c) { f.one() } else { c })
        .collect();
    let mut pre = vec![f.one(); n + 1];
    for i in 0..n {
        pre[i + 1] = f.mul(pre[i], w[i]);
    }
    let mut suf = vec![f.one(); n + 1];
    for i in (0..n).rev() {
        suf[i] = f.mul(suf[i + 1], w[i]);
    }
    (0..n)
        .map(|i| {
            if f.is_zero(v[i]) {
                f.zero()
            } else {
                f.mul(pre[i], suf[i + 1])
            }
        })
        .collect()
}

/// A (2,...,2)-isogeny with kernel K2 of the domain structure.
#[derive(Clone, Debug)]
pub struct IsoG<E> {
    pub dinv: Vec<E>,
    pub zeros: Vec<usize>,
    pub codomain: Vec<E>,
}

/// The isogeny from the coordinates of T''_j (order 8, 2 T''_j = e_j / 4 + K2), j = 0..g.
/// None if inconsistent (two paths of the hypercube give different D beyond a sign, or D^2 is
/// not proportional to H(S(null))).
pub fn isogeny_g<F: Field>(f: &F, null: &[F::E], t8: &[Vec<F::E>]) -> Option<IsoG<F::E>> {
    let masked: Vec<(usize, Vec<F::E>)> = t8
        .iter()
        .enumerate()
        .map(|(j, t)| (1usize << j, t.clone()))
        .collect();
    isogeny_g_masks(f, null, &masked)
}

/// As `isogeny_g` with points T''_m = sum_{j in m} T''_j for arbitrary masks m (each gives
/// D_{chi ^ m} / D_chi = x_m[chi ^ m] / x_m[chi]); on products the single-bit relations can
/// leave coordinates of D unreachable.
pub fn isogeny_g_masks<F: Field>(
    f: &F,
    null: &[F::E],
    t8: &[(usize, Vec<F::E>)],
) -> Option<IsoG<F::E>> {
    let n = null.len();
    let xs: Vec<(usize, Vec<F::E>)> = t8
        .iter()
        .map(|(m, t)| (*m, wht(f, &sq_all(f, t))))
        .collect();
    let x0 = wht(f, &sq_all(f, null));
    // D_chi = 0 exactly where H(S(null)) vanishes
    let b = (0..n).find(|&c| !f.is_zero(x0[c]) && xs.iter().any(|(_, x)| !f.is_zero(x[c])))?;
    let mut d: Vec<Option<F::E>> = (0..n)
        .map(|c| {
            if f.is_zero(x0[c]) {
                Some(f.zero())
            } else {
                None
            }
        })
        .collect();
    d[b] = Some(f.one());
    let mut queue = vec![b];
    while let Some(c) = queue.pop() {
        let dc = d[c].unwrap();
        if f.is_zero(dc) {
            continue;
        }
        for (m, x) in xs.iter() {
            if f.is_zero(x[c]) {
                continue;
            }
            let nb = c ^ m;
            let val = f.div(f.mul(dc, x[nb]), x[c]);
            match d[nb] {
                Some(v) => {
                    // a sign-only conflict is tolerated (it happens on the last step into a
                    // product, where the dual coordinates of f(T'') have complementary zeros)
                    if v != val && v != f.neg(val) {
                        return None;
                    }
                }
                None => {
                    d[nb] = Some(val);
                    queue.push(nb);
                }
            }
        }
    }
    let d: Vec<F::E> = d.into_iter().collect::<Option<Vec<_>>>()?;
    if !proj_eq_v(f, &sq_all(f, &d), &x0) {
        return None;
    }
    let zeros: Vec<usize> = (0..n).filter(|&i| f.is_zero(d[i])).collect();
    if zeros.len() == n {
        return None;
    }
    Some(IsoG {
        dinv: proj_inv_skip(f, &d),
        zeros,
        codomain: wht(f, &d),
    })
}

impl<E: Copy> IsoG<E> {
    fn dual<F: Field<E = E>>(&self, f: &F, x: &[E]) -> Vec<E> {
        let y = wht(f, &sq_all(f, x));
        y.iter()
            .zip(self.dinv.iter())
            .map(|(&a, &b)| f.mul(a, b))
            .collect()
    }
    pub fn eval<F: Field<E = E>>(&self, f: &F, x: &[E]) -> Vec<E> {
        wht(f, &self.dual(f, x))
    }
    /// Image of x given images-to-be of x + T'_m for the masks m in `trans` (needed when some
    /// dual coordinates vanish).
    pub fn eval_translates<F: Field<E = E>>(
        &self,
        f: &F,
        x: &[E],
        trans: &[(usize, Vec<E>)],
    ) -> Option<Vec<E>> {
        let mut y = self.dual(f, x);
        if self.zeros.is_empty() {
            return Some(wht(f, &y));
        }
        let duals: Vec<(usize, Vec<E>)> =
            trans.iter().map(|(m, xt)| (*m, self.dual(f, xt))).collect();
        let y0 = y.clone();
        for &z in &self.zeros {
            let mut done = false;
            for (m, yt) in &duals {
                if self.zeros.contains(&(z ^ m)) {
                    continue;
                }
                // normalisation coordinate c: Y(f(x+T'_m))_{c^m} = lambda Y(f x)_c
                let Some(c) = (0..y0.len()).find(|&c| {
                    !self.zeros.contains(&c)
                        && !self.zeros.contains(&(c ^ m))
                        && !f.is_zero(y0[c])
                        && !f.is_zero(yt[c ^ m])
                }) else {
                    continue;
                };
                y[z] = f.div(f.mul(yt[z ^ m], y0[c]), yt[c ^ m]);
                done = true;
                break;
            }
            if !done {
                return None;
            }
        }
        Some(wht(f, &y))
    }
}

/// Doubling on the Kummer variety (no zero in null or H(S(null))).
#[derive(Clone, Debug)]
pub struct DoublerG<E> {
    inv_x0: Vec<E>,
    inv_null: Vec<E>,
}

impl<E: Copy> DoublerG<E> {
    pub fn new<F: Field<E = E>>(f: &F, null: &[E]) -> Option<Self> {
        let x0 = wht(f, &sq_all(f, null));
        if x0.iter().chain(null.iter()).any(|&v| f.is_zero(v)) {
            return None;
        }
        Some(DoublerG {
            inv_x0: proj_inv_skip(f, &x0),
            inv_null: proj_inv_skip(f, null),
        })
    }
    pub fn double<F: Field<E = E>>(&self, f: &F, x: &[E]) -> Vec<E> {
        let xx = sq_all(f, &wht(f, &sq_all(f, x)));
        let y: Vec<E> = xx
            .iter()
            .zip(self.inv_x0.iter())
            .map(|(&a, &b)| f.mul(a, b))
            .collect();
        let z = wht(f, &y);
        z.iter()
            .zip(self.inv_null.iter())
            .map(|(&a, &b)| f.mul(a, b))
            .collect()
    }
}

/// theta[a/2, b/2](0)^2 for all (a, b) in (Z/2)^g x (Z/2)^g, index a * 2^g + b.
pub fn fundamental_squares_g<F: Field>(f: &F, n: &[F::E]) -> Vec<F::E> {
    let m = n.len();
    let mut out = vec![f.zero(); m * m];
    for a in 0..m {
        for b in 0..m {
            let mut s = f.zero();
            for c in 0..m {
                let t = f.mul(n[c], n[c ^ a]);
                s = if (b & c).count_ones() % 2 == 1 {
                    f.sub(s, t)
                } else {
                    f.add(s, t)
                };
            }
            out[a * m + b] = s;
        }
    }
    out
}

/// Evidence that the codomain B of a chain is a product: the level-2 null coordinates of the
/// last step's domain A are theta constants theta[c/2, 0](0, Omega_B) of B, so their zeros are
/// vanishing even theta constants of B (none for a generic principally polarised B), plus the
/// vanishing even constants of the computed final null.
pub fn split_score<F: Field>(f: &F, res: &ChainResultG<F::E>, first_null: &[F::E]) -> usize {
    let k = res.nulls.len();
    let dom = if k >= 2 {
        &res.nulls[k - 2]
    } else {
        first_null
    };
    dom.iter().filter(|&&v| f.is_zero(v)).count() + even_zero_count(f, &res.nulls[k - 1])
}

/// Number of vanishing even theta constants (0 for a generic ppav; > 0 for products).
pub fn even_zero_count<F: Field>(f: &F, n: &[F::E]) -> usize {
    let m = n.len();
    let t = fundamental_squares_g(f, n);
    (0..m * m)
        .filter(|&i| ((i / m) & (i % m)).count_ones() % 2 == 0 && f.is_zero(t[i]))
        .count()
}

// ------------------------------------------------------------------ the starting structure

/// A 2-torsion point of a product as per-factor (z, x) bits: z = component on the "1/2" point
/// T0 = 2 P4 (sign change), x = on the "tau/2" point T1 = 2 Q4 (swap).
type Tor2 = Vec<(bool, bool)>;

fn symp(a: &Tor2, b: &Tor2) -> bool {
    a.iter()
        .zip(b.iter())
        .fold(false, |acc, (u, v)| acc ^ (u.0 & v.1) ^ (u.1 & v.0))
}

/// Apply the operator of a 2-torsion point to a coordinate vector: per factor, Z then X.
fn apply_op<F: Field>(f: &F, v: &Tor2, scal: F::E, u: &[F::E]) -> Vec<F::E> {
    let n = u.len();
    let mut out = u.to_vec();
    for (k, &(z, x)) in v.iter().enumerate() {
        if z {
            for (c, val) in out.iter_mut().enumerate() {
                if (c >> k) & 1 == 1 {
                    *val = f.neg(*val);
                }
            }
        }
        if x {
            let prev = out.clone();
            for c in 0..n {
                out[c] = prev[c ^ (1 << k)];
            }
        }
    }
    out.iter().map(|&a| f.mul(a, scal)).collect()
}

/// Gaussian elimination inverse of a square matrix (rows).
fn mat_inv<F: Field>(f: &F, m: &[Vec<F::E>]) -> Option<Vec<Vec<F::E>>> {
    let n = m.len();
    let mut a: Vec<Vec<F::E>> = m.to_vec();
    let mut inv: Vec<Vec<F::E>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { f.one() } else { f.zero() })
                .collect()
        })
        .collect();
    for col in 0..n {
        let piv = (col..n).find(|&r| !f.is_zero(a[r][col]))?;
        a.swap(col, piv);
        inv.swap(col, piv);
        let pi = f.inv(a[col][col]);
        for j in 0..n {
            a[col][j] = f.mul(a[col][j], pi);
            inv[col][j] = f.mul(inv[col][j], pi);
        }
        for r in 0..n {
            if r != col && !f.is_zero(a[r][col]) {
                let fac = a[r][col];
                for j in 0..n {
                    a[r][j] = f.sub(a[r][j], f.mul(fac, a[col][j]));
                    inv[r][j] = f.sub(inv[r][j], f.mul(fac, inv[col][j]));
                }
            }
        }
    }
    Some(inv)
}

/// Theta structure on E_1 x ... x E_g in which a given kernel is K2.
#[derive(Clone, Debug)]
pub struct ProductStructure<E> {
    pub t: Vec<Theta1<E>>,
    /// change of coordinates (rows), new = m * product
    pub m: Vec<Vec<E>>,
}

impl<E: Copy> ProductStructure<E> {
    pub fn product_coords<F: Field<E = E>>(&self, f: &F, pts: &[Pt<E>]) -> Vec<E> {
        let g = pts.len();
        let per: Vec<[E; 2]> = (0..g).map(|k| self.t[k].point(f, &pts[k])).collect();
        (0..1usize << g)
            .map(|c| (0..g).fold(f.one(), |acc, k| f.mul(acc, per[k][(c >> k) & 1])))
            .collect()
    }
    pub fn point<F: Field<E = E>>(&self, f: &F, pts: &[Pt<E>]) -> Vec<E> {
        let p = self.product_coords(f, pts);
        self.m
            .iter()
            .map(|row| {
                row.iter()
                    .zip(p.iter())
                    .fold(f.zero(), |acc, (&a, &b)| f.add(acc, f.mul(a, b)))
            })
            .collect()
    }
    pub fn null<F: Field<E = E>>(&self, f: &F, g: usize) -> Vec<E> {
        self.point(f, &vec![Pt::Inf; g])
    }
}

/// Candidate structures (one per sign choice) for which the kernel generated by 2 T'_j is K2 and
/// T'_j has the zero pattern of e_j / 4. `iota` is a square root of -1.
pub fn product_structures<F: Field>(
    f: &F,
    curves: &[Curve<F::E>],
    t4: &[Vec<Pt<F::E>>],
    iota: F::E,
) -> Vec<ProductStructure<F::E>> {
    let g = curves.len();
    let n = 1usize << g;
    let dbl = |k: usize, p: &Pt<F::E>| padd(f, &curves[k], p, p);
    // per-factor 4-torsion basis taken from the components of the T'_j
    let mut t1s = Vec::with_capacity(g);
    for k in 0..g {
        let comps: Vec<Pt<F::E>> = t4
            .iter()
            .map(|t| t[k])
            .filter(|p| *p != Pt::Inf && dbl(k, p) != Pt::Inf)
            .collect();
        let mut found = None;
        'outer: for a in 0..comps.len() {
            for b in 0..comps.len() {
                if dbl(k, &comps[a]) != dbl(k, &comps[b]) {
                    if let Some(t) = theta1_structure(f, &curves[k], &comps[a], &comps[b]) {
                        found = Some((t, comps[a], comps[b]));
                        break 'outer;
                    }
                }
            }
        }
        let Some(fd) = found else { return vec![] };
        t1s.push(fd);
    }
    // 2-torsion classification per factor
    let classify = |k: usize, p: &Pt<F::E>| -> (bool, bool) {
        let t0 = dbl(k, &t1s[k].1);
        let t1 = dbl(k, &t1s[k].2);
        let x = |q: &Pt<F::E>| match q {
            Pt::Aff(x, _) => Some(*x),
            Pt::Inf => None,
        };
        if *p == Pt::Inf {
            (false, false)
        } else if x(p) == x(&t0) {
            (true, false)
        } else if x(p) == x(&t1) {
            (false, true)
        } else {
            (true, true)
        }
    };
    let s: Vec<Tor2> = t4
        .iter()
        .map(|t| (0..g).map(|k| classify(k, &dbl(k, &t[k]))).collect())
        .collect();
    // isotropy and independence
    for i in 0..g {
        for j in 0..g {
            if symp(&s[i], &s[j]) {
                return vec![];
            }
        }
    }
    // symplectic complement: R_j with <S_i, R_j> = delta_ij, <R_i, R_j> = 0
    let all: Vec<Tor2> = (0..1usize << (2 * g))
        .map(|bits| {
            (0..g)
                .map(|k| ((bits >> (2 * k)) & 1 == 1, (bits >> (2 * k + 1)) & 1 == 1))
                .collect()
        })
        .collect();
    let mut r: Vec<Tor2> = Vec::with_capacity(g);
    for j in 0..g {
        let cand = all
            .iter()
            .find(|v| (0..g).all(|i| symp(&s[i], v) == (i == j)));
        let Some(c) = cand else { return vec![] };
        let mut v = c.clone();
        for k in 0..j {
            if symp(&v, &r[k]) {
                v = v
                    .iter()
                    .zip(s[k].iter())
                    .map(|(a, b)| (a.0 ^ b.0, a.1 ^ b.1))
                    .collect();
            }
        }
        r.push(v);
    }
    // operators normalised to square to the identity
    let sq_sign = |v: &Tor2| v.iter().filter(|(z, x)| *z && *x).count() % 2 == 1;
    let scal = |v: &Tor2| if sq_sign(v) { iota } else { f.one() };
    let mut out = vec![];
    let mut rng = crate::field::Rng::new(0x7e7a);
    for eps in 0..1usize << g {
        for bsg in 0..1usize << g {
            let a_ops: Vec<(Tor2, F::E)> = (0..g)
                .map(|j| {
                    (
                        s[j].clone(),
                        if (eps >> j) & 1 == 1 {
                            f.neg(scal(&s[j]))
                        } else {
                            scal(&s[j])
                        },
                    )
                })
                .collect();
            let b_ops: Vec<(Tor2, F::E)> = (0..g)
                .map(|j| {
                    (
                        r[j].clone(),
                        if (bsg >> j) & 1 == 1 {
                            f.neg(scal(&r[j]))
                        } else {
                            scal(&r[j])
                        },
                    )
                })
                .collect();
            // common +1 eigenvector
            let mut w0: Vec<F::E> = (0..n).map(|_| f.random(&mut rng)).collect();
            for (v, sc) in &a_ops {
                let aw = apply_op(f, v, *sc, &w0);
                w0 = w0
                    .iter()
                    .zip(aw.iter())
                    .map(|(&x, &y)| f.add(x, y))
                    .collect();
            }
            if w0.iter().all(|&c| f.is_zero(c)) {
                continue;
            }
            let cols: Vec<Vec<F::E>> = (0..n)
                .map(|chi| {
                    let mut w = w0.clone();
                    for (j, (v, sc)) in b_ops.iter().enumerate() {
                        if (chi >> j) & 1 == 1 {
                            w = apply_op(f, v, *sc, &w);
                        }
                    }
                    w
                })
                .collect();
            // W has columns w_chi; new coordinates = W^-1 * product coordinates
            let wmat: Vec<Vec<F::E>> = (0..n)
                .map(|i| (0..n).map(|chi| cols[chi][i]).collect())
                .collect();
            let Some(m) = mat_inv(f, &wmat) else { continue };
            let st = ProductStructure {
                t: t1s.iter().map(|x| x.0).collect(),
                m,
            };
            // zero pattern of T'_j
            let ok = (0..g).all(|j| {
                let c = st.point(f, &t4[j]);
                (0..n)
                    .filter(|chi| (chi >> j) & 1 == 1)
                    .all(|chi| f.is_zero(c[chi]))
            });
            if ok {
                out.push(st);
            }
        }
    }
    out
}

// ------------------------------------------------------------------ chains

#[derive(Clone, Debug)]
pub struct ChainResultG<E> {
    /// theta null of every codomain. The last one is only reliable up to the choice of signs of
    /// its dual null: the 8-torsion points do not fix them on the last step into a product (the
    /// relations around a square of the hypercube multiply to -1 there), and for g > 2 a wrong
    /// sign pattern need not be a theta null at all. Use `split_score` to decide splitting.
    pub nulls: Vec<Vec<E>>,
    pub images: Vec<Vec<E>>,
    /// number of leading steps whose dual null had vanishing coordinates
    pub gluing_steps: usize,
}

struct GlueStep<E> {
    iso: IsoG<E>,
    masks: Vec<usize>,
}

/// Masks m so that every vanishing dual coordinate z has z ^ m non-vanishing for some m.
fn choose_masks(zeros: &[usize], g: usize) -> Vec<usize> {
    let mut masks = vec![];
    for &z in zeros {
        if masks.iter().any(|m| !zeros.contains(&(z ^ m))) {
            continue;
        }
        if let Some(m) = (1..1usize << g).find(|m| !zeros.contains(&(z ^ m))) {
            masks.push(m);
        }
    }
    masks
}

/// The (2^n, ..., 2^n)-isogeny of E_1 x ... x E_g with kernel generated by [4] K_j, where the
/// K_j (g points of the product, order 2^(n+2)) generate a maximal isotropic subgroup. Steps
/// whose dual null has zeros (gluings) are evaluated with translates computed on the product;
/// the rest by an optimal strategy over theta doublings. `iota` = sqrt(-1).
pub fn chain_g<F: Field>(
    f: &F,
    curves: &[Curve<F::E>],
    k: &[Vec<Pt<F::E>>],
    n: u32,
    extra: &[Vec<Pt<F::E>>],
    iota: F::E,
) -> Option<ChainResultG<F::E>> {
    chain_g_with(f, curves, k, n, extra, iota, None)
}

/// As `chain_g`, optionally with a given starting structure.
pub fn chain_g_with<F: Field>(
    f: &F,
    curves: &[Curve<F::E>],
    k: &[Vec<Pt<F::E>>],
    n: u32,
    extra: &[Vec<Pt<F::E>>],
    iota: F::E,
    fixed: Option<&ProductStructure<F::E>>,
) -> Option<ChainResultG<F::E>> {
    let g = curves.len();
    let n = n as usize;
    assert!(n >= 1 && k.len() == g);
    let addp = |p: &[Pt<F::E>], q: &[Pt<F::E>]| -> Vec<Pt<F::E>> {
        (0..g).map(|i| padd(f, &curves[i], &p[i], &q[i])).collect()
    };
    let mults: Vec<Vec<Vec<Pt<F::E>>>> = k
        .iter()
        .map(|kj| {
            let mut v = vec![kj.clone()];
            for t in 0..=n {
                let nx = addp(&v[t], &v[t]);
                v.push(nx);
            }
            v
        })
        .collect();
    let t4: Vec<Vec<Pt<F::E>>> = (0..g).map(|j| mults[j][n].clone()).collect();
    // offsets on the product for the translates of step s (1-based): sum over bits of m of
    // [2^(n + 1 - s)] K_j
    let offset = |s: usize, m: usize| -> Vec<Pt<F::E>> {
        let mut acc = vec![Pt::Inf; g];
        for j in 0..g {
            if (m >> j) & 1 == 1 {
                acc = addp(&acc, &mults[j][n + 1 - s]);
            }
        }
        acc
    };
    let sts = match fixed {
        Some(st) => vec![st.clone()],
        None => product_structures(f, curves, &t4, iota),
    };
    for st in sts {
        if let Some(r) = chain_with_structure(f, g, n, &st, &mults, extra, &offset, &addp) {
            return Some(r);
        }
    }
    None
}

/// Push a product point through the recorded gluing steps 1..=s.
fn push<F: Field>(
    f: &F,
    st: &ProductStructure<F::E>,
    steps: &[GlueStep<F::E>],
    s: usize,
    x: &[Pt<F::E>],
    offset: &dyn Fn(usize, usize) -> Vec<Pt<F::E>>,
    addp: &dyn Fn(&[Pt<F::E>], &[Pt<F::E>]) -> Vec<Pt<F::E>>,
) -> Option<Vec<F::E>> {
    if s == 0 {
        return Some(st.point(f, x));
    }
    let y = push(f, st, steps, s - 1, x, offset, addp)?;
    let step = &steps[s - 1];
    if step.iso.zeros.is_empty() {
        return Some(step.iso.eval(f, &y));
    }
    let mut trans = vec![];
    for &m in &step.masks {
        let xt = addp(x, &offset(s, m));
        trans.push((m, push(f, st, steps, s - 1, &xt, offset, addp)?));
    }
    step.iso.eval_translates(f, &y, &trans)
}

#[allow(clippy::too_many_arguments)]
fn chain_with_structure<F: Field>(
    f: &F,
    g: usize,
    n: usize,
    st: &ProductStructure<F::E>,
    mults: &[Vec<Vec<Pt<F::E>>>],
    extra: &[Vec<Pt<F::E>>],
    offset: &dyn Fn(usize, usize) -> Vec<Pt<F::E>>,
    addp: &dyn Fn(&[Pt<F::E>], &[Pt<F::E>]) -> Vec<Pt<F::E>>,
) -> Option<ChainResultG<F::E>> {
    let null = st.null(f, g);
    let mut steps: Vec<GlueStep<F::E>> = vec![];
    let mut nulls = vec![];
    // gluing phase: compute steps until one has no vanishing dual coordinate
    let mut s = 1;
    loop {
        // T''_m = sum_{j in m} [2^(n - s)] K_j for every non-zero mask (computed on the product)
        let mut t8 = vec![];
        for m in 1..1usize << g {
            let mut acc = vec![Pt::Inf; g];
            for j in 0..g {
                if (m >> j) & 1 == 1 {
                    acc = addp(&acc, &mults[j][n - s]);
                }
            }
            t8.push((m, push(f, st, &steps, s - 1, &acc, offset, addp)?));
        }
        let dom = if s == 1 {
            null.clone()
        } else {
            steps[s - 2].iso.codomain.clone()
        };
        let iso = isogeny_g_masks(f, &dom, &t8)?;
        let masks = choose_masks(&iso.zeros, g);
        let zero = !iso.zeros.is_empty();
        nulls.push(iso.codomain.clone());
        steps.push(GlueStep { iso, masks });
        if !zero || s == n {
            break;
        }
        s += 1;
    }
    let gl = steps.iter().filter(|x| !x.iso.zeros.is_empty()).count();
    let done = steps.len();
    let mut images: Vec<Vec<F::E>> = vec![];
    for x in extra {
        images.push(push(f, st, &steps, done, x, offset, addp)?);
    }
    if done == n {
        return Some(ChainResultG {
            nulls,
            images,
            gluing_steps: gl,
        });
    }
    // generic phase: generators of order 2^(n + 2 - done), m = n - done steps
    let gens: Vec<Vec<F::E>> = (0..g)
        .map(|j| push(f, st, &steps, done, &mults[j][0], offset, addp))
        .collect::<Option<_>>()?;
    let m = n - done;
    let splits = crate::kernel::two_power::optimal_splits(m, 1.0, 0.6);
    let mut cur = steps.last().unwrap().iso.codomain.clone();
    rec_strategy(f, &mut cur, gens, m, &mut images, &splits, &mut nulls)?;
    Some(ChainResultG {
        nulls,
        images,
        gluing_steps: gl,
    })
}

fn rec_strategy<F: Field>(
    f: &F,
    cur: &mut Vec<F::E>,
    gens: Vec<Vec<F::E>>,
    m: usize,
    stack: &mut Vec<Vec<F::E>>,
    splits: &[usize],
    nulls: &mut Vec<Vec<F::E>>,
) -> Option<()> {
    let g = gens.len();
    if m == 1 {
        let iso = isogeny_g(f, cur, &gens)?;
        if !iso.zeros.is_empty() {
            return None;
        }
        for x in stack.iter_mut() {
            *x = iso.eval(f, x);
        }
        *cur = iso.codomain.clone();
        nulls.push(cur.clone());
        return Some(());
    }
    let i = splits[m];
    let dbl = DoublerG::new(f, cur)?;
    let mut h = gens.clone();
    for _ in 0..(m - i) {
        h = h.iter().map(|x| dbl.double(f, x)).collect();
    }
    for x in gens {
        stack.push(x);
    }
    rec_strategy(f, cur, h, i, stack, splits, nulls)?;
    let mut back = vec![];
    for _ in 0..g {
        back.push(stack.pop().unwrap());
    }
    back.reverse();
    rec_strategy(f, cur, back, m - i, stack, splits, nulls)
}
