//! (2,2)-isogenies of abelian surfaces in level-2 theta coordinates, gluing E1 x E2 -> J(C),
//! detection of split codomains, and chains of (2,2)-isogenies with kernel of type (2^n, 2^n)
//! as used by Kani-lemma methods (Castryck–Decru, Maino–Martindale–Panny–Pope–Wesolowski,
//! Robert 2022; algorithms after Dartois–Maino–Pope–Robert 2023).
//!
//! Notation. A point x of an abelian surface A = C^2/(Z^2 + Omega Z^2) has level-2 theta
//! coordinates theta_c(x) = theta[c/2, 0](2x, 2 Omega), c in (Z/2)^2, stored at index 2 c1 + c2.
//! H is the Hadamard transform H(u)_chi = sum_c (-1)^{chi.c} u_c and S squares coordinatewise.
//! The duplication formula
//!     theta[a/2,b/2](z1, W) theta[a/2,b/2](z2, W)
//!         = sum_c (-1)^{b.c} theta[c/2,0](z1+z2, 2W) theta[(c+a)/2,0](z1-z2, 2W)
//! applied with W = 2 Omega, z1 = z2 = 2x gives, for f: A -> B = C^2/(Z^2 + 2 Omega Z^2),
//! z -> 2z (kernel K2 = (1/2)Z^2, the points acting on coordinates by sign changes),
//!     H(S(theta^A(x))) = H(theta^B(f x)) * D   (pointwise),   D = H(theta^B(0)),
//! so f(x) = H(H(S(theta^A(x))) / D) and theta^B(0) = H(D). D needs a square root of
//! H(S(theta^A(0))) unless points T'' of order 8 with 4 T'' generating the kernel are known:
//! f(T''_i) = e_i / 4, whose theta coordinates vanish where c_i = 1, which makes
//! H(theta^B(f T''_i)) invariant under chi -> chi + e_i and gives D_{chi + e_i} / D_chi as a
//! ratio of coordinates of H(S(theta^A(T''_i))) (no square root). The structure on B obtained
//! this way is compatible with the images of the next kernel, so chains need no basis changes.
//!
//! Gluing. On E1 x E2 with the product of one-dimensional structures, the kernel
//! {(T, psi(T))} of a gluing is not of type K2; the change of coordinates
//! u = (x00 + x11, x00 - x11, x01 + x10, x10 - x01) makes it so (it intertwines the action of
//! the 2-torsion). For a product one coordinate of D vanishes; the missing coordinate of an
//! image is recovered from the image of x + T' (T' of order 4 above the kernel, f(T') a
//! sign-type 2-torsion point), which permutes the dual coordinates.
use crate::curve::{padd, pmul, Curve, Pt};
use crate::field::Field;

/// Level-2 theta coordinates of a point of an abelian surface (index 2 c1 + c2).
pub type Th<E> = [E; 4];

pub fn hadamard<F: Field>(f: &F, x: &Th<F::E>) -> Th<F::E> {
    let (s01, d01) = (f.add(x[0], x[1]), f.sub(x[0], x[1]));
    let (s23, d23) = (f.add(x[2], x[3]), f.sub(x[2], x[3]));
    [
        f.add(s01, s23),
        f.add(d01, d23),
        f.sub(s01, s23),
        f.sub(d01, d23),
    ]
}

pub fn squared<F: Field>(f: &F, x: &Th<F::E>) -> Th<F::E> {
    [f.sq(x[0]), f.sq(x[1]), f.sq(x[2]), f.sq(x[3])]
}

/// Projective equality of two theta points.
pub fn proj_eq<F: Field>(f: &F, x: &Th<F::E>, y: &Th<F::E>) -> bool {
    (0..4).all(|i| (0..4).all(|j| f.mul(x[i], y[j]) == f.mul(x[j], y[i])))
}

// ------------------------------------------------------------------ dimension one

/// Level-2 theta structure on y^2 = x^3 + a x + b from a basis (P4, Q4) of E[4]:
/// x_m = kappa (x - e0) is a Montgomery x-coordinate with 2 P4 at 0 and P4 at 1; the theta
/// coordinates of x are (x_m + 1 : s (x_m - 1)), with null (1 : s), s = (x_m(Q4) + 1)/(x_m(Q4) - 1).
/// Then P4 = 1/4 and Q4 = tau/4: theta(P4) = (* : 0), theta(Q4) = (1 : 1).
#[derive(Clone, Copy, Debug)]
pub struct Theta1<E> {
    pub e0: E,
    pub kappa: E,
    pub s: E,
}

pub fn theta1_structure<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    p4: &Pt<F::E>,
    q4: &Pt<F::E>,
) -> Option<Theta1<F::E>> {
    let Pt::Aff(xp, _) = *p4 else { return None };
    let Pt::Aff(xq, _) = *q4 else { return None };
    let Pt::Aff(e0, _) = pmul(f, e, p4, 2) else {
        return None;
    };
    let d = f.sub(xp, e0);
    if f.is_zero(d) {
        return None;
    }
    let kappa = f.inv(d);
    let xm = f.mul(kappa, f.sub(xq, e0));
    let den = f.sub(xm, f.one());
    if f.is_zero(den) {
        return None;
    }
    let s = f.div(f.add(xm, f.one()), den);
    if f.is_zero(s) {
        return None;
    }
    Some(Theta1 { e0, kappa, s })
}

impl<E: Copy> Theta1<E> {
    pub fn null<F: Field<E = E>>(&self, f: &F) -> [E; 2] {
        [f.one(), self.s]
    }
    pub fn point<F: Field<E = E>>(&self, f: &F, p: &Pt<E>) -> [E; 2] {
        match *p {
            Pt::Inf => self.null(f),
            Pt::Aff(x, _) => {
                let xm = f.mul(self.kappa, f.sub(x, self.e0));
                [f.add(xm, f.one()), f.mul(self.s, f.sub(xm, f.one()))]
            }
        }
    }
}

impl<E: Copy> Theta1<E> {
    /// x-coordinate of a point from its theta coordinates [u : v] ~ [x_m + 1 : s(x_m - 1)]:
    /// x_m = (s u + v)/(s u - v), x = e0 + x_m / kappa. None at the identity (u = v = 0) or when
    /// s u = v (the 2-torsion point 2 Q4, x_m = infinity: returns None).
    pub fn x_of<F: Field<E = E>>(&self, f: &F, t: &[E; 2]) -> Option<E> {
        let su = f.mul(self.s, t[0]);
        let den = f.sub(su, t[1]);
        if f.is_zero(den) {
            return None;
        }
        let xm = f.div(f.add(su, t[1]), den);
        Some(f.add(self.e0, f.div(xm, self.kappa)))
    }
}

/// Factor the product theta coordinates [a0 b0 : a0 b1 : a1 b0 : a1 b1] (rank-1 2x2) into
/// ([a0 : a1], [b0 : b1]) projectively. None if all zero.
pub fn factor_product<F: Field>(f: &F, p: &Th<F::E>) -> Option<([F::E; 2], [F::E; 2])> {
    // rows (a0 b*, a1 b*): [a0 : a1] from a non-zero column; [b0 : b1] from a non-zero row
    let a = if !f.is_zero(p[0]) || !f.is_zero(p[2]) {
        [p[0], p[2]]
    } else {
        [p[1], p[3]]
    };
    let b = if !f.is_zero(p[0]) || !f.is_zero(p[1]) {
        [p[0], p[1]]
    } else {
        [p[2], p[3]]
    };
    if a.iter().chain(b.iter()).all(|&v| f.is_zero(v)) {
        return None;
    }
    Some((a, b))
}

/// Legendre lambda and j-invariant of an elliptic curve from its level-2 theta null (a : b):
/// theta[0,0]^2 = a^2 + b^2, theta[1/2,0]^2 = 2ab, lambda = theta[1/2,0]^4 / theta[0,0]^4.
pub fn j_from_theta1<F: Field>(f: &F, n: &[F::E; 2]) -> F::E {
    let (a, b) = (n[0], n[1]);
    let t3 = f.add(f.sq(a), f.sq(b));
    let t2 = f.mul(f.from_u64(2), f.mul(a, b));
    j_from_lambda(f, f.div(f.sq(t2), f.sq(t3)))
}

pub fn j_from_lambda<F: Field>(f: &F, l: F::E) -> F::E {
    let one = f.one();
    let num = f.sub(f.add(f.sq(l), one), l);
    let den = f.mul(f.sq(l), f.sq(f.sub(l, one)));
    f.div(f.mul(f.from_u64(256), f.mul(num, f.sq(num))), den)
}

// ------------------------------------------------------------------ products and gluing

/// Coordinates on E1 x E2 in the structure where {(T, psi(T))} (psi matching 2P4 -> 2P4',
/// 2Q4 -> 2Q4') is the K2 subgroup: the product coordinates followed by the gluing change.
#[derive(Clone, Copy, Debug)]
pub struct ProductTheta<E> {
    pub t1: Theta1<E>,
    pub t2: Theta1<E>,
}

fn glue_change<F: Field>(f: &F, x: &Th<F::E>) -> Th<F::E> {
    [
        f.add(x[0], x[3]),
        f.sub(x[0], x[3]),
        f.add(x[1], x[2]),
        f.sub(x[2], x[1]),
    ]
}

impl<E: Copy> ProductTheta<E> {
    /// Product coordinates (no change of basis).
    pub fn product_point<F: Field<E = E>>(&self, f: &F, p1: &Pt<E>, p2: &Pt<E>) -> Th<E> {
        let a = self.t1.point(f, p1);
        let b = self.t2.point(f, p2);
        [
            f.mul(a[0], b[0]),
            f.mul(a[0], b[1]),
            f.mul(a[1], b[0]),
            f.mul(a[1], b[1]),
        ]
    }
    pub fn point<F: Field<E = E>>(&self, f: &F, p1: &Pt<E>, p2: &Pt<E>) -> Th<E> {
        glue_change(f, &self.product_point(f, p1, p2))
    }
    pub fn null<F: Field<E = E>>(&self, f: &F) -> Th<E> {
        self.point(f, &Pt::Inf, &Pt::Inf)
    }
    /// Inverse of glue_change: product coordinates from the K2-structure coordinates.
    pub fn unglue<F: Field<E = E>>(f: &F, g: &Th<E>) -> Th<E> {
        [
            f.add(g[0], g[1]),
            f.sub(g[2], g[3]),
            f.add(g[2], g[3]),
            f.sub(g[0], g[1]),
        ]
    }
    /// The (x1, x2) of a point from its K2-structure coordinates (x-coordinates on the two
    /// factors; None where a factor coordinate is the identity or a troublesome 2-torsion point).
    pub fn factor_point<F: Field<E = E>>(&self, f: &F, g: &Th<E>) -> (Option<E>, Option<E>) {
        let prod = Self::unglue(f, g);
        match factor_product(f, &prod) {
            Some((a, b)) => (self.t1.x_of(f, &a), self.t2.x_of(f, &b)),
            None => (None, None),
        }
    }
}

// ------------------------------------------------------------------ the (2,2)-isogeny

/// A (2,2)-isogeny with kernel K2 of the domain structure.
#[derive(Clone, Debug)]
pub struct ThetaIso<E> {
    /// 1 / D_chi (zero where D_chi = 0)
    pub dinv: Th<E>,
    /// index of the vanishing dual coordinate (gluing), if any
    pub zero: Option<usize>,
    pub domain: Th<E>,
    pub codomain: Th<E>,
}

/// The isogeny with kernel K2 from the theta coordinates of T''_1, T''_2 (order 8, with
/// 2 T''_i = e_i / 4 + K2). None if the data is inconsistent (wrong structure or torsion):
/// D^2 must be proportional to H(S(null)). Projective throughout: no field inversion.
pub fn theta_isogeny<F: Field>(
    f: &F,
    null: &Th<F::E>,
    t1: &Th<F::E>,
    t2: &Th<F::E>,
) -> Option<ThetaIso<F::E>> {
    let x1 = hadamard(f, &squared(f, t1));
    let x2 = hadamard(f, &squared(f, t2));
    // D_{chi ^ 2} / D_chi = x1[chi ^ 2] / x1[chi],  D_{chi ^ 1} / D_chi = x2[chi ^ 1] / x2[chi];
    // with D_b = x1[b] x2[b] x1[b ^ 1] all four are products (b: a base index with the needed
    // coordinates non-zero; for a gluing one D vanishes and b avoids it)
    let b =
        (0..4usize).find(|&b| !f.is_zero(x1[b]) && !f.is_zero(x2[b]) && !f.is_zero(x1[b ^ 1]))?;
    let mut d = [f.zero(); 4];
    let u = f.mul(x1[b], x1[b ^ 1]);
    d[b] = f.mul(u, x2[b]);
    d[b ^ 1] = f.mul(u, x2[b ^ 1]);
    d[b ^ 2] = f.mul(f.mul(x1[b ^ 2], x1[b ^ 1]), x2[b]);
    d[b ^ 3] = f.mul(f.mul(x1[b ^ 3], x1[b]), x2[b ^ 1]);
    // consistency: D^2 proportional to H(S(null)), checked against the base index
    let x0 = hadamard(f, &squared(f, null));
    let d2 = squared(f, &d);
    if (0..4).any(|i| f.mul(d2[i], x0[b]) != f.mul(d2[b], x0[i])) {
        return None;
    }
    let zeros: Vec<usize> = (0..4).filter(|&i| f.is_zero(d[i])).collect();
    if zeros.len() > 1 {
        return None;
    }
    // projective 1/D: dinv_i = prod_{j != i} D_j (skipping the zero coordinate, whose dinv is 0)
    let w = d.map(|v| if f.is_zero(v) { f.one() } else { v });
    let (p01, p23) = (f.mul(w[0], w[1]), f.mul(w[2], w[3]));
    let mut dinv = [
        f.mul(w[1], p23),
        f.mul(w[0], p23),
        f.mul(p01, w[3]),
        f.mul(p01, w[2]),
    ];
    if let Some(&z) = zeros.first() {
        dinv[z] = f.zero();
    }
    Some(ThetaIso {
        dinv,
        zero: zeros.first().copied(),
        domain: *null,
        codomain: hadamard(f, &d),
    })
}

impl<E: Copy> ThetaIso<E> {
    /// Dual coordinates H(theta^B(f x)) (the zero index, if any, left at zero).
    fn dual<F: Field<E = E>>(&self, f: &F, x: &Th<E>) -> Th<E> {
        let y = hadamard(f, &squared(f, x));
        [0, 1, 2, 3].map(|i| f.mul(y[i], self.dinv[i]))
    }
    /// Image of a point (no vanishing dual coordinate).
    pub fn eval<F: Field<E = E>>(&self, f: &F, x: &Th<E>) -> Th<E> {
        debug_assert!(self.zero.is_none());
        hadamard(f, &self.dual(f, x))
    }
    /// Image of x through a gluing isogeny, given also the coordinates of x + T'_i where
    /// T'_i has order 4 and f(T'_i) = e_i / 2 (i = 1: `shift` = 2, i = 2: `shift` = 1).
    pub fn eval_glue<F: Field<E = E>>(
        &self,
        f: &F,
        x: &Th<E>,
        xt: &Th<E>,
        shift: usize,
    ) -> Option<Th<E>> {
        let Some(z) = self.zero else {
            return Some(self.eval(f, x));
        };
        let y = self.dual(f, x);
        let yt = self.dual(f, xt);
        // Y(f(x + T'))_chi = lambda Y(f x)_{chi ^ shift}: y[z] = yt[z ^ shift] y[c] / yt[c ^ shift],
        // written projectively (all coordinates times yt[c ^ shift])
        for c in 0..4usize {
            if c == z || c ^ shift == z || f.is_zero(y[c]) || f.is_zero(yt[c ^ shift]) {
                continue;
            }
            let s = yt[c ^ shift];
            let mut out = y.map(|v| f.mul(v, s));
            out[z] = f.mul(yt[z ^ shift], y[c]);
            return Some(hadamard(f, &out));
        }
        None
    }
}

// ------------------------------------------------------------------ invariants of the codomain

/// theta[a/2, b/2](0, Omega)^2 for the 16 characteristics, index b1 + 2 b2 + 4 a1 + 8 a2
/// (Dupont's numbering), from the level-2 null theta[c/2, 0](0, 2 Omega). Odd ones are zero.
pub fn fundamental_squares<F: Field>(f: &F, n: &Th<F::E>) -> [F::E; 16] {
    let mut out = [f.zero(); 16];
    for a in 0..4usize {
        for b in 0..4usize {
            let (a1, a2, b1, b2) = (a & 1, a >> 1, b & 1, b >> 1);
            let mut s = f.zero();
            for c1 in 0..2usize {
                for c2 in 0..2usize {
                    let ci = 2 * c1 + c2;
                    let cai = 2 * (c1 ^ a1) + (c2 ^ a2);
                    let t = f.mul(n[ci], n[cai]);
                    s = if (b1 * c1 + b2 * c2) % 2 == 1 {
                        f.sub(s, t)
                    } else {
                        f.add(s, t)
                    };
                }
            }
            out[b1 + 2 * b2 + 4 * a1 + 8 * a2] = s;
        }
    }
    out
}

pub const EVEN: [usize; 10] = [0, 1, 2, 3, 4, 6, 8, 9, 12, 15];

/// Rosenhain invariants (lambda, mu, nu) of y^2 = x (x - 1)(x - lambda)(x - mu)(x - nu):
/// lambda = t0 t2 / (t1 t3), mu = t2 t12 / (t1 t15), nu = t0 t12 / (t3 t15), t_i = theta_i^2
/// in the numbering of `fundamental_squares`. (Cited from memory with t1 and t3 exchanged in
/// mu and nu; the version here is the one that matches the Howe–Leprévost–Poonen gluing in
/// `tests/theta.rs`.)
pub fn rosenhain<F: Field>(f: &F, n: &Th<F::E>) -> Option<[F::E; 3]> {
    let t = fundamental_squares(f, n);
    if [1, 3, 15].iter().any(|&i| f.is_zero(t[i])) {
        return None;
    }
    let lam = f.div(f.mul(t[0], t[2]), f.mul(t[1], t[3]));
    let mu = f.div(f.mul(t[2], t[12]), f.mul(t[1], t[15]));
    let nu = f.div(f.mul(t[0], t[12]), f.mul(t[3], t[15]));
    Some([lam, mu, nu])
}

/// Parity a.b of a characteristic in Dupont's numbering (b1 + 2 b2 + 4 a1 + 8 a2).
fn parity(i: usize) -> usize {
    ((i >> 2) & i & 1) ^ ((i >> 3) & (i >> 1) & 1)
}

/// If the null point is that of a product of elliptic curves (exactly one even theta constant
/// vanishes, at z0), the j-invariants of the two factors. The nine other even characteristics
/// form a 3 x 3 grid (theta^2 = u_r v_c up to roots of unity: factor 1's constants times factor
/// 2's): m and n share a row or a column iff m + n + z0 is odd (azygetic; invariant under the
/// symplectic group, and true for the product structure, where z0 = 15).
pub fn split_j<F: Field>(f: &F, n: &Th<F::E>) -> Option<(F::E, F::E)> {
    let t = fundamental_squares(f, n);
    let zeros: Vec<usize> = EVEN.iter().copied().filter(|&i| f.is_zero(t[i])).collect();
    if zeros.len() != 1 {
        return None;
    }
    let z0 = zeros[0];
    let ch: Vec<usize> = EVEN.iter().copied().filter(|&i| i != z0).collect();
    let adj = |m: usize, k: usize| m != k && parity(m ^ k ^ z0) == 1;
    // the two triangles (row and column) through ch[0]
    let nb: Vec<usize> = ch.iter().copied().filter(|&k| adj(ch[0], k)).collect();
    if nb.len() != 4 {
        return None;
    }
    let mate = nb.iter().copied().find(|&k| adj(nb[0], k))?;
    let row0 = [ch[0], nb[0], mate];
    let col0: Vec<usize> = std::iter::once(ch[0])
        .chain(nb.iter().copied().filter(|&k| k != nb[0] && k != mate))
        .collect();
    if col0.len() != 3 {
        return None;
    }
    let mut grid = [[0usize; 3]; 3];
    grid[0] = row0;
    for r in 1..3 {
        let c0 = col0[r];
        grid[r][0] = c0;
        // the rest of c0's row: its neighbours outside column 0, placed under the row-0 element
        // of their column (each shares a column with exactly one of row0[1], row0[2])
        for o in ch
            .iter()
            .copied()
            .filter(|&k| adj(c0, k) && !col0.contains(&k))
        {
            let c = (1..3).find(|&c| adj(row0[c], o))?;
            grid[r][c] = o;
        }
    }
    // theta^4 = t^2 up to sign, theta^8 exact: check rank one on theta^8
    let v4 = |i: usize| f.sq(t[i]);
    let v8 = |i: usize| f.sq(v4(i));
    for r in 0..3 {
        for c in 0..3 {
            if f.mul(v8(grid[0][0]), v8(grid[r][c])) != f.mul(v8(grid[0][c]), v8(grid[r][0])) {
                return None;
            }
        }
    }
    let col = [v4(grid[0][0]), v4(grid[1][0]), v4(grid[2][0])];
    let row = [v4(grid[0][0]), v4(grid[0][1]), v4(grid[0][2])];
    Some((
        j_from_jacobi_triple(f, &col)?,
        j_from_jacobi_triple(f, &row)?,
    ))
}

/// j from (theta_3^4, theta_4^4, theta_2^4) of one elliptic factor, in unknown order, up to a
/// common factor and up to individual signs: Jacobi's theta_3^4 = theta_2^4 + theta_4^4 with
/// some choice of signs identifies a Legendre lambda (any solution gives the same j).
fn j_from_jacobi_triple<F: Field>(f: &F, s: &[F::E; 3]) -> Option<F::E> {
    for h in 0..3 {
        let (i, k) = ((h + 1) % 3, (h + 2) % 3);
        if f.is_zero(s[h]) {
            continue;
        }
        for (ei, ek) in [(false, false), (false, true), (true, false), (true, true)] {
            let si = if ei { f.neg(s[i]) } else { s[i] };
            let sk = if ek { f.neg(s[k]) } else { s[k] };
            if s[h] == f.add(si, sk) {
                return Some(j_from_lambda(f, f.div(si, s[h])));
            }
        }
    }
    None
}

// ------------------------------------------------------------------ chains (Kani)

#[derive(Clone, Debug)]
pub struct ChainResult<E> {
    /// theta null of every codomain (the last one first split-tested)
    pub nulls: Vec<Th<E>>,
    /// j-invariants of the factors if the final codomain is a product
    pub split: Option<(E, E)>,
    /// images of the extra points under the whole chain
    pub images: Vec<Th<E>>,
}

/// Projective inverse (v1 v2 v3, v0 v2 v3, v0 v1 v3, v0 v1 v2); None if a coordinate is zero.
fn proj_inv<F: Field>(f: &F, v: &Th<F::E>) -> Option<Th<F::E>> {
    if v.iter().any(|&c| f.is_zero(c)) {
        return None;
    }
    let (p01, p23) = (f.mul(v[0], v[1]), f.mul(v[2], v[3]));
    Some([
        f.mul(v[1], p23),
        f.mul(v[0], p23),
        f.mul(p01, v[3]),
        f.mul(p01, v[2]),
    ])
}

/// Doubling on the Kummer surface of A in level-2 theta coordinates: with X = H(S(x)),
/// theta(2x) = H(S(X) / H(S(null))) / null. (Composite of the isogeny A -> B, z -> 2z, whose
/// dual coordinates are X / D with D^2 = H(S(null)), and B -> A, z -> z, given by
/// theta^A(z) * theta^A(0) = H(S(H(theta^B(z)))); the signs of D cancel in S.) Needs a Jacobian
/// (no zero coordinate in null or H(S(null))).
#[derive(Clone, Debug)]
pub struct ThetaDoubler<E> {
    inv_x0: Th<E>,
    inv_null: Th<E>,
}

impl<E: Copy> ThetaDoubler<E> {
    pub fn new<F: Field<E = E>>(f: &F, null: &Th<E>) -> Option<Self> {
        Some(ThetaDoubler {
            inv_x0: proj_inv(f, &hadamard(f, &squared(f, null)))?,
            inv_null: proj_inv(f, null)?,
        })
    }
    pub fn double<F: Field<E = E>>(&self, f: &F, x: &Th<E>) -> Th<E> {
        let xx = squared(f, &hadamard(f, &squared(f, x)));
        let y = [0, 1, 2, 3].map(|i| f.mul(xx[i], self.inv_x0[i]));
        let z = hadamard(f, &y);
        [0, 1, 2, 3].map(|i| f.mul(z[i], self.inv_null[i]))
    }
}

/// The (2^n, 2^n)-isogeny of E1 x E2 with kernel generated by [4] K1, [4] K2, where
/// K_i = (K_i^(1), K_i^(2)) have order 2^(n+2) and 2^(n+1) K_1, 2^(n+1) K_2 are a
/// gluing kernel {(T, psi(T))} (i.e. not contained in E1 x 0 or 0 x E2). The chain is a gluing
/// followed by n - 1 (2,2)-isogenies between Jacobians (or a final split). `extra` are points of
/// E1 x E2 to push through. None if some step's data is inconsistent.
///
/// After the gluing only K1, K2 (and `extra`) are pushed; the Jacobian part finds its 8-torsion
/// points by an optimal strategy over theta doublings and evaluations (O(n log n) operations).
pub fn chain<F: Field>(
    f: &F,
    e1: &Curve<F::E>,
    e2: &Curve<F::E>,
    k: [(Pt<F::E>, Pt<F::E>); 2],
    n: u32,
    extra: &[(Pt<F::E>, Pt<F::E>)],
) -> Option<ChainResult<F::E>> {
    chain_opts(f, e1, e2, k, n, extra, true)
}

/// `strategy = false`: push all the multiples [2^j] K_i through every step instead (O(n^2)
/// evaluations; kept for comparison).
pub fn chain_opts<F: Field>(
    f: &F,
    e1: &Curve<F::E>,
    e2: &Curve<F::E>,
    k: [(Pt<F::E>, Pt<F::E>); 2],
    n: u32,
    extra: &[(Pt<F::E>, Pt<F::E>)],
    strategy: bool,
) -> Option<ChainResult<F::E>> {
    assert!(n >= 1);
    let dbl1 = |p: &(Pt<F::E>, Pt<F::E>)| (padd(f, e1, &p.0, &p.0), padd(f, e2, &p.1, &p.1));
    let add = |p: &(Pt<F::E>, Pt<F::E>), q: &(Pt<F::E>, Pt<F::E>)| {
        (padd(f, e1, &p.0, &q.0), padd(f, e2, &p.1, &q.1))
    };
    // multiples [2^j] K_i for j = 0..=n, incrementally
    let mults: Vec<Vec<(Pt<F::E>, Pt<F::E>)>> = (0..2)
        .map(|i| {
            let mut v = vec![k[i]];
            for j in 0..n as usize {
                let nx = dbl1(&v[j]);
                v.push(nx);
            }
            v
        })
        .collect();
    // 4-torsion above the gluing kernel and the one-dimensional structures
    let t4 = [mults[0][n as usize], mults[1][n as usize]];
    let th = ProductTheta {
        t1: theta1_structure(f, e1, &t4[0].0, &t4[1].0)?,
        t2: theta1_structure(f, e2, &t4[0].1, &t4[1].1)?,
    };
    let null = th.null(f);
    let t8 = [mults[0][n as usize - 1], mults[1][n as usize - 1]];
    let glue = theta_isogeny(
        f,
        &null,
        &th.point(f, &t8[0].0, &t8[0].1),
        &th.point(f, &t8[1].0, &t8[1].1),
    )?;
    let push_glue = |p: &(Pt<F::E>, Pt<F::E>)| -> Option<Th<F::E>> {
        let x = th.point(f, &p.0, &p.1);
        for (ti, shift) in [(0usize, 2usize), (1, 1)] {
            let q = add(p, &t4[ti]);
            let xt = th.point(f, &q.0, &q.1);
            if let Some(y) = glue.eval_glue(f, &x, &xt, shift) {
                return Some(y);
            }
        }
        None
    };
    let mut images: Vec<Th<F::E>> = Vec::new();
    for p in extra {
        images.push(push_glue(p)?);
    }
    let mut nulls = vec![glue.codomain];
    if n == 1 {
        let split = split_j(f, &glue.codomain);
        return Some(ChainResult {
            nulls,
            split,
            images,
        });
    }
    let mut cur = glue.codomain;
    if strategy {
        // generators of order 2^(n+1) on the first Jacobian; m = n - 1 steps remain
        let g = [push_glue(&mults[0][0])?, push_glue(&mults[1][0])?];
        let m = (n - 1) as usize;
        let splits = crate::kernel::two_power::optimal_splits(m, 1.0, 0.6);
        fn rec<F: Field>(
            f: &F,
            cur: &mut Th<F::E>,
            g: [Th<F::E>; 2],
            m: usize,
            stack: &mut Vec<Th<F::E>>,
            splits: &[usize],
            nulls: &mut Vec<Th<F::E>>,
        ) -> Option<()> {
            if m == 1 {
                let iso = theta_isogeny(f, cur, &g[0], &g[1])?;
                if iso.zero.is_some() {
                    return None;
                }
                for x in stack.iter_mut() {
                    *x = iso.eval(f, x);
                }
                *cur = iso.codomain;
                nulls.push(*cur);
                return Some(());
            }
            let i = splits[m];
            let dbl = ThetaDoubler::new(f, cur)?;
            let mut h = g;
            for _ in 0..(m - i) {
                h = [dbl.double(f, &h[0]), dbl.double(f, &h[1])];
            }
            stack.push(g[0]);
            stack.push(g[1]);
            rec(f, cur, h, i, stack, splits, nulls)?;
            let g1 = stack.pop().unwrap();
            let g0 = stack.pop().unwrap();
            rec(f, cur, [g0, g1], m - i, stack, splits, nulls)
        }
        rec(f, &mut cur, g, m, &mut images, &splits, &mut nulls)?;
    } else {
        // after the gluing, [2^j] K_i is needed for j <= n - 2
        let mut gens: Vec<Vec<Th<F::E>>> = Vec::new();
        for i in 0..2 {
            let mut v = Vec::new();
            for j in 0..n as usize - 1 {
                v.push(push_glue(&mults[i][j])?);
            }
            gens.push(v);
        }
        for step in 2..=n {
            // kernel 2^(n - step + 1) G_i; T''_i = 2^(n - step) G_i = gens[i][n - step]
            let jj = (n - step) as usize;
            let iso = theta_isogeny(f, &cur, &gens[0][jj], &gens[1][jj])?;
            if iso.zero.is_some() {
                return None;
            }
            for g in gens.iter_mut() {
                g.truncate(jj);
                for x in g.iter_mut() {
                    *x = iso.eval(f, x);
                }
            }
            for x in images.iter_mut() {
                *x = iso.eval(f, x);
            }
            cur = iso.codomain;
            nulls.push(cur);
        }
    }
    let split = split_j(f, &cur);
    Some(ChainResult {
        nulls,
        split,
        images,
    })
}
