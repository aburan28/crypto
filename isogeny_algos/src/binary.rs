//! Ordinary elliptic curves over GF(2^n): E: y^2 + xy = x^3 + a2 x^2 + a4 x + a6 (a1 = 1, a3 = 0),
//! j = 1/(a6 + a4^2). Group law, point counting, Vélu and Kohel in characteristic 2, division
//! polynomials, kernel polynomials of odd prime degree, and the j-line neighbour oracles used by
//! the path-finding algorithms (Galbraith, GHS) on binary curves.
//!
//! Vélu in characteristic 2 (odd kernel G, S = one point of each pair +-Q): with a1 = 1, a3 = 0,
//! g^y_Q = x_Q, t_Q = x_Q, u_Q = x_Q^2, so w = sum(u_Q + x_Q t_Q) = 0 and the codomain is
//!   E': (a2, a4 + t, a6 + t),  t = sum_{Q in S} x_Q  (= the x^{d-1} coefficient of h),
//! with x-map  x + sum x_Q x/(x + x_Q)^2 = x (h^2 + x h'^2 + h h')/h^2  for h = prod (x + x_Q)
//! (char 2: h'' = 0, (h'/h)' = h'^2/h^2).
use crate::curve::Pt;
use crate::field::{Field, Rng};
use crate::gf2n::GF2n;
use crate::poly::{self, Poly};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BinCurve {
    pub a2: u64,
    pub a4: u64,
    pub a6: u64,
}

pub type BPt = Pt<u64>;

impl BinCurve {
    pub fn new(a2: u64, a4: u64, a6: u64) -> Self {
        BinCurve { a2, a4, a6 }
    }
    /// a6 + a4^2 (non-zero for a non-singular curve).
    pub fn disc(&self, f: &GF2n) -> u64 {
        self.a6 ^ f.sq(self.a4)
    }
    pub fn j(&self, f: &GF2n) -> u64 {
        f.inv(self.disc(f))
    }
    /// The curve with a2 = a4 = 0 and the given j != 0.
    pub fn from_j(f: &GF2n, j: u64) -> Self {
        BinCurve::new(0, 0, f.inv(j))
    }
    /// Isomorphic model with a4 = 0 via y -> y + a4 (x unchanged).
    pub fn normalized(&self, f: &GF2n) -> Self {
        BinCurve::new(self.a2, 0, self.disc(f))
    }
    pub fn rhs(&self, f: &GF2n, x: u64) -> u64 {
        // x^3 + a2 x^2 + a4 x + a6
        let x2 = f.sq(x);
        f.mul(x2, x) ^ f.mul(self.a2, x2) ^ f.mul(self.a4, x) ^ self.a6
    }
    pub fn on_curve(&self, f: &GF2n, p: &BPt) -> bool {
        match *p {
            Pt::Inf => true,
            Pt::Aff(x, y) => f.sq(y) ^ f.mul(x, y) == self.rhs(f, x),
        }
    }
    /// Is x the abscissa of an F_q-rational point?
    pub fn has_x(&self, f: &GF2n, x: u64) -> bool {
        if x == 0 {
            return true;
        }
        f.trace(f.mul(self.rhs(f, x), f.inv(f.sq(x)))) == 0
    }
    pub fn neg(&self, p: &BPt) -> BPt {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => Pt::Aff(x, y ^ x),
        }
    }
    pub fn add(&self, f: &GF2n, p: &BPt, q: &BPt) -> BPt {
        let (x1, y1, x2, y2) = match (*p, *q) {
            (Pt::Inf, _) => return *q,
            (_, Pt::Inf) => return *p,
            (Pt::Aff(a, b), Pt::Aff(c, d)) => (a, b, c, d),
        };
        let lam = if x1 != x2 {
            f.div(y1 ^ y2, x1 ^ x2)
        } else if y2 == y1 ^ x1 {
            // q = -p (covers doubling the 2-torsion point x = 0)
            return Pt::Inf;
        } else {
            // doubling: (3x^2 + 2 a2 x + a4 - y)/(2y + x) = x + (y + a4)/x
            x1 ^ f.div(y1 ^ self.a4, x1)
        };
        let x3 = f.sq(lam) ^ lam ^ self.a2 ^ x1 ^ x2;
        let y3 = f.mul(lam ^ 1, x3) ^ y1 ^ f.mul(lam, x1);
        Pt::Aff(x3, y3)
    }
    pub fn mul(&self, f: &GF2n, p: &BPt, mut k: u128) -> BPt {
        let mut r = Pt::Inf;
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
    pub fn random_point(&self, f: &GF2n, rng: &mut Rng) -> BPt {
        loop {
            let x = f.random(rng);
            if x == 0 {
                continue;
            }
            // y = x z, z^2 + z = rhs/x^2
            let c = f.mul(self.rhs(f, x), f.inv(f.sq(x)));
            if let Some(z) = f.solve_quadratic(c) {
                let z = if rng.next() & 1 == 1 { z ^ 1 } else { z };
                return Pt::Aff(x, f.mul(x, z));
            }
        }
    }

    /// #E(F_q) by the trace criterion, O(q): for testing.
    pub fn order_naive(&self, f: &GF2n) -> u128 {
        let q = f.size();
        let mut n = 2u128; // infinity and (0, sqrt(a6 + ...))
        for x in 1..q as u64 {
            if self.has_x(f, x) {
                n += 2;
            }
        }
        n
    }

    /// #E(F_q) by baby-step giant-step in the Hasse interval: for random points, a multiple of
    /// the order found by BSGS is reduced to the exact order (factoring), and the lcm L of the
    /// orders narrows the candidates to the even multiples of L in the interval.
    pub fn order(&self, f: &GF2n, rng: &mut Rng) -> u128 {
        let q = f.size();
        if q < (1 << 14) {
            return self.order_naive(f);
        }
        let w = 2 * ((q as f64).sqrt().ceil() as u128) + 2;
        let (lo, hi) = (q + 1 - w, q + 1 + w);
        let m = (((hi - lo + 1) as f64).sqrt().ceil() as u128).max(1);
        let mut l = 1u128;
        for _ in 0..64 {
            let p = self.random_point(f, rng);
            let mut baby: HashMap<BPt, u128> = HashMap::new();
            let mut r = Pt::Inf;
            for j in 0..m {
                baby.entry(r).or_insert(j);
                r = self.add(f, &r, &p);
            }
            let step = self.mul(f, &p, m);
            let mut g = self.mul(f, &p, lo);
            let mut k = None;
            let mut i = 0u128;
            while lo + i * m <= hi + m {
                if let Some(&j) = baby.get(&self.neg(&g)) {
                    k = Some(lo + i * m + j);
                    break;
                }
                g = self.add(f, &g, &step);
                i += 1;
            }
            let k = k.expect("BSGS: no multiple of the point order in the Hasse interval");
            let mut ord = k as u64;
            for (r, _) in crate::field::factor_u64(k as u64) {
                while ord.is_multiple_of(r) && self.mul(f, &p, (ord / r) as u128) == Pt::Inf {
                    ord /= r;
                }
            }
            l = lcm(l, ord as u128);
            let first = lo.div_ceil(l) * l;
            let cands: Vec<u128> = (0..)
                .map(|t| first + t * l)
                .take_while(|&c| c <= hi)
                .filter(|c| c % 2 == 0)
                .collect();
            if cands.len() == 1 {
                return cands[0];
            }
        }
        panic!("order: group exponent too small to isolate #E in the Hasse interval");
    }
}

fn lcm(a: u128, b: u128) -> u128 {
    let (mut x, mut y) = (a, b);
    while y != 0 {
        let t = x % y;
        x = y;
        y = t;
    }
    a / x * b
}

/// The x-map of the normalised isogeny with kernel polynomial h: (num, den) with
/// phi_x = num/den = x (h^2 + x h'^2 + h h') / h^2.
pub fn x_map(f: &GF2n, h: &Poly<GF2n>) -> (Poly<GF2n>, Poly<GF2n>) {
    let hp = poly::derivative(f, h);
    let h2 = poly::mul(f, h, h);
    let x = poly::x_poly(f);
    let t = poly::add(
        f,
        &poly::add(f, &h2, &poly::mul(f, &x, &poly::mul(f, &hp, &hp))),
        &poly::mul(f, h, &hp),
    );
    (poly::mul(f, &x, &t), h2)
}

/// Codomain of the isogeny with kernel polynomial h of degree d = (l-1)/2 (Kohel in char 2):
/// (a2, a4 + t, a6 + t) with t the sum of the roots of h.
pub fn kohel_codomain(_f: &GF2n, e: &BinCurve, h: &Poly<GF2n>) -> BinCurve {
    let d = h.len() - 1;
    let t = h[d - 1]; // monic h: x^d + h_{d-1} x^{d-1} + ..., sum of roots = h_{d-1} in char 2
    BinCurve::new(e.a2, e.a4 ^ t, e.a6 ^ t)
}

/// Full Vélu isogeny from the points of S (one of each +-Q in the kernel of odd order).
pub struct BinVelu {
    pub dom: BinCurve,
    pub cod: BinCurve,
    pub s: Vec<(u64, u64)>,
}

pub fn velu(f: &GF2n, e: &BinCurve, gen: &BPt, ell: u64) -> BinVelu {
    let d = (ell - 1) / 2;
    let mut s = vec![];
    let mut cur = *gen;
    for _ in 0..d {
        match cur {
            Pt::Aff(x, y) => s.push((x, y)),
            Pt::Inf => panic!("generator order too small"),
        }
        cur = e.add(f, &cur, gen);
    }
    let t = s.iter().fold(0u64, |acc, q| acc ^ q.0);
    BinVelu {
        dom: *e,
        cod: BinCurve::new(e.a2, e.a4 ^ t, e.a6 ^ t),
        s,
    }
}

impl BinVelu {
    pub fn eval(&self, f: &GF2n, p: &BPt) -> BPt {
        let (x, y) = match *p {
            Pt::Inf => return Pt::Inf,
            Pt::Aff(x, y) => (x, y),
        };
        if self.s.iter().any(|q| q.0 == x) {
            return Pt::Inf;
        }
        let a4 = self.dom.a4;
        let (mut xs, mut ys) = (x, y);
        for &(xq, yq) in &self.s {
            let dx = x ^ xq;
            let i1 = f.inv(dx);
            let i2 = f.sq(i1);
            let i3 = f.mul(i2, i1);
            // X += x_Q x / (x + x_Q)^2
            xs ^= f.mul(f.mul(xq, x), i2);
            // Y += x_Q^2 x/(dx)^3 + x_Q (dx + y + y_Q)/(dx)^2 + (x_Q^2 + (x_Q^2 + a4 + y_Q) x_Q)/(dx)^2
            let gx = f.sq(xq) ^ a4 ^ yq;
            let term = f.mul(f.mul(f.sq(xq), x), i3)
                ^ f.mul(f.mul(xq, dx ^ y ^ yq), i2)
                ^ f.mul(f.sq(xq) ^ f.mul(gx, xq), i2);
            ys ^= term;
        }
        Pt::Aff(xs, ys)
    }
}

/// Division polynomials f_n of the normalised curve y^2 + xy = x^3 + a2 x^2 + a6 (a4 = 0):
/// f_0 = 0, f_1 = 1, f_2 = x, f_3 = x^4 + x^3 + a6, f_4 = x^6 + a6 x^2,
/// f_{2m+1} = f_{m+2} f_m^3 + f_{m-1} f_{m+1}^3,  f_{2m} = (f_{m+2} f_{m-1}^2 + f_{m-2} f_{m+1}^2) f_m / x.
/// For odd n the roots of f_n are the abscissas of the non-zero n-torsion points.
pub fn division_poly(f: &GF2n, e: &BinCurve, n: usize) -> Poly<GF2n> {
    assert!(e.a4 == 0, "normalise first");
    let mut memo: HashMap<usize, Poly<GF2n>> = HashMap::new();
    fn get(f: &GF2n, a6: u64, n: usize, memo: &mut HashMap<usize, Poly<GF2n>>) -> Poly<GF2n> {
        if let Some(p) = memo.get(&n) {
            return p.clone();
        }
        let r: Poly<GF2n> = match n {
            0 => vec![],
            1 => vec![1],
            2 => vec![0, 1],
            3 => vec![a6, 0, 0, 1, 1],
            4 => vec![0, 0, a6, 0, 0, 0, 1],
            _ => {
                let m = n / 2;
                if n % 2 == 1 {
                    let (a, b, c, d) = (
                        get(f, a6, m + 2, memo),
                        get(f, a6, m, memo),
                        get(f, a6, m - 1, memo),
                        get(f, a6, m + 1, memo),
                    );
                    let b3 = poly::mul(f, &poly::mul(f, &b, &b), &b);
                    let d3 = poly::mul(f, &poly::mul(f, &d, &d), &d);
                    poly::add(f, &poly::mul(f, &a, &b3), &poly::mul(f, &c, &d3))
                } else {
                    let (a, c, g, d, b) = (
                        get(f, a6, m + 2, memo),
                        get(f, a6, m - 1, memo),
                        get(f, a6, m - 2, memo),
                        get(f, a6, m + 1, memo),
                        get(f, a6, m, memo),
                    );
                    let t = poly::add(
                        f,
                        &poly::mul(f, &a, &poly::mul(f, &c, &c)),
                        &poly::mul(f, &g, &poly::mul(f, &d, &d)),
                    );
                    let num = poly::mul(f, &t, &b);
                    let (q, r) = poly::divrem(f, &num, &poly::x_poly(f));
                    debug_assert!(r.is_empty());
                    q
                }
            }
        };
        memo.insert(n, r.clone());
        r
    }
    get(f, e.a6, n, &mut memo)
}

/// x-only doubling on a normalised curve: x(2P) = x^2 + a6/x^2.
fn xdbl_norm(f: &GF2n, a6: u64, x: u64) -> Option<u64> {
    if x == 0 {
        return None;
    }
    let x2 = f.sq(x);
    Some(x2 ^ f.mul(a6, f.inv(x2)))
}

/// Does h (monic, degree (l-1)/2, dividing f_l) define a subgroup? Checked by the x-map commuting
/// with doubling on a few random points: phi_x(x(2P)) = x(2 phi(P)) on the normalised codomain.
pub fn is_kernel(f: &GF2n, e: &BinCurve, h: &Poly<GF2n>, rng: &mut Rng) -> bool {
    let en = e.normalized(f);
    let cod = kohel_codomain(f, &en, h).normalized(f);
    let (num, den) = x_map(f, h);
    let phi = |x: u64| -> Option<u64> {
        let dv = poly::eval(f, &den, x);
        if dv == 0 {
            None
        } else {
            Some(f.div(poly::eval(f, &num, x), dv))
        }
    };
    let mut checked = 0;
    for _ in 0..20 {
        let p = en.random_point(f, rng);
        let Pt::Aff(x, _) = p else { continue };
        let (Some(x2), Some(px)) = (xdbl_norm(f, en.a6, x), phi(x)) else {
            continue;
        };
        let (Some(lhs), Some(rhs)) = (phi(x2), xdbl_norm(f, cod.a6, px)) else {
            continue;
        };
        if lhs != rhs {
            return false;
        }
        checked += 1;
        if checked == 4 {
            return true;
        }
    }
    checked > 0
}

/// All F_q-rational kernel polynomials of degree-l isogenies from e (l odd prime).
pub fn kernel_polys(f: &GF2n, e: &BinCurve, ell: u64, rng: &mut Rng) -> Vec<Poly<GF2n>> {
    let en = e.normalized(f);
    let psi = poly::monic(f, &division_poly(f, &en, ell as usize));
    let facs = poly::factor_squarefree(f, &psi, rng);
    let d = ((ell - 1) / 2) as usize;
    let degs: Vec<usize> = facs.iter().map(|p| p.len() - 1).collect();
    let mut subsets = vec![];
    fn rec(i: usize, left: usize, degs: &[usize], cur: &mut Vec<usize>, res: &mut Vec<Vec<usize>>) {
        if left == 0 {
            res.push(cur.clone());
            return;
        }
        if i == degs.len() {
            return;
        }
        if degs[i] <= left {
            cur.push(i);
            rec(i + 1, left - degs[i], degs, cur, res);
            cur.pop();
        }
        rec(i + 1, left, degs, cur, res);
    }
    rec(0, d, &degs, &mut vec![], &mut subsets);
    let mut out = vec![];
    for s in subsets {
        let mut h = vec![1u64];
        for &i in &s {
            h = poly::mul(f, &h, &facs[i]);
        }
        if is_kernel(f, e, &h, rng) {
            out.push(h);
        }
    }
    out
}

/// Kernel-polynomial neighbour oracle on the j-line of ordinary binary curves: for each l, the
/// j-invariants of the codomains of the rational l-isogenies (division polynomial factoring).
pub struct BinKernelOracle<'a> {
    pub f: &'a GF2n,
    pub ells: Vec<u64>,
}

impl crate::path::graph::Oracle<u64> for BinKernelOracle<'_> {
    fn neighbors(&self, j: u64, rng: &mut Rng) -> Vec<(usize, u64)> {
        let e = BinCurve::from_j(self.f, j);
        let mut out = vec![];
        for &l in &self.ells {
            for h in kernel_polys(self.f, &e, l, rng) {
                let jn = kohel_codomain(self.f, &e, &h).j(self.f);
                if !out.contains(&(l as usize, jn)) {
                    out.push((l as usize, jn));
                }
            }
        }
        out
    }
}

/// Phi_l modulo 2 (integer coefficients by CRT over 62-bit primes), over GF(2^n).
pub fn phi_mod2(f: &GF2n, ell: usize) -> crate::find::modpoly::Phi<GF2n> {
    crate::find::modpoly::Phi::via_crt(f, ell)
}
