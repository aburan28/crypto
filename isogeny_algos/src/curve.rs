//! Short Weierstrass curves y^2 = x^3 + a x + b (char > 3), affine arithmetic,
//! point counting over F_p, and the isogeny interface shared by all algorithms.
use crate::field::{Field, Rng, Zp};
use crate::poly::{self, Poly};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Curve<E> {
    pub a: E,
    pub b: E,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Pt<E> {
    Inf,
    Aff(E, E),
}

impl<E: Copy> Curve<E> {
    pub fn new(a: E, b: E) -> Self {
        Curve { a, b }
    }
}

pub fn rhs<F: Field>(f: &F, c: &Curve<F::E>, x: F::E) -> F::E {
    f.add(f.mul(f.add(f.mul(x, x), c.a), x), c.b)
}
pub fn on_curve<F: Field>(f: &F, c: &Curve<F::E>, p: &Pt<F::E>) -> bool {
    match *p {
        Pt::Inf => true,
        Pt::Aff(x, y) => f.mul(y, y) == rhs(f, c, x),
    }
}
pub fn is_smooth<F: Field>(f: &F, c: &Curve<F::E>) -> bool {
    let d = f.add(
        f.mul(f.from_u64(4), f.mul(c.a, f.mul(c.a, c.a))),
        f.mul(f.from_u64(27), f.mul(c.b, c.b)),
    );
    !f.is_zero(d)
}
/// j = 1728 * 4a^3 / (4a^3 + 27 b^2)
pub fn jinv<F: Field>(f: &F, c: &Curve<F::E>) -> F::E {
    let a3 = f.mul(f.from_u64(4), f.mul(c.a, f.mul(c.a, c.a)));
    let d = f.add(a3, f.mul(f.from_u64(27), f.mul(c.b, c.b)));
    f.mul(f.from_u64(1728), f.div(a3, d))
}
/// A curve with the given j (j != 0, 1728).
pub fn from_j<F: Field>(f: &F, j: F::E) -> Curve<F::E> {
    let k = f.sub(f.from_u64(1728), j);
    let a = f.mul(f.from_u64(3), f.mul(j, k));
    let b = f.mul(f.from_u64(2), f.mul(j, f.mul(k, k)));
    Curve { a, b }
}
pub fn neg<F: Field>(f: &F, p: &Pt<F::E>) -> Pt<F::E> {
    match *p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(x, f.neg(y)),
    }
}
pub fn padd<F: Field>(f: &F, c: &Curve<F::E>, p: &Pt<F::E>, q: &Pt<F::E>) -> Pt<F::E> {
    match (*p, *q) {
        (Pt::Inf, r) | (r, Pt::Inf) => r,
        (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
            let lam = if x1 == x2 {
                if f.is_zero(f.add(y1, y2)) {
                    return Pt::Inf;
                }
                let num = f.add(f.mul(f.from_u64(3), f.mul(x1, x1)), c.a);
                f.div(num, f.add(y1, y1))
            } else {
                f.div(f.sub(y2, y1), f.sub(x2, x1))
            };
            let x3 = f.sub(f.sub(f.mul(lam, lam), x1), x2);
            let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
            Pt::Aff(x3, y3)
        }
    }
}
pub fn pmul<F: Field>(f: &F, c: &Curve<F::E>, p: &Pt<F::E>, mut k: u128) -> Pt<F::E> {
    let mut r = Pt::Inf;
    let mut b = *p;
    while k > 0 {
        if k & 1 == 1 {
            r = padd(f, c, &r, &b);
        }
        b = padd(f, c, &b, &b);
        k >>= 1;
    }
    r
}

/// [k]P for a big scalar.
pub fn pmul_big<F: Field>(
    f: &F,
    c: &Curve<F::E>,
    p: &Pt<F::E>,
    k: &crate::bigint::Big,
) -> Pt<F::E> {
    let mut r = Pt::Inf;
    for i in (0..k.bits()).rev() {
        r = padd(f, c, &r, &r);
        if k.bit(i) {
            r = padd(f, c, &r, p);
        }
    }
    r
}

/// Random affine point over any field of odd size.
pub fn random_point_f<F: Field>(f: &F, c: &Curve<F::E>, rng: &mut Rng) -> Pt<F::E> {
    loop {
        let x = f.random(rng);
        if let Some(y) = f.sqrt(rhs(f, c, x)) {
            return Pt::Aff(x, y);
        }
    }
}

impl Curve<u64> {
    pub fn random_point(&self, fp: &Zp, rng: &mut Rng) -> Pt<u64> {
        loop {
            let x = fp.random(rng);
            if let Some(y) = fp.sqrt(rhs(fp, self, x)) {
                return Pt::Aff(x, y);
            }
        }
    }
}

/// #E(F_p): table of squares for small p, Mestre-style BSGS in the Hasse interval otherwise.
pub fn order(fp: &Zp, c: &Curve<u64>, rng: &mut Rng) -> u64 {
    let p = fp.p;
    if p < (1 << 22) {
        let mut sq = vec![false; p as usize];
        for i in 0..p {
            sq[((i as u128 * i as u128) % p as u128) as usize] = true;
        }
        let mut n = 1i64 + p as i64;
        for x in 0..p {
            let r = rhs(fp, c, x);
            if r == 0 {
            } else if sq[r as usize] {
                n += 1;
            } else {
                n -= 1;
            }
        }
        return n as u64;
    }
    let w = 2 * ((p as f64).sqrt() as u64) + 3;
    let mid = p + 1;
    let m = ((2.0 * w as f64).sqrt() as u64) + 1;
    let mut cands: Option<Vec<u64>> = None;
    loop {
        let pt = c.random_point(fp, rng);
        let mut baby: HashMap<u64, Vec<(u64, bool)>> = HashMap::new();
        let mut cur = Pt::Inf;
        for j in 0..=m {
            if let Pt::Aff(x, y) = cur {
                baby.entry(x).or_default().push((j, y == 0 || y < fp.p / 2));
            } else {
                baby.entry(u64::MAX).or_default().push((0, true));
            }
            cur = padd(fp, c, &cur, &pt);
        }
        let step = pmul(fp, c, &pt, (2 * m + 1) as u128);
        let lo = mid - w;
        let mut r = pmul(fp, c, &pt, lo as u128);
        let mut found = vec![];
        let reps = (2 * w) / (2 * m + 1) + 2;
        for a in 0..reps {
            let base = lo + a * (2 * m + 1);
            // need base + b' = N with R_a + b' P = O, b' in [0, 2m]; use b' = m + b, b in [-m, m]
            // equivalently (base + m) P = -b P.
            let rm = padd(fp, c, &r, &pmul(fp, c, &pt, m as u128));
            let key = match rm {
                Pt::Inf => u64::MAX,
                Pt::Aff(x, _) => x,
            };
            if let Some(v) = baby.get(&key) {
                for &(j, _) in v {
                    // rm = +/- jP
                    let jp = pmul(fp, c, &pt, j as u128);
                    if rm == neg(fp, &jp) {
                        // (base+m) P = -jP  => N = base + m + j
                        found.push(base + m + j);
                    }
                    if rm == jp {
                        found.push(base + m - j);
                    }
                }
            }
            r = padd(fp, c, &r, &step);
        }
        found.retain(|&n| n + w >= mid && n <= mid + w);
        found.sort();
        found.dedup();
        let cs: Vec<u64> = match cands.take() {
            None => found,
            Some(old) => old.into_iter().filter(|n| found.contains(n)).collect(),
        };
        if cs.len() == 1 {
            return cs[0];
        }
        if cs.is_empty() {
            cands = None;
        } else {
            cands = Some(cs);
        }
    }
}

// ---------------------------------------------------------------- isogenies

pub trait Isogeny<F: Field> {
    fn domain(&self) -> &Curve<F::E>;
    fn codomain(&self) -> &Curve<F::E>;
    fn degree(&self) -> u64;
    /// x-coordinate map; None means the point is in the kernel.
    fn eval_x(&self, f: &F, x: F::E) -> Option<F::E>;
    fn eval(&self, f: &F, p: &Pt<F::E>) -> Pt<F::E>;
}

/// Normalised isogeny given by x -> num(x)/den(x), den = ker^2.
#[derive(Clone, Debug)]
pub struct RatIsogeny<F: Field> {
    pub dom: Curve<F::E>,
    pub cod: Curve<F::E>,
    pub deg: u64,
    pub ker: Poly<F>,
    pub num: Poly<F>,
    pub den: Poly<F>,
}

impl<F: Field> Isogeny<F> for RatIsogeny<F> {
    fn domain(&self) -> &Curve<F::E> {
        &self.dom
    }
    fn codomain(&self) -> &Curve<F::E> {
        &self.cod
    }
    fn degree(&self) -> u64 {
        self.deg
    }
    fn eval_x(&self, f: &F, x: F::E) -> Option<F::E> {
        let d = poly::eval(f, &self.den, x);
        if f.is_zero(d) {
            None
        } else {
            Some(f.div(poly::eval(f, &self.num, x), d))
        }
    }
    fn eval(&self, f: &F, p: &Pt<F::E>) -> Pt<F::E> {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => {
                let (d, dd) = poly::eval_with_derivative(f, &self.den, x);
                if f.is_zero(d) {
                    return Pt::Inf;
                }
                let (n, dn) = poly::eval_with_derivative(f, &self.num, x);
                // f = n/d, f' = (n' d - n d') / d^2: one inversion
                let di = f.inv(d);
                let fx = f.mul(n, di);
                let fp = f.mul(f.sub(f.mul(dn, d), f.mul(n, dd)), f.sq(di));
                Pt::Aff(fx, f.mul(y, fp))
            }
        }
    }
}

/// Scale of curve model: map E(a,b) -> E(u^-4 a, u^-6 b) is an isomorphism over F(u).
pub fn iso_exists_over<F: Field>(f: &F, c1: &Curve<F::E>, c2: &Curve<F::E>) -> Option<F::E> {
    // find u with u^4 a2 = a1 ... brute-force is unavailable; use sqrt-free test:
    // u^2 = (b1 a2)/(b2 a1) when a,b != 0.
    if f.is_zero(c1.a) || f.is_zero(c2.a) || f.is_zero(c1.b) || f.is_zero(c2.b) {
        return None;
    }
    let u2 = f.div(f.mul(c1.b, c2.a), f.mul(c2.b, c1.a));
    Some(u2)
}

/// Verify the isogeny as an algebraic identity: (x^3+ax+b) f'(x)^2 == f^3 + a' f + b'
/// at `trials` random x, and that deg matches.
pub fn check_x_identity<F: Field, I: Isogeny<F> + ?Sized>(
    f: &F,
    iso: &I,
    rng: &mut Rng,
    trials: usize,
) -> bool {
    let (c1, c2) = (iso.domain(), iso.codomain());
    let h = f.from_u64(1);
    let _ = h;
    for _ in 0..trials {
        let x = f.random(rng);
        // finite-difference-free derivative: use the point map with a formal y
        // choose y^2 = rhs(x); the image of (x,y) must lie on c2.  Over F_p y may not exist,
        // so compare through y^2: y'^2 = y^2 f'^2.
        let (fx, fpx) = match deriv_x(f, iso, x) {
            Some(v) => v,
            None => continue,
        };
        let lhs = f.mul(rhs(f, c1, x), f.mul(fpx, fpx));
        let r = rhs(f, c2, fx);
        if lhs != r {
            return false;
        }
    }
    true
}

/// (f(x), f'(x)) from the point map using the formal second coordinate y=1 trick:
/// eval at (x, y) returns (f(x), y f'(x)); we call it with a symbolic y of 1.
fn deriv_x<F: Field, I: Isogeny<F> + ?Sized>(f: &F, iso: &I, x: F::E) -> Option<(F::E, F::E)> {
    match iso.eval(f, &Pt::Aff(x, f.one())) {
        Pt::Inf => None,
        Pt::Aff(fx, fpx) => Some((fx, fpx)),
    }
}

/// Homomorphism check on F_p-points: phi(P+Q) = phi(P)+phi(Q) and image on the codomain.
pub fn check_homomorphism(fp: &Zp, iso: &dyn Isogeny<Zp>, rng: &mut Rng, trials: usize) -> bool {
    let (c1, c2) = (*iso.domain(), *iso.codomain());
    for _ in 0..trials {
        let p = c1.random_point(fp, rng);
        let q = c1.random_point(fp, rng);
        let (ip, iq) = (iso.eval(fp, &p), iso.eval(fp, &q));
        if !on_curve(fp, &c2, &ip) || !on_curve(fp, &c2, &iq) {
            return false;
        }
        let s = padd(fp, &c1, &p, &q);
        if iso.eval(fp, &s) != padd(fp, &c2, &ip, &iq) {
            return false;
        }
    }
    true
}

/// Kernel polynomial from a generator: prod_{k=1}^{(l-1)/2} (x - x(kP)).
pub fn kernel_poly_from_point(
    fp: &Zp,
    c: &Curve<u64>,
    p: &Pt<u64>,
    ell: u64,
) -> (Poly<Zp>, Vec<Pt<u64>>) {
    kernel_poly_from_point_f(fp, c, p, ell)
}

/// `kernel_poly_from_point` over any field.
pub fn kernel_poly_from_point_f<F: Field>(
    fp: &F,
    c: &Curve<F::E>,
    p: &Pt<F::E>,
    ell: u64,
) -> (Poly<F>, Vec<Pt<F::E>>) {
    let n = ((ell - 1) / 2) as usize;
    let mut pts = Vec::with_capacity(n);
    let mut cur = *p;
    for _ in 0..n {
        pts.push(cur);
        cur = padd(fp, c, &cur, p);
    }
    let xs: Vec<F::E> = pts
        .iter()
        .map(|q| match q {
            Pt::Aff(x, _) => *x,
            Pt::Inf => panic!("point order too small"),
        })
        .collect();
    (poly::from_roots(fp, &xs), pts)
}
