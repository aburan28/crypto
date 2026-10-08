//! General Weierstrass curves y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 over any field, and
//! Vélu's formulas in the characteristic-independent form (Vélu 1971, as in Silverman III.4 and
//! Kohel's thesis), so that isogenies can be computed in characteristic 2 and 3 where the
//! short-form coefficients (which divide by 2, 3, 6, 24, 27) are unavailable.
use crate::curve::Pt;
use crate::field::Field;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GWCurve<E> {
    pub a1: E,
    pub a2: E,
    pub a3: E,
    pub a4: E,
    pub a6: E,
}

impl<E: Copy + PartialEq> GWCurve<E> {
    pub fn new(a: [E; 5]) -> Self {
        GWCurve { a1: a[0], a2: a[1], a3: a[2], a4: a[3], a6: a[4] }
    }
    pub fn on_curve<F: Field<E = E>>(&self, f: &F, p: &Pt<E>) -> bool {
        match *p {
            Pt::Inf => true,
            Pt::Aff(x, y) => {
                let l = f.add(f.add(f.sq(y), f.mul(self.a1, f.mul(x, y))), f.mul(self.a3, y));
                let r = f.add(f.add(f.add(f.mul(f.sq(x), x), f.mul(self.a2, f.sq(x))), f.mul(self.a4, x)), self.a6);
                l == r
            }
        }
    }
    pub fn neg<F: Field<E = E>>(&self, f: &F, p: &Pt<E>) -> Pt<E> {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => Pt::Aff(x, f.sub(f.sub(f.neg(y), f.mul(self.a1, x)), self.a3)),
        }
    }
    pub fn add<F: Field<E = E>>(&self, f: &F, p: &Pt<E>, q: &Pt<E>) -> Pt<E> {
        let (x1, y1) = match *p {
            Pt::Inf => return *q,
            Pt::Aff(x, y) => (x, y),
        };
        let (x2, y2) = match *q {
            Pt::Inf => return *p,
            Pt::Aff(x, y) => (x, y),
        };
        if x1 == x2 && y2 == self.neg_y(f, x1, y1) {
            return Pt::Inf;
        }
        let lam = if x1 == x2 {
            // doubling: (3x^2 + 2 a2 x + a4 - a1 y) / (2y + a1 x + a3)
            let num = f.sub(f.add(f.add(f.mul(f.from_u64(3), f.sq(x1)), f.mul(f.mul(f.from_u64(2), self.a2), x1)), self.a4), f.mul(self.a1, y1));
            let den = f.add(f.add(f.mul(f.from_u64(2), y1), f.mul(self.a1, x1)), self.a3);
            f.div(num, den)
        } else {
            f.div(f.sub(y2, y1), f.sub(x2, x1))
        };
        // x3 = lam^2 + a1 lam - a2 - x1 - x2;  y3 = -(lam (x3 - x1) + y1) - a1 x3 - a3
        let x3 = f.sub(f.sub(f.sub(f.add(f.sq(lam), f.mul(self.a1, lam)), self.a2), x1), x2);
        let y3 = f.sub(f.sub(f.neg(f.add(f.mul(lam, f.sub(x3, x1)), y1)), f.mul(self.a1, x3)), self.a3);
        Pt::Aff(x3, y3)
    }
    fn neg_y<F: Field<E = E>>(&self, f: &F, x: E, y: E) -> E {
        f.sub(f.sub(f.neg(y), f.mul(self.a1, x)), self.a3)
    }
    pub fn mul<F: Field<E = E>>(&self, f: &F, p: &Pt<E>, mut k: u64) -> Pt<E> {
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
    /// b2, b4, b6, b8, c4, discriminant, j (standard formulas; integer constants reduce in F).
    pub fn disc<F: Field<E = E>>(&self, f: &F) -> E {
        let (b2, b4, b6, b8) = self.b_invariants(f);
        let c = |v: u64| f.from_u64(v);
        // Delta = -b2^2 b8 - 8 b4^3 - 27 b6^2 + 9 b2 b4 b6
        let t1 = f.neg(f.mul(f.sq(b2), b8));
        let t2 = f.mul(c(8), f.mul(b4, f.sq(b4)));
        let t3 = f.mul(c(27), f.sq(b6));
        let t4 = f.mul(c(9), f.mul(b2, f.mul(b4, b6)));
        f.add(f.sub(f.sub(t1, t2), t3), t4)
    }
    fn b_invariants<F: Field<E = E>>(&self, f: &F) -> (E, E, E, E) {
        let c = |v: u64| f.from_u64(v);
        let b2 = f.add(f.sq(self.a1), f.mul(c(4), self.a2));
        let b4 = f.add(f.mul(c(2), self.a4), f.mul(self.a1, self.a3));
        let b6 = f.add(f.sq(self.a3), f.mul(c(4), self.a6));
        let b8 = f.add(
            f.sub(f.add(f.mul(f.sq(self.a1), self.a6), f.mul(c(4), f.mul(self.a2, self.a6))), f.mul(self.a1, f.mul(self.a3, self.a4))),
            f.sub(f.mul(self.a2, f.sq(self.a3)), f.sq(self.a4)),
        );
        (b2, b4, b6, b8)
    }
    pub fn j<F: Field<E = E>>(&self, f: &F) -> E {
        let (b2, b4, _, _) = self.b_invariants(f);
        let c4 = f.sub(f.sq(b2), f.mul(f.from_u64(24), b4));
        f.div(f.mul(f.sq(c4), c4), self.disc(f))
    }
}

/// An odd-degree Vélu isogeny with cyclic kernel, in the general Weierstrass form.
pub struct GWVelu<E> {
    pub dom: GWCurve<E>,
    pub cod: GWCurve<E>,
    pub deg: u64,
    /// per kernel pair {Q, -Q}: (xQ, yQ, tQ, uQ)
    pub data: Vec<(E, E, E, E)>,
}

/// Vélu from a generator P of a cyclic subgroup of odd order ell (coprime to char), using the
/// representatives P, 2P, ..., ((ell-1)/2) P (no 2-torsion in an odd-order group).
pub fn gw_velu_cyclic<F: Field>(f: &F, e: &GWCurve<F::E>, p: &Pt<F::E>, ell: u64) -> GWVelu<F::E> {
    let mut q = *p;
    let mut data = Vec::with_capacity((ell / 2) as usize);
    let (mut t, mut w) = (f.zero(), f.zero());
    for i in 0..ell / 2 {
        if i > 0 {
            q = e.add(f, &q, p);
        }
        let (xq, yq) = match q {
            Pt::Aff(x, y) => (x, y),
            Pt::Inf => panic!("kernel generator order mismatch"),
        };
        // gxQ = 3 xQ^2 + 2 a2 xQ + a4 - a1 yQ;  gyQ = -2 yQ - a1 xQ - a3
        let gx = f.sub(f.add(f.add(f.mul(f.from_u64(3), f.sq(xq)), f.mul(f.mul(f.from_u64(2), e.a2), xq)), e.a4), f.mul(e.a1, yq));
        let gy = f.sub(f.sub(f.neg(f.mul(f.from_u64(2), yq)), f.mul(e.a1, xq)), e.a3);
        // Q not 2-torsion: tQ = 2 gxQ - a1 gyQ, uQ = gyQ^2
        let tq = f.sub(f.mul(f.from_u64(2), gx), f.mul(e.a1, gy));
        let uq = f.sq(gy);
        t = f.add(t, tq);
        w = f.add(w, f.add(uq, f.mul(xq, tq)));
        let _ = gx;
        data.push((xq, yq, tq, uq));
    }
    let c = |v: u64| f.from_u64(v);
    let b2 = f.add(f.sq(e.a1), f.mul(c(4), e.a2));
    let cod = GWCurve {
        a1: e.a1,
        a2: e.a2,
        a3: e.a3,
        a4: f.sub(e.a4, f.mul(c(5), t)),
        a6: f.sub(f.sub(e.a6, f.mul(b2, t)), f.mul(c(7), w)),
    };
    GWVelu { dom: *e, cod, deg: ell, data }
}

impl<E: Copy + PartialEq> GWVelu<E> {
    pub fn eval<F: Field<E = E>>(&self, f: &F, p: &Pt<E>) -> Pt<E> {
        let (x, y) = match *p {
            Pt::Inf => return Pt::Inf,
            Pt::Aff(x, y) => (x, y),
        };
        let mut xx = x;
        // X'(x) = 1 - sum tQ/(x-xQ)^2 + 2 uQ/(x-xQ)^3
        let mut dxdx = f.one();
        for &(xq, _yq, tq, uq) in &self.data {
            let dx = f.sub(x, xq);
            if f.is_zero(dx) {
                return Pt::Inf; // in the kernel
            }
            let di = f.inv(dx);
            let di2 = f.sq(di);
            let di3 = f.mul(di2, di);
            xx = f.add(xx, f.add(f.mul(tq, di), f.mul(uq, di2)));
            dxdx = f.sub(dxdx, f.add(f.mul(tq, di2), f.mul(f.mul(f.from_u64(2), uq), di3)));
        }
        // the Vélu isogeny is normalised: dX/(2Y + A1 X + A3) = dx/(2y + a1 x + a3), so
        // 2Y + A1 X + A3 = (2y + a1 x + a3) X'(x), hence Y = ((2y+a1x+a3) X'(x) - A1 X - A3)/2
        let two_y = f.add(f.add(f.mul(f.from_u64(2), y), f.mul(self.dom.a1, x)), self.dom.a3);
        let yy = f.div(
            f.sub(f.sub(f.mul(two_y, dxdx), f.mul(self.cod.a1, xx)), self.cod.a3),
            f.from_u64(2),
        );
        Pt::Aff(xx, yy)
    }
}
