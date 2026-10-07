//! Velu (1971): enumerate the kernel, O(l) work to build the codomain and O(l) per evaluation.
use crate::curve::*;
use crate::field::Field;

pub struct VeluIso<F: Field> {
    pub dom: Curve<F::E>,
    pub cod: Curve<F::E>,
    pub deg: u64,
    /// (x_Q, v_Q, u_Q) for representatives Q of (G\{O})/{+-1}
    pub reps: Vec<(F::E, F::E, F::E)>,
}

/// Odd-degree Velu from the affine points {P, 2P, ..., ((l-1)/2) P}.
pub fn velu<F: Field>(f: &F, e: &Curve<F::E>, reps: &[Pt<F::E>], ell: u64) -> VeluIso<F> {
    let mut t = f.zero();
    let mut w = f.zero();
    let mut out = Vec::with_capacity(reps.len());
    for q in reps {
        let (x, y) = match *q {
            Pt::Aff(x, y) => (x, y),
            Pt::Inf => panic!("kernel point at infinity"),
        };
        let x2 = f.mul(x, x);
        let v = f.add(f.mul(f.from_u64(6), x2), f.mul(f.from_u64(2), e.a));
        let u = f.mul(f.from_u64(4), f.mul(y, y));
        t = f.add(t, v);
        w = f.add(w, f.add(u, f.mul(x, v)));
        out.push((x, v, u));
    }
    let a2 = f.sub(e.a, f.mul(f.from_u64(5), t));
    let b2 = f.sub(e.b, f.mul(f.from_u64(7), w));
    VeluIso { dom: *e, cod: Curve::new(a2, b2), deg: ell, reps: out }
}

impl<F: Field> Isogeny<F> for VeluIso<F> {
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
        self.eval(f, &Pt::Aff(x, f.one())).into_x()
    }
    fn eval(&self, f: &F, p: &Pt<F::E>) -> Pt<F::E> {
        let (x, y) = match *p {
            Pt::Inf => return Pt::Inf,
            Pt::Aff(x, y) => (x, y),
        };
        let mut fx = x;
        let mut fd = f.one();
        for &(xq, v, u) in &self.reps {
            let d = f.sub(x, xq);
            if f.is_zero(d) {
                return Pt::Inf;
            }
            let di = f.inv(d);
            let di2 = f.mul(di, di);
            let di3 = f.mul(di2, di);
            fx = f.add(fx, f.add(f.mul(v, di), f.mul(u, di2)));
            fd = f.sub(fd, f.add(f.mul(v, di2), f.mul(f.from_u64(2), f.mul(u, di3))));
        }
        Pt::Aff(fx, f.mul(y, fd))
    }
}

pub trait IntoX<E> {
    fn into_x(self) -> Option<E>;
}
impl<E> IntoX<E> for Pt<E> {
    fn into_x(self) -> Option<E> {
        match self {
            Pt::Inf => None,
            Pt::Aff(x, _) => Some(x),
        }
    }
}

/// Velu from a generator: enumerate kP, k = 1..(l-1)/2, then apply the formulas.
pub fn velu_from_point<F: Field>(f: &F, e: &Curve<F::E>, p: &Pt<F::E>, ell: u64) -> VeluIso<F> {
    let n = ((ell - 1) / 2) as usize;
    let mut pts = Vec::with_capacity(n);
    let mut cur = *p;
    for _ in 0..n {
        pts.push(cur);
        cur = padd(f, e, &cur, p);
    }
    velu(f, e, &pts, ell)
}
