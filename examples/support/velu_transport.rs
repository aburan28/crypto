//! Shared checked Vélu point evaluation; extracted unchanged from the prime61 pilot.
use crypto_lib::cryptanalysis::isogeny_walk::{
    curve::Model,
    field::{Fe, Field},
    kernel,
    poly::{self, Poly},
};
use num_traits::ToPrimitive;
fn word(f: &Field, x: &Fe) -> u64 {
    f.to_big(x).to_u64().unwrap()
}
/// Known Vélu formula, with sums evaluated as traces in F_p[T]/h(T).
/// This handles Frobenius-stable kernels without rational kernel points.
pub(crate) struct KernelMap<'a> {
    f: &'a Field,
    h: Poly,
    v: Poly,
    u: Poly,
    sums: Vec<Fe>,
    pub(crate) u2: Fe,
    pub(crate) u3: Fe,
}

impl<'a> KernelMap<'a> {
    pub(crate) fn new(f: &'a Field, source: &Model, target: &Model, ell: u64, h: Poly) -> Self {
        let codomain = kernel::verify_kernel(f, source, ell, &h).expect("kernel certificate");
        let u2 = kernel::link_to_target(f, &codomain, target).expect("codomain isomorphism");
        let u_root = f.sqrt(&u2).expect("u2 is square");
        let u3 = f.mul(&u2, &u_root);
        let v = poly::rem(
            f,
            &vec![f.mul(&f.from_u64(2), &source.a), f.zero(), f.from_u64(6)],
            &h,
        );
        let u = poly::rem(
            f,
            &vec![
                f.mul(&f.from_u64(4), &source.b),
                f.mul(&f.from_u64(4), &source.a),
                f.zero(),
                f.from_u64(4),
            ],
            &h,
        );
        let d = h.len() - 1;
        let mut sums = vec![f.from_u64(d as u64)];
        sums.extend(poly::power_sums(f, &h, d.saturating_sub(1)));
        Self {
            f,
            h,
            v,
            u,
            sums,
            u2,
            u3,
        }
    }

    fn trace(&self, a: &Poly) -> Fe {
        a.iter().zip(&self.sums).fold(self.f.zero(), |acc, (a, s)| {
            self.f.add(&acc, &self.f.mul(a, s))
        })
    }

    pub(crate) fn eval(&self, point: Option<(u64, u64)>) -> Option<(u64, u64)> {
        let (x, y) = point?;
        let f = self.f;
        let x = f.from_u64(x);
        let y = f.from_u64(y);
        let z = poly::invmod(f, &vec![x, f.neg(&f.one())], &self.h)?;
        let z2 = poly::mulmod(f, &z, &z, &self.h);
        let z3 = poly::mulmod(f, &z2, &z, &self.h);
        let dx = poly::add(
            f,
            &poly::mulmod(f, &self.v, &z, &self.h),
            &poly::mulmod(f, &self.u, &z2, &self.h),
        );
        let dy = poly::add(
            f,
            &poly::mulmod(f, &self.v, &z2, &self.h),
            &poly::scale(f, &poly::mulmod(f, &self.u, &z3, &self.h), &f.from_u64(2)),
        );
        let xx = f.mul(&self.u2, &f.add(&x, &self.trace(&dx)));
        let yy = f.mul(&self.u3, &f.mul(&y, &f.sub(&f.one(), &self.trace(&dy))));
        Some((word(f, &xx), word(f, &yy)))
    }
}
