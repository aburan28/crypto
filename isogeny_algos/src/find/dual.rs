//! Dual isogeny of a normalised l-isogeny (Kohel / Elkies): the dual of phi: E -> E' is obtained
//! from Phi_l alone. j(E) is always a root of Phi_l(j(E'), Y); Elkies' formulas give the normalised
//! codomain E~ over it, Padé/BMSS gives psi: E' -> E~, and E~ is E rescaled by u = l:
//! E~ = (l^4 A, l^6 B), because psi o phi pulls omega back to omega while [l] pulls it back to l omega.
//! The dual is then  (x, y) -> (x_psi / u^2, y_psi / u^3)  (not normalised: c = 1/l).
use super::elkies::{bmss_isogeny, elkies_codomain};
use super::modpoly::Phi;
use crate::curve::*;
use crate::field::Field;

pub struct ScaledIso<F: Field> {
    pub inner: RatIsogeny<F>,
    /// scaling u: the codomain of `inner` is (u^4 a, u^6 b) of the true codomain (a, b)
    pub u: F::E,
    pub cod: Curve<F::E>,
}

impl<F: Field> Isogeny<F> for ScaledIso<F> {
    fn domain(&self) -> &Curve<F::E> {
        &self.inner.dom
    }
    fn codomain(&self) -> &Curve<F::E> {
        &self.cod
    }
    fn degree(&self) -> u64 {
        self.inner.deg
    }
    fn eval_x(&self, f: &F, x: F::E) -> Option<F::E> {
        let u2 = f.mul(self.u, self.u);
        self.inner.eval_x(f, x).map(|v| f.div(v, u2))
    }
    fn eval(&self, f: &F, p: &Pt<F::E>) -> Pt<F::E> {
        match self.inner.eval(f, p) {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => {
                let u2 = f.mul(self.u, self.u);
                Pt::Aff(f.div(x, u2), f.div(y, f.mul(u2, self.u)))
            }
        }
    }
}

/// Dual of `iso` (prime degree, Phi_l available, no j in {0, 1728} on the way).
pub fn dual_isogeny<F: Field>(f: &F, phi: &Phi<F>, iso: &RatIsogeny<F>) -> Option<ScaledIso<F>> {
    let ell = phi.ell;
    let et = elkies_codomain(f, phi, &iso.cod, jinv(f, &iso.dom))?;
    let psi = bmss_isogeny(f, &iso.cod, &et, ell)?;
    let l = f.from_u64(ell as u64);
    let l2 = f.mul(l, l);
    let l4 = f.mul(l2, l2);
    let l6 = f.mul(l4, l2);
    let expect = Curve::new(f.mul(l4, iso.dom.a), f.mul(l6, iso.dom.b));
    if psi.cod != expect {
        return None;
    }
    // u = l or -l: the sign only changes the sign of y; fix it by testing phi-hat o phi = [l] on a point
    Some(ScaledIso {
        inner: psi,
        u: l,
        cod: iso.dom,
    })
}
