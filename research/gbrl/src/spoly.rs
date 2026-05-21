use crate::poly::Poly;

// S-polynomial: lcm(LM(f),LM(g))/LT(f) * f  -  lcm/LT(g) * g
// (with leading-term cancellation arranged via coefficient cross-multiplication).
pub fn spoly(f: &Poly, g: &Poly) -> Poly {
    let ft = f.lt().expect("spoly: f must be nonzero");
    let gt = g.lt().expect("spoly: g must be nonzero");
    let lcm = ft.mono.lcm(&gt.mono);
    let mf = lcm.div(&ft.mono).expect("lcm/LM(f) divides cleanly");
    let mg = lcm.div(&gt.mono).expect("lcm/LM(g) divides cleanly");
    // Cross-multiply by the other leading coefficient to cancel LTs without inversion.
    let p1 = f.mul_term(gt.coef, &mf);
    let p2 = g.mul_term(ft.coef, &mg);
    p1.sub(&p2)
}
