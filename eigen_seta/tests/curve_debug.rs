mod common;
use common::*;
use eigen_seta::cm::{mul_vartime, p_plus_1, random_point};
use eigen_seta::weier::{Isog2, WIso};
use rand_core::OsRng;

#[test]
fn debug_curves() {
    let mut rng = OsRng;
    let ed = ed_curve();
    let p = random_point(&ed, &mut rng);
    assert!(bool::from(ed.on_curve(&p)), "random point not on curve");
    let q = random_point(&ed, &mut rng);
    let s = ed.add(&p, &q);
    assert!(bool::from(ed.on_curve(&s)), "sum not on curve");
    let d = ed.add(&p, &p);
    assert!(bool::from(ed.on_curve(&d)), "double not on curve");
    let n = p_plus_1::<F>();
    let z = mul_vartime(&ed, &p, &n);
    assert!(bool::from(z.inf), "[p+1]P != O");
    let half = &n >> 1;
    let t = mul_vartime(&ed, &p, &half);
    eprintln!("2-torsion inf? {}  y zero? {}", bool::from(t.inf), bool::from(t.y.is_zero()));
    if !bool::from(t.inf) {
        let phi = Isog2::new(&ed, &t.x);
        let im = phi.eval(&q);
        eprintln!("image on codomain: {}", bool::from(phi.codomain.on_curve(&im)));
        let iso = WIso::find(&phi.codomain, &phi.codomain);
        eprintln!("self-iso u2 = {:?}", iso.u2);
        eprintln!("j(ed) = {:?}\nj(cod) = {:?}", ed.j_invariant(), phi.codomain.j_invariant());
    }
}
