use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{dual::*, elkies, modpoly::Phi};
use isogeny_algos::testdata::*;

#[test]
fn dual_composes_to_multiplication_by_l() {
    let mut rng = Rng::new(61);
    let fp = Zp::new(next_prime(1u64 << 30));
    for &ell in &[3usize, 5, 7, 11, 13] {
        let phi = Phi::compute(&fp, ell);
        let mut checked = 0;
        while checked < 3 {
            let (e, _) = curve_with_point(&fp, ell as u64, &mut rng);
            let j = jinv(&fp, &e);
            if j == 0 || j == 1728 {
                continue;
            }
            for iso in elkies::elkies_isogenies(&fp, &phi, &e, &mut rng) {
                let Some(d) = dual_isogeny(&fp, &phi, &iso) else {
                    continue;
                };
                assert_eq!(*d.codomain(), e);
                for _ in 0..5 {
                    let p = e.random_point(&fp, &mut rng);
                    let back = d.eval(&fp, &iso.eval(&fp, &p));
                    let lp = pmul(&fp, &e, &p, ell as u128);
                    assert_eq!(back, lp, "dual(phi(P)) must equal [l]P, l = {ell}");
                }
                // and in the other direction: phi(dual(Q)) = [l]Q on E'
                for _ in 0..3 {
                    let q = iso.cod.random_point(&fp, &mut rng);
                    let fwd = iso.eval(&fp, &d.eval(&fp, &q));
                    assert_eq!(fwd, pmul(&fp, &iso.cod, &q, ell as u128));
                }
                checked += 1;
            }
        }
    }
}
