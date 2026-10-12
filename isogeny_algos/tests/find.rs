use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{divpoly, elkies, modpoly::Phi};
use isogeny_algos::kernel::kohel::kohel;
use isogeny_algos::testdata::*;

#[test]
fn phi2_matches_known() {
    let fp = Zp::new(next_prime(1u64 << 30));
    let phi = Phi::compute(&fp, 2);
    let m = |x: i128| -> u64 { (((x % fp.p as i128) + fp.p as i128) % fp.p as i128) as u64 };
    assert_eq!(phi.c[3][0], 1);
    assert_eq!(phi.c[2][2], m(-1));
    assert_eq!(phi.c[2][1], m(1488));
    assert_eq!(phi.c[2][0], m(-162000));
    assert_eq!(phi.c[1][1], m(40773375));
    assert_eq!(phi.c[1][0], m(8748000000));
    assert_eq!(phi.c[0][0], m(-157464000000000));
}

#[test]
fn divpoly_kernels_are_isogenies_and_match_velu() {
    let mut rng = Rng::new(5);
    for &ell in &[3u64, 5, 7, 11] {
        let fp = Zp::new(next_prime(1u64 << 28));
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let (h, _) = kernel_poly_from_point(&fp, &e, &p, ell);
        let ks = divpoly::kernel_polys(&fp, &e, ell, &mut rng);
        assert!(ks.contains(&h), "known kernel poly missing, l={ell}");
        for k in &ks {
            let iso = kohel(&fp, &e, k, ell);
            assert!(check_x_identity(&fp, &iso, &mut rng, 5));
        }
    }
}

#[test]
fn elkies_bmss_matches_kohel() {
    let mut rng = Rng::new(6);
    for &ell in &[3u64, 5, 7, 11, 13] {
        let fp = Zp::new(next_prime(1u64 << 28));
        let phi = Phi::compute(&fp, ell as usize);
        let mut checked = 0;
        for _ in 0..40 {
            let (e, p) = curve_with_point(&fp, ell, &mut rng);
            let j = jinv(&fp, &e);
            if j == 0 || j == 1728 {
                continue;
            }
            // Phi vanishes at (j, j(codomain))
            let (h, _) = kernel_poly_from_point(&fp, &e, &p, ell);
            let kh = kohel(&fp, &e, &h, ell);
            assert_eq!(
                phi.eval(&fp, j, jinv(&fp, &kh.cod)),
                0,
                "Phi_l(j, j') != 0, l={ell}"
            );
            let isos = elkies::elkies_isogenies(&fp, &phi, &e, &mut rng);
            let m = isos.iter().find(|i| i.ker == h);
            let m = m.expect("BMSS kernel for known subgroup not found");
            assert_eq!(m.cod, kh.cod);
            assert_eq!(m.num, kh.num);
            assert!(check_x_identity(&fp, m, &mut rng, 5));
            checked += 1;
            if checked == 3 {
                break;
            }
        }
        assert!(checked > 0);
    }
}
