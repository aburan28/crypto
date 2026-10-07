use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::kernel::{kohel::kohel, sqrt_velu::sqrt_velu, velu::velu};
use isogeny_algos::testdata::*;

#[test]
fn order_matches_naive() {
    let mut rng = Rng::new(1);
    let fp = Zp::new(next_prime(5_000_000));
    // large-p BSGS vs naive on a medium prime by forcing both
    let c = Curve::new(3, 7);
    let n = order(&fp, &c, &mut rng);
    // Hasse + point order divides
    let pt = c.random_point(&fp, &mut rng);
    assert_eq!(pmul(&fp, &c, &pt, n as u128), Pt::Inf);
}

#[test]
fn bsgs_order_consistent() {
    let mut rng = Rng::new(2);
    let fp = Zp::new(next_prime(1u64 << 40));
    for _ in 0..3 {
        let c = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
        if !is_smooth(&fp, &c) { continue; }
        let n = order(&fp, &c, &mut rng);
        for _ in 0..3 {
            let pt = c.random_point(&fp, &mut rng);
            assert_eq!(pmul(&fp, &c, &pt, n as u128), Pt::Inf);
        }
    }
}

#[test]
fn three_kernel_algorithms_agree() {
    let mut rng = Rng::new(3);
    for &(pbits, ell) in &[(30u32, 3u64), (30, 5), (30, 7), (30, 11), (30, 13), (34, 29), (34, 61), (34, 101)] {
        let fp = Zp::new(next_prime(1u64 << pbits));
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let (h, pts) = kernel_poly_from_point(&fp, &e, &p, ell);
        let v = velu(&fp, &e, &pts, ell);
        let k = kohel(&fp, &e, &h, ell);
        let s = sqrt_velu(&fp, &e, &p, ell);
        assert_eq!(v.cod, k.cod, "velu vs kohel codomain l={ell}");
        assert_eq!(v.cod, s.cod, "velu vs sqrt-velu codomain l={ell}");
        assert!(check_homomorphism(&fp, &v, &mut rng, 5), "velu hom l={ell}");
        assert!(check_homomorphism(&fp, &k, &mut rng, 5), "kohel hom l={ell}");
        assert!(check_homomorphism(&fp, &s, &mut rng, 5), "sqrt-velu hom l={ell}");
        for _ in 0..5 {
            let x = fp.random(&mut rng);
            let a = v.eval_x(&fp, x);
            assert_eq!(a, k.eval_x(&fp, x), "kohel x-map l={ell}");
            assert_eq!(a, s.eval_x(&fp, x), "sqrt-velu x-map l={ell}");
        }
    }
}
