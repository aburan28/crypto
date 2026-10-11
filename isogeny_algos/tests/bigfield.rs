use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::bmss;
use isogeny_algos::kernel::{
    kohel::kohel, sqrt_velu::sqrt_velu_fast, velu::velu_cyclic, xonly::velu_xonly_fast,
};
use isogeny_algos::testdata::supersingular_workload;

fn kernel_algorithms_agree<const N: usize>(bits: usize, ells: &[u64], seed: u64) {
    let mut rng = Rng::new(seed);
    let w = supersingular_workload::<N>(ells, bits, &mut rng);
    let f = &w.f;
    assert_eq!(f.modulus().bits(), bits);
    for &(ell, p) in &w.pts {
        let Pt::Aff(x0, _) = p else { unreachable!() };
        let r = velu_xonly_fast(f, &w.e, x0, ell).unwrap();
        // the codomain is supersingular too: a random point times p+1 is the identity
        let q = random_point_f(f, &r.cod, &mut rng);
        assert_eq!(
            pmul_big(f, &r.cod, &q, &f.modulus().add_small(1)),
            Pt::Inf,
            "l = {ell}"
        );
        let sv = sqrt_velu_fast(f, &w.e, x0, ell).unwrap();
        assert_eq!(sv.cod, r.cod, "sqrt-Velu, l = {ell}");
        let t = random_point_f(f, &w.e, &mut rng);
        assert_eq!(sv.eval(f, &t), r.eval(f, &t), "sqrt-Velu image, l = {ell}");
        if ell <= 101 {
            assert_eq!(velu_cyclic(f, &w.e, &p, ell).cod, r.cod);
            let (g, _) = kernel_poly_from_point_f(f, &w.e, &p, ell);
            let kh = kohel(f, &w.e, &g, ell);
            assert_eq!(kh.cod, r.cod, "Kohel, l = {ell}");
            let iso = bmss::isogeny(
                f,
                bmss::Method::FastElkiesPrime,
                &w.e,
                &kh.cod,
                ell as usize,
                None,
            )
            .expect("BMSS");
            assert_eq!(iso.ker, g, "fastElkies', l = {ell}");
        }
    }
}

#[test]
fn kernel_algorithms_agree_at_256_bits() {
    kernel_algorithms_agree::<4>(256, &[3, 5, 7, 31, 101, 1009], 501);
}

#[test]
fn kernel_algorithms_agree_at_512_bits() {
    kernel_algorithms_agree::<8>(511, &[3, 13, 61, 401], 502);
}
