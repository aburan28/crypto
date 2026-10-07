use isogeny_algos::bigint::Big;
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::montgomery::*;
use isogeny_algos::kernel::sqrt_velu_mont::SqrtVeluMont;
use isogeny_algos::path::csidh::Csidh;

/// On CSIDH-512 curves (p + 1 = 4 prod l): kernel point of order l from a random point, then
/// sqrt-Velu codomain and x-map against x-only Velu.
#[test]
fn sqrt_velu_mont_matches_velu_on_csidh512() {
    let cs = Csidh::<FpM<8>>::csidh512();
    let f = &cs.fp;
    let mut rng = Rng::new(601);
    // a random curve in the isogeny class
    let mut e = vec![0i32; 74];
    e[0] = 1;
    e[3] = -1;
    let a = cs.action_fast(f.zero(), &e, &mut rng);
    for &ell in &[3u64, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587] {
        let cof = cs.p1.divrem_small(ell).0;
        let k = loop {
            let x = f.random(&mut rng);
            let k = ladder_xz(f, a24(f, a), (x, f.one()), &cof);
            if !f.is_zero(k.1) {
                break k;
            }
        };
        let d = ((ell - 1) / 2) as usize;
        let ms = multiples(f, a24(f, a), k, d);
        let a_ref = velu_codomain(f, a, &ms, ell);
        let sv = SqrtVeluMont::new(f, a, k, ell).expect("kps");
        assert_eq!(sv.codomain(f), a_ref, "codomain, l = {ell} (b = {}, b' = {})", sv.b, sv.bp);
        let u = f.random(&mut rng);
        assert_eq!(sv.eval(f, u), isog_x(f, &ms, u), "x-map, l = {ell}");
    }
    let _ = Big::from_u64(1);
}

#[test]
fn sqrt_velu_mont_large_degree() {
    // p = 4 * 3 * 1009 * 10007 * c - 1 over 256 bits: kernels of degree 1009 and 10007
    let mut rng = Rng::new(602);
    let w = isogeny_algos::testdata::supersingular_workload::<4>(&[3, 1009, 10007], 256, &mut rng);
    let f = &w.f;
    // Montgomery A = 0 curve y^2 = x^3 + x over the same p has p + 1 points
    let a = f.zero();
    let p1 = f.modulus().add_small(1);
    for &ell in &[1009u64, 10007] {
        let cof = p1.divrem_small(ell).0;
        let k = loop {
            let x = f.random(&mut rng);
            let k = ladder_xz(f, a24(f, a), (x, f.one()), &cof);
            if !f.is_zero(k.1) {
                break k;
            }
        };
        let ms = multiples(f, a24(f, a), k, ((ell - 1) / 2) as usize);
        let sv = SqrtVeluMont::new(f, a, k, ell).unwrap();
        assert_eq!(sv.codomain(f), velu_codomain(f, a, &ms, ell), "l = {ell}");
        let u = f.random(&mut rng);
        assert_eq!(sv.eval(f, u), isog_x(f, &ms, u));
    }
}
