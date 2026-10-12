use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::fpm::FpM;

#[test]
fn hecke_construction_matches_linear_algebra() {
    let fp = Zp::new(next_prime(1u64 << 61));
    for ell in [2usize, 3, 5, 7, 11, 13, 17, 19, 23] {
        let a = Phi::compute_hecke(&fp, ell);
        let b = Phi::compute_linear_algebra(&fp, ell);
        assert_eq!(a.c, b.c, "l = {ell}");
    }
    let f4 = FpM::<4>::from_dec(
        "115792089210356248762697446949407573530086143415290314195533631308867097853951",
    );
    for ell in [3usize, 7, 13] {
        assert_eq!(
            Phi::compute_hecke(&f4, ell).c,
            Phi::compute_linear_algebra(&f4, ell).c
        );
    }
    // small characteristic (p > l + 1): agrees with the integer coefficients reduced mod p
    let small = Zp::new(1259);
    for ell in [3usize, 5, 7] {
        assert_eq!(
            Phi::compute_hecke(&small, ell).c,
            Phi::via_crt(&small, ell).c,
            "p = 1259, l = {ell}"
        );
    }
}
