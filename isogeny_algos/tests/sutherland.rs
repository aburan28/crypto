use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::{integer_coeffs, Phi};
use isogeny_algos::find::sutherland::{phi_crt, phi_mod_p};

/// Sutherland's isogeny-graph Phi_l mod p equals the q-expansion (Hecke) Phi_l mod p exactly.
/// The two methods share nothing: one walks the l-isogeny graph (Velu codomains), the other
/// solves q-series relations. l = 3, 5, 7 start from X_1(l) samples, l = 11 from uniformly
/// random curves behind the l^2 filter; all use the prime-to-l walk with transported torsion.
#[test]
fn sutherland_phi_mod_p_matches_qexpansion() {
    for &ell in &[3usize, 5, 7, 11] {
        let mut found = 0;
        let mut cand = (1u64 << 16) + 1;
        while found < 2 {
            if cand % ell as u64 == 1 && is_prime(cand) {
                let fp = Zp::new(cand);
                let mut rng = Rng::new(0x5107 ^ cand);
                let c = phi_mod_p(cand, ell, &mut rng).expect("enough full-structure curves");
                assert_eq!(
                    c,
                    Phi::compute(&fp, ell).c,
                    "l={ell} p={cand}: Sutherland != Hecke mod p"
                );
                found += 1;
            }
            cand += 2;
        }
    }
}

/// Sutherland's CRT reconstruction gives the exact integer Phi_l (equal to the crate's
/// q-expansion + CRT `integer_coeffs`).
#[test]
fn sutherland_phi_crt_matches_integer_phi() {
    for &ell in &[3usize, 5] {
        let c = phi_crt(ell).expect("phi_crt");
        assert_eq!(
            c,
            integer_coeffs(ell),
            "l={ell}: Sutherland CRT != integer Phi_l"
        );
    }
}
