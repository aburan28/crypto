//! Standalone exact-equality acceptance run; no timing or speedup claim.
//! Build against the release isogeny_algos rlib, as recorded in commands.sh.
use isogeny_algos::field::Rng;
use isogeny_algos::find::{modpoly::integer_coeffs, sutherland::cm::CmPlan};
fn main() {
    for ell in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31] {
        let plan = CmPlan::new(ell).unwrap();
        let primes = plan.crt_primes().unwrap();
        eprintln!(
            "checking ell={ell} D={} CRT primes={}",
            plan.discriminant(),
            primes.len()
        );
        let actual = plan.phi_integer(&mut Rng::new(0x5107)).unwrap();
        let expected = integer_coeffs(ell);
        assert_eq!(
            actual, expected,
            "integer coefficient mismatch at ell={ell}"
        );
        let (h_top, h_floor) = plan.class_numbers();
        println!("{{\"ell\":{ell},\"D\":{},\"surface_classes\":{h_top},\"floor_classes\":{h_floor},\"crt_primes\":{},\"largest_prime\":{},\"matrix_entries\":{},\"seed\":20743,\"exact_integer_equality\":true}}",
            plan.discriminant(),primes.len(),primes.last().unwrap().0,(ell+2)*(ell+2));
    }
}
