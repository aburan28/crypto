mod common;
use common::*;
use eigen_seta::cm::{mul_vartime, torsion_basis};
use eigen_seta::pairing::{dlog_3k, full_order_3k, weil_sq};
use eigen_seta::params::{ParamSet, Toy16};
use num_bigint::BigUint;
use rand_core::{OsRng, RngCore};

#[test]
fn weil_bilinear_and_dlog() {
    let mut rng = OsRng;
    let ed = ed_curve();
    let k = Toy16::F3;
    let n = BigUint::from(3u32).pow(k as u32);
    let (p, q) = torsion_basis(&ed, 3, k, &mut rng);
    let e = weil_sq(&ed, &p, &q, &n);
    assert!(
        bool::from(full_order_3k(&e, k)),
        "pairing of a basis not of full order"
    );
    assert_eq!(cube_n(&e, k), F::one(), "pairing not in mu_n");
    for _ in 0..5 {
        let a = BigUint::from(rng.next_u64()) % &n;
        let b = BigUint::from(rng.next_u64()) % &n;
        let ap = mul_vartime(&ed, &p, &a);
        let bq = mul_vartime(&ed, &q, &b);
        let lhs = weil_sq(&ed, &ap, &bq, &n);
        let rhs = e.pow_public(&(&a * &b % &n).to_u64_digits());
        assert_eq!(lhs, rhs, "bilinearity");
        let s = dlog_3k(&e, &lhs, k);
        assert_eq!(s.to_biguint() % &n, &a * &b % &n, "dlog");
        // alternating
        let pq = ed.add(&ap, &bq);
        let e1 = weil_sq(&ed, &pq, &q, &n);
        let e2 = weil_sq(&ed, &ap, &q, &n);
        assert_eq!(e1, e2, "e(aP + bQ, Q) = e(aP, Q)");
    }
}

fn cube_n(x: &F, k: usize) -> F {
    eigen_seta::pairing::cube_iter(x, k)
}
