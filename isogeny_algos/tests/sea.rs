use isogeny_algos::bigint::Big;
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::sea::*;
use isogeny_algos::fpm::FpM;
use std::collections::HashMap;

#[test]
fn atkin_candidate_sets_partition() {
    // every t mod l with a non-zero ratio falls in exactly one order class
    for l in [3u64, 5, 7, 11, 13] {
        for q in 1..l {
            let mut total = 0;
            for r in 1..=(l as usize + 1) {
                total += atkin_candidates(l, q, r).len();
            }
            assert!(total <= l as usize);
        }
    }
}

#[test]
fn sea_matches_bsgs_at_40_and_61_bits() {
    let mut rng = Rng::new(1100);
    for bits in [40u32, 61] {
        let fp = Zp::new(next_prime((1u64 << (bits - 1)) + 12345));
        let mut phis = HashMap::new();
        for _ in 0..4 {
            let e = loop {
                let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
                let j = jinv(&fp, &e);
                if is_smooth(&fp, &e) && j != 0 && j != 1728 {
                    break e;
                }
            };
            let (n, st) = sea(&fp, &e, 43, &mut phis, &mut rng).expect("SEA");
            let reference = order(&fp, &e, &mut rng);
            assert_eq!(n, Big::from_u64(reference), "bits = {bits}, stats {st:?}");
        }
    }
}

#[test]
fn sea_at_128_bits() {
    let mut rng = Rng::new(1101);
    // 2^127 - 1
    let f = FpM::<2>::from_dec("170141183460469231731687303715884105727");
    let mut phis = HashMap::new();
    let e = loop {
        let e = Curve::new(f.random(&mut rng), f.random(&mut rng));
        let j = jinv(&f, &e);
        if is_smooth(&f, &e) && j != f.zero() && j != f.from_u64(1728) {
            break e;
        }
    };
    let (n, st) = sea(&f, &e, 61, &mut phis, &mut rng).expect("SEA");
    for _ in 0..5 {
        let p = random_point_f(&f, &e, &mut rng);
        assert_eq!(pmul_big(&f, &e, &p, &n), Pt::Inf, "stats {st:?}");
    }
    eprintln!("128-bit SEA: elkies {:?}, atkin {:?}, candidates {}", st.elkies, st.atkin, st.candidates);
}
