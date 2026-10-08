use isogeny_algos::curve::*;
use isogeny_algos::field::{Field, Zp2};
use isogeny_algos::testdata::kani_instance;
use isogeny_algos::theta::split_j;
use isogeny_algos::theta_g::*;

/// The general-dimension code with g = 2 (algorithmic change of basis) reproduces the
/// dimension-2 Kani results: split at E0 x X with the Velu-computed j-invariants.
#[test]
fn theta_g2_matches_dimension_two_kani() {
    for &(a, b) in &[(8u32, 4u32), (10, 6), (12, 7), (16, 10)] {
        let inst = kani_instance(a, b, 21_000 + a as u64);
        let f: &Zp2 = &inst.f;
        let iota = (0u64, 1u64);
        assert_eq!(f.sq(iota), f.neg(f.one()));
        let k: Vec<Vec<Pt<(u64, u64)>>> = inst.k.iter().map(|(x, y)| vec![*x, *y]).collect();
        let res = chain_g(f, &[inst.c, inst.e], &k, a, &[], iota).expect("chain");
        let last = res.nulls.last().unwrap();
        let arr = [last[0], last[1], last[2], last[3]];
        let (j1, j2) = split_j(f, &arr).expect("split");
        let want = (jinv(f, &inst.e0), jinv(f, &inst.x));
        assert!((j1, j2) == want || (j1, j2) == (want.1, want.0), "a = {a}");
        assert_eq!(res.gluing_steps, 1);
        let kt: Vec<Vec<Pt<(u64, u64)>>> = inst.k_twisted.iter().map(|(x, y)| vec![*x, *y]).collect();
        let res2 = chain_g(f, &[inst.c, inst.e], &kt, a, &[], iota).expect("chain");
        assert_eq!(even_zero_count(f, res2.nulls.last().unwrap()), 0);
        assert_eq!(even_zero_count(f, last), 1);
    }
}

/// Dimension-4 Kani embedding: the chain splits (vanishing even theta constants) for the true
/// kernel and not for the twisted one.
#[test]
fn theta_g4_kani_embedding_splits() {
    use isogeny_algos::bigint::Big;
    use isogeny_algos::field::{is_prime, Rng};
    use isogeny_algos::testdata::kani4_instance;
    for &(n, b) in &[(12u32, 7u32), (16, 9), (20, 12)] {
        let base = (1u64 << (n + 2)) * 3u64.pow(b);
        let p = (1..).map(|c| base * c - 1).find(|&p| is_prime(p)).unwrap();
        let f = Zp2::new(p);
        let iota = (0u64, 1u64);
        let mut rng = Rng::new(400 + n as u64);
        let inst = kani4_instance(&f, &Big::from_u64(p + 1), n, b, iota, &mut rng);
        let res = chain_g(&f, &inst.curves, &inst.k, n, &[], iota).expect("chain");
        let z = split_score(&f, &res, &[]);
        assert!(z > 0, "n = {n}: no split");
        let res2 = chain_g(&f, &inst.curves, &inst.k_twisted, n, &[], iota).expect("chain");
        assert_eq!(split_score(&f, &res2, &[]), 0, "n = {n}");
        eprintln!("n = {n}, b = {b}: zeros {z}, gluing steps {}", res.gluing_steps);
    }
}
