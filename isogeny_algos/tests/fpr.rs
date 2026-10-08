use isogeny_algos::bigint::Big;
use isogeny_algos::field::*;
use isogeny_algos::fpr::*;

/// F_{p^r} satisfies the field axioms, Frobenius is additive, and x^(p^r) = x for all x.
#[test]
fn fpr_is_a_field() {
    for &(pp, r) in &[
        (97u64, 2usize),
        (97, 3),
        (101, 4),
        (1_000_003, 2),
        (1_000_003, 3),
        ((1u64 << 40) + 15, 3),
    ] {
        let mut p = pp;
        while !is_prime(p) {
            p += 1;
        }
        let f = FpR::new(p, r);
        let mut rng = Rng::new(p ^ r as u64);
        for _ in 0..40 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            let c = f.random(&mut rng);
            assert_eq!(f.add(a, b), f.add(b, a));
            assert_eq!(f.mul(a, b), f.mul(b, a));
            assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a)), f.one(), "inverse p={p} r={r}");
            }
            // Frobenius x -> x^p is additive
            let fp = |x: ER| f.pow_big(x, &Big::from_u64(p));
            assert_eq!(
                fp(f.add(a, b)),
                f.add(fp(a), fp(b)),
                "frobenius additive p={p} r={r}"
            );
            // every element is fixed by the full Frobenius x -> x^(p^r)
            assert_eq!(f.pow_big(a, &f.q()), a, "a^(p^r)=a p={p} r={r}");
        }
        // there is an element outside the prime field (the extension is non-trivial)
        let gen = {
            let mut g = [0u64; RMAX];
            g[1] = 1;
            g
        };
        assert!(
            f.project(gen).is_none(),
            "generator z should not lie in F_p"
        );
    }
}
