use isogeny_algos::field::{Field, Rng};
use isogeny_algos::gf3n::GF3n;

#[test]
fn gf3n_field_axioms() {
    for n in 1..=12u32 {
        let f = GF3n::new(n);
        let q = 3u64.pow(n);
        assert_eq!(f.size() as u64, q);
        let mut rng = Rng::new(n as u64 * 7 + 1);
        // characteristic 3
        assert_eq!(f.add(f.add(f.one(), f.one()), f.one()), f.zero());
        let mut nonsquares = 0;
        for _ in 0..60 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            let c = f.random(&mut rng);
            // ring axioms
            assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
            assert_eq!(f.mul(a, b), f.mul(b, a));
            assert_eq!(f.sub(f.add(a, b), b), a);
            if a != 0 {
                assert_eq!(f.mul(a, f.inv(a)), f.one());
                // Fermat: a^(q-1) = 1
                let mut p = f.one();
                let mut base = a;
                let mut e = q - 1;
                while e > 0 {
                    if e & 1 == 1 {
                        p = f.mul(p, base);
                    }
                    base = f.mul(base, base);
                    e >>= 1;
                }
                assert_eq!(p, f.one());
            }
            // sqrt
            let s = f.sq(a);
            let r = f.sqrt(s).expect("square has a root");
            assert!(r == a || r == f.neg(a) || f.sq(r) == s);
            if a != 0 && f.sqrt(a).is_none() {
                nonsquares += 1;
            }
        }
        if n >= 1 && q > 3 {
            assert!(nonsquares > 0, "n = {n}: no non-squares seen");
        }
    }
}
