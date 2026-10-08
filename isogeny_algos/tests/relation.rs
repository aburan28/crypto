use isogeny_algos::field::*;
use isogeny_algos::path::csidh::*;
use isogeny_algos::path::relation::*;

#[test]
fn form_group_and_class_number() {
    for p in [1019i128, 1031, 2003, 10007, 100003] {
        if p % 4 != 3 {
            continue;
        }
        let d = -4 * p;
        let h = class_number(d as i64) as i128;
        // group laws
        let f = Form { a: 3, b: 2, c: (4 - d) / 12 };
        if (4 - d) % 12 == 0 {
            let f = f.reduce();
            assert_eq!(f.compose(&f.inverse()), Form::identity(d));
            assert_eq!(f.pow(h), Form::identity(d), "f^h = 1, p = {p}");
        }
        // BSGS exponent around the analytic estimate finds a multiple of every order
        let est = class_number_estimate(d, 100_000);
        assert!((est - h as f64).abs() / (h as f64) < 0.1, "estimate {est} vs h {h}");
    }
}

#[test]
fn relations_act_trivially_on_toy_csidh() {
    let mut rng = Rng::new(1300);
    let mut tested = 0;
    for n in [6usize, 8, 10] {
        let cs = Csidh::with_n_primes(n).unwrap();
        let p = cs.p() as i128;
        let Some(rl) = relation_lattice(p, &cs.primes) else { continue };
        if n == 6 {
            assert_eq!(rl.h, class_number(-4 * p as i64) as i128);
        }
        // every reduced basis vector is a relation: [prod l_i^{v_i}] E_0 = E_0
        for v in &rl.basis {
            let e: Vec<i32> = v.iter().map(|&x| x as i32).collect();
            assert_eq!(cs.action_fast(0, &e, &mut rng), 0, "n = {n}, relation {v:?}");
        }
        // reduction mod L does not change the action
        for _ in 0..5 {
            let e: Vec<i64> = (0..n).map(|_| rng.below(41) as i64 - 20).collect();
            let r = rl.reduce(&e);
            let l1 = |v: &[i64]| v.iter().map(|x| x.abs()).sum::<i64>();
            assert!(l1(&r) <= l1(&e));
            let ei: Vec<i32> = e.iter().map(|&x| x as i32).collect();
            let ri: Vec<i32> = r.iter().map(|&x| x as i32).collect();
            assert_eq!(cs.action_fast(0, &ei, &mut rng), cs.action_fast(0, &ri, &mut rng));
        }
        // classes [l_1]^a for random a
        for _ in 0..3 {
            let a = rng.below(rl.h as u64) as i128;
            let v = rl.vector_of(a);
            let mut e1 = vec![0i32; n];
            e1[0] = a as i32;
            let vi: Vec<i32> = v.iter().map(|&x| x as i32).collect();
            assert_eq!(cs.action_fast(0, &vi, &mut rng), cs.action(0, &e1, &mut rng));
        }
        tested += 1;
    }
    assert!(tested >= 2, "lattices computed for {tested} parameter sets");
}

#[test]
fn explicit_lattice_handles_non_cyclic_groups() {
    let mut rng = Rng::new(1301);
    for n in [6usize, 8, 10] {
        let cs = Csidh::with_n_primes(n).unwrap();
        let p = cs.p() as i128;
        let rl = relation_lattice_explicit(p, &cs.primes, 20_000_000).expect("explicit lattice");
        if let Some(c) = relation_lattice(p, &cs.primes) {
            assert_eq!(c.h, rl.h, "same class number, n = {n}");
        }
        let est = class_number_estimate(-4 * p, 200_000);
        assert!((rl.h as f64 - est).abs() / est < 0.1, "h = {} vs estimate {est}", rl.h);
        for v in &rl.basis {
            let e: Vec<i32> = v.iter().map(|&x| x as i32).collect();
            assert_eq!(cs.action_fast(0, &e, &mut rng), 0, "n = {n}, relation {v:?}");
        }
    }
}
