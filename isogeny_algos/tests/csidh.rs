use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::kernel::montgomery::j_invariant;
use isogeny_algos::path::csidh::*;

#[test]
fn class_number_known_values() {
    for (d, h) in [
        (-3i64, 1u64),
        (-4, 1),
        (-15, 2),
        (-20, 2),
        (-23, 3),
        (-47, 5),
        (-71, 7),
        (-84, 4),
        (-163, 1),
        (-296, 10),
        (-4 * 23, 3),
    ] {
        assert_eq!(class_number(d), h, "h({d})");
    }
}

#[test]
fn csidh_toy_action_is_a_group_action() {
    let cs = Csidh::with_n_primes(6).expect("params");
    assert_eq!(cs.p(), 1021019);
    let mut rng = Rng::new(71);
    let fp = cs.fp;
    // every positive step is an l-isogeny: Phi_l(j, j') = 0
    let phis: Vec<Phi<Zp>> = [3usize, 5, 7, 11, 13]
        .iter()
        .map(|&l| Phi::compute(&fp, l))
        .collect();
    let mut a = 0u64;
    for step in 0..20 {
        let idx = step % 5;
        let a2 = cs.step(a, idx, step % 3 != 0, &mut rng);
        assert_eq!(
            phis[idx].eval(&fp, j_invariant(&fp, a), j_invariant(&fp, a2)),
            0,
            "step {step}"
        );
        a = a2;
    }
    // commutativity and inverse
    for _ in 0..5 {
        let e1: Vec<i32> = (0..6).map(|_| rng.below(7) as i32 - 3).collect();
        let e2: Vec<i32> = (0..6).map(|_| rng.below(7) as i32 - 3).collect();
        let ab = cs.action(cs.action(0, &e1, &mut rng), &e2, &mut rng);
        let ba = cs.action(cs.action(0, &e2, &mut rng), &e1, &mut rng);
        assert_eq!(ab, ba, "class group action must commute");
        let sum: Vec<i32> = e1.iter().zip(&e2).map(|(x, y)| x + y).collect();
        assert_eq!(ab, cs.action(0, &sum, &mut rng));
        let neg: Vec<i32> = e1.iter().map(|x| -x).collect();
        assert_eq!(cs.action(cs.action(0, &e1, &mut rng), &neg, &mut rng), 0);
    }
}

#[test]
fn ideal_orders_divide_class_number() {
    let cs = Csidh::with_n_primes(6).unwrap();
    let mut rng = Rng::new(72);
    let h = class_number(-4 * cs.p() as i64);
    for idx in 0..3 {
        let ord = cs
            .ideal_order(idx, &mut rng, h as usize + 1)
            .expect("cycle");
        assert_eq!(
            h % ord as u64,
            0,
            "order {ord} of l_{idx} must divide h = {h}"
        );
    }
}

#[test]
fn mitm_recovers_exponents() {
    let cs = Csidh::with_n_primes(6).unwrap();
    let mut rng = Rng::new(73);
    for _ in 0..3 {
        let e: Vec<i32> = (0..6).map(|_| rng.below(5) as i32 - 2).collect();
        let target = cs.action(0, &e, &mut rng);
        let (found, nodes) = cs.mitm(0, target, 2, &mut rng).expect("collision");
        assert!(found.iter().all(|&x| x.abs() <= 2));
        assert_eq!(cs.action(0, &found, &mut rng), target);
        assert_eq!(nodes, 2 * 3usize.pow(6));
    }
}

#[test]
fn batched_action_equals_stepwise() {
    for n in [6usize, 8, 10] {
        let cs = Csidh::with_n_primes(n).unwrap();
        let mut rng = Rng::new(74 + n as u64);
        for _ in 0..4 {
            let e: Vec<i32> = (0..n).map(|_| rng.below(9) as i32 - 4).collect();
            let a1 = cs.action(0, &e, &mut rng);
            let a2 = cs.action_batched(0, &e, &mut rng);
            assert_eq!(a1, a2, "n = {n}, e = {e:?}");
        }
    }
}
