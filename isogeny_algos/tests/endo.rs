use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::path::endo::*;
use isogeny_algos::path::graph::*;
use isogeny_algos::path::volcano::*;
use isogeny_algos::testdata::*;

#[test]
fn conductor_drops_by_l_along_ascent() {
    let mut rng = Rng::new(81);
    let mut p = next_prime(1u64 << 17);
    while ![1u64, 4, 7].contains(&(p % 9)) {
        p = next_prime(p + 1);
    }
    let fp = Zp::new(p);
    let cache = PhiCache::new(&Zp::new(p), &[3, 5, 7]);
    let mut verified_vertical = 0;
    for _ in 0..6 {
        let (e, t) = curve_with_trace(&fp, &mut rng, |t| {
            let d = (t * t - 4 * p as i64).abs();
            t % 2 != 0 && d % 27 == 0 && d % 25 != 0 && d % 49 != 0
        });
        let j = jinv(&fp, &e);
        let end = endomorphism_ring(&fp, &cache, p, t, j, &mut rng).expect("end ring");
        let (_, h, level) = end.per_prime[0];
        assert_eq!(end.per_prime[0].0, 3);
        assert!(level <= h);
        assert_eq!(
            end.disc,
            end.fundamental_disc * (end.conductor * end.conductor) as i128
        );
        // ascend: the conductor loses exactly one factor of 3 per step
        let up = ascend_to_crater(&fp, &cache, 3, h, j, &mut rng);
        assert_eq!(up.len() as u32, level, "ascent length must equal the level");
        let mut prev_level = level;
        for jj in up.js.iter().skip(1) {
            let e2 = endomorphism_ring(&fp, &cache, p, t, *jj, &mut rng).unwrap();
            let l2 = e2.per_prime[0].2;
            assert_eq!(l2, prev_level - 1);
            prev_level = l2;
            verified_vertical += 1;
        }
        // the crater has maximal order at 3
        assert_eq!(prev_level, 0);
    }
    assert!(verified_vertical > 0, "need at least one non-crater curve");
}

#[test]
fn conductor_one_when_discriminant_is_squarefree() {
    let mut rng = Rng::new(82);
    let p = next_prime(1u64 << 20);
    let fp = Zp::new(p);
    let cache = PhiCache::new(&Zp::new(p), &[3, 5]);
    let (e, t) = curve_with_trace(&fp, &mut rng, |t| {
        let d = (t * t - 4 * p as i64).abs();
        d % 9 != 0 && d % 25 != 0 && d % 4 != 0 && t % 2 != 0
    });
    let r = endomorphism_ring(&fp, &cache, p, t, jinv(&fp, &e), &mut rng).unwrap();
    assert_eq!(r.conductor, 1);
    assert_eq!(
        r.disc,
        t as i128 * t as i128
            - 4 * p as i128 / ((r.conductor_frobenius * r.conductor_frobenius) as i128)
    );
}
