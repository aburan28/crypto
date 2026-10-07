use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::path::graph::*;
use isogeny_algos::path::{couveignes, delfs_galbraith, galbraith, ghs, volcano};
use isogeny_algos::testdata::*;

fn check_chain(fp: &Zp, cache: &PhiCache, e: &Curve<u64>, path: &Path<u64>, rng: &mut Rng) {
    assert!(verify_path(fp, cache, path));
    assert_eq!(*path.js.first().unwrap(), jinv(fp, e));
    if let Some(chain) = explicit_chain(fp, cache, e, path) {
        let last = chain.last().map(|c| c.cod).unwrap_or(*e);
        assert_eq!(jinv(fp, &last), *path.js.last().unwrap());
        for _ in 0..3 {
            let p = e.random_point(fp, rng);
            let q = apply_chain(fp, &chain, &p);
            assert!(on_curve(fp, &last, &q));
        }
    }
}

#[test]
fn galbraith_ghs_find_paths() {
    let mut rng = Rng::new(11);
    let p = next_prime(1u64 << 24);
    let fp = Zp::new(p);
    let ells = [3usize, 5, 7];
    let cache = PhiCache::new(p, &ells);
    // heights zero for all l: trace with no l^2 | t^2-4p
    let (e1, t) = curve_with_trace(&fp, &mut rng, |t| {
        let d = (t * t - 4 * p as i64).abs();
        ells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
    });
    let j1 = jinv(&fp, &e1);
    let walk = random_walk(&fp, &cache, &ells, j1, 40, &mut rng);
    let j2 = *walk.js.last().unwrap();
    let (r, _) = galbraith::galbraith(&fp, &cache, &ells, j1, j2, 1_000_000, &mut rng);
    check_chain(&fp, &cache, &e1, &r.expect("galbraith failed"), &mut rng);
    let (r, _) = ghs::ghs(&fp, &cache, &ells, j1, j2, 1_000_000, &mut rng);
    check_chain(&fp, &cache, &e1, &r.expect("ghs failed"), &mut rng);
    let (r, _) = ghs::ghs_volcano(&fp, &cache, &ells, p, t, j1, j2, 1_000_000, &mut rng);
    check_chain(&fp, &cache, &e1, &r.expect("ghs_volcano failed"), &mut rng);
}

#[test]
fn kohel_volcano_and_ghs_volcano_with_conductor() {
    let mut rng = Rng::new(12);
    let p = next_prime(1u64 << 17);
    let fp = Zp::new(p);
    let ells = [3usize, 5];
    let cache = PhiCache::new(p, &ells);
    let (e1, t) = curve_with_trace(&fp, &mut rng, |t| {
        let d = (t * t - 4 * p as i64).abs();
        d % 9 == 0 && d % 25 != 0
    });
    let h = volcano::volcano_height(p, t, 3);
    assert!(h >= 1);
    let j1 = jinv(&fp, &e1);
    let walk = random_walk(&fp, &cache, &[3], j1, 25, &mut rng);
    let j2 = *walk.js.last().unwrap();
    let r = volcano::kohel_volcano_path(&fp, &cache, 3, h, j1, j2, &mut rng).expect("kohel failed");
    check_chain(&fp, &cache, &e1, &r, &mut rng);
    let (r, _) = ghs::ghs_volcano(&fp, &cache, &ells, p, t, j1, j2, 1_000_000, &mut rng);
    check_chain(&fp, &cache, &e1, &r.expect("ghs_volcano failed"), &mut rng);
}

#[test]
fn couveignes_recovers_exponents() {
    let mut rng = Rng::new(13);
    let p = next_prime(1u64 << 22);
    let fp = Zp::new(p);
    let ells = [3usize, 5, 7, 11, 13, 17, 19, 23];
    let cache = PhiCache::new(p, &ells);
    let (e1, _t) = curve_with_trace(&fp, &mut rng, |t| {
        let d = (t * t - 4 * p as i64).abs();
        ells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
    });
    let j1 = jinv(&fp, &e1);
    let mut act = couveignes::Action {
        fp: &fp,
        cache: &cache,
        plus_class: Default::default(),
    };
    let primes = couveignes::select_primes(&mut act, j1, &mut rng);
    assert!(primes.len() >= 2, "need >=2 split primes, got {:?}", primes);
    let primes = &primes[..primes.len().min(3)];
    let secret: Vec<i64> = primes.iter().map(|_| rng.below(7) as i64 - 3).collect();
    let pa = act.apply(j1, primes, &secret, &mut rng).expect("apply");
    let j2 = *pa.js.last().unwrap();
    let (r, _) = couveignes::couveignes(&act, primes, 3, j1, j2, &mut rng);
    let (e, path) = r.expect("couveignes failed");
    assert!(verify_path(&fp, &cache, &path));
    assert_eq!(*path.js.last().unwrap(), j2);
    // exponent vector must reproduce j2
    let pb = act.apply(j1, primes, &e, &mut rng).unwrap();
    assert_eq!(*pb.js.last().unwrap(), j2);
}

#[test]
fn delfs_galbraith_supersingular() {
    let mut rng = Rng::new(14);
    let mut p = next_prime(1u64 << 22);
    while p % 8 != 7 {
        p = next_prime(p + 1);
    }
    let fp = Zp::new(p);
    let f2 = Zp2::new(p);
    let cache = PhiCache::new(p, &[2, 3]);
    let start = (1728u64, 0u64);
    let w1 = random_walk(&f2, &cache, &[2], start, 60, &mut rng);
    let w2 = random_walk(&f2, &cache, &[2], start, 60, &mut rng);
    let (j1, j2) = (*w1.js.last().unwrap(), *w2.js.last().unwrap());
    let (r, st) =
        delfs_galbraith::delfs_galbraith(&f2, &fp, &cache, j1, j2, 1_000_000, 1_000_000, &mut rng);
    let path = r.expect("delfs-galbraith failed");
    assert_eq!(path.js[0], j1);
    assert_eq!(*path.js.last().unwrap(), j2);
    assert!(verify_path(&f2, &cache, &path));
    eprintln!(
        "dg: path len {} walk {} bfs {}",
        path.len(),
        st.walk_steps,
        st.bfs_nodes
    );
}

#[test]
fn galbraith_stolbunov_weighted_walk() {
    use isogeny_algos::path::ghs::galbraith_stolbunov;
    let mut rng = Rng::new(15);
    let p = next_prime(1u64 << 24);
    let fp = Zp::new(p);
    let ells = [3usize, 5, 7, 11];
    let cache = PhiCache::new(p, &ells);
    let (e1, _t) = curve_with_trace(&fp, &mut rng, |t| {
        let d = (t * t - 4 * p as i64).abs();
        ells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
    });
    let j1 = jinv(&fp, &e1);
    let walk = random_walk(&fp, &cache, &ells, j1, 60, &mut rng);
    let j2 = *walk.js.last().unwrap();
    for weights in [[1u32, 1, 1, 1], [8, 4, 2, 1], [1, 1, 1, 8]] {
        let (r, _) = galbraith_stolbunov(&fp, &cache, &ells, &weights, j1, j2, 1_000_000, &mut rng);
        check_chain(
            &fp,
            &cache,
            &e1,
            &r.expect("galbraith-stolbunov failed"),
            &mut rng,
        );
    }
}
