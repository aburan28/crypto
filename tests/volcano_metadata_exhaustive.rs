//! Bounded correctness census, with independent arithmetic oracles.
use crypto_lib::isogeny::{cm::*, volcano::*, SmallCurve};
use std::collections::{BTreeMap, BTreeSet};

const PRIMES: &[u64] = &[5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43];
fn curve(p: u64, a: u64, b: u64) -> SmallCurve {
    SmallCurve {
        name: "exhaustive",
        p,
        a,
        b,
    }
}
fn nonsingular(c: &SmallCurve) -> bool {
    (4 * c.a.pow(3) + 27 * c.b.pow(2)) % c.p != 0
}
// Counts square roots by enumerating y, without Legendre symbols or EC arithmetic.
fn count(c: &SmallCurve) -> i64 {
    let mut roots = vec![0i64; c.p as usize];
    for y in 0..c.p {
        roots[(y * y % c.p) as usize] += 1;
    }
    1 + (0..c.p)
        .map(|x| roots[((x * x % c.p * x + c.a * x + c.b) % c.p) as usize])
        .sum::<i64>()
}
fn squarefree(n: i64) -> bool {
    let n = n.abs();
    !(2..).take_while(|d| d * d <= n).any(|d| n % (d * d) == 0)
}
fn fundamental(d: i64) -> bool {
    d < 0
        && ((d.rem_euclid(4) == 1 && squarefree(d))
            || (d % 4 == 0 && [2, 3].contains(&(d / 4).rem_euclid(4)) && squarefree(d / 4)))
}
fn v(mut n: i64, ell: u64) -> u32 {
    let mut out = 0;
    while n % ell as i64 == 0 {
        n /= ell as i64;
        out += 1;
    }
    out
}

#[test]
fn exhaustive_order_discriminants() {
    let mut tested = 0;
    for n in -100_000i64..0 {
        if ![0, 1].contains(&n.rem_euclid(4)) {
            continue;
        }
        let (d, f) = fundamental_discriminant_and_conductor(n);
        assert!(fundamental(d), "n={n}, D={d}, f={f}");
        assert!(f > 0);
        assert_eq!(d * f * f, n);
        tested += 1;
    }
    println!("valid negative discriminants checked: {tested}");
}

#[test]
fn exhaustive_small_prime_metadata() {
    let (mut curves, mut ss, mut unknown, mut positions) = (0, 0, 0, 0);
    let mut exact_bsgs = 0;
    for &p in PRIMES {
        for a in 0..p {
            for b in 0..p {
                let c = curve(p, a, b);
                if !nonsingular(&c) {
                    continue;
                }
                let m = cm_discriminant(&c);
                let n = count(&c);
                let t = p as i64 + 1 - n;
                if let Some(candidate) = frobenius_trace_bsgs(&c) {
                    assert_eq!(candidate, t, "BSGS {c:?}");
                    exact_bsgs += 1;
                }
                assert_eq!((m.order, m.trace), (n, t), "{c:?}");
                assert!(t * t <= 4 * p as i64);
                assert_eq!(m.frobenius_disc, t * t - 4 * p as i64);
                assert!(fundamental(m.fundamental_disc));
                assert_eq!(
                    m.fundamental_disc * m.frobenius_order_conductor.pow(2),
                    m.frobenius_disc
                );
                assert_eq!(m.anomalous, n == p as i64);
                assert_eq!(m.supersingular, t == 0);
                if t == 0 {
                    ss += 1;
                    assert_eq!(m.endomorphism_conductor, None);
                    assert_eq!(m.endomorphism_disc, None);
                    assert_eq!(m.endomorphism_evidence, None);
                } else if let Some(f) = m.endomorphism_conductor {
                    assert_eq!(f, 1);
                    assert_eq!(m.endomorphism_disc, Some(m.fundamental_disc));
                    assert!(m.frobenius_order_conductor == 1 || a == 0 || b == 0);
                    if a == 0 {
                        assert_eq!(m.fundamental_disc, -3);
                    }
                    if b == 0 {
                        assert_eq!(m.fundamental_disc, -4);
                    }
                    assert!(verify_cm(&c, m.fundamental_disc));
                    assert!(!verify_cm(&c, m.fundamental_disc * 4));
                } else {
                    unknown += 1;
                    assert!(m.frobenius_order_conductor > 1 && a != 0 && b != 0);
                    assert!(!verify_cm(&c, m.frobenius_disc));
                }
                for ell in [2, 3, 5, 7].into_iter().filter(|&ell| ell != p) {
                    let pos = position(&c, ell);
                    let expected_max = (t != 0).then(|| v(m.frobenius_order_conductor, ell));
                    let expected_depth = if t == 0 {
                        None
                    } else if m.endomorphism_conductor == Some(1)
                        || m.frobenius_order_conductor % ell as i64 != 0
                    {
                        Some(0)
                    } else {
                        None
                    };
                    assert_eq!(pos.max_depth, expected_max, "{c:?}, ell={ell}");
                    assert_eq!(pos.depth, expected_depth);
                    assert_eq!(pos.on_crater, expected_depth.map(|d| d == 0));
                    for cap in 0..=3 {
                        assert_eq!(
                            volcano_depth(&c, ell, cap),
                            expected_max
                                .filter(|&d| d as usize <= cap)
                                .map(|d| d as usize)
                        );
                    }
                    positions += 1;
                }
                curves += 1;
            }
        }
    }
    println!("metadata: curves={curves}, supersingular={ss}, unknown_ordinary={unknown}, positions={positions}, exact_bsgs_certificates={exact_bsgs}");
}

#[test]
fn exhaustive_two_kernel_maps() {
    let (mut curves, mut kernels, mut images) = (0, 0, 0);
    for &p in PRIMES.iter().filter(|&&p| p <= 19) {
        for a in 0..p {
            for b in 0..p {
                let c = curve(p, a, b);
                if !nonsingular(&c) {
                    continue;
                }
                let expected: BTreeSet<_> = (0..p)
                    .filter(|&x| (x * x * x + a * x + b) % p == 0)
                    .collect();
                let isos = neighbors_ell(&c, 2);
                assert_eq!(isos.len(), expected.len());
                assert_eq!(
                    isos.iter()
                        .map(|i| i.kernel_half[0].0)
                        .collect::<BTreeSet<_>>(),
                    expected
                );
                for iso in isos {
                    assert!(nonsingular(&iso.codomain));
                    assert_eq!(count(&iso.codomain), count(&c));
                    for x in 0..p {
                        for y in 0..p {
                            if y * y % p != (x * x * x + a * x + b) % p {
                                continue;
                            }
                            match iso.evaluate(x, y) {
                                None => assert!(
                                    y == 0
                                        && expected.contains(&x)
                                        && iso.kernel_half.contains(&(x, y))
                                ),
                                Some((u, w)) => {
                                    assert_eq!(
                                        w * w % p,
                                        (u * u * u + iso.codomain.a * u + iso.codomain.b) % p,
                                        "{c:?}, ({x},{y})"
                                    );
                                    images += 1;
                                }
                            }
                        }
                    }
                    kernels += 1;
                }
                curves += 1;
            }
        }
    }
    println!("2-kernel census: curves={curves}, kernels={kernels}, nonkernel_images={images}");
}

#[test]
fn exhaustive_capped_neighborhood_edges() {
    let mut maps = 0;
    for &p in PRIMES.iter().filter(|&&p| p <= 19) {
        for a in 0..p {
            for b in 0..p {
                let c = curve(p, a, b);
                if !nonsingular(&c) || count(&c) == p as i64 + 1 {
                    continue;
                }
                for cap in 1..=4 {
                    for depth in 0..=2 {
                        let map = map_volcano(&c, 2, depth, cap);
                        let mut retained = BTreeMap::new();
                        for level in &map.levels {
                            assert!(level.bfs_distance <= depth);
                            for (&j, &(p, a, b)) in level.j_invariants.iter().zip(&level.curves) {
                                assert!(retained
                                    .insert(j, (curve(p, a, b), level.bfs_distance))
                                    .is_none());
                            }
                        }
                        assert!(retained.len() <= cap);
                        // Every retained vertex inside the depth bound must retain all
                        // kernel edges to retained vertices, even after reaching the cap.
                        let mut expected = vec![];
                        for (&j, (source, d)) in &retained {
                            if *d >= depth {
                                continue;
                            }
                            for iso in neighbors_ell(source, 2) {
                                let target = j_invariant(&iso.codomain);
                                if retained.contains_key(&target) {
                                    expected.push((j, target, iso.kernel_half));
                                }
                            }
                        }
                        let mut actual: Vec<_> = map
                            .edges
                            .iter()
                            .map(|e| (e.source_j, e.target_j, e.kernel_half.clone()))
                            .collect();
                        expected.sort();
                        actual.sort();
                        assert_eq!(actual, expected, "{c:?}, cap={cap}, depth={depth}");
                        maps += 1;
                    }
                }
            }
        }
    }
    println!("capped neighborhoods checked: {maps}");
}

#[test]
fn point_count_branch_boundary_against_square_histogram() {
    let mut checked = 0;
    let mut bsgs_certificates = 0;
    // Fixed inputs on either side of the 2^14 branch threshold.
    for p in [16381, 16411, 65537] {
        for (a, b) in [(0, 1), (1, 0), (1, 1), (2, 3), (7, 11), (19, 23)] {
            let c = curve(p, a, b);
            if !nonsingular(&c) {
                continue;
            }
            let expected = p as i64 + 1 - count(&c);
            assert_eq!(frobenius_trace(&c), expected, "{c:?}");
            if let Some(t) = frobenius_trace_bsgs(&c) {
                assert_eq!(t, expected, "{c:?}");
                bsgs_certificates += 1;
            }
            checked += 1;
        }
    }
    assert!(bsgs_certificates > 0);
    println!("branch-boundary point counts: {checked}; BSGS certificates: {bsgs_certificates}");
}
