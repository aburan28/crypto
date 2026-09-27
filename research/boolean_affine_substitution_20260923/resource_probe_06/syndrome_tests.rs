#[test]
fn every_gray_block_matches_direct_quadratic_evaluation() {
    let mut rng = 17731;
    for n in 4..=12 {
        let mut system = Vec::new();
        for _ in 0..32 {
            let mut p = Vec::new();
            for b in 0..n {
                for a in 0..b {
                    if next(&mut rng) % 3 == 0 {
                        p.push((1u64 << a) | (1u64 << b));
                    }
                }
                if next(&mut rng) & 1 != 0 {
                    p.push(1u64 << b);
                }
            }
            if next(&mut rng) & 1 != 0 {
                p.push(0);
            }
            algebra::canonical(&mut p);
            system.push(p);
        }
        let form = SyndromeForm::from_system(&system, n).unwrap();
        let offsets = form.low_offsets();
        let mut cursor = GrayCursor::new(&form);
        for step in 0..1u64 << (n - 4) {
            cursor.advance(&form, step);
            for low in 0..16 {
                let point = ((step ^ (step >> 1)) << 4) | low;
                let expected = system.iter().enumerate().fold(0u32, |syndrome, (e, p)| {
                    syndrome
                        ^ if p.iter().fold(false, |v, &m| v ^ (point & m == m)) {
                            1u32 << e
                        } else {
                            0
                        }
                });
                let mut actual = cursor.constant ^ offsets[low as usize / 4][low as usize % 4];
                for i in 0..4 {
                    if low & (1 << i) != 0 {
                        actual ^= cursor.linear[i];
                    }
                }
                assert_eq!(actual, expected);
                assert_eq!(form.value(point), expected);
            }
            assert_eq!(
                scalar_syndrome_block(cursor.constant, &cursor.linear, &offsets),
                simd_syndrome_block(cursor.constant, &cursor.linear, &offsets)
            );
        }
    }
}

#[test]
fn simd_zero_detection_and_wrapping_checksum_match_every_forced_lane() {
    let mut rng = 57109;
    for _ in 0..256 {
        let constant = next(&mut rng) as u32;
        let linear = std::array::from_fn(|_| next(&mut rng) as u32);
        let offsets: [[u32; 4]; 4] =
            std::array::from_fn(|_| std::array::from_fn(|_| next(&mut rng) as u32));
        for lane in 0..16 {
            let mut forced = offsets;
            let mut value = constant;
            for i in 0..4 {
                if lane & (1 << i) != 0 {
                    value ^= linear[i];
                }
            }
            forced[lane / 4][lane % 4] = value;
            let expected = scalar_syndrome_block(constant, &linear, &forced);
            assert!(expected.1.is_some());
            assert_eq!(simd_syndrome_block(constant, &linear, &forced), expected);
        }
        assert_eq!(
            simd_syndrome_block(constant, &linear, &offsets),
            scalar_syndrome_block(constant, &linear, &offsets)
        );
    }
}

#[test]
fn enumeration_caps_are_censored_and_do_not_overcharge_points() {
    let form = SyndromeForm::from_system(&vec![vec![0]], 6).unwrap();
    for cap in [0, 1, 15, 16, 17, 31, 32, 63, 64, 100] {
        let a = enumerate_syndromes::<false>(&form, cap);
        let b = enumerate_syndromes::<true>(&form, cap);
        assert_eq!(
            (a.model, a.complete, a.points, a.batches, a.checksum),
            (b.model, b.complete, b.points, b.batches, b.checksum)
        );
        assert!(a.points <= cap);
        assert_eq!(a.complete, cap >= 64);
        assert!(a.model.is_none());
    }
}

#[test]
fn tiny_enumeration_and_embedded_leaves_match_independent_truth_tables() {
    let terms: Vec<_> = (0u64..8).filter(|m| m.count_ones() <= 2).collect();
    let embed = |m: u64| {
        (0..3).fold(0, |v, i| {
            v | if m & (1 << i) != 0 {
                1u64 << [1, 7, 15][i]
            } else {
                0
            }
        })
    };
    for chosen in 0..128 {
        let mut p: Vec<_> = terms
            .iter()
            .enumerate()
            .filter(|(i, _)| chosen & (1 << i) != 0)
            .map(|(_, &m)| m)
            .collect();
        algebra::canonical(&mut p);
        for q in [vec![], vec![0], vec![1, 0], vec![6, 2]] {
            let mut q = q;
            algebra::canonical(&mut q);
            let original = vec![p.clone(), q];
            let possible = (0..8).any(|m| satisfies(&original, m));
            let a = solve_gray(&original, 3, 200000, false);
            let b = solve_gray(&original, 3, 200000, true);
            assert_eq!(a.outcome, b.outcome);
            assert_eq!(a.logical, b.logical);
            assert_eq!(a.trace, b.trace);
            assert_eq!(matches!(a.outcome, Outcome::Sat(_)), possible);
            let embedded: System = original
                .iter()
                .map(|p| p.iter().map(|&m| embed(m)).collect())
                .collect();
            for cut in [8, 12] {
                let x = solve_leaf(&embedded, 16, 200000, cut, false);
                let y = solve_leaf(&embedded, 16, 200000, cut, true);
                assert_eq!(x.outcome, y.outcome);
                assert_eq!(x.logical, y.logical);
                assert_eq!(x.trace, y.trace);
                assert_eq!(matches!(x.outcome, Outcome::Sat(_)), possible);
                if let Outcome::Sat(m) = x.outcome {
                    assert!(satisfies(&embedded, m));
                }
            }
        }
    }
}

#[test]
fn hybrid_disabled_matches_packed_and_active_policies_match_exactly() {
    for n in [12, 16, 20] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let reference = solve_packed(&system, n, 200000);
                let disabled = solve_leaf_backend::<false, 0>(&system, n, 200000);
                assert_eq!(reference.outcome, disabled.outcome);
                assert_eq!(reference.logical, disabled.logical);
                assert_eq!(reference.trace, disabled.trace);
                for cut in [8, 12] {
                    let a = solve_leaf(&system, n, 200000, cut, false);
                    let b = solve_leaf(&system, n, 200000, cut, true);
                    assert_eq!(a.outcome, b.outcome);
                    assert_eq!(a.logical, b.logical);
                    assert_eq!(a.trace, b.trace);
                    assert_eq!(status(&a.outcome), status(&reference.outcome));
                    if let Outcome::Sat(m) = a.outcome {
                        assert!(satisfies(&system, m));
                    }
                }
                let a = solve_gray(&system, n, 200000, false);
                let b = solve_gray(&system, n, 200000, true);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(status(&a.outcome), status(&reference.outcome));
                if let Outcome::Sat(m) = a.outcome {
                    assert!(satisfies(&system, m));
                }
            }
        }
    }
}
