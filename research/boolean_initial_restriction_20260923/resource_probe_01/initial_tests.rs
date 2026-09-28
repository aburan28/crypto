fn enumeration_signature(e: &Enumeration) -> (Option<u64>, bool, u64, u64, u64) {
    (e.model, e.complete, e.points, e.batches, e.checksum)
}
#[test]
fn transported_blocks_match_every_direct_point_through_twelve_variables() {
    let mut rng = 937193;
    for n in 4..=12 {
        let mut form = SyndromeForm::new(n);
        form.constant = next(&mut rng) as u32;
        for i in 0..n {
            form.linear[i] = next(&mut rng) as u32;
            for j in i + 1..n {
                form.quadratic[i][j] = next(&mut rng) as u32;
                form.quadratic[j][i] = form.quadratic[i][j];
            }
        }
        let mut a = DeltaCursor::<ScalarDeltaBlock>::new(&form);
        let mut b = DeltaCursor::<NativeDeltaBlock>::new(&form);
        for step in 0..1 << n - 4 {
            a.advance(&form, step);
            b.advance(&form, step);
            let mut expected = [[0; 4]; 4];
            for y in 0..16 {
                expected[y / 4][y % 4] = form.value(((step ^ (step >> 1)) << 4) | y as u64);
            }
            assert_eq!(a.block.values(), expected);
            assert_eq!(b.block.values(), expected);
            assert_eq!(a.block.inspect(), inspect_delta_values(&expected));
            assert_eq!(b.block.inspect(), inspect_delta_values(&expected));
        }
        assert_eq!(
            enumeration_signature(&enumerate_syndromes::<true>(&form, 1 << 24)),
            enumeration_signature(&enumerate_delta::<NativeDeltaBlock>(&form, 1 << 24))
        );
    }
}
#[test]
fn transported_blocks_preserve_every_hit_lane_and_wrapping_sum() {
    let mut rng = 179;
    for _ in 0..128 {
        let original = std::array::from_fn(|_| std::array::from_fn(|_| next(&mut rng) as u32));
        for hit in 0..16 {
            let mut values = original;
            values[hit / 4][hit % 4] = 0;
            assert_eq!(
                NativeDeltaBlock::from_values(&values).inspect(),
                inspect_delta_values(&values)
            );
        }
    }
    assert_eq!(
        NativeDeltaBlock::from_values(&[[u32::MAX; 4]; 4]).inspect(),
        (u32::MAX - 15, None)
    );
}
#[test]
fn transported_enumeration_matches_censored_and_tiny_controls() {
    for n in 0..=12 {
        for system in [vec![], vec![vec![0]]] {
            let form = SyndromeForm::from_system(&system, n).unwrap();
            for cap in [0, 1, 15, 16, 17, 31, 32, (1 << n) - 1, 1 << n, (1 << n) + 1] {
                let old = enumerate_syndromes::<true>(&form, cap);
                for got in [
                    enumerate_delta::<ScalarDeltaBlock>(&form, cap),
                    enumerate_delta::<NativeDeltaBlock>(&form, cap),
                ] {
                    assert_eq!(enumeration_signature(&old), enumeration_signature(&got));
                    assert!(got.points <= cap);
                }
            }
        }
    }
}
#[test]
fn every_small_row_space_has_exact_initial_restriction_and_recovery() {
    let monomials = [0, 1, 2, 4, 3, 5, 6];
    let polynomials: Vec<Vec<u64>> = (0..128)
        .map(|mask| {
            let mut p: Vec<_> = monomials
                .iter()
                .enumerate()
                .filter_map(|(j, &m)| if mask & (1 << j) != 0 { Some(m) } else { None })
                .collect();
            algebra::canonical(&mut p);
            p
        })
        .collect();
    for first in &polynomials {
        for second in &polynomials {
            let system = vec![first.clone(), second.clone()];
            let list = initial_affine_rows(TailList, &system, 3);
            let wide = initial_affine_rows(WideBasis::tail(3), &system, 3);
            assert_eq!(list, wide);
            let old: Vec<_> = (0..8).filter(|&p| satisfies(&system, p)).collect();
            match InitialMap::from_rows(&wide, 3) {
                None => assert!(old.is_empty()),
                Some(map) => {
                    let original = SyndromeForm::from_system(&system, 3).unwrap();
                    let transformed = transform_syndrome(&original, &map);
                    let mut recovered = Vec::new();
                    let mut image = std::collections::BTreeSet::new();
                    for y in 0..1 << map.free {
                        let x = map.recover(y, 3);
                        assert!(image.insert(x));
                        assert_eq!(transformed.value(y), original.value(x));
                        if transformed.value(y) == 0 {
                            recovered.push(x);
                        }
                    }
                    recovered.sort_unstable();
                    assert_eq!(old, recovered);
                }
            }
        }
    }
    // The affine row requires a combination of THREE nonlinear input rows.
    let mut system = vec![vec![3, 1], vec![6, 2], vec![3, 6, 4]];
    for p in &mut system {
        algebra::canonical(p);
    }
    assert_eq!(initial_affine_rows(WideBasis::tail(3), &system, 3), vec![7]);
}
#[test]
fn packed_affine_images_match_all_small_maps_and_all_quadratics() {
    let monomials = [0, 1, 2, 4, 3, 5, 6];
    for batch in 0..4 {
        let system: System = (batch * 32..(batch + 1) * 32)
            .map(|mask| {
                monomials
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &m)| if mask & (1 << j) != 0 { Some(m) } else { None })
                    .collect()
            })
            .collect();
        let original = SyndromeForm::from_system(&system, 3).unwrap();
        for code in 0..4096 {
            let mut images = [0; MAX];
            for j in 0..3 {
                images[j] = (code >> (4 * j)) & 15;
            }
            let map = InitialMap {
                images,
                free: 3,
                rank: 0,
            };
            let actual = transform_syndrome(&original, &map);
            for y in 0..8 {
                assert_eq!(actual.value(y), original.value(map.recover(y, 3)));
            }
        }
    }
}
#[test]
fn dense_maps_handle_dimension_changes_boolean_diagonals_and_full_equation_word() {
    let mut rng = 9019283;
    for n in [4, 8, 16, 24] {
        for free in [0, 1, 3, 8, 16, 24] {
            let mut original = SyndromeForm::new(n);
            original.constant = next(&mut rng) as u32;
            let mut images = [0; MAX];
            for i in 0..n {
                images[i] = next(&mut rng) & ((1u64 << (free + 1)) - 1);
                original.linear[i] = next(&mut rng) as u32;
                for j in i + 1..n {
                    original.quadratic[i][j] = next(&mut rng) as u32;
                    original.quadratic[j][i] = original.quadratic[i][j];
                }
            }
            let map = InitialMap {
                images,
                free,
                rank: 0,
            };
            let transformed = transform_syndrome(&original, &map);
            for _ in 0..64 {
                let y = next(&mut rng) & ((1 << free) - 1);
                assert_eq!(
                    transformed.value(y),
                    original.value(map.recover(y, n as u8))
                );
            }
        }
    }
}
#[test]
fn initial_and_delta_complete_policies_match_models_work_and_trace() {
    for n in [12, 16] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let old = solve_gray(&system, n, 200000, true);
                for got in [
                    solve_gray_delta(&system, n, 200000, false),
                    solve_gray_delta(&system, n, 200000, true),
                ] {
                    assert_eq!(old.outcome, got.outcome);
                    assert_eq!(old.logical, got.logical);
                    assert_eq!(old.trace, got.trace);
                }
                let a = solve_initial(&system, n, 200000, false);
                let b = solve_initial(&system, n, 200000, true);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(status(&a.outcome), status(&old.outcome));
                if let Outcome::Sat(x) = a.outcome {
                    assert!(satisfies(&system, x));
                }
                let old = solve_leaf_backend::<true, 16>(&system, n, 200000);
                for got in [
                    solve_leaf_delta::<false>(&system, n, 200000),
                    solve_leaf_delta::<true>(&system, n, 200000),
                ] {
                    assert_eq!(old.outcome, got.outcome);
                    assert_eq!(old.logical, got.logical);
                    assert_eq!(old.trace, got.trace);
                }
            }
        }
    }
    for system in [
        vec![],
        vec![vec![0]],
        vec![vec![1, 1, 0, 0]],
        vec![vec![3, 1], vec![3, 2, 0]],
        vec![vec![7]],
    ] {
        for n in [3, 25] {
            let expected = solve_packed(&system, n, 200000);
            for optimized in [false, true] {
                let got = solve_initial(&system, n, 200000, optimized);
                assert_eq!(status(&got.outcome), status(&expected.outcome));
                if let Outcome::Sat(x) = got.outcome {
                    assert!(satisfies(&system, x));
                }
                assert!(matches!(
                    solve_initial(&system, n, 0, optimized).outcome,
                    Outcome::Unknown(_)
                ));
            }
        }
    }
}
