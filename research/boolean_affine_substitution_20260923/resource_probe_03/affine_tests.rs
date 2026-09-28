fn evaluate_map(map: &RecoveryMap, n: u8, point: u64) -> u64 {
    let mut value = 0;
    for i in 0..n as usize {
        let bit = (map[i] & point).count_ones() % 2 ^ ((map[i] >> n) & 1) as u32;
        value |= u64::from(bit) << i;
    }
    value
}
fn small_polynomials() -> System {
    let terms: Vec<_> = (0u64..8).filter(|m| m.count_ones() <= 2).collect();
    (0usize..128)
        .map(|chosen| {
            let mut p: Vec<_> = terms
                .iter()
                .enumerate()
                .filter(|(i, _)| chosen & (1 << i) != 0)
                .map(|(_, &m)| m)
                .collect();
            algebra::canonical(&mut p);
            p
        })
        .collect()
}
fn decode_small_map(code: u64) -> RecoveryMap {
    let mut out = [0; MAX];
    for i in 0..3 {
        out[i] = (code >> (4 * i)) & 15;
    }
    out
}

#[test]
fn all_small_affine_substitutions_match_truth_table_mobius_oracle() {
    let polys = small_polynomials();
    let wide = WideBasis::tail(3);
    let encoded = wide.build(&polys);
    let compact = CompactBackend::new(3);
    let compact_input = compact.build(&polys);
    let truths: Vec<Vec<u8>> = polys
        .iter()
        .map(|p| {
            (0..8)
                .map(|x| p.iter().fold(0u8, |v, &m| v ^ u8::from(x & m == m)))
                .collect()
        })
        .collect();
    let mut cases = 0;
    for code in 0..4096 {
        let map = decode_small_map(code);
        let (got, _) = wide.transform_pair(encoded.clone(), Vec::new(), &map, 3);
        let actual = wide.materialize(&got);
        let (short, _) =
            compact.transform_pair(compact_input.clone(), compact.build(&Vec::new()), &map, 3);
        assert_eq!(compact.materialize(&short), actual);
        let inputs: Vec<_> = (0..8).map(|y| evaluate_map(&map, 3, y) as usize).collect();
        for (i, poly) in polys.iter().enumerate() {
            let mut anf: Vec<_> = inputs.iter().map(|&x| truths[i][x]).collect();
            for bit in [1, 2, 4] {
                for subset in 0..8 {
                    if subset & bit != 0 {
                        anf[subset] ^= anf[subset ^ bit];
                    }
                }
            }
            assert_eq!(anf[7], 0, "Affine substitution must preserve degree <= 2");
            let mut expected: Vec<_> = (0..8).filter(|&m| anf[m] != 0).map(|m| m as u64).collect();
            algebra::canonical(&mut expected);
            assert_eq!(actual[i], expected);
            assert_eq!(substitute_polynomial(poly, &map, 3), expected);
            cases += 1;
        }
    }
    assert_eq!(cases, 524288);
}

#[test]
fn recovery_composition_matches_direct_affine_composition() {
    let mut projectors = 0;
    for eliminated in 0..8 {
        for code in 0..4096 {
            let map = decode_small_map(code);
            if (0..3).any(|j| {
                if eliminated & (1 << j) != 0 {
                    map[j] & eliminated != 0
                } else {
                    map[j] != (1 << j)
                }
            }) {
                continue;
            }
            projectors += 1;
            for old_code in 0..4096 {
                let original = decode_small_map(old_code);
                let mut expected = [0u64; MAX];
                for i in 0..3 {
                    expected[i] = original[i] & 8;
                    for j in 0..3 {
                        if original[i] & (1 << j) != 0 {
                            expected[i] ^= map[j];
                        }
                    }
                }
                let mut actual = original;
                compose_recovery(&mut actual, &map, eliminated, 3);
                assert_eq!(actual, expected);
                assert_eq!(
                    recover_zero_free(&actual, 3),
                    evaluate_map(&original, 3, evaluate_map(&map, 3, 0))
                );
            }
        }
    }
    assert_eq!(projectors, 81);
}

#[test]
fn dense_affine_images_cross_all_coefficient_word_boundaries() {
    let mut rng = 90125;
    for n in [12, 16, 24, 36] {
        let wide = WideBasis::tail(n);
        let all = (1u64 << n) - 1;
        for round in 0..6 {
            let input: System = (0..4)
                .map(|_| {
                    wide.monomials
                        .iter()
                        .copied()
                        .filter(|_| next(&mut rng) % 3 == 0)
                        .collect()
                })
                .collect();
            let mut map = [0u64; MAX];
            for expr in &mut map[..n as usize] {
                *expr = next(&mut rng) & ((1u64 << (n + 1)) - 1);
            }
            if round % 2 == 1 {
                let allowed = (((1u64 << 9) - 1) << (n - 9)) | (1u64 << n);
                for row in &mut map[..n as usize] {
                    *row &= allowed;
                }
            }
            let (got, _) = wide.transform_pair(wide.build(&input), Vec::new(), &map, n);
            let compact = CompactBackend::new(n);
            let (short, _) =
                compact.transform_pair(compact.build(&input), compact.build(&Vec::new()), &map, n);
            assert_eq!(compact.materialize(&short), wide.materialize(&got));
            if round % 2 == 1 {
                assert!(matches!(short.rows, CompactRows::Word(_)));
            }

            let actual = wide.materialize(&got);
            for (original, transformed) in input.iter().zip(&actual) {
                assert_eq!(*transformed, substitute_polynomial(original, &map, n));
                for _ in 0..8 {
                    let y = next(&mut rng) & all;
                    let x = evaluate_map(&map, n, y);
                    let eval = |p: &Vec<u64>, v| p.iter().fold(false, |a, &m| a ^ (v & m == m));
                    assert_eq!(eval(original, x), eval(transformed, y));
                }
            }
        }
    }
}

#[test]
fn complete_affine_solves_match_all_small_solution_sets() {
    let polys = small_polynomials();
    for a in &polys {
        for b in &polys {
            let system = vec![a.clone(), b.clone()];
            let possible = (0..8).any(|m| satisfies(&system, m));
            let x = solve_affine(&system, 3, 100, false);
            let y = solve_affine(&system, 3, 100, true);
            let z = solve_affine_compact(&system, 3, 100);
            assert_eq!(x.outcome, z.outcome);
            assert_eq!(x.logical, z.logical);
            assert_eq!(x.trace, z.trace);
            assert_eq!(x.outcome, y.outcome);
            assert_eq!(x.logical, y.logical);
            assert_eq!(x.trace, y.trace);
            match x.outcome {
                Outcome::Sat(model) => {
                    assert!(possible);
                    assert!(satisfies(&system, model));
                }
                Outcome::Unsat => assert!(!possible),
                Outcome::Unknown(_) => panic!("small exhaustive solve must complete"),
            }
        }
    }
}

#[test]
fn affine_policy_preserves_models_trace_and_censoring_on_larger_systems() {
    for n in [12, 16, 20] {
        for seed in [17, 937, 20261015, 1966087] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let reference = solve_packed(&system, n, 200000);
                let a = solve_affine(&system, n, 200000, false);
                let b = solve_affine(&system, n, 200000, true);
                let c = solve_affine_compact(&system, n, 200000);
                assert_eq!(a.outcome, c.outcome);
                assert_eq!(a.logical, c.logical);
                assert_eq!(a.trace, c.trace);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(status(&a.outcome), status(&reference.outcome));
                if let Outcome::Sat(model) = a.outcome {
                    assert!(satisfies(&system, model));
                }
            }
        }
    }
    let (system, _) = fixture(16, 17, "unplanted");
    for cap in [0, 1, 2, 7, 31] {
        let a = solve_affine(&system, 16, cap, false);
        let b = solve_affine(&system, 16, cap, true);
        let c = solve_affine_compact(&system, 16, cap);
        assert_eq!(a.outcome, c.outcome);
        assert_eq!(a.logical, c.logical);
        assert_eq!(a.trace, c.trace);
        assert_eq!(a.outcome, b.outcome);
        assert_eq!(a.logical, b.logical);
        assert_eq!(a.trace, b.trace);
        assert!(a.logical.nodes <= cap);
        if cap == 0 {
            assert_eq!(a.outcome, Outcome::Unknown("NODE_CAP"));
        }
    }
}

#[test]
fn compact_point_schedules_are_exact_across_label_embeddings_and_word_transition() {
    let compact = CompactBackend::new(16);
    for labels in [[0u8, 1, 2], [1, 3, 7], [4, 8, 11]] {
        let active = labels.iter().fold(0u64, |v, &j| v | (1u64 << j));
        let embed = |m: u64| {
            (0..3).fold(0, |v, j| {
                v | if m & (1 << j) != 0 {
                    1u64 << labels[j]
                } else {
                    0
                }
            })
        };
        let polys: System = small_polynomials()
            .iter()
            .map(|p| p.iter().map(|&m| embed(m)).collect())
            .collect();
        let state = compact.build_at(&polys, active);
        assert!(state
            .layout
            .monomials
            .windows(2)
            .all(|w| algebra::mono_order(w[0], w[1]) == std::cmp::Ordering::Less));
        for mask in 0..8 {
            for values in (0..8).filter(|v| v & !mask == 0) {
                let a = compact.specialize(state.clone(), embed(mask), embed(values), active);
                let b = TailList.specialize(polys.clone(), embed(mask), embed(values), active);
                assert_eq!(compact.materialize(&a), b);
            }
        }
    }
    for active in [(1u64 << 11) - 1, ((1u64 << 10) - 1) | (1u64 << 15)] {
        let layout = compact.layout(active);
        let input: System = layout.monomials.iter().map(|&m| vec![m]).collect();
        let initial = compact.build_at(&input, active);
        assert!(matches!(initial.rows, CompactRows::Wide(_)));
        for &j in &layout.labels {
            for one in [false, true] {
                let bit = 1u64 << j;
                let value = if one { bit } else { 0 };
                let got = compact.specialize(initial.clone(), bit, value, active);
                assert!(matches!(got.rows, CompactRows::Word(_)));
                assert_eq!(
                    compact.materialize(&got),
                    TailList.specialize(input.clone(), bit, value, active)
                );
            }
        }
    }
    assert!(compact.points.borrow().len() <= 256);
    assert!(compact.layouts.borrow().len() <= 128);
}
