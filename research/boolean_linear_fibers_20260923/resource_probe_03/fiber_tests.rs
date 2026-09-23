#[test]
fn maximum_sets_match_every_simple_graph_through_five_vertices() {
    for n in 0..=5 {
        let mut edges = Vec::new();
        for j in 0..n {
            for i in 0..j {
                edges.push((i, j));
            }
        }
        for code in 0..1u64 << edges.len() {
            let mut graph = [0; MAX];
            for (p, &(i, j)) in edges.iter().enumerate() {
                if code & (1 << p) != 0 {
                    graph[i] |= 1 << j;
                    graph[j] |= 1 << i;
                }
            }
            let mut expected = 0;
            for mask in 0u32..1 << n {
                if (0..n).all(|i| mask & (1 << i) == 0 || mask & graph[i] == 0) {
                    expected = better_independent(expected, mask);
                }
            }
            assert_eq!(maximum_independent(&graph, n).0, expected);
        }
    }
}
#[test]
fn every_small_linear_system_matches_its_complete_solution_set() {
    let mut checked = 0u64;
    for equations in 0..=4 {
        for k in 0..=4 {
            for code in 0u64..1u64 << (equations * k) {
                let columns: Vec<u32> = (0..k)
                    .map(|j| ((code >> (equations * j)) & ((1 << equations) - 1)) as u32)
                    .collect();
                let mut first = [None; 16];
                for point in 0u32..1 << k {
                    let value = columns.iter().enumerate().fold(0, |v, (j, &a)| {
                        v ^ if point & (1 << j) != 0 { a } else { 0 }
                    });
                    if first[value as usize].is_none() {
                        first[value as usize] = Some(point);
                    }
                }
                let size = first.iter().filter(|v| v.is_some()).count();
                assert!(size.is_power_of_two());
                let rank = size.trailing_zeros();
                for rhs in 0..1 << equations {
                    let expected = LinearAnswer {
                        model: first[rhs as usize],
                        rank,
                    };
                    assert_eq!(fiber_rows(&columns, rhs, equations), expected);
                    assert_eq!(fiber_columns(&columns, rhs), expected);
                    checked += 1;
                }
            }
        }
    }
    assert_eq!(
        checked,
        (0..=4)
            .flat_map(|m| (0..=4).map(move |k| 1u64 << (m * k + m)))
            .sum::<u64>()
    );
}
#[test]
fn linear_membership_handles_full_words_dependent_columns_and_rank_changes() {
    let mut rng = 71583;
    for k in 0..=24 {
        for _ in 0..64 {
            let columns: Vec<_> = (0..k).map(|_| next(&mut rng) as u32).collect();
            for rhs in [
                0,
                next(&mut rng) as u32,
                columns.iter().fold(0, |v, a| v ^ a),
            ] {
                let a = fiber_rows(&columns, rhs, 32);
                let b = fiber_columns(&columns, rhs);
                assert_eq!(a, b);
                if let Some(y) = a.model {
                    assert!(y < 1 << k);
                    assert_eq!(
                        columns
                            .iter()
                            .enumerate()
                            .fold(0, |v, (j, &a)| v ^ if y & (1 << j) != 0 { a } else { 0 }),
                        rhs
                    );
                }
            }
        }
    }
    for columns in [
        vec![0, 0, 0],
        vec![1 << 31, 1 << 31, 0],
        vec![1, 2, 3],
        vec![3, 5, 6],
    ] {
        for rhs in [0, 1, 2, 3, 1 << 31, u32::MAX] {
            assert_eq!(fiber_rows(&columns, rhs, 32), fiber_columns(&columns, rhs));
        }
    }
}
#[test]
fn every_fiber_coefficient_and_screen_matches_direct_original_evaluation() {
    let mut rng = 762037;
    for n in 4..=12 {
        let mut system = Vec::new();
        for _ in 0..32 {
            let mut p = Vec::new();
            for i in 0..n {
                if next(&mut rng) & 1 != 0 {
                    p.push(1u64 << i);
                }
                for j in i + 1..n {
                    if next(&mut rng) % 7 == 0 {
                        p.push((1u64 << i) | (1u64 << j));
                    }
                }
            }
            if next(&mut rng) & 1 != 0 {
                p.push(0);
            }
            algebra::canonical(&mut p);
            system.push(p);
        }
        let original = SyndromeForm::from_system(&system, n).unwrap();
        let selected = maximum_independent(&fiber_graph(&original), n as usize).0;
        let plan = FiberPlan::new(&original, 32, selected);
        let h = plan.outside.len();
        let size = 1usize << h.min(4);
        let mut cursor = GrayCursor::new(&plan.form);
        let offsets = plan.form.low_offsets();
        let mut columns = plan.columns.clone();
        for step in 0..1u64 << h.saturating_sub(4) {
            cursor.advance(&plan.form, step);
            if step != 0 {
                let j = step.trailing_zeros() as usize;
                for i in 0..columns.len() {
                    columns[i] ^= plan.cross_high[j][i];
                }
            }
            let screen =
                screen_fibers_scalar(&cursor, &offsets, &columns, &plan.column_offsets, size);
            if size == 16 {
                let native =
                    screen_fibers_native(&cursor, &offsets, &columns, &plan.column_offsets);
                assert_eq!(screen.survivors, native.survivors);
                assert_eq!(screen.values, native.values);
                assert_eq!(screen.rounds, native.rounds);
            }
            for y in 0..size {
                let outside = ((step ^ (step >> 1)) << 4) | y as u64;
                let b = original.value(plan.recover(outside, 0));
                assert_eq!(screen.values[y / 4][y % 4], b);
                let actual: Vec<_> = (0..columns.len())
                    .map(|i| columns[i] ^ plan.column_offsets[i][y / 4][y % 4])
                    .collect();
                for (i, &a) in actual.iter().enumerate() {
                    assert_eq!(a, original.value(plan.recover(outside, 1 << i)) ^ b);
                }
                let union = actual.iter().fold(0, |a, b| a | b);
                let mut contradiction = b & !union;
                for shift in FIBER_SHIFTS {
                    let transformed = actual.iter().fold(0, |v, &a| v | (a ^ (a >> shift)));
                    contradiction |= (b ^ (b >> shift)) & !transformed;
                }
                assert_eq!(screen.survivors & (1 << y) != 0, contradiction == 0);
                let mut first = None;
                for inside in 0..1u32 << columns.len() {
                    let v = actual.iter().enumerate().fold(b, |v, (j, &a)| {
                        v ^ if inside & (1 << j) != 0 { a } else { 0 }
                    });
                    assert_eq!(v, original.value(plan.recover(outside, inside)));
                    if v == 0 && first.is_none() {
                        first = Some(inside);
                    }
                }
                assert_eq!(fiber_columns(&actual, b).model, first);
            }
        }
    }
}
#[test]
fn simd_screen_preserves_all_hit_lanes_and_high_equation_bits() {
    let mut rng = 819371;
    for k in 0..=6 {
        for _ in 0..64 {
            let cursor = GrayCursor {
                constant: next(&mut rng) as u32,
                linear: std::array::from_fn(|_| next(&mut rng) as u32),
                differences: [0; MAX],
            };
            let columns: Vec<_> = (0..k).map(|_| next(&mut rng) as u32).collect();
            let offsets = std::array::from_fn(|_| std::array::from_fn(|_| next(&mut rng) as u32));
            let images: Vec<[[u32; 4]; 4]> = (0..k)
                .map(|_| std::array::from_fn(|_| std::array::from_fn(|_| next(&mut rng) as u32)))
                .collect();
            for target in 0..16 {
                let mut forced = offsets;
                let mut v = cursor.constant;
                for j in 0..4 {
                    if target & (1 << j) != 0 {
                        v ^= cursor.linear[j];
                    }
                }
                forced[target / 4][target % 4] = v;
                let a = screen_fibers_scalar(&cursor, &forced, &columns, &images, 16);
                let b = screen_fibers_native(&cursor, &forced, &columns, &images);
                assert_ne!(a.survivors & (1 << target), 0);
                assert_eq!(a.survivors, b.survivors);
                assert_eq!(a.values, b.values);
                assert_eq!(a.rounds, b.rounds);
            }
        }
    }
}
#[test]
fn complete_fiber_policies_match_small_truth_tables_and_each_other() {
    let monomials = [0, 1, 2, 4, 3, 5, 6];
    for code in 0..128 {
        let p: Vec<_> = monomials
            .iter()
            .enumerate()
            .filter_map(|(j, &m)| if code & (1 << j) != 0 { Some(m) } else { None })
            .collect();
        for other in [
            vec![],
            vec![0],
            vec![1],
            vec![1, 0],
            vec![2, 4],
            vec![3, 1, 0],
        ] {
            let mut system = vec![p.clone(), other];
            for p in &mut system {
                algebra::canonical(p);
            }
            let expected = brute(&system, 3);
            let mut reference = None;
            for mode in 0..=2 {
                let got = solve_fiber(&system, 3, 200000, mode);
                assert_eq!(status(&expected), status(&got.outcome));
                if let Outcome::Sat(point) = got.outcome {
                    assert!(satisfies(&system, point));
                }
                let signature = (got.outcome, got.logical, got.trace);
                if let Some(ref r) = reference {
                    assert_eq!(r, &signature);
                } else {
                    reference = Some(signature);
                }
            }
        }
    }
    for n in [12, 16] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let expected = solve_gray(&system, n, 200000, true);
                let mut reference = None;
                for mode in 0..=2 {
                    let got = solve_fiber(&system, n, 200000, mode);
                    assert_eq!(status(&expected.outcome), status(&got.outcome));
                    if let Outcome::Sat(point) = got.outcome {
                        assert!(satisfies(&system, point));
                    }
                    let signature = (got.outcome, got.logical, got.trace);
                    if let Some(ref r) = reference {
                        assert_eq!(r, &signature);
                    } else {
                        reference = Some(signature);
                    }
                }
            }
        }
    }
}
#[test]
fn fiber_caps_and_rejected_extensions_have_exact_accounting() {
    for n in 0..=8 {
        let system = vec![vec![0]];
        let original = SyndromeForm::from_system(&system, n).unwrap();
        for selected in [0, (1u32 << n) - 1] {
            let plan = FiberPlan::new(&original, 1, selected);
            let h = plan.outside.len();
            let size = 1u64 << h.min(4);
            for cap in [0, 1, size - 1, size, size + 1, (1 << h) - 1, 1 << h] {
                let mut reference = None;
                for mode in 0..=2 {
                    let got = enumerate_fibers(&plan, cap, mode);
                    assert!(got.model.is_none());
                    assert!(got.logical.fiber_prefixes <= cap);
                    assert_eq!(got.complete, got.logical.fiber_prefixes == 1 << h);
                    assert_eq!(got.logical.fiber_prefixes, got.logical.fiber_zero_rejected);
                    assert_eq!(got.logical.fiber_queries, 0);
                    assert_eq!(
                        got.logical.fiber_extensions_rejected,
                        got.logical.fiber_prefixes << plan.low.len()
                    );
                    let signature = (got.model, got.complete, got.logical, got.trace);
                    if let Some(ref r) = reference {
                        assert_eq!(r, &signature);
                    } else {
                        reference = Some(signature);
                    }
                }
            }
        }
    }
}
