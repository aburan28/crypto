#[test]
fn every_small_row_span_commutes_with_every_assignment() {
    let monomials: Vec<_> = (0u64..8).filter(|m| m.count_ones() <= 2).collect();
    let polys: Vec<Vec<u64>> = (0..128)
        .map(|chosen| {
            let mut p: Vec<_> = monomials
                .iter()
                .enumerate()
                .filter(|(i, _)| chosen & (1 << i) != 0)
                .map(|(_, &m)| m)
                .collect();
            algebra::canonical(&mut p);
            p
        })
        .collect();
    let wide = WideBasis::new(3);
    for a in &polys {
        for b in &polys {
            let original = vec![a.clone(), b.clone()];
            let list_basis = ListBasis.reduce(original.clone());
            let word_basis = wide.reduce(wide.build(&original));
            assert_eq!(wide.materialize(&word_basis), list_basis);
            for model in 0..8 {
                assert_eq!(satisfies(&original, model), satisfies(&list_basis, model));
            }
            for mask in 0..8 {
                for values in (0..8).filter(|v| v & !mask == 0) {
                    let rebuilt =
                        ListBasis.reduce(ListBasis.specialize(original.clone(), mask, values, 7));
                    let transported =
                        ListBasis.reduce(ListBasis.specialize(list_basis.clone(), mask, values, 7));
                    let actual = wide.reduce(wide.specialize(word_basis.clone(), mask, values, 7));
                    assert_eq!(rebuilt, transported);
                    assert_eq!(wide.materialize(&actual), rebuilt);
                }
            }
        }
    }
}

#[test]
fn dense_transport_handles_word_boundaries_and_repeated_degree_drops() {
    for n in [12, 16, 24, 36] {
        let wide = WideBasis::new(n);
        assert!(wide
            .monomials
            .windows(2)
            .all(|w| algebra::mono_order(w[0], w[1]) == std::cmp::Ordering::Less));
        let mut rng = 31337;
        let original: System = (0..12)
            .map(|_| {
                wide.monomials
                    .iter()
                    .copied()
                    .filter(|_| next(&mut rng) % 4 == 0)
                    .collect()
            })
            .collect();
        let mut carried = wide.reduce(wide.build(&original));
        let (mut known, mut values) = (0, 0);
        let all = (1u64 << n) - 1;
        for step in 0..n {
            let j = if step % 2 == 0 {
                step / 2
            } else {
                n - 1 - step / 2
            };
            let bit = 1u64 << j;
            let one = if next(&mut rng) & 1 != 0 { bit } else { 0 };
            carried = wide.reduce(wide.specialize(carried, bit, one, all & !known));
            known |= bit;
            values |= one;
            let rebuilt =
                ListBasis.reduce(ListBasis.specialize(original.clone(), known, values, all));
            assert_eq!(wide.materialize(&carried), rebuilt);
            assert!(carried.len() <= original.len());
        }
    }
}

#[test]
fn basis_backends_match_complete_models_counters_and_trace() {
    for n in [12, 16, 20] {
        for seed in [17, 937, 20261014, 1835009] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let reference = solve_packed(&system, n, 200000);
                let a = solve_basis(&system, n, 200000, false);
                let b = solve_basis(&system, n, 200000, true);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(
                    a.profile.kernel_calls_by_active,
                    b.profile.kernel_calls_by_active
                );
                assert_eq!(status(&a.outcome), status(&reference.outcome));
                if let Outcome::Sat(model) = a.outcome {
                    assert!(satisfies(&system, model));
                }
            }
        }
    }
}

#[test]
fn basis_node_caps_are_censored_and_backend_independent() {
    let (system, _) = fixture(16, 17, "unplanted");
    for limit in [0, 1, 2, 3, 7, 29] {
        let a = solve_basis(&system, 16, limit, false);
        let b = solve_basis(&system, 16, limit, true);
        assert_eq!(a.outcome, b.outcome);
        assert_eq!(a.logical, b.logical);
        assert_eq!(a.trace, b.trace);
        assert!(a.logical.nodes <= limit);
        if limit == 0 {
            assert_eq!(a.outcome, Outcome::Unknown("NODE_CAP"));
        }
    }
}
