#[test]
fn quadratic_product_admission_matches_direct_boolean_multiplication() {
    let backend = CompactBackend::new(3);
    let layout = backend.layout(7);
    for chosen in 0u64..128 {
        let poly: Vec<_> = layout
            .monomials
            .iter()
            .enumerate()
            .filter(|(i, _)| chosen & (1 << i) != 0)
            .map(|(_, &m)| m)
            .collect();
        let common = quadratic_common_variables(chosen, &layout);
        for j in 0..3 {
            let mut product: Vec<_> = poly.iter().map(|&m| m | (1u64 << j)).collect();
            algebra::canonical(&mut product);
            let admitted = product.iter().all(|m| m.count_ones() <= 2);
            assert_eq!(common & (1 << j) != 0, admitted);
            if admitted {
                let actual = multiply_word_by_variable(chosen, &layout, j);
                let state = CompactState {
                    layout: layout.clone(),
                    rows: CompactRows::Word(vec![actual]),
                };
                assert_eq!(backend.materialize(&state), vec![product]);
            }
        }
    }
}

#[test]
fn degree_two_products_preserve_every_small_solution_set_and_match_references() {
    let layout = CompactLayout::new(7);
    let polys: System = (0u64..128)
        .map(|chosen| {
            layout
                .monomials
                .iter()
                .enumerate()
                .filter(|(i, _)| chosen & (1 << i) != 0)
                .map(|(_, &m)| m)
                .collect()
        })
        .collect();
    let reference = LinearizedList::new(3);
    let compact = LinearizedCompact::new(3);
    for a in &polys {
        for b in &polys {
            let original = vec![a.clone(), b.clone()];
            let expected = reference.reduce(original.clone());
            let ref_work = reference.take_extra_reduction();
            let got = compact.reduce(compact.build(&original));
            let work = compact.take_extra_reduction();
            assert_eq!(compact.inner.materialize(&got), expected);
            assert_eq!(work, ref_work);
            for model in 0..8 {
                assert_eq!(satisfies(&original, model), satisfies(&expected, model));
            }
            let possible = (0..8).any(|m| satisfies(&original, m));
            let x = solve_affine_linearized(&original, 3, 100, false);
            let y = solve_affine_linearized(&original, 3, 100, true);
            assert_eq!(x.outcome, y.outcome);
            assert_eq!(x.logical, y.logical);
            assert_eq!(x.trace, y.trace);
            match x.outcome {
                Outcome::Sat(model) => {
                    assert!(possible);
                    assert!(satisfies(&original, model));
                }
                Outcome::Unsat => assert!(!possible),
                Outcome::Unknown(_) => panic!("small system must complete"),
            }
        }
    }
}

#[test]
fn products_derive_a_forced_boolean_factor_and_preserve_the_width_gate() {
    let system = vec![vec![3, 1, 0]];
    let got = solve_affine_linearized(&system, 3, 100, true);
    assert_eq!(got.outcome, Outcome::Sat(1));
    assert_eq!(got.logical.decisions, 0);
    assert!(got.logical.derived_rows > 0);
    let wide = solve_affine_linearized(&system, 12, 100, true);
    let prior = solve_affine_compact(&system, 12, 100);
    assert_eq!(wide.logical.derived_rows, 0);
    assert_eq!(wide.outcome, prior.outcome);
    assert_eq!(wide.trace, prior.trace);
}

#[test]
fn linearized_affine_solvers_match_complete_controls_and_node_caps() {
    for n in [12, 16, 20] {
        for seed in [17, 937, 20261015, 1966087] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let oracle = solve_packed(&system, n, 200000);
                let a = solve_affine_linearized(&system, n, 200000, false);
                let b = solve_affine_linearized(&system, n, 200000, true);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(status(&a.outcome), status(&oracle.outcome));
                if let Outcome::Sat(model) = a.outcome {
                    assert!(satisfies(&system, model));
                }
            }
        }
    }
    let (system, _) = fixture(16, 17, "unplanted");
    for limit in [0, 1, 2, 7, 31] {
        let a = solve_affine_linearized(&system, 16, limit, false);
        let b = solve_affine_linearized(&system, 16, limit, true);
        assert_eq!(a.outcome, b.outcome);
        assert_eq!(a.logical, b.logical);
        assert_eq!(a.trace, b.trace);
        if limit == 0 {
            assert_eq!(a.outcome, Outcome::Unknown("NODE_CAP"));
        }
    }
}
