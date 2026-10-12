#[test]
fn byte_affine_images_and_every_zero_lane_match_scalar() {
    let mut rng = 712983;
    for constant in 0..=255u8 {
        for _ in 0..32 {
            let linear = std::array::from_fn(|_| next(&mut rng) as u8);
            let a = ScalarByte64::affine(constant, &linear);
            let b = NativeByte64::affine(constant, &linear);
            let planes = NativePlane64::affine(constant, &linear);
            assert_eq!(a.values(), planes.values());
            assert_eq!(a.zero_masks(), planes.zero_masks());
            assert_eq!(a.values(), b.values());
            assert_eq!(a.zero_masks(), b.zero_masks());
        }
    }
    for value in [1, 2, 127, 128, 255] {
        for lane in 0..64 {
            let mut values = [value; 64];
            values[lane] = 0;
            let mut expected = [0; 4];
            expected[lane / 16] = 1 << (lane % 16);
            assert_eq!(NativeByte64::from_values(&values).zero_masks(), expected);
            assert_eq!(NativePlane64::from_values(&values).zero_masks(), expected);
        }
    }
    assert_eq!(
        NativeByte64::from_values(&[0; 64]).zero_masks(),
        [u16::MAX; 4]
    );
    assert_eq!(NativeByte64::from_values(&[255; 64]).zero_masks(), [0; 4]);
}
#[test]
fn sixty_four_point_order_is_the_exact_legacy_order() {
    for step in 0..1u64 << 18 {
        for r in 0..4 {
            let old = step * 4 + r as u64;
            assert_eq!(
                ((step ^ (step >> 1)) << 2) | legacy_group64(step, r) as u64,
                old ^ (old >> 1)
            );
        }
    }
}
#[test]
fn full_and_projected_transport_match_every_point_through_twelve_variables() {
    let mut rng = 761923;
    for n in 6..=12 {
        let mut form = SyndromeForm::new(n);
        form.constant = next(&mut rng) as u32;
        for i in 0..n {
            form.linear[i] = next(&mut rng) as u32;
            for j in i + 1..n {
                form.quadratic[i][j] = next(&mut rng) as u32;
                form.quadratic[j][i] = form.quadratic[i][j];
            }
        }
        let offsets = low_quadratic_offsets64(&form);
        let mut cursor = SixCursor::new(&form);
        let mut scalar = ScalarByte64::affine(form.constant as u8, &cursor.linear.map(|v| v as u8));
        let mut native = NativeByte64::affine(form.constant as u8, &cursor.linear.map(|v| v as u8));
        scalar.xor_block(&ScalarByte64::from_values(&offsets.map(|v| v as u8)));
        native.xor_block(&NativeByte64::from_values(&offsets.map(|v| v as u8)));
        let mut word = Word64::affine(form.constant, &cursor.linear);
        word.xor_block(&Word64::from_values(&offsets));
        let crosses: Vec<_> = cursor
            .cross
            .iter()
            .map(|a| {
                (
                    ScalarByte64::affine(0, &a.map(|v| v as u8)),
                    NativeByte64::affine(0, &a.map(|v| v as u8)),
                    Word64::affine(0, a),
                )
            })
            .collect();
        for step in 0..1u64 << (n - 6) {
            if let Some((j, d)) = cursor.advance::<true>(&form, step) {
                scalar.xor_uniform(d as u8);
                native.xor_uniform(d as u8);
                word.xor_uniform(d);
                scalar.xor_block(&crosses[j].0);
                native.xor_block(&crosses[j].1);
                word.xor_block(&crosses[j].2);
            }
            let a = scalar.values();
            let b = native.values();
            assert_eq!(a, b);
            assert_eq!(scalar.zero_masks(), native.zero_masks());
            for low in 0..64 {
                let point = ((step ^ (step >> 1)) << 6) | low as u64;
                let expected = form.value(point);
                assert_eq!(cursor.value(&offsets, low), expected);
                assert_eq!(u32::from(a[low]), expected & 255);
                assert_eq!(word.0[low / 16].values()[(low % 16) / 4][low % 4], expected);
            }
        }
    }
}
#[test]
fn projected_false_positives_require_full_equation_checks() {
    for n in 6..=10 {
        let mut form = SyndromeForm::new(n);
        form.constant = 1 << 31;
        for got in [
            enumerate_bytes::<ScalarByte64>(&form, 1 << n),
            enumerate_bytes::<NativeByte64>(&form, 1 << n),
        ] {
            assert!(got.complete && got.model.is_none());
            assert_eq!(got.logical.screen_points, 1 << n);
            assert_eq!(got.logical.screen_passes, 1 << n);
            assert_eq!(got.logical.screen_full_checks, 1 << n);
            assert_eq!(got.logical.screen_full_rejected, 1 << n);
        }
        form.constant = 0;
        let got = enumerate_bytes::<NativeByte64>(&form, 1 << n);
        assert_eq!(got.model, Some(0));
        assert_eq!(got.logical.screen_points, 64);
        assert_eq!(got.logical.screen_passes, 64);
        assert_eq!(got.logical.screen_full_checks, 1);
        form.constant = 255;
        let got = enumerate_bytes::<NativeByte64>(&form, 1 << n);
        assert!(got.complete && got.model.is_none());
        assert_eq!(got.logical.screen_full_checks, 0);
        form.constant = 1 << 8;
        let got = enumerate_tiered::<NativeByte64, false>(&form, 1 << n);
        assert!(got.complete && got.model.is_none());
        assert_eq!(got.logical.screen_second_checks, 1 << n);
        assert_eq!(got.logical.screen_second_rejected, 1 << n);
        assert_eq!(got.logical.screen_full_checks, 0);
    }
}

#[test]
fn packed_secondary_coefficients_and_quiet_modes_preserve_solutions() {
    let mut rng = 571938;
    for _ in 0..4096 {
        let linear = std::array::from_fn(|_| next(&mut rng) as u32);
        let constant = next(&mut rng) as u32;
        let offset = next(&mut rng) as u32;
        let packed = packed_secondary(&linear);
        for low in 0..64 {
            let expected = linear
                .iter()
                .enumerate()
                .fold(constant ^ offset, |v, (i, &c)| {
                    v ^ if low & (1 << i) != 0 { c } else { 0 }
                });
            assert_eq!(
                secondary_value(constant, offset, packed, low),
                (expected >> 8) as u8
            );
        }
    }
    for n in [8, 12, 16] {
        for family in ["planted", "cross_planted", "unplanted"] {
            let (system, _) = fixture(n, 17, family);
            let form = SyndromeForm::from_system(&system, n).unwrap();
            let reference = enumerate_quiet(&form, 1 << n);
            let traced = enumerate_tiered::<ScalarByte64, true>(&form, 1 << n);
            let a = enumerate_tiered::<ScalarByte64, false>(&form, 1 << n);
            let b = enumerate_tiered::<NativeByte64, false>(&form, 1 << n);
            let planes = enumerate_tiered::<NativePlane64, false>(&form, 1 << n);
            assert_eq!(
                (a.model, a.complete, a.logical.clone()),
                (planes.model, planes.complete, planes.logical)
            );
            assert_eq!(reference.model, a.model);
            assert_eq!(
                (a.model, a.complete, a.logical.clone()),
                (b.model, b.complete, b.logical)
            );
            assert_eq!(a.logical, traced.logical);
            assert_eq!(a.model, traced.model);
            assert_eq!((a.trace, b.trace), (0, 0));
            let single = enumerate_single::<NativeByte64, false, false>(&form, 1 << n);
            assert_eq!(single.model, a.model);
            assert_eq!(single.logical.screen_points, a.logical.screen_points);
            assert_eq!(single.logical.screen_passes, a.logical.screen_passes);
            assert!(a.logical.screen_full_checks <= single.logical.screen_full_checks);
        }
    }
}

#[test]
fn compiled_gray_schedule_matches_looped_models_work_traces_and_every_small_cap() {
    for n in 6..=12 {
        let mut form = SyndromeForm::new(n);
        form.constant = 1 << 31;
        for cap in 0..=(1u64 << n) {
            let a = enumerate_tiered::<NativeByte64, true>(&form, cap);
            let b = enumerate_tiered_unrolled::<NativeByte64, true>(&form, cap);
            assert_eq!(
                (a.model, a.complete, a.logical, a.trace, a.secondary_updates),
                (b.model, b.complete, b.logical, b.trace, b.secondary_updates)
            );
            let a = enumerate_word64(&form, cap);
            let b = enumerate_word64_unrolled(&form, cap);
            assert_eq!(
                (a.model, a.complete, a.points, a.batches),
                (b.model, b.complete, b.points, b.batches)
            );
        }
    }
    for desired in 0..1024u32 {
        let mut form = SyndromeForm::new(10);
        form.constant = desired;
        for i in 0..10 {
            form.linear[i] = 1 << i;
        }
        let a = enumerate_tiered_unrolled::<NativeByte64, false>(&form, 1024);
        assert_eq!(a.model, Some(u64::from(desired)));
        assert_eq!(a.model, enumerate_quiet(&form, 1024).model);
        assert_eq!(a.model, enumerate_word64_unrolled(&form, 1024).model);
    }
    for family in ["planted", "cross_planted", "unplanted"] {
        let (system, _) = fixture(16, 937, family);
        let form = SyndromeForm::from_system(&system, 16).unwrap();
        let a = enumerate_tiered::<NativeByte64, true>(&form, 1 << 16);
        let b = enumerate_tiered_unrolled::<NativeByte64, true>(&form, 1 << 16);
        assert_eq!(
            (a.model, a.complete, a.logical, a.trace),
            (b.model, b.complete, b.logical, b.trace)
        );
        let a = solve_leaf_scan_engine::<true, 16, false, 3>(&system, 16, 200000);
        let b = solve_leaf_scan_engine::<true, 16, false, 6>(&system, 16, 200000);
        assert_eq!(
            (a.outcome, a.logical, a.trace),
            (b.outcome, b.logical, b.trace)
        );
    }
}
#[test]
fn byte_and_word_caps_are_censored_with_actual_block_accounting() {
    for n in 0..=12 {
        let mut form = SyndromeForm::new(n);
        form.constant = 1 << 31;
        let size = if n < 4 {
            1
        } else if n < 6 {
            16
        } else {
            64
        };
        for cap in [
            0,
            1,
            size - 1,
            size,
            size + 1,
            (1 << n) - 1,
            1 << n,
            (1 << n) + 1,
        ] {
            let a = enumerate_bytes::<ScalarByte64>(&form, cap);
            let b = enumerate_bytes::<NativeByte64>(&form, cap);
            assert_eq!(
                (a.model, a.complete, a.logical.clone(), a.trace),
                (b.model, b.complete, b.logical.clone(), b.trace)
            );
            let points = a.logical.screen_points + a.logical.screen_fallback_points;
            assert!(points <= cap);
            assert_eq!(a.complete, points == 1 << n);
            assert!(a.model.is_none());
            let word = enumerate_word64(&form, cap);
            assert_eq!(word.points, points);
            assert_eq!(word.complete, a.complete);
        }
    }
}

#[test]
fn deferred_coefficients_match_eager_updates_after_arbitrary_gaps() {
    let mut rng = 743987;
    for n in 6..=12 {
        let mut form = SyndromeForm::new(n);
        form.constant = next(&mut rng) as u32;
        for i in 0..n {
            form.linear[i] = next(&mut rng) as u32;
            for j in i + 1..n {
                form.quadratic[i][j] = next(&mut rng) as u32;
                form.quadratic[j][i] = form.quadratic[i][j];
            }
        }
        for period in [1, 2, 3, 7, 13] {
            let mut eager = SixCursor::new(&form);
            let mut lazy = SixCursor::new(&form);
            let mut prior = 0;
            for step in 0..1u64 << (n - 6) {
                eager.advance::<true>(&form, step);
                lazy.advance::<false>(&form, step);
                if step % period == 0 || step + 1 == 1 << (n - 6) {
                    let high = step ^ (step >> 1);
                    let mut changed = high ^ prior;
                    while changed != 0 {
                        let j = changed.trailing_zeros() as usize;
                        changed &= changed - 1;
                        for i in 0..6 {
                            lazy.linear[i] ^= lazy.cross[j][i];
                        }
                    }
                    prior = high;
                    assert_eq!(eager.constant, lazy.constant);
                    assert_eq!(eager.linear, lazy.linear);
                }
            }
        }
        for cap in [0, 1, 63, 64, 65, (1 << n) - 1, 1 << n] {
            let a = enumerate_byte_engine::<ScalarByte64, false>(&form, cap);
            let b = enumerate_byte_engine::<ScalarByte64, true>(&form, cap);
            assert_eq!(
                (a.model, a.complete, a.logical, a.trace),
                (b.model, b.complete, b.logical, b.trace)
            );
            assert!(b.linear_updates <= a.linear_updates);
        }
    }
}
#[test]
fn quiet_and_byte_solvers_preserve_complete_models_and_domains() {
    for n in [8, 12, 16] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let old = solve_gray(&system, n, 200000, true);
                let quiet = solve_gray_control(&system, n, 200000, false);
                assert_eq!(old.outcome, quiet.outcome);
                assert_eq!(old.logical, quiet.logical);
                let wide = solve_gray_control(&system, n, 200000, true);
                assert_eq!(old.outcome, wide.outcome);
                let a = solve_byte(&system, n, 200000, false);
                let b = solve_byte(&system, n, 200000, true);
                assert_eq!(old.outcome, a.outcome);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                if let Outcome::Sat(point) = a.outcome {
                    assert!(satisfies(&system, point));
                }
                if a.outcome == Outcome::Unsat {
                    assert_eq!(a.logical.screen_points, 1 << n);
                    assert_eq!(a.logical.screen_second_checks, a.logical.screen_passes);
                    assert_eq!(a.logical.screen_full_rejected, a.logical.screen_full_checks);
                    assert_eq!(
                        a.logical.screen_second_rejected + a.logical.screen_full_checks,
                        a.logical.screen_second_checks
                    );
                }
                assert!(matches!(
                    solve_byte(&system, n, 0, true).outcome,
                    Outcome::Unknown("SCREEN_CAP")
                ));
            }
        }
    }
    for n in 0..6 {
        for system in [vec![], vec![vec![0]]] {
            let a = solve_byte(&system, n, 200000, false);
            let b = solve_byte(&system, n, 200000, true);
            assert_eq!(a.outcome, brute(&system, n));
            assert_eq!(a.outcome, b.outcome);
            assert_eq!(a.logical, b.logical);
            assert_eq!(a.logical.screen_points, 0);
        }
    }
}
#[test]
fn partial_leaves_preserve_original_prefixes_and_recovered_models() {
    for n in [12, 16, 20] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let old = solve_leaf_backend::<true, 16>(&system, n, 200000);
                let quiet = solve_leaf_scan_engine::<true, 16, false, 1>(&system, n, 200000);
                assert_eq!(old.outcome, quiet.outcome);
                assert_eq!(old.logical, quiet.logical);
                let a = solve_leaf_scan_engine::<true, 16, false, 2>(&system, n, 200000);
                let b = solve_leaf_scan_engine::<true, 16, false, 3>(&system, n, 200000);
                assert_eq!(old.outcome, a.outcome);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
                assert_eq!(
                    (
                        old.logical.nodes,
                        old.logical.decisions,
                        old.logical.forced,
                        old.logical.specialized_terms,
                        old.logical.max_depth,
                        old.logical.enumeration_leaves
                    ),
                    (
                        a.logical.nodes,
                        a.logical.decisions,
                        a.logical.forced,
                        a.logical.specialized_terms,
                        a.logical.max_depth,
                        a.logical.enumeration_leaves
                    )
                );
                assert_eq!(a.logical.screen_calls, a.logical.enumeration_leaves);
                if let Outcome::Sat(point) = a.outcome {
                    assert!(satisfies(&system, point));
                }
            }
        }
    }
}
