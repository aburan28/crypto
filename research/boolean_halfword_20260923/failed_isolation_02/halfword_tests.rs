#[cfg(test)]
mod halfword_tests {
    use super::*;
    fn random_system(n: u8, equations: usize, seed: &mut u64) -> System {
        (0..equations)
            .map(|_| {
                (0..1u64 << n)
                    .filter(|m| m.count_ones() <= 2 && next(seed) & 3 == 0)
                    .collect()
            })
            .collect()
    }
    fn direct(system: &System, n: u8, width: usize, cap: u64) -> HalfEnumeration {
        let full = SyndromeForm::from_system(system, n).unwrap();
        let low = if n < 4 {
            0
        } else if width == 64 && n >= 6 {
            6
        } else {
            4
        };
        let size = 1usize << low;
        let mut out = HalfEnumeration {
            model: None,
            complete: false,
            work: HalfWork {
                block_width: size,
                ..HalfWork::default()
            },
        };
        for step in 0..1u64 << (usize::from(n) - low) {
            if cap.saturating_sub(out.work.points) < size as u64 {
                return out;
            }
            out.work.points += size as u64;
            out.work.batches += 1;
            for y in 0..size {
                // Express the 64-point order as four original 16-point Gray
                // blocks, independently of legacy_group64 and the new schedule.
                let point = if low == 0 {
                    step
                } else {
                    let block = if low == 6 {
                        4 * step + (y / 16) as u64
                    } else {
                        step
                    };
                    ((block ^ (block >> 1)) << 4) | (y % 16) as u64
                };
                let value = full.value(point);
                if value as u16 == 0 {
                    out.work.projected_hits += 1;
                    out.work.original_checks += 1;
                    if value == 0 {
                        out.work.verified_hits += 1;
                        out.model = Some(point);
                        out.complete = true;
                        return out;
                    }
                    out.work.rejected_hits += 1;
                }
            }
        }
        out.complete = true;
        out
    }
    fn equal(a: &HalfEnumeration, b: &HalfEnumeration) {
        assert_eq!(
            (a.model, a.complete, &a.work),
            (b.model, b.complete, &b.work)
        );
    }
    fn block_check<B: HalfBlock>() {
        for lane in 0..B::POINTS {
            let mut values = [0x8000; 64];
            assert!(!B::from_values(&values).any_zero());
            values[lane] = 0;
            let mut block = B::from_values(&values);
            assert!(block.any_zero());
            assert_eq!(&block.values()[..B::POINTS], &values[..B::POINTS]);
            block.xor_uniform(0x8000);
            for v in &mut values[..B::POINTS] {
                *v ^= 0x8000;
            }
            assert_eq!(&block.values()[..B::POINTS], &values[..B::POINTS]);
            let other = std::array::from_fn(|i| (i as u16).wrapping_mul(1973));
            block.xor_block(&B::from_values(&other));
            for i in 0..B::POINTS {
                values[i] ^= other[i];
            }
            assert_eq!(&block.values()[..B::POINTS], &values[..B::POINTS]);
            assert_eq!(block.any_zero(), values[..B::POINTS].contains(&0));
        }
    }
    #[test]
    fn native_and_scalar_blocks_preserve_every_lane_and_unsigned_high_bit() {
        block_check::<ScalarHalf<16>>();
        block_check::<ScalarHalf<64>>();
        block_check::<NativeHalf<2>>();
        block_check::<NativeHalf<8>>();
    }
    #[test]
    fn fused_encoding_matches_full_projection_and_rejects_unsupported_omitted_rows() {
        let mut seed = 93701;
        for n in 0..=8 {
            for equations in [0, 1, 15, 16, 17, 32] {
                let system = random_system(n, equations, &mut seed);
                let full = SyndromeForm::from_system(&system, n).unwrap();
                let projected = HalfForm::from_full(&full);
                let fused = HalfForm::from_system(&system, n).unwrap();
                for point in 0..1u64 << n {
                    assert_eq!(fused.value(point), full.value(point) as u16);
                    assert_eq!(fused.value(point), projected.value(point));
                }
            }
        }
        let mut unsupported = vec![vec![]; 17];
        unsupported[16] = vec![7];
        assert!(HalfForm::from_system(&unsupported, 3).is_none());
        unsupported[16] = vec![8];
        assert!(HalfForm::from_system(&unsupported, 3).is_none());
        assert!(HalfForm::from_system(&vec![vec![]; 33], 3).is_none());
        assert!(HalfForm::from_system(&vec![], 25).is_none());
    }
    fn schedule_check<B: HalfBlock, const LOW: usize>(system: &System, n: u8) {
        let full = SyndromeForm::from_system(system, n).unwrap();
        let half = HalfForm::from_system(system, n).unwrap();
        let mut schedule = HalfSchedule::<B, LOW>::new(&half);
        for step in 0..1u64 << (usize::from(n) - LOW) {
            schedule.advance(&half, step);
            let values = schedule.block.values();
            for (low, &value) in values[..B::POINTS].iter().enumerate() {
                let point = ((step ^ (step >> 1)) << LOW) | low as u64;
                assert_eq!(value, full.value(point) as u16);
            }
        }
    }
    #[test]
    fn every_initial_block_and_gray_transition_matches_direct_equations() {
        let mut seed = 19791;
        for n in [8, 10, 12] {
            for equations in [16, 17, 32] {
                let system = random_system(n, equations, &mut seed);
                schedule_check::<NativeHalf<2>, 4>(&system, n);
                schedule_check::<NativeHalf<8>, 6>(&system, n);
                schedule_check::<ScalarHalf<16>, 4>(&system, n);
                schedule_check::<ScalarHalf<64>, 6>(&system, n);
            }
        }
    }
    #[test]
    fn sixteen_seventeen_equation_boundary_rejects_false_and_multiple_hits() {
        let mut system = vec![vec![]; 17];
        system[16] = vec![1, 0];
        for arm in [
            "half16_scalar",
            "half64_scalar",
            "half16_native",
            "half64_native",
            "half_dispatch",
        ] {
            let (got, work) = solve_half(&system, 8, 1, arm);
            assert_eq!(got.outcome, Outcome::Sat(1));
            assert_eq!(
                (
                    work.projected_hits,
                    work.original_checks,
                    work.rejected_hits,
                    work.verified_hits
                ),
                (2, 2, 1, 1)
            );
        }
        system[16] = vec![0];
        for arm in ["half16_native", "half64_native"] {
            let (got, work) = solve_half(&system, 8, 1, arm);
            assert_eq!(got.outcome, Outcome::Unsat);
            assert_eq!(work.rejected_hits, 256);
            assert_eq!(work.verified_hits, 0);
        }
        system.truncate(16);
        system[15] = vec![0];
        let (got, work) = solve_half(&system, 8, 1, "half64_native");
        assert_eq!(got.outcome, Outcome::Unsat);
        assert_eq!(work.projected_hits, 0);
    }
    #[test]
    fn every_small_cap_is_censored_with_exact_work_and_no_skipped_false_hits() {
        let mut system = vec![vec![]; 17];
        system[16] = vec![0];
        for n in 0..=7 {
            let form = HalfForm::from_system(&system, n).unwrap();
            for cap in 0..=(1u64 << n) + 1 {
                let a = direct(&system, n, 16, cap);
                let b = direct(&system, n, 64, cap);
                equal(&enumerate_half::<NativeHalf<2>, 4>(&form, &system, cap), &a);
                equal(
                    &enumerate_half::<ScalarHalf<16>, 4>(&form, &system, cap),
                    &a,
                );
                equal(&enumerate_half::<NativeHalf<8>, 6>(&form, &system, cap), &b);
                equal(
                    &enumerate_half::<ScalarHalf<64>, 6>(&form, &system, cap),
                    &b,
                );
            }
        }
    }
    #[test]
    fn every_late_model_preserves_the_original_reflected_assignment_order() {
        for model in 0..256u64 {
            let mut system = vec![vec![]; 16];
            for i in 0..8 {
                system.push(if model & (1 << i) == 0 {
                    vec![1 << i]
                } else {
                    vec![1 << i, 0]
                });
            }
            let form = HalfForm::from_system(&system, 8).unwrap();
            let got16 = enumerate_half::<NativeHalf<2>, 4>(&form, &system, 256);
            let got64 = enumerate_half::<NativeHalf<8>, 6>(&form, &system, 256);
            equal(&got16, &direct(&system, 8, 16, 256));
            equal(&got64, &direct(&system, 8, 64, 256));
            assert_eq!(got16.model, Some(model));
            assert_eq!(got64.model, Some(model));
        }
    }
    #[test]
    fn every_pair_of_small_polynomials_matches_complete_truth_tables() {
        let monomials = [0, 1, 2, 4, 3, 5, 6];
        let poly = |code: usize| {
            monomials
                .iter()
                .enumerate()
                .filter_map(|(i, &m)| (code & (1 << i) != 0).then_some(m))
                .collect::<Vec<_>>()
        };
        for a in 0..128 {
            for b in 0..128 {
                let system = vec![poly(a), poly(b)];
                let form = HalfForm::from_system(&system, 3).unwrap();
                let expected = direct(&system, 3, 16, 8);
                equal(
                    &enumerate_half::<NativeHalf<2>, 4>(&form, &system, 8),
                    &expected,
                );
                equal(
                    &enumerate_half::<NativeHalf<8>, 6>(&form, &system, 8),
                    &expected,
                );
            }
        }
    }
    #[test]
    fn compiled_scalar_and_native_paths_match_original_models_and_block_counts() {
        for seed in [17, 937, 8177] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(12, seed, family);
                for arm in [
                    "half16_scalar",
                    "half64_scalar",
                    "half16_native",
                    "half64_native",
                    "half_dispatch",
                ] {
                    let (got, work) = solve_half(&system, 12, 200000, arm);
                    let reference = solve(
                        &system,
                        12,
                        if arm.contains("64") {
                            "wide64_unrolled"
                        } else {
                            "word16_unrolled"
                        },
                        200000,
                    );
                    assert_eq!(got.outcome, reference.outcome);
                    assert_eq!(got.logical, reference.logical);
                    assert_eq!(work.projected_hits, work.original_checks);
                    assert_eq!(
                        work.original_checks,
                        work.rejected_hits + work.verified_hits
                    );
                }
            }
        }
    }
}
