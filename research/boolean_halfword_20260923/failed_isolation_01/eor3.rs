// Three-input XOR changes only the implementation of block ^= delta ^ cross.
// Feature detection is outside the unsafe kernels. Unsupported hosts use the
// unchanged two-XOR policy; the producer records the selected capability.
fn eor3_available() -> bool {
    #[cfg(target_arch = "aarch64")]
    {
        std::arch::is_aarch64_feature_detected!("sha3")
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        false
    }
}
fn half_eor3<const G: usize, const LOW: usize>(
    form: &HalfForm,
    system: &System,
    cap: u64,
) -> HalfEnumeration {
    #[cfg(target_arch = "aarch64")]
    if eor3_available() {
        // The runtime flag is required before entering this target-feature body.
        return unsafe { half_eor3_kernel::<G, LOW>(form, system, cap) };
    }
    enumerate_half::<NativeHalf<G>, LOW>(form, system, cap)
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "sha3")]
unsafe fn half_eor3_kernel<const G: usize, const LOW: usize>(
    form: &HalfForm,
    system: &System,
    cap: u64,
) -> HalfEnumeration {
    use std::arch::aarch64::*;
    if form.n < LOW + 4 {
        return enumerate_half::<NativeHalf<G>, LOW>(form, system, cap);
    }
    let mut schedule = HalfSchedule::<NativeHalf<G>, LOW>::new(form);
    let low: [NativeHalf<G>; 4] = std::array::from_fn(|j| schedule.cross[j]);
    let mut out = HalfEnumeration {
        model: None,
        complete: false,
        work: HalfWork {
            block_width: 1 << LOW,
            ..HalfWork::default()
        },
    };
    for base in (0..1u64 << (form.n - LOW)).step_by(16) {
        if cap - out.work.points < 1u64 << LOW {
            return out;
        }
        if base != 0 {
            let j = base.trailing_zeros() as usize;
            let delta = vdupq_n_u16(schedule.differences[j]);
            for g in 0..G {
                schedule.block.0[g] =
                    veor3q_u16(schedule.block.0[g], schedule.cross[j].0[g], delta);
            }
            for i in 0..j {
                schedule.differences[i] ^= form.quadratic[i + LOW][j + LOW];
            }
        }
        if half_process::<NativeHalf<G>, LOW>(&schedule.block, system, &mut out, base) {
            return out;
        }
        macro_rules! step {
            ($j:literal,$r:literal) => {{
                if cap - out.work.points < 1u64 << LOW {
                    return out;
                }
                let delta = vdupq_n_u16(schedule.differences[$j]);
                for g in 0..G {
                    schedule.block.0[g] = veor3q_u16(schedule.block.0[g], low[$j].0[g], delta);
                }
                for i in 0..$j {
                    schedule.differences[i] ^= form.quadratic[i + LOW][$j + LOW];
                }
                if half_process::<NativeHalf<G>, LOW>(&schedule.block, system, &mut out, base + $r)
                {
                    return out;
                }
            }};
        }
        step!(0, 1);
        step!(1, 2);
        step!(0, 3);
        step!(2, 4);
        step!(0, 5);
        step!(1, 6);
        step!(0, 7);
        step!(3, 8);
        step!(0, 9);
        step!(1, 10);
        step!(0, 11);
        step!(2, 12);
        step!(0, 13);
        step!(1, 14);
        step!(0, 15);
    }
    out.complete = true;
    out
}
fn word_eor3(form: &SyndromeForm, cap: u64, wide: bool) -> Enumeration {
    #[cfg(target_arch = "aarch64")]
    if eor3_available() {
        return unsafe {
            if wide {
                word64_eor3_kernel(form, cap)
            } else {
                word16_eor3_kernel(form, cap)
            }
        };
    }
    if wide {
        enumerate_word64_unrolled(form, cap)
    } else {
        enumerate_word16_unrolled(form, cap)
    }
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "sha3")]
unsafe fn word16_eor3_kernel(form: &SyndromeForm, cap: u64) -> Enumeration {
    use std::arch::aarch64::*;
    if form.n < 8 {
        return enumerate_word16_unrolled(form, cap);
    }
    let mut cursor = Word16Schedule::new(form);
    let low: [NativeDeltaBlock; 4] = std::array::from_fn(|j| cursor.cross[j]);
    let mut out = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0,
    };
    for base in (0..1u64 << (form.n - 4)).step_by(16) {
        if cap - out.points < 16 {
            return out;
        }
        if base != 0 {
            let j = base.trailing_zeros() as usize;
            let delta = vdupq_n_u32(cursor.differences[j]);
            for g in 0..4 {
                cursor.block.0[g] = veor3q_u32(cursor.block.0[g], cursor.cross[j].0[g], delta);
            }
            for i in 0..j {
                cursor.differences[i] ^= form.quadratic[i + 4][j + 4];
            }
        }
        if process_word16(&cursor.block, &mut out, base) {
            return out;
        }
        macro_rules! step {
            ($j:literal,$r:literal) => {{
                if cap - out.points < 16 {
                    return out;
                }
                let delta = vdupq_n_u32(cursor.differences[$j]);
                for g in 0..4 {
                    cursor.block.0[g] = veor3q_u32(cursor.block.0[g], low[$j].0[g], delta);
                }
                for i in 0..$j {
                    cursor.differences[i] ^= form.quadratic[i + 4][$j + 4];
                }
                if process_word16(&cursor.block, &mut out, base + $r) {
                    return out;
                }
            }};
        }
        step!(0, 1);
        step!(1, 2);
        step!(0, 3);
        step!(2, 4);
        step!(0, 5);
        step!(1, 6);
        step!(0, 7);
        step!(3, 8);
        step!(0, 9);
        step!(1, 10);
        step!(0, 11);
        step!(2, 12);
        step!(0, 13);
        step!(1, 14);
        step!(0, 15);
    }
    out.complete = true;
    out
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "sha3")]
unsafe fn word64_eor3_kernel(form: &SyndromeForm, cap: u64) -> Enumeration {
    use std::arch::aarch64::*;
    if form.n < 10 {
        return enumerate_word64_unrolled(form, cap);
    }
    let mut cursor = SixCursor::new(form);
    let offsets = low_quadratic_offsets64(form);
    let mut block = Word64::affine(form.constant, &cursor.linear);
    block.xor_block(&Word64::from_values(&offsets));
    let cross: Vec<_> = cursor.cross.iter().map(|a| Word64::affine(0, a)).collect();
    let low: [Word64; 4] = std::array::from_fn(|j| cross[j]);
    let mut out = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0,
    };
    for base in (0..1u64 << (form.n - 6)).step_by(16) {
        if cap - out.points < 64 {
            return out;
        }
        if let Some((j, d)) = cursor.advance::<false>(form, base) {
            let delta = vdupq_n_u32(d);
            for g in 0..4 {
                for v in 0..4 {
                    block.0[g].0[v] = veor3q_u32(block.0[g].0[v], cross[j].0[g].0[v], delta);
                }
            }
        }
        if process_word64(&block, &mut out, base) {
            return out;
        }
        macro_rules! step {
            ($j:literal,$r:literal) => {{
                if cap - out.points < 64 {
                    return out;
                }
                let delta = vdupq_n_u32(cursor.advance_fixed::<$j>(form));
                for g in 0..4 {
                    for v in 0..4 {
                        block.0[g].0[v] = veor3q_u32(block.0[g].0[v], low[$j].0[g].0[v], delta);
                    }
                }
                if process_word64(&block, &mut out, base + $r) {
                    return out;
                }
            }};
        }
        step!(0, 1);
        step!(1, 2);
        step!(0, 3);
        step!(2, 4);
        step!(0, 5);
        step!(1, 6);
        step!(0, 7);
        step!(3, 8);
        step!(0, 9);
        step!(1, 10);
        step!(0, 11);
        step!(2, 12);
        step!(0, 13);
        step!(1, 14);
        step!(0, 15);
    }
    out.complete = true;
    out
}
fn solve_word_eor3(system: &System, n: u8, limit: u64, wide: bool) -> Solved {
    if n > 24 {
        return Solved {
            outcome: Outcome::Unknown("EOR3_DOMAIN"),
            logical: Logical::default(),
            profile: Profile::default(),
            trace: 0,
        };
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return Solved {
            outcome: Outcome::Unknown("EOR3_DOMAIN"),
            logical: Logical::default(),
            profile: Profile::default(),
            trace: 0,
        };
    };
    let got = word_eor3(&form, if limit == 0 { 0 } else { 1 << 24 }, wide);
    let outcome = if !got.complete {
        Outcome::Unknown("ENUM_CAP")
    } else {
        got.model.map_or(Outcome::Unsat, Outcome::Sat)
    };
    Solved {
        outcome,
        logical: Logical {
            enumeration_points: got.points,
            enumeration_batches: got.batches,
            enumeration_leaves: 1,
            ..Logical::default()
        },
        profile: Profile {
            enumeration_ns: tick.elapsed().as_nanos(),
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: 0,
    }
}

#[cfg(test)]
mod eor3_tests {
    use super::*;
    #[test]
    fn fused_and_portable_kernels_have_identical_models_work_and_caps() {
        let mut seed = 41779;
        for n in [4, 6, 8, 10, 12] {
            let system: System = (0..20)
                .map(|_| {
                    (0..1u64 << n)
                        .filter(|m| m.count_ones() <= 2 && next(&mut seed) & 7 == 0)
                        .collect()
                })
                .collect();
            let full = SyndromeForm::from_system(&system, n).unwrap();
            let half = HalfForm::from_system(&system, n).unwrap();
            for cap in [
                0,
                1,
                15,
                16,
                17,
                63,
                64,
                65,
                127,
                128,
                255,
                256,
                (1u64 << n) - 1,
                1u64 << n,
            ] {
                for wide in [false, true] {
                    let a = word_eor3(&full, cap, wide);
                    let b = if wide {
                        enumerate_word64_unrolled(&full, cap)
                    } else {
                        enumerate_word16_unrolled(&full, cap)
                    };
                    assert_eq!(
                        (a.model, a.complete, a.points, a.batches),
                        (b.model, b.complete, b.points, b.batches)
                    );
                }
                let a = half_eor3::<2, 4>(&half, &system, cap);
                let b = enumerate_half::<NativeHalf<2>, 4>(&half, &system, cap);
                assert_eq!((a.model, a.complete, a.work), (b.model, b.complete, b.work));
                let a = half_eor3::<8, 6>(&half, &system, cap);
                let b = enumerate_half::<NativeHalf<8>, 6>(&half, &system, cap);
                assert_eq!((a.model, a.complete, a.work), (b.model, b.complete, b.work));
            }
        }
    }
    #[test]
    fn fused_projection_keeps_all_false_hits_and_late_true_hits() {
        for target in [0, 1, 47, 64, 127, 255, 1023] {
            let mut system = vec![vec![]; 16];
            for i in 0..10 {
                system.push(if target & (1 << i) == 0 {
                    vec![1 << i]
                } else {
                    vec![1 << i, 0]
                });
            }
            let half = HalfForm::from_system(&system, 10).unwrap();
            let got = half_eor3::<8, 6>(&half, &system, 1024);
            let old = enumerate_half::<NativeHalf<8>, 6>(&half, &system, 1024);
            assert_eq!(got.model, Some(target));
            assert_eq!(got.work, old.work);
        }
        let mut system = vec![vec![]; 16];
        system.push(vec![0]);
        let half = HalfForm::from_system(&system, 10).unwrap();
        let got = half_eor3::<8, 6>(&half, &system, 1024);
        assert!(got.complete && got.model.is_none());
        assert_eq!(got.work.rejected_hits, 1024);
    }
}
