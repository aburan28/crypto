// Discovery-only cost partition. Builders and scans are shared by both modes.
struct Prepared16 {
    cursor: Word16Schedule,
    low: [NativeDeltaBlock; 4],
}
struct Prepared64 {
    cursor: SixCursor,
    block: Word64,
    cross: Vec<Word64>,
    low_cross: [Word64; 4],
}
#[derive(Default, Clone, Debug)]
struct ConstructionPhases {
    encoding_ns: u128,
    construction_ns: u128,
    scan_ns: u128,
    release_ns: u128,
}
#[inline(never)]
fn prepare16(form: &SyndromeForm) -> Prepared16 {
    let cursor = Word16Schedule::new(form);
    let low = std::array::from_fn(|j| cursor.cross[j]);
    Prepared16 { cursor, low }
}
#[inline(never)]
fn prepare64(form: &SyndromeForm) -> Prepared64 {
    let cursor = SixCursor::new(form);
    let offsets = low_quadratic_offsets64(form);
    let mut block = Word64::affine(form.constant, &cursor.linear);
    block.xor_block(&Word64::from_values(&offsets));
    let cross: Vec<_> = cursor.cross.iter().map(|a| Word64::affine(0, a)).collect();
    let low_cross = std::array::from_fn(|j| cross[j]);
    Prepared64 {
        cursor,
        block,
        cross,
        low_cross,
    }
}
fn construction_solve<const PROFILE: bool>(
    system: &System,
    n: u8,
    limit: u64,
    width: usize,
) -> (Solved, Option<ConstructionPhases>) {
    assert!((12..=24).contains(&n) && [16, 64].contains(&width));
    let mut phases = ConstructionPhases::default();
    let mut tick = if PROFILE { Some(Instant::now()) } else { None };
    let form = SyndromeForm::from_system(system, n)
        .expect("Diagnostic domain is quadratic and has at most 32 equations");
    if let Some(start) = tick {
        let now = Instant::now();
        phases.encoding_ns = now.duration_since(start).as_nanos();
        tick = Some(now);
    }
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let answer;
    if width == 16 {
        let mut prepared = prepare16(&form);
        if let Some(start) = tick {
            let now = Instant::now();
            phases.construction_ns = now.duration_since(start).as_nanos();
            tick = Some(now);
        }
        answer = scan16_shared(&form, &mut prepared, cap);
        if let Some(start) = tick {
            let now = Instant::now();
            phases.scan_ns = now.duration_since(start).as_nanos();
            tick = Some(now);
        }
        drop(prepared);
    } else {
        let mut prepared = prepare64(&form);
        if let Some(start) = tick {
            let now = Instant::now();
            phases.construction_ns = now.duration_since(start).as_nanos();
            tick = Some(now);
        }
        answer = scan64_shared(&form, &mut prepared, cap);
        if let Some(start) = tick {
            let now = Instant::now();
            phases.scan_ns = now.duration_since(start).as_nanos();
            tick = Some(now);
        }
        drop(prepared);
    }
    drop(form);
    if let Some(start) = tick {
        phases.release_ns = start.elapsed().as_nanos();
    }
    let solved = Solved {
        outcome: if !answer.complete {
            Outcome::Unknown("ENUM_CAP")
        } else {
            answer.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
        logical: Logical {
            enumeration_points: answer.points,
            enumeration_batches: answer.batches,
            enumeration_leaves: 1,
            ..Logical::default()
        },
        profile: Profile::default(),
        trace: 0,
    };
    (solved, if PROFILE { Some(phases) } else { None })
}

#[inline(never)]
fn scan16_shared(form: &SyndromeForm, prepared: &mut Prepared16, cap: u64) -> Enumeration {
    let cursor = &mut prepared.cursor;
    let low = &prepared.low;
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
        cursor.advance(form, base);
        if process_word16(&cursor.block, &mut out, base) {
            return out;
        }
        macro_rules! step {
            ($j:literal,$r:literal) => {{
                if cap - out.points < 16 {
                    return out;
                }
                cursor.fixed::<$j>(form, &low[$j]);
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

#[inline(never)]
fn scan64_shared(form: &SyndromeForm, prepared: &mut Prepared64, cap: u64) -> Enumeration {
    let cursor = &mut prepared.cursor;
    let block = &mut prepared.block;
    let cross = &prepared.cross;
    let low_cross = &prepared.low_cross;
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
            block.xor_uniform(d);
            block.xor_block(&cross[j]);
        }
        if process_word64(&block, &mut out, base) {
            return out;
        }
        macro_rules! next_step {
            ($j:literal,$r:literal) => {{
                if cap - out.points < 64 {
                    return out;
                }
                let d = cursor.advance_fixed::<$j>(form);
                block.xor_uniform(d);
                block.xor_block(&low_cross[$j]);
                if process_word64(&block, &mut out, base + $r) {
                    return out;
                }
            }};
        }
        next_step!(0, 1);
        next_step!(1, 2);
        next_step!(0, 3);
        next_step!(2, 4);
        next_step!(0, 5);
        next_step!(1, 6);
        next_step!(0, 7);
        next_step!(3, 8);
        next_step!(0, 9);
        next_step!(1, 10);
        next_step!(0, 11);
        next_step!(2, 12);
        next_step!(0, 13);
        next_step!(1, 14);
        next_step!(0, 15);
    }
    out.complete = true;
    out
}

#[cfg(test)]
mod construction_tests {
    use super::*;
    fn signature(v: &Enumeration) -> (Option<u64>, bool, u64, u64, u64) {
        (v.model, v.complete, v.points, v.batches, v.checksum)
    }
    #[test]
    fn extracted_scans_preserve_every_small_cap() {
        let mut form = SyndromeForm::new(12);
        form.constant = 1 << 31;
        for cap in 0..=4096 {
            let mut a = prepare16(&form);
            let mut b = prepare64(&form);
            assert_eq!(
                signature(&scan16_shared(&form, &mut a, cap)),
                signature(&enumerate_word16_unrolled(&form, cap))
            );
            assert_eq!(
                signature(&scan64_shared(&form, &mut b, cap)),
                signature(&enumerate_word64_unrolled(&form, cap))
            );
        }
    }
    #[test]
    fn profiled_and_rebound_paths_preserve_complete_models_and_work() {
        for n in [12, 16, 20, 24] {
            for seed in [17, 937] {
                for family in ["planted", "cross_planted", "unplanted"] {
                    let (system, _) = fixture(n, seed, family);
                    for width in [16, 64] {
                        let old = solve(
                            &system,
                            n,
                            if width == 16 {
                                "word16_unrolled"
                            } else {
                                "wide64_unrolled"
                            },
                            200000,
                        );
                        let (plain, absent) =
                            construction_solve::<false>(&system, n, 200000, width);
                        assert!(absent.is_none());
                        let tick = Instant::now();
                        let (measured, phases) =
                            construction_solve::<true>(&system, n, 200000, width);
                        let elapsed = tick.elapsed().as_nanos();
                        let phases = phases.unwrap();
                        assert!(
                            phases.encoding_ns
                                + phases.construction_ns
                                + phases.scan_ns
                                + phases.release_ns
                                <= elapsed
                        );
                        assert_eq!(old.outcome, plain.outcome);
                        assert_eq!(old.logical, plain.logical);
                        assert_eq!(old.trace, plain.trace);
                        assert_eq!(plain.outcome, measured.outcome);
                        assert_eq!(plain.logical, measured.logical);
                        assert_eq!(plain.trace, measured.trace);
                        if let Outcome::Sat(point) = measured.outcome {
                            assert!(satisfies(&system, point));
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn prepared_initial_blocks_match_direct_equation_values() {
        let mut rng = 612971;
        for n in [12, 16, 20, 24] {
            let mut form = SyndromeForm::new(n);
            form.constant = next(&mut rng) as u32;
            for i in 0..n {
                form.linear[i] = next(&mut rng) as u32;
                for j in i + 1..n {
                    form.quadratic[i][j] = next(&mut rng) as u32;
                    form.quadratic[j][i] = form.quadratic[i][j];
                }
            }
            let a = prepare16(&form);
            let b = prepare64(&form);
            for y in 0..16 {
                assert_eq!(a.cursor.block.values()[y / 4][y % 4], form.value(y as u64));
            }
            for y in 0..64 {
                assert_eq!(
                    b.block.0[y / 16].values()[(y % 16) / 4][y % 4],
                    form.value(y as u64)
                );
            }
        }
    }
}
