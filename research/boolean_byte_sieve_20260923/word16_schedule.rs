// Compiled four-bit Gray transition schedule for the retained 16-point domain.
struct Word16Schedule {
    block: NativeDeltaBlock,
    differences: [u32; MAX],
    cross: Vec<NativeDeltaBlock>,
}
impl Word16Schedule {
    fn new(form: &SyndromeForm) -> Self {
        let offsets = form.low_offsets();
        let mut values = offsets;
        for y in 0..16 {
            values[y / 4][y % 4] ^= form.constant;
            for i in 0..4 {
                if y & (1 << i) != 0 {
                    values[y / 4][y % 4] ^= form.linear[i];
                }
            }
        }
        let mut differences = [0; MAX];
        let cross = (0..form.n - 4)
            .map(|j| {
                differences[j] = form.linear[j + 4]
                    ^ if j == 0 {
                        0
                    } else {
                        form.quadratic[j + 4][j + 3]
                    };
                let mut flat = [0u32; 16];
                for y in 1usize..16 {
                    let bit = y & y.wrapping_neg();
                    flat[y] = flat[y ^ bit] ^ form.quadratic[bit.trailing_zeros() as usize][j + 4];
                }
                NativeDeltaBlock::from_values(&std::array::from_fn(|r| {
                    std::array::from_fn(|i| flat[r * 4 + i])
                }))
            })
            .collect();
        Self {
            block: NativeDeltaBlock::from_values(&values),
            differences,
            cross,
        }
    }
    #[inline(always)]
    fn advance(&mut self, form: &SyndromeForm, step: u64) {
        if step == 0 {
            return;
        }
        let j = step.trailing_zeros() as usize;
        self.block.xor_uniform(self.differences[j]);
        self.block.xor_block(&self.cross[j]);
        for i in 0..j {
            self.differences[i] ^= form.quadratic[i + 4][j + 4];
        }
    }
    #[inline(always)]
    fn fixed<const J: usize>(&mut self, form: &SyndromeForm, cross: &NativeDeltaBlock) {
        self.block.xor_uniform(self.differences[J]);
        self.block.xor_block(cross);
        for i in 0..J {
            self.differences[i] ^= form.quadratic[i + 4][J + 4];
        }
    }
}
#[inline(always)]
fn process_word16(block: &NativeDeltaBlock, out: &mut Enumeration, step: u64) -> bool {
    out.points += 16;
    out.batches += 1;
    if let Some(low) = native_first_zero16(block) {
        out.model = Some(((step ^ (step >> 1)) << 4) | u64::from(low));
        out.complete = true;
        true
    } else {
        false
    }
}
fn enumerate_word16_unrolled(form: &SyndromeForm, cap: u64) -> Enumeration {
    if form.n < 8 {
        return enumerate_quiet(form, cap);
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
fn solve_word_schedule(system: &System, n: u8, limit: u64, dispatch: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let got = if dispatch && n > 20 {
        enumerate_word64_unrolled(&form, cap)
    } else {
        enumerate_word16_unrolled(&form, cap)
    };
    Solved {
        outcome: if !got.complete {
            Outcome::Unknown("ENUM_CAP")
        } else {
            got.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
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
