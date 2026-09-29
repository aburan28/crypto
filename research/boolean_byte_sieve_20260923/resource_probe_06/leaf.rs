// Retained packed search with a bounded exact enumeration leaf.
struct LeafSolver<const SIMD: bool, const CUT: usize, const DELTA: bool = false, const SCAN: u8 = 0>
{
    model: PackedModel,
    limit: u64,
    logical: Logical,
    profile: Profile,
    trace: u64,
}
impl<const SIMD: bool, const CUT: usize, const DELTA: bool, const SCAN: u8>
    LeafSolver<SIMD, CUT, DELTA, SCAN>
{
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn visit(
        &mut self,
        mut state: PackedState,
        mut known: u64,
        mut values: u64,
        depth: usize,
    ) -> Outcome {
        if self.logical.nodes == self.limit {
            self.event(0xcaf);
            return Outcome::Unknown("NODE_CAP");
        }
        self.logical.nodes += 1;
        self.logical.max_depth = self.logical.max_depth.max(depth);
        self.event(0x100);
        self.event(known);
        self.event(values);
        self.event(u64::from(state.live.count_ones()));
        let mut live = state.live;
        while live != 0 {
            let e = live.trailing_zeros() as usize;
            self.event(u64::from(state.rows[e].count_ones()));
            let mut hash = self.trace;
            self.model
                .terms(&state, e, |m| hash = (hash ^ m).wrapping_mul(0x100000001b3));
            self.trace = hash;
            live &= live - 1;
        }
        loop {
            if state.live == 0 {
                self.event(0x501);
                self.event(values);
                return Outcome::Sat(values);
            }
            let mut live = state.live;
            let mut affine = StackAffine::new();
            let mut all_affine = true;
            let constant = 1u64 << self.model.n;
            while live != 0 {
                let e = live.trailing_zeros() as usize;
                let equation = &self.model.equations[e];
                let row = state.rows[e];
                if row & equation.q_mask == 0 {
                    let low = row >> equation.q_len;
                    if low == constant {
                        self.event(0x500);
                        return Outcome::Unsat;
                    }
                    affine.insert(low);
                } else {
                    all_affine = false;
                }
                live &= live - 1;
            }
            affine.reduce();
            let (mut mask, mut ones) = (0, 0);
            let mut rows = affine.present;
            while rows != 0 {
                let p = rows.trailing_zeros() as usize;
                let row = affine.rows[p];
                if row == constant {
                    self.event(0x502);
                    return Outcome::Unsat;
                }
                let variables = row & (constant - 1);
                if variables.count_ones() == 1 {
                    mask |= variables;
                    if row & constant != 0 {
                        ones |= variables;
                    }
                }
                rows &= rows - 1;
            }
            if all_affine {
                let mut rows = affine.present;
                while rows != 0 {
                    let p = rows.trailing_zeros() as usize;
                    let row = affine.rows[p];
                    if row & constant != 0 {
                        values |= 1u64 << p;
                    }
                    rows &= rows - 1;
                }
                self.event(0x503);
                self.event(values);
                return Outcome::Sat(values);
            }
            if mask != 0 {
                self.logical.forced += u64::from(mask.count_ones());
                self.event(0x200);
                self.event(mask);
                self.event(ones);
                known |= mask;
                values |= ones;
                state = self.model.specialize(state, mask, ones, &mut self.logical);
                continue;
            }
            break;
        }
        let mut counts = [0usize; MAX];
        let mut live = state.live;
        while live != 0 {
            let e = live.trailing_zeros() as usize;
            self.model.terms(&state, e, |mut m| {
                while m != 0 {
                    counts[m.trailing_zeros() as usize] += 1;
                    m &= m - 1;
                }
            });
            live &= live - 1;
        }
        let relevant = (0..self.model.n as usize)
            .filter(|&i| counts[i] > 0)
            .fold(0u64, |m, i| m | (1u64 << i));
        if relevant.count_ones() as usize <= CUT && state.live.count_ones() <= 32 {
            let tick = Instant::now();
            let labels: Vec<_> = (0..self.model.n as usize)
                .filter(|&i| relevant & (1u64 << i) != 0)
                .collect();
            let mut positions = [0u8; MAX];
            for (p, &j) in labels.iter().enumerate() {
                positions[j] = p as u8;
            }
            let mut form = SyndromeForm::new(labels.len());
            let mut live = state.live;
            let mut e = 0;
            while live != 0 {
                let row = live.trailing_zeros() as usize;
                self.model.terms(&state, row, |mut monomial| {
                    let mut local = 0u64;
                    while monomial != 0 {
                        local |= 1u64 << positions[monomial.trailing_zeros() as usize];
                        monomial &= monomial - 1;
                    }
                    form.add_term(e, local);
                });
                e += 1;
                live &= live - 1;
            }
            let got = if SCAN >= 2 {
                let got = if SCAN == 6 {
                    enumerate_tiered_unrolled::<NativeByte64, false>(&form, 1u64 << labels.len())
                } else if SCAN == 5 {
                    enumerate_tiered::<NativePlane64, false>(&form, 1u64 << labels.len())
                } else if SCAN == 4 {
                    enumerate_single::<NativeByte64, false, false>(&form, 1u64 << labels.len())
                } else if SCAN == 3 {
                    enumerate_tiered::<NativeByte64, false>(&form, 1u64 << labels.len())
                } else {
                    enumerate_tiered::<ScalarByte64, false>(&form, 1u64 << labels.len())
                };
                self.profile.full_update_words += got.linear_updates;
                self.profile.secondary_update_words += got.secondary_updates;
                add_screen_work(&mut self.logical, &got.logical);
                Enumeration {
                    model: got.model,
                    complete: got.complete,
                    points: 0,
                    batches: 0,
                    checksum: got.trace,
                }
            } else if SCAN == 1 {
                enumerate_quiet(&form, 1u64 << labels.len())
            } else if DELTA {
                if SIMD {
                    enumerate_delta::<NativeDeltaBlock>(&form, 1u64 << labels.len())
                } else {
                    enumerate_delta::<ScalarDeltaBlock>(&form, 1u64 << labels.len())
                }
            } else {
                enumerate_syndromes::<SIMD>(&form, 1u64 << labels.len())
            };
            assert!(got.complete);
            self.logical.enumeration_leaves += 1;
            self.logical.enumeration_points += got.points;
            self.logical.enumeration_batches += got.batches;
            self.event(0x850);
            self.event(relevant);
            self.event(got.points);
            self.event(got.checksum);
            let outcome = if let Some(local) = got.model {
                let mut model = values;
                for (p, &j) in labels.iter().enumerate() {
                    if local & (1u64 << p) != 0 {
                        model |= 1u64 << j;
                    }
                }
                self.event(0x851);
                self.event(model);
                Outcome::Sat(model)
            } else {
                self.event(0x852);
                Outcome::Unsat
            };
            if SCAN >= 2 {
                self.profile.partial_ns += tick.elapsed().as_nanos();
            } else {
                self.profile.enumeration_ns += tick.elapsed().as_nanos();
            }
            return outcome;
        }
        let variable = (0..self.model.n as usize)
            .filter(|&i| counts[i] > 0)
            .max_by_key(|&i| (counts[i], std::cmp::Reverse(i)))
            .unwrap();
        let bit = 1u64 << variable;
        self.logical.decisions += 1;
        self.event(0x400);
        self.event(variable as u64);
        let left = self
            .model
            .specialize(state.clone(), bit, 0, &mut self.logical);
        match self.visit(left, known | bit, values, depth + 1) {
            Outcome::Unsat => {
                let right = self.model.specialize(state, bit, bit, &mut self.logical);
                self.visit(right, known | bit, values | bit, depth + 1)
            }
            answer => answer,
        }
    }
}

fn solve_leaf_backend<const SIMD: bool, const CUT: usize>(
    system: &System,
    n: u8,
    limit: u64,
) -> Solved {
    solve_leaf_engine::<SIMD, CUT, false>(system, n, limit)
}
fn solve_leaf_delta<const SIMD: bool>(system: &System, n: u8, limit: u64) -> Solved {
    solve_leaf_engine::<SIMD, 16, true>(system, n, limit)
}
fn solve_leaf_engine<const SIMD: bool, const CUT: usize, const DELTA: bool>(
    system: &System,
    n: u8,
    limit: u64,
) -> Solved {
    solve_leaf_scan_engine::<SIMD, CUT, DELTA, 0>(system, n, limit)
}
fn solve_leaf_scan_engine<const SIMD: bool, const CUT: usize, const DELTA: bool, const SCAN: u8>(
    system: &System,
    n: u8,
    limit: u64,
) -> Solved {
    let Some((model, state)) = PackedModel::compile(system, n) else {
        return solve_packed(system, n, limit);
    };
    let mut solver = LeafSolver::<SIMD, CUT, DELTA, SCAN> {
        model,
        limit,
        logical: Logical::default(),
        profile: Profile::default(),
        trace: 0xcbf29ce484222325,
    };
    let outcome = solver.visit(state, 0, 0, 0);
    Solved {
        outcome,
        logical: solver.logical,
        profile: Profile {
            enumeration_ns: solver.profile.enumeration_ns,
            partial_ns: solver.profile.partial_ns,
            full_update_words: solver.profile.full_update_words,
            secondary_update_words: solver.profile.secondary_update_words,
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: solver.trace,
    }
}
fn solve_leaf(system: &System, n: u8, limit: u64, cut: usize, simd: bool) -> Solved {
    match (cut, simd) {
        (8, false) => solve_leaf_backend::<false, 8>(system, n, limit),
        (8, true) => solve_leaf_backend::<true, 8>(system, n, limit),
        (12, false) => solve_leaf_backend::<false, 12>(system, n, limit),
        (12, true) => solve_leaf_backend::<true, 12>(system, n, limit),
        (16, false) => solve_leaf_backend::<false, 16>(system, n, limit),
        (16, true) => solve_leaf_backend::<true, 16>(system, n, limit),
        _ => unreachable!(),
    }
}
