// Accounting control: the retained packed algorithm without diagnostic hashing.
struct UntracedPackedSolver {
    model: PackedModel,
    limit: u64,
    logical: Logical,
    trace: u64,
}
impl UntracedPackedSolver {
    #[inline(always)]
    fn event(&mut self, _value: u64) {}
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

fn solve_packed_untraced(system: &System, n: u8, limit: u64) -> Solved {
    let Some((model, state)) = PackedModel::compile(system, n) else {
        return solve_packed(system, n, limit);
    };
    let mut solver = UntracedPackedSolver {
        model,
        limit,
        logical: Logical::default(),
        trace: 0,
    };
    let outcome = solver.visit(state, 0, 0, 0);
    Solved {
        outcome,
        logical: solver.logical,
        profile: Profile {
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: solver.trace,
    }
}
