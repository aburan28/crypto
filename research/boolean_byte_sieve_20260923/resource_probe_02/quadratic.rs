// Exact residual representation for degree-at-most-two Boolean systems.
// The quadratic coefficients between unassigned variables never change.
const MAX: usize = 36;
struct QuadraticModel {
    n: u8,
    neighbours: Vec<[u64; MAX]>,
}
#[derive(Clone)]
struct QuadraticState {
    active: u64,
    live: u64,
    linear: [u64; MAX],
    constants: u64,
    edges: [u16; MAX],
}
impl QuadraticModel {
    fn compile(system: &System, n: u8) -> Option<(Self, QuadraticState)> {
        if n as usize > MAX || system.len() > MAX {
            return None;
        }
        let mut state = QuadraticState {
            active: (1u64 << n) - 1,
            live: (1u64 << system.len()) - 1,
            linear: [0; MAX],
            constants: 0,
            edges: [0; MAX],
        };
        let mut neighbours = vec![[0u64; MAX]; system.len()];
        for (e, poly) in system.iter().enumerate() {
            if poly
                .windows(2)
                .any(|w| algebra::mono_order(w[0], w[1]) != std::cmp::Ordering::Less)
            {
                return None;
            }
            for &m in poly {
                if m >= 1u64 << n || m.count_ones() > 2 {
                    return None;
                }
                match m.count_ones() {
                    0 => state.constants ^= 1u64 << e,
                    1 => state.linear[e] ^= m,
                    2 => {
                        let a = m.trailing_zeros() as usize;
                        let b = (m & (m - 1)).trailing_zeros() as usize;
                        neighbours[e][a] |= 1u64 << b;
                        neighbours[e][b] |= 1u64 << a;
                        state.edges[e] += 1;
                    }
                    _ => unreachable!(),
                }
            }
        }
        Some((Self { n, neighbours }, state))
    }
    fn terms(&self, state: &QuadraticState, e: usize, mut visit: impl FnMut(u64)) {
        // Numeric order within degree 2 is colex order: larger variable first
        // as the outer index, then smaller variables in increasing order.
        let mut higher = state.active;
        while higher != 0 {
            let bit = higher & higher.wrapping_neg();
            let j = bit.trailing_zeros() as usize;
            let mut lower = self.neighbours[e][j] & state.active & (bit - 1);
            while lower != 0 {
                let i = lower & lower.wrapping_neg();
                visit(bit | i);
                lower &= lower - 1;
            }
            higher &= higher - 1;
        }
        let mut linear = state.linear[e];
        while linear != 0 {
            let bit = linear & linear.wrapping_neg();
            visit(bit);
            linear &= linear - 1;
        }
        if state.constants & (1u64 << e) != 0 {
            visit(0);
        }
    }
    fn same_equation(&self, state: &QuadraticState, a: usize, b: usize) -> bool {
        if state.linear[a] != state.linear[b]
            || state.edges[a] != state.edges[b]
            || ((state.constants >> a) ^ (state.constants >> b)) & 1 != 0
        {
            return false;
        }
        if state.edges[a] == 0 {
            return true;
        }
        let mut active = state.active;
        while active != 0 {
            let bit = active & active.wrapping_neg();
            let j = bit.trailing_zeros() as usize;
            if (self.neighbours[a][j] ^ self.neighbours[b][j]) & state.active & (bit - 1) != 0 {
                return false;
            }
            active &= active - 1;
        }
        true
    }
    fn specialize(
        &self,
        mut state: QuadraticState,
        mask: u64,
        values: u64,
        logical: &mut Logical,
    ) -> QuadraticState {
        let mut equations = state.live;
        while equations != 0 {
            let e = equations.trailing_zeros() as usize;
            logical.specialized_terms += u64::from(state.edges[e])
                + u64::from(state.linear[e].count_ones())
                + ((state.constants >> e) & 1);
            equations &= equations - 1;
        }
        let mut assigned = mask & state.active;
        while assigned != 0 {
            let bit = assigned & assigned.wrapping_neg();
            let j = bit.trailing_zeros() as usize;
            let remaining = state.active & !bit;
            let mut equations = state.live;
            while equations != 0 {
                let e = equations.trailing_zeros() as usize;
                let neighbours = self.neighbours[e][j] & remaining;
                state.edges[e] -= neighbours.count_ones() as u16;
                if values & bit != 0 {
                    if state.linear[e] & bit != 0 {
                        state.constants ^= 1u64 << e;
                    }
                    state.linear[e] ^= neighbours;
                }
                state.linear[e] &= !bit;
                equations &= equations - 1;
            }
            state.active = remaining;
            assigned &= assigned - 1;
        }
        // Remove zeros and exact duplicate residuals in original equation order.
        // No hash-only acceptance and no representative-order change.
        let mut kept = 0u64;
        let mut equations = state.live;
        while equations != 0 {
            let bit = equations & equations.wrapping_neg();
            let e = bit.trailing_zeros() as usize;
            equations &= equations - 1;
            if state.edges[e] == 0 && state.linear[e] == 0 && state.constants & bit == 0 {
                continue;
            }
            let mut prior = kept;
            let mut duplicate = false;
            while prior != 0 {
                let p = prior.trailing_zeros() as usize;
                if self.same_equation(&state, e, p) {
                    duplicate = true;
                    break;
                }
                prior &= prior - 1;
            }
            if !duplicate {
                kept |= bit;
            }
        }
        state.live = kept;
        state
    }
    #[cfg(test)]
    fn materialize(&self, state: &QuadraticState) -> System {
        let mut out = Vec::new();
        let mut live = state.live;
        while live != 0 {
            let e = live.trailing_zeros() as usize;
            let mut row = Vec::new();
            self.terms(state, e, |m| row.push(m));
            out.push(row);
            live &= live - 1;
        }
        out
    }
}
struct QuadraticSolver {
    model: QuadraticModel,
    limit: u64,
    logical: Logical,
    trace: u64,
}
impl QuadraticSolver {
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn visit(
        &mut self,
        mut state: QuadraticState,
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
            self.event(
                u64::from(state.edges[e])
                    + u64::from(state.linear[e].count_ones())
                    + ((state.constants >> e) & 1),
            );
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
            let mut affine = Vec::new();
            let mut all_affine = true;
            while live != 0 {
                let e = live.trailing_zeros() as usize;
                let constant = (state.constants >> e) & 1;
                if state.edges[e] == 0 {
                    if state.linear[e] == 0 && constant != 0 {
                        self.event(0x500);
                        return Outcome::Unsat;
                    }
                    affine.push(state.linear[e] | (constant << self.model.n));
                } else {
                    all_affine = false;
                }
                live &= live - 1;
            }
            let reduced = algebra::affine_rref(self.model.n, affine.into_iter());
            let constant = 1u64 << self.model.n;
            let (mut mask, mut ones) = (0, 0);
            for &row in &reduced {
                if row == constant {
                    self.event(0x502);
                    return Outcome::Unsat;
                }
                let variable = row & (constant - 1);
                if variable.count_ones() == 1 {
                    mask |= variable;
                    if row & constant != 0 {
                        ones |= variable;
                    }
                }
            }
            if all_affine {
                for row in reduced {
                    let variable = row & (constant - 1);
                    if variable != 0 && row & constant != 0 {
                        values |= variable & variable.wrapping_neg();
                    }
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
            let mut variables = state.active;
            while variables != 0 {
                let i = variables.trailing_zeros() as usize;
                counts[i] += (self.model.neighbours[e][i] & state.active).count_ones() as usize
                    + usize::from(state.linear[e] & (1u64 << i) != 0);
                variables &= variables - 1;
            }
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
fn solve_quadratic(system: &System, n: u8, limit: u64) -> Solved {
    let Some((model, state)) = QuadraticModel::compile(system, n) else {
        return solve(system, n, "merge_search", limit);
    };
    let mut solver = QuadraticSolver {
        model,
        limit,
        logical: Logical::default(),
        trace: 0xcbf29ce484222325,
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
