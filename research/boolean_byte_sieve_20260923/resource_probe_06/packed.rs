// One coefficient word per residual equation. Original quadratic terms occupy
// local low bits; shared variable coordinates and the constant follow them.
struct PackedEquation {
    quadratic: [u64; 32],
    neighbours: [u64; MAX],
    touches: [u64; MAX],
    q_len: u8,
    q_mask: u64,
}
struct PackedModel {
    n: u8,
    equations: Vec<PackedEquation>,
}
#[derive(Clone)]
struct PackedState {
    active: u64,
    live: u64,
    rows: [u64; MAX],
}
impl PackedModel {
    fn compile(system: &System, n: u8) -> Option<(Self, PackedState)> {
        if n as usize > MAX || system.len() > MAX {
            return None;
        }
        let mut state = PackedState {
            active: (1u64 << n) - 1,
            live: (1u64 << system.len()) - 1,
            rows: [0; MAX],
        };
        let mut equations = Vec::with_capacity(system.len());
        for (e, poly) in system.iter().enumerate() {
            if poly
                .windows(2)
                .any(|w| algebra::mono_order(w[0], w[1]) != std::cmp::Ordering::Less)
                || poly.iter().any(|&m| m >= 1u64 << n || m.count_ones() > 2)
            {
                return None;
            }
            let q_len = poly.iter().take_while(|m| m.count_ones() == 2).count();
            if q_len > 32 || q_len + usize::from(n) >= 64 {
                return None;
            }
            let mut equation = PackedEquation {
                quadratic: [0; 32],
                neighbours: [0; MAX],
                touches: [0; MAX],
                q_len: q_len as u8,
                q_mask: (1u64 << q_len) - 1,
            };
            for (q, &monomial) in poly.iter().take(q_len).enumerate() {
                equation.quadratic[q] = monomial;
                let a = monomial.trailing_zeros() as usize;
                let b = (monomial & (monomial - 1)).trailing_zeros() as usize;
                equation.neighbours[a] |= 1u64 << b;
                equation.neighbours[b] |= 1u64 << a;
                equation.touches[a] |= 1u64 << q;
                equation.touches[b] |= 1u64 << q;
                state.rows[e] |= 1u64 << q;
            }
            for j in 0..n as usize {
                equation.touches[j] |= 1u64 << (q_len + j);
            }
            for &monomial in &poly[q_len..] {
                state.rows[e] |= if monomial == 0 {
                    1u64 << (q_len + usize::from(n))
                } else {
                    monomial << q_len
                };
            }
            equations.push(equation);
        }
        Some((Self { n, equations }, state))
    }
    fn terms(&self, state: &PackedState, e: usize, mut visit: impl FnMut(u64)) {
        let equation = &self.equations[e];
        let row = state.rows[e];
        let mut quad = row & equation.q_mask;
        while quad != 0 {
            visit(equation.quadratic[quad.trailing_zeros() as usize]);
            quad &= quad - 1;
        }
        let affine = row >> equation.q_len;
        let mut linear = affine & ((1u64 << self.n) - 1);
        while linear != 0 {
            let bit = linear & linear.wrapping_neg();
            visit(bit);
            linear &= linear - 1;
        }
        if affine & (1u64 << self.n) != 0 {
            visit(0);
        }
    }
    fn same_equation(&self, state: &PackedState, a: usize, b: usize) -> bool {
        let (ma, mb) = (&self.equations[a], &self.equations[b]);
        let (ra, rb) = (state.rows[a], state.rows[b]);
        if ra >> ma.q_len != rb >> mb.q_len {
            return false;
        }
        let (mut qa, mut qb) = (ra & ma.q_mask, rb & mb.q_mask);
        if qa.count_ones() != qb.count_ones() {
            return false;
        }
        while qa != 0 {
            if ma.quadratic[qa.trailing_zeros() as usize]
                != mb.quadratic[qb.trailing_zeros() as usize]
            {
                return false;
            }
            qa &= qa - 1;
            qb &= qb - 1;
        }
        true
    }
    fn specialize(
        &self,
        mut state: PackedState,
        mask: u64,
        values: u64,
        logical: &mut Logical,
    ) -> PackedState {
        let mut live = state.live;
        while live != 0 {
            let e = live.trailing_zeros() as usize;
            logical.specialized_terms += u64::from(state.rows[e].count_ones());
            live &= live - 1;
        }
        let mut assigned = mask & state.active;
        while assigned != 0 {
            let bit = assigned & assigned.wrapping_neg();
            let j = bit.trailing_zeros() as usize;
            let remaining = state.active & !bit;
            let mut live = state.live;
            while live != 0 {
                let e = live.trailing_zeros() as usize;
                let equation = &self.equations[e];
                let row = &mut state.rows[e];
                if values & bit != 0 {
                    if *row & (bit << equation.q_len) != 0 {
                        *row ^= 1u64 << (self.n + equation.q_len);
                    }
                    *row ^= (equation.neighbours[j] & remaining) << equation.q_len;
                }
                *row &= !equation.touches[j];
                live &= live - 1;
            }
            state.active = remaining;
            assigned &= assigned - 1;
        }
        let mut buckets = [0u64; 64];
        let mut kept = 0u64;
        let mut live = state.live;
        while live != 0 {
            let bit = live & live.wrapping_neg();
            let e = bit.trailing_zeros() as usize;
            live &= live - 1;
            let row = state.rows[e];
            if row == 0 {
                continue;
            }
            let equation = &self.equations[e];
            let affine = row >> equation.q_len;
            let degree_two_terms = u64::from((row & equation.q_mask).count_ones());
            let key = affine.wrapping_mul(0x9e3779b97f4a7c15)
                ^ degree_two_terms.wrapping_mul(0xbf58476d1ce4e5b9);
            let bucket = ((key ^ (key >> 32)) & 63) as usize;
            let mut prior = buckets[bucket];
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
                buckets[bucket] |= bit;
            }
        }
        state.live = kept;
        state
    }
    #[cfg(test)]
    fn materialize(&self, state: &PackedState) -> System {
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
struct StackAffine {
    rows: [u64; 37],
    present: u64,
}
impl StackAffine {
    fn new() -> Self {
        Self {
            rows: [0; 37],
            present: 0,
        }
    }
    fn insert(&mut self, mut row: u64) {
        while row != 0 {
            let p = row.trailing_zeros() as usize;
            if self.rows[p] == 0 {
                self.rows[p] = row;
                self.present |= 1u64 << p;
                return;
            }
            row ^= self.rows[p];
        }
    }
    fn reduce(&mut self) {
        let mut higher = self.present;
        while higher != 0 {
            let p = 63 - higher.leading_zeros();
            let bit = 1u64 << p;
            let mut lower = self.present & (bit - 1);
            while lower != 0 {
                let q = lower.trailing_zeros() as usize;
                if self.rows[q] & bit != 0 {
                    self.rows[q] ^= self.rows[p as usize];
                }
                lower &= lower - 1;
            }
            higher &= !bit;
        }
    }
}
struct PackedSolver {
    model: PackedModel,
    limit: u64,
    logical: Logical,
    trace: u64,
}
impl PackedSolver {
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
fn solve_packed(system: &System, n: u8, limit: u64) -> Solved {
    let Some((model, state)) = PackedModel::compile(system, n) else {
        return solve_quadratic(system, n, limit);
    };
    let mut solver = PackedSolver {
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
