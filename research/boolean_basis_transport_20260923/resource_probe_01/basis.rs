// Canonical coefficient-row-space transport. This is linear row reduction over
// quadratic Boolean polynomials, not a Macaulay multiplication closure.
const BASIS_WORDS: usize = 11;
const BASIS_COLUMNS: usize = 667; // C(36, 2) + 36 + 1
type CoefficientRow = [u64; BASIS_WORDS];

trait BasisBackend {
    type State: Clone;
    fn build(&self, system: &System) -> Self::State;
    fn row_count(&self, state: &Self::State) -> usize;
    fn term_count(&self, state: &Self::State, row: usize) -> usize;
    fn terms(&self, state: &Self::State, row: usize, visit: impl FnMut(u64));
    fn column_count(&self, state: &Self::State) -> usize;
    fn reduce(&self, state: Self::State) -> Self::State;
    fn specialize(&self, state: Self::State, mask: u64, values: u64, active: u64) -> Self::State;
}

fn xor_polynomials(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match algebra::mono_order(a[i], b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}
struct ListBasis;
impl BasisBackend for ListBasis {
    type State = System;
    fn build(&self, system: &System) -> System {
        system.clone()
    }
    fn row_count(&self, state: &System) -> usize {
        state.len()
    }
    fn term_count(&self, state: &System, row: usize) -> usize {
        state[row].len()
    }
    fn terms(&self, state: &System, row: usize, mut visit: impl FnMut(u64)) {
        for &term in &state[row] {
            visit(term);
        }
    }
    fn column_count(&self, state: &System) -> usize {
        state
            .iter()
            .flatten()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    }
    fn reduce(&self, state: System) -> System {
        let mut basis: System = Vec::new();
        for mut row in state {
            while !row.is_empty() {
                if let Some(pivot) = basis.iter().find(|p| p[0] == row[0]) {
                    row = xor_polynomials(&row, pivot);
                } else {
                    basis.push(row);
                    break;
                }
            }
        }
        basis.sort_unstable_by(|a, b| algebra::mono_order(a[0], b[0]));
        for p in (0..basis.len()).rev() {
            let pivot = basis[p][0];
            for q in 0..p {
                if basis[q].contains(&pivot) {
                    basis[q] = xor_polynomials(&basis[q], &basis[p]);
                }
            }
        }
        basis
    }
    fn specialize(&self, state: System, mask: u64, values: u64, _active: u64) -> System {
        let zero = mask & !values;
        state
            .into_iter()
            .map(|row| {
                let mut out: Vec<_> = row
                    .into_iter()
                    .filter(|m| m & zero == 0)
                    .map(|m| m & !mask)
                    .collect();
                algebra::canonical(&mut out);
                out
            })
            .collect()
    }
}

struct WideBasis {
    n: u8,
    quadratic_columns: usize,
    words: usize,
    monomials: Vec<u64>,
    kill: Vec<CoefficientRow>,
}
fn row_bit(row: &CoefficientRow, column: usize) -> bool {
    row[column / 64] & (1u64 << (column % 64)) != 0
}
fn toggle_bit(row: &mut CoefficientRow, column: usize) {
    row[column / 64] ^= 1u64 << (column % 64);
}
fn read_span(row: &CoefficientRow, start: usize, length: usize) -> u64 {
    if length == 0 {
        return 0;
    }
    let (word, offset) = (start / 64, start % 64);
    let mut value = row[word] >> offset;
    if offset + length > 64 {
        value |= row[word + 1] << (64 - offset);
    }
    value & ((1u64 << length) - 1)
}
fn xor_span(row: &mut CoefficientRow, start: usize, value: u64, length: usize) {
    let (word, offset) = (start / 64, start % 64);
    row[word] ^= value << offset;
    if offset + length > 64 {
        row[word + 1] ^= value >> (64 - offset);
    }
}
fn first_column(row: &CoefficientRow, words: usize) -> Option<usize> {
    (0..words)
        .find(|&w| row[w] != 0)
        .map(|w| w * 64 + row[w].trailing_zeros() as usize)
}
impl WideBasis {
    fn new(n: u8) -> Self {
        let q = usize::from(n) * usize::from(n.saturating_sub(1)) / 2;
        let mut monomials = Vec::with_capacity(q + usize::from(n) + 1);
        let mut kill = vec![[0u64; BASIS_WORDS]; n as usize];
        for b in 1..n as usize {
            for a in 0..b {
                let column = monomials.len();
                monomials.push((1u64 << a) | (1u64 << b));
                toggle_bit(&mut kill[a], column);
                toggle_bit(&mut kill[b], column);
            }
        }
        for j in 0..n as usize {
            monomials.push(1u64 << j);
            toggle_bit(&mut kill[j], q + j);
        }
        monomials.push(0);
        Self {
            n,
            quadratic_columns: q,
            words: monomials.len().div_ceil(64),
            monomials,
            kill,
        }
    }
    fn coordinate(&self, monomial: u64) -> usize {
        match monomial.count_ones() {
            0 => self.quadratic_columns + usize::from(self.n),
            1 => self.quadratic_columns + monomial.trailing_zeros() as usize,
            2 => {
                let a = monomial.trailing_zeros() as usize;
                let b = (monomial & (monomial - 1)).trailing_zeros() as usize;
                b * (b - 1) / 2 + a
            }
            _ => unreachable!(),
        }
    }
    #[cfg(test)]
    fn materialize(&self, state: &[CoefficientRow]) -> System {
        state
            .iter()
            .map(|row| {
                let mut out = Vec::new();
                for w in 0..self.words {
                    let mut bits = row[w];
                    while bits != 0 {
                        out.push(self.monomials[w * 64 + bits.trailing_zeros() as usize]);
                        bits &= bits - 1;
                    }
                }
                out
            })
            .collect()
    }
}
impl BasisBackend for WideBasis {
    type State = Vec<CoefficientRow>;
    fn build(&self, system: &System) -> Self::State {
        system
            .iter()
            .map(|poly| {
                let mut row = [0u64; BASIS_WORDS];
                for &monomial in poly {
                    toggle_bit(&mut row, self.coordinate(monomial));
                }
                row
            })
            .collect()
    }
    fn row_count(&self, state: &Self::State) -> usize {
        state.len()
    }
    fn term_count(&self, state: &Self::State, row: usize) -> usize {
        state[row][..self.words]
            .iter()
            .map(|w| w.count_ones() as usize)
            .sum()
    }
    fn terms(&self, state: &Self::State, row: usize, mut visit: impl FnMut(u64)) {
        for w in 0..self.words {
            let mut bits = state[row][w];
            while bits != 0 {
                visit(self.monomials[w * 64 + bits.trailing_zeros() as usize]);
                bits &= bits - 1;
            }
        }
    }
    fn column_count(&self, state: &Self::State) -> usize {
        let mut union = [0u64; BASIS_WORDS];
        for row in state {
            for w in 0..self.words {
                union[w] |= row[w];
            }
        }
        union[..self.words]
            .iter()
            .map(|v| v.count_ones() as usize)
            .sum()
    }
    fn reduce(&self, state: Self::State) -> Self::State {
        let mut pivots = [-1i8; BASIS_COLUMNS];
        let mut basis: Vec<CoefficientRow> = Vec::with_capacity(state.len());
        for mut row in state {
            while let Some(p) = first_column(&row, self.words) {
                if pivots[p] < 0 {
                    pivots[p] = basis.len() as i8;
                    basis.push(row);
                    break;
                }
                let pivot = &basis[pivots[p] as usize];
                for w in p / 64..self.words {
                    row[w] ^= pivot[w];
                }
            }
        }
        basis.sort_unstable_by_key(|r| first_column(r, self.words).unwrap());
        for p in (0..basis.len()).rev() {
            let column = first_column(&basis[p], self.words).unwrap();
            let pivot = basis[p];
            for q in 0..p {
                if row_bit(&basis[q], column) {
                    for w in column / 64..self.words {
                        basis[q][w] ^= pivot[w];
                    }
                }
            }
        }
        basis
    }
    fn specialize(
        &self,
        mut state: Self::State,
        mask: u64,
        values: u64,
        mut active: u64,
    ) -> Self::State {
        let mut assigned = mask & active;
        while assigned != 0 {
            let bit = assigned & assigned.wrapping_neg();
            let j = bit.trailing_zeros() as usize;
            for row in &mut state {
                let mut delta = 0;
                if values & bit != 0 {
                    delta = read_span(row, j * j.saturating_sub(1) / 2, j);
                    let mut higher = active & !((bit << 1) - 1);
                    while higher != 0 {
                        let k = higher.trailing_zeros() as usize;
                        if row_bit(row, k * (k - 1) / 2 + j) {
                            delta ^= 1u64 << k;
                        }
                        higher &= higher - 1;
                    }
                    if row_bit(row, self.quadratic_columns + j) {
                        delta ^= 1u64 << self.n;
                    }
                }
                for w in 0..self.words {
                    row[w] &= !self.kill[j][w];
                }
                xor_span(row, self.quadratic_columns, delta, usize::from(self.n) + 1);
            }
            active &= !bit;
            assigned &= assigned - 1;
        }
        state
    }
}

struct BasisSolver<B: BasisBackend> {
    backend: B,
    n: u8,
    limit: u64,
    logical: Logical,
    profile: Profile,
    trace: u64,
}
impl<B: BasisBackend> BasisSolver<B> {
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn specialize(&mut self, state: B::State, mask: u64, values: u64, active: u64) -> B::State {
        for r in 0..self.backend.row_count(&state) {
            self.logical.specialized_terms += self.backend.term_count(&state, r) as u64;
        }
        self.backend.specialize(state, mask, values, active)
    }
    fn visit(
        &mut self,
        mut state: B::State,
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
        self.event(0x600);
        self.event(known);
        self.event(values);
        let all_variables = (1u64 << self.n) - 1;
        let constant = 1u64 << self.n;
        loop {
            let tick = Instant::now();
            self.logical.kernel_calls += 1;
            self.profile.kernel_calls_by_active[(all_variables & !known).count_ones() as usize] +=
                1;
            self.logical.source_rows += self.backend.row_count(&state) as u64;
            self.logical.source_columns += self.backend.column_count(&state) as u64;
            state = self.backend.reduce(state);
            self.profile.kernel_ns += tick.elapsed().as_nanos();
            self.event(0x610);
            self.event(self.backend.row_count(&state) as u64);
            let mut trace = self.trace;
            let (mut mask, mut ones, mut model_ones) = (0, 0, 0);
            let mut all_affine = true;
            let mut contradiction = false;
            let mut counts = [0usize; MAX];
            for r in 0..self.backend.row_count(&state) {
                trace =
                    (trace ^ self.backend.term_count(&state, r) as u64).wrapping_mul(0x100000001b3);
                let (mut has_quadratic, mut affine) = (false, 0u64);
                self.backend.terms(&state, r, |m| {
                    trace = (trace ^ m).wrapping_mul(0x100000001b3);
                    if m.count_ones() > 1 {
                        has_quadratic = true;
                    } else {
                        affine |= if m == 0 { constant } else { m };
                    }
                    let mut bits = m;
                    while bits != 0 {
                        counts[bits.trailing_zeros() as usize] += 1;
                        bits &= bits - 1;
                    }
                });
                all_affine &= !has_quadratic;
                if !has_quadratic {
                    contradiction |= affine == constant;
                    let variables = affine & all_variables;
                    if variables.count_ones() == 1 {
                        mask |= variables;
                        if affine & constant != 0 {
                            ones |= variables;
                        }
                    }
                    if variables != 0 && affine & constant != 0 {
                        model_ones |= variables & variables.wrapping_neg();
                    }
                }
            }
            self.trace = trace;
            if contradiction {
                self.event(0x650);
                return Outcome::Unsat;
            }
            if all_affine {
                values |= model_ones;
                self.event(0x651);
                self.event(values);
                return Outcome::Sat(values);
            }
            if mask != 0 {
                self.logical.forced += u64::from(mask.count_ones());
                self.event(0x620);
                self.event(mask);
                self.event(ones);
                state = self.specialize(state, mask, ones, all_variables & !known);
                known |= mask;
                values |= ones;
                continue;
            }
            let variable = (0..self.n as usize)
                .filter(|&i| counts[i] > 0)
                .max_by_key(|&i| (counts[i], std::cmp::Reverse(i)))
                .unwrap();
            let bit = 1u64 << variable;
            self.logical.decisions += 1;
            self.event(0x640);
            self.event(variable as u64);
            let left = self.specialize(state.clone(), bit, 0, all_variables & !known);
            return match self.visit(left, known | bit, values, depth + 1) {
                Outcome::Unsat => {
                    let right = self.specialize(state, bit, bit, all_variables & !known);
                    self.visit(right, known | bit, values | bit, depth + 1)
                }
                outcome => outcome,
            };
        }
    }
}
fn solve_basis_backend<B: BasisBackend>(system: &System, n: u8, limit: u64, backend: B) -> Solved {
    let state = backend.build(system);
    let mut solver = BasisSolver {
        backend,
        n,
        limit,
        logical: Logical::default(),
        profile: Profile {
            kernel_calls_by_active: vec![0; usize::from(n) + 1],
            ..Profile::default()
        },
        trace: 0xcbf29ce484222325,
    };
    let outcome = solver.visit(state, 0, 0, 0);
    Solved {
        outcome,
        logical: solver.logical,
        profile: solver.profile,
        trace: solver.trace,
    }
}
fn solve_basis(system: &System, n: u8, limit: u64, wide: bool) -> Solved {
    if n as usize > MAX
        || system.len() > MAX
        || system.iter().any(|p| {
            p.iter().any(|&m| m >= 1u64 << n || m.count_ones() > 2)
                || p.windows(2)
                    .any(|w| algebra::mono_order(w[0], w[1]) != std::cmp::Ordering::Less)
        })
    {
        return solve_packed(system, n, limit);
    }
    if wide {
        solve_basis_backend(system, n, limit, WideBasis::new(n))
    } else {
        solve_basis_backend(system, n, limit, ListBasis)
    }
}
