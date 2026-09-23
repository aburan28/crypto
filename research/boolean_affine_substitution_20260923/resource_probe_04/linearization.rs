// Sound Boolean products are admitted only when canonical degree stays <= 2.
// This is a bounded inference rule, not a claim of complete ideal closure.
struct LinearizedList {
    active: std::cell::Cell<u64>,
    extra: std::cell::RefCell<ExtraReduction>,
}
impl LinearizedList {
    fn new(n: u8) -> Self {
        Self {
            active: std::cell::Cell::new((1u64 << n) - 1),
            extra: Default::default(),
        }
    }
}
impl BasisBackend for LinearizedList {
    type State = System;
    fn build(&self, s: &System) -> System {
        TailList.build(s)
    }
    fn row_count(&self, s: &System) -> usize {
        s.len()
    }
    fn term_count(&self, s: &System, r: usize) -> usize {
        s[r].len()
    }
    fn terms(&self, s: &System, r: usize, f: impl FnMut(u64)) {
        TailList.terms(s, r, f)
    }
    fn column_count(&self, s: &System) -> usize {
        TailList.column_count(s)
    }
    fn specialize(&self, s: System, m: u64, v: u64, a: u64) -> System {
        TailList.specialize(s, m, v, a)
    }
    fn reduce(&self, s: System) -> System {
        *self.extra.borrow_mut() = ExtraReduction::default();
        let mut state = TailList.reduce(s);
        if self.active.get().count_ones() > 10 {
            return state;
        }
        loop {
            if state.is_empty() || state.iter().any(|p| p[0].count_ones() <= 1) {
                return state;
            }
            let mut products = Vec::new();
            for row in &state {
                let mut variables = self.active.get();
                while variables != 0 {
                    let bit = variables & variables.wrapping_neg();
                    variables &= variables - 1;
                    let mut product: Vec<_> = row.iter().map(|&m| m | bit).collect();
                    algebra::canonical(&mut product);
                    if !product.is_empty()
                        && product.iter().all(|m| m.count_ones() <= 2)
                        && !state.contains(&product)
                        && !products.contains(&product)
                    {
                        products.push(product);
                    }
                }
            }
            if products.is_empty() {
                return state;
            }
            let old_rank = state.len();
            self.extra.borrow_mut().derived += products.len() as u64;
            state.extend(products);
            {
                let mut extra = self.extra.borrow_mut();
                extra.calls += 1;
                extra.rows += state.len() as u64;
                extra.columns += TailList.column_count(&state) as u64;
            }
            state = TailList.reduce(state);
            assert!(state.len() >= old_rank && state.len() <= 56);
            if state.len() == old_rank {
                return state;
            }
        }
    }
}
impl AffineBackend for LinearizedList {
    fn set_active(&self, active: u64) {
        self.active.set(active);
    }
    fn take_extra_reduction(&self) -> ExtraReduction {
        std::mem::take(&mut *self.extra.borrow_mut())
    }
    fn transform_pair(&self, a: System, b: System, m: &RecoveryMap, n: u8) -> (System, System) {
        TailList.transform_pair(a, b, m, n)
    }
    fn clean_source(&self, s: System) -> System {
        TailList.clean_source(s)
    }
}

fn quadratic_common_variables(row: u64, layout: &CompactLayout) -> u64 {
    let mut common = (1u64 << layout.labels.len()) - 1;
    let mut terms = row & ((1u64 << layout.quadratic) - 1);
    while terms != 0 {
        let (a, b) = layout.factors[terms.trailing_zeros() as usize];
        common &= (1u64 << a) | (1u64 << b);
        terms &= terms - 1;
    }
    common
}
fn multiply_word_by_variable(row: u64, layout: &CompactLayout, variable: usize) -> u64 {
    let k = layout.labels.len();
    let q = layout.quadratic;
    let mut out = row & ((1u64 << q) - 1);
    let mut linear = (row >> q) & ((1u64 << k) - 1);
    while linear != 0 {
        let j = linear.trailing_zeros() as usize;
        out ^= 1u64 << local_coordinate((1u64 << j) | (1u64 << variable), k);
        linear &= linear - 1;
    }
    if row & (1u64 << (q + k)) != 0 {
        out ^= 1u64 << (q + variable);
    }
    out
}
struct LinearizedCompact {
    inner: CompactBackend,
    active: std::cell::Cell<u64>,
    extra: std::cell::RefCell<ExtraReduction>,
}
impl LinearizedCompact {
    fn new(n: u8) -> Self {
        Self {
            inner: CompactBackend::new(n),
            active: std::cell::Cell::new((1u64 << n) - 1),
            extra: Default::default(),
        }
    }
    fn new_fast(n: u8) -> Self {
        let mut result = Self::new(n);
        result.inner = CompactBackend::new_fast(n);
        result
    }
}
impl BasisBackend for LinearizedCompact {
    type State = CompactState;
    fn build(&self, s: &System) -> Self::State {
        self.inner.build(s)
    }
    fn row_count(&self, s: &Self::State) -> usize {
        self.inner.row_count(s)
    }
    fn term_count(&self, s: &Self::State, r: usize) -> usize {
        self.inner.term_count(s, r)
    }
    fn terms(&self, s: &Self::State, r: usize, f: impl FnMut(u64)) {
        self.inner.terms(s, r, f)
    }
    fn column_count(&self, s: &Self::State) -> usize {
        self.inner.column_count(s)
    }
    fn specialize(&self, s: Self::State, m: u64, v: u64, a: u64) -> Self::State {
        self.inner.specialize(s, m, v, a)
    }
    fn reduce(&self, s: Self::State) -> Self::State {
        *self.extra.borrow_mut() = ExtraReduction::default();
        assert_eq!(s.layout.active, self.active.get());
        let mut state = self.inner.reduce(s);
        if state.layout.labels.len() > 10 {
            return state;
        }
        loop {
            let CompactRows::Word(rows) = &state.rows else {
                unreachable!()
            };
            if rows.is_empty()
                || rows.last().unwrap().trailing_zeros() as usize >= state.layout.quadratic
            {
                return state;
            }
            let mut products = Vec::new();
            for &row in rows {
                let mut common = quadratic_common_variables(row, &state.layout);
                while common != 0 {
                    let j = common.trailing_zeros() as usize;
                    common &= common - 1;
                    let product = multiply_word_by_variable(row, &state.layout, j);
                    if product != 0 && !rows.contains(&product) && !products.contains(&product) {
                        products.push(product);
                    }
                }
            }
            if products.is_empty() {
                return state;
            }
            let old_rank = rows.len();
            let layout = state.layout;
            let CompactRows::Word(mut expanded) = state.rows else {
                unreachable!()
            };
            self.extra.borrow_mut().derived += products.len() as u64;
            expanded.extend(products);
            {
                let mut extra = self.extra.borrow_mut();
                extra.calls += 1;
                extra.rows += expanded.len() as u64;
                extra.columns += u64::from(expanded.iter().fold(0u64, |a, b| a | b).count_ones());
            }
            state = self.inner.reduce(CompactState {
                layout,
                rows: CompactRows::Word(expanded),
            });
            let rank = self.inner.row_count(&state);
            assert!(rank >= old_rank && rank <= 56);
            if rank == old_rank {
                return state;
            }
        }
    }
}
impl AffineBackend for LinearizedCompact {
    fn set_active(&self, active: u64) {
        self.active.set(active);
    }
    fn take_extra_reduction(&self) -> ExtraReduction {
        std::mem::take(&mut *self.extra.borrow_mut())
    }
    fn transform_pair(
        &self,
        a: Self::State,
        b: Self::State,
        m: &RecoveryMap,
        n: u8,
    ) -> (Self::State, Self::State) {
        self.inner.transform_pair(a, b, m, n)
    }
    fn clean_source(&self, s: Self::State) -> Self::State {
        self.inner.clean_source(s)
    }
}
fn solve_affine_linearized(system: &System, n: u8, limit: u64, compact: bool) -> Solved {
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
    if compact {
        solve_affine_backend(system, n, limit, LinearizedCompact::new(n))
    } else {
        solve_affine_backend(system, n, limit, LinearizedList::new(n))
    }
}

fn solve_affine_linearized_fast(system: &System, n: u8, limit: u64) -> Solved {
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
    solve_affine_backend(system, n, limit, LinearizedCompact::new_fast(n))
}
