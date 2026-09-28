// Exact general affine substitution in the Boolean quotient. Coordinates keep
// their original labels; eliminated variables are represented by recovery maps.
type RecoveryMap = [u64; MAX];

fn identity_affine(active: u64, n: u8) -> RecoveryMap {
    let mut map = [0u64; MAX];
    for j in 0..n as usize {
        if active & (1u64 << j) != 0 {
            map[j] = 1u64 << j;
        }
    }
    map
}
fn affine_terms(expression: u64, n: u8) -> Vec<u64> {
    let mut terms = Vec::new();
    let mut variables = expression & ((1u64 << n) - 1);
    while variables != 0 {
        terms.push(variables & variables.wrapping_neg());
        variables &= variables - 1;
    }
    if expression & (1u64 << n) != 0 {
        terms.push(0);
    }
    terms
}
fn substitute_polynomial(poly: &[u64], map: &RecoveryMap, n: u8) -> Vec<u64> {
    let mut out = Vec::new();
    for &monomial in poly {
        match monomial.count_ones() {
            0 => out.push(0),
            1 => out.extend(affine_terms(map[monomial.trailing_zeros() as usize], n)),
            2 => {
                let a = monomial.trailing_zeros() as usize;
                let b = (monomial & (monomial - 1)).trailing_zeros() as usize;
                for x in affine_terms(map[a], n) {
                    for y in affine_terms(map[b], n) {
                        out.push(x | y);
                    }
                }
            }
            _ => unreachable!(),
        }
    }
    algebra::canonical(&mut out);
    out
}
fn compose_recovery(recovery: &mut RecoveryMap, map: &RecoveryMap, eliminated: u64, n: u8) {
    // Every eliminated coordinate is expressed solely in surviving coordinates.
    let mut pivots = eliminated;
    while pivots != 0 {
        let bit = pivots & pivots.wrapping_neg();
        let j = bit.trailing_zeros() as usize;
        assert_eq!(map[j] & eliminated, 0);
        let delta = bit ^ map[j];
        for expression in &mut recovery[..n as usize] {
            if *expression & bit != 0 {
                *expression ^= delta;
            }
        }
        pivots &= pivots - 1;
    }
}
fn recover_zero_free(recovery: &RecoveryMap, n: u8) -> u64 {
    let mut model = 0;
    for (j, &expression) in recovery[..n as usize].iter().enumerate() {
        if expression & (1u64 << n) != 0 {
            model |= 1u64 << j;
        }
    }
    model
}

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq)]
struct ExtraReduction {
    calls: u64,
    rows: u64,
    columns: u64,
    derived: u64,
}
trait AffineBackend: BasisBackend {
    fn set_active(&self, _active: u64) {}
    fn take_extra_reduction(&self) -> ExtraReduction {
        ExtraReduction::default()
    }
    fn transform_pair(
        &self,
        basis: Self::State,
        source: Self::State,
        map: &RecoveryMap,
        n: u8,
    ) -> (Self::State, Self::State);
    fn clean_source(&self, state: Self::State) -> Self::State;
}
impl AffineBackend for TailList {
    fn transform_pair(
        &self,
        basis: System,
        source: System,
        map: &RecoveryMap,
        n: u8,
    ) -> (System, System) {
        let transform = |s: System| {
            s.into_iter()
                .map(|p| substitute_polynomial(&p, map, n))
                .collect()
        };
        (transform(basis), transform(source))
    }
    fn clean_source(&self, state: System) -> System {
        let mut out = Vec::new();
        for row in state {
            if !row.is_empty() && !out.contains(&row) {
                out.push(row);
            }
        }
        out
    }
}
impl WideBasis {
    fn affine_product(&self, mut a: u64, b: u64) -> CoefficientRow {
        let mut out = [0u64; BASIS_WORDS];
        while a != 0 {
            let i = a.trailing_zeros() as usize;
            let x = if i == self.n as usize { 0 } else { 1u64 << i };
            let mut rhs = b;
            while rhs != 0 {
                let j = rhs.trailing_zeros() as usize;
                let y = if j == self.n as usize { 0 } else { 1u64 << j };
                toggle_bit(&mut out, self.coordinate(x | y));
                rhs &= rhs - 1;
            }
            a &= a - 1;
        }
        out
    }
    fn image(&self, monomial: u64, map: &RecoveryMap) -> CoefficientRow {
        let mut image = [0u64; BASIS_WORDS];
        match monomial.count_ones() {
            0 => toggle_bit(&mut image, self.quadratic_columns + self.n as usize),
            1 => xor_span(
                &mut image,
                self.quadratic_columns,
                map[monomial.trailing_zeros() as usize],
                self.n as usize + 1,
            ),
            2 => {
                let i = monomial.trailing_zeros() as usize;
                let j = (monomial & (monomial - 1)).trailing_zeros() as usize;
                image = self.affine_product(map[i], map[j]);
            }
            _ => unreachable!(),
        }
        image
    }
}
impl AffineBackend for WideBasis {
    fn transform_pair(
        &self,
        basis: Self::State,
        source: Self::State,
        map: &RecoveryMap,
        _n: u8,
    ) -> (Self::State, Self::State) {
        // Build each used monomial image once and share it across both carried
        // representations. Unused global coordinates allocate no image payload.
        let mut used = [0u64; BASIS_WORDS];
        for row in basis.iter().chain(&source) {
            for w in 0..self.words {
                used[w] |= row[w];
            }
        }
        let mut index = [-1i16; BASIS_COLUMNS];
        let mut images = Vec::new();
        for w in 0..self.words {
            let mut bits = used[w];
            while bits != 0 {
                let column = w * 64 + bits.trailing_zeros() as usize;
                index[column] = images.len() as i16;
                images.push(self.image(self.monomials[column], map));
                bits &= bits - 1;
            }
        }
        let transform = |state: Vec<CoefficientRow>| {
            state
                .into_iter()
                .map(|row| {
                    let mut out = [0u64; BASIS_WORDS];
                    for w in 0..self.words {
                        let mut bits = row[w];
                        while bits != 0 {
                            let column = w * 64 + bits.trailing_zeros() as usize;
                            let image = &images[index[column] as usize];
                            for k in 0..self.words {
                                out[k] ^= image[k];
                            }
                            bits &= bits - 1;
                        }
                    }
                    out
                })
                .collect()
        };
        (transform(basis), transform(source))
    }
    fn clean_source(&self, state: Self::State) -> Self::State {
        let mut out = Vec::new();
        for row in state {
            if row[..self.words].iter().any(|&x| x != 0) && !out.contains(&row) {
                out.push(row);
            }
        }
        out
    }
}

struct AffineSolver<B: AffineBackend> {
    backend: B,
    n: u8,
    limit: u64,
    logical: Logical,
    profile: Profile,
    trace: u64,
}
impl<B: AffineBackend> AffineSolver<B> {
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn recovery_event(&mut self, recovery: &RecoveryMap) {
        self.event(0x730);
        for &row in &recovery[..self.n as usize] {
            self.event(row);
        }
    }
    fn charge_terms(&mut self, state: &B::State) {
        for r in 0..self.backend.row_count(state) {
            self.logical.specialized_terms += self.backend.term_count(state, r) as u64;
        }
    }
    fn fixed(
        &mut self,
        basis: B::State,
        source: B::State,
        mut recovery: RecoveryMap,
        active: u64,
        bit: u64,
        one: bool,
    ) -> (B::State, B::State, RecoveryMap) {
        let tick = Instant::now();
        self.charge_terms(&basis);
        self.charge_terms(&source);
        let values = if one { bit } else { 0 };
        let basis = self.backend.specialize(basis, bit, values, active);
        let source = self
            .backend
            .clean_source(self.backend.specialize(source, bit, values, active));
        let delta = bit ^ if one { 1u64 << self.n } else { 0 };
        for row in &mut recovery[..self.n as usize] {
            if *row & bit != 0 {
                *row ^= delta;
            }
        }
        self.profile.substitution_ns += tick.elapsed().as_nanos();
        (basis, source, recovery)
    }
    fn visit(
        &mut self,
        mut basis: B::State,
        mut source: B::State,
        mut recovery: RecoveryMap,
        mut active: u64,
        depth: usize,
    ) -> Outcome {
        if self.logical.nodes == self.limit {
            self.event(0xcaf);
            return Outcome::Unknown("NODE_CAP");
        }
        self.logical.nodes += 1;
        self.logical.max_depth = self.logical.max_depth.max(depth);
        self.event(0x700);
        self.event(active);
        self.recovery_event(&recovery);
        let constant = 1u64 << self.n;
        let variables = constant - 1;
        loop {
            let tick = Instant::now();
            self.backend.set_active(active);
            self.logical.kernel_calls += 1;
            self.profile.kernel_calls_by_active[active.count_ones() as usize] += 1;
            self.logical.source_rows += self.backend.row_count(&basis) as u64;
            self.logical.source_columns += self.backend.column_count(&basis) as u64;
            basis = self.backend.reduce(basis);
            let extra = self.backend.take_extra_reduction();
            self.logical.kernel_calls += extra.calls;
            self.profile.kernel_calls_by_active[active.count_ones() as usize] += extra.calls;
            self.logical.source_rows += extra.rows;
            self.logical.source_columns += extra.columns;
            self.logical.derived_rows += extra.derived;
            self.profile.kernel_ns += tick.elapsed().as_nanos();
            self.event(0x710);
            self.event(self.backend.row_count(&basis) as u64);
            let mut map = identity_affine(active, self.n);
            let (mut eliminated, mut unit_count) = (0u64, 0u64);
            let (mut all_affine, mut contradiction) = (true, false);
            let mut trace = self.trace;
            for r in 0..self.backend.row_count(&basis) {
                trace =
                    (trace ^ self.backend.term_count(&basis, r) as u64).wrapping_mul(0x100000001b3);
                let (mut high, mut low) = (false, 0u64);
                self.backend.terms(&basis, r, |m| {
                    trace = (trace ^ m).wrapping_mul(0x100000001b3);
                    if m.count_ones() > 1 {
                        high = true;
                    } else {
                        low |= if m == 0 { constant } else { m };
                    }
                });
                all_affine &= !high;
                if !high {
                    let vars = low & variables;
                    if vars == 0 {
                        contradiction |= low == constant;
                    } else {
                        let pivot = vars & vars.wrapping_neg();
                        assert_eq!(pivot & !active, 0);
                        eliminated |= pivot;
                        map[pivot.trailing_zeros() as usize] = low ^ pivot;
                        unit_count += u64::from(vars.count_ones() == 1);
                    }
                }
            }
            self.trace = trace;
            if contradiction {
                self.event(0x750);
                return Outcome::Unsat;
            }
            if eliminated != 0 {
                let tick = Instant::now();
                self.logical.affine_eliminated += u64::from(eliminated.count_ones());
                self.logical.forced += unit_count;
                self.event(0x720);
                self.event(eliminated);
                let mut pivots = eliminated;
                while pivots != 0 {
                    let p = pivots.trailing_zeros() as usize;
                    self.event(map[p]);
                    pivots &= pivots - 1;
                }
                compose_recovery(&mut recovery, &map, eliminated, self.n);
                active &= !eliminated;
                self.recovery_event(&recovery);
                if !all_affine {
                    self.charge_terms(&basis);
                    self.charge_terms(&source);
                    (basis, source) = self.backend.transform_pair(basis, source, &map, self.n);
                    source = self.backend.clean_source(source);
                }
                self.profile.substitution_ns += tick.elapsed().as_nanos();
            }
            if all_affine {
                let model = recover_zero_free(&recovery, self.n);
                self.event(0x751);
                self.event(model);
                return Outcome::Sat(model);
            }
            if eliminated != 0 {
                continue;
            }
            let mut counts = [0usize; MAX];
            for r in 0..self.backend.row_count(&source) {
                self.backend.terms(&source, r, |mut m| {
                    while m != 0 {
                        counts[m.trailing_zeros() as usize] += 1;
                        m &= m - 1;
                    }
                });
            }
            let variable = (0..self.n as usize)
                .filter(|&i| active & (1u64 << i) != 0 && counts[i] > 0)
                .max_by_key(|&i| (counts[i], std::cmp::Reverse(i)))
                .unwrap();
            let bit = 1u64 << variable;
            self.logical.decisions += 1;
            self.event(0x740);
            self.event(variable as u64);
            let (left, left_source, left_map) =
                self.fixed(basis.clone(), source.clone(), recovery, active, bit, false);
            return match self.visit(left, left_source, left_map, active & !bit, depth + 1) {
                Outcome::Unsat => {
                    let (right, right_source, right_map) =
                        self.fixed(basis, source, recovery, active, bit, true);
                    self.visit(right, right_source, right_map, active & !bit, depth + 1)
                }
                answer => answer,
            };
        }
    }
}
fn solve_affine_backend<B: AffineBackend>(
    system: &System,
    n: u8,
    limit: u64,
    backend: B,
) -> Solved {
    let basis = backend.build(system);
    let source = backend.clean_source(backend.build(system));
    let active = (1u64 << n) - 1;
    let recovery = identity_affine(active, n);
    let mut solver = AffineSolver {
        backend,
        n,
        limit,
        logical: Logical::default(),
        profile: Profile {
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: 0xcbf29ce484222325,
    };
    let outcome = solver.visit(basis, source, recovery, active, 0);
    Solved {
        outcome,
        logical: solver.logical,
        profile: solver.profile,
        trace: solver.trace,
    }
}
fn solve_affine(system: &System, n: u8, limit: u64, wide: bool) -> Solved {
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
        solve_affine_backend(system, n, limit, WideBasis::tail(n))
    } else {
        solve_affine_backend(system, n, limit, TailList)
    }
}
