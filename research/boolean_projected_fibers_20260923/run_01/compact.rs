// Current free coordinates are compacted, while term visitors decode original
// variable labels so the affine policy's trace remains representation-invariant.
#[derive(Clone)]
struct CompactLayout {
    active: u64,
    labels: Vec<u8>,
    positions: [u8; MAX],
    quadratic: usize,
    words: usize,
    monomials: Vec<u64>,
    factors: Vec<(u8, u8)>,
}
impl CompactLayout {
    fn new(active: u64) -> Self {
        let labels: Vec<_> = (0..MAX)
            .filter(|&i| active & (1u64 << i) != 0)
            .map(|i| i as u8)
            .collect();
        let k = labels.len();
        let quadratic = k * k.saturating_sub(1) / 2;
        let mut positions = [u8::MAX; MAX];
        for (i, &label) in labels.iter().enumerate() {
            positions[label as usize] = i as u8;
        }
        let mut monomials = Vec::new();
        let mut factors = Vec::new();
        for b in 1..k {
            for a in 0..b {
                monomials.push((1u64 << labels[a]) | (1u64 << labels[b]));
                factors.push((a as u8, b as u8));
            }
        }
        for (i, &label) in labels.iter().enumerate() {
            monomials.push(1u64 << label);
            factors.push((i as u8, k as u8));
        }
        monomials.push(0);
        factors.push((k as u8, k as u8));
        Self {
            active,
            labels,
            positions,
            quadratic,
            words: monomials.len().div_ceil(64),
            monomials,
            factors,
        }
    }
    fn coordinate(&self, monomial: u64) -> usize {
        match monomial.count_ones() {
            0 => self.quadratic + self.labels.len(),
            1 => self.quadratic + self.positions[monomial.trailing_zeros() as usize] as usize,
            2 => {
                let a = self.positions[monomial.trailing_zeros() as usize] as usize;
                let b =
                    self.positions[(monomial & (monomial - 1)).trailing_zeros() as usize] as usize;
                b * (b - 1) / 2 + a
            }
            _ => unreachable!(),
        }
    }
    fn expression(&self, value: u64, original_n: u8) -> u64 {
        let mut out = ((value >> original_n) & 1) << self.labels.len();
        let mut bits = value & ((1u64 << original_n) - 1);
        while bits != 0 {
            let j = bits.trailing_zeros() as usize;
            assert!(self.positions[j] != u8::MAX);
            out ^= 1u64 << self.positions[j];
            bits &= bits - 1;
        }
        out
    }
}
#[derive(Clone)]
enum CompactRows {
    Word(Vec<u64>),
    Wide(Vec<CoefficientRow>),
}
#[derive(Clone)]
struct CompactState {
    layout: std::rc::Rc<CompactLayout>,
    rows: CompactRows,
}
struct WordPlan {
    tables: Vec<[u64; 16]>,
}
impl WordPlan {
    fn new(images: &[u64]) -> Self {
        let mut tables = Vec::with_capacity(images.len().div_ceil(4));
        for start in (0..images.len()).step_by(4) {
            let mut table = [0u64; 16];
            for v in 1usize..16 {
                let bit = v.trailing_zeros() as usize;
                table[v] = table[v & (v - 1)] ^ images.get(start + bit).copied().unwrap_or(0);
            }
            tables.push(table);
        }
        Self { tables }
    }
    fn word(&self, row: u64) -> u64 {
        let mut result = 0;
        for (j, table) in self.tables.iter().enumerate() {
            result ^= table[((row >> (4 * j)) & 15) as usize];
        }
        result
    }
    fn wide(&self, row: &CoefficientRow) -> u64 {
        let mut result = 0;
        for (j, table) in self.tables.iter().enumerate() {
            result ^= table[((row[j / 16] >> (4 * (j % 16))) & 15) as usize];
        }
        result
    }
    fn apply(&self, rows: CompactRows) -> CompactRows {
        CompactRows::Word(match rows {
            CompactRows::Word(rows) => rows.into_iter().map(|r| self.word(r)).collect(),
            CompactRows::Wide(rows) => rows.iter().map(|r| self.wide(r)).collect(),
        })
    }
}
struct CompactBackend {
    original_n: u8,
    fast_products: bool,
    layouts: std::cell::RefCell<std::collections::BTreeMap<u64, std::rc::Rc<CompactLayout>>>,
    points: std::cell::RefCell<std::collections::BTreeMap<(u8, u8, bool), std::rc::Rc<WordPlan>>>,
}
fn local_coordinate(monomial: u64, k: usize) -> usize {
    let q = k * k.saturating_sub(1) / 2;
    match monomial.count_ones() {
        0 => q + k,
        1 => q + monomial.trailing_zeros() as usize,
        2 => {
            let a = monomial.trailing_zeros() as usize;
            let b = (monomial & (monomial - 1)).trailing_zeros() as usize;
            b * (b - 1) / 2 + a
        }
        _ => unreachable!(),
    }
}
fn compact_image(a: u64, b: u64, k: usize) -> CoefficientRow {
    let mut out = [0u64; BASIS_WORDS];
    let mut left = a;
    while left != 0 {
        let i = left.trailing_zeros() as usize;
        let x = if i == k { 0 } else { 1u64 << i };
        let mut right = b;
        while right != 0 {
            let j = right.trailing_zeros() as usize;
            let y = if j == k { 0 } else { 1u64 << j };
            toggle_bit(&mut out, local_coordinate(x | y, k));
            right &= right - 1;
        }
        left &= left - 1;
    }
    out
}
fn compact_product_word(a: u64, b: u64, k: usize) -> u64 {
    let q = k * k.saturating_sub(1) / 2;
    let mask = (1u64 << k) - 1;
    let av = a & mask;
    let bv = b & mask;
    let ac = (a >> k) & 1;
    let bc = (b >> k) & 1;
    if av == 0 {
        return if ac != 0 { b << q } else { 0 };
    }
    if bv == 0 {
        return if bc != 0 { a << q } else { 0 };
    }
    let linear = (av & bv) ^ (bv & 0u64.wrapping_sub(ac)) ^ (av & 0u64.wrapping_sub(bc));
    let mut quadratic = 0;
    if av & (av - 1) == 0 && bv & (bv - 1) == 0 {
        if av != bv {
            quadratic = 1u64 << local_coordinate(av | bv, k);
        }
    } else {
        let mut start = 0;
        for j in 1..k {
            let low =
                (bv & 0u64.wrapping_sub((av >> j) & 1)) ^ (av & 0u64.wrapping_sub((bv >> j) & 1));
            quadratic |= (low & ((1u64 << j) - 1)) << start;
            start += j;
        }
    }
    quadratic | (linear << q) | ((ac & bc) << (q + k))
}
impl CompactBackend {
    fn new(n: u8) -> Self {
        Self {
            original_n: n,
            fast_products: false,
            layouts: Default::default(),
            points: Default::default(),
        }
    }
    fn new_fast(n: u8) -> Self {
        let mut result = Self::new(n);
        result.fast_products = true;
        result
    }
    fn layout(&self, active: u64) -> std::rc::Rc<CompactLayout> {
        {
            if let Some(layout) = self.layouts.borrow().get(&active) {
                return layout.clone();
            }
        }
        let layout = std::rc::Rc::new(CompactLayout::new(active));
        let mut cache = self.layouts.borrow_mut();
        if cache.len() >= 128 {
            cache.clear();
        }
        cache.insert(active, layout.clone());
        layout
    }
    fn build_at(&self, system: &System, active: u64) -> CompactState {
        let layout = self.layout(active);
        let rows = if layout.words == 1 {
            CompactRows::Word(
                system
                    .iter()
                    .map(|p| p.iter().fold(0, |v, &m| v ^ (1u64 << layout.coordinate(m))))
                    .collect(),
            )
        } else {
            CompactRows::Wide(
                system
                    .iter()
                    .map(|p| {
                        let mut row = [0; BASIS_WORDS];
                        for &m in p {
                            toggle_bit(&mut row, layout.coordinate(m));
                        }
                        row
                    })
                    .collect(),
            )
        };
        CompactState { layout, rows }
    }
    fn images(
        &self,
        input: &CompactLayout,
        target: &CompactLayout,
        map: &RecoveryMap,
    ) -> Vec<CoefficientRow> {
        let mut expressions = Vec::with_capacity(input.labels.len() + 1);
        for &label in &input.labels {
            expressions.push(target.expression(map[label as usize], self.original_n));
        }
        expressions.push(1u64 << target.labels.len());
        input
            .factors
            .iter()
            .map(|&(a, b)| {
                compact_image(
                    expressions[a as usize],
                    expressions[b as usize],
                    target.labels.len(),
                )
            })
            .collect()
    }
    fn word_images(
        &self,
        input: &CompactLayout,
        target: &CompactLayout,
        map: &RecoveryMap,
    ) -> Vec<u64> {
        let mut expressions = Vec::with_capacity(input.labels.len() + 1);
        for &j in &input.labels {
            expressions.push(target.expression(map[j as usize], self.original_n));
        }
        expressions.push(1u64 << target.labels.len());
        input
            .factors
            .iter()
            .map(|&(a, b)| {
                compact_product_word(
                    expressions[a as usize],
                    expressions[b as usize],
                    target.labels.len(),
                )
            })
            .collect()
    }
    fn to_wide(&self, input: CompactRows, images: &[CoefficientRow], words: usize) -> CompactRows {
        let convert = |row: &[u64]| {
            let mut out = [0u64; BASIS_WORDS];
            for (w, &word) in row.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let image = &images[w * 64 + bits.trailing_zeros() as usize];
                    for j in 0..words {
                        out[j] ^= image[j];
                    }
                    bits &= bits - 1;
                }
            }
            out
        };
        CompactRows::Wide(match input {
            CompactRows::Word(rows) => rows.into_iter().map(|r| convert(&[r])).collect(),
            CompactRows::Wide(rows) => rows.iter().map(|r| convert(r)).collect(),
        })
    }
    fn point_plan(&self, k: usize, j: usize, one: bool) -> std::rc::Rc<WordPlan> {
        let key = (k as u8, j as u8, one);
        {
            if let Some(plan) = self.points.borrow().get(&key) {
                return plan.clone();
            }
        }
        assert!((1..=11).contains(&k));
        let target = k - 1;
        let expressions: Vec<_> = (0..k)
            .map(|i| {
                if i == j {
                    if one {
                        1u64 << target
                    } else {
                        0
                    }
                } else {
                    1u64 << if i < j { i } else { i - 1 }
                }
            })
            .chain(std::iter::once(1u64 << target))
            .collect();
        let source = CompactLayout::new((1u64 << k) - 1);
        let images: Vec<_> = source
            .factors
            .iter()
            .map(|&(a, b)| {
                if self.fast_products {
                    compact_product_word(expressions[a as usize], expressions[b as usize], target)
                } else {
                    compact_image(expressions[a as usize], expressions[b as usize], target)[0]
                }
            })
            .collect();
        let plan = std::rc::Rc::new(WordPlan::new(&images));
        let mut cache = self.points.borrow_mut();
        if cache.len() >= 256 {
            cache.clear();
        }
        cache.insert(key, plan.clone());
        plan
    }
    #[cfg(test)]
    fn materialize(&self, state: &CompactState) -> System {
        (0..self.row_count(state))
            .map(|r| {
                let mut p = Vec::new();
                self.terms(state, r, |m| p.push(m));
                p
            })
            .collect()
    }
}
impl BasisBackend for CompactBackend {
    type State = CompactState;
    fn build(&self, system: &System) -> Self::State {
        self.build_at(system, (1u64 << self.original_n) - 1)
    }
    fn row_count(&self, state: &Self::State) -> usize {
        match &state.rows {
            CompactRows::Word(r) => r.len(),
            CompactRows::Wide(r) => r.len(),
        }
    }
    fn term_count(&self, state: &Self::State, r: usize) -> usize {
        match &state.rows {
            CompactRows::Word(rows) => rows[r].count_ones() as usize,
            CompactRows::Wide(rows) => rows[r][..state.layout.words]
                .iter()
                .map(|v| v.count_ones() as usize)
                .sum(),
        }
    }
    fn terms(&self, state: &Self::State, r: usize, mut visit: impl FnMut(u64)) {
        let mut word = |w: usize, mut bits: u64| {
            while bits != 0 {
                visit(state.layout.monomials[w * 64 + bits.trailing_zeros() as usize]);
                bits &= bits - 1;
            }
        };
        match &state.rows {
            CompactRows::Word(rows) => word(0, rows[r]),
            CompactRows::Wide(rows) => {
                for w in 0..state.layout.words {
                    word(w, rows[r][w]);
                }
            }
        }
    }
    fn column_count(&self, state: &Self::State) -> usize {
        match &state.rows {
            CompactRows::Word(rows) => rows.iter().fold(0u64, |a, b| a | b).count_ones() as usize,
            CompactRows::Wide(rows) => {
                let mut union = [0u64; BASIS_WORDS];
                for r in rows {
                    for w in 0..state.layout.words {
                        union[w] |= r[w];
                    }
                }
                union.iter().map(|v| v.count_ones() as usize).sum()
            }
        }
    }
    fn reduce(&self, state: Self::State) -> Self::State {
        let layout = state.layout;
        let rows = match state.rows {
            CompactRows::Word(rows) => {
                let mut pivots = [0u64; 64];
                let mut present = 0u64;
                for mut row in rows {
                    while row != 0 {
                        let p = row.trailing_zeros() as usize;
                        if pivots[p] == 0 {
                            pivots[p] = row;
                            present |= 1u64 << p;
                            break;
                        }
                        row ^= pivots[p];
                    }
                }
                let affine = present & !((1u64 << layout.quadratic) - 1);
                let mut upper = affine;
                while upper != 0 {
                    let p = 63 - upper.leading_zeros();
                    let bit = 1u64 << p;
                    let mut lower = affine & (bit - 1);
                    while lower != 0 {
                        let q = lower.trailing_zeros() as usize;
                        if pivots[q] & bit != 0 {
                            pivots[q] ^= pivots[p as usize];
                        }
                        lower &= lower - 1;
                    }
                    upper &= !bit;
                }
                let mut out = Vec::with_capacity(present.count_ones() as usize);
                while present != 0 {
                    out.push(pivots[present.trailing_zeros() as usize]);
                    present &= present - 1;
                }
                CompactRows::Word(out)
            }
            CompactRows::Wide(rows) => {
                let mut pivots = [-1i8; BASIS_COLUMNS];
                let mut out: Vec<CoefficientRow> = Vec::with_capacity(rows.len());
                for mut row in rows {
                    while let Some(p) = first_column(&row, layout.words) {
                        if pivots[p] < 0 {
                            pivots[p] = out.len() as i8;
                            out.push(row);
                            break;
                        }
                        let pivot = &out[pivots[p] as usize];
                        for w in p / 64..layout.words {
                            row[w] ^= pivot[w];
                        }
                    }
                }
                out.sort_unstable_by_key(|r| first_column(r, layout.words).unwrap());
                for p in (0..out.len()).rev() {
                    let column = first_column(&out[p], layout.words).unwrap();
                    if column < layout.quadratic {
                        continue;
                    }
                    let pivot = out[p];
                    for q in 0..p {
                        if first_column(&out[q], layout.words).unwrap() >= layout.quadratic
                            && row_bit(&out[q], column)
                        {
                            for w in column / 64..layout.words {
                                out[q][w] ^= pivot[w];
                            }
                        }
                    }
                }
                CompactRows::Wide(out)
            }
        };
        CompactState { layout, rows }
    }
    fn specialize(
        &self,
        mut state: Self::State,
        mask: u64,
        values: u64,
        active: u64,
    ) -> Self::State {
        let mut assigned = mask & active & state.layout.active;
        while assigned != 0 {
            let bit = assigned & assigned.wrapping_neg();
            let global = bit.trailing_zeros() as usize;
            let j = state.layout.positions[global] as usize;
            let k = state.layout.labels.len();
            let one = values & bit != 0;
            let target = self.layout(state.layout.active & !bit);
            let rows = if target.words == 1 {
                self.point_plan(k, j, one).apply(state.rows)
            } else {
                let mut map = identity_affine(state.layout.active, self.original_n);
                map[global] = if one { 1u64 << self.original_n } else { 0 };
                let images = self.images(&state.layout, &target, &map);
                self.to_wide(state.rows, &images, target.words)
            };
            state = CompactState {
                layout: target,
                rows,
            };
            assigned &= assigned - 1;
        }
        state
    }
}
impl AffineBackend for CompactBackend {
    fn transform_pair(
        &self,
        basis: Self::State,
        source: Self::State,
        map: &RecoveryMap,
        _n: u8,
    ) -> (Self::State, Self::State) {
        assert_eq!(basis.layout.active, source.layout.active);
        let mut active = 0;
        for &j in &basis.layout.labels {
            active |= map[j as usize] & ((1u64 << self.original_n) - 1);
        }
        let target = self.layout(active);
        let (a, b) = if target.words == 1 {
            let words: Vec<_> = if self.fast_products {
                self.word_images(&basis.layout, &target, map)
            } else {
                self.images(&basis.layout, &target, map)
                    .iter()
                    .map(|r| r[0])
                    .collect()
            };
            let plan = WordPlan::new(&words);
            (plan.apply(basis.rows), plan.apply(source.rows))
        } else {
            let images = self.images(&basis.layout, &target, map);
            (
                self.to_wide(basis.rows, &images, target.words),
                self.to_wide(source.rows, &images, target.words),
            )
        };
        (
            CompactState {
                layout: target.clone(),
                rows: a,
            },
            CompactState {
                layout: target,
                rows: b,
            },
        )
    }
    fn clean_source(&self, state: Self::State) -> Self::State {
        let rows = match state.rows {
            CompactRows::Word(rows) => {
                let mut out = Vec::new();
                for row in rows {
                    if row != 0 && !out.contains(&row) {
                        out.push(row);
                    }
                }
                CompactRows::Word(out)
            }
            CompactRows::Wide(rows) => {
                let mut out = Vec::new();
                for row in rows {
                    if row[..state.layout.words].iter().any(|&v| v != 0) && !out.contains(&row) {
                        out.push(row);
                    }
                }
                CompactRows::Wide(out)
            }
        };
        CompactState {
            layout: state.layout,
            rows,
        }
    }
}
fn solve_affine_compact(system: &System, n: u8, limit: u64) -> Solved {
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
    solve_affine_backend(system, n, limit, CompactBackend::new(n))
}

fn solve_affine_fast(system: &System, n: u8, limit: u64) -> Solved {
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
    solve_affine_backend(system, n, limit, CompactBackend::new_fast(n))
}
