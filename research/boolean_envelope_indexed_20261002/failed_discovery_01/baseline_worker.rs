//! Standalone, bounded parameterized Boolean support-envelope experiment. Uses only std.
use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;
use std::mem::size_of;
use std::time::Instant;

type Poly = Vec<u16>;
#[derive(Clone, Debug, PartialEq, Eq)]
struct Input {
    n: u8,
    degree: u8,
    active: u16,
    polys: Vec<Poly>,
}
#[derive(Clone, Copy)]
struct Caps {
    rows: usize,
    columns: usize,
}
const CAPS: Caps = Caps {
    rows: 512,
    columns: 512,
};
#[derive(Clone, Debug, PartialEq, Eq)]
struct Matrix {
    columns: Vec<u16>,
    rows: Vec<Vec<u64>>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
enum Error {
    InvalidInput,
    RowCap,
    ColumnCap,
    PlanCap,
}
type Built = Result<Matrix, Error>;

impl Input {
    fn validate(&self) -> Result<(), Error> {
        if !(1..=12).contains(&self.n)
            || !(1..=4).contains(&self.degree)
            || self.polys.len() > 12
            || self.active >= (1 << self.n)
        {
            return Err(Error::InvalidInput);
        }
        for poly in &self.polys {
            if poly.len() > 32
                || poly.iter().any(|&m| m >= (1 << self.n))
                || poly.windows(2).any(|w| w[0] >= w[1])
            {
                return Err(Error::InvalidInput);
            }
        }
        Ok(())
    }
    fn retained(&self) -> usize {
        size_of::<Self>()
            + self.polys.capacity() * size_of::<Poly>()
            + self
                .polys
                .iter()
                .map(|p| p.capacity() * size_of::<u16>())
                .sum::<usize>()
    }
}
impl Matrix {
    fn retained(&self) -> usize {
        size_of::<Self>()
            + self.columns.capacity() * size_of::<u16>()
            + self.rows.capacity() * size_of::<Vec<u64>>()
            + self
                .rows
                .iter()
                .map(|r| r.capacity() * size_of::<u64>())
                .sum::<usize>()
    }
    fn payload(&self) -> usize {
        self.columns.len() * 2 + self.rows.iter().map(|r| r.len() * 8).sum::<usize>()
    }
    fn check_caps(&self, caps: Caps) -> Result<(), Error> {
        check_shape(self.rows.len(), self.columns.len(), caps)
    }
}
fn check_shape(rows: usize, cols: usize, caps: Caps) -> Result<(), Error> {
    if rows > caps.rows {
        return Err(Error::RowCap);
    }
    if cols > caps.columns {
        return Err(Error::ColumnCap);
    }
    Ok(())
}
fn multipliers(active: u16, gap: u32) -> Vec<u16> {
    let mut out = vec![0u16];
    for bit in 0..12 {
        if active & (1 << bit) == 0 {
            continue;
        }
        let old_len = out.len();
        for i in 0..old_len {
            if out[i].count_ones() < gap {
                out.push(out[i] | (1 << bit));
            }
        }
    }
    out.sort_unstable();
    out
}

// The measured constructor sorts products and cancels even multiplicities.
fn products(input: &Input, caps: Caps) -> Result<Vec<Poly>, Error> {
    input.validate()?;
    let mut schedules = BTreeMap::new();
    let mut rows = Vec::new();
    for poly in &input.polys {
        let degree = poly.iter().map(|m| m.count_ones()).max().unwrap_or(0);
        if degree > u32::from(input.degree) || poly.is_empty() {
            continue;
        }
        let gap = u32::from(input.degree) - degree;
        let schedule = schedules
            .entry(gap)
            .or_insert_with(|| multipliers(input.active, gap));
        for &mult in schedule.iter() {
            let mut terms: Vec<_> = poly.iter().map(|&term| term | mult).collect();
            terms.sort_unstable();
            let mut row = Vec::new();
            let mut i = 0;
            while i < terms.len() {
                let mut end = i + 1;
                while end < terms.len() && terms[end] == terms[i] {
                    end += 1;
                }
                if (end - i) % 2 != 0 {
                    row.push(terms[i]);
                }
                i = end;
            }
            if !row.is_empty() {
                if rows.len() == caps.rows {
                    return Err(Error::RowCap);
                }
                rows.push(row);
            }
        }
    }
    Ok(rows)
}
fn columns(rows: &[Poly], caps: Caps) -> Result<Vec<u16>, Error> {
    let mut columns: Vec<_> = rows.iter().flatten().copied().collect();
    columns.sort_unstable();
    columns.dedup();
    check_shape(rows.len(), columns.len(), caps)?;
    Ok(columns)
}
fn pack(rows: &[Poly], columns: &[u16]) -> Matrix {
    let words = columns.len().div_ceil(64);
    let rows = rows
        .iter()
        .map(|terms| {
            let mut row = vec![0; words];
            for term in terms {
                let index = columns.binary_search(term).unwrap();
                row[index / 64] |= 1 << (index % 64);
            }
            row
        })
        .collect();
    Matrix {
        columns: columns.to_vec(),
        rows,
    }
}
fn direct(input: &Input, caps: Caps) -> Built {
    let rows = products(input, caps)?;
    Ok(pack(&rows, &columns(&rows, caps)?))
}

// Independent correctness oracle: enumerate submasks, use set symmetric
// difference for parity, and test each matrix entry instead of setting bits.
fn oracle(input: &Input, caps: Caps) -> Built {
    input.validate()?;
    let mut rows = Vec::new();
    let mut support = BTreeSet::new();
    for poly in &input.polys {
        let d = poly.iter().map(|m| m.count_ones()).max().unwrap_or(0);
        if d > u32::from(input.degree) || poly.is_empty() {
            continue;
        }
        for mult in 0u16..(1 << input.n) {
            if mult & !input.active != 0 || mult.count_ones() > u32::from(input.degree) - d {
                continue;
            }
            let mut row = BTreeSet::new();
            for &term in poly {
                if !row.insert(term | mult) {
                    row.remove(&(term | mult));
                }
            }
            if !row.is_empty() {
                support.extend(row.iter().copied());
                rows.push(row);
            }
        }
    }
    check_shape(rows.len(), support.len(), caps)?;
    let columns: Vec<_> = support.into_iter().collect();
    let rows = rows
        .iter()
        .map(|terms| {
            columns
                .chunks(64)
                .map(|chunk| {
                    chunk
                        .iter()
                        .enumerate()
                        .fold(0, |word, (i, m)| word | (u64::from(terms.contains(m)) << i))
                })
                .collect()
        })
        .collect();
    Ok(Matrix { columns, rows })
}

struct Schedule {
    signature: Input,
    columns: Vec<u16>,
    rows: Vec<Vec<u16>>,
}
impl Schedule {
    fn compile(input: &Input, caps: Caps) -> Result<Self, Error> {
        let products = products(input, caps)?;
        let mut columns = columns(&products, caps)?;
        // Keep the compact support, not the product-discovery capacity.
        // The matrix-cache control's pack() already copies compact columns.
        columns.shrink_to_fit();
        let rows = products
            .iter()
            .map(|row| {
                row.iter()
                    .map(|m| columns.binary_search(m).unwrap() as u16)
                    .collect()
            })
            .collect();
        Ok(Self {
            signature: input.clone(),
            columns,
            rows,
        })
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        if input != &self.signature {
            return (direct(input, caps), false);
        }
        if let Err(error) = check_shape(self.rows.len(), self.columns.len(), caps) {
            return (Err(error), true);
        }
        let words = self.columns.len().div_ceil(64);
        let rows = self
            .rows
            .iter()
            .map(|slots| {
                let mut row = vec![0u64; words];
                for &slot in slots {
                    row[slot as usize / 64] |= 1 << (slot % 64);
                }
                row
            })
            .collect();
        (
            Ok(Matrix {
                columns: self.columns.clone(),
                rows,
            }),
            true,
        )
    }
    fn retained(&self) -> usize {
        size_of::<Self>() - size_of::<Input>()
            + self.signature.retained()
            + self.columns.capacity() * 2
            + self.rows.capacity() * size_of::<Vec<u16>>()
            + self.rows.iter().map(|r| r.capacity() * 2).sum::<usize>()
    }
}
#[derive(Clone, Copy)]
struct PlanCaps {
    rows: usize,
    groups: usize,
    retained_bytes: usize,
}
const PLAN_CAPS: PlanCaps = PlanCaps {
    rows: 8192,
    groups: 131072,
    retained_bytes: 2 * 1024 * 1024,
};

// A group is the linear form in generator coefficients contributing to one
// product monomial. Boolean multiplication uses union; addition uses parity.
struct Group {
    monomial: u16,
    coefficients: u32,
}
struct ProductPlan {
    multiplier_degree: u8,
    groups: Vec<Group>,
}
struct Envelope {
    support: Input,
    plans: Vec<Vec<ProductPlan>>,
}
impl Envelope {
    fn compile(support: &Input, caps: PlanCaps) -> Result<Self, Error> {
        support.validate()?;
        let mut plans = Vec::new();
        let (mut row_count, mut group_count) = (0, 0);
        for terms in &support.polys {
            let mut generator = Vec::new();
            if let Some(min_degree) = terms.iter().map(|m| m.count_ones()).min() {
                if min_degree <= u32::from(support.degree) {
                    for multiplier in
                        multipliers(support.active, u32::from(support.degree) - min_degree)
                    {
                        if row_count == caps.rows {
                            return Err(Error::PlanCap);
                        }
                        row_count += 1;
                        let k = multiplier.count_ones();
                        let mut grouped = BTreeMap::<u16, u32>::new();
                        for (slot, &term) in terms.iter().enumerate() {
                            // This term cannot be active in any generator for
                            // which this multiplier is eligible otherwise.
                            if term.count_ones() <= u32::from(support.degree) - k {
                                *grouped.entry(term | multiplier).or_default() ^= 1u32 << slot;
                            }
                        }
                        group_count += grouped.len();
                        if group_count > caps.groups {
                            return Err(Error::PlanCap);
                        }
                        generator.push(ProductPlan {
                            multiplier_degree: k as u8,
                            groups: grouped
                                .into_iter()
                                .map(|(monomial, coefficients)| Group {
                                    monomial,
                                    coefficients,
                                })
                                .collect(),
                        });
                    }
                }
            }
            plans.push(generator);
        }
        let result = Self {
            support: support.clone(),
            plans,
        };
        if result.retained() > caps.retained_bytes {
            return Err(Error::PlanCap);
        }
        Ok(result)
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        if let Err(error) = input.validate() {
            return (Err(error), false);
        }
        if input.n != self.support.n
            || input.degree != self.support.degree
            || input.active != self.support.active
            || input.polys.len() != self.support.polys.len()
        {
            return (direct(input, caps), false);
        }
        // Validate every generator before producing anything. Generator slots
        // are ordered; each has its own independently declared support envelope.
        let mut coefficients = Vec::with_capacity(input.polys.len());
        for (poly, support) in input.polys.iter().zip(&self.support.polys) {
            let mut bits = 0u32;
            for term in poly {
                let Ok(slot) = support.binary_search(term) else {
                    return (direct(input, caps), false);
                };
                bits |= 1u32 << slot;
            }
            coefficients.push(bits);
        }
        let mut rows = Vec::new();
        for ((poly, &bits), plans) in input.polys.iter().zip(&coefficients).zip(&self.plans) {
            let Some(degree) = poly.iter().map(|m| m.count_ones()).max() else {
                continue;
            };
            if degree > u32::from(input.degree) {
                continue;
            }
            for plan in plans {
                if u32::from(plan.multiplier_degree) > u32::from(input.degree) - degree {
                    continue;
                }
                let row: Vec<_> = plan
                    .groups
                    .iter()
                    .filter(|g| (g.coefficients & bits).count_ones() % 2 == 1)
                    .map(|g| g.monomial)
                    .collect();
                if !row.is_empty() {
                    if rows.len() == caps.rows {
                        return (Err(Error::RowCap), true);
                    }
                    rows.push(row);
                }
            }
        }
        // The envelope is not the actual output support. Recompact it on every
        // application, including after cancellation or disappearance of rows.
        (columns(&rows, caps).map(|cols| pack(&rows, &cols)), true)
    }
    fn retained(&self) -> usize {
        size_of::<Self>() - size_of::<Input>()
            + self.support.retained()
            + self.plans.capacity() * size_of::<Vec<ProductPlan>>()
            + self
                .plans
                .iter()
                .map(|plans| {
                    plans.capacity() * size_of::<ProductPlan>()
                        + plans
                            .iter()
                            .map(|p| p.groups.capacity() * size_of::<Group>())
                            .sum::<usize>()
                })
                .sum::<usize>()
    }
}

enum Cache {
    Direct,
    Layout(Vec<u16>),
    Schedule(Schedule),
    Matrix(Input, Matrix),
    Envelope(Envelope),
}
impl Cache {
    fn compile(variant: usize, input: &Input, envelope: &Input) -> Self {
        match variant {
            0 => Self::Direct,
            1 => Self::Layout(columns(&products(input, CAPS).unwrap(), CAPS).unwrap()),
            2 => Self::Schedule(Schedule::compile(input, CAPS).unwrap()),
            3 => Self::Matrix(input.clone(), direct(input, CAPS).unwrap()),
            4 => Self::Envelope(Envelope::compile(envelope, PLAN_CAPS).unwrap()),
            _ => unreachable!(),
        }
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        match self {
            Self::Direct => (direct(input, caps), false),
            Self::Schedule(schedule) => schedule.apply(input, caps),
            Self::Envelope(schedule) => schedule.apply(input, caps),
            Self::Matrix(signature, matrix) => {
                if input != signature {
                    (direct(input, caps), false)
                } else {
                    (matrix.check_caps(caps).map(|_| matrix.clone()), true)
                }
            }
            Self::Layout(layout) => {
                let rows = match products(input, caps) {
                    Ok(rows) => rows,
                    Err(e) => return (Err(e), false),
                };
                // Validate exact support, including absent cached columns.
                let mut seen = vec![false; layout.len()];
                let fits = layout.len() <= caps.columns
                    && rows.iter().flatten().all(|m| {
                        if let Ok(i) = layout.binary_search(m) {
                            seen[i] = true;
                            true
                        } else {
                            false
                        }
                    })
                    && seen.iter().all(|x| *x);
                if fits {
                    (Ok(pack(&rows, layout)), true)
                } else {
                    (columns(&rows, caps).map(|cols| pack(&rows, &cols)), false)
                }
            }
        }
    }
    fn retained(&self) -> usize {
        match self {
            Self::Direct => 0,
            Self::Layout(cols) => size_of::<Vec<u16>>() + cols.capacity() * 2,
            Self::Schedule(s) => s.retained(),
            Self::Matrix(input, matrix) => input.retained() + matrix.retained(),
            Self::Envelope(s) => s.retained(),
        }
    }
}
fn next(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *seed;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn fixture(n: u8, mut seed: u64) -> Input {
    let mut polys = Vec::new();
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < n as usize + 4 {
            let a = (next(&mut seed) % u64::from(n)) as u8;
            let b = (next(&mut seed) % u64::from(n)) as u8;
            terms.insert((1 << a) | (1 << b));
        }
        polys.push(terms.into_iter().collect());
    }
    Input {
        n,
        degree: 3,
        active: (1 << n) - 1,
        polys,
    }
}

// Envelope construction depends only on the base fixture, never on a later
// coefficient assignment. Every slot permits a constant and all linear terms.
fn declared_envelope(base: &Input) -> Input {
    let mut out = base.clone();
    for poly in &mut out.polys {
        poly.extend(std::iter::once(0).chain((0..base.n).map(|i| 1 << i)));
        poly.sort_unstable();
        poly.dedup();
    }
    out.validate().unwrap();
    out
}
fn assignments(
    base: &Input,
    envelope: &Input,
    seed: u64,
    batch: usize,
    family: &str,
) -> Vec<Input> {
    assert!(["repeat", "coefficients", "degree_cycle", "escape"].contains(&family));
    let mut state = seed ^ 0x6334be35a74c2189;
    (0..batch)
        .map(|i| {
            if i == 0 || family == "repeat" {
                return base.clone();
            }
            let mut input = base.clone();
            for (poly, support) in input.polys.iter_mut().zip(&envelope.polys) {
                *poly = support
                    .iter()
                    .copied()
                    .filter(|_| next(&mut state) & 1 != 0)
                    .collect();
                // Keep every generator quadratic unless this family explicitly
                // drops the first one. This gives a bounded matrix at all sizes.
                if !poly.iter().any(|m| m.count_ones() == 2) {
                    poly.push(*support.iter().find(|m| m.count_ones() == 2).unwrap());
                    poly.sort_unstable();
                }
            }
            if family == "degree_cycle" {
                input.polys[0] = match i % 4 {
                    1 => (0..base.n)
                        .filter(|_| next(&mut state) & 1 != 0)
                        .map(|j| 1 << j)
                        .chain(std::iter::once(1))
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                    2 => vec![0],
                    3 => vec![],
                    _ => input.polys[0].clone(),
                };
            }
            if family == "escape" && i % 4 == 3 {
                input.polys[0].push(7); // a cubic term excluded from every envelope
                input.polys[0].sort_unstable();
            }
            assert_ne!(input, *base);
            input
        })
        .collect()
}
fn eligible(input: &Input, poly: &Poly) -> Vec<u16> {
    match poly.iter().map(|m| m.count_ones()).max() {
        Some(d) if d <= u32::from(input.degree) => {
            multipliers(input.active, u32::from(input.degree) - d)
        }
        _ => vec![],
    }
}
fn newly_required(base: &Input, input: &Input) -> usize {
    base.polys
        .iter()
        .zip(&input.polys)
        .map(|(before, after)| {
            let old = eligible(base, before);
            eligible(input, after)
                .iter()
                .filter(|m| old.binary_search(m).is_err())
                .count()
        })
        .sum()
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 6, "worker N SEED BATCH FAMILY REPETITIONS");
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let batch = args[3].parse::<usize>().unwrap();
    let family = &args[4];
    let repetitions = args[5].parse::<usize>().unwrap();
    assert!([6, 8, 10, 12].contains(&n));
    assert!((1..=64).contains(&batch) && (1..=32).contains(&repetitions));
    let base = fixture(n, seed);
    let envelope = declared_envelope(&base);
    let inputs = assignments(&base, &envelope, seed, batch, family);
    if family != "repeat" {
        let distinct: BTreeSet<_> = inputs.iter().map(|x| &x.polys).collect();
        assert_eq!(distinct.len(), batch);
    }
    let expected: Vec<_> = inputs.iter().map(|x| oracle(x, CAPS).unwrap()).collect();
    let new_multipliers: usize = inputs.iter().map(|x| newly_required(&base, x)).sum();
    let polys: Vec<_> = inputs.iter().map(|x| &x.polys).collect();
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"batch\":{batch},\"family\":\"{family}\",\"degree\":3,\"active\":{},\"base\":{:?},\"envelope\":{:?},\"inputs\":{:?},\"newly_required_multipliers\":{new_multipliers}}}", base.active, base.polys, envelope.polys, polys);
    let names = ["direct", "layout", "schedule", "matrix_cache", "envelope"];
    for rep in 0..repetitions {
        for order in 0..names.len() {
            let variant = (rep + order) % names.len();
            let start = Instant::now();
            let cache = Cache::compile(variant, black_box(&base), black_box(&envelope));
            let setup_ns = start.elapsed().as_nanos();
            let retained_bytes = cache.retained();
            let (mut hits, mut fallbacks, mut changed_hits) = (0, 0, 0);
            let (mut apply_ns, mut validation_ns, mut output_bytes) = (0u128, 0u128, 0usize);
            for (input, expected) in inputs.iter().zip(&expected) {
                let tick = Instant::now();
                let (got, hit) = cache.apply(black_box(input), CAPS);
                let got = got.unwrap();
                apply_ns += tick.elapsed().as_nanos();
                hits += usize::from(hit);
                fallbacks += usize::from(variant != 0 && !hit);
                changed_hits += usize::from(hit && input != &base);
                let tick = Instant::now();
                assert_eq!(black_box(&got), expected);
                output_bytes += got.payload();
                validation_ns += tick.elapsed().as_nanos();
            }
            drop(cache);
            let total_ns = start.elapsed().as_nanos();
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{}\",\"setup_ns\":{setup_ns},\"apply_ns\":{apply_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"retained_bytes\":{retained_bytes},\"output_bytes\":{output_bytes},\"hits\":{hits},\"fallbacks\":{fallbacks},\"changed_hits\":{changed_hits},\"verified_outputs\":{batch}}}", names[variant]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn small(polys: Vec<Poly>, degree: u8) -> Input {
        Input {
            n: 3,
            active: 7,
            degree,
            polys,
        }
    }
    #[test]
    fn exhaustive_coefficients_degrees_and_active_masks() {
        // Compile once per context, then reuse on all 256 assignments. Includes
        // zero, constants, degree drops, and all Boolean product cancellations.
        for degree in 1..=4 {
            for active in 0..8 {
                let support = Input {
                    active,
                    ..small(vec![(0..8).collect()], degree)
                };
                let plan = Envelope::compile(&support, PLAN_CAPS).unwrap();
                for coefficients in 0..256u16 {
                    let input = Input {
                        polys: vec![(0..8).filter(|m| coefficients & (1 << m) != 0).collect()],
                        ..support.clone()
                    };
                    let expected = oracle(&input, CAPS);
                    assert_eq!(direct(&input, CAPS), expected);
                    assert_eq!(plan.apply(&input, CAPS), (expected, true));
                }
            }
        }
    }
    #[test]
    fn cancellations_degree_drops_and_zero_rows() {
        let support = small(vec![vec![0, 1, 2, 3]], 3);
        let plan = Envelope::compile(&support, PLAN_CAPS).unwrap();
        let quadratic = small(vec![vec![1, 3]], 3);
        let linear = small(vec![vec![0, 1]], 3);
        let constant = small(vec![vec![0]], 3);
        let zero = small(vec![vec![]], 3);
        assert_eq!(newly_required(&quadratic, &linear), 3);
        assert_eq!(newly_required(&quadratic, &constant), 4);
        assert!(products(&linear, CAPS).unwrap().len() < eligible(&linear, &linear.polys[0]).len());
        for input in [&quadratic, &linear, &constant, &zero] {
            assert_eq!(plan.apply(input, CAPS), (oracle(input, CAPS), true));
        }
        assert_eq!(
            plan.apply(&zero, CAPS).0.unwrap(),
            Matrix {
                columns: vec![],
                rows: vec![]
            }
        );
    }
    #[test]
    fn support_escape_and_context_changes_use_direct_fallback() {
        let support = small(vec![vec![0, 1, 3], vec![2]], 3);
        let plan = Envelope::compile(&support, PLAN_CAPS).unwrap();
        let changes = [
            small(vec![vec![7], vec![2]], 3),
            small(vec![vec![2], vec![0, 1, 3]], 3),
            small(vec![vec![0, 1, 3]], 3),
            Input {
                degree: 2,
                ..support.clone()
            },
            Input {
                n: 4,
                ..support.clone()
            },
            Input {
                active: 3,
                ..support.clone()
            },
        ];
        for input in changes {
            assert_eq!(plan.apply(&input, CAPS), (oracle(&input, CAPS), false));
        }
        assert_eq!(plan.apply(&support, CAPS), (oracle(&support, CAPS), true));
        let invalid = small(vec![vec![1, 1]], 3);
        assert_eq!(
            plan.apply(&invalid, CAPS),
            (Err(Error::InvalidInput), false)
        );
    }
    #[test]
    fn current_output_caps_apply_after_parity_and_compaction() {
        let support = small(vec![(0..8).collect()], 3);
        let empty = small(vec![vec![]], 3);
        let constant = small(vec![vec![0]], 3);
        for variant in 0..5 {
            let cache = Cache::compile(variant, &support, &support);
            for input in [&support, &empty, &constant] {
                for caps in [
                    CAPS,
                    Caps {
                        rows: 0,
                        columns: 0,
                    },
                    Caps {
                        rows: 512,
                        columns: 0,
                    },
                    Caps {
                        rows: 1,
                        columns: 8,
                    },
                ] {
                    assert_eq!(cache.apply(input, caps).0, oracle(input, caps));
                }
            }
        }
    }
    #[test]
    fn compilation_limits_are_explicit_refusals() {
        let support = small(vec![(0..8).collect()], 3);
        for caps in [
            PlanCaps {
                rows: 0,
                ..PLAN_CAPS
            },
            PlanCaps {
                groups: 0,
                ..PLAN_CAPS
            },
            PlanCaps {
                retained_bytes: 0,
                ..PLAN_CAPS
            },
        ] {
            assert!(matches!(
                Envelope::compile(&support, caps),
                Err(Error::PlanCap)
            ));
        }
    }
    #[test]
    fn coefficient_slot_31_and_empty_or_high_degree_envelopes() {
        let support = Input {
            n: 5,
            active: 31,
            degree: 4,
            polys: vec![(0..32).collect(), vec![], vec![31]],
        };
        let plan = Envelope::compile(&support, PLAN_CAPS).unwrap();
        for poly in [vec![31], vec![0, 31], vec![0], vec![]] {
            let input = Input {
                polys: vec![poly, vec![], vec![31]],
                ..support.clone()
            };
            assert_eq!(plan.apply(&input, CAPS), (oracle(&input, CAPS), true));
        }
    }
    #[test]
    fn changed_assignment_is_recomputed_and_all_control_arms_agree() {
        let base = fixture(6, 17);
        let envelope = declared_envelope(&base);
        for family in ["repeat", "coefficients", "degree_cycle", "escape"] {
            let inputs = assignments(&base, &envelope, 17, 8, family);
            for variant in 0..5 {
                let cache = Cache::compile(variant, &base, &envelope);
                for (i, input) in inputs.iter().enumerate() {
                    let (got, hit) = cache.apply(input, CAPS);
                    assert_eq!(got, oracle(input, CAPS));
                    if variant == 4 {
                        assert_eq!(hit, family != "escape" || i % 4 != 3);
                    }
                    if [2, 3].contains(&variant) {
                        assert_eq!(hit, i == 0 || family == "repeat");
                    }
                }
            }
            if family != "repeat" {
                assert_ne!(direct(&base, CAPS), direct(&inputs[1], CAPS));
            }
        }
    }
}
