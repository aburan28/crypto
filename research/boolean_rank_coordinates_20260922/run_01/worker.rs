//! Standalone, bounded combinatorial-coordinate Boolean matrix construction experiment. Uses only std.
use std::collections::{BTreeMap, BTreeSet};
use std::hint::black_box;
use std::mem::size_of;
use std::time::Instant;

type Poly = Vec<u64>;
#[derive(Clone, Debug, PartialEq, Eq)]
struct Input {
    n: u8,
    degree: u8,
    active: u64,
    polys: Vec<Poly>,
}
#[derive(Clone, Copy)]
struct Caps {
    rows: usize,
    columns: usize,
}
const CAPS: Caps = Caps {
    rows: 4096,
    columns: 8192,
};
#[derive(Clone, Debug, PartialEq, Eq)]
struct Matrix {
    columns: Vec<u64>,
    rows: Vec<Vec<u64>>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
enum Error {
    InvalidInput,
    RowCap,
    ColumnCap,
    DenseLimit,
}
type Built = Result<Matrix, Error>;

impl Input {
    fn validate(&self) -> Result<(), Error> {
        if !(1..=36).contains(&self.n)
            || !(1..=3).contains(&self.degree)
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
}
impl Matrix {
    fn payload(&self) -> usize {
        self.columns.len() * 8 + self.rows.iter().map(|r| r.len() * 8).sum::<usize>()
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
fn multipliers(active: u64, gap: u32) -> Vec<u64> {
    let mut out = vec![0u64];
    for bit in 0..36 {
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
fn columns(rows: &[Poly], caps: Caps) -> Result<Vec<u64>, Error> {
    let mut columns: Vec<_> = rows.iter().flatten().copied().collect();
    columns.sort_unstable();
    columns.dedup();
    check_shape(rows.len(), columns.len(), caps)?;
    Ok(columns)
}
fn pack(rows: &[Poly], columns: &[u64]) -> Matrix {
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

fn retain_row(
    rows: &mut Vec<Vec<u64>>,
    row: Vec<u64>,
    occupied: &mut [u64],
    caps: Caps,
) -> Result<(), Error> {
    if row.iter().all(|&w| w == 0) {
        return Ok(());
    }
    if rows.len() == caps.rows {
        return Err(Error::RowCap);
    }
    for (seen, &word) in occupied.iter_mut().zip(&row) {
        *seen |= word;
    }
    rows.push(row);
    Ok(())
}
fn compact_packed(rows: Vec<Vec<u64>>, basis: &[u64], occupied: &[u64], caps: Caps) -> Built {
    let count = occupied
        .iter()
        .map(|w| w.count_ones() as usize)
        .sum::<usize>();
    check_shape(rows.len(), count, caps)?;
    if count == basis.len() {
        return Ok(Matrix {
            columns: basis.to_vec(),
            rows,
        });
    }
    let mut columns = Vec::with_capacity(count);
    let mut remap = vec![u32::MAX; basis.len()];
    for (block, &word) in occupied.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let source = block * 64 + bits.trailing_zeros() as usize;
            remap[source] = columns.len() as u32;
            columns.push(basis[source]);
            bits &= bits - 1;
        }
    }
    let words = count.div_ceil(64);
    let rows = rows
        .into_iter()
        .map(|source| {
            let mut out = vec![0u64; words];
            for (block, mut bits) in source.into_iter().enumerate() {
                while bits != 0 {
                    let to = remap[block * 64 + bits.trailing_zeros() as usize] as usize;
                    out[to / 64] |= 1u64 << (to % 64);
                    bits &= bits - 1;
                }
            }
            out
        })
        .collect();
    Ok(Matrix { columns, rows })
}

// Independent oracle: enumerate combinations recursively, cancel in ordered
// sets, and fill the result column-by-column using monomial-to-row incidence.
fn reference_multipliers(active: u64, gap: u32) -> Vec<u64> {
    fn visit(bits: &[u32], start: usize, left: u32, mask: u64, out: &mut Vec<u64>) {
        out.push(mask);
        if left == 0 {
            return;
        }
        for i in start..bits.len() {
            visit(bits, i + 1, left - 1, mask | (1 << bits[i]), out);
        }
    }
    let bits: Vec<_> = (0..64).filter(|i| active & (1 << i) != 0).collect();
    let mut out = Vec::new();
    visit(&bits, 0, gap, 0, &mut out);
    out.sort_unstable();
    out
}
fn oracle(input: &Input, caps: Caps) -> Built {
    input.validate()?;
    let mut incidence = BTreeMap::<u64, Vec<usize>>::new();
    let mut row_count = 0;
    for poly in &input.polys {
        let Some(d) = poly.iter().map(|m| m.count_ones()).max() else {
            continue;
        };
        if d > u32::from(input.degree) {
            continue;
        }
        for t in reference_multipliers(input.active, u32::from(input.degree) - d) {
            let mut row = BTreeSet::new();
            for m in poly {
                if !row.insert(m | t) {
                    row.remove(&(m | t));
                }
            }
            if row.is_empty() {
                continue;
            }
            for m in row {
                incidence.entry(m).or_default().push(row_count);
            }
            row_count += 1;
        }
    }
    check_shape(row_count, incidence.len(), caps)?;
    let mut rows = vec![vec![0u64; incidence.len().div_ceil(64)]; row_count];
    let mut columns = Vec::new();
    for (column, (monomial, members)) in incidence.into_iter().enumerate() {
        columns.push(monomial);
        for row in members {
            rows[row][column / 64] |= 1 << (column % 64);
        }
    }
    Ok(Matrix { columns, rows })
}

// prefix[b][r] = sum_{j=0}^r C(b,j). For each set bit b from high to low,
// add all smaller prefixes with that bit zero, then spend one unit of degree.
// No table is indexed by a 2^n monomial mask.
struct CombinatorialRank {
    prefix: Vec<[usize; 4]>,
    degree: u8,
}
impl CombinatorialRank {
    fn new(n: u8, degree: u8) -> Self {
        let mut choose = vec![[0usize; 4]; n as usize + 1];
        for b in 0..=n as usize {
            choose[b][0] = 1;
            if b > 0 {
                for j in 1..=3 {
                    choose[b][j] = choose[b - 1][j] + choose[b - 1][j - 1];
                }
            }
        }
        for row in &mut choose {
            for j in 1..=3 {
                row[j] += row[j - 1];
            }
        }
        Self {
            prefix: choose,
            degree,
        }
    }
    fn index(&self, mut monomial: u64) -> usize {
        let mut remaining = self.degree as usize;
        let mut rank = 0;
        while monomial != 0 {
            let bit = 63 - monomial.leading_zeros() as usize;
            rank += self.prefix[bit][remaining];
            remaining -= 1;
            monomial ^= 1u64 << bit;
        }
        rank
    }
    fn retained(&self) -> usize {
        size_of::<Self>() + self.prefix.capacity() * size_of::<[usize; 4]>()
    }
}
enum Lookup {
    Binary,
    Ranked(CombinatorialRank),
    Dense(Vec<u32>),
}
struct Packed {
    n: u8,
    degree: u8,
    active: u64,
    columns: Vec<u64>,
    schedules: Vec<Vec<u64>>,
    lookup: Lookup,
}
impl Packed {
    fn compile(input: &Input, mode: &str) -> Result<Self, Error> {
        input.validate()?;
        if mode == "dense" && input.n > 12 {
            return Err(Error::DenseLimit);
        }
        let columns = multipliers((1 << input.n) - 1, u32::from(input.degree));
        let lookup = match mode {
            "binary" => Lookup::Binary,
            "ranked" => Lookup::Ranked(CombinatorialRank::new(input.n, input.degree)),
            "dense" => {
                let mut slots = vec![u32::MAX; 1usize << input.n];
                for (i, &m) in columns.iter().enumerate() {
                    slots[m as usize] = i as u32;
                }
                Lookup::Dense(slots)
            }
            _ => panic!("unknown lookup"),
        };
        Ok(Self {
            n: input.n,
            degree: input.degree,
            active: input.active,
            columns,
            schedules: (0..=input.degree)
                .map(|gap| multipliers(input.active, u32::from(gap)))
                .collect(),
            lookup,
        })
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        if let Err(error) = input.validate() {
            return (Err(error), false);
        }
        if (input.n, input.degree, input.active) != (self.n, self.degree, self.active) {
            return (direct(input, caps), false);
        }
        let words = self.columns.len().div_ceil(64);
        let mut seen = vec![0u64; words];
        let mut rows = Vec::new();
        for poly in &input.polys {
            let Some(d) = poly.iter().map(|m| m.count_ones() as u8).max() else {
                continue;
            };
            if d > input.degree {
                continue;
            }
            for &t in &self.schedules[(input.degree - d) as usize] {
                let mut row = vec![0u64; words];
                for &m in poly {
                    let product = m | t;
                    let slot = match &self.lookup {
                        Lookup::Binary => self.columns.binary_search(&product).unwrap(),
                        Lookup::Ranked(rank) => rank.index(product),
                        Lookup::Dense(slots) => slots[product as usize] as usize,
                    };
                    row[slot / 64] ^= 1u64 << (slot % 64);
                }
                if let Err(error) = retain_row(&mut rows, row, &mut seen, caps) {
                    return (Err(error), true);
                }
            }
        }
        (compact_packed(rows, &self.columns, &seen, caps), true)
    }
    fn retained(&self) -> usize {
        size_of::<Self>()
            + self.columns.capacity() * 8
            + self.schedules.capacity() * size_of::<Vec<u64>>()
            + self
                .schedules
                .iter()
                .map(|s| s.capacity() * 8)
                .sum::<usize>()
            + match &self.lookup {
                Lookup::Binary => 0,
                Lookup::Ranked(r) => r.retained() - size_of::<CombinatorialRank>(),
                Lookup::Dense(s) => s.capacity() * 4,
            }
    }
    fn lookup_entries(&self) -> usize {
        match &self.lookup {
            Lookup::Binary => 0,
            Lookup::Ranked(r) => r.prefix.len() * 4,
            Lookup::Dense(s) => s.len(),
        }
    }
}
enum Arm {
    Sorted,
    Packed(Packed),
}
impl Arm {
    fn compile(name: &str, input: &Input) -> Self {
        if name == "sorted" {
            Self::Sorted
        } else {
            Self::Packed(Packed::compile(input, name).unwrap())
        }
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        match self {
            Self::Sorted => (direct(input, caps), false),
            Self::Packed(p) => p.apply(input, caps),
        }
    }
    fn retained(&self) -> usize {
        match self {
            Self::Sorted => 0,
            Self::Packed(p) => p.retained(),
        }
    }
    fn lookup_entries(&self) -> usize {
        match self {
            Self::Sorted => 0,
            Self::Packed(p) => p.lookup_entries(),
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
fn fixtures(n: u8, seed: u64, batch: usize, family: &str) -> Vec<Input> {
    assert!(["quadratic", "linear_drop", "restricted_cycle"].contains(&family));
    let mut state = seed;
    let mut supports = Vec::new();
    for _ in 0..8 {
        let mut terms = BTreeSet::new();
        while terms.len() < 16 {
            let a = next(&mut state) % u64::from(n);
            let b = next(&mut state) % u64::from(n);
            if a != b {
                terms.insert((1 << a) | (1 << b));
            }
        }
        terms.insert(0);
        terms.extend((0..n.min(10)).map(|i| 1 << i));
        supports.push(terms.into_iter().collect::<Poly>());
    }
    let active = if family == "restricted_cycle" {
        (0..8).fold(0, |m, i| m | (1 << (i * n / 8)))
    } else {
        (1 << n) - 1
    };
    (0..batch)
        .map(|i| {
            let mut polys: Vec<Poly> = supports
                .iter()
                .map(|support| {
                    let mut poly: Poly = support
                        .iter()
                        .copied()
                        .filter(|_| next(&mut state) & 1 != 0)
                        .collect();
                    if !poly.iter().any(|m| m.count_ones() == 2) {
                        poly.push(*support.iter().find(|m| m.count_ones() == 2).unwrap());
                        poly.sort_unstable();
                    }
                    poly
                })
                .collect();
            if family != "quadratic" {
                match i % 4 {
                    1 => {
                        polys[0].retain(|m| m.count_ones() == 1);
                        polys[0].push(1);
                        polys[0].sort_unstable();
                        polys[0].dedup();
                    }
                    2 if family == "restricted_cycle" => polys[0] = vec![0],
                    3 => polys[0].clear(),
                    _ => (),
                }
            }
            Input {
                n,
                degree: 3,
                active,
                polys,
            }
        })
        .collect()
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 6, "worker N SEED BATCH FAMILY REPETITIONS");
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let batch = args[3].parse::<usize>().unwrap();
    let family = &args[4];
    let repetitions = args[5].parse::<usize>().unwrap();
    assert!(
        [12, 20, 28, 36].contains(&n)
            && [1, 8, 32].contains(&batch)
            && (1..=24).contains(&repetitions)
    );
    let inputs = fixtures(n, seed, batch, family);
    assert_eq!(
        inputs
            .iter()
            .map(|x| &x.polys)
            .collect::<BTreeSet<_>>()
            .len(),
        batch
    );
    let expected: Vec<_> = inputs.iter().map(|x| oracle(x, CAPS).unwrap()).collect();
    let polys: Vec<_> = inputs.iter().map(|x| &x.polys).collect();
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"batch\":{batch},\"family\":\"{family}\",\"degree\":3,\"active\":{},\"inputs\":{:?}}}", inputs[0].active, polys);
    let names: &[&str] = if n == 12 {
        &["sorted", "binary", "ranked", "dense"]
    } else {
        &["sorted", "binary", "ranked"]
    };
    for rep in 0..repetitions {
        for order in 0..names.len() {
            let variant = (rep + order) % names.len();
            let start = Instant::now();
            let arm = Arm::compile(names[variant], black_box(&inputs[0]));
            let setup_ns = start.elapsed().as_nanos();
            let retained_bytes = arm.retained();
            let lookup_entries = arm.lookup_entries();
            let (mut apply_ns, mut validation_ns) = (0u128, 0u128);
            let (mut hits, mut output_bytes, mut total_rows, mut total_columns) = (0, 0, 0, 0);
            for (input, expected) in inputs.iter().zip(&expected) {
                let tick = Instant::now();
                let (got, hit) = arm.apply(black_box(input), CAPS);
                let got = got.unwrap();
                apply_ns += tick.elapsed().as_nanos();
                hits += usize::from(hit);
                let tick = Instant::now();
                assert_eq!(black_box(&got), expected);
                output_bytes += got.payload();
                total_rows += got.rows.len();
                total_columns += got.columns.len();
                validation_ns += tick.elapsed().as_nanos();
            }
            drop(arm);
            let total_ns = start.elapsed().as_nanos();
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{}\",\"setup_ns\":{setup_ns},\"apply_ns\":{apply_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"retained_bytes\":{retained_bytes},\"lookup_entries\":{lookup_entries},\"output_bytes\":{output_bytes},\"total_rows\":{total_rows},\"total_columns\":{total_columns},\"hits\":{hits},\"verified_outputs\":{batch}}}", names[variant]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn small(poly: Poly, degree: u8, active: u64) -> Input {
        Input {
            n: 3,
            degree,
            active,
            polys: vec![poly],
        }
    }
    #[test]
    fn ranking_matches_numeric_enumeration_and_large_combinations() {
        for n in 1..=12 {
            for d in 0..=3 {
                let rank = CombinatorialRank::new(n, d);
                let monomials: Vec<_> = (0u64..(1 << n))
                    .filter(|m| m.count_ones() <= u32::from(d))
                    .collect();
                for (index, &m) in monomials.iter().enumerate() {
                    assert_eq!(rank.index(m), index, "n={n}, d={d}, m={m}");
                }
                assert_eq!(rank.prefix[n as usize][d as usize], monomials.len());
            }
        }
        for n in [20, 28, 36] {
            let rank = CombinatorialRank::new(n, 3);
            let monomials = reference_multipliers((1 << n) - 1, 3);
            assert_eq!(rank.prefix[n as usize][3], monomials.len());
            for (i, m) in monomials.into_iter().enumerate() {
                assert_eq!(rank.index(m), i);
            }
        }
    }
    #[test]
    fn exhaustive_coefficients_degrees_masks_and_all_constructors() {
        for d in 1..=3 {
            for active in 0..8 {
                let context = small(vec![], d, active);
                let arms: Vec<_> = ["binary", "ranked", "dense"]
                    .iter()
                    .map(|name| Packed::compile(&context, name).unwrap())
                    .collect();
                for selected in 0..256u16 {
                    let input = small(
                        (0..8).filter(|m| selected & (1 << m) != 0).collect(),
                        d,
                        active,
                    );
                    let expected = oracle(&input, CAPS);
                    assert_eq!(direct(&input, CAPS), expected);
                    for arm in &arms {
                        assert_eq!(arm.apply(&input, CAPS), (expected.clone(), true));
                    }
                }
            }
        }
    }
    #[test]
    fn fixture_families_and_high_variable_bits_match() {
        for n in [12, 20, 28, 36] {
            for family in ["quadratic", "linear_drop", "restricted_cycle"] {
                let inputs = fixtures(n, 31337, 4, family);
                let arms: Vec<_> = ["binary", "ranked"]
                    .iter()
                    .map(|name| Packed::compile(&inputs[0], name).unwrap())
                    .collect();
                assert!(inputs
                    .iter()
                    .flat_map(|x| &x.polys)
                    .flatten()
                    .any(|m| m & (1 << (n - 1)) != 0));
                for input in inputs {
                    let expected = oracle(&input, CAPS);
                    assert_eq!(direct(&input, CAPS), expected);
                    for arm in &arms {
                        assert_eq!(arm.apply(&input, CAPS), (expected.clone(), true));
                    }
                }
            }
        }
    }
    #[test]
    fn current_caps_and_empty_actual_support() {
        let base = small(vec![0, 1, 2, 3], 3, 7);
        let inputs = [
            base.clone(),
            small(vec![], 3, 7),
            small(vec![0], 3, 7),
            small(vec![0, 1], 3, 7),
        ];
        for name in ["sorted", "binary", "ranked", "dense"] {
            let arm = Arm::compile(name, &base);
            for input in &inputs {
                for caps in [
                    CAPS,
                    Caps {
                        rows: 0,
                        columns: 0,
                    },
                    Caps {
                        rows: 1,
                        columns: 1,
                    },
                    Caps {
                        rows: 4096,
                        columns: 0,
                    },
                ] {
                    assert_eq!(arm.apply(input, caps).0, oracle(input, caps));
                }
            }
        }
        let narrow = Input {
            n: 36,
            degree: 3,
            active: 0,
            polys: vec![vec![1u64 << 35]],
        };
        let caps = Caps {
            rows: 1,
            columns: 1,
        };
        assert_eq!(
            Packed::compile(&narrow, "ranked")
                .unwrap()
                .apply(&narrow, caps)
                .0,
            oracle(&narrow, caps)
        );
    }
    #[test]
    fn changed_context_falls_back_and_invalid_input_is_rejected() {
        let base = small(vec![0, 1, 3], 3, 7);
        for name in ["binary", "ranked", "dense"] {
            let arm = Packed::compile(&base, name).unwrap();
            for input in [
                Input {
                    n: 4,
                    ..base.clone()
                },
                small(vec![0, 1, 3], 2, 7),
                small(vec![0, 1, 3], 3, 3),
            ] {
                assert_eq!(arm.apply(&input, CAPS), (oracle(&input, CAPS), false));
            }
            for input in [
                small(vec![1, 1], 3, 7),
                small(vec![8], 3, 7),
                small(vec![0], 3, 8),
            ] {
                assert_eq!(arm.apply(&input, CAPS), (Err(Error::InvalidInput), false));
            }
        }
    }
    #[test]
    fn dense_lookup_refuses_large_dimensions_before_allocation() {
        let input = Input {
            n: 36,
            degree: 3,
            active: 0,
            polys: vec![],
        };
        assert!(matches!(
            Packed::compile(&input, "dense"),
            Err(Error::DenseLimit)
        ));
        let rank = Packed::compile(&input, "ranked").unwrap();
        assert_eq!(rank.lookup_entries(), 37 * 4);
        assert_eq!(rank.columns.len(), 7807);
    }
    #[test]
    fn full_width_constant_can_be_resource_censored() {
        let input = Input {
            n: 36,
            degree: 3,
            active: (1 << 36) - 1,
            polys: vec![vec![0]],
        };
        assert_eq!(
            Packed::compile(&input, "ranked")
                .unwrap()
                .apply(&input, CAPS)
                .0,
            Err(Error::RowCap)
        );
        assert_eq!(direct(&input, CAPS), Err(Error::RowCap));
    }
}
