//! Standalone, bounded Boolean product-schedule experiment. Uses only std.
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
    rows: 256,
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
        let columns = columns(&products, caps)?;
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
enum Cache {
    Direct,
    Layout(Vec<u16>),
    Schedule(Schedule),
    Matrix(Input, Matrix),
}
impl Cache {
    fn compile(variant: usize, input: &Input) -> Self {
        match variant {
            0 => Self::Direct,
            1 => Self::Layout(columns(&products(input, CAPS).unwrap(), CAPS).unwrap()),
            2 => Self::Schedule(Schedule::compile(input, CAPS).unwrap()),
            3 => Self::Matrix(input.clone(), direct(input, CAPS).unwrap()),
            _ => unreachable!(),
        }
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        match self {
            Self::Direct => (direct(input, caps), false),
            Self::Schedule(schedule) => schedule.apply(input, caps),
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
fn changed(input: &Input, family: &str) -> Input {
    let mut out = input.clone();
    match family {
        "repeat" => (),
        "constant_toggle" => out.polys[0].insert(0, 0), // generated fixtures have no constant
        "term_toggle" => {
            out.polys[0].remove(0);
        }
        "degree_drop" => out.polys[0] = vec![1, 2],
        _ => panic!("unknown family"),
    }
    out
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
    let alternate = changed(&base, family);
    let inputs: Vec<_> = (0..batch)
        .map(|i| {
            if i % 2 == 0 {
                base.clone()
            } else {
                alternate.clone()
            }
        })
        .collect();
    let expected: Vec<_> = inputs.iter().map(|x| oracle(x, CAPS).unwrap()).collect();
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"batch\":{batch},\"family\":\"{family}\",\"degree\":3,\"active\":{},\"base\":{:?},\"alternate\":{:?}}}", base.active, base.polys, alternate.polys);
    let names = ["direct", "layout", "schedule", "matrix_cache"];
    for rep in 0..repetitions {
        for order in 0..4 {
            let variant = (rep + order) % 4;
            let start = Instant::now();
            let cache = Cache::compile(variant, black_box(&base));
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
    #[test]
    fn exhaustive_small_supports_and_evaluations() {
        for n in 1..=3u8 {
            for selected in 0..(1u16 << (1 << n)) {
                let poly: Vec<_> = (0u16..(1 << n))
                    .filter(|m| selected & (1 << m) != 0)
                    .collect();
                for degree in 1..=3 {
                    let input = Input {
                        n,
                        degree,
                        active: (1 << n) - 1,
                        polys: vec![poly.clone()],
                    };
                    let expected = oracle(&input, CAPS).unwrap();
                    assert_eq!(direct(&input, CAPS).unwrap(), expected);
                    let plan = Schedule::compile(&input, CAPS).unwrap();
                    assert_eq!(plan.apply(&input, CAPS), (Ok(expected.clone()), true));
                    // Every output row must vanish wherever its generator vanishes.
                    for point in 0u16..(1 << n) {
                        let eval = |m: u16| u64::from(point & m == m);
                        let p = poly.iter().fold(0, |s, &m| s ^ eval(m));
                        if p == 0 {
                            for row in &expected.rows {
                                let value =
                                    expected.columns.iter().enumerate().fold(0, |s, (i, &m)| {
                                        s ^ (((row[i / 64] >> (i % 64)) & 1) & eval(m))
                                    });
                                assert_eq!(value, 0);
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn changed_supports_and_degrees_fall_back() {
        let base = fixture(6, 17);
        for family in ["constant_toggle", "term_toggle", "degree_drop"] {
            let input = changed(&base, family);
            for variant in [2, 3] {
                let cache = Cache::compile(variant, &base);
                let (got, hit) = cache.apply(&input, CAPS);
                assert!(!hit);
                assert_eq!(got, oracle(&input, CAPS));
            }
        }
        let mut variants = Vec::new();
        let mut input = base.clone();
        input.polys.swap(0, 1);
        variants.push(input);
        let mut input = base.clone();
        input.degree = 2;
        variants.push(input);
        let mut input = base.clone();
        input.active >>= 1;
        variants.push(input);
        let plan = Schedule::compile(&base, CAPS).unwrap();
        for input in variants {
            assert!(!plan.apply(&input, CAPS).1);
            assert_eq!(plan.apply(&input, CAPS).0, oracle(&input, CAPS));
        }
    }
    #[test]
    fn all_arms_honor_lowered_caps_and_smaller_fallback() {
        let base = Input {
            n: 2,
            degree: 2,
            active: 3,
            polys: vec![vec![0, 1, 2, 3]],
        };
        let small = Input {
            polys: vec![vec![1, 2, 3]],
            ..base.clone()
        };
        for variant in 0..4 {
            let cache = Cache::compile(variant, &base);
            for caps in [
                Caps {
                    rows: 0,
                    columns: 4,
                },
                Caps {
                    rows: 1,
                    columns: 3,
                },
            ] {
                assert_eq!(cache.apply(&base, caps).0, oracle(&base, caps));
                assert_eq!(cache.apply(&small, caps).0, oracle(&small, caps));
            }
            assert_eq!(cache.apply(&base, CAPS).0, oracle(&base, CAPS));
        }
    }
    #[test]
    fn exact_support_means_identical_matrix() {
        for n in [6, 8, 10, 12] {
            let base = fixture(n, 937);
            assert_eq!(
                Schedule::compile(&base, CAPS)
                    .unwrap()
                    .apply(&base.clone(), CAPS)
                    .0,
                direct(&base, CAPS)
            );
        }
    }
}
