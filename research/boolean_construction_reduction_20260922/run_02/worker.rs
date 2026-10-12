//! Standalone, bounded Boolean construction and canonical linear reduction experiment. Uses only std.
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

fn compact_sparse(
    rows: Vec<Vec<(usize, u64)>>,
    basis: &[u64],
    occupied: &[u64],
    caps: Caps,
) -> Built {
    let count = occupied
        .iter()
        .map(|w| w.count_ones() as usize)
        .sum::<usize>();
    check_shape(rows.len(), count, caps)?;
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
    let identity = count == basis.len();
    let rows = rows
        .into_iter()
        .map(|source| {
            let mut out = vec![0u64; words];
            for (block, mut bits) in source {
                if identity {
                    out[block] = bits;
                    continue;
                }
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
    sparse: bool,
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
            "ranked" | "sparse_rank" => {
                Lookup::Ranked(CombinatorialRank::new(input.n, input.degree))
            }
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
            sparse: mode == "sparse_rank",
        })
    }
    fn apply(&self, input: &Input, caps: Caps) -> (Built, bool) {
        if let Err(error) = input.validate() {
            return (Err(error), false);
        }
        if (input.n, input.degree, input.active) != (self.n, self.degree, self.active) {
            return (direct(input, caps), false);
        }
        if self.sparse {
            return (self.apply_sparse(input, caps), true);
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
    fn apply_sparse(&self, input: &Input, caps: Caps) -> Built {
        let Lookup::Ranked(rank) = &self.lookup else {
            unreachable!()
        };
        let words = self.columns.len().div_ceil(64);
        let mut scratch = vec![0u64; words];
        let mut marked = vec![false; words];
        let mut touched = Vec::with_capacity(32);
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
                for &m in poly {
                    let slot = rank.index(m | t);
                    let word = slot / 64;
                    if !marked[word] {
                        marked[word] = true;
                        touched.push(word);
                    }
                    scratch[word] ^= 1u64 << (slot % 64);
                }
                let mut row = Vec::with_capacity(touched.len());
                for word in touched.drain(..) {
                    let value = scratch[word];
                    scratch[word] = 0;
                    marked[word] = false;
                    if value != 0 {
                        row.push((word, value));
                        seen[word] |= value;
                    }
                }
                if !row.is_empty() {
                    if rows.len() == caps.rows {
                        return Err(Error::RowCap);
                    }
                    rows.push(row);
                }
            }
        }
        compact_sparse(rows, &self.columns, &seen, caps)
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

#[derive(Clone, Copy, Debug, Default)]
struct ReduceStats {
    row_xors: u64,
    word_xors: u64,
}
struct Basis {
    by_column: Vec<u32>,
    pivots: Vec<usize>,
    rows: Vec<Vec<u64>>,
    stats: ReduceStats,
}
impl Basis {
    fn new(columns: usize) -> Self {
        Self {
            by_column: vec![u32::MAX; columns],
            pivots: vec![],
            rows: vec![],
            stats: ReduceStats::default(),
        }
    }
    fn insert(&mut self, mut row: Vec<u64>) {
        let mut start = 0;
        loop {
            let Some(offset) = row[start..].iter().position(|&w| w != 0) else {
                return;
            };
            start += offset;
            let pivot = start * 64 + row[start].trailing_zeros() as usize;
            let index = self.by_column[pivot];
            if index == u32::MAX {
                self.by_column[pivot] = self.rows.len() as u32;
                self.pivots.push(pivot);
                self.rows.push(row);
                return;
            }
            let base = &self.rows[index as usize];
            self.stats.row_xors += 1;
            self.stats.word_xors += (row.len() - start) as u64;
            for (to, &from) in row[start..].iter_mut().zip(&base[start..]) {
                *to ^= from;
            }
        }
    }
    fn storage_bytes(&self) -> usize {
        size_of::<Self>()
            + self.by_column.capacity() * 4
            + self.pivots.capacity() * size_of::<usize>()
            + self.rows.capacity() * size_of::<Vec<u64>>()
            + self.rows.iter().map(|r| r.capacity() * 8).sum::<usize>()
    }
    fn finish(self) -> (Vec<Vec<u64>>, ReduceStats, usize) {
        let storage = self.storage_bytes();
        let mut stats = self.stats;
        let mut ordered: Vec<_> = self.pivots.into_iter().zip(self.rows).collect();
        ordered.sort_unstable_by_key(|(pivot, _)| *pivot);
        for i in (0..ordered.len()).rev() {
            let (above, tail) = ordered.split_at_mut(i);
            let (pivot, row) = &tail[0];
            let word = pivot / 64;
            for (_, target) in above {
                if target[word] & (1u64 << (pivot % 64)) != 0 {
                    stats.row_xors += 1;
                    stats.word_xors += (row.len() - word) as u64;
                    for (to, &from) in target[word..].iter_mut().zip(&row[word..]) {
                        *to ^= from;
                    }
                }
            }
        }
        (
            ordered.into_iter().map(|(_, row)| row).collect(),
            stats,
            storage,
        )
    }
}
struct Answer {
    matrix: Matrix,
    source_rows: usize,
    stats: ReduceStats,
    basis_bytes: usize,
    auxiliary_bytes: usize,
}
fn reduce(matrix: Matrix) -> Answer {
    let source_rows = matrix.rows.len();
    let mut basis = Basis::new(matrix.columns.len());
    for row in matrix.rows {
        basis.insert(row);
    }
    let (rows, stats, basis_bytes) = basis.finish();
    Answer {
        matrix: Matrix {
            columns: matrix.columns,
            rows,
        },
        source_rows,
        stats,
        basis_bytes,
        auxiliary_bytes: 0,
    }
}

// Independent column-oriented Gauss-Jordan oracle. It searches matrix rows by
// column and clears above and below immediately, unlike pivot-table insertion.
fn oracle_reduce(mut matrix: Matrix) -> Matrix {
    let mut rank = 0;
    for column in 0..matrix.columns.len() {
        let Some(found) = (rank..matrix.rows.len())
            .find(|&r| matrix.rows[r][column / 64] & (1u64 << (column % 64)) != 0)
        else {
            continue;
        };
        matrix.rows.swap(rank, found);
        let pivot = matrix.rows[rank].clone();
        for (i, row) in matrix.rows.iter_mut().enumerate() {
            if i != rank && row[column / 64] & (1u64 << (column % 64)) != 0 {
                for (to, &from) in row.iter_mut().zip(&pivot) {
                    *to ^= from;
                }
            }
        }
        rank += 1;
        if rank == matrix.rows.len() {
            break;
        }
    }
    matrix.rows.truncate(rank);
    matrix
}
fn is_rref(matrix: &Matrix) -> bool {
    let mut previous = None;
    let words = matrix.columns.len().div_ceil(64);
    for (i, row) in matrix.rows.iter().enumerate() {
        if row.len() != words {
            return false;
        }
        let Some(word) = row.iter().position(|&w| w != 0) else {
            return false;
        };
        let pivot = word * 64 + row[word].trailing_zeros() as usize;
        if pivot >= matrix.columns.len() || previous.is_some_and(|p| p >= pivot) {
            return false;
        }
        if matrix
            .rows
            .iter()
            .enumerate()
            .any(|(j, other)| j != i && other[word] & (1u64 << (pivot % 64)) != 0)
        {
            return false;
        }
        previous = Some(pivot);
    }
    true
}
fn stream(context: &Packed, input: &Input, caps: Caps) -> Result<Answer, Error> {
    input.validate()?;
    if (input.n, input.degree, input.active) != (context.n, context.degree, context.active) {
        return direct(input, caps).map(reduce);
    }
    let Lookup::Ranked(rank) = &context.lookup else {
        unreachable!()
    };
    let words = context.columns.len().div_ceil(64);
    let mut seen = vec![0u64; words];
    let mut basis = Basis::new(context.columns.len());
    let mut source_rows = 0;
    for poly in &input.polys {
        let Some(d) = poly.iter().map(|m| m.count_ones() as u8).max() else {
            continue;
        };
        if d > input.degree {
            continue;
        }
        for &t in &context.schedules[(input.degree - d) as usize] {
            let mut row = vec![0u64; words];
            for &m in poly {
                let slot = rank.index(m | t);
                row[slot / 64] ^= 1u64 << (slot % 64);
            }
            if row.iter().all(|&w| w == 0) {
                continue;
            }
            if source_rows == caps.rows {
                return Err(Error::RowCap);
            }
            source_rows += 1;
            for (to, &from) in seen.iter_mut().zip(&row) {
                *to |= from;
            }
            basis.insert(row);
        }
    }
    let actual_columns = seen.iter().map(|w| w.count_ones() as usize).sum();
    check_shape(source_rows, actual_columns, caps)?;
    let (rows, stats, basis_bytes) = basis.finish();
    let matrix = compact_packed(rows, &context.columns, &seen, caps)?;
    Ok(Answer {
        matrix,
        source_rows,
        stats,
        basis_bytes,
        auxiliary_bytes: 0,
    })
}

struct WordRow {
    pivot: usize,
    words: Vec<u64>,
    nonzero: Vec<usize>,
}
fn nonzero_words(row: &[u64]) -> Vec<usize> {
    row.iter()
        .enumerate()
        .filter(|(_, w)| **w != 0)
        .map(|(i, _)| i)
        .collect()
}
fn reduce_incidence(matrix: Matrix) -> Answer {
    let source_rows = matrix.rows.len();
    let mut by_column = vec![u32::MAX; matrix.columns.len()];
    let mut rows = Vec::<WordRow>::new();
    let mut stats = ReduceStats::default();
    for mut row in matrix.rows {
        let mut start = 0;
        loop {
            let Some(offset) = row[start..].iter().position(|&w| w != 0) else {
                break;
            };
            start += offset;
            let pivot = start * 64 + row[start].trailing_zeros() as usize;
            let index = by_column[pivot];
            if index == u32::MAX {
                by_column[pivot] = rows.len() as u32;
                let nonzero = nonzero_words(&row);
                rows.push(WordRow {
                    pivot,
                    words: row,
                    nonzero,
                });
                break;
            }
            let base = &rows[index as usize];
            stats.row_xors += 1;
            stats.word_xors += base.nonzero.len() as u64;
            for &word in &base.nonzero {
                row[word] ^= base.words[word];
            }
        }
    }
    let basis_bytes = size_of::<Vec<u32>>()
        + size_of::<Vec<WordRow>>()
        + size_of::<ReduceStats>()
        + by_column.capacity() * 4
        + rows.capacity() * size_of::<WordRow>()
        + rows
            .iter()
            .map(|r| r.words.capacity() * 8 + r.nonzero.capacity() * size_of::<usize>())
            .sum::<usize>();
    rows.sort_unstable_by_key(|r| r.pivot);
    let mut pivot_mask = vec![0u64; matrix.columns.len().div_ceil(64)];
    for (i, row) in rows.iter().enumerate() {
        by_column[row.pivot] = i as u32;
        pivot_mask[row.pivot / 64] |= 1u64 << (row.pivot % 64);
    }
    // Higher-pivot elimination cannot alter a lower-pivot column. Therefore
    // incidence in the forward-echelon rows is valid for the entire backward pass.
    let mut targets = vec![Vec::<usize>::new(); rows.len()];
    for (i, row) in rows.iter().enumerate() {
        for &word in &row.nonzero {
            let mut bits = row.words[word] & pivot_mask[word];
            while bits != 0 {
                let column = word * 64 + bits.trailing_zeros() as usize;
                let j = by_column[column] as usize;
                debug_assert!(j >= i);
                if j > i {
                    targets[j].push(i);
                }
                bits &= bits - 1;
            }
        }
    }
    let auxiliary_bytes = size_of::<Vec<u64>>()
        + size_of::<Vec<Vec<usize>>>()
        + pivot_mask.capacity() * 8
        + targets.capacity() * size_of::<Vec<usize>>()
        + targets
            .iter()
            .map(|t| t.capacity() * size_of::<usize>())
            .sum::<usize>();
    for j in (0..rows.len()).rev() {
        // This row has now received all higher-pivot eliminations. Refresh its
        // sparse word list once before using it as a backward pivot.
        let (above, tail) = rows.split_at_mut(j);
        let base = &mut tail[0];
        base.nonzero.clear();
        base.nonzero.extend(
            base.words
                .iter()
                .enumerate()
                .filter(|(_, w)| **w != 0)
                .map(|(w, _)| w),
        );
        for &i in &targets[j] {
            let to = &mut above[i].words;
            stats.row_xors += 1;
            stats.word_xors += base.nonzero.len() as u64;
            for &word in &base.nonzero {
                to[word] ^= base.words[word];
            }
        }
    }
    Answer {
        matrix: Matrix {
            columns: matrix.columns,
            rows: rows.into_iter().map(|r| r.words).collect(),
        },
        source_rows,
        stats,
        basis_bytes,
        auxiliary_bytes,
    }
}
enum Pipeline {
    Sorted,
    Staged(Packed),
    Stream(Packed),
    Incidence(Packed),
}
impl Pipeline {
    fn compile(name: &str, input: &Input) -> Self {
        match name {
            "sorted_reduce" => Self::Sorted,
            "ranked_reduce" => Self::Staged(Packed::compile(input, "ranked").unwrap()),
            "sparse_reduce" => Self::Staged(Packed::compile(input, "sparse_rank").unwrap()),
            "dense_reduce" => Self::Staged(Packed::compile(input, "dense").unwrap()),
            "stream_reduce" => Self::Stream(Packed::compile(input, "ranked").unwrap()),
            "incidence_reduce" => Self::Incidence(Packed::compile(input, "sparse_rank").unwrap()),
            _ => panic!("unknown pipeline"),
        }
    }
    fn retained(&self) -> usize {
        match self {
            Self::Sorted => 0,
            Self::Staged(c) | Self::Stream(c) | Self::Incidence(c) => c.retained(),
        }
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        [6, 7].contains(&args.len()),
        "worker N SEED BATCH FAMILY REPETITIONS [with-incidence]"
    );
    let include_incidence = args.len() == 7;
    if include_incidence {
        assert_eq!(args[6], "with-incidence");
    }
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let batch = args[3].parse::<usize>().unwrap();
    let family = &args[4];
    let repetitions = args[5].parse::<usize>().unwrap();
    assert!(
        [12, 20, 28, 36].contains(&n)
            && [1, 4, 8].contains(&batch)
            && (1..=30).contains(&repetitions)
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
    let expected: Vec<_> = inputs
        .iter()
        .map(|x| {
            let matrix = oracle(x, CAPS).unwrap();
            let source_rows = matrix.rows.len();
            let matrix = oracle_reduce(matrix);
            assert!(is_rref(&matrix));
            (matrix, source_rows)
        })
        .collect();
    let polys: Vec<_> = inputs.iter().map(|x| &x.polys).collect();
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"batch\":{batch},\"family\":\"{family}\",\"degree\":3,\"active\":{},\"inputs\":{:?}}}",inputs[0].active,polys);
    let names: &[&str] = match (n == 12, include_incidence) {
        (true, false) => &[
            "sorted_reduce",
            "ranked_reduce",
            "sparse_reduce",
            "stream_reduce",
            "dense_reduce",
        ],
        (false, false) => &[
            "sorted_reduce",
            "ranked_reduce",
            "sparse_reduce",
            "stream_reduce",
        ],
        (true, true) => &[
            "sorted_reduce",
            "ranked_reduce",
            "sparse_reduce",
            "stream_reduce",
            "incidence_reduce",
            "dense_reduce",
        ],
        (false, true) => &[
            "sorted_reduce",
            "ranked_reduce",
            "sparse_reduce",
            "stream_reduce",
            "incidence_reduce",
        ],
    };
    for rep in 0..repetitions {
        for order in 0..names.len() {
            let name = names[(rep + order) % names.len()];
            let start = Instant::now();
            let pipeline = Pipeline::compile(name, black_box(&inputs[0]));
            let setup_ns = start.elapsed().as_nanos();
            let retained_bytes = pipeline.retained();
            let (mut construct_ns, mut reduce_ns, mut fused_ns, mut validation_ns) =
                (0u128, 0u128, 0u128, 0u128);
            let (
                mut output_bytes,
                mut source_rows,
                mut total_rank,
                mut total_columns,
                mut basis_max_bytes,
            ) = (0, 0, 0, 0, 0);
            let (mut row_xors, mut word_xors) = (0u64, 0u64);
            let mut auxiliary_max_bytes = 0;
            for (input, (expected, expected_source_rows)) in inputs.iter().zip(&expected) {
                let answer = match &pipeline {
                    Pipeline::Stream(context) => {
                        let tick = Instant::now();
                        let out = stream(context, black_box(input), CAPS).unwrap();
                        fused_ns += tick.elapsed().as_nanos();
                        out
                    }
                    _ => {
                        let tick = Instant::now();
                        let matrix = match &pipeline {
                            Pipeline::Sorted => direct(black_box(input), CAPS),
                            Pipeline::Staged(c) | Pipeline::Incidence(c) => {
                                c.apply(black_box(input), CAPS).0
                            }
                            _ => unreachable!(),
                        }
                        .unwrap();
                        construct_ns += tick.elapsed().as_nanos();
                        let tick = Instant::now();
                        let out = if matches!(&pipeline, Pipeline::Incidence(_)) {
                            reduce_incidence(matrix)
                        } else {
                            reduce(matrix)
                        };
                        reduce_ns += tick.elapsed().as_nanos();
                        out
                    }
                };
                let tick = Instant::now();
                assert_eq!(black_box(&answer.matrix), expected);
                assert_eq!(answer.source_rows, *expected_source_rows);
                output_bytes += answer.matrix.payload();
                source_rows += answer.source_rows;
                total_rank += answer.matrix.rows.len();
                total_columns += answer.matrix.columns.len();
                basis_max_bytes = basis_max_bytes.max(answer.basis_bytes);
                auxiliary_max_bytes = auxiliary_max_bytes.max(answer.auxiliary_bytes);
                row_xors += answer.stats.row_xors;
                word_xors += answer.stats.word_xors;
                validation_ns += tick.elapsed().as_nanos();
            }
            drop(pipeline);
            let total_ns = start.elapsed().as_nanos();
            let (construction, reduction, fused) = if name == "stream_reduce" {
                ("null".to_owned(), "null".to_owned(), fused_ns.to_string())
            } else {
                (
                    construct_ns.to_string(),
                    reduce_ns.to_string(),
                    "null".to_owned(),
                )
            };
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{name}\",\"setup_ns\":{setup_ns},\"construction_ns\":{construction},\"reduction_ns\":{reduction},\"fused_ns\":{fused},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"retained_bytes\":{retained_bytes},\"basis_max_bytes\":{basis_max_bytes},\"auxiliary_max_bytes\":{auxiliary_max_bytes},\"row_xors\":{row_xors},\"word_xors\":{word_xors},\"output_bytes\":{output_bytes},\"source_rows\":{source_rows},\"total_rank\":{total_rank},\"total_columns\":{total_columns},\"verified_outputs\":{batch}}}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn span(rows: &[Vec<u64>]) -> BTreeSet<Vec<u64>> {
        (0..(1usize << rows.len()))
            .map(|selected| {
                let mut row = vec![0u64; rows.first().map_or(0, Vec::len)];
                for (i, other) in rows.iter().enumerate() {
                    if selected & (1 << i) != 0 {
                        for (a, b) in row.iter_mut().zip(other) {
                            *a ^= b;
                        }
                    }
                }
                row
            })
            .collect()
    }
    #[test]
    fn exhaustive_small_matrices_match_oracle_and_row_span() {
        for cols in 1..=4 {
            for count in 1..=3 {
                for code in 0..(1u64 << (cols * count)) {
                    let rows = (0..count)
                        .map(|i| vec![(code >> (i * cols)) & ((1 << cols) - 1)])
                        .collect::<Vec<_>>();
                    let matrix = Matrix {
                        columns: (0..cols).collect(),
                        rows,
                    };
                    let got = reduce(matrix.clone());
                    let optimized = reduce_incidence(matrix.clone());
                    assert_eq!(optimized.matrix, got.matrix);
                    assert_eq!(optimized.stats.row_xors, got.stats.row_xors);
                    assert_eq!(got.matrix, oracle_reduce(matrix.clone()));
                    assert!(is_rref(&got.matrix));
                    if got.matrix.rows.is_empty() {
                        assert!(matrix.rows.iter().all(|r| r[0] == 0));
                    } else {
                        assert_eq!(span(&matrix.rows), span(&got.matrix.rows));
                    }
                    assert_eq!(got.source_rows, count as usize);
                }
            }
        }
    }
    #[test]
    fn exhaustive_boolean_coefficients_all_construction_paths() {
        for d in 1..=3 {
            for active in 0..8 {
                let context = Input {
                    n: 3,
                    degree: d,
                    active,
                    polys: vec![],
                };
                let constructors: Vec<_> = ["ranked", "sparse_rank", "dense"]
                    .iter()
                    .map(|mode| Packed::compile(&context, mode).unwrap())
                    .collect();
                for selected in 0..256u16 {
                    let input = Input {
                        polys: vec![(0..8).filter(|m| selected & (1 << m) != 0).collect()],
                        ..context.clone()
                    };
                    let original = oracle(&input, CAPS).unwrap();
                    let source_rows = original.rows.len();
                    let expected = oracle_reduce(original);
                    let baseline = reduce(direct(&input, CAPS).unwrap());
                    let optimized = reduce_incidence(direct(&input, CAPS).unwrap());
                    assert_eq!(optimized.matrix, baseline.matrix);
                    assert_eq!(optimized.stats.row_xors, baseline.stats.row_xors);
                    assert_eq!(baseline.matrix, expected);
                    for constructor in &constructors {
                        let answer = reduce(constructor.apply(&input, CAPS).0.unwrap());
                        assert_eq!(answer.matrix, expected);
                        assert_eq!(answer.source_rows, source_rows);
                        assert_eq!(answer.stats.row_xors, baseline.stats.row_xors);
                    }
                    let streamed = stream(&constructors[0], &input, CAPS).unwrap();
                    assert_eq!(streamed.matrix, expected);
                    assert_eq!(streamed.source_rows, source_rows);
                    assert_eq!(streamed.stats.row_xors, baseline.stats.row_xors);
                }
            }
        }
    }
    #[test]
    fn larger_fixtures_match_canonical_oracle_and_logical_eliminations() {
        for n in [12, 20, 28, 36] {
            for family in ["quadratic", "linear_drop", "restricted_cycle"] {
                let inputs = fixtures(n, 31337, 4, family);
                let dense = Packed::compile(&inputs[0], "ranked").unwrap();
                let sparse = Packed::compile(&inputs[0], "sparse_rank").unwrap();
                for input in inputs {
                    let source = oracle(&input, CAPS).unwrap();
                    let source_rows = source.rows.len();
                    let expected = oracle_reduce(source);
                    let baseline = reduce(direct(&input, CAPS).unwrap());
                    let optimized = reduce_incidence(direct(&input, CAPS).unwrap());
                    assert_eq!(optimized.matrix, baseline.matrix);
                    assert_eq!(optimized.stats.row_xors, baseline.stats.row_xors);
                    assert_eq!(baseline.matrix, expected);
                    assert!(is_rref(&expected));
                    for answer in [
                        reduce(dense.apply(&input, CAPS).0.unwrap()),
                        reduce(sparse.apply(&input, CAPS).0.unwrap()),
                        stream(&dense, &input, CAPS).unwrap(),
                    ] {
                        assert_eq!(answer.matrix, expected);
                        assert_eq!(answer.source_rows, source_rows);
                        assert_eq!(answer.stats.row_xors, baseline.stats.row_xors);
                    }
                }
            }
        }
    }
    #[test]
    fn dependent_source_rows_still_consume_the_row_budget() {
        let input = Input {
            n: 3,
            degree: 1,
            active: 7,
            polys: vec![vec![0], vec![0], vec![0]],
        };
        let context = Packed::compile(&input, "ranked").unwrap();
        let caps = Caps {
            rows: 5,
            columns: 8,
        };
        assert_eq!(direct(&input, caps), Err(Error::RowCap));
        assert!(matches!(stream(&context, &input, caps), Err(Error::RowCap)));
        let answer = stream(&context, &input, CAPS).unwrap();
        assert_eq!(answer.source_rows, 12);
        assert_eq!(answer.matrix.rows.len(), 4);
    }
    #[test]
    fn current_column_caps_and_narrow_outputs_in_wide_contexts() {
        let input = Input {
            n: 36,
            degree: 3,
            active: 0,
            polys: vec![vec![1u64 << 35]],
        };
        let context = Packed::compile(&input, "ranked").unwrap();
        let one = Caps {
            rows: 1,
            columns: 1,
        };
        assert_eq!(
            stream(&context, &input, one).unwrap().matrix,
            oracle_reduce(oracle(&input, one).unwrap())
        );
        assert!(matches!(
            stream(
                &context,
                &input,
                Caps {
                    rows: 1,
                    columns: 0
                }
            ),
            Err(Error::ColumnCap)
        ));
        let empty = Input {
            polys: vec![],
            ..input
        };
        let zero = Caps {
            rows: 0,
            columns: 0,
        };
        assert_eq!(
            stream(&context, &empty, zero).unwrap().matrix,
            Matrix {
                columns: vec![],
                rows: vec![]
            }
        );
    }
    #[test]
    fn changed_context_and_invalid_input_are_checked() {
        let input = Input {
            n: 3,
            degree: 3,
            active: 7,
            polys: vec![vec![0, 1, 2]],
        };
        let context = Packed::compile(&input, "ranked").unwrap();
        for changed in [
            Input {
                n: 4,
                ..input.clone()
            },
            Input {
                degree: 2,
                ..input.clone()
            },
            Input {
                active: 3,
                ..input.clone()
            },
        ] {
            assert_eq!(
                stream(&context, &changed, CAPS).unwrap().matrix,
                oracle_reduce(oracle(&changed, CAPS).unwrap())
            );
        }
        let invalid = Input {
            polys: vec![vec![1, 1]],
            ..input
        };
        assert!(matches!(
            stream(&context, &invalid, CAPS),
            Err(Error::InvalidInput)
        ));
    }
    #[test]
    fn pivot_order_and_word_boundaries_are_canonical() {
        let matrix = Matrix {
            columns: (0..70).collect(),
            rows: vec![vec![0, 2], vec![1, 2], vec![1u64 << 63, 1], vec![0, 2]],
        };
        let expected = oracle_reduce(matrix.clone());
        assert_eq!(reduce(matrix.clone()).matrix, expected);
        assert_eq!(reduce_incidence(matrix.clone()).matrix, expected);
        assert!(is_rref(&expected));
        let mut permuted = matrix;
        permuted.rows.reverse();
        assert_eq!(reduce(permuted).matrix, expected);
    }
    #[test]
    fn small_boolean_solution_sets_are_preserved() {
        for degree in 1..=3 {
            let terms: Vec<_> = (0u64..8)
                .filter(|m| m.count_ones() <= u32::from(degree))
                .collect();
            for selected in 0..(1usize << terms.len()) {
                let poly: Poly = terms
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| selected & (1 << i) != 0)
                    .map(|(_, &m)| m)
                    .collect();
                let input = Input {
                    n: 3,
                    degree,
                    active: 7,
                    polys: vec![poly.clone()],
                };
                let context = Packed::compile(&input, "ranked").unwrap();
                let answer = stream(&context, &input, CAPS).unwrap();
                assert_eq!(
                    reduce_incidence(direct(&input, CAPS).unwrap()).matrix,
                    answer.matrix
                );
                for point in 0u64..8 {
                    let value = |m: u64| u64::from(point & m == m);
                    let source_zero = poly.iter().fold(0, |s, &m| s ^ value(m)) == 0;
                    let reduced_zero = answer.matrix.rows.iter().all(|row| {
                        answer
                            .matrix
                            .columns
                            .iter()
                            .enumerate()
                            .fold(0, |s, (i, &m)| {
                                s ^ (((row[i / 64] >> (i % 64)) & 1) & value(m))
                            })
                            == 0
                    });
                    assert_eq!(source_zero, reduced_zero);
                }
            }
        }
    }
}
