//! Standalone, bounded Boolean linear-tail construction and elimination experiment. Uses only std.
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::hash::{BuildHasherDefault, Hasher};
use std::hint::black_box;
use std::rc::Rc;
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
}
type Built = Result<Matrix, Error>;

impl Input {
    fn canonical(mut self) -> Self {
        for poly in &mut self.polys {
            poly.sort_by(|&a, &b| source_mono_cmp(a, b));
        }
        self
    }
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
                || poly
                    .windows(2)
                    .any(|w| source_mono_cmp(w[0], w[1]) != std::cmp::Ordering::Less)
            {
                return Err(Error::InvalidInput);
            }
        }
        Ok(())
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
fn multipliers(active: u64, degree: u32) -> Vec<u64> {
    let variables: Vec<u64> = (0..64)
        .filter(|bit| active & (1u64 << bit) != 0)
        .map(|bit| 1u64 << bit)
        .collect();
    let mut out = vec![0u64];
    let mut level = vec![(0u64, 0usize)];
    for _ in 0..degree {
        let mut next = Vec::new();
        for &(monomial, start) in &level {
            for (index, &variable) in variables.iter().enumerate().skip(start) {
                next.push((monomial | variable, index + 1));
            }
        }
        out.extend(next.iter().map(|(m, _)| *m));
        level = next;
        if level.is_empty() {
            break;
        }
    }
    out
}

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
    let mut ordered: Vec<_> = incidence.keys().copied().collect();
    ordered.sort_by(|&a, &b| source_mono_cmp(a, b));
    for (column, monomial) in ordered.into_iter().enumerate() {
        columns.push(monomial);
        for &row in &incidence[&monomial] {
            rows[row][column / 64] |= 1 << (column % 64);
        }
    }
    Ok(Matrix { columns, rows })
}

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

fn source_mono_cmp(a: u64, b: u64) -> std::cmp::Ordering {
    let degrees = b.count_ones().cmp(&a.count_ones());
    if degrees != std::cmp::Ordering::Equal {
        return degrees;
    }
    let diff = a ^ b;
    if diff == 0 {
        return std::cmp::Ordering::Equal;
    }
    let bit = 63 - diff.leading_zeros();
    if a & (1u64 << bit) != 0 {
        std::cmp::Ordering::Greater
    } else {
        std::cmp::Ordering::Less
    }
}
fn order_key(m: u64) -> u64 {
    ((3 - m.count_ones()) as u64) << 56 | m
}
fn affine_bit(m: u64, n: u8) -> u64 {
    if m == 0 {
        1u64 << n
    } else {
        m
    }
}

// Same deterministic u64 hashing strategy as the retained source control.
#[derive(Default)]
struct FastHasher(u64);
impl Hasher for FastHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        let mut v = 0xcbf29ce484222325u64;
        for &b in bytes {
            v = (v ^ u64::from(b)).wrapping_mul(0x100000001b3);
        }
        self.write_u64(v);
    }
    fn write_u64(&mut self, value: u64) {
        let mut z = value.wrapping_add(0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        self.0 = z ^ (z >> 31);
    }
}
type FastMap<V> = HashMap<u64, V, BuildHasherDefault<FastHasher>>;
#[derive(Default)]
struct Context {
    schedules: FastMap<Rc<[u64]>>,
    layouts: FastMap<Rc<Layout>>,
    censuses: FastMap<Rc<Census>>,
}
impl Context {
    fn schedule(&mut self, active: u64, gap: u8) -> Rc<[u64]> {
        let key = active | (u64::from(gap) << 56);
        self.schedules
            .entry(key)
            .or_insert_with(|| multipliers(active, u32::from(gap)).into())
            .clone()
    }
    fn retained_counts(&self) -> (usize, usize, usize) {
        (
            self.schedules.values().map(|s| s.len()).sum(),
            self.layouts.values().map(|l| l.columns.len()).sum(),
            self.censuses.values().map(|c| c.choose.len() * 4).sum(),
        )
    }
    fn census(&mut self, n: u8, degree: u8) -> Rc<Census> {
        self.censuses
            .entry((u64::from(n) << 8) | u64::from(degree))
            .or_insert_with(|| Rc::new(Census::new(n, degree)))
            .clone()
    }
}

// A collision-free index for the exact column census, not a matrix layout.
// Within weight k, rank(m)=sum_i C(b_i,i), with set bits b_i in ascending order.
struct Census {
    choose: Vec<[usize; 4]>,
    offset: [usize; 4],
    columns: usize,
}
impl Census {
    fn new(n: u8, degree: u8) -> Self {
        let mut choose = vec![[0usize; 4]; n as usize + 1];
        for b in 0..=n as usize {
            choose[b][0] = 1;
            if b > 0 {
                for k in 1..=3 {
                    choose[b][k] = choose[b - 1][k] + choose[b - 1][k - 1];
                }
            }
        }
        let mut offset = [0usize; 4];
        let mut columns = 0;
        for k in (0..=degree as usize).rev() {
            offset[k] = columns;
            columns += choose[n as usize][k];
        }
        Self {
            choose,
            offset,
            columns,
        }
    }
    fn index(&self, mut m: u64) -> usize {
        let mut index = self.offset[m.count_ones() as usize];
        let mut ordinal = 1;
        while m != 0 {
            index += self.choose[m.trailing_zeros() as usize][ordinal];
            ordinal += 1;
            m &= m - 1;
        }
        index
    }
}
enum Observed {
    Hash(FastMap<()>),
    Bitmap { index: Rc<Census>, seen: Vec<u64> },
}
impl Observed {
    fn insert(&mut self, m: u64) {
        match self {
            Self::Hash(set) => {
                set.insert(m, ());
            }
            Self::Bitmap { index, seen } => {
                let c = index.index(m);
                seen[c / 64] |= 1u64 << (c % 64);
            }
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Hash(set) => set.len(),
            Self::Bitmap { seen, .. } => seen.iter().map(|w| w.count_ones() as usize).sum(),
        }
    }
}
#[cfg(test)]
fn monomial_bound(variables: u32, degree: u8) -> usize {
    let (mut choose, mut total) = (1usize, 1usize);
    for k in 1..=variables.min(u32::from(degree)) {
        choose = choose * (variables + 1 - k) as usize / k as usize;
        total += choose;
    }
    total
}
fn prefer_flat(active: u64, degree: u8) -> bool {
    // Exactly equivalent to sum_{j<=degree} C(popcount(active),j) <= 512.
    // At degree 2: C(31,<=2)=497, C(32,<=2)=529.
    // At degree 3: C(14,<=3)=470, C(15,<=3)=576.
    match degree {
        0 | 1 => true,
        2 => active.count_ones() <= 31,
        3 => active.count_ones() <= 14,
        _ => false,
    }
}
fn product_into(poly: &[u64], mult: u64, out: &mut Vec<u64>) {
    out.clear();
    out.extend(poly.iter().map(|&m| m | mult));
    out.sort_unstable();
    let (mut read, mut write) = (0, 0);
    while read < out.len() {
        let mut end = read + 1;
        while end < out.len() && out[end] == out[read] {
            end += 1;
        }
        if (end - read) % 2 == 1 {
            out[write] = out[read];
            write += 1;
        }
        read = end;
    }
    out.truncate(write);
}
fn source_products(input: &Input, context: &mut Context, caps: Caps) -> Result<Vec<Poly>, Error> {
    input.validate()?;
    let mut rows = Vec::new();
    let mut scratch = Vec::new();
    for poly in &input.polys {
        let Some(d) = poly.iter().map(|m| m.count_ones() as u8).max() else {
            continue;
        };
        if d > input.degree {
            continue;
        }
        for &t in context.schedule(input.active, input.degree - d).iter() {
            product_into(poly, t, &mut scratch);
            if !scratch.is_empty() {
                if rows.len() == caps.rows {
                    return Err(Error::RowCap);
                }
                rows.push(scratch.clone());
            }
        }
    }
    Ok(rows)
}
struct Layout {
    columns: Vec<u64>,
    index: FastMap<usize>,
    low_start: usize,
}
impl Layout {
    fn from_rows(rows: &[Poly], caps: Caps) -> Result<Self, Error> {
        let mut columns: Vec<_> = rows.iter().flatten().copied().collect();
        columns.sort_unstable();
        columns.dedup();
        check_shape(rows.len(), columns.len(), caps)?;
        columns.sort_unstable_by_key(|&m| order_key(m));
        let low_start = columns
            .iter()
            .position(|m| m.count_ones() <= 1)
            .unwrap_or(columns.len());
        let index = columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        Ok(Self {
            columns,
            index,
            low_start,
        })
    }
}
struct Flat {
    data: Vec<u64>,
    rows: usize,
    words: usize,
}
fn pack_flat(rows: &[Poly], layout: &Layout) -> Flat {
    let words = layout.columns.len().div_ceil(64);
    let mut data = vec![0u64; rows.len() * words];
    for (i, row) in rows.iter().enumerate() {
        for m in row {
            let c = layout.index[m];
            data[i * words + c / 64] |= 1u64 << (c % 64);
        }
    }
    Flat {
        data,
        rows: rows.len(),
        words,
    }
}
fn fused_cached(
    input: &Input,
    context: &mut Context,
    layout: &Layout,
    caps: Caps,
) -> Result<Option<Flat>, Error> {
    if layout.columns.len() > caps.columns {
        return Ok(None);
    }
    let words = layout.columns.len().div_ceil(64);
    let mut seen = vec![false; layout.columns.len()];
    let mut data = Vec::new();
    let mut rows = 0;
    let mut product = Vec::new();
    for poly in &input.polys {
        let Some(d) = poly.iter().map(|m| m.count_ones() as u8).max() else {
            continue;
        };
        if d > input.degree {
            continue;
        }
        for &t in context.schedule(input.active, input.degree - d).iter() {
            product_into(poly, t, &mut product);
            if product.is_empty() {
                continue;
            }
            if rows == caps.rows {
                return Err(Error::RowCap);
            }
            let start = data.len();
            data.resize(start + words, 0u64);
            for m in &product {
                let Some(&c) = layout.index.get(m) else {
                    return Ok(None);
                };
                data[start + c / 64] |= 1u64 << (c % 64);
                seen[c] = true;
            }
            rows += 1;
        }
    }
    if seen.iter().any(|v| !*v) {
        return Ok(None);
    }
    Ok(Some(Flat { data, rows, words }))
}
fn build_flat(
    input: &Input,
    context: &mut Context,
    caps: Caps,
) -> Result<(Flat, Rc<Layout>, bool), Error> {
    input.validate()?;
    let key = input.active | (u64::from(input.degree) << 40) | (u64::from(input.n) << 48);
    if let Some(layout) = context.layouts.get(&key).cloned() {
        if let Some(matrix) = fused_cached(input, context, &layout, caps)? {
            return Ok((matrix, layout, true));
        }
    }
    let rows = source_products(input, context, caps)?;
    let layout = Rc::new(Layout::from_rows(&rows, caps)?);
    let matrix = pack_flat(&rows, &layout);
    context.layouts.insert(key, layout.clone());
    Ok((matrix, layout, false))
}

struct LowBasis {
    pivots: [u64; 37],
}
impl Default for LowBasis {
    fn default() -> Self {
        Self { pivots: [0; 37] }
    }
}
impl LowBasis {
    fn insert(&mut self, mut row: u64) {
        while row != 0 {
            let p = row.trailing_zeros() as usize;
            if self.pivots[p] == 0 {
                self.pivots[p] = row;
                return;
            }
            row ^= self.pivots[p];
        }
    }
    fn finish(mut self, n: u8) -> Vec<u64> {
        for p in (0..=n as usize).rev() {
            if self.pivots[p] == 0 {
                continue;
            }
            for q in 0..p {
                if self.pivots[q] & (1u64 << p) != 0 {
                    self.pivots[q] ^= self.pivots[p];
                }
            }
        }
        self.pivots[..=n as usize]
            .iter()
            .copied()
            .filter(|&v| v != 0)
            .collect()
    }
}
#[derive(Debug, PartialEq, Eq)]
struct Answer {
    tail: Vec<u64>,
    source_rows: usize,
    source_columns: usize,
    high_rank: usize,
}
#[derive(Default)]
struct Work {
    high_xors: u64,
    word_xors: u64,
    merge_items: u64,
    pivot_exchanges: usize,
    max_sparse_terms: usize,
    layout_hit: bool,
    flat_dispatch: bool,
}

// Faithful generic port of the retained flat high-column elimination loop.
fn flat_tail(mut matrix: Flat, layout: &Layout, n: u8, work: &mut Work) -> Answer {
    let source_rows = matrix.rows;
    let mut pivot_row = 0;
    for column in 0..layout.low_start {
        let (word, bit) = (column / 64, 1u64 << (column % 64));
        let Some(pivot) =
            (pivot_row..matrix.rows).find(|&r| matrix.data[r * matrix.words + word] & bit != 0)
        else {
            continue;
        };
        if pivot != pivot_row {
            for i in 0..matrix.words {
                matrix
                    .data
                    .swap(pivot_row * matrix.words + i, pivot * matrix.words + i);
            }
        }
        let split = (pivot_row + 1) * matrix.words;
        let (head, tail) = matrix.data.split_at_mut(split);
        let pivot = &head[pivot_row * matrix.words..(pivot_row + 1) * matrix.words];
        for row in tail.chunks_exact_mut(matrix.words) {
            if row[word] & bit != 0 {
                for (to, &from) in row[word..].iter_mut().zip(&pivot[word..]) {
                    *to ^= from;
                }
                work.high_xors += 1;
                work.word_xors += (matrix.words - word) as u64;
            }
        }
        pivot_row += 1;
        if pivot_row == matrix.rows {
            break;
        }
    }
    let mut low = LowBasis::default();
    if matrix.words != 0 {
        for row in matrix.data[pivot_row * matrix.words..].chunks_exact(matrix.words) {
            let mut affine = 0;
            for c in layout.low_start..layout.columns.len() {
                if row[c / 64] & (1u64 << (c % 64)) != 0 {
                    affine ^= affine_bit(layout.columns[c], n);
                }
            }
            low.insert(affine);
        }
    }
    Answer {
        tail: low.finish(n),
        source_rows,
        source_columns: layout.columns.len(),
        high_rank: pivot_row,
    }
}
fn xor_indices(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
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
// Retained lowest-weight sparse-bucket strategy from sparse_macaulay.rs.
fn bucket_tail(products: Vec<Poly>, layout: &Layout, n: u8, work: &mut Work) -> Answer {
    let source_rows = products.len();
    let mut rows: Vec<Vec<u32>> = products
        .into_iter()
        .map(|r| {
            let mut row: Vec<_> = r.iter().map(|m| layout.index[m] as u32).collect();
            row.sort_unstable();
            row
        })
        .collect();
    work.max_sparse_terms = rows.iter().map(Vec::len).max().unwrap_or(0);
    let mut buckets = vec![Vec::<usize>::new(); layout.columns.len().max(1)];
    for (i, row) in rows.iter().enumerate() {
        if let Some(&c) = row.first() {
            buckets[c as usize].push(i);
        }
    }
    let mut pivot_flags = vec![false; rows.len()];
    let mut high_rank = 0;
    for c in 0..layout.low_start {
        let bucket = std::mem::take(&mut buckets[c]);
        if bucket.is_empty() {
            continue;
        }
        let pivot = *bucket.iter().min_by_key(|&&i| rows[i].len()).unwrap();
        pivot_flags[pivot] = true;
        high_rank += 1;
        let base = rows[pivot].clone();
        for i in bucket {
            if i == pivot {
                continue;
            }
            work.high_xors += 1;
            work.merge_items += (rows[i].len() + base.len()) as u64;
            let row = xor_indices(&rows[i], &base);
            work.max_sparse_terms = work.max_sparse_terms.max(row.len());
            if let Some(&lead) = row.first() {
                buckets[lead as usize].push(i);
            }
            rows[i] = row;
        }
    }
    let mut low = LowBasis::default();
    for (i, row) in rows.into_iter().enumerate() {
        if pivot_flags[i] || row.is_empty() || (row[0] as usize) < layout.low_start {
            continue;
        }
        low.insert(
            row.iter()
                .fold(0, |a, &c| a ^ affine_bit(layout.columns[c as usize], n)),
        );
    }
    Answer {
        tail: low.finish(n),
        source_rows,
        source_columns: layout.columns.len(),
        high_rank,
    }
}

struct HighRow {
    high: Vec<u64>,
    low: u64,
}
fn xor_keys(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
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
fn sparse_stream(
    input: &Input,
    context: &mut Context,
    caps: Caps,
    exchange: bool,
    census: bool,
    work: &mut Work,
) -> Result<Answer, Error> {
    input.validate()?;
    let mut observed = if census {
        let index = context.census(input.n, input.degree);
        let seen = vec![0u64; index.columns.div_ceil(64)];
        Observed::Bitmap { index, seen }
    } else {
        Observed::Hash(FastMap::default())
    };
    let mut pivots = FastMap::<usize>::default();
    let mut bases = Vec::<HighRow>::new();
    let mut tail = LowBasis::default();
    let mut source_rows = 0;
    for poly in &input.polys {
        let Some(d) = poly.iter().map(|m| m.count_ones() as u8).max() else {
            continue;
        };
        if d > input.degree {
            continue;
        }
        for &t in context.schedule(input.active, input.degree - d).iter() {
            let mut high = Vec::with_capacity(poly.len());
            let mut low = 0u64;
            for &m in poly {
                let product = m | t;
                if product.count_ones() > 1 {
                    high.push(order_key(product));
                } else {
                    low ^= affine_bit(product, input.n);
                }
            }
            high.sort_unstable();
            let (mut read, mut write) = (0, 0);
            while read < high.len() {
                let mut end = read + 1;
                while end < high.len() && high[end] == high[read] {
                    end += 1;
                }
                if (end - read) % 2 == 1 {
                    high[write] = high[read];
                    write += 1;
                }
                read = end;
            }
            high.truncate(write);
            if high.is_empty() && low == 0 {
                continue;
            }
            if source_rows == caps.rows {
                return Err(Error::RowCap);
            }
            source_rows += 1;
            for &m in &high {
                observed.insert(m & ((1u64 << 56) - 1));
            }
            let mut bits = low;
            while bits != 0 {
                let bit = bits.trailing_zeros();
                observed.insert(if bit == u32::from(input.n) {
                    0
                } else {
                    1u64 << bit
                });
                bits &= bits - 1;
            }
            work.max_sparse_terms = work
                .max_sparse_terms
                .max(high.len() + low.count_ones() as usize);
            loop {
                let Some(&lead) = high.first() else {
                    tail.insert(low);
                    break;
                };
                if let Some(&i) = pivots.get(&lead) {
                    let base = &mut bases[i];
                    if exchange
                        && (high.len() + low.count_ones() as usize)
                            < base.high.len() + base.low.count_ones() as usize
                    {
                        std::mem::swap(&mut high, &mut base.high);
                        std::mem::swap(&mut low, &mut base.low);
                        work.pivot_exchanges += 1;
                    }
                    work.high_xors += 1;
                    work.merge_items += (high.len() + base.high.len()) as u64;
                    high = xor_keys(&high, &base.high);
                    low ^= base.low;
                    work.max_sparse_terms = work
                        .max_sparse_terms
                        .max(high.len() + low.count_ones() as usize);
                } else {
                    pivots.insert(lead, bases.len());
                    bases.push(HighRow { high, low });
                    break;
                }
            }
        }
    }
    let source_columns = observed.len();
    check_shape(source_rows, source_columns, caps)?;
    Ok(Answer {
        tail: tail.finish(input.n),
        source_rows,
        source_columns,
        high_rank: bases.len(),
    })
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kernel {
    Flat,
    Bucket,
    Plain,
    Exchange,
    Census,
    Hybrid,
}
fn select_kernel(name: &str) -> Kernel {
    match name {
        "flat_cached" | "flat_alias" => Kernel::Flat,
        "sparse_bucket" => Kernel::Bucket,
        "stream_plain" => Kernel::Plain,
        "stream_exchange" => Kernel::Exchange,
        "stream_census" => Kernel::Census,
        "hybrid_census" => Kernel::Hybrid,
        _ => panic!("unknown arm"),
    }
}
#[cfg(test)]
fn apply(
    name: &str,
    input: &Input,
    context: &mut Context,
    caps: Caps,
) -> Result<(Answer, Work), Error> {
    run_kernel(select_kernel(name), input, context, caps)
}
#[inline(never)]
fn flat_leaf(input: &Input, context: &mut Context, caps: Caps) -> Result<(Answer, Work), Error> {
    let mut work = Work::default();
    work.flat_dispatch = true;
    let (matrix, layout, hit) = build_flat(input, context, caps)?;
    work.layout_hit = hit;
    Ok((flat_tail(matrix, &layout, input.n, &mut work), work))
}
fn run_kernel(
    kernel: Kernel,
    input: &Input,
    context: &mut Context,
    caps: Caps,
) -> Result<(Answer, Work), Error> {
    let kernel = if kernel == Kernel::Hybrid {
        if prefer_flat(input.active, input.degree) {
            Kernel::Flat
        } else {
            Kernel::Census
        }
    } else {
        kernel
    };
    let mut work = Work::default();
    let answer = match kernel {
        Kernel::Flat => return flat_leaf(input, context, caps),
        Kernel::Bucket => {
            let rows = source_products(input, context, caps)?;
            let layout = Layout::from_rows(&rows, caps)?;
            bucket_tail(rows, &layout, input.n, &mut work)
        }
        Kernel::Plain => sparse_stream(input, context, caps, false, false, &mut work)?,
        Kernel::Exchange => sparse_stream(input, context, caps, true, false, &mut work)?,
        Kernel::Census => sparse_stream(input, context, caps, true, true, &mut work)?,
        Kernel::Hybrid => unreachable!(),
    };
    Ok((answer, work))
}
fn reference(input: &Input, caps: Caps) -> Result<Answer, Error> {
    let matrix = oracle(input, caps)?;
    let source_rows = matrix.rows.len();
    let source_columns = matrix.columns.len();
    let reduced = oracle_reduce(matrix);
    let high_start = reduced
        .columns
        .iter()
        .position(|m| m.count_ones() <= 1)
        .unwrap_or(source_columns);
    let mut tail = Vec::new();
    let mut high_rank = 0;
    for row in &reduced.rows {
        if (0..high_start).any(|c| row[c / 64] & (1u64 << (c % 64)) != 0) {
            high_rank += 1;
            continue;
        }
        let mut low = 0;
        for c in high_start..source_columns {
            if row[c / 64] & (1u64 << (c % 64)) != 0 {
                low ^= affine_bit(reduced.columns[c], input.n);
            }
        }
        if low != 0 {
            tail.push(low);
        }
    }
    Ok(Answer {
        tail,
        source_rows,
        source_columns,
        high_rank,
    })
}

fn next(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *seed;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn fixtures(n: u8, seed: u64, batch: usize, family: &str) -> Vec<Input> {
    assert!([
        "quadratic",
        "linear_drop",
        "restricted_cycle",
        "cross_cancel"
    ]
    .contains(&family));
    let positions: Vec<u8> = if family == "restricted_cycle" {
        (0..8).map(|i| i * n / 8).collect()
    } else {
        (0..n).collect()
    };
    let mut state = seed;
    let mut supports = Vec::new();
    for _ in 0..8 {
        let mut terms = BTreeSet::new();
        while terms.len() < 16 {
            let a = positions[next(&mut state) as usize % positions.len()];
            let b = positions[next(&mut state) as usize % positions.len()];
            if a != b {
                terms.insert((1u64 << a) | (1u64 << b));
            }
        }
        terms.insert(0);
        terms.extend(positions.iter().take(10).map(|&i| 1u64 << i));
        supports.push(terms.into_iter().collect::<Poly>());
    }
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
            if family == "linear_drop" || family == "restricted_cycle" {
                match i % 4 {
                    1 => {
                        polys[0].retain(|m| m.count_ones() == 1);
                        polys[0].push(1u64 << positions[0]);
                        polys[0].sort_unstable();
                        polys[0].dedup();
                    }
                    2 if family == "restricted_cycle" => polys[0] = vec![0],
                    3 => polys[0].clear(),
                    _ => (),
                }
            }
            if family == "cross_cancel" {
                let mut second: BTreeSet<_> = polys[0].iter().copied().collect();
                let variable = 1u64 << positions[i % positions.len()];
                if !second.insert(variable) {
                    second.remove(&variable);
                }
                if i % 2 == 1 && !second.insert(0) {
                    second.remove(&0);
                }
                polys[1] = second.into_iter().collect();
            }
            let active = polys.iter().flatten().fold(0, |mask, &m| mask | m);
            Input {
                n,
                degree: 3,
                active,
                polys,
            }
            .canonical()
        })
        .collect()
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        [6, 7, 8, 9].contains(&args.len()),
        "worker N SEED BATCH FAMILY REPETITIONS [with-census [balanced-order]]"
    );
    let include_census = args.len() >= 7;
    if include_census {
        assert_eq!(args[6], "with-census");
    }
    let balanced = args.len() >= 8;
    let matched = balanced && args[7] == "matched-order";
    if balanced {
        assert!(["balanced-order", "matched-order"].contains(&args[7].as_str()));
    }
    let inner_batches = if args.len() == 9 {
        assert!(matched);
        args[8].parse::<usize>().unwrap()
    } else {
        1
    };
    assert!([1, 16].contains(&inner_batches));
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let batch = args[3].parse::<usize>().unwrap();
    let family = &args[4];
    let repetitions = args[5].parse::<usize>().unwrap();
    assert!(
        [12, 20, 28, 36].contains(&n)
            && [1, 4, 8].contains(&batch)
            && (1..=42).contains(&repetitions)
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
        .map(|input| reference(input, CAPS).unwrap())
        .collect();
    let serialized: Vec<_> = inputs
        .iter()
        .map(|x| format!("{{\"active\":{},\"polys\":{:?}}}", x.active, x.polys))
        .collect();
    let expected_tails: Vec<_> = expected.iter().map(|r| &r.tail).collect();
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"batch\":{batch},\"family\":\"{family}\",\"degree\":3,\"inputs\":[{}],\"reference_tails\":{:?}}}",serialized.join(","),expected_tails);
    let names: &[&str] = if balanced {
        &[
            "flat_cached",
            "sparse_bucket",
            "stream_plain",
            "stream_exchange",
            "stream_census",
            "hybrid_census",
            "flat_alias",
        ]
    } else if include_census {
        &[
            "flat_cached",
            "sparse_bucket",
            "stream_plain",
            "stream_exchange",
            "stream_census",
            "hybrid_census",
        ]
    } else {
        &[
            "flat_cached",
            "sparse_bucket",
            "stream_plain",
            "stream_exchange",
        ]
    };
    for (rep, order, variant, measured) in
        measurement_order(names.len(), repetitions, balanced, matched)
    {
        let name = names[variant];
        let kernel = select_kernel(name); // experimental label decoding is harness work
        let start = Instant::now();
        let mut setup_ns = 0u128;
        let (mut apply_ns, mut validation_ns) = (0u128, 0u128);
        let (
            mut source_rows,
            mut source_columns,
            mut high_rank,
            mut tail_rank,
            mut nonempty_tails,
            mut layout_hits,
            mut max_sparse_terms,
            mut pivot_exchanges,
        ) = (0, 0, 0, 0, 0, 0, 0, 0);
        let (mut high_xors, mut word_xors, mut merge_items) = (0u64, 0u64, 0u64);
        let mut flat_dispatches = 0;
        let mut retained_schedule_masks = 0usize;
        let mut retained_layout_columns = 0usize;
        let mut retained_census_entries = 0usize;
        let repetitions_inside = inner_batches;
        for _ in 0..repetitions_inside {
            let tick = Instant::now();
            let mut context = Context::default();
            setup_ns += tick.elapsed().as_nanos();
            for (input, expected) in inputs.iter().zip(&expected) {
                let tick = Instant::now();
                let (answer, work) =
                    run_kernel(kernel, black_box(input), &mut context, CAPS).unwrap();
                apply_ns += tick.elapsed().as_nanos();
                let tick = Instant::now();
                assert_eq!(black_box(&answer), expected);
                source_rows += answer.source_rows;
                source_columns += answer.source_columns;
                high_rank += answer.high_rank;
                tail_rank += answer.tail.len();
                nonempty_tails += usize::from(!answer.tail.is_empty());
                layout_hits += usize::from(work.layout_hit);
                flat_dispatches += usize::from(work.flat_dispatch);
                max_sparse_terms = max_sparse_terms.max(work.max_sparse_terms);
                pivot_exchanges += work.pivot_exchanges;
                high_xors += work.high_xors;
                word_xors += work.word_xors;
                merge_items += work.merge_items;
                validation_ns += tick.elapsed().as_nanos();
            }
            let (s, l, c) = context.retained_counts();
            retained_schedule_masks = retained_schedule_masks.max(s);
            retained_layout_columns = retained_layout_columns.max(l);
            retained_census_entries = retained_census_entries.max(c);
            drop(context);
        }
        let sample_outputs = batch * repetitions_inside;
        let total_ns = start.elapsed().as_nanos();
        let word_count = if flat_dispatches > 0 {
            word_xors.to_string()
        } else {
            "null".into()
        };
        let merge_count = if flat_dispatches == sample_outputs {
            "null".into()
        } else {
            merge_items.to_string()
        };
        let fill = if flat_dispatches == sample_outputs {
            "null".into()
        } else {
            max_sparse_terms.to_string()
        };
        if measured {
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{name}\",\"setup_ns\":{setup_ns},\"apply_ns\":{apply_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"source_rows\":{source_rows},\"source_columns\":{source_columns},\"high_rank\":{high_rank},\"tail_rank\":{tail_rank},\"nonempty_tails\":{nonempty_tails},\"layout_hits\":{layout_hits},\"high_xors\":{high_xors},\"high_word_xors\":{word_count},\"merge_items\":{merge_count},\"max_sparse_terms\":{fill},\"pivot_exchanges\":{pivot_exchanges},\"retained_schedule_masks\":{retained_schedule_masks},\"retained_layout_columns\":{retained_layout_columns},\"flat_dispatches\":{flat_dispatches},\"retained_census_entries\":{retained_census_entries},\"verified_outputs\":{sample_outputs},\"inner_batches\":{repetitions_inside}}}");
        }
    }
}

// An Euler cycle on the complete directed arm graph, including self edges.
// Each method follows every method exactly once per cycle. The unmeasured
// primer supplies the first predecessor; all measured contexts remain cold.
fn measurement_order(
    arms: usize,
    repetitions: usize,
    balanced: bool,
    matched: bool,
) -> Vec<(usize, usize, usize, bool)> {
    if !balanced {
        return (0..repetitions)
            .flat_map(|rep| (0..arms).map(move |order| (rep, order, (rep + order) % arms, true)))
            .collect();
    }
    assert_eq!(repetitions % arms, 0);
    let mut next = vec![0usize; arms];
    let mut stack = vec![0usize];
    let mut circuit = Vec::new();
    while let Some(&u) = stack.last() {
        if next[u] < arms {
            let v = next[u];
            next[u] += 1;
            stack.push(v);
        } else {
            circuit.push(stack.pop().unwrap());
        }
    }
    circuit.reverse();
    assert_eq!(circuit.len(), arms * arms + 1);
    assert_eq!(circuit[0], 0);
    assert_eq!(*circuit.last().unwrap(), 0);
    let mut visits = vec![0usize; arms];
    let mut out = vec![(0, 0, 0, false)];
    let mut previous = 0;
    for cycle in 0..repetitions / arms {
        for &arm in &circuit[1..] {
            let rep = if matched {
                cycle * arms + previous
            } else {
                visits[arm]
            };
            out.push((rep, out.len() - 1, arm, true));
            visits[arm] += 1;
            previous = arm;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn span(rows: &[u64]) -> BTreeSet<u64> {
        (0..(1usize << rows.len()))
            .map(|selected| {
                rows.iter().enumerate().fold(0, |a, (i, &r)| {
                    a ^ if selected & (1 << i) != 0 { r } else { 0 }
                })
            })
            .collect()
    }
    #[test]
    fn encoded_order_matches_source_degrevlex() {
        for n in 1..=12 {
            let mut source: Vec<_> = (0u64..(1 << n)).filter(|m| m.count_ones() <= 3).collect();
            let mut keyed = source.clone();
            source.sort_by(|&a, &b| source_mono_cmp(a, b));
            keyed.sort_unstable_by_key(|&m| order_key(m));
            assert_eq!(source, keyed);
        }
    }
    #[test]
    fn source_multiplier_layers_and_input_order_are_exact() {
        assert_eq!(multipliers(7, 2), vec![0, 1, 2, 4, 3, 5, 6]);
        assert_eq!(multipliers(0b10101, 2), vec![0, 1, 4, 16, 5, 17, 20]);
        let ordered = Input {
            n: 3,
            degree: 3,
            active: 7,
            polys: vec![vec![0, 4, 3, 1, 7]],
        }
        .canonical();
        assert_eq!(ordered.polys[0], vec![7, 3, 1, 4, 0]);
        assert_eq!(ordered.validate(), Ok(()));
        let wrong = Input {
            polys: vec![vec![0, 1, 3]],
            ..ordered
        };
        assert_eq!(wrong.validate(), Err(Error::InvalidInput));
    }
    #[test]
    fn predecessor_schedule_balances_all_directed_pairs_and_cold_samples() {
        for arms in [4, 6, 7] {
            let schedule = measurement_order(arms, arms * 2, true, false);
            let mut counts = vec![vec![0usize; arms]; arms];
            let mut per_arm = vec![0usize; arms];
            let mut previous = 0;
            assert!(!schedule[0].3);
            for &(rep, order, arm, measured) in &schedule[1..] {
                assert!(measured);
                assert_eq!(rep, per_arm[arm]);
                assert_eq!(order, per_arm.iter().sum());
                counts[previous][arm] += 1;
                per_arm[arm] += 1;
                previous = arm;
            }
            assert!(counts.iter().flatten().all(|&n| n == 2));
            assert!(per_arm.iter().all(|&n| n == arms * 2));
        }
    }
    #[test]
    fn flat_alias_is_the_identical_mathematical_control() {
        assert_eq!(select_kernel("flat_cached"), select_kernel("flat_alias"));
        let inputs = fixtures(20, 31337, 4, "restricted_cycle");
        let (mut a, mut b) = (Context::default(), Context::default());
        for input in inputs {
            let (x, wx) = apply("flat_cached", &input, &mut a, CAPS).unwrap();
            let (y, wy) = apply("flat_alias", &input, &mut b, CAPS).unwrap();
            assert_eq!(x, y);
            assert_eq!(wx.word_xors, wy.word_xors);
            assert_eq!(wx.layout_hit, wy.layout_hit);
        }
    }
    #[test]
    fn matched_pairs_have_the_same_predecessor_and_cycle() {
        let arms = 7;
        let schedule = measurement_order(arms, 14, true, true);
        let mut previous = 0;
        let mut pairs = BTreeMap::new();
        for &(rep, event, arm, _) in &schedule[1..] {
            assert_eq!(rep, (event / (arms * arms)) * arms + previous);
            assert!(pairs.insert((rep, arm), previous).is_none());
            previous = arm;
        }
        for rep in 0..14 {
            for arm in 0..arms {
                assert_eq!(pairs[&(rep, arm)], rep % arms);
            }
        }
    }
    #[test]
    fn census_is_injective_and_counts_actual_columns() {
        for n in [3, 8, 12, 20, 28, 36] {
            for degree in 1..=3 {
                let c = Census::new(n, degree);
                let mut monomials = reference_multipliers((1u64 << n) - 1, u32::from(degree));
                monomials.sort_unstable_by_key(|&m| order_key(m));
                assert_eq!(c.columns, monomials.len());
                for (i, &m) in monomials.iter().enumerate() {
                    assert_eq!(c.index(m), i);
                }
                let mut observed = Observed::Bitmap {
                    seen: vec![0; c.columns.div_ceil(64)],
                    index: Rc::new(c),
                };
                for &m in monomials.iter().step_by(3) {
                    observed.insert(m);
                    observed.insert(m);
                }
                assert_eq!(observed.len(), monomials.iter().step_by(3).count());
            }
        }
        assert_eq!(monomial_bound(8, 3), 93);
        assert_eq!(monomial_bound(14, 3), 470);
        assert_eq!(monomial_bound(15, 3), 576);
    }
    #[test]
    fn hybrid_dispatch_and_bitmap_caps_preserve_the_same_contract() {
        for n in [12, 20, 36] {
            let input = fixtures(n, 31337, 1, "quadratic").pop().unwrap();
            let (answer, work) =
                apply("hybrid_census", &input, &mut Context::default(), CAPS).unwrap();
            assert_eq!(answer, reference(&input, CAPS).unwrap());
            assert_eq!(
                work.flat_dispatch,
                monomial_bound(input.active.count_ones(), 3) <= 512
            );
        }
        let input = Input {
            n: 3,
            degree: 3,
            active: 7,
            polys: vec![vec![1, 2]],
        }
        .canonical();
        for caps in [
            Caps {
                rows: 1,
                columns: 8,
            },
            Caps {
                rows: 4096,
                columns: 0,
            },
        ] {
            let expected = reference(&input, caps).unwrap_err();
            assert_eq!(
                apply("stream_census", &input, &mut Context::default(), caps).err(),
                Some(expected)
            );
        }
    }
    #[test]
    fn constant_cutoffs_are_exactly_the_original_dispatch_rule() {
        for variables in 0..=64u32 {
            let mask = if variables == 64 {
                u64::MAX
            } else {
                (1u64 << variables) - 1
            };
            for degree in 0..=3 {
                assert_eq!(
                    prefer_flat(mask, degree),
                    monomial_bound(variables, degree) <= 512
                );
            }
        }
        assert!(!prefer_flat(0, 4));
    }
    #[test]
    fn small_split_matrices_match_explicit_intersection() {
        for cols in 1..=4 {
            for count in 1..=3 {
                for code in 0..(1u64 << (cols * count)) {
                    for high in 0..=cols {
                        let columns: Vec<u64> = [7, 3, 5, 6]
                            .into_iter()
                            .take(high)
                            .chain([1, 2, 4, 0].into_iter().take(cols - high))
                            .collect();
                        let rows: Vec<_> = (0..count)
                            .map(|i| vec![(code >> (i * cols)) & ((1 << cols) - 1)])
                            .collect();
                        let full = oracle_reduce(Matrix {
                            columns: columns.clone(),
                            rows: rows.clone(),
                        });
                        let mut expected = Vec::new();
                        for row in full.rows {
                            if row[0] & ((1 << high) - 1) == 0 {
                                let low = (high..cols).fold(0, |a, c| {
                                    a ^ if row[0] & (1 << c) != 0 {
                                        affine_bit(columns[c], 3)
                                    } else {
                                        0
                                    }
                                });
                                if low != 0 {
                                    expected.push(low);
                                }
                            }
                        }
                        let layout = Layout {
                            index: columns.iter().enumerate().map(|(i, &m)| (m, i)).collect(),
                            columns,
                            low_start: high,
                        };
                        let flat = Flat {
                            data: rows.iter().map(|r| r[0]).collect(),
                            rows: count,
                            words: 1,
                        };
                        let got = flat_tail(flat, &layout, 3, &mut Work::default());
                        assert_eq!(got.tail, expected);
                        let all = span(&rows.iter().map(|r| r[0]).collect::<Vec<_>>());
                        let intersection: BTreeSet<_> = all
                            .into_iter()
                            .filter(|r| r & ((1 << high) - 1) == 0)
                            .map(|r| {
                                (high..cols).fold(0, |a, c| {
                                    a ^ if r & (1 << c) != 0 {
                                        affine_bit(layout.columns[c], 3)
                                    } else {
                                        0
                                    }
                                })
                            })
                            .collect();
                        assert_eq!(span(&got.tail), intersection);
                    }
                }
            }
        }
    }
    #[test]
    fn exhaustive_boolean_coefficients_all_arms() {
        for degree in 1..=3 {
            for active in 0..8 {
                let mut contexts: [Context; 6] = std::array::from_fn(|_| Context::default());
                for selected in 0..256u16 {
                    let input = Input {
                        n: 3,
                        degree,
                        active,
                        polys: vec![(0..8).filter(|m| selected & (1 << m) != 0).collect()],
                    }
                    .canonical();
                    let expected = reference(&input, CAPS).unwrap();
                    for (name, context) in [
                        "flat_cached",
                        "sparse_bucket",
                        "stream_plain",
                        "stream_exchange",
                        "stream_census",
                        "hybrid_census",
                    ]
                    .iter()
                    .zip(&mut contexts)
                    {
                        assert_eq!(apply(name, &input, context, CAPS).unwrap().0, expected);
                    }
                }
            }
        }
    }
    #[test]
    fn all_larger_families_and_hidden_affine_consequences_match() {
        for n in [12, 20, 28, 36] {
            for family in [
                "quadratic",
                "linear_drop",
                "restricted_cycle",
                "cross_cancel",
            ] {
                let inputs = fixtures(n, 31337, 4, family);
                let mut contexts: [Context; 6] = std::array::from_fn(|_| Context::default());
                for (i, input) in inputs.iter().enumerate() {
                    let expected = reference(input, CAPS).unwrap();
                    for (name, context) in [
                        "flat_cached",
                        "sparse_bucket",
                        "stream_plain",
                        "stream_exchange",
                        "stream_census",
                        "hybrid_census",
                    ]
                    .iter()
                    .zip(&mut contexts)
                    {
                        assert_eq!(apply(name, input, context, CAPS).unwrap().0, expected);
                    }
                    if family == "cross_cancel" {
                        let mut value =
                            (1u64 << (i % n as usize)) | if i % 2 == 1 { 1u64 << n } else { 0 };
                        for &row in &expected.tail {
                            let pivot = row.trailing_zeros();
                            if value & (1u64 << pivot) != 0 {
                                value ^= row;
                            }
                        }
                        assert_eq!(value, 0);
                    }
                }
            }
        }
    }
    #[test]
    fn caps_count_sources_and_cache_revalidates_smaller_support() {
        let input = Input {
            n: 3,
            degree: 1,
            active: 7,
            polys: vec![vec![0]; 3],
        }
        .canonical();
        let caps = Caps {
            rows: 5,
            columns: 8,
        };
        for name in [
            "flat_cached",
            "sparse_bucket",
            "stream_plain",
            "stream_exchange",
            "stream_census",
            "hybrid_census",
        ] {
            assert!(matches!(
                apply(name, &input, &mut Context::default(), caps),
                Err(Error::RowCap)
            ));
        }
        let input = Input {
            n: 3,
            degree: 2,
            active: 0,
            polys: vec![vec![1, 2, 3]],
        }
        .canonical();
        let mut c = Context::default();
        assert!(
            !apply("flat_cached", &input, &mut c, CAPS)
                .unwrap()
                .1
                .layout_hit
        );
        assert!(
            apply("flat_cached", &input, &mut c, CAPS)
                .unwrap()
                .1
                .layout_hit
        );
        let smaller = Input {
            polys: vec![vec![1]],
            ..input.clone()
        };
        let caps = Caps {
            rows: 1,
            columns: 1,
        };
        assert_eq!(
            apply("flat_cached", &smaller, &mut c, caps).unwrap().0,
            reference(&smaller, caps).unwrap()
        );
        assert!(matches!(
            apply("flat_cached", &input, &mut c, caps),
            Err(Error::ColumnCap)
        ));
    }
    #[test]
    fn cancellation_constant_empty_and_invalid_inputs() {
        for poly in [vec![], vec![0], vec![0, 1], vec![0, 1, 2, 3]] {
            let input = Input {
                n: 3,
                degree: 3,
                active: 7,
                polys: vec![poly],
            }
            .canonical();
            let expected = reference(&input, CAPS).unwrap();
            for name in [
                "flat_cached",
                "sparse_bucket",
                "stream_plain",
                "stream_exchange",
                "stream_census",
                "hybrid_census",
            ] {
                assert_eq!(
                    apply(name, &input, &mut Context::default(), CAPS)
                        .unwrap()
                        .0,
                    expected
                );
            }
        }
        let bad = Input {
            n: 3,
            degree: 3,
            active: 7,
            polys: vec![vec![1, 1]],
        };
        for name in [
            "flat_cached",
            "sparse_bucket",
            "stream_plain",
            "stream_exchange",
            "stream_census",
            "hybrid_census",
        ] {
            assert!(matches!(
                apply(name, &bad, &mut Context::default(), CAPS),
                Err(Error::InvalidInput)
            ));
        }
    }
}
