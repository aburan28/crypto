//! Exact graded reuse of a degree-three Boolean Macaulay high block.
//! The benchmark producer is added only after the algebraic kernel is checked.

use std::collections::{BTreeMap, BTreeSet};
use std::mem::size_of;
use std::time::Instant;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

const PROTOCOL_BYTES: &[u8] = include_bytes!("protocol.json");
const SOURCE_BYTES: &[u8] = include_bytes!("worker.rs");
const VERIFY_BYTES: &[u8] = include_bytes!("verify.rs");
const FROZEN_PROTOCOL_SHA256: &str =
    "7c25cfa898c690c755c8d6a6377dada44314404774f0c30e037de457b8077862";

mod verify;

#[derive(Deserialize)]
struct Protocol {
    variables: Vec<u8>,
    batches: Vec<usize>,
    families: Vec<String>,
    discovery_seeds: Vec<u64>,
    holdout_seeds: Vec<u64>,
    repetitions: usize,
    row_cap: usize,
    column_cap: usize,
    retained_context_cap_bytes: usize,
    campaign_seconds: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct System {
    n: u8,
    quadratic: Vec<Vec<u32>>,
    // Bit 0 is constant; bit v+1 is the coefficient of variable v.
    affine: Vec<u32>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Row {
    high: Vec<u64>,
    low: Vec<u64>,
}

impl Row {
    fn zero(high_words: usize, low_words: usize) -> Self {
        Self {
            high: vec![0; high_words],
            low: vec![0; low_words],
        }
    }
    fn is_zero(&self) -> bool {
        self.high.iter().chain(&self.low).all(|&w| w == 0)
    }
}

#[derive(Clone)]
struct Basis {
    high: Vec<u32>,
    low: Vec<u32>,
    high_slots: BTreeMap<u32, usize>,
    low_slots: BTreeMap<u32, usize>,
}

impl Basis {
    fn new(n: u8) -> Self {
        let mut high = monomials(n, 3)
            .into_iter()
            .filter(|m| m.count_ones() == 3)
            .collect::<Vec<_>>();
        let mut low = monomials(n, 2);
        high.reverse();
        low.reverse();
        let high_slots = high.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let low_slots = low.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        Self {
            high,
            low,
            high_slots,
            low_slots,
        }
    }
    fn high_words(&self) -> usize {
        self.high.len().div_ceil(64)
    }
    fn low_words(&self) -> usize {
        self.low.len().div_ceil(64)
    }
    fn toggle(&self, row: &mut Row, monomial: u32) {
        if monomial.count_ones() == 3 {
            let slot = self.high_slots[&monomial];
            row.high[slot / 64] ^= 1u64 << (slot % 64);
        } else {
            let slot = self.low_slots[&monomial];
            row.low[slot / 64] ^= 1u64 << (slot % 64);
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Operation {
    Swap(usize, usize),
    Xor(usize, usize),
}

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, Serialize)]
struct Work {
    high_word_xors: u64,
    low_word_xors: u64,
    tail_low_word_xors: u64,
    row_swaps: u64,
}

struct Compiled {
    core: Vec<Vec<u32>>,
    n: u8,
    basis: Basis,
    high_echelon: Vec<Vec<u64>>,
    fixed_low_echelon: Vec<Vec<u64>>,
    high_rank: usize,
    schedule: Vec<Operation>,
    setup_work: Work,
}

fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn monomials(n: u8, max_degree: u8) -> Vec<u32> {
    fn choose(out: &mut Vec<u32>, n: u8, left: u8, next_bit: u8, mask: u32) {
        if left == 0 {
            out.push(mask);
            return;
        }
        for bit in next_bit..n {
            choose(out, n, left - 1, bit + 1, mask | (1u32 << bit));
        }
    }
    let mut out = Vec::new();
    for degree in 0..=max_degree.min(n) {
        choose(&mut out, n, degree, 0, 0);
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn fixture_with_state(n: u8, seed: u64) -> (System, u64) {
    assert!([8, 12, 16, 20, 24].contains(&n));
    let mut state = seed;
    let mut quadratic = Vec::new();
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < usize::from(2 * n) {
            let a = (next(&mut state) % u64::from(n)) as u8;
            let b = (next(&mut state) % u64::from(n)) as u8;
            if a != b {
                terms.insert((1u32 << a) | (1u32 << b));
            }
        }
        quadratic.push(terms.into_iter().collect());
    }
    (
        System {
            n,
            quadratic,
            affine: vec![0; n as usize],
        },
        state,
    )
}

#[cfg(test)]
fn fixture(n: u8, seed: u64) -> System {
    fixture_with_state(n, seed).0
}

fn assignments(n: u8, seed: u64, batch: usize, family: &str) -> Vec<System> {
    assert!((1..=64).contains(&batch));
    assert!([
        "independent_affine",
        "walk_affine",
        "repeat",
        "support_escape"
    ]
    .contains(&family));
    let (base, mut state) = fixture_with_state(n, seed);
    let mut current = base.clone();
    let mut out = Vec::with_capacity(batch);
    out.push(base.clone());
    for index in 1..batch {
        let mut system = base.clone();
        match family {
            "repeat" => {}
            "independent_affine" | "support_escape" => {
                for generator in &mut system.affine {
                    for slot in 0..=n {
                        if next(&mut state) & 1 != 0 {
                            *generator |= 1u32 << slot;
                        }
                    }
                }
            }
            "walk_affine" => {
                for generator in &mut current.affine {
                    let slot = (next(&mut state) % u64::from(n + 1)) as u8;
                    *generator ^= 1u32 << slot;
                }
                system.affine.clone_from(&current.affine);
            }
            _ => unreachable!(),
        }
        if family == "support_escape" && index % 4 == 3 {
            system.quadratic[0].remove(0);
        }
        out.push(system);
    }
    out
}

fn validate(system: &System) {
    assert!(system.n <= 24 && system.n >= 2);
    assert_eq!(system.quadratic.len(), system.affine.len());
    assert!(system
        .quadratic
        .iter()
        .all(|poly| !poly.is_empty() && poly.windows(2).all(|w| w[0] < w[1])));
    assert!(system
        .quadratic
        .iter()
        .flatten()
        .all(|m| { m.count_ones() == 2 && *m < (1u32 << system.n) }));
    assert!(system
        .affine
        .iter()
        .all(|mask| *mask < (1u32 << (system.n + 1))));
}

fn row_count(system: &System) -> usize {
    system.quadratic.len() * (usize::from(system.n) + 1)
}

fn toggle_quadratic(row: &mut Row, poly: &[u32], multiplier: u32, basis: &Basis) {
    for &term in poly {
        basis.toggle(row, term | multiplier);
    }
}

fn toggle_affine(row: &mut Row, mask: u32, multiplier: u32, basis: &Basis, n: u8) {
    if mask & 1 != 0 {
        basis.toggle(row, multiplier);
    }
    for variable in 0..n {
        if mask & (1u32 << (variable + 1)) != 0 {
            basis.toggle(row, (1u32 << variable) | multiplier);
        }
    }
}

fn construct_full(system: &System, basis: &Basis) -> Vec<Row> {
    validate(system);
    let multipliers = monomials(system.n, 1);
    let mut rows = Vec::with_capacity(row_count(system));
    for (poly, &affine) in system.quadratic.iter().zip(&system.affine) {
        for &multiplier in &multipliers {
            let mut row = Row::zero(basis.high_words(), basis.low_words());
            toggle_quadratic(&mut row, poly, multiplier, basis);
            toggle_affine(&mut row, affine, multiplier, basis, system.n);
            rows.push(row);
        }
    }
    rows
}

fn binomial(n: u8, k: u8) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 1..=k {
        result = result * usize::from(n - k + i) / usize::from(i);
    }
    result
}

fn rank_fixed_three(mask: u32) -> usize {
    assert_eq!(mask.count_ones(), 3);
    let mut rank = 0;
    let mut position = 0u8;
    for bit_index in 0..32 {
        if mask & (1u32 << bit_index) != 0 {
            position += 1;
            rank += binomial(bit_index as u8, position);
        }
    }
    rank
}

fn rank_up_to_two(mask: u32) -> usize {
    assert!(mask.count_ones() <= 2);
    let mut rank = 0;
    let mut seen = 0u8;
    for bit_index in (0..32).rev() {
        if mask & (1u32 << bit_index) != 0 {
            let remaining = 2 - seen;
            for degree in 0..=remaining {
                rank += binomial(bit_index as u8, degree);
            }
            seen += 1;
        }
    }
    rank
}

fn ranked_slot(basis: &Basis, monomial: u32) -> (bool, usize) {
    if monomial.count_ones() == 3 {
        let slot = basis.high.len() - 1 - rank_fixed_three(monomial);
        debug_assert_eq!(basis.high[slot], monomial);
        (true, slot)
    } else {
        let slot = basis.low.len() - 1 - rank_up_to_two(monomial);
        debug_assert_eq!(basis.low[slot], monomial);
        (false, slot)
    }
}

fn affine_terms(mask: u32, n: u8) -> impl Iterator<Item = u32> {
    std::iter::once(0)
        .chain((0..n).map(|v| 1u32 << v))
        .enumerate()
        .filter_map(move |(slot, term)| ((mask >> slot) & 1 != 0).then_some(term))
}

fn construct_ranked_dense(system: &System, basis: &Basis) -> Vec<Row> {
    validate(system);
    let multipliers = monomials(system.n, 1);
    let mut rows = Vec::with_capacity(row_count(system));
    for (poly, &affine) in system.quadratic.iter().zip(&system.affine) {
        for &multiplier in &multipliers {
            let mut row = Row::zero(basis.high_words(), basis.low_words());
            for product in poly
                .iter()
                .copied()
                .chain(affine_terms(affine, system.n))
                .map(|term| term | multiplier)
            {
                let (high, slot) = ranked_slot(basis, product);
                let target = if high { &mut row.high } else { &mut row.low };
                target[slot / 64] ^= 1u64 << (slot % 64);
            }
            rows.push(row);
        }
    }
    rows
}

fn construct_ranked_sparse(system: &System, basis: &Basis) -> Vec<Row> {
    validate(system);
    let multipliers = monomials(system.n, 1);
    let mut rows = Vec::with_capacity(row_count(system));
    for (poly, &affine) in system.quadratic.iter().zip(&system.affine) {
        for &multiplier in &multipliers {
            let mut sparse = BTreeMap::<(bool, usize), u64>::new();
            for product in poly
                .iter()
                .copied()
                .chain(affine_terms(affine, system.n))
                .map(|term| term | multiplier)
            {
                let (high, slot) = ranked_slot(basis, product);
                *sparse.entry((high, slot / 64)).or_default() ^= 1u64 << (slot % 64);
            }
            let mut row = Row::zero(basis.high_words(), basis.low_words());
            for ((high, word), bits) in sparse {
                if high {
                    row.high[word] = bits;
                } else {
                    row.low[word] = bits;
                }
            }
            rows.push(row);
        }
    }
    rows
}

struct DenseSlots {
    n: u8,
    table: Vec<u16>,
}

impl DenseSlots {
    fn new(basis: &Basis, n: u8) -> Option<Self> {
        let entries = 1usize << n;
        let table_bytes = entries.checked_mul(size_of::<u16>())?;
        let basis_bytes = (basis.high.len() + basis.low.len()) * size_of::<u32>();
        if table_bytes + basis_bytes > 64 * 1024 * 1024 {
            return None;
        }
        assert!(basis.high.len() < 0x8000 && basis.low.len() < 0x8000);
        let mut table = vec![u16::MAX; entries];
        for (slot, &monomial) in basis.high.iter().enumerate() {
            table[monomial as usize] = (slot as u16) | 0x8000;
        }
        for (slot, &monomial) in basis.low.iter().enumerate() {
            table[monomial as usize] = slot as u16;
        }
        Some(Self { n, table })
    }
    fn slot(&self, monomial: u32) -> (bool, usize) {
        assert!(monomial < (1u32 << self.n));
        let encoded = self.table[monomial as usize];
        assert_ne!(encoded, u16::MAX);
        (encoded & 0x8000 != 0, (encoded & 0x7fff) as usize)
    }
    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.table.capacity() * size_of::<u16>()
    }
}

fn construct_dense_lookup(system: &System, basis: &Basis, slots: &DenseSlots) -> Vec<Row> {
    validate(system);
    let multipliers = monomials(system.n, 1);
    let mut rows = Vec::with_capacity(row_count(system));
    for (poly, &affine) in system.quadratic.iter().zip(&system.affine) {
        for &multiplier in &multipliers {
            let mut row = Row::zero(basis.high_words(), basis.low_words());
            for product in poly
                .iter()
                .copied()
                .chain(affine_terms(affine, system.n))
                .map(|term| term | multiplier)
            {
                let (high, slot) = slots.slot(product);
                let target = if high { &mut row.high } else { &mut row.low };
                target[slot / 64] ^= 1u64 << (slot % 64);
            }
            rows.push(row);
        }
    }
    rows
}

fn construct_affine_low(system: &System, basis: &Basis) -> Vec<Vec<u64>> {
    let multipliers = monomials(system.n, 1);
    let mut rows = Vec::with_capacity(row_count(system));
    for &affine in &system.affine {
        for &multiplier in &multipliers {
            let mut row = Row::zero(basis.high_words(), basis.low_words());
            toggle_affine(&mut row, affine, multiplier, basis, system.n);
            debug_assert!(row.high.iter().all(|&word| word == 0));
            rows.push(row.low);
        }
    }
    rows
}

fn bit(row: &[u64], column: usize) -> bool {
    row[column / 64] & (1u64 << (column % 64)) != 0
}

fn xor_words(target: &mut [u64], source: &[u64]) {
    for (left, right) in target.iter_mut().zip(source) {
        *left ^= right;
    }
}

fn xor_row(rows: &mut [Row], source: usize, target: usize, work: &mut Work) {
    assert!(source < target);
    let (prefix, suffix) = rows.split_at_mut(target);
    let pivot = &prefix[source];
    let destination = &mut suffix[0];
    xor_words(&mut destination.high, &pivot.high);
    xor_words(&mut destination.low, &pivot.low);
    work.high_word_xors += pivot.high.len() as u64;
    work.low_word_xors += pivot.low.len() as u64;
}

fn eliminate_full(rows: &mut [Row], basis: &Basis) -> (usize, Work) {
    let mut work = Work::default();
    let mut pivot = 0;
    for column in 0..basis.high.len() {
        let selected = (pivot..rows.len()).find(|&r| bit(&rows[r].high, column));
        if let Some(source) = selected {
            if source != pivot {
                rows.swap(source, pivot);
                work.row_swaps += 1;
            }
            for later in pivot + 1..rows.len() {
                if bit(&rows[later].high, column) {
                    xor_row(rows, pivot, later, &mut work);
                }
            }
            pivot += 1;
        }
    }
    tail_eliminate(rows, pivot, basis.low.len(), &mut work);
    let rank = rows.iter().filter(|row| !row.is_zero()).count();
    (rank, work)
}

fn tail_eliminate(rows: &mut [Row], mut pivot: usize, low_columns: usize, work: &mut Work) {
    for column in 0..low_columns {
        let selected = (pivot..rows.len()).find(|&r| bit(&rows[r].low, column));
        if let Some(source) = selected {
            if source != pivot {
                rows.swap(source, pivot);
                work.row_swaps += 1;
            }
            for later in pivot + 1..rows.len() {
                if bit(&rows[later].low, column) {
                    let (prefix, suffix) = rows.split_at_mut(later);
                    let pivot_row = &prefix[pivot].low;
                    xor_words(&mut suffix[0].low, pivot_row);
                    work.tail_low_word_xors += pivot_row.len() as u64;
                }
            }
            pivot += 1;
        }
    }
}

fn fresh(system: &System, basis: &Basis) -> (Vec<Row>, usize, Work) {
    let mut rows = construct_full(system, basis);
    let (rank, work) = eliminate_full(&mut rows, basis);
    rows.retain(|r| !r.is_zero());
    (rows, rank, work)
}

fn reduce_constructed(mut rows: Vec<Row>, basis: &Basis) -> (Vec<Row>, usize, Work) {
    let (rank, work) = eliminate_full(&mut rows, basis);
    rows.retain(|r| !r.is_zero());
    (rows, rank, work)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Arm {
    PackedDense,
    RankedDense,
    RankedSparse,
    MatrixCache,
    Graded,
}

impl Arm {
    const ALL: [Self; 5] = [
        Self::PackedDense,
        Self::RankedDense,
        Self::RankedSparse,
        Self::MatrixCache,
        Self::Graded,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::PackedDense => "packed_dense",
            Self::RankedDense => "ranked_dense",
            Self::RankedSparse => "ranked_sparse",
            Self::MatrixCache => "matrix_cache",
            Self::Graded => "graded",
        }
    }
}

enum Engine {
    PackedDense(Basis, DenseSlots),
    RankedDense(Basis),
    RankedSparse(Basis),
    MatrixCache(Basis, System, Vec<Row>, usize, Work),
    Graded(Compiled),
}

impl Engine {
    fn compile(arm: Arm, base: &System) -> Option<Self> {
        let basis = Basis::new(base.n);
        match arm {
            Arm::PackedDense => {
                DenseSlots::new(&basis, base.n).map(|slots| Self::PackedDense(basis, slots))
            }
            Arm::RankedDense => Some(Self::RankedDense(basis)),
            Arm::RankedSparse => Some(Self::RankedSparse(basis)),
            Arm::MatrixCache => {
                let (rows, rank, setup_work) =
                    reduce_constructed(construct_ranked_dense(base, &basis), &basis);
                Some(Self::MatrixCache(
                    basis,
                    base.clone(),
                    rows,
                    rank,
                    setup_work,
                ))
            }
            Arm::Graded => Some(Self::Graded(Compiled::new(base))),
        }
    }

    fn apply(&mut self, input: &System) -> (Vec<Row>, usize, Work, bool) {
        match self {
            Self::PackedDense(basis, slots) => {
                let (rows, rank, work) =
                    reduce_constructed(construct_dense_lookup(input, basis, slots), basis);
                (rows, rank, work, true)
            }
            Self::RankedDense(basis) => {
                let (rows, rank, work) =
                    reduce_constructed(construct_ranked_dense(input, basis), basis);
                (rows, rank, work, true)
            }
            Self::RankedSparse(basis) => {
                let (rows, rank, work) =
                    reduce_constructed(construct_ranked_sparse(input, basis), basis);
                (rows, rank, work, true)
            }
            Self::MatrixCache(basis, base, cached, rank, _) => {
                if input == base {
                    (cached.clone(), *rank, Work::default(), true)
                } else {
                    let (rows, rank, work) =
                        reduce_constructed(construct_ranked_dense(input, basis), basis);
                    (rows, rank, work, false)
                }
            }
            Self::Graded(cache) => cache.apply(input),
        }
    }

    fn retained_bytes(&self) -> usize {
        fn basis_bytes(basis: &Basis) -> usize {
            basis.high.capacity() * size_of::<u32>()
                + basis.low.capacity() * size_of::<u32>()
                + (basis.high_slots.len() + basis.low_slots.len())
                    * (size_of::<u32>() + size_of::<usize>())
        }
        fn row_bytes(rows: &Vec<Row>) -> usize {
            rows.iter()
                .map(|r| (r.high.capacity() + r.low.capacity()) * size_of::<u64>())
                .sum::<usize>()
                + rows.capacity() * size_of::<Row>()
        }
        match self {
            Self::PackedDense(basis, slots) => basis_bytes(basis) + slots.retained_bytes(),
            Self::RankedDense(basis) | Self::RankedSparse(basis) => basis_bytes(basis),
            Self::MatrixCache(basis, base, cached, _, _) => {
                basis_bytes(basis)
                    + row_bytes(cached)
                    + base
                        .quadratic
                        .iter()
                        .map(|q| q.capacity() * 4)
                        .sum::<usize>()
                    + base.affine.capacity() * 4
            }
            Self::Graded(cache) => cache.retained_bytes(),
        }
    }

    fn setup_work(&self) -> Work {
        match self {
            Self::MatrixCache(_, _, _, _, work) => *work,
            Self::Graded(cache) => cache.setup_work,
            _ => Work::default(),
        }
    }
}

#[derive(Serialize)]
struct Sample {
    kind: &'static str,
    cell: String,
    repetition: usize,
    order: usize,
    arm: &'static str,
    setup_ns: u128,
    apply_ns: u128,
    validate_ns: u128,
    total_ns: u128,
    retained_bytes: usize,
    setup_work: Work,
    output_bytes: usize,
    output_digest: String,
    hits: usize,
    fallbacks: usize,
    work: Work,
}

struct SampleContext<'a> {
    cell: &'a str,
    inputs: &'a [System],
    expected: &'a [(Vec<Row>, usize)],
    cap: usize,
}

struct CellSpec<'a> {
    n: u8,
    seed: u64,
    batch: usize,
    family: &'a str,
    repetitions: usize,
    cap: usize,
    row_cap: usize,
    column_cap: usize,
    split: &'a str,
}

fn hash_output(hasher: &mut Sha256, rows: &[Row], rank: usize) -> usize {
    hasher.update((rank as u64).to_le_bytes());
    hasher.update((rows.len() as u64).to_le_bytes());
    let mut bytes = 0;
    for row in rows {
        for &word in row.high.iter().chain(&row.low) {
            hasher.update(word.to_le_bytes());
            bytes += 8;
        }
    }
    bytes
}

fn sample(
    arm: Arm,
    name: &'static str,
    repetition: usize,
    order: usize,
    context: &SampleContext<'_>,
) -> Option<Sample> {
    let cell = context.cell;
    let inputs = context.inputs;
    let expected = context.expected;
    let cap = context.cap;
    let start = Instant::now();
    let mut engine = Engine::compile(arm, &inputs[0])?;
    let setup_ns = start.elapsed().as_nanos();
    let retained_bytes = engine.retained_bytes();
    let setup_work = engine.setup_work();
    assert!(retained_bytes <= cap, "retained context cap");
    let mut apply_ns = 0;
    let mut validate_ns = 0;
    let mut output_bytes = 0;
    let mut hits = 0;
    let mut work = Work::default();
    let mut hasher = Sha256::new();
    hasher.update(b"graded-macaulay-output-v1");
    for (input, (reference, expected_rank)) in inputs.iter().zip(expected) {
        let tick = Instant::now();
        let (got, rank, observed_work, hit) = engine.apply(input);
        apply_ns += tick.elapsed().as_nanos();
        hits += usize::from(hit);
        work.high_word_xors += observed_work.high_word_xors;
        work.low_word_xors += observed_work.low_word_xors;
        work.tail_low_word_xors += observed_work.tail_low_word_xors;
        work.row_swaps += observed_work.row_swaps;
        let tick = Instant::now();
        assert_eq!(rank, *expected_rank);
        assert_eq!(&got, reference);
        output_bytes += hash_output(&mut hasher, &got, rank);
        validate_ns += tick.elapsed().as_nanos();
    }
    let tick = Instant::now();
    let output_digest = format!("{:x}", hasher.finalize());
    validate_ns += tick.elapsed().as_nanos();
    drop(engine);
    let total_ns = start.elapsed().as_nanos();
    Some(Sample {
        kind: "sample",
        cell: cell.to_string(),
        repetition,
        order,
        arm: name,
        setup_ns,
        apply_ns,
        validate_ns,
        total_ns,
        retained_bytes,
        setup_work,
        output_bytes,
        output_digest,
        hits,
        fallbacks: inputs.len() - hits,
        work,
    })
}

fn run_cell(spec: &CellSpec<'_>) {
    let n = spec.n;
    let seed = spec.seed;
    let batch = spec.batch;
    let family = spec.family;
    let repetitions = spec.repetitions;
    let cap = spec.cap;
    let row_cap = spec.row_cap;
    let column_cap = spec.column_cap;
    let split = spec.split;
    let inputs = assignments(n, seed, batch, family);
    let basis = Basis::new(n);
    assert!(row_count(&inputs[0]) <= row_cap, "row cap; censored");
    assert!(
        basis.high.len() + basis.low.len() <= column_cap,
        "column cap; censored"
    );
    let mut expected = Vec::with_capacity(batch);
    for input in &inputs {
        let raw = construct_full(input, &basis);
        oracle_check_input(input, &raw, &basis);
        let (rows, rank, _) = reduce_constructed(raw, &basis);
        expected.push((rows, rank));
    }
    let cell = format!("n{n}-{split}-{seed}-{family}-b{batch}");
    println!(
        "{}",
        serde_json::json!({
            "kind":"fixture", "cell":cell, "n":n, "seed":seed,
            "batch":batch, "family":family,
            "output_ranks":expected.iter().map(|(_,r)|*r).collect::<Vec<_>>(),
            "quadratic_support":inputs[0].quadratic,
            "affine_masks":inputs.iter().map(|input|&input.affine).collect::<Vec<_>>(),
            "escape_indices":inputs.iter().enumerate().filter_map(|(i,input)|
                (input.quadratic!=inputs[0].quadratic).then_some(i)).collect::<Vec<_>>()
        })
    );
    let context = SampleContext {
        cell: &cell,
        inputs: &inputs,
        expected: &expected,
        cap,
    };
    for repetition in 0..repetitions {
        for order in 0..2 {
            let aa = sample(
                Arm::RankedDense,
                if order == 0 { "aa_a" } else { "aa_b" },
                repetition,
                order,
                &context,
            )
            .unwrap();
            println!("{}", serde_json::to_string(&aa).unwrap());
        }
        for order in 0..Arm::ALL.len() {
            let arm = Arm::ALL[(repetition + order) % Arm::ALL.len()];
            if let Some(observation) = sample(arm, arm.name(), repetition, order + 2, &context) {
                println!("{}", serde_json::to_string(&observation).unwrap());
            }
        }
    }
}

fn campaign(phase: &str, protocol_path: &str) {
    assert_eq!(
        format!("{:x}", Sha256::digest(PROTOCOL_BYTES)),
        FROZEN_PROTOCOL_SHA256
    );
    let supplied = std::fs::read(protocol_path).expect("read protocol");
    assert_eq!(
        supplied, PROTOCOL_BYTES,
        "protocol differs from compiled source"
    );
    let protocol: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("protocol JSON");
    assert_eq!(protocol.variables, [8, 12, 16, 20, 24]);
    assert_eq!(protocol.batches, [1, 2, 8, 32, 64]);
    assert_eq!(protocol.repetitions, 10);
    assert_eq!(protocol.row_cap, 4096);
    assert_eq!(protocol.column_cap, 8192);
    assert_eq!(protocol.retained_context_cap_bytes, 64 * 1024 * 1024);
    assert_eq!(protocol.campaign_seconds, 1200);
    assert_eq!(
        protocol.families,
        [
            "independent_affine",
            "walk_affine",
            "repeat",
            "support_escape"
        ]
    );
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "full" => ("holdout", protocol.holdout_seeds),
        _ => panic!("expected discovery or full phase"),
    };
    let mut hasher = Sha256::new();
    hasher.update(PROTOCOL_BYTES);
    println!(
        "{}",
        serde_json::json!({
            "kind":"campaign", "phase":phase,
            "protocol_sha256":format!("{:x}",hasher.finalize()),
            "source_sha256":format!("{:x}",Sha256::digest(SOURCE_BYTES)),
            "verifier_sha256":format!("{:x}",Sha256::digest(VERIFY_BYTES))
        })
    );
    let started = Instant::now();
    for n in protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    run_cell(&CellSpec {
                        n,
                        seed,
                        batch,
                        family,
                        repetitions: protocol.repetitions,
                        cap: protocol.retained_context_cap_bytes,
                        row_cap: protocol.row_cap,
                        column_cap: protocol.column_cap,
                        split,
                    });
                    assert!(
                        started.elapsed().as_secs() <= protocol.campaign_seconds,
                        "campaign cap; incomplete run is censored"
                    );
                }
            }
        }
    }
}

impl Compiled {
    fn new(base: &System) -> Self {
        validate(base);
        let basis = Basis::new(base.n);
        let fixed = System {
            n: base.n,
            quadratic: base.quadratic.clone(),
            affine: vec![0; base.affine.len()],
        };
        let rows = construct_full(&fixed, &basis);
        let mut high = rows.iter().map(|r| r.high.clone()).collect::<Vec<_>>();
        let mut low = rows.iter().map(|r| r.low.clone()).collect::<Vec<_>>();
        let mut schedule = Vec::new();
        let mut work = Work::default();
        let mut pivot = 0;
        for column in 0..basis.high.len() {
            let selected = (pivot..high.len()).find(|&r| bit(&high[r], column));
            if let Some(source) = selected {
                if source != pivot {
                    high.swap(pivot, source);
                    low.swap(pivot, source);
                    schedule.push(Operation::Swap(pivot, source));
                    work.row_swaps += 1;
                }
                for later in pivot + 1..high.len() {
                    if bit(&high[later], column) {
                        let (high_prefix, high_suffix) = high.split_at_mut(later);
                        xor_words(&mut high_suffix[0], &high_prefix[pivot]);
                        let (low_prefix, low_suffix) = low.split_at_mut(later);
                        xor_words(&mut low_suffix[0], &low_prefix[pivot]);
                        schedule.push(Operation::Xor(pivot, later));
                        work.high_word_xors += basis.high_words() as u64;
                        work.low_word_xors += basis.low_words() as u64;
                    }
                }
                pivot += 1;
            }
        }
        Self {
            core: base.quadratic.clone(),
            n: base.n,
            basis,
            high_echelon: high,
            fixed_low_echelon: low,
            high_rank: pivot,
            schedule,
            setup_work: work,
        }
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>()
            + self.schedule.capacity() * size_of::<Operation>()
            + self.core.capacity() * size_of::<Vec<u32>>()
            + self.core.iter().map(|r| r.capacity() * 4).sum::<usize>()
            + self.high_echelon.capacity() * size_of::<Vec<u64>>()
            + self
                .high_echelon
                .iter()
                .map(|r| r.capacity() * 8)
                .sum::<usize>()
            + self.fixed_low_echelon.capacity() * size_of::<Vec<u64>>()
            + self
                .fixed_low_echelon
                .iter()
                .map(|r| r.capacity() * 8)
                .sum::<usize>()
            + self.basis.high.capacity() * 4
            + self.basis.low.capacity() * 4
            + (self.basis.high_slots.len() + self.basis.low_slots.len())
                * (size_of::<u32>() + size_of::<usize>())
    }

    fn apply(&self, input: &System) -> (Vec<Row>, usize, Work, bool) {
        validate(input);
        if input.n != self.n || input.quadratic != self.core {
            let (rows, rank, work) = fresh(input, &Basis::new(input.n));
            return (rows, rank, work, false);
        }
        let mut low = construct_affine_low(input, &self.basis);
        let mut work = Work::default();
        for &op in &self.schedule {
            match op {
                Operation::Swap(a, b) => {
                    low.swap(a, b);
                    work.row_swaps += 1;
                }
                Operation::Xor(source, target) => {
                    let (prefix, suffix) = low.split_at_mut(target);
                    xor_words(&mut suffix[0], &prefix[source]);
                    work.low_word_xors += self.basis.low_words() as u64;
                }
            }
        }
        let mut rows = Vec::with_capacity(low.len());
        for ((high, fixed), affine) in self
            .high_echelon
            .iter()
            .zip(&self.fixed_low_echelon)
            .zip(low)
        {
            let mut result_low = fixed.clone();
            xor_words(&mut result_low, &affine);
            work.low_word_xors += self.basis.low_words() as u64;
            rows.push(Row {
                high: high.clone(),
                low: result_low,
            });
        }
        tail_eliminate(&mut rows, self.high_rank, self.basis.low.len(), &mut work);
        rows.retain(|r| !r.is_zero());
        let rank = rows.len();
        (rows, rank, work, true)
    }
}

fn independent_product_row(system: &System, generator: usize, multiplier: u32) -> BTreeSet<u32> {
    let mut parity = BTreeSet::new();
    let toggle = |set: &mut BTreeSet<u32>, monomial| {
        if !set.insert(monomial) {
            set.remove(&monomial);
        }
    };
    for &term in &system.quadratic[generator] {
        toggle(&mut parity, term | multiplier);
    }
    let affine = system.affine[generator];
    if affine & 1 != 0 {
        toggle(&mut parity, multiplier);
    }
    for variable in 0..system.n {
        if affine & (1u32 << (variable + 1)) != 0 {
            toggle(&mut parity, (1u32 << variable) | multiplier);
        }
    }
    parity
}

fn oracle_check_input(system: &System, rows: &[Row], basis: &Basis) {
    let multipliers = monomials(system.n, 1);
    assert_eq!(rows.len(), row_count(system));
    for (index, row) in rows.iter().enumerate() {
        let generator = index / multipliers.len();
        let multiplier = multipliers[index % multipliers.len()];
        let expected = independent_product_row(system, generator, multiplier);
        let actual = basis
            .high
            .iter()
            .enumerate()
            .filter(|(i, _)| bit(&row.high, *i))
            .chain(
                basis
                    .low
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| bit(&row.low, *i)),
            )
            .map(|(_, &m)| m)
            .collect::<BTreeSet<_>>();
        assert_eq!(actual, expected);
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() == 3 && args[1] == "--seal" {
        verify::seal(&args[2]);
        return;
    }
    if args.len() == 3 && args[1] == "--verify-bundle" {
        verify::verify_bundle(&args[2]);
        return;
    }
    if args.len() == 3 && args[1] == "--wait-quiet" {
        verify::wait_quiet(&args[2]);
        return;
    }
    if args.len() == 4 && args[1] == "--failure" {
        verify::failure(&args[2], &args[3]);
        return;
    }
    if args.len() == 4 && args[1] == "--check-discovery" {
        verify::check_discovery(&args[2], &args[3]);
        return;
    }
    if args.len() == 6 && args[1] == "--verify" {
        verify::verify_run(&args[2], &args[3], &args[4], &args[5]);
        return;
    }
    if args.len() == 4 && args[1] == "--campaign" {
        campaign(&args[2], &args[3]);
        return;
    }
    panic!("usage: --campaign PHASE PROTOCOL | --verify PHASE RAW CONDITIONS RESULTS | --check-discovery BUNDLE BINDING | --wait-quiet OUTPUT | --failure BUNDLE REASON | --seal BUNDLE | --verify-bundle BUNDLE");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn with_affine(mut system: System, mut bits: u64) -> System {
        for a in &mut system.affine {
            bits = next(&mut bits);
            *a = (bits as u32) & ((1u32 << (system.n + 1)) - 1);
        }
        system
    }

    #[test]
    fn high_block_is_invariant_and_full_echelon_matches() {
        for n in [8, 12, 16, 20, 24] {
            let base = fixture(n, 17);
            let basis = Basis::new(n);
            let compiled = Compiled::new(&base);
            assert!(compiled.retained_bytes() < 64 * 1024 * 1024);
            assert!(compiled.setup_work.high_word_xors > 0);
            for seed in [1, 19, 97] {
                let system = with_affine(base.clone(), seed);
                let original = construct_full(&system, &basis);
                oracle_check_input(&system, &original, &basis);
                let fixed = construct_full(&base, &basis);
                for (a, b) in original.iter().zip(&fixed) {
                    assert_eq!(a.high, b.high);
                }
                let (expected, rank, _) = fresh(&system, &basis);
                let (actual, actual_rank, _, hit) = compiled.apply(&system);
                assert!(hit);
                assert_eq!(rank, actual_rank);
                assert_eq!(expected, actual);
            }
        }
    }

    #[test]
    fn quadratic_escape_uses_full_fallback() {
        let base = fixture(8, 23);
        let compiled = Compiled::new(&base);
        let mut changed = with_affine(base.clone(), 11);
        changed.quadratic[0].remove(0);
        let expected = fresh(&changed, &Basis::new(8));
        let actual = compiled.apply(&changed);
        assert!(!actual.3);
        assert_eq!(actual.0, expected.0);
        assert_eq!(actual.1, expected.1);
    }

    #[test]
    fn boolean_low_tail_cancellations_are_exact() {
        let base = fixture(8, 29);
        let compiled = Compiled::new(&base);
        let mut system = base.clone();
        for affine in &mut system.affine {
            *affine = (1u32 << 9) - 1;
        }
        let (expected, rank, _) = fresh(&system, &Basis::new(8));
        let (actual, actual_rank, _, hit) = compiled.apply(&system);
        assert!(hit);
        assert_eq!(rank, actual_rank);
        assert_eq!(expected, actual);
    }

    #[test]
    fn exhaustive_small_affine_assignments_preserve_exact_echelon() {
        let base = System {
            n: 4,
            quadratic: vec![vec![0b0011, 0b0101], vec![0b0110, 0b1001]],
            affine: vec![0, 0],
        };
        let basis = Basis::new(4);
        let compiled = Compiled::new(&base);
        for left in 0..32 {
            for right in 0..32 {
                let mut system = base.clone();
                system.affine = vec![left, right];
                let input = construct_full(&system, &basis);
                oracle_check_input(&system, &input, &basis);
                let (expected, rank, _) = fresh(&system, &basis);
                let (actual, actual_rank, _, hit) = compiled.apply(&system);
                assert!(hit);
                assert_eq!(actual_rank, rank);
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn ranked_and_dense_slots_match_every_ambient_column() {
        for n in [8, 12, 16, 20, 24] {
            let basis = Basis::new(n);
            let dense = DenseSlots::new(&basis, n);
            assert!(dense.as_ref().unwrap().retained_bytes() < 64 * 1024 * 1024);
            for (slot, &monomial) in basis.high.iter().enumerate() {
                assert_eq!(ranked_slot(&basis, monomial), (true, slot));
                if let Some(dense) = &dense {
                    assert_eq!(dense.slot(monomial), (true, slot));
                }
            }
            for (slot, &monomial) in basis.low.iter().enumerate() {
                assert_eq!(ranked_slot(&basis, monomial), (false, slot));
                if let Some(dense) = &dense {
                    assert_eq!(dense.slot(monomial), (false, slot));
                }
            }
        }
    }

    #[test]
    fn every_declared_family_matches_all_fresh_constructors() {
        for n in [8, 12, 16, 20, 24] {
            let basis = Basis::new(n);
            let dense = DenseSlots::new(&basis, n);
            let base = fixture(n, 41);
            let compiled = Compiled::new(&base);
            for family in [
                "independent_affine",
                "walk_affine",
                "repeat",
                "support_escape",
            ] {
                let batch = assignments(n, 41, 16, family);
                assert_eq!(batch.len(), 16);
                for (index, system) in batch.iter().enumerate() {
                    let expected_input = construct_full(system, &basis);
                    oracle_check_input(system, &expected_input, &basis);
                    assert_eq!(construct_ranked_dense(system, &basis), expected_input);
                    assert_eq!(construct_ranked_sparse(system, &basis), expected_input);
                    if let Some(slots) = &dense {
                        assert_eq!(
                            construct_dense_lookup(system, &basis, slots),
                            expected_input
                        );
                    }
                    let (expected, rank, _) = fresh(system, &basis);
                    let (actual, actual_rank, _, hit) = compiled.apply(system);
                    assert_eq!(actual_rank, rank);
                    assert_eq!(actual, expected);
                    assert_eq!(hit, family != "support_escape" || index % 4 != 3);
                }
            }
        }
    }

    #[test]
    fn campaign_arms_share_one_exact_digest_and_guarded_hits() {
        for n in [8, 24] {
            let inputs = assignments(n, 41, 8, "support_escape");
            let basis = Basis::new(n);
            let expected = inputs
                .iter()
                .map(|input| {
                    let raw = construct_full(input, &basis);
                    oracle_check_input(input, &raw, &basis);
                    let (rows, rank, _) = reduce_constructed(raw, &basis);
                    (rows, rank)
                })
                .collect::<Vec<_>>();
            let mut digests = Vec::new();
            let context = SampleContext {
                cell: "development",
                inputs: &inputs,
                expected: &expected,
                cap: 64 * 1024 * 1024,
            };
            for arm in Arm::ALL {
                let got = sample(arm, arm.name(), 0, 0, &context);
                let got = got.unwrap();
                assert_eq!(
                    got.fallbacks,
                    if arm == Arm::Graded {
                        2
                    } else if arm == Arm::MatrixCache {
                        7
                    } else {
                        0
                    }
                );
                digests.push(got.output_digest);
            }
            assert!(digests.windows(2).all(|pair| pair[0] == pair[1]));
        }
    }
}
