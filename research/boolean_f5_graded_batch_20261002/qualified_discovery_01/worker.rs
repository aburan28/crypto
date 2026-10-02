//! Experimental exact high-degree Boolean F5 batch cache.
//! Only deterministic generated public quadratic systems are accepted.

#[cfg(test)]
use crypto_lib::cryptanalysis::gf2_elim::echelon_counted;
use crypto_lib::cryptanalysis::gf2_elim::{echelon_prefix_counted, echelon_resume_counted};
use crypto_lib::cryptanalysis::matrix_f5_f2::{
    canonical_row_space_fingerprint, matrix_f5_f2_with_form_timed, F5Criterion, F5OutputForm,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use std::collections::BTreeSet;
use std::mem::size_of;

mod campaign;

#[derive(Clone, Debug, PartialEq, Eq)]
struct System {
    n: u8,
    quadratic: Vec<Vec<u64>>,
    // Bit zero is constant; bit v+1 is the coefficient of variable v.
    affine: Vec<u64>,
}

#[derive(Clone, Debug)]
struct Basis {
    n: u8,
    cols: Vec<u64>,
    offsets: [usize; 5],
    choose: [[usize; 5]; 25],
    high_cols: usize,
}

impl Basis {
    fn new(n: u8) -> Self {
        assert!((1..=24).contains(&n));
        let mut choose = [[0usize; 5]; 25];
        for width in 0..=usize::from(n) {
            choose[width][0] = 1;
            for degree in 1..=4.min(width) {
                choose[width][degree] = if width == 0 {
                    0
                } else {
                    choose[width - 1][degree - 1] + choose[width - 1][degree]
                };
            }
        }
        let mut offsets = [0usize; 5];
        let mut count = 0;
        for degree in (0..=4).rev() {
            offsets[degree] = count;
            count += choose[n as usize][degree];
        }
        let mut cols = monomials_up_to(n, 4);
        cols.sort_unstable_by(|&a, &b| b.count_ones().cmp(&a.count_ones()).then_with(|| a.cmp(&b)));
        assert_eq!(cols.len(), count);
        let result = Self {
            n,
            cols,
            offsets,
            choose,
            high_cols: choose[n as usize][4],
        };
        for (slot, &m) in result.cols.iter().enumerate() {
            assert_eq!(result.slot(m), slot);
        }
        result
    }

    fn slot(&self, mask: u64) -> usize {
        assert!(mask < (1u64 << self.n));
        let degree = mask.count_ones() as usize;
        assert!(degree <= 4);
        let mut index = self.offsets[degree];
        let mut bits = mask;
        let mut position = 1;
        while bits != 0 {
            let bit = bits.trailing_zeros() as usize;
            index += self.choose[bit][position];
            bits &= bits - 1;
            position += 1;
        }
        index
    }

    fn words(&self) -> usize {
        self.cols.len().div_ceil(64)
    }
}

fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn monomials_up_to(n: u8, max_degree: u8) -> Vec<u64> {
    fn choose(out: &mut Vec<u64>, n: u8, left: u8, start: u8, mask: u64) {
        if left == 0 {
            out.push(mask);
            return;
        }
        for bit in start..n {
            choose(out, n, left - 1, bit + 1, mask | (1u64 << bit));
        }
    }
    let mut out = Vec::new();
    for degree in 0..=max_degree.min(n) {
        choose(&mut out, n, degree, 0, 0);
    }
    out.sort_unstable();
    out
}

// Match the inherited Macaulay builder's degree layers and combination
// traversal exactly. Numeric sorting changes row tie breaks in M4RI and can
// return a different, though row-space-equivalent, Echelon basis.
fn source_multipliers(n: u8, max_degree: u8) -> Vec<u64> {
    let variables = (0..n).map(|bit| 1u64 << bit).collect::<Vec<_>>();
    let mut out = vec![0u64];
    let mut level = vec![(0u64, 0usize)];
    for _ in 0..max_degree {
        let mut next = Vec::new();
        for &(monomial, start) in &level {
            for (index, &variable) in variables.iter().enumerate().skip(start) {
                next.push((monomial | variable, index + 1));
            }
        }
        out.extend(next.iter().map(|(mask, _)| *mask));
        level = next;
    }
    out
}

#[cfg(test)]
fn fixture(n: u8, seed: u64) -> System {
    assert!([12, 16, 20, 24].contains(&n));
    let mut state = seed;
    let mut quadratic = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < usize::from(2 * n) {
            let a = next(&mut state) % u64::from(n);
            let b = next(&mut state) % u64::from(n);
            if a != b {
                terms.insert((1u64 << a) | (1u64 << b));
            }
        }
        quadratic.push(terms.into_iter().collect());
    }
    System {
        n,
        quadratic,
        affine: vec![0; n as usize],
    }
}

#[cfg(test)]
fn with_affine(mut system: System, mut seed: u64) -> System {
    for affine in &mut system.affine {
        for slot in 0..=system.n {
            if next(&mut seed) & 1 != 0 {
                *affine |= 1u64 << slot;
            }
        }
    }
    system
}

fn polynomials(system: &System) -> Vec<F2BoolPoly> {
    system
        .quadratic
        .iter()
        .zip(&system.affine)
        .map(|(quadratic, &affine)| {
            let mut terms = quadratic
                .iter()
                .copied()
                .map(F2BoolMono::from_mask)
                .collect::<Vec<_>>();
            if affine & 1 != 0 {
                terms.push(F2BoolMono::from_mask(0));
            }
            for variable in 0..system.n {
                if affine & (1u64 << (variable + 1)) != 0 {
                    terms.push(F2BoolMono::from_mask(1u64 << variable));
                }
            }
            F2BoolPoly::from_monos(terms, system.n as usize)
        })
        .collect()
}

fn baseline(system: &System) -> Vec<F2BoolPoly> {
    matrix_f5_f2_with_form_timed(
        &polynomials(system),
        system.n as usize,
        4,
        F5OutputForm::Echelon,
    )
    .expect("baseline F5")
    .0
}

#[derive(Clone)]
struct Selected {
    labels: Vec<(usize, u64)>,
    rows: Vec<Vec<u64>>,
    full_columns: bool,
    rows_f4: u64,
    criterion_rows: u64,
    criterion_word_xors: u64,
}

fn selected_matrix(system: &System, basis: &Basis) -> Selected {
    let polys = polynomials(system);
    let mask = (1u64 << system.n) - 1;
    let criterion = F5Criterion::new(&polys, system.n as usize, 4, mask);
    let multipliers = source_multipliers(system.n, 2);
    let mut labels = Vec::new();
    let mut rows = Vec::new();
    let mut occupied = vec![0u64; basis.words()];
    let mut rows_f4 = 0u64;
    for (generator, poly) in polys.iter().enumerate() {
        for &multiplier in &multipliers {
            let mut row = vec![0u64; basis.words()];
            for term in &poly.terms {
                let slot = basis.slot(term.mask | multiplier);
                row[slot / 64] ^= 1u64 << (slot % 64);
            }
            if row.iter().all(|&word| word == 0) {
                continue;
            }
            rows_f4 += 1;
            if criterion.prunes(generator, multiplier) {
                continue;
            }
            for (seen, &word) in occupied.iter_mut().zip(&row) {
                *seen |= word;
            }
            labels.push((generator, multiplier));
            rows.push(row);
        }
    }
    let last_bits = basis.cols.len() % 64;
    let last_mask = if last_bits == 0 {
        u64::MAX
    } else {
        (1u64 << last_bits) - 1
    };
    let full_columns = occupied[..occupied.len() - 1]
        .iter()
        .all(|&word| word == u64::MAX)
        && occupied[occupied.len() - 1] == last_mask;
    Selected {
        labels,
        rows,
        full_columns,
        rows_f4,
        criterion_rows: criterion.lower_level_rows().0,
        criterion_word_xors: criterion.word_ops(),
    }
}

#[cfg(test)]
fn bit(row: &[u64], slot: usize) -> bool {
    row[slot / 64] & (1u64 << (slot % 64)) != 0
}

fn xor_words(target: &mut [u64], source: &[u64]) -> u64 {
    for (left, &right) in target.iter_mut().zip(source) {
        *left ^= right;
    }
    target.len() as u64
}

// Multiply the sparse changing suffix by the fixed high-prefix row transform.
// A changing row has few low-degree terms; traversing those terms and XORing
// one packed transform column avoids scanning every transform-row bit for
// every assignment. The result is transposed back into the inherited row
// layout before the exact lower-column M4RI continuation.
fn apply_transform_columns(
    columns: &[Vec<u64>],
    low: &[Vec<u64>],
    low_bits: usize,
) -> (Vec<Vec<u64>>, u64) {
    let count = low.len();
    let output_words = count.div_ceil(64);
    let low_words = low_bits.div_ceil(64);
    assert_eq!(columns.len(), count);
    let mut by_column = vec![vec![0u64; output_words]; low_bits];
    let mut word_xors = 0;
    for (source, row) in low.iter().enumerate() {
        assert_eq!(row.len(), low_words);
        for (block, &word) in row.iter().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let column = block * 64 + remaining.trailing_zeros() as usize;
                assert!(column < low_bits);
                word_xors += xor_words(&mut by_column[column], &columns[source]);
                remaining &= remaining - 1;
            }
        }
    }
    let mut output = vec![vec![0u64; low_words]; count];
    for column_word in 0..low_words {
        for row_start in (0..count).step_by(64) {
            let mut tile = [0u64; 64];
            for bit in 0..64 {
                let column = column_word * 64 + bit;
                if column >= low_bits {
                    break;
                }
                let mut rows = by_column[column][row_start / 64];
                while rows != 0 {
                    let row = rows.trailing_zeros() as usize;
                    tile[row] |= 1u64 << bit;
                    rows &= rows - 1;
                }
            }
            for row in 0..64.min(count - row_start) {
                output[row_start + row][column_word] = tile[row];
            }
        }
    }
    (output, word_xors)
}

#[cfg(test)]
fn high_rank(selected: &Selected, basis: &Basis) -> (usize, u64) {
    let mut rows = selected
        .rows
        .iter()
        .map(|source| {
            let mut high = vec![0u64; basis.high_cols.div_ceil(64)];
            for slot in 0..basis.high_cols {
                if bit(source, slot) {
                    high[slot / 64] |= 1u64 << (slot % 64);
                }
            }
            high
        })
        .collect::<Vec<_>>();
    let mut word_ops = 0;
    let rank = echelon_counted(&mut rows, basis.high_cols, &mut word_ops);
    (rank, word_ops)
}

struct HighCache {
    base: System,
    basis: Basis,
    labels: Vec<(usize, u64)>,
    base_rows: Vec<Vec<u64>>,
    fixed_rows: Vec<Vec<u64>>,
    transform_columns: Vec<Vec<u64>>,
    prefix_word: usize,
    prefix_rank: usize,
    compile_word_xors: u64,
    context_bytes: usize,
}

impl HighCache {
    fn compile(base: &System) -> Result<Self, &'static str> {
        #[cfg(test)]
        let profile_start = std::time::Instant::now();
        let basis = Basis::new(base.n);
        let selected = selected_matrix(base, &basis);
        #[cfg(test)]
        let selected_at = profile_start.elapsed();
        if selected.rows.is_empty() {
            return Err("base-empty-selected-matrix");
        }
        let count = selected.rows.len();
        let full_cols = basis.cols.len();
        let prefix_cols = basis.high_cols / 64 * 64;
        let augmented_cols = full_cols + count;
        let mut augmented = vec![vec![0u64; augmented_cols.div_ceil(64)]; count];
        for (row_index, source) in selected.rows.iter().enumerate() {
            augmented[row_index][..basis.words()].copy_from_slice(source);
            let identity = full_cols + row_index;
            augmented[row_index][identity / 64] |= 1u64 << (identity % 64);
        }
        #[cfg(test)]
        let augmented_at = profile_start.elapsed();
        let mut compile_word_xors = 0;
        let prefix_rank = echelon_prefix_counted(
            &mut augmented,
            prefix_cols,
            augmented_cols,
            &mut compile_word_xors,
        );
        if prefix_rank == 0 {
            return Err("empty-prefix-rank");
        }
        #[cfg(test)]
        let reduced_at = profile_start.elapsed();
        let mut fixed_rows = Vec::with_capacity(count);
        let mut transform_columns = vec![vec![0u64; count.div_ceil(64)]; count];
        for (output_row, source) in augmented.iter().enumerate() {
            let mut fixed = source[..basis.words()].to_vec();
            if full_cols % 64 != 0 {
                *fixed.last_mut().unwrap() &= (1u64 << (full_cols % 64)) - 1;
            }
            for input_word in 0..count.div_ceil(64) {
                let offset = full_cols + input_word * 64;
                let first = offset / 64;
                let shift = offset % 64;
                let mut bits = source[first] >> shift;
                if shift != 0 && first + 1 < source.len() {
                    bits |= source[first + 1] << (64 - shift);
                }
                while bits != 0 {
                    let input_row = input_word * 64 + bits.trailing_zeros() as usize;
                    if input_row < count {
                        transform_columns[input_row][output_row / 64] |= 1u64 << (output_row % 64);
                    }
                    bits &= bits - 1;
                }
            }
            fixed_rows.push(fixed);
        }
        let context_bytes = size_of::<Self>()
            + transform_columns.capacity() * size_of::<Vec<u64>>()
            + selected.rows.capacity() * size_of::<Vec<u64>>()
            + fixed_rows.capacity() * size_of::<Vec<u64>>()
            + base.quadratic.capacity() * size_of::<Vec<u64>>()
            + base
                .quadratic
                .iter()
                .map(|r| r.capacity() * 8)
                .sum::<usize>()
            + base.affine.capacity() * size_of::<u64>()
            + transform_columns
                .iter()
                .map(|r| r.capacity() * 8)
                .sum::<usize>()
            + selected
                .rows
                .iter()
                .map(|r| r.capacity() * 8)
                .sum::<usize>()
            + fixed_rows.iter().map(|r| r.capacity() * 8).sum::<usize>()
            + basis.cols.capacity() * size_of::<u64>()
            + selected.labels.capacity() * size_of::<(usize, u64)>();
        if context_bytes > 128 * 1024 * 1024 {
            return Err("context-cap");
        }
        #[cfg(test)]
        if std::env::var("KIC_DEV_PROFILE").as_deref() == Ok("1") {
            eprintln!(
                "development n={} compile_selected_ns={} compile_augment_ns={} compile_reduce_ns={} compile_extract_ns={} compile_word_xors={}",
                base.n,
                selected_at.as_nanos(),
                (augmented_at - selected_at).as_nanos(),
                (reduced_at - augmented_at).as_nanos(),
                (profile_start.elapsed() - reduced_at).as_nanos(),
                compile_word_xors
            );
        }
        Ok(Self {
            base: base.clone(),
            basis,
            labels: selected.labels,
            base_rows: selected.rows,
            fixed_rows,
            transform_columns,
            prefix_word: prefix_cols / 64,
            prefix_rank,
            compile_word_xors,
            context_bytes,
        })
    }

    fn apply(&self, input: &System) -> (Vec<F2BoolPoly>, bool, &'static str, u64, u64) {
        #[cfg(test)]
        let profile_start = std::time::Instant::now();
        if input.n != self.base.n || input.quadratic != self.base.quadratic {
            return (baseline(input), false, "context-change", 0, 0);
        }
        let selected = selected_matrix(input, &self.basis);
        #[cfg(test)]
        let selected_at = profile_start.elapsed();
        if selected.labels != self.labels || !selected.full_columns {
            return (baseline(input), false, "signature-or-support-change", 0, 0);
        }
        let mut changed_rows = Vec::with_capacity(selected.rows.len());
        for (row_index, current) in selected.rows.iter().enumerate() {
            let base = &self.base_rows[row_index];
            if current[..self.prefix_word] != base[..self.prefix_word] {
                return (baseline(input), false, "high-prefix-changed", 0, 0);
            }
            changed_rows.push(
                current[self.prefix_word..]
                    .iter()
                    .zip(&base[self.prefix_word..])
                    .map(|(&a, &b)| a ^ b)
                    .collect(),
            );
        }
        let (changed, mut word_xors) = apply_transform_columns(
            &self.transform_columns,
            &changed_rows,
            self.basis.cols.len() - self.prefix_word * 64,
        );
        let mut packed_rows = self.fixed_rows.clone();
        for (row, delta) in packed_rows.iter_mut().zip(&changed) {
            word_xors += xor_words(&mut row[self.prefix_word..], delta);
        }
        #[cfg(test)]
        let transform_at = profile_start.elapsed();
        #[cfg(test)]
        {
            if std::env::var("KIC_DEV_ASSERT_PREFIX").as_deref() == Ok("1") {
                let mut fresh_prefix = selected.rows.clone();
                let mut prefix_ops = 0;
                let fresh_rank = echelon_prefix_counted(
                    &mut fresh_prefix,
                    self.prefix_word * 64,
                    self.basis.cols.len(),
                    &mut prefix_ops,
                );
                assert_eq!(fresh_rank, self.prefix_rank);
                if packed_rows != fresh_prefix {
                    let row = packed_rows
                        .iter()
                        .zip(&fresh_prefix)
                        .position(|(a, b)| a != b);
                    panic!("cached prefix differs from fresh prefix at row {row:?}");
                }
            }
        }
        let mut tail_word_xors = 0;
        let rank = echelon_resume_counted(
            &mut packed_rows,
            self.basis.cols.len(),
            self.prefix_word,
            self.prefix_rank,
            &mut tail_word_xors,
        );
        word_xors += tail_word_xors;
        #[cfg(test)]
        let reduce_at = profile_start.elapsed();
        let mut rows = Vec::with_capacity(rank);
        for full in packed_rows.iter().take(rank) {
            // `Basis::cols` already follows the same DegRevLex order as the
            // inherited direct unpacker. Every bit is unique, so sorting and
            // cancellation would only redo work after the packed scan.
            let terms = full.iter().map(|word| word.count_ones() as usize).sum();
            let mut monos = Vec::with_capacity(terms);
            for (word_index, &word) in full.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let slot = word_index * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    monos.push(F2BoolMono::from_mask(self.basis.cols[slot]));
                }
            }
            rows.push(F2BoolPoly {
                terms: monos,
                n_vars: input.n as usize,
            });
        }
        #[cfg(test)]
        if std::env::var("KIC_DEV_PROFILE").as_deref() == Ok("1") {
            eprintln!(
                "development n={} apply_selected_ns={} apply_transform_ns={} apply_resume_ns={} apply_unpack_ns={} transform_word_xors={} resume_word_xors={}",
                input.n,
                selected_at.as_nanos(),
                (transform_at - selected_at).as_nanos(),
                (reduce_at - transform_at).as_nanos(),
                (profile_start.elapsed() - reduce_at).as_nanos(),
                word_xors - tail_word_xors,
                tail_word_xors
            );
        }
        (
            rows,
            true,
            "hit",
            word_xors + selected.criterion_word_xors,
            tail_word_xors,
        )
    }
}

fn main() {
    campaign::main_cli();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combinadic_basis_matches_numeric_column_order() {
        for n in [12, 16, 20, 24] {
            let basis = Basis::new(n);
            assert_eq!(
                basis.high_cols,
                basis.cols.iter().filter(|m| m.count_ones() == 4).count()
            );
        }
    }

    #[test]
    fn fixed_quartic_projection_ignores_affine_tail() {
        let base = fixture(12, 17);
        let changed = with_affine(base.clone(), 29);
        let basis = Basis::new(12);
        let first = selected_matrix(&base, &basis);
        let second = selected_matrix(&changed, &basis);
        assert_eq!(first.labels, second.labels);
        for (a, b) in first.rows.iter().zip(&second.rows) {
            for col in 0..basis.high_cols {
                assert_eq!(bit(a, col), bit(b, col));
            }
        }
    }

    #[test]
    fn candidate_matches_baseline_or_uses_exact_fallback() {
        let mut compiled = 0;
        let mut hits = 0;
        for n in [12, 16, 20, 24] {
            for seed in 17..=20 {
                let base = fixture(n, seed);
                let basis = Basis::new(n);
                let selected = selected_matrix(&base, &basis);
                let (high_rank, _) = high_rank(&selected, &basis);
                eprintln!(
                    "development n={n} seed={seed}: selected={} high_rank={} deficiency={}",
                    selected.rows.len(),
                    high_rank,
                    selected.rows.len() - high_rank
                );
                let changed = with_affine(base.clone(), 29);
                match HighCache::compile(&base) {
                    Ok(cache) => {
                        compiled += 1;
                        let (got, hit, reason, _, _) = cache.apply(&changed);
                        hits += usize::from(hit);
                        let expected = baseline(&changed);
                        if got != expected {
                            let first = got.iter().zip(&expected).position(|(a, b)| a != b);
                            let same_span = canonical_row_space_fingerprint(&got)
                                == canonical_row_space_fingerprint(&expected);
                            panic!("n={n} seed={seed} hit={hit} reason={reason} got_rows={} expected_rows={} first_difference={first:?} same_span={same_span}",
                                got.len(), expected.len());
                        }
                        assert!(cache.context_bytes <= 128 * 1024 * 1024);
                        assert!(cache.compile_word_xors > 0);
                    }
                    Err(reason) => eprintln!("development n={n} seed={seed}: {reason}"),
                }
            }
        }
        assert!(compiled > 0, "no cache compiled on development cores");
        assert!(hits > 0, "no cache application hit on development cores");
    }

    #[test]
    #[ignore = "development timing only; registered seeds must remain untouched"]
    fn development_n24_batch_timing() {
        use std::time::Instant;
        let base = fixture(24, 17);
        let mut inputs = vec![base.clone()];
        inputs.extend((29..60).map(|seed| with_affine(base.clone(), seed)));
        let start = Instant::now();
        let cache = HighCache::compile(&base).expect("development cache");
        let compile_ns = start.elapsed().as_nanos();
        let mut fresh_ns = 0;
        let mut cached_ns = compile_ns;
        let mut hits = 0;
        let mut word_xors = 0;
        for input in &inputs {
            let start = Instant::now();
            let expected = baseline(input);
            fresh_ns += start.elapsed().as_nanos();
            let start = Instant::now();
            let (actual, hit, reason, ops, _) = cache.apply(input);
            cached_ns += start.elapsed().as_nanos();
            hits += usize::from(hit);
            word_xors += ops;
            assert_eq!(actual, expected, "reason={reason}");
        }
        eprintln!("development n24 B32 fresh_ns={fresh_ns} cached_ns={cached_ns} compile_ns={compile_ns} ratio={:.3} hits={hits} word_xors={word_xors}", fresh_ns as f64 / cached_ns as f64);
    }
}
