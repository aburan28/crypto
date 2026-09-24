//! # F4 over the boolean ring `F_2[v_0, …, v_{n−1}] / (v_i² − v_i)`.
//!
//! Faugère's F4 (1999) is Buchberger's algorithm with the reductions
//! batched: at each step it takes **every** critical pair of the lowest
//! degree, writes the two halves of each S-polynomial as rows of one
//! matrix, adds a reducer row for every monomial some basis element's
//! leading monomial divides (*symbolic preprocessing*), and reduces the
//! whole matrix at once.  Rows whose pivot is a monomial no basis
//! element leads join the basis.  Over `F_2` the matrix is bit-packed,
//! sixty-four columns to a word, so a reduction step is a run of 64-bit
//! XORs — the unit this module counts.
//!
//! [`pq_groebner_f2`](crate::cryptanalysis::pq_groebner_f2) reduces the
//! same pairs one at a time on sorted monomial lists.  This module is the
//! batched version of the same computation, and a test holds it to that:
//! wherever the Buchberger output is certified to be a Gröbner basis of
//! the boolean ideal, the two reduced bases are identical.
//!
//! ## The field pairs, which the boolean ring needs
//!
//! In `F_2[v]/(v² − v)` a Gröbner basis must also be closed under the
//! products `v_i · g` for every variable `v_i` of `LM(g)`: those are the
//! S-polynomials of `g` against the field equations `v_i² + v_i`, and
//! multiplying by a variable already in the leading monomial is not
//! order-compatible there — `v_0·(v_0v_1 + v_0 + v_1) = v_0`, a smaller
//! leading monomial.  An engine that forms only the S-polynomials of its
//! own basis can stop early: on `{v_0v_1 + v_0 + v_1}` it returns the
//! generator itself, three standard monomials for an ideal with one
//! solution.  This engine queues a **field pair** `(g, v_i)` for every
//! basis element and every variable of its leading monomial, at sugar
//! degree `deg LM(g) + 1`, and reduces it in the same matrices as the
//! critical pairs.
//!
//! A field product is a single row, not half of a pair, so Faugère's
//! "the pivot is not a leading monomial of an input row" test for new
//! elements does not apply to it.  The test used here is the one that
//! holds for both: a reduced row is **new** iff its pivot is divisible by
//! no active basis leading monomial.  Symbolic preprocessing makes that
//! exact, because every such divisible monomial other than an S-pair's
//! common leading monomial receives a reducer row led by it.
//!
//! ## Pair criteria
//!
//! Critical pairs go through the Gebauer–Möller installation (Becker and
//! Weispfenning's `UPDATE`): the chain criterion on the new pairs and on
//! the queued ones, and Buchberger's product criterion.  These hold for
//! the ideal `I + ⟨v_i² + v_i⟩` of the polynomial ring, whose basis
//! elements other than the field equations are exactly the boolean
//! polynomials here, with square-free leading monomials; the field pairs
//! are never pruned and never used to prune.
//!
//! ## What is counted
//!
//! `word_xors`: 64-bit XORs in the eliminations — the forward pass of
//! every step and of the initial echelon, M4RI combination-table formation,
//! and the backward pass of the final inter-reduction — counting only the
//! words actually XORed. The
//! construction of the matrices (the products, the column map, the
//! packing) is **not** in that count; its wall time is reported beside
//! it as `build_ns`, so a reader can see what share of the run the unit
//! covers.  Solution extraction counts its monomial tests separately.
//!
//! ## Solutions, without enumeration
//!
//! The ideal of a boolean system is radical and zero-dimensional, so its
//! reduced basis has exactly as many standard monomials as the system
//! has solutions.  A variable is either the leading monomial of a linear
//! basis element or itself standard — *free*.  Every linear element
//! reads `v_p + (free variables) + c`, so an assignment of the free
//! variables determines the point, and there are at most
//! `#solutions − 1` of them.  [`solutions_from_reduced_basis`]
//! enumerates those assignments and keeps the points every basis element
//! vanishes at; the cost is `2^{#free}`, not `2^n`.
//!
//! ## References
//!
//! - J.-C. Faugère, *A new efficient algorithm for computing Gröbner
//!   bases (F4)*, J. Pure Appl. Algebra 139 (1999).
//! - T. Becker, V. Weispfenning, *Gröbner Bases*, Springer 1993 —
//!   `UPDATE`, p. 230.
//! - M. Brickenstein, *Boolean Gröbner bases*, PhD thesis, TU
//!   Kaiserslautern 2010 — the boolean ring and its field pairs.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cryptanalysis::pq_groebner_f2::{cmp_mono, F2BoolMono, F2BoolPoly};

/// What one F4 run cost, and how far it got.
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct F4Stats {
    /// F4 steps: one matrix each, after the initial echelon.
    pub steps: u64,
    /// Critical pairs reduced, and field pairs reduced.
    pub pairs_reduced: u64,
    pub field_pairs_reduced: u64,
    /// Pairs Buchberger's product criterion dropped.
    pub pairs_product_skipped: u64,
    /// Pairs the Gebauer–Möller chain criterion dropped.
    pub pairs_chain_skipped: u64,
    /// Pending critical-pair survival tests and full queue passes made while
    /// installing new basis elements.
    pub pair_prune_tests: u64,
    pub pair_prune_submask_lookups: u64,
    pub pair_prune_linear_tests: u64,
    pub pair_prune_passes: u64,
    /// Multi-element installation groups and elements handled by them.
    pub batch_insert_groups: u64,
    pub batch_insert_elements: u64,
    /// Active basis entries visited while generating new pairs and removing
    /// leaders dominated by a newly installed element.
    pub active_candidate_visits: u64,
    pub active_deactivation_tests: u64,
    /// Pair-selection backend use and exact low-degree cover probes.
    pub pair_dense_select_calls: u64,
    pub pair_sorted_select_calls: u64,
    pub pair_lcm_groups: u64,
    pub pair_cover_lookups: u64,
    /// Exact heap bytes reserved by the dense pair-selection arrays.
    pub pair_dense_scratch_bytes_max: u64,
    /// Matrices packed through a dense monomial-to-column index and its
    /// largest exact allocation.
    pub dense_column_matrices: u64,
    pub dense_column_bytes_max: u64,
    /// Rows symbolic preprocessing added as reducers.
    pub reducer_rows: u64,
    /// Largest matrix built, and the sum of rows over every matrix.
    pub matrix_rows_max: u64,
    pub matrix_cols_max: u64,
    pub matrix_rows_sum: u64,
    /// The unit: 64-bit XORs performed by elimination and M4RI tables.
    pub word_xors: u64,
    /// Matrices routed through block-4 Method of Four Russians elimination.
    pub m4ri_matrices: u64,
    /// Largest M4RI block width used in this solve.
    pub m4ri_block_width_max: u64,
    /// XORs used to construct M4RI combination tables (also in `word_xors`).
    pub m4ri_table_word_xors: u64,
    pub m4ri_blocks: u64,
    pub m4ri_consecutive_blocks: u64,
    pub m4ri_trimmed_word_xors_avoided: u64,
    /// Largest reusable M4RI pivot/table scratch allocation.
    pub m4ri_scratch_bytes_max: u64,
    /// Matrices packed directly into one contiguous row arena for M4RI and
    /// the largest row-order allocation used by that representation.
    pub flat_m4ri_matrices: u64,
    pub flat_m4ri_order_bytes_max: u64,
    /// Exact divisor candidates tested by symbolic preprocessing, either one
    /// active leading-monomial word test or one indexed submask lookup;
    /// reported, not in `word_xors`.
    pub divisor_tests: u64,
    pub divisor_submask_lookups: u64,
    pub divisor_linear_tests: u64,
    /// Dense boolean-domain multiplication by a monomial.  For small
    /// systems this toggles equal output masks in an epoch-stamped array,
    /// then sorts only the surviving terms instead of sorting every input
    /// term before duplicate cancellation.
    pub dense_mul_calls: u64,
    pub dense_mul_input_terms: u64,
    pub dense_mul_output_terms: u64,
    pub dense_mul_cancelled_terms: u64,
    pub dense_mul_scratch_bytes_max: u64,
    /// Symbolic-preprocessing monomial sets backed by the complete Boolean
    /// mask domain rather than hash tables.
    pub dense_symbolic_set_steps: u64,
    pub dense_symbolic_set_bytes_max: u64,
    /// Elements added to the basis after the initial echelon.
    pub new_elements: u64,
    /// Highest step degree processed.
    pub degree_reached: u32,
    /// Highest step degree at which the basis gained an element: the
    /// solving degree, comparable to `GbStats::solving_degree`.
    pub solving_degree: u32,
    /// Highest degree of any row in any matrix.
    pub max_poly_degree: u32,
    /// Elements of the reduced basis returned.
    pub basis_len: u64,
    /// Largest matrix held, in bytes of packed words.
    pub peak_matrix_bytes: u64,
    /// Wall time spent building matrices (products, columns, packing)
    /// and eliminating them.
    pub build_ns: u64,
    pub eliminate_ns: u64,
    pub pair_update_ns: u64,
    pub wall_ns: u64,
    /// The budget ran out; the returned set is not a Gröbner basis.
    pub timed_out: bool,
    /// A matrix would have exceeded [`MAX_MATRIX_WORDS`].
    pub oversize: bool,
    /// Conservative peak for the sparse symbolic-preprocessing objects
    /// held before a packed matrix exists.
    pub symbolic_bytes_estimate_max: u64,
    /// Symbolic preprocessing crossed [`MAX_SYMBOLIC_BYTES`] and stopped
    /// before allocating the packed matrix.
    pub symbolic_cap_hit: bool,
    pub pairs_left: u64,
}

/// The largest matrix a step may build, in 64-bit words (1 GiB).
pub const MAX_MATRIX_WORDS: u64 = 1 << 27;

/// The symbolic row and monomial sets can be materially larger than the
/// packed matrix they are preparing.  Stop at the same one-GiB envelope
/// before allocator overhead turns a nominally bounded F4 call into a
/// multi-gigabyte process.  The estimate below deliberately overprices hash
/// entries and sparse terms; reaching it is a censored resource terminal.
pub const MAX_SYMBOLIC_BYTES: u64 = 1 << 30;

fn symbolic_bytes_estimate(
    examined: usize,
    queued: usize,
    row_terms: usize,
    rows: usize,
    seen_rows: usize,
) -> u64 {
    (examined as u64)
        .saturating_mul(64)
        .saturating_add((queued as u64).saturating_mul(8))
        .saturating_add((row_terms as u64).saturating_mul(16))
        .saturating_add((rows as u64).saturating_mul(32))
        .saturating_add((seen_rows as u64).saturating_mul(32))
}

fn symbolic_cap_exceeded(
    st: &mut F4Stats,
    examined: usize,
    queued: usize,
    row_terms: usize,
    rows: usize,
    seen_rows: usize,
) -> bool {
    let estimate = symbolic_bytes_estimate(examined, queued, row_terms, rows, seen_rows);
    st.symbolic_bytes_estimate_max = st.symbolic_bytes_estimate_max.max(estimate);
    if estimate > MAX_SYMBOLIC_BYTES {
        st.oversize = true;
        st.symbolic_cap_hit = true;
        true
    } else {
        false
    }
}

const NONE: u32 = u32::MAX;

/// Deterministic SplitMix-style hashing for trusted internal monomial masks.
/// Hash-map key equality still resolves collisions exactly; this only avoids
/// SipHash's adversarial-input cost on masks constructed inside the solver.
#[derive(Default)]
struct FastU64Hasher(u64);

impl Hasher for FastU64Hasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut value = 0xcbf29ce484222325u64;
        for &byte in bytes {
            value = (value ^ u64::from(byte)).wrapping_mul(0x100000001b3);
        }
        self.write_u64(value);
    }

    fn write_u64(&mut self, value: u64) {
        let mut mixed = value.wrapping_add(0x9e3779b97f4a7c15);
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d049bb133111eb);
        self.0 = mixed ^ (mixed >> 31);
    }
}

type FastU64Map<V> = HashMap<u64, V, BuildHasherDefault<FastU64Hasher>>;
type FastU64Set = HashSet<u64, BuildHasherDefault<FastU64Hasher>>;

enum MonomialMembership {
    Dense(Vec<u64>),
    Sparse(FastU64Set),
}

struct MonomialSet {
    masks: Vec<u64>,
    membership: MonomialMembership,
}

impl MonomialSet {
    fn new(n_vars: usize) -> Self {
        let dense =
            n_vars <= 20 && std::env::var("PQ_F4_DISABLE_DENSE_SYMBOLIC_SET").as_deref() != Ok("1");
        Self {
            masks: Vec::new(),
            membership: if dense {
                MonomialMembership::Dense(vec![0; (1usize << n_vars).div_ceil(64)])
            } else {
                MonomialMembership::Sparse(FastU64Set::default())
            },
        }
    }

    fn clear(&mut self) {
        self.masks.clear();
        match &mut self.membership {
            MonomialMembership::Dense(bits) => bits.fill(0),
            MonomialMembership::Sparse(set) => set.clear(),
        }
    }

    fn insert(&mut self, mask: u64) -> bool {
        let inserted = match &mut self.membership {
            MonomialMembership::Dense(bits) => {
                let (word, bit) = (mask as usize / 64, 1u64 << (mask % 64));
                let inserted = bits[word] & bit == 0;
                bits[word] |= bit;
                inserted
            }
            MonomialMembership::Sparse(set) => set.insert(mask),
        };
        if inserted {
            self.masks.push(mask);
        }
        inserted
    }

    fn contains(&self, mask: &u64) -> bool {
        match &self.membership {
            MonomialMembership::Dense(bits) => {
                bits[*mask as usize / 64] & (1u64 << (*mask % 64)) != 0
            }
            MonomialMembership::Sparse(set) => set.contains(mask),
        }
    }

    fn copy_from(&mut self, source: &Self) {
        self.clear();
        self.masks.reserve(source.len());
        for &mask in &source.masks {
            self.insert(mask);
        }
    }

    fn len(&self) -> usize {
        self.masks.len()
    }

    fn iter(&self) -> impl Iterator<Item = &u64> {
        self.masks.iter()
    }

    fn dense_bytes(&self) -> Option<u64> {
        let MonomialMembership::Dense(bits) = &self.membership else {
            return None;
        };
        Some(
            (bits.capacity() * std::mem::size_of::<u64>()
                + self.masks.capacity() * std::mem::size_of::<u64>()) as u64,
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PairKind {
    /// The S-polynomial of basis elements `i` and `j`.
    Critical(usize, usize),
    /// `v · g_i` for a variable `v` of `LM(g_i)`: the S-polynomial of
    /// `g_i` against the field equation `v² + v`.
    Field(usize, u32),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Pair {
    kind: PairKind,
    /// The lcm of the two leading monomials (critical pairs only).
    lcm: u64,
    /// Selection degree: `deg lcm`, or `deg LM(g) + 1` for a field pair.
    deg: u32,
}

/// A bit-packed row.  Words outside `[start, end)` are zero.
struct Row {
    bits: Vec<u64>,
    start: usize,
    end: usize,
}

impl Row {
    fn lead(&mut self) -> Option<usize> {
        while self.start < self.end {
            let w = self.bits[self.start];
            if w != 0 {
                return Some(self.start * 64 + w.trailing_zeros() as usize);
            }
            self.start += 1;
        }
        None
    }
}

struct FlatRows {
    data: Vec<u64>,
    order: Vec<u32>,
    rows: usize,
    words: usize,
}

impl FlatRows {
    fn new(rows: usize, words: usize) -> Self {
        assert!(rows <= u32::MAX as usize);
        let entries = rows
            .checked_mul(words)
            .expect("flat F4 matrix size fits usize");
        Self {
            data: vec![0; entries],
            order: (0..rows as u32).collect(),
            rows,
            words,
        }
    }

    #[inline(always)]
    fn row(&self, logical: usize) -> &[u64] {
        let physical = self.order[logical] as usize;
        &self.data[physical * self.words..(physical + 1) * self.words]
    }

    #[inline(always)]
    fn row_mut(&mut self, logical: usize) -> &mut [u64] {
        let physical = self.order[logical] as usize;
        &mut self.data[physical * self.words..(physical + 1) * self.words]
    }

    fn swap(&mut self, left: usize, right: usize) {
        self.order.swap(left, right);
    }

    fn order_bytes(&self) -> u64 {
        (self.order.capacity() * std::mem::size_of::<u32>()) as u64
    }
}

enum PackedRows {
    Nested(Vec<Row>),
    Flat(FlatRows),
}

enum EchelonRows {
    Nested(Vec<Row>),
    Flat(FlatRows),
}

/// The columns of one matrix: monomials in descending order, so the
/// first set bit of a row is its leading monomial.
struct Columns {
    monos: Vec<u64>,
    index: FastU64Map<usize>,
    dense_index: Option<Vec<u32>>,
}

fn cmp_mono_mask_descending(a: u64, b: u64) -> std::cmp::Ordering {
    b.count_ones().cmp(&a.count_ones()).then_with(|| a.cmp(&b))
}

impl Columns {
    fn from_monomials(set: impl IntoIterator<Item = u64>, n_vars: usize) -> Self {
        let mut monos: Vec<u64> = set.into_iter().collect();
        monos.sort_unstable_by(|&a, &b| cmp_mono_mask_descending(a, b));
        monos.dedup();
        if n_vars <= 20 {
            let mut dense_index = vec![NONE; 1usize << n_vars];
            for (index, &monomial) in monos.iter().enumerate() {
                debug_assert!(index < u32::MAX as usize);
                dense_index[monomial as usize] = index as u32;
            }
            Self {
                monos,
                index: FastU64Map::default(),
                dense_index: Some(dense_index),
            }
        } else {
            let index = monos.iter().enumerate().map(|(i, &m)| (m, i)).collect();
            Self {
                monos,
                index,
                dense_index: None,
            }
        }
    }

    fn words(&self) -> usize {
        self.monos.len().div_ceil(64).max(1)
    }

    fn record_dense_cost(&self, st: &mut F4Stats) {
        if let Some(index) = &self.dense_index {
            st.dense_column_matrices += 1;
            st.dense_column_bytes_max = st
                .dense_column_bytes_max
                .max((index.capacity() * std::mem::size_of::<u32>()) as u64);
        }
    }

    fn pack(&self, p: &F2BoolPoly) -> Row {
        let mut bits = vec![0u64; self.words()];
        let (mut start, mut end) = (usize::MAX, 0usize);
        for t in &p.terms {
            let c = self.dense_index.as_ref().map_or_else(
                || self.index[&t.mask],
                |index| {
                    let column = index[t.mask as usize];
                    debug_assert_ne!(column, NONE);
                    column as usize
                },
            );
            bits[c / 64] ^= 1u64 << (c % 64);
            start = start.min(c / 64);
            end = end.max(c / 64 + 1);
        }
        if start == usize::MAX {
            start = 0;
        }
        Row { bits, start, end }
    }

    fn pack_into(&self, p: &F2BoolPoly, bits: &mut [u64]) {
        debug_assert_eq!(bits.len(), self.words());
        for t in &p.terms {
            let c = self.dense_index.as_ref().map_or_else(
                || self.index[&t.mask],
                |index| {
                    let column = index[t.mask as usize];
                    debug_assert_ne!(column, NONE);
                    column as usize
                },
            );
            bits[c / 64] ^= 1u64 << (c % 64);
        }
    }

    fn unpack(&self, row: &Row, n_vars: usize) -> F2BoolPoly {
        self.unpack_bits(&row.bits, row.start, row.end, n_vars)
    }

    fn unpack_bits(&self, bits: &[u64], start: usize, end: usize, n_vars: usize) -> F2BoolPoly {
        let mut monos = Vec::new();
        for (w, &packed) in bits.iter().enumerate().take(end).skip(start) {
            let mut word = packed;
            while word != 0 {
                let c = w * 64 + word.trailing_zeros() as usize;
                word &= word - 1;
                monos.push(F2BoolMono::from_mask(self.monos[c]));
            }
        }
        F2BoolPoly::from_monos(monos, n_vars)
    }
}

struct DenseMulScratch {
    generation: u32,
    n_vars: u32,
    mask_bits: u32,
    state: Vec<u32>,
    touched: Vec<u32>,
}

impl DenseMulScratch {
    fn new(n_vars: usize) -> Self {
        debug_assert!(n_vars <= 20);
        let size = 1usize << n_vars;
        Self {
            generation: 0,
            n_vars: n_vars as u32,
            mask_bits: (size - 1) as u32,
            state: vec![0; size],
            touched: Vec::new(),
        }
    }

    fn bytes(&self) -> u64 {
        (self.state.capacity() * std::mem::size_of::<u32>()
            + self.touched.capacity() * std::mem::size_of::<u32>()) as u64
    }

    fn multiply(&mut self, p: &F2BoolPoly, multiplier: u64) -> F2BoolPoly {
        debug_assert_eq!(p.n_vars, self.n_vars as usize);
        debug_assert_eq!(multiplier & !u64::from(self.mask_bits), 0);
        self.generation = self.generation.wrapping_add(2);
        if self.generation == 0 {
            self.state.fill(0);
            self.generation = 2;
        }
        self.touched.clear();
        self.touched.reserve(p.terms.len());
        for term in &p.terms {
            let mask = (term.mask | multiplier) as u32;
            let index = mask as usize;
            if self.state[index] & !1 != self.generation {
                self.state[index] = self.generation | 1;
                let order_key = ((self.n_vars - mask.count_ones()) << self.n_vars) | mask;
                self.touched.push(order_key);
            } else {
                self.state[index] ^= 1;
            }
        }
        let state = &self.state;
        let mask_bits = self.mask_bits;
        self.touched
            .retain(|&key| state[(key & mask_bits) as usize] & 1 != 0);
        self.touched.sort_unstable();
        F2BoolPoly {
            terms: self
                .touched
                .iter()
                .map(|&key| F2BoolMono::from_mask((key & mask_bits) as u64))
                .collect(),
            n_vars: p.n_vars,
        }
    }
}

fn dense_mul_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("PQ_F4_DISABLE_DENSE_MUL").as_deref() != Ok("1"))
}

fn multiply_for_f4(
    p: &F2BoolPoly,
    multiplier: u64,
    scratch: &mut Option<DenseMulScratch>,
    st: &mut F4Stats,
) -> F2BoolPoly {
    let Some(scratch) = scratch else {
        return p.mul_mono(F2BoolMono::from_mask(multiplier));
    };
    st.dense_mul_calls += 1;
    st.dense_mul_input_terms += p.terms.len() as u64;
    let product = scratch.multiply(p, multiplier);
    st.dense_mul_output_terms += product.terms.len() as u64;
    st.dense_mul_cancelled_terms += (p.terms.len().saturating_sub(product.terms.len())) as u64;
    st.dense_mul_scratch_bytes_max = st.dense_mul_scratch_bytes_max.max(scratch.bytes());
    product
}

struct EchelonOutput {
    pivot_leads: Vec<usize>,
    pivot_rows: EchelonRows,
    m4ri: bool,
    m4ri_block_width: u64,
    table_word_xors: u64,
    scratch_bytes: u64,
    m4ri_blocks: u64,
    consecutive_blocks: u64,
    trimmed_word_xors_avoided: u64,
}

impl EchelonOutput {
    fn len(&self) -> usize {
        self.pivot_leads.len()
    }

    fn lead(&self, index: usize) -> usize {
        self.pivot_leads[index]
    }

    fn flat_order_bytes(&self) -> Option<u64> {
        match &self.pivot_rows {
            EchelonRows::Nested(_) => None,
            EchelonRows::Flat(rows) => Some(rows.order_bytes()),
        }
    }

    fn unpack(&self, index: usize, columns: &Columns, n_vars: usize) -> F2BoolPoly {
        match &self.pivot_rows {
            EchelonRows::Nested(rows) => columns.unpack(&rows[index], n_vars),
            EchelonRows::Flat(rows) => {
                let bits = rows.row(index);
                let start = self.pivot_leads[index] / 64;
                let mut end = bits.len();
                while end > start && bits[end - 1] == 0 {
                    end -= 1;
                }
                columns.unpack_bits(bits, start, end, n_vars)
            }
        }
    }

    #[cfg(test)]
    fn bit_rows(&self) -> Vec<Vec<u64>> {
        match &self.pivot_rows {
            EchelonRows::Nested(rows) => rows.iter().map(|row| row.bits.clone()).collect(),
            EchelonRows::Flat(rows) => (0..self.len())
                .map(|index| rows.row(index).to_vec())
                .collect(),
        }
    }
}

/// Rows times words of a matrix's non-leading rows below which reductions by
/// the leading block stay on one thread.
const PAR_WORDS: usize = 1 << 16;

/// Streaming forward elimination. The initial block of rows with distinct
/// leads becomes pivots unchanged; reductions of the remaining rows by only
/// that fixed block are independent and may run in parallel. The serial pass
/// then resumes at exactly the point reached by the original row-at-a-time
/// loop, preserving pivots and the charged XOR count.
fn echelon_streaming(
    mut rows: Vec<Row>,
    n_cols: usize,
    word_xors: &mut u64,
    deadline: Option<Instant>,
) -> Option<EchelonOutput> {
    let mut pivot_of = vec![NONE; n_cols];
    let mut pivots: Vec<(usize, Row)> = Vec::with_capacity(rows.len());
    let mut fixed = 0usize;
    while fixed < rows.len() {
        match rows[fixed].lead() {
            Some(lead) if pivot_of[lead] != NONE => break,
            Some(lead) => pivot_of[lead] = 0,
            None => {}
        }
        fixed += 1;
    }
    let rest = rows.split_off(fixed);
    for mut row in rows {
        if let Some(lead) = row.lead() {
            pivot_of[lead] = pivots.len() as u32;
            pivots.push((lead, row));
        }
    }
    let mut rest = rest;
    {
        let (pivot_of, fixed_pivots) = (&pivot_of, &pivots);
        let expired = std::sync::atomic::AtomicBool::new(false);
        let by_fixed = |row: &mut Row| -> u64 {
            if deadline.is_some_and(|d| Instant::now() >= d) {
                expired.store(true, std::sync::atomic::Ordering::Relaxed);
                return 0;
            }
            let mut xors = 0u64;
            while let Some(lead) = row.lead() {
                let p = pivot_of[lead];
                if p == NONE {
                    break;
                }
                let pivot = &fixed_pivots[p as usize].1;
                let (from, to) = (lead / 64, pivot.end);
                for (a, b) in row.bits[from..to].iter_mut().zip(&pivot.bits[from..to]) {
                    *a ^= *b;
                }
                xors += (to - from) as u64;
                row.end = row.end.max(to);
            }
            xors
        };
        let words = n_cols.div_ceil(64).max(1);
        *word_xors += if rest.len() > 1 && rest.len() * words > PAR_WORDS {
            rest.par_iter_mut().map(by_fixed).sum::<u64>()
        } else {
            rest.iter_mut().map(by_fixed).sum::<u64>()
        };
        if expired.load(std::sync::atomic::Ordering::Relaxed) {
            return None;
        }
    }
    for (k, mut row) in rest.into_iter().enumerate() {
        if k % 128 == 0 && deadline.is_some_and(|d| Instant::now() >= d) {
            return None;
        }
        while let Some(lead) = row.lead() {
            match pivot_of[lead] {
                NONE => {
                    pivot_of[lead] = pivots.len() as u32;
                    pivots.push((lead, row));
                    break;
                }
                p => {
                    let pivot = &pivots[p as usize].1;
                    let (from, to) = (lead / 64, pivot.end);
                    for (a, b) in row.bits[from..to].iter_mut().zip(&pivot.bits[from..to]) {
                        *a ^= *b;
                    }
                    *word_xors += (to - from) as u64;
                    row.end = row.end.max(to);
                }
            }
        }
    }
    let (pivot_leads, pivot_rows) = pivots.into_iter().unzip();
    Some(EchelonOutput {
        pivot_leads,
        pivot_rows: EchelonRows::Nested(pivot_rows),
        m4ri: false,
        m4ri_block_width: 0,
        table_word_xors: 0,
        scratch_bytes: 0,
        m4ri_blocks: 0,
        consecutive_blocks: 0,
        trimmed_word_xors_avoided: 0,
    })
}

fn pq_m4ri_block_width() -> usize {
    static BLOCK_WIDTH: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *BLOCK_WIDTH.get_or_init(|| {
        std::env::var("PQ_F4_M4RI_BLOCK")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(8)
            .clamp(2, 10)
    })
}

#[derive(Default)]
struct PqM4riScratch {
    pivot_columns: Vec<usize>,
    block_pivots: Vec<Vec<u64>>,
    table: Vec<u64>,
}

thread_local! {
    static PQ_M4RI_SCRATCH: std::cell::RefCell<PqM4riScratch> =
        std::cell::RefCell::new(PqM4riScratch::default());
}

/// Method of Four Russians row echelon form. Pivot combinations are
/// materialized once per block, then each remaining row clears all block
/// columns with one suffix XOR. The returned rows span the same space and use
/// the same pivot columns as ordinary left-to-right Gaussian elimination. The
/// measured default block width is eight; `PQ_F4_M4RI_BLOCK` supports controlled
/// width ablations from two through ten.
fn echelon_m4ri(
    rows: Vec<Row>,
    n_cols: usize,
    word_xors: &mut u64,
    deadline: Option<Instant>,
) -> Option<EchelonOutput> {
    let mut matrix: Vec<Vec<u64>> = rows.into_iter().map(|row| row.bits).collect();
    let n_rows = matrix.len();
    let words = n_cols.div_ceil(64).max(1);
    let block_width = pq_m4ri_block_width();
    PQ_M4RI_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        scratch.pivot_columns.resize(block_width, 0);
        scratch.block_pivots.resize_with(block_width, Vec::new);
        for pivot in &mut scratch.block_pivots {
            pivot.resize(words, 0);
        }
        scratch.table.resize((1usize << block_width) * words, 0);
        let scratch_bytes = (scratch.pivot_columns.capacity() * std::mem::size_of::<usize>()
            + scratch
                .block_pivots
                .iter()
                .map(|row| row.capacity() * std::mem::size_of::<u64>())
                .sum::<usize>()
            + scratch.table.capacity() * std::mem::size_of::<u64>())
            as u64;
        let PqM4riScratch {
            pivot_columns,
            block_pivots,
            table,
        } = &mut *scratch;
        let mut all_pivots = Vec::with_capacity(n_rows.min(n_cols));
        let mut pivot_row = 0usize;
        let mut column = 0usize;
        let mut table_word_xors = 0u64;
        let mut m4ri_blocks = 0u64;
        let mut consecutive_blocks = 0u64;
        let mut trimmed_word_xors_avoided = 0u64;
        while pivot_row < n_rows && column < n_cols {
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            let block_start = pivot_row;
            let mut block_rows = 0usize;
            while block_rows < block_width && pivot_row < n_rows && column < n_cols {
                let next_pivot = block_start + block_rows;
                let (word, bit) = (column / 64, 1u64 << (column % 64));
                let mut found = None;
                for row in next_pivot..n_rows {
                    if row % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                        return None;
                    }
                    for (index, &pivot_column) in pivot_columns[..block_rows].iter().enumerate() {
                        let (pivot_word, pivot_bit) =
                            (pivot_column / 64, 1u64 << (pivot_column % 64));
                        if matrix[row][pivot_word] & pivot_bit != 0 {
                            for (target, &source) in matrix[row][pivot_word..words]
                                .iter_mut()
                                .zip(&block_pivots[index][pivot_word..words])
                            {
                                *target ^= source;
                            }
                            *word_xors += (words - pivot_word) as u64;
                        }
                    }
                    if matrix[row][word] & bit != 0 {
                        found = Some(row);
                        break;
                    }
                }
                if let Some(found) = found {
                    matrix.swap(next_pivot, found);
                    block_pivots[block_rows].copy_from_slice(&matrix[next_pivot]);
                    let (previous_pivots, current_pivots) = block_pivots.split_at_mut(block_rows);
                    let pivot = &current_pivots[0];
                    for previous in block_start..next_pivot {
                        if matrix[previous][word] & bit != 0 {
                            let block_index = previous - block_start;
                            for (target, &source) in matrix[previous][word..words]
                                .iter_mut()
                                .zip(&pivot[word..words])
                            {
                                *target ^= source;
                            }
                            for (target, &source) in previous_pivots[block_index][word..words]
                                .iter_mut()
                                .zip(&pivot[word..words])
                            {
                                *target ^= source;
                            }
                            *word_xors += 2 * (words - word) as u64;
                        }
                    }
                    pivot_columns[block_rows] = column;
                    block_rows += 1;
                }
                column += 1;
            }
            if block_rows == 0 {
                break;
            }
            m4ri_blocks += 1;

            let first_word = pivot_columns[0] / 64;
            let mut block_end = first_word + 1;
            for pivot in &block_pivots[..block_rows] {
                let mut end = words;
                while end > block_end && pivot[end - 1] == 0 {
                    end -= 1;
                }
                block_end = block_end.max(end);
            }
            let suffix_words = block_end - first_word;
            let combinations = 1usize << block_rows;
            trimmed_word_xors_avoided += ((combinations - 1) * (words - block_end)) as u64;
            table[..suffix_words].fill(0);
            let mut filled = 1usize;
            for pivot in block_pivots[..block_rows]
                .iter()
                .map(|row| &row[first_word..block_end])
            {
                let split = filled * suffix_words;
                let (source_tables, target_tables) =
                    table[..combinations * suffix_words].split_at_mut(split);
                for mask in 0..filled {
                    let offset = mask * suffix_words;
                    for index in 0..suffix_words {
                        target_tables[offset + index] =
                            source_tables[offset + index] ^ pivot[index];
                    }
                }
                let added_word_xors = (filled * suffix_words) as u64;
                *word_xors += added_word_xors;
                table_word_xors += added_word_xors;
                filled *= 2;
            }
            let consecutive = pivot_columns[..block_rows]
                .windows(2)
                .all(|pair| pair[1] == pair[0] + 1);
            consecutive_blocks += u64::from(consecutive);
            for row in block_start + block_rows..n_rows {
                if row % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                    return None;
                }
                let pattern = if consecutive {
                    let first_column = pivot_columns[0];
                    let (word, offset) = (first_column / 64, first_column % 64);
                    let mut packed = matrix[row][word] >> offset;
                    if offset + block_rows > 64 {
                        packed |= matrix[row][word + 1] << (64 - offset);
                    }
                    packed as usize & ((1usize << block_rows) - 1)
                } else {
                    let mut pattern = 0usize;
                    for (index, &pivot_column) in pivot_columns[..block_rows].iter().enumerate() {
                        if matrix[row][pivot_column / 64] & (1u64 << (pivot_column % 64)) != 0 {
                            pattern |= 1usize << index;
                        }
                    }
                    pattern
                };
                if pattern != 0 {
                    let offset = pattern * suffix_words;
                    for (target, &source) in matrix[row][first_word..block_end]
                        .iter_mut()
                        .zip(&table[offset..offset + suffix_words])
                    {
                        *target ^= source;
                    }
                    *word_xors += suffix_words as u64;
                    trimmed_word_xors_avoided += (words - block_end) as u64;
                }
            }
            all_pivots.extend_from_slice(&pivot_columns[..block_rows]);
            pivot_row += block_rows;
        }

        let mut pivot_rows = Vec::with_capacity(all_pivots.len());
        for (bits, &lead) in matrix.iter_mut().take(all_pivots.len()).zip(&all_pivots) {
            let start = lead / 64;
            let mut end = bits.len();
            while end > start && bits[end - 1] == 0 {
                end -= 1;
            }
            pivot_rows.push(Row {
                bits: std::mem::take(bits),
                start,
                end,
            });
        }
        Some(EchelonOutput {
            pivot_leads: all_pivots,
            pivot_rows: EchelonRows::Nested(pivot_rows),
            m4ri: true,
            m4ri_block_width: block_width as u64,
            table_word_xors,
            scratch_bytes,
            m4ri_blocks,
            consecutive_blocks,
            trimmed_word_xors_avoided,
        })
    })
}

fn echelon_m4ri_flat(
    mut matrix: FlatRows,
    n_cols: usize,
    word_xors: &mut u64,
    deadline: Option<Instant>,
) -> Option<EchelonOutput> {
    let n_rows = matrix.rows;
    let words = matrix.words;
    debug_assert_eq!(words, n_cols.div_ceil(64).max(1));
    let block_width = pq_m4ri_block_width();
    PQ_M4RI_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        scratch.pivot_columns.resize(block_width, 0);
        scratch.block_pivots.resize_with(block_width, Vec::new);
        for pivot in &mut scratch.block_pivots {
            pivot.resize(words, 0);
        }
        scratch.table.resize((1usize << block_width) * words, 0);
        let scratch_bytes = (scratch.pivot_columns.capacity() * std::mem::size_of::<usize>()
            + scratch
                .block_pivots
                .iter()
                .map(|row| row.capacity() * std::mem::size_of::<u64>())
                .sum::<usize>()
            + scratch.table.capacity() * std::mem::size_of::<u64>())
            as u64;
        let PqM4riScratch {
            pivot_columns,
            block_pivots,
            table,
        } = &mut *scratch;
        let mut all_pivots = Vec::with_capacity(n_rows.min(n_cols));
        let mut pivot_row = 0usize;
        let mut column = 0usize;
        let mut table_word_xors = 0u64;
        let mut m4ri_blocks = 0u64;
        let mut consecutive_blocks = 0u64;
        let mut trimmed_word_xors_avoided = 0u64;
        while pivot_row < n_rows && column < n_cols {
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            let block_start = pivot_row;
            let mut block_rows = 0usize;
            while block_rows < block_width && pivot_row < n_rows && column < n_cols {
                let next_pivot = block_start + block_rows;
                let (word, bit) = (column / 64, 1u64 << (column % 64));
                let mut found = None;
                for row in next_pivot..n_rows {
                    if row % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                        return None;
                    }
                    for (index, &pivot_column) in pivot_columns[..block_rows].iter().enumerate() {
                        let (pivot_word, pivot_bit) =
                            (pivot_column / 64, 1u64 << (pivot_column % 64));
                        let target = matrix.row_mut(row);
                        if target[pivot_word] & pivot_bit != 0 {
                            for (target, &source) in target[pivot_word..words]
                                .iter_mut()
                                .zip(&block_pivots[index][pivot_word..words])
                            {
                                *target ^= source;
                            }
                            *word_xors += (words - pivot_word) as u64;
                        }
                    }
                    if matrix.row(row)[word] & bit != 0 {
                        found = Some(row);
                        break;
                    }
                }
                if let Some(found) = found {
                    matrix.swap(next_pivot, found);
                    block_pivots[block_rows].copy_from_slice(matrix.row(next_pivot));
                    let (previous_pivots, current_pivots) = block_pivots.split_at_mut(block_rows);
                    let pivot = &current_pivots[0];
                    for previous in block_start..next_pivot {
                        let target = matrix.row_mut(previous);
                        if target[word] & bit != 0 {
                            let block_index = previous - block_start;
                            for (target, &source) in
                                target[word..words].iter_mut().zip(&pivot[word..words])
                            {
                                *target ^= source;
                            }
                            for (target, &source) in previous_pivots[block_index][word..words]
                                .iter_mut()
                                .zip(&pivot[word..words])
                            {
                                *target ^= source;
                            }
                            *word_xors += 2 * (words - word) as u64;
                        }
                    }
                    pivot_columns[block_rows] = column;
                    block_rows += 1;
                }
                column += 1;
            }
            if block_rows == 0 {
                break;
            }
            m4ri_blocks += 1;

            let first_word = pivot_columns[0] / 64;
            let mut block_end = first_word + 1;
            for pivot in &block_pivots[..block_rows] {
                let mut end = words;
                while end > block_end && pivot[end - 1] == 0 {
                    end -= 1;
                }
                block_end = block_end.max(end);
            }
            let suffix_words = block_end - first_word;
            let combinations = 1usize << block_rows;
            trimmed_word_xors_avoided += ((combinations - 1) * (words - block_end)) as u64;
            table[..suffix_words].fill(0);
            let mut filled = 1usize;
            for pivot in block_pivots[..block_rows]
                .iter()
                .map(|row| &row[first_word..block_end])
            {
                let split = filled * suffix_words;
                let (source_tables, target_tables) =
                    table[..combinations * suffix_words].split_at_mut(split);
                for mask in 0..filled {
                    let offset = mask * suffix_words;
                    for index in 0..suffix_words {
                        target_tables[offset + index] =
                            source_tables[offset + index] ^ pivot[index];
                    }
                }
                let added_word_xors = (filled * suffix_words) as u64;
                *word_xors += added_word_xors;
                table_word_xors += added_word_xors;
                filled *= 2;
            }
            let consecutive = pivot_columns[..block_rows]
                .windows(2)
                .all(|pair| pair[1] == pair[0] + 1);
            consecutive_blocks += u64::from(consecutive);
            for row in block_start + block_rows..n_rows {
                if row % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                    return None;
                }
                let pattern = if consecutive {
                    let first_column = pivot_columns[0];
                    let (word, offset) = (first_column / 64, first_column % 64);
                    let target = matrix.row(row);
                    let mut packed = target[word] >> offset;
                    if offset + block_rows > 64 {
                        packed |= target[word + 1] << (64 - offset);
                    }
                    packed as usize & ((1usize << block_rows) - 1)
                } else {
                    let target = matrix.row(row);
                    let mut pattern = 0usize;
                    for (index, &pivot_column) in pivot_columns[..block_rows].iter().enumerate() {
                        if target[pivot_column / 64] & (1u64 << (pivot_column % 64)) != 0 {
                            pattern |= 1usize << index;
                        }
                    }
                    pattern
                };
                if pattern != 0 {
                    let offset = pattern * suffix_words;
                    for (target, &source) in matrix.row_mut(row)[first_word..block_end]
                        .iter_mut()
                        .zip(&table[offset..offset + suffix_words])
                    {
                        *target ^= source;
                    }
                    *word_xors += suffix_words as u64;
                    trimmed_word_xors_avoided += (words - block_end) as u64;
                }
            }
            all_pivots.extend_from_slice(&pivot_columns[..block_rows]);
            pivot_row += block_rows;
        }

        matrix.order.truncate(all_pivots.len());
        matrix.rows = all_pivots.len();
        Some(EchelonOutput {
            pivot_leads: all_pivots,
            pivot_rows: EchelonRows::Flat(matrix),
            m4ri: true,
            m4ri_block_width: block_width as u64,
            table_word_xors,
            scratch_bytes,
            m4ri_blocks,
            consecutive_blocks,
            trimmed_word_xors_avoided,
        })
    })
}

fn m4ri_shape(n_rows: usize, n_cols: usize) -> bool {
    n_rows >= 128
        && n_cols >= 256
        && n_cols <= 4 * n_rows
        && std::env::var("PQ_F4_DISABLE_M4RI").as_deref() != Ok("1")
}

fn flat_m4ri_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("PQ_F4_FLAT_M4RI").as_deref() == Ok("1"))
}

fn flat_m4ri_min_words() -> usize {
    static MIN_WORDS: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *MIN_WORDS.get_or_init(|| {
        std::env::var("PQ_F4_FLAT_M4RI_MIN_WORDS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(0)
    })
}

fn pack_rows_for_echelon<'a>(
    columns: &Columns,
    n_rows: usize,
    rows: impl Iterator<Item = &'a F2BoolPoly>,
    deadline: Option<Instant>,
) -> Option<PackedRows> {
    if flat_m4ri_enabled()
        && m4ri_shape(n_rows, columns.monos.len())
        && n_rows.saturating_mul(columns.words()) >= flat_m4ri_min_words()
    {
        let mut packed = FlatRows::new(n_rows, columns.words());
        let mut packed_rows = 0usize;
        for (row_index, polynomial) in rows.enumerate() {
            if row_index % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            columns.pack_into(polynomial, packed.row_mut(row_index));
            packed_rows += 1;
        }
        debug_assert_eq!(packed_rows, n_rows);
        Some(PackedRows::Flat(packed))
    } else {
        let mut packed = Vec::with_capacity(n_rows);
        for (row_index, polynomial) in rows.enumerate() {
            if row_index % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            packed.push(columns.pack(polynomial));
        }
        Some(PackedRows::Nested(packed))
    }
}

fn echelon(
    rows: PackedRows,
    n_cols: usize,
    word_xors: &mut u64,
    deadline: Option<Instant>,
) -> Option<EchelonOutput> {
    match rows {
        PackedRows::Flat(rows) => {
            debug_assert!(m4ri_shape(rows.rows, n_cols));
            echelon_m4ri_flat(rows, n_cols, word_xors, deadline)
        }
        PackedRows::Nested(rows) if m4ri_shape(rows.len(), n_cols) => {
            echelon_m4ri(rows, n_cols, word_xors, deadline)
        }
        PackedRows::Nested(rows) => echelon_streaming(rows, n_cols, word_xors, deadline),
    }
}

fn is_one(p: &F2BoolPoly) -> bool {
    p.terms.len() == 1 && p.terms[0].mask == 0
}

fn degree(p: &F2BoolPoly) -> u32 {
    p.terms.iter().map(|t| t.degree()).max().unwrap_or(0)
}

struct State {
    n_vars: usize,
    polys: Vec<F2BoolPoly>,
    lm: Vec<u64>,
    active: Vec<bool>,
    active_indices: Vec<usize>,
    pairs: Vec<Pair>,
    pair_select_scratch: PairSelectScratch,
}

#[derive(Default)]
struct PairSelectScratch {
    grouped: Vec<(usize, u64)>,
    lcms: FastU64Set,
    survivor: FastU64Map<usize>,
    dense: Option<DensePairSelectScratch>,
}

struct DensePairSelectScratch {
    epoch: u32,
    stamp: Vec<u32>,
    first_noncoprime: Vec<u32>,
    noncoprime_count: Vec<u32>,
    has_coprime: Vec<u8>,
    survivor: Vec<u32>,
    touched: Vec<u32>,
}

impl DensePairSelectScratch {
    fn new(size: usize) -> Self {
        Self {
            epoch: 0,
            stamp: vec![0; size],
            first_noncoprime: vec![NONE; size],
            noncoprime_count: vec![0; size],
            has_coprime: vec![0; size],
            survivor: vec![NONE; size],
            touched: Vec::new(),
        }
    }

    fn bytes(&self) -> u64 {
        (self.stamp.capacity() * std::mem::size_of::<u32>()
            + self.first_noncoprime.capacity() * std::mem::size_of::<u32>()
            + self.noncoprime_count.capacity() * std::mem::size_of::<u32>()
            + self.has_coprime.capacity() * std::mem::size_of::<u8>()
            + self.survivor.capacity() * std::mem::size_of::<u32>()
            + self.touched.capacity() * std::mem::size_of::<u32>()) as u64
    }

    fn begin(&mut self) {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.stamp.fill(0);
            self.epoch = 1;
        }
        self.touched.clear();
    }
}

fn active_leading_monomial_index(s: &State, active: &[usize]) -> FastU64Map<usize> {
    let mut active_by_lm = FastU64Map::default();
    for &g in active {
        active_by_lm
            .entry(s.lm[g])
            .and_modify(|best: &mut usize| {
                if (s.polys[g].terms.len(), g) < (s.polys[*best].terms.len(), *best) {
                    *best = g;
                }
            })
            .or_insert(g);
    }
    active_by_lm
}

/// Becker--Weispfenning `UPDATE` selection for the pairs made by one new
/// leading monomial.  The direct formulation checks every candidate against
/// every other candidate.  Here equal LCMs are grouped and proper divisor
/// LCMs are found by exact submask lookup.  Iterating the original candidates
/// in reverse at the end preserves the direct algorithm's pair order and its
/// lowest-index representative for duplicate non-coprime LCMs.
fn select_new_pairs(
    n_vars: usize,
    lh: u64,
    lm: &[u64],
    active: &[usize],
    scratch: &mut PairSelectScratch,
    st: &mut F4Stats,
) -> Vec<(usize, u64)> {
    st.active_candidate_visits += active.len() as u64;
    if active.is_empty() {
        return Vec::new();
    }
    if n_vars <= 20 {
        let size = 1usize << n_vars;
        if scratch.dense.as_ref().map(|dense| dense.stamp.len()) != Some(size) {
            scratch.dense = Some(DensePairSelectScratch::new(size));
        }
        let dense = scratch.dense.as_mut().unwrap();
        st.pair_dense_select_calls += 1;
        dense.begin();
        dense.touched.reserve(active.len());
        for &g in active {
            debug_assert!(g < u32::MAX as usize);
            let lcm = (lh | lm[g]) as usize;
            debug_assert!(lcm < size);
            if dense.stamp[lcm] != dense.epoch {
                dense.stamp[lcm] = dense.epoch;
                dense.first_noncoprime[lcm] = NONE;
                dense.noncoprime_count[lcm] = 0;
                dense.has_coprime[lcm] = 0;
                dense.survivor[lcm] = NONE;
                dense.touched.push(lcm as u32);
            }
            if lh & lm[g] == 0 {
                dense.has_coprime[lcm] = 1;
            } else {
                dense.noncoprime_count[lcm] += 1;
                let first = dense.first_noncoprime[lcm];
                if first == NONE || g < first as usize {
                    dense.first_noncoprime[lcm] = g as u32;
                }
            }
        }
        st.pair_lcm_groups += dense.touched.len() as u64;
        st.pair_dense_scratch_bytes_max = st.pair_dense_scratch_bytes_max.max(dense.bytes());
        for &lcm_u32 in &dense.touched {
            let lcm = lcm_u32 as usize;
            let mut proper_cover = false;
            let remainder = lcm & !(lh as usize);
            if remainder != 0 {
                let mut submask = (remainder - 1) & remainder;
                loop {
                    st.pair_cover_lookups += 1;
                    if dense.stamp[(lh as usize) | submask] == dense.epoch {
                        proper_cover = true;
                        break;
                    }
                    if submask == 0 {
                        break;
                    }
                    submask = (submask - 1) & remainder;
                }
            }
            let count = dense.noncoprime_count[lcm] as u64;
            if proper_cover || dense.has_coprime[lcm] != 0 {
                st.pairs_chain_skipped += count;
            } else {
                let first = dense.first_noncoprime[lcm];
                if first != NONE {
                    dense.survivor[lcm] = first;
                    st.pairs_chain_skipped += count.saturating_sub(1);
                }
            }
        }
        let mut selected = Vec::new();
        for &g in active.iter().rev() {
            let lcm = (lh | lm[g]) as usize;
            if lh & lm[g] == 0 {
                st.pairs_product_skipped += 1;
            } else if dense.survivor[lcm] == g as u32 {
                selected.push((g, lcm as u64));
            }
        }
        return selected;
    }

    st.pair_sorted_select_calls += 1;
    let PairSelectScratch {
        grouped,
        lcms,
        survivor,
        ..
    } = scratch;
    grouped.clear();
    grouped.extend(active.iter().map(|&g| (g, lh | lm[g])));
    grouped.sort_unstable_by_key(|&(_, lcm)| lcm);
    lcms.clear();
    lcms.reserve(grouped.len());
    let mut previous_lcm = None;
    for &(_, lcm) in grouped.iter() {
        if previous_lcm != Some(lcm) {
            lcms.insert(lcm);
            previous_lcm = Some(lcm);
        }
    }
    st.pair_lcm_groups += lcms.len() as u64;
    survivor.clear();
    survivor.reserve(lcms.len());
    let mut start = 0usize;
    while start < grouped.len() {
        let lcm = grouped[start].1;
        let mut end = start + 1;
        while end < grouped.len() && grouped[end].1 == lcm {
            end += 1;
        }
        let group = &grouped[start..end];
        let mut has_coprime = false;
        let mut first = None;
        let mut count = 0usize;
        for &(g, _) in group {
            if lh & lm[g] == 0 {
                has_coprime = true;
            } else {
                count += 1;
                first = Some(first.map_or(g, |current: usize| current.min(g)));
            }
        }
        let mut proper_cover = false;
        let remainder = lcm & !lh;
        if remainder != 0 {
            let mut submask = (remainder - 1) & remainder;
            loop {
                st.pair_cover_lookups += 1;
                if lcms.contains(&(lh | submask)) {
                    proper_cover = true;
                    break;
                }
                if submask == 0 {
                    break;
                }
                submask = (submask - 1) & remainder;
            }
        }
        if proper_cover || has_coprime {
            st.pairs_chain_skipped += count as u64;
        } else if let Some(g) = first {
            survivor.insert(lcm, g);
            st.pairs_chain_skipped += count.saturating_sub(1) as u64;
        }
        start = end;
    }

    let mut selected = Vec::with_capacity(survivor.len());
    for &g in active.iter().rev() {
        let lcm = lh | lm[g];
        if lh & lm[g] == 0 {
            st.pairs_product_skipped += 1;
        } else if survivor.get(&lcm) == Some(&g) {
            selected.push((g, lcm));
        }
    }
    selected
}

fn critical_pair_survives(p: &Pair, lh: u64, lm: &[u64]) -> bool {
    match p.kind {
        PairKind::Field(..) => true,
        PairKind::Critical(i, j) => {
            lh & !p.lcm != 0 || (lm[i] | lh) == p.lcm || (lm[j] | lh) == p.lcm
        }
    }
}

fn batch_basis_insert_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("PQ_F4_DISABLE_BATCH_INSERT").as_deref() != Ok("1"))
}

fn critical_pair_survives_batch(
    p: &Pair,
    start: usize,
    batch_lms: &[u64],
    latest_batch_index: &FastU64Map<usize>,
    lm: &[u64],
    st: &mut F4Stats,
) -> bool {
    if matches!(p.kind, PairKind::Field(..)) || start >= batch_lms.len() {
        return true;
    }
    let degree = p.lcm.count_ones();
    let submask_count = if degree < usize::BITS {
        (1usize << degree) - 1
    } else {
        usize::MAX
    };
    if submask_count <= batch_lms.len() - start {
        let mut lh = p.lcm;
        while lh != 0 {
            if !critical_pair_survives(p, lh, lm) {
                st.pair_prune_tests += 1;
                st.pair_prune_submask_lookups += 1;
                if latest_batch_index
                    .get(&lh)
                    .is_some_and(|&latest| latest >= start)
                {
                    return false;
                }
            }
            lh = (lh - 1) & p.lcm;
        }
    } else {
        for &lh in &batch_lms[start..] {
            st.pair_prune_tests += 1;
            st.pair_prune_linear_tests += 1;
            if !critical_pair_survives(p, lh, lm) {
                return false;
            }
        }
    }
    true
}

impl State {
    /// Becker–Weispfenning `UPDATE` for the critical pairs, plus the
    /// field pairs of the new element.
    fn insert(&mut self, h_poly: F2BoolPoly, st: &mut F4Stats) {
        let h = self.polys.len();
        let lh = h_poly.lt().expect("a new element is non-zero").mask;
        self.polys.push(h_poly);
        self.lm.push(lh);
        self.active.push(false);

        let dh = lh.count_ones();
        let mut bits = lh;
        while bits != 0 {
            let v = bits.trailing_zeros();
            bits &= bits - 1;
            self.pairs.push(Pair {
                kind: PairKind::Field(h, v),
                lcm: lh,
                deg: dh + 1,
            });
        }

        // New pairs (h, g): keep one per minimal LCM (chain criterion).
        let selected = select_new_pairs(
            self.n_vars,
            lh,
            &self.lm,
            &self.active_indices,
            &mut self.pair_select_scratch,
            st,
        );
        // Old pairs whose lcm `h` divides strictly on both sides.
        let lm = &self.lm;
        let mut dropped = 0u64;
        let mut tests = 0u64;
        self.pairs.retain(|p| {
            if matches!(p.kind, PairKind::Critical(..)) {
                tests += 1;
            }
            let keep = critical_pair_survives(p, lh, lm);
            dropped += u64::from(!keep);
            keep
        });
        st.pairs_chain_skipped += dropped;
        st.pair_prune_tests += tests;
        st.pair_prune_linear_tests += tests;
        st.pair_prune_passes += 1;
        for (g, l) in selected {
            self.pairs.push(Pair {
                kind: PairKind::Critical(g, h),
                lcm: l,
                deg: l.count_ones(),
            });
        }
        st.active_deactivation_tests += self.active_indices.len() as u64;
        let lm = &self.lm;
        let active = &mut self.active;
        self.active_indices.retain(|&g| {
            let keep = lh & !lm[g] != 0;
            if !keep {
                active[g] = false;
            }
            keep
        });
        self.active[h] = true;
        self.active_indices.push(h);
    }

    /// Install one matrix's new basis elements in the same order as repeated
    /// [`State::insert`], but defer pending-pair pruning until the whole batch
    /// is known. New-pair selection and active-basis updates remain sequential;
    /// an old pair is tested against every batch leader until the first one
    /// that would have removed it, and a pair born in the batch starts with the
    /// following leader. This preserves the final pair order and all UPDATE
    /// decisions while avoiding one full queue scan per new element.
    fn insert_batch(&mut self, elements: Vec<F2BoolPoly>, st: &mut F4Stats) {
        if elements.len() <= 1 {
            for element in elements {
                self.insert(element, st);
            }
            return;
        }
        st.batch_insert_groups += 1;
        st.batch_insert_elements += elements.len() as u64;
        let existing = std::mem::take(&mut self.pairs);
        let mut added: Vec<(Pair, Option<usize>)> = Vec::new();
        let mut batch_lms = Vec::with_capacity(elements.len());

        for (batch_index, h_poly) in elements.into_iter().enumerate() {
            let h = self.polys.len();
            let lh = h_poly.lt().expect("a new element is non-zero").mask;
            self.polys.push(h_poly);
            self.lm.push(lh);
            self.active.push(false);
            batch_lms.push(lh);

            let dh = lh.count_ones();
            let mut bits = lh;
            while bits != 0 {
                let v = bits.trailing_zeros();
                bits &= bits - 1;
                added.push((
                    Pair {
                        kind: PairKind::Field(h, v),
                        lcm: lh,
                        deg: dh + 1,
                    },
                    None,
                ));
            }

            let selected = select_new_pairs(
                self.n_vars,
                lh,
                &self.lm,
                &self.active_indices,
                &mut self.pair_select_scratch,
                st,
            );
            for (g, l) in selected {
                added.push((
                    Pair {
                        kind: PairKind::Critical(g, h),
                        lcm: l,
                        deg: l.count_ones(),
                    },
                    Some(batch_index),
                ));
            }
            st.active_deactivation_tests += self.active_indices.len() as u64;
            let lm = &self.lm;
            let active = &mut self.active;
            self.active_indices.retain(|&g| {
                let keep = lh & !lm[g] != 0;
                if !keep {
                    active[g] = false;
                }
                keep
            });
            self.active[h] = true;
            self.active_indices.push(h);
        }

        let mut pairs = Vec::with_capacity(existing.len() + added.len());
        let mut dropped = 0u64;
        let mut latest_batch_index = FastU64Map::default();
        for (index, &lh) in batch_lms.iter().enumerate() {
            latest_batch_index.insert(lh, index);
        }
        for p in existing {
            let keep =
                critical_pair_survives_batch(&p, 0, &batch_lms, &latest_batch_index, &self.lm, st);
            if keep {
                pairs.push(p);
            } else {
                dropped += 1;
            }
        }
        for (p, birth) in added {
            let keep = birth.is_none_or(|birth| {
                critical_pair_survives_batch(
                    &p,
                    birth + 1,
                    &batch_lms,
                    &latest_batch_index,
                    &self.lm,
                    st,
                )
            });
            if keep {
                pairs.push(p);
            } else {
                dropped += 1;
            }
        }
        self.pairs = pairs;
        st.pairs_chain_skipped += dropped;
        st.pair_prune_passes += 1;
    }

    /// Active basis elements whose leading monomial divides `m`, the
    /// shortest first.
    fn reducer_for(
        &self,
        m: u64,
        active: &[usize],
        active_by_lm: &FastU64Map<usize>,
        st: &mut F4Stats,
        deadline: Option<Instant>,
    ) -> Result<Option<usize>, ()> {
        let mut best: Option<usize> = None;
        let degree = m.count_ones();
        let submask_count = if degree < usize::BITS {
            (1usize << degree) - 1
        } else {
            usize::MAX
        };
        if submask_count <= active.len() {
            let mut divisor = m;
            let mut index = 0usize;
            while divisor != 0 {
                if index % 1024 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                    return Err(());
                }
                st.divisor_tests += 1;
                st.divisor_submask_lookups += 1;
                if let Some(&g) = active_by_lm.get(&divisor) {
                    let candidate = (self.polys[g].terms.len(), g);
                    if best.is_none_or(|b| candidate < (self.polys[b].terms.len(), b)) {
                        best = Some(g);
                    }
                }
                divisor = (divisor - 1) & m;
                index += 1;
            }
        } else {
            for (index, &g) in active.iter().enumerate() {
                if index % 1024 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                    return Err(());
                }
                st.divisor_tests += 1;
                st.divisor_linear_tests += 1;
                if self.lm[g] & !m == 0
                    && best.is_none_or(|b| self.polys[g].terms.len() < self.polys[b].terms.len())
                {
                    best = Some(g);
                }
            }
        }
        Ok(best)
    }
}

/// **Compute the reduced Gröbner basis** of the ideal `initial` generates
/// in `F_2[v]/(v² − v)` by F4, within an optional wall-clock budget.
///
/// On a budget or size overrun the returned set generates the same ideal
/// but is not a Gröbner basis, and the stats say so.
pub fn groebner_basis_f4(
    initial: Vec<F2BoolPoly>,
    n_vars: usize,
    budget: Option<Duration>,
) -> (Vec<F2BoolPoly>, F4Stats) {
    let started = Instant::now();
    let deadline = budget.map(|b| started + b);
    let mut st = F4Stats::default();
    let mut dense_mul_scratch = if n_vars <= 20 && dense_mul_enabled() {
        Some(DenseMulScratch::new(n_vars))
    } else {
        None
    };
    let one = || vec![F2BoolPoly::one(n_vars)];
    let finish = |basis: Vec<F2BoolPoly>, mut st: F4Stats| {
        st.basis_len = basis.len() as u64;
        st.wall_ns = started.elapsed().as_nanos() as u64;
        (basis, st)
    };

    let inputs: Vec<F2BoolPoly> = initial.into_iter().filter(|p| !p.is_zero()).collect();
    if inputs.iter().any(is_one) {
        return finish(one(), st);
    }
    for p in &inputs {
        st.max_poly_degree = st.max_poly_degree.max(degree(p));
    }

    // The initial echelon: distinct leading monomials to start from.
    let t = Instant::now();
    let cols = Columns::from_monomials(
        inputs.iter().flat_map(|p| p.terms.iter().map(|t| t.mask)),
        n_vars,
    );
    cols.record_dense_cost(&mut st);
    let Some(rows) = pack_rows_for_echelon(&cols, inputs.len(), inputs.iter(), deadline) else {
        st.timed_out = true;
        st.build_ns += t.elapsed().as_nanos() as u64;
        return finish(inputs, st);
    };
    st.build_ns += t.elapsed().as_nanos() as u64;
    let t = Instant::now();
    let Some(initial_echelon) = echelon(rows, cols.monos.len(), &mut st.word_xors, deadline) else {
        st.timed_out = true;
        st.eliminate_ns += t.elapsed().as_nanos() as u64;
        return finish(inputs, st);
    };
    st.m4ri_matrices += u64::from(initial_echelon.m4ri);
    st.m4ri_block_width_max = st
        .m4ri_block_width_max
        .max(initial_echelon.m4ri_block_width);
    st.m4ri_table_word_xors += initial_echelon.table_word_xors;
    st.m4ri_blocks += initial_echelon.m4ri_blocks;
    st.m4ri_consecutive_blocks += initial_echelon.consecutive_blocks;
    st.m4ri_trimmed_word_xors_avoided += initial_echelon.trimmed_word_xors_avoided;
    st.m4ri_scratch_bytes_max = st.m4ri_scratch_bytes_max.max(initial_echelon.scratch_bytes);
    if let Some(order_bytes) = initial_echelon.flat_order_bytes() {
        st.flat_m4ri_matrices += 1;
        st.flat_m4ri_order_bytes_max = st.flat_m4ri_order_bytes_max.max(order_bytes);
    }
    st.eliminate_ns += t.elapsed().as_nanos() as u64;
    let mut start: Vec<F2BoolPoly> = (0..initial_echelon.len())
        .map(|index| initial_echelon.unpack(index, &cols, n_vars))
        .collect();
    if start.iter().any(is_one) {
        return finish(one(), st);
    }
    start.sort_by(|a, b| cmp_mono(a.lt().unwrap(), b.lt().unwrap()));

    let mut s = State {
        n_vars,
        polys: Vec::new(),
        lm: Vec::new(),
        active: Vec::new(),
        active_indices: Vec::new(),
        pairs: Vec::new(),
        pair_select_scratch: PairSelectScratch::default(),
    };
    let mut lcm_columns = MonomialSet::new(n_vars);
    let mut examined = MonomialSet::new(n_vars);
    let mut no_divisor = MonomialSet::new(n_vars);
    let t = Instant::now();
    if batch_basis_insert_enabled() {
        s.insert_batch(start, &mut st);
    } else {
        for p in start {
            s.insert(p, &mut st);
        }
    }
    st.pair_update_ns += t.elapsed().as_nanos() as u64;

    while !s.pairs.is_empty() {
        if deadline.is_some_and(|d| Instant::now() >= d) {
            st.timed_out = true;
            break;
        }
        let d = s.pairs.iter().map(|p| p.deg).min().unwrap();
        let (selected, rest): (Vec<Pair>, Vec<Pair>) = s.pairs.drain(..).partition(|p| p.deg == d);
        s.pairs = rest;
        st.steps += 1;
        st.degree_reached = st.degree_reached.max(d);

        let t = Instant::now();
        // The rows the step reduces: both halves of every critical pair
        // (their common leading monomial is the lcm, which one half will
        // pivot and the other lose), and every field product.
        let mut seen_rows: HashSet<(u64, usize)> = HashSet::new();
        let mut half_rows: Vec<F2BoolPoly> = Vec::new();
        let mut field_rows: Vec<F2BoolPoly> = Vec::new();
        lcm_columns.clear();
        let mut row_terms = 0usize;
        'selected_pairs: for (selected_index, p) in selected.iter().enumerate() {
            if selected_index % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                st.timed_out = true;
                break;
            }
            match p.kind {
                PairKind::Critical(i, j) => {
                    st.pairs_reduced += 1;
                    lcm_columns.insert(p.lcm);
                    for g in [i, j] {
                        let mult = p.lcm & !s.lm[g];
                        if seen_rows.insert((mult, g)) {
                            let source_terms = s.polys[g].terms.len();
                            if symbolic_cap_exceeded(
                                &mut st,
                                lcm_columns.len(),
                                0,
                                row_terms.saturating_add(source_terms),
                                half_rows.len() + field_rows.len() + 1,
                                seen_rows.len(),
                            ) {
                                break 'selected_pairs;
                            }
                            let row =
                                multiply_for_f4(&s.polys[g], mult, &mut dense_mul_scratch, &mut st);
                            row_terms = row_terms.saturating_add(row.terms.len());
                            half_rows.push(row);
                        }
                    }
                }
                PairKind::Field(g, v) => {
                    st.field_pairs_reduced += 1;
                    // A field multiplier lies inside `LM(g)` and a half's is
                    // disjoint from it, so the two kinds of key never meet.
                    let mult = 1u64 << v;
                    if seen_rows.insert((mult, g)) {
                        let source_terms = s.polys[g].terms.len();
                        if symbolic_cap_exceeded(
                            &mut st,
                            lcm_columns.len(),
                            0,
                            row_terms.saturating_add(source_terms),
                            half_rows.len() + field_rows.len() + 1,
                            seen_rows.len(),
                        ) {
                            break 'selected_pairs;
                        }
                        let prod =
                            multiply_for_f4(&s.polys[g], mult, &mut dense_mul_scratch, &mut st);
                        if !prod.is_zero() {
                            row_terms = row_terms.saturating_add(prod.terms.len());
                            field_rows.push(prod);
                        }
                    }
                }
            }
        }
        if st.timed_out || st.oversize {
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }

        // Symbolic preprocessing.  Every monomial other than an S-pair's
        // lcm is examined once; a divisible one gets a reducer led by it.
        let active = s.active_indices.clone();
        let active_by_lm = active_leading_monomial_index(&s, &active);
        examined.copy_from(&lcm_columns);
        no_divisor.clear();
        let mut queue: Vec<u64> = Vec::new();
        'seed_queue: for p in half_rows.iter().chain(field_rows.iter()) {
            for t in &p.terms {
                if examined.insert(t.mask) {
                    queue.push(t.mask);
                    if queue.len() % 1024 == 0
                        && symbolic_cap_exceeded(
                            &mut st,
                            examined.len(),
                            queue.len(),
                            row_terms,
                            half_rows.len() + field_rows.len(),
                            seen_rows.len(),
                        )
                    {
                        break 'seed_queue;
                    }
                }
            }
        }
        if st.oversize {
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }
        let mut reducers: Vec<F2BoolPoly> = Vec::new();
        let mut symbolic_steps = 0usize;
        while let Some(m) = queue.pop() {
            if symbolic_steps % 128 == 0 && deadline.is_some_and(|limit| Instant::now() >= limit) {
                st.timed_out = true;
                break;
            }
            symbolic_steps += 1;
            match s.reducer_for(m, &active, &active_by_lm, &mut st, deadline) {
                Err(()) => {
                    st.timed_out = true;
                    break;
                }
                Ok(Some(g)) => {
                    let source_terms = s.polys[g].terms.len();
                    if symbolic_cap_exceeded(
                        &mut st,
                        examined.len(),
                        queue.len(),
                        row_terms.saturating_add(source_terms),
                        half_rows.len() + field_rows.len() + reducers.len() + 1,
                        seen_rows.len(),
                    ) {
                        break;
                    }
                    let r =
                        multiply_for_f4(&s.polys[g], m & !s.lm[g], &mut dense_mul_scratch, &mut st);
                    debug_assert_eq!(
                        r.lt().map(|t| t.mask),
                        Some(m),
                        "a reducer must lead with its monomial"
                    );
                    row_terms = row_terms.saturating_add(r.terms.len());
                    for (term_index, t) in r.terms.iter().enumerate() {
                        if examined.insert(t.mask) {
                            queue.push(t.mask);
                        }
                        if term_index % 1024 == 0
                            && symbolic_cap_exceeded(
                                &mut st,
                                examined.len(),
                                queue.len(),
                                row_terms,
                                half_rows.len() + field_rows.len() + reducers.len() + 1,
                                seen_rows.len(),
                            )
                        {
                            break;
                        }
                    }
                    if st.oversize {
                        break;
                    }
                    reducers.push(r);
                }
                Ok(None) => {
                    no_divisor.insert(m);
                }
            }
        }
        if st.timed_out || st.oversize {
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }
        if symbolic_cap_exceeded(
            &mut st,
            examined.len(),
            queue.len(),
            row_terms,
            half_rows.len() + field_rows.len() + reducers.len(),
            seen_rows.len(),
        ) {
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }
        if let (Some(lcm_bytes), Some(examined_bytes), Some(no_divisor_bytes)) = (
            lcm_columns.dense_bytes(),
            examined.dense_bytes(),
            no_divisor.dense_bytes(),
        ) {
            st.dense_symbolic_set_steps += 1;
            st.dense_symbolic_set_bytes_max = st
                .dense_symbolic_set_bytes_max
                .max(lcm_bytes + examined_bytes + no_divisor_bytes);
        }
        st.reducer_rows += reducers.len() as u64;

        let n_rows = reducers.len() + half_rows.len() + field_rows.len();
        let cols = Columns::from_monomials(examined.iter().copied(), n_vars);
        cols.record_dense_cost(&mut st);
        if deadline.is_some_and(|limit| Instant::now() >= limit) {
            st.timed_out = true;
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }
        let words = cols.words() as u64;
        if words * n_rows as u64 > MAX_MATRIX_WORDS {
            st.oversize = true;
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        }
        st.matrix_rows_max = st.matrix_rows_max.max(n_rows as u64);
        st.matrix_cols_max = st.matrix_cols_max.max(cols.monos.len() as u64);
        st.matrix_rows_sum += n_rows as u64;
        st.peak_matrix_bytes = st.peak_matrix_bytes.max(words * n_rows as u64 * 8);
        // Reducers first, each on a column of its own; then the S-rows by
        // leading monomial, so the second half of a pair meets the first.
        let mut s_rows: Vec<&F2BoolPoly> = half_rows.iter().chain(field_rows.iter()).collect();
        s_rows.sort_by(|a, b| cmp_mono(b.lt().unwrap(), a.lt().unwrap()));
        for p in reducers.iter().chain(s_rows.iter().copied()) {
            st.max_poly_degree = st.max_poly_degree.max(degree(p));
        }
        let Some(rows) = pack_rows_for_echelon(
            &cols,
            n_rows,
            reducers.iter().chain(s_rows.into_iter()),
            deadline,
        ) else {
            st.timed_out = true;
            st.build_ns += t.elapsed().as_nanos() as u64;
            s.pairs.extend(selected);
            break;
        };
        st.build_ns += t.elapsed().as_nanos() as u64;

        let t = Instant::now();
        let Some(step_echelon) = echelon(rows, cols.monos.len(), &mut st.word_xors, deadline)
        else {
            st.timed_out = true;
            st.eliminate_ns += t.elapsed().as_nanos() as u64;
            break;
        };
        st.m4ri_matrices += u64::from(step_echelon.m4ri);
        st.m4ri_block_width_max = st.m4ri_block_width_max.max(step_echelon.m4ri_block_width);
        st.m4ri_table_word_xors += step_echelon.table_word_xors;
        st.m4ri_blocks += step_echelon.m4ri_blocks;
        st.m4ri_consecutive_blocks += step_echelon.consecutive_blocks;
        st.m4ri_trimmed_word_xors_avoided += step_echelon.trimmed_word_xors_avoided;
        st.m4ri_scratch_bytes_max = st.m4ri_scratch_bytes_max.max(step_echelon.scratch_bytes);
        if let Some(order_bytes) = step_echelon.flat_order_bytes() {
            st.flat_m4ri_matrices += 1;
            st.flat_m4ri_order_bytes_max = st.flat_m4ri_order_bytes_max.max(order_bytes);
        }
        st.eliminate_ns += t.elapsed().as_nanos() as u64;

        let mut fresh: Vec<F2BoolPoly> = (0..step_echelon.len())
            .filter(|&index| no_divisor.contains(&cols.monos[step_echelon.lead(index)]))
            .map(|index| step_echelon.unpack(index, &cols, n_vars))
            .collect();
        if fresh.iter().any(is_one) {
            return finish(one(), st);
        }
        if !fresh.is_empty() {
            st.solving_degree = st.solving_degree.max(d);
        }
        fresh.sort_by(|a, b| cmp_mono(a.lt().unwrap(), b.lt().unwrap()));
        st.new_elements += fresh.len() as u64;
        let t = Instant::now();
        if batch_basis_insert_enabled() {
            s.insert_batch(fresh, &mut st);
        } else {
            for p in fresh {
                s.insert(p, &mut st);
            }
        }
        st.pair_update_ns += t.elapsed().as_nanos() as u64;
    }
    st.pairs_left = s.pairs.len() as u64;
    if st.timed_out || st.oversize {
        let basis: Vec<F2BoolPoly> = (0..s.polys.len())
            .filter(|&g| s.active[g])
            .map(|g| s.polys[g].clone())
            .collect();
        return finish(basis, st);
    }
    let active: Vec<F2BoolPoly> = (0..s.polys.len())
        .filter(|&g| s.active[g])
        .map(|g| s.polys[g].clone())
        .collect();
    let reduced = interreduce(active, s.n_vars, &mut dense_mul_scratch, &mut st);
    finish(reduced, st)
}

/// Keep one element per minimal leading monomial, then reduce every tail
/// by Gauss–Jordan on the symbolic-preprocessing matrix of that minimal
/// basis.  Returns the reduced Gröbner basis, sorted by leading monomial,
/// largest first.
fn interreduce(
    mut elements: Vec<F2BoolPoly>,
    n_vars: usize,
    dense_mul_scratch: &mut Option<DenseMulScratch>,
    st: &mut F4Stats,
) -> Vec<F2BoolPoly> {
    let t = Instant::now();
    elements.sort_by(|a, b| cmp_mono(a.lt().unwrap(), b.lt().unwrap()));
    let mut minimal: Vec<F2BoolPoly> = Vec::new();
    for p in elements {
        let l = p.lt().unwrap().mask;
        if !minimal.iter().any(|q| q.lt().unwrap().mask & !l == 0) {
            minimal.push(p);
        }
    }
    let lms: Vec<u64> = minimal.iter().map(|p| p.lt().unwrap().mask).collect();
    let mut examined: FastU64Set = lms.iter().copied().collect();
    let mut queue: Vec<u64> = Vec::new();
    for p in &minimal {
        for t in &p.terms[1..] {
            if examined.insert(t.mask) {
                queue.push(t.mask);
            }
        }
    }
    let mut reducers: Vec<F2BoolPoly> = Vec::new();
    while let Some(m) = queue.pop() {
        let mut best: Option<usize> = None;
        for (k, &l) in lms.iter().enumerate() {
            st.divisor_tests += 1;
            st.divisor_linear_tests += 1;
            if l & !m == 0 && best.is_none_or(|b| minimal[k].terms.len() < minimal[b].terms.len()) {
                best = Some(k);
            }
        }
        if let Some(k) = best {
            let r = multiply_for_f4(&minimal[k], m & !lms[k], dense_mul_scratch, st);
            for t in &r.terms {
                if examined.insert(t.mask) {
                    queue.push(t.mask);
                }
            }
            reducers.push(r);
        }
    }
    let cols = Columns::from_monomials(examined.iter().copied(), n_vars);
    cols.record_dense_cost(st);
    let mut rows: Vec<(usize, Row)> = minimal
        .iter()
        .chain(reducers.iter())
        .map(|p| {
            let mut r = cols.pack(p);
            let lead = r.lead().expect("non-zero");
            (lead, r)
        })
        .collect();
    let n_rows = rows.len() as u64;
    st.matrix_rows_max = st.matrix_rows_max.max(n_rows);
    st.matrix_cols_max = st.matrix_cols_max.max(cols.monos.len() as u64);
    st.matrix_rows_sum += n_rows;
    st.peak_matrix_bytes = st.peak_matrix_bytes.max(cols.words() as u64 * n_rows * 8);
    st.build_ns += t.elapsed().as_nanos() as u64;

    // Every row already leads on a column of its own, so the forward pass
    // is empty; the backward pass clears each pivot column from every row
    // above it, rightmost pivot first.
    let t = Instant::now();
    rows.sort_by_key(|(lead, _)| *lead);
    for k in (0..rows.len()).rev() {
        let (lead, from) = (rows[k].0, rows[k].0 / 64);
        let (bit_word, bit) = (lead / 64, 1u64 << (lead % 64));
        let (above, rest) = rows.split_at_mut(k);
        let pivot = &rest[0].1;
        for (_, row) in above.iter_mut() {
            if row.bits[bit_word] & bit != 0 {
                let to = pivot.end;
                for (a, b) in row.bits[from..to].iter_mut().zip(&pivot.bits[from..to]) {
                    *a ^= *b;
                }
                st.word_xors += (to - from) as u64;
                row.end = row.end.max(to);
            }
        }
    }
    st.eliminate_ns += t.elapsed().as_nanos() as u64;
    let lm_set: FastU64Set = lms.iter().copied().collect();
    let mut out: Vec<F2BoolPoly> = rows
        .iter()
        .filter(|(lead, _)| lm_set.contains(&cols.monos[*lead]))
        .map(|(_, r)| cols.unpack(r, n_vars))
        .collect();
    out.sort_by(|a, b| cmp_mono(b.lt().unwrap(), a.lt().unwrap()));
    out
}

/// **The solutions of a reduced boolean Gröbner basis**, without
/// enumerating `{0,1}^n`.
///
/// Every variable is either the leading monomial of a linear element
/// `v_p + Σ c_q v_q + c` (with every `v_q` free) or free itself, so an
/// assignment of the free variables fixes the point.  The number of free
/// variables is at most the number of solutions less one, since the
/// basis of a radical zero-dimensional ideal has as many standard
/// monomials as the ideal has points.  Returns `None` above `max_free`
/// free variables, which on a reduced basis means the ideal has more than
/// `2^max_free` solutions, or the input was not a reduced basis.
///
/// `tests` counts monomial tests, one per term evaluated.
pub fn solutions_from_reduced_basis(
    gb: &[F2BoolPoly],
    n_vars: usize,
    max_free: usize,
    tests: &mut u64,
) -> Option<Vec<u64>> {
    if gb.iter().any(is_one) {
        return Some(Vec::new());
    }
    let mut linear_of = vec![None; n_vars];
    for (k, g) in gb.iter().enumerate() {
        let l = g.lt()?.mask;
        if l.count_ones() == 1 {
            linear_of[l.trailing_zeros() as usize] = Some(k);
        }
    }
    let free: Vec<usize> = (0..n_vars).filter(|&v| linear_of[v].is_none()).collect();
    if free.len() > max_free {
        return None;
    }
    let pivots: Vec<(usize, usize)> = (0..n_vars)
        .filter_map(|v| linear_of[v].map(|k| (v, k)))
        .collect();
    let mut out = Vec::new();
    for a in 0u64..(1u64 << free.len()) {
        let mut point = 0u64;
        for (bit, &v) in free.iter().enumerate() {
            if (a >> bit) & 1 == 1 {
                point |= 1u64 << v;
            }
        }
        for &(v, k) in &pivots {
            // The tail involves free variables and the constant only.
            *tests += gb[k].terms.len() as u64;
            if gb[k].eval(point) == 1 {
                point |= 1u64 << v;
            }
        }
        let mut ok = true;
        for g in gb {
            *tests += g.terms.len() as u64;
            if g.eval(point) != 0 {
                ok = false;
                break;
            }
        }
        if ok {
            out.push(point);
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::pq_groebner_f2::{groebner_basis_f2, reduce, spoly};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn poly(masks: &[u64], n: usize) -> F2BoolPoly {
        F2BoolPoly::from_monos(masks.iter().map(|&m| F2BoolMono::from_mask(m)).collect(), n)
    }

    fn brute_force(eqs: &[F2BoolPoly], n: usize) -> Vec<u64> {
        (0u64..1 << n)
            .filter(|v| eqs.iter().all(|e| e.eval(*v) == 0))
            .collect()
    }

    #[test]
    fn direct_mask_order_matches_degrevlex() {
        let mut rng = StdRng::seed_from_u64(0x451f_c3a8_97d2_6be0);
        for _ in 0..100_000 {
            let a = rng.gen::<u64>();
            let b = rng.gen::<u64>();
            assert_eq!(
                cmp_mono_mask_descending(a, b),
                cmp_mono(F2BoolMono::from_mask(b), F2BoolMono::from_mask(a))
            );
        }
    }

    #[test]
    fn dense_monomial_multiply_matches_sorted_reference() {
        let mut rng = StdRng::seed_from_u64(0x37c8_d51a_24e6_90bf);
        for n_vars in 0..=20usize {
            let cap = (1u64 << n_vars) - 1;
            let mut scratch = DenseMulScratch::new(n_vars);
            for _ in 0..1_000 {
                let terms = (0..rng.gen_range(0..=96usize))
                    .map(|_| F2BoolMono::from_mask(rng.gen::<u64>() & cap))
                    .collect();
                let p = F2BoolPoly::from_monos(terms, n_vars);
                let multiplier = rng.gen::<u64>() & cap;
                assert_eq!(
                    scratch.multiply(&p, multiplier),
                    p.mul_mono(F2BoolMono::from_mask(multiplier))
                );
            }
            scratch.generation = u32::MAX - 1;
            let p = poly(&[0, 1 & cap, 2 & cap, 3 & cap, cap], n_vars);
            assert_eq!(
                scratch.multiply(&p, cap >> 1),
                p.mul_mono(F2BoolMono::from_mask(cap >> 1))
            );
        }
    }

    #[test]
    fn dense_symbolic_monomial_sets_match_hash_sets() {
        let mut rng = StdRng::seed_from_u64(0x728e_4bf1_93a5_c06d);
        let mut dense = MonomialSet::new(18);
        let mut copied = MonomialSet::new(18);
        let mut reference = FastU64Set::default();
        assert!(dense.dense_bytes().is_some());
        for _ in 0..200 {
            dense.clear();
            reference.clear();
            for _ in 0..2_000 {
                let mask = rng.gen::<u64>() & ((1 << 18) - 1);
                assert_eq!(dense.insert(mask), reference.insert(mask));
            }
            assert_eq!(dense.len(), reference.len());
            assert!(dense.iter().all(|mask| reference.contains(mask)));
            for _ in 0..200 {
                let mask = rng.gen::<u64>() & ((1 << 18) - 1);
                assert_eq!(dense.contains(&mask), reference.contains(&mask));
            }
            copied.copy_from(&dense);
            assert_eq!(copied.len(), dense.len());
            assert!(dense.iter().all(|mask| copied.contains(mask)));
        }
    }

    fn reference_new_pairs(lh: u64, lm: &[u64], active: &[bool]) -> (Vec<(usize, u64)>, u64, u64) {
        let mut c: Vec<(usize, u64)> = (0..lm.len())
            .filter(|&g| active[g])
            .map(|g| (g, lh | lm[g]))
            .collect();
        let mut d = Vec::with_capacity(c.len());
        let mut chain = 0u64;
        while let Some((g1, l1)) = c.pop() {
            let coprime = lh & lm[g1] == 0;
            let covered = c.iter().chain(d.iter()).any(|&(_, l2)| l2 & !l1 == 0);
            if coprime || !covered {
                d.push((g1, l1));
            } else {
                chain += 1;
            }
        }
        let mut product = 0u64;
        let mut selected = Vec::new();
        for (g, lcm) in d {
            if lh & lm[g] == 0 {
                product += 1;
            } else {
                selected.push((g, lcm));
            }
        }
        (selected, chain, product)
    }

    #[test]
    fn grouped_pair_selection_matches_quadratic_update() {
        let mut state = 0x7c3a_59d1_a641_2f0bu64;
        let mut scratch = PairSelectScratch::default();
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for n_vars in 1..=18usize {
            let cap = (1u64 << n_vars) - 1;
            for _ in 0..200 {
                let len = 1 + (next() as usize % 96);
                let lh = (next() & cap).max(1);
                let lm: Vec<u64> = (0..len).map(|_| (next() & cap).max(1)).collect();
                let active: Vec<bool> = (0..len).map(|_| next() & 3 != 0).collect();
                let active_indices: Vec<usize> = (0..len).filter(|&g| active[g]).collect();
                let expected = reference_new_pairs(lh, &lm, &active);
                let mut stats = F4Stats::default();
                let actual =
                    select_new_pairs(n_vars, lh, &lm, &active_indices, &mut scratch, &mut stats);
                assert_eq!(actual, expected.0);
                assert_eq!(stats.pairs_chain_skipped, expected.1);
                assert_eq!(stats.pairs_product_skipped, expected.2);
            }
        }
    }

    #[test]
    fn indexed_reducer_selection_matches_linear_scan() {
        let mut rng = StdRng::seed_from_u64(0x6bf4_1d9e_a320_57c8);
        for n_vars in 1..=18usize {
            let cap = (1u64 << n_vars) - 1;
            for _ in 0..200 {
                let len = 1 + rng.gen_range(0..192usize);
                let lm: Vec<u64> = (0..len).map(|_| (rng.gen::<u64>() & cap).max(1)).collect();
                let polys: Vec<F2BoolPoly> = lm
                    .iter()
                    .enumerate()
                    .map(|(g, &lead)| {
                        if g % 3 == 0 {
                            poly(&[lead], n_vars)
                        } else {
                            poly(&[lead, 0], n_vars)
                        }
                    })
                    .collect();
                let active_flags: Vec<bool> = (0..len).map(|_| rng.gen_ratio(3, 4)).collect();
                let active: Vec<usize> = (0..len).filter(|&g| active_flags[g]).collect();
                let state = State {
                    n_vars,
                    polys,
                    lm,
                    active: active_flags,
                    active_indices: active.clone(),
                    pairs: Vec::new(),
                    pair_select_scratch: PairSelectScratch::default(),
                };
                let active_by_lm = active_leading_monomial_index(&state, &active);
                for _ in 0..64 {
                    let monomial = rng.gen::<u64>() & cap;
                    let expected = active
                        .iter()
                        .copied()
                        .filter(|&g| state.lm[g] & !monomial == 0)
                        .min_by_key(|&g| (state.polys[g].terms.len(), g));
                    let mut stats = F4Stats::default();
                    let actual = state
                        .reducer_for(monomial, &active, &active_by_lm, &mut stats, None)
                        .unwrap();
                    assert_eq!(actual, expected);
                    assert_eq!(
                        stats.divisor_tests,
                        stats.divisor_submask_lookups + stats.divisor_linear_tests
                    );
                }
            }
        }
    }

    #[test]
    fn batched_basis_install_matches_repeated_update() {
        let n_vars = 10usize;
        let cap = (1u64 << n_vars) - 1;
        let empty_state = || State {
            n_vars,
            polys: Vec::new(),
            lm: Vec::new(),
            active: Vec::new(),
            active_indices: Vec::new(),
            pairs: Vec::new(),
            pair_select_scratch: PairSelectScratch::default(),
        };
        let mut reference = empty_state();
        let mut batched = empty_state();
        let mut reference_stats = F4Stats::default();
        let mut batched_stats = F4Stats::default();
        let mut rng = StdRng::seed_from_u64(0x918e_2d40_b73a_65cf);

        for _ in 0..40 {
            let batch_len = rng.gen_range(2..=8usize);
            let elements: Vec<F2BoolPoly> = (0..batch_len)
                .map(|_| poly(&[(rng.gen::<u64>() & cap).max(1)], n_vars))
                .collect();
            for element in elements.iter().cloned() {
                reference.insert(element, &mut reference_stats);
            }
            batched.insert_batch(elements, &mut batched_stats);
            assert_eq!(batched.lm, reference.lm);
            assert_eq!(batched.active, reference.active);
            assert_eq!(batched.active_indices, reference.active_indices);
            assert_eq!(batched.pairs, reference.pairs);
            assert_eq!(
                batched_stats.pairs_product_skipped,
                reference_stats.pairs_product_skipped
            );
            assert_eq!(
                batched_stats.pairs_chain_skipped,
                reference_stats.pairs_chain_skipped
            );
        }
    }

    fn canonical_rref(mut matrix: Vec<Vec<u64>>, n_cols: usize) -> Vec<Vec<u64>> {
        let words = n_cols.div_ceil(64);
        let mut pivot_row = 0usize;
        for column in 0..n_cols {
            let (word, bit) = (column / 64, 1u64 << (column % 64));
            let Some(found) = (pivot_row..matrix.len()).find(|&row| matrix[row][word] & bit != 0)
            else {
                continue;
            };
            matrix.swap(pivot_row, found);
            let pivot = matrix[pivot_row].clone();
            for row in 0..matrix.len() {
                if row != pivot_row && matrix[row][word] & bit != 0 {
                    for index in word..words {
                        matrix[row][index] ^= pivot[index];
                    }
                }
            }
            pivot_row += 1;
            if pivot_row == matrix.len() {
                break;
            }
        }
        matrix.truncate(pivot_row);
        matrix
    }

    #[test]
    fn m4ri_and_streaming_have_the_same_row_space() {
        let mut seed = 0x9a71_d05c_3f28_44e1u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for (n_rows, n_cols, keep_mask) in [
            (128usize, 257usize, 1u64),
            (173, 319, 3),
            (211, 385, 7),
            (257, 511, 15),
        ] {
            let words = n_cols.div_ceil(64);
            let mut raw = Vec::with_capacity(n_rows);
            for _ in 0..n_rows {
                let mut bits: Vec<u64> = (0..words)
                    .map(|_| {
                        let value = next();
                        if keep_mask == 1 {
                            value
                        } else {
                            value & next() & keep_mask.wrapping_neg()
                        }
                    })
                    .collect();
                if n_cols % 64 != 0 {
                    *bits.last_mut().unwrap() &= (1u64 << (n_cols % 64)) - 1;
                }
                raw.push(bits);
            }
            let rows = |source: &[Vec<u64>]| {
                source
                    .iter()
                    .cloned()
                    .map(|bits| Row {
                        bits,
                        start: 0,
                        end: words,
                    })
                    .collect::<Vec<_>>()
            };
            let mut streaming_ops = 0u64;
            let streaming =
                echelon_streaming(rows(&raw), n_cols, &mut streaming_ops, None).unwrap();
            let mut m4ri_ops = 0u64;
            let m4ri = echelon_m4ri(rows(&raw), n_cols, &mut m4ri_ops, None).unwrap();
            let mut flat_rows = FlatRows::new(n_rows, words);
            for (row_index, source) in raw.iter().enumerate() {
                flat_rows.row_mut(row_index).copy_from_slice(source);
            }
            let mut flat_ops = 0u64;
            let flat = echelon_m4ri_flat(flat_rows, n_cols, &mut flat_ops, None).unwrap();
            let mut streaming_columns = streaming.pivot_leads.clone();
            let mut m4ri_columns = m4ri.pivot_leads.clone();
            let mut flat_columns = flat.pivot_leads.clone();
            streaming_columns.sort_unstable();
            m4ri_columns.sort_unstable();
            flat_columns.sort_unstable();
            assert_eq!(m4ri_columns, streaming_columns, "shape {n_rows}x{n_cols}");
            assert_eq!(
                flat_columns, streaming_columns,
                "flat shape {n_rows}x{n_cols}"
            );
            assert_eq!(flat_ops, m4ri_ops, "flat charged work at {n_rows}x{n_cols}");
            let streaming_space = canonical_rref(streaming.bit_rows(), n_cols);
            let m4ri_space = canonical_rref(m4ri.bit_rows(), n_cols);
            let flat_space = canonical_rref(flat.bit_rows(), n_cols);
            assert_eq!(m4ri_space, streaming_space, "shape {n_rows}x{n_cols}");
            assert_eq!(flat_space, streaming_space, "flat shape {n_rows}x{n_cols}");
        }
    }

    fn standard_monomials(gb: &[F2BoolPoly], n: usize) -> u64 {
        let lms: Vec<u64> = gb.iter().filter_map(|p| p.lt()).map(|m| m.mask).collect();
        (0u64..1 << n)
            .filter(|&m| !lms.iter().any(|&l| l & !m == 0))
            .count() as u64
    }

    /// The Gröbner-basis certificate, checked directly: every
    /// S-polynomial and every field product reduces to zero.
    fn is_boolean_groebner_basis(gb: &[F2BoolPoly]) -> bool {
        for (i, g) in gb.iter().enumerate() {
            let l = g.lt().unwrap().mask;
            let mut bits = l;
            while bits != 0 {
                let v = bits.trailing_zeros();
                bits &= bits - 1;
                if !reduce(&g.mul_mono(F2BoolMono::from_mask(1 << v)), gb).is_zero() {
                    return false;
                }
            }
            for h in &gb[i + 1..] {
                if !reduce(&spoly(g, h), gb).is_zero() {
                    return false;
                }
            }
        }
        true
    }

    fn random_quadratic(n: usize, m: usize, rng: &mut StdRng) -> Vec<F2BoolPoly> {
        (0..m)
            .map(|_| {
                let mut monos = Vec::new();
                for i in 0..n {
                    for j in 0..i {
                        if rng.gen_bool(0.5) {
                            monos.push((1u64 << i) | (1 << j));
                        }
                    }
                    if rng.gen_bool(0.5) {
                        monos.push(1u64 << i);
                    }
                }
                if rng.gen_bool(0.5) {
                    monos.push(0);
                }
                poly(&monos, n)
            })
            .collect()
    }

    /// **The field pairs are what makes it a boolean basis.**  On
    /// `v_0v_1 + v_0 + v_1`, whose only root is `(0, 0)`, the basis is
    /// `{v_0, v_1}`: `v_0·g = v_0` and `v_1·g = v_1`.  The repository's
    /// Buchberger used to return the generator itself — three standard
    /// monomials for one solution — and was closed under the field
    /// equations after this engine was written (34154ed9); the two now
    /// agree here, and the ledger's §14–§17 rows record which engine
    /// they measured.
    #[test]
    fn the_field_pairs_close_the_toy_ideal() {
        let g = poly(&[3, 1, 2], 2);
        let (gb, st) = groebner_basis_f4(vec![g.clone()], 2, None);
        let mut lts: Vec<u64> = gb.iter().map(|p| p.lt().unwrap().mask).collect();
        lts.sort();
        assert_eq!(lts, vec![1, 2]);
        assert!(gb.iter().all(|p| p.terms.len() == 1), "{gb:?}");
        assert_eq!(standard_monomials(&gb, 2), 1);
        assert!(st.field_pairs_reduced > 0);
        let buchberger = groebner_basis_f2(vec![g], 2);
        assert_eq!(
            standard_monomials(&buchberger, 2),
            1,
            "Buchberger closes under the field equations too"
        );
    }

    /// **Certified on random systems**: the output is a Gröbner basis of
    /// the boolean ideal (every S-polynomial and field product reduces to
    /// zero), it has exactly as many standard monomials as the system has
    /// solutions, an inconsistent system gives `{1}`, and the extraction
    /// recovers every solution and nothing else.
    #[test]
    fn the_basis_is_a_certified_boolean_groebner_basis() {
        let mut rng = StdRng::seed_from_u64(20260922);
        let (mut consistent, mut inconsistent) = (0, 0);
        for trial in 0..60 {
            let n = 4 + trial % 7;
            let m = n - 2 + trial % 5;
            let eqs = random_quadratic(n, m, &mut rng);
            let (gb, st) = groebner_basis_f4(eqs.clone(), n, None);
            assert!(!st.timed_out && !st.oversize);
            let sols = brute_force(&eqs, n);
            assert!(is_boolean_groebner_basis(&gb), "trial {trial}: not closed");
            assert_eq!(
                standard_monomials(&gb, n),
                sols.len() as u64,
                "trial {trial}"
            );
            if sols.is_empty() {
                inconsistent += 1;
                assert_eq!(gb.len(), 1);
                assert!(is_one(&gb[0]));
            } else {
                consistent += 1;
            }
            let mut tests = 0;
            let mut got = solutions_from_reduced_basis(&gb, n, 24, &mut tests).unwrap();
            got.sort_unstable();
            assert_eq!(got, sols, "trial {trial}: extraction");
            // Every generator lies in the ideal the basis generates.
            for e in &eqs {
                assert!(
                    reduce(e, &gb).is_zero(),
                    "trial {trial}: a generator escaped"
                );
            }
        }
        assert!(
            consistent > 5 && inconsistent > 5,
            "{consistent} / {inconsistent}"
        );
    }

    /// **Batched Buchberger is Buchberger.**  Wherever the Buchberger
    /// output is certified a boolean basis — everywhere, since it closes
    /// under the field equations — the reduced bases must coincide
    /// polynomial for polynomial: the reduced Gröbner basis of an ideal
    /// is unique.
    #[test]
    fn f4_and_buchberger_agree_where_buchberger_is_certified() {
        use crate::cryptanalysis::ic_boundary::random_binary_instance;
        use crate::cryptanalysis::pq_descent_symbolic::descend;
        let mut rng = StdRng::seed_from_u64(5);
        let mut compared = 0;
        for &(n, np, seed) in &[(7u32, 4u32, 1u64), (9, 5, 2), (11, 6, 3)] {
            let inst = random_binary_instance(n, seed, 1 << 20).unwrap();
            let words: Vec<u64> = (0..np).map(|k| 1u64 << k).collect();
            for _ in 0..4 {
                let x_r = rng.gen::<u64>() & ((1u64 << n) - 1);
                let sys = descend(&inst.gf, inst.b, x_r, &words, 2).unwrap();
                let bb = groebner_basis_f2(sys.equations.clone(), sys.n_vars);
                if !is_boolean_groebner_basis(&bb) {
                    continue;
                }
                let (f4, _) = groebner_basis_f4(sys.equations.clone(), sys.n_vars, None);
                assert_eq!(f4, bb, "n = {n}, n' = {np}: the reduced bases differ");
                compared += 1;
            }
        }
        assert!(compared >= 8, "only {compared} systems were comparable");
    }

    /// The extraction enumerates free variables only.  The standard
    /// monomials include `1` and every free variable, and there are as
    /// many of them as solutions, so `1 + #free ≤ #solutions`.
    #[test]
    fn extraction_enumerates_only_the_free_variables() {
        let mut rng = StdRng::seed_from_u64(99);
        let mut consistent = 0;
        for _ in 0..20 {
            let eqs = random_quadratic(9, 9, &mut rng);
            let (gb, _) = groebner_basis_f4(eqs.clone(), 9, None);
            let sols = brute_force(&eqs, 9).len() as u64;
            let mut tests = 0;
            let got = solutions_from_reduced_basis(&gb, 9, 24, &mut tests).unwrap();
            assert_eq!(got.len() as u64, sols);
            if sols > 0 {
                consistent += 1;
                let free = (0..9)
                    .filter(|&v| !gb.iter().any(|g| g.lt().unwrap().mask == 1u64 << v))
                    .count() as u64;
                assert!(
                    1 + free <= sols,
                    "{free} free variables for {sols} solutions"
                );
            }
        }
        assert!(consistent > 3);
    }

    /// A budget is honoured and reported, not silently turned into a
    /// wrong basis.
    #[test]
    fn the_budget_is_reported() {
        let mut rng = StdRng::seed_from_u64(3);
        let eqs = random_quadratic(18, 17, &mut rng);
        let (_, st) = groebner_basis_f4(eqs, 18, Some(Duration::from_nanos(1)));
        assert!(st.timed_out);
    }

    #[test]
    fn symbolic_memory_cap_is_reported_as_oversize() {
        let mut st = F4Stats::default();
        assert!(!symbolic_cap_exceeded(&mut st, 100, 100, 100, 10, 10));
        assert!(symbolic_cap_exceeded(
            &mut st,
            (MAX_SYMBOLIC_BYTES / 64 + 1) as usize,
            0,
            0,
            0,
            0,
        ));
        assert!(st.oversize);
        assert!(st.symbolic_cap_hit);
        assert!(st.symbolic_bytes_estimate_max > MAX_SYMBOLIC_BYTES);
    }
}
