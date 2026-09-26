//! # Structured sparse elimination for binary Macaulay matrices.
//!
//! The solving-degree measurement in
//! [`super::koblitz_groebner::solving_profile`] is bounded by one dense
//! `F_2` row reduction.  At `n = 5, m = 3` that reduction is 8 340 rows
//! over 21 778 columns and takes about eight minutes; at 27 unknowns —
//! the `n = 9` and `n = 15` cells, which are the ones that would turn a
//! single measured gap into a scaling claim — degree 6 needs 397 594
//! columns and the same dense pass runs to roughly a fortnight per draw.
//! That is the wall, and it is a representation problem rather than a
//! mathematical one.
//!
//! ## Why not Wiedemann or Lanczos
//!
//! Those solve a sparse system or find a kernel vector.  The measurement
//! needs neither: it asks whether the reduced row space contains the
//! constant `1` (a refutation), and whether it contains `v` or `v + 1`
//! for each variable (a pinned unknown).  Those are questions about a
//! *reduced basis*, which a Krylov method does not produce.  Plain
//! Lanczos is doubly wrong here — its normal-equations step is unsound
//! over `F_2`, where a nonzero vector can be self-orthogonal, which is
//! why production GF(2) linear algebra uses *block* Wiedemann.
//!
//! So this is structured Gaussian elimination, the same family as
//! [`super::koblitz_sparse_la`]'s filtering stage, over a different
//! field and aimed at a different question.
//!
//! ## The structure that makes it cheap
//!
//! [`super::pq_groebner_f2::cmp_mono`] orders monomials by total degree
//! first, and the Macaulay builder sorts columns *descending*.  So the
//! high-degree monomials occupy the leading columns and the degree-≤1
//! monomials and the constant occupy the tail.  Rows are stored as
//! ascending column indices, so a row's leading index is its minimum —
//! and therefore
//!
//! ```text
//!     leading index ≥ the degree-≤1 boundary
//!         ⟺  the row's whole support lies in the tail
//!         ⟺  the row is a linear consequence of the system.
//! ```
//!
//! That equivalence is what lets the expensive half of the work be
//! restricted to the high columns and done sparsely:
//!
//! 1. Echelonise by leading column across the high columns only,
//!    choosing the lowest-weight row in each column's bucket as pivot to
//!    hold down fill-in.
//! 2. Whatever is left with a leading index past the boundary is exactly
//!    the set of linear consequences, and spans at most `n_vars + 1`
//!    columns.
//! 3. Reduce *that* block densely.  Echelon form alone is not enough —
//!    `{v + w, w}` has to reduce to `{v, w}` before `v` reads as pinned —
//!    but the full reduction now runs on a 65-column block instead of a
//!    65 000-column one.

use std::mem::take;

/// Outcome of eliminating the high-degree columns.
#[derive(Clone, Debug, Default)]
pub struct SparseElimination {
    /// Pivots found among the high columns.
    pub high_rank: usize,
    /// Rows whose support lies entirely in the low columns, as ascending
    /// column indices.  These are the linear consequences of the system
    /// at this degree.
    pub linear_rows: Vec<Vec<u32>>,
    /// Rows that reduced to nothing: syzygies among the Macaulay rows.
    pub vanished: usize,
    /// Largest number of nonzeros any row reached during elimination.
    ///
    /// Fill-in is the one thing that can make this approach lose to the
    /// dense path, so it is measured rather than assumed.
    pub max_weight: usize,
}

/// Symmetric difference of two ascending, duplicate-free index lists —
/// addition over `F_2`.
pub fn xor_sorted(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
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
            // Present in both: cancels.
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

/// Echelonise `rows` by leading column across columns `0 .. low_start`,
/// and return the rows left spanning only `low_start ..`.
///
/// `rows` must be ascending and duplicate-free per row; empty rows are
/// allowed and counted as vanished.
pub fn eliminate_high_columns(
    mut rows: Vec<Vec<u32>>,
    n_cols: usize,
    low_start: usize,
) -> SparseElimination {
    let mut out = SparseElimination {
        max_weight: rows.iter().map(|r| r.len()).max().unwrap_or(0),
        ..Default::default()
    };

    // Bucket every row under its leading (minimum) column.
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n_cols.max(1)];
    for (i, r) in rows.iter().enumerate() {
        match r.first() {
            Some(&lead) => buckets[lead as usize].push(i),
            None => out.vanished += 1,
        }
    }

    let mut is_pivot = vec![false; rows.len()];
    for c in 0..low_start.min(n_cols) {
        let bucket = take(&mut buckets[c]);
        if bucket.is_empty() {
            continue;
        }
        // Lowest weight wins: the cheapest row to add into the others is
        // the one that introduces the fewest new nonzeros.
        let piv = *bucket
            .iter()
            .min_by_key(|&&i| rows[i].len())
            .expect("bucket is non-empty");
        is_pivot[piv] = true;
        out.high_rank += 1;

        let pivot_row = rows[piv].clone();
        for i in bucket {
            if i == piv {
                continue;
            }
            let reduced = xor_sorted(&rows[i], &pivot_row);
            out.max_weight = out.max_weight.max(reduced.len());
            match reduced.first() {
                Some(&lead) => {
                    debug_assert!(
                        lead as usize > c,
                        "reduction must advance the leading column"
                    );
                    buckets[lead as usize].push(i);
                }
                None => out.vanished += 1,
            }
            rows[i] = reduced;
        }
    }

    // Everything still bucketed at or past the boundary spans the low
    // columns alone; so does any non-pivot row whose leading index was
    // already there.
    for (i, r) in rows.into_iter().enumerate() {
        if is_pivot[i] || r.is_empty() {
            continue;
        }
        if r[0] as usize >= low_start {
            out.linear_rows.push(r);
        }
    }
    out
}

/// [`eliminate_high_columns`] with a dense finish: the columns before
/// `sparse_until` are eliminated sparsely, exactly as there, and the rows
/// that survive them are packed as bit rows over the remaining columns and
/// put in row echelon form by
/// [`crate::cryptanalysis::koblitz_groebner::echelon_f2_counted`].
///
/// Why: in a degree-`D` Macaulay matrix the leading band (the degree-`D`
/// columns) is where rows are short, and the sparse merge wins there.
/// Once that band is gone the survivors are thousands of entries long —
/// a ladder draw at 20 unknowns spends 93 of its 100 seconds merging index
/// lists over the degree-5 and degree-4 columns — and a bit row does the
/// same row addition at a fraction of the cost.
///
/// The row space is the same, so the linear rows span the same space; the
/// refutation and the pinned variables read off it, and `high_rank`, are
/// unchanged.  The dense finish also reduces the linear rows among
/// themselves, so they come back independent; `vanished` keeps the identity
/// `rows = high_rank + vanished + linear rows` and therefore counts the
/// dependent linear rows that the sparse path lists.  `max_weight` covers
/// the sparse phase only.  The linear rows themselves
/// are a different basis of that space, which is all
/// [`crate::cryptanalysis::koblitz_groebner::solving_profile_sparse`] reads.
pub fn eliminate_high_columns_dense_finish(
    mut rows: Vec<Vec<u32>>,
    n_cols: usize,
    low_start: usize,
    sparse_until: usize,
) -> SparseElimination {
    use crate::cryptanalysis::koblitz_groebner::echelon_f2_counted;

    let n_rows = rows.len();
    let sparse_until = sparse_until.min(low_start).min(n_cols);
    let mut out = SparseElimination {
        max_weight: rows.iter().map(|r| r.len()).max().unwrap_or(0),
        ..Default::default()
    };

    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n_cols.max(1)];
    for (i, r) in rows.iter().enumerate() {
        if let Some(&lead) = r.first() {
            buckets[lead as usize].push(i);
        }
    }
    let mut is_pivot = vec![false; n_rows];
    for c in 0..sparse_until {
        let bucket = take(&mut buckets[c]);
        if bucket.is_empty() {
            continue;
        }
        let piv = *bucket
            .iter()
            .min_by_key(|&&i| rows[i].len())
            .expect("bucket is non-empty");
        is_pivot[piv] = true;
        out.high_rank += 1;
        let pivot_row = rows[piv].clone();
        for i in bucket {
            if i == piv {
                continue;
            }
            let reduced = xor_sorted(&rows[i], &pivot_row);
            out.max_weight = out.max_weight.max(reduced.len());
            if let Some(&lead) = reduced.first() {
                buckets[lead as usize].push(i);
            }
            rows[i] = reduced;
        }
    }
    drop(buckets);

    // Survivors span columns `sparse_until ..` only.  Pack them, freeing
    // each sparse row as it goes.
    let width = n_cols - sparse_until;
    let words = width.div_ceil(64).max(1);
    let mut dense: Vec<Vec<u64>> = Vec::new();
    for (i, r) in rows.iter_mut().enumerate() {
        if is_pivot[i] || r.is_empty() {
            continue;
        }
        let mut bits = vec![0u64; words];
        for &c in r.iter() {
            let k = c as usize - sparse_until;
            bits[k / 64] |= 1 << (k % 64);
        }
        *r = Vec::new();
        dense.push(bits);
    }
    drop(rows);

    let mut word_ops = 0u64;
    let rank = if dense.is_empty() {
        0
    } else {
        echelon_f2_counted(&mut dense, width, &mut word_ops)
    };
    let low_offset = low_start - sparse_until;
    for row in dense.iter().take(rank) {
        let lead = row
            .iter()
            .enumerate()
            .find(|(_, &w)| w != 0)
            .map(|(k, &w)| k * 64 + w.trailing_zeros() as usize)
            .expect("an echelon row within the rank is nonzero");
        if lead < low_offset {
            out.high_rank += 1;
            continue;
        }
        let mut support = Vec::new();
        for (k, &w) in row.iter().enumerate() {
            let mut w = w;
            while w != 0 {
                let b = w.trailing_zeros() as usize;
                support.push((sparse_until + k * 64 + b) as u32);
                w &= w - 1;
            }
        }
        out.linear_rows.push(support);
    }
    out.vanished = n_rows - out.high_rank - out.linear_rows.len();
    out
}

/// Where [`eliminate_high_columns_dense_finish`] should stop eliminating
/// sparsely: the end of the leading degree band.  Columns are in
/// descending, degree-first monomial order, so that band is a prefix.
pub fn leading_band_end(cols: &[u64]) -> usize {
    match cols.first() {
        Some(first) => {
            let top = first.count_ones();
            cols.iter()
                .position(|m| m.count_ones() < top)
                .unwrap_or(cols.len())
        }
        None => 0,
    }
}

/// Index of the first column whose monomial has degree at most 1.
///
/// Columns are in descending monomial order and [`cmp_mono`] is
/// degree-first, so the degree-≤1 monomials form a suffix and this is
/// the boundary the elimination targets.
///
/// [`cmp_mono`]: super::pq_groebner_f2::cmp_mono
pub fn low_column_start(cols: &[u64]) -> usize {
    cols.iter()
        .position(|m| m.count_ones() <= 1)
        .unwrap_or(cols.len())
}
