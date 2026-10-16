//! Per-thread tallies of the `bitop/v1` primitives an engine performs
//! beyond the count it reports natively.
//!
//! An algebraic engine's native unit (`word XORs (elimination only)`, say)
//! deliberately counts one kind of work so that committed records stay
//! comparable.  Everything else it does — building and packing Macaulay
//! matrices, sorting and hashing monomials, searching for pivots, reading
//! rows back, the splitting tree's substitutions, building and
//! specialising a summation polynomial — is tallied here, in the
//! primitives of `docs/units/BITOPS.md`, and reported beside the native
//! count through the `bitop_*` keys of `SolverCost::extra`
//! ([`crate::cryptanalysis::ecbench::work::SOLVER_BITOP_KEYS`]).
//!
//! The granularity is one counted unit per inner-loop iteration that
//! touches a matrix word or a monomial mask:
//!
//! * **word op** — one 64-bit test, XOR, OR, compare or move on a packed row
//!   or a monomial mask (a monomial over at most 64 variables is one word);
//!   a comparison of a sort counts one, and a sort of `k` keys is charged
//!   `k · ⌈log₂ k⌉` ([`sort_cost`]);
//! * **hash probe** — one insertion into or lookup in a hashed monomial
//!   set or map;
//! * **field add / mul / sqr / inv** — one operation in `F_{2^n}`.
//!
//! The tallies are thread-local, like the koblitz engine's word-op counter,
//! so a reader takes the difference of two [`snapshot`]s on its own thread.
//! Work shared out to rayon workers is charged by the calling thread from
//! the sizes it hands out, never from inside the workers, so every count is
//! a deterministic function of the input.

use std::cell::Cell;
use std::collections::BTreeMap;

/// What one thread has tallied so far.  Counts only grow; readers subtract.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    pub word_ops: u64,
    pub hash_probes: u64,
    pub field_adds: u64,
    pub field_muls: u64,
    pub field_sqrs: u64,
    pub field_invs: u64,
}

impl Tally {
    /// The work done between `before` and `self`.
    pub fn since(self, before: Tally) -> Tally {
        Tally {
            word_ops: self.word_ops.wrapping_sub(before.word_ops),
            hash_probes: self.hash_probes.wrapping_sub(before.hash_probes),
            field_adds: self.field_adds.wrapping_sub(before.field_adds),
            field_muls: self.field_muls.wrapping_sub(before.field_muls),
            field_sqrs: self.field_sqrs.wrapping_sub(before.field_sqrs),
            field_invs: self.field_invs.wrapping_sub(before.field_invs),
        }
    }

    /// Write the whole run into `extra` as the `bitop_*` keys: `native`
    /// word operations (the engine's own count, whatever its unit's name)
    /// plus this tally, and `bitop_complete = 1`.
    pub fn put_complete(self, extra: &mut BTreeMap<String, u64>, native_word_ops: u64) {
        let mut put = |k: &str, v: u64| {
            extra.insert(k.to_string(), v);
        };
        put(
            "bitop_word_ops",
            native_word_ops.wrapping_add(self.word_ops),
        );
        put("bitop_hash_probes", self.hash_probes);
        put("bitop_field_adds", self.field_adds);
        put("bitop_field_muls", self.field_muls);
        put("bitop_field_sqrs", self.field_sqrs);
        put("bitop_field_invs", self.field_invs);
        put("bitop_complete", 1);
    }
}

thread_local! {
    static WORD_OPS: Cell<u64> = const { Cell::new(0) };
    static HASH_PROBES: Cell<u64> = const { Cell::new(0) };
    static FIELD_ADDS: Cell<u64> = const { Cell::new(0) };
    static FIELD_MULS: Cell<u64> = const { Cell::new(0) };
    static FIELD_SQRS: Cell<u64> = const { Cell::new(0) };
    static FIELD_INVS: Cell<u64> = const { Cell::new(0) };
}

#[inline]
fn bump(c: &'static std::thread::LocalKey<Cell<u64>>, n: u64) {
    c.with(|c| c.set(c.get().wrapping_add(n)));
}

/// Charge `n` word operations to this thread.
#[inline]
pub fn word_ops(n: u64) {
    bump(&WORD_OPS, n);
}

/// Charge `n` hash probes to this thread.
#[inline]
pub fn hash_probes(n: u64) {
    bump(&HASH_PROBES, n);
}

/// Charge `F_{2^n}` operations to this thread.
#[inline]
pub fn field_ops(adds: u64, muls: u64, sqrs: u64, invs: u64) {
    bump(&FIELD_ADDS, adds);
    bump(&FIELD_MULS, muls);
    bump(&FIELD_SQRS, sqrs);
    bump(&FIELD_INVS, invs);
}

/// This thread's tallies so far.
pub fn snapshot() -> Tally {
    Tally {
        word_ops: WORD_OPS.with(Cell::get),
        hash_probes: HASH_PROBES.with(Cell::get),
        field_adds: FIELD_ADDS.with(Cell::get),
        field_muls: FIELD_MULS.with(Cell::get),
        field_sqrs: FIELD_SQRS.with(Cell::get),
        field_invs: FIELD_INVS.with(Cell::get),
    }
}

/// Comparisons charged for sorting `k` keys: `k · ⌈log₂ k⌉`, the
/// comparison count of a balanced merge sort (a pattern-defeating
/// quicksort makes about as many on these inputs).
#[inline]
pub fn sort_cost(k: usize) -> u64 {
    if k < 2 {
        return k as u64;
    }
    let k = k as u64;
    k * (64 - (k - 1).leading_zeros()) as u64
}

/// Counts gathered in a hot loop and charged to this thread when dropped,
/// so every early return of the loop's function is charged too.
#[derive(Debug, Default)]
pub struct Pending {
    pub word_ops: u64,
    pub hash_probes: u64,
}

impl Pending {
    /// One lookup of a monomial in a column index: a word operation on a
    /// direct table, a hash probe on a hashed one.
    #[inline]
    pub fn lookup(&mut self, hashed: bool, n: u64) {
        if hashed {
            self.hash_probes += n;
        } else {
            self.word_ops += n;
        }
    }
}

impl Drop for Pending {
    fn drop(&mut self) {
        word_ops(self.word_ops);
        hash_probes(self.hash_probes);
    }
}

/// Word operations of the Boolean-polynomial primitives of
/// [`crate::cryptanalysis::pq_groebner_f2::F2BoolPoly`], from the term
/// counts, so a caller charges them without that module counting anything.
pub mod poly {
    use super::sort_cost;

    /// `p · m` on `k` terms: one OR a term, the re-sort, one compare a term
    /// for the cancelling scan.
    #[inline]
    pub fn mul_mono(k: usize) -> u64 {
        2 * k as u64 + sort_cost(k)
    }

    /// `from_monos` on `k` monomials: the sort and the cancelling scan.
    #[inline]
    pub fn from_monos(k: usize) -> u64 {
        k as u64 + sort_cost(k)
    }

    /// `substitute` on `k` terms: a test a term and a compare a term of the
    /// merge of the kept and folded runs.
    #[inline]
    pub fn substitute(k: usize) -> u64 {
        2 * k as u64
    }

    /// `p + q`: one compare a term of the merge.
    #[inline]
    pub fn add(a: usize, b: usize) -> u64 {
        (a + b) as u64
    }

    /// `eval` (or any single scan) on `k` terms: one AND-and-compare a term.
    #[inline]
    pub fn scan(k: usize) -> u64 {
        k as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tallies_are_per_thread_differences() {
        let before = snapshot();
        word_ops(5);
        hash_probes(2);
        field_ops(1, 2, 3, 4);
        let d = snapshot().since(before);
        assert_eq!(
            d,
            Tally {
                word_ops: 5,
                hash_probes: 2,
                field_adds: 1,
                field_muls: 2,
                field_sqrs: 3,
                field_invs: 4,
            }
        );
        std::thread::spawn(|| assert_eq!(snapshot(), Tally::default()))
            .join()
            .unwrap();
        let mut extra = BTreeMap::new();
        d.put_complete(&mut extra, 10);
        assert_eq!(extra["bitop_word_ops"], 15);
        assert_eq!(extra["bitop_complete"], 1);
    }

    #[test]
    fn sort_cost_is_k_log_k() {
        assert_eq!(sort_cost(0), 0);
        assert_eq!(sort_cost(1), 1);
        assert_eq!(sort_cost(2), 2);
        assert_eq!(sort_cost(8), 24);
        assert_eq!(sort_cost(9), 36);
    }
}
