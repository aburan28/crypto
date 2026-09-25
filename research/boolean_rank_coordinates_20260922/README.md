# Combinatorial monomial coordinates without an exponential lookup

The previous turn made verified progress: packed direct construction passed a
2x cold-construction gate twice, but used a table indexed by every n-bit mask.
This experiment removes that exponential table and tests larger generic Boolean
matrices. The previous files and measurements remain unchanged. The complete
algorithm's speed remains unproven; this is another construction-stage test.

## Rank formula and contract

Let `B(b,r) = sum(binomial(b,j), j=0..r)`. This counts b-bit masks having at most
r set bits. For a monomial mask m of degree k <= D, list its set-bit positions
in decreasing order `h_1 > ... > h_k`. Its zero-based index in ascending numeric
order among all monomials of degree at most D is

`rank_D(m) = sum(B(h_i, D-i+1), i=1..k)`.

Proof: partition the smaller masks by the highest bit at which they differ
from m. That bit must be one of the h_i, changed from one to zero. The i-1 higher
ones remain fixed, leaving at most D-i+1 ones in h_i lower positions. These
classes are disjoint and contain every smaller allowed mask exactly once.
The constant monomial has rank zero. Pascal's recurrence constructs B using
O(nD) table entries; evaluation visits at most D set bits.

Every product supplied to the ranker satisfies the precondition: the selected
multiplier t has degree at most D-deg(p), so every `m OR t` has degree <= D.
Input validation and the exact `(n,D,active)` context guard remain mandatory.
Changed contexts rebuild directly. Duplicate or out-of-range input monomials
are rejected. Cancellation is XOR, empty rows are omitted, actual occupied
columns are compacted in numeric order, and current output caps apply after
compaction. Row order stays generator order then numeric multiplier order.

This removes only the `2^n` lookup. The ambient basis still has
`sum(binomial(n,j),j=0..D)` columns, and temporary packed rows use that width.
The experiment charges these costs and tests sparse restricted-mask cases that
could expose their overhead.

## Fixed experiment and falsification

The source and `protocol.json` are frozen before the run. There are 144 cells:
n=12/20/28/36, three families, batches 1/8/32, two discovery seeds, two new
holdout seeds and twelve repetitions. Each system has eight generators, each
with at most 27 possible terms. All batch inputs are distinct.

The families are quadratic coefficients, full-mask linear degree drops, and
a restricted eight-variable multiplier mask cycling through quadratic, linear,
constant and zero generators. A full-mask degree-3 constant at n=36 requires
7,807 rows by itself, exceeding the common 4,096-row budget; its explicit
resource refusal is tested. It is not silently interpreted as a zero matrix.

Three scalable arms share each input:

- **sorted:** direct product sorting, parity cancellation and compact packing;
- **binary:** packed rows with binary search in the ambient basis;
- **ranked:** the same packed constructor with combinatorial coordinate ranks.

A **dense** lookup control runs only at the n=12 bridge. It is a u64-monomial
port of the earlier approach, under the new matched workload. The old timing
numbers are not cross-run baselines. Dense controls above n=12 are deliberately
NOT_EXECUTED with null costs; an omitted run is not a measured failure.

For dramatic scalable-construction promotion, at batch32 on **all three
families** and n=20/28/36, the 95% paired-bootstrap lower bound must exceed
**2.0 against sorted** and **1.05 against binary**. All 18 comparisons must pass.
The n=12 dense comparisons are separately reported. Failed families cannot be
removed or retitled after timing. No tuning uses the new holdouts.

The reference boundary is the best correct measured constructor of the same
matrices. Returning the output words remains necessary work. All cold setup,
applications, allocations, support compaction, exact validation and destruction
are timed. Fixture/oracle generation is common supplied-input work outside arm
totals and inside process receipts. Capacity bytes exclude allocator metadata;
whole-worker RSS includes all arms and references. Intervals describe fixed
seeds and repetitions, not an independently sampled population.

The seven Rust tests independently check the rank against numeric enumeration
through n=12 and combinations at n=20/28/36, all 6,144 small coefficient/degree/
mask cases across every constructor, high variable bits, current caps, context
changes and resource refusals. The oracle uses recursive multiplier generation,
set parity and column-major monomial-to-row incidence, separately from the
candidate's iterative multipliers, rank calculation and row-major XOR.

## Run and replay

```sh
python3 research/boolean_rank_coordinates_20260922/run.py --out research/boolean_rank_coordinates_20260922/run_01
python3 -m unittest discover -s research/boolean_rank_coordinates_20260922 -p 'test_*.py'
```

The runner retains source/configuration hashes, executable hashes, raw fixtures,
measurements, process receipts, test output and a manifest. Replays operate on
temporary copies, preserving frozen runs. This worker constructs generic
matrices only. No curve inputs, solver integration or key-related work occurs;
the WDSat/full IC suite is inapplicable. Full-solving and rho costs stay null.

## Additive sparse-intermediate follow-up

The first run removed the exponential lookup and improved every larger-size
comparison against binary search, but only 11/18 promotion gates passed. The
remaining ambient-width intermediate rows were a specific untested cost. The
frozen `run_01` evidence is retained without changes.

`sparse_protocol.json` predeclares a successor with holdouts 20260928/155921.
It retains the sorted, binary and ranked constructors as fresh controls and
adds `sparse_rank`. Each application uses one scratch bitmap and a marked-word
list. Only touched nonzero words survive into intermediate rows. The scratch
and markers reset even after cancellation, and only the final output matrix is
allocated at the actual compact width. Context, rank, parity, row ordering and
cap rules are unchanged. The final dense-output cost is still charged.

There are twenty balanced repetitions of four arms, or five at the n=12 bridge.
All **27** larger-size gates must pass: a lower bound above 2.0 against sorted,
and above 1.05 against each of binary and ranked, for every family and size at
batch32. No tuning uses these new holdouts. The eighth Rust test specifically
covers cancellation/reactivation inside a word and resetting scratch state.

```sh
python3 research/boolean_rank_coordinates_20260922/run.py --protocol research/boolean_rank_coordinates_20260922/sparse_protocol.json --out research/boolean_rank_coordinates_20260922/run_02
```

The current worker accepts an optional `with-sparse` argument that the runner
adds only for the successor protocol. Each frozen run retains the exact worker,
runner, analyzer and protocol it executed. The successor changes the algorithm;
it is a fresh experiment, not a confirmation of the original implementation.
