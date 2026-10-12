# Specialized linear-tail workload with retained controls

The preceding goal turn produced useful negative evidence: full-RREF streaming
and incidence caches lost once complete reduction was charged. Its source audit
also established that the specialized consumer already computes a smaller
mathematical output. This study targets that output directly and retains the
existing flat-tail and sparse-bucket strategies as controls. It does not claim
that merely restricting the output is a new optimization.

## Exact output and invariants

Let W be the row space of every nonzero degree-D product admitted by the current
multiplier mask and caps. Let L contain only the variable and constant columns.
Return the canonical affine basis of **W intersect L**, encoding x_i in bit i
and the constant in bit n. This is the full linear-tail kernel output; the
consumer's later contradiction/forced-variable filtering and recursive solving
are outside this experiment. Source row count, source column count and high-
block rank must also agree. No arm may stop early when it finds a constant or
fills the affine space. All nonzero source rows consume the row cap, even when
they are dependent.

The controls use descending degree-reverse-lex order, putting nonlinear columns
before the linear/constant tail. The flat control retains exact-support column
caching, fused cached packing, complete fallback and the existing flat high-
column elimination loop. Both controls cache multiplier schedules within a
cold batch; doing so for every mask strengthens the source control. The second
control retains the existing minimum-weight sparse-bucket elimination strategy.
`SOURCE_PORTS.json` pins the audited source files and describes port differences.

The streamed candidates generate parity into two parts: sorted nonlinear
monomial keys, and one affine word. High pivots are indexed by their leading
key. Eliminating high terms updates both parts. A row with no high terms enters
a small affine reducer. Distinct high pivots are independent in the high-block
projection, so they cannot form a nonzero combination supported wholly in L;
the high-zero residues therefore span exactly W intersect L.

The plain candidate keeps the first pivot. The exchange candidate may replace
it with a lighter same-leading row, then continues with the old XOR new residual.
Exchanging **both** the high vector and affine word and retaining that residual
preserves the row space. Neither implementation needs a globally indexed matrix,
but both track original support exactly to honor resource caps.

## Inputs and predeclared gate

`protocol.json` is frozen before timing. There are 192 cells: n=12/20/28/36,
four families, batches 1/4/8, two discovery seeds and two new holdout seeds.
Twelve repetitions rotate four arms evenly. Every batch contains distinct
coefficient assignments, with no performance screening.

The active mask equals variables actually occurring in each input. Quadratic
and linear-drop families use the nominal variable range. Restricted-cycle
fixtures use eight evenly spaced variables within that range and cycle through
quadratic, linear, constant and zero generators. Cross-cancel fixtures make two
quadratic generators differ by a variable and optional constant, ensuring a
nontrivial affine consequence without supplying a linear generator. These
fixtures address empty and nonempty tails and changing cache contexts.

For **each candidate separately**, at batch8 on every family and n=20/28/36,
the 95% paired-bootstrap lower bound for the faster retained control / candidate
must exceed **2.0** in all twelve comparisons. Both existing controls participate
in the reference minimum. No failed family may be removed after observation,
and no tuning uses the new holdouts. This is a finite kernel gate, not an
end-to-end polynomial-solving or cryptanalytic claim.

## Accounting and verification

Cold totals include schedule/layout creation, exact cache guards and misses,
products and cancellations, source-support accounting, high elimination, affine
canonicalization, output validation and destruction. Cache state survives only
within the batch. Reference preparation is common supplied-input work outside
arm totals and inside fresh-worker process receipts. Retained schedule-mask
and layout-column counts are counts, not invented byte measurements. Packed
word XORs and sparse merge items are different phase diagnostics and are not
treated as a common total-operation unit. Inapplicable counters stay null.

An independent recursive/set-parity constructor forms the full matrix. A full
Gauss-Jordan oracle then extracts its canonical affine intersection. Every
answer, source dimension and high rank must match. Seven current Rust tests additionally
cover all 24,456 small binary matrix/split cases by explicit span enumeration,
6,144 small Boolean coefficient/degree/mask cases, ordering equivalence, larger
fixtures, hidden affine consequences, cache/cap changes, cancellations and
invalid inputs. The evidence analyzer independently reconstructs exact source
support and cached-layout hits from the retained polynomial inputs.

The reference boundary is the fastest correct retained kernel strategy on the
same complete workload. Full polynomial solving, degree-ladder iteration, root
extraction, IC costs and rho ratios remain unmeasured. No production solver
path is changed or run, no curve/target input is accepted, and the WDSat/full
IC suite is inapplicable to this standalone mathematical kernel study.

## Run and replay

```sh
python3 research/boolean_linear_tail_20260923/run.py --out research/boolean_linear_tail_20260923/run_01
python3 -m unittest discover -s research/boolean_linear_tail_20260923 -p 'test_*.py'
```

Each run uses a new directory with frozen source/configuration, executable
hashes, raw fixtures and measurements, process receipts, test output and a
manifest. Replays use temporary copies. All prior experiments remain unchanged.

## Ordering correction and source-matched run

The first run is retained as a numeric-order diagnostic. Its inherited worker
sorted polynomial terms and multiplier masks numerically. The source orders
terms by descending DegRevLex and multipliers by degree layers, then combination
order within each layer. The row space is identical, but pivot selection and
timing can differ. The first run is therefore not promoted as source-order
timing evidence, regardless of its gate results.

The current worker and protocol match both source orders. A seventh regression
pins explicit multiplier sequences and rejects the wrong input-term order.
`run_02` uses new holdouts 20261002/393241 and the unchanged performance gate.
It is an additive corrected-port experiment, not a same-worker confirmation.
The predecessor worker is hash-bound in the new protocol, and both frozen
source bundles remain replayable. To rerun either historical configuration,
use that run's frozen `run.py` with a new output directory. The top-level
`run.py` uses the current source-order protocol.

```sh
python3 research/boolean_linear_tail_20260923/run.py --out research/boolean_linear_tail_20260923/run_02
```

## Exact census and fixed density dispatch

The source-order run passes seven of twelve 2x cells for pivot exchange, but
sparse streaming loses badly on eight-active-variable systems. A third,
additive experiment keeps all four earlier arms and introduces two changes
under `census_protocol.json`, with new holdouts 20261003/458879.

`stream_census` replaces only the original-support hash set with an exact
bitmap. For a weight-k monomial whose set bits are `b_1 < ... < b_k`, its
within-weight index is `sum(C(b_i,i), i=1..k)`. Add the sizes of the higher-degree
layers to obtain a unique degree-ordered index. Cached binomial entries use
O(nD) storage; one temporary bitset uses `sum(C(n,j), j=0..D)` bits. Population
count gives the **exact** original column count. This does not weaken the
resource or output contract, and no full ambient matrix is constructed.

`hybrid_census` uses the identical cached flat path when the possible degree-D
monomial universe over the active variables has at most **512 columns**
(eight packed u64 words per potential row). Otherwise it uses `stream_census`.
The threshold is fixed before the new holdouts, not selected per measured cell.
On these fixtures active variables equal occurring variables. Both branches
remain mathematically valid on every supported input.

The original per-cell 2x gate remains reported and unchanged. A **separate
portfolio gate** measures one cold batch8 for every one of the twelve
n>=20/family cells, summed for each paired seed/repetition. It requires a 95%
lower bound above 2.0 against the pointwise best retained flat/bucket controls
and above 1.05 against the pointwise best of all four prior arms. It also
requires every size/family cell, including n=12, to have a lower bound at least
0.95 against the best prior arm. No cell is omitted. These fixed benchmark
weights are not measured solver-call frequencies and cannot establish whole-
solver performance.

Nine current Rust tests include exact census indexing/cardinality and both
dispatch branches. All earlier coefficient, source-cap and output checks apply
to both added arms. The current worker accepts `with-census` only when the
runner is given the corresponding protocol. Prior frozen runs remain intact.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/census_protocol.json --out research/boolean_linear_tail_20260923/run_03
```

The attempted `run_03` failed before its first measurement because the launcher
appended the optional mode twice. `LAUNCH_FAILURE.json` binds the untouched
failure files; `run_04` uses the identical worker/protocol with the corrected
command builder. A regression pins the exact command argument list.

The completed `run_04` passes both aggregate ratio thresholds but fails several
small-path non-regression checks, so its portfolio gate remains REJECTED. The
dispatcher made an extra recursive call before executing the selected method.
`dispatch_protocol.json` fixes a follow-up that normalizes the method once,
keeps the same algorithms, 512-column threshold and every gate, and uses fresh
holdouts 20261004/524309. The earlier result is not overwritten or relaxed.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/dispatch_protocol.json --out research/boolean_linear_tail_20260923/run_05
```

## Predecessor-balanced timing and identical-control check

The nonrecursive dispatcher still failed the small-path guard in `run_05`.
Rotating start positions preserved fixed neighbors: on small cases the flat
control usually followed the identical hybrid-flat path, while the hybrid
followed a sparse method. This creates a possible first-order carryover bias.
The earlier data are retained as observations under that schedule.

`balanced_protocol.json` keeps the algorithms, threshold and performance gates,
uses fresh holdouts 20261005/655373, and adds an identical `flat_alias` control.
Six Euler circuits over the seven directed arm labels, including self edges,
give every method every predecessor exactly six times. One validated, untimed
primer establishes the first predecessor. Its cost is in process receipts.
Every measured invocation still creates and destroys a fresh algorithm context.
This is cold algorithmic setup, not a claim of flushed hardware caches.

There are 42 observations per arm per cell. The occurrence index of each arm
defines paired repetitions. The alias must match exact work/cache counters and
have a 95% timing interval wholly within [0.95,1.05] on all sixteen batch8
size/family cells. This additional gate is required for portfolio promotion.
It tests the timing setup; it is not an additional candidate or a second
opportunity to pick a faster baseline. First-order balance does not guarantee
the absence of every higher-order hardware effect.

Eleven current Rust tests include schedule balance and alias identity. The
worker accepts the additional `balanced-order` argument only with the census
arm set. The command builder's exact-argument regressions cover both modes.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/balanced_protocol.json --out research/boolean_linear_tail_20260923/run_06
```

The `run_06` alias-equivalence gate also fails. `typed_protocol.json` declares a
follow-up with fresh holdouts 20261006/786433. Experimental labels now resolve
to a typed kernel before timing: both flat labels select the identical enum
value and timed entry point. The hybrid's real density decision remains timed.
The paired repetition index is `cycle * arm_count + predecessor`, matching
control and candidate under the same predecessor and cycle. Algorithms,
workloads, dispatch threshold and every gate are unchanged.

The current twelve Rust tests include typed alias identity and exact matched-
predecessor pairing. The command builder distinguishes `matched-order` from
the earlier occurrence-paired `balanced-order` protocol.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/typed_protocol.json --out research/boolean_linear_tail_20260923/run_07
```

The short-window `run_07` still does not satisfy all alias and non-regression
intervals. `long_window_protocol.json` retains every threshold and uses fresh
holdouts 20261007/917519. Each observation now contains sixteen independent
cold batches, each with a new context, full validation and destruction. Raw
times/work counts sum those batches; reports normalize by sixteen. The primer
observation uses the same sixteen-batch structure and stays outside measured
sample counts. This increases timing duration without warming algorithmic
caches across batches. Both flat labels also share one non-inlined leaf.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/long_window_protocol.json --out research/boolean_linear_tail_20260923/run_08
```

`run_08` passes every identical-control check and both aggregate ratios, but one
small-case non-regression interval still fails. `cutoff_protocol.json` uses
fresh holdouts 20261008/1048583 and replaces dispatch arithmetic with equivalent
integer comparisons. For degree 2 the boundary is 31 active variables
(`497 <= 512 < 529`); for degree 3 it is 14 (`470 <= 512 < 576`). This preserves
the original 512-column policy exactly. A thirteenth Rust test checks every
active-variable count 0–64 and degree 0–3. All timing and promotion thresholds
remain unchanged.

```sh
python3 research/boolean_linear_tail_20260923/run.py --protocol research/boolean_linear_tail_20260923/cutoff_protocol.json --out research/boolean_linear_tail_20260923/run_09
```
