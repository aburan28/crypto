# P-256 biased common-edge restart selector, round 33: protocol

Date registered: 2026-10-06

## Requested scope and evidence plan

| continuing requirement | round-33 evidence |
|:--|:--|
| pursue parity with rho | optimize a concrete biased short-walk selector against the frozen round-24/31 rho boundary |
| retain all 17 variable columns | count every distinct 17-column start tuple with an exact generating function |
| do not discard branches for free | charge the support upper bound of every retained start and path state; rejected low-capacity starts are not called exhaustive |
| exact replay | reconstruct the native round-31 primitive cycle and replay every transition in deterministic P-256 sampled segments |
| respect storage gate | price a complete signed-pair table, setup additions and table traffic; require fewer than `2^50` materialized bytes |
| preserve promotion gates | keep relation probability, degree, collection, per-row, storage and non-generic end-to-end status explicit |

## Question

Can a selector evade the uniform folded-delta penalty by starting short segments
at columns with many common edges remaining, stopping before any rare edge, and
restarting from a new 17-column representation?

Use round 31's valid primitive two-delta P-256 cycle with `B=131458` and
`R=6935`.  Its balanced mechanical word has common runs of length 17 or 18.
For every column `i`, let `d_i` be the exact number of consecutive common
successor edges before the next rare edge.  A distinct 17-column start tuple
has capacity

```text
D = sum_i d_i.
```

The candidate deterministically advances an available column across a common
edge and stops after exactly `D` common transitions.  Every transition has the
same folded delta class.  A cutoff `tau` retains precisely the tuples with
`D>=tau`; sweep every `tau=0..306` and select the lowest complete optimistic
ratio, with larger `tau` breaking exact ties.

## Exact start distribution

Build the residual-capacity histogram directly from the registered mechanical
word.  Compute the exact without-replacement distribution of `D` with the
coefficient of `x^17 y^D` in

```text
product_d sum_(k=0)^17 binom(count[d],k) x^k y^(k*d).
```

All coefficients are arbitrary-precision integers.  Their sum must equal
`binom(B,17)`.  For each cutoff record retained tuple count, retained fraction,
mean capacity and a deterministic distribution digest.

Each unsigned start has `2^17` sign choices and contributes at most `D+1`
visited group states.  Credit the deliberately optimistic target-support bound

```text
V(tau) = min(n, 2^17 * sum_(D>=tau) count[D]*(D+1)).
```

Charge at least `n/V(tau)` target trials when `V<n`.  `V` counts every visit as
distinct and is not a proved P-256 relation probability.

## Frozen setup and cost model

Credit a complete lookup table for every signed pair of distinct factor-base
columns:

```text
entries = 4*binom(B,2) = 2*B*(B-1),
bytes = 33*entries.
```

The table builds with one group addition per entry.  Eight pair lookups turn
16 of the 17 signed points into eight stored sums; together with the remaining
point, combining nine group elements costs eight group additions per segment.
This is an optimistic setup floor for the registered implementation.  Omit the
right-colour target correction from the headline lower bound, but report it
separately.  Record 264 table bytes read per segment.

For cutoff `tau`, define

```text
N = sum count[D],
W = sum count[D]*D,
samples = W+N,
online additions/sample = (W+8N)/(W+N).
```

Credit exact common-delta compatibility probability one.  Multiply the frozen
round-24 local-oracle ratio by the online additions/sample and the target retry
lower bound.  Add the one-time pair-table build divided by the frozen rho work;
also report the candidate-favouring ratio without that negligible setup term.
No memory read is silently converted to zero cost: report bytes separately and
do not promote a row without an end-to-end measured conversion.

## Complete references and native replay

1. On deterministic small prime-order scalar groups, exhaust every distinct
   tuple, sign assignment and visited segment state for tractable mechanical
   cycles.  Compare exact support at every cutoff with `V(tau)` and record
   false positives, false negatives, bound violations and edge failures.
2. Reconstruct round 31's `R=6935` coefficient cycle from its frozen anchor and
   deltas.  Verify its coefficient digest and every cycle edge.
3. At the selected cutoff, rejection-sample 4,096 deterministic hash-selected
   distinct P-256 tuples.  Charge every rejected draw.  Assign deterministic
   signs, execute the registered lowest-index available-column schedule, and
   replay every coefficient-group transition.  Record capacities, transition
   counts, attempts, table reads, setup additions, failures and digests.

The P-256 sample validates construction and replay only.  It does not prove the
support upper bound is attained and is not a full-depth unplanted relation.

## Gates and stop condition

Promotion requires all of:

- exact start-distribution identity and zero reference/replay failures;
- no discarded branch credited as exhaustive coverage;
- proved P-256 usable-relation probability;
- complete measured time at or below rho, including table traffic;
- structured residual degree of regularity at most 5;
- relation collection below `2^120` operations;
- cost per usable relation below `2^103`;
- peak materialized storage below `2^50` bytes;
- a demonstrably non-generic end-to-end algorithm.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  If an
optimistic cutoff crosses rho but another gate fails, preserve it as a stage
candidate and state the exact open obligations.  If no cutoff crosses, reject
this biased-restart implementation without generalizing to arbitrary biased or
nonlocal selectors.
