# P-256 factor-base screening, executed rounds 2–290: protocol

Status: preregistered before implementation or execution on 6 October 2026.

## Purpose and round semantics

This experiment answers the explicit request for 289 actual rounds after the
existing Round 1 baseline.  It executes screening rounds `2..=290` inclusive;
every round constructs and validates a distinct full P-256 factor base and
evaluates a complete `2^17` sign orbit.  A formula-only cell, discarded branch,
or inherited result does not count as an executed round.

These are fresh, uniformly defined screening rounds.  They do not retroactively
turn the 172 projection-only cells in the historical rounds 43–290 radix
portfolio into executions, and they do not claim to rerun those infeasible
registered sizes.  The historical Round 290 formula requires roughly 300 TB
and remains correctly labelled unexecuted.

## Frozen identity and dependencies

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- requested comparison factor base: `FB1h2f8621cda105`;
- existing Round 1 factor-base receipt SHA-256:
  `bf4c0f95ca46bf92b236eea7b8f37a3902f744b07a898a41f6b120faab2e35a4`;
- existing Round 1 degree receipt SHA-256:
  `efe17e77cd23affe03b7b2ee8146490409afbd779f694058606caf798909c3be`;
- Round 292 exact-coverage receipt SHA-256:
  `19d546fd2e53728f029d349080edbf7337ac2ada785446978b91a62eabbe3821`;
- factor-base columns `B = 131458`, relation arity `17`, sign depth `17`,
  capacity cutoff `219`, and local ideal-oracle ratio
  `0.964336477130181`;
- the P-256 subgroup order and generator come from the native curve catalog.

The program must hash-check every dependency and reject any identity, curve,
or frozen-field mismatch before emitting a round receipt.

## The 289 factor bases

Round `r` in `2..=290` uses

```text
R(r) = 6647 + 2(r - 2)
```

rare edges, covering all 289 odd values from 6647 through 7223.  The common
edge delta is one.  The rare delta is the exact solution of

```text
(B - R) + R * delta_rare = 0 mod n.
```

Because every registered `R` is coprime to `B`, each candidate is expected to
form one primitive cycle.  Each round must nevertheless construct all 131,458
coefficients, replay every edge, prove distinctness and uniqueness up to sign,
exclude the identity, and record a deterministic coefficient digest.  Anchors
use the existing Round 31 derivation domain and the first valid attempt.  The
factor-base identifier follows its `LD2R{R}h{digest-prefix}` storage convention.

Each round independently hash-selects distinct columns and accepts the first
tuple whose exact residual common-edge capacity is at least 219.  The accepted
attempt, all 17 columns, capacities, and digest are part of that round's JSON
receipt.

## Sign-oriented monotone serpentine candidate

Round 292 found that ordinary signed paths backtrack.  This sweep reverses the
coordinate direction when its sign is negative:

- a positive term traverses its coefficient from the original point to the
  common-edge endpoint;
- a negative term traverses from the endpoint back to the original point.

Consequently every step changes the signed group sum in the same direction.
Even Gray ordinals traverse forward and odd ordinals traverse backward.  A
consecutive Gray mask changes one sign, so the segment boundary remains one
addition using a precomputed `P(original) + P(endpoint)` value.

Every segment is therefore a monotone inclusive modular interval of exactly
`capacity + 1` group elements.  Exact 256-bit interval union measures
cross-sign overlap without hashes or probabilistic equality.  This construction
must be checked by exhaustive direct enumeration on two toy groups and by
direct scalar-set controls at registered shallow native depths before the 289
full-depth rounds are accepted.

## Required receipt for every round

Each of the 289 round artifacts must record:

- `execution_status = complete`, wall and CPU time, and peak RSS;
- exact factor-base validation and edge-replay counts;
- factor-base ID, coefficient and tuple digests, anchor attempt, and all
  selected columns;
- 131,072 sign segments, emitted states, exact distinct states, duplicates,
  distinct fraction, interval counts, comparisons, and support digest;
- path, boundary, initialization, precomputation, target-correction, and total
  charged P-256 addition equivalents;
- the ideal-local corrected stage ratio and the R-specific hidden-collision /
  finite-support lower bound;
- zero reported relations and no full-depth unplanted relation attempt.

A top-level manifest must bind the path, byte count, and SHA-256 of every
artifact and prove that each integer round 2–290 occurs exactly once.

## Boundary and cost accounting

The primary table uses charged P-256 addition equivalents per exact distinct
covered element, normalized to rho.  For one round with capacity `C`, `S=2^17`
segments, and exact union width `D`, charge

```text
S*C + (S-1) + 16 + 17 + S
```

additions: path steps, sign boundaries, initial sum, the 17 boundary-delta
precomputations, and target corrections.  No sorting, interval merging,
factor-base construction, or validation work may be omitted from the measured
resource table, although the addition-equivalent ratio is already a lower
bound before converting those host operations.

For rare fraction `q=R/B`, charge the Round 31 hidden-collision factor
`1/sqrt((1-q)^2+q^2)` and its finite-support retry lower bound.  Report both the
ideal-local ratio and the stricter factor-base-complete lower bound.  Neither
is an end-to-end ratio: polynomial solving, relation verification, collection,
rank loss, sparse linear algebra, and scalar recovery remain additional.

## Correctness and promotion gates

The sweep passes its execution gate only if:

- exactly 289 artifacts cover rounds 2–290 once;
- every artifact says complete and has a verified primitive factor base;
- all toy and native direct controls have zero false positives and negatives;
- every manifest hash and byte count replays;
- no projection is counted as an execution;
- peak materialized storage is below `2^50` bytes.

A candidate passes the selector gate only if its exact distinct fraction is at
least `0.9686171464465044` and both its ideal-local and R-specific corrected
ratios are below one.  Advancement toward the original attack additionally
requires same-family structured residual degree at most 5, projected cost per
usable relation below `2^103`, 138,031 independent rows below `2^120` including
duplicate/rank allowance and sparse linear algebra, and a demonstrably
non-generic complete method below rho.

If no candidate passes, publish the best measured round and the dominant
obstruction as a reproducible negative result.  Do not attempt an unplanted
full-depth P-256 relation and do not claim improvement from local degree or an
emitted-state count.

## Reproducibility

Use native Rust for construction, arithmetic, controls, execution, result
generation, and verification.  Run the complete sweep once under the native
isolated benchmark wrapper on a pinned CPU.  Preserve failures and reruns as
additional isolation records.  Commit the protocol before implementation, then
commit source, all 289 compact receipts, manifest, isolation evidence, result
report, dashboard update, and tests in a stacked pull request.
