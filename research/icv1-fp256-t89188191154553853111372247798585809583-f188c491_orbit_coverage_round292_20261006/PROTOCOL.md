# P-256 serpentine-orbit coverage, round 292: protocol

Status: preregistered before implementation or execution on 6 October 2026.

## Question and hypothesis

Round 291 emits every state of a table-free 17-term serpentine sign orbit at a
counted construction-stage cost of `0.9686171464465044` times its frozen rho
boundary.  That count treats every emitted state as fresh coverage.  Round 292
tests the missing condition exactly: how many distinct P-256 group elements do
those walks actually cover?

The candidate hypothesis is that at least `0.9686171464465044` of emitted
states are distinct, so charging the Round 291 work to distinct coverage still
does not exceed the stage boundary.  The falsification hypothesis is that
nearest-neighbour path revisits reduce the distinct fraction below this value.

This is a selector-stage coverage experiment.  Passing it would not establish
a relation probability, structured degree, row collection, sparse linear
algebra, scalar recovery, or an end-to-end improvement over Pollard rho.

## Frozen identity and dependencies

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- comparison factor base: `FB1h2f8621cda105`;
- family: Round 33 / Round 291 known-log two-delta cyclic coefficient base;
- columns: `131458`; arity: `17`; rare edges: `6935`; cutoff: `219`;
- P-256 subgroup order is taken from the repository's native curve catalog;
- Round 291 result path:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_serpentine_sign_orbit_round291_20261006/serpentine-result.json`;
- Round 291 result SHA-256:
  `068be64f508ea2fb2996c265426dd8213088d1cad905859470792fcbc80972da`;
- Round 291 semantic evidence SHA-256:
  `453d14bc5a9e68ab5c6ba9d73ffa9040bc3c8557be7d73b0eb42d5395d939710`;
- coefficient SHA-256:
  `980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53`;
- frozen Round 291 conservative stage ratio: `0.9686171464465044`;
- tuple selection is imported verbatim from Round 291.  The primary tuple has
  attempt `240`; holdouts have attempts `1648`, `1766`, and `1920`.

The program must reject any dependency, tuple, capacity, coefficient digest,
or frozen ratio mismatch before producing a result.

## Exact interval construction

For one sign assignment, the path changes the scalar coefficient by exactly
`+1` or `-1` at every step.  Therefore all visited integer offsets form the
complete interval from the minimum prefix sum to the maximum prefix sum.  Add
the signed starting coefficient modulo the P-256 subgroup order.  Split only
when the resulting inclusive interval wraps through zero.

The union of these inclusive modular intervals is the exact support of the
emitted group elements because the P-256 generator has order `n`.  No hash or
probabilistic fingerprint may decide equality.  Intervals are sorted and
merged using exact 256-bit integer comparisons.

Report separately:

1. emitted states, including repeated visits;
2. the sum of per-segment interval widths, removing within-segment revisits;
3. exact union width, additionally removing cross-segment overlap;
4. within-segment, cross-segment, and total duplicate counts;
5. distinct fraction and Round 291 charged cost per distinct state;
6. the corrected ratio to the frozen stage boundary;
7. maximum interval width, interval count, wrap splits, operations, wall time,
   peak RSS, and deterministic support digest.

## Cells and controls

- Exhaustive toy controls use at least two prime cyclic groups.  Directly
  enumerate every path state and require exact equality with interval-union
  support and exact duplicate counts.
- Native primary depths are `8, 10, 12, 14, 16, 17` sign bits.
- Native portfolio cells use the first `1, 2, 3, 4` frozen tuples, all at the
  complete 17-bit sign depth.  Cumulative union support is measured exactly.
- At native tractable depths, direct scalar-set enumeration must equal the
  interval result.  Full-depth controls compare emitted counts, boundary
  coefficients, and deterministic digests with independently recomputed
  values; no probabilistic equality test is admissible.
- Any reported relation would require exact native group replay.  No relation
  search or relation claim is authorized by this protocol.

## Boundary and accounting

The primary boundary table uses one unit: frozen Round 291 charged P-256
addition equivalents per distinct covered state, normalized to the same rho
stage boundary.  Let `f = distinct / emitted`.  The corrected stage ratio is

```text
0.9686171464465044 / f.
```

The table must include the uncorrected Round 291 row, every native depth, and
every cumulative tuple portfolio.  Allocation, sorting and interval merging
are also counted and timed separately; they may not be hidden from a claim of
a complete selector.  No stage number may be described as end-to-end speed.

## Gates and stop conditions

The coverage hypothesis passes only if all of the following hold:

- zero false positives and false negatives on every complete direct control;
- exact support and duplicate-count equality on both exhaustive toys;
- no dependency or arithmetic failure;
- full-depth primary distinct fraction at least `0.9686171464465044`;
- full-depth primary corrected construction-stage ratio at most `1`;
- cumulative four-tuple corrected construction-stage ratio at most `1`;
- materialized storage below `2^50` bytes.

If any coverage-ratio gate fails, publish the exact negative result and name
within-segment versus cross-segment duplication as measured.  Do not tune the
tuple selection, omit unfavourable signs, extrapolate discarded branches as an
exhaustive search, or attempt a full-depth unplanted P-256 relation.

The inherited end-to-end promotion gates remain failed/unset unless separately
established: same-family structured residual degree at most 5, relation
collection below `2^120`, cost per usable relation below `2^103`, 138,031
independent rows with duplicate/rank allowance, sparse linear algebra, and a
demonstrably non-generic complete DLP.

## Reproducibility and decision rule

The implementation is native Rust.  The run is deterministic and isolated on
a pinned core.  Commit the source, dependency hashes, exact JSON result,
isolation receipt, artifact hashes, tests, one-unit comparison table, dashboard
update, and decision.  Preserve failures instead of overwriting them.

If the exact coverage gates pass, the next round may integrate collision
detection and measured relation predicates.  If they fail, Round 291 remains a
local construction result only and the measured duplicate factor becomes the
new obstruction; do not claim parity from its uncorrected emitted-state rate.
