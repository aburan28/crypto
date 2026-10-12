# N83 recorded-query sampler contract mismatch

- **Status:** `open`
- **Kind:** implementation or test-contract mismatch; cause unknown
- **Affected work:** [conductor-gap validation](../research/ecbench_conductor_gap_20261012/VALIDATION.md)
- **Next action owner:** library sampler maintainer; unassigned
- **Last updated:** 2026-10-12

## Problem

A green whole-library validation cannot be established on the pinned baseline
because the N83 query sampler's recorded first coefficient differs from the
reference stream expected by its regression test.

The fixture is `KoblitzCurve::known_n83_k0()`, its signed generator orbit,
target `[7]G`, seed 2026100901, m = 2, enumeration, eight maximum trials and
`allow_direct_relation = false`. The acceptance condition is full-width,
zero-inclusive additive query draws matching the specified deterministic
stream, with the recorded point equal to `[a]G + [b]Q`.

## Evidence and uncertainty

- **Observed:** the individual release test fails at
  `src/cryptanalysis/koblitz_index_calculus.rs:15515` on original baseline
  e5a8c7b8cb66cc1d3c332f6ac80c7492482fde41. Recorded first `a` is
  899991076855513638720967; the test expects 1842027667052925756179639.
  [Individual reproducer log](../research/ecbench_conductor_gap_20261012/library-n83-baseline-retry.log).
- **Derived or proved:** no mathematical obstruction established.
- **Suspected explanation:** unknown; RNG consumption order, sampler version
  and test expectations require tracing rather than assumption.
- **Missing capability:** none identified for reproducing the assertion.
- **Limits:** this says nothing about N83 hardness or the validity of the
  separately audited prime-field conductor-gap sessions. No experiment kernel
  was changed by the conductor-gap work.

## What would unblock it

- [ ] Trace the reference and producer streams for the exact eight-query
  fixture, retaining seeds, draw order, coefficients and point checks.
- [ ] Resolve the implementation/test contract without removing the full-width
  or zero-inclusive requirement; version any intentional stream change and
  preserve historical replay behavior.
- [ ] Pass the named regression and relevant query-law tests with a retained
  receipt. This resolves the sampler validation issue, not broader N83 research.

## Next discriminating step

- **Action:** trace the first differing draw and any RNG use before relation
  generation; compare with the test's stated stream.
- **Controls:** unchanged fixture and seed; point-check each recorded query;
  do not replace a failed expectation with a copied producer value.
- **Limits:** one fixture, at most eight attempts, one worker, 60-second
  watchdog; retain any infrastructure refusal or timeout.
- **Outcomes:** unintended extra draws require a producer correction with
  replay compatibility; an intentionally different law needs an explicit
  versioned contract and independent verification. An unexplained discrepancy
  leaves this record open.
- **Recheck trigger:** a sampler correction or explicit stream-version record.

## Resolution and history

2026-10-12: opened from the full-library run and an individual reproducer.
Unresolved. No sampling-law or scientific status transition is made here.
