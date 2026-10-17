# Conductor-gap CM orbit experiments

## Request and boundary

Execute the controls, 100-pair orbit-folded comparison, replication and
representation/floor/family controls, and larger-prime-gap construction work
proposed in the 15 October research summary. Preserve the historical results;
land native sources, frozen inputs, unsuccessful runs, independent replays,
visual reports and a PR. A planning PR alone does not discharge execution.

The primary population is the existing 100 distinct prime gaps
33,555,419–35,282,963, with two explicit 53-bit endpoints per gap. Each
execution solves one public target, with no shared target work. The generic
reference is `rho.negation` on that same point. The negation generic reference
floor is sqrt(pi*r/4) group additions; hash-table and non-group arithmetic are
also retained in their native units. An unpriced cost leaves a complete
operation-price or speedup unknown, rather than zero. No exponent, sub-rho,
production-target, or intrinsic-hardness claim is a success requirement.

## Frozen experiments

1. **Orbit construction.** All 100 pairs, four unknown seed logs, 2,048
   levels, three arms: crater/Velu, the identical crater/scalar-lambda, and
   floor/scalar-lambda. Twenty shuffled paired repetitions plus A/A controls.
   Require identical point/coefficient digests in the two crater arms.
   Measure native arithmetic counts, cold setup and repeated construction;
   record earned isolation and do not upgrade contended wall observations.
2. **Complete single-target solves.** Retain all 200 original public points
   as a historical-regression cohort. Four seed columns, 8,192 abscissae,
   signed MITM2, rejection of scalar-proportional column parts, exact modular
   elimination and independent scalar replay. Add crater/scalar-lambda and
   same-point rho controls. Cap each target at 1,000,000,000 point queries
   and 1,800 seconds; retain all unsuccessful rows. No target can read its
   validation scalar. Reusable setup and the five exclusive target-dependent
   phases are recorded separately; setup is also included in supplementary
   cold totals. Code, inputs, target law, method order and executable hashes
   are pinned before the campaign, in a versioned run manifest.
3. **Replication.** Use the first ten pairs in the frozen shuffled cohort,
   two fresh public targets per curve, and five alternating query-seed
   rounds on the original pilot. Include degree-one coordinate changes,
   two additional certified floor vertices per selected gap, and a second
   admissible class-number-one CM family. Add same-point rho references.
   These are independent single-target workloads, not shared-table averages.
4. **Structural versus generic control.** The Velu map and scalar action
   must agree on the crater. Compare its two implementations before
   comparing different models or conductors. Keep attempted/dependent rows,
   useful relations, rank, queries, lookups and per-query arithmetic costs.
   MITM results are not prime-field F4/SAT complexity measurements. A new
   algebraic-backend claim is conditional on a separate cross-checked backend;
   none is presumed by this protocol.
5. **Larger gaps.** Implement/test a native short-Weierstrass square-root
   Velu codomain construction against all 100 known codomains and a direct
   small-kernel reference, then attempt rational-torsion-family constructions
   through 96 field bits. Independently certify primes, point orders,
   nonsingularity, j, trace and conductor classification for each successful
   endpoint. The older near-limit, non-rational-kernel candidate remains a
   separate requirement; a new rational-kernel family does not silently
   supply its floor.

## Mathematics and accounting

Canonical scalar action is lambda=(1+(2-t)/f_pi)/2 mod r. Select the sign of
the rational self-isogeny so that its actual subgroup action equals this
value, not just a root of its polynomial. This differs from the historical
pilot's sign choice and is a versioned native candidate, not a rewrite.
Seed selection is target-independent. Scalar-generated orbit points remain
available on the floor, so a map-stage advantage alone is not an ECDLP gain.

The main family is D_K=-19. D_K=-11 cannot give both prime field order and
prime group order in this odd-conductor construction (one is divisible by
3); use D_K=-43 as the second admissible family and retain the obstruction.
Attach checked Lean lemmas where feasible, with exact hypotheses and a proof
map; arithmetic lemmas alone do not formalize the entire CM classification.

Counts and target recovery are the primary reproducibility requirements.
Wall/CPU values are admitted only within their earned measurement tier.
Freeze A/A observations, affinity, environment, budgets and order; record
native independent replay outside the arithmetic window when separate.
Report per-target online IC/rho costs, supplementary cold costs, and paired
curve-level intervals. Never drop censored targets to construct a speedup.
The 100-pair population has one CM family and initially one floor per gap.

## Coordination and implementation

Base revision: e15e29d (crypto/main, checked before changes). New execution
and research analysis are native Rust, with an existing compiled native
backend where required; historical Python sources/results are immutable
provenance only. Python may render already-derived plots/PDFs as allowed by
the repository's research-visuals skill, not perform research computation.

`conductor check` was attempted before code changes: the CLI was initially
absent; the verified release was installed and returned `not logged in`.
There is no configured Conductor credential in this executor, so no shared
lease is claimed. Work is on a new isolated branch and confined to the new
study, its native plug-in, tests, and the required identity/report indexes.
The preflight is an access/tool limitation, not a mathematical obstruction.

## Acceptance and reporting

- Every executed target has a status, exact inputs/seed, accounting boundary,
  native-work counters, phase costs, and recovery certificate or stop reason.
- Same-crater implementations produce identical orbit records.
- New floor constructions pass independent native arithmetic certificates;
  a stratum or discriminant alone is not an explicit endpoint.
- Results, controls, negatives, native-source hashes and verification receipts
  land beside a source-linked report, SVG/PNG and PDF. Mirror performance
  visuals under `benchmarks/` and update affected indexes and dashboards.
- Track every A–E requirement as complete, partial or blocked with its exact
  evidence. Passing tests, opening or merging a PR is not research completion.

Inference is conditional on the executed configuration and population.
Bounded failure, a faster known-scalar map, or a single favorable timing
does not establish a general conductor-dependent ECDLP hardness difference.
