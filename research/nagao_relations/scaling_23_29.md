# Hybrid decomposition beyond eleven bits

The image-space hybrid resolved all 48 target/mode slots at 11, 23 and 29 bits, including every full enumeration, within five seconds. The slowest took 3.391237 seconds. Across 240 new trials no validation errors were recorded. The matched S3/S4 controls resolved 0/64 larger-field slots under the same budget. This is a finite-panel engineering improvement, not an ECDLP breakthrough.

## What changed

1. Cache L_V per instance; move the fixed-b norm outside the inner loop; batch inversions of root-dependent constants.
2. Replace the cubic remainder test by exact membership in the image of w²+uw on V. Retain preimages during Gaussian elimination so the support test also recovers both residual roots.
3. Charge every image-space setup, candidate, duplicate and verification. Check actual deadline compliance and preserve timeouts as unknown.

Proofs and frozen source contracts are in solver_07, solver_08 and solver_09. The image-space successor was selected after profiling the cache variant, so this is an adaptive development comparison; a later held-out panel is still necessary.

## Matched completion

Every table cell is resolutions within five seconds / four attempted target slots. U is uniform affine sampling. D is sampled from random signed factor triples and is biased toward supported targets with more decompositions. Strata may overlap. First mode includes verified UNSAT; enumerate mode requires the full projected set.

### first

| Variant | n11,d5 U | n11,d5 D | n23,d6 U | n23,d6 D | n29,d6 U | n29,d6 D | Class | S / rho |
|---|---:|---:|---:|---:|---:|---:|---|---|
| quadratic-original | 4/4 | 4/4 | 0/4 | 2/4 | 0/4 | 2/4 | Engineering | Unmeasured |
| quadratic-optimized | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 4/4 | Engineering | Unmeasured |
| quadratic-image | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | Engineering | Unmeasured |
| s4-symmetric | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| chained-s3 | 0/4 | 1/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |

### enumerate

| Variant | n11,d5 U | n11,d5 D | n23,d6 U | n23,d6 D | n29,d6 U | n29,d6 D | Class | S / rho |
|---|---:|---:|---:|---:|---:|---:|---|---|
| quadratic-original | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| quadratic-optimized | 4/4 | 4/4 | 4/4 | 3/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| quadratic-image | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | Engineering | Unmeasured |
| s4-symmetric | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| chained-s3 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |

## Counted work on matched complete enumerations

This table uses one unit: field API operations (additions + multiplications + squarings, inversions expanded). All four targets in each eleven-bit stratum completed, so work is comparable. This unit excludes Python control and coordinate conversion, and assigns uncalibrated equal weight to different primitives. Raw files retain the separate operation vectors.

| Variant | Uniform operations | Ratio to original | Supported operations | Ratio to original |
|---|---:|---:|---:|---:|
| quadratic-original | 1,037,611 | 1.000 | 1,043,058 | 1.000 |
| quadratic-optimized | 535,750 | 0.516 | 540,318 | 0.518 |
| quadratic-image | 301,229 | 0.290 | 302,090 | 0.290 |

The image solver uses about 71% fewer counted field API operations than the original on these completed panels. Ratios for larger incomplete baseline enumerations are deliberately absent; comparing a full run with a timeout would change the amount of work being solved.

## Preliminary timing for verified relations

Seconds per verified unique projected relation, including every attempted target and all its setup/failures. The larger uniform strata have zero relations, so their per-relation cost and ratio are undefined. Timeouts make these budgeted diagnostics, not uncapped expected costs.

| Field | Stratum | Mode | Image solver seconds / relation | Symmetric S4 | Chained S3 |
|---:|---|---|---:|---:|---:|
| 11 | uniform | first | 0.140421 | 1.076124 | No verified relation |
| 11 | uniform | enumerate | 0.205191 | 2.354431 | No verified relation |
| 11 | known_decomposable | first | 0.097652 | 0.617215 | 15.362349 |
| 11 | known_decomposable | enumerate | 0.126566 | 1.369089 | 20.079616 |
| 23 | uniform | first | No verified relation | No verified relation | No verified relation |
| 23 | uniform | enumerate | No verified relation | No verified relation | No verified relation |
| 23 | known_decomposable | first | 0.822202 | No verified relation | No verified relation |
| 23 | known_decomposable | enumerate | 2.252873 | No verified relation | No verified relation |
| 29 | uniform | first | No verified relation | No verified relation | No verified relation |
| 29 | uniform | enumerate | No verified relation | No verified relation | No verified relation |
| 29 | known_decomposable | first | 1.167391 | No verified relation | No verified relation |
| 29 | known_decomposable | enumerate | 3.167307 | No verified relation | No verified relation |

## Correctness and sampling limits

The cached candidate sequence and complete projected sets match the original
on all 43 affine five-bit targets. The image-space predicate matches the
previous cubic predicate on all 2,798 generated candidate occurrences across
those targets, including rejections. The chained control also passes complete
oracle equality on all 43. Larger outputs are checked against an independent
curve-group pair table. That pair oracle itself agrees with the earlier
signed-triple oracle on the full five-bit panel. All source hashes, 240 raw
trial aggregates and cross-variant target/mode lists were verified.

All eight uniform larger-field targets have no restricted decomposition.
Their resolution is verified nonexistence, not successful relation yield.
The supported samples demonstrate recovery but cannot estimate natural yield.
There are only 64 possible factor abscissas at d=6 before curve membership
and exclusions, while the ambient field grows to 2^29. This controlled test
shows field-degree scaling with a small base; it does not establish practical
relation collection with factor bases sized for a full ECDLP attack.

The current hybrid still enumerates O(2^(2d)) conditioned branches. A crucial
next test increases d to 7 and 8, rather than only increasing n, on fresh
seeds with preregistered budgets and identical Semaev controls. The image
setup must remain charged. Repeated runs, calibrated SAT operation counts,
relation dependence, matrix cost and final recovery remain outside this
campaign. The original three-size 20% all-cost goal against the strongest
Semaev baseline is therefore still open; no generic counting floor moved.

## Publication

The original PR #316 merged before several later experiment commits. This
follow-up includes those already-published research dependencies and their
unchanged frozen evidence, plus solver_07/08/09. Current main was merged;
its unrelated work is preserved. Publication records map original run IDs
to exact source snapshot trees before that integration.
