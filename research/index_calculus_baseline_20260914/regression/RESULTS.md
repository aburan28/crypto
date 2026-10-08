# Full-corpus benchmark results

Class: **accounting**. Frozen target: [frozen_targets.json](frozen_targets.json),
derived from [baseline_v2/summary.json](results/baseline_v2/summary.json) and
[raw.jsonl](results/baseline_v2/raw.jsonl). Both configurations use the identical
pinned WDSat binary and input corpus. No solver algorithm was changed.

| Variant | Full-DLP S | Cost / rho | Cost / floor | Correctness per 60-input sweep | Class |
|---|---|---|---|---|---|
| Symmetry on | Unmeasured | Unmeasured | Unmeasured | 30 relations, 29 UNSAT, 1 algebraic-only rejection | accounting |
| Symmetry off | Unmeasured | Unmeasured | Unmeasured | 30 relations, 29 UNSAT, 1 algebraic-only rejection | accounting |
| Matched rho reference | Not run | Unmeasured | Unmeasured | No full-DLP recovery in this suite | reference required |

The normalized operation metric and both boundary ratios remain unmeasured.
The support-counting floor in the parent report applies to uniform-target
attempt counts, not to this selected corpus; no counter-to-operation conversion
or matched rho calibration is supplied here. The existing attack verdict and
exponents are unchanged.

The following are **supplementary solver diagnostics**. Each entry sums the
per-input median of three repetitions, rather than summing all repetitions.

| Configuration | n | l | Inputs | Sum of median conflicts | Sum of median process seconds |
|---|---:|---:|---:|---:|---:|
| symmetry_on | 15 | 5 | 20 | 80,314 | 0.589849 |
| symmetry_on | 17 | 6 | 20 | 614,275 | 5.545287 |
| symmetry_on | 19 | 6 | 20 | 585,692 | 5.362547 |
| symmetry_off | 15 | 5 | 20 | 386,438 | 2.600279 |
| symmetry_off | 17 | 6 | 20 | 3,036,969 | 28.848458 |
| symmetry_off | 19 | 6 | 20 | 2,823,074 | 27.480456 |

| Configuration | All-input conflict target | All-input process seconds (secondary) |
|---|---:|---:|
| symmetry_on | 1,280,281 | 11.497683 |
| symmetry_off | 6,246,481 | 58.929193 |

The symmetry-off/on conflict ratio is 4.878992
on this fixed corpus. It is not a ratio of full-attack operations. Three repeated
processes per input all produced identical statuses and conflict counts within
each configuration. The timing host reports AMD EPYC 9V74; it was shared, with
no CPU pinning. Complete per-case medians, p90 values and raw repetitions remain
in the linked files. No runtime distribution or exponent was fitted.

The v2 batch completed 360 processes in 244.837 seconds of loop
wall time: 211.688 seconds in solver processes,
0.599 seconds in SAT certification, and
32.432 seconds in independent exhaustive UNSAT checks.
The remaining loop time is orchestration. This excludes compilation, downloads,
and initial input preparation. The oracle checked 929,280 unordered triples,
with repeated summands allowed. No whole-DLP time is implied.

## Admission correction retained

The first batch used a strict point-only admission rule and completed 360
processes, with 354 passing that rule and six rejected. Its complete failed
record remains under [strict_point_v1](results/strict_point_v1). The six rejected
runs all concern `n19l6-19-U`: ANF and S4 are satisfiable, but the target does not
lift to this curve over the stated base field. The absolute trace certificate
is in [negative_control.json](negative_control.json).

The v2 contract designated this case as an algebraic-only rejection control
**before rerunning all 360 processes**. V2 certifies 180 SAT point-relation
results, 174 UNSAT results, and six algebraic-only rejections. The latter always
contribute zero point relations. Neither batch had a process timeout or nonzero
solver exit. These inputs do not estimate natural decomposition probability.

## Target for the next iteration

Run the commands in [README.md](README.md#run-every-iteration-against-this-target).
Keep both baseline configurations immutable and compare every candidate on the
same 60 inputs. The optional diagnostic target is at least 20% fewer summed
median conflicts in each configuration, with no input more than 10% worse.
Passing it does not establish an operation-normalized improvement; that requires
the calibrated all-phase accounting and independent holdouts in the parent
contract. Algebraic-only solutions must never become claimed point relations.
