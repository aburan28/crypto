# Function-first goal: a completed eleven-bit panel, not yet a breakthrough

Across 288 new trials, the direct quadratic function solver resolved every first-relation and full-enumeration slot at 5, 9 and 11 bits. All answers were checked against exhaustive group arithmetic. Exact same-base S4 controls now show why completion alone is insufficient: elementary S4 also resolves every first-relation slot and beats the function solver at smaller sizes. The three-size 20% total-cost acceptance criterion remains unmet.

## Goal status

1. Modular-squaring support circuit: implemented and validated in solver_02.
2. Eleven-bit completion: all eight development target slots resolved in both modes by the hybrid quadratic solver. This is finite-panel success, not population reliability or success for the pure root-free SAT circuit. Two uniform eleven-bit targets have no restricted relation.
3. Scaling and strongest baseline: matched S3 and two S4 implementations now measured at three sizes. Calibrated cross-solver operation accounting, wider held-out samples and larger fields remain open.

## Completion tables

Each entry is resolved / four attempted target slots. First mode means a verified relation or verified restricted UNSAT. Enumerate mode requires complete projected-set equality; an incomplete run is not promoted because it appears to have found all tuples. Uniform and known-decomposable strata remain separate.

### first

| Variant | n5 U | n5 D | n9 U | n9 D | n11 U | n11 D | Class | S / rho |
|---|---:|---:|---:|---:|---:|---:|---|---|
| chained-s3 | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 2/4 | Engineering | Unmeasured |
| rr-norm | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 3/4 | Engineering | Unmeasured |
| fixed-u | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 1/4 | Engineering | Unmeasured |
| quadratic-function | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | Engineering | Unmeasured |
| s4-symmetric | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | Engineering | Unmeasured |
| s4-transformed | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 0/4 | Engineering | Unmeasured |

### enumerate

| Variant | n5 U | n5 D | n9 U | n9 D | n11 U | n11 D | Class | S / rho |
|---|---:|---:|---:|---:|---:|---:|---|---|
| chained-s3 | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| rr-norm | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 3/4 | Engineering | Unmeasured |
| fixed-u | 4/4 | 4/4 | 2/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| quadratic-function | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | Engineering | Unmeasured |
| s4-symmetric | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | Engineering | Unmeasured |
| s4-transformed | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 0/4 | Engineering | Unmeasured |

## Preliminary cost comparison

All timings below are seconds per verified unique projected relation, including every attempted target in the stratum and its setup, failures, rejected/duplicate models, extraction and verification. The ratio uses the lowest measured Semaev control value in the same stratum and mode. This is a censored-budget diagnostic, not calibrated operation cost or an estimate of uncapped completion cost. A zero-yield cell has no ratio.

| Mode | Field | Stratum | Quadratic seconds / relation | Best measured Semaev seconds / relation | Diagnostic ratio |
|---|---:|---|---:|---:|---:|
| first | 5 | uniform | 0.008493 | 0.005403 | 1.572 |
| first | 5 | known_decomposable | 0.006553 | 0.004838 | 1.354 |
| first | 9 | uniform | 0.018291 | 0.030613 | 0.597 |
| first | 9 | known_decomposable | 0.074534 | 0.044991 | 1.657 |
| first | 11 | uniform | 1.076956 | 2.965193 | 0.363 |
| first | 11 | known_decomposable | 0.365710 | 0.714990 | 0.511 |
| enumerate | 5 | uniform | 0.011132 | 0.003312 | 3.361 |
| enumerate | 5 | known_decomposable | 0.010112 | 0.003879 | 2.607 |
| enumerate | 9 | uniform | 0.047644 | 0.044641 | 1.067 |
| enumerate | 9 | known_decomposable | 0.111963 | 0.082941 | 1.350 |
| enumerate | 11 | uniform | 1.679022 | 4.339943 | 0.387 |
| enumerate | 11 | known_decomposable | 0.767995 | 1.794750 | 0.428 |

The direct function solver is slower than the best measured Semaev control
at five bits; it also loses at nine bits in complete enumeration and in the
supported first-relation stratum. Its eleven-bit diagnostic ratios improve,
but that cannot discharge the three-size criterion. No full attack cost,
new counting-floor ratio, asymptotic exponent or rho speedup is established.

## Evidence and proof trail

* solver_04: 144 trials and 6,944 coefficient/branch vector checks. Conditioned
  quadratic membership improves nine-bit first solving but still fails eleven
  bits. Its single exhaustive circuit-validation target has no admitted
  functions; positive coverage comes from checked SAT campaign models.
* solver_05: 48 trials plus complete oracle equality on all 43 affine five-bit
  targets. Two Artin-Schreier solves replace branch SAT. The full derivation,
  completeness proof and finite-search bound are in its README. The slowest
  benchmark cell took 1.092205 seconds, including field-operation counting.
* solver_06: 96 trials. Exact polynomial identities validate the elementary
  S4 expression and the denominator-cleared transformed S4 adapter. Both
  preserve the x-subspace base and exact target; no two-torsion-shift shortcut.

Source hashes, all 288 trial aggregates, exact cross-variant target lists,
and exclusion of the old solver_02/03 targets were verified. Zero validation
errors occurred. Full raw evidence is retained, including timeouts. Field
API operation vectors are available for the direct solver; comparable SAT
operation counts and field-to-Boolean calibration are not. Instrumentation
and Python overhead differ, so even the timing ratios remain preliminary.

## Next decision gates

Cache L_V once per instance and move fixed-b norm work out of the inner
root loop; charge that setup explicitly. Freeze the optimized solver before
selecting a new held-out panel. Reuse the same amortization policy for the
controls. Extend across larger supported odd field degrees and several
subspace dimensions; repeat seeds and budgets to measure censored outcomes.
Instrument SAT solving before asserting an operation-count improvement.

The pure coefficient SAT route remains an open alternative. The successful
candidate is a hybrid: it conditions on a root abscissa, generates function
coefficients algebraically, enforces full H-divides-L_V support, and only then
extracts and verifies curve points. This distinction is part of the result.
