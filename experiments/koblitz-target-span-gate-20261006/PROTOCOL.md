# Frozen target-span gate for compact-orbit rank preparation

Status: preregistered before reading the target-span results. This is a
counterfactual **linear-algebra diagnostic** on the frozen native transcripts
from the [base-size holdout](../koblitz-base-size-cold-panel-20261004/RESULT.md),
not a new IC producer, timing run, `IC1` candidate, or speedup claim. The
source Koblitz curves are n=41 and n=53; no isogeny is involved. The fixed
comparison is the same actual base, point, relation order, and rank seed as
each archived run. Strong same-Q rho remains the reference, but no new rho
measurement or operation conversion is made in this diagnostic.

## Hypothesis and frozen input

The hypothesis is that a fixed target-independent prefix of the guided rank
stream might already determine the held-out target log, allowing some
precomputation to stop before all `K` orbit-column logs are known. The null
is that target coefficients are outside every proper prefix row span. Analyze
**all 24** successful held-out IC/rho directories: `n=41, K=64/85` and
`n=53, K=160/220`, repeats 1 through 6. There are only **four distinct
point/base cells**, since repeats intentionally share one Q per n and K;
do not treat 24 runs as 24 independent target samples. The input files in
each directory are `base.jsonl`, `rank.jsonl`, `ic_targets.jsonl`,
`replay.jsonl`, and `status.jsonl`. Verify the parent
`HOLDOUT_SHA256SUMS` before analysis; do not modify those files or select
additional targets after seeing their rows.

## Exact analysis and decision rule

Independently build the target coefficient vector modulo the recorded prime
subgroup order by summing the four point labels referenced by its accepted
relation. Insert the recorded rank equations in their original order into a
modular row basis, retaining each equation's right-hand scalar. After each
prefix, reduce the target vector against that basis. The earliest prefix
whose coefficient remainder is zero is `first_target_span_rank`; its reduced
right-hand side must equal the producer's recovered scalar. Check the
recorded row lengths, modular values, rank transitions, base identity,
success statuses and passing independent full-point replay. At full rank,
the target must be in span with the archived verified scalar, or the analysis
fails. Preserve a row for every input directory and report the four unique
cells separately from repeats.

The prespecified admission threshold for a *new producer* is a proper
prefix of at most `floor(0.8 K)` in **all six** repeats of at least one
cell, with the same verified scalar and no input/replay failure. If no cell
passes, close the fixed-prefix shortcut on these four frozen cells and
prioritize a fixed-base probe-reducing index/query experiment instead. A
positive cell would only justify a fresh disjoint-Q coverage and complete
producer test. Rank rows saved, rank probes removed, and any phase-cost
subtraction are counterfactual stage diagnostics, never measured wall-time
gains. A target-conditioned stopping rule would make subsequent rank work
target-dependent and must charge it to the **online** interval; the archived
3–12 ms online timings cannot be reused for that altered method. The
host-isolation gate and a complete `IC1` manifest remain required for a
controlled one-target IC/rho speedup claim.
