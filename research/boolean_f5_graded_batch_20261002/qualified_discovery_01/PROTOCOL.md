# Exact graded high-block reuse inside complete Boolean F5 batches

## Hypothesis and source boundary

The [fixed-core selected-row study](https://github.com/aburan28/crypto/pull/1214)
found identical Boolean-safe F5 prune signatures across changing affine
tails on fresh generated fixtures. The [call-only ceiling study](https://github.com/aburan28/crypto/pull/1221)
then found an optimistic n=24 build/reduction-only ceiling above 4x on fresh
fixtures under the inherited polynomial-output API. Neither measured an
actual cached F5 implementation. This protocol tests one: cache and replay
the elimination of the high-degree quartic block shared by a fixed quadratic
core, then construct and reduce each changed lower-degree block. The F5
criterion is still recomputed and charged **for every assignment**.

For quadratic `f_j=q_j+a_j` and degree-at-most-two multiplier `t`, the
degree-four projection of `t*f_j` depends only on `q_j`. The cache compiles
the exact selected generator–multiplier row labels and quartic packed
columns for the first input of each cold batch, plus a verified ordered
high-block row-operation schedule. A later input may reuse that schedule
only after the native criterion yields the same complete selected-row
bitstream, n, degree, generator order, quadratic support, multiplier mask,
row labels and column layout. Replay must verify each scheduled pivot and
row operation; any missing pivot, changed signature, support escape, cap or
semantic mismatch takes the charged fresh-F5 fallback. A hash may shortlist
a match, but full signature and context equality are required.

## Exact returned-output contract

The reference is the repository's current
`matrix_f5_f2_with_form_timed(..., F5OutputForm::Echelon)` with
`KIC_F5_DIRECT_PACK=1`, `KIC_F5_UNPACK_DIRECT=1`,
`KIC_GF2_TABLES=4`, `KIC_F5_AVX512_UNPACK=0`,
`KIC_GF2_REUSE_TABLE=0`, `RAYON_NUM_THREADS=1`, and all other F5/GF2
option variables unset. The source-defined direct-pack full-column guard
may take its sorted-row fallback; record and retain those cases. Compile
one native binary containing both arms, launch fresh cold processes per
cell and rotate A/B order. Do not divide in historical absolute times.

Both arms must return **byte-identical ordered `Vec<F2BoolPoly>` output**,
not only equal rank or canonical row space, on every timed input. The
reference's row-space fingerprint and small F4 cross-check are additional
controls. Record criterion prunes and word work, selected/built row counts,
rank, output terms/digest, high/low GF(2) word operations, schedule hits,
fallback reasons, retained context bytes and peak process RSS. The cache
may reduce measured reduction operations, but it may not change the F5
criterion, silently omit a syzygy, or relabel UNKNOWN/caps as UNSAT.

The primary timing unit is one **complete cold batch**. Candidate time
includes fixed-core cache compilation, all criterion recomputations,
changing-low-block construction, schedule application, low elimination,
every full-build fallback, returned polynomial unpacking, output allocation
and context/output destruction. Baseline time includes fresh complete F5
calls and output destruction for the same inputs. Independent exact output
verification is outside both arm clocks but inside process receipts, as in
the confirmed call-only ceiling study. Neither arm retains state between
separately timed batches. Correctness and the same output form take
precedence over an apparent timing win.

## Frozen workloads, resources and gates

`protocol.json` fixes n=12/16/20/24, m=n, 2n distinct quadratic monomials
per generator, two public changing-affine families (`independent_affine`
and `walk_affine`), batches 2/8/32, and seven balanced repetitions.
Discovery seeds are 20261020/3141691; untouched holdout seeds are
20261027/4242547. The generator recurrence, affine coefficient order and
input validation are identical to
`research/boolean_f5_call_only_ceiling_20261002/PROTOCOL.md`.
The n=12 cells are correctness/noise guards. Primary performance groups are
n=24, batch 32, both families and both seeds: **four groups**.

Under a single native Linux x86-64 AVX2 worker, reserve the whole physical
core named by logical CPU 2's SMT-sibling list. Keep one worker thread;
require the repository isolation controller's uncontended receipt with at
most 10% other-process CPU and PSI some avg10 at most 5.0. Record an A/A
baseline/baseline pair in every fixed cell, all refused preparations and
the exact source/binary/protocol hashes. Cap each arm at 128 MiB retained
context, the campaign at 1,200 worker seconds and raw evidence at 64 MiB;
cap, timeout, OOM, failed cell or verifier disagreement is CENSORED with
null performance, never a positive result. Seal complete and failed runs
without overwriting earlier artifacts. The native Rust verifier must
rebuild every fixture, compare exact output and counters, enforce source
and resource binding, and compute 4,000 fixed-seed paired bootstrap samples.

For a dramatic **complete-F5-batch engineering gate**, every primary group
must have a paired reference/candidate median and 95% lower bound above
**2.0**, above its own A/A 97.5th-percentile floor, with zero correctness
or resource failures. Also require a 95% lower bound above **0.95** on the
eight n=16/20, batch-32 seed/family non-regression groups. Report all n=12
and smaller-batch guardrails without selecting favourable cases. Discovery
must satisfy every gate before the untouched holdouts are run. The unchanged
source and thresholds must then satisfy the same gates on holdouts; a
discovery-only success is not confirmation. No source tuning after holdout
timing, seed exclusion or reference-arm change is allowed.

This is a Boolean matrix-F5 API/batch study on generated public systems.
Even a verified 2x complete batch is an engineering result for related
systems, not a changed exponent, natural Semaev decomposition probability
or Pollard-rho crossover. A later full index-calculus comparison would
still need independent relation yield and rank, factor-base setup, memory,
all online targets and automorphism-discounted rho cost. Those fields stay
null here. No curve, scalar or key input is admitted. Use Rust for all
new algorithm, producer, verifier and analysis code; thin shell may call
the repository's existing CPU-isolation controller.
