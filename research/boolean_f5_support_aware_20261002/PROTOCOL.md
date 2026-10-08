# Preregistered support-aware Boolean F5 batch continuation

## Boundary and hypothesis

The [exact graded batch discovery](https://github.com/aburan28/crypto/pull/1236)
finished on a reserved Linux x86-64 core with complete-batch paired medians
of 1.784–1.838×, below its universal >2× gate. It retained 48 cells and
18,816 timed calls at [run 37045285547](https://github.com/aburan28/crypto/actions/runs/37045285547);
0/4 n=24 primary groups passed, 8/8 n=16/20 nonregression groups passed,
and its holdout stayed unused. In each primary independent-affine batch one
assignment took the fresh-F5 fallback; in each walk-affine batch two did.
Cold compilation cost about 70–72 ms per n=24 batch. A roughly 5.6×
reduction in counted GF(2) word XORs yielded only about 1.8× wall speed
because materialisation, setup and fallbacks remained charged.

The *stage reference* here is the inherited, same-binary complete Boolean
matrix-F5 `Echelon` API at degree four on the same generated public systems.
The *anchor* is the exact graded-prefix cache of that rejected discovery,
rebuilt as an unchanged arm in the same native binary. The mathematical
output-size floor remains the number of terms in the required ordered
`Vec<F2BoolPoly>`: every arm must emit them and pay allocation/destruction.
The generic-group and automorphism-discounted Pollard-rho boundary for a
complete ECDLP remains separate and unmeasured. This experiment can cross
only a complete-F5-batch engineering gate; it cannot establish a new
index-calculus exponent, natural Semaev relation yield or a rho crossover.

The primary hypothesis is that the high prefix remains fixed even when the
source F5 builder's *lower-degree column support* changes. Let `M_i` be the
selected packed matrix for affine assignment `i` and `k` the largest
multiple of 64 at or below `binom(n,4)`. Its ordered selected-row labels
and first `k` columns `H = M_i[:,0..k]` must be exactly the same as the
compiled base. Let `U_i` be the sorted source-column list actually occupied
outside that prefix, and `P_i` delete only columns absent from `U_i`. If
every prefix column is present, the row transform `E` computed from `H`
commutes with this column projection:

`E · P_i(M_i) = P_i(E · M_i)`.

The deterministic M4RI pivot trace through the first `k` columns is the
same for every such assignment. The candidate can therefore apply `E` to
the changed suffix, project to the exact source fallback column layout,
then resume elimination at word `k/64` and the retained prefix rank. A
missing high-prefix column, changed criterion row-label bitstream,
quadratic core, generator order, multiplier schedule, or layout takes a
charged fresh-F5 fallback. The projection theorem is conditional on these
guards; it is not a claim that arbitrary sparse-support systems share a
pivot schedule.

## Frozen arms and attribution

One native Rust binary contains these arms on identical input order:

1. `fresh`: inherited complete F5 reference.
2. `anchor`: the exact graded cache semantics frozen by PR #1236.
3. `support`: `anchor` plus the guarded lower-column projection and exact
   word-aligned M4RI continuation; no unpack change.
4. `support_unpack`: `support` plus a direct preallocated unpack equivalent
   to the inherited direct scalar unpacker; it may not change the ordered
   output or omit destruction.
5. `support_unpack_layout`: `support_unpack` plus cached static
   quadratic-product row layout. It still recomputes and charges the native
   Boolean F5 criterion on *every* assignment, rebuilds each changing
   affine tail, checks the complete selected-row labels and source support,
   and takes a charged fallback on any mismatch.

The last arm is the registered primary candidate; the preceding arms are
fixed ablations and may not be selected after seeing results. For each
repetition, an A/A fresh/fresh pair measures the noise floor. The even
ordering is `aa_a, aa_b, fresh, anchor, support, support_unpack,
support_unpack_layout`; odd repetitions reverse it. Every arm rotates the
same batch input order by repetition. One child process per cell starts
with no retained context; every timed arm constructs and destroys its own
context. The complete charged batch time is compilation plus every call,
criterion, packed construction, projection, low elimination, unpack,
fallback, returned-output destruction and context destruction. Direct
ordered-vector preflight is outside arm timers; every timed output is
checked against the same ordered-output SHA-256 digest outside arm timers.
The verifier independently rebuilds every unique fixture, compares actual
ordered vectors, output terms and criterion counters, and checks all
resource and source bindings.

Before registered timing, development-only seeds 17–20 may be used for
implementation, exact-output tests and profiling. Candidate source,
thresholds, input generator, arm ordering and analysis must be committed
before discovery seed execution. The parent discovery/holdout seeds may
not be used to tune the new candidate. No changes to the registered source
are allowed after inspecting this study's holdout timing; a correction
requires additive evidence and a new protocol with new seeds.

## Workloads, gates and stopping rule

`protocol.json` fixes n=12/16/20/24, m=n, 2n distinct quadratic terms per
generator, two changing-affine families, batches 2/8/32, discovery seeds
20261103/3141753, untouched holdouts 20261110/4242601, seven balanced
repetitions and 4,000 fixed-seed paired bootstrap resamples. Fixture
generation and the F5 route match the previous graded study. The n=12
cells are correctness/noise guards. The four primary groups are n=24,
batch 32, both families and both discovery seeds. Eight n=16/20 batch-32
groups are the nonregression controls. Report every smaller-batch and n=12
cell, including failures and fallbacks.

First require exact output, complete selected-row signature agreement,
source-column-layout replay, and a verified cache recovery for every
primary sparse assignment whose first `k` columns are fully occupied.
If that structural hypothesis fails, the corresponding arm falls back and
the structural gate rejects; it must not be relabelled as a timing win.
For a **complete-F5-batch engineering success**, all four primary groups
must have `support_unpack_layout`/`fresh` paired median and 95% bootstrap
lower bound above 2.0 and above their own A/A 97.5th-percentile floor.
All eight n=16/20 nonregression groups require a 95% lower bound above
0.95. The `support`/`anchor` ablation's 95% lower bound must exceed its
own A/A floor in each primary group, with the exact support-recovery counts
reported. The `support_unpack` and `support_unpack_layout` ablations are both reported
without selecting the better arm after timing. Discovery must pass every
gate before the untouched holdouts are run under identical source,
thresholds and analysis. A discovery-only pass is not confirmation.

Use Linux x86-64 AVX2, `RAYON_NUM_THREADS=1`, and reserve the entire
physical core identified by logical CPU 2's SMT-sibling list. The
repository's existing isolation controller must record uncontended
conditions, at most 10% other-process CPU and PSI some avg10 at most 5.0.
Cap retained context at 128 MiB, worker wall time at 2,100 seconds and
raw evidence at 64 MiB. Cap, timeout, OOM, failed cell, incomplete
receipt, changed source, or verifier disagreement is **CENSORED** with
null performance, never a positive result. A complete run below any
registered gate is **REJECTED**, and its fresh holdouts remain unused.
Retain full source/binary/protocol hashes, matched A/A and all arm samples,
fallback reasons, prefix/support checks, GF(2) word work, output digests,
peak RSS, failed attempts and native replay in sealed, additive bundles.

No curve, point, scalar or key input is admitted. Use Rust for all new
algorithm, producer, verifier and analysis code; thin shell may call the
repository's existing CPU-isolation controller. Full IC cost and rho ratio
remain null regardless of this stage's outcome.
