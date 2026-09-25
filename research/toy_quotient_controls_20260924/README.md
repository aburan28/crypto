# Finite symmetry quotient and representation controls

Follow-on to [PR #701](https://github.com/aburan28/crypto/pull/701), restricted to
its GF(256) curves, four x-coordinates and at most six Boolean selector variables.
This experiment measures complete **toy decomposition calls** from supplied
supports and targets. It accepts no external cryptographic inputs and performs
no discrete-log recovery.

```sh
python3 -m pip install -r research/toy_quotient_controls_20260924/requirements.txt
python3 -m unittest discover -s research/toy_quotient_controls_20260924 -v
python3 research/toy_quotient_controls_20260924/verify_replay.py
python3 research/toy_quotient_controls_20260924/run.py --output /tmp/quotient-replay
python3 research/toy_quotient_controls_20260924/verify_replay.py /tmp/quotient-replay
```

Use CPython 3.12 for opcode comparisons. SymPy 1.14.0 and mpmath 1.3.0 are pinned.
A run refuses to overwrite existing evidence and journals each completed case
before producing gzip-compressed JSON. The CI job performs an entire replay.
Uninstrumented phase times are retained as secondary observations; they are
excluded from deterministic replay comparisons. No timing speedup is claimed.

## Protocol and hypotheses

`protocol-v1.json` preserves the initial protocol, committed before implementation.
Correctness preflights showed that some six-variable SymPy F5B calls exceeded
15 seconds. Before the production run, version 2 added a deterministic 100,000
Python-call budget for F5B, preserving its incompletions instead of allowing an
unbounded audit. Both initial versions are retained. No production cases or seeds were
selected using measured outcomes.

Run-001 is retained as superseded because its profiling-hook exception could
fire during generator cleanup. Version 3 fixes only that audit instrumentation:
the same call budget is checked at line boundaries in the main F5B routine;
one internal step can overshoot the nominal cap. Garbage collection is disabled
only during audit counting and restored afterward. Run-002 is the canonical
corrected replay. Original independent-audit source and contract snapshots are
stored alongside run-001. No numerical success threshold or corpus input changed.

The production contract fixes four support seeds: 101 and 503 from the previous
experiment, and fresh seeds 809 and 1601. It selects eight declared subgroup
target indices, both summand counts, all five curve models and six variants:
1,920 systems, three uninstrumented repetitions and one instrumented replay each.
The full earlier 1,280-system corpus is also replayed as a dependency check.
The new eight-target subset is explicitly a different frozen suite, not a claim
to have added these controls to every target in the preceding corpus.

The hypothesis is that a finite invariant-coordinate quotient can lower the
complete-call Python opcode proxy after setup and lifting. The frozen threshold
is at least 20% lower aggregate proxy cost on **each** fresh seed, model and
summand count, with identical verified signed-point multisets and no correctness
failure. All source-neighbor deltas and regressions are retained. We do not
select a winning neighbor or representation after seeing the holdouts.

The stage reference is the ordered variant on the same supplied workload.
The equivalent-suite exception is appropriate here: the WDSat full-DLP benchmark
cannot run this finite-domain reference experiment. Full-DLP cost, S, rho/floor
ratios and speedup remain null. A Python proxy reduction is neither a native
runtime speedup nor evidence of an asymptotic or cryptographic advance.

## The six variants

| Variant | Change | What is held fixed |
|---|---|---|
| ordered | Original ordered selector encoding | Reference inputs and output convention |
| canonical | Require nondecreasing selector labels | Same Boolean variables, all signed lifts |
| relabel | Fixed permutation [2,0,3,1] of the four support labels | Actual support, target and decompositions |
| scale | x'=u²x, y'=u³y, u=2; evaluate equations with a1'=u, a6'=u⁶b | Isomorphic curve and transported x-values |
| reverse_equations | Reverse the field-bit equation list | Exactly the same ideal and degree filtration |
| invariant | Elementary symmetric coefficient tuple; finite image addressed by rank | Same unordered signed-point decompositions after root recovery |

We test one fixed permutation of the four selector values, not all 24 permutations. Every permutation on four values is affine over
F2². Consequently this relabeling preserves polynomial degree and the complete
Macaulay filtration, though elimination work can change. Equation reordering
also preserves the filtration; it can change the F5 criterion's prefix cost.
These are correctness controls as well as cost measurements.

Rescaling is an actual isomorphism to a generalized binary model. The S3 value
scales by u⁸ and the S4 value by u²⁴; infinity targets have their corresponding
lower-polynomial weights. Multiplication by a nonzero field constant mixes the
bit equations invertibly. The zero set is unchanged, but the presentation by
individually degree-bounded generators can change. Tests verify the transformed
curve equations and summation-polynomial identities, not merely root counts.
All five unscaled models have a1=1, so point verification uses their shared
normalized addition formulas on the supplied original points.

## Genuine finite quotient, with its enumeration cost exposed

For an unordered x-multiset, form

    H(T) = product_i (T + x_i) = T^m + e1*T^(m-1) + ... + em.

The coefficients e1,...,em are elementary symmetric functions. This tuple
separates multisets, including multiplicities, in characteristic two. We retain
only the coefficients, not a table of original root tuples. There are 10 possible
coefficient tuples for two summands and 20 for three, given four support values.
They are sorted lexicographically and addressed by four or five Boolean bits;
a guard excludes unused bit patterns. Original selectors use four or six bits.

This is a finite invariant image with a rank encoding, **not** a scalable system
with freely varying field-valued invariant coefficients. The coefficients alone
would use 16 or 24 field-basis bits. Enumerating the permitted image is deliberate
preprocessing and is fully charged on every decomposition call; no target-to-
target cache or precomputed inverse table is used. Rank labeling can itself
change Boolean equation presentation, so no representation-independent degree
claim follows from this experiment.

The Semaev equation is evaluated directly in e1,...,em before roots are recovered.
For three summands, `symmetric_value` implements the exact characteristic-two
symmetric rewrite of the quadratic S3/S3 resultant. Tests compare that identity
on arbitrary field values, including repeated roots and infinity targets.
For each surviving coefficient tuple, `recover` divides H repeatedly by T+x for
each support x. Repeated division recovers multiplicities; a non-splitting tuple
is rejected. We then enumerate signed lifts and verify point sums.

Every variant outputs the same convention: **all unordered signed-point tuples**
whose sum is the target. This avoids comparing ordered counts against quotient
counts or rewarding an implementation for returning less information. An
independent exhaustive point oracle checks exact set equality, not only counts.

## Accounting and degree labels

Each decomposition call has five exclusive measured phases:

1. Construct selector metadata or the finite invariant image.
2. Evaluate field equations and interpolate their Boolean ANF.
3. Run the matrix-F5 criterion, certify ideal completion, and extract all roots.
4. Decode selectors or recover polynomial roots with multiplicities.
5. Lift signs, verify point sums, and deduplicate the output multisets.

An additional traced replay counts executed Python opcode dispatches in each
phase and all of its Python callees. Their sum is the common complete-call
**proxy unit**. CPython 3.12 requires an opcode-enabled frame before `settrace`;
a regression test catches silently missing the first phase. Three uninstrumented
repetitions must agree on every deterministic result and matrix counter.

C-internal instruction costs, hardware calibration, interpreter startup,
arithmetic-table generation, fixture construction, support/target generation,
map certification and independent audit work are outside this call metric.
The supplied arithmetic tables and curve inputs are shared fixtures. Table and
fixture setup are not silently amortized into a full-DLP claim: that total stays
null. The proxy counts some common measurement bookkeeping in the matrix path
and is specific to this Python reference, not native F4/F5 performance.

The Boolean solver imports the previously verified bounded matrix implementation.
We independently repeat its F4/F5 row-space comparison outside the measured path.
The reported completion degree is the first Macaulay budget whose span contains
all generators and is closed under multiplication by all variables. This is not
intrinsic degree of regularity. F4 and the F5 row criterion have identical degree
on the same input; pruning affects work, not which row space is obtained.

## Independent algebra audit

A predeclared 240-case subset covers seeds 101 and 809, targets 0 and 7, all five
models, both summand counts and every variant. We build an ordinary polynomial
ring over GF(2), add x_i²+x_i explicitly, and use grevlex order. Boolean field
equations precede the recorded input generators. This construction is independent
of the Boolean bit-matrix elimination code.

SymPy's Buchberger computation must complete, and its reduced basis must pass
input-membership, all S-polynomial reductions, and exhaustive Boolean root checks.
SymPy's **actual signature-based F5B** is also run with its declared call cap.
Whenever it completes, its reduced basis must equal Buchberger's exactly.
Budget-exceeded F5B runs remain in the evidence; they are not passes and their
partial degree observations are not complete solving-degree measurements.

The audit records polynomial degrees at entry/exit to F5B signature reductions,
plus final reduced-basis degree. It does not observe all intermediate work or
preprocessing; these fields are explicitly not intrinsic regularity or a global
maximum F5 solving degree. The observer restores SymPy's function and Python's
profiling state even after a cap, timeout or verification error. Negative tests
check the cap and reject deliberately wrong roots.

References: [SymPy polynomial algorithms](https://docs.sympy.org/1.14.0/modules/polys/internals.html#groebner-basis-algorithms),
[SymPy F5B source](https://github.com/sympy/sympy/blob/1.14/sympy/polys/groebnertools.py),
and the [preceding experiment](../toy_f5_neighbors_20260924/README.md).
