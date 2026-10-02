# Complete Boolean F5 call: optimistic ceiling for graded batch reuse

## Question and accounting boundary

The [fixed-core F5 signature study](https://github.com/aburan28/crypto/pull/1214)
verified that the inherited Boolean-safe F5 criterion selects the same
generator–multiplier rows as affine tails change on its frozen generated
systems. It did not time a complete F5 call. Before implementing a batch
cache, measure whether reusing only the high-degree matrix construction and
elimination could **possibly** double a complete F5 call while the criterion,
returned polynomial representation and public API stay unchanged.

For each call, use the outer cold time `T` from immediately before
`matrix_f5_f2_with_form_timed` through destruction of its returned rows and
report. Exact returned-row digest validation occurs within that outer timer;
the independent small F4 cross-check is supplied-reference work outside it.
The source already reports exclusive `criterion_ns`, `build_ns`,
`reduce_ns` and `unpack_ns` inside the call. Define the deliberately
optimistic ceiling

`U = T / (T - build_ns - reduce_ns)`.

The denominator keeps criterion, unpacking, destruction and any uninstrumented
call overhead, while pretending that **all** matrix building and reduction
cost zero. A graded cache that changes only those two phases cannot beat U
under this fixed output contract. This is a conditional Amdahl bound on one
F5 API route, not a speedup measurement or a bound for methods that also
change criterion evaluation or output representation. Require all phase times
nonnegative and their sum at most T; refuse any sample that violates the
accounting identity. Do not subtract fixture generation or failed calls.

## Fixed route and public systems

Call `matrix_f5_f2_with_form_timed(..., degree=4,
F5OutputForm::Echelon)` on the repository's current source. Launch a fresh
process for each cell with `KIC_F5_DIRECT_PACK=1`,
`KIC_F5_UNPACK_DIRECT=1`, `KIC_GF2_TABLES=4`,
`KIC_F5_AVX512_UNPACK=0`, `KIC_GF2_REUSE_TABLE=0` and
`RAYON_NUM_THREADS=1`; leave all other F5/GF2 option variables unset.
The worker must confirm direct unpacking and record the actual direct-pack
flag in every returned `F5Timings`. The source's full-column guard can
fall back to its sorted packed builder when selected rows leave ambient
columns unused; such cases remain in the fixed grid and are reported,
not censored or selected away. This option set is the comparison's named
reference; do not silently substitute a historical option. The
candidate batch cache is not implemented in this study.

Use the same generated public fixed-quadratic-core fixture construction and
two changing-affine families as
`research/boolean_f5_affine_signature_20261002/PROTOCOL.md`: n=12/16/20/24,
m=n, 2n distinct quadratic monomials per generator, and full n-variable
multiplier mask. Batches are 2/8/32. Each timed F5 call starts with fresh
internal matrix and criterion state on its immutable supplied polynomial
input; input construction is common fixture work outside the call clock.
No cached F5 state leaks across repetitions. The
source-pinned Boolean F5 result must complete, report unchanged selected and
built row counts within a core, and have an exact stable row-space fingerprint
across repeated calls on the same input. An independent small-system F4
cross-check and existing repository F5 tests remain correctness controls;
none of their costs is hidden in the timed F5 call.

## Frozen measurement and stop rule

`protocol.json` fixes discovery seeds 20261008/3141627 and untouched
holdout seeds 20261015/4242463, n=12/16/20/24, two affine-tail families,
batches 2/8/32 and seven balanced repetitions. The n=12/16/20 cells are
reported guardrails. The primary cells are n=24, batch 32, both families
and both seeds: **four groups**. Each group has its own paired A/A timing
noise, source and binary hash, CPU feature record, peak RSS, phase times,
F5 report counters, row-space fingerprint and cap status.

Use a reserved Linux x86-64 AVX2 CPU through the existing isolation
controller, with one worker thread, every thread left on the reserved CPU
recorded, other-process CPU fraction at most 10%, and CPU PSI some avg10 at
most 5.0. Retain every refused preparation. A 900-second worker cap, 64 MiB raw
evidence cap and exact source/protocol hashes apply. A timeout, OOM, false
route flag, incorrect result, incomplete cell or contended receipt is
CENSORED and supplies no ceiling. The native Rust verifier independently
reconstructs fixtures, validates every phase/counter/route record, checks
hashes and resource admission, and computes deterministic 4,000-resample
bootstrap intervals. Its result is sealed before interpretation.

Advance to fresh holdouts only if every primary discovery group has a 95%
bootstrap **upper** bound for the median optimistic ceiling strictly above
2.0 and all correctness/resource checks pass. If any primary upper bound is
at most 2.0, this fixed-contract high-block-only route cannot support a
universal 2x complete-call claim on the measured grid; retain the failure
and stop without running holdouts. A passing ceiling merely permits a later
separate implementation and matched timing of the actual cache; it is not
evidence that any speedup occurred. Do not choose a subset, alter the output
form, change seeds, or reclassify a censored case after seeing phase times.

This screen has no curve/key input, complete index-calculus relation yield,
calibrated group operation cost or Pollard-rho comparison. Those fields stay
null. Use Rust for producer, verifier and analysis; thin shell may invoke
the repository's required CPU-isolation controller. No Python research
execution path is introduced.
