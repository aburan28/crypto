# Boolean F5 API call ceiling with output validation outside the timer

## Why this is a new experiment

The merged [validation-inclusive screen](https://github.com/aburan28/crypto/pull/1218)
reported an optimistic n=24 ceiling near 1.61x for a metric that spent about
half its outer time SHA-256 hashing returned polynomials. That registered
negative result remains correct for its own metric, but it cannot bound the
F5 API call alone. Its raw samples, source and unused holdouts are not
relabelled. This protocol fixes the cost boundary **before** collecting new
data on disjoint discovery and holdout seeds.

For each generated public Boolean system, let `C` be the time from
immediately before `matrix_f5_f2_with_form_timed` until immediately after
it returns, with its complete `Vec<F2BoolPoly>` materialized. Then compute
and independently compare an exact returned-row digest **outside** the
timed interval. Let `D` be a separate timer around destruction of those
returned rows and the report, after validation. Define charged API cost
`T = C + D`; retain `C`, `D` and untimed validation work separately. The
library's exclusive `criterion_ns`, `build_ns`, `reduce_ns` and `unpack_ns`
must sum to at most `C` on every call. The deliberately optimistic ceiling
for a cache that changes **only** matrix build and reduction is

`U_call = T / (T - build_ns - reduce_ns)`.

This is a conditional Amdahl bound under the unchanged polynomial-output
API, not a measured candidate speedup. It pretends both affected phases
become free; a real cache adds setup, fallback and memory costs. It does not
bound a method that also changes criterion evaluation, output representation
or relation collection. Exact digest work is outside `T` but inside the
whole-process receipt, so it is visible and cannot be mistaken for a free
full-pipeline audit.

## Frozen route, fixtures and correctness

Use the same source route and public synthetic fixture contract as
`research/boolean_f5_amortized_ceiling_20261002/PROTOCOL.md`:
`F5OutputForm::Echelon`, degree four, n=12/16/20/24, m=n, 2n distinct
quadratic monomials per generator, fixed quadratic core, independent-affine
and walk-affine tails, full multiplier mask, batches 2/8/32, and one-thread
execution. Launch a fresh cell process with
`KIC_F5_DIRECT_PACK=1`, `KIC_F5_UNPACK_DIRECT=1`,
`KIC_GF2_TABLES=4`, `KIC_F5_AVX512_UNPACK=0`,
`KIC_GF2_REUSE_TABLE=0`, and `RAYON_NUM_THREADS=1`; all other F5/GF2
option variables are unset. Record the inherited exact full-column
direct-pack hit or sorted-row fallback on every call; require direct
unpacking. No route, input or unfavorable cell may be selected after timing.

For every timed call, check F5 report counters and exact returned-row digest
against an untimed reference for the same immutable input. Validate the
small n=12/16 row space against the inherited F4 routine once per cell
outside `C` and `D`. A changed result, false route flag, phase overrun or
failure to materialize the returned polynomials censors the cell. Fingerprints
are diagnostics for this bounded experiment, not an unaffiliated proof of
the F5 criterion. No curve, key or external target is accepted.

## Fixed grid and gate

`protocol.json` fixes discovery seeds 20261011/3141667 and untouched
holdout seeds 20261018/4242509. Each phase includes the four n values,
two families, batches 2/8/32 and seven balanced A/A repetitions, with one
fresh process per cell. The four primary groups are n=24, batch 32,
each seed/family. Bootstrap 4,000 resamples of complete cold batch ceilings
with a fixed seed; retain each median, 95% **upper** bound and A/A noise
floor. Also report the development-size guardrails, source/binary/protocol
hashes, direct-pack fallback count, peak RSS, all exclusive phases, digest
validation work and total process CPU/wall.

Reserve one Linux x86-64 AVX2 physical core by reading logical CPU 2's full
SMT sibling list from sysfs and passing all siblings to the repository's
isolation controller. Keep one worker thread, record all remaining threads
and require an uncontended receipt: other-process CPU at most 10% and CPU
PSI some avg10 at most 5.0. The worker cap is 900 seconds and raw evidence
cap is 64 MiB. A build failure, resource refusal, timeout, OOM, incomplete
cell or verifier error is CENSORED with null ceiling, never a mathematical
rejection. Seal and retain each attempt.

Advance to fresh holdouts only if **all four** primary discovery groups have
a 95% bootstrap upper bound of median `U_call` strictly above 2.0, with
correct outputs and qualified resources. If one is at most 2.0, retain the
negative and stop before holdouts. Passing means only that a later actual
graded F5 batch cache might have room to reach 2x; that cache requires a
separate same-binary complete-call comparison, then natural relation-yield,
independent rank and automorphism-aware rho accounting before any
cryptanalytic claim.

Implement the producer, verifier and analysis in native Rust. Thin shell
may invoke the repository's existing `tools/isolated_bench.py` controller.
Do not introduce Python research algorithms or alter the prior frozen
evidence. Full IC cost and rho ratio remain null in this screen.
