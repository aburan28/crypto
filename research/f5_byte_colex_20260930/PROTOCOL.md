# Byte-table colex indexing for Boolean matrix-F5 row packing

## Hypothesis and frozen reference

The direct full-column F5 builder computes a combinatorial column index for
every product term by repeatedly finding and removing its lowest set bit.
At n24 degree 4 that loop is executed millions of times. Replace the loop
with three small byte lookups whose entries sum the same binomial coefficients.
This is a strict implementation change: packed rows, returned polynomials,
all counters, fingerprints, and ranks must remain identical. The intended
incremental gain is in matrix construction; its upper bound on the complete
call is limited by that phase's measured share. The user's further 2×
complete-call goal is retained as a separate threshold, not inferred from
build-only timing.

The baseline is merged main `d40843d50e13e541943a89786e009165a702fc5d`.
Source SHA-256 values: `koblitz_groebner.rs`
`cc50fa7dee14d6c606883713fe842139a1d25f636c57be99817ca63607c2b816`,
`matrix_f5_f2.rs`
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`,
`f4_f2_bench.rs`
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`,
and isolation wrapper
`ff59e53fd5566078ad493c8be0d06168699c4fd194bc881d99df19fa20e007ef`.
PR #801 merged the historical frozen-source snapshot; this experiment does
not alter that snapshot or carry its old measurements forward.

## Candidate and frozen workload

`KIC_F5_COLEX_BYTES=0` retains the existing loop and is the reference;
`=1` enables the byte lookup only within the direct full-column packed
builder. Split each at-most-24-variable mask into low, middle and high bytes.
For each byte position and number of prior set bits `0..=4`, precompute
the sum of `C(bit_position,j)` for set bits in that byte. A 256-entry byte
popcount table supplies the prior counts and total degree. The lookup table
has at most `3 × 5 × 256` 16-bit sums plus 256 byte counts. Build it once per
F5 call when opted in and charge that setup to the build phase. Retain all
fallback paths. Record an actual route indicator in the benchmark JSON.

The seven benchmark cases in `examples/f4_f2_bench.rs` use seed XORs `0` and
`badc0de1`; the primary case is `f5_n24_m24_d4`. Fix one Rayon thread,
selective echelon, fused counting, direct packed build, scalar direct unpack,
four Gray-code tables with reuse, AVX2 row XOR and branchless strip where
available; disable deferred-above, word batching, AVX-512 unpack and other
experimental table construction. Use separate processes per arm to isolate
OnceLock environment choices. The reference and candidate are one binary.

## Gates, accounting and stop conditions

Commit this protocol and open its draft PR before candidate code or timing.
Build and test through `tools/isolated_bench.py busy`. Unit tests must cover
all 8-bit patterns at prior counts `0..4` against the scalar colex sum for
valid degree-at-most-four masks and compare direct packed matrices with
and without the lookup on small random and sparse systems. The release GF(2)
and F5 test suites must pass.

First run two seeds of the full seven-case suite in separate reference and
candidate processes, retaining complete JSON stdout, stderr, exit status,
host, binary and source hashes. Every case must match exactly in `rows_fp`,
`row_space_fp`, output terms, rank, all reported row/column and word-operation
counts, and criterion counters. The route indicator must be true only when
the direct builder runs. Any mismatch or failure rejects the candidate
before timing.

If exactness passes, run Apple ARM64 same-binary local pairs under the
repository's isolated benchmark reservation when available: one warmup per
arm, five A/A pairs and five alternating A/B pairs per seed, preserving each
complete-call and build phase time, full output, resource receipt and load.
Qualified local advancement requires both seeds' complete-call median at
least 1.03, no seven-case regression below the A/A minimum, and build phase
median at least 1.10. If local isolation is unavailable, treat its ratios as
exploratory and advance to Linux unless both seeds' complete-call medians
are below 0.95 or exactness fails. A failed gate archives the source patch
and every receipt, then removes the runtime option.

A passing local candidate requires Linux x86-64 AVX2/BMI2 one- and two-thread
paired A/A and A/B runs on the same two seeds plus `5eed2026` and `f5c02a28`,
with five pairs per mode and seed. Use the existing isolated CI benchmark
method: first-clean preflight, zero eligible user threads and zero contended
samples. Incremental promotion requires a one-thread full-call paired median
at least 1.05, a 3,125-resample 95% bootstrap lower bound above 1.02, no
seed or smaller case below its A/A minimum, and no two-thread regression.
The user's further 2× target requires the full-call median **and** lower
bound above 2.0 with those guards. Keep any accepted option opt-in until
the exact-output tests and cross-architecture checks justify a default.

These are Boolean matrix-F5 solver-stage diagnostics. They do not measure
one-target IC online wall time, a recovered DLP, or a paired rho solve.
