# Projected full-rank certificate for sparse F5 output

## Hypothesis, output contract, and proof boundary

On the frozen n24 degree-4 workload, the F5 criterion produces 6,924 rows
and the current eliminator reports rank 6,924. Its returned echelon basis
has about 13.7 million terms, and output materialization is about 46 ms of
an 86 ms one-thread F5 call on the latest qualified x86 host. The original
F5-surviving Macaulay rows are much sparser. This experiment tests an
explicit **independent-row basis** output form, separate from `Reduced`,
`Echelon`, and `SelectiveEchelon`; it does not promise the same ordered rows
or leading terms. Callers requiring an echelon/Gröbner basis must continue
to use the existing forms. No default API changes in this experiment.

Let `M` be the built `r × c` binary matrix. Fix a deterministic linear map
`P: F₂^c → F₂^k` that maps each original column to one projected column by
XOR. If the projected matrix `MP` has rank `r`, then `rank(M) ≥ r`, hence
`rank(M)=r`. The original `r` rows are therefore independent and form a
basis of exactly the same row space as the echelon output. This implication
is exact for any projection; hash collisions only cause a false negative.
If projected rank is below `r`, use the established full elimination and
return its echelon rows. The result remains exact for every input. Record
whether the certificate succeeded and never infer full rank from a prior
benchmark, a fingerprint, or a probabilistic assumption.

The objective is a further **2×** on the complete one-thread n24 degree-4
call for this explicitly different output form, including projection,
certificate elimination, fallback if needed, and materialization of returned
polynomials. A ratio for different output forms is labelled a solver-stage
alternative-contract comparison; it is not an exact-row speedup or an
end-to-end IC/DLP result. Promotion requires a separate decision about
caller semantics and downstream work even if the timing threshold passes.

## Frozen reference and candidate

The baseline is merged main
`d4d46d59bd1190d105af7678163ca5d74f04b505`.
SHA-256: `matrix_f5_f2.rs`
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`,
`gf2_elim.rs`
`6fd3f10ef8e681490784c624e703eeb8fa887c5817ec9bdf34844e3f943c7f5c`,
benchmark `f4_f2_bench.rs`
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`,
isolation wrapper
`ff59e53fd5566078ad493c8be0d06168699c4fd194bc881d99df19fa20e007ef`.
Build and retain the unmodified release binary and host record before code.
No timing from a different CPU is used as the paired reference.

Add `F5OutputForm::IndependentRows`, selected in the benchmark only by
`KIC_F5_ECHELON=3`. It may attempt the certificate only for degree 4,
`n_vars ≥ 24`, at least 4,096 rows and at least 8,192 columns. All other
cases follow `SelectiveEchelon` exactly. Build `M` with the existing F5
criterion and direct packed builder. For `k=min(c, r+s)`, with frozen slack
`s=512` or `1024`, map column `j` to
`splitmix64(0x9e3779b97f4a7c15 XOR j) mod k` and XOR every set entry of
each original row into its mapped projected column. This map is independent
of row values, seeds and timings. Here `splitmix64(x)` adds
`0x9e3779b97f4a7c15` modulo 2⁶⁴, then applies XOR-shift 30 and multiply
by `0xbf58476d1ce4e5b9`, XOR-shift 27 and multiply by
`0x94d049bb133111eb`, and final XOR-shift 31, all modulo 2⁶⁴.
Use the existing exact GF(2) echelon
kernel on the projected matrix, charging projection and its elimination
to the reduction phase. If rank is `r`, decode and return the retained
original packed rows; otherwise run the full reference eliminator and
return its echelon rows. Cap peak packed-row memory at 2× the original
matrix plus small maps. Keep a portable scalar CPU path and unchanged
default behavior. Emit an actual route/certificate indicator.

## Frozen workloads and progression

Use one Rayon thread, the seven F5 cases in `examples/f4_f2_bench.rs`,
seed XORs `0` and `badc0de1`, and primary `f5_n24_m24_d4`. Fix selective
echelon on the reference, fused row counting, direct packed build, scalar
direct unpack, four Gray-code tables with table reuse, AVX2 row XOR and
branchless strip where available; disable deferred-above, word batching,
AVX-512 unpack and experimental table building. Use separate processes
for reference and each candidate slack so env caches do not cross arms.

Commit this protocol and open its draft PR before candidate code or timing.
First, release tests must prove the projected-full-rank implication on
small full-rank, rank-deficient, sparse and dense matrices, and show fallback
equals the reference output exactly. The F5 and GF(2) focused release tests
must pass. Then run an untimed-for-decision structural screen: one reference
and one candidate process per slack and seed, preserving full stdout,
stderr, source/binary hashes, host and failures. Correctness requires equal
rank, canonical row-space fingerprint, criterion work and built/pruned row
and column counts in all seven cases. Six inactive cases must retain exact
raw-row fingerprints, term counts and counted reduction work. The primary
must have a successful full-rank projection and at most 25% of reference
output terms on both seeds. Select slack 512 if it passes; otherwise 1024
if it passes; otherwise reject without paired timing.

If selected, run one Apple ARM64 local screen: warmup each arm, five A/A
and five alternating A/B pairs per seed, preserving every process. Use
`tools/isolated_bench.py reserve` when available; this macOS host currently
lacks its Linux CPU reservation path, so its ratios are exploratory. Reject
locally only if both primary complete-call medians are below 1.10 or any
exactness check fails; otherwise advance to qualified Linux x86-64.

On Linux x86-64 AVX2/BMI2, run same-binary one- and two-thread five-pair
A/A and A/B blocks on seeds `0`, `badc0de1`, `5eed2026`, `f5c02a28`.
Select the first clean block for each seed without reading timings; retain
all preflight failures, timeouts, OOMs and contended blocks. Qualification
requires zero eligible user threads and zero contended samples. The 2×
target requires the one-thread paired complete-call median and its exact
3,125-resample 95% lower bound above 2.0 on all four seeds, no inactive
case below its A/A minimum, and no two-thread regression. A promising
alternative-contract result below 2× is recorded but does not meet the
user's requested threshold. No output-form promotion follows from a
benchmark alone.

These are Boolean matrix-F5 solver-stage measurements, not one-target IC
online wall time, DLP recovery, or paired Pollard-rho speedup.
