# AVX2 four-column compaction for matrix-F5 row unpacking

## Frozen hypothesis and reference

The accepted one-thread matrix-F5 fast path spends about 57 ms unpacking
13.7 million terms after n24 degree-4 elimination on an AMD EPYC 7763
(`research/gf2_avx2_table_build_20260929/RESULT.md`). The direct scalar
loop finds one set bit at a time. For each four-column nibble, load its
four monomial masks into one AVX2 register, permute the selected 64-bit
lanes to the front, and write one overlapping four-lane vector. Advance
the output by the nibble popcount. Allocate three extra slots so the
last overlapping write remains in capacity. Handle the final partial
nibble scalar. The hypothesis is that this reduces unpack cost enough to
improve the complete call; returned full `Vec<F2BoolPoly>` rows, term
order, raw fingerprints and all counted work must be bit for bit exact.
Keep scalar and AVX-512 paths available, with runtime AVX2 detection.

The earlier AVX-512 compress-store experiment was exact but had no gain
on an EPYC 9V45 (`research/f5_avx512_unpack_20260929/RESULT.md`). This
is a distinct AVX2 decoder for an AVX2-only host. The frozen source
reference is main commit `b82b87aacb897c15e0e84701f18b4b9169d73c75`, with
SHA-256 `bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`,
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen workloads, cost and gates

Use the seven F5 cases from `examples/f4_f2_bench.rs`, primary
`f5_n24_m24_d4`, and seed XOR values `0`, `badc0de1`, `5eed2026`,
`f5c02a28`. Use one release binary, separate arm processes, and one
Rayon thread. Both arms enable selective echelon output, fused row
counting, direct packed rows, AVX2 row XOR, Gray-code table reuse,
direct scalar unpack as the prior arm, and disable deferred-above,
word-batch and AVX-512 unpack. Only `KIC_F5_AVX2_NIBBLE_UNPACK=0|1`
differs; the new arm reports whether it actually selected AVX2 unpack.
The entire criterion, build, reduction and full unpack is charged to
`wall_ms`. Process launch, fixture construction and fingerprints are
outside it.

First check exact output on random sparse, dense and partial-word rows,
including output capacity safety under sanitizable bounds. Then run on
Linux x86-64 with AVX2/BMI2, one pinned allowed CPU, one warmup per arm,
five prior/prior A/A pairs and five alternating prior/new pairs on each
seed. Preserve full outputs and statuses for successes, failures,
timeouts and OOMs; source/binary hashes, CPU model/features, affinity,
load and Rust version; phase costs, term counts, row signatures and
counted operations. Report exact five-pair bootstrap 95% intervals,
A/A ranges and every smaller case. If the host lacks AVX2/BMI2, retain
an `unsupported_host` zero-call receipt.

The requested further 2× passes only if the frozen prior/new complete-
call median and its lower interval bound both reach 2.00, each holdout
primary beats its A/A maximum, no smaller case falls below its own A/A
minimum, and all exactness checks pass. An incremental opt-in may remain
only if the frozen median and lower interval bound exceed 1.05 with the
same correctness and guards. Otherwise archive the negative result and
remove the runtime option. This is a matrix-F5 solver-stage diagnostic,
not an IC online-time or DLP speedup.
