# AVX2 construction of GF(2) Gray-code tables

## Frozen hypothesis and reference

The accepted one-thread matrix-F5 fast path uses AVX2 for row updates on
the EPYC 7763 host but still constructs each Gray-code table entry with a
generic per-word XOR loop. Table entries are rebuilt for every pivot block.
Compile that exact `dst = previous_entry XOR pivot_row` loop for AVX2 and
dispatch only after a runtime AVX2 check. The hypothesis is that faster
table construction reduces the complete F5 call. Pivot order, table
contents, rank, full `Vec<F2BoolPoly>` output, raw row fingerprints and
counted word XORs must remain identical. Retain the portable scalar path.

The frozen source reference is main commit
`bd9815732611edd886075902865ffc12e1971810`. SHA-256 hashes are
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`. The prior eligible x86 direct-unpack
complete-call marginal median was 122.96 ms on its host in
`research/f5_further_2x_20260929/RESULT.md`; use only ratios from the
new paired run for this decision.

## Workloads, arms and accounting

Run the seven F5 cases in `examples/f4_f2_bench.rs`, primary
`f5_n24_m24_d4`, with seed XOR values `0`, `badc0de1`, `5eed2026` and
`f5c02a28`. Use one release binary and separate prior and new arm
processes. Both arms enable selective echelon output, fused row counting,
direct packed rows, forced AVX2 row XOR, table reuse, direct scalar
unpack, one Rayon thread, and disable deferred-above and word-batch
options. Only `KIC_GF2_AVX2_TABLE_BUILD=0|1` differs. The arm-1 route
must report that AVX2 table construction was actually selected on the
primary case. Charge table construction, elimination, criterion, build
and full unpack to the complete call. Fixture generation, process launch
and fingerprinting are outside `wall_ms`.

First run the portable and x86 exactness tests. Then use a Linux x86-64
host reporting AVX2 and BMI2, pin one allowed CPU, warm each arm once,
take five prior/prior A/A pairs and five alternating prior/new pairs per
seed. If the host lacks the features, preserve an `unsupported_host`
zero-call receipt. Record every process status, failure, timeout and OOM,
full output, phase times, binary/source hashes, CPU model/features,
affinity, load, Rust version and counted work. Report five-pair medians,
exact five-pair bootstrap 95% intervals, A/A ranges, all holdouts and
smaller cases. Preserve the CI artifact and a compact raw receipt in
the PR.

The requested further 2× passes only if the frozen prior/new full-call
paired median and lower interval bound both reach 2.00, every holdout
primary beats its own A/A maximum, no smaller-case median falls below
its own A/A minimum, and all exactness checks pass. An incremental
opt-in may remain only if the frozen median and lower bound exceed 1.03
with the same correctness and guard checks. Otherwise archive the
negative result and remove the experimental runtime option. The
measurement is a matrix-F5 solver-stage diagnostic, not an IC online-time
or DLP speedup.
