# Seven-node inverse with five/six shared X slots on G7

Preregistered engineering screen on the local AWS g7.2xlarge / RTX PRO 4500,
GPU UUID GPU-827e82b5-4739-d339-d95e-b7214337553e, 165 W. The objective is
12 billion complete scalar ECC2K-130 rho iterations/s on one GPU. The selected
binary is `build/ecc2k130-local-packed`, SHA-256
`8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02`;
the preceding seven-node screen measured a fresh selected median of
6.182635 B/s. Generic work remains `sqrt(n/262)`, ratio 1, and full-DLP S is
null. Every iteration retains the original eight-way jump partition.

## Hypothesis

The independently verified seven-node block inverse reduces inverse scratch
from 10,624 to 4,480 bytes. Spend that 6,144-byte saving on persistent X state:
cache five or six of the eight physical slots per thread, versus four selected.
Predicted kernel shared allocations are 27,776 and 32,128 bytes, respectively,
including the unchanged 1,536-byte Frobenius table. With the driver's observed
1,024-byte reserve, both must retain three 256-thread blocks per SM. The
seven-node/four-slot implementation is the arithmetic control. The candidate
changes no field operation, jump, state format, DP rule, replay, restart or
checkpoint contract.

Build `xcache-control` with logical-pair inverse mode 2 and four cached slots,
`xcache-five` with five, and `xcache-six` with six. All use the selected CUDA
13.3 front end / 13.4 assembler, B16, split2, 256 threads and minblocks3. The
selected executable is a fourth timing row. Compile in an isolated source copy;
do not alter the selected executable or prior seven-node source archive.

Before timing, require the prior 4,096 zero-mask algebra identities and guarded
seven-node inverse proof, the independent packed field suite, actual-walk state
comparison for each new cache size, and full client state/DP/CPU replay checks
on IDs 0, 139 and fresh 48103, including partial populations and bidirectional
resume. At minimum run memcheck and shared-memory initcheck on the changed walk.
A correctness failure blocks timing.

Time three interleaved repetitions of all distinct verified rows at 524,288
workers, B16, 1,024 steps and four launches: 34,359,738,368 complete updates per
sample. Preserve every sample and hardware identity. A candidate advances only
if its paired log-ratio Student-t 95% interval against selected is wholly above
1. Promotion still requires the full validation suite and five fresh alternating
pairs in both benchmark and DP34 workloads, both confidence intervals above 1,
identical DP multisets and zero drops. Otherwise retain the selected runtime and
the negative evidence. Update this note and the canonical local scoreboard when
the round closes. No hardware, power, driver, service or cloud-resource changes.

## Result

The screen closed without a qualifying candidate. All twelve preregistered
samples completed and none was excluded. The selected median was 6.262430 B/s.
The four-slot arithmetic control measured 6.243798 B/s, paired ratio 0.995699
with 95% CI [0.987072, 1.004402]. Five slots measured 6.295066 B/s, ratio
1.000250 [0.989535, 1.011081]. Six slots measured 6.374050 B/s, ratio 1.009801
[0.991428, 1.028515]. Every interval includes 1, so the full confirmation gate
did not open and the selected executable remains unchanged. The best observed
median is 0.531171 of the 12 B/s objective.

Both new cache sizes passed an independent actual-walk/cache oracle (136
scenarios and 1,088 launches each), memcheck, shared-memory initcheck, complete
state comparison, DP multiset comparison, CPU endpoint replay and bidirectional
resume. Across selected/candidate comparisons each size covered 64 CPU endpoint
replays and 38,835 DP records with zero drops or mismatches. The binaries use
80 registers; their kernel shared allocations are 27,776 bytes for five slots
and 32,128 bytes for six, retaining three blocks per SM after the observed
1,024-byte driver reserve. Generic work remains `sqrt(n/262)`, ratio 1, and
full-DLP S remains null.
