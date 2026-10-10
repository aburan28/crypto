# Independent N53 C1200 combined-kernel confirmation

This panel uses T197–T204, separate from the T189–T196 selection screen.
T197–T200 were already frozen in the original public-point law. Checked Sage
extends that exact SHA-256 scalar law to T201–T204 and verifies the existing
T200 fixture before writing the four new public points and workload records.
`GENERATED_TARGETS.json` and `sage_runtime_info.json` bind that step.

Each point receives two independent processes per arm in ABBA or BAAB order,
alternating by point. Both arms use 14 preparation workers, one target worker,
the same 32 GiB RSS envelope, and identical one-target timing boundaries.
The reference and combined candidate keep their already screened source,
binary, factor-base, and algorithm identities. The runner stops on the first
invalid or over-budget run and retains all completed rows and raw reports.

The point, rather than a process repeat, is the independent sampling unit.
For each point, the online ratio is the geometric mean of the two reference
times divided by the geometric mean of the two combined times. The frozen
local confirmation gate requires all 32 native and independent Sage scalar
replays, the same first relation in every pair, a geometric online ratio of
at least **2×** across eight points, and a target-bootstrap 95% lower
endpoint above 1×. A 95% lower endpoint above 2× is reported separately as
evidence that the full 2× threshold is resolved rather than just reached by
the point estimate. Per-arm A/A spread and setup-time drift are independent
host-noise checks. A controlled wall-time claim still requires host-level
CPU isolation, affinity, frequency, NUMA, IRQ, and contention receipts.

## Complete confirmation result

All 32 one-target runs passed native group verification, scalar fixture
checks, complete rank 1200, exclusive phase accounting, and the 32 GiB RSS
gate. Checked Sage independently replayed every scalar and the generator
order. The first relation witness was identical across all four runs for
each point. Maximum peak RSS was 12,067,045,376 bytes.

| Point | Reference online ms, R1/R2 | Combined online ms, R1/R2 | Target-level ratio |
| --- | ---: | ---: | ---: |
| T197 | 0.512 / 0.532 | 0.278 / 0.230 | 2.065× |
| T198 | 0.321 / 0.353 | 0.211 / 0.181 | 1.721× |
| T199 | 0.405 / 0.461 | 0.187 / 0.190 | 2.292× |
| T200 | 0.602 / 0.613 | 0.295 / 0.291 | 2.072× |
| T201 | 0.876 / 1.161 | 0.439 / 0.449 | 2.271× |
| T202 | 2.159 / 2.202 | 2.127 / 1.139 | 1.401× |
| T203 | 3.814 / 2.921 | 1.518 / 1.711 | 2.072× |
| T204 | 2.573 / 2.964 | 1.375 / 1.447 | 1.958× |

The independent point-level geometric online ratio is **1.960×** with a
target-bootstrap 95% interval **[1.748×, 2.148×]**. The frozen local
confirmation gate requires a point estimate of at least 2×, so it failed.
The 95% lower endpoint also falls below 2×. The geometric reference/combined
S3-call ratio is 0.992×; the candidate charges its final-batch look-ahead
calls. The T202 combined repeats differ by 1.87× despite identical S3 calls
and witness, making their timing variation visible in the complete ledger.

The 95th percentile A/A absolute log ratio corresponds to a 1.319× repeat
spread for the reference and 1.603× for the combined candidate. Setup times
also exceeded the 10% noise gate in both arms, reaching 1.410× and 1.585×
their respective medians. The host has no exclusive CPU partition receipt;
the controlled speedup remains unknown. `analysis.json` verifies every
receipt, raw report, phase sum, run ID, scalar, first witness, and independent
Sage replay. Reanalysis of the archived directory reproduced the result
byte for byte; `local/archive_receipt.json` binds all 133 raw files.
