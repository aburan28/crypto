# N53 C1200 one-versus-four target-worker screen

This experiment keeps the same exact N53 C1200 factor base, root index,
relation collector, and complete rank-1200 precomputation. The frozen source
and release binary are identical in both arms. The reference uses one target
query worker; the parallel candidate starts four independent rotating target
searches and stops them cooperatively at the first verified relation. The
prepared index and base-log system remain target-independent. Each solver
process receives one public point, and all target workers' completed work is
charged inside that point's online interval.

`SCREEN_FREEZE.json` binds a correctness calibration pair on T025, followed
by eight previously unused public points T181–T188. The T025 pair is excluded
from the held-out ratio. The run order alternates AB/BA, with one repetition
per arm and point, 14 preparation workers, at most four target workers, a
32 GiB RSS envelope, and a 120-second process timeout. A fresh helper process
provides per-solve RSS and resource counters. The runner preserves failures
and stops at the first invalid or over-budget row. After the T025 pair, it
pauses at a frozen verified prefix so the report and resource gates can be
checked before the held-out panel starts.

The predeclared local advancement threshold is an eight-target geometric
mean of at least **2×** for one-worker/four-worker online wall time, with all
16 held-out scalars and phase sums verified. Target S3 calls, state probes,
preparation times, and phase costs remain diagnostics. The same physical CPU
resources are available to both arms; this ordinary macOS host has no
exclusive CPU partition or thread-affinity receipt. Any wall-time result is
therefore exploratory under the repository's CPU isolation rule.

## Complete local result

The T025 calibration pair and all 16 held-out runs verified the native group
equation, fixture scalar, rank 1200, exclusive target phase sum, and 32 GiB
RSS gate. Checked Sage independently replayed all 18 scalars and the
generator order. The maximum process RSS was 12,056,477,696 bytes.

| Public point | One worker online ms | Four workers online ms | One/four online | One/four S3 calls |
| --- | ---: | ---: | ---: | ---: |
| T181 | 0.462 | 0.835 | 0.554× | 618 / 3,074 |
| T182 | 0.146 | 0.241 | 0.605× | 8 / 12 |
| T183 | 0.171 | 0.257 | 0.668× | 58 / 234 |
| T184 | 1.223 | 0.449 | 2.725× | 2,368 / 794 |
| T185 | 0.699 | 0.715 | 0.978× | 1,008 / 2,474 |
| T186 | 0.563 | 0.634 | 0.888× | 818 / 1,842 |
| T187 | 0.562 | 0.352 | 1.597× | 806 / 256 |
| T188 | 0.453 | 0.436 | 1.039× | 674 / 1,120 |

The held-out geometric online ratio is **0.984×**, with a target-bootstrap
95% interval **[0.715×, 1.429×]**. The geometric ratio of total charged S3
calls is **0.655×**, interval **[0.349×, 1.318×]**. The 2× online advancement
gate failed. T184 shows the benefit on a longer search, while the short T182
search and extra work on T181/T183 show the fixed cost of four concurrent
searches. The T025 calibration ratio was 3.409×; it is excluded from these
held-out summaries.

The unchanged one-worker setup stayed within 10% of its median on all eight
held-out points. The four-worker arm ranged from 0.899× to 1.105× its setup
median, failing the 10% noise gate. The host has no exclusive partition
receipt, so `controlled_speedup` remains `null`. `analysis.json` audits each
run receipt and compressed report; `local/archive_receipt.json` binds all 77
files under `local/paired_screen/` by SHA-256. Reanalysis from the archived
directory reproduced the original JSON byte for byte.
