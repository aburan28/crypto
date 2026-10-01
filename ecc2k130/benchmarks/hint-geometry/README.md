# Reconverged-hints geometry sweep

This bounded RTX PRO 6000 experiment holds the reconverged table-v3 walk and
the per-block slot footprint fixed while changing batch geometry:

| arm | batch | block threads | runtime threads | live slots | 32-launch slot populations |
|---|---:|---:|---:|---:|---:|
| `b16` | 16 | 512 | 96,256 | 1,540,096 | 49,283,072 |
| `b32` | 32 | 256 | 48,128 | 1,540,096 | 49,283,072 |
| `b64` | 64 | 128 | 24,064 | 1,540,096 | 49,283,072 |

Every arm uses `TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 MINBLOCKS=1` and
therefore runs the same v3 step. Before timing, each arm must re-walk 300
reports with zero drops/mismatches, identify the expected runtime geometry,
and produce the same header-aware table-v3 corpus record multiset.

After one excluded warmup per arm, two equal-launch screens run in opposite
orders (`b16,b32,b64`, then `b64,b32,b16`). The faster of `b32` and `b64`
qualifies for exactly three alternating confirmation pairs against `b16` only
if its screened median is at least 1.005 times the `b16` median. An unqualified
screen stops. This is same-walk kernel engineering and does not prove 26 B/s.

## Measured result

All three arms replayed 300/300 reports with zero drops and produced the same
1,480,278-record v3 corpus. B16 compiled with 128 registers; B32 and B64 used
160. All had a 400-byte frame and zero reported spills.

The two 32-launch screens gave medians of 2.931747 B/s (B16), 2.950038 B/s
(B32) and 2.333948 B/s (B64). B32 crossed the 0.5% admission threshold, but
the longer confirmation rejected it in all three pairs:

| Pair | B16 B/s | B32 B/s | B32 / B16 |
|---|---:|---:|---:|
| 1 | 2.449343 | 2.432091 | 0.992956 |
| 2 | 2.449763 | 2.431911 | 0.992713 |
| 3 | 2.449612 | 2.431759 | 0.992712 |

B16/T512 remains selected. [`result.json`](result.json) retains the complete
screen and confirmation. The raw archive is 92,074,189 bytes with SHA-256
`06ccc4323b8d8ca4dc70ab5d3d78a2481235651af780da9f6589702c3d7d4354`.
[`independent-audit.json`](independent-audit.json) reopens every timed log,
checks its digest and exact work, and independently reaches the B16 decision.

Prepared launch; do not run without admission:

```sh
modal run --detach modal_job.py \
  --job benchmarks/hint-geometry/gpujob.sh \
  --out /tmp/ecc2k-hint-geometry --gpu RTX-PRO-6000
```
