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

Prepared launch; do not run without admission:

```sh
modal run --detach modal_job.py \
  --job benchmarks/hint-geometry/gpujob.sh \
  --out /tmp/ecc2k-hint-geometry --gpu RTX-PRO-6000
```
