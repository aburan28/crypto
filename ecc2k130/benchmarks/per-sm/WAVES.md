# Per-SM occupancy — oversubscribed waves

The packed walk on an RTX PRO 6000 is already against the integer pipe
([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md)): 74.1 lane-instructions per
SM-clock, 95% of the best mix the hardware-limits probe ever measured.
Instruction cuts are the only way past that ceiling. Occupancy is a
different axis: the CUDA occupancy API reports **two** 256-thread blocks
resident per SM, and `autoThreads` launches exactly that many workers.
The shipping preset launches **four times** that grid (385,024 on 188 SMs).
GOAL28 measured 1× / 2× / 4× on G7 and kept 4×. **6× and 8× were never
measured**, and no other Modal SKU was measured above 1×.

This note is that sweep. It does not change the campaign default. Live
6000 searchers stay on 385,024 workers.

`MINBLOCKS=3` is a different lever (ptxas register budget). GOAL28 already
rejected it: three resident blocks, 80 registers, −2% rate. This sweep
keeps `MINBLOCKS=2` and only changes the launch grid.

## 1. How to measure

```bash
make bench-waves-modal                         # 6000, waves 1,4,6,8
WAVES_GPU=L40S make bench-waves-modal          # Ada, same multipliers
bash benchmarks/per-sm/run.sh                  # log to /tmp, then freeze
```

One Modal allocation, one binary, four grids. Automatic occupancy is
probed with `--steps 1`, then each wave is `--threads $(auto×wave)` with
the shipping `--bench` geometry: batch 16, 256-thread blocks, 1024 steps,
32 launches, 3 repeats. Logs stay under `/tmp` until Modal returns.

## 2. Boundaries, written before the run

Unit: **M complete scalar updates per SM per second**. That makes the
shipping 6000 a flat ~80.4 at every wave count that saturates the pipe,
and turns "did occupancy move the SM" into reading one column. B/s stays
in the table as the campaign unit.

| | value | kind |
|---|---|---|
| Floor (integer-pipe remainder on a 6000 SM) | ≈ 85 M it/s/SM | 78.3/74.1 × 80.4; occupancy cannot cross the pipe |
| One-add walk floor, per 6000 SM | ≈ 122–133 M it/s/SM | 23–25 B/s / 188; different iteration function |
| Reference, 4 waves | **80.40 M it/s/SM** (15.115792 B/s) | [top-clmad](../top-clmad/summary.json), 385,024 workers |
| Automatic, 1 wave | 77.0 M it/s/SM (14.470 B/s) | survey, same recipe, 96,256 workers |
| GOAL28 2 waves | slower than 4 | occupancy-summary, G7; not re-cited as a 6000 number |

Class: engineering (same walk, 5.3125 products/update). Occupancy that
raises B/s at a flat ratio to the pipe floor is latency hiding, not an
advance.

Success is a verified median on the requested GPU, all four waves
complete, identity matching the shipping packed build, automatic
occupancy probed on that same allocation. A wave beats 4× only if its
per-SM median is **> 1%** above the 4-wave row **on that allocation**
(not against the historical 80.40, which is a different card-day).
Inadmissible: `MINBLOCKS=3`, changing batch or the walk, 385,024 workers
on a non-188-SM part, quoting a different SKU as occupancy, skipping
identity, or promoting the campaign from `--bench`.

Falsify "more waves help the 6000" if 6× and 8× are both ≤ 4× per-SM.
Falsify "the survey starved other SMs" if L40S 4× / 1× is within 1%.

Do not move the live searchers.

## 3. Receipt

Frozen in this directory. Cite those files; do not recompute.

| GPU | waves | workers | median B/s | M it/s/SM | / 4-wave | / pipe ~85 | class | correctness |
|---|---:|---:|---:|---:|---:|---:|---|---|
| RTX PRO 6000 shipping | 4 | 385024 | 15.115792 | 80.40 | 1.000 | 0.946 | reference | top-clmad, other allocation |
| RTX-PRO-6000 | 1 | *not yet run* | | | | | | |
| RTX-PRO-6000 | 4 | *not yet run* | | | | | | |
| RTX-PRO-6000 | 6 | *not yet run* | | | | | | |
| RTX-PRO-6000 | 8 | *not yet run* | | | | | | |
