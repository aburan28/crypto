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

One table, one unit. The 85 M it/s/SM floor is a 6000 integer-pipe remainder;
L40S rows keep the column because the receipt writes it, not because Ada can
be asked to cross a 6000 pipe. Occupancy on Ada is the `/ 4-wave` column on
that allocation.

| GPU | waves | workers | median B/s | M it/s/SM | / 4-wave | / pipe ~85 | class | correctness |
|---|---:|---:|---:|---:|---:|---:|---|---|
| RTX PRO 6000 shipping | 4 | 385024 | 15.115792 | 80.40 | — | 0.946 | reference | top-clmad, other allocation |
| RTX-PRO-6000 | 1 | 96256 | 14.479876 | 77.02 | 0.970 | 0.906 | engineering | 3/3, 104 regs, identity shipping packed, clmad 1 |
| RTX-PRO-6000 | 4 | 385024 | 14.932231 | 79.43 | 1.000 | 0.934 | engineering | 3/3, same binary, 188 SMs, 2×256 resident |
| RTX-PRO-6000 | 6 | 577536 | 14.861568 | 79.05 | 0.995 | 0.930 | engineering | 3/3, same allocation |
| RTX-PRO-6000 | 8 | 770048 | 14.844084 | 78.96 | 0.994 | 0.929 | engineering | 3/3, same allocation |
| L40S | 1 | 72704 | 8.618576 | 60.69 | 1.115 | 0.714 | engineering | 3/3, 94 regs, identity shipping packed, clmad 1 |
| L40S | 4 | 290816 | 7.727937 | 54.42 | 1.000 | 0.640 | engineering | 3/3, 142 SMs, 2×256 resident |
| L40S | 6 | 436224 | 7.464325 | 52.57 | 0.966 | 0.618 | engineering | 3/3, same allocation |
| L40S | 8 | 581632 | 7.268352 | 51.19 | 0.941 | 0.602 | engineering | 3/3, same allocation |

6000 rates from [rtx-pro-6000-waves.json](rtx-pro-6000-waves.json)
(`ap-CreQoLEbhXsEhB1glybxbo`): automatic occupancy 96,256 threads on 188 SMs.
Best wave is **4**. 6× is 0.995 of that row, 8× is 0.994. Neither is >1%
above 4×, so `fourBeaten` is false. This 4-wave median is 0.988 of the
historical 80.40 on a different card-day; that gap is not occupancy.

L40S rates from [l40s-waves.json](l40s-waves.json)
(`ap-aNsvFhZH5TKB6fETEVUfLP`): automatic occupancy 72,704 threads on 142 SMs.
Best wave is **1**, 1.115× the 4-wave row on that allocation, so
`fourBeaten` is true — 4× lost, it was not beaten from above by 6× or 8×.

Every row is **engineering**: same walk, 5.3125 products/update, grid
multiplier only. The ratio to the 6000 pipe did not fall. 6× and 8× on the
6000 slightly lowered `S` at a flat pipe ratio (~0.93). L40S 1-wave raised
B/s relative to 4-wave by cutting oversubscription, which is still
engineering.

## 4. Verdict

Falsify "more waves help the 6000": **yes**. 6× and 8× are both ≤ 4× per-SM
(0.995 and 0.994). The shipping 385,024-worker grid stays the 6000 default.

Falsify "the survey starved other SMs": **yes** on L40S. 4× / 1× is 0.897,
not within 1%, and 4× is slower. Ada wants automatic occupancy. Putting
385,024 workers on a 142-SM part is still inadmissible; so is promoting
L40S 1-wave onto the campaign 6000.

The occupancy axis on the campaign part is spent. Remaining 6000 headroom
is instruction cuts against the pipe ([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md)).
Do not move the live searchers. This is `--bench`; it does not promote a
collecting rate.
