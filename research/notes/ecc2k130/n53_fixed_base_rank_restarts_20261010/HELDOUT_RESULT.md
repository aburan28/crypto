# n53 fixed-base rank restarts: six-seed held-out result

The 400,000-probe rank-restart policy reduced median paired total rank probes
by **5.189%** on the new public point, with a median paired complete-cold
cost ratio of **1.020** relative to the unbounded control. It did not pass
the preregistered advancement gate of at least 15% probe reduction without
higher paired cold cost. Both IC arms reached rank 220/220 and independently
verified the same recovered scalar in all six pairs. The same-point strong-rho
arm also independently verified that scalar in all six runs.

The point is `[2939726529610565,1384157319972946]` on the n53 `a=0`
Koblitz curve (`r=21044858204113`), derived from `hash:53261307` at counter 1
after the generation rule was pushed. Independent group arithmetic verified
the point and recovered scalar `11184763905218`; the IC and rho processes
received Q only. The fixed base held 23,320 subgroup-usable points in 220
signed-Frobenius columns, digest
`7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973`.

| Rank seed / rho seed | Control probes | 400k probes | Paired probe reduction | Control cold ms | 400k cold ms | Control target online ms | 400k target online ms | Rho online ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 531053 / 530153 | 19,737,459 | 18,637,430 | 5.573% | 10,088.708 | 13,125.130 | 2.342 | 3.337 | 188.496 |
| 531054 / 530154 | 20,084,341 | 16,486,220 | 17.915% | 8,914.042 | 9,470.037 | 2.712 | 2.462 | 258.650 |
| 531055 / 530155 | 17,998,946 | 18,314,535 | -1.753% | 8,869.200 | 7,886.839 | 3.535 | 4.131 | 147.745 |
| 531056 / 530156 | 18,312,430 | 16,456,095 | 10.137% | 9,116.698 | 7,653.340 | 2.672 | 2.387 | 215.297 |
| 531057 / 530157 | 19,080,823 | 18,163,947 | 4.805% | 9,655.798 | 9,844.086 | 2.696 | 2.616 | 132.403 |
| 531058 / 530158 | 17,248,790 | 17,248,790 | 0.000% | 8,347.347 | 8,518.139 | 2.836 | 2.652 | 173.012 |

The cap abandoned 2, 3, 2, 1, 2, and 0 rank searches respectively, then
drew the next scalar from the unchanged stream; every attempted support probe
and failed search is charged. Across the six rank streams, unbounded search
used 112,462,789 probes and the capped search used 105,307,017. Paired probe
reductions ranged from -1.753% to 17.915%. A descriptive exact six-draw
bootstrap interval for the median is -0.877% to 14.026%. Paired cold-cost
ratios ranged from 0.839 to 1.301, with a descriptive bootstrap interval for
their median of 0.864 to 1.182. These intervals describe the six frozen seed
streams; they are not a controlled host-performance interval.

The primary one-target online interval starts when target-dependent IC query
work begins after the reusable base/index/log system is ready, and stops after
the recovered scalar is checked. Its five exclusive phases are target query,
PDP, relation check, descent, and recovery check. The matched rho interval
starts at its first target-dependent walk computation and stops after scalar
verification; it comprises walk, collision, and recovery-check time. The
median measured online intervals were 2.704 ms (unbounded IC), 2.634 ms
(400k IC), and 180.754 ms (rho). Median *paired observed* rho/IC online
ratios were 70.753 and 60.865 respectively. These are exploratory timings on
a shared Apple T6041 host, without the required isolation receipt; the
controlled online-speedup field is `null` in `HELDOUT_ANALYSIS.json`.
Complete cold IC accounting separately includes base/index construction,
every rank attempt, relation checks, matrix construction and final LA, target
work, and scalar check. The median unpaired complete-cold costs were
9,015.370 ms (control) and 8,994.088 ms (400k); the preregistered gate uses
the paired ratio shown above.

`analyze_heldout.py` reproduces `HELDOUT_ANALYSIS.json` from all 18 raw
process records and verifies source/binary/input hashes, the fixed base,
probe and cap accounting, rank transitions, exclusive phase sums, and every
independent rank, IC-target, and rho-target replay receipt. Each cell has an
immutable `intent.json`, raw JSONL, stderr, exit, external wall, wait4 peak
RSS, and `status.json` under `runs/heldout/<rank-seed>/<arm>/`. All 18 cells
exited zero below the 60-second wall and 16-GiB observed RSS caps. The exact
source, public Q, verifier-only fixture, generator amendment, and six workload
IDs were committed and pushed before the first timed cell.

The result directs the next fixed-base experiment toward reducing cost per
support probe or changing the pair-root index policy. The 400k restart acts
on only ten of the 1,330 guided rank attempts across the six capped runs, so
its effect is bounded by the rare long searches while extra attempts can
erase their savings. A new candidate should freeze a matched full-cold and
target-online comparison before promotion. The degree-263 transported-base
question remains an equal-useful-size, same-target experiment with explicit
maps and charged transport.
