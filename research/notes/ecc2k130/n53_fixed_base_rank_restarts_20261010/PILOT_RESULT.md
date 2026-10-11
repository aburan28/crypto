# n53 fixed-base rank-restart pilot

All nine preregistered cells reached rank 220/220 and independently replayed
the same public target, recovering scalar `17385600002971`. The actual base
contained 23,320 subgroup-usable points in 220 signed-Frobenius columns, with
digest `7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973`.
The 400,000-probe cap had the lowest median complete cold cost of the two
capped variants, so it is the preregistered held-out candidate. The unbounded
control's median cold cost was lower in this pilot.

| Rank seed | Arm | Total rank probes | Rank attempts | Capped attempts | Complete cold ms | Rank PDP ms | Target online ms |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 530053 | unbounded | 18,413,613 | 220 | 0 | 5,912.214 | 4,526.345 | 6.335 |
| 530053 | 200,000 | 16,582,773 | 237 | 17 | 10,881.242 | 8,358.424 | 6.514 |
| 530053 | 400,000 | 17,913,901 | 224 | 4 | 6,940.423 | 5,315.598 | 14.793 |
| 530054 | unbounded | 18,661,879 | 220 | 0 | 6,144.949 | 4,807.309 | 6.753 |
| 530054 | 200,000 | 16,993,325 | 240 | 20 | 6,704.194 | 5,020.686 | 7.857 |
| 530054 | 400,000 | 17,104,035 | 222 | 2 | 5,936.396 | 4,371.030 | 6.043 |
| 530055 | unbounded | 16,807,264 | 220 | 0 | 6,466.231 | 4,831.133 | 6.518 |
| 530055 | 200,000 | 17,763,255 | 240 | 20 | 6,910.439 | 5,326.872 | 7.004 |
| 530055 | 400,000 | 16,807,264 | 220 | 0 | 6,575.210 | 4,937.373 | 6.490 |

Across the three seeds, median complete cold costs were 6,144.949 ms
(unbounded), 6,910.439 ms (200,000), and 6,575.210 ms (400,000).
Median paired probe reductions against the unbounded arm were 8.941% for
200,000 (seed range -5.688% to 9.942%) and 2.714% for 400,000 (0% to
8.348%). The corresponding median paired cold-cost ratios were approximately
1.091 and 1.017; the full per-seed ratios are in `PILOT_ANALYSIS.json`.
These three-seed ranges show the uncertainty directly. Timing is exploratory
on this shared Apple T6041 host; no CPU isolation receipt was available.
The 200,000 arm on seed 530053, for example, had fewer probes but higher
index and rank-PDP times, so probe savings do not determine wall cost alone.

`analyze_pilot.py` recomputes `PILOT_ANALYSIS.json` from each cell's raw
summary, base, target, rank trace, independent replay receipts, hashes, and
process status. Every cell used the committed source SHA-256
`4dd3c5e7823b453ad1fa9046d9be68a04bddc79343ea29576dd1f32b0bc66299`
and executable SHA-256
`3c03013e4a12d5eac74e1737c5549558ddf734aa9a324426914778765da0f51e`.
The producer had a 60-second wall cap and 16-GiB observed RSS cap per cell;
the nine recorded peaks were 332.06–334.31 MiB. Each cell has its own
`intent.json`, raw JSONL, stderr, two replay receipts, and `status.json` under
`runs/pilot/<seed>/<arm>/`. The source, input, protocol, and candidate/workload
identities were committed and pushed before the first timed cell.

The next preregistered step generates the new point with `hash:53261307`,
then pairs the 400,000-probe cap and unbounded arm on six independent rank
seeds and same-point rho seeds. Its advancement gate is a median paired total
rank-probe reduction of at least 15% without higher paired complete cold cost,
with full rank and target verification in every cell.
