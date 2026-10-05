# The `ℓ = 6` rung with a dense in-tree engine: results

Written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since its
registration commit `e215db702`. The readout is
[runs/registered/readout.txt](runs/registered/readout.txt), printed by the unchanged
`examples/gb_ladder_analyze`.

## Registered verdict: growing

By §6's rule the verdict is **growing**. Both arms resolve at `ℓ = 6` on both `n = 19`
curves, on every draw, exactly at the registered law: the direct `S₄` descent (`x4`)
refutes at **11** and the symmetric norm form (`rr`) at **8**. No reading on either curve
is below them. The `x4` excess over the sharp first-fall bound 7 is now **4**, one more
than at `ℓ = 5`, on the canonical object. The first-fall assumption stays four degrees
behind the data, and the stop decision's item 3 stays closed.

| `ℓ` | cell | `rr` | `x4` | `x4 − 7` | control |
|--:|:--|:--|:--|:--|:--|
| 2 | `K₁/2¹⁷` | 4 triv 4 4 | 5 5 5 6 | −2 | ≥8 ≥8 ≥8 ≥8 |
| 3 | `K₁/2¹⁷` | 5 5 5 5 | 8 8 8 8 | 1 | ≥8 ≥8 ≥8 ≥8 |
| 4 | `K₁/2¹⁷` | 6 6 6 6 | 9 9 9 9 | 2 | ≥8 ≥8 ≥8 ≥8 |
| 5 | `K₁/2¹⁷` | 8 triv 8 8 | 10 10 10 10 | 3 | ≥8 ≥8 ≥8 ≥8 |
| **6** | `K₁/2¹⁹` | **8 8 8 8** | **11 11 11 11** | **4** | ≥8 ≥8 ≥8 ≥8 |
| **6** | `K₀/2¹⁹` | **8 8 8 8** | **11 11 11 11** | **4** | ≥8 ≥8 ≥8 ≥8 |

`triv` marks the two trace-constant draws, whose system contains the constant `1` and
refutes at degree 1; they are excluded from medians, as registered. A control reading
`≥8` means the scan reached the cap of 7 without refuting; this engine does not detect
pinning, so that is all a control reading says.

Slopes of the exact medians on `K₁` over `ℓ = 2…6`: `rr` 1.10, `x4` 1.40. From `ℓ = 3` on,
`x4` rises by exactly one degree per rung: 8, 9, 10, 11.

## Predictions

1. **P1, calibration: holds.** On all 32 `rr` and `x4` draws of the four `K₁/2¹⁷`
   cells, this engine's reading equals the external engine's
   ([ic_gb_ladder_20261003](../ic_gb_ladder_20261003/RESULTS.md)), including the two
   trivial draws. Against the in-tree ladder, the readout reports 20 agreements and 3
   disagreements. All 3 are the in-tree floor of 6 at `x4 ℓ = 2`, which the earlier round
   already recorded. Paired `rr − x4` over the 22 draws exact on both arms is −2.64
   (range −3 to −1).
2. **P2, `x4` at `ℓ = 6`: holds.** The median is exact and equals 11 on each `n = 19`
   curve, on all eight draws. Degree 10 is unrefuted on every draw, and degree 11
   refutes on every draw.
3. **P3, `rr` at `ℓ = 6`: holds.** The median is exact and equals 8 on each curve, on all
   eight draws.
4. **P4, controls: holds.** No control is refuted at or below 7 on any of the six cells,
   24 draws in all.
5. **P5, pairs: fails as registered, through a slip in the registration.** `rr` is below
   `x4` on every draw exact on both arms. The gap is 2 or 3 at every `ℓ ≥ 3`, including
   3 on all draws at `ℓ = 6`. At `ℓ = 2` it is 1 on two of the three paired draws. The registration said
   "by 2 to 4", but the earlier round had already measured a gap of 1 at `ℓ = 2` and
   registered its own P5 as "by 1 to 3". The failure reflects how the prediction was
   written. It carries no new information about the systems.

## Run

| | |
|:--|:--|
| engine | `examples/macaulay_dense.rs` at commit `bafb36a4f`, unchanged through the run |
| systems | the earlier round's exports, regenerated from `rr_degree_ladder --dump-dir` at seed 20260930 and checked against [dump.sha256](../ic_gb_ladder_20261003/runs/registered/dump.sha256) before use, on every host |
| serial phase | 2026-10-04 14:35Z to 2026-10-05 02:08Z, [run.sh](run.sh), one process at a time, four threads: calibration cells, the `ℓ = 6` `rr` arms, the `K₁/2¹⁹` `x4` arm |
| parallel phase | 2026-10-05 01:48Z to 05:22Z, [lanes.sh](lanes.sh): the controls in a one-thread lane, `K₀/2¹⁹` `x4` draws 0 and 2 in two two-thread lanes on this host, draws 3 and 4 on two other 4-core hosts |
| remote records | [runs/registered-remote-k0d3](runs/registered-remote-k0d3/) and [runs/registered-remote-k0d4](runs/registered-remote-k0d4/), as pushed by those hosts, with their `HOST.txt`; their lines are folded verbatim into `runs/registered/K0n19l6.jsonl` for the readout |

### Departures from the registration, disclosed

- **Two host restarts.** The container restarted at about 22:50Z and again at about
  00:05Z. Both times the `K₁/2¹⁹` draw 3 `x4` degree-11 process was the one running. No
  line was lost; the runner resumed from its log and recomputed that degree from scratch.
  The interrupted processes left no line, so they are not counted as censored. A restart
  is not a budget kill. After each restart, `progress.txt` repeats the "finished" lines
  of the cells it skipped.
- **Parallel lanes.** The registration says "one process at a time". From 01:48Z the
  remaining units ran concurrently to finish sooner. Each unit still ran as one process
  per degree, with the same `measure` function, log schema and per-process memory cap.
  This affects only wall-clock time; the degree readings are deterministic.
  [lanes.sh](lanes.sh) claims each unit with a lock. In a test with three competing
  lanes, every unit ran exactly once and the readings matched the serial runner.
- **CPU budget.** The registration says 86,400 CPU-s per process. `run.sh` set 43,200 and
  `lanes.sh` sets 86,400. No process came near either limit. The largest, a degree-11
  `x4` process, used about 20,000 CPU-s.

## What the engine reached, and what it cost

| arm, degree | columns | rows | rank at refutation | wall, four threads |
|:--|--:|--:|--:|:--|
| `rr` 8, `ℓ = 6` | 169,766 | 285,552 | 98,037–98,038 | 15–26 min |
| `x4` 10, `ℓ = 6` (unrefuted) | 199,140 | 76,912 | 76,912 (full) | 5–7 min |
| `x4` 11, `ℓ = 6` | 230,964 | 239,704 | 188,495–192,255 | 80–122 min |
| `x4` 10, `ℓ = 5` | 30,827 | 32,997 | 26,586–26,907 | 18–24 s |

The two `K₀/2¹⁹` draws run locally on two threads each, beside each other, took 3.0 h
at degree 11. Peak resident memory for a degree-11 process was about 5.5 GB. Singular's truncated
`slimgb` died at 12 GB at degree 8 (`rr`) and 9 (`x4`) on these systems. The "rank at
refutation" is the rank when the constant row first appeared; the scan stops there.

## What this says

- **The law reaches one more rung.** On both `n = 19` curves the direct descent needs
  degree `ℓ + 5 = 11`, four above the sharp first-fall bound. Its excess over that bound
  is −2, 1, 2, 3, 4 at `ℓ = 2…6`: it rises by one per rung from `ℓ = 3` and shows no sign
  of flattening. The symmetric form stays three below the direct descent at `ℓ = 5` and
  `ℓ = 6`.
- **It is the same number on both curves.** `K₀` and `K₁` give identical readings
  draw for draw at `ℓ = 6` on both arms, and identical matrix sizes.
- **Scope.** This is a stage diagnostic (AGENTS.md §2, §5): two Koblitz curves, `n` 17–19,
  `ℓ` 2–6, `m = 3`, four rootless draws per cell, one engine. It says nothing about
  end-to-end cost, nothing at `n ≈ 83` or 131, and nothing asymptotic. At `m = 3`
  nothing here could tie rho even with a free oracle.

## Next

The registration's next step after "growing" is the `ℓ = 7` rung. At degree 9 that is
about 1.2 million columns, a dense basis near 180 GB, which no host available here can
hold. Going further needs a sparse or structured elimination, not a larger dense one.
