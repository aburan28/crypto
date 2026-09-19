# Modal GPU survey — one packed walk, every SKU

The campaign default remains the **RTX PRO 6000** Blackwell workstation part,
`sm_120`, 188 SMs: 15.115792 B/s on the matched TOP_CLMAD control
([benchmarks/top-clmad](../top-clmad/summary.json)) and ~14.1 B/s collecting.
Modal also rents T4, L4, A10, L40S, A100-40GB, A100-80GB, H100, H200, B200
and B300. This note is one table of those parts under one recipe. It does
not change the campaign default.

The walk is integer/carryless-bound, not tensor-core-bound
([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md)). LLM SKU rankings do not
transfer. SM count, clock, and whether `clmad` exists on the part do.

## 1. How to measure

```bash
bash benchmarks/modal-gpus/run-all.sh          # every catalog GPU, sequential
bash benchmarks/modal-gpus/run-one.sh L40S     # one SKU
make bench-modal-gpu SURVEY_GPU=H100!          # same recipe; automatic workers
make bench-modal-gpu SURVEY_GPU=T4 SURVEY_CLMAD=0
```

Same geometry as Ada ([ADA-L4-L40S.md](../../ADA-L4-L40S.md)): batch 16,
256-thread blocks, two resident blocks, compact storage, shared masks,
CUDA 13.3.1. Workers stay **automatic**. The 6000 preset's 385,024 is four
waves on a 188-SM part and has no claim on any other SM count;
`autoThreads` sizes the grid from `multiProcessorCount`.

`PACKED_CLMAD` is on except T4 (`sm_75` has no `clmad`). Ampere and Hopper
are sm_80+, so the device path emits `clmad`
([include/packed131.h](../../include/packed131.h)).

The catalog is [catalog.json](catalog.json). Aliases that Modal may upgrade
(`A100`, `H100`, `B200+`) are skipped so a row is one SKU: `A100-40GB` /
`A100-80GB`, `H100!`, `B200` / `B300`. `freeze.py` refuses a log whose
nvidia-smi name or compute capability does not match the requested type.

Writes under [benchmarks/modal-gpus](.). Logs go to `/tmp` until Modal
returns: the image mount is this tree, and writing into it during
`modal run` raises ExecutionError. The page cites the frozen JSON and
does not recompute.

Live RTX PRO 6000 searchers are not this survey. Do not stop them.

## 2. Boundaries, written before the run

Unit: billions of completed scalar updates per second. Class: engineering
(same walk, 5.3125 products/update).

SM-count prior: `15.115792 × publishedSMs / 188`, equal clock and per-SM
issue to the 6000. It is a prior, not a prediction that will hold: Ada
CLMAD already beat it on L40S, B200 fell short of it on clock.

Modal list prices from [modal.com/pricing](https://modal.com/pricing)
(2026-09-19), used only for the it/$ practicality column. Break-even is
the `--bench` rate that matches a 6000 collecting at 14.1 B/s and
$3.0312/h.

| Modal type | SMs (published) | SM prior B/s | $/h | it/$ break-even vs 14.1 B/s 6000 |
|---|---:|---:|---:|---:|
| T4 | 40 | 3.216 | 0.590 | 2.75 |
| L4 | 58 | 4.663 | 0.799 | 3.72 |
| A10 | 72 | 5.789 | 1.102 | 5.12 |
| L40S | 142 | 11.417 | 1.951 | 9.08 |
| A100-40GB | 108 | 8.684 | 2.099 | 9.76 |
| A100-80GB | 108 | 8.684 | 2.498 | 11.62 |
| RTX-PRO-6000 | 188 | 15.116 | 3.031 | 14.1 |
| H100! | 132 | 10.613 | 3.949 | 18.37 |
| H200 | 132 | 10.613 | 4.540 | 21.11 |
| B200 | 148 | 11.900 | 6.250 | 29.07 |
| B300 | 160 | 12.865 | 7.099 | 33.02 |

The Hopper and Blackwell-datacenter break-evens sit at or above the
6000's own integer-pipe ceiling (~23–25 B/s on 188 SMs). A part that
merely matches 15.116 B/s still loses on Modal iterations per dollar
once its list price exceeds the 6000's. The rates still have to be
measured: a per-SM surprise would show up against the SM prior, and
would still be engineering unless it also beat 15.116.

| | value | kind |
|---|---|---|
| Floor (one-add walk on a 188-SM 6000) | ≈ 23–25 B/s | 6000 number; other SM counts cannot be asked to cross it |
| Reference | 15.115792 B/s on RTX PRO 6000, 385,024 workers | [top-clmad](../top-clmad/summary.json) |
| Collecting reference | ~14.1 B/s on the same 6000 | campaign |
| Modal list, 6000 | $3.03 / h | $0.000842 / s |

Success is a verified median on a GPU whose nvidia-smi name and compute
capability match the requested catalog entry, identity matching the
requested build (including `packed native carryless multiply` equal to
the catalog `clmad` bit), automatic workers, three valid `--bench`
repeats. Inadmissible: filing an upgraded SKU under a pinned name
(H200 as H100!, A100 80 GB as A100-40GB, B300 as B200), using 385,024
workers on a non-188-SM part, skipping identity, quoting tensor-core
marketing, or promoting the campaign onto another SKU from `--bench`
without a collecting rate.

Falsify "GPU X is faster than the 6000" if that row's median is below
15.115792 B/s. Falsify "GPU X is cheaper per iteration than a collecting
6000" if its it/$ is below 14.1e9 × 3600 / 3.0312. Record unavailable
rather than invent a rate.

## 3. Receipt

Frozen in this directory. Cite those files; do not recompute. The 6000
shipping row is the reference, not a survey number; the survey also
re-measures a 6000 with **automatic** workers so occupancy matches the
other rows.

| variant | median B/s | / SM prior | / 6000 15.116 | T it/$ | class | correctness |
|---|---:|---:|---:|---:|---|---|
| RTX PRO 6000 shipping, 385k workers | 15.115792 | 1.000 | 1.000 | 17.96 | reference | top-clmad control |
| T4 CLMAD=0 | *not yet run* | | | | | |
| L4 CLMAD=1 | *not yet run* | | | | | |
| A10 CLMAD=1 | *not yet run* | | | | | |
| L40S CLMAD=1 | *not yet run* | | | | | |
| A100-40GB CLMAD=1 | *not yet run* | | | | | |
| A100-80GB CLMAD=1 | *not yet run* | | | | | |
| RTX-PRO-6000 automatic CLMAD=1 | *not yet run* | | | | | |
| H100! CLMAD=1 | *not yet run* | | | | | |
| H200 CLMAD=1 | *not yet run* | | | | | |
| B200 CLMAD=1 | *not yet run* | | | | | |
| B300 CLMAD=1 | *not yet run* | | | | | |

Collecting rates on every non-6000 SKU are unmeasured. Do not move the
live searchers.
