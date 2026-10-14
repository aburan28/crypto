# DP32 ordinary-sigma table screen

Run: `20261012T1815Z`

Provider: RunPod, pod `l1of60e0sm1spt`

GPU: NVIDIA RTX PRO 6000 Blackwell Server Edition

Driver: 580.82.07

CUDA: 13.3.73

Source revision: `615381ba5c5fa89236a228cd75338d636fa81903` with the experiment changes uncommitted on the pod

## Decision

The exact polynomial-coordinate sigma table is correct for the ordinary witness-carrying DP32 walk, but this screen does not justify a production rollout or the 20 billion iterations/s target. The candidate median was **9.254351 B/s** on this server-class RTX PRO 6000: **+0.7947%** over the matched 512-thread control, but **−1.0330%** against the live 256-thread geometry. The live arm was fastest at **9.350950 B/s** on this card.

This is a same-GPU result. It must not be compared directly with the existing workstation result or with the separate table-walk benchmark: those use different hardware and/or a different campaign iteration function.

## Arms

All arms use the ordinary sigma walk (`WALK_TABLE=0`, `SIGMA_FUSED=0`, `WITNESS=1`), DP32, 1,024 steps/launch, and 385,024 benchmark threads (`16` slots/thread = 6,160,384 walks). The only candidate code change is using the exact polynomial-coordinate sigma table for the ordinary weighted-prefix path.

| Arm | Geometry | Table | Median | Samples |
| --- | --- | --- | ---: | ---: |
| live | 256 threads/block, min 2 blocks/SM | off | 9.350950 B/s | 9.348547–9.354615 |
| geometry | 512 threads/block, min 1 block/SM | off | 9.181391 B/s | 9.179768–9.181455 |
| candidate | 512 threads/block, min 1 block/SM | on | 9.254351 B/s | 9.252145–9.255851 |

The timing panel used one warmup per arm followed by three rotating A/B panels; each arm therefore has three recorded timing samples. The raw per-launch logs, hashes, GPU clock/power/temperature samples, and derived summary are in [`raw/run11`](raw/run11). The graph source is [`data.tsv`](data.tsv), and the rendered graph is [`dp32-rates.svg`](dp32-rates.svg).

## Correctness gate

- CUDA arithmetic suite: passed for all three arms.
- Each arm: 5 DP32 distinguished points, 5/5 CPU-reference verified, 0 dropped.
- Sorted v2 corpus identity: live is the reference; geometry and candidate are identical to it.
- Candidate resource receipt: 126 registers/thread, 68,352 dynamic shared bytes, 1 block/SM at 512 threads/block.
- Matched control resource receipt: 91 registers/thread, 1,792 shared bytes/block, 1 block/SM at 512 threads/block.

The gate is deliberately a screen, not a full campaign qualification. A prior 300-report collection attempt timed out at the explicit 45-minute limit because CPU replay/verification became the bottleneck at DP32; that failed receipt is preserved in [`raw/failed-run9`](raw/failed-run9). The earlier short gate is preserved in [`raw/failed-run8`](raw/failed-run8).

## Reproduction

The isolated RunPod job is [`gpujob.sh`](../../../../ecc2k130/benchmarks/dp32-sigma-table/gpujob.sh). Its source receipt records the exact file hashes and the RunPod host metadata. The code can be built locally with:

```sh
make -B -C ecc2k130 gpu-rtx-pro6000-dp32-sigma-table
```

The target is intentionally separate from the production presets. No live worker or campaign configuration was changed by this experiment.

## Evidence limits

This run establishes a correct, reproducible optimization candidate and a measured same-GPU comparison. It does **not** establish 20 B/s, a production speedup over the current live geometry, or a full-length DP32 campaign qualification. The next useful optimization should reduce the table candidate's register/shared-memory cost or improve the existing 256-thread geometry; the table itself is not the limiting win on this run.
