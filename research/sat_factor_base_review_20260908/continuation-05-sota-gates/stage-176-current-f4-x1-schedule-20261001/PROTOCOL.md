# Stage 176 protocol: current F4 fixed-X1 schedule on one target

## Hypothesis

Stage 175 ran the current repository `BlockTables` F4 engine with all 242
rational fixed-X1 systems in one outer Rayon batch. Its 4.03 GB median RSS and
1.81x Stage 174 CPU regression suggest that simultaneous independent F4 calls
are contending for cache and memory while the modern engine also tries to
parallelise large matrices internally.

Reducing `PQ_F4_X1_BATCH` may lower total core-seconds and wall time without
changing a single equation or F4 algorithm. This stage tests only that schedule
on the same already-opened target. It is implementation engineering, not a new
target or cryptanalytic claim.

## Frozen source and target

- Binary SHA-256:
  `b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6`.
- Binary source commit:
  `c38626a5e2674dbfdff796474ea83af0c863917a`.
- Input: Stage 175 `input/manifest.json`, SHA-256
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Cell and target: the same `n=59, ell=9, m=3` blind true-negative.
- Expected equations BLAKE3:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- `RAYON_NUM_THREADS=12`; all BLAS/OpenMP thread controls remain one.

## Screen

One complete fresh process is run at each X1 batch size, in this fixed order:

```text
1, 4, 12, 24, 39, 64, 128, 256
```

Stage 175 supplies the three-repeat batch-512 baseline. The screen is a
diagnostic; no one-repeat result is a selected improvement. Every screen must
remain exhaustive UNSAT, visit all 512 masks, complete all 242 rational
systems, reproduce the equation fingerprint, and repeat the exact structural
and operation counts from Stage 175.

If at least one arm improves both wall and total core-seconds over the Stage
175 medians, select the arm minimizing total core-seconds (wall breaks ties)
and run three interleaved candidate/batch-512 pairs. Otherwise stop and reject
schedule tuning.

## Confirmation rule

The selected schedule passes only if the median paired candidate/control wall
ratio and total-core ratio are both below `0.97`, every correctness gate passes,
and RSS is reported. It is stronger if its absolute medians also beat the Stage
174 medians (`26.358217916989815` wall and `147.841771` core-seconds), but that
is reported separately and is not inferred from a one-target schedule win.

Every attempted process, including rejected screen arms and controls, is
charged in the stage total. Single-core time remains unmeasured in this
multi-worker schedule experiment and must remain null.
