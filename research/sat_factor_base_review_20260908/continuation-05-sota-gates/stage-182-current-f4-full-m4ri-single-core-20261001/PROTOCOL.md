# Stage 182 protocol: full-matrix M4RI single-core adjudication

## Reason for this adjudication

Stage 181 full M4RI reduced actual XORs by 32.89%, median paired CPU by 10.07%,
and median paired RSS by 8.62%, but uneven host descheduling made its median
paired wall ratio 1.0549 and blocked multi-worker promotion. A one-worker run
removes nested/outer scheduling ambiguity and supplies the current-engine
single-core row requested by the larger benchmark gate.

## Frozen inputs and modes

- Exact Stage 181 binary SHA-256:
  `fe17f002fedb2e4ef3037d100a4c541ffe0d2c7dcbfa6b2121a540804cc9050c`.
- Binary source commit:
  `850182efba8d8a9499e82aa6f11600941d3cb9d6`.
- Same opened `n=59, ell=9, m=3` target, algebraic factor base, and equation
  fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- `RAYON_NUM_THREADS=1`, `PQ_F4_X1_BATCH=1`, and every external thread
  control set to one.
- Current control: `F4_F2_FULL_M4RI=0`.
- Candidate: `F4_F2_FULL_M4RI=1`.

Run one current/full screen. Continue only if candidate wall, total core-seconds,
and performed XORs all improve. Confirmation uses three fresh interleaved pairs
in fixed order `current, full, full, current, current, full`.

## Correctness and accounting

Every process must authenticate the same source, visit all 512 masks, complete
all 242 rational systems, return exhaustive UNSAT, and reproduce the equation
fingerprint. The candidate must repeat Stage 181's 723 matrices, 438,923
blocks, logical and performed XOR counts. The control must repeat the current
logical/performed counts. Record valid single-core seconds, total core-seconds,
wall, RSS, matrix/table memory, and conflicts (`null` for F4).

Carry the Stage 181 exact build and validation costs into cumulative reporting;
charge every new screen and confirmation process.

## Decision and runtime policy

The candidate passes if the median paired full/current wall and total-core
ratios are both below `0.97`, every correctness gate passes, and every process
has non-null single-core seconds equal to total core-seconds.

If it passes, a follow-up code commit may select full M4RI automatically only
when Rayon has one worker, while retaining explicit `F4_F2_FULL_M4RI=0` and
`=1` controls and leaving the multi-worker default on current `BlockTables`.
That policy must pass the same unit/backend tests and a selected-mode replay.

This is still one public solver target. It cannot establish natural relation
yield, an unknown-scalar index-calculus recovery, rho crossover, independent
reproduction, novelty, or SOTA.
