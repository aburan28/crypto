# Engine suite: run `hybrids_v3`

Problems: none

## Declared targets (ledger §17.2)

- **1 reach**: no data in this run — f4-f2: 0 cells, all decided False, agree True, matrix-f4: 0 cells, all decided False, agree True, matrix-f5: 0 cells, all decided False, agree True, inherited-f4: 0 cells, all decided False, agree True
- **2 gate**: no data in this run — 
- **3 question**: not met — 0 crossing cell(s)
- **abandon**: no data in this run — gaining 16→26: {'buchberger-f2': None, 'f4-f2': None, 'matrix-f4': None, 'matrix-f5': None, 'inherited-f4': None, 'crossbred-f2': None}

## Stage: wall / the cell's reference (fes-f2-wide), pooled over families, median [95% CI] (pairs)

| engine | 34 |
|:--|--:|
| buchberger-f2 | — |
| f4-f2 | — |
| matrix-f4 | 2.45 [2.37, 2.56] (k=16) |
| matrix-f5 | 2.48 [2.41, 2.6] (k=16) |
| inherited-f4 | 1.2 [1.2, 1.23] (k=14) |
| crossbred-f2 | — |

## Extrapolation (marked: not a measurement)


## Stage cells

| part | E | n:n':m | vars | engine | finds | decided | budget | D_learn | D_reach | D_sr | ms | vs reference | vs baseline | agrees |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|
| Y1-main-34 | R | 33:17:2 | 34 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3192.49 | 1.21 [1.2, 1.24] (k=8) | — | yes |
| Y1-main-34 | R | 33:17:2 | 34 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 6793.36 | 2.53 [2.48, 2.8] (k=8) | — | yes |
| Y1-main-34 | R | 33:17:2 | 34 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 6806.85 | 2.55 [2.47, 2.73] (k=8) | — | yes |
| Y1-main-34 | R | 33:17:2 | 34 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 20.43 | — | — | — |
| Y1-main-34 | R | 33:17:2 | 34 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 2628.85 | ref | — | yes |
| Y2-holdout-34 | R | 33:17:2 | 34 | inherited-f4 | all | 6/8 | 2 | — | 3.00 | 7 | 3198.79 | 1.2 [1.16, 1.22] (k=6) | — | yes |
| Y2-holdout-34 | R | 33:17:2 | 34 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 6378.45 | 2.35 [2.32, 2.4] (k=8) | — | yes |
| Y2-holdout-34 | R | 33:17:2 | 34 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 6544.80 | 2.4 [2.36, 2.48] (k=8) | — | yes |
| Y2-holdout-34 | R | 33:17:2 | 34 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 20.51 | — | — | — |
| Y2-holdout-34 | R | 33:17:2 | 34 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 2668.05 | ref | — | yes |

## Whole method (S; per-relation GAE; paired over curve × planted target)

