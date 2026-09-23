# Engine suite: run `hybrids_v2`

Problems: none

## Declared targets (ledger §17.2)

- **1 reach**: no data in this run — f4-f2: 0 cells, all decided False, agree True, matrix-f4: 0 cells, all decided False, agree True, matrix-f5: 0 cells, all decided False, agree True, inherited-f4: 0 cells, all decided False, agree True
- **2 gate**: no data in this run — 
- **3 question**: not met — 0 crossing cell(s)
- **abandon**: no data in this run — gaining 16→26: {'buchberger-f2': None, 'f4-f2': None, 'matrix-f4': None, 'matrix-f5': None, 'inherited-f4': None, 'crossbred-f2': None}

## Stage: wall / the cell's reference (fes-f2-wide), pooled over families, median [95% CI] (pairs)

| engine | 32 |
|:--|--:|
| buchberger-f2 | — |
| f4-f2 | — |
| matrix-f4 | 5.4 [5, 5.82] (k=32) |
| matrix-f5 | 5.42 [5.17, 5.77] (k=32) |
| inherited-f4 | 2.09 [2.03, 2.11] (k=32) |
| crossbred-f2 | — |

## Extrapolation (marked: not a measurement)


## Stage cells

| part | E | n:n':m | vars | engine | finds | decided | budget | D_learn | D_reach | D_sr | ms | vs reference | vs baseline | agrees |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|
| X1-main | K | 31:16:2 | 32 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 1416.81 | 2.11 [2.09, 2.17] (k=8) | — | yes |
| X1-main | K | 31:16:2 | 32 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3751.85 | 5.73 [4.95, 6.37] (k=8) | — | yes |
| X1-main | K | 31:16:2 | 32 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 3731.85 | 5.75 [5.17, 6.04] (k=8) | — | yes |
| X1-main | K | 31:16:2 | 32 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 15.21 | — | — | — |
| X1-main | K | 31:16:2 | 32 | fes-f2 | all | 8/8 | 0 | — | — | 7 | 2413.93 | 3.63 [3.52, 3.77] (k=8) | — | yes |
| X1-main | K | 31:16:2 | 32 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 663.32 | ref | — | yes |
| X1-main | R | 31:16:2 | 32 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 1374.48 | 2.03 [1.97, 2.09] (k=8) | — | yes |
| X1-main | R | 31:16:2 | 32 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3527.99 | 5.12 [4.65, 5.83] (k=8) | — | yes |
| X1-main | R | 31:16:2 | 32 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 3623.74 | 5.36 [4.73, 6.03] (k=8) | — | yes |
| X1-main | R | 31:16:2 | 32 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 16.18 | — | — | — |
| X1-main | R | 31:16:2 | 32 | fes-f2 | all | 8/8 | 0 | — | — | 7 | 2462.19 | 3.67 [3.5, 3.8] (k=8) | — | yes |
| X1-main | R | 31:16:2 | 32 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 673.45 | ref | — | yes |
| X1-main | R | 33:17:2 | 34 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3206.43 | ref | — | yes |
| X1-main | R | 33:17:2 | 34 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 6863.92 | 2.15 [1.97, 2.29] (k=8) | — | yes |
| X1-main | R | 33:17:2 | 34 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 6967.83 | 2.16 [2.03, 2.3] (k=8) | — | yes |
| X1-main | R | 33:17:2 | 34 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 19.94 | — | — | — |
| X2-holdout | K | 31:16:2 | 32 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 1354.49 | 2.09 [1.93, 2.17] (k=8) | — | yes |
| X2-holdout | K | 31:16:2 | 32 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3737.70 | 5.84 [5.44, 5.95] (k=8) | — | yes |
| X2-holdout | K | 31:16:2 | 32 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 3731.04 | 5.69 [5.29, 6.03] (k=8) | — | yes |
| X2-holdout | K | 31:16:2 | 32 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 14.58 | — | — | — |
| X2-holdout | K | 31:16:2 | 32 | fes-f2 | all | 8/8 | 0 | — | — | 7 | 2574.11 | 3.9 [3.8, 4.09] (k=8) | — | yes |
| X2-holdout | K | 31:16:2 | 32 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 653.61 | ref | — | yes |
| X2-holdout | R | 31:16:2 | 32 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 1371.97 | 2.04 [1.88, 2.11] (k=8) | — | yes |
| X2-holdout | R | 31:16:2 | 32 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 3332.94 | 4.93 [4.39, 5.36] (k=8) | — | yes |
| X2-holdout | R | 31:16:2 | 32 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 3398.10 | 4.88 [4.51, 5.63] (k=8) | — | yes |
| X2-holdout | R | 31:16:2 | 32 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 14.84 | — | — | — |
| X2-holdout | R | 31:16:2 | 32 | fes-f2 | all | 8/8 | 0 | — | — | 7 | 2634.71 | 3.88 [3.72, 3.96] (k=8) | — | yes |
| X2-holdout | R | 31:16:2 | 32 | fes-f2-wide | all | 8/8 | 0 | — | — | 7 | 680.28 | ref | — | yes |
| X2-holdout | R | 33:17:2 | 34 | inherited-f4 | all | 6/8 | 2 | — | 3.00 | 7 | 3117.09 | ref | — | yes |
| X2-holdout | R | 33:17:2 | 34 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 7 | 6376.07 | 2.04 [1.9, 2.23] (k=6) | — | yes |
| X2-holdout | R | 33:17:2 | 34 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 7 | 6452.58 | 2.05 [1.94, 2.25] (k=6) | — | yes |
| X2-holdout | R | 33:17:2 | 34 | crossbred-f2 | all | 0/8 | 8 | — | — | 7 | 20.23 | — | — | — |

## Whole method (S; per-relation GAE; paired over curve × planted target)

