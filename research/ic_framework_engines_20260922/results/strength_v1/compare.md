# Engine suite: run `strength_v1`

Problems: none

## Declared targets (ledger §17.2)

- **1 reach**: no data in this run — f4-f2: 0 cells, all decided False, agree True, matrix-f4: 0 cells, all decided False, agree True, matrix-f5: 0 cells, all decided False, agree True, inherited-f4: 0 cells, all decided False, agree True
- **2 gate**: no data in this run — 
- **3 question**: not met — 0 crossing cell(s)
- **abandon**: no data in this run — gaining 16→26: {'buchberger-f2': None, 'f4-f2': None, 'matrix-f4': None, 'matrix-f5': None, 'inherited-f4': None, 'crossbred-f2': None}

## Stage: wall / the cell's reference (fes-f2-wide), pooled over families, median [95% CI] (pairs)

| engine | 24 | 26 | 28 | 30 |
|:--|--:|--:|--:|--:|
| buchberger-f2 | — | — | — | — |
| f4-f2 | — | — | — | — |
| matrix-f4 | — | — | — | — |
| matrix-f5 | — | — | — | — |
| inherited-f4 | 22.6 [21.7, 23.2] (k=32) | 15.4 [14.7, 15.8] (k=16) | 14.8 [14.5, 15.8] (k=16) | 4.55 [4.47, 4.62] (k=32) |
| crossbred-f2 | 3 [2.81, 3.1] (k=32) | 2.54 [2.46, 2.59] (k=16) | 2.45 [2.36, 2.55] (k=16) | — |

## Extrapolation (marked: not a measurement)

- inherited-f4: ×0.62 per two unknowns over [24, 26, 28, 30]; crossover at 37.4 unknowns

## Stage cells

| part | E | n:n':m | vars | engine | finds | decided | budget | D_learn | D_reach | D_sr | ms | vs reference | vs baseline | agrees |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|
| S1-main | K | 23:12:2 | 24 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 10.18 | 3.08 [2.72, 3.25] (k=8) | — | yes |
| S1-main | K | 23:12:2 | 24 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 10.67 | 3.14 [2.95, 3.19] (k=8) | — | yes |
| S1-main | K | 23:12:2 | 24 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 3.54 | ref | — | yes |
| S1-main | K | 23:12:2 | 24 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 81.80 | 22.3 [20.7, 23.4] (k=8) | — | yes |
| S1-main | K | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 13.91 | — | — | — |
| S1-main | K | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 680.10 | 4.06 [3.83, 4.12] (k=8) | — | yes |
| S1-main | K | 29:15:2 | 30 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 169.73 | ref | — | yes |
| S1-main | K | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 779.35 | 4.6 [4.47, 4.74] (k=8) | — | yes |
| S1-main | R | 23:12:2 | 24 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 10.32 | 3.04 [2.75, 3.15] (k=8) | — | yes |
| S1-main | R | 23:12:2 | 24 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 10.66 | 3.09 [3.04, 3.25] (k=8) | — | yes |
| S1-main | R | 23:12:2 | 24 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 3.43 | ref | — | yes |
| S1-main | R | 23:12:2 | 24 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 74.91 | 21.8 [20.9, 23.1] (k=8) | — | yes |
| S1-main | R | 25:13:2 | 26 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 29.54 | 2.57 [2.53, 2.65] (k=8) | — | yes |
| S1-main | R | 25:13:2 | 26 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 40.67 | 3.55 [3.5, 3.63] (k=8) | — | yes |
| S1-main | R | 25:13:2 | 26 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 11.43 | ref | — | yes |
| S1-main | R | 25:13:2 | 26 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 178.63 | 15.6 [15, 16.1] (k=8) | — | yes |
| S1-main | R | 27:14:2 | 28 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 109.28 | 2.47 [2.27, 2.61] (k=8) | — | yes |
| S1-main | R | 27:14:2 | 28 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 166.24 | 3.87 [3.75, 3.93] (k=8) | — | yes |
| S1-main | R | 27:14:2 | 28 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 42.94 | ref | — | yes |
| S1-main | R | 27:14:2 | 28 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 675.27 | 14.8 [14.1, 17.2] (k=8) | — | yes |
| S1-main | R | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 11.72 | — | — | — |
| S1-main | R | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 661.04 | 3.87 [3.8, 4.07] (k=8) | — | yes |
| S1-main | R | 29:15:2 | 30 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 168.84 | ref | — | yes |
| S1-main | R | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 757.29 | 4.47 [4.28, 4.72] (k=8) | — | yes |
| S2-holdout | K | 23:12:2 | 24 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 10.00 | 2.87 [2.71, 3.15] (k=8) | — | yes |
| S2-holdout | K | 23:12:2 | 24 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 10.50 | 3.09 [2.96, 3.18] (k=8) | — | yes |
| S2-holdout | K | 23:12:2 | 24 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 3.39 | ref | — | yes |
| S2-holdout | K | 23:12:2 | 24 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 77.94 | 23.4 [21.7, 23.9] (k=8) | — | yes |
| S2-holdout | K | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 12.06 | — | — | — |
| S2-holdout | K | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 653.85 | 3.92 [3.79, 3.94] (k=8) | — | yes |
| S2-holdout | K | 29:15:2 | 30 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 169.06 | ref | — | yes |
| S2-holdout | K | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 763.86 | 4.57 [4.35, 4.75] (k=8) | — | yes |
| S2-holdout | R | 23:12:2 | 24 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 10.25 | 2.92 [2.68, 3.15] (k=8) | — | yes |
| S2-holdout | R | 23:12:2 | 24 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 10.57 | 3.02 [2.95, 3.13] (k=8) | — | yes |
| S2-holdout | R | 23:12:2 | 24 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 3.50 | ref | — | yes |
| S2-holdout | R | 23:12:2 | 24 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 78.96 | 23.2 [21.3, 24] (k=8) | — | yes |
| S2-holdout | R | 25:13:2 | 26 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 28.28 | 2.44 [2.33, 2.58] (k=8) | — | yes |
| S2-holdout | R | 25:13:2 | 26 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 42.18 | 3.63 [3.58, 3.76] (k=8) | — | yes |
| S2-holdout | R | 25:13:2 | 26 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 11.56 | ref | — | yes |
| S2-holdout | R | 25:13:2 | 26 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 174.38 | 14.9 [14.5, 15.8] (k=8) | — | yes |
| S2-holdout | R | 27:14:2 | 28 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 108.84 | 2.44 [2.36, 2.55] (k=8) | — | yes |
| S2-holdout | R | 27:14:2 | 28 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 171.66 | 3.86 [3.77, 3.91] (k=8) | — | yes |
| S2-holdout | R | 27:14:2 | 28 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 44.54 | ref | — | yes |
| S2-holdout | R | 27:14:2 | 28 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 667.24 | 14.7 [14.3, 15.8] (k=8) | — | yes |
| S2-holdout | R | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 12.60 | — | — | — |
| S2-holdout | R | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 668.39 | 3.95 [3.91, 4.12] (k=8) | — | yes |
| S2-holdout | R | 29:15:2 | 30 | fes-f2-wide | all | 8/8 | 0 | — | — | 6 | 167.70 | ref | — | yes |
| S2-holdout | R | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 759.69 | 4.54 [4.42, 4.64] (k=8) | — | yes |

## Whole method (S; per-relation GAE; paired over curve × planted target)

