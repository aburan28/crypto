# Engine suite: run `confirm_v1`

Problems: none

## Declared targets (ledger §17.2)

- **1 reach**: no data in this run — f4-f2: 0 cells, all decided False, agree True, matrix-f4: 0 cells, all decided False, agree True, matrix-f5: 0 cells, all decided False, agree True, inherited-f4: 0 cells, all decided False, agree True
- **2 gate**: no data in this run — 
- **3 question**: met — 3 crossing cell(s); crossbred-f2 at 26 unknowns (R, C1-holdout-26): 0.735 [0.693, 0.767] (k=8); crossbred-f2 at 28 unknowns (R, C2-main-28-30): 0.686 [0.661, 0.738] (k=8); crossbred-f2 at 28 unknowns (R, C3-holdout-28): 0.67 [0.651, 0.695] (k=8)
- **abandon**: no data in this run — gaining 16→26: {'buchberger-f2': None, 'f4-f2': None, 'matrix-f4': None, 'matrix-f5': None, 'inherited-f4': None, 'crossbred-f2': None}

## Stage: wall / the cell's reference (fes-f2), pooled over families, median [95% CI] (pairs)

| engine | 26 | 28 | 30 |
|:--|--:|--:|--:|
| buchberger-f2 | — | — | — |
| f4-f2 | 206 [191, 210] (k=8) | 114 [110, 121] (k=16) | 49.2 [31.7, 49.9] (k=5) |
| matrix-f4 | 7.78 [7.26, 8.22] (k=8) | 4.56 [4.4, 4.84] (k=16) | 2.75 [2.65, 2.82] (k=16) |
| matrix-f5 | 8.09 [7.68, 8.5] (k=8) | 4.81 [4.75, 4.87] (k=16) | 2.87 [2.76, 2.97] (k=16) |
| inherited-f4 | 4.5 [4.31, 4.64] (k=8) | 4.1 [3.95, 4.37] (k=16) | 1.24 [1.2, 1.25] (k=16) |
| crossbred-f2 | 0.735 [0.693, 0.767] (k=8) | 0.678 [0.666, 0.71] (k=16) | — |

## Extrapolation (marked: not a measurement)


## Stage cells

| part | E | n:n':m | vars | engine | finds | decided | budget | D_learn | D_reach | D_sr | ms | vs reference | vs baseline | agrees |
|:--|:--|:--|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|
| C1-holdout-26 | R | 25:13:2 | 26 | f4-f2 | all | 8/8 | 0 | 3.50 | 4.50 | 6 | 7526.35 | 206 [191, 210] (k=8) | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 299.26 | 7.78 [7.26, 8.22] (k=8) | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 6 | 309.24 | 8.09 [7.68, 8.43] (k=8) | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 173.37 | 4.5 [4.31, 4.64] (k=8) | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 28.04 | 0.735 [0.693, 0.767] (k=8) | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 38.28 | ref | — | yes |
| C1-holdout-26 | R | 25:13:2 | 26 | exhaustive | all | 8/8 | 0 | — | — | 6 | 4191.29 | 115 [90.7, 122] (k=8) | — | yes |
| C2-main-28-30 | K | 29:15:2 | 30 | f4-f2 | all | 3/8 | 5 | 3.00 | 4.00 | 6 | 27261.37 | 49.6 [49.2, 49.9] (k=3) | — | yes |
| C2-main-28-30 | K | 29:15:2 | 30 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 1702.61 | 2.72 [2.55, 2.91] (k=8) | — | yes |
| C2-main-28-30 | K | 29:15:2 | 30 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 6 | 1758.05 | 2.91 [2.6, 2.98] (k=8) | — | yes |
| C2-main-28-30 | K | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 769.79 | 1.24 [1.17, 1.26] (k=8) | — | yes |
| C2-main-28-30 | K | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 11.39 | — | — | — |
| C2-main-28-30 | K | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 623.13 | ref | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | f4-f2 | all | 8/8 | 0 | 3.88 | 4.88 | 6 | 17218.53 | 111 [106, 116] (k=8) | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 711.78 | 4.45 [4.38, 4.95] (k=8) | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 6 | 741.57 | 4.8 [4.42, 5.13] (k=8) | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 655.59 | 4.29 [3.89, 4.52] (k=8) | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 109.02 | 0.686 [0.661, 0.738] (k=8) | — | yes |
| C2-main-28-30 | R | 27:14:2 | 28 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 155.26 | ref | — | yes |
| C2-main-28-30 | R | 29:15:2 | 30 | f4-f2 | all | 2/8 | 6 | 3.00 | 4.00 | 6 | 26886.77 | 39.7 [31.7, 47.7] (k=2) | — | yes |
| C2-main-28-30 | R | 29:15:2 | 30 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 1703.77 | 2.77 [2.65, 2.83] (k=8) | — | yes |
| C2-main-28-30 | R | 29:15:2 | 30 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 6 | 1767.51 | 2.86 [2.78, 2.95] (k=8) | — | yes |
| C2-main-28-30 | R | 29:15:2 | 30 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 762.20 | 1.23 [1.17, 1.31] (k=8) | — | yes |
| C2-main-28-30 | R | 29:15:2 | 30 | crossbred-f2 | all | 0/8 | 8 | — | — | 6 | 11.77 | — | — | — |
| C2-main-28-30 | R | 29:15:2 | 30 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 619.33 | ref | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | f4-f2 | all | 8/8 | 0 | 4.00 | 5.00 | 6 | 17971.72 | 120 [111, 123] (k=8) | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | matrix-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 707.10 | 4.6 [4.48, 4.88] (k=8) | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | matrix-f5 | all | 8/8 | 0 | — | 3.00 | 6 | 745.60 | 4.81 [4.79, 4.85] (k=8) | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | inherited-f4 | all | 8/8 | 0 | — | 3.00 | 6 | 625.99 | 4.04 [3.95, 4.33] (k=8) | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | crossbred-f2 | all | 8/8 | 0 | — | 3.00 | 6 | 102.89 | 0.67 [0.651, 0.695] (k=8) | — | yes |
| C3-holdout-28 | R | 27:14:2 | 28 | fes-f2 | all | 8/8 | 0 | — | — | 6 | 151.91 | ref | — | yes |

## Whole method (S; per-relation GAE; paired over curve × planted target)

