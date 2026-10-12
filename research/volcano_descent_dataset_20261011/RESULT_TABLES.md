### Campaign

|   ell |   volcanoes |   max_depth |   vertices |   stored |   max_p_bits |
|------:|------------:|------------:|-----------:|---------:|-------------:|
|     2 |          36 |          15 |     391133 |    92541 |           68 |
|     3 |          22 |           9 |      73991 |    35879 |           62 |
|     5 |          20 |           5 |      59374 |    35974 |           63 |
|     7 |          12 |           4 |      26412 |    26412 |           62 |
|    11 |           6 |           3 |       5326 |     5326 |           60 |
|    13 |           6 |           3 |       8790 |     8790 |           62 |

volcanoes: 102; vertices in the exact volcanoes: 565026; rows stored: 204922; deepest volcano: 15 levels.

### Verification (fail-closed)

| check              |   volcanoes passing |   of |
|:-------------------|--------------------:|-----:|
| V1_level_counts    |                 102 |  102 |
| V2_degrees         |                 102 |  102 |
| V3_group_structure |                 102 |  102 |
| V4_point_count     |                 102 |  102 |
| V5_blind_descent   |                 102 |  102 |

vertex rows whose measured group structure equals the O_d/(pi-1) prediction: 204922 of 204922.

blind descents (start vertex only j, p, t, ell): 2448 of 2448 reach a floor curve with the correct level and path length; median 0.0001 s, max 0.059 s.

volcanoes with a non-cyclic group somewhere (level-dependent group structure): 80 of 102.

### ECDLP E1

4864 solves, 0 incorrect, 6 volcanoes, 152 vertices, largest-prime sizes 28-32 bits.

| mode   | layer   |   solves |   vertices |   mean ops_rho / sqrt(pi r/2) |   ci_lo |   ci_hi |
|:-------|:--------|---------:|-----------:|------------------------------:|--------:|--------:|
| neg    | crater  |      144 |          9 |                         0.624 |   0.566 |   0.684 |
| neg    | floor   |      192 |         12 |                         0.712 |   0.662 |   0.764 |
| neg    | middle  |     2096 |        131 |                         0.675 |   0.66  |   0.69  |
| plain  | crater  |      144 |          9 |                         0.859 |   0.787 |   0.934 |
| plain  | floor   |      192 |         12 |                         0.956 |   0.887 |   1.03  |
| plain  | middle  |     2096 |        131 |                         0.914 |   0.895 |   0.934 |

| layer   | ratio of mean ops   |   ratio |   ci_lo |   ci_hi |
|:--------|:--------------------|--------:|--------:|--------:|
| crater  | plain / neg         |    1.38 |    1.21 |    1.56 |
| floor   | plain / neg         |    1.34 |    1.21 |    1.49 |
| middle  | plain / neg         |    1.35 |    1.31 |    1.4  |

| volcano                       | mode   |   m |   slope of norm. cost per level |   mean norm. cost |
|:------------------------------|:-------|----:|--------------------------------:|------------------:|
| v_l2_D11_c1_k12_m6_s20261011  | neg    |   6 |                        0.00087  |             0.766 |
| v_l2_D11_c1_k12_m6_s20261011  | plain  |   6 |                        0.00348  |             1.08  |
| v_l2_D19_c1_k12_m6_s20261011  | neg    |   6 |                        0.00669  |             0.653 |
| v_l2_D19_c1_k12_m6_s20261011  | plain  |   6 |                        0.0038   |             0.923 |
| v_l2_D23_c1_k12_m12_s20261011 | neg    |  12 |                       -0.000892 |             0.685 |
| v_l2_D23_c1_k12_m12_s20261011 | plain  |  12 |                        0.000561 |             0.885 |
| v_l2_D24_c1_k12_m0_s20261011  | neg    |   0 |                        0.0109   |             0.645 |
| v_l2_D24_c1_k12_m0_s20261011  | plain  |   0 |                       -0.00343  |             0.853 |
| v_l2_D24_c1_k12_m12_s20261011 | neg    |  12 |                       -0.000284 |             0.616 |
| v_l2_D24_c1_k12_m12_s20261011 | plain  |  12 |                        0.00456  |             0.809 |
| v_l2_D7_c1_k12_m0_s20261011   | neg    |   0 |                       -0.00185  |             0.693 |
| v_l2_D7_c1_k12_m0_s20261011   | plain  |   0 |                        0.00917  |             0.953 |

|   m |   mean ops_total (PH-aware) |   mean ops_rho |   mean sqrt(ord G) (no PH) |
|----:|----------------------------:|---------------:|---------------------------:|
|   0 |                    3.23e+04 |       2.24e+04 |                   5.7e+05  |
|   6 |                    6.77e+04 |       5.76e+04 |                   1.81e+07 |
|  12 |                    6.18e+04 |       5.14e+04 |                   8.21e+08 |

### ECDLP E2

8352 solves, 0 incorrect, 4 volcanoes, 96 vertices, largest-prime sizes 28-29 bits.

| mode   | layer   |   solves |   vertices |   mean ops_rho / sqrt(pi r/2) |   ci_lo |   ci_hi |
|:-------|:--------|---------:|-----------:|------------------------------:|--------:|--------:|
| aut    | crater  |      576 |          4 |                         0.461 |   0.437 |   0.487 |
| aut    | floor   |      192 |          8 |                         0.685 |   0.636 |   0.74  |
| aut    | middle  |     2016 |         84 |                         0.644 |   0.628 |   0.659 |
| neg    | crater  |      576 |          4 |                         0.64  |   0.612 |   0.67  |
| neg    | floor   |      192 |          8 |                         0.6   |   0.556 |   0.646 |
| neg    | middle  |     2016 |         84 |                         0.657 |   0.642 |   0.672 |
| plain  | crater  |      576 |          4 |                         0.884 |   0.846 |   0.922 |
| plain  | floor   |      192 |          8 |                         0.884 |   0.818 |   0.946 |
| plain  | middle  |     2016 |         84 |                         0.876 |   0.856 |   0.896 |

| layer   | ratio of mean ops   |   ratio |   ci_lo |   ci_hi |
|:--------|:--------------------|--------:|--------:|--------:|
| crater  | plain / neg         |   1.38  |   1.3   |   1.47  |
| crater  | neg / aut           |   1.39  |   1.3   |   1.49  |
| floor   | plain / neg         |   1.47  |   1.32  |   1.64  |
| floor   | neg / aut           |   0.875 |   0.79  |   0.969 |
| middle  | plain / neg         |   1.33  |   1.29  |   1.38  |
| middle  | neg / aut           |   1.02  |   0.986 |   1.05  |

Per volcano and layer (mean of ops_rho / sqrt(pi r / 2); ratios of means, bootstrap 95% intervals):

| volcano                      | layer   |   solves |   vertices |   plain |   neg |   aut | plain/neg         | neg/aut           | plain/aut         |
|:-----------------------------|:--------|---------:|-----------:|--------:|------:|------:|:------------------|:------------------|:------------------|
| v_l2_D3_c1_k12_m0_s20261011  | crater  |      432 |          1 |   0.937 | 0.705 | 0.493 | 1.33 [1.18, 1.49] | 1.43 [1.21, 1.68] | 1.90 [1.64, 2.22] |
| v_l2_D3_c1_k12_m0_s20261011  | middle  |     1512 |         21 |   0.910 | 0.661 | 0.645 | 1.38 [1.29, 1.47] | 1.02 [0.96, 1.09] | 1.41 [1.32, 1.50] |
| v_l2_D3_c1_k12_m0_s20261011  | floor   |      144 |          2 |   0.916 | 0.573 | 0.663 | 1.60 [1.26, 2.03] | 0.86 [0.68, 1.12] | 1.38 [1.09, 1.75] |
| v_l2_D3_c1_k12_m6_s20261011  | crater  |      432 |          1 |   0.785 | 0.559 | 0.363 | 1.41 [1.26, 1.59] | 1.54 [1.37, 1.74] | 2.16 [1.93, 2.42] |
| v_l2_D3_c1_k12_m6_s20261011  | middle  |     1512 |         21 |   0.791 | 0.583 | 0.585 | 1.36 [1.28, 1.44] | 1.00 [0.94, 1.06] | 1.35 [1.27, 1.44] |
| v_l2_D3_c1_k12_m6_s20261011  | floor   |      144 |          2 |   0.851 | 0.584 | 0.617 | 1.46 [1.17, 1.86] | 0.95 [0.75, 1.15] | 1.38 [1.12, 1.71] |
| v_l2_D4_c1_k12_m0_s20261011  | crater  |      432 |          1 |   0.783 | 0.623 | 0.428 | 1.26 [1.11, 1.41] | 1.46 [1.30, 1.63] | 1.83 [1.62, 2.06] |
| v_l2_D4_c1_k12_m0_s20261011  | middle  |     1512 |         21 |   0.812 | 0.594 | 0.587 | 1.37 [1.28, 1.46] | 1.01 [0.95, 1.08] | 1.38 [1.30, 1.48] |
| v_l2_D4_c1_k12_m0_s20261011  | floor   |      144 |          2 |   0.714 | 0.535 | 0.669 | 1.33 [1.11, 1.62] | 0.80 [0.68, 0.95] | 1.07 [0.87, 1.30] |
| v_l2_D4_c1_k12_m12_s20261011 | crater  |      432 |          1 |   1.029 | 0.674 | 0.560 | 1.53 [1.34, 1.73] | 1.20 [1.04, 1.38] | 1.84 [1.62, 2.08] |
| v_l2_D4_c1_k12_m12_s20261011 | middle  |     1512 |         21 |   0.990 | 0.788 | 0.759 | 1.26 [1.18, 1.34] | 1.04 [0.98, 1.11] | 1.30 [1.23, 1.40] |
| v_l2_D4_c1_k12_m12_s20261011 | floor   |      144 |          2 |   1.054 | 0.707 | 0.792 | 1.49 [1.26, 1.75] | 0.89 [0.75, 1.05] | 1.33 [1.12, 1.58] |

Pooled by D_K (ideal gains: j = 0, mu_6: neg/aut 1.73, plain/aut 2.45; j = 1728, mu_4: neg/aut 1.41, plain/aut 2.00; elsewhere neg/aut 1.00):

|   D_K | layer      |   solves per mode | neg / aut         | plain / aut       |
|------:|:-----------|------------------:|:------------------|:------------------|
|    -4 | crater     |               288 | 1.31 [1.20, 1.44] | 1.84 [1.68, 2.02] |
|    -4 | not crater |              1104 | 1.01 [0.97, 1.05] | 1.33 [1.27, 1.38] |
|    -3 | crater     |               288 | 1.48 [1.33, 1.64] | 2.01 [1.82, 2.24] |
|    -3 | not crater |              1104 | 1.00 [0.96, 1.05] | 1.38 [1.33, 1.44] |

### ML baseline: predicting the level (distance to the crater) from observables

| volcanoes                                     | feature set                                 |   train rows |   test rows |   balanced accuracy |   constant-prediction baseline |
|:----------------------------------------------|:--------------------------------------------|-------------:|------------:|--------------------:|-------------------------------:|
| every curve in the volcano has a cyclic group | raw curve only (j,a,b stats)                |        23865 |       13310 |               0.111 |                         0.111  |
| every curve in the volcano has a cyclic group | + degrees for other ell (no own-ell degree) |        23865 |       13310 |               0.11  |                         0.111  |
| every curve in the volcano has a cyclic group | + group structure (n1)                      |        23865 |       13310 |               0.11  |                         0.111  |
| every curve in the volcano has a cyclic group | + own-ell degree (graph-derived)            |        23865 |       13310 |               0.279 |                         0.111  |
| group structure varies with level             | raw curve only (j,a,b stats)                |       120323 |       47424 |               0.13  |                         0.0909 |
| group structure varies with level             | + degrees for other ell (no own-ell degree) |       120323 |       47424 |               0.146 |                         0.0909 |
| group structure varies with level             | + group structure (n1)                      |       120323 |       47424 |               0.576 |                         0.0909 |
| group structure varies with level             | + own-ell degree (graph-derived)            |       120323 |       47424 |               0.54  |                         0.0909 |

Target: round(10 * height / k), 11 classes (0 = crater side ... 10 = at the floor is height 0), held-out volcanoes (grouped split), so no volcano contributes to both sides.
