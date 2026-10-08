# Full-cost index-calculus iterations against rho

Rho parity is not established: the frozen all-cases gates are not met. 5 engineering hypotheses were implemented and compared against their predecessors. The strongest structural change compresses the rational support of the degree-eleven ordinary seed from three coordinates to one while retaining every factor-base point. All claims below concern the measured small fields, with degree-131 solving and asymptotic parity unestablished. On fresh confirmation targets, the candidate uses 1.06–1.59 times rho's instructions. Its full instruction cost is 72.6–94.0% below the original chart implementation.

The public algebraic visitor still returns every original S3 root. Group-witness search may use the smaller rational support. Point lifts, Frobenius coefficients and fixed-base point powers are reused; small projected systems use packed checks and fixed-capacity row elimination. All setup is charged, direct-relation shortcuts and caches remain disabled.

| Case | n | k | Seed dimension | Family | Subgroup order N | Role |
|---|---|---|---|---|---|---|
| 0 | 7 | 1 | 2 | power_span | 71 | useful |
| 1 | 9 | 1 | 3 | power_span | 127 | zero-column control |
| 2 | 9 | 1 | 3 | scaled_subfield | 127 | useful |
| 3 | 9 | 1 | 3 | subfield_control | 127 | zero-column control |
| 4 | 11 | 1 | 3 | power_span | 991 | useful |
| 5 | 11 | 1 | 3 | random | 991 | useful |
| 6 | 13 | 1 | 3 | power_span | 2003 | useful |
| 7 | 15 | 3 | 3 | power_span | 661 | useful |
| 8 | 15 | 3 | 3 | scaled_subfield | 661 | useful |
| 9 | 15 | 3 | 5 | subfield_control | 661 | useful |

## Correctness and archived workload

Every completed full-DLP oracle input is checked against independent pair truth, every reported logarithm verifies `[k]G=Q`, and matched coefficient/target prefixes agree. Incomplete zero-column controls remain in every suite; they never count as speedups.

| Comparison | Processes | Matched native targets | Complete algebraic root sets | Oracle checks |
|---|---|---|---|---|
| confirmation | 2004 | 912 | 1920 | 88956 |
| iteration-1 | 2004 | 912 | 1920 | 89682 |
| iteration-2 | 2004 | 912 | 1920 | 89682 |
| iteration-3 | 2004 | 912 | 1920 | 89682 |
| iteration-4 | 2004 | 912 | 1920 | 89682 |
| iteration-5 | 2004 | 912 | 1920 | 89682 |

## Twenty-percent engineering targets

These lists identify useful cases meeting a candidate/predecessor cost ratio at most 0.8, with the paired 95% upper limit below 1. All other useful cases miss that numeric target or its confidence gate. The original four new inputs are reused validation inputs after iteration 1. These engineering gates are separate from the all-cases rho parity gate, which remains unmet.

| Iteration | Cases meeting instruction target | Cases meeting native runtime target |
|---|---|---|
| iteration-1 | 0, 2, 4, 5, 6, 7, 8, 9 | 0, 2, 4, 5, 6, 7, 8, 9 |
| iteration-2 | 4, 5, 6 | 4, 5, 6 |
| iteration-3 | 4, 5, 9 | 4, 5 |
| iteration-4 | 0, 9 | 0, 8, 9 |
| iteration-5 | 6 | 5 |

## Complete instruction accounting

The table uses one machine-specific instruction model. `S_I = mean cold Ir / (measured Ir per addition × sqrt(N))`. The conversion uses 4,096 public random point additions minus a loop control on the same curve; all variants share the candidate calibration. `Ir/rho` is the ratio of arithmetic mean cold instruction counts on identical targets. This is distinct from the historical algebraic-operation S; it does not establish an attack exponent or a calibrated mathematical floor. Kernel execution and physical memory latency are outside this model. Empty controls have no verified-cost result.

| Cohort | Case | Variant | Verified | Mean cold million Ir | S_I | Ir/rho | Original/variant Ir | Floor ratio | Class |
|---|---|---|---|---|---|---|---|---|---|
| validation | 0 | original charts | 8/8 | 6.580 | 901.40 | 3.883 | 1.000 | null | engineering |
| validation | 0 | iteration-1 charts | 8/8 | 3.118 | 427.15 | 1.840 | 2.110 | null | engineering |
| validation | 0 | iteration-2 charts | 8/8 | 2.886 | 395.32 | 1.703 | 2.280 | null | engineering |
| validation | 0 | iteration-3 charts | 8/8 | 2.388 | 327.13 | 1.409 | 2.755 | null | engineering |
| validation | 0 | iteration-4 charts | 8/8 | 1.917 | 262.54 | 1.131 | 3.433 | null | engineering |
| validation | 0 | iteration-5 charts | 8/8 | 1.878 | 257.30 | 1.108 | 3.503 | null | engineering |
| validation | 0 | original pair table | 8/8 | 6.257 | 857.19 | 3.692 | 1.000 | null | engineering |
| validation | 0 | candidate pair table | 8/8 | 2.127 | 291.32 | 1.255 | 2.942 | null | engineering |
| validation | 0 | rho | 8/8 | 1.695 | 232.16 | 1.000 | 1.080 | null | engineering |
| validation | 1 | original charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | iteration-1 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | iteration-2 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | iteration-3 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | iteration-4 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | iteration-5 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | original pair table | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | candidate pair table | 0/8 | null | null | null | null | null | engineering |
| validation | 1 | rho | 8/8 | 2.483 | 202.24 | 1.000 | 1.226 | null | engineering |
| validation | 2 | original charts | 8/8 | 14.763 | 1202.32 | 5.945 | 1.000 | null | engineering |
| validation | 2 | iteration-1 charts | 8/8 | 5.377 | 437.89 | 2.165 | 2.746 | null | engineering |
| validation | 2 | iteration-2 charts | 8/8 | 5.053 | 411.51 | 2.035 | 2.922 | null | engineering |
| validation | 2 | iteration-3 charts | 8/8 | 4.068 | 331.30 | 1.638 | 3.629 | null | engineering |
| validation | 2 | iteration-4 charts | 8/8 | 3.396 | 276.56 | 1.367 | 4.347 | null | engineering |
| validation | 2 | iteration-5 charts | 8/8 | 3.297 | 268.50 | 1.328 | 4.478 | null | engineering |
| validation | 2 | original pair table | 8/8 | 14.375 | 1170.76 | 5.789 | 1.000 | null | engineering |
| validation | 2 | candidate pair table | 8/8 | 4.166 | 339.29 | 1.678 | 3.451 | null | engineering |
| validation | 2 | rho | 8/8 | 2.483 | 202.24 | 1.000 | null | null | engineering |
| validation | 3 | original charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | iteration-1 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | iteration-2 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | iteration-3 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | iteration-4 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | iteration-5 charts | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | original pair table | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | candidate pair table | 0/8 | null | null | null | null | null | engineering |
| validation | 3 | rho | 8/8 | 2.483 | 202.24 | 1.000 | null | null | engineering |
| validation | 4 | original charts | 8/8 | 94.593 | 1991.21 | 27.347 | 1.000 | null | engineering |
| validation | 4 | iteration-1 charts | 8/8 | 71.893 | 1513.36 | 20.785 | 1.316 | null | engineering |
| validation | 4 | iteration-2 charts | 8/8 | 12.256 | 258.00 | 3.543 | 7.718 | null | engineering |
| validation | 4 | iteration-3 charts | 8/8 | 6.677 | 140.54 | 1.930 | 14.168 | null | engineering |
| validation | 4 | iteration-4 charts | 8/8 | 5.621 | 118.33 | 1.625 | 16.828 | null | engineering |
| validation | 4 | iteration-5 charts | 8/8 | 4.774 | 100.50 | 1.380 | 19.812 | null | engineering |
| validation | 4 | original pair table | 8/8 | 26.481 | 557.43 | 7.656 | 1.000 | null | engineering |
| validation | 4 | candidate pair table | 8/8 | 3.659 | 77.02 | 1.058 | 7.238 | null | engineering |
| validation | 4 | rho | 8/8 | 3.459 | 72.81 | 1.000 | 1.122 | null | engineering |
| validation | 5 | original charts | 8/8 | 48.298 | 1016.68 | 13.963 | 1.000 | null | engineering |
| validation | 5 | iteration-1 charts | 8/8 | 31.033 | 653.25 | 8.972 | 1.556 | null | engineering |
| validation | 5 | iteration-2 charts | 8/8 | 8.445 | 177.76 | 2.441 | 5.719 | null | engineering |
| validation | 5 | iteration-3 charts | 8/8 | 6.474 | 136.28 | 1.872 | 7.460 | null | engineering |
| validation | 5 | iteration-4 charts | 8/8 | 5.506 | 115.91 | 1.592 | 8.771 | null | engineering |
| validation | 5 | iteration-5 charts | 8/8 | 4.942 | 104.02 | 1.429 | 9.774 | null | engineering |
| validation | 5 | original pair table | 8/8 | 24.649 | 518.86 | 7.126 | 1.000 | null | engineering |
| validation | 5 | candidate pair table | 8/8 | 5.793 | 121.95 | 1.675 | 4.255 | null | engineering |
| validation | 5 | rho | 8/8 | 3.459 | 72.81 | 1.000 | null | null | engineering |
| validation | 6 | original charts | 8/8 | 58.299 | 739.63 | 14.883 | 1.000 | null | engineering |
| validation | 6 | iteration-1 charts | 8/8 | 32.104 | 407.29 | 8.196 | 1.816 | null | engineering |
| validation | 6 | iteration-2 charts | 8/8 | 12.548 | 159.20 | 3.203 | 4.646 | null | engineering |
| validation | 6 | iteration-3 charts | 8/8 | 10.937 | 138.75 | 2.792 | 5.331 | null | engineering |
| validation | 6 | iteration-4 charts | 8/8 | 9.621 | 122.06 | 2.456 | 6.060 | null | engineering |
| validation | 6 | iteration-5 charts | 8/8 | 7.339 | 93.10 | 1.873 | 7.944 | null | engineering |
| validation | 6 | original pair table | 8/8 | 34.813 | 441.66 | 8.887 | 1.000 | null | engineering |
| validation | 6 | candidate pair table | 8/8 | 7.107 | 90.17 | 1.814 | 4.898 | null | engineering |
| validation | 6 | rho | 8/8 | 3.917 | 49.70 | 1.000 | 1.097 | null | engineering |
| validation | 7 | original charts | 8/8 | 47.124 | 911.47 | 7.649 | 1.000 | null | engineering |
| validation | 7 | iteration-1 charts | 8/8 | 13.615 | 263.34 | 2.210 | 3.461 | null | engineering |
| validation | 7 | iteration-2 charts | 8/8 | 12.231 | 236.58 | 1.985 | 3.853 | null | engineering |
| validation | 7 | iteration-3 charts | 8/8 | 9.851 | 190.55 | 1.599 | 4.783 | null | engineering |
| validation | 7 | iteration-4 charts | 8/8 | 8.265 | 159.86 | 1.341 | 5.702 | null | engineering |
| validation | 7 | iteration-5 charts | 8/8 | 7.470 | 144.49 | 1.212 | 6.308 | null | engineering |
| validation | 7 | original pair table | 8/8 | 43.473 | 840.86 | 7.056 | 1.000 | null | engineering |
| validation | 7 | candidate pair table | 8/8 | 6.989 | 135.18 | 1.134 | 6.220 | null | engineering |
| validation | 7 | rho | 8/8 | 6.161 | 119.17 | 1.000 | 1.157 | null | engineering |
| validation | 8 | original charts | 8/8 | 48.263 | 933.52 | 7.834 | 1.000 | null | engineering |
| validation | 8 | iteration-1 charts | 8/8 | 12.795 | 247.48 | 2.077 | 3.772 | null | engineering |
| validation | 8 | iteration-2 charts | 8/8 | 12.693 | 245.51 | 2.060 | 3.802 | null | engineering |
| validation | 8 | iteration-3 charts | 8/8 | 10.107 | 195.49 | 1.640 | 4.775 | null | engineering |
| validation | 8 | iteration-4 charts | 8/8 | 8.508 | 164.56 | 1.381 | 5.673 | null | engineering |
| validation | 8 | iteration-5 charts | 8/8 | 7.833 | 151.51 | 1.271 | 6.162 | null | engineering |
| validation | 8 | original pair table | 8/8 | 46.257 | 894.72 | 7.508 | 1.000 | null | engineering |
| validation | 8 | candidate pair table | 8/8 | 7.677 | 148.48 | 1.246 | 6.026 | null | engineering |
| validation | 8 | rho | 8/8 | 6.161 | 119.17 | 1.000 | null | null | engineering |
| validation | 9 | original charts | 8/8 | 102.260 | 1977.93 | 16.598 | 1.000 | null | engineering |
| validation | 9 | iteration-1 charts | 8/8 | 15.848 | 306.53 | 2.572 | 6.453 | null | engineering |
| validation | 9 | iteration-2 charts | 8/8 | 14.344 | 277.44 | 2.328 | 7.129 | null | engineering |
| validation | 9 | iteration-3 charts | 8/8 | 10.931 | 211.44 | 1.774 | 9.355 | null | engineering |
| validation | 9 | iteration-4 charts | 8/8 | 8.255 | 159.67 | 1.340 | 12.387 | null | engineering |
| validation | 9 | iteration-5 charts | 8/8 | 8.001 | 154.76 | 1.299 | 12.780 | null | engineering |
| validation | 9 | original pair table | 8/8 | 98.937 | 1913.67 | 16.059 | 1.000 | null | engineering |
| validation | 9 | candidate pair table | 8/8 | 8.005 | 154.84 | 1.299 | 12.359 | null | engineering |
| validation | 9 | rho | 8/8 | 6.161 | 119.17 | 1.000 | null | null | engineering |
| confirmation | 0 | original charts | 8/8 | 6.758 | 925.76 | 4.048 | 1.000 | null | engineering |
| confirmation | 0 | candidate charts | 8/8 | 1.868 | 255.94 | 1.119 | 3.617 | null | engineering |
| confirmation | 0 | original pair table | 8/8 | 6.076 | 832.28 | 3.639 | 1.000 | null | engineering |
| confirmation | 0 | candidate pair table | 8/8 | 2.101 | 287.80 | 1.258 | 2.892 | null | engineering |
| confirmation | 0 | rho | 8/8 | 1.670 | 228.72 | 1.000 | 1.081 | null | engineering |
| confirmation | 1 | original charts | 0/8 | null | null | null | null | null | engineering |
| confirmation | 1 | candidate charts | 0/8 | null | null | null | null | null | engineering |
| confirmation | 1 | original pair table | 0/8 | null | null | null | null | null | engineering |
| confirmation | 1 | candidate pair table | 0/8 | null | null | null | null | null | engineering |
| confirmation | 1 | rho | 8/8 | 2.583 | 210.34 | 1.000 | 1.217 | null | engineering |
| confirmation | 2 | original charts | 8/8 | 14.716 | 1198.51 | 5.698 | 1.000 | null | engineering |
| confirmation | 2 | candidate charts | 8/8 | 3.339 | 271.96 | 1.293 | 4.407 | null | engineering |
| confirmation | 2 | original pair table | 8/8 | 14.305 | 1164.99 | 5.539 | 1.000 | null | engineering |
| confirmation | 2 | candidate pair table | 8/8 | 4.225 | 344.08 | 1.636 | 3.386 | null | engineering |
| confirmation | 2 | rho | 8/8 | 2.583 | 210.34 | 1.000 | null | null | engineering |
| confirmation | 3 | original charts | 0/8 | null | null | null | null | null | engineering |
| confirmation | 3 | candidate charts | 0/8 | null | null | null | null | null | engineering |
| confirmation | 3 | original pair table | 0/8 | null | null | null | null | null | engineering |
| confirmation | 3 | candidate pair table | 0/8 | null | null | null | null | null | engineering |
| confirmation | 3 | rho | 8/8 | 2.583 | 210.34 | 1.000 | null | null | engineering |
| confirmation | 4 | original charts | 8/8 | 84.574 | 1780.30 | 21.058 | 1.000 | null | engineering |
| confirmation | 4 | candidate charts | 8/8 | 4.997 | 105.18 | 1.244 | 16.926 | null | engineering |
| confirmation | 4 | original pair table | 8/8 | 24.929 | 524.77 | 6.207 | 1.000 | null | engineering |
| confirmation | 4 | candidate pair table | 8/8 | 4.014 | 84.50 | 0.999 | 6.210 | null | engineering |
| confirmation | 4 | rho | 8/8 | 4.016 | 84.54 | 1.000 | 1.105 | null | engineering |
| confirmation | 5 | original charts | 8/8 | 36.400 | 766.22 | 9.063 | 1.000 | null | engineering |
| confirmation | 5 | candidate charts | 8/8 | 5.226 | 110.01 | 1.301 | 6.965 | null | engineering |
| confirmation | 5 | original pair table | 8/8 | 25.908 | 545.37 | 6.451 | 1.000 | null | engineering |
| confirmation | 5 | candidate pair table | 8/8 | 7.255 | 152.72 | 1.806 | 3.571 | null | engineering |
| confirmation | 5 | rho | 8/8 | 4.016 | 84.54 | 1.000 | null | null | engineering |
| confirmation | 6 | original charts | 8/8 | 56.769 | 720.21 | 12.898 | 1.000 | null | engineering |
| confirmation | 6 | candidate charts | 8/8 | 7.478 | 94.87 | 1.699 | 7.592 | null | engineering |
| confirmation | 6 | original pair table | 8/8 | 34.889 | 442.63 | 7.927 | 1.000 | null | engineering |
| confirmation | 6 | candidate pair table | 8/8 | 7.437 | 94.35 | 1.690 | 4.692 | null | engineering |
| confirmation | 6 | rho | 8/8 | 4.401 | 55.84 | 1.000 | 1.086 | null | engineering |
| confirmation | 7 | original charts | 8/8 | 42.131 | 814.91 | 6.346 | 1.000 | null | engineering |
| confirmation | 7 | candidate charts | 8/8 | 7.542 | 145.89 | 1.136 | 5.586 | null | engineering |
| confirmation | 7 | original pair table | 8/8 | 39.405 | 762.17 | 5.935 | 1.000 | null | engineering |
| confirmation | 7 | candidate pair table | 8/8 | 7.279 | 140.80 | 1.096 | 5.413 | null | engineering |
| confirmation | 7 | rho | 8/8 | 6.639 | 128.41 | 1.000 | 1.146 | null | engineering |
| confirmation | 8 | original charts | 8/8 | 43.539 | 842.15 | 6.558 | 1.000 | null | engineering |
| confirmation | 8 | candidate charts | 8/8 | 7.924 | 153.26 | 1.194 | 5.495 | null | engineering |
| confirmation | 8 | original pair table | 8/8 | 42.167 | 815.59 | 6.351 | 1.000 | null | engineering |
| confirmation | 8 | candidate pair table | 8/8 | 7.958 | 153.92 | 1.199 | 5.299 | null | engineering |
| confirmation | 8 | rho | 8/8 | 6.639 | 128.41 | 1.000 | null | null | engineering |
| confirmation | 9 | original charts | 8/8 | 108.839 | 2105.19 | 16.394 | 1.000 | null | engineering |
| confirmation | 9 | candidate charts | 8/8 | 8.458 | 163.60 | 1.274 | 12.868 | null | engineering |
| confirmation | 9 | original pair table | 8/8 | 105.497 | 2040.55 | 15.891 | 1.000 | null | engineering |
| confirmation | 9 | candidate pair table | 8/8 | 8.433 | 163.12 | 1.270 | 12.510 | null | engineering |
| confirmation | 9 | rho | 8/8 | 6.639 | 128.41 | 1.000 | null | null | engineering |

The validation cohort includes the four frozen targets and four reused validation targets. Every retained candidate is shown on those identical inputs; earlier costs are preserved. The confirmation cohort uses the four frozen targets and four new uniformly sampled nonzero scalars. The fresh-only parity gate uses the latter four, and does not pool them with tuning inputs.

### Exclusive candidate phases on fresh confirmation targets

| Case | Curve setup million Ir | Target generation million Ir | Base/plan/index million Ir | Driver/verification million Ir |
|---|---|---|---|---|
| 0 | 0.933 | 0.189 | 0.371 | 0.431 |
| 2 | 1.531 | 0.306 | 0.875 | 0.779 |
| 4 | 1.994 | 0.665 | 0.894 | 2.026 |
| 5 | 1.994 | 0.665 | 1.216 | 1.917 |
| 6 | 2.111 | 0.662 | 1.099 | 4.146 |
| 7 | 4.299 | 0.716 | 0.536 | 2.253 |
| 8 | 4.299 | 0.716 | 0.805 | 2.386 |
| 9 | 4.299 | 0.716 | 0.851 | 3.500 |

### Remaining gap

The original factor base and counting bound are unchanged. Rational-support compression helps the degree-eleven ordinary seed but does not reduce the degree-thirteen seed dimension. The latter still pays substantial work for failed chart queries and relation generation. The final relation-matrix timer is a small part of cold native runtime; setting just that measured timer to zero gives the following Amdahl limits. These are runtime bounds for that scope, not operation-count floors and not bounds on all polynomial preprocessing or chart projection.

| Case | Timed relation-matrix fraction | Maximum speedup from removing that timer |
|---|---|---|
| 0 | 0.76% | 1.008 |
| 2 | 0.42% | 1.004 |
| 4 | 0.42% | 1.004 |
| 5 | 0.84% | 1.008 |
| 6 | 0.33% | 1.003 |
| 7 | 0.36% | 1.004 |
| 8 | 0.36% | 1.004 |
| 9 | 0.37% | 1.004 |

Scopes include curve setup, target generation, factor-base/plan/index setup, and the entire driver plus final verification. Diagnostic JSON and independent exhaustive oracle replay are outside those scopes. Raw profiles, exclusive phase totals, compiler/binary hashes, and calibration controls are archived in `counts/`. No Valgrind runtime is used as a native timing result.

## Secondary native runtime comparisons

Ratios below one favor the candidate. Intervals are paired 95% seed-cluster bootstrap intervals; three repetitions share a seed and are not treated as three independent instances. Iteration rows compare against their immediate predecessor. The confirmation compares directly against the original revision on newly sampled uniform targets. The original four new seed/log pairs became reused validation inputs after iteration 1. No cross-round timing ratios are multiplied.

| Comparison | Case | Candidate/predecessor | Candidate/rho | Rho ratio 95% interval |
|---|---|---|---|---|
| iteration-1 | 0 | 0.544 | 1.701 | [1.555, 1.883] |
| iteration-1 | 1 | — | — | — |
| iteration-1 | 2 | 0.420 | 2.042 | [1.827, 2.190] |
| iteration-1 | 3 | — | — | — |
| iteration-1 | 4 | 0.756 | 16.895 | [12.454, 26.814] |
| iteration-1 | 5 | 0.689 | 8.665 | [6.291, 11.933] |
| iteration-1 | 6 | 0.568 | 5.956 | [3.897, 9.638] |
| iteration-1 | 7 | 0.315 | 1.992 | [1.400, 2.833] |
| iteration-1 | 8 | 0.310 | 2.075 | [1.408, 3.059] |
| iteration-1 | 9 | 0.179 | 2.291 | [1.754, 2.994] |
| iteration-2 | 0 | 0.918 | 1.350 | [1.300, 1.396] |
| iteration-2 | 1 | — | — | — |
| iteration-2 | 2 | 0.953 | 1.732 | [1.386, 2.081] |
| iteration-2 | 3 | — | — | — |
| iteration-2 | 4 | 0.213 | 4.213 | [3.202, 6.139] |
| iteration-2 | 5 | 0.276 | 2.919 | [2.318, 3.807] |
| iteration-2 | 6 | 0.452 | 2.612 | [2.146, 3.179] |
| iteration-2 | 7 | 0.861 | 1.661 | [1.405, 1.956] |
| iteration-2 | 8 | 0.908 | 1.812 | [1.623, 2.029] |
| iteration-2 | 9 | 0.883 | 1.819 | [1.556, 2.127] |
| iteration-3 | 0 | 0.748 | 1.298 | [1.107, 1.524] |
| iteration-3 | 1 | — | — | — |
| iteration-3 | 2 | 0.911 | 1.558 | [1.457, 1.666] |
| iteration-3 | 3 | — | — | — |
| iteration-3 | 4 | 0.429 | 1.812 | [1.543, 2.191] |
| iteration-3 | 5 | 0.696 | 1.809 | [1.589, 2.065] |
| iteration-3 | 6 | 0.883 | 3.262 | [2.422, 4.611] |
| iteration-3 | 7 | 0.881 | 1.571 | [1.072, 2.303] |
| iteration-3 | 8 | 0.804 | 1.516 | [1.038, 2.213] |
| iteration-3 | 9 | 0.897 | 1.745 | [1.213, 2.523] |
| iteration-4 | 0 | 0.791 | 1.081 | [0.810, 1.414] |
| iteration-4 | 1 | — | — | — |
| iteration-4 | 2 | 0.909 | 1.543 | [1.235, 2.243] |
| iteration-4 | 3 | — | — | — |
| iteration-4 | 4 | 0.848 | 1.732 | [1.333, 2.338] |
| iteration-4 | 5 | 0.921 | 1.783 | [1.395, 2.280] |
| iteration-4 | 6 | 0.880 | 2.209 | [2.029, 2.354] |
| iteration-4 | 7 | 0.865 | 1.263 | [1.119, 1.424] |
| iteration-4 | 8 | 0.794 | 1.239 | [1.082, 1.419] |
| iteration-4 | 9 | 0.735 | 1.209 | [1.027, 1.424] |
| iteration-5 | 0 | 0.873 | 1.025 | [0.930, 1.174] |
| iteration-5 | 1 | — | — | — |
| iteration-5 | 2 | 0.859 | 1.535 | [1.397, 1.753] |
| iteration-5 | 3 | — | — | — |
| iteration-5 | 4 | 1.010 | 1.725 | [1.563, 1.929] |
| iteration-5 | 5 | 0.796 | 1.773 | [1.626, 1.928] |
| iteration-5 | 6 | 0.813 | 2.502 | [1.973, 3.301] |
| iteration-5 | 7 | 0.897 | 1.513 | [1.260, 1.817] |
| iteration-5 | 8 | 0.922 | 1.535 | [1.426, 1.701] |
| iteration-5 | 9 | 0.960 | 1.577 | [1.349, 1.814] |
| confirmation | 0 | 0.354 | 1.126 | [0.914, 1.363] |
| confirmation | 1 | — | — | — |
| confirmation | 2 | 0.322 | 1.689 | [1.259, 2.825] |
| confirmation | 3 | — | — | — |
| confirmation | 4 | 0.079 | 1.284 | [0.974, 1.543] |
| confirmation | 5 | 0.173 | 1.229 | [1.028, 1.437] |
| confirmation | 6 | 0.222 | 2.178 | [1.799, 3.108] |
| confirmation | 7 | 0.262 | 1.264 | [1.143, 1.381] |
| confirmation | 8 | 0.241 | 1.315 | [1.219, 1.471] |
| confirmation | 9 | 0.088 | 1.380 | [1.209, 1.578] |

## Limits and reproduction

The two zero-column controls cannot recover a logarithm. Cold single-target costs are reported; there is no hidden warm amortization, precomputed target logarithm, or rho fallback. The small-coordinate adapter remains opt-in, supports at most seven seed coordinates, and uses n<=63 field arithmetic (packed point arithmetic through n<=62). No Redis/network latency or degree-131 solve is inferred from this cache-off local experiment.

The instruction audit checks every non-timing field emitted by the immutable cost harness against the matching native run, including F4 reductions and every oracle witness. A stale shared-target reference build and two truncated raw-profile archives were rejected and retained with reasons. Accepted replacements use isolated package builds and validated, completed profile archives. The first iteration-5 native run overlapped profiling and is retained as a correctness diagnostic; its timing values are excluded from the tables. All native timings remain shared-host diagnostics.

See [README.md](README.md) for hypotheses frozen before each change, [COST_PROTOCOL.md](COST_PROTOCOL.md) for the instruction model, and [support-census.json](support-census.json) for the independent rational-support census. Use `build.py` with a new target directory per revision. Run `run.py`, then `compare.py` (the witness-order flag is required for the support/visitation experiments); `count.py` requires the identically generated `rho_parity_cost` executable and Valgrind. Run `audit.py` before `report.py`, which renders this note and the canonical scoreboard from archived comparisons.
