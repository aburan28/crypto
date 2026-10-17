# endo-chain-bench (counting run): arkworks/ed_on_bls12_381_bandersnatch

curve: arkworks/ed_on_bls12_381_bandersnatch (255 -bit p, 253-bit n, cofactor 4, D_K = -8); doubling general a: dbl-2007-bl 1M + 8S + 1M for a; inversion 164 M + 255 S
chain set: research/endo_native_bench_20261016/frozen_inputs/bandersnatch.chains.constants.json (sha256 12100c87b6ffe5a3469f3c8ab625f2e8f23411df025a01bf3c06136cd646114b); M = 256 pairs, seed 20261016

| arm | M | S | inversions | M_eq = M + S (mean) | min | max | model M_eq | counted / model − 1 | correct |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| baseline/affine/w3 | 1127.28 | 2528.69 | 1.000 | 3655.97 | 3531 | 3736 | 3655.28 | +0.019 % | yes |
| baseline/affine/w4 | 1073.10 | 2487.22 | 1.000 | 3560.32 | 3478 | 3638 | 3558.45 | +0.053 % | yes |
| baseline/affine/w5 | 1080.60 | 2474.23 | 1.000 | 3554.83 | 3473 | 3631 | 3551.39 | +0.097 % | yes |
| baseline/affine/w6 | 1173.62 | 2494.97 | 1.000 | 3668.59 | 3583 | 3728 | 3668.18 | +0.011 % | yes |
| baseline/affine/w7 | 1413.00 | 2567.86 | 1.000 | 3980.86 | 3896 | 4051 | 3980.70 | +0.004 % | yes |
| baseline/jacobian/w3 | 1203.69 | 2334.29 | 0.000 | 3537.98 | 3370 | 3650 | 3536.66 | +0.037 % | yes |
| baseline/jacobian/w4 | 1087.63 | 2278.35 | 0.000 | 3365.98 | 3258 | 3468 | 3363.07 | +0.087 % | yes |
| baseline/jacobian/w5 | 1036.96 | 2252.82 | 0.000 | 3289.78 | 3190 | 3388 | 3285.34 | +0.135 % | yes |
| baseline/jacobian/w6 | 1057.87 | 2259.53 | 0.000 | 3317.40 | 3224 | 3394 | 3317.33 | +0.002 % | yes |
| baseline/jacobian/w7 | 1183.36 | 2311.95 | 0.000 | 3495.31 | 3400 | 3580 | 3495.09 | +0.006 % | yes |
| glv/0+1w/2/generic/affine/w3 | 907.62 | 1535.12 | 1.000 | 2442.75 | 2350 | 2547 | 2441.23 | +0.062 % | yes |
| glv/0+1w/2/generic/affine/w4 | 887.06 | 1506.08 | 1.000 | 2393.14 | 2313 | 2464 | 2393.42 | -0.012 % | yes |
| glv/0+1w/2/generic/affine/w5 | 963.44 | 1518.38 | 1.000 | 2481.81 | 2420 | 2537 | 2481.70 | +0.005 % | yes |
| glv/0+1w/2/generic/affine/w6 | 1192.64 | 1587.08 | 1.000 | 2779.71 | 2703 | 2829 | 2780.73 | -0.037 % | yes |
| glv/0+1w/2/generic/affine/w7 | 1704.87 | 1759.34 | 1.000 | 3464.21 | 3395 | 3509 | 3467.04 | -0.082 % | yes |
| glv/0+1w/2/generic/jacobian/w3 | 972.69 | 1338.89 | 0.000 | 2311.58 | 2175 | 2457 | 2309.12 | +0.106 % | yes |
| glv/0+1w/2/generic/jacobian/w4 | 877.92 | 1293.29 | 0.000 | 2171.21 | 2065 | 2271 | 2171.27 | -0.003 % | yes |
| glv/0+1w/2/generic/jacobian/w5 | 872.50 | 1289.14 | 0.000 | 2161.64 | 2081 | 2233 | 2161.59 | +0.002 % | yes |
| glv/0+1w/2/generic/jacobian/w6 | 981.71 | 1335.85 | 0.000 | 2317.56 | 2227 | 2383 | 2318.51 | -0.041 % | yes |
| glv/0+1w/2/generic/jacobian/w7 | 1284.05 | 1471.64 | 0.000 | 2755.70 | 2675 | 2809 | 2758.97 | -0.119 % | yes |
| glv/0+1w/2/optimised/affine/w3 | 905.62 | 1534.12 | 1.000 | 2439.75 | 2347 | 2544 | 2437.23 | +0.103 % | yes |
| glv/0+1w/2/optimised/affine/w4 | 885.06 | 1505.08 | 1.000 | 2390.14 | 2310 | 2461 | 2389.42 | +0.030 % | yes |
| glv/0+1w/2/optimised/affine/w5 | 961.44 | 1517.38 | 1.000 | 2478.81 | 2417 | 2534 | 2477.70 | +0.045 % | yes |
| glv/0+1w/2/optimised/affine/w6 | 1190.64 | 1586.08 | 1.000 | 2776.71 | 2700 | 2826 | 2776.73 | -0.001 % | yes |
| glv/0+1w/2/optimised/affine/w7 | 1702.87 | 1758.34 | 1.000 | 3461.21 | 3392 | 3506 | 3463.04 | -0.053 % | yes |
| glv/0+1w/2/optimised/jacobian/w3 | 970.69 | 1337.89 | 0.000 | 2308.58 | 2172 | 2454 | 2305.12 | +0.150 % | yes |
| glv/0+1w/2/optimised/jacobian/w4 | 875.92 | 1292.29 | 0.000 | 2168.21 | 2062 | 2268 | 2167.27 | +0.043 % | yes |
| glv/0+1w/2/optimised/jacobian/w5 | 870.50 | 1288.14 | 0.000 | 2158.64 | 2078 | 2230 | 2157.59 | +0.049 % | yes |
| glv/0+1w/2/optimised/jacobian/w6 | 979.71 | 1334.85 | 0.000 | 2314.56 | 2224 | 2380 | 2314.51 | +0.002 % | yes |
| glv/0+1w/2/optimised/jacobian/w7 | 1282.05 | 1470.64 | 0.000 | 2752.70 | 2672 | 2806 | 2754.97 | -0.083 % | yes |
| glv/1+1w/3/optimised/jacobian/w5 | 870.70 | 1286.09 | 0.000 | 2156.79 | 2062 | 2246 | 2157.18 | -0.018 % | yes |
| phi/0+1w/2/generic | 7.00 | 2.00 | 0.000 | 9.00 | 9 | 9 | 10.00 | -10.000 % | yes |
| phi/0+1w/2/optimised | 5.00 | 1.00 | 0.000 | 6.00 | 6 | 6 | 6.00 | +0.000 % | yes |
| phi/1+1w/3/generic | 24.00 | 2.00 | 0.000 | 26.00 | 26 | 26 | 26.00 | +0.000 % | yes |
| phi/1+1w/3/optimised | 6.00 | 0.00 | 0.000 | 6.00 | 6 | 6 | 6.00 | +0.000 % | yes |
| decompose/0+1w/2 | 0.00 | 0.00 | 0.000 | 0.00 | 0 | 0 | 0 | — | yes |
| decompose/1+1w/3 | 0.00 | 0.00 | 0.000 | 0.00 | 0 | 0 | 0 | — | yes |
| aa/jacobian/w5 | 1036.96 | 2252.82 | 0.000 | 3289.78 | 3190 | 3388 | 3285.34 | +0.135 % | yes |

all correctness checks passed: yes
json record written to research/endo_native_bench_20261016/bandersnatch.counts.json
