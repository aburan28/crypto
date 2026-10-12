# Controlled re-test: results (v1)

v1, `PROTOCOL.md` exactly as frozen. Its N fail the power check, so no rank test here can reach its α (see the power column). Seed base 20261011. Protocol file hashes are recorded in `results_v1.json`. Every verdict below is the decision rule applied mechanically.

## E1 — small-factor content (M1 = log2 of 2^32-smooth part found by factor(x, 2^18); M2 = # distinct such primes)

Null: N = 400 random primes n' in the Hasse interval per curve. Bonferroni α = 0.05/18 = 0.0028. Power check: min attainable p = 2/(N+1) = 0.0050, **FAILS (no target can be declared unusual)**.

| curve | value | metric | target | control median [2.5%, 97.5%] | rank / N | p (two-sided) | verdict |
|:--|:--|:--|--:|:--|--:|--:|:--|
| P-192 | n − 1 | M1 | 17.5 | 29.1 [2.6, 89.9] | 117/400 | 0.5885 | within control distribution |
| P-192 | n − 1 | M2 | 3.0 | 4.0 [2.0, 8.0] | 70/400 | 0.8429 | within control distribution |
| P-192 | n + 1 | M1 | 7.9 | 30.6 [2.0, 83.9] | 46/400 | 0.2344 | within control distribution |
| P-192 | n + 1 | M2 | 3.0 | 4.0 [1.0, 7.0] | 61/400 | 0.7282 | within control distribution |
| P-192 | twist | M1 | 4.5 | 23.9 [0.0, 74.9] | 52/400 | 0.2793 | within control distribution |
| P-192 | twist | M2 | 1.0 | 2.0 [0.0, 5.0] | 39/400 | 0.6983 | within control distribution |
| P-224 | n − 1 | M1 | 29.0 | 27.9 [2.0, 89.7] | 206/400 | 0.9726 | within control distribution |
| P-224 | n − 1 | M2 | 5.0 | 4.0 [1.0, 7.0] | 262/400 | 0.6933 | within control distribution |
| P-224 | n + 1 | M1 | 8.1 | 28.7 [2.6, 78.0] | 40/400 | 0.2045 | within control distribution |
| P-224 | n + 1 | M2 | 3.0 | 4.0 [2.0, 8.0] | 60/400 | 0.7581 | within control distribution |
| P-224 | twist | M1 | 59.0 | 26.7 [0.0, 76.4] | 360/400 | 0.2045 | within control distribution |
| P-224 | twist | M2 | 5.0 | 3.0 [0.0, 6.0] | 336/400 | 0.3242 | within control distribution |
| P-256 | n − 1 | M1 | 125.0 | 28.7 [2.6, 84.3] | 399/400 | 0.0100 | within control distribution |
| P-256 | n − 1 | M2 | 10.0 | 4.0 [2.0, 8.0] | 399/400 | 0.0100 | within control distribution |
| P-256 | n + 1 | M1 | 41.6 | 31.4 [3.6, 88.2] | 258/400 | 0.7132 | within control distribution |
| P-256 | n + 1 | M2 | 4.0 | 4.0 [1.0, 8.0] | 149/400 | 1.0000 | within control distribution |
| P-256 | twist | M1 | 15.1 | 27.9 [0.0, 84.3] | 118/400 | 0.5935 | within control distribution |
| P-256 | twist | M2 | 4.0 | 3.0 [0.0, 7.0] | 246/400 | 0.7731 | within control distribution |

Positive control (planted 2^20-smooth part of ≈120 bits):

| curve | detection |
|:--|:--|
| P-192 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (89.9) |
| P-224 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (89.7) |
| P-256 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (84.3) |

## E2 — census power (|x| ≤ 10^6; x = u/d², 2 ≤ d ≤ 40, |u| ≤ 10^5)

| curve | planted integer found | planted d=2 found | planted d=3 found | random-b curves with any find | NIST b finds |
|:--|--:|--:|--:|--:|--:|
| P-192 | 10/10 | 10/10 | 10/10 | 0/10 | 0 |
| P-224 | 10/10 | 10/10 | 10/10 | 0/10 | 0 |
| P-256 | 10/10 | 10/10 | 10/10 | 0/10 | 0 |

## E3 — ĥ of the one-point lift: generator vs random points (nats; ĥ := lim h_x([2^k]P)/4^k)

Bonferroni α = 0.05/6 = 0.0083. Power check: min attainable p = 2/(N+1) = 0.0100, **FAILS (no target can be declared unusual)**.

| curve | metric | G̃ | random points median [2.5%, 97.5%] (N) | rank | p | verdict | small-x positive control median [min, max] |
|:--|:--|--:|:--|--:|--:|:--|:--|
| P-192 | ĥ/ln p | 0.9865 | 0.9922 [0.9616, 1.0023] (200) | 59/200 | 0.597 | within control distribution | 0.6613 [0.6438, 0.6694] |
| P-192 | h_x − ĥ | -0.5493 | 0.1438 [-0.5493, 2.4836] (200) | 11/200 | 0.498 | within control distribution | -66.5123 [-68.9485, -64.2055] |
| P-224 | ĥ/ln p | 0.9858 | 0.9943 [0.9751, 1.0032] (200) | 40/200 | 0.408 | within control distribution | 0.6611 [0.6453, 0.6691] |
| P-224 | h_x − ĥ | 1.8637 | 0.1438 [-0.5493, 2.0999] (200) | 191/200 | 0.100 | within control distribution | -81.0740 [-83.1827, -79.1396] |
| P-256 | ĥ/ln p | 0.9982 | 0.9947 [0.9733, 1.0024] (200) | 145/200 | 0.557 | within control distribution | 0.6625 [0.6508, 0.6683] |
| P-256 | h_x − ĥ | -0.5493 | 0.0000 [-0.5493, 2.8883] (200) | 0/200 | 0.129 | within control distribution | -96.0648 [-99.8508, -93.5027] |

## E4 — two-point-lift pairing: (G̃, R̃) vs random pairs

Bonferroni α = 0.05/6 = 0.0083. Power check: min attainable p = 2/(N+1) = 0.0198, **FAILS (no target can be declared unusual)**.

| curve | metric | target | random pairs median [2.5%, 97.5%] (N) | rank | p | verdict | dependent-pair positive control ([2]G̃; [3]G̃) |
|:--|:--|--:|:--|--:|--:|:--|:--|
| P-192 | |cos| | 0.00071 | 0.00130 [0.00128, 0.00688] (100) | 0/100 | 0.020 | within control distribution | 1.002460; 1.001076 |
| P-192 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99995, 1.00000] (100) | 100/100 | 0.020 | within control distribution | -4.93e-3; -2.15e-3 |
| P-224 | |cos| | 0.00129 | 0.00130 [0.00070, 0.01129] (100) | 22/100 | 0.455 | within control distribution | 1.002456; 1.001071 |
| P-224 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99987, 1.00000] (100) | 78/100 | 0.455 | within control distribution | -4.92e-3; -2.14e-3 |
| P-256 | |cos| | 0.00173 | 0.00130 [0.00046, 0.00795] (100) | 69/100 | 0.634 | within control distribution | 1.002447; 1.001081 |
| P-256 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99994, 1.00000] (100) | 31/100 | 0.634 | within control distribution | -4.90e-3; -2.16e-3 |

---
Run time 2804 s. Raw per-sample rows: `rows_v1.jsonl`.
