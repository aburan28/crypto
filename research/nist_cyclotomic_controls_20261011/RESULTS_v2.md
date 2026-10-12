# Controlled re-test: results (v2)

v2, `PROTOCOL.md` as amended by `AMENDMENT_1.md`: larger N only. E2 is not re-run. Seed base 20261011. Protocol file hashes are recorded in `results_v2.json`. Every verdict below is the decision rule applied mechanically.

## E1 — small-factor content (M1 = log2 of 2^32-smooth part found by factor(x, 2^18); M2 = # distinct such primes)

Null: N = 2000 random primes n' in the Hasse interval per curve. Bonferroni α = 0.05/18 = 0.0028. Power check: min attainable p = 2/(N+1) = 0.0010, passes.

| curve | value | metric | target | control median [2.5%, 97.5%] | rank / N | p (two-sided) | verdict |
|:--|:--|:--|--:|:--|--:|--:|:--|
| P-192 | n − 1 | M1 | 17.5 | 29.4 [2.6, 87.2] | 593/2000 | 0.5937 | within control distribution |
| P-192 | n − 1 | M2 | 3.0 | 4.0 [1.0, 8.0] | 343/2000 | 0.8056 | within control distribution |
| P-192 | n + 1 | M1 | 7.9 | 29.9 [2.0, 84.0] | 240/2000 | 0.2409 | within control distribution |
| P-192 | n + 1 | M2 | 3.0 | 4.0 [1.0, 8.0] | 303/2000 | 0.7326 | within control distribution |
| P-192 | twist | M1 | 4.5 | 25.4 [0.0, 76.1] | 292/2000 | 0.3058 | within control distribution |
| P-192 | twist | M2 | 1.0 | 2.0 [0.0, 5.0] | 206/2000 | 0.6557 | within control distribution |
| P-224 | n − 1 | M1 | 29.0 | 29.1 [2.0, 88.8] | 991/2000 | 0.9915 | within control distribution |
| P-224 | n − 1 | M2 | 5.0 | 4.0 [1.0, 8.0] | 1278/2000 | 0.7226 | within control distribution |
| P-224 | n + 1 | M1 | 8.1 | 28.9 [2.0, 85.7] | 238/2000 | 0.2399 | within control distribution |
| P-224 | n + 1 | M2 | 3.0 | 4.0 [1.0, 8.0] | 337/2000 | 0.7816 | within control distribution |
| P-224 | twist | M1 | 59.0 | 27.7 [0.0, 83.1] | 1791/2000 | 0.2099 | within control distribution |
| P-224 | twist | M2 | 5.0 | 3.0 [0.0, 7.0] | 1623/2000 | 0.3778 | within control distribution |
| P-256 | n − 1 | M1 | 125.0 | 28.7 [2.6, 90.2] | 1992/2000 | 0.0090 | within control distribution |
| P-256 | n − 1 | M2 | 10.0 | 4.0 [1.0, 8.0] | 1991/2000 | 0.0100 | within control distribution |
| P-256 | n + 1 | M1 | 41.6 | 28.9 [2.0, 82.8] | 1387/2000 | 0.6137 | within control distribution |
| P-256 | n + 1 | M2 | 4.0 | 4.0 [1.0, 8.0] | 791/2000 | 1.0000 | within control distribution |
| P-256 | twist | M1 | 15.1 | 27.3 [0.0, 85.0] | 566/2000 | 0.5667 | within control distribution |
| P-256 | twist | M2 | 4.0 | 3.0 [0.0, 7.0] | 1308/2000 | 0.6927 | within control distribution |

Positive control (planted 2^20-smooth part of ≈120 bits):

| curve | detection |
|:--|:--|
| P-192 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (87.2) |
| P-224 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (88.8) |
| P-256 | 20/20 planted-smooth values above the control 97.5th percentile of M1 (90.2) |


## E3 — ĥ of the one-point lift: generator vs random points (nats; ĥ := lim h_x([2^k]P)/4^k)

Bonferroni α = 0.05/6 = 0.0083. Power check: min attainable p = 2/(N+1) = 0.0020, passes.

| curve | metric | G̃ | random points median [2.5%, 97.5%] (N) | rank | p | verdict | small-x positive control median [min, max] |
|:--|:--|--:|:--|--:|--:|:--|:--|
| P-192 | ĥ/ln p | 0.9865 | 0.9923 [0.9662, 1.0034] (1000) | 281/1000 | 0.563 | within control distribution | 0.6600 [0.6463, 0.6653] |
| P-192 | h_x − ĥ | -0.5493 | 0.1438 [-0.5493, 2.1972] (1000) | 43/1000 | 0.438 | within control distribution | -65.9053 [-69.2966, -65.0771] |
| P-224 | ĥ/ln p | 0.9858 | 0.9938 [0.9729, 1.0031] (1000) | 188/1000 | 0.378 | within control distribution | 0.6615 [0.6471, 0.6681] |
| P-224 | h_x − ĥ | 1.8637 | 0.1438 [-0.5493, 2.2034] (1000) | 956/1000 | 0.090 | within control distribution | -80.8559 [-83.0641, -78.5936] |
| P-256 | ĥ/ln p | 0.9982 | 0.9941 [0.9746, 1.0025] (1000) | 737/1000 | 0.527 | within control distribution | 0.6625 [0.6444, 0.6684] |
| P-256 | h_x − ĥ | -0.5493 | 0.1438 [-0.5493, 2.1715] (1000) | 0/1000 | 0.110 | within control distribution | -96.1205 [-101.4126, -92.3954] |

## E4 — two-point-lift pairing: (G̃, R̃) vs random pairs

Bonferroni α = 0.05/6 = 0.0083. Power check: min attainable p = 2/(N+1) = 0.0040, passes.

| curve | metric | target | random pairs median [2.5%, 97.5%] (N) | rank | p | verdict | dependent-pair positive control ([2]G̃; [3]G̃) |
|:--|:--|--:|:--|--:|--:|:--|:--|
| P-192 | |cos| | 0.00071 | 0.00130 [0.00071, 0.00904] (500) | 14/500 | 0.060 | within control distribution | 1.002460; 1.001076 |
| P-192 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99992, 1.00000] (500) | 486/500 | 0.060 | within control distribution | -4.93e-3; -2.15e-3 |
| P-224 | |cos| | 0.00129 | 0.00130 [0.00071, 0.00740] (500) | 132/500 | 0.531 | within control distribution | 1.002456; 1.001071 |
| P-224 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99995, 1.00000] (500) | 368/500 | 0.531 | within control distribution | -4.92e-3; -2.14e-3 |
| P-256 | |cos| | 0.00173 | 0.00130 [0.00051, 0.00592] (500) | 363/500 | 0.551 | within control distribution | 1.002447; 1.001081 |
| P-256 | det/(ĥĥ) | 1.00000 | 1.00000 [0.99996, 1.00000] (500) | 137/500 | 0.551 | within control distribution | -4.90e-3; -2.16e-3 |

---
Run time 8966 s. Raw per-sample rows: `rows_v2.jsonl`.
