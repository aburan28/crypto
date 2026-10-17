# External comparison: isogeny_algos vs FLINT and PARI/GP

Best time over all rounds per tool, in microseconds; ratio = ours / theirs (below 1: ours is faster). `same` = every tool and round printed the same digest.

## mul

| p | n | ours | FLINT | PARI | ours/FLINT | ours/PARI | outputs |
|:--|--:|--:|--:|--:|--:|--:|:--|

## sqr

| p | n | ours | FLINT | PARI | ours/FLINT | ours/PARI | outputs |
|:--|--:|--:|--:|--:|--:|--:|:--|

## powx

| p | n | ours | FLINT | PARI | ours/FLINT | ours/PARI | outputs |
|:--|--:|--:|--:|--:|--:|--:|:--|

## Modular polynomials Phi_l mod p, all primes 3 <= l <= 89, p = 2^127 - 1

| computation | best total over rounds |
|:--|--:|
| PARI: Phi_l(j, Y) mod p for one j (polmodular(l, 0, Mod(j, p))) | 9089 ms |
| PARI: full Phi_l over Z, reduced mod p (polmodular(l)) | 9523 ms |
| ours: full Phi_l mod p (Hecke, compute_many) | 1311 ms |

Digests Phi_l(987654321, 123456789) mod p: 23 of 23 levels agree across tools.

## SEA point counting, p61

| curve k | ours ms | PARI (seadata) ms | PARI (no seadata) ms | orders |
|--:|--:|--:|--:|:--|
| 0 | 1.1 | 22.0 | 133.0 | same |
| 1 | 1.0 | 3.0 | 122.0 | same |
| 2 | 1.0 | 3.0 | 128.0 | same |
| 3 | 0.9 | 3.0 | 119.0 | same |
| 4 | 0.9 | 3.0 | 99.0 | same |
| 5 | 0.7 | 3.0 | 102.0 | same |
| 6 | 0.9 | 3.0 | 107.0 | same |
| 7 | 1.1 | 3.0 | 106.0 | same |
| 8 | 1.0 | 3.0 | 99.0 | same |
| 9 | 1.1 | 3.0 | 110.0 | same |
| 10 | 0.9 | 3.0 | 112.0 | same |
| 11 | 1.0 | 3.0 | 121.0 | same |
| 12 | 0.7 | 3.0 | 114.0 | same |
| 13 | 1.0 | 3.0 | 107.0 | same |
| 14 | 0.9 | 3.0 | 101.0 | same |
| 15 | 1.1 | 3.0 | 95.0 | same |
| 16 | 1.0 | 3.0 | 104.0 | same |
| 17 | 0.7 | 3.0 | 106.0 | same |
| 18 | 0.6 | 3.0 | 120.0 | same |
| 19 | 0.6 | 3.0 | 137.0 | same |
| **sum** | **18** | **79** | **2242** | 20/20 agree |

Ours needs Phi_l for l up to 61 computed once per field first: 75 ms (PARI with seadata reads precomputed polynomials from disk; without it, PARI computes what it needs per curve, inside the times above).

## SEA point counting, p127

| curve k | ours ms | PARI (seadata) ms | PARI (no seadata) ms | orders |
|--:|--:|--:|--:|:--|
| 0 | 67.1 | 90.0 | 565.0 | same |
| 1 | 55.9 | 61.0 | 549.0 | same |
| 2 | 48.6 | 67.0 | 627.0 | same |
| 3 | 30.1 | 67.0 | 526.0 | same |
| 4 | 47.3 | 57.0 | 447.0 | same |
| 5 | 40.3 | 58.0 | 426.0 | same |
| 6 | 56.6 | 58.0 | 446.0 | same |
| 7 | 82.1 | 63.0 | 520.0 | same |
| 8 | 59.7 | 64.0 | 467.0 | same |
| 9 | 57.6 | 62.0 | 479.0 | same |
| **sum** | **545** | **647** | **5052** | 10/10 agree |

Ours needs Phi_l for l up to 89 computed once per field first: 1296 ms (PARI with seadata reads precomputed polynomials from disk; without it, PARI computes what it needs per curve, inside the times above).

**Output mismatches: 0.**
