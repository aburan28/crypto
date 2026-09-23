# n=53: shared-factor-log IC against batched signed-Frobenius rho

The working paper's n=53 amortization result reported IC 25.26x faster than rho on
1,024 targets. That comparison used **independent** rho: every target got a fresh
table, and its jump points depended on that target's Q. This lab reruns the
comparison against rho that also shares work across targets.

When rho shares work across targets, IC does not win at any batch size tested
(32–8,192). It uses 10 GB against rho's 9–395 MB. IC keeps one advantage: at
equal precomputation time, its warm per-target queries are about 2x cheaper.

## Arms

- **IC**: the frozen producer `../autolab_n53_eta_sweep_20260912/frozen/current/koblitz_rank_fixture`
  (sha256 6884b57d…), unchanged. It was rerun here with the frozen batch1024
  settings and reproduced the base hash, the 1,244 relations and the 103,361,536
  support queries exactly.
- **Batched rho**: `frozen/v*/koblitz_rho_batch_ks`. Its field arithmetic,
  canonicalization, 32-way partition and signed-Frobenius quotient (size 106) are
  copied verbatim from the packed independent comparator.
  - Kuhn–Struik batching: jump points are multiples of G only, and one
    distinguished-point table persists across targets.
  - v2: successive walk starts step by a fixed stride, one addition per walk.
  - v3: optional Bernstein–Lange precomputation, where G-only walks are stored
    with known logs, and an optional frozen table.
- **Targets** use the same derivations as the frozen runs (`KIC_BATCH_CORPUS` /
  `KIC_RHO_BATCH_CORPUS`). In every run, every fixture's scalar, public point and
  recovered scalar agree between the arms.

## 1. Rerun of the paper's own corpora (producer v1, d=8)

3 blocks per size, alternating which arm runs first.

| L | IC wall (median) | Batched rho wall | IC / rho (range) | Peak RSS IC / rho |
|---:|---:|---:|---:|---:|
| 32 | 38.1 s | 9.0 s | 4.36 (3.88–4.39) | 10.05 GB / 9.8 MB |
| 1,024 | 55.9 s | 43.6 s | 1.25 (1.23–1.28) | 10.64 GB / 19 MB |

Batched rho took 20,208,130 group additions for 1,024 targets, against 558,876,665
for independent rho on the same corpus. That matches Kuhn–Struik's ≈√(2Lℓ/|Aut|).

## 2. Scaling panel on fresh corpora (producer v2, d=4, `protocol.json`)

The distinguished-point density d=4 was chosen on a disjoint tuning corpus
(`n53-ks-dp-tuning-v1`) by counting group additions, not wall time.

| L | IC wall (median) | Rho wall | IC / rho wall (range) | IC / rho charged (range) | Verdict (all 3 blocks) |
|---:|---:|---:|---:|---:|---|
| 1,024 | 50.1 s | 38.0 s | 1.32 (1.32–1.45) | 1.15 (1.14–1.26) | rho faster |
| 2,048 | 78.9 s | 56.5 s | 1.40 (1.24–1.48) | 1.14 (1.01–1.21) | rho faster |
| 4,096 | 119.1 s | 87.3 s | 1.29 (1.26–1.94) | 1.02 (0.98–1.60) | wall: rho; charged: unresolved |
| 8,192 | 227.2 s | 132.7 s | 1.74 (1.71–1.82) | 1.35 (1.35–1.36) | rho faster |

How each charged figure is built:
- IC charged is `full_algorithm_charged_total_ms`. It excludes reference validation (~4 ms/target) and fixture generation (~2 ms/target).
- Rho charged is setup plus the sum of per-target `total_ms`. It includes target generation and one verification scalar multiplication per target.

The two cost curves have different shapes:
- IC's cost per extra target is flat, averaging about 13 ms. The paper's 7.49 ms is the median.
- Batched rho's cost per extra target falls roughly as 1/√i. Its median over the second half of targets was 17.8 ms at L=1,024 and 6.7 ms at L=8,192.

So the gap widens with L, and this panel shows no crossover.

## 3. Warm queries after precomputation (producer v3, d=6, frozen table)

Setup: `bl_warm.sh`, corpus `n53-bl-warm-1024-v1`, 1,024 targets, the same targets for both arms.

| Arm | Precompute | Warm median | Mean | p95 | p99 | Max | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| IC (shared factor logs) | 35.2 s | 8.18 ms | 11.18 | 30.40 | 49.35 | 80.09 | 10.03 GB |
| BL rho, 0.78x IC precompute | 27.3 s | 17.43 ms | 24.48 | 71.29 | 115.68 | 140.18 | 57 MB |
| BL rho, 3.11x IC precompute | 109.6 s | 4.18 ms | 6.34 | 20.12 | 30.74 | 41.70 | 190 MB |
| BL rho, 13.2x IC precompute | 465.9 s | 1.58 ms | 2.15 | 5.86 | 9.20 | 30.78 | 395 MB |

In the BL rows, every target was solved from the precomputed table.

Bernstein–Lange warm cost falls roughly as 1/precompute. Interpolating, BL rho
matches IC's warm median at about 1.6x IC's precompute time, and its mean at about
1.7x. At equal precompute time, IC's warm queries are about 2x cheaper. That is the
one measured regime where IC is better, and it uses 10 GB to get there, against 57–395 MB
for the Bernstein–Lange rows (25–175x less).

## Admissible conclusions

- The 25.26x n=53 batch advantage holds only against unbatched rho. It should not
  be cited as an IC advantage.
- Against Kuhn–Struik batched rho on this host, IC does not win on process wall at
  L ∈ {32, 1,024, 2,048, 4,096, 8,192}. It wins on charged time in only one block
  of twelve (L=4,096, by 1.8%).
- With precomputation excluded, IC's warm query beats a Bernstein–Lange table built
  with equal time. A table built with ~1.7x the time matches IC, and ~3x beats it.

## Limits

- The host was shared with other sessions. Load averages were 7–27 and are
  recorded per run. IC's L=8,192 runs showed 20+ s of system time.
- Rho group-addition counts are deterministic and identical across blocks.
- The three arms are one implementation each on one host, at n=53 only.
- Rho is single-threaded, as is IC.
- Bernstein–Lange here stores the first DP of each walk and does not select the
  most useful DPs, so a tuned table would do better.
- None of this is an asymptotic statement.
- `runs/ks1024_cyclebug.jsonl` is a superseded v1 run with a 2-cycle detection bug
  that wasted 72% of its steps. It is kept only as a record.
