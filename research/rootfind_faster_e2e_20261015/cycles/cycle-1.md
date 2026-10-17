# Cycle 1 — single-hash probe path + multi-target e2e

**Change.** Walking `m = 2` descent and blocked `m = 3`/`m = 4` scans used to
call `pair_filter_hash` independently for prefetch, presence-filter admit, and
compact bucket/rest. They now hash **once** per key and reuse it
(`contains_key_hash`, `prefetch_key_hash`, `pairs_for_key_hash`).

Also: `semaev_rootfind_e2e_bench` gains a **real** multi-target section
(shared FB logs + `PairSumTable` via `IndividualLogSolver`). Knobs:
`E2E_MULTI_N`, `E2E_BINARY_N`, `E2E_SKIP_MULTI`, `E2E_SKIP_PRIME`.

**Hypothesis.** Miss-dominated probe loops get cheaper CPU without changing
answers; multi-target e2e shows amortized setup.

**Verify.** `contains_key_hash_matches_contains_key`; every multi-target log
checked `[d]G = Q`.

## Measure

| artifact | result |
|---|---|
| unit test | ok (hash ≡ contains_key) |
| miss probe (`logs/cycle1-probe-bench.txt`) | shared **1.38 ns** vs doubled **4.08 ns** → **2.95×** (n=13, 8e6 misses, 0 hits) |
| e2e multi (`logs/cycle1-e2e-multi.txt`) | n=9 N=8: batch 1.89 ms vs cold×8 5.54 ms → **3×**; n=13 N=8: 1.18 s vs 5.35 s → **5×** |
| e2e binary n=9 | root mean 602 µs vs F4 2.31 ms → **4×** (CPU) |

## Outcome

**Win.** Probe miss path ~3× cheaper; multi-target amortization 3–5× on
n∈{9,13}. No answer change. Classification: **engineering**. Does not move
cold single-target vs rho (~1.54× hard stop still stands).
