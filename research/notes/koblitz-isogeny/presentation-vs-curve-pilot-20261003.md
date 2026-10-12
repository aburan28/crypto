# Pilot: presentation-vs-curve experiment — feasibility and budgets

These are pilot timings for
[`presentation-vs-curve-design-20261003.md`](presentation-vs-curve-design-20261003.md).
They use the Koblitz curve only, with a handful of probes per setting. **No
hypothesis is tested here.** Raw rows are in
`experiments/koblitz_ic_pilot.jsonl`.

Runs: the default subspace V = ⟨1, z, …, z^{l−1}⟩ ran alone. Random V ran
three jobs at a time on 4 cores, so its wall-clock numbers may include
contention. Its reductions-per-call figures do not.

| n | m | l | V | \|F\| | µs/call | reductions/call | first fall | relations / probes | core-h per 100 relations |
|---|---|---|---|---|---|---|---|---|---|
| 17 | 3 | 6 | default | 32 | **2.9 × 10⁸** | 472 | 2 | 1 / 2 | 24.6 |
| 19 | 2 | 9 | default | 249 | 22 927 | 1.83 | 3 | 69 / 400 | 0.003 |
| 19 | 2 | 9 | random ×4 | 246–276 | 61 000–64 000 | 3.3–3.4 | 3 | 88–95 / 400 | 0.006–0.008 |
| 19 | 2 | 10 | default | 522 | 89 443 | 6.06 | 3 | 137 / 200 | 0.002 |
| 23 | 2 | 11 | default | 1 009 | 83 696 | 1.76 | 3 | 16 / 100 | 0.010 |
| 23 | 2 | 12 | default | 2 021 | 690 176 | 7.7 | 3 | 24 / 50 | 0.020 |
| 31 | 2 | 10 | default | 527 | **510** | 1.0 | 2 | 0 / 2000 | 0.055 |
| 31 | 2 | 10 | random ×2 | 489–518 | **75 000–92 000** | 1.0 | 2 | 0 / 40 | 9.4–10.2 |
| 31 | 2 | 11 | default | 1 051 | 2 879 | 1.0 | 2 | 0 / 20 | 0.078 |
| 31 | 2 | 12 | default | 2 123 | 399 492 | 1.0 | 2 | 0 / 20 | 2.6 |
| 31 | 2 | 14 | default | — | > 2.2 × 10⁸ (8 probes timed out at 30 min) | — | — | — | — |
| 31 | 2 | 16 | default | 32 944 | 4.1 × 10⁷ | **20.5** | — | 2 / 4 | 1.1 (4 probes; noisy) |
| 31 | 4 | 6 | default | 35 | system not built (`unknown_probes = 1`) | — | — | — | — |

Exact and full-group probes at n = 19 match within noise (88–101 relations
against 88–95). That is expected, since the cofactor there is 2. They serve
as a negative control for the new instrument.

## What the pilot settles

1. **m ≥ 3 is not feasible with this solver.**
   - m = 3 takes about 295 s per call even at n = 17.
   - m = 4 is not supported: the decomposition system is not built.
   - The design's m ∈ {3, 4} arms are dropped under the pre-registered
     budget rule (≤ 2 core-h per cell). This is a solver limit, not a
     finding about the curves.
2. **At n = 31, Gröbner does no real work at any feasible l ≤ 12.** Each
   call is 1.0 reduction at first fall degree 2. The Gröbner step becomes
   nontrivial only at the balanced point l = 16 (20.5 reductions/call,
   41 s/call). At n = 19 and 23 it is already nontrivial at l ≈ n/2
   (first fall degree 3, 2–8 reductions/call).
3. **The choice of subspace changes cost by orders of magnitude.**
   - At n = 31, l = 10, a random V costs about 150× more per call than the
     monomial V.
   - At n = 19, l = 9, a random V costs about 2.7× more per call and about
     1.8× more reductions per call. Unlike wall clock, the reduction count
     is free of contention.
   - The monomial subspace is therefore a special, cheap setup, and any
     curve-level claim must be made across subspaces, as the design
     requires.
4. **The yield heuristic** (2|F|)^m / (m!·#E) is about 1.3–1.5× optimistic
   where it can be checked: observed/predicted is 0.66–0.77 at n = 19 and
   23.

## Budget consequences, applying the pre-registered rule

| class | setting | cells (curves × V) | projected core-h | verdict |
|---|---|---|---|---|
| n = 16, 17, 19, 23 (full classes) | m = 2, l ≈ n/2 | all members × 4 | ≲ 0.01 per cell, about 10–40 total | **run** |
| n = 31, random V, l = 10 | m = 2 | 64 × 4 | about 10 per cell → about 2 500 | **exceeds rule** |
| n = 31, l = 16 | m = 2 | 64 × 4 | ≥ 1.1 per cell for 100 relations, random V untimed | **needs one more pilot** (random V at l = 16) before committing |
| any n, m ≥ 3 | — | — | — | **dropped (solver)** |

Recommendation: run the full design at n = 16–23 now. Gate n = 31 on one
random-V timing at l = 16. If that fits, use a reduced n = 31 arm of
16 curves × 2 V. Otherwise report n = 31 as out of budget for this solver.
Either way, the Gröbner channel at n = 31 is only informative at l = 16.
