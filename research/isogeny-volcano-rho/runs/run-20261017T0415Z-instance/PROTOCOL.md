# Instance study: does rho's cost depend on the generator or on k?

Frozen before any measured trial (smoke run of 8 trials per arm discarded).

- Curve C of this study: p = 146833024135393, y² = x³ + 14391545271884x + 141956518099848,
  #E = 146833030777920 = 960 · n, n = 152951073727 (prime).
- Instrument: `volcano-rho --instance-study`, its native single-walker 32-adding DP
  rho (no negation map), the same binary for every arm. 1,000 solves per arm,
  `--seed 20261018`, 4 threads; iterations (group additions) are the metric.
- Arms (trial t draws what its arm varies from its own seed):
  1. `walk`: P and k fixed, walk seed redrawn (what the volcano runs measured).
  2. `k`: P and walk seed fixed, k uniform in [1, n).
  3. `generator`: k and walk seed fixed, P a fresh point of order n.
  4. `all`: walk seed, k and P all redrawn.
  5. `small_k`: P and walk seed fixed, k uniform in [1, 2^16].
- Null: every arm draws from one distribution, mean ≈ √(πn/2) · 1/√(1 − 1/32)
  ≈ 490,158 · 1.016 ≈ 498,000 and sd/mean ≈ 0.52.
- Decision: an arm differs from `walk` only if its Welch z exceeds 2.58 in size
  (four comparisons, ≈ Bonferroni at 0.05) or its KS p is below 0.0125.
  Otherwise report "no detectable difference" with the interval on the mean
  difference as the size of effect excluded.
- Every solve must return its own k; the binary exits non-zero otherwise.
- Replay: re-run with `--trials 5` and match the first five iteration counts of every arm.
