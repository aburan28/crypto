# Cycle 2 — m=2 hit path reuses known key/hash

**Change.** After the walking descent admits an `m = 2` probe via
`contains_key_hash`, recover summands with `decompose_fast_m2_known`
(known key + hash + reusable pair buffer) instead of
`decompose_fast` (which re-ran `key_of` / hash / admit).

**Hypothesis.** Hits are rare, but each hit was paying a full second
canonicalisation on folded tables; removing it shortens time-to-first
relation slightly without changing answers.

**Verify.** `decompose_fast_m2_known_matches_decompose_fast` on n∈{9,13}.

## Measure

| artifact | result |
|---|---|
| unit test | ok |
| e2e multi N=16 (`logs/cycle2-e2e.txt`) | n=9: **4×** vs cold×16; n=13: **8×** vs cold×16 |
| e2e binary n=9 | root mean 544 µs (CPU); F4 noisy at tiny n |

## Outcome

**Win (modest on cold, strong on multi).** Hit path no longer re-keys.
Multi-target amortization scales with N as expected. Classification:
**engineering**. Cold vs rho unchanged.
