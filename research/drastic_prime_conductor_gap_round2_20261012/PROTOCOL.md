# Protocol: the j = 0 crater — the biggest provable generic discrepancy across a prime conductor gap

Frozen 2026-10-12 UTC, before any search or session of this round ran.
Round 1 is
[`../drastic_prime_conductor_gap_20261012/`](../drastic_prime_conductor_gap_20261012/README.md):
rho floor/crater quotients sat at 1 on five of six pairs (one 95%-level
lead, 0.528596 at `d43-ell1016839`), and the prior MITM IC configuration
stopped completing at `q ≥ 2^41`.

**Result class.** Boundary control and counted-work accounting, as round
1. A quotient of counted lower bounds is a diagnostic; no wall-clock
speedup, general hardness, or deployed-curve claim follows. The wall-time
figures are descriptive (this host records at L0).

## Question

Round 1's curves agree for every method ecbench ran because both ends of
each pair offer rho the same automorphism group (negation only, `A = 2`).
The one CM family where that is false is `D = −3`: the crater of a
`D = −3` class-number-one volcano is a `j = 0` curve `y² = x³ + b`,
whose rational automorphism group has order **6** (the order-3 map
`ζ₃: (x, y) ↦ (ωx, y)` for a primitive cube root `ω` of `1` in `F_p`,
together with negation), while every floor of the same class has
`j ≠ 0` and only negation.  A rho that walks on canonical orbit
representatives searches `r/A` classes, so the crater's best generic
rho is shorter by `√3` — the **largest generic discrepancy possible
between two isogenous curves**, and this round's target:

    S(rho.negation, floor) / S(rho.cubic, crater) = √3 ≈ 1.7321

at a single-prime conductor gap `ℓ` of 20 to 21 bits (1 to 2 million),
with **32 independent one-target workloads per curve** so the bootstrap
intervals are about twice as tight as round 1's — the higher-confidence
comparison round 1's eight workloads could not deliver.

## New method

`rho.cubic` (`ic_boundary::CubicClasses`, `A = 6`): the tuned walk of
`rho.negation` on the classes `{±ζ₃^i(P)}`.  It exists exactly when
`a = 0` (j = 0), `p ≡ 1 (mod 3)` (so `ω ∈ F_p`), `3 | r − 1` (so `ζ₃`
acts on the order-`r` subgroup as `[κ]` for a primitive cube root `κ` of
`1 mod r`), and `ζ₃(G) ≠ G`; construction verifies `[κ]G = ζ₃(G)` with
two counted scalar multiplications, and every canonicalisation's two
`ω`-rotations are counted field multiplications.  It refuses any curve
without the structure — the floor cannot run it, by design.

## Pre-registered predictions

- **P1 (same-method control).**  `rho.negation` floor/crater `S`
  quotient `= 1`, 99% interval covering 1, both pairs.  Generic
  theorem; a deviation is a harness fault.
- **P2 (the discrepancy).**  Best-method quotient
  `S(rho.negation, floor) / S(rho.cubic, crater)` `= √3 ≈ 1.7321`,
  99% interval excluding 1, same sign on both pairs: the repeated-lead
  rule met with a theorem-backed mechanism.  A value materially below
  `√3` (say below 1.55) would mean `rho.cubic` pays more than its
  rotations' counted cost — an implementation finding, reported as such.
- **P3 (within-curve validation).**  On the crater itself,
  `S(rho.negation) / S(rho.cubic) = √3` within the same intervals.
- **P4 (round-1 lead replication).**  `d43-ell1016839` re-measured with
  32 fresh targets: either the 0.5286 quotient persists (99% interval
  excluding 1, same sign — the lead becomes a result to investigate) or
  it regresses to 1 (round 1's reading was noise).  Either outcome is
  reported, never averaged away.
- **P5 (search).**  Certified `D = −3` pairs with `q, N` prime,
  `q ≡ N ≡ 1 (mod 3)`, exist at `ℓ` near `2^20` and `2^21`; every empty
  window or failed twist scan is recorded as a failure.

## Decision rule (frozen)

P2 is a *result* only if both pairs' 99% intervals exclude 1 with the
same sign and the A/A control arms' intervals include 0 everywhere.
P4's lead is a result only if its 99% interval excludes 1.  Nothing in
this round promotes a frontier, leaderboard row or hypothesis status.

## E1.  Search (certified `D = −3` pairs)

Runner: `examples/drastic_gap_round2_pairs.rs`, deterministic from its
arguments, failures kept.

- `ℓ` prime near `2^k`, `k ∈ {20, 21}`, ±3% window; odd `t`,
  `3 ≤ |t| ≤ 3999`, `t ≠ 1`; `q = (t² + 3ℓ²)/4` prime with `q ≡ 1 (mod 3)`
  and `q < 2^44`; `N = q + 1 − t` prime with `N ≡ 1 (mod 3)`; `q^i ≢ 1
  (mod N)` for `i ≤ 32`.
- **Crater**: scan `b` for `y² = x³ + b` (j = 0 by construction) whose
  order is `N` (64 random `[N]P = O` per candidate); the sextic-twist
  scan runs to 400 candidates under a 900 s bound.
- **Floor**: random `(a, b)` with order `N` and
  `j ∉ {0, 1728}` (1728 excluded so no quartic automorphism hides in
  the floor), same trial accounting as round 1.
- Certification: `t² − 4q = −3ℓ²` with `D_K = −3` fundamental,
  `h(−3) = 1`, `j(O_K) = 0`; crater `j = 0` gives `End = O_K`; floor
  `j ≠ 0` gives `f_E = ℓ`.  The automorphism structure is certified by
  `rho.cubic`'s own construction checks (`p ≡ 1 (mod 3)`,
  `3 | r − 1`, `[κ]G = ζ₃(G)`), which refuse the floor by `j ≠ 0`.

## E2.  Paired measurement

One ecbench spec per pair, crater first, `targets_per_curve = 32`,
fresh `target_seed`.  Arms: `rho-neg` and `rho-neg-aa` on both curves,
`rho-cubic` on the crater only (it refuses the floor).  Plus one
replication spec for `d43-ell1016839` (round 1's exact curves) with
32 fresh targets and the rho arms only.  `rounds = 1`, `warmup = 1`,
order alternate, isolation L2 requested, timeout 600 s.

Analysis `examples/drastic_gap_round2_analyze.rs`: round 1's
statistical core (20 000 resamples, frozen seed 20261012,
independent-workload bootstrap and label permutation) with two
quotients per pair — same-method and best-method — and, new in this
round, descriptive wall-time means with the same bootstrap intervals
(L0 host; never interpreted as a speedup).

## What is not measured here

The degree-`ℓ` isogeny; `D = −4` (a `√2` discrepancy, strictly smaller);
index calculus (round 1 showed the prior configuration exhausts at these
sizes); external-host replay; intrinsic or cryptographic-scale hardness.

## Amendments

None so far.
