# Protocol: round 3 — the paired design, the √2 rung, a 30× gap, and the cubic fold

Frozen 2026-10-12 UTC, before any search or session of this round ran.
Rounds 1 and 2 are
[`../drastic_prime_conductor_gap_20261012/`](../drastic_prime_conductor_gap_20261012/README.md)
and
[`../drastic_prime_conductor_gap_round2_20261012/`](../drastic_prime_conductor_gap_round2_20261012/README.md).

**Result class.** Boundary control and counted-work accounting, as both
prior rounds. No wall-clock speedup, general hardness, or
deployed-curve claim follows.

## E1. The paired planted-scalar design (higher confidence)

Round 2's quotients carry the full between-target variance of rho: the
two curves drew unrelated targets. This round adds
`target_kind: planted_shared`: the scalar `k` derives from
`[r, target_seed, index]` alone (law `uniform_scalar_sha256_shared_v1`,
no curve-model hash), so both curves of a pair — same `r` — solve the
**same** scalars. The paired statistic divides per-index costs and
resamples indices jointly, which removes the between-target variance
from the quotient entirely. It is a different statistic from the
independent law, reported beside it, never replacing it.

The two D = −3 pairs of round 2 re-measure under this design with
`targets_per_curve = 64`: rho-neg on both curves, rho-cubic on the
crater (its own spec, same seed, so the crater solves the same `k`s in
both specs). Predictions as frozen in round 2 (√3 best-method, 1
same-method), now with the tighter intervals the paired design buys.

## E2. The √2 rung: D = −4 (j = 1728) and `rho.quartic`

`QuarticClasses` (A = 4): the orbit of the rational order-4 automorphism
`ι: (x, y) ↦ (−x, iy)` of a `j = 1728` curve `y² = x³ + ax` over
`p ≡ 1 (mod 4)` with `4 | r − 1`, which already contains negation.  A
2-adic obstruction forces this family to cofactor 2: with `t ≡ 0
(mod 4)` (integrality of `q = (t/2)² + ℓ²`) and `q ≡ 1 (mod 4)`, the
order `N = q + 1 − t` is `2 (mod 4)`, so the pair runs on the prime
subgroup `r = N/2` with `r ≡ 1 (mod 4)` (`N ≡ 2 (mod 8)`).

- **P-E2a (the rung).**  Best-method quotient
  `S(rho.negation, floor) / S(rho.quartic, crater) = √2 ≈ 1.4142`,
  99% interval excluding 1.  A value materially below 1.27 would mean
  `rho.quartic` pays more than its one counted rotation per
  canonicalisation.
- **P-E2b (control).**  Same-method rho floor/crater `= 1`.
- Search: `ℓ` prime near `2^20` (±3%), `t ≡ 0 (mod 4)`,
  `3 ≤ t ≤ 3999`, `q = (t/2)² + ℓ²` prime `≡ 1 (mod 4)`, `q < 2^45`,
  `N = q + 1 − t ≡ 2 (mod 8)`, `r = N/2` prime `≡ 1 (mod 4)`,
  `q^i ≢ 1 (mod r)` for `i ≤ 32`; crater by quartic-twist scan
  (order `N`, verified by `[N]P = O` with `[r]P ≠ O` and `y ≠ 0`, so
  the point has order `2r` and `[2]P` generates the subgroup); floor
  random with `j ∉ {0, 1728}`.  Failures kept.

## E3. The 30× gap: D = −19 at ℓ ≈ 2^23

One pair at the round-1 construction (`q = (t² + 19ℓ²)/4`, cofactor 1)
with `ℓ` prime near `2^22` — roughly 2× round 1's largest gap and 32×
the 131101 the programme started from.  (Frozen amendment before the
search ran: the `2^23` window of the pre-registered draft would have
cost hours of floor scanning per candidate on this host; `2^22` keeps
the 30× reach at a measurable budget.)  rho-neg both curves,
planted-shared targets, 32 per curve.  **P-E3:** same-method
floor/crater `= 1` (the generic theorem; a deviation is a harness
fault, and this cell exists to show the theorem holding as the gap
grows, not to re-test it).

## E4 (deferred lead, frozen for the next session): the cubic fold of the MITM table

Round 2's √3 is the generic ceiling.  The algebraic lever — a
`PairTable` keyed by the canonical `ζ₃`-orbit of the pair sum, one
entry per 6-orbit instead of per ±-orbit, with multiplier-carrying
relation rows — is frozen here as the round-4 lead: at matched raw
base budget it should cut the crater's table build about ninefold
against the floor's negation fold and change the relation-row algebra
from ±1 entries to `Z_r` coefficients.  Implementation and
measurement are the next session's work; nothing in E1–E3 depends on
it.

## Measurement and analysis

All sessions: `rounds = 1`, `warmup = 1`, order alternate, isolation L2
requested, timeout 600–900 s, ecbench with exact replay audits.  The
analyzer computes, per pair: the **paired** quotient (per-index
differences of per-index means, joint bootstrap over indices) and the
independent-law quotient of the prior rounds, both at 20 000
resamples, frozen seed 20261012, plus A/A exact-count checks and
descriptive wall means.

## Decision rule (frozen)

A rung is a result only if its 99% interval excludes 1 with the same
sign as its prediction on every pair of its family, and every A/A
interval includes 0.  E3's quotient must include 1.  No frontier,
leaderboard row or hypothesis status is promoted.

## What is not measured here

The degree-`ℓ` isogenies; index calculus (E4 deferred); the
root-finding PDP of #1729 (an example harness, not yet an ecbench arm
— wiring it is recorded as the step-5 follow-on); Tier-1 Linux wall
times (this host records at L0; the crypto repository's
`tools/isolated_bench.py` needs a Linux host with PSI and `taskset`,
not available in this environment — recorded as the step-3 blocker);
external-host replay (same blocker).

## Amendments

None so far.
