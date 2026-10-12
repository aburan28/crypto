# Protocol: drastic prime conductor gaps, ECDLP measured on both ends

Frozen 2026-10-12 UTC, before any search ran and before any session was
sealed.  This extends
[`../ecbench_conductor_gap_20261012/`](../ecbench_conductor_gap_20261012/README.md),
which measured one-target ECDLP across the archived `ℓ = 8209` and
`ℓ = 131101` pairs (`D = −19`), to **more drastic single-prime gaps** and a
**second discriminant**.

**Result class.**  Boundary control and counted-work accounting, like the
prior round: a quotient of counted lower bounds is a diagnostic, not a
bound on complete attack cost, and no wall-clock speedup, general hardness,
or deployed-curve claim follows.

## Question

Two curves over the same prime field `GF(q)` with the same prime order
`N = q + 1 − t`, the same Frobenius discriminant
`t² − 4q = −D·ℓ²` (`D ∈ {19, 43}`, both fundamental and class number one),
differing in exactly one invariant: the crater curve has
`End(E) = O_K` (conductor `f_E = 1`) and the floor curve has `f_E = ℓ`,
where `ℓ` is a prime of **19 to 21 bits** — a single-prime conductor gap
between `2^19 ≈ 5.2·10⁵` and `2^21 ≈ 2.1·10⁶`, i.e. 4 to 16 times the
largest gap measured before (131101).  Does the paired counted cost

    S(method) = cost(floor) / cost(crater),  counted group-addition
    equivalents per `√r`,

differ from 1 for any method ecbench runs?

## Why these discriminants

`D = −19` continues the existing line at larger `ℓ`, so the `ℓ`-trend of
the prior round extends.  `D = −43` (`j(O_K) = −884736000`) is a second
class-number-one CM field with no prior conductor-gap measurement: an
independent family, so a repeated lead cannot be an artefact of one CM
field.  Both `D ≡ 1 (mod 4)` are odd, so odd `t` and odd `ℓ` keep
`q = (t² + Dℓ²)/4` integral.

## Pre-registered predictions

- **Generic methods** (`rho.negation`, and its A/A control): `S = 1`
  within the bootstrap interval at every `ℓ` and both `D`.  This is the
  same theorem the prior round cites (generic methods see only the group);
  a deviation is a harness fault, not a discovery.
- **`ic.pipeline` (prime-abscissa base, `mitm` m=2 oracle)**: prior
  rounds found `S ≈ 1` at `ℓ = 263, 8209, 131101`; prediction `S = 1` here
  too.  A confirmed `S ≠ 1` with the same sign on two independent pairs
  would be the first evidence an implemented algebraic method sees the
  conductor at a drastic gap.
- **Search prediction**: in a window of ±3% around `2^k` for
  `k ∈ {19, 20, 21}` and `|t| < 4000`, at least one certified pair per
  `(D, k)` exists with `q, N` both prime and `N` not anomalous; failures
  and empty windows are recorded, never silently narrowed.

## Decision rule (frozen before running)

A method "sees the gap" only if, on **two independent pairs** — different
`ℓ` or different `D` — the paired bootstrap interval of `log₂ S` excludes
0 with the same sign, and the A/A control arm's interval includes 0.  One
pair alone is a lead, not a result.  Intervals include 1 ⇒ the repeated-lead
criterion is not met; this is not an equivalence proof.

## E1.  Search (certified pair construction)

Runner: `examples/drastic_prime_gap_pairs.rs`, deterministic from its
arguments, one JSON output, failures kept.

- For `D ∈ {19, 43}`, `k ∈ {19, 20, 21}`: scan odd `ℓ` in
  `[2^k·0.97, 2^k·1.03]` prime; for odd `t` with `3 ≤ |t| ≤ 3999` (t ≠ 1:
  exclude the anomalous trace) form `q = (t² + Dℓ²)/4`; require `q` prime,
  `N = q + 1 − t` prime (`N ∤ q^i − 1` for `i ≤ 32`: no small embedding
  degree), and `q < 2^46` so u64 arithmetic and rho both stay exact.
- **Crater**: `j₀ = j(O_K) mod q`; curve `y² = x³ + 3c·x + 2c` with
  `c = j₀/(1728 − j₀)`; whichever of it and its quadratic twist has order
  `N` (one scalar multiplication of a random point) is the crater.  If
  neither does, record a failure.
- **Floor**: random `(a, b)`, kept when a random point satisfies
  `[N]P = O` and `j ≠ j₀`.  At most 200 trials per pair; exhausted scans
  are recorded as failures.
- **Certification without Sage** (same class-number-one argument the
  prior round reused): `t² − 4q = −Dℓ²` with `D` fundamental fixes
  `f_π = ℓ`; `j = j(O_K)` with `h(D) = 1` gives `End = O_K`; `j ≠ j(O_K)`
  gives `f_E | ℓ`, `f_E ≠ 1`, hence `f_E = ℓ`.  Both curves' orders are
  verified by scalar multiplication, never assumed.
- Output `pairs/pairs.json`: full model, generator (first
  non-identity point of order `N`), j-invariants, certification notes,
  timings, seeds, and every failure.

## E2.  Paired one-target ECDLP measurement

Specs: one ecbench spec per pair, crater first, `targets_per_curve = 8`,
fresh `target_seed`, arms `rho-neg` (reference), `rho-neg-aa` (A/A
control) on every pair, plus `ic-mitm` (the prior round's exact
configuration: `prime-abscissa:size=2048`, `mitm:m=2,negation_folded=1`,
`structured-gauss`, walk targets) on every pair where it completes inside
the timeout.  `rounds = 1`, `warmup = 1`, order alternate, isolation L2,
timeout 3600 s.  The target law hashes the curve model, so cross-curve
samples are independent and no target is transported; the isogeny itself
is not walked (its cost is not assigned zero, it is out of scope, as in
the prior round).

Runner: `run_sessions.sh`; audit `ecbench verify --replay-all`;
analysis `examples/drastic_gap_analyze.rs` (same statistical core as the
prior round: 20 000 resamples, frozen seed, independent-workload
bootstrap and label permutation; the IC arm is reported per pair and
marked inadmissible when it does not complete).

## What is not measured here

Walking the degree-`ℓ` isogeny (for `ℓ > 10⁶` the modular polynomial is
far out of reach); wall-clock speedup; external-host replay; intrinsic or
cryptographic-scale hardness.  Eight targets on a fixed curve are eight
workloads on one curve, not eight curves.

## Amendments

None so far.
