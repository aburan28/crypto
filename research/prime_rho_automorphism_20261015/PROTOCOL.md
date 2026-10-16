# Protocol: rho on the full automorphism group of the round-2 crater curves

**Date.** 2026-10-15. Frozen before any session ran.

## Question

Round 2 ([`../isogeny_conductor_gap_round2_20261011/`](../isogeny_conductor_gap_round2_20261011/README.md),
finding 4) folded the crater's order-6 and order-4 automorphism groups into
index calculus, but its matched rho was `rho.negation` (A = 2) only:
"Rho gains √3 from the same automorphisms, which this repository does not yet
run on prime curves." This round measures that rho.

## Arms and workloads

| spec | reference | candidate | curves |
|:--|:--|:--|:--|
| `specs/cubic.json` | `rho.negation` | `rho.cubic` (A = 6) | the six j = 0 crater curves of round-2 E2 |
| `specs/quartic.json` | `rho.negation` | `rho.quartic` (A = 4) | the six j = 1728 crater curves of round-2 E2 |

Curves, target seed (20261011), four targets per curve, three alternating
rounds, warm-up, seed 1, 900 s timeout and L2 isolation are copied unchanged
from `e2_glv_crater_floor.json`. All twelve curves satisfy the methods'
preconditions (p ≡ 1 mod 3 and 3 | r − 1, or p ≡ 1 mod 4 and 4 | r − 1).

## Prediction

Expected steps scale as √(πr / 2A). The step ratio against `rho.negation`
should be **√3 ≈ 1.732** for `rho.cubic` and **√2 ≈ 1.414** for
`rho.quartic`, before any per-step overhead for computing the class
representative. Under S (group-addition equivalents, the unit of round 2),
the gain is expected to be somewhat below these values.

## Reading

- Report the paired ratio (`ecbench compare`) with its interval, per spec.
- Recompute round 2's "best IC cell vs matched rho" ratio (12.8×) against the
  automorphism rho, for the crater rows only.
- Every collision must be verified (`ecbench verify --replay 12`). A failed
  replay voids the session.
- This is a rho baseline on synthetic 22–28-bit curves. It makes no claim
  about deployed curves.

## Host

One AMD EPYC 9R05 VM, 2 logical cores, 7 GiB RAM, no GPU. This differs from
round 2's host, so absolute times are not compared across rounds. Ratios are
compared within this round only.
