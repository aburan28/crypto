# Protocol: graph-cycle index calculus with the linear (A, P) pair oracle

Frozen 2026-10-07, before the confirmatory run. **This is a full ECDLP measurement
on toy curves, not a stage diagnostic.** Each run solves a discrete log end to
end and verifies `kP = Q`. At toy sizes it is compared against a plain Pollard
rho on the same instance. No claim is made at `n = 131`.

This is Idea 4 from the screen in
[`../geometric_v_linear_20261006`](../geometric_v_linear_20261006/README.md),
proposed by an idea-generator subagent.

## Method (`cyc.py`)

The curve is `y² + xy = x³ + 1` over `F_{2^n}` with `#E = 4r`, `r` prime:
`n ∈ {19, 23}`, with `r = 130873` and `2095853`. `P` generates the order-`r`
subgroup and `Q = kP`.

1. Each test asks one linear `(A, P)` solve whether a walk point
   `W = aP + bQ` equals `s₂A₂ + s₃A₃`, with `A₂, A₃` lifts of abscissae in a
   geometric `V` of dimension `l`.
2. A success is an edge with the label `s₂λ₂ + s₃λ₃ = 4(a + bk) mod r`, where
   `λ = log_P(4A)`.
3. A union-find tracks each node's value as an affine function `σ·root + α + βk`
   over `Z/r`. The first usable cycle yields `k`: either an even cycle, or two
   constraints on one component.
4. There is no linear algebra.

## Two modes

- **walk:** an r-adding walk drives the tests, and the same walk is watched for a
  rho collision. The run records which of the two solves first.

  *Derivation:* a walk longer than about `√r` collides. The cycle method needs
  about `2^{n−l}` tests, with `l ≤ (n+1)/3 < n/2`. So it should almost always lose
  to the rho embedded in its own walk.
- **fresh:** independent uniform `W` per test, with the scalar multiplications
  counted separately. This isolates the graph law: about `N/2` edges to a usable
  cycle, where `N ≈ 2^{l−1}` nodes. Each test succeeds with probability about
  `2^{2l−n−1}`, so tests come to about `c · 2^{n−l}`.

## Reference

A separate plain Pollard rho on the same instance: an r-adding walk, collisions
found by storing visited points, no negation or Frobenius speedup. Cost is
reported as `S = operations / √r`.

The fresh mode's `S_tests` charges each oracle solve as **one** group operation
and leaves the scalar multiplications uncharged. Both choices favour the cycle
method, so a loss under them is robust.

## Smoke tests (disclosed, not cited)

- **Seed 1, `n = 13`, `l = 3, 4`, walk only.** The graph had 3–11 nodes, and the
  walk cycled before gaining edges. That exposed the domination point and led to
  the two modes above. `n = 13` was dropped as too small.
- **Seed 1, `n = 19`, `l = 5, 6`, 2 instances.** Walk mode: rho solved first in
  4 of 4. Fresh mode: solved in 4 of 4, with 5,256–13,622 tests, 0.43–0.73 edges
  per node, and `k` verified each time. Rho's `S` was 1.03–1.65; fresh `S_tests`
  was 14.5–37.7.
- **Before testing:** a sign error in the union-find merge branch was fixed. It
  occurred when the absorbed root already had a known value.

## Predictions (pass/fail)

Grid: `n = 19` with `l ∈ {5, 6}` (6 instances), and `n = 23` with `l ∈ {7, 8}`
(4 instances). The test cap is `2^{n−l+5}`. Seeds are `2026100719` (`n = 19`)
and `2026100723` (`n = 23`).

- **Y1 (correctness).** Every solve, by any method, verifies `kP = Q`.
- **Y2 (domination).** In walk mode, rho solves first in at least 90% of instances
  in every cell.
- **Y3 (graph law, fresh).** The median edges per node at the solve lies in
  `[0.2, 2]`.
- **Y4 (test law, fresh).** The median of `log₂(tests) − (n − l)` lies in
  `[−2, +3]`.
- **Y5 (cost, fresh).** In every cell, the median `S_tests` exceeds the median `S`
  of the reference rho. This holds even with the oracle charged as one group
  operation.

Commands:

    python3 cyc.py 19 6 2026100719 5,6 results/n19.json
    python3 cyc.py 23 4 2026100723 7,8 results/n23.json
