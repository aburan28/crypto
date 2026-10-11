# Graph-cycle index calculus with the linear (A, P) pair oracle

This is a **full ECDLP measurement on toy Koblitz curves**. Every run solves the
discrete log end to end and verifies `kP = Q`. **Class: boundary (negative
result).** The method works, but it is dominated by Pollard rho, both at the
sizes measured and asymptotically. No claim is made at `n = 131` beyond the
arithmetic below.

## Answer

Dropping linear algebra lifts the factor-base cap from `l ≈ 30` to the oracle's
reach, `l ≤ (n+1)/3`, which is 44 at `n = 131`. But the method then needs about
`2^{n−l}` oracle tests, and that is more than `√r` whenever `l < n/2`. So:

1. **Its own walk beats it.** A walk long enough to collect the edges collides
   first, and the collision is a rho solve. Measured: rho solved first in **20 of
   20** walk-mode runs.
2. **Without a walk it is still slower.** Driven by fresh independent points, with
   each oracle solve charged as just one group operation, its median `S` is
   **8.9–35.9** against rho's **1.0–2.0**.
3. **At `n = 131`:** about `2^{n−l−1.3} ≈ 2^{85.7}` oracle solves against
   `2^{60.81}` for rho, with about `2^{43}` stored nodes. Asymptotically it costs
   `2^{2n/3}` against `2^{n/2}`.

## Method (`cyc.py`)

The curve is `y² + xy = x³ + 1` over `F_{2^n}` with `#E = 4r`, `r` prime. That
holds for `n = 13, 19, 23` in the toy range, found by the Frobenius-trace
recurrence. `P` generates the order-`r` subgroup and `Q = kP`.

1. Each test is one `n × (3l − 1)` linear `(A, P)` solve over a geometric `V`, as
   in [`../geometric_v_linear_20261006`](../geometric_v_linear_20261006/README.md).
2. A success `W = aP + bQ = s₂A₂ + s₃A₃` is an edge with the label
   `s₂λ₂ + s₃λ₃ = 4(a + bk) mod r`, where `λ = log_P(4A)`.
3. A union-find tracks each node's value as `σ·root + α + βk` over `Z/r`.
   - An even cycle gives `α + βk = 0`.
   - An odd cycle fixes the root, and a second constraint on that component
     gives `k`.
4. There is no linear algebra.

The reference is a plain Pollard rho on the same instance: an r-adding walk,
collisions found by storing visited points, no negation or Frobenius speedup.
`S = operations / √r`.

## Results

The protocol is `PROTOCOL.md`, frozen before the run. Seeds were `2026100719`
(`n = 19`, 6 instances) and `2026100723` (`n = 23`, 4 instances). All of Y1–Y5
passed (`results/score.json`).

| n | l | walk: rho solved first | fresh: solved | median edges / node | median `log₂ tests − (n−l)` | median `S_tests`, fresh | median `S`, rho |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 19 | 5 | 6/6 | 6/6 | 0.565 | −0.36 | 35.87 | 1.00 |
| 19 | 6 | 6/6 | 6/6 | 0.418 | −1.31 | 9.41 | 1.00 |
| 23 | 7 | 4/4 | 4/4 | 0.447 | −1.23 | 19.61 | 1.96 |
| 23 | 8 | 4/4 | 4/4 | 0.366 | −1.38 | 8.92 | 1.96 |

- **Y1 (correctness):** every solve, by rho, by cycle, or in either mode, verified
  `kP = Q`.
- **Y2 (domination):** in walk mode, rho solved first in 20 of 20 runs.
- **Y3 (graph law):** about 0.4–0.6 edges per node close a usable cycle,
  consistent with a random graph `G(N, M)` at `M = Θ(N)`.
- **Y4 (test law):** tests were within 1.4 bits below `2^{n−l}`.
- **Y5 (cost):** fresh-mode `S_tests` exceeds rho's `S` in every cell.

That holds under two charges that favour the cycle method: each oracle solve
(a linear solve of about 80 unknowns) is charged as one group operation, and the
fresh mode's two scalar multiplications per test are not charged at all. Those
uncharged scalar multiplications are counted in `results/n19.json` and
`results/n23.json`.

## Reading

This closes the graph-cycle variant of the outer-algorithm branch of the kill test from #1442.

Removing linear algebra removes the `c ≤ 30` cap, but it exposes a harder one:
a method that touches about `2^{n−l}` group elements cannot beat a collision
search that needs about `2^{n/2}` of them, unless `l > n/2`. The pair oracle's
reach, `(n+1)/3`, is below that. The cycle method would need a pair test reaching
`l ≥ n/2 + 1`, and there the linear `(A, P)` system has about `3n/2` unknowns in
`n` equations.

The smoke tests (`n = 13, 19`, disclosed in `PROTOCOL.md`, not cited) are what
first exposed the walk-domination point. They also caught a sign error in the
union-find merge branch, which was fixed before the protocol was frozen.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md` | frozen before the run: Y1–Y5 and the disclosed smoke tests |
| `cyc.py` | graph-cycle IC, walk and fresh modes, and the rho reference |
| `score.py` | scores Y1–Y5 |
| `results/` | `n19.json`, `n23.json` and their logs, plus `score.json` |

```sh
python3 cyc.py 19 6 2026100719 5,6 results/n19.json
python3 cyc.py 23 4 2026100723 7,8 results/n23.json
python3 score.py
```

The runs wrote into the session scratchpad, and their outputs were copied into
`results/` unchanged. Pure Python 3, no dependencies. It reuses the frozen
`lr.py` and `glin.py`.
