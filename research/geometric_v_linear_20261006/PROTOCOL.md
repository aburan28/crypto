# Protocol: the linear (A, P) oracle over a geometric factor base

Frozen 2026-10-06, before the confirmatory run. **Stage diagnostic.** `S`,
end-to-end cost and speedup are **unset**. No work below `2^61` is claimed.
This corrects the frontier stated in
[`../linearization_reach_20260930`](../linearization_reach_20260930/README.md)
and [`../xl_bilinear_core_20261006`](../xl_bilinear_core_20261006/README.md).

## Derivation (stated before measuring)

On `y² + xy = x³ + 1`, `S₃ = e₂² + e₃ + 1`. Fix `x₁ = t` and put `A = X₂ + X₃`
and `P = X₂X₃`. Then:

    S₃(t, X₂, X₃) = P² + t²A² + tP + 1.

This is `F₂`-linear in the bits of `(A, P)`, because squaring is `F₂`-linear.
`X₂` and `X₃` are the roots of `Z² + AZ + P`.

For `X₂, X₃ ∈ V`:

- `A ∈ V`, which is `l` bits.
- `P ∈ span(V·V)`. That span has dimension at most `2l − 1` when
  `V = θ⟨1, g, …, g^{l−1}⟩` is geometric, and `l(l+1)/2` when `V` is random.

So one trial is a single linear solve with `N` unknowns:

- geometric `V`: `N = 3l − 1`, determined while `3l − 1 ≤ n`;
- random `V`: `N = l(l+3)/2`, the same count as the earlier threads at `b = l`.

The earlier threads used random `V` and coordinates over `W < V`. There the reach
ended near `l = 14` at `n = 131`. With geometric `V`, both summands are resolved
(`b = l`) up to `l = ⌊(n+1)/3⌋`, which is 44 at `n = 131`.

The repository already measured `dim span(V·V) ≤ 2l − 1` and used it inside a
Gröbner solve
([`../notes/koblitz-isogeny/subspace-structure-20261004.md`](../notes/koblitz-isogeny/subspace-structure-20261004.md)).
What is tested here is the purely linear `(A, P)` solve, and what it does to the
chain cost map. The derivation was proposed by an idea-generator subagent and
checked by hand before freezing.

**Law.** A random target abscissa `t` is hit by a pair from `V` with probability
about `2^{2l − n − 1}`: there are about `2^{2l−2}` liftable unordered pairs, each
giving 2 abscissae, over about `2^{n−1}` point abscissae. This holds for either
`V`. The two bases differ only in the solve: past the reach, the solution family
has `2^{N − rank}` points to enumerate.

## Instrument

`glin.py` reuses the frozen `lr.py`. Each trial does the following:

1. Draw a fresh `V` (geometric: random `θ`, `g`).
2. Draw a random curve point `T`.
3. Solve the `n × N` system.
4. Enumerate the solution family, up to `2^6` points; a larger family counts as
   `too_big` and is unresolved.
5. Accept `(A, P)` only if:
   - both roots of `Z² + AZ + P` lie in `V`;
   - both lift to curve points;
   - `±P₂ ± P₃ = ±T` holds on the curve.

Everything else counts as `fake`. Planted cells use 30 pairs per cell, with
`T = P₂ + P₃`.

## Smoke test (disclosed, not cited)

Seed 1: `n = 19` with `l = 5, 6, 7`, and `n = 23` with `l = 8`, 2,000–4,000
trials each.

- Success rates were within 0.6 bits of `2l − n − 1` for geometric `V`, and
  within 0.5 bits for random `V` at `n = 19, l = 6`.
- Planted recovery was 50 of 50 in every cell.

## Predictions (pass/fail)

Grid: `n ∈ {19, 23, 29, 31}`, `l = 4 … ⌊(n+1)/3⌋ + 2`, both bases. Trials are
sized for about 24 expected successes, capped at `2^15` for geometric and `2^13`
for random. Seed `20261008`.

- **G1 (law).** Take the geometric cells within reach (`3l − 1 ≤ n`) that have
  at least 8 successes. In each, `|log₂ p̂ − (2l − n − 1)| ≤ 1.5`. Per `n`, with
  at least 3 cells, the least-squares slope of `log₂ p̂` on `l` lies in
  `[1.6, 2.4]`.
- **G2 (reach).**
  - Geometric cells within reach have a mean rank defect `N − rank` of at
    most 2.0.
  - Every cell, of either base, with `N > n` has a mean defect within
    `[0, +1.5]` of `N − n`. That is, both bases sit at generic rank, so the
    family size is `2^{N − n}` and is set by `N`.
- **G3 (recovery).** Geometric cells within reach recover the planted pair in
  at least 29 of 30 instances.

## n = 131 arithmetic (`gcost.py`)

The heuristics and budget convention are the same as `costmap_a1.py`. Per
relation the oracle needs `2^{n − 2c}` trials. Each trial is charged
`ω · log₂(3c − 1)` bits, plus the family size past the reach. With
RHO = 60.8090, the log₂ of matched rho:

| | bits-only gap | charged gap, ω = 2 | charged gap, ω = 2.81 |
|---|---:|---:|---:|
| previous frontier (random `V`, `l = b = 14`) | 56.19 | 71.82 | 78.15 |
| geometric `V`, subspace, `m = 5`, `c = 28` | **42.47** | **55.22** | **60.39** |

The orbit-union rows agree within 0.1 bit (`m = 4`, `c = 28`: 42.41 / 55.16).
They are not headlined, because an orbit-union trial must place `X₂` and `X₃`
in one Frobenius conjugate for `span(V·V)` to stay small, and that case is not
measured here.

Both corrected gaps remain positive, so nothing here comes in under rho. If G1–G3
pass, these figures replace 56.19 and 71.82 as this thread's frontier.

## Stop condition

The run is bounded. If G1 or G3 fails, the derivation is wrong and the earlier
frontier stands. If G2 fails, the reach is not set by `N`, and the n = 131 table
must be recomputed from measured families.

Command: `python3 run_confirm.py results 20261008 && python3 score_confirm.py`.
