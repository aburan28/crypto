# Protocol: split and hybrid point-decomposition systems against matched rho

Registered 2026-10-11, **before any measured run of the grid below.** The
smoke runs that sized the grid are listed under *Smoke runs* and are not
results. Conductor task `T-1`.

## Question

Can one point decomposition `R = P₁ + … + P_m` be cut into pieces joined by
auxiliary abscissae ("fan out" the summation polynomial) so that each piece
has a lower solving degree, and does that make the decomposition cheaper?
Does guessing `g` Boolean unknowns (Bettale–Faugère–Perret hybrid) do the
same? And what does each variant cost against matched rho on the same curve?

The repository already measured the chained `S₃` system on its own
(`research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md`, Results 1–10:
refutation degree 6–7 at `m = 3`, one higher at `m = 4`). It never measured
the unsplit `S₄` system, a two-way split with an `S₄` link, or the hybrid on
the same targets, and never put them side by side against rho.

## Systems (native, `src/cryptanalysis/pdp_split.rs`)

Curve `K_a / GF(2^n)` (`b = 1`). Factor base: every point whose abscissa lies
in the standard subspace `V = span(1, z, …, z^{ℓ−1})` (ecbench's
`koblitz-standard-subspace:dimension=ℓ`). A plan is a list of block sizes:

| plan | `m` | links | Boolean unknowns |
|---|---|---|---|
| `unsplit` `[3]` | 3 | `S₄(x₁,x₂,x₃,x_R)` | `3ℓ` |
| `split-2+1` `[2,1]` | 3 | `S₃(x₁,x₂,u)`, `S₃(u,x₃,x_R)` | `3ℓ + n` |
| `split-2+2` `[2,2]` (two-way) | 4 | `S₃(x₁,x₂,u)`, `S₄(u,x₃,x₄,x_R)` | `4ℓ + n` |
| `split-2+1+1` `[2,1,1]` (three-way) | 4 | `S₃(x₁,x₂,u₁)`, `S₃(u₁,x₃,u₂)`, `S₃(u₂,x₄,x_R)` | `4ℓ + 2n` |

The unsplit `m = 4` system needs a symbolic binary `S₅`, which the repository
does not have. It is **not measured**, and the code refuses it rather than
substituting anything. At `m = 3`, the two-way split and the full chain are
the same system.

**Hybrid.** For `g ∈ {0, 2, 4}`, fix the first `g` unknowns (bits of `x₁`) to
each of `2^g` values. Substitute, adjoin the linear pins, and solve every
branch. All branches are charged. Nothing stops at the first solution.

**Solver.** `pq_f4_f2::groebner_basis_f4` (F4 over `F₂[v]/(v²−v)`).

- *Solving degree* means its `F4Stats::solving_degree`: the highest step
  degree at which the basis gained an element. This is **not** the Macaulay
  refutation degree of `RESEARCH_DREG_MEASUREMENT.md`, and the two are not
  compared numerically.
- *Cost* is `F4Stats::word_xors` (64-bit XORs, deterministic) and F4 wall
  time.
- Per hybrid run, the maximum solving degree over branches and the sums of
  XORs and wall time are reported.

## Correctness gates (each run)

1. Every Boolean solution is lifted to curve points. Some signed sum must
   have abscissa `x_R`. Solutions that do not lift, or false ones, are
   counted and kept.
2. The verified decomposition set must equal an **exhaustive enumeration**
   of all multisets of `m` liftable abscissae in `V` summing to `±R`
   (`complete_and_exact`).
3. A run with a timed-out branch or an unresolved basis has
   `complete_and_exact = null` and is reported as censored, never dropped.

## Frozen grid

Seed `20261011`. Four targets per cell (subgroup points `[k]G`, `k` from the
cell's own seeded stream). `--guesses 0,2,4`. Budget 300 s per branch,
`--max-free 24`. 256 extra targets per cell estimate the yield `λ`
(decompositions per target) by exhaustive enumeration, independently of the
solver. Each cell runs as its own process under `timeout 7200`. A killed cell
keeps its flushed lines and is reported as censored.

| cell `m:a:n:ℓ` | curve | surplus `n − mℓ` |
|---|---|--:|
| `3:1:11:4` | `icv1-f2m11-t67-05f5aa36` | −1 |
| `3:0:13:4` | `icv1-f2m13-t181-515ee569` | +1 |
| `3:1:17:4` | `icv1-f2m17-tm101-00378d4e` | +5 |
| `3:1:19:4` | `icv1-f2m19-tm797-9c54981b` | +7 |
| `3:0:13:5` | `icv1-f2m13-t181-515ee569` | −2 |
| `3:1:17:5` | `icv1-f2m17-tm101-00378d4e` | +2 |
| `4:1:11:3` | `icv1-f2m11-t67-05f5aa36` | −1 |
| `4:0:13:3` | `icv1-f2m13-t181-515ee569` | +1 |

Command, per cell:

```sh
timeout 7200 target/release/pdp_split_dreg --cells CELL --guesses 0,2,4 \
  --targets 4 --seed 20261011 --budget-secs 300 --max-free 24 \
  --lambda-targets 256 --out research/pdp_split_dreg_20261011/runs/CELL.jsonl
```

## Matched rho (ecbench)

`SPEC-rho.json` (spec `ECS1h10c3ed7d3628`) covers the same four curves:

- eight public one-target workloads per curve, five rounds;
- `rho.signed_frobenius_strong` (the f6 ladder's frozen parameters) as
  reference, and `rho.signed_frobenius` as baseline;
- `ecbench run`, then `ecbench verify --replay-all`, then `ecbench table`.

The unit is `S = charged group-addition equivalents / √r`.

## Pricing the PDP variants in `S` (derived; a projection)

F4 work is not a group operation, so ecbench leaves it unpriced (README §7).
To put the variants in `S`, this study prices word XORs with the pinned ratio
`ns_per_word_xor` of `docs/ic/calibration.json` for that curve. Those ratios
were calibrated on another host, as the f6 ladder's "pinned-ratio
sensitivity" row was. Per variant and cell:

```
relations needed   ≈ F + 1                     (F = liftable abscissae in V, up to sign)
PDP calls needed   ≈ (F + 1) / λ               (λ from the 256-target exhaustive count)
S_relations        ≈ (F + 1) / λ · X̄ · ratio / √r
```

`X̄` is the variant's mean word XORs per target, all branches included.

This is the **relation phase only**. It omits linear algebra, descent,
factor-base construction and verification, so it is a lower bound on full
IC cost. It is reported against the strong-rho `S` on the same curve and
labelled a projection, never a measured end-to-end `S`. A cell where any
of the four targets is censored gets no projection.

## Predictions (written first)

- **P1.** At every `m = 3` cell, `split-2+1` has a lower solving degree than
  `unsplit`.
- **P2.** At every `m = 3` cell, `split-2+1` costs more than `unsplit` in word
  XORs at `g = 0`, and the ratio grows with `n` at fixed `ℓ = 4`. The cut
  adds `n` unknowns.
- **P3.** At both `m = 4` cells, `split-2+1+1` (three-way) has a lower
  solving degree than `split-2+2` (two-way).
- **P4.** Guessing lowers the per-branch solving degree. For `unsplit`,
  total XORs at `g = 2` or `g = 4` are below `g = 0`. For `split-2+1`,
  `g = 4` is above `g = 0`.
- **P5.** Every complete run passes gate 2.
- **P6.** No variant's projected `S_relations` is below the strong-rho `S`
  on the same curve.

**Falsification target:** any complete cell with a variant whose projected
`S_relations` is below the same curve's strong-rho `S`.

## Inadmissible

The following are inadmissible after any measured run:

- changing the grid, seeds, budgets, solver or pricing ratio;
- dropping censored, failed or regressing runs;
- quoting the projection as a measured `S`, or quoting wall time beyond the
  isolation level ecbench earned;
- any claim about deployed curves, asymptotic sub-rho behaviour, or
  ECC2K-130.

## Smoke runs (sizing only, not results)

Smoke runs on seed 1, kept in the scratchpad and not committed:

- `3:1:11:4` and `4:1:11:3` at `g ∈ {0,2}`;
- `3:0:13:4` and `3:1:17:6`, which showed that `ℓ = 6` at `n = 17` exceeds
  the budget (unsplit times out at 120 s, and the 35-unknown chain overran
  its budget inside one F4 step);
- `3:1:17:4` and `3:1:19:4`.

They fixed `max_free = 24`: at 16, one `[2,2]` basis with no solutions read
as unresolved. They also capped `ℓ` at 5 for `m = 3`.
