# Pre-registration: a triple-sum table for the tournament's IC arm

Written and committed **before the collector exists and before any
measurement**, per `AGENTS.md` §4.  Nothing below is a result.

## Why this collector

The tournament's IC arm (`koblitz_tiny_ic.rs`, round-0020 `both` +
`round21-wide-pair-table.patch` + `round23-scaled-base.patch`, i.e. `scaled`)
builds a Frobenius-canonicalised **pair** table from `t` orbit representatives
and, for each walk target, scans the base and looks `target − q` up in it.  Its
own model states the boundary and scopes it:

> At its own optimal factor-base size this method is `r^(2/3)` where rho is
> `r^(1/2)` … It is a statement about this pair-table collector, not about index
> calculus.  (`campaign_20260916/round19_model.py`)

Tabulating `v`-subset sums and scanning one summand per target costs
`#E^{v/(2v−1)}` at the balanced optimum, so IC/rho grows as
`r^{1/(2(2v−1))}`: `r^{1/6}` for pairs (`v = 2`, measured `r^[0.13, 0.19]`),
**`r^{1/10}` for triples**.  This is the first lever in this campaign aimed at
the rate rather than a constant.  `RESEARCH_ECC2K130_DECOMPOSITION.md` §6.2
derives the same law from the other side.

## The prediction, in the solver's own units

`model.py` mirrors the solver's published cost model exactly —
`W(t) = t·size + (K+3)/(λ·cov(t))`, minimised over `t` at run time — and writes
the triple-table analogue with the same units:

```text
W₃(t) = t·size(size+1)/2 + (K+3)/(λ₃·cov₃(t)),   λ₃ = C(size+2, 3)/r,   cov₃ = 1 − ((K−t)/K)³
```

**Calibration first.**  `W₂` puts its own optimum at 6, 12 and 26 orbits
against measured best bases of 8, 16 and 24, and its value at the measured base
is within 3.5% of its own minimum at every cell (8 is the batch floor).  The
continuous family model, by contrast, predicts 30 and 83.  So `W` is the model
to predict with.

| cell | `scaled` base | `W₃*/W₂` | `K₃*`, `t₃` | whole job, triple/`scaled` | IC/rho, measured → predicted |
|---|--:|--:|--:|--:|--:|
| `n23a1` | 8 orbits | 0.924 | 1, 1 | **0.938–0.951** | 0.831 → 0.78–0.79 |
| `n37a0` | 16 | 0.773 | 2, 1 | **0.816–0.852** | 1.396 → 1.14–1.19 |
| `n43a1` | 24 | **0.460** | 2, 1 | **0.563–0.649** | 2.207 → 1.24–1.43 |

The whole-job ratio applies `W₃/W₂` to the collection share only, taken as
0.65–0.81 from round 0023's phase profile, and holds every other phase fixed.
That is conservative: a one- or two-orbit base is also cheaper to build and to
solve.  IC/rho is round 0022's merged measurement times that ratio.

## What must be true for the prediction to hold

1. **A triple-sum entry costs what a pair-sum entry costs**: one addition, one
   canonicalisation, one insert.  This is the assumption most likely to fail,
   and `n37a0` is the cell most exposed to it (its `W₃` is 82% table).
2. **The base builder can make one or two orbits.**  Today it samples in
   batches of eight and cannot go below one batch.
3. **Four-term relations verify.**  `oracle.py` must accept a witness of four
   base points.

## Decision rule, fixed now

Measured as paired instruction counts, triple against `scaled`, on identical
fixtures and seeds, every report oracle-verified.

- **Correctness gate, absolute:** every triple-arm report verifies and recovers
  the same logarithm as `scaled`.  One failure is a defect, not a measurement.
- **Supported:** the whole-job ratio at `n43a1` falls inside or below
  `[0.563, 0.649]`, and the ratio decreases from `n23a1` to `n37a0` to `n43a1`,
  which is the signature of a rate change rather than a constant.
- **Refuted:** the ratio at `n43a1` is `0.85` or above, or the gain does not grow
  with `r`.  Refuted means assumption 1 failed, and the note will say by how much.
- **Otherwise:** partial, reported as that.

A whole-job ratio below one at the crossover cells would still not beat rho
there: the best predicted cell, `n43a1`, stays at 1.24–1.43 × rho.  What this
would move is the crossover and the rate, not the verdict at 131 bits.  This is
an **engineering** candidate for the tournament, not a claim against rho, and it
is offered to that lane's harness rather than scored outside it.
