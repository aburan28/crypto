# Two cheap closures: torsion/model invariants and additive energy

**Stage diagnostic.** `S`, end-to-end cost and speedup are unset. No work below
`2^61` is claimed. Both ideas are closed: C1 and C2 passed as frozen in
`PROTOCOL.md`, and the decisive "idea survives" outcome did not occur
(`results/score.json`).

## Answer

1. **Torsion/model invariants give at most 2 bits per summand (Idea 2).** On
   `y² + xy = x³ + 1`, the invariants used by translation-invariant factor bases
   are pullbacks by the separable endomorphisms `1 + τ`, `1 − τ` and `[2]`:
   - `x((1+τ)P) = x + 1/x`;
   - `x(P + T₂) = 1/x`;
   - `(1−τ)(P+T) = (1−τ)P` for `T ∈ E(F₂)`;
   - `λ² + λ = x(2P)`.

   All four held on 100% of points at `n = 7, 11, 13`, which is 10,238 points
   and 71,666 checks. Each kernel has at most 4 points, so replacing `x` by any
   of these invariants changes the decomposition problem by at most
   `log₂ 4 = 2` bits per summand. That cannot move a gap of 42–55 bits.
2. **A geometric factor base has no excess additive energy (Idea 6).** Its
   nontrivial energy matches a size-matched random point set: the pooled ratio
   `R` was 0.94–1.01 in 6 of 7 cells and 1.34 in one, with no growth in `n`.
   There are no free 4-term relations, and the gap formulas stand.

## Results (pooled over 3 seeds per cell; seed `20261007`)

| n | l | R = geometric / random set | geometric / random-V | geometric / `m⁴/#E` | random set / `m⁴/#E` |
|---:|---:|---:|---:|---:|---:|
| 13 | 6 | 0.936 | 1.984 | 0.886 | 0.947 |
| 17 | 7 | 0.944 | 1.121 | 1.003 | 1.062 |
| 17 | 8 | 1.011 | 1.006 | 0.997 | 0.986 |
| 19 | 8 | 0.983 | 1.157 | 0.965 | 0.981 |
| 19 | 9 | 1.339 | 1.789 | 1.329 | 0.993 |
| 23 | 10 | 1.001 | 1.060 | 0.990 | 0.989 |
| 23 | 11 | 0.999 | 1.057 | 0.998 | 0.999 |

- **C1 (identities):** passed.
- **C2 (`R ∈ [0.67, 1.5]` in every cell):** passed.
- **Decisive (`R ≥ 1.5` everywhere, growing with `n`):** not met.

Two cells need comment:

- **`(19, 9)`, `R = 1.34`.** This comes from one seed. Geometric seed 0 has
  nontrivial energy 1.97× expectation, while the other two seeds give 0.99 and
  1.02. The cell stays inside the frozen band, and the effect does not recur at
  `n = 23`, where both cells are within 0.2% of 1. It was not investigated
  further: a single draw of `(θ, g)` that is about 2× expectation is the expected
  tail for these counts. It is not a law that would give free relations.
- **The random-V column at `(13, 6)`.** This reflects the random-V control's
  small, noisy sets (seed 1 had 54 points and a ratio of 0.34). It does not
  reflect the geometric arm, which sits at 0.89× expectation.

## Run and disclosure

- The first launch of the frozen command (`python3 closures.py <dir> 20261007`)
  was stopped by the session's background time limit during the last cell, after
  58 of 66 log lines. Nothing from it is cited.
- The same command was rerun unchanged, detached, in about 2.8 hours. Its first
  58 lines are identical to the stopped run's log, which shows the run is
  deterministic. `results/` holds the complete rerun output, copied unchanged.
- Smoke tests, the pre-freeze metric correction and the dropped `l = 5` cells are
  disclosed in `PROTOCOL.md`.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md` | frozen before the run: C1, C2, the decisive outcome and the smoke tests |
| `closures.py` | identity checks (Part I) and additive energy (Part II) |
| `score.py` | scores C1, C2 and the decisive outcome |
| `results/` | `closures.json`, `closures.log`, `score.json` |

```sh
python3 closures.py results 20261007   # ~3 h, dominated by the n = 23 point enumeration
python3 score.py
```

Pure Python 3, no dependencies. It reuses the frozen `lr.py` and `glin.py`.
