# Free-relation census of the K = 600 pair-root tables: results

[PREDICTION.md](PREDICTION.md) was committed (`e51993fc`) before the fixture runs. The
census harness is
[`examples/koblitz_free_relation_census.rs`](../../examples/koblitz_free_relation_census.rs),
and the run used binary sha256 `2e183196…b203`. Records: `runs/census_n{61,73,83}.json` and
`runs/calibration.jsonl`. The run report, graph and PDF are in
[`reports/run-20261011T2236Z/`](reports/run-20261011T2236Z/).

## Verdict: the note's falsifier fails at every rung. The rank stage is not free.

| fixture | `log₂ r` | states | note `n⁴K⁴/r` | classes `K⁴n³/r` | nontrivial rows (measured) | rank mod `r` (measured) | needed |
|:--|--:|--:|--:|--:|--:|--:|--:|
| n = 61, a = 0 (`e94f7ffb…`) | 47.2 | 21,960,000 | 11,016 | 181 | **147** | **49** | 599 |
| n = 73 (`37f22aed…`) | 56.3 | 26,280,000 | 43 | 0.6 | **3** | **1** | 599 |
| n = 83, a = 1 (`d4a4da26…`) | 52.9 | 29,880,000 | 718 | 8.6 | **18** | **6** | 599 |

- **Every collision checks out.** Each one verified as a group identity
  `S₁ = η·π^t S₂` (identity failures: 0). Each kept row verified as `Σ [c]P = O` by
  scalar multiplication (row failures: 0).
- **Structural collisions.** 1,800 = 3K collisions per fixture give rows that vanish mod
  `r`. They are counted as trivial and dropped.
- **Wall time** (logged, not a measurement): 10 s, 36 s and 46 s.

**What the numbers say:**

1. **The note's count is high by a factor of about `n`.** It binned Frobenius classes into
   `r/2` x-values. The subgroup has about `r/(2n)` Frobenius classes of x-values.
   - With the class count, measured and predicted agree within Poisson scatter on the
     fixtures: 147 against 181.
   - They also agree on the unsaturated synthetic bases: 2,698 against 2,847 at n = 37, and
     8,585 against 8,667 at n = 43 (see `PREDICTION.md`).
   - At n = 73 and n = 83 the counts are too small to separate the two formulas from
     scatter. But both are far below the note's figure: 3 against 43, and 18 against 718.
2. **Rank is exactly one third of the rows, at every rung** (147→49, 18→6, 3→1).
   - A 4-term zero-sum `A + B + C + D = O` pairs up three ways:
     `(A+B) = −(C+D)`, `(A+C) = −(B+D)` and `(A+D) = −(B+C)`.
   - Each of the three is a separate collision in the table. So the distinct relations
     number `≈ K⁴n³/(3r)`: about 60 at n = 61, about 3 at n = 83, and about 0.2 at n = 73.
   - The measured distinct relations are 49, 6 and 1.
3. **The ledger's precompute charge stands.** The free relations cover 49/599 of the rank at
   n = 61, and essentially none of it at n = 73 or 83. The rank stage still needs its
   targeted relations. The "online after precompute" rows keep their precompute.
4. **Making it free would need** `K⁴n³/(3r) ≳ K`, i.e. `K ≳ (3r/n³)^{1/3}`:
   - about 1,300 at n = 61, so a table about 5× larger;
   - about 8,600 at n = 73;
   - about 3,600 at n = 83.

   This is the same generic time–memory line the note itself identifies, at a larger `K`.

## Scope

- **Covered.** The three committed `point_defined_factor_base` fixtures. Only pair sums of
  orbit representatives were used, the table `koblitz_orbit_dlp_fast.rs` builds.
- **Independence.** The census is computed from the group law (library point arithmetic),
  not from that example's normal-basis root solver. So it is an independent census of the
  same set. It does not instrument the production table, and it counts unordered states.
- **No claims.** No speedup or exponent claim. This changes no ledger row; it confirms the
  existing precompute charge.

## Graph

[`reports/run-20261011T2236Z/census.svg`](reports/run-20261011T2236Z/census.svg): the two
predictions, the measured rows and the rank per rung, against the needed 599.
