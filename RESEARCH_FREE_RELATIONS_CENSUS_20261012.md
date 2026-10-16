# Free relations by pair-root self-collision: measured

**Date:** 2026-10-12. **Status:** measurement closing Lead B of
`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §3. Every number below is a
measurement from `experiments/koblitz-free-relations-20261012/` (receipts
and verified rows in that directory); no number here is a calculation
except where marked. No security claim about a deployed curve is made or
withdrawn.

## What was measured

The Lead B proposal: two pair-states of the compact-orbit S3 root table
that share a root x-coordinate up to Frobenius rotation are a four-term
relation among factor-base points, so sorting the table should fix the
orbit logs without the guided-rank probes. The prediction was ≈11,000
free relations at n=61/K=600 and ≈730 at n=83, with the rank stage "free"
at the rungs the ladder runs.

Two new census examples (`examples/koblitz_free_relations.rs` for n ≤ 63,
`examples/koblitz_free_relations_wide.rs` for 64 < n ≤ 127) enumerate the
same states and roots as the production `build_index`, sort every
`(canonical root, state, rotation)` entry, verify each collision as an
exact group identity, and insert the verified rows into the rank stage's
own echelon mod r.

## Results (measured)

| rung | K | collisions | verified nontrivial | homogeneous rank | K−1 needed |
|---|---:|---:|---:|---:|---:|
| n=61 | 600 (frozen) | 10,188 | 588 | **49** | 599 |
| n=61 | 1200 | 12,732 | 4,800 | **746** | 1,199 |
| n=61 | 1500 | 10,924 | 6,000 | **1,197** | 1,499 |
| n=83 | 600 (frozen) | 9,672 | 72 | **6** | 599 |

- The mechanism is real: every nontrivial collision verified as four
  factor-base points with a sign choice summing to O.
- The yield is not: the rows are extremely correlated (the 588 rows at
  n=61 touch 168 of 600 columns; the 72 rows at n=83 span 6 dimensions),
  and 9,600 collisions at every rung are structured and trivial.
- rank(K) grows ≈ K⁴ with a small constant, so making the rank free needs
  K ≈ 1,500–1,900 at n=61 and K ≈ 2,600 at n=83 — where the census scan
  plus the K²-linear-algebra growth already exceeds the guided rank wall
  (13.1 s at n=61; 4.8 min with 12 threads at n=83), before charging
  memory (36 GB at n=83/K=2600).

## Verdict

Lead B is closed in both directions: the relations exist (positive
evidence, zero probes), but the rank stage does **not** collapse at any
rung the ladder runs, and forcing it to costs more end-to-end than the
guided rank it would replace. Any future free-relation claim must report
the homogeneous rank, not the relation count. The remaining "do" items of
the leads note (preprocessing-rho control, the (m,ℓ) = (4,4..6)
solving-degree cells, equivariant F4) are unaffected.

Boundary: stage census only. No target solved, no DLP recovered, no
vs_rho comparison, no ledger promotion.
