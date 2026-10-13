# Shape replay for the `k = 4` Macaulay solve: a 2–5 % lever, not a 10× one

2026-10-13.  Step 4 of [`../ic_speed_and_rho_scope_20261013/`](../ic_speed_and_rho_scope_20261013/README.md):
the first cut at a compiled `S₅` solver for Gaudry's index calculus over
`F_{p⁴}`, whose per-decomposition cost `C₄` (§11.17 of
`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`: `1.227·10¹²` `F_p`
multiplications, 47 % echelon, 45 % normal forms, 6 % characteristic
polynomial) sets the crossover with rho near a group of `2^{150}`.  The
scope note's decision rule: continue only if `C₄` falls at least 10×.

## What was built

`gaudry_quartic::solve_system4_traced` (the public `solve_system4` is
unchanged and delegates to it).  Every trial solves a system of the same
monomial support, so the Macaulay matrix has the same shape each time.  The
trace records, on a first system, the leading monomials after
triangularisation, the rows built, the input rows that reduce to zero, the
pivot count and the quotient dimension.  A replay does not build the
zero rows; if its pivots, dimension, leads or row count differ from the
trace, the full solve runs and both costs are charged.  This is
F4Remake's record-and-replay applied to the Macaulay design.

`examples/gaudry_quartic_replay.rs` runs the full solve (the counter §11.17
measured `C₄` with) and the replay on the same coefficients and compares
counts, wall time and solution sets.

## Result (generic planted systems of degree `d`, 4 systems each, `p = 269`)

Raw: [`raw/generic_d3-5.jsonl`](raw/generic_d3-5.jsonl).

| `d` | Macaulay degree | quotient dim | zero rows skipped | full `F_p` muls | replay ÷ full (muls) | replay ÷ full (wall, 4 systems) | solutions agree, planted found |
|--:|--:|--:|--:|--:|--:|--:|:--|
| 3 | 9 | 81 | 41 of 675 | `6.4·10⁶` | 0.980–0.981 | 0.70–1.45 (noise) | 4 / 4 |
| 4 | 13 | 256 | 190 of 2,314 | `2.3·10⁸` | 0.965 | 0.91–1.10 | 4 / 4 |
| 5 | 17 | 625 | 570 of 5,930 | `3.7·10⁹` | 0.952–0.953 | 0.95–1.02 | 4 / 4 |

Every replay matched its trace; every solution set equals the full solve's.
The saving grows slowly with `d` (2.0, 3.5, 4.7 %); on the trend, about 7 %
at the real `d = 8`.  Zero rows are a small share of the elimination, and the
rows that become pivots cannot be skipped: their reductions are the
solve.  Wall time moves within the noise of these short runs (shared host).

**Decision (the scope note's rule): record-and-replay does not cut `C₄`
10×.**  It stays in the tree as the harness for the next lever: a replay
that matches its trace is also the precondition for compiling the
elimination into a fixed straight-line kernel, which can cut wall time but
not the counted multiplications that set the crossover.

## Where 10× has to come from

The counted cost scales with the quotient dimension `D` (the Bézout bound,
4,096 at `k = 4`): from the §11.17 panel, `C` grows about as `D^{3.1}`
between `d = 5` (`D = 625`, `3.7·10⁹`) and the `S₅` solve (`D = 4,096`,
`1.2·10¹²`).  The requirements document's candidate C8 (the
`D₄ = (Z/2)³ ⋊ S₄` invariants at `m = 4`, Bézout 4,096 → 512 on paper,
after FGHR's use of 2-torsion) would put `C₄` near `D = 512`: on that
trend about `2·10⁹`, some 600× below today, which would move the
crossover by about `2^{36}`.  That needs a curve model over `F_{p⁴}` with
the right rational torsion and the invariant rewriting of `S₅`; the
`k = 3` analogue exists (`fghr_full`, E17 of
`research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`).
It is the next step this programme should take.

Not yet run: the replay on the real `S₅` systems (`--s5`, about an hour a
pair at `p = 269` and 3.4 GB a matrix), deferred until the ecbench ladder
in `../pdp_root_finder_ecbench_20261013/` has finished, since it would
compete with its timed runs.

## Reproduction

```sh
cargo build --release --example gaudry_quartic_replay
target/release/examples/gaudry_quartic_replay --generic 3,4,5 --systems 4 --json OUT.jsonl
target/release/examples/gaudry_quartic_replay --s5 --p 269 --seed 1 --systems 3 --json OUT.jsonl
```
