# Uniform-query support bounds for retained N83 pilot pairs

The [frozen screen](screen.json) applies exact integer counting to all **648**
ordered cross-stream pairs in the independently replayed [archive overlap
receipt](../pilot-overlap-census-20261013/overlap.json). Ninety pairs have fewer
than two residual signed-Frobenius orbits and remain unavailable to the current
worker. For each of the other **558** pairs, the screen covers the four current
five/six-summand, one/two-residual shapes. These are 2,232 conditional bounds,
not SAT searches or measured relation yields.

For a small domain of `B` points, disjoint residual domain of `L` points,
`u` small summands and `v` residual summands, the number of point multisets,
including repetitions, is exactly
`M = C(B+u-1,u) C(L+v-1,v)`. The archived signed orbits each contain 166
points, and the exact overlap receipt supplies the residual orbit count.
Each multiset has one group sum. If a query target is independently selected
uniformly from the **nonidentity** elements of the primary subgroup of order
`r = 2417851639230796216685689`, the number of nonidentity sums represented
is at most `M`; hence `Pr(at least one decomposition) <= min(1,M/(r-1))`.
This argument needs no independence among the sums. For multiple queries with
uniform marginals, the union bound implies that at least
`ceil((r-1)/(2M))` queries are necessary before its hit ceiling can reach
one half. These are ceilings and necessary query-count floors, never predicted
success probabilities or sufficient trial budgets.

| Current-worker shape | Admissible pairs | Nonvacuous hit ceilings | Largest `M/(r-1)`, approximately |
| --- | ---: | ---: | ---: |
| Five summands, one residual | 558 | 558 | 0.168920347 |
| Five summands, two residual | 558 | 558 | 0.337833909 |
| Six summands, one residual | 558 | 209 | Vacuous for the largest supports |
| Six summands, two residual | 558 | 198 | Vacuous for the largest supports |

Thus even an ideal exhaustive oracle has a hit probability below 0.338 for
one uniformly chosen nonidentity query in **any** retained-pilot five-summand
pair under the 3+2 split. In particular, a single bounded five-summand trial
with no model cannot by itself establish poor solver access or poor base
quality. The previously screened K1182/K2048 3+2 count was larger than
`r-1`, so its corresponding ceiling was vacuous; this new result concerns the
smaller, exactly replayed primary pilot archives.

A concrete future comparison can keep the small base fixed at K600 sequential
seed `2026100802` and the nominal pool size fixed at K600. The Gray-prefix
seed `2026100801` pool shares 503 orbits and leaves 97 residual orbits; its
3+2 single-query ceiling is approximately **0.0088301243**, with a 57-query half-ceiling
floor. The hash seed `2026100801` pool shares none and leaves 600 residual
orbits; its ceiling is approximately **0.3378339093**, with a two-query floor. The receipt
pins all three archive hashes and exact fractions. Pool policy changes
between arms, so any future timing or yield difference would need further
controls before being attributed to overlap alone. These figures are useful
for allocating trial effort, not for selecting a total-runtime winner.

The proof concerns a fresh uniform query distribution. It does **not** apply
to the frozen public target, a nonuniform deterministic coefficient schedule,
or a target chosen after inspecting the bases. It does not measure SAT
construction, solver conflicts, large-prime graph closure, useful matrix
rank, linear algebra, individual logarithm, or complete cold cost. The
current local charge remains **22,200 seconds** and no N83 solver call was
made for this screen. S3 objects were not reread: the input is the saved
archive-backed overlap receipt, whose SHA-256 and audit chain are pinned in
the new JSON; the separately retained [S3 content audit](../s3-stream-integrity-20261014/RESULT.md)
covered all 57 objects earlier that day.

Reproduce the saved result and focused checks from the repository root:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-support-screen-20261014/screen.py
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-support-screen-20261014 -p 'test_*.py' -v
```

The tests independently enumerate a small cyclic group to check the
multiset ceiling, pin a same-base pair contrast, and reject a changed
reciprocal overlap. The screen's `--write` mode creates the receipt once;
normal invocation checks its exact bytes, including the source and input
hashes.

Validation on this branch: the saved receipt rechecked byte-for-byte; its
four focused tests, the three prior overlap tests, the 63-test N83 study
suite, and the 43-test boundary autolab suite passed. The overlap suite used
`uv --with blake3` because the system Python lacked that dependency. No Rust
source or example changed in this screen. The prior release-library attempt
on this branch is [recorded as interrupted](../known-witness-stage-preflight-20261014/RESULT.md);
this Python and documentation change does not turn it into a full-suite pass.
