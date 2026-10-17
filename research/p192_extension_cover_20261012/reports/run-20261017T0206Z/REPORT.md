# P-192 all-pivot residual geometry and fixed-base support bound

**Result.** The complete 141-pivot line screen yields 700,993 split cuts and
1,401,712 distinct rational residual points on the descended genus-three
quartic. Every one of the 133 cross-fan repeated points is an original
selected point on the two pivots' source lines. Each pivot's single
within-fan repeat is its source partner. The resulting geometry has no
additional repeated residual point across the declared 141 fans. Even if
every selected and residual geometric point supplied a distinct usable
quotient class, the fixed base would have at most 1,411,637 labels, and its
unordered three-sum support would have at most 468,833,646,178,255,739
targets. Under the certified accepted-scalar lower bound, an ideal uniform
accepted query would hit that support with probability below \(2^{-131}\).

| Exact P-192 screen | Value |
| --- | ---: |
| Field prime \(p\) | 6277101735386680763835789423207666416083908700390324961279 |
| Selected quartic points | 9,925 |
| Residual pivots | 141 |
| Split / nonsplit selected-pivot cuts | 700,993 / 698,432 |
| Distinct residual geometric points | 1,401,712 |
| Within-fan / cross-fan repeats | 141 / 133 |
| Cross-fan repeats outside the original source lines | 0 |
| Generous quotient-base cardinality cap | 1,411,637 |
| Unordered three-sum support cap | 468,833,646,178,255,739 |

For each pivot, the [producer](../../line_all_pivot_geometry_screen.py)
reconstructs the quartic line through that pivot and each of the 9,925
selected points. It divides the quartic intersection polynomial by the
known pivot and selected roots, tests the remaining quadratic's
discriminant, and enumerates its rational roots. One batched inversion per
pivot supplies the line slopes. The [full receipt](../../type_i_pivot_geometry_screen_141.json)
binds the input hashes, all 141 per-pivot count rows, the ordered residual
stream digest, and every cross-fan reuse witness. The 284.071 s interval is
a shared-host producer diagnostic. The interrupted inversion-per-line
attempt and its ten-pivot [log](../../type_i_all_pivot_geometry_screen.log)
remain recorded; the completed batch-inversion source produced the frozen
[32-pivot prefix](../../type_i_pivot_geometry_screen_032.json) and then the
full screen.

The [independent verifier](../../verify_pivot_geometry_screen.py) expands
the plane quartic directly as a polynomial along each sampled line,
checks 1,409 stratified line cases and every one of the 133 cross-fan
witnesses, and verifies the source-pair geometry for every pivot. The
[verification receipt](../../type_i_pivot_geometry_screen_141_verification.json)
also binds source and input hashes. Its classification shows that the
cross-fan repeated points are selected indices 141--145, each from an
original source line. The 141 within-fan repeats are accounted for by two
presentations of each source partner. The shorter 32-pivot receipt has its
own [replay](../../type_i_pivot_geometry_screen_032_verification.json).

The quotient projection maps each geometric quartic point to at most one
subgroup class. Thus the union of 9,925 selected points and 1,401,712
residual points supplies at most 1,411,637 quotient labels; this addition
deliberately counts known geometric overlaps twice. The
[capacity producer](../../line_all_pivot_capacity.py) combines that cap
with the certified accepted nonzero scalar-set lower bound
\(L=1569275433846670190958947353987675545133315652035359680560\).
Its [receipt](../../type_i_all_pivot_capacity.json) and
[independent integer replay](../../type_i_all_pivot_capacity_verification.json)
give \(\binom{1411637+2}{3}=468833646178255739\) and
\(468833646178255739\cdot2^{131}<L\). The first unordered arity whose
*capacity* reaches \(L\) is 11: the ten-sum capacity is
8659267211365270960962016489109601808335080570368176291, whereas the
eleven-sum capacity is
1111257143738377332384467971891100187630071952901684086969207.
These are upper capacities for the declared fixed-base policy; collisions
between sums, subgroup filtering, and query admissibility can only shrink
support.

The Lean [UnorderedSupport.lean](../../lean/UnorderedSupport.lean) theorem
`unordered_policy_support_card_le_multichoose` proves the structural
bound for the image of **any** map from unordered \(m\)-tuples of a finite
base. The same module checks the three-sum capacity, the \(2^{-131}\)
integer inequality, and the ten/eleven arity threshold as separate
theorems. The pinned [`lake build UnorderedSupport`](../../lean/unordered-support-build-20261017.log)
completed with 935 jobs. `#print axioms` for all four declarations lists
only `propext`, `Classical.choice`, and `Quot.sound`; `lean4checker --fresh`
was unavailable on this host. The concrete P-192 point counts and
source-line incidences are established by the finite-field receipts, not
by the generic Lean theorem.

| Mathematical step | Lean declaration | Concrete bridge |
| --- | --- | --- |
| Unordered policy image has at most \(\binom{B+m-1}{m}\) outputs | `unordered_policy_support_card_le_multichoose` | The collector uses an unordered fixed-base tuple policy |
| Generous P-192 three-sum capacity | `all_pivot_three_sum_capacity` | Complete geometry screen and independent source-overlap replay |
| Accepted-query probability bound | `all_pivot_three_sum_hit_bound` | Frozen accepted scalar lower bound and independent integer replay |
| Ten/eleven capacity threshold | `all_pivot_arity_capacity_threshold` | Same frozen base and accepted-set sizes |

The earlier [two-pivot reduction and guided collector](../run-20261016T2324Z/REPORT.md)
give fixed-origin degree-three Jacobian representatives for 64 exact
quotient classes and a 48-class, known-scalar three-sum test under uniform
and character-guided query laws. This full geometric census identifies the
next discriminating experiment: check whether *distinct* residual quartic
points from separated pivots project to the same subgroup class, then
measure useful relation rank under a frozen ordinary-query law. Geometric
uniqueness alone is insufficient for the subgroup-column condition in
`line_private_pairs_extend`; a projected collision audit is the required
bridge. The exact tuple-capacity result also directs higher-arity
collector designs toward arities whose support can cover a material
fraction of the accepted scalar set, with point-decomposition and cost
measured separately.

Reproduce from the `crypto` worktree root:

```sh
python3 -B research/p192_extension_cover_20261012/line_all_pivot_geometry_screen.py --pivots 141
python3 -B research/p192_extension_cover_20261012/verify_pivot_geometry_screen.py --pivots 141
python3 -B research/p192_extension_cover_20261012/line_all_pivot_capacity.py
python3 -B research/p192_extension_cover_20261012/verify_all_pivot_capacity.py
cd research/p192_extension_cover_20261012/lean && lake build UnorderedSupport
```
