# P-192 two-pivot genus-three representatives and guided three-sum coverage

**Result.** Sixty-four exact quotient classes selected from the two new
residual-pivot fans now have fixed-origin degree-three representatives on
the smooth P-192 genus-three plane quartic. The independent PARI/GP replay
checks all 64 two-cubic reductions. Adding 16 classes from each fan to the
original 16-class base gives 19,600 distinct three-sum targets, of which
4,926 satisfy the frozen two-character geometric query predicate. Both
frozen 64-query known-scalar laws have zero hits on this 48-class base.

| Frozen base | Usable classes | Unordered triples | Distinct targets | Geometrically admitted `++` targets | Triple-collision rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original | 16 | 816 | 816 | 224 | 0 |
| Original + pivot 000 | 32 | 5,984 | 5,984 | 1,535 | 0 |
| Original + pivots 000 and 001 | 48 | 19,600 | 19,600 | 4,926 | 0 |

The selection takes the first 32 distinct exact projected quotient points
from each frozen pivot fan, scanning fan rows and residual sides in recorded
order. For every selected point `Q`, the checked fiber constructor takes
`H = Q/2` in the odd-order subgroup, builds the two inverse-cover fibers,
and cuts their six Frobenius support points plus three fixed rational
anchors with cubics. The second pair of cubic cuts yields a degree-three
divisor `D` relative to the same origin `T` for all classes. The section
identity is `[D-T]=[B-A]`, with `A` and `B` the positive and negative input
residual cubics. [The reduction receipt](../../type_i_two_pivot_reduced64.json)
contains all 64 exact quotient coordinates, cubic coefficients, residual
polynomials, and representatives. [The independent GP log](../../type_i_two_pivot_reduced64_native.log)
records `TWO_PIVOT_FIXED_ORIGIN_NATIVE_PASS cases=64 prime_bits=192` after
checking both resultants, prescribed intersections, squarefreeness,
residual coordinates, rank-nine interpolation, and a noncollinear final
degree-three divisor for every class. [The selection replay](../../type_i_two_pivot_reduced64_replay.json)
binds those cuts to the original fan receipts. The first Sage attempt
stopped after all 64 cubic inputs because the wrapper omitted the reducer's
character-cell field; its [failed log](../../type_i_two_pivot_reduced64_attempt1.log)
is retained. The corrected source recomputes that cell from the frozen
sextic predicate.

The [coverage producer](../../line_pivot_guided_coverage.py) enumerates
every unordered triple with repetition over the 16-, 32-, and 48-class
bases. It counts a target as admitted exactly when both characters at
`P=-Q/2` and `2P` are positive, the same predicate used by the earlier
known-scalar guided collector. The two frozen scalar lists come from the
existing line-base collector: 64 uniform-law scalars and 64 accepted
character-guided scalars. Each query is checked as
`Q=-2[a]G=[a](-2G)`, so a hit would give an exact logarithm equation in
the three base classes. The [coverage receipt](../../type_i_two_pivot_guided_coverage.json)
records 11 `++` cases in the uniform list and all 64 in the guided list;
neither list hits any of the three declared supports. Its local 165.954 s
support-and-query interval is a shared-host diagnostic.

The certified accepted-scalar lower bound from the earlier character
density analysis is
`L=1569275433846670190958947353987675545133315652035359680560`.
For an ideal uniform scalar drawn from the accepted set, the measured
48-class support therefore gives a per-query hit bound `4926/L < 2^-177`.
This bound concerns the declared fixed base and query law. Adding these
two fan slices increases admitted support by `4926/224`, but their
three-sum map is injective on the tested 48-class base, so this slice
does not produce internal collision rows. The full two-pivot fans contain
20,014 distinct signed columns; their fixed-origin reduction and
relations with a larger usable base remain the next gates.

The [Lean module](../../lean/PivotGuided.lean) proves that any hit at the
known-scalar target yields the stated logarithm equation and checks the
48-class capacity and integer probability certificate. The abstract
fixed-origin cubic class identity is already proved in
`P192ExtensionCover.fixed_origin_cubic_class_identity`; the concrete
P-192 sections are checked by GP. The finite-field cell counts are
checked by the producer and [independent integer-curve replay](../../type_i_two_pivot_guided_coverage_replay.json),
and are inputs to the Lean arithmetic implication. The pinned
`lake build PivotGuided` [completed](../../lean/pivot-guided-build-20261016.log)
with 1,213 jobs; `#print axioms` reports only Lean's `propext`,
`Classical.choice`, and `Quot.sound` foundations. The first direct Lean
invocation lacked compiled local imports and is retained as a
[failed build](../../lean/pivot-guided-build-20261016.failed.log).

| Result | Lean declaration | Concrete verification |
| --- | --- | --- |
| Known-scalar three-sum equation | `guided_scalar_three_sum_log_equation` | Exact query point and scalar replay in the coverage receipt and independent verifier |
| Fixed-origin divisor identity | `fixed_origin_cubic_class_identity` | 64 independent GP cubic-cut certificates |
| 48-class capacity and accepted-set bound | `two_pivot_48_class_three_sum_capacity`, `two_pivot_guided_hit_bound_certificate` | Exact support/cell census plus prior certified accepted-set lower bound |

The proposal retains `candidate_id: null`: the full-pipeline configuration,
ordinary-query relation yield at a scalable base, useful rank, final
linear algebra, and unseen-target descent are not yet measured for this
P-192 collector. The current exact support counts show what the new
geometry contributes and how far this fixed-base query policy must grow.

Reproduce from the `crypto` worktree root with the checked local Sage
launcher:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_pivot_reduced_classes.py
python3 -B research/p192_extension_cover_20261012/prepare_two_pivot_gp.py
gp -s 256000000 -q research/p192_extension_cover_20261012/verify_two_pivot_reduced_native.gp
python3 -B research/p192_extension_cover_20261012/verify_two_pivot_reduced_receipt.py
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_pivot_guided_coverage.py
python3 -B research/p192_extension_cover_20261012/verify_two_pivot_guided_coverage.py
cd research/p192_extension_cover_20261012/lean && lake build PivotGuided
```
