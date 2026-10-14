# Source-only capacity census of the v4 cross-stream recipes

The [saved census](census.json) covers every one of the 169,884 ordered,
distinct-stream, same-arm public-X pair recipes in the frozen
[v4 design](../../paired-frontier-v4.json). The current signed-Frobenius
formula applies to 17,424 recipes whose two closures match that worker. It
evaluates four shapes for those recipes: five or six summands with one or two
residual summands, giving 69,696 pair/shape capacity calculations. The
remaining 152,460 recipes have a cross-closure or non-signed adapter
disposition; no signed-orbit cap verdict is assigned to them. All 679,536
pair/shape cells have an explicit capacity or adapter disposition. These
are **not** SAT instances or an execution of the full 26,939,278,080-case
hyperparameter grid. The other solver axes retain their v4 source
dispositions.

For a base of (B) signed orbits and a separate pool of (P), the unknown
nonempty residual domain has between `max(1, P-B)` and `P` orbits. The screen
evaluates both endpoints with the exact first-pair variable and orbit-domain
clause formulas already [checked against two retained nested trials](../paired-capacity-screen-20261013/RESULT.md).
Both counts increase with residual size. A pair that fits at `P` fits the
current 150,000-variable and 12-million-domain-clause caps for every nonempty
exact import; a pair above a cap even at the lower endpoint cannot fit this
worker configuration with a nonempty residual domain. The middle class
requires an exact overlap measurement. A fully absorbed pool has no residual
domain for these positive-residual shapes.

| Shape, total/residual summands | Fits even at zero overlap | Depends on exact overlap | Over cap for every nonempty residual |
| --- | ---: | ---: | ---: |
| 5 / 1 | 14,112 | 144 | 3,168 |
| 5 / 2 | 14,112 | 288 | 3,024 |
| 6 / 1 | 11,952 | 144 | 5,328 |
| 6 / 2 | 11,808 | 720 | 4,896 |

Each row sums to 17,424 signed-Frobenius ordered pairs, with 152,460 other
pairs separately marked adapter-unavailable in the saved JSON. Some signed
pairs use the other curve arm or lack exact pilot archives; the table is
conditional arithmetic for those recipes, not an admission to the current
worker. The frozen v4 source classification gives 648 primary signed-Frobenius
pilot pairs with archive identities and an exact-import path. **All 648 fit
both current caps at zero overlap in all four shapes.** That removes these
two caps as a reason to exclude a pilot pair before import. It does not
predict memory peak, conflict count, relation yield, or complete runtime.

The sole retained cross-stream import in this panel is v4 pair index 36,953,
K600 public-X hash versus K600 public-X sequential. The copied
[raw search result](retained-first-result.json) matches the SHA-256 in its
[independent audit](../cross-stream-first-gate-source-20261013/audit.json).
Its exact overlap was zero. The observed 97,080 variables, 269,830
domain clauses and 4,980 XOR rows match the zero-overlap formula for the
five-summand/two-residual shape. Search reached `UNKNOWN_search_cap` after
5,000 conflicts on its first identity pattern. That is a censored observation,
not a yield estimate. The frozen v4 source-class label still says
`pilot_exact_import_gate_available_unexecuted` because the design preceded
this one retained run; the census reports its outcome separately.

`census.py` pins the v4 JSON and generator, earlier capacity formula source,
audit and copied raw result by SHA-256. `test_census.py` checks complete pair
and shape counts, endpoint monotonicity, the retained counter replay, and
rejection of altered raw evidence. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-capacity-census-20261013/census.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-capacity-census-20261013 -p 'test_*.py' -v
```

This source-only census launched no container, import or SAT search. It adds
zero local seconds to the already charged 19,800-second ledger, creates no
new S3 object, and leaves the minimum complete cold IC runtime unresolved.
