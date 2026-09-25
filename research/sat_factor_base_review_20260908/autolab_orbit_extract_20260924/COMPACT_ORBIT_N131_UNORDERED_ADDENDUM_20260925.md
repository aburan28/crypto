# Exact unordered four-sum coverage addendum

This addendum sharpens the necessary coverage bound in the [original compact-orbit n=131 transfer note](COMPACT_ORBIT_N131_TRANSFER.md). The original ordered-tuple bound and its arithmetic were valid but loose. The present result is an exact combinatorial **upper bound on coverage**, not a measured hit probability, memory lower bound, or solver forecast. The integer calculation is reproduced by [`unordered_four_sum_bound.py`](../../notes/ecc2k130/compact_base_sweep_20260925/unordered_four_sum_bound.py).

For a fixed factor base of `F` subgroup points, addition ignores the order of four selected points. Repeated selections are allowed in this optimistic count. There are exactly `C(F+3,4)` unordered four-multisets, and each has one group sum. Collisions can only reduce the number of distinct covered targets. Therefore a uniformly chosen target from the order-`q` subgroup has

```text
p <= min(1, C(F+3,4)/q).
```

For uniform nonidentity targets, use `q-1` in the denominator. The bound applies to any four-sum extractor over those same `F` points, regardless of its S3 implementation or SAT/linear-algebra engine. It does not imply that the current compact producer finds every covered target. Nor may a finite sample frequency be compared to the bound without sampling uncertainty.

At n=131, `q=680564733841876926932320129493409985129` and full signed-Frobenius columns give `F=262R`. The table keeps the original model's optimistic 24-byte record for each of `2nR²` potential roots. That is a **conditional materialized-index projection**; a different index may use different memory. Decimal EB means `10^18` bytes.

| Condition | `R` | `F` | Exact unordered coverage ceiling | Conditional packed-root bytes |
|---|---:|---:|---:|---:|
| Packed roots fit 1 TiB | 13,223 | 3,464,426 | `8.819525633e-15` (`2^-46.69`) | 1,099,442,519,952 |
| First orbit-rounded base with ceiling at least 1% | 13,644,854 | 3,574,951,748 | 0.01000000129 | 1,170,712,671,804,115,008 (1.171 EB) |
| First orbit-rounded base with ceiling at least 50% | 36,283,685 | 9,506,325,470 | 0.50000003516 | 8,278,188,452,662,966,800 (8.278 EB) |

Before orbit rounding, the exact smallest `F` thresholds are 3,574,951,633 for 1% and 9,506,325,303 for 50%. The original ordered bound put these at 1,615,166,730 and 4,294,967,297 respectively. Thus the original transfer note's 239 PB and 1.690 EB conditional projections are superseded by 1.171 EB and 8.278 EB for the **same** materialized-root design. At 1 TiB, the original ceiling `2.116682486e-13` is superseded by `8.819525633e-15`. None of these ceilings estimates the actual attainable coverage.

Changing the number of summands is a separate architecture hypothesis. The same necessary unordered-count calculation would first permit 1% coverage at `F=60,591,280` with five summands and `F=4,121,293` with six, compared with `F=3,574,951,633` with four. Those numbers do not transfer the four-sum index memory, S3 query work, relation rank, or solver cost to five or six summands. They motivate a newly preregistered small-rung solver and group-law test only after its own representation, cost model, and correctness gate are defined.
