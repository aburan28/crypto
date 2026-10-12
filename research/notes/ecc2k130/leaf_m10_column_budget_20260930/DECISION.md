# Decision: equal-column thinning loses m10 capacity on the measured leaf slots

The two fixed degree-263 leaves in the [exact m10 support census](../leaf_m10_support_outcome_20260930/DECISION.md)
have 26–32% more **raw** ordered physical tuples than the source, but they
also have about 41,000–45,000 distinct uncompressed `[4]` sign classes across
the ten measured slots. The source uses cheap same-curve Frobenius to normalize
the balanced and unequal slots to, respectively, **3,988** and **8,062**
nonzero log columns in the pinned m10 screens. Coordinate squaring does not
act on either fixed leaf. This note asks how much raw tuple capacity remains
if a native leaf policy is restricted to the **same number of log columns**
and only selects from the archived, mutually disjoint rotated slots.

Let `c_i` be selected nonzero projected sign classes in slot `i`. The full
raw census proves there are no `[4]` collisions within or between these
slots, so each selected class has exactly two physical sign lifts and each
slot can also use its single x=0 point. With `Σ c_i ≤ B`, the physical
ordered tuple count is `N = ∏_{i=1}^{10}(1+2c_i)`. Moving one class from a
larger slot to a smaller slot increases this product whenever their counts
differ by at least two. Thus the exact integer maximum has `B mod 10` slots
with `⌈B/10⌉` classes and the rest with `⌊B/10⌋`; [bound.py](bound.py)
computes it directly from the pinned source budgets and the replayed leaf
census. For a uniformly distributed subgroup target, at most `N/q` of the
targets can be represented by these tuples, regardless of the solver.

| Frozen arm | Source normalized columns `B` | Maximum tuples with leaf `B` | Necessary uniform-target support ceiling `N/q` | Minimum disjoint-slot columns for `N/q ≥ 1` | Minimum / source columns |
| --- | ---: | ---: | ---: | ---: | ---: |
| Balanced m10 | 3,988 | 105,509,332,938,070,832,394,432,577,609 | 1.55032031034663e−10 | 38,213 | 9.581996× |
| Unequal m10 | 8,062 | 119,514,333,106,613,002,600,816,201,540,225 | 1.75610529261396e−7 | 38,213 | 4.739891× |

The shared threshold 38,213 is only where the **counting ceiling** first
reaches one. It is not a predicted factor-base size, a lower bound on every
possible leaf design, or evidence that one target can be solved there. The
full measured leaf domains have 40,875/40,994 balanced or 44,932/45,075
unequal uncompressed columns, so they clear this necessary capacity threshold
only by paying far more log columns than the source's normalized budgets.

**Gate decision.** The particular strategy “thin the measured disjoint leaf
slots to the source's normalized column budget” cannot have broad uniform
target support: its rigorous ceilings are `1.55e−10` and `1.76e−7` per
target. This is a conditional counting no-go for that factor-base policy,
not a native-leaf or ECC2K-130 no-go in general. Reusing classes across
overlapping slots, finding a low-cost leaf log action, or changing the
factor-domain geometry can break the disjoint-slot premise. Such a proposal
must freeze its actual distinct projected columns, account for repeated
factor/permutation collisions and action evaluation, independently measure
held-out Q support and rank, and compare all four policies at equal useful
size with full cold cost. The separate review-gated m10 representation
[PR #937](https://github.com/aburan28/crypto/pull/937) is not released by
this arithmetic result. Natural-target PDP yield, full ECDLP cost, and
matched-rho crossover remain unset.

[BOUND.json](BOUND.json) preserves the exact products, subgroup order,
input hashes, threshold, and the original full-domain leaf counts.
