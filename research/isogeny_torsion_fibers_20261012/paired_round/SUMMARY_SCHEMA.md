# Summary floor interpretation

The `mean_s_over_floor` field in `summary.tsv`, `summary-v2.tsv`, and
`summary-v4.tsv` is the mean counted-group `S` divided by
`sqrt(pi/4)`, the two-class **negation** walk-step floor selected in the
frozen protocol. It is not the floor for every rational automorphism
available on these `j=0` curves.

For both registered pairs, `p ≡ 1 (mod 3)` and the action
`(x,y) ↦ (ζx,y)` has order three. Together with negation, it gives six
rational automorphisms. The corresponding curve walk-step floor is
`sqrt(pi/12)`. Divide a mean `S` by that floor, or multiply the stored
`mean_s_over_floor` by `sqrt(3)`, to obtain the six-class ratio. The
independent verifier checked the cubic action's group law on all 89,764
ordered point pairs of the four curve models. See the [result report](REPORT.md).

This note corrects the interpretation of a derived column. It changes
neither the frozen protocol nor any raw or derived numeric record.

The raw and phase-summary `row_mul_subs` counter is also broader than
its short name: it counts modular products in elimination,
back substitution, and pivot/RHS normalization. `row_inversions` is
separate. These native operations have no pinned group-operation
conversion in this round.
