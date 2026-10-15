# Next discriminating gates after full-image setup reuse

**Execution update, 2026-10-15:** the bounded Gate D candidate screen
is [complete](../different_j_order_screen_20261015/RESULTS.md). It
selected the first eligible pair at `p=211`, `r=229`, with different
`j=77` and `160`; independent replay checked 192 map witnesses and a
direct `F_(p²)` census matched both curve orders at each of three
accepted primes. A later separately frozen
[full map/fiber audit](../different_j_order_screen_20261015/FIBER_REPORT.md)
checked all 44,655 source points, 14,885 three-point fibers and
52,441 subgroup pairs with independent native replay. The full
64-target relation/cost comparison and
Gate A's direct χ-coordinate resultant remain pending. The frozen
protocol text below is retained as the pre-screen decision rule.

These are **protocols with the execution status above**. The completed
[60-target result](RESULTS.md) exhausts every nonidentity target in
the frozen order-61 subgroup. Repeating that same group with duplicate
targets would not test further genuine setup reuse.

## Gate D: increasing-order different-j pair

First run a bounded native Rust screen over primes `p≤271` for the
same source family and rational degree-three kernel, requiring a
different destination `j`, an explicit checked Vélu map, and a prime
subgroup order `r≥211` coprime to three. Record every candidate and
rejection with exact field/curve identity and cost. If no candidate
passes, retain that bounded result and preregister the next prime
bound before widening the search. A bounded miss says nothing about
all curve families.

For the first accepted pair, commit a public source generator, target
point order formed by a seeded hash of point coordinates, separate
checker scalars, source/projected-image/destination-direct/full-image
factor bases, and direct source BSGS, rho and rho A/A controls.
Preregister prefixes `N=1,4,16,64` distinct targets, relation and
rho caps, arm rotation, 512 MiB RSS and per-arm wall ceilings, source
and binary hashes before timing. Enumerate the rational image once
per arm process and charge all failed attempts, pair-table work,
linear recovery, verification, process startup and output. Independently
replay every scalar and relation witness. Report cold and per-target
amortized totals in one counted unit `S=total operations/√r`, with
measured field/group conversions and a matched direct rho reference;
retain process CPU/wall and memory as separate units. A direct BSGS
control is required at each prefix because it was cheaper than rho
throughout the order-61 panel.

If the larger pair's full-image support differs, test whether its
relation-yield change survives same-stream, equal-memory paired
collection and whether the *whole* method crosses the cheaper direct
control. Equal support, capped rank, or no eligible pair narrows that
specific gate only. Isolated interleaving and second-host replay are
required before promoting a performance result.

## Gate A remainder: direct χ-coordinate resultant

The bounded nonfactorable χ finite-scan and direct source/quotient
`S3` roots are already recorded in the stacked predecessor work.
A direct χ-coordinate resultant/root solver has not been executed.
Implement it natively, prove the exact finite decomposition contract
against independent pair enumeration on the smallest frozen cell,
then run complete relation collection and scalar replay on the same
four targets with direct source and quotient roots, projected-index
and rho controls. Charge polynomial construction, factorization,
all failed roots, memory, linear recovery and verification. Its
potential algebraic saving remains a hypothesis until that full gate
and a portable counted-unit comparison complete.
