# Prospective successor after the rejected primary

No candidate below is implemented, timed or accepted by this record. A successor
must use its own source snapshot, protocol and unused holdouts. All primary
failures remain regression inputs; no successor changes this primary's verdict.

## Size-based dispatch is a selection experiment

A fixed rule could choose direct SIMD enumeration when n <= 20 and the
16-variable-leaf SIMD hybrid otherwise, with retained fallbacks outside each
method's supported domain. Those thresholds are suggested by the completed
primary data. They therefore cannot be tested as if selected before that data.

The contract must include the dispatcher and domain-check costs, exact equality
with the selected method, boundary cases n=20/21 and unsupported row/degree
domains, outcome/model checks, all enumeration work, and censored limits. A
comparison must include both constituent methods and all other retained arms.

Selection alone cannot beat the pointwise minimum of its constituent costs in
an ideal cost model: T_dispatch >= min(T_direct, T_leaf), before dispatch overhead.
It may make a convenient single policy, but it does not establish a new gain
against the current fastest-method boundary. Any historical pre-SIMD comparison
must be explicitly separate and cannot be labelled a current-frontier gain.

## A different mechanism: one initial affine restriction

The recursive affine policy rebuilds coefficient images at many nodes. A separate
mathematical hypothesis is to extract the original system's linear consequences
once, compile an injective parametrization x=A y+b once, and then run exact syndrome
evaluation in the remaining coordinates. This could remove an exponential factor
from enumeration when the initial affine rank is positive while avoiding repeated
per-node substitutions. Rank zero offers no assignment-count benefit.

The extraction must either span all affine rows obtainable from the stated row
space or clearly document a sound incomplete rule. It must never assume a
quadratic-support collision is sufficient without exact comparison. Inconsistent
affine rows imply UNSAT; otherwise Boolean diagonal and parity cancellations,
new quadratic coefficients, degree drops, and original-variable recovery must all
be handled. Transforming every equation is required. Costs include row reduction,
map and syndrome construction, enumeration, model recovery and verification.

The exhaustive domain changes from n variables to n-r free coordinates for affine
rank r, but early SAT stopping and assignment order also change. A 2^r ratio of
domain sizes is not a measured solve-time ratio. Tests should compare solution
sets through explicit recovery, not require models or traces to equal those from
a different ordering. A same-policy direct transformation should supply the
coefficient and enumeration-order reference.

Before any successor holdout timing, freeze its mechanism, budgets, fixtures and
acceptance gate against all retained complete methods, including direct SIMD and
both existing SIMD hybrids. Charge rank-zero failures and retain UNKNOWN values.
Do not remove the current SIMD methods from the reference to manufacture a gain.
The predecessor's universal 2x rejection remains unchanged regardless of a later
result. Full index-calculus, production and rho costs remain null.
