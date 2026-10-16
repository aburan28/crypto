#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= Decomposition programme on ECC2K-130: plan report

*Date:* 2026-10-16. *Study:* `research/decomposition_programme_20261016/` (source: `PLAN.md`).
*Status:* plan only. Nothing in this report is a measurement, and no work below $2^61$ is claimed.

== Question

What is left for summation-polynomial decomposition on $y^2 + x y = x^3 + 1$ over $FF_(2^131)$, with $\#E = 4r$ and matched rho at $2^(60.8090)$?

== Baseline and closed branches

The baseline is the linear $(A, P)$ pair oracle over a geometric factor base (frozen, PR \#1442). Its gap to rho is 42.47 bits, or 55.22 bits charged with linear algebra at $omega = 2$.

These branches are closed by frozen experiments:
- higher-degree XL on the bilinear core (\#1420);
- XL on the trilinear system with $t$ unknown (\#1429);
- the three-summand linearization (\#1442);
- graph-cycle IC (\#1510);
- torsion and model invariants, and additive energy (\#1509).

== Scoping facts (proposed; not cited)

These came from throwaway checks in the planning session. They are re-established in Rust under frozen protocols in W1–W3.
+ The order of 2 modulo 131 is 130. So $FF_(2^131)$ has no small Frobenius-invariant subspace, and the same holds at $n = 13, 19$. At $n = 23$ the order is 11, which gives an 11-dimensional invariant subspace and makes it a positive control.
+ $dim(V V) >= 2 l - 1$ for every sampled $V$, and geometric $V$ attains it. This is the linear Kneser bound (Hou–Leung–Xiang). If it holds, the pair-oracle reach $l <= (n+1)/3$ is a ceiling for every subspace factor base.
+ Under the `gcost.py` heuristics, a one-shot $m = 4$ test fits the rho budget at about $2^11$ group operations per attempt when it handles same-conjugate decompositions only. It fits at about $2^34$ when it also handles mixed Frobenius conjugates. For geometric $V$, $dim(V dot V^(2^k))$ was 13–19 against the bound 9 at $n = 19, l = 5$.

== Workstreams and gates

#figure(image("workstreams.svg", width: 100%), caption: [Planned workstreams, their order and their outcomes. Source: `workstreams.svg`, hand-written, in this directory. No box is a measurement.])

- *W1, ceiling note.* Statements C1–C3, their proofs and scope, a Rust verifier for the $2 l - 1$ bound, and an identity certificate. Do this first.
- *W2, Frobenius-orbit factor base.* Measure the excess $e(l)$ of $dim(V dot V^(2^k))$ over $2 l - 1$. Kill if $e(l) >= l - 2$; reopen if $e(l) <= 2$ with no trend in $l$.
- *W3, four- or five-summand test.* Freeze the budget first. Then run a bounded search over at most 6 candidates, with a glued-pair-oracle negative control and an analytic-then-toy gate.
- *W4, rho constants.* Optional and independent. Uses ecbench measure and bounds, and claims nothing beyond constants.

== Success and kill rules

Any verified improvement over the baseline counts as a success, reported with its distance to rho (owner criterion, 2026-10-15). A branch closes only when no candidate can improve the baseline.

== Graphs checked

No canonical graph changes, because this PR lands no measured value. No scoreboard, leaderboard or timeline entry is added or altered.
