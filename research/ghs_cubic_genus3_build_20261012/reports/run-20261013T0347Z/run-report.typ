#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)
= Run report: order-332 class at q = 7 (addendum), run 20261013T0347Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`.

== What ran
- Command: `rustc --edition 2021 -O genus3_build.rs && ./genus3_build 7 3 2 332`.
- Purpose: build the order-332 (= 4·83) isogeny class at q = 7. The frozen first-in-index selection of run 20261013T0320Z put all five of its q = 7 cases in the order-356 class; a target-order argument (the 4th) was added to `genus3_build.rs` in this commit to reach the other 4·prime class near q³ + 1 = 344.
- Regression: with no target argument, `./genus3_build 7 3 2` reproduces the q = 7 records of run 20261013T0320Z byte for byte (checked before this run).
- Spec: the construction and checks are those of `PROTOCOL.md` (unchanged); the target-order selection is the one addition, defaulting to off.
- Isolation: none; shared host. No timing is a result. Values are exact counts, deterministic.
- Raw records: `results/run-20261013T0347Z.jsonl` (`rec` = case, search).

#raw(read("../../results/run-20261013T0347Z.env.txt"), block: true)

== Recorded values
Requested 3 split and 2 irreducible of order 332; the (b, c) space was exhausted. Found: *2 irreducible, 0 split*. `candidates_examined` 65432, `linearly_dependent_skipped` 11232.

#table(columns: 9, stroke: 0.4pt,
  [*u (idx b, c)*], [*roots in K*], [*\#E0(K)*], [*t_E*], [*in F_q*], [*smooth*], [*N1,N2,N3*], [*\#Jac = \#E0*], [*divisor points on quartic / tested*],
  [49, 8], [0], [332 = 4·83], [12], [yes], [yes], [8, 50, 308], [yes], [912 / 912],
  [49, 9], [0], [332 = 4·83], [12], [yes], [yes], [8, 50, 308], [yes], [912 / 912],
)

Both cases: `N` equals `N_predicted = (q+1, q²+1, q³+1−3t_E)`, `N3` also counted on the K model and equal, every transfer divisor F_q-rational with matching Frobenius (912 of 912 points each). P1–P4 hold.

== Observation
Exhausting the parametrization `u = t² + bt + c`, the order-332 class has exactly two irreducible-u representatives meeting the genus-3 conditions and no split-u representative. (The order-356 run was not exhaustive — it stopped at the requested 3 split + 2 irreducible — so this is not a comparison of multiplicities between the two classes.)

== Missing data
- As in the parent run: no divisor-class arithmetic, so the transfer is not executed as a homomorphism; points of E0 whose preimages are not K-rational are not tested (178 per case).
