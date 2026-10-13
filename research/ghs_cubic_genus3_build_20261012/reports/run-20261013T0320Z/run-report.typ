#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)
= Run report: genus-3 build, run 20261013T0320Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`.

== What ran
- Command: `rustc --edition 2021 -O genus3_build.rs && ./genus3_build 5 1 1 && ./genus3_build 7 3 2 && ./genus3_build 11 1 1`.
- Spec: `PROTOCOL.md`, frozen before the run. Deterministic: no randomness, cases taken in index order.
- Isolation: none; shared host under heavy load. No timing is a result. All recorded values are exact counts.
- Raw records: `results/run-20261013T0320Z.jsonl` (`rec` = case, search).

#raw(read("../../results/run-20261013T0320Z.env.txt"), block: true)

== Runs requested and completed
#table(columns: 4, stroke: 0.4pt,
  [*Arm*], [*Requested*], [*Completed*], [*Candidates examined / skipped as linearly dependent*],
  [q = 5], [1 split, 1 irreducible], [2], [634 / 620],
  [q = 7], [3 split, 2 irreducible], [5], [2443 / 2394],
  [q = 11], [1 split, 1 irreducible], [2], [14657 / 14630],
)

== Recorded values
Unit: number of points (N1, N2, N3, divisor checks) or group elements (orders). Exact counts, no uncertainty. In all 9 cases `descends`, `smooth` and `jac_equals_E0` are true, `N` equals `N_predicted`, and `N3_K_model` equals `N[2]`.

#figure(image("orders.svg", width: 100%), caption: [Run 20261013T0320Z, 9 cases. Ring: \#E0(K) counted on the elliptic curve. Dot: \#Jac(C)(F_q) derived from the three measured point counts of the F_q quartic through the L-polynomial. Log axis.])

#figure(image("../construction.svg", width: 100%), caption: [What was built and which steps were measured.])

== Missing data
- No divisor-class arithmetic on the quartic: homomorphism and injectivity of the transfer were not executed.
- Points of E0 whose preimages on C are not K-rational were not tested (counts in `points_skipped`).
