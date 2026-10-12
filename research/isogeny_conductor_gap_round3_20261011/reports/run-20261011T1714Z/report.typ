#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Run report: E4, two-torsion quotients of the S#sub[3] system
_2026-10-11, 17:14–17:15 UTC. research/isogeny_conductor_gap_round3_20261011, lead E4._

*Command.* `isogeny_gap_two_torsion_f4 --sizes 32,64 --budget-ms 20000` on the round-1 prime pairs (`research/isogeny_conductor_gap_ic_20261010/pairs/prime_pairs.json`). Runner `examples/isogeny_gap_two_torsion_f4.rs` (commit fb8185437 plus the orbit fix committed with the records); solver `cryptanalysis::f4_fp::f4`, grevlex, degree bound 2·deg + 8, 20 s budget. Host: cloud container, not an isolated ecbench host; wall times are observations only.

*Workloads.* 24 even-cofactor prime pairs (22 and 26 bits) × crater, floor, isomorphic floor model = 72 curves; base sizes B = 32, 64; subgroups 1 (direct), C2, C4, and E[2] where rational; two random and one planted target each.

*Completed.* 1,128 of 1,128 systems (752 random, 376 planted). No timeout; F4 agrees with enumeration on every system; every planted solution annihilates its basis; the Vélu map passed 32 random-sum checks per (curve, subgroup). A first launch looped on a C4-stable base whose orbit repeated an abscissa; it was stopped, fixed, and restarted from an empty record file.

*Rational 2-torsion.* Rank 2 on the crater and rank 1 on floor and model for all 12 D = −7 pairs and the D = −11 and D = −3 pairs with ℓ = 2; rank 1 everywhere else.

*Measured (exact counts; range over curves, roles, targets).*
#table(columns: 5, align: right,
  [base-poly degree], [solving degree], [widest matrix], [mults, random], [mults, planted],
  [8], [10], [45], [7,408–7,414], [13,903–13,929],
  [16], [18], [93], [34,152–34,158], [81,128–81,178],
  [32], [34], [189], [145,630], [483,552–483,561],
  [64], [66], [381], [600,510], [3,123,348–3,123,354])

*Derived.* Direct ÷ quotient at matched B: solving degree 1.89–3.67; multiplications 4.12–4.26 (|H| = 2) and 17.6–19.6 (|H| = 4). Crater E[2] against floor and model C4: same degree and matrix, multiplications within 6 of 34,158 at B = 64.

#figure(image("e4_ops_vs_degree.svg", width: 92%), caption: [Every E4 system: F#sub[p] row-reduction multiplications against the factor-base polynomial degree B/|H|. 1,128 systems, 72 curves; records `e4/records.jsonl`.])

#figure(image("e4_gap_b64.svg", width: 92%), caption: [B = 64, |H| = 4: crater with rational E[2] against floor and floor model with C4, mean of 2 random targets per bar; the bars are equal to within 6 multiplications.])

Raw records: `research/isogeny_conductor_gap_round3_20261011/e4/records.jsonl`; tables `TABLES.md`, `tables.json`; plot source `plot.py`.
