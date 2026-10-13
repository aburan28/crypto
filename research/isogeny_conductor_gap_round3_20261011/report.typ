#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Round 3 of the conductor-gap study: two-torsion summation polynomials, and the combined descent audited
_2026-10-11. research/isogeny_conductor_gap_round3_20261011. Protocol frozen before both leads ran (PROTOCOL.md). Result classes: E4 stage diagnostic (no whole-method or rho comparison); E5 engineering, an audited negative._

== Findings
+ *Two-torsion symmetrisation is a quotient isogeny and pays as predicted.* Writing the S#sub[3] system for an H-stable factor base on E/H (Vélu) divides the factor-base polynomial's degree by |H|. At matched base B = 32 or 64 on E, the F4 solving degree falls 1.9–3.7× and F#sub[p] row-reduction work 4.1–4.3× (|H| = 2) and 17.6–19.6× (|H| = 4). F4's cost on these systems depends only on that degree.
+ *The crater's richer 2-torsion changes nothing at matched |H|.* On 12 pairs (every D = −7 pair and the ℓ = 2 pairs of D = −11 and D = −3) the crater has rational E[2] and the floor rank-1 2-torsion; crater E[2] and the floor's and floor model's C4 give the same solving degree and widest matrix, with work within 6 of 34,158 at B = 64.
+ *Round 1's combined Weil descent is negative in audited sessions.* On the a = 1 Koblitz curve over F#sub[2#super[17]] it verifies in 4 of 4 runs at 324–613 s online (blind 0.9–2.0 s), maximum solving degree 4 (blind 3), 234–529× the blind F4 work paired by run, for at most h = 2 times the relations. The predicted 600–700× fails; round 1's unaudited 629 and 654× lie above the audited range.

== E4: 1,128 F4 systems on 72 curves
24 even-cofactor round-1 prime pairs × crater, floor, floor model; B = 32, 64; H = 1, C2, C4, E[2] where rational; two random and one planted target. Every system completed; F4 agrees with enumeration on all; every planted solution annihilates its basis; the Vélu map passed 32 random-sum checks per curve and subgroup. Run report: reports/run-20261011T1714Z.

#table(columns: 5, align: right,
  [base-poly degree], [solving degree], [widest matrix], [mults, random], [mults, planted],
  [8], [10], [45], [7,408–7,414], [13,903–13,929],
  [16], [18], [93], [34,152–34,158], [81,128–81,178],
  [32], [34], [189], [145,630], [483,552–483,561],
  [64], [66], [381], [600,510], [3,123,348–3,123,354])

#figure(image("reports/run-20261011T1714Z/e4_ops_vs_degree.svg", width: 88%), caption: [Every E4 system: F#sub[p] row-reduction multiplications against the factor-base polynomial degree B/|H|; records e4/records.jsonl.])

#figure(image("reports/run-20261011T1714Z/e4_gap_b64.svg", width: 88%), caption: [B = 64, |H| = 4: crater E[2] against floor and floor-model C4 on the 12 pairs whose crater has rational E[2].])

The protocol expected rational E[2] on the crater only from ℓ = 2 volcanoes; the run found it on every D = −7 crater too (their floors' Frobenius conductors are even). The predictions did not depend on this. Scope: two summands, bivariate systems, split factor-base polynomials; higher m is untested.

== E5: two audited ecbench sessions
#table(columns: 6, align: right,
  [session, round], [blind wall, F4 ops], [combined wall, F4 ops], [ratio], [max degree], [peak RSS],
  [t0, 0], [0.9 s, 65,072,207], [324.5 s, 32,012,560,285], [492], [4], [0.66 GB],
  [t0, 1], [1.8 s, 133,408,080], [613.2 s, 70,613,284,746], [529], [4], [0.64 GB],
  [t1, 0], [2.0 s, 164,172,136], [453.5 s, 38,387,092,945], [234], [4], [0.65 GB],
  [t1, 1], [1.3 s, 84,877,444], [464.5 s, 34,853,985,256], [411], [4], [0.99 GB])
Sessions ECBS1h05f404b1f89b and ECBS1he4c2c379c257, 8 records each, all verified, all replays reproduced; isolation earned L1 on 6 and L2 on 2 records each. Solver work depends on each run's algorithm seed: round 1's two identical counts (59,916,380,641) came from two warm-ups with the same seed. S leaves solver XORs unpriced and is a loose lower bound for the combined arm. The h = 6 random-curve timeout of round 1 was not rerun.

#figure(image("reports/run-20261011T1727Z/e5_arms.svg", width: 88%), caption: [Session ECBS1h05f404b1f89b (target 0, seed 20261010): online wall and F4 operations per arm.])
#figure(image("reports/run-20261011T1758Z/e5_arms.svg", width: 88%), caption: [Session ECBS1he4c2c379c257 (target 0, seed 20261011).])

== Graphs checked
Scoreboard: a round-3 panel added and round 1's combined-descent row updated to the audited figures. Progress timeline: unchanged (E4 makes no rho comparison; E5's curve was already plotted by round 1). Built from the round-2 head 75b95656a because main at 1613b90fe does not compile (union-merge damage in koblitz_index_calculus.rs and koblitz_factor_base_search.rs).
