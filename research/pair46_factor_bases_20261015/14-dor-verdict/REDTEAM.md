# 14-dor-verdict: red-team review

Reviewer: independent red-team subagent, 2026-10-15. Only this file was added; no other file in the directory was changed.

## Re-run

- I rebuilt `dorv/` from a copy in a scratch directory, then replayed `m2_n13_d5_s1`, `m2_n11_d5_s2` and `m2_n9_d3_s1` against the raw campaign files. The outputs are identical to `out/replay/*.json` except for the `source` path string. Each took less than 1 s.
- I wrote my own comparison of `out/replay/m2_*.json` against `raw/m2_*.json`. Target, `d_star` and `rank_profile` match on 4638/4638 targets. **The P1 replay gate is confirmed.**

## What holds

- **Protocol.** PROTOCOL.md (mtime 16:57Z) predates `out/started` (17:03Z) and `dorv/` (17:01Z). The deviations are disclosed. The null scope was restricted to n<=15 as preregistered.
- **Instrument.** One binary is used on every replayed arm.
- **Positive controls.** P2 is +1.92 [1.83, 2.00]. P2b, the planted-q control, is monotone. Both work.
- **E2.** The subspace counts are the correct Gaussian binomials:
  - [7,2]=2667, [7,3]=11811, [11,2]=698027, [11,3]=50955971, [9,3]=788035.
  - The 73 subspaces at n=9 with dim(V·V)=3 are exactly the 73 = 511/7 scalings of F_8.
  - The theorem (Hou–Leung–Xiang, applied with no intermediate subfields for prime n) is applied correctly.
- **Numbers.** The H1, H5, N1 and N2 numbers in RESULTS.md match results.json:
  - 40/45 cells within ±0.25;
  - nb_window − random has |diff| ≤ 0.0625 in all 31 cells;
  - standard −1.53 at n=23, n'=8.

## Problems

1. **E1 is false as stated (an EXACT claim).**
   - **The claim.** PROTOCOL and RESULTS say the *degree-2 Macaulay row space* lies in span{1}+L+Q, so `rank_2 <= 1+l+q`, "checked on 4638/4638".
   - **What the code checks.** It checks `rank_eq`, the rank of the equation list itself (`main.rs:76`), not `rank_profile["2"]`.
   - **Why they differ.** `refutation_degree` multiplies each equation by monomials up to `d − actual degree`. Coordinate equations that happen to be purely linear are therefore multiplied by variables at D=2.
   - **What the data show.** In the replay outputs, `rank_profile["2"] > 1+l+q` on **937/4638** targets, and `rank_profile["2"] > rank_eq` on 940. The violations are concentrated at small n' and in subfield arms. Examples:
     - n=11 n'=2: 45/72 standard targets;
     - n=15 n'=3 subfield: 61/71 targets.
   - **What survives.** The l and q identities, and the bound on the equation span, hold.
   - **Consequence.** The statement "the structure label enters degree 2 only through l and r" is not established for the instrument actually used. Which coordinate equations are linear depends on the coordinates, not only on (l, r). The N2 null almost never produces individually linear equations.
2. **The headline slope range is misstated.**
   - **The claim.** The headline says "1.07 to 1.36 degrees per doubling for generic bases", and the PDF says "at least one degree per doubling".
   - **The data.** n=23 random is 0.87 [0.86, 0.88] and n=23 nb_window is 0.86 [0.85, 0.88].
   - **Correction.** The generic range is 0.86 to 1.36, and the rule fires in 14/16 groups. RESULTS.md's table discloses this, but the headline and the report do not.
   - Standard never fires: its slope is 0.17 to 0.98.
3. **The step from "r is bounded" to "no structure beats standard" is empirical, not exact.**
   - The verdict `exact_obstruction` overstates it. E2 bounds r exactly. That D* is a function of (l, r) is the N2 test, and the N2 test has these weaknesses:
     - a point-estimate ±0.25 rule with no equivalence test;
     - about 17 of the 45 cells sit at the D*=2 floor (both means ≤ 2.1), where agreement is nearly automatic;
     - all 5 failures are at n'=5, the n'≈n/2 regime that relation yield needs.
   - In the failures the Semaev system is easier than its invariant-matched twin, by up to 0.9 degree. So D* has a Semaev-specific component that the model does not capture. A base not tested here could in principle exploit it.
   - RESULTS.md's "it only helps the incumbent standard window" is wrong. random and nb_window at n=13 n'=5 are also 0.88 and 0.90 below the null, although standard remains the lowest in absolute D*.
4. **A post-hoc predictor became the headline.**
   - r(D*, (l+r)/n) = 0.88 [0.85, 0.91] was not preregistered. H5 named nvars and n. Yet it is presented in the headline as confirmatory.
   - Its CI overlaps that of r(D*, n') = 0.82 [0.75, 0.88], so "explains more than n' alone" is not shown.
   - The cell-level r uses the maximum over seeds.
5. **Scope.** Toy fields with n ≤ 23 and m=2 only. The non-transfer to pair 46 is correctly stated.

## Corrected verdict

The verdict should be **empirical negative for dimension and structure levers in the m=2 toy model, with an exact lower bound on r for prime n**, not an exact obstruction. The supported parts are:

- The abandon rule fires for generic bases in 14/16 groups, with slopes 0.86 to 1.36.
- The standard window attains the E2 minimum of r.
- An (l, q)-matched null reproduces 40/45 cells. Its failures are one-sided, in the relevant n'≈n/2 regime, and favour Semaev.

E1 should be restated as a bound on the rank of the equation span, not on the degree-2 Macaulay rank.
