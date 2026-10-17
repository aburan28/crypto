# REDTEAM — 02-deuring-lift

Reviewer: independent red-team subagent, 2026-10-15. No other file in this directory was altered.

## Verdict: supported (not refuted), with minor caveats

## Checks performed

- **Preregistration.** PROTOCOL.md mtime 15:17:58; run_deuring_lift.py 15:21:56; results.json 15:21:59. The protocol predates the reported data. The precision change (realprecision 80) is an instrument fix, applied to every arm, and the frozen 1e-20 tolerance is unchanged. It is disclosed.
- **Reproduction.** I re-ran `python3 run_deuring_lift.py` into a scratch file (2.2 s, 1 core). Apart from the per-twist `sec` timing fields, the output is **identical** to the committed results.json.
- **Numbers vs claims.** All of these match results.json/summary.json:
  - a_p = 34489701
  - N = 5947697996161643 (= p+1-t)
  - isomorphism checks pass
  - ellrank [1,1] with generator (0,9), trivial torsion, h = 0.181292
  - ord(G mod p) = N
  - 115 twists, rank census 51/54/8/2, no analytic-rank fallbacks, no undetermined twists. The rank-3 twists are d = 69 and 104.
  - A2/NC/PC medians, bootstrap CI and Mann-Whitney U = 321 (p = 0.50)
  - Xedni 0/96 (all p-integral and on the curve), PC Xedni 24/24
- **Math.** These are correct:
  - (0,9) lies on y^2+y = x^3-38x+90.
  - E(Q(sqrt-19)) has Z-rank 1+1 = 2, because E^(-19) is isogenous to E.
  - Points on twists with distinct characters chi_d lie in distinct eigenspaces of E(L)⊗Q, so they are independent.
  - The lifts of T in E(Q) are (k+jN)G, so producing one requires log T.

## Caveats (none overturns the verdict)

1. **Expected lift height overstated by about 4x.** "h·N^2/3" takes k uniform on [0,N). The minimal-height lift uses |k| <= N/2, which gives about h·N^2/12 ≈ 5.3e29 (log2 ≈ 98.8) rather than 2.1e30. The conclusion does not change.
2. **A2 "0 Q-relations" comes from the derivation, not from I-Q.** The script checks that each twist height is nonzero and relies on the block-diagonal argument. It never builds an ellheightmatrix for A2. The NC arm draws uniform logs directly and does not call elllog on points. That is equivalent for prime N, but "same instrument on every arm" is slightly stronger than what was done.
3. **The positive controls are gross and not independent.** PC_lifted_relations is one prime-independent Q computation (the Gram matrix of G, 2G, 3G) repeated 24 times. "48 found" is therefore one check counted 24 times. The I-p PC plants G+2G-3G = 0, a vector of length sqrt3. At 16 bits it gives rho 0.46, against the threshold 0.5.
4. **Power and equivalence.** The PC does not show that the instrument could detect a modest bias. The CI half-width (about 0.22) is larger than the 0.1 equivalence margin. "Indistinguishable" passes only through the "CI contains 0" clause, so it shows absence of evidence, not equivalence. The point estimate (+0.14) points away from the hypothesis (A2 relations are longer). The obstruction rests on the exact eigenspace and DLP-circularity arguments, not on this statistic.
5. **Out of scope.** The Xedni arm is one trial per (prime, r), with r <= 5. The GLV / Q(sqrt-19) item is a derivation, and RESULTS.md labels it correctly. Nothing about the toy results is extrapolated to pair 46. The pair-46 statements are exact.
