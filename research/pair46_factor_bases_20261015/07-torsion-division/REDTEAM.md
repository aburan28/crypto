# 07-torsion-division: red-team review

2026-10-15. Reviewer: independent red-team subagent. No other file in this directory was changed.

**Outcome: not refuted.** The verdict (exact obstruction) and the headline numbers are supported. The issues below are wording and scope problems. None of them changes the verdict.

## Checks that pass

- **Preregistration.** PROTOCOL.md mtime is 16:04:44, src/ 16:07:13 and results.json 16:07:44. The sha256 matches PROTOCOL.sha256 (97f57630...8f62b6), and results.json records the same hash. All arms, sizes and seeds ran as preregistered, and every prediction held.
- **Reproduction.** Rebuilt in a scratch target dir and re-ran on 1 core (`taskset -c 3`) in 10.3 s. results.json is byte-identical in content once the `secs`/`wall_secs` timing fields are removed.
- **Controls and instrument.** Every arm is scored with the same exact coverage enumeration (`coverage()` over discrete logs). RND is the right null for SX, and RND - SX has CIs that contain 0. The KL positive control works: coverage falls to about 0.037 against 0.41. TOR is a working positive control for the root-count detector: 27/27 curves matched the prediction, and q-subgroup coverage was 1/q on all 27.
- **Pair-46 arithmetic.** I checked independently that the twist order is 2(p+1) - N = 5947698065141045 = 5 v^2. Also v mod 5 = 3, so the 5-part is cyclic Z/5 and psi_5 has exactly 2 roots. Neither 5 v^2 nor N is divisible by 2, 3, 7, 11 or 13. The 20/20 killed-by-5v result on the crater forces Z/5v x Z/v: exponent | 5v and order 5v^2 leave no other two-generator group. This holds with high probability, not with certainty.
- **E1-E4.** The math is correct. In E4 the polynomial phi_m - s psi_m^2 can have no F_p-root on the twist, because s is on E. This agrees with off_E_hits = 0.

## Issues (minor)

1. **Number mismatch.** RESULTS.md (rho table) says the level sets have "mean size 0.58-0.62 points", and the summary says "about 0.6". results.json gives per-size means over 5 curves of 0.47-0.56: 12 bits 0.50/0.50/0.53, 14 bits 0.56/0.47/0.50, 16 bits 0.54/0.53/0.52 for m = 3/5/7. The 0.58-0.62 range matches only curve (12, seed 0). The rho values quoted (0.124-0.140, etc.) are correct, and the true sizes sit even closer to the predicted 0.5.
2. **Wrong justification in E2.** RESULTS.md says m in {2,3,7,11,13} have 0 roots "since f has no root because both group orders are odd". That argument covers only m = 2. For odd m the reason is that m divides neither N nor 5 v^2. The checked numbers are right either way.
3. **"Consequence" sentence is trivial or overstated.** rho <= 1 holds for any membership polynomial, since a polynomial with |V| distinct roots has degree >= |V|. So "no base has membership degree below its size" is not specific to division polynomials. Also, "the best case is the transported base, at ratio 1/m^2" contradicts the paper's own table: TOR_3 has rho = 0.25 (equal to PRE_2), and LVL_3 has rho of about 0.13, above PRE_3's 0.111. The finding that holds is narrower: no tested torsion or division-polynomial arm beats an explicit list on degree or coverage.
4. **Scope of the metric.** rho counts degree only. It does not see sparsity or the low-degree Weil-restriction structure that makes summation-polynomial bases tractable (for example Diem/Petit subspace bases, which also have rho = 1). "Give no useful factor base" should be scoped to the arms and metrics M1-M4. The exact parts (E1-E4) carry over to pair 46. The level-set rho and coverage figures are toy 12-16-bit only, and RESULTS.md labels them that way correctly.
5. **M4 is a consistency check, not an independent solve.** roots_e comes from evaluating g at every x in F_p (exhaustive), and the log table is the ground truth. 400/400 agreement shows that g and the univariate route are consistent. It is not Groebner evidence, and RESULTS.md does not claim that it is.
