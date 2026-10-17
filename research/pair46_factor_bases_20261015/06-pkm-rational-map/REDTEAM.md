# 06-pkm-rational-map: red-team review (2026-10-15)

Verdict: **not refuted**. The negative verdict is supported, with the minor corrections and caveats below.
No other file in this directory was changed.

## Re-run
- `analyze.py` on the committed `results.json`, run in a scratch copy, gives a `summary.json` byte-identical to the committed one.
- `taskset -c 3 python3 run_pkm.py 4` (n1 = 2..4, 95 s), also in a scratch copy, gives the same results as `results.json` for every arm, A8 included: valid counts, nbsols, sugar D_corr and lcm Dl_corr are identical target by target. Mean T is within about 25%.

## Checks that passed
- PROTOCOL.md's hash matches PROTOCOL.sha256 (baa68336…). The timestamps (16:03) come before results.json (16:29). Amendment 1 is disclosed.
- The numbers in RESULTS.md match summary.json:
  - slopes 2.02 [2.00, 2.04], 1.91, 2.11 and flat 1.28-1.31;
  - rho 0.52; PKM Table 4 slope 1.51;
  - IC/rho: A4 at n1 = 5 is 4.0e6 [2.4e6, 7.6e6], A1 is 7.0e6 [4.1e6, 1.7e7], A7 at n1 = 7 is 5.8e6;
  - 42 yield intervals, and only A7 at n1 = 5 excludes 1;
  - the A8 ratios; and 168 = 7 x 2 x 12 sympy checks, all matching.
- X1: A1 = A2 and A4 = A5 target by target, at every level of the log.
- X3: the factorisations of p-1 and p+1 are correct (sympy).
  - The degree <= 7 subgroup is 1134 ≈ 2^10.15, so m ≈ 5.2.
  - The 21061-smooth part of p+1 is ≈ 2^21.8.
  - The extrapolation (2^30.4)^1.51 ≈ 2^45.9 checks out.
- The cost model is sound. The yield per system is ~2|F|^2/N = Θ(1) and about p^(1/2) relations are needed, so IC ≈ p^(1/2)·T. To beat rho, T must stop growing (slope ≤ 0).
- The last-3-level slopes (1.97, 1.92, 2.01) agree with the 4-point fits, so n1 = 2 overhead is not driving the slope.

## Corrections and caveats (do not change the verdict)
1. **X3 wording is wrong in a small way.** "The only subgroup of size >= p^(1/2) has order 2·3^4·7·241·631" is false: 3^4·7·241·631 ≈ 2^26.36 also exceeds 2^26.2, and so do subgroups that contain 34489703. The conclusion still holds: every such subgroup needs a chain step of degree 241 and 631, or of degree 34489703.
2. **The headline "p^2.0 (95% CI 2.00-2.04)" is a local slope over p = 17..1153 (4 points).** The bootstrap CI covers only resampling of targets, not model or range uncertainty. Along this family the chained system has 2·n1 + 1 variables and p ≈ 4^n1, so the "exponent in p" is really cost growth per level in n1. It is not an asymptotic exponent. RESULTS labels it empirical and toy, which is adequate, but the headline CI reads as more precise than it is.
3. **The degree positive control that "passed" is post hoc.** The lcm metric was added after a smoke test at n1 <= 3 showed the preregistered sugar metric missing PKM (this is disclosed). The metric also has limited power: degree-2 chain equations give an lcm of about 4 almost by construction. It confirms the system construction, not solver equivalence with Magma F4.
4. **The yield detector's preregistered power (>=10x) is shown only at n1 = 6, 7, where the chained arms were not run.** "No yield effect" for chained bases rests on X1, which is exact, not on M1.
5. **The PKM Table 4 times (1.14 … 4975 s) were not re-checked against the paper in this review.** They are consistent with the "about 8x per n1" written into the protocol before data.
6. **The absolute IC/rho ratios (e.g. 4e6x) are Python-constant-laden, as RESULTS itself says.** Only the exponent gap carries over. The headline's "4e6x rho" should be read the same way.
