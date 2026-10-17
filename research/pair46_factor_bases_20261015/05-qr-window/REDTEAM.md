# 05-qr-window: red-team review (2026-10-15)

**Verdict: the result stands (not refuted).** The `constant_only` verdict and the numbers are supported.
Two caveats about how strong the positive control is are listed below.

## Checks done

- **Protocol before data.** PROTOCOL.md mtime 15:53:30, results/run.log 15:57:27. Its sha256 matches PROTOCOL.sha256
  (3ee2d891...8928). Arms A0-A4, sizes, seeds, bootstrap count and decision rule all match `src/main.rs`.
  No amendment was needed.
- **Reproduction.** I rebuilt in a scratch copy and ran `taskset -c 3 qr_window . 40 100 100 2000`.
  results.json matches the committed file exactly once `wall_s` is excluded. run.log differs only in timing lines.
  The input CSV sha256 is 88c26329...7c12, as stated.
- **Numbers.** Every figure in RESULTS.md and in the claim matches run.log and results.json: G values, CIs, counts,
  toy ratios, 36/36 E1 and 8000/8000 E2.
- **Math.** p = 5t^2+19t+19 with t = 34489701, and N = p+1-t. The twist order 2p+2-N = 5947698065141045 = 5*34489703^2.
  (N-1)/2 + (N'-1)/2 = p, so f has no root. The global gain 2p/(N-1) = 2 + 2(p-N+1)/(N-1) = 2.0000000116.
  The Polya-Vinogradov scale sqrt(p) ln p is about 2.8e9, which is more than B (about 1.09e8), so calling the window
  gain empirical is correct.
- **Side finding.** Confirmed: `research/factor_base_lab.md` line 155 reads `V = {x < B : Legendre(x) in {0,1}}`.
  `factor_base.py` `small_x` keeps only x where `lift_x` succeeds, so "gain 1 against existing liftable factor bases"
  is correct.
- **Overstatement check.** The verdict is `constant_only`, a factor of 2 in candidate count, and no exponent is claimed.
  The pair-46 density was measured directly, not extrapolated from toy curves. The relation identity rests on exact E2,
  not on toy extrapolation.

## Caveats (they do not refute the result)

1. **The positive control is true by construction.** `count_window` and `toy_measure` decide "is an E-point" with
   `chif != -1`. That is the same `jacobi(E.rhs(x), p)` value the A1 and A4 membership predicates use. So
   G(A1) = 1/M1(A0), M1(A4) = 0 and the per-target E2 check `c[0] == c[1]` all follow from the code itself.
   A broken `jacobi` would not show up in these checks. Two checks are independent: E1 full-field counts against N,
   with N*P = O verified, and the A2/A3 halving controls. The first validates `jacobi` against an externally fixed
   group order on 36 toy curves. The second shows that the M2/D relation instrument registers lost FB points.
   "Positive control detected" therefore means the bookkeeping is consistent and the relation counter responds to
   thinning. It does not show that the instrument can detect a non-trivial density effect.
2. **The resolution sentence goes slightly too far.** RESULTS.md says "resolves a 2x density effect to about ±0.02%...
   later items can read a null as a real null". The ±0.02% is the Wilson sampling width of the liftable fraction on a
   1.09e8 window. It is not a calibrated sensitivity to a non-trivial factor-base effect, so later items should cite
   the A2/A3 relation-loss controls as their sensitivity evidence, not A1.

Minor: the D count adds both R-P and R+P with `x2 >= x1`, so the x1 == x2 case can count twice. This is the same in
every arm and does not affect the ratios.
