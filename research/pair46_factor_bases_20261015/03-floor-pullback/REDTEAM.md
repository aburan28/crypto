# REDTEAM — 03-floor-pullback

Reviewer: independent red-team pass, 2026-10-15. No other file in this
directory was modified.

## Verdict

Not refuted. The exact obstruction (phi: E(F_p) -> E'(F_p) is a group
isomorphism when N is prime, so Y_floor(phi(S)) = Y_crater(S)) is correct, and
the empirical numbers match the records. The issues below are minor.

## Reproduction

Both runs were re-executed on one core (`taskset -c 0`) into the scratchpad,
not this directory. Main run (v = 359, 409, 443, about 12 s) and extension E1
(v = 1213 to 1613, about 3 min) gave instance records identical to
`results.json` and `results_ext.json` in every field except the timing fields.
The current script, edited at 15:49 after the main run, therefore reproduces
the main run too.

## Checks

- Protocol order. File mtimes: PROTOCOL.md 15:26, then results.json 15:38,
  then PROTOCOL_EXT.md 15:39, then the crash log 15:41, then
  results_ext.json 15:52. The directory is untracked in git, so this order
  rests on mtimes only. Arms, metric, seeds, m, the C0 null and the decision
  rule in the code match PROTOCOL.md. D1, [5]P in place of [5v]P, is a
  necessary fix and is disclosed.
- Instrument. One `floor_yield` (one log table built from G' = phi(G)) is
  used on A1, A2, C1, C0 and P+.
- Numbers. The point ratios are A1 0.99778 to 1.00144 and A2 0.99964 to
  1.00226, so the headline "0.9978 to 1.0023" is right. Max |z| is 1.40 for
  A1, 1.28 for A2 and 1.41 for C1, and min |z| for P+ is 562. All match.
  `summary.json` `ratio_ranges` holds the ends of the 95% CIs, not the point
  ratios. RESULTS.md says so.
- Math. With v = t + 2 and p = 5v^2 - v + 1: N = 5v^2 - 2v + 4, the twist
  order is 5v^2, and disc Z[pi] = -19 v^2 (conductor v). The crater twist is
  O_K/(v w), which is Z/5v x Z/v. pi = v w - 1 = -1 on E[v], so all v+1
  v-isogenies are rational. "Some [5v]P != O on a twist of order 5v^2" does
  imply that the floor twist is cyclic. All correct.

## Minor issues (none changes the verdict)

1. The X0 check is not independent evidence. The floor table is generated
   from phi(G), and phi was already checked to be a homomorphism, so
   dlog_floor(phi(P)) = dlog_crater(P) holds by construction. "Two
   independent log tables, bit-identical" only confirms that the
   homomorphism holds on S. It does not corroborate the theorem separately.
2. The cost figure "about 6e3 field ops" for the dual isogeny at pair 46 is
   sqrt(v) with the constant and log factors dropped. A practical sqrt-Velu
   evaluation at l = 3.4e7 costs roughly an order of magnitude more. This
   strengthens the "costs more" conclusion, but the number is not a measured
   or derived constant.
3. The positive control (an interval in log space, ratio about 0.002 to
   0.006) shows that the detector sees gross structure. It does not show
   sensitivity to small effects. The minimum detectable shift at 3 sd0 is
   about 0.5% of mu0 at v about 400 and about 0.1% at v about 1600. Because the
   result is exact, this matters only for C1, the native floor base, which is
   empirical only.
4. The script never asserts that [5]P has order exactly v. This is checked
   only implicitly: the Velu homomorphism check on 20 pairs would fail
   otherwise. "[5v]P = O" was sampled only at v = 359, which RESULTS.md
   states.
5. RESULTS.md says the runs used `taskset -c 3`, but the reproduce line says
   `-c 0`. This is cosmetic.
6. Scope. The exact claim covers transported bases only. Native floor bases
   are tested only for C1, with 2-term decompositions at m about sqrt(N), and
   RESULTS.md states that limit. Nothing toy-scale is extrapolated as an
   exponent: the pair-46 statements are either the exact identity or derived
   costs, and both are labelled as such.
