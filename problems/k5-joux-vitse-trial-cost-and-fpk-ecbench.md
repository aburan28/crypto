# k = 5 Joux–Vitse decomposition: the trial constant puts rho parity near q ≈ 2^55, and ecbench has no `F_{p^k}` instance to price it

- **Status:** `open`
- **Kind:** empirical bottleneck (the per-trial Gröbner constant); implementation gap (no `F_{p^k}` curve family in ecbench on `main`; extension-field curves not harvested by the curve registry)
- **Affected work:** [`research/joux_vitse_n1_20261015/`](../research/joux_vitse_n1_20261015/README.md); [`RESEARCH_K5_TORSION_JOUX_VITSE.md`](../research/notes/index-calculus/RESEARCH_K5_TORSION_JOUX_VITSE.md) §5–6; PR #1863 / #1792 (the `F_{p⁴}` ecbench instance, unmerged)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-15

## Problem

Joux–Vitse `(n − 1)`-decomposition on `E(F_{q⁵})` cannot approach rho at
sizes this repository runs, because one decomposition trial costs
`≈ 5–9·10⁹` `F_p` multiplications (measured, flat in `q`) while about
`24 q` trials buy one relation and `q/2` relations are needed: the method
costs `≈ 12 q² · C_trial` against rho's `≈ 0.89 q^{5/2} c_add`, a ratio
`√(q*/q)` with `q* ≈ 2^{55}`.  Separately, the requested ecbench run with
full pricing cannot be made: `main`'s ecbench `Instance` has prime, binary
and wide-Koblitz variants only (`Instance::Fp4` lives in the unmerged
#1792), and `scripts/build_curve_registry.py` harvests no extension-field
JSON, so the `F_{p⁵}` curves of the complete DLPs have computed but
unregistered ICV1 identities (`research/joux_vitse_n1_20261015/curves.json`).

## Evidence and uncertainty

- **Observed:** complete DLPs at `p = 11` (three seeds) and `p = 31` (one
  seed), every relation re-added in the group and every trial's
  decompositions equal to an independent meet-in-the-middle oracle, at
  IC / rho `≈ 4·10⁷–7·10⁷` against rho run on the same target in the same
  unit (`research/joux_vitse_n1_20261015/summary.json`); the trial-cost
  ladder `p = 11 … 241` (`data/ladder/`); the committed `jv_quintic` test
  on today's engine, `C″ = 1.104·10¹⁰` (`data/jv_quintic_rerun/`).
- **Derived or proved:** the counting bound on the rate
  (`16·C(|F| + 3, 4)/#E`) and the crossover formula
  `q* = (κ C_trial / (2·√(π/4)·c_add))²` (MODEL.md).  No lower bound on
  `C_trial` is proved.
- **Suspected explanation:** the descended system's certificate of
  inconsistency needs F4 to degree 19 on a `≈ 2,000 × 3,700` matrix for
  five generic octics in four unknowns; nothing about the curve lowers it
  without a symmetry (the 2-torsion Edwards variant reaches `7.3·10⁷`,
  RESEARCH_K5 §3).
- **Missing capability:** an ecbench `F_{p^k}` instance with a counted
  group and a full-price table; registry harvesting for extension curves.
- **What this evidence does not establish:** it does not establish that
  no `k = 5` route reaches parity at feasible sizes (covers do on weak
  classes; symmetries change the constant); the crossover is an
  extrapolation on measured constants and derived exponents.

## What would unblock it

- [ ] A trial cost `C_trial ≤ 10⁵` `F_p` multiplications on the plain
  system (parity near `q ≈ 2^{21}`, group `≈ 2^{105}`), or a route with a
  different exponent — measured with the same MITM cross-check.
- [ ] An ecbench `F_{p^k}` instance (`Instance` variant, counted group,
  `full_prices.json` row) on `main`, so that a `k = 5` complete DLP can be
  run as an `ecbench.spec/v1` session against `rho.strong` and audited
  with replays.
- [ ] Extension-field harvesting in `build_curve_registry.py`, so the
  study's curves can be cited by slug.

The first is what parity needs; the second and third only make the
existing measurement auditable in the harness.

## Next discriminating step

- **Action:** run the 2-torsion Edwards test (`jv_quintic_edwards`) under
  F4 trace replay and price its end-to-end DLP at `p ≈ 271` (RESEARCH_K5
  §6 sized it at 27 h on four cores without replay).
- **Inputs and controls:** the same pair-table oracle; rho on the same
  target; the plain route of this study as the baseline.
- **Outcomes:** a replay factor of 3–10× moves the extrapolated crossover
  by 16–33 bits (RESEARCH_K5 §5), still past `2^{200}`; a larger factor
  would reopen the route.
- **Recheck trigger:** #1792 merging (ecbench `F_{p⁴}`), or a new symmetry
  for `E(F_{p⁵})`.

## Resolution and history

- 2026-10-15: recorded with the study's complete DLPs and ladder.
  Unresolved.
