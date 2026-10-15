# 14-dor-verdict: results

**Campaign judged:** `research/factor_base_dor_campaign_20261013`, a binary Koblitz toy study (`y^2+xy=x^3+1` over `F_{2^n}`, m = 2 Semaev `S_3`, Weil descent, GF(2) Macaulay refutation degree `D*` on non-decomposable targets).
**Protocol:** `PROTOCOL.md`, frozen 2026-10-15T17:00Z before any data.
**Run:** 2026-10-15, 17:03Z to 17:15Z, on one core (`taskset -c 0`). Measured compute was about 2 min. A null stage that ran past the preregistered scope was killed after about 9 min (see Deviations).
**Reproduce:** `bash research/pair46_factor_bases_20261015/14-dor-verdict/run.sh && python3 research/pair46_factor_bases_20261015/14-dor-verdict/plot.py`. This needs a Rust toolchain only. Sage is not used.

## Verdict

**Factor-base dimension changes D\*.** For generic bases the increase is at least one degree per doubling of n'. **Structure changes D\* at matched n', but almost entirely through two exact invariants of the base:**

- `l = dim(V1+V2)`;
- `r = dim(V1·V2)`, the dimension of the product space.

When n is prime (the ECC2K-130 / n = 131 analogue), `r >= min(n, 2n'-1)` holds for every subspace. The plain polynomial window ("standard") already reaches this bound. So in the model the campaign measured, no factor-base structure lowers D\* below what the standard window already gets. The subfield closure that gives `r = n'` cannot exist when n is prime. This does not give a route to a drastic speedup.

Some structure is left over beyond (l, r): in 5 of 45 cells the standard window and n = 13 / n' = 5 sit 0.4 to 0.9 degree below the matched null. It only helps the incumbent standard window, and its cause is not established.

## What ran (replay gate P1: PASS)

The sage engine cannot run on this host, so its m = 2 path was ported line for line to Rust (`dorv/`, no dependencies). The port includes CPython's MT19937 so the recorded targets and bases can be regenerated. Sage's Conway modulus was confirmed with the recorded n = 9 subfield basis (`w = z^73`, giving `[0x1, 0x14c, 0x150]`).

| check | measured |
|---|---|
| recorded m=2 cells replayed (R1 n≤15, R3 n=17, R4 n=19 n'≤7) | 62 cells, 239 structure arms |
| bases, target lists, decomposable fraction, factor-base point count identical | 239 / 239 arms |
| per-target `d_star` and `rank_profile` identical | 4638 / 4638 targets |
| per-target `macaulay_shape` identical | 4255 / 4255 where recorded (the n=19 raw files do not record shapes) |
| δ-sweeps n=11, 13 replayed: mean D\* and fraction identical | 6 / 6 rows |

All other numbers come from the recorded raw files, which this replay validated. Cells that were not replayed (n = 19 with n' = 8, 9; n = 23; the δ-sweeps for n ≥ 15; sym4) are recorded values from the same engine.

## Exact statements (C-exact)

- **E1 (proved, and checked on 4638/4638 targets).** In the m = 2 system the linear part comes only from `x3^2 (X1+X2)^2` and the quadratic part only from `x3·B + B^2` with `B = X1X2`. Every quadratic monomial has the bilinear form `a_i b_j`.
  - The linear forms have rank `l = dim(V1+V2)`.
  - The quadratic forms have rank `q = dim M_x3(V1·V2) = r - [x3 ∈ V1·V2]`. Here `M_x3(y) = x3·y + y^2` has kernel `{0, x3}`.
  - So the degree-2 row space satisfies `rank_2 <= 1 + l + q`.
  - The structure label enters degree 2 only through l and r.
- **E2 (theorem, checked).** For prime n, `dim(V·V) >= min(n, 2d-1)` for every d-dimensional F2-subspace. This is the linear Cauchy–Davenport bound (Hou–Leung–Xiang 2002; Bachoc–Serra–Zémor). `<1, z, …, z^{d-1}>` attains it.
  - Exhaustive checks: n = 7 (d = 2, 3: 2667 and 11811 subspaces), n = 11 (d = 2, 3: 698 027 and 50 955 971 subspaces). In every case the minimum is exactly 2d−1 and the standard window attains it.
  - n = 13 (d = 3 to 6) was checked by sampling 200 000 subspaces per d; no sample fell below 2d−1.
  - Contrast: n = 9 is composite, and d = 3 reaches 3 < 5 through the subfield F_8.
- Consequence: for symmetric bases with n prime, `l + r >= d + min(n, 2d-1)`. With 2d ≈ n, which relation yield at m = 2 requires, this is above n.

## Hypotheses (C-empirical)

Every interval below is a 95% bootstrap with 10 000 resamples, using seed clusters.

| | rule (PROTOCOL) | result | adjudication |
|---|---|---|---|
| H1 closure | subfield − random < 0 | n=9 n'=3: **−0.57 [−0.73, −0.41]**; n=15 n'=5: **−0.99 [−1.00, −0.94]**; n=15 n'=3: −0.01 [−0.06, 0] | confirmed where `r` drops (r = n' vs 6/15) |
| H1 windows | nb_window − random within ±0.25 | all 31 matched cells: diff ∈ [−0.06, +0.06] | confirmed (equal l, r) |
| (nb_break) | — | +0.0 to +1.0 above random (l doubles, r up to n'^2) | worse |
| (standard) | — | down to **−1.53 [−1.72, −1.34]** (n=23 n'=8); r = 2n'−1 vs generic | best, and at the E2 minimum |
| N1 null | random seed-vs-seed | 58 contrasts, 95% of \|diff\| ≤ 0.15 | structure effects above exceed this |
| H2 dimension | slope per doubling of n' > 1 means "dimension is not a lever" (abandon) | random / nb_window / nb_break: 1.07 to 1.36, CI lower > 1 in 14/16 (n, structure) groups (n=23 random/nb_window 0.86/0.87: grid capped at n'=8 of 11); standard 0.17 to 0.98 | **abandon rule fires for generic bases**; does not fire for standard |
| P2 detector | n=11 random n'=5 − n'=2, CI excluding 0 | **+1.92 [1.83, 2.00]** | detector sees the dimension effect |
| H3 δ | ≈ 1 degree per +2δ; ≥10× W/rel only at δ ≤ 2 | D\* drops 0.55 to 0.67 degree per +2δ (n = 11 to 19); W/rel gain at δ=+3 vs δ=−1: 6.5×, 11×, 14.6×, 7.9×, 19.7× | refuted as stated: weaker D\* drop, and large gains at δ=3. The gains rest on fractions of 1 to 2 of 16 to 24 targets, so this is not a reliable lever estimate |
| H4 sym4 | geometric = standard | identical D\* in all 5 cells (n=9,15,15,19,21); same nvars | confirmed; random higher (19 vars), subfield lower (9 vars) |
| H5 | r(D\*, nvars) high, r(D\*, n) ≈ 0 | r(D\*, n') = **0.82 [0.75, 0.88]**; r(D\*, n) = **−0.03 [−0.18, 0.12]**; r(D\*, (l+r)/n) = **0.88 [0.85, 0.91]** | confirmed; the invariant explains more than n' alone |

### Invariant null N2 (negative control) and Hypothesis I

N2 generates random Boolean systems matched on (n, d1, d2, l, q), with bilinear monomials only, 200 unsatisfiable systems per cell. They go through the same `refutation_degree`.

- **40 / 45 cells** have Semaev − null within ±0.25, which meets the frozen rule (≥ 80%). **Hypothesis I passes**: in this grid, structure affects D\* through (n, l, q) and not through the structure label.
- The 5 exceptions all go the same way, with Semaev easier:
  - standard n=11 n'=5: −0.48 [−0.65, −0.29];
  - standard n=13 n'=5: −0.62;
  - standard n=15 n'=5: −0.37;
  - random and nb_window at n=13 n'=5: −0.88 / −0.90 [−0.97, −0.83].
- Constraining the null's linear forms to `a_i + b_i` (N2s) leaves the same cells unexplained. This variant is exploratory and was not preregistered.
- Planted-q control P2b (n=13, n'=4, l=4): q = 4, 7, 10, 13 gives D\* = 2.03, 2.21, 2.76, 2.94. The null's own CIs are disjoint, so the detector sees r.
- **Conjecture, not tested:** the remaining effect is degree-3 structure, such as `dim(V·V·V)`, which is 3n'−2 for the standard window.

## Figures

- `fig_dstar_vs_dim.svg/.png`: the recorded D\* against n' for each structure and field, with 95% CIs.
- `fig_dstar_vs_invariant.svg/.png`: cell-mean D\* against `(l + r)/n`. This is the exact invariant; all structures collapse onto one curve, which rises at `l + r ≈ n`.
- `fig_null_vs_semaev.svg/.png`: Semaev against the N2 null, with a ±0.25 band.

## Cuts and deviations

- **Cut** (the port does not support them): chained m = 3 (the campaign never completed it either), symmetric-S4 replay, stacking, the n = 131 Frobenius-window check. H4 is judged from the recorded sym4 files only.
- **Cut from the replay:** n = 19 with n' = 8, 9, all of n = 23, and the δ-sweeps for n = 15, 17, 19. These use recorded values. Wall clock was not the binding limit; R5 was simply never preregistered.
- **Deviation 1.** The replay matcher first marked targets as non-matching when the raw file had no `macaulay_shape` (the n = 19 files). I changed the matcher to compare only recorded fields and counted shapes separately (table above). No D\* or rank value changed.
- **Deviation 2.** The first N2 run included n ≥ 17 cells, which went against the protocol ("N2 on every R1 cell"), and stalled on high-degree null systems. I killed it after about 9 minutes and reran it restricted to n ≤ 15 as preregistered. Nothing from the killed run was used.
- **Deviation 3.** N2s (symmetric linear forms) is exploratory and was added after seeing N2.
- Censoring: there are no censored targets in the m = 2 grid. In the sym4 file `sym4_n15_l3`, the random arm has 0 measured targets (D\* > 5 or not sampled).
- Scope: these are toy fields (n ≤ 23), D\* only on non-decomposable targets, and stage diagnostics. They support no claim about deployed curves, speed below rho, or end-to-end cost. Pair 46 is a prime-field pair and this binary study does not transfer to it. The repo's prime-field conclusion still stands.

## Files

`PROTOCOL.md`, `run.sh`, `dorv/` (Rust: replay, nulls, E2, analysis), `out/` (replay outputs, null, null_sym, e2), `results.json` (a copy of `out/results.json`), `plot.py` (layout only), the figures, `report.typ`, `report.pdf`.
