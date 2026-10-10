# IC lead ladder, 2026-10-10 (Conductor T-248)

Working directory for `docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md`.

| file | what |
|---|---|
| `sweeps/baseline_l1_l4.json` | the baseline configuration and the L1 (signed Frobenius folding) and L4 (invariant base) toggles as a matrix (paired by seed) |
| `run_baseline.sh` | runs `ic bench` on `K_0` per rung, operation counts, one report per rung under `runs/<stamp>/` with command, commit and binary hash |
| `summarise_runs.py` | turns the reports into the per-phase ratio table of the methodology §3 |
| `e7_vertical_step_costs.py`, `e7_vertical_step_costs.json` | L8 / E7: order-of-magnitude cost of one vertical isogeny step per `(n, ell)` in `F_q`-multiplications, three routes, against rho |
| `RESULTS.md`, `RESULTS_TABLES.md` | the first baseline run (n = 13, 17, 23): L1 and L4 measured, readings, what was not measured |
| `runs/` | raw reports; wall time in them is not evidence (host load 90–150) |
| `e2_floor_curves_polclass.sage` | L8 / E2: floor curves one level down from K_0 by the CM route (ring class polynomial mod 2) |

## E7 result (2026-10-10, arithmetic model, no Sage)

Model per route, in `F_q`-multiplications with `M(d) = d^1.585` for one `F_{q^d}`
multiplication: kernel route = point of order `ell` over `F_{q^d}` + Vélu
(`(ell−1)/2` kernel multiples) or √élu (`√ell·log ell`) + descent of `j`;
modular route = evaluate `Φ_ell(X, j)` (`ell²` coefficients) + root-find a
degree-`(ell+1)` polynomial over `F_q`; transfer = the isogeny on two points.
These are models, not measurements; the constants (4, 12, 6 field
multiplications per step) are generous round numbers.

| n | ell | χ(−7) | d | cheapest descent | route | rho (A = 2n) | descent < rho |
|--:|--:|--:|--:|--:|:--|--:|:--|
| 23 | 967 | +1 | 21 | 2^19.4 | kernel | 2^8.1 | no (toy: rho is trivial) |
| 37 | 73 | −1 | 18 | 2^13.4 | modular | 2^14.7 | yes |
| 37 | 2663 | −1 | 2662 | 2^22.8 | modular | 2^14.7 | no |
| 41 | 409 | −1 | 408 | 2^17.6 | modular | 2^16.7 | no, within 2× |
| 41 | 1721 | −1 | 215 | 2^21.6 | modular | 2^16.7 | no |
| 53 | 68476319 | +1 | 646003 | 2^52.1 | modular | 2^22.5 | no |
| 61 | 1951 | −1 | 975 | 2^22.0 | modular | 2^26.4 | yes |
| 61 | 587551 | −1 | 117510 | 2^38.3 | modular | 2^26.4 | no |
| 71 | 3.2·10¹⁰ | −1 | 3.2·10¹⁰ | 2^69.8 | modular | 2^31.3 | no |
| 73 | 2.3·10¹⁰ | −1 | 1.2·10¹⁰ | 2^68.9 | modular | 2^32.2 | no |
| 83 | 6473 | −1 | 6472 | 2^25.4 | modular | 2^37.1 | yes |
| 83 | 53676929 | −1 | 1.3·10⁷ | 2^51.4 | modular | 2^37.1 | no |
| 97 | 751943 | −1 | 375971 | 2^39.0 | modular | 2^44.0 | yes |
| 97 | 352124743 | −1 | 5.9·10⁷ | 2^56.8 | modular | 2^44.0 | no |
| **131** | **263** | +1 | **2** | **2^13.0** | kernel | 2^60.8 | **yes** |
| **131** | 1.47·10¹⁷ | −1 | 1.2·10¹⁶ | 2^114 | modular | 2^60.8 | **no** |

Readings:

1. The two components of the ECC2K-130 class are separated by at least
   `2^114 / 2^60.8 ≈ 2^53` under every route modelled; the 263-level costs
   `2^13` to reach (the parallel workspace measured `2^15.6`, within the
   model's slack). This is the operational content of "conductor gap".
2. On the toy ladder the kernel route is almost never the cheapest; the
   modular route wins whenever `ell ≤ 10⁴`, and above `ell ≈ 10⁶` nothing is
   cheaper than solving the DLP by rho. The constructible floor levels for
   hardness-difference experiments (E2) are therefore `n = 61 / 1951`,
   `n = 83 / 6473`, `n = 97 / 751943` and (marginally) `n = 41 / 409`; at
   `n = 23` everything is trivial.
3. "descent < rho" is the only honest reachability criterion: a level whose
   cheapest descent exceeds rho's cost is unreachable *for the purpose of
   changing the DLP's hardness*, whatever its mathematical status.
