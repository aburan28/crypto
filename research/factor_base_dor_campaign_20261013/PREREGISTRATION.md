# Factor base → degree of regularity: yield and correlation campaign

**Date:** 2026-10-13
**Producer:** `sage/factor_base_dor.sage` (selftest required before trusting a number)
**Evidence scope:** synthetic toy binary Koblitz fields only; no deployed-curve security claim, no sub-ρ claim. Falsification target is a correlation, not a speedup claim.

## Question

Does the factor base (dimension `n'`, algebraic structure, and Shantz–Teske slack `δ`) change the degree of regularity `D*` and the decomposable yield of the Weil-descended Semaev decomposition system enough to constitute a lever toward a drastic index-calculus speedup on the ECC2K-130 family (`y² + xy = x³ + 1`, `x^n − 1 = (x+1)·irreducible`, no proper subfields)?

## Boundaries

- **Reference:** decomposition search by brute enumeration over the base, cost `|V|^m` per target — the engine already reports decomposable fraction and `trials = (unknowns+1)/coverage` where applicable.
- **Floor:** yield per target is bounded by counting: for base `V` (|V| = 2^n') and m summands the expected decomposable fraction is `~|V|^m / 2^n` (up to orbit constants). `D*` is bounded below by the refutation degree of the descendent system.

## Hypotheses (falsification targets)

- H1 (structure): bases closed under multiplication (subfield) refute at lower `D*` than standard/random bases at the same `n'` (prior: `D* = 3` vs ≥4). Frobenius/normal-basis windows are NOT closed; prediction: `D*` equal to random.
- H2 (dimension): `D*` is flat or slowly rising in `n'` (prior result at n=8–9 suggests mean `D*` 2–3 at small dims); the yield `≈ |V|^m/2^n` rises exponentially. If `D*` is flat, the work-per-relation proxy `W(n') = rows(D*)×cols / fraction` has an interior optimum — the drastic-speedup candidate.
- H3 (δ, Shantz–Teske): raising `δ` (smaller base, `m·n' = n − δ`) lowers `D*` ≈ one degree per +2δ but multiplies yield by `2^{-δ}`; net work-per-relation improves 10× only at `δ ≤ 2`.
- H4 (m=3, symmetric S4): the geometric special subspace does not lower `D*` vs the chain; only multiplicatively closed bases do.
- H5 (correlation): `D*` correlates with the number of variables (Macaulay cols) and with the multiplicative-closure of the base, not with the field degree n per se; quantified by Pearson r over the grid.

Success: H1–H5 each confirmed or refuted with numbers, ≥3 seeds where random structures are involved, selftest PASS before any cited number. Abandon: if `D*` rises faster than one degree per doubling of `n'`, dimension is not a lever and the campaign reports the negative.

## Frozen inputs

- Fields: n = 9, 11, 13, 15, 21 (prime-shape controls 11, 13 mirror 131; 9, 15, 21 carry subfields as closure contrasts; all odd).
- m = 2 (chained S3) and m = 3 (chained and `--symmetric-s4`).
- Structures: standard, random, subfield, frobenius, nb_window, nb_break, geometric (where defined).
- Dims: n' = 2..(n−1)/m; δ-sweep −1..+4 on standard.
- Targets: 48 per cell (single-target D* sampling `--dstar-sample 6` where chained m=3); seeds 1,2,3.
- `dmax = 6` (censor above and report).
- n=131 applicability checks: `--frobenius-windows --field 131`, weight-class `--stacking` at n = 11, 13.

## Accounting

Per cell: `d_star_mean/min/max`, `censored_above_dmax`, `decomposable_fraction`, `unsat_targets`, Macaulay shape at `D*`, work proxy `rows×cols/fraction`. Analysis computes correlations (D* vs n', vs nvars, vs closure flag; fraction vs |V|^m/2^n) and the work-per-relation optimum per field/summands. All numbers are stage diagnostics of the decomposition oracle; any end-to-end speedup remains an extrapolation until a matched m=83 pipeline run prices it (repo AGENTS.md §8a).

## Commands

See `run_campaign.sh` (thin shell orchestration of the sage engine only; no Python analysis beyond the engine).
