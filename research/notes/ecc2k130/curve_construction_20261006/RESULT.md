# Constructing a curve whose Jacobian carries the ECC2K-130 subgroup — result

Run 2026-10-06 against `PROTOCOL.md` (frozen before the run). Native Rust
(`crypto_lib::cryptanalysis::curve_construction`, binary
`ecc2k130_curve_construction`); report in `results/curve_construction.json`;
analysis in `research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md`.

## Against the falsification target

- **Success (curve of genus ≤ 300 carrying `A` with an evaluable
  correspondence): not met.**
- **Boundary refinement (construction below genus `2^129`): met.** The modular
  curve `X_H(3²·7²·263²)` mod 2 carries `A` at genus 508,799,809 (`2^28.92`).
  Modelled index-calculus cost ≈ `2^169,981` (extrapolation) — not an attack.

## Per route

| route | outcome |
|---|---|
| 1, lift to characteristic 0 | explicit curve at genus `2^28.92`; at `n = 3` the same construction is the Klein quartic (`X_H(49)`, genus 3) |
| 2, cyclic covers | abandoned by derived bound: genus ≥ 1300 for every geometrically cyclic degree-131 cover of a curve of genus ≤ 1; two-branch family excluded; Fermat quotients of exponent 1841 excluded |
| 3, direct search | toy: `A_n` is not a Jacobian at `n = 3, 5`; genus-3 curves carry `A` at `n = 3` for both signs; nothing at `n = 5` through genus 4 (complete) or genus-5 hyperelliptic. `n = 131`: `2^395` hyperelliptic models |
| 4, no curve | cited: bounded above rho or closed (scoreboard) |

## Checks

- Legacy replay of the Python-era boundary file: 20 values, max |diff|
  0.0048 (tolerance 0.006), every exact field equal (`r`, smoothness bounds,
  window crossover 290/300, `#A(F_2) = r`).
- GHS census: 400 random `b` (seed 20261006), magic number ∈ {1, 130, 131} with
  the trace correspondence on every sample (196 at 130, 204 at 131).
- Klein quartic: counts over `F_{2^k}`, `k ≤ 8`, equal `Res_{F_8/F_2}(E_1)`;
  over `Q`, 93 primes ≤ 499 match the Weil restriction of the conductor-49
  curve from `Q(ζ_7)^+`, 0 mismatches.
- Modular genera: `Γ_0` formulas equal enumeration at five levels; every cusp
  splits at all five.
- Kummer toys: 44 covers, all Weil polynomials, `E_0` splits off every one,
  every Prym in `Z[T^d]`.

## Correction carried in the same change

`RESEARCH_ECC2K130_HYPERELLIPTIC.md` claimed the genus any construction reaches
is 1, `2^129` or `2^130`. That is true of GHS/Hess descent only; corrected in
place (accounting).

## Changes after the first native run, disclosed

The protocol's route 2 named cyclic covers of `E_0` and Kummer toys at
`n = 5`. After the first native run showed genus-3 toy curves whose elliptic
factor is `E_1`, not `E_0`, the cyclic bound was extended to every base of
genus ≤ 1 over `F_2` (`P¹` and the five elliptic isogeny classes; minimum 1300
against 1301 for `E_0`), and the Kummer toy to `ℓ = 3` as well. Both are
extensions of the frozen route, not changes to it. The Picard-part condition
was tightened in the same edit, from "`n` divides the group order" to the
equivariant eigenvalue condition. That changes no `E_0` figure, since neither
version ever fires for `E_0`. The report was regenerated from the final code.
Contrary to the repository's rule that no run is overwritten, the first run's
file was overwritten during development and cannot be recovered. As far as can be
reconstructed, it differed only in the cyclic-bound section (`E_0` base only,
minimum 1301), the wording of two construction-table rows, and its runtime;
every other section was produced by unchanged code.
