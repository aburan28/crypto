# 02-deuring-lift: results

Item 02 of the pair-46 factor-base campaign. Protocol: `PROTOCOL.md`, frozen
before any data was generated. Raw data: `results.json`. Statistics:
`summary.json`. Log: `run.log`. Figures: `fig_relations.{svg,png}` and
`fig_twist_ranks.{svg,png}`. Report: `report.pdf`.

**Verdict: exact obstruction. The prediction (failure by rank obstruction) holds.**
The canonical lift gives no characteristic-0 relation that involves the
target. The lifted generators of the canonical-lift family behave exactly like
random points mod p. A generic Xedni lift also gives no relations.

## Exact statements (derivation plus exact PARI computation)

1. **The canonical lift is 361a1 itself.** E: y^2 + y = x^3 - 38x + 90 has
   j = -884736 = -96^3 and conductor 361. For pair 46, a_p(E) = 34489701 = t,
   so #E(F_p) = 5947697996161643 = N. The isomorphism to the crater is given
   explicitly by u^2 = (B_c/B)/(A_c/A), which is a square mod p, and both
   checks u^4 A = A_c and u^6 B = B_c pass. No quadratic twist is needed.
2. **E(Q) = Z·G.** The rank is 1 unconditionally (`ellrank` bounds [1,1]),
   and the analytic rank is also 1. The generator is G = (0, 9) on the
   minimal model, the torsion subgroup is trivial, and ĥ(G) = 0.181292 in
   PARI's normalization.
3. **Image and circularity.** G mod p has order N, so the image of E(Q) is all
   of E(F_p). Every target T therefore has a canonical lift, namely kG with
   k = log_G T, and that lift has height k^2 ĥ(G). For uniform k the expected
   height is ĥ N^2/3 ≈ 2.1·10^30 nats (log2 ≈ 100.8).
   - Writing down the lift of T requires k, so lifting is the DLP itself.
   - Every lifted factor-base point is a known multiple mG, so its log is free
     and tells us nothing about T.
   - The number of lifted relations that involve the target is 0 unless
     log T is already known.
4. **Twists (one per {d, -19d} isogeny class, 1 < |d| ≤ 200, 115 twists, all
   `ellrank` bounds closed).**
   - Rank census: 51 twists have rank 0, 54 rank 1, 8 rank 2 and 2 rank 3. The
     maximum is 3.
   - The twists with χ_d(p) = +1 for pair 46 are listed in
     `results.json → exact_pair46.twists_chi_plus_at_pair46`. Examples: d = -67
     has rank 2, d = 69 rank 3, and 8, 21, 37, ... have rank 1.
   - The points of E_d(Q) for different d lie in distinct Gal(L/Q)-eigenspaces
     of E(L) ⊗ Q, where L is the multiquadratic field. They are therefore
     Z-independent. A factor base built from twist generators has rank
     equal to its size and **zero** Q-relations.
   - Every relation among these points exists only mod p. Such relations are
     points of the kernel lattice of Z^R → Z/N, which has determinant N, and
     finding them is the R-dimensional DLP.
5. **Over K = Q(√-19)** (derivation, not computed here). E(K) has Z[w]-rank 1
   and Z-rank 2. The only extra structure is wG ≡ λG mod p. That is the
   known GLV/endomorphism relation, which the glv-orbit item in the inventory
   already covers. It is not a new relation source.

## Empirical results (24 toy primes: 8 each at 16, 20 and 24 bits, each with prime #E(F_p); factor-base size R = 6)

| metric | A2 canonical twists | NC random points | PC planted | effect A2 − NC (95% bootstrap CI) |
|---|---|---|---|---|
| λ1/GH, median | 1.393 | 1.250 | 0.290 (max 0.461) | +0.143 [−0.180, +0.262]; Mann–Whitney p = 0.50 |
| coverage, \|c_i\| ≤ 2 (mean diff) | — | — | — | +0.0008 [−0.052, +0.056] |
| lifted Q-relations, total | 0 / 24 primes | n/a | 48 (= 2 × 24, planted) | — |
| Xedni generic lift, r = 2..5 | 0 of 96 trials (all lifts p-integral and on the curve) | n/a | 24/24 (P, −P planted) | — |

- **Positive controls.** Both detectors see planted structure. The mod-p
  λ1/GH was below 0.5 at every prime, and every planted Q-relation was found.
- **Same instrument on every arm.** Every arm was measured with the same
  instrument (PARI `elllog`/`qflll`/`qfminim`, `ellheightmatrix`).
- **Effect size.** The A2 arm cannot be distinguished from random points.
  The confidence interval contains 0, which meets the preregistered
  "indistinguishable" criterion.

## Deviations and cuts

- **Nothing in the preregistered design was cut.** The whole run took 2.4 s
  of CPU on 1 core.
- **Instrument precision (fixed after the smoke run, protocol unchanged).**
  cypari2's default precision is 15 digits, and at that precision the
  planted G, 2G, 3G Gram matrix showed eigenvalues around 1e-17. That was
  above the frozen 1e-20 tolerance, so the positive control failed. I raised
  `realprecision` to 80 for every arm. The tolerance itself was not changed.
- **Scope.**
  - One generator per twist (the protocol says one each).
  - Twists only up to |d| ≤ 200.
  - No DLP was attempted at pair-46 size. The pair-46 statements are exact
    and need no DLP.
- **Repository note.** AGENTS.md prefers non-Python research code. All the
  arithmetic here is done by PARI/GP (through cypari2); Python only
  orchestrates the calls and computes the bootstrap statistics.

## Reproduce

```sh
pip install cypari2 matplotlib typst
cd research/pair46_factor_bases_20261015/02-deuring-lift
python3 run_deuring_lift.py results.json && python3 analyze.py \
  && python3 -c "import typst; typst.compile('report.typ', output='report.pdf')"
```
