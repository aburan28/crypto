# Item 3: the PKM tower's solving degree past the memory wall (protocol, frozen before data)

Campaign: `research/ic_next10_20261014/RULES.md`. Design: "### 3." of
`reports/Index calculus leads next steps.md`. Background:
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md` §§3, 10–16.

Written 2026-10-14, before any PKM system was handed to msolve in this item.

## 0. Check first (what already exists)

- `git log` and `research/` hold the PKM rounds 1–7 (`research/pkm_tower_*`),
  all measured with the in-house `f4_fp` / `f4_fp_tower`. The Kummer `m = 3`
  line at `p₁ = 2013265921` reads `D = 6, 6, 7, 7` at `N = 9, 12, 15, 18`; the
  `m = 4` line `D = 6, 7, 7` at `N = 8, 12, 16`; `m = 2` reads 5 at `N = 12–18`
  and 6 at `N = 20, 22`.
- `research/magma_f4_comparison_20261013/` exported every committed Kummer
  tower system as a Magma script (raw presentation), but no Magma result is
  committed. `research/f4_fp_msolve_parity_20261011/` compared `f4_fp` and
  msolve on other systems (timing only), not on PKM systems.
- No msolve run on a PKM system, no Hilbert-series computation and no
  Macaulay (XL) degree measurement exist in the repository. The note's §3.3
  bound `D_XL ≥ d + δ` is derived but `δ` was never computed.
- Caveat recorded in the note (§1.2, §4.4): the packaged msolve 0.6.5 was
  reported to fail a known-answer battery above `2^16` in another project. So
  msolve's verdict and degree are admitted only where the positive control
  (below) passes at the run's prime.

## 1. Hypotheses

- **H-plateau.** On the Kummer tower family at `p₁`, msolve's F4 solving
  degree `D` stays at 7 for `m = 3` at `N ≥ 21` (and at `m = 4` at `N = 20`).
- **H0 (growth resumes; the design's kill criterion).** `D ≥ 8` at `N = 21`
  and `D ≥ 9` at `N = 24` for `m = 3`.
- **H-XL (the conflict with `d + δ`).** The Macaulay/XL solving degree
  `D_XL = d + δ` grows like `N`, while F4's `D` does not: the §3.3 bound binds
  XL only, and F4 escapes it by degree falls. This item computes `δ` exactly
  where it can.

## 2. Arms

| arm | what | role |
|:--|:--|:--|
| PC | msolve 0.6.5 on committed Kummer tower systems of rounds 2–7 with known `f4_fp_tower` `D` and verdict: `m = 2`, `N = 12–22`; `m = 3`, `N = 9, 12, 15` (and 18 if it fits); `m = 4`, `N = 8, 12` (and 16 if it fits). Both targets each. | positive control: msolve's max round degree and verdict must equal the committed `D` and `inconsistent` |
| T | msolve on new systems drawn by the same generator (`examples/pkm_tower_pilot.rs`, same seed and flags, targets 0 and 1): `m = 3` at `N = 21, 24, 27` (`t = 7, 8, 9`), `m = 4` at `N = 20` (`t = 5`) | the measurement |
| PL | planted-solution systems at small `N` (pilot's `--planted`) | positive control for msolve's correctness at `p₁` (a planted point must be among msolve's solutions) |
| XL | exact Macaulay/XL degree via `d + δ`: for `m = 2` (`full`) `δ = deg(1/S̄)` by a 2-D number-theoretic transform on the Kummer grid; for `m = 3` (`chain`) a direct Macaulay rank computation in the tower quotient at small `N` | the comparison arm for the lower bound |
| HS | Hilbert series: the top-degree (homogenised-at-infinity) ideal's series and its degree of regularity, and the affine Hilbert function of the homogenised generators from the XL rank computation | does a Hilbert-series predictor reproduce `D` at `N ≤ 18`? |

There is no comparative performance claim, so no negative-control
population is needed; the shape-matched null was already run in rounds 1–2.

## 3. Sizes and budgets

- Prime `p₁ = 2013265921` throughout (the prime of the committed plateau; the
  design text says `p ≈ 2^24`, but the plateau was measured at `p₁ ≈ 2^31`, so
  `p₁` is used for comparability — a stated deviation from the design text).
- Two targets per cell (the generator's targets 0 and 1), as in rounds 2–7.
- Each msolve job: `nice`, `-t 1`, `-v 2`, under `prlimit --as=5000000000`
  (5 GB, this item's cap) and `timeout 2700` (45 min).
- Order: controls smallest first, then `N = 21`, then larger only if `N = 21`
  completes.

## 4. Metric

- `D_msolve` = the maximum `deg` over msolve's F4 rounds (`-v 2` table) on a
  run that ends with a verdict. A run killed by the memory or time cap gives
  `D ≥ d*`, where `d*` is the largest degree of a completed round that added
  new basis elements (the note's §16.1 convention).
- Peak memory: `/usr/bin/time -v` maximum resident set size; wall time.
- Verdict: `[-1]` (no solution) vs a parametrisation.
- `δ`, `D_XL = d + δ`, `d_reg` of the top-part ideal.

## 5. Decision rules (fixed now)

- **Instrument first.** If msolve's verdict disagrees with a committed verdict
  or with exhaustive search, or a planted solution is missing, msolve is not
  admitted at `p₁` and the item is BLOCKED on the instrument. If verdicts agree
  but `D_msolve ≠ D_f4` on the controls, msolve's degree is reported as its own
  line (engine-dependent `D`), and plateau conclusions at new `N` are drawn
  only from msolve's own line.
- **Kill criterion (from the report).** `D ≥ 8` at `N = 21` and `≥ 9` at
  `N = 24` (m = 3) closes the PKM route (NEGATIVE for the plateau). `D = 7` at
  `N = 24` supports the plateau. `D = 7` at `N = 21` alone is one more level:
  INCONCLUSIVE under the note's A1 rule unless `N = 24` also completes.
- A run that hits the 5 GB / 45 min cap is reported as a wall, with its lower
  bound, never extrapolated to a value of `D`.
- **H-XL** is supported if `D_XL − D_F4` grows with `N` across the measured
  sizes; refuted if `D_XL` stays within one of `D_F4`.

## 6. Predictions (before data; low confidence)

1. msolve reproduces the committed `D` on every control that fits (both are
   normal-strategy F4 on the same ideal, though msolve takes the raw and
   `f4_fp_tower` the reduced presentation).
2. `m = 3`, `N = 15` fits in 5 GB; `N = 18` does not (f4_fp_tower needed
   13 GB); `N = 21` does not. So the memory wall at 5 GB falls between
   `N = 15` and `N = 18` for `m = 3`, and between 12 and 16 for `m = 4`.
3. `δ` for `m = 2` grows like `N − O(1)`, so `D_XL ≈ N + 4` while `D_F4 ≤ 6`.
4. The top-part ideal's Hilbert series gives `d_reg = N + 1` (full) or no
   finite `d_reg` (chain), so it does not predict `D`.
