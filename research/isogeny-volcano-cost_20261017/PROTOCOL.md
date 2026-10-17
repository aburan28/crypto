# PROTOCOL v1: vertical, horizontal and coordinate effects on ECDLP cost across isogeny volcanoes

**Frozen 2026-10-17, before any measured run.** Everything after this
commit that changes N, the classes, the arms, the metrics, the margin, the
decision rules or the correction is a new, labelled protocol version, and
both versions are reported. This file and the exported specs are pinned by
the commit that adds them. The digest of
`sha256sum classes/*/spec*.json classes/*/curves.jsonl` is
`f8f169d62dda91e110a8950b0e721027bbaef0539b240469e8d010e113b7047c`.
The design and its reasoning are in [`PLAN.md`](PLAN.md).

## 1. Question and licensed claims

Within one ℓ-isogeny volcano over `F_p`, every curve has the same `#E`, and
one ECDLP instance is carried to all of them. Does solving cost differ:

- **vertically**: between levels (crater against floor)?
- **horizontally**: between curves of the same level? Such curves have one
  endomorphism ring, so they are horizontally isogenous through the class
  group.
- **by coordinates**: between isomorphic models of one curve? This is the
  null an isogeny effect has to beat.

The methods are Pollard rho (`rho.negation`), BSGS (`bsgs.negation`) and
index calculus (`ic.pipeline`, `prime-abscissa:size=32`,
`mitm:negation_folded=1`).

Wording is fixed by §6. "No difference" is never written; only an
equivalence bound or a measured effect is.

## 2. Population (exported by `volcano-rho --ecbench`, 24-bit `p`, depth 3)

| class | p | root a, b | curves per level (crater → floor) | curves | `E(F_p)` structure (count) | crater |
|:--|--:|:--|:--|--:|:--|:--|
| `l2-ramified` | 16517689 | 1675631, 136892 | 2 / 4 / 8 / 12 | 26 | Z/16522852: 12, Z/8261426 × Z/2: 14 | ramified |
| `l2-split` | 13151689 | 668241, 11547039 | 11 / 9 / 12 / 12 | 44 | Z/13146004: 12, Z/6573002 × Z/2: 32 | split |
| `l2-inert` | 9298073 | 1597478, 3640773 | 1 / 3 / 6 / 12 | 22 | Z/4647098 × Z/2: 10, Z/9294196: 12 | inert |
| `l3-ramified` | 10315087 | 4986021, 8775343 | 2 / 6 / 12 / 12 | 32 | cyclic: 32 | ramified |
| `l3-split` | 16759147 | 69845, 9777429 | 11 / 12 / 12 / 12 | 47 | cyclic: 47 | split |
| `l3-inert` | 12450847 | 11354569, 10721363 | 1 / 4 / 12 / 12 | 29 | cyclic: 29 | inert |
| `j0` (positive control) | 13151689 | 0, 9554511 | 1 / 1 / 2 / 4 | 8 | Z/6577962 × Z/2: 4, Z/13155924: 4 | j = 0 |
| `j1728` (positive control) | 9420017 | 7905883, 0 | 1 / 1 / 2 / 4 | 8 | Z/4710578 × Z/2: 4, Z/9421156: 4 | j = 1728 |

- **Commands.** Each class is
  `volcano-rho --bits 24 --depth 3 --max-cofactor 4096 <class args> --controls 0 --per-level 12 --ecbench classes/<class> --models 4 --targets 16`.
  The class args and seeds are in `classes/<class>.log`.
- **Selection.** Up to 12 curves per level in BFS order, with the root and
  every crater curve always kept. Every walk edge is a verified Vélu
  isogeny.
- **Models.** Each curve gets 4 isomorphic models (model 0 as walked,
  models 1–3 with random u), so 216 curves give 864 workload curves.
- **Pooled analysis.** The six ordinary classes form the pooled analysis.
  `j0` and `j1728` are positive controls, analysed apart.
- **Structure and level.** Group structure is reported, not adjusted for:
  for ℓ = 2 it is determined by level (floor curves have one rational
  2-torsion point).

## 3. Sessions (all at L0, counted units only)

| spec | arms | workloads | role |
|:--|:--|:--|:--|
| `spec.json` | `rho-neg`, `ic-mitm` | all models, 16 `planted_shared` targets | primary |
| `spec-aa.json` | `rho-neg-aa`, another measurement seed | same | independent A/A |
| `spec-bsgs.json` | `bsgs-neg` | all models, 2 targets | exactness |
| `spec-icctl.json` | `ic-mitm-plus1` (`size=33`) | model 0, 16 targets | planted IC control |
| `spec-crater-aut.json` (`j0`, `j1728`) | `rho-neg-crater`, `rho-cubic` / `rho-quartic` | crater, 4 models, 16 targets | automorphism control |

Total: 46,912 runs. Expected time is about 70 min serial, from the dry
runs' ≈ 64 ms per-run process overhead plus solve time.

- **Run:** each session with
  `ecbench run --wait --allow-busy --cpus inherit` into `sessions/<class>/<spec>/`.
- **Audit:** `ecbench verify --replay 12` on every session; the receipts
  are committed.
- **Storage:** sessions are committed with `records.jsonl` gzipped (about
  280 MB raw, too large for git) and a SHA-256 manifest.
  `ecbench verify` runs on the decompressed files.

## 4. Metrics

Each run is a ratio `num/den`, analysed by `volcano-stats` (native).

| family | metric | per run |
|:--|:--|:--|
| F1 | rho cost | `S` (GAE/√r) |
| F2 | BSGS exactness | number of distinct `S` per paired target across all curves and models |
| F3 | IC relation yield | `relations_found / targets_tried` |
| F3s | IC cost | `S` (reported with the §5a band through `ecbench band`, against `rho-neg`; the gap to rho is distance remaining) |

## 5. Statistics and decision rules

Tested families: F1 and F3 × {vertical, horizontal}, so m = 4. α = 0.05,
Bonferroni α/m = 0.0125. Equivalence margin M = ±5% on a 90% interval.
F3s, the model tests and F4 are reported but not tested.

- **Vertical:** V = floor/crater ratio of class-normalised pooled means.
  - Interval: a 90% (and 95%) two-stage bootstrap with 2,000 replicates.
    Each replicate draws one set of paired target indices per class,
    shared by every curve, then resamples curves within class × level.
  - p: 10,000 permutations of level labels among curves within class.
  - **Equivalent within ±5%** if the 90% interval lies in [0.95, 1.05].
    **Vertical effect** if p < 0.0125 and the effect is outside the
    interval of the A/A arm (`rho-neg-aa`, same analysis).
- **Horizontal:** cells are curves, strata are (class, level).
  - Statistics: Cochran Q and the largest standardised deviation max|z|.
  - Null: 2,000 permutations reassigning runs among the curves of a
    stratum.
  - **Horizontal effect** if either permutation p is below 0.0125/2,
    splitting the family's α over its two statistics.
  - Effect size: the DerSimonian–Laird τ (relative), reported beside the
    null's 95th percentile.
- **Model (coordinate) null:** the same test with cells = models and
  strata = curves. A horizontal effect counts as an isogeny effect only if
  the model test is not significant, or the horizontal τ exceeds the model
  τ beyond its null band. Otherwise it is reported as a coordinate effect.
- **F2 BSGS:** prediction is exactly 1 distinct cost per paired target over
  all 864 workload curves. Any other count falsifies the transport or the
  instrument, with no statistics.
- **Positive controls,** each of which must pass:
  - `ic-mitm-plus1 : ic-mitm` yield, matched by workload: predicted
    ≈ (33/32)² = 1.063, so the 95% CI must exclude 1 and contain the
    prediction's direction.
  - `rho-cubic : rho-neg-crater` on `j0`: predicted √(1/3) ≈ 0.58 in S, so
    the 95% CI must exclude 1.
  - `rho-quartic : rho-neg-crater` on `j1728`: predicted ≈ 0.71. It is
    reported, but not required to pass, because the earlier sessions found
    `rho.quartic` below its √2 fold.
  - The analysis's own controls: `volcano-stats --plant-level 3:1.05` and
    `--plant-curve <a floor vertex>:1.10` on the real records must be
    detected, which shows the power at the margins.

  If the first two or the analysis controls fail, the run is
  **inconclusive**.
- **Power:** the smallest permutation p is 1/10,001 against α/m = 0.0125,
  and N + 1 ≥ 4m/α = 320, so the check passes.
  - Expected interval half-widths from the pre-flight variances:
    - F1 vertical: about ±2% (crater ≈ 28 curves × 64 runs, rho sd/mean
      0.46).
    - F3 vertical: about ±1% (≥ 35,000 crater relations).
    - Per curve: rho ±5.6% and IC yield ±3.5% (64 runs each).

## 6. Fixed verdict wording

- "level changes *method* cost by less than ±5% (90% CI [lo, hi]) on these
  six classes at 24 bits"
- "level changes *method* cost by d% [95% CI], permutation p = …"
- "curves of one level differ by τ = …% (null 95th percentile …%),
  p = …: an isogeny / coordinate effect"
- "positive control failed: inconclusive"

Scope: 24-bit `p`, depth-3 volcanoes, ℓ ∈ {2, 3}, these three methods, at
these parameters. This is not a security statement, and it says nothing
about larger fields or other factor bases.

## 7. Dry-run disclosures (scratch output, not data)

- **Pre-flight** (`preflight/`): the instance was carried to all 32 curves
  of one 24-bit volcano. BSGS was identical per target; rho sd/mean was
  0.46; IC found 20 relations per run.
- **Dry run on `l2-ramified`'s volcano, 2 models × 6 targets.**
  - BSGS was exact.
  - The planted IC control read 1.068 [1.030, 1.108].
  - A planted ×1.5 level effect was detected (p = 0.0005).
  - A planted ×1.5 single-curve yield effect was detected by max|z|
    (p = 0.005) but not by Q, which is why max|z| is in the rule.
  - Unplanted, rho's vertical ratio read 1.26 with p = 0.019–0.035
    (1–2 crater curves × 12 runs). It was not significant at 0.0125, and
    that level of noise is what this protocol's N addresses.
- **Dry run on `j0`, 2 targets.** All five spec types ran, and the five
  sessions ran in parallel without conflict.
- **Tooling bugs found and fixed before freezing:**
  - The bootstrap resampled paired targets independently per curve, which
    gave BSGS a spurious ±5%.
  - Classes were keyed by p, but `j0` and `l2-split` share p.
  - The crater control's `rho-neg` shared a name and seeds with
    `spec.json`'s, so the two would have been pooled.
  - The planted pair was pooled across classes instead of matched by
    workload.

## 8. Does not count

- Changing N, the margin, the metrics, classes, arms or the correction
  after seeing results.
- Dropping a curve, a model or a failed run (`incomplete` stays
  `incomplete`).
- Re-running a session because of its result.
- Selecting a fast curve without fresh-target confirmation.
- Citing wall time from this L0 host.
