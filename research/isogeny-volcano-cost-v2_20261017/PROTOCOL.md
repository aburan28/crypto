# PROTOCOL v2: more crater curves, a confirmatory test of the ℓ = 2 IC-yield lead

**Frozen 2026-10-17, before any measured run.** The digest of
`sha256sum classes/*/spec*.json classes/*/curves.jsonl` is
`0f7c63764c796c02256cbf6072519fd0ef48dfdda193d1aa9b19f9d906b0e12b`.
Any later change is a new, labelled version, and both are reported.

The predecessor is [`../isogeny-volcano-cost_20261017`](../isogeny-volcano-cost_20261017/README.md)
(PROTOCOL v1 and v1.1). Its limits were:

- the rho vertical bound was ±6.2%;
- the IC vertical test was inconclusive;
- rho single-curve outliers below about ±15% were not excluded;
- an exploratory lead: ℓ = 2 floors gave about 6.5% fewer IC relations
  than craters, while ℓ = 3 showed nothing.

v1 found the fix for all of these to be more crater curves, through more
classes and more models per crater curve. v2 does that on fresh classes.

## 1. Hypotheses and licensed claims

Holm correction over H1–H3 at family α = 0.05.

| id | hypothesis | test | licensed wording |
|:--|:--|:--|:--|
| **H1** (primary) | in ℓ = 2 volcanoes, IC relation yield is lower on the floor than on the crater | `volcano-stats --vertical stratified --alternative less` on the four ℓ = 2 classes, IC yield, permutation p (10,000) | if Holm-significant: "ℓ = 2 floors give d% [95% CI] fewer IC relations per trial than craters (4 fresh classes)". Otherwise: "the v1 ℓ = 2 lead is not confirmed" |
| **H2** | in ℓ = 3 volcanoes, IC yield is equal across levels | the same analysis on the four ℓ = 3 classes, two-sided; equivalence if the 90% CI ⊂ [0.95, 1.05] | "within ±5%" if equivalent; a significant result is reported as an effect |
| **H3** | rho cost is equal across levels | all 8 classes, rho S, two-sided; equivalence if the 90% CI ⊂ [0.95, 1.05] | "within ±5%" |
| H4 (replication, α = 0.05) | the v1 coordinate effect: IC yield varies more between isomorphic models of one curve than between same-level curves | model Q permutation p < 0.05 **and** model τ > horizontal τ, all 8 classes | "replicated" / "not replicated" |
| H5 (reported, untested) | IC cost S floor/crater per ℓ group, and rho against IC across the §5a band | `ecbench band` per class | descriptive |

**Mechanism contrast, inside H2.** `l3-s1003` has `E(F_p) ≅ Z/n × Z/2` at
every level, so its 2-power torsion does not change between levels. If
the ℓ = 2 effect comes from 2-torsion structure, this class should show
none. Its per-class contrast is reported but not tested separately.

## 2. Population (24-bit `p`, depth 3, split craters)

**Selection rule (applied before freezing):** for ℓ ∈ {2, 3}, run seeds
1001, 1002, … with `volcano-rho --bits 24 --depth 3 --max-cofactor 4096
--ell ℓ --crater split`. Accept a seed if its walk reaches at least 8
crater curves and its p is not one of v1's. Take the first four accepted
seeds per ℓ. Seeds 1001–1004 were accepted for both ℓ, and none were
rejected.

| class | p | root a, b | curves per level (crater → floor) | workload curves | `E(F_p)` (count of curves) |
|:--|--:|:--|:--|--:|:--|
| `l2-s1001` | 15808601 | 2021557, 15763640 | 15 / 4 / 4 / 12 | 320 | Z/7905498 × Z/2: 23, Z/15810996: 12 |
| `l2-s1002` | 13897337 | 8685960, 9251477 | 11 / 4 / 4 / 12 | 256 | Z/13904116: 12, Z/6952058 × Z/2: 19 |
| `l2-s1003` | 11073121 | 438735, 8479511 | 13 / 4 / 4 / 12 | 288 | Z/5539082 × Z/2: 21, Z/11078164: 12 |
| `l2-s1004` | 10962617 | 7780919, 5486135 | 15 / 4 / 4 / 12 | 320 | Z/5483642 × Z/2: 23, Z/10967284: 12 |
| `l3-s1001` | 10102303 | 5274709, 3331729 | 15 / 4 / 4 / 12 | 320 | Z/10101316: 35 |
| `l3-s1002` | 10308967 | 3471117, 9312808 | 13 / 4 / 4 / 12 | 288 | Z/10314484: 33 |
| `l3-s1003` | 11845507 | 5796691, 8048851 | 11 / 4 / 4 / 12 | 256 | Z/5922716 × Z/2: 31 |
| `l3-s1004` | 9887107 | 8527579, 8217853 | 9 / 4 / 4 / 12 | 224 | Z/9881434: 29 |

**Export command:** `volcano-rho … --level-caps 16,4,4,12 --models 4
--crater-models 16 --aa-model0 --targets 8`.

- Every crater curve is kept, in **16 isomorphic models**.
- The two middle levels keep 4 curves each and the floor keeps 12, all
  in 4 models.
- Every model carries the isogeny-transported generator and 8
  `planted_shared` targets.

## 3. Sessions

The arms are those of v1, with these differences:

- the A/A arm covers model 0 only;
- BSGS uses 2 targets;
- the planted IC control (`size = 33`) covers model 0 only.

**Size:** 45,088 runs. They are run in class groups, each under the
host's 2-hour background limit. An interrupted session is kept and re-run
into a new directory, and the analysis uses the complete one.

**Audit:** `ecbench verify --replay 12` on every session.

## 4. Analysis

- **Tool:** `volcano-stats`, the class-stratified vertical statistic
  (v1.1), bootstrap 2,000 and permutations 10,000. The ℓ = 2 and ℓ = 3
  groups are analysed separately for H1 and H2, and all 8 classes
  together for H3 and H4.
- **Positive controls:** if any of these fails, the dependent hypothesis
  is **inconclusive**.
  - **Planted IC:** `ic-mitm-plus1 : ic-mitm` yield, matched by workload;
    the 95% CI must exclude 1. Required for H1 and H2.
  - **BSGS exactness:** 1 distinct cost per paired target. Required for
    every hypothesis, since it checks transport and pairing.
  - **H1's power:** `--plant-level 3:0.95` on the ℓ = 3 group's IC yield
    must give a one-sided `less` p < 0.0167, Holm's smallest threshold.
  - **H3's power:** `--plant-level 3:1.05` on all 8 classes' rho must
    give p < 0.0167.
  - **Single-curve power, reported:** `--plant-curve <a crater vertex>:1.15`
    must be detected by max|z| for rho and IC yield. It is reported, not
    required.
- **Expected precision,** from v1's variance components, so it is a
  prediction:
  - IC-yield floor/crater: about ±1.8% per class, about ±0.9% for the
    four ℓ = 2 classes together.
  - Rho: about ±2.6% per class, about ±0.9% across all eight.
  - Against H1's predicted −6.5% effect, that is a ratio of about 7
    standard errors.

## 5. Does not count

- Changing classes, N, arms, metrics, margins, the direction of H1 or
  the correction after data.
- Dropping a run, curve or model.
- Re-running a session because of its result.
- Citing wall time from this L0 host.

Scope: 24-bit `p`, depth 3, split craters, ℓ ∈ {2, 3}, the three methods
at v1's parameters. No security statement. An H1 effect is an IC-yield
statement; it becomes an improvement only if it also shows in IC cost S
and survives the §5a band (CLAUDE.md success criterion).
