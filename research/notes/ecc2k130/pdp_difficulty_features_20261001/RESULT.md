# Public-point PDP difficulty features: negative

**Decision: `NEGATIVE`.** Class: diagnostic. These three coordinate features do not predict compact four-point PDP probe cost on the frozen holdouts. This is not a solver change, not a factor-base improvement, and not a speedup. Common `S`, n=83, and n=131 stay unset. The sealed W64 cold CPU no-gos are unchanged.

The protocol and `CONFIG.json` were already on main (`c1a2e5f8d`, pull requests #1114 and #1113). The producer hashed both sealed tarballs and every manifest member, then wrote `TRAIN_CHOICE.json` from n41 blocks 0–2 only. Holdout target files were parsed after that file existed. An independent checker, which does not import the producer, recomputed the choice, both cohorts, and both 2,000-draw permutation hashes.

## Training choice

n41 blocks 0–2, 3,072 unique Q. Average-tie Spearman of each feature with integer A-trace probes. A zero-variance correlation is defined as 0. The largest absolute value is `signed_y_weight`, and its sign is positive, so lower weight was treated as easier.

| Feature | Training Spearman |
| --- | ---: |
| `x_weight` | −0.006680109605771484 |
| `signed_y_weight` | 0.01054714743356203 |
| `frobenius_x_min_weight` | −0.002242890351345753 |

## Holdouts

The same feature and sign, with no refit. Median is the average of the two central values. The transfer gate requires signed Spearman at least 0.10, permutation p at most 0.01, and a selected-quarter median probe ratio at most 0.80, on both cohorts.

| Cohort | Q | Signed Spearman | Permutation p | Selected-quarter median probe ratio | Gate |
| --- | ---: | ---: | ---: | ---: | --- |
| n41 blocks 3–4 | 2,048 | −0.01623388189548811 | 0.769615192403798 | 1.0131137541933517 | fail |
| n53 blocks 0–4 | 5,120 | −0.005156678549466277 | 0.6471764117941029 | 1.0258028379387603 | fail |

Both cohorts miss every numeric bar. The predicted-easy quarter is slightly more expensive than the cohort median, not cheaper. These ratios are probe counts inside an already sealed panel. They are not an index-calculus speedup.

The first draft of this replay sorted the opposite quarter: high signed-y weight first. Those withdrawn ratios were 0.9777 and 0.9864. The corrected order matches the merged producer in #1118, and the quartile ratios match the result recorded in #1120. The permutation p-values differ from that record by one draw. Both are far above 0.01, so the decision does not change. #1120 is the archive of that merged producer.

## Scope

A/B traces matched on public Q, base hash, probes, and recovered scalar, with `group_verified` true and 1,024 solved rows in every consumed block. Feature evaluations: 16,384. Permutations use Python 3.12 `random.Random(seed).shuffle` for 2,000 successive Fisher–Yates shuffles of the archive-order probe list. Seeds are the frozen per-cohort values 2026100141 and 2026100153.

The corrected producer run took 12.886850003007567 s wall and 12.885518794 s CPU on Linux x86_64, Python 3.12.8. An independent checker of that outcome passed in 13.017912187002366 s with peak resident set 137,944 KiB, inside the 120 s and 512 MiB limits. The withdrawn high-weight sort had quoted 12.64 s and 137,868 KiB; those figures are not this run. No input mismatch and no censoring.

This rejects these three public-coordinate predictors under this split. It does not reject every PDP difficulty discriminator. A charged relation-selection pilot on new Q is not authorized.
# Held-out public-point difficulty diagnostic: three coordinate features fail

**Decision: `NO_TRANSFER_WORTHY_PREDICTOR`.** This is the one-shot outcome of the [merged, frozen protocol](PROTOCOL.md) and [merged source lock](SOURCE_LOCK.json). The producer completed under its limits, and the separately written checker recomputed every feature, split, rank, permutation result, archive hash, A/B same-Q equality and decision. Complete observations and the two receipts are in [the committed evidence directory](evidence_run_20261001_source3a8c4699/RESULT.json). This is a predictor diagnostic on already solved public points, not a new factor base, a method speedup, or an ECC2K-130 logarithm.

The n41 training blocks contained 3,072 distinct Q. Their average-tie Spearman correlations with integer compact-PDP target probes were `x_weight = -0.006680109606`, `signed_y_weight = +0.010547147434`, and `frobenius_x_min_weight = -0.002242890351`. The frozen largest-absolute rule therefore selected **signed_y_weight, positive direction**, and wrote [train_choice.json](evidence_run_20261001_source3a8c4699/train_choice.json) before either evaluation cohort was read for feature/probe scoring. Within each fixed-degree cohort the source records the raw bit-count numerator; division by the positive constant n in the protocol leaves ranks, ties, selected Q, correlations and permutation p-values exactly unchanged.

| Independent evaluation cohort | Q | Signed Spearman (gate ≥ 0.10) | One-sided p (gate ≤ 0.01) | Selected-quarter median / all median (gate ≤ 0.80) | Probe medians selected / all |
|:--|--:|--:|--:|--:|--:|
| n41 blocks 3–4, K255/L1,024 | 2,048 | −0.016233881895 | 0.770114942529 | **1.013113754193** | 1,661 / 1,639.5 |
| n53 blocks 0–4, K440/L1,024 | 5,120 | −0.005156678549 | 0.646676661669 | **1.025802837939** | 13,735.5 / 13,390 |

The selected quartiles were slightly *harder* by the frozen median-probe measure on both holdouts. Neither cohort passes any of the three thresholds. All 2,000 seeded permutation scores per cohort, 10,240 individual Q observations, and both A/B target traces' upstream hash custody are retained. A post-run shape audit found 10,240 unique `(cell,x,y)` points and positive integer probes in every observation. The frozen source and raw archive manifest hashes appear in [RESULT.json](evidence_run_20261001_source3a8c4699/RESULT.json); the independent [VERIFY.json](evidence_run_20261001_source3a8c4699/VERIFY.json) binds that result by SHA-256 `808a422c9c292355681205c463afd21ed5560255582ec6008b15c226f4f95315`.

The producer ran on merged source commit `3a8c4699729cbc0aaf08ccb7d3bbc14cde590c04` with:

```sh
python3 -B research/notes/ecc2k130/pdp_difficulty_features_20261001/produce.py \
  --out research/notes/ecc2k130/pdp_difficulty_features_20261001/evidence_run_20261001_source3a8c4699
python3 -B research/notes/ecc2k130/pdp_difficulty_features_20261001/verify.py \
  --evidence research/notes/ecc2k130/pdp_difficulty_features_20261001/evidence_run_20261001_source3a8c4699 \
  --out research/notes/ecc2k130/pdp_difficulty_features_20261001/evidence_run_20261001_source3a8c4699/VERIFY.json
```

The producer reported 13.87 s wall, 9.19 s CPU, 73.1 MB peak RSS; 10,240 feature evaluations and 481,280 field squarings used 2.72 s feature CPU. The checker reported 11.56 s wall and 71.8 MB peak RSS. Both were below the prospectively fixed 120 s / 512 MiB limits. The producer receipt records the macOS host and Python 3.13.1 interpreter. The independent checker wrote `PASS` with the same narrow decision. The original W64 compact/rho n41 and n53 CPU ratios remain [1.314× and 1.605×](../disjoint_cold_v2_outcome_20261001/RESULT.md), respectively; those are prior complete-cost boundaries, not speed ratios for this feature diagnostic.

This rejects **only these three cheap coordinate-only predictors on this split and the frozen compact W64 source**. It does not prove that PDP difficulty is unlearnable, that a relation/factor-base structural score cannot help, or that a fixed Q becomes cheaper by abstaining. The [phase admission](../selector_phase_admission_20261001/DECISION.md) points to a target-blind, one-index n41 factor-base pilot: score candidate bases *before* constructing an index, include all sample and selection cost, then compare full cold rank plus target recovery on new orbit-disjoint Q against matched rho. A later natural high-arity descendant-native PDP pilot remains separate; this result supplies no n83/n131 or degree-263 transfer claim.
