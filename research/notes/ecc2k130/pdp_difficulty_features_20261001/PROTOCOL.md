# Frozen public-point feature test for compact PDP difficulty

**Status: protocol only.** No x/y feature correlation, feature choice, selected-quartile probe ratio, or permutation outcome has been computed. The archived per-orbit witness frequencies and aggregate probe quantiles were inspected before this freeze; they are not inputs to the feature choice below. Commit and merge this protocol and its CONFIG before running the scored producer. Preserve any run failure or changed archive in a follow-on outcome PR.

## Question, boundary, and scope

Can a feature computable from a public subgroup point Q before the compact four-point PDP distinguish targets with fewer solver probes on genuinely held-out Q? The reference is the frozen n41/K255/L1,024 and n53/K440/L1,024 W64 compact source process in the [six-cell cold outcome](../disjoint_cold_v2_outcome_20261001/RESULT.md). That outcome recovered every Q but cost 1.3137× and 1.6054× matched rho CPU, respectively. Its full CPU gaps and the [selector phase admission](../selector_phase_admission_20261001/DECISION.md) are boundaries for any later charged policy. The present test is a **difficulty-prediction diagnostic**, not a new solver, a factor-base improvement, a speedup or a native degree-263 PDP result.

Use only the two sealed v2 batch tarballs, their SHA-256 values, the archived manifest and exact source head in CONFIG. Rehash each full tar and every consumed target/base/run member against the manifest before scoring. For each block, read both `ic_a` and `ic_b` traces, require identical public Q, base hash, integer `probes`, independently verified scalar and success for each fixture index, and require 1,024 complete point-only rows. Use the A trace once for statistics; B is a same-Q deterministic replication, not another independent sample. The old scalar labels and recovered scalar must not enter any predictor. The v2 input freeze already proves Q orbit disjointness; this experiment does not resample or replace Q.

The three candidate features, in fixed tie order, are:

1. `x_weight`: `bit_count(x)/n`.
2. `signed_y_weight`: `min(bit_count(y), bit_count(x XOR y))/n`, invariant under point negation.
3. `frobenius_x_min_weight`: the minimum `bit_count` over `x, x², ..., x^(2^(n−1))`, divided by n, using the pinned bit-polynomial moduli in CONFIG. This is invariant under the source's geometric Frobenius. Verify the next square returns to x after n steps.

All features use the public affine coordinates exactly as archived, with zero-based polynomial-basis bits. They require no group log, rank row, oracle output, fixture scalar, or target probe. Record feature evaluation counts and elapsed CPU as diagnostics. A later charged relation selector would pay for them; this study does not compare its end-to-end cost.

## Fixed split and decision

Train on **n41 blocks 0, 1, 2** (3,072 unique Q). Compute the average-tie Spearman correlation of each feature with integer `probes`; choose the largest absolute correlation, resolving exact ties by the feature order above. Freeze its sign as the direction: a positive correlation selects lower feature values as easier, a negative one selects higher values. If all three correlations are zero, record `NO_PREDICTOR` and still evaluate the first feature by the fixed tie rule. Do not inspect n41 blocks 3 or 4 or any n53 feature/probe pairs before this choice is recorded in the raw output.

Evaluate that *one* feature and sign on n41 blocks 3, 4 (2,048 unique Q) and all five n53 blocks (5,120 unique Q) without refitting. For each evaluation cohort, compute signed Spearman rho (`training_sign × observed_rho`), and order Q by predicted ease, breaking feature ties by public `(x,y)` integers only. Select the first floor(N/4) points and report `median(probes_selected)/median(probes_all)` plus every raw observation. Compute a one-sided permutation p-value for the signed rho from 2,000 independent label shuffles using the fixed per-cohort RNG seeds in CONFIG and `p=(1+number of shuffled signed rhos ≥ observed)/(2001)`. Retain all 2,000 scores or their SHA-256 plus exact seed and algorithm. Use a separately written checker to recompute every feature, split, tie, rank, p-value and A/B equality.

A **transfer-worthy predictor** requires *both* evaluation cohorts to have signed rho at least 0.10, permutation p at most 0.01, and selected-quarter median probe ratio at most 0.80. Otherwise the result is negative or mixed for these three coordinate features under this fixed split. Preserve the n41 training correlations, both holdout cohorts and failures; no alternate threshold, feature, split or later seed may replace them. A positive result only permits a separately frozen, charged relation-selection pilot on new orbit-disjoint Q and matched rho. It cannot lower the cost of a fixed public Q by abstaining, and it cannot establish n83 or n131 transfer. A negative result rejects these cheap coordinate-only predictors, not every PDP difficulty discriminator.

The producer and independent checker must each finish within 120 seconds and 512 MiB on a normal CPU host. An input/hash/replication mismatch is `FAIL`; a resource stop is `CENSORED`; neither is a negative correlation result. Archive source hashes, exact commands, host/interpreter, complete observations, train choice, holdout metrics, elapsed CPU, independent receipt and the narrow decision in a linked outcome PR. The canonical scoreboard changes only when the verified outcome exists.

The non-scoring input check is reproducible before the outcome:

```sh
python3 research/notes/ecc2k130/pdp_difficulty_features_20261001/check_config.py \
  --out /tmp/ecc2k130-pdp-difficulty-input-check.json
cmp /tmp/ecc2k130-pdp-difficulty-input-check.json \
  research/notes/ecc2k130/pdp_difficulty_features_20261001/INPUT_CHECK.json
```
