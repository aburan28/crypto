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
| n41 blocks 3–4 | 2,048 | −0.01623388189548811 | 0.769615192403798 | 0.9777371149740774 | fail |
| n53 blocks 0–4 | 5,120 | −0.005156678549466277 | 0.6471764117941029 | 0.9864451082897685 | fail |

Both cohorts miss every numeric bar. The quartile ratios near 1 are probe-count descriptions inside an already sealed panel. They are not an index-calculus speedup.

## Scope

A/B traces matched on public Q, base hash, probes, and recovered scalar, with `group_verified` true and 1,024 solved rows in every consumed block. Feature evaluations: 16,384. Permutations use Python 3.12 `random.Random(seed).shuffle` for 2,000 successive Fisher–Yates shuffles of the archive-order probe list. Seeds are the frozen per-cohort values 2026100141 and 2026100153.

Producer wall time 12.64 s and CPU 12.64 s on Linux x86_64, Python 3.12.8. The checker finished in 13.17 s with peak resident set 137,868 KiB, inside the 120 s and 512 MiB limits. No input mismatch and no censoring.

This rejects these three public-coordinate predictors under this split. It does not reject every PDP difficulty discriminator. A charged relation-selection pilot on new Q is not authorized.
