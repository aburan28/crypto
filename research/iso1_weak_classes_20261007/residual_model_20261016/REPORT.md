# ISO-1 continuation (2026-10-16): certified labels, the 8th-power mass, and why the residual zeros have no compact criterion

This round answers the two open items in the [main report](../REPORT.md):
the zero labels were probabilistic, and the high-depth classes with no weak
representative had no explanation beyond an association with class number.

1. **Every census label through p = 83 is now certified.** A new counter
   returns a group order only when it is the unique Hasse-interval candidate
   compatible with the tested points on the curve and its twist. Reruns at
   p = 7, 13, 37, 41, 43, 47, 53 and 59 reproduce every retained CSV byte for
   byte. p = 61, 67, 71, 73, 79 and 83 are new. The queue continues through p = 199.
2. **The number of 8th-power Legendre parameters in a class is an exact
   2-adic formula times a Hurwitz class number.** The formula reproduces
   exhaustive censuses at p = 7, 11, 13, 17 and the registered holdouts
   p = 19 and 23 (0 mismatches in 13,146 rows).
3. **Beyond that 2-adic structure, the norm-one torus samples a class at
   random.** A parameter-free binomial thinning predicts the residual zero
   count. It passes at the pristine holdouts p = 53 and 59, and at p = 61,
   67, 71, 73, 79 and 83, which were predicted before their censuses finished. The dispersion is
   1.00 ± 0.01 throughout.

The consequence for the original request: past the proved obstruction
`v₂(f_π) ≥ 2`, there is no compact deterministic criterion to find. Whether a
class is weak is a draw. Its probability `1 − (1 − 8/(p²−1))^(N₈(t)/12)` is
computed from `t` alone, and the evidence does not support anything sharper.

## 1. Certified point counts

The original counter accepted a Hasse-interval multiple `m` of one random
point's order once a second random point was also killed by `m`. That is
not a proof: when the group exponent is below the interval width, several
candidates survive. It mislabeled two representatives at p = 7 and one each
at p = 11 and 17.

`certified_order` in
[`point_count_kernel.rs`](../census_runtime/src/point_count_kernel.rs)
keeps the whole candidate set. It works as follows:

- Because `Q = p⁶`, the Hasse interval is exactly `[Q+1−2p³, Q+1+2p³]`.
- Every curve here has full rational 2-torsion, so `4 | #E`, and the twist
  has order `2Q+2−#E`.
- A full baby-step/giant-step scan of `4P` lists every `N ≡ 0 (mod 4)` in
  the interval with `N·P = O`, not just the first hit.
- Each further point, alternating between the curve and its twist, filters
  that list.
- The order is returned only when exactly one candidate remains. The true
  order passes every filter, so a returned order is proved.

The census aborts if 16 points on each side leave two or more candidates.
That never happened. The scan uses stride 4 and the absolute-Frobenius
quotient, so the new runs are also much faster: p = 37 takes 32 s on four
threads, against 1,107 s for the original direct run.

| p | point counts | extra isolation points | wall s (4 threads, shared host) | result |
|---:|---:|---:|---:|:--|
| 7 | 213 | 5 | 0.01 | identical to the exact square-table census |
| 13 | 2,423 | 8 | 0.2 | identical to the retained CSV |
| 37 | 156,523 | 16 | 32 | identical |
| 41 | 235,901 | 20 | 56 | identical |
| 43 | 285,363 | 27 | 110 | identical |
| 47 | 407,193 | 35 | 118 | identical |
| 53 | 658,243 | 27 | 347 | identical |
| 59 | 1,010,651 | 23 | 941 | identical to the retained p59 run |
| **61** | **1,154,751** | **41** | **1,265** | **new**: 111,084 / 223,260 ordinary rows weak; 0 / 111,630 depth-1 weak; 546 / 111,630 high-depth zero |
| **67** | **1,680,383** | **30** | **1,448** | **new**: 147,496 / 296,274 ordinary rows weak; 0 / 148,138 depth-1 weak; 640 / 148,136 high-depth zero |
| **71** | **2,118,901** | **34** | resumed\* | **new**: 175,800 / 352,870 ordinary rows weak; 0 / 176,436 depth-1 weak; 634 / 176,434 high-depth zero |
| **73** | **2,367,853** | **27** | resumed\* | **new**: 191,222 / 383,688 ordinary rows weak; 0 / 191,844 depth-1 weak; 622 / 191,844 high-depth zero |
| **79** | **3,247,401** | **35** | resumed\* | **new**: 242,652 / 486,798 ordinary rows weak; 0 / 243,400 depth-1 weak; 746 / 243,398 high-depth zero |
| **83** | **3,956,583** | **33** | resumed\* | **new**: 281,650 / 564,898 ordinary rows weak; 0 / 282,450 depth-1 weak; 798 / 282,448 high-depth zero |

"Extra isolation points" counts the points needed beyond the first, summed
over all point counts. Almost every order is isolated by one scan.

Evidence: [comparison receipt](../certified_20261016/comparison_receipt.txt),
per-prime stderr receipts and CSV hashes, the
[queue script](../certified_20261016/queue.sh) and its
[status log](../certified_20261016/queue_status.tsv), and the
[p = 61](../certified_20261016/p61_cert.csv.gz) and
[p = 67](../certified_20261016/p67_cert.csv.gz) and
[p = 71](../certified_20261016/p71_cert.csv.gz) CSVs. \*p = 71, 73, 79 and 83 ran in
several container sessions with `ISO1_CHECKPOINT`; its job log is
[p71_cert.ckpt.gz](../certified_20261016/p71_cert.ckpt.gz), and its
wall time is not a single measurement. The unit test
`certified_order_matches_exhaustive_count_at_p7` compares 12 random p = 7
curves with brute-force counts.

The retained p = 37–59 zero labels are therefore exact. They no longer
depend on the randomized counter. The remaining assumption is the
orbit-quotient enumeration, which is proved in [THEOREM.pdf](../THEOREM.pdf),
Theorem 6.

## 2. An exhaustive λ census at small p

`iso1_lambda_census` ([source](../census_runtime/src/lambda_census.rs))
takes every `λ ∈ F_Q \ {0,1}`, or every element of a subgroup
`(F_Q^*)^(2^k)`, and records three things:

- the certified trace of the Legendre curve `y² = x(x−1)(x−λ)`;
- whether `λ` is a `2^k`-th power;
- whether `λ` lies in the norm-one torus `T` of order `q²+q+1`, whose
  elements give exactly the weak curves.

The results:

* At p = 11 the census agrees with the orbit census on every weak count
  (666 rows). This is an independent certification by a different
  enumeration.
* The total number of `λ` with `|trace| = t` is exactly `6·H(|t²−4Q|/4)`.
  This checks the Hurwitz computation.
* `T` has odd order, so it lies in `Ω = (F_Q^*)^(2^e)` with `e = v₂(p²−1)`.
  At p = 11, membership in `T` thins `Ω` at random inside each class: the
  observed rate is 0.0667 against `|T|/|Ω|` = 0.0667, the binomial
  dispersion is 1.012, and there are 31 zero classes against 32.4
  expected. At p = 7 the same test gives dispersion 0.65 and 7 zeros
  against 9.3.
* `Ω` itself is not 2-adically local when `e ≥ 4`. At p = 17 (`e = 5`),
  32nd-power counts scatter inside every 2-adic cell. **8th powers are
  local** at every prime tested (`e` = 3, 4, 5). This matches the fact that
  the Kummer cover `λ = μᴺ` of the λ-line is a congruence modular curve
  only for `N | 8`. The fact is quoted, not reproved here.

The model below therefore uses 8th powers, `|T| / |(F_Q^*)^8| = 8/(q−1)`.

## 3. The 8th-power mass formula (conjecture, exact on six primes)

Write `Δ = t² − 4Q = f_π² D_K`, `k = v₂(f_π)`, and `s = (D_K/2)`. Let
`H_odd(t) = (h_K/(w_K/2)) · Σ_{d | f_odd} d ∏_{ℓ|d}(1 − (D_K/ℓ)/ℓ)` be the
Hurwitz sum with the 2-part of the conductor removed. For `k ≥ 2`:

```
N₈(t) = H_odd(t) · (2 − s) · 2^(k−2)                      if t ≢ ±(Q+1) (mod 32)
N₈(t) = H_odd(t) · [5 (2 − s) · 2^(k−2) − 6 (1 − s)]      if t ≡ ±(Q+1) (mod 32)
```

and `N₈(t) = 0` for `k = 1`, which is the proved theorem. Here `N₈(t)` is
the number of 8th powers `λ` with `|trace(L_λ)| = t`.

The weighted mass of curves at 2-volcano level `j` in the class is
`m_j = H_odd · 2^j (1 − s/2)` for `j ≥ 1`, and `m_0 = H_odd`. In those
terms the formula equals:

- for `v = 4`, one 8th power per curve at level `k−1`, and none elsewhere;
- for `v ≥ 5`, six per curve at levels `≤ k−2`, 2.5 per curve at level
  `k−1`, and none on the floor.

Either way the floor contributes nothing, as the earlier curve-level bound
requires. The per-level reading is consistent with the class totals but was
not checked curve by curve.

| test | primes | rows | mismatches |
|:--|:--|---:|---:|
| derivation (exhaustive 8th-power census) | 7, 11, 13, 17 (e = 4, 3, 3, 5) | 4,078 | 0 |
| registered exactness holdout | 19 | 3,249 | 0 |
| registered exactness holdout | 23 (e = 4) | 5,819 | 0 |
| consistency at large p: `(8/(q−1))·Σ N₈` against the exact weak total `q²+q` | 37–83 | all | off by 1 at each prime, from the `λ = 1` and special-orbit terms |

The formula was read off the exhaustive tables; it is not proved. A proof
would count `F_Q`-points on the congruence cover `λ = μ⁸` by Deuring's
correspondence, class by class.

## 4. Residual zero classes: a parameter-free thinning model

[PROTOCOL.md](PROTOCOL.md) registered **V1**: weak orbits are Poisson with
mean proportional to `H(|Δ|/16)`. V1 **failed** on its training primes
37–47, with dispersion about 5 and roughly twice the predicted zero count.
The diagnosis was the 2-adic factor in Section 3.

[PROTOCOL_V2.md](PROTOCOL_V2.md) registered **V2**, which uses only
small-prime information: weak orbits follow `Binomial(N₈(t)/12, 8/(q−1))`.
Nothing is fitted. Primes 37–47 were seen during the V1 diagnosis, so they
are reported but are not pristine.

| p | role | high-depth rows (t > 0) | zero rows observed | predicted [95%] | verdict | dispersion | AUC |
|---:|:--|---:|---:|:--|:--|---:|---:|
| 37 | seen | 12,321 | 145 | 168.8 [146.5, 191.1] | outside, by 1.5 | 0.981 | 0.959 |
| 41 | seen | 16,810 | 198 | 186.2 [162.8, 209.6] | inside | 0.998 | 0.969 |
| 43 | seen | 19,414 | 188 | 195.9 [171.9, 219.9] | inside | 1.003 | 0.972 |
| 47 | seen | 25,403 | 212 | 212.6 [187.6, 237.6] | inside | 0.997 | 0.978 |
| **53** | pristine holdout | 36,517 | **247** | 240.4 [213.8, 266.9] | **pass** | 1.005 | 0.982 |
| **59** | pristine holdout | 50,474 | **270** | 265.9 [238.0, 293.8] | **pass** | 1.001 | 0.986 |
| **61** | forward prediction | 55,815 | **273** | 275.3 [247.0, 303.6] | **pass** | 0.992 | 0.984 |
| **67** | forward prediction | 74,068 | **320** | 301.1 [271.5, 330.8] | **pass** | 0.996 | 0.987 |
| **71** | forward prediction (e = 4) | 88,217 | **317** | 317.9 [287.5, 348.3] | **pass** | 1.000 | 0.989 |
| **73** | forward prediction (e = 4) | 95,922 | **311** | 327.1 [296.3, 357.9] | **pass** | 1.011 | 0.989 |
| **79** | forward prediction (e = 5) | 121,699 | **373** | 353.5 [321.5, 385.5] | **pass** | 0.998 | 0.993 |
| **83** | forward prediction (e = 3) | 141,224 | **399** | 372.6 [339.7, 405.4] | **pass** | 1.001 | 0.993 |

**Pooled check.** Over the eight scored primes 53–83, observed minus
predicted zero rows is +56.2 with standard deviation 43.1 (z = 1.30, not
significant). Three of the four largest scored primes (67, 79, 83) lie above
their means and one (73) below; the pooled excess is recorded as a watch
item for the next forward predictions.
| 89, 97, 101, 103 | forward predictions | | *census running* | see [predictions_v2.txt](predictions_v2.txt) | | | |

![Observed residual zero classes against the registered V2 predictions](figures/zero_rows.svg)

The figure plots the recorded values in [figures/zero_rows.csv](figures/zero_rows.csv)
(source [plot_zero_rows.py](figures/plot_zero_rows.py)). Filled markers are
scored censuses; hollow markers are forward predictions whose census is
still running.

The predictions for p = 61, 67 and 71 were committed in `edc0c3ddf`, and
those for 73, 79 and 83 in `788e11c8c`. Both commits precede the census
output for those primes. Dispersion is the mean of
`(N_t − μ_t)²/Var` over rows with mean ≥ 5 weak orbits; a value of 1 means
binomial. AUC measures how well `P(zero)` ranks the zero rows. For
comparison, `h_K` alone gives 0.71–0.78.

The calibration bins (in [model_v2_results.txt](model_v2_results.txt)) show
no excess zeros in any band of predicted probability. Every observed zero
lies in a row with predicted `P(zero) > 10⁻⁶`. The earlier counterexample
pairs, each one weak and one zero class with equal `h_K` or equal trace mod
2¹⁷, are what this model expects:

| p | zero class t | P(zero) | weak class t | P(zero) | N₈ ratio |
|---:|---:|---:|---:|---:|---:|
| 37 | −92,218 | 0.168 | 38,854 | 0.0020 | 3.5 |
| 41 | 136,542 | 0.343 | 5,470 | 1.1·10⁻⁵ | 10.7 |
| 47 | 204,238 | 0.071 | 73,166 | 0.00016 | 3.3 |
| 41 | 133,826 | 0.263 | 134,082 | 0.263 | 1 |

In the first three pairs `H_odd` differs by a large factor while `h_K`
agrees. The near-edge p = 41 pair has identical predictors, so one zero
and one weak class is simply a draw with probability about 0.39.

**Scaling.** The mean number of weak orbits in a central class grows like
`p`, so only near-edge and small-`H_odd` classes can be empty. Predicted
zero rows (t > 0) are 275, 301, 318, 327, 354 and 373 at p = 61, 67, 71, 73,
79 and 83: roughly `4.5p`, a vanishing fraction of the `≈ p³/2` high-depth
rows.

## 5. What this settles and what it does not

* Proved (THEOREM.pdf): no ordinary class with `v₂(f_π) = 1` is weak, at
  every odd `p`.
* Exact for p ≤ 83: the class labels, certified.
* Conjectured, exact on six primes: the `N₈` formula.
* Measured, three pristine or forward holdouts passed: inside high-depth
  classes, weakness is a binomial thinning of the 8th-power mass. No
  further deterministic invariant is needed to explain the data, and none
  of the tested invariants (2-splitting, depth, trace residues to 2¹⁷,
  `h_K`) separates the residual zeros.
* Not shown: that no exotic arithmetic criterion exists. A future failure of
  the forward predictions would be evidence for one.
* Census status: p ≤ 83 complete and certified; the queue runs the remaining
  primes through 199 in order. At the measured rate, p = 199 alone needs
  about 4 days on 4 threads.

## Reproduction

```sh
cargo build --release --locked --manifest-path research/iso1_weak_classes_20261007/census_runtime/Cargo.toml --bins
cargo test  --release --locked --manifest-path research/iso1_weak_classes_20261007/census_runtime/Cargo.toml --lib
B=research/iso1_weak_classes_20261007/census_runtime/target/release
RAYON_NUM_THREADS=4 $B/iso1_class_census 61 p61.csv --derive-twists --absolute-orbit-quotient --certified
ISO1_SKIP_SYLOW=1 $B/iso1_lambda_census 11 lam11.csv                       # all lambda
ISO1_OMEGA_ONLY=1 ISO1_OMEGA_EXP=3 $B/iso1_lambda_census 17 o8_17.csv     # 8th powers only
python3 -c "p=19; print('\n'.join(str(t) for t in range(2,2*p**3,4) if t%p))" > t19.txt
ISO1_P=19 ISO1_TRACES=t19.txt gp -q -f research/iso1_weak_classes_20261007/residual_model_20261016/class_levels.gp > L19.csv
$B/iso1_residual_model exact om19.csv L19.csv
$B/iso1_residual_model holdout 61 p61.csv L61.csv
$B/iso1_residual_model predict 71 L71.csv
```

The model is evaluated natively by `iso1_residual_model`
([source](../census_runtime/src/residual_model.rs)). The first evaluations and
the committed predictions came from Python scripts, now kept unchanged in
[legacy_python/](legacy_python/README.md) as provenance; the native binary
reproduces every number they printed. A PDF of this report with the figure is
[REPORT.pdf](REPORT.pdf).

`class_levels.gp` uses `quadclassunit`, so its large-`|D_K|` class
numbers are GRH-conditional. That affects the model's predictions, not the
census labels. The small-prime census CSVs and the gzipped class-level
tables are in [data/](data/). This round used PARI/GP 2.15.4 on a 4-core
Linux container.
