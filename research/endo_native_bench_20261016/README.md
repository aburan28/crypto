# Native GLV-chain benchmark: Bandersnatch, Tom-384 and CP6-782 G1

Constructive measurement of GLV-2 through a cheap isogeny-chain endomorphism,
natively, on the three curves the autoresearcher's cross-curve sweep named as
the next ones to measure (`research/endosweep_curves_20261006`, finding 5:
**Tom-384 modelled 1.525, Bandersnatch 1.522, CP6-782 G1 1.494**), with a new
curve-generic harness (`endo-chain-bench`) modelled on the CryptoPro-B chain
sweep (`research/cryptopro_b_chain_sweep_20261006`, 1.408× counted, 1.417×
median).

The hypothesis, arms, units, predictions, success and stop conditions and
commands were committed in [`PROTOCOL.md`](PROTOCOL.md) (commit `33e2dd10b`,
harness in `85793758e`) before any counting or timing run on these curves,
and that file is unchanged.  Everything below comes from one counting run
and one timing run per curve at commit `33e2dd10b`.

* **Bandersnatch is measured on its short Weierstrass model**, as the model
  priced it.  Its deployed twisted-Edwards arithmetic is **not** measured.
  Nothing here says how GLV performs in an Edwards implementation.
* **Variable-time harness**: for measurement only, never for secret scalars.
* Class **engineering**: `k·P` on the same curve, computed faster.  This is
  not an attack and makes no security claim.
* Counted field operations (`M_eq = M + S`) are the primary result.  Wall
  time comes second, and its ratios hold within one curve only.  No claim is
  made about CPU cycles per field operation.

## Verdict

| preregistered condition | Bandersnatch | Tom-384 | CP6-782 G1 | |
|:--|:--|:--|:--|:--|
| 1. every correctness check passes in both runs | yes: every arm equals textbook `BigUint` `k·P` on 256/256 pairs, `φ(P) = λ·P` on 256 points + test vectors, decompositions valid | yes | yes | **met** |
| 2. `φ(P)` counts exact (the one declared exception at its declared 7M + 2S); totals within 1 % | 4/4 exact (`phi/0+1w/2/generic` = 7M + 2S as declared; the model charges 8M + 2S); worst total +0.15 % | 44/44 exact; worst −0.24 % | 114/114 exact; worst +0.31 % | **met** |
| 3. counted ratio ≥ 99 % of the model's | **1.525** ≥ 1.508 (model 1.523; sweep 1.522) | **1.521** ≥ 1.509 (model 1.524; sweep 1.525) | **1.492** ≥ 1.479 (model 1.494; sweep 1.494) | **met** |
| 4. declared GLV config beats the fastest baseline by median > 1.40, and timed/counted within ±5 % | **1.504** median, 1.528 min; −1.3 % | **1.504** median, 1.530 min; −1.2 % | **1.570** median, 1.520 min; +4.7 % | **met** (CP6-782 +4.7 % is near the 5 % limit, and see the noise note) |
| 5. the lowest `φ(P)` median is an optimised arm at the model's least count | `0+1w/2` optimised, 226 ns | `4+1w/7.5.5` optimised, 7 652 ns | `4+1w/7.5.3` optimised, 28 440 ns (its tie `4+1w/7.3.5`: 28 798) | **met** |

The declared GLV configuration is the grid chain at
`optimised/jacobian/w5`: `0+1w/2` (`√−2`, one 2-isogeny) on Bandersnatch,
`4+1w/7.5.5` on Tom-384 and `4+1w/7.3.5` on CP6-782.  In condition 4 the
fastest baseline by median is the reference itself on Bandersnatch
(`jacobian/w5`) and Tom-384 (`jacobian/w6`).  On CP6-782 it is `jacobian/w5`
(2 145 µs), not the counted cheapest `jacobian/w6` (2 175 µs).  The two
differ by 1.4 %, which is inside that run's 3.8 % A/A floor.

**Positive control** (run before the protocol and recorded in it): on the
frozen CryptoPro-B chain set, with the inputs of the CryptoPro-B sweep, the
generic harness's counting record equals that sweep's committed `counts.json`
on all 97 arms and every counted field: `mul`, `sqr`, `inv`, `add`, mean,
min and max, with 0 differences.  See
[`positive_control/`](positive_control/) (`compare_counts.py` and its output).
So the generic field, formulas, evaluators, recoding, tables and
decomposition reproduce the specialised harness exactly before they are used
on any new curve.

## Headline (counted primary, timed secondary)

| curve | `p` / `n` bits | best baseline (counted) | best GLV-2 (counted) | counted ratio | model ratio | declared GLV, median ratio vs fastest baseline | min ratio | A/A median ratio |
|:--|:--|:--|:--|--:|--:|--:|--:|--:|
| Bandersnatch (short Weierstrass) | 255 / 253, h = 4 | `jacobian/w5` 3 289.8 | `1+1w/3` opt/jac/w5 2 156.8 | **1.525** | 1.523 | **1.504** | 1.528 | 0.997 |
| Tom-384 | 384 / 384, h = 1 | `jacobian/w6` 4 929.6 | `4+1w/7.5.5` opt/jac/w5 3 240.4 | **1.521** | 1.524 | **1.504** | 1.530 | 1.001 |
| CP6-782 G1 | 782 / 377, h of 405 bits | `jacobian/w6` 4 462.4 | `4+1w/7.3.5` opt/jac/w5 2 991.5 | **1.492** | 1.494 | **1.570** | 1.520 | 1.038 |

For comparison, the CryptoPro-B sweep measured 1.408 counted and 1.417 by
median.  The ordering the model predicts holds when counted: GLV saves more
on the general-`a` curves Tom-384 and Bandersnatch, whose doubling costs
2M + 8S, than on CryptoPro-B (`a = −3`, 3M + 5S).  On Bandersnatch the
cheapest GLV arm by count is `1+1w/3` (2 156.8).  The declared `√−2` arm
`0+1w/2` counts 2 158.6, and the two medians are 93.8 and 94.1 µs.  The
difference between these two elements is within noise, in counts and in
time.

## Results

Full per-arm tables (every arm, counted and timed) are in
`<curve>.timing.stdout.md`, with per-round times in `<curve>.timing.json`.
The counting runs are `<curve>.counts.stdout.md` and `.json`, and their tables
show the model column and the deviation.  The rows below are verbatim from
those files.  Units: counted M_eq per `k·P` (mean over 256 pairs); ns per
`k·P` (median and min over 20 rounds); ratios are reference / arm.

### Bandersnatch (short Weierstrass model; reference `jacobian/w5`)

| arm | counted M_eq | model | ns median | ns min | ratio, median | ratio, counted |
|:--|--:|--:|--:|--:|--:|--:|
| baseline jacobian/w5 (reference) | 3 289.8 | 3 285.3 | 141 443 | 136 457 | 1.000 | 1.000 |
| baseline jacobian/w6 | 3 317.4 | 3 317.3 | 143 020 | 138 086 | 0.989 | 0.992 |
| baseline affine/w5 | 3 554.8 | 3 551.4 | 156 618 | 151 167 | 0.903 | 0.925 |
| GLV `0+1w/2` optimised/jacobian/w5 (declared) | 2 158.6 | 2 157.6 | 94 059 | 89 303 | 1.504 | 1.524 |
| GLV `0+1w/2` optimised/jacobian/w4 | 2 168.2 | 2 167.3 | 94 103 | 89 850 | 1.503 | 1.517 |
| GLV `0+1w/2` generic/jacobian/w5 | 2 161.6 | 2 161.6 (model) / 2 160.6 (declared) | 94 432 | 89 860 | 1.498 | 1.522 |
| GLV `0+1w/2` optimised/affine/w4 | 2 390.1 | 2 389.4 | 106 325 | 102 892 | 1.330 | 1.376 |
| GLV `1+1w/3` optimised/jacobian/w5 | 2 156.8 | 2 157.2 | 93 757 | 89 312 | 1.509 | 1.525 |
| A/A control | 3 289.8 | 3 285.3 | 141 821 | 136 115 | 0.997 | 1.000 |

Note on the `generic/jacobian/w5` row: its mean count, 2 161.6, equals the
model's value and not the declared model − 1.  This is not a contradiction.
The arm's chain costs exactly 1 M_eq less than the model's chain, which is
the declared exception, and the remaining difference is sampling: the model
averages 400 other scalars.  The same arm and the optimised one differ by
exactly 3.0 M_eq here, and 3 is the declared chain difference (9 against 6).
All 30 other totals lie within 0.15 % of the model.

### Tom-384 (reference `jacobian/w6`)

| arm | counted M_eq | model | ns median | ns min | ratio, median | ratio, counted |
|:--|--:|--:|--:|--:|--:|--:|
| baseline jacobian/w6 (reference) | 4 929.6 | 4 930.5 | 466 490 | 450 769 | 1.000 | 1.000 |
| baseline jacobian/w5 | 4 955.6 | 4 946.0 | 479 320 | 453 846 | 0.973 | 0.995 |
| baseline affine/w5 | 5 375.4 | 5 367.9 | 518 847 | 492 562 | 0.899 | 0.917 |
| GLV `4+1w/7.5.5` optimised/jacobian/w5 (declared) | 3 240.4 | 3 234.6 | 310 273 | 294 598 | 1.503 | 1.521 |
| GLV `4+1w/7.5.5` optimised/jacobian/w4 | 3 317.1 | 3 313.7 | 319 141 | 303 809 | 1.462 | 1.486 |
| GLV `4+1w/7.5.5` generic/jacobian/w5 | 3 288.4 | 3 282.6 | 314 279 | 297 399 | 1.484 | 1.499 |
| GLV `4+1w/7.5.5` optimised/affine/w4 | 3 672.5 | 3 670.8 | 354 598 | 336 754 | 1.316 | 1.342 |
| GLV `9+1w/7.5.7` optimised/jacobian/w5 | 3 248.4 | 3 251.6 | 308 878 | 293 613 | 1.510 | 1.518 |
| GLV `37+3w/23.5.5.5` optimised/jacobian/w5 (dearest element) | 3 309.3 | 3 317.1 | 320 911 | 300 420 | 1.454 | 1.490 |
| A/A control | 4 929.6 | 4 930.5 | 466 007 | 453 229 | 1.001 | 1.000 |

Across the 11 elements at the declared configuration, the counted totals
span 3 240–3 309 M_eq (2.1 %).  The medians span 308.9–333.2 µs.  Two
element arms (`15+2w/7.5.5.5` 332.2 µs and `25+1w/23.5.7` 333.2 µs) sit
about 6 % above their counted position, but their minima (295.8 and 302.6
µs) fall in the counted order.  This reads as round-to-round contention, not
as cost.  The order inside the pack is not established by wall time.

### CP6-782 G1 (reference `jacobian/w6`; inputs `r·G` with `G` of order `n`, checked)

| arm | counted M_eq | model | ns median | ns min | ratio, median | ratio, counted |
|:--|--:|--:|--:|--:|--:|--:|
| baseline jacobian/w6 (reference) | 4 462.4 | 4 463.9 | 2 175 424 | 1 990 214 | 1.000 | 1.000 |
| baseline jacobian/w5 (fastest baseline by median) | 4 483.3 | 4 479.5 | 2 145 010 | 2 024 679 | 1.014 | 0.995 |
| baseline affine/w5 | 5 398.7 | 5 396.2 | 2 618 756 | 2 416 705 | 0.831 | 0.827 |
| GLV `4+1w/7.3.5` optimised/jacobian/w5 (declared) | 2 991.5 | 2 987.6 | 1 366 596 | 1 332 371 | 1.592 | 1.492 |
| GLV `4+1w/7.3.5` optimised/jacobian/w4 | 3 069.0 | 3 068.2 | 1 438 455 | 1 371 838 | 1.512 | 1.454 |
| GLV `4+1w/7.3.5` generic/jacobian/w5 | 3 039.5 | 3 035.6 | 1 439 407 | 1 342 397 | 1.511 | 1.468 |
| GLV `4+1w/7.3.5` optimised/affine/w4 | 3 919.0 | 3 917.8 | 1 901 864 | 1 736 209 | 1.144 | 1.139 |
| GLV `9+1w/7.5.5` optimised/jacobian/w5 | 2 993.0 | 2 993.1 | 1 401 138 | 1 322 906 | 1.553 | 1.491 |
| A/A control | 4 462.4 | 4 463.9 | 2 096 275 | 1 999 904 | 1.038 | 1.000 |

**CP6-782's timing run was the noisiest of the three.**  The A/A median ratio
is 1.038 (a 3.8 % floor; the other curves' floors are 0.3 % and 0.1 %).  The
per-round spread is 49 % on the reference and 54 % on the A/A arm.  The
isolation record shows 11.1 CPU-seconds of other activity during the
570-second run, 10.3 s of it from the agent runtime.  The run was still
admitted as uncontended under the declared 0.5-CPU budget.  The declared
arm's median ratio, 1.570 against the fastest baseline, exceeds its counted
ratio of 1.499 by 4.7 %.  That excess is inside twice the floor, so per the
protocol it is not called a difference.  The minimum ratio is 1.520 against
`jacobian/w5`, and 1.494 against the counted reference `jacobian/w6`, which
is the counted value.  Read CP6-782's wall time as consistent with its
counts, not as a larger gain.  Across its 18 elements the counted totals
span 2 991–3 055 M_eq.

### The grid (counted M_eq / median µs), declared chain of each curve

| curve, evaluator, table | w = 3 | w = 4 | w = 5 | w = 6 | w = 7 |
|:--|--:|--:|--:|--:|--:|
| Bandersnatch `0+1w/2`, optimised, jacobian | 2 309 / 98.3 | 2 168 / 94.1 | 2 159 / 94.1 | 2 315 / 99.2 | 2 753 / 114.7 |
| Bandersnatch `0+1w/2`, optimised, affine | 2 440 / 108.8 | 2 390 / 106.3 | 2 479 / 109.0 | 2 777 / 120.2 | 3 461 / 144.5 |
| Bandersnatch baseline, jacobian | 3 538 / 153.5 | 3 366 / 148.4 | 3 290 / 141.4 | 3 317 / 143.0 | 3 495 / 149.8 |
| Bandersnatch baseline, affine | 3 656 / 161.4 | 3 560 / 156.6 | 3 555 / 156.6 | 3 669 / 160.8 | 3 981 / 171.7 |
| Tom-384 `4+1w/7.5.5`, optimised, jacobian | 3 562 / 374.0 | 3 317 / 319.1 | 3 240 / 310.3 | 3 346 / 320.1 | 3 741 / 356.9 |
| Tom-384 `4+1w/7.5.5`, optimised, affine | 3 795 / 367.4 | 3 673 / 354.6 | 3 715 / 364.7 | 3 979 / 378.2 | 4 633 / 446.2 |
| Tom-384 baseline, jacobian | 5 379 / 510.1 | 5 100 / 499.1 | 4 956 / 479.3 | 4 930 / 466.5 | 5 073 / 533.9 |
| Tom-384 baseline, affine | 5 597 / 533.9 | 5 427 / 527.4 | 5 375 / 518.8 | 5 452 / 515.5 | 5 741 / 563.1 |
| CP6-782 `4+1w/7.3.5`, optimised, jacobian | 3 311 / 1 534 | 3 069 / 1 438 | 2 992 / 1 367 | 3 100 / 1 461 | 3 503 / 1 677 |
| CP6-782 `4+1w/7.3.5`, optimised, affine | 4 038 / 1 985 | 3 919 / 1 902 | 3 961 / 1 895 | 4 226 / 2 024 | 4 887 / 2 241 |
| CP6-782 baseline, jacobian | 4 894 / 2 312 | 4 625 / 2 200 | 4 483 / 2 145 | 4 462 / 2 175 | 4 608 / 2 303 |
| CP6-782 baseline, affine | 5 612 / 2 664 | 5 450 / 2 617 | 5 399 / 2 619 | 5 479 / 2 628 | 5 771 / 2 736 |

* On all three curves **Jacobian tables beat affine tables** for both arms
  at every width ≥ 4.  The reason is the Fermat inversion: it costs 419,
  684 and 1 173 M_eq here (255 S + 164 M, 384 S + 300 M, 782 S + 391 M).  On
  CP6-782, with its 782-bit field, the affine tables lose 0.9–1.4 k M_eq.
  The conclusion depends on that inversion, as on CryptoPro-B.
* GLV's optimum is at width 5 everywhere (4 is close behind).  The
  baseline's optimum is at width 5 on Bandersnatch and at 6 on Tom-384 and
  CP6-782 by count, with 5 and 6 within the A/A floor by time on CP6-782.

### `φ(P)` alone (stage diagnostics, never a speed)

Bandersnatch, both chains:

| chain | generic: counted M + S | generic: ns median / min | optimised: counted M + S | optimised: ns median / min |
|:--|--:|--:|--:|--:|
| `0+1w/2` (`√−2`, 2-isogeny) | 7 + 2 | 303 / 289 | 5 + 1 | 226 / 221 |
| `1+1w/3` (3-isogeny) | 24 + 2 | 925 / 857 | 6 + 0 | 265 / 256 |

Tom-384 (22 chains, the same set and counts as CryptoPro-B, which has the
same `D_K = −619`), the five lowest by optimised median:

| chain | generic: counted M + S | generic: ns median / min | optimised: counted M + S | optimised: ns median / min |
|:--|--:|--:|--:|--:|
| `4+1w/7.5.5` | 124 + 4 | 11 737 / 11 325 | 78 + 2 | 7 652 / 7 343 |
| `4+1w/5.7.5` | 124 + 4 | 11 597 / 11 348 | 87 + 2 | 8 311 / 8 070 |
| `4+1w/5.5.7` | 124 + 4 | 11 859 / 11 298 | 87 + 2 | 8 579 / 8 114 |
| `9+1w/7.7.5` | 139 + 4 | 13 398 / 12 633 | 93 + 2 | 9 022 / 8 788 |
| `9+1w/7.5.7` | 139 + 4 | 13 160 / 12 796 | 93 + 2 | 9 119 / 8 701 |

CP6-782 (57 chains), the five lowest by optimised median:

| chain | generic: counted M + S | generic: ns median / min | optimised: counted M + S | optimised: ns median / min |
|:--|--:|--:|--:|--:|
| `4+1w/7.5.3` | 109 + 4 | 49 642 / 47 669 | 63 + 2 | 28 440 / 27 722 |
| `4+1w/7.3.5` | 109 + 4 | 50 105 / 47 836 | 63 + 2 | 28 798 / 27 808 |
| `4+1w/5.7.3` | 109 + 4 | 49 483 / 47 667 | 72 + 2 | 32 785 / 31 411 |
| `4+1w/5.3.7` | 109 + 4 | 49 966 / 47 972 | 72 + 2 | 33 208 / 31 466 |
| `0+1w/17.5` | 164 + 3 | 72 133 / 70 614 | 78 + 1 | 35 056 / 34 170 |

* Every `φ(P)` arm's counted `M` and `S` equals the chain set's model, 162
  of 162, apart from the one exception the protocol declared in advance.
  That exception is Bandersnatch `0+1w/2` under the generic evaluator.  The
  model charges it 8M + 2S: a Jacobian 2-step followed by the generic
  evaluator's projective tail.  That combination would not compute the right
  point.  The harness's projective 2-step costs 3M + 1S, so the arm counts
  7M + 2S, exactly as declared.  The optimised 2-step, 5M + 1S, matches the
  model.  Every count is identical on every input (min = max).
* The medians follow the counts.  The Spearman rank correlation between
  counted M_eq and median over the `φ(P)` arms is 0.949 on Bandersnatch (4
  arms), 0.984 on Tom-384 (44) and 0.994 on CP6-782 (114).
* The chain is a small part of the total: 0.2–0.3 % of the declared GLV
  arm on Bandersnatch, 2.5 % on Tom-384 and 2.2 % on CP6-782, by count.
* Decomposition stage (`num-bigint` Babai rounding; no field operations, so
  absent from the counts): 0.82–0.83 µs on Bandersnatch, 1.3–1.6 µs on
  Tom-384, 0.9–1.0 µs on CP6-782 (377-bit `n`).  That is under 1 % of the
  GLV total on every curve.

## Model against measurement

The frozen model files predicted every total arm's counted M_eq from the
same recoding and lattice bases, averaged over 400 other scalars.  The
counting runs land within 0.31 % of them on every arm of the three curves
(Bandersnatch worst +0.15 %, Tom-384 −0.24 %, CP6-782 +0.31 %).  They are
exact on every chain except the declared one.  The model's counted ratios
(1.523 / 1.524 / 1.494) and the counted ones (1.525 / 1.521 / 1.492) agree
to within 0.2 %.  The sweep README's rounded figures 1.522 / 1.525 / 1.494
come from an earlier sampling of the same model and agree to the same
precision.  The timed median ratios are 1.504, 1.504 and 1.570.  They lie
within 1.3 %, 1.2 % and 4.7 % of the counted ratio against the same
baseline, and the last of these is inside twice its A/A floor.

## Host, build and isolation

* Host: `Intel(R) Xeon(R) Processor @ 2.10GHz`, a cloud VM with 4 logical
  CPUs (`pclmulqdq sse4_2 popcnt aes avx bmi1 avx2 bmi2 avx512f adx`),
  kernel `6.18.44-fc-v139`, Linux x86-64.  The bench's host line says
  `logical cpus: 1` because the bench ran pinned to CPU 3.
* Build: `rustc 1.97.0 (2d8144b78 2026-07-07)`, `release` profile, no
  `target-cpu` flags.  The counting binary is the same source built with
  `--features endo-chain-opcount` into `target/opcount`.  The field has a
  run-time modulus over the crate's generic `Uint<L>::mont_mul` (4, 6 and 13
  limbs), with no dedicated squaring and no assembly.
* Commit at run time: `33e2dd10bd8b75a13122683576bf5a81f1583577` for all six
  runs.  The records say `git_dirty: true` only because the runs' own output
  files were untracked at that moment.  No tracked file differed from the
  commit.  The positive-control record was made from the same source before
  it was committed (its `git_commit` is the parent, `d5382a530`, and dirty).
* Timing ran with the host otherwise idle.  Before the timing runs, a
  low-priority background job of the orchestrator was paused through its
  pause file and left to settle for 95 s.  `uptime` then showed a load
  average of 0.46, and `top` showed 100 % idle with no user process busy.
  No build or other work of this session ran during timing.  The pause file
  was removed afterwards.
* Isolation: the three timing runs went through the native
  `isolated_bench run` on CPU 3 with `--max-other-cpu 0.5`, as declared.
  Each was admitted at the first attempt, with no retries.  The records
  (`<curve>.isolated_bench_record.jsonl`) show:

  | curve | wall s | `contended` | other CPU-s (agent runtime) | involuntary switches | major faults |
  |:--|--:|:--|--:|--:|--:|
  | Bandersnatch | 23.3 | false | 0.48 (0.45) | 55 | 0 |
  | Tom-384 | 98.9 | false | 1.79 (1.65) | 136 | 0 |
  | CP6-782 | 570.2 | false | 11.11 (10.28) | 621 | 0 |

  The counting runs (3 s, 14 s and 94 s) are deterministic and were not
  isolated.  Counts do not depend on contention.  One earlier attempt to
  launch the counting runs under `/usr/bin/time` failed because that tool is
  not installed on the host, so the bench never started.  That attempt is
  recorded here.
* Noise: A/A median ratios 0.997, 1.001 and 1.038 (minimum ratios 1.003,
  0.995 and 0.995).  Per-round spread (max/min − 1 over the 20 rounds): on
  the reference 38 %, 32 % and 49 %; on the declared GLV arm 85 %, 55 % and
  66 %.  Medians are the robust statistic.  The minimum ratios are reported
  as the least-contended ones.
* ns per counted M_eq across the total arms: 41.6–45.0 on Bandersnatch (4
  limbs), 94.4–105.3 on Tom-384 (6 limbs), 457–515 on CP6-782 (13 limbs).
  These are measured facts for this harness only, not a cost model.

## Honesty notes

* **Bandersnatch: short Weierstrass model only.**  Its deployed
  twisted-Edwards arithmetic (and arkworks' own GLV on it) is not measured.
  The 1.525 / 1.504 here compares two scalar multiplications on the
  Weierstrass model with general `a`, where a doubling costs 2M + 8S.
* **The baseline is a plain width-`w` NAF** with the same Jacobian
  arithmetic, not any production library's.  A faster field, a dedicated
  squaring or a cheaper inversion would move both arms.  The table
  conclusion depends on the Fermat inversion.
* **The chain sets are bounded by their generator.**  Chains outside them are
  not measured.  The cost-completeness argument (every chain within 2× of
  the cheapest is in the set) is the autoresearcher's, and it was not
  re-derived here.
* **CP6-782**: the GLV acts on the order-`n` subgroup (377-bit `n`, cofactor
  of 405 bits).  Every input point is `r·G` with `G` of order `n`, checked
  by textbook arithmetic.  Points outside the subgroup are not handled.
* **Counts are of field operations**, not instructions.  Additions are
  recorded in the JSON (`add`) but not in M_eq.
* **The generic evaluator's 2-step** differs from the model's pricing by
  1 M, as declared in advance.  The model's figure for it does not
  correspond to a correct evaluator.  This affects only the generic arms on
  Bandersnatch's `√−2` chain.
* **Lean**: the constants loaded here are certified by the autoresearcher's
  `bandersnatch_native_chain_certificate`, `tom384_native_chain_certificate`
  and `cp6_782_native_chain_certificate` (`NativeChains.lean`; primality is
  a hypothesis there).  Those theorems were not re-checked in this
  repository.  The native verification here is the loader's checks plus the
  correctness gate.  Timings and op counts are not Lean claims.

## Commands

Exactly what was run, in this order (harness `85793758e`, protocol
`33e2dd10b`):

```sh
cargo test --release --lib endo_generic                                  # 13 passed
cargo test --release --features endo-chain-opcount --lib endo_generic   # 16 passed
cargo clippy --release --lib --bins --tests [--features endo-chain-opcount]  # no warnings in the new files
cargo build --release --bin endo-chain-bench --bin isolated_bench
CARGO_TARGET_DIR=target/opcount cargo build --release \
  --features endo-chain-opcount --bin endo-chain-bench

# positive control (before the protocol)
C=research/cryptopro_b_chain_sweep_20261006
target/opcount/release/endo-chain-bench --count-only \
  --chains $C/frozen_inputs/chains.constants.json \
  --expect-sha256 683f92d0e56f5606f78b6c685f747a7b67106a7b5e5883e81869c8d6199d7336 \
  --zoo-generator gost-cryptopro-b --pairs 256 --seed 20261006 \
  --continuity-chain 4+1w/5.5.7 --continuity-config generic/affine/5 \
  --json research/endo_native_bench_20261016/positive_control/cryptopro_b.counts.json \
  > research/endo_native_bench_20261016/positive_control/cryptopro_b.counts.stdout.md
python3 research/endo_native_bench_20261016/positive_control/compare_counts.py

D=research/endo_native_bench_20261016; F=$D/frozen_inputs
# per curve C, SHA, REF:
#   bandersnatch 12100c87b6ffe5a3469f3c8ab625f2e8f23411df025a01bf3c06136cd646114b jacobian/5
#   tom384       c93b4b4752a3bad0331d6bdc6a7056c7a8608df16db626cb968c4e013a7ed00f jacobian/6
#   cp6_782      9e08fe6a45643b75c513c75164ecb51d39393d65c078dc07dadd9461bfa8ac90 jacobian/6
for each curve:
  target/opcount/release/endo-chain-bench --count-only \
    --chains $F/$C.chains.constants.json --expect-sha256 $SHA --model $F/$C.model.json \
    --reference $REF --pairs 256 --seed 20261016 \
    --json $D/$C.counts.json > $D/$C.counts.stdout.md
# then, host paused and idle:
for each curve:
  target/release/isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5 \
    --label endo-chain-bench-$C --out $D/$C.isolated_bench_record.jsonl -- \
    target/release/endo-chain-bench \
    --chains $F/$C.chains.constants.json --expect-sha256 $SHA --model $F/$C.model.json \
    --reference $REF --pairs 256 --rounds 20 --seed 20261016 --counts $D/$C.counts.json \
    --json $D/$C.timing.json > $D/$C.timing.stdout.md
```

## Files

| file | what |
|:--|:--|
| `PROTOCOL.md` | the protocol, committed before the runs |
| `README.md` | this record |
| `frozen_inputs/*.chains.constants.json`, `*.model.json`, `SHA256SUMS` | the frozen chain sets and model predictions (Python provenance, aburan28/crypto-autoresearcher `b2f0ac3095`) |
| `<curve>.counts.json`, `<curve>.counts.stdout.md` | the counting runs (with the model column) |
| `<curve>.timing.json`, `<curve>.timing.stdout.md` | the timing runs (every per-round time, host, build, correctness) |
| `<curve>.isolated_bench_record.jsonl` | the isolation records |
| `positive_control/` | the CryptoPro-B counting record of the generic harness and its arm-by-arm comparison |
| `src/ecc/endo_generic_field.rs`, `endo_generic_point.rs`, `endo_generic_chain_set.rs`, `src/bin/endo_chain_bench.rs` | the harness and its tests |
