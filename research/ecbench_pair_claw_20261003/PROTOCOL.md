# The pair claw of cryptanalysis#175 against strong rho, on the same targets

Registered 2026-10-03, before the measured session (`spec.json`) ran.
The scratch smoke and calibration runs that preceded it are described
under "Before this protocol" and are not evidence for any prediction.

## Question

[aburan28/cryptanalysis#175](https://github.com/aburan28/cryptanalysis/pull/175)
solved a fresh public target on `EC1N83Ckb1h876c2921cb64` with a
four-point signed-Frobenius quotient pair-table search over a
**known-log** orbit base (candidate
`IC1N83Ckb1fb8000204PDP4qtableRCdirectLAnoneTDdirectISO0h49b47d79e9e3`).
It did not run same-target rho in one unit and claims no speedup. This
protocol puts a native port of that search, `claw.pair_table`, into
`ecbench` beside `rho.signed_frobenius_strong`, on the same one-target
workloads, in `S = group-addition equivalents / √r`.

## The method, and why it is generic

`src/cryptanalysis/ecbench/claw.rs`. `S` seed scalars `a_i`, seed points
`[a_i]G` and all `2n` signed Frobenius conjugates (`B = 2nS` points, all
logs known); a target-independent table of `M` zero-pair quotient
descriptors `P_i + σφ^t P_j`, keyed by the least normal-basis rotation of
the abscissa; unique unordered query pairs `Q − (F_k + F_l)`, probed by the
same key. A hit is a natural four-point relation whose every term has a
known log, so the scalar follows without linear algebra (`LAnone`).

Every base log is known, so each table entry is a known multiple of `G`
and the method is a randomised baby-step giant-step on signed Frobenius
classes. The generic-group bound applies to it.

## Boundaries, derived before measuring

Per curve, with `A = 2n` automorphisms:

- **Floor:** `√(π/2A) = √(π/4n)`, the signed-Frobenius rho constant.
- **Reference:** `rho.signed_frobenius_strong` measured in this session on
  the same targets. Its A/A control is `rho-strong-aa`.
- **Claw expectation:** a query lands in the table with probability
  `2nM/r` and costs two additions; with `M = c·√(r/n)` the mean is
  `(c + 1/c)/√n` plus `O(S log r)` base set-up. Finite base support and
  duplicate classes are not in the formula.

| curve | log₂ r | n | floor | claw `c = 1` | / floor | claw PR shape | / floor |
|---|---:|---:|---:|---:|---:|---:|---:|
| `icv1-f2m43-tm998717-e2e742b0` | 32.11 | 43 | 0.1351 | 0.3050 | 2.26 | 4.885 | 36.1 |
| `icv1-f2m47-t22705043-f4e44623` | 36.64 | 47 | 0.1293 | 0.2917 | 2.26 | 4.672 | 36.1 |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.00 | 41 | 0.1384 | 0.3123 | 2.26 | 5.002 | 36.1 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.26 | 53 | 0.1217 | 0.2747 | 2.26 | 4.400 | 36.1 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.52 | 59 | 0.1154 | 0.2604 | 2.26 | 4.170 | 36.1 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.21 | 61 | 0.1135 | 0.2561 | 2.26 | 4.101 | 36.1 |

## Arms

| arm | role | method |
|---|---|---|
| `rho-strong` | reference | `rho.signed_frobenius_strong` (lanes 32, dp 4) |
| `claw-c1` | candidate | `claw.pair_table`, `base_scale = 4`, `table_scale = 1` (cost-minimising) |
| `claw-pr` | candidate | `claw.pair_table`, `base_scale = 1.06`, `table_scale = 1/32` (the PR's n=83 shape: 48,194 seeds, `2^31`–`2^33` tables) |
| `bsgs-neg` | baseline | `bsgs.negation`, the generic table method without the Frobenius fold |
| `rho-strong-aa` | control | the reference again |

Six Koblitz curves, eight public targets each (`target_seed 20261003175`),
three interleaved rounds after one warm-up, `--cpus none` on the macOS
development host. That host earns isolation level L0, so **operation
counts are the result and wall time is descriptive only**.

## Predictions

- **P1, correctness.** Every claw run that returns a candidate verifies
  as `[k]G = Q`; no wrong answers.
- **P2, the claw at `c = 1` sits where theory puts it.** Mean `S` of
  `claw-c1` is within `[0.8, 1.25] × 2/√n` on every curve.
- **P3, it does not beat rho.** On every curve, the mean `S` of `claw-c1`
  is at least `1.5×` that of `rho-strong`. **Falsifier:** any curve where
  the claw's paired ratio to `rho-strong` has a 95 % interval wholly
  below 1. That would contradict the generic bound, and the first
  suspect would be a charging bug.
- **P4, the PR's shape is an order of magnitude worse.** On every curve
  the mean `S` of verified `claw-pr` runs is at least `10×` that of
  `rho-strong`. Some `claw-pr` runs may exhaust the query-pair domain
  without a relation (finite base support; at this shape the expected
  query count is 12–17 % of the domain). Exhausted runs are counted and
  reported, never dropped, and any comparison that has one is marked
  incomplete.
- **P5, flat in r.** Across the six curves, a least-squares fit of
  `log₂(S·√n)` for `claw-c1` against `log₂ r` has slope within
  `[−0.05, 0.05]`. Total operations then grow as `r^{1/2}`, like rho.

## Inadmissible

Changing a parameter, curve, target seed or round count after the
session starts; dropping exhausted, failed or timed-out runs; quoting a
stage (table or queries alone) as the method's cost; reporting a `vs_rho`
claim. The claw is not an index-calculus candidate, and the
`ecbench claim` path is for `ic.pipeline` only.

## Before this protocol

Two scratch sessions ran during development and are not evidence:

- A 12-run smoke on `n = 53, 61`. Every run that finished verified. One
  `claw-pr` run on `n = 53` exhausted the query domain. Two `claw-pr`
  runs on `n = 61` hit at about 0.1 of the expected query count, which
  prompted the calibration below.
- A calibration on `icv1-f2m47-t22705043-f4e44623` (24 targets per arm)
  and `n = 61`. On `n = 47`, the mean of queries over expected queries,
  `r/(2nM)`, was 0.994 for `claw-c1` and 0.885 for `claw-pr`, with 24 of
  24 runs verified in each arm. The early `n = 61` hits were chance.

## n = 83 is not covered

`ecbench` handles `GF(2^m)` only up to `m ≤ 62` (one-word field elements),
so it cannot run the PR's own curve. Any statement this session makes
about n = 83 is an extrapolation from the six curves above under P5.
