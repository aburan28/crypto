# Does the compact-orbit lead survive a size and batch sweep against the strong rho?

Follow-on to `RESEARCH_STRONG_RHO_LADDER_20260929.md`. Written and committed
**before any sweep cell is run**. The ladder note fixed the reference at rung R3
(`examples/koblitz_rho_batch_ks_strong.rs`, `KIC_RHO_RUNG=3`, 32 lanes) and
scored one cell; at n = 53, L = 1,024, K = 440 IC/rho was measured before this
note was written (see the ladder note's results, added by the commit that
follows this one). This note asks the question that one cell cannot answer.

## Why a sweep

A constant-factor loss at one cell says nothing about scaling. Two facts keep
the lead alive as a *scaling* claim until they are tested against a fair
reference:

1. `RESULT.md` §3 of the PR #830 evidence reports in-process slopes on log₂ r of
   0.467 (IC) against 0.531 (batched rho) over n = 41→53 (three points, tuned K,
   wall clock, **unmatched** rho, "no interval"). A slope gap of about 0.06 would
   not be a constant factor: extrapolated over the 86 doublings from r ≈ 2^44 to
   ECC2K-130's r ≈ 2^130 it would be a factor near 2^5.
2. The cost model behind "both arms ~ √(L·r/n)" makes the *exponents* equal at
   optimal K and leaves the constants to decide. The IC constant is not known to
   be flat in n.

`AGENTS.md` §5 says an exponent claim needs at least four sizes, priced
whole-process, with the fit next to the prediction. This is that test, with the
reference the ladder produced.

## Frozen design

All IC runs are `koblitz_orbit_dlp_fast construct:<n>:0:<K> <scalars> 7 <out>`;
all rho runs are `KIC_RHO_RUNG=3 KIC_RHO_LANES=32 KIC_RHO_DP_BITS=4
KIC_RHO_BATCH_CORPUS=<corpus> koblitz_rho_batch_ks_strong <n> 0 signed_frobenius
<L> 531310`. The scalars fed to IC are the `published_fixture_scalar` values the
rho run derives for the same corpus, so both arms solve the same targets. The
binaries are the ones whose hashes the ladder note records; the IC source is
unchanged (sha256 `8dea682a…`).

- **Size sweep (A).** L = 1,024; n ∈ {37, 41, 53, 61}; IC's K is the value the
  PR #830 evidence chose for that n by lowest wall clock on a disjoint tuning
  corpus: K = 42, 255, 440, 600. Corpus names `n<n>-strong-sweep-L1024-v1`, except
  n = 53, which reuses the ladder's `n53-ks-growing-1024-v1` cell. K was tuned for
  *wall clock against an unmatched rho*, not for instructions; that favours
  neither arm in a way this note can quantify, and is disclosed rather than
  re-tuned.
- **Batch sweep (B).** n = 53; (L, K) ∈ {(1,024, 440), (4,096, 660), (16,384, 880)},
  the K the PR #830 evidence chose per L. Corpora `n53-strong-sweep-L<L>-v1`,
  except L = 1,024, which is the ladder's cell.
- **Metric.** Whole-process retired instructions (`valgrind --tool=callgrind`,
  `Ir`) for both arms in every cell, plus native wall clock for reference only.
- **Gates.** Every target in every cell verified and equal to the planted scalar;
  otherwise the cell is reported as failed and excluded from fits (not repaired
  or re-seeded).

## Quantities and decision rules (fixed now)

- `ρ*(cell) = Ir(IC) / Ir(rho R3)`, `< 1` meaning IC is cheaper.
- **Slope fit (A).** Ordinary least squares of `log₂ Ir` on `log₂ r`, separately
  for each arm, over the gate-passing size cells, with the standard error of each
  slope and of their difference computed from the residuals (n − 2 degrees of
  freedom, Student-t 95 % interval). `r` is the subgroup order of the cell's
  curve.
- **The lead survives as a scaling result** iff *either* (i) `ρ*(cell) < 0.8` in
  at least one measured cell with the whole-process cost, *or* (ii) the fitted
  slope difference `s_IC − s_rho` is `≤ −0.05` **and** its 95 % interval excludes
  0. A crossover size implied by (ii) is an **extrapolation** and is reported
  with the slopes it rests on; it is not a measurement.
- **The lead dies as a scaling result** iff neither holds. In particular, four
  points with a slope-difference interval that includes 0 do not establish a
  gap in either direction, and are reported as "not resolved".
- **Batch sweep (B)** is descriptive: the trend of `ρ*` in L is reported, and
  survival under (i) also applies to B cells.

**Inadmissible, stated in advance:** re-tuning K after seeing Ir; dropping or
adding cells after seeing results; changing rho's rung or lanes; excluding a
cost from either arm; scoring on wall clock; fitting a subset of the size cells
chosen after the fact.

## Scope and honesty limits

- The reference is the best rho built in this repository, not the best rho
  possible: after R3 the canonicalization scan is still about half of its
  instructions, and the DP table still uses the default hasher. IC is as merged.
  Both are one implementation each, on one x86-64 host class with `pclmulqdq`.
- Sizes 37, 41, 53, 61 are prime-degree fields of one Koblitz family with
  different cofactors and subgroup structure; a slope over four such sizes is a
  fit over a small, irregular set, and the scatter around the line is reported.
- The `n` in the cost model `√(L·r/n)` and the automorphism count both change
  with the field size; nothing here corrects for that, because the metric is the
  measured cost, not a modelled one.
- Nothing here speaks to m = 83 or to ECC2K-130 (`AGENTS.md` §8a, §8b); any
  crossover extrapolation is a statement about a fit, and the m = 83 gate is
  separate.
- Class (`AGENTS.md` §3): **accounting** for the reference change; the sweep
  itself produces a measured ratio-to-boundary trend, and is classified after it
  is read.

---

## Results (2026-09-29, additive — nothing above this line was edited after seeing it)

Generated by `strong_rho_sweep_20260929_run/sweep_analysis.py` from the raw logs;
host, binaries and build as in the ladder note (x86-64 Xeon 2.1 GHz with
`pclmulqdq`, `rustc 1.94.1`, `valgrind 3.22.0`; `koblitz_rho_batch_ks_strong` sha256
`3a4eb7b4…`, `koblitz_orbit_dlp_fast` sha256 `8dea682a…`).

**Gates.** All six cells passed G1 for both arms: rho recovered and verified every
target and equal to the planted scalar; IC solved every target with zero failures
and full rank K. The rho records of every cell also pass the independent
pure-Python replay (`independent_replay_rho.json` per cell: 1,024 / 1,024 / 4,096 /
16,384 / 1,024, plus the ladder cell). Peak RSS: IC 0.01 GiB (n = 37) … 2.5 GiB
(n = 61) … 4.94 GiB (n = 53, L = 16,384); rho 11–591 MiB.

**Whole-process retired instructions.**

| n | L | K | log₂ r | IC Ir | rho R3 Ir | IC / rho | IC wall s | rho wall s |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 37 | 1,024 | 42 | 27.78 | 560,613,246 | 705,849,698 | **0.794** | 0.2 | 0.1 |
| 41 | 1,024 | 255 | 39.00 | 12,441,253,531 | 8,184,739,666 | 1.520 | 4.2 | 1.3 |
| 53 | 1,024 | 440 | 44.26 | 89,190,346,806 | 46,384,789,571 | 1.923 | (30.8, PR #955) | 5.3 |
| 61 | 1,024 | 600 | 47.21 | 292,412,736,909 | 127,304,678,714 | 2.297 | 89.3 | 15.6 |
| 53 | 4,096 | 660 | 44.26 | 153,201,588,995 | 95,528,133,576 | 1.604 | 52.6 | 11.9 |
| 53 | 16,384 | 880 | 44.26 | 298,529,760,485 | 199,039,001,579 | 1.500 | 111.4 | 24.6 |

**Size fits (L = 1,024, four sizes, as registered).** `log₂ Ir` on `log₂ r`: IC
slope 0.459 ± 0.029, rho R3 0.381 ± 0.031 (95 % intervals 0.335–0.583 and
0.249–0.513, df = 2; both far from the asymptotic ½ because the two smallest
cells are dominated by fixed costs). The registered quantity, the slope of
`log₂(IC/rho)`, is **+0.078 ± 0.003, 95 % interval 0.067–0.089**: IC's instruction
count grows *faster* than the strong rho's, not slower.

**The registered decision, applied literally.** Clause (i) — some measured cell
with `IC/rho < 0.8` — is met, **by one cell only: n = 37, L = 1,024, 0.794**.
Clause (ii) — slope gap `≤ −0.05` with an interval excluding 0 — is not met. By
the rule as written the outcome is therefore **"survives as a scaling result"**,
and that is what the analysis script prints and what is recorded.

**Why that label should not be read as the lead surviving, and what was wrong
with the rule.** I wrote clause (i) so that any single cell could satisfy it,
with no condition on its size or on the other cells. That was a design flaw I
should have caught before registering, and I am recording it as such rather than
correcting the rule after the fact. The evidence around the cell that triggers
it:

1. It is the smallest cell (a 2^27.8-order group) and both runs are tiny (0.56 B
   and 0.71 B instructions in total), so fixed costs — curve construction, jump
   table, 1,024 target generations and verifications — are a large share of each.
2. The margin is 0.794 against a threshold of 0.8, inside the roughly 10 %
   build-to-build variation the ladder measured for an identical algorithm.
3. At the same cell IC is *slower* in wall time (0.182 s against 0.113 s, 1.61×),
   although it retires fewer instructions.
4. **Every cell with n ≥ 41 is at or above 1.5**, and the ratio rises with n:
   0.794, 1.520, 1.923, 2.297 at n = 37, 41, 53, 61.
5. The unfavourable slope is not produced by the small cell. Dropping any one
   size gives a slope difference of +0.072, +0.078, +0.079 or +0.078, and the local
   slope between adjacent sizes is +0.083 (37→41), +0.064 (41→53) and +0.087
   (53→61) (post hoc, not registered).

**What the sweep does show.** (a) Against the strong reference, IC has no
instruction-count advantage at any size from n = 41 up, at any tested batch size.
(b) Its cost ratio to the reference grows with the field size (about `r^0.08`)
over the four sizes measured. If that slope persisted — an **extrapolation** from
four points with K as chosen at each size, not a measurement — the ratio at n = 61
(2.30) would grow by about 7× by r ≈ 2^83 and about 89× by r ≈ 2^130. It is in
any case not a reason to expect an eventual crossover.
(c) At n = 53 the ratio **falls with the batch size**, 1.923 → 1.604 → 1.500 at
L = 1,024 / 4,096 / 16,384, by −0.131 then −0.048 per doubling of L (post hoc): a
real but flattening trend, still 1.5 at 16× the targets, with IC then holding
4.94 GiB against rho's 0.58 GiB and taking 4.5× rho's wall time. Three points do
not support extrapolating it to a crossover.

**Exploratory operation counts (unregistered; computed after the native runs were
seen, and not entering the decision above).** The runs record IC's probes and
rho's canonicalizations, which are the dominant primitive of both arms (see the
ladder note's profile). IC / rho probes at L = 1,024 is 0.242, 0.763, 1.418, 1.932
at n = 37, 41, 53, 61 (slope of its log₂ on log₂ r: +0.155 ± 0.003), and at n = 53
it is 1.418, 0.935, 0.946 for L = 1,024, 4,096, 16,384. So at n = 53 and L ≥ 4,096
IC needs slightly *fewer* probes than rho needs canonicalizations, yet retires
1.5–1.6× the instructions: the difference is the index build (about K² entries,
10.3 M at K = 440 and 41 M at K = 880) and the poorer locality of a multi-GiB
table. `strong_rho_sweep_20260929_run/sweep_analysis.txt` has the table.

**Scope and limits.** Four sizes of prime-degree fields of one Koblitz family with
different cofactors, K as chosen for IC's own wall-clock optimum in PR #830 and not
re-tuned for instructions, one x86-64 host class, one implementation of each arm, a
reference (R3) that is not the best rho possible, and IC as merged. The fitted
intervals have two degrees of freedom. The n = 83 gate and ECC2K-130 transfer
(`AGENTS.md` §8a, §8b) are untouched; the extrapolations above are statements
about a fit.

**Classification (`AGENTS.md` §3).** Not an advance for IC: its ratio to the
reference rose with n over the measured range. The reference change is accounting
(hardware-matched rung) and engineering (stronger rungs); the sweep itself is a
measured ratio-to-boundary trend and is reported as such.
