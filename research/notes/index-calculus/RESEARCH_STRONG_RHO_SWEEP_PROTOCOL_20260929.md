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
