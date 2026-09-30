# Ledger §22 protocol, v1: the descent's fixed cost per target

Declared 2026-09-30, before any candidate code existed and before any
timed comparison.

The only measurements made first are a pin check and a probe, both on the
unmodified binary (`main` at `57e7ce3a`):

- **Pin check.** `ic price` at `n = 19` and `n = 41`, `M1`, reproduces
  §21's counts and recovered logarithms exactly, so `main` still runs
  §20's recipes as §21 measured them.
- **Probe.** `examples/koblitz_descent_prices.rs` priced each target's
  descent part by part at six sizes, on §20's `M1` files, with the log
  table from an `ic workflow` run of the same file. Its trials match the
  pricer's at every size. The outputs are in `probe/`, and they chose the
  design below.

## Where §21 left it

After §21 the descent is the largest phase at six of nine sizes, 37–45%
of `S` below `2^37`. It is paid per target, so no number of targets
amortises it (§20.5). Below `2^37` it is not probes that cost: at
`n = 19` a target takes eight trials and about 1,300–1,500 units.

The probe, in units per target (one batched affine addition, measured in
the same process):

| `M1` | trials | descent | walk start | of which 63 additions | relation assembly | recovery check |
|:--|--:|--:|--:|--:|--:|--:|
| `K_1/GF(2^19)` | 8.2 | 1,481 | 635 | 549 | 391 | 118 |
| `K_1/GF(2^23)` | 105.5 | 1,984 | 752 | 658 | 429 | 146 |
| `K_1/GF(2^45)` | 160.8 | 2,449 | 1,162 | 1,045 | 345 | 155 |
| `K_0/GF(2^37)` | 1,619 | 4,408 | 1,060 | 918 | 387 | 177 |
| `K_1/GF(2^43)` | 5,115 | 9,530 | 1,120 | 970 | 371 | 190 |
| `K_0/GF(2^41)` | 24,134 | 38,789 | 1,156 | 986 | 488 | 210 |

- **The walk start** is three scalar multiplications (21–29 units each)
  and 63 additions made one at a time, each paying its own inversion. The
  same 63 points, as one batched addition of the start to precomputed
  multiples of the stride, cost 78–90 units.
- **The relation assembly** turns the one decomposition found into `d`,
  in big-integer arithmetic: a `modpow` per summand, a dot product with
  the column logarithms and a modular inverse.
- **The recovery check** is `[d]G = Q` in the general arithmetic. The
  single-word arithmetic makes the same check in 61–117 units.
- **What is left** is stepping and lookups: about 1.5 units a probe, one
  batched addition and one table lookup. That is near its floor, and this
  round does not touch it.

## The change (candidate)

Every walk, every trial and every recovered logarithm stays as it is.

1. **The walk start, batched.**
   - The solver precomputes the 63 multiples `[j·stride]G`, `j = 1..63`,
     once. The stride depends only on `r`.
   - Each target then gets its 64 starts from one batched addition of its
     first start to those multiples, instead of 63 additions one at a
     time.
   - The first start is `[a₀]G + [b]Q` as now, from the same random
     draws, so every walk starts where it did.
2. **The relation's arithmetic in single words.**
   - The fast path already requires `r < 2⁶²`. There, the solver keeps
     the column logarithms, `λᵏ mod r` and `h mod r` as machine words.
   - It computes `d = (Σ ±λᵏ·x_o − h·a)·(h·b)⁻¹ mod r` in `u128`
     arithmetic, with the inverse by the extended Euclidean algorithm.
   - `d` is unique modulo `r`, so it is the same `d`.
3. **The solver's recovery check in single words.**
   - `[d]G = Q` is checked with the fast curve's scalar multiplication,
     which is tested equal to the general one.
   - The pipeline's own final verification (`verify_final`) stays in the
     general arithmetic, as an independent check.

**Tests** check that:
- the batched starts equal the sequential ones, point for point;
- the single-word `d` equals the big-integer `d` on every relation of a
  full descent at several degrees;
- a descent recovers the same logarithm, with the same trials, as before.

## Frozen inputs

§20's 36 measurement parameter files, the same as §21's: nine sizes,
sets `M1`–`M4`. They are listed with their sha256 in `inputs.sha256`.

## Baseline and candidate

- **Baseline:** `ic` built at `57e7ce3a` (`main`), rustc 1.94.1. Kept
  outside the tree, sha256
  `07e527a4d9339e54020db55b5477a95f472b6abb6e6969ed1c003d48e7ee2612`.
  The probe binary (the same library) has sha256 `75e59a86…`.
- **Candidate:** `ic` built at the candidate commit, from a clean tree.
  Its sha256 is recorded when it is built.

## Procedure

**The main comparison**, as in §21. For every size and every set, five
rounds of baseline then candidate: 360 processes.

- Each process is `ic price --repeats 3 --repeats-fast 15`, with
  `RAYON_NUM_THREADS=1` under `taskset -c 2`. Nothing else runs.

**Controls:**

- **Pin:** every candidate report's `counts` (trials per target
  included) and `recovered` must equal the baseline's in all 180 pairs.
- **Control 1:** `ic workflow` and `ic price` agree field by field on the
  candidate, on `M1` at every size.
- **Rho:** batch rho re-priced on the baseline at `n = 41`, `M1`, with
  §20's seeds. Its counts must equal §20's, and its price is reported
  beside §20's.

**The probe, after.** `koblitz_descent_prices` on the candidate's library
at the same six sizes, with the same log tables.

**Many threads.** At `n = 41`, `M1`, three ABAB rounds with
`RAYON_NUM_THREADS=4` under `taskset -c 0-3`.

## Measures

§21.4 found that the unit itself runs up to 11% faster in one binary than
in another, with its code unchanged. A ratio of totals each converted by
its own process's unit therefore mixes two calibrations. So this round
converts both arms at one unit.

- **Speedup per pair (primary):** baseline total nanoseconds over
  candidate total nanoseconds. This is `baseline_total_operations /
  candidate_total_operations` in one calibrated unit (AGENTS.md §8),
  whichever arm's unit it is.
  - Per size: the geometric mean over 20 pairs, with a 95% interval (`t`,
    19 degrees of freedom, on the logarithms).
- **Speedup in each process's own unit:** reported beside it, as §21's
  declared measure was.
- **`S`**, before and after:
  - after = `S` before (the baseline arm) over the primary speedup, that
    is at the baseline arm's unit;
  - `S` in the candidate's own unit is reported too;
  - the ratio to §20's batch rho per size, as in §21.
- **Stage diagnostics:** the descent phase per target, as a paired
  ratio; the probe's parts before and after.
- **A/A spread:** max over min of the five baseline totals, per set.

## Targets

1. **Identical outputs** in all 180 pairs, and Control 1 on the candidate
   at every size.
2. **The relation and its check fall at least 3× at every probed size.**
   That is the probe's relation assembly plus recovery check, the
   session's `target_descent + recovery_check`.
3. **The descent per target falls at least 1.8×** (the probe's `solve`)
   at `n = 19`, `23` and `45`.
4. **A whole-pipeline speedup (primary) whose 95% interval excludes 1**
   at the four sizes where the probe put the fixed part at 25% or more of
   the descent. The fixed part is the walk start, the relation assembly
   and the recovery check. The four sizes are `K_1/GF(2^19)`,
   `K_1/GF(2^23)`, `K_1/GF(2^45)` and `K_0/GF(2^37)`.
5. **No regression.** No size's primary interval lies wholly below 1, at
   one thread or at four.

## Stop

- **Any mismatch in counts or recovered logarithms:** stop, the candidate
  is wrong. Fix it, and restart the comparison from nothing.
- **An A/A spread above 1.25 on a set:** that set's five rounds are
  rerun once with double the repetitions, and both are reported.

## Suite and class

**Suite.** The AGENTS.md §8 frozen WDSat suite does not apply: no
decomposition-solver path changes. The matched suite is §20's frozen
inputs, priced as §20 and §21 priced them.

**Class: engineering** (AGENTS.md §3).

- The counts, the relations, the counting floor and the ratio to that
  floor do not move.
- `S` falls where the descent's fixed part weighs.
- It is not an advance: nothing the method finds changes.
