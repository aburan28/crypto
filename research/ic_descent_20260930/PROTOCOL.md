# Ledger §22 protocol, v1: the descent's fixed cost per target

Declared 2026-09-30, before any candidate code existed and before any
timed comparison.  [v2](#v2-the-timed-runs-again-isolated) amends the
procedure after v1's runs, before any isolated run; everything v2 does
not name stands as v1 declared it.

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

## v2: the timed runs again, isolated

Declared 2026-09-30, after v1's runs and before any isolated run.

**Why.**
- AGENTS.md §10 makes `tools/isolated_bench.py` mandatory for every
  wall-clock or native-time number, and says `taskset` alone is not
  enough. The rule came with `aa677e4c`, which the baseline `57e7ce3a`
  already contains, so it was in force when v1 ran.
- v1 ran every process under `taskset -c 2` alone.
- Every figure this round quotes is timed:
  - the primary speedup is a ratio of times;
  - `S` is time converted by a unit that is itself timed;
  - the descent and the probe's parts are times.
- So v1's figures are not evidence under §10. They stay in `runs/` as the
  first run, labelled as not isolated. Nothing in them is deleted or
  pooled with v2.

**What does not change:**
- the change and the four binaries, by sha256;
- the frozen inputs;
- the pins;
- the measures and the analysis;
- the five targets and the stop rules;
- the class.

**What changes:**

1. **Isolation.** Every timed process runs as
   `tools/isolated_bench.py run --wait --cpus 2 -- …` with
   `RAYON_NUM_THREADS=1`, instead of under `taskset -c 2`.
   - The tool's defaults stand: a 2 s settle, other processes at most
     0.10 CPUs on average, PSI `some avg10` at most 5.
   - A refused start (a busy machine) is logged to `refusals.log` and
     tried again after 15 s.
   - Each process's isolation record is kept beside its report
     (`*.isolation.jsonl`).
2. **Contention.** A pair is excluded from every figure, but kept, if
   either process is marked `contended` by the tool or exits non-zero.
   - It is run again in the same slot, up to twice.
   - Each slot's figure is its first clean pair. A slot without one is
     reported as missing.
   - The contended processes are counted. Contended and clean pairs are
     never pooled.
3. **A/A.** Before the comparison, the baseline runs against a
   byte-identical copy of itself: `M1` at all nine sizes, five rounds,
   isolated the same way. Its paired ratio is the noise floor, reported
   beside every speedup.
4. **Many threads: three, not four.**
   - The tool will not reserve every CPU of this four-CPU host.
   - So the thread check runs at three threads on CPUs 1–3
     (`--cpus 1-3`, `RAYON_NUM_THREADS=3`), instead of four on 0–3.
   - Target 5's "at four" reads "at three".
5. **Controls.**
   - The pin check runs on every pair, retries and A/A included.
   - Control 1 compares counts only (`control.py`: base, collection,
     logs, trials, recovered logarithms), not times. v1's stands and is
     not rerun.
   - The rho control and the probe are timed, so they run again,
     isolated.
6. **Nothing else runs.** No build, test, lint or other job runs on the
   machine while the timed steps run.
7. **Order and place.**
   - The order is: manifest, inputs, A/A, main comparison, rho, threads,
     probe.
   - The outputs go to `runs-isolated/`, never over `runs/`. The analysis
     is `analysis-isolated.json`.

**Recorded beside it, not measured:** the nine curves' EC1 identities
(`curve_ids.json`, docs/curve-identities.md). They come from
`examples/koblitz_curve_records.rs`, which reads the defining
polynomial, subgroup and generator off `KoblitzCurve::new(a, n)`.

**Scope, per AGENTS.md §8a and §8b.**
- The `K_0` rows are `E_0: y² + xy = x³ + 1`, the ECC2K-130 family. The
  `K_1` rows are `E_1`, a different Koblitz model.
- The largest field here is `GF(2^61)`, and there is no `m = 83` run. So
  this round establishes nothing about ECC2K-130 at high fidelity.
- `GF(2^45)` has proper intermediate subfields over `GF(2)`: `GF(2^3)`,
  `GF(2^5)`, `GF(2^9)` and `GF(2^15)`. The change uses none of them.
