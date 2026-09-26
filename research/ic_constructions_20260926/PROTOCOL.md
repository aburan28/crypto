# Ledger §21 protocol, v1: the workflow's constructions, built once

Declared 2026-09-26, before any candidate code existed and before any
timed comparison.

The only measurements made first are two runs of
`examples/koblitz_construction_prices.rs` on the unmodified binary. They
priced each constructor alone, at `n = 41` and `n = 61`, on §20's `M1`
parameter files, and they are the probe that chose the design (below).

## Where §20 left it

Ledger §20.5 found that the constructions the Koblitz workflow builds
around its work cost:

- 7–70% of `S` with the curve's setup;
- 107–433 units per base point without it.

Three constructions are rebuilt in every run. The probe priced each one
at `n = 41` (`|F|` = 6,560) and `n = 61` (39,040).

| construction | builds per run | units per base point, `n = 41` | `n = 61` |
|:--|--:|--:|--:|
| **projected signed-orbit map** | 5: selection's column count, each coverage (two when extension units run), the log solver, the descent solver | 23.3 | 59.1 |
| point index (`index_map`, keyed by big integers) | 3: each collector (two with extension units), the descent solver | 5.7 | 5.6 |
| decomposition check (`m_can_decompose`, `[r]P` per signed orbit in big-integer arithmetic) | 2: each collector | 10.8 | 8.5 |

The map dominates, and its single build is dominated by one
cofactor multiplication per base point. The base's points come in
Frobenius orbits and `[h]` commutes with Frobenius, so one multiplication
per orbit is enough.

The curve's own construction (`setup`: 72.7 and 161 units per point in
the probe) factorises the group order. It is out of scope here and stays
in `S` as declared in §20.

## The change (candidate)

1. **Build each construction once per base.** `FrobeniusFactorBase`
   keeps the projected map, the point index and the cofactor classes of
   the decomposition check in `OnceLock`s. They are built on first use
   and shared by every consumer of that base.
   - The map and the classes are keyed by the curve's identity (degree,
     subfield, `a`, `b`, cofactor) and the base's size.
   - A call with a different curve computes afresh and is not cached.
   - The consumers borrow instead of rebuilding: selection's column
     count, `ColumnCoverage`, `RelationCollector`, the log solver and
     `IndividualLogSolver`. No public signature changes.
2. **One cofactor multiplication per Frobenius orbit.**
   - The fast map lifts every point once.
   - It walks each recorded orbit `i₀, π(i₀), π²(i₀), …` by Frobenius to
     confirm the order, and projects only `i₀`.
   - Every other point's projection follows by Frobenius.
   - An orbit whose recorded order does not match falls back to one
     multiplication per point.

**The same answers.**

- Every map, index and class set is the one the old code built.
- Every count, every relation and every recovered logarithm is
  unchanged.
- Tests compare the new map with the per-point construction on several
  bases, and compare the cached constructions with fresh ones.

## Frozen inputs

§20's 36 measurement parameter files: nine sizes, sets `M1`–`M4`. They
are listed with their sha256 in `inputs.sha256`.

## Baseline and candidate

- **Baseline:** `ic` built at `e7022b75` (`main`, the merge of #824),
  rustc 1.94.1. Kept outside the tree, sha256
  `e51e9431777a33572dc56ec4d4cca6ccf8d6d8419dbbe1e329c0446adef9e612`.
- **Candidate:** `ic` built at the candidate commit, from a clean tree.
  Its sha256 is recorded when it is built.

## Procedure

**The main comparison.** For every size and every set, five rounds of
baseline then candidate: ten processes, ABAB….

- Each process is `ic price --repeats 3 --repeats-fast 15`, with
  `RAYON_NUM_THREADS=1` under `taskset -c 2`. Nothing else runs.
- No reference runs in these. The batch-rho code does not change, so
  each set's §20 priced rho stands.
- One control run re-prices batch rho on the candidate at `n = 41`
  `M1`, and its counts must equal §20's.

**Controls:**

- **Pin:** every candidate report's `counts` and `recovered` fields must
  equal the baseline's for the same parameter file, in all 180 pairs.
- **Control 1:** on the candidate, `ic workflow` and `ic price` must agree
  field by field (`research/ic_exponent_20260926/control.py`) on `M1` at
  every size.

**Many threads.** At `n = 53` and `n = 61`, `M1`, three ABAB rounds with
`RAYON_NUM_THREADS=4` under `taskset -c 0-3`.

**Construction prices.** `koblitz_construction_prices` on both binaries at
`n = 41` and `n = 61`, before and after. This is the per-constructor
diagnostic.

## Measures

- **Speedup per pair:** baseline total over candidate total, in the
  unit each process measured for itself. This is AGENTS.md §8's
  `baseline_total_operations / candidate_total_operations`.
- **Per size:** the geometric mean over its 20 pairs, with a 95% interval
  (`t`, 19 degrees of freedom, on the logarithms).
- **The constructions group** (`select_projection`, `collect_setup`,
  `logs_setup`, `descent_setup`) gets the same treatment, as a stage
  diagnostic.
- **A/A spread:** max over min of the five baseline totals, per set.
- **Reference:** `S` and its ratio to §20's batch rho per size, before
  and after.

## Targets

1. **Identical outputs** in all 180 pairs, and Control 1 on the candidate
   at every size.
2. **The constructions group falls at least 3×** at every size.
3. **A whole-pipeline speedup** whose 95% interval excludes 1, at every
   size where §20 measured the constructions at 15% or more of `S`. That
   is seven sizes: all but `K_0/GF(2^53)` and `K_0/GF(2^61)`.
4. **No regression.** No size's speedup interval lies wholly below 1, at
   one thread or at four.

## Stop

- **Any mismatch in counts or recovered logarithms:** stop, the candidate
  is wrong. Fix it, and restart the comparison from nothing.
- **An A/A spread above 1.25 on a set:** that set's five rounds are
  rerun once with double the repetitions, and both are reported.

## Suite and class

**Suite.** The AGENTS.md §8 frozen WDSat suite does not apply: no
decomposition-solver path changes. The matched suite is §20's frozen
inputs, priced as §20 priced them.

**Class: engineering** (AGENTS.md §3).

- The counts, the relations and the counting floor do not move, and
  neither does the ratio to that floor.
- `S` falls, so the ratio to the generic floor and to batch rho falls
  with it.
- It is not an advance: nothing the method finds changes.
