# Ledger §23 protocol, v1: the top-end exponent, in range

> **Withdrawn before anything ran.** AGENTS.md's IC measurement rules
> (#1076) made one previously unseen target the primary comparison, and
> bar a ratio to batch rho at `k = 32` as a headline. This declaration's
> primary figure was that ratio. Nothing below ran. §23 was declared
> again as a single-target round in
> [`../ic_single_target_20260930/`](../ic_single_target_20260930/PROTOCOL.md).
> This file is kept unchanged below as the record.

Declared 2026-09-30, before anything below ran. The only computations
made first are `predict.py`, which carries §20's frozen model to the six
sizes below (`prediction.json`), and a pin check. The pin check ran the
round's binary on §20's `K_0/GF(2^41)` `M1` file, and its counts and
recovered logarithms equal §22's.

## Question

On current `main`, how does the Koblitz collection thread's ratio to
batch rho grow across the top of the range the library can build? And
does the growth resolve from zero?

§22 left the declared fit over the four largest sizes at `r^0.137`
[−0.059, 0.332]. That interval includes zero, and it separates neither
the law's `1/6` nor §20's model (0.141 there) from it. §22.12 planned new
rungs at `2^49`–`2^56`, but they need `n > MAX_N = 63`, which the library
refuses. That correction is §22.12's, merged as #1077.

This round works within that limit. It adds the only two other curves
with `n ≤ 63` and `r ≥ 2^36`, and doubles the sets per size. It can
narrow the fit. It cannot lengthen the fit's lever arm past `2^47.2`.

## Boundaries (unchanged from §20)

- **Floor**, per target at `k = 32`: `L(32)·√(π/4n)`, with
  `L(32) = 0.19869`.
- **Reference:** batch rho at `k = 32` on the index-calculus run's own 32
  targets, in the same process. Its counted operations are priced at the
  canonical step, measured against the same unit (§20).

## Sizes

Every Koblitz curve with `n ≤ 63` and `r ≥ 2^36`:

| curve | `log₂ r` | `r` | in §20 | proper intermediate subfields over `GF(2)` |
|:--|--:|--:|:--|:--|
| `K_1/GF(2^47)` | 36.6 | 106,781,081,677 | yes | none |
| `K_0/GF(2^57)` | 38.0 | 275,295,876,199 | **new** | `GF(2^3)`, `GF(2^19)` |
| `K_0/GF(2^41)` | 39.0 | 549,756,390,943 | yes | none |
| `K_0/GF(2^53)` | 44.3 | 21,044,858,204,113 | yes | none |
| `K_1/GF(2^59)` | 44.5 | 25,179,555,920,633 | **new** | none |
| `K_0/GF(2^61)` | 47.2 | 162,888,033,982,417 | yes | none |

- **Where `r` comes from.** It is the subgroup order that
  `KoblitzCurve::new(a, n)` builds, read by
  `examples/koblitz_curve_records.rs` (`curve_records.json`).
- **Curve families.** The `K_0` rows are `E_0`, the ECC2K-130 family.
  The `K_1` rows are `E_1`, a different Koblitz model (AGENTS.md §8b).
- **Subfields.** `GF(2^57)`'s subfields are disclosed. The method's base
  is Frobenius orbits of subgroup points over `GF(2)`, and it uses no
  subfield.

## Recipe

§20's rules, unchanged: `research/ic_exponent_20260926/make_params.py`.
- `summands = 3`, pair table, tier `auto`, collection aimed at the
  least-mentioned columns.
- The window, unit and descent-cap rules.
- 32 targets.

**The free parameter is re-swept at every size, on this round's binary,**
by §20's procedure:
- the column count on `prediction.json`'s grid (the model's optimum times
  `2^{j/2}`, `j = −2 … 4`, read in actual columns of `8⌈c/8⌉`), and the
  descent's summands, 2 or 3;
- edge extension by §20's Amendment 1 rule, at most three steps;
- on sweep set `W` (seed 101, targets 10100–10131);
- the least `S` per target is chosen before any measurement set runs.

Why re-sweep:
- The fit must compare every size's best recipe on one binary.
- §20's choices were made on an older binary. `main`'s n59 stack has
  since made the top end's build and collection 12–34% cheaper.
- A recipe left off its optimum at the large sizes would bias the
  exponent upward.
- Each of the four old sizes' choices is reported beside §20's.

## Sets

- **`M1`–`M8`:** seeds 201–208, targets `100·seed` to `100·seed + 31`.
  `M1`–`M4` are §20's sets. `M5`–`M8` are new.
- **Rho:** batch rho on each set's targets, seed `0x200000 + seed`.
- **Cold rho:** `M1` only, seeds `0x210000 + i`, as in §20.

## Accounting

§20's, unchanged:
- one batched addition as the unit, measured around each repetition;
- every phase timed and divided by it;
- three repetitions per set (fifteen under 50 ms), the median as the
  figure;
- counts identical across repetitions.

`S` per target is total units over `k·√r`.

**The ratio at a size** is the mean over `M1`–`M8` of the
index-calculus `S`, over the mean of the priced rho `S`. Its 95% interval
comes from the eight per-set ratios (`t`, 7 degrees of freedom).

## Isolation (AGENTS.md §10)

- **Every timed process runs through the isolation tool.** That is every
  `ic price`, sweep and sets alike, run as
  `tools/isolated_bench.py run --wait --cpus 2`, with
  `RAYON_NUM_THREADS=1`.
- **The tool's defaults stand:** a 2 s settle, other processes at most
  0.10 CPUs, PSI `some avg10` at most 5.
- **Waiting for pressure to fall.** Before each call the harness waits
  until PSI `some avg10` (CPU and memory) is below 4.0, polling each
  second for at most 120 s. The tool's own gate stays the authority.
  - A refused start is logged to `refusals.log` and tried again after
    15 s.
  - This wait answers §22.11's cost, where most starts were refused once
    on the pressure the benchmark before had left.
- **Contention.** A run the tool marks contended, or one that exits
  non-zero, is kept, excluded, and run again up to twice. The first clean
  run is the figure. Contended and clean runs are never pooled.
- **Control 1** (`ic workflow` against `ic price`, counts only) is not
  timed. It runs on `M1` at every size, after that size's sets, and
  nothing timed runs beside it.
- **The A/A.** The eight sets' repetitions are the A/A spread, as in §20.
  No second binary is compared, so no between-binary effect enters.

## Predictions (from `prediction.json`)

§20's model and law, carried unchanged:

| curve | `log₂ r` | model ratio | law ratio | model's `(c, m)` |
|:--|--:|--:|--:|:--|
| `K_1/GF(2^47)` | 36.6 | 3.56 | 4.53 | 55, 2 |
| `K_0/GF(2^57)` | 38.0 | 3.70 | 4.82 | 61, 2 |
| `K_0/GF(2^41)` | 39.0 | 4.51 | 6.38 | 104, 2 |
| `K_0/GF(2^53)` | 44.3 | 6.76 | 10.3 | 263, 2 |
| `K_1/GF(2^59)` | 44.5 | 6.61 | 10.1 | 251, 2 |
| `K_0/GF(2^61)` | 47.2 | 8.68 | 13.5 | 448, 2 |

- **The model's local exponent** (`ratio·√n` on `r`) is 0.139 over the
  six sizes, and 0.149 over the four largest.
- **The law's** is `1/6`.

## Targets

1. **Correct.**
   - Every target of every run is recovered and verified, index calculus
     and rho alike.
   - No relation is rejected.
   - Counts are identical across repetitions, and every declared run
     completes.
2. **Control 1** holds on `M1` at all six sizes.
3. **Exponent, primary.** Fit `ln(ratio·√n)` against `ln r` over the six
   sizes. Give the slope `β` with its 95% interval (`t`, 4 degrees of
   freedom).
   - The rise is *resolved from zero* if the interval excludes 0.
   - The law's `1/6` and the model's 0.139 are each *consistent* if the
     interval contains them, and *falsified at these sizes* if it does
     not.
   - **The round's aim is met if the interval excludes 0 or excludes
     `1/6`.** Otherwise the round reports that the sizes the library can
     build cannot separate them.
4. **Exponent, secondary, for continuity:** §20's four-largest fit
   (`2^39.0`–`2^47.2`), read the same way against the model's 0.149.
5. **Crossing.** A size whose ratio interval lies wholly below one is a
   crossing. None is expected.

## Inadmissible

- Choosing `(c, m)` on a measurement set.
- Tuning the window, unit or cap rules per size.
- Dropping a size, set, repetition, failed or contended run from the
  record.
- Pooling contended with clean runs.
- Changing `k`.
- Fitting per-set points as if they were independent sizes. That fit is
  reported only as a diagnostic.

## Stop and abandon

- **A verification failure:** stop that size and report it. It gives no
  figure until resolved.
- **Control 1 fails:** the pricer is wrong. Stop before any figure is
  quoted.
- **The repetition spread** above 1.25 on a set: that set is rerun once
  with double the repetitions, and both are reported. The rerun is the
  figure (§20).
- **The table budget refuses a grid point:** record it and continue.

## Cost, as an estimate

- 84 sweep pricings, 48 set pricings and 6 workflow runs.
- About 2–2½ hours on this host, isolated. The grid's upper points at
  `2^47.2` cost minutes each.

## Host

- **Binary:** `ic` at the declaration commit, built from a clean tree
  with rustc 1.94.1, kept outside the tree; its sha256 is recorded in
  `runs/host.json`.
- **Host class:** one x86-64 cloud container. Nothing is claimed for
  Arm64, GPUs or other hosts.

## Suite and class

**Suite.** The AGENTS.md §8 WDSat suite does not apply: no code changes.

**Class: accounting** (a measurement, AGENTS.md §3).
- No algorithm changes.
- The round re-measures the thread on current `main` and adds two sizes.
- It re-derives the exponent panel and the boundary facts, and rewrites
  the verdict if the reading changes.

**Scope, per AGENTS.md §8a.** There is no `m = 83` run: `n = 83` is past
`MAX_N`. High-fidelity ECC2K-130 improvement stays unestablished.
