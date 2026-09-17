# Koblitz index calculus vs rho, 2026-09-12: no crossover at n=13, 37, 41, 53

This directory is the committed evidence behind retracting the
`regimes.koblitz.records.vs_rho` record in `docs/ic/boundary_targets.json`. The
retracted row claimed a charged-time crossover at n=37; it pointed at a run
directory that is absent from the tree and appears in no commit, so it was never
independently replayable. Everything quoted in the ledger now resolves to a file
here.

Public synthetic Koblitz fixtures only. Not key recovery, not asymptotic
sub-rho, not imported points. A negative result for this implementation at these
rungs, not a proof that no crossover exists.

## What is here

| Path | Contents |
| --- | --- |
| `target_mode_sweep.json` | The headline comparison: n=37, 1024 targets, all three target modes against rho. |
| `summaries/` | Per-arm producer batch summaries and aggregated per-target cost totals. |
| `sweep_target_modes.py` | The script that produced both. Re-runnable. |
| `autolab_runs/` | Five `boundary_autolab.py` run bundles, one per rung quoted in the ledger. |
| `rung_breakdown.json` | Cost components for those five rungs, reduced out of the gitignored producer logs. |
| `extract_rung_breakdown.py` | The reduction. Runnable only while the source logs are still on disk. |

## Headline: n=37, 1024 targets, charged

`rho_over_ic` above 1.0 would be an index-calculus win.

| Arm | Charged ms/target | rho_over_ic |
| --- | --- | --- |
| rho (`signed_frobenius`, packed), 1024/1024 verified | 12.822 | 1.000 |
| direct, `partition_walk` | 14.778 | 0.868 |
| direct, `coefficient_walk` | 15.330 | 0.836 |
| direct, `independent` | 20.270 | 0.633 |

Index calculus is 1.15x behind rho in its best target mode. The registered
autolab beat `koblitz.vs_rho.n37_wall` runs `independent`, the worst of the
three, which costs 37% more per target than `partition_walk`. Any ledger row
sourced from that beat alone therefore reads worse than the method's best
measured mode.

## How charged cost is defined

The direct arm's charged cost is `full_algorithm_charged_total_ms` from the
producer's own batch summary row: curve setup plus support-table setup charged
**once** for the batch, plus per-target fixture setup and relation collection.

Do not derive it by summing the per-target `charged_total_ms` field. That field
re-adds the shared `setup_ms` for every target. At n=37 that inflates the batch
by about 14 ms per target and turns 20.27 into 34.08. The retracted ledger row's
successor draft made exactly this mistake.

The rho arm's charged cost is the sum of per-target `total_ms`, which is
`setup_ms + walk_ms + validation_ms`.

## Two asymmetries, both recorded, neither corrected for

Ratios above charge each arm its full cost as its producer defines it. Nothing
is subtracted. Two things are nevertheless worth knowing when reading them, and
they point in opposite directions:

- **rho is charged for building its own instance.** Its `setup_ms` (3.976
  ms/target of the 12.822) covers computing the target `Q = d0*G` and the jump
  table. Constructing `Q` is instance generation; the direct arm reports the
  equivalent separately as `fixture_generation_ms` and does not charge it. This
  runs against rho.
- **The direct arm's largest component is an assertion.**
  `solution_validation_ms` is 7.780 of `partition_walk`'s 14.778 ms/target, or
  53%. It re-derives the discrete log of every factor-base representative with a
  scalar multiplication and replays every collected relation against the solved
  vector, under `assert_eq!`. It confirms an answer the linear solve has already
  produced. This runs against index calculus.

Setup amortization is not a lever in either direction: curve plus support-table
setup is under 0.02 ms/target at 1024 targets. Whatever closes the remaining gap
has to come out of per-target relation collection.

## The other rungs

From `autolab_runs/`, one bundle each. Ratios are `rho_cost / ic_cost` in the
bundle's own timing class.

| Rung | Bundle | Timing class | ratio |
| --- | --- | --- | --- |
| n=13 | `20260912T161059Z-07471fc673` | whole_process_wall | 0.437 |
| n=37, 1 target | `20260912T161100Z-5f1d599743` | whole_process_wall | 0.146 |
| n=41 | `20260912T161101Z-d038475a6b` | algorithmic_charged | 0.032 |
| n=53 | `20260912T161208Z-d6cb39c596` | whole_process_wall | 0.006 |
| n=37, 1024 targets | `20260912T163230Z-3f7e563ca0` | whole_process_wall | 0.155 |

n=53 constructs and both arms recover, but index calculus is 180x rho there. It
measures construction feasibility, not competitiveness.

All four of those are single-fixture runs, so the direct arm carries its whole
support-table build with nothing to amortize it against. `rung_breakdown.json`
separates the components, and the separation is what rules n=41 out as a next
target: discount the 6877 ms support build entirely and per-target relation
collection is still 821.7 ms against rho's 247.5 ms. No number of targets
sharing one base closes a 3.3x gap in the per-target term.

## Reproducing

The two producers must hash to the values in `target_mode_sweep.json`
(`executable_sha256`):

```bash
cargo build --release --example koblitz_rank_fixture --example koblitz_rho_fixture
python3 research/sat_factor_base_review_20260908/autolab/evidence/20260912-koblitz-vs-rho-no-crossover/sweep_target_modes.py
```

Roughly four and a half minutes: three direct arms at about 80 s each plus a 13 s
rho arm. `--fixtures 4` gives a fast smoke run, though at 4 targets the batch
setup is not amortized and the charged figures are correspondingly high. Pass
`--from-summaries` to rebuild `target_mode_sweep.json` from `summaries/` without
re-running the producers.

The `independent` arm reuses the seeds of autolab run
`20260912T163230Z-3f7e563ca0`, so it doubles as a cross-check against that
bundle: 20.270 ms/target here against 20.094 ms/target there, and rho 12.822
against 12.801.

## Why the producer logs are not here

`boundary_autolab.py` keeps each run's raw producer stdout under
`autolab/runs/<id>/logs/`, which `autolab/.gitignore` excludes. Those files are
56 MB for these five runs, and the per-relation records carry factor-base point
coordinates, target point keys and walk coefficients. The bundles promoted here
are the aggregate layers only: `artifacts/`, `inputs/`, `receipts/` and
`state.json`.
