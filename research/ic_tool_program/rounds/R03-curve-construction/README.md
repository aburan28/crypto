# R03: the curve's construction at composite degrees

The `ic` tool programme's third round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md`, Track A, A9). It
tests one lever on the set-up's curve construction:
- `factorise_u64` ran a deterministic Miller–Rabin test before every
  trial division;
- the candidate tests primality only when a division changes what is
  left, and skips even divisors after 2.

Its class is **engineering**. The declaration is
[`PROTOCOL.md`](PROTOCOL.md), merged in #1119 before anything ran.

## Status: accepted, as baseline v1

Every declared condition holds:

| condition | declared | measured | met |
|:--|:--|:--|:--|
| the pin | every pinned output of the candidate equals v0's, on all 90 rows | 90 of 90 rows, names agree | yes |
| `icv1-f2m57-tm747311035-c1f545af`, suite | cold interval above 1.3 | 1.771 [1.694, 1.851] | yes |
| `icv1-f2m57-tm747311035-c1f545af`, holdouts | cold interval above 1.3 | 1.922 [1.756, 2.104] | yes |
| no regression | no size's interval wholly below its A/A band | none | yes |

**What follows from the acceptance:**
- **v1 is the newest baseline:** `30f6c153` on v0′ `c1a2e5f8`. Its
  change reaches `src/` with these results, and `baselines.json`
  carries it.
- **Its class is engineering.** The factors are identical, divisor for
  divisor, so nothing the pipeline decides changes, and the ratio to the
  floor is flat.
- **R05 runs on v1**, as R05's protocol fixed before any R05 run.

## What ran

All steps ran on one host, in R01's container session. Its manifest
matches R01's (`aa-source.json`), so the A/A bands are R01's. Each timed
process ran under `tools/isolated_bench.py run --wait --cpus 2`, with
`RAYON_NUM_THREADS=1`, from the run tree `f0ee5da6`, whose R03 scripts
and harness are byte for byte main's.

| step | processes | contended | failed |
|:--|--:|--:|--:|
| pin | 90, untimed | — | 0 |
| comparison, five rounds | 382 | 42, each rerun | 0 |
| holdouts, five rounds | 42 | 2, each rerun | 0 |

Every row and round has a clean pair: the analysis found no missing
pair at any size.

## Refused attempts: none was timed

**The isolation tool refused 105 attempts before they ran.** A refused
attempt runs no benchmark and writes no figure; the harness waited and
tried the same row again until it ran. Every attempt is in
`refusals.log` beside its step's runs, and every wait for pressure to
settle in `psi-waits.log`.

| cause | attempts | rows | when (UTC) |
|:--|--:|--:|:--|
| the benchmark CPU outside the process's mask (below) | 44 | 5 | 15:54–16:00 and 17:14–17:22 |
| the machine busy: other processes above 0.20 CPU s in 2 s, or PSI above 5 | 61 | 61 | 16:00–18:37 |

**The busy refusals are the tool working as designed.** In every one,
this session's own agent process was among the busiest neighbours. The
host's environment manager was among them in 19, and a `git` or `pip` in
4.

### The CPU-mask incident

**The cause.** An isolated run moves every other movable thread off the
benchmark CPU while it runs, and restores them afterwards.
- A process forked during that time inherits the narrowed mask.
- When R03's runner started its next isolated child in that window, the
  child's mask lacked CPU 2.
- The isolation tool then refused the child, correctly: it pins the
  benchmark to a CPU the process may use.

**When it happened.**
- **15:54–16:00 UTC, 25 refusals** of one row
  (`base/icv1-f2m19-tm797-9c54981b/M1-T01/r1`). The runner itself had
  inherited the narrowed mask. The mask was widened by hand
  (`taskset -p f`) on the runner and on one other process, and the row
  then ran.
- **17:14–17:22 UTC, 19 refusals** of four rows at
  `icv1-f2m57-tm747311035-c1f545af`, round 3. An exploration was
  interleaving its own isolated runs under the same lock. Each refusal
  cleared once the other run restored the runner's mask.

**What it did and did not change.**
- Every figure in this round comes from a run that started with the
  benchmark CPU in its mask, and the contended-run rule applied as usual.
- The refusals cost wall time only.

**The fix.** `tools/isolated_bench.py` now widens an inherited mask to
include the benchmark CPU and records the mask it found as
`affinity_widened_from` (#1162). R03 ran with the tool as it was before
that fix, from its own run tree's commit. R05 and later rounds run with
the fix.

## The comparison: v0′ against the candidate

The ratio is v0′ over the candidate, so above 1 means the candidate is
faster. The two composite sizes have eight rows each (seeds 201–204 ×
two targets); the nine prime sizes have `M1`'s two rows. Every row ran
five rounds. The intervals are 95% `t` intervals of the geometric mean
of the paired ratios. The set-up and the construction stage are stage
diagnostics, not speedups. The construction is the pricer's `setup`
phase, which builds the curve; its milliseconds are medians. `S` is the
median cold cost in v0's unit (R01's `unit_ns`), each arm from this
round's own runs.

| curve | log₂ r | rows | cold [95%] | set-up | construction stage | construction, ms | `S` cold | R01's A/A, cold |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 2 | 0.995 [0.960, 1.032] | 0.994 [0.959, 1.030] | 1.020 [0.961, 1.082] | 0.24 → 0.23 | 132.01 → 132.33 | 1.020 [0.930, 1.118] |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 2 | 0.976 [0.938, 1.016] | 0.976 [0.937, 1.017] | 0.985 [0.931, 1.041] | 0.22 → 0.22 | 34.26 → 35.13 | 1.021 [0.984, 1.059] |
| `icv1-f2m45-tm6236725-40939294`, composite | 24.8 | 8 | 1.007 [0.985, 1.030] | 1.007 [0.985, 1.029] | 1.047 [1.017, 1.077] | 0.82 → 0.77 | 22.10 → 22.09 | 0.985 [0.916, 1.060] |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 2 | 1.049 [1.014, 1.084] | 1.050 [1.014, 1.087] | 1.082 [1.037, 1.130] | 0.67 → 0.63 | 10.83 → 10.53 | 0.976 [0.903, 1.055] |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 2 | 1.009 [0.985, 1.034] | 1.009 [0.986, 1.034] | 1.315 [1.280, 1.352] | 1.20 → 0.91 | 6.73 → 6.69 | 1.028 [0.936, 1.129] |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 2 | 0.971 [0.873, 1.079] | 0.971 [0.874, 1.079] | 1.299 [1.212, 1.393] | 1.08 → 0.79 | 4.45 → 4.39 | 1.019 [0.922, 1.127] |
| `icv1-f2m57-tm747311035-c1f545af`, composite | 38.0 | 8 | **1.771 [1.694, 1.851]** | 1.779 [1.702, 1.860] | 42.886 [41.025, 44.831] | **56.00 → 1.28** | 8.00 → 4.45 | 1.012 [0.935, 1.095] |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 2 | 1.034 [0.944, 1.133] | 1.035 [0.945, 1.133] | 1.101 [1.014, 1.195] | 0.55 → 0.52 | 5.90 → 5.83 | 1.075 [1.005, 1.149] |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 2 | 0.970 [0.914, 1.029] | 0.970 [0.914, 1.029] | 0.991 [0.909, 1.080] | 1.40 → 1.34 | 8.04 → 8.16 | 1.017 [0.951, 1.088] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 2 | 0.989 [0.950, 1.031] | 0.990 [0.950, 1.031] | 3.266 [3.039, 3.510] | 7.16 → 2.15 | 9.17 → 9.43 | 0.990 [0.934, 1.049] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 2 | 1.010 [0.987, 1.032] | 1.010 [0.988, 1.033] | 1.930 [1.696, 2.195] | 3.60 → 1.77 | 15.82 → 15.61 | 0.986 [0.926, 1.049] |

**What the table says.**
- **At `icv1-f2m57-tm747311035-c1f545af` the lever is the cost.** The
  construction falls from 56.0 ms to 1.28 ms, 42.9 times less, and the
  cold time is 1.771 [1.694, 1.851] times faster. `S` falls from 8.00
  to 4.45.
- **At `icv1-f2m45-tm6236725-40939294` it is not.** The construction
  moves from 0.82 ms to 0.77 ms, and the cold time reads 1.007
  [0.985, 1.030].
- **The prime sizes read 0.970–1.049**, and every interval overlaps its
  size's A/A band. Their construction stage is faster too, by up to 3.3 times at
  `icv1-f2m59-tm943548413-98844ecc`, but it is a few milliseconds of a
  set-up of seconds. The one interval that excludes 1,
  `icv1-f2m37-tm534059-32aad96b` at 1.049 [1.014, 1.084], overlaps its
  A/A band, [0.903, 1.055].

## Against the prediction

| curve | predicted | suite | holdouts |
|:--|--:|--:|--:|
| `icv1-f2m57-tm747311035-c1f545af` | 1.75× [1.5, 2.0] | 1.771 [1.694, 1.851] | 1.922 [1.756, 2.104] |
| `icv1-f2m45-tm6236725-40939294` | 1.3× [1.1, 1.5] | 1.007 [0.985, 1.030] | 1.001 [0.913, 1.098] |
| the nine prime sizes | 1.00× | 0.970–1.049 | — |

**`icv1-f2m57-tm747311035-c1f545af` landed on the prediction.** The
protocol predicted 1.75× from callgrind at that size, where `factorise`
was 71% of the instructions up to factor-base selection.

**`icv1-f2m45-tm6236725-40939294` missed it, and the order says why.**
- The prediction assumed the loop dominated the construction there as
  at `n = 57`. Callgrind had run at `n = 57` only.
- `#K_1(GF(2^45))` is `2 · 7 · 11 · 37 · 211 · 29264761`. The trial
  division stops at 211, so the old loop tested primality about two
  hundred times, and the candidate tests it six.
- `#K_0(GF(2^57))` is `4 · 130873 · 275295876199`. There the old loop
  ran to 130,873 and tested at every divisor; the candidate tests three
  times.
- So R01's 0.8 ms of construction at `n = 45` is mostly elsewhere, and
  this lever does not reach it.

## The holdouts

R02's frozen files for seed 205, targets `T101` and `T102`, at the two
composite sizes, five rounds.

| curve | pairs | cold [95%] | set-up | construction stage | construction, ms | `S` cold |
|:--|--:|--:|--:|--:|--:|--:|
| `icv1-f2m45-tm6236725-40939294` | 10 | 1.001 [0.913, 1.098] | 1.000 [0.914, 1.095] | 1.042 [0.938, 1.158] | 0.81 → 0.77 | 22.90 → 22.40 |
| `icv1-f2m57-tm747311035-c1f545af` | 10 | 1.922 [1.756, 2.104] | 1.930 [1.765, 2.111] | 47.915 [43.312, 53.007] | 56.91 → 1.21 | 8.30 → 4.29 |

The holdouts agree with the suite at both sizes: 1.922 against 1.771
at `n = 57`, and 1.001 against 1.007 at `n = 45`.

## Files

- `runs.tar.xz` (`runs.tar.xz.sha256`): every process's report and
  isolation record, the pin, the host manifest, and the refusal and PSI
  logs. `tar -xJf runs.tar.xz && python3 analyse.py > analysis.json`
  reproduces `analysis.json` byte for byte.
- `analysis.json`: `analyse.py`'s output from the archive.
- The candidate: `30f6c153` on v0′ `c1a2e5f8`, one function in
  `src/cryptanalysis/koblitz_index_calculus.rs` and its test. The binary
  measured has SHA-256 `23c83a47…` (`host.json`).
