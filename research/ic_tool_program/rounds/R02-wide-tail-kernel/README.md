# R02: the AVX-512 batched addition for fields with a wide reduction tail

The `ic` tool programme's second round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md`, Track A, A2). It
tests one lever on the collection scan:
- the 8-lane AVX-512 batched addition refuses the two fields whose
  polynomial has `n + deg t = 66`;
- the candidate carries `H·t` past `z^63`, so the kernel serves them too.

Its class would be **engineering**. The declaration is
[`PROTOCOL.md`](PROTOCOL.md), merged in #1117 with amendments 1–3 before
anything ran.

## Status: rejected

**Two declared conditions are missed.**
- **The holdouts.** At `icv1-f2m59-tm943548413-98844ecc` they read
  1.161 [1.089, 1.239], and the lower end is not above 1.10. This is the
  success condition, so it decides the round.
- **The callgrind control.** The candidate's portable path is not
  instruction-identical to v0′'s: the difference is 0.1%, in the
  compiler's code for two scan functions (below).

Every declared condition, met or not:

| condition | declared | measured | met |
|:--|:--|:--|:--|
| the control: both AVX-512 scan kernels off, at `icv1-f2m53-tm56619371-dac20a85` (amendment 3) | collection per summand, scalar over kernel, not below 1.3, or R02 stops | 1.890 [1.738, 2.055] | yes: R02 goes on. The two kernels together clear the bar, which says nothing of the addition alone; the comparison tests that |
| the pin | every output of both arms equals v0's, and the names agree | 90 of 90 rows, names agree | yes |
| `icv1-f2m59-tm943548413-98844ecc`, suite | cold interval above 1.10 | 1.182 [1.151, 1.213] | yes |
| `icv1-f2m59-tm943548413-98844ecc`, holdouts | cold interval above 1.10 | 1.161 [1.089, 1.239] | **no** |
| `icv1-f2m61-t158598901-ab42b6c5`, suite | cold interval above 1.10 | 1.267 [1.242, 1.292] | yes |
| `icv1-f2m61-t158598901-ab42b6c5`, holdouts | cold interval above 1.10 | 1.262 [1.216, 1.310] | yes |
| no regression | no size's interval wholly below its A/A band | none | yes |
| callgrind control (amendment 3) | the portable paths instruction-identical (`analyse.py`: within `10⁻³`) | 0.085–0.102% fewer instructions in the candidate, all in two scan functions' codegen | **no** |

**What follows from the rejection:**
- **No new baseline.** v0 stays the newest, and `baselines.json` is
  unchanged.
- **The code is kept on record.** It is not merged into `src/`:
  - [`candidate.patch`](candidate.patch) is `a1645ac6` against
    `c1a2e5f8`, one file, `src/cryptanalysis/koblitz_fast.rs`;
  - the binary measured has SHA-256 `19240609…` (the run tree's
    `host.json`).
- **The scan has one failed round.** The plan sets a phase aside after
  two consecutive declared rounds on it fail (§11). The next scan round
  is chosen from R04's answer.

## What ran

All steps ran on one host, in R01's container session. Its manifest
matches R01's (`aa-source.json`), so the A/A bands are R01's. Each timed
process ran under `tools/isolated_bench.py run --wait --cpus 2`, with
`RAYON_NUM_THREADS=1`.

| step | processes | contended | failed |
|:--|--:|--:|--:|
| control | 20 | 0 | 0 |
| pin | 180, untimed | — | 0 |
| comparison, rounds 1–5, then rounds 6–10 at two sizes | 1,076 | 36, each rerun | 0 |
| holdouts | 232 | 13, each rerun | 0 |
| v0 against v0′ | 238 | 18, each rerun | 0 |
| callgrind | 6, untimed | — | 0 |

- **The comparison's extension.** Amendment 1 gives five more rounds to
  a size whose interval's half-width exceeds 5%. Two sizes got them:
  `icv1-f2m47-t22705043-f4e44623` and `icv1-f2m57-tm747311035-c1f545af`.
  Both read about 1 either way.
- **One holdout pair is missing.** `icv1-f2m19-tm797-9c54981b` `M5-T101`
  round 2 was contended on its first attempt and on both reruns, in the
  base arm. That size's holdout interval rests on nine pairs. Every
  attempt is in the run tree.

## The comparison: v0′ against the candidate

The ratio is v0′ over the candidate, so above 1 means the candidate is
faster. Each size has eight rows (seeds 201–204 × two targets) and five
rounds, or ten where marked. The intervals are 95% `t` intervals of the
geometric mean of the paired ratios. Collection and build are stage
diagnostics, not speedups. Units per summand and `S` are medians in v0's
unit (R01's `unit_ns`), each arm from this round's own runs.

| curve | `n + deg t` | cold [95%] | collection | build | online | units per scanned summand | `S` cold | R01's A/A, cold |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | — | 1.001 [0.977, 1.025] | 0.969 | 0.999 | 1.014 | 27.47 → 28.07 | 125.93 → 126.48 | 1.020 [0.930, 1.118] |
| `icv1-f2m23-tm5197-1f85e9e1` | — | 0.998 [0.962, 1.035] | 0.987 | 1.002 | 0.995 | 10.41 → 10.67 | 35.29 → 35.68 | 1.021 [0.984, 1.059] |
| `icv1-f2m45-tm6236725-40939294` | — | 1.027 [1.002, 1.054] | 1.017 | 1.039 | 1.037 | 8.28 → 7.82 | 22.61 → 22.02 | 0.985 [0.916, 1.060] |
| `icv1-f2m37-tm534059-32aad96b` | — | 1.024 [0.984, 1.066] | 1.016 | 1.017 | 1.011 | 4.89 → 4.86 | 9.96 → 10.01 | 0.976 [0.903, 1.055] |
| `icv1-f2m43-tm998717-e2e742b0` | — | 1.000 [0.969, 1.033] | 0.991 | 0.988 | 0.989 | 2.63 → 2.63 | 6.22 → 6.28 | 1.028 [0.936, 1.129] |
| `icv1-f2m47-t22705043-f4e44623` (ten rounds) | 52 | 0.980 [0.942, 1.019] | 0.976 | 0.970 | 1.026 | 1.50 → 1.51 | 4.61 → 4.75 | 1.019 [0.922, 1.127] |
| `icv1-f2m57-tm747311035-c1f545af` (ten rounds) | 61 | 0.973 [0.939, 1.008] | 0.959 | 0.972 | 0.958 | 1.58 → 1.62 | 8.08 → 8.40 | 1.012 [0.935, 1.095] |
| `icv1-f2m41-tm2308219-7f48b14a` | 44 | 1.006 [0.961, 1.052] | 1.013 | 0.987 | 0.993 | 1.45 → 1.44 | 5.97 → 6.12 | 1.075 [1.005, 1.149] |
| `icv1-f2m53-tm56619371-dac20a85` | 59 | 1.001 [0.969, 1.034] | 0.994 | 1.020 | 0.966 | 1.76 → 1.80 | 8.77 → 8.86 | 1.017 [0.951, 1.088] |
| `icv1-f2m59-tm943548413-98844ecc` | **66** | **1.182** [1.151, 1.213] | **1.288** [1.250, 1.328] | 0.977 [0.944, 1.011] | 1.167 | **2.51 → 1.95** | 9.91 → 8.53 | 0.990 [0.934, 1.049] |
| `icv1-f2m61-t158598901-ab42b6c5` | **66** | **1.267** [1.242, 1.292] | **1.312** [1.283, 1.341] | 1.003 [0.972, 1.035] | 1.159 | **2.47 → 1.90** | 16.07 → 12.91 | 0.986 [0.926, 1.049] |

**What the table says.**
- **The mechanism is right.** At the two wide sizes a scanned summand
  falls from about 2.5 units to 1.9, and collection runs 1.29–1.31×
  faster. The online interval, whose descent walks through the same
  kernel, runs 1.16–1.17× faster.
- **The cold ratio is what collection's share gives.** Collection is 72%
  of the set-up at `icv1-f2m59-tm943548413-98844ecc` and 89% at
  `icv1-f2m61-t158598901-ab42b6c5` (R01). A 1.29× and a 1.31× collection
  give 1.19× and 1.27× cold: what was measured.
- **The build did not move** (0.98 and 1.00), although it takes the
  kernel too. The batched addition is not what the build waits on at
  these sizes.
- **The nine other sizes read 0.973–1.027,** each inside its A/A band.
  Their kernels did not change.

## Against the prediction

| curve | predicted | suite | holdouts |
|:--|--:|--:|--:|
| `icv1-f2m59-tm943548413-98844ecc` | 1.31× [1.15, 1.45] | 1.182 [1.151, 1.213] | 1.161 [1.089, 1.239] |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.37× [1.18, 1.53] | 1.267 [1.242, 1.292] | 1.262 [1.216, 1.310] |

**The prediction was high.** It assumed a wide field's summand would
cost what `icv1-f2m53-tm56619371-dac20a85`'s does, with a range of 1.4
to 1.9 units. The candidate's 1.90–1.95 units sit at the top of that
range, 8–11% above `icv1-f2m53-tm56619371-dac20a85`'s 1.76–1.80 in the
same runs. Two things differ between those sizes, and this round does
not separate them:
- the carry adds a shift and an XOR per tail term to every
  multiplication;
- the two wide sizes' tables are larger: 3.9 M and 6.3 M stored pairs,
  against 3.7 M (plan §8).

R04 prices the scan's stages inside the scan.

**Why the holdout missed.** Each size's holdout set has two rows, so ten
pairs against the suite's forty. At `icv1-f2m59-tm943548413-98844ecc`
the holdouts' point estimate, 1.161, is close to the suite's 1.182, but
their interval is twice as wide. The rule asks for a lower end above
1.10 on each set separately, and the ten-pair set missed it by 0.011.
The rule was declared before the runs, so it decides: the round is
rejected.

## The holdouts

Two rows per size, on seed 205 with targets `T101` and `T102`, five
rounds.

| curve | cold [95%] | collection | build | units per scanned summand | `S` cold |
|:--|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` (nine pairs) | 1.034 [0.934, 1.145] | 1.016 | 1.056 | 25.27 → 24.76 | 134.19 → 126.98 |
| `icv1-f2m23-tm5197-1f85e9e1` | 0.918 [0.832, 1.014] | 0.924 | 0.951 | 10.48 → 10.81 | 35.01 → 35.78 |
| `icv1-f2m45-tm6236725-40939294` | 0.940 [0.883, 1.000] | 0.955 | 0.951 | 7.56 → 7.83 | 22.85 → 23.57 |
| `icv1-f2m37-tm534059-32aad96b` | 1.009 [0.983, 1.037] | 0.995 | 0.970 | 4.80 → 4.79 | 9.89 → 9.84 |
| `icv1-f2m43-tm998717-e2e742b0` | 0.989 [0.954, 1.026] | 1.004 | 0.992 | 2.68 → 2.64 | 5.87 → 5.86 |
| `icv1-f2m47-t22705043-f4e44623` | 1.016 [0.925, 1.118] | 1.008 | 0.995 | 1.48 → 1.49 | 5.01 → 4.91 |
| `icv1-f2m57-tm747311035-c1f545af` | 0.956 [0.889, 1.028] | 0.947 | 0.957 | 1.60 → 1.63 | 8.18 → 8.40 |
| `icv1-f2m41-tm2308219-7f48b14a` | 0.973 [0.855, 1.107] | 0.974 | 0.950 | 1.47 → 1.42 | 6.03 → 5.92 |
| `icv1-f2m53-tm56619371-dac20a85` | 0.979 [0.932, 1.027] | 0.990 | 0.953 | 1.84 → 1.84 | 8.94 → 8.91 |
| `icv1-f2m59-tm943548413-98844ecc` | **1.161** [1.089, 1.239] | 1.245 | 0.997 | 2.49 → 1.96 | 9.89 → 8.43 |
| `icv1-f2m61-t158598901-ab42b6c5` | **1.262** [1.216, 1.310] | 1.303 | 1.020 | 2.50 → 1.90 | 15.87 → 12.61 |

The no-regression rule applies to the suite rows, not the holdouts. Two
holdout intervals sit below 1 at their point estimates, at
`icv1-f2m23-tm5197-1f85e9e1` and `icv1-f2m45-tm6236725-40939294`, but
both include 1 or touch it, on ten pairs, at sizes whose kernel did not
change.

## v0 against v0′ (accounting)

v0 (R01's binary, SHA-256 `c776cc04…`) ran against v0′ (`c1a2e5f8`,
SHA-256 `a89e1be0…`) on `M1`'s 22 rows, five rounds ABAB, interleaved in
this round's session: ten pairs per size. The ratio is v0 over v0′. It
gates nothing (amendment 3).

| curve | cold, v0 over v0′ [95%] | R01's A/A, cold |
|:--|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 0.978 [0.885, 1.081] | 1.020 [0.930, 1.118] |
| `icv1-f2m23-tm5197-1f85e9e1` | 1.031 [0.994, 1.071] | 1.021 [0.984, 1.059] |
| `icv1-f2m45-tm6236725-40939294` | 0.970 [0.940, 1.001] | 0.985 [0.916, 1.060] |
| `icv1-f2m37-tm534059-32aad96b` | 1.002 [0.948, 1.060] | 0.976 [0.903, 1.055] |
| `icv1-f2m43-tm998717-e2e742b0` | 1.001 [0.970, 1.033] | 1.028 [0.936, 1.129] |
| `icv1-f2m47-t22705043-f4e44623` | 1.016 [0.914, 1.130] | 1.019 [0.922, 1.127] |
| `icv1-f2m57-tm747311035-c1f545af` | 0.970 [0.874, 1.077] | 1.012 [0.935, 1.095] |
| `icv1-f2m41-tm2308219-7f48b14a` | 1.065 [0.985, 1.151] | 1.075 [1.005, 1.149] |
| `icv1-f2m53-tm56619371-dac20a85` | 1.023 [0.984, 1.064] | 1.017 [0.951, 1.088] |
| `icv1-f2m59-tm943548413-98844ecc` | 0.966 [0.911, 1.024] | 0.990 [0.934, 1.049] |
| `icv1-f2m61-t158598901-ab42b6c5` | 0.988 [0.955, 1.021] | 0.986 [0.926, 1.049] |

**The renaming moved nothing measurable.** Every interval includes 1,
none lies outside its size's A/A band, and the point estimates run from
0.966 to 1.065.

**Runs in different sessions do differ.** In this round's runs, v0′'s
`S` reads 4–17% above v0's in R01's profile, size by size. The two
binaries read alike when interleaved, so the difference lies between
the sessions, not in the code. Each process also measures its own unit,
a compute-bound batched addition. That rose only 0–9% (median 3%), while
`S` rose by up to 16.5% (at `icv1-f2m53-tm56619371-dac20a85`, against 2%
for its unit). So most of the drift lies in what the unit does not
measure, plausibly the big tables' memory accesses. This round does not
test that.

What follows for the programme: only a paired ratio compares across
rounds on this host. An `S` from one session against an `S` from
another carries the drift, so the ledger quotes each round's `S` beside
its own base arm's.

## Callgrind (a control)

Amendment 3 made this step a control. Valgrind hides AVX-512, so both
arms take the portable paths, and those must be instruction-identical:
a total ratio of 1.000 per row and the same logarithm. `analyse.py`
reads "identical" as a ratio within `10⁻³` of 1.

| row | v0′ instructions | candidate instructions | v0′ over candidate | same logarithm |
|:--|--:|--:|--:|:--|
| `icv1-f2m41-tm2308219-7f48b14a` `M1-T01` | 2,141,041,508 | 2,138,855,867 | 1.00102 | yes |
| `icv1-f2m59-tm943548413-98844ecc` `M1-T01` | 18,967,935,307 | 18,951,740,727 | 1.00085 | yes |
| `icv1-f2m61-t158598901-ab42b6c5` `M1-T01` | 88,670,745,219 | 88,584,288,309 | 1.00098 | yes |

**The control fails by its declared rule.**
- The candidate's portable path runs 0.085–0.102% fewer instructions at
  all three rows.
- At `icv1-f2m41-tm2308219-7f48b14a` the ratio, 1.00102, falls outside
  the tolerance.

**Where the difference lies.**
- **Two functions hold all of it:** `witnesses_fast_scan` and
  `pairs_for_key`. The candidate does not change their source.
- **Every arithmetic kernel matches to the instruction:** the keys
  (`keys_of`), the batched inversion, the lazy addition and the field
  operations.
- **One unnamed function moved address between the two binaries,** and
  its counts cancel.

So editing `koblitz_fast.rs` changed the code the compiler generated for
those two callers, not the work they do. By the protocol's reading, the
candidate changed more than its kernel.

**What that does to the reading.** The effect is a tenth of a per cent of
the instructions. The nine sizes whose kernel did not change run the
same two functions, with the same codegen, and read 0.973–1.027 in
time, so it is not what moved the two wide sizes. It adds a second
failed condition to a decision the holdout had already made.

## Files

- [`PROTOCOL.md`](PROTOCOL.md): the declaration and amendments 1–3.
- [`run.py`](run.py) and [`analyse.py`](analyse.py): the runner and the
  analysis, as merged before the runs.
- [`analysis.json`](analysis.json): every figure above.
- `runs.tar.xz` and `runs.tar.xz.sha256`: the run tree, every attempt
  kept, contended ones included.
- [`candidate.patch`](candidate.patch): the rejected kernel, on record.
- [`holdouts/`](holdouts/): the holdouts' parameter files.

To reproduce the figures:

```
sha256sum -c runs.tar.xz.sha256 && tar -xJf runs.tar.xz && python3 analyse.py > analysis.json
```
