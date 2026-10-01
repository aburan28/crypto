# R01: baseline v0

The `ic` tool programme's first round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §12). It:
- freezes suite v1;
- pins v0 against ledger §23;
- profiles every phase;
- measures the host's noise floor.

It has no candidate and claims no speedup. Its class is **accounting**.
The declaration is [`PROTOCOL.md`](PROTOCOL.md), merged in #1104 before
anything ran.

## Status: complete

- **The pin held.** v0's outputs equal §23's on all 12 shared rows.
- **Every timed row was verified.** Every process exited cleanly; seven
  contended attempts are kept, and their reruns are the figures.
- **All declared steps ran**: the profile (88 rows), the A/A (220
  processes), the huge-page probe (60), the memory calibration (30
  points) and callgrind (two rows).

| step | processes | contended | failed |
|:--|--:|--:|--:|
| profile pass | 95 | 7, each rerun clean | 0 |
| A/A | 220 | 0 | 0 |
| huge-page probe | 60 | 0 | 0 |
| memory calibration | 30 | 0 | 0 |

## Curves

Each curve is named by its compact alias, which links to the full
identity (`docs/curve-identities.md`). The fields are the ones
`find_irreducible_sparse` picks, and the registry's moduli match them.

| alias | EC1 | curve UID | `n + deg t` |
|:--|:--|:--|--:|
| `icv1-f2m19-tm797-9c54981b` (K_1, n = 19) | `EC1N19Ce1h7e63aa468dcc` | `urn:ec-record:1:sha256:7e63aa468dcc9ee8813f014acabe3395ac3e3d9f533d00361a821770038f4d23` | — |
| `icv1-f2m23-tm5197-1f85e9e1` (K_1, 23) | `EC1N23Ce1heec6a3619e3f` | `urn:ec-record:1:sha256:eec6a3619e3fdf05035e0ef5d3ec8cf3af507268cf760c77d47518d62141b00f` | — |
| `icv1-f2m45-tm6236725-40939294` (K_1, 45) | `EC1N45Ce1h97bc367649de` | `urn:ec-record:1:sha256:97bc367649ded9e1309ece1e66c630f558ee953b017378c53b2c3cf2896747df` | — |
| `icv1-f2m37-tm534059-32aad96b` (K_0, 37) | `EC1N37Ce0h0c51aa4aa7c3` | `urn:ec-record:1:sha256:0c51aa4aa7c3c76c45885e230d717aedcbefad98e2604d1664ac38b57299f180` | — |
| `icv1-f2m43-tm998717-e2e742b0` (K_1, 43) | `EC1N43Ce1h7d34f0f22c70` | `urn:ec-record:1:sha256:7d34f0f22c707e3756957364281079ea2835960013fd13a1e648b7ff5e5fea8c` | — |
| `icv1-f2m47-t22705043-f4e44623` (K_1, 47) | `EC1N47Ce1h9669b6d5105b` | `urn:ec-record:1:sha256:9669b6d5105b035b2ff286a4e32c52725841e8546ab03de8a8077c5114e98b0c` | 52 |
| `icv1-f2m57-tm747311035-c1f545af` (K_0, 57) | `EC1N57Ce0hb5884395c7f2` | `urn:ec-record:1:sha256:b5884395c7f28d0d22eb1acc451c34ff979b0a40f93629e44b30107592b76695` | 61 |
| `icv1-f2m41-tm2308219-7f48b14a` (K_0, 41) | `EC1N41Ce0he09550ab560a` | `urn:ec-record:1:sha256:e09550ab560a9bb6de7981f4b47c48235d8553626dc43dbe92bf73262a95a8e1` | 44 |
| `icv1-f2m53-tm56619371-dac20a85` (K_0, 53) | `EC1N53Ce0hb097de99be9a` | `urn:ec-record:1:sha256:b097de99be9a3af970966707e65d24bf63302ba5324f05b714f916a680651ae6` | 59 |
| `icv1-f2m59-tm943548413-98844ecc` (K_1, 59) | `EC1N59Ce1h20858f64459b` | `urn:ec-record:1:sha256:20858f64459bbef8a18a8b75dbe579ce4cb0901d5c1b70d236ac154427d86980` | **66** |
| `icv1-f2m61-t158598901-ab42b6c5` (K_0, 61) | `EC1N61Ce0hd8c15448081e` | `urn:ec-record:1:sha256:d8c15448081eb349784c3d5d9d3144eaa4a5252d115f9049e10bbadf373db689` | **66** |

`n + deg t` is the field degree plus the degree of the field
polynomial's tail. It decides whether the AVX-512 batched addition runs
(≤ 65) or the scalar path does. It is given for the sizes where the scan
dominates.

## The profile

These are v0's medians over each size's eight rows (seeds 201–204 × two
targets). One thread, isolated. `S` is units over `√r` in v0's own unit,
which is this round's pinned unit (IC_TOOL_PROGRAM.md §5). Collection and
build are shares of the set-up.

| curve | log₂ r | set-up | online | rho online | `S` cold | `S` rho online | collection | build | units per scanned summand | units per stored pair |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 1.34 ms | 0.01 ms | 0.05 ms | 120.56 | 4.03 | 7% | 13% | 23.79 | 5.75 |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 1.51 ms | 0.02 ms | 0.22 ms | 33.48 | 4.78 | 11% | 12% | 9.84 | 4.98 |
| `icv1-f2m45-tm6236725-40939294` | 24.8 | 3.14 ms | 0.02 ms | 0.56 ms | 21.78 | 3.82 | 6% | 11% | 7.76 | 3.70 |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 3.11 ms | 0.06 ms | 0.57 ms | 8.53 | 1.54 | 28% | 8% | 4.68 | 4.25 |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 11 ms | 0.22 ms | 1.79 ms | 5.87 | 0.97 | 53% | 7% | 2.54 | 2.49 |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 39 ms | 0.69 ms | 6.99 ms | 4.43 | 0.77 | 41% | 27% | 1.41 | 2.48 |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 118 ms | 0.48 ms | 7.74 ms | 7.65 | 0.50 | 26% | 14% | 1.47 | 2.31 |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 96 ms | 1.29 ms | 17 ms | 5.42 | 0.92 | 69% | 17% | 1.36 | 2.60 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.03 s | 5.05 ms | 75 ms | 7.53 | 0.55 | 66% | 27% | 1.49 | 2.56 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.35 s | 3.01 ms | 58 ms | 8.77 | 0.37 | 72% | 22% | **2.18** | 2.48 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 5.70 s | 6.32 ms | 188 ms | 14.49 | 0.48 | 89% | 9% | **2.17** | 2.69 |

**What the profile says.**
- **The cold cost is the set-up at every size.** The online interval is
  0.1–2% of it.
- **Below `2^28` the set-up is fixed cost.** At those sizes:
  - verification, the curve, selection and the collection's set-up
    dominate;
  - `S` is 8.5–121 because a few milliseconds of fixed work sit over a
    tiny `√r`.
- **From `2^32` up, collection is the largest phase**: 41–89% of the
  set-up, except at `icv1-f2m57-tm747311035-c1f545af`. There the curve's
  construction takes 45%, the artefact §23 found.
- **The scan costs 2.17–2.18 units a summand at exactly the two sizes
  with `n + deg t = 66`**: `icv1-f2m59-tm943548413-98844ecc` and
  `icv1-f2m61-t158598901-ab42b6c5`. Every other size from `2^36.6` up
  costs 1.36–1.49.
  - That is §23's split, measured again on v0.
  - It is what R02 tests: the AVX-512 kernel refuses those two fields.
- **The build costs 2.3–2.7 units a stored pair from `2^32` up,** with
  no step at the two wide sizes.

## The A/A: v0 against a byte-identical copy

Each size has ten pairs (`M1`'s two rows × five rounds). The ratio is the
geometric mean of `A` over `A2`, with a 95% `t` interval.

| curve | cold [95%] | online | rho online |
|:--|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 1.020 [0.930, 1.118] | 1.032 [0.970, 1.099] | 0.983 [0.848, 1.139] |
| `icv1-f2m23-tm5197-1f85e9e1` | 1.021 [0.984, 1.059] | 1.026 [0.973, 1.082] | 1.010 [0.973, 1.049] |
| `icv1-f2m45-tm6236725-40939294` | 0.985 [0.916, 1.060] | 0.977 [0.901, 1.061] | 0.985 [0.923, 1.050] |
| `icv1-f2m37-tm534059-32aad96b` | 0.976 [0.903, 1.055] | 0.993 [0.925, 1.067] | 0.992 [0.922, 1.068] |
| `icv1-f2m43-tm998717-e2e742b0` | 1.028 [0.936, 1.129] | 1.016 [0.940, 1.098] | 1.035 [0.944, 1.135] |
| `icv1-f2m47-t22705043-f4e44623` | 1.019 [0.922, 1.127] | 0.989 [0.895, 1.093] | 1.017 [0.924, 1.119] |
| `icv1-f2m57-tm747311035-c1f545af` | 1.012 [0.935, 1.095] | 1.008 [0.930, 1.092] | 1.011 [0.908, 1.127] |
| `icv1-f2m41-tm2308219-7f48b14a` | 1.075 [1.005, 1.149] | 1.027 [0.891, 1.184] | 1.025 [0.929, 1.130] |
| `icv1-f2m53-tm56619371-dac20a85` | 1.017 [0.951, 1.088] | 1.039 [0.913, 1.183] | 0.995 [0.921, 1.074] |
| `icv1-f2m59-tm943548413-98844ecc` | 0.990 [0.934, 1.049] | 0.923 [0.826, 1.032] | 0.923 [0.819, 1.040] |
| `icv1-f2m61-t158598901-ab42b6c5` | 0.986 [0.926, 1.049] | 0.959 [0.856, 1.075] | 1.018 [0.958, 1.080] |

- **The cold band is ±4–11% at ten pairs.**
- **One A/A interval excludes 1.** At `icv1-f2m41-tm2308219-7f48b14a` the
  cold ratio reads 1.075 [1.005, 1.149]. Identical binaries differ there
  by about 7% on this host. This is the noise floor, and it is why a
  ratio inside the band is not a gain.
- **What this means for the rounds after it.** R01's protocol tied R02's
  round count to this band. R02's amendment 1, declared before R02 ran,
  explains why its 40-pair design meets that rule's intent at five
  rounds.

## The huge-page probe (a diagnostic, not a candidate)

v0 ran with `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` (`thp`) and without
it (`4k`), on `M1`'s rows at the three largest sizes, for five rounds.
The ratio is `4k` over `thp`, so a ratio above 1 means huge pages were
faster.

| curve | cold | collection | build | process user s, 4k / thp | process system s, 4k / thp |
|:--|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 0.963 [0.877, 1.059] | 1.120 [1.045, 1.199] | 0.668 [0.562, 0.794] | 4.46 / 4.38 | 0.12 / 0.55 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.000 [0.922, 1.085] | 1.102 [1.010, 1.202] | 0.970 [0.879, 1.070] | 5.30 / 5.21 | 0.18 / 0.46 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.004 [0.956, 1.054] | 1.053 [1.004, 1.105] | 0.697 [0.658, 0.738] | 18.08 / 17.81 | 0.28 / 0.92 |

**Huge pages are not a lever, at least not through glibc's tunable.**
- Collection gets 5–12% faster.
- The build gets up to 1.5× slower, because faulting in huge pages costs
  kernel time when the table is first touched. System time rises 3–5×.
- Cold time is unchanged within the band.

A candidate that asks for huge pages for the table alone, after the
build, would be a later round's question. This probe does not support
one.

## The memory calibration

The tool is the earlier study's `mem_calib.c`. It measures a random
pointer chase, one load at a time, with the median over runs, on this
host.

| working set | 16 KiB | 256 KiB | 1 MiB | 2 MiB | 4 MiB | 8 MiB | 16 MiB | 32 MiB | 64 MiB | 256 MiB | 1 GiB |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| ns, 4 KiB pages | 1.3 | 5.1 | 7.6 | 44.0 | 58.2 | 63.4 | 70.2 | 65.7 | 137.4 | 221.1 | 269.7 |
| ns, huge pages | 1.5 | 5.1 | 8.4 | 41.0 | 58.3 | 58.2 | 65.1 | 92.2 | 151.9 | 180.5 | 215.6 |

- **Latency jumps past 1 MiB.** It is 44–70 ns up to 32 MiB, and
  221–270 ns at 256 MiB–1 GiB.
- **Huge pages help only the largest working sets,** by 18–20% at 256
  MiB and above.
- **Memory-level parallelism is high.** At 256 MiB, 8 independent chases
  take 222 ns a step, against 207 ns for one. That is about 7–8 loads in
  flight.
- **The pair tables sit in the middle range.** At the two largest sizes
  they are 3.7–6.3 M entries. Their filters are 2–4 MiB, so a filter
  probe lands at 44–58 ns.

## Callgrind

The declared step ran `ic price --repeats 1` on `M1-T01` at
`icv1-f2m41-tm2308219-7f48b14a` and `icv1-f2m61-t158598901-ab42b6c5`.
The simulated I1 and D1 are 32 KiB each and the LL is 2 MiB.

**Each run is split into parts, one per phase boundary.**
- `ic_measurement` asks valgrind to dump its statistics at every phase
  boundary, and labels each part with the phase that just ended
  (`DUMP_STATS_AT`).
- So the run leaves one file per part. `../../harness/callgrind_phases.py`
  sums them by label.
- The step's own annotation call read only the empty base file. The
  per-phase split below re-reads every part.

**What the split shows at `icv1-f2m41-tm2308219-7f48b14a`:**

| phase | parts | instructions | D1 read misses | LL read misses |
|:--|--:|--:|--:|--:|
| set-up, and everything else outside the intervals | 4 | 3.83 × 10⁹ | 15.4 M | 283 k |
| rho online | 2 | 2.08 × 10⁸ | 3.90 M | 48 |
| target PDP | 780 | 2.36 × 10⁷ | 122 k | 6.3 k |
| target query | 779 | 2.35 × 10⁷ | 64 k | 13 |
| recovery check, relation check, descent | 5 | 9.0 × 10⁴ | 350 | 8 |

**What this changes: the declared profile mostly measures the pricer,
not the index calculus.**
- The first set-up part runs from process start to the first online
  boundary, 3.82 × 10⁹ instructions in all.
- Of that, the pricer's unit calibration (`UnitBench::measure`) is
  1.56 × 10⁹. The index-calculus pass up to that boundary is
  1.2 × 10⁸.
- So the declared profile measures the pricer more than the pipeline.

**What R02 does about it.** R02's amendment 2 profiles `ic workflow`,
which has no calibration, for its cross-check.

**At `icv1-f2m61-t158598901-ab42b6c5` the set-up dominates the profile.**
- The run has 6,696 parts and 9.34 × 10¹⁰ instructions.
- The first set-up part alone is 9.07 × 10¹⁰. The calibration is the same
  size as at `icv1-f2m41-tm2308219-7f48b14a`, so here the set-up
  dominates.
- The online phases add 2.65 × 10⁸, and rho's walk 2.45 × 10⁹.

| function (first set-up part, exclusive) | instructions | share | LL read misses |
|:--|--:|--:|--:|
| `PairSumTable::keys_of` | 42.05 × 10⁹ | 46.4% | 0 |
| `Gf2::batch_inv` | 16.28 × 10⁹ | 18.0% | 0.03 M |
| `FastCurve::add_many_lazy` | 8.31 × 10⁹ | 9.2% | 1.84 M |
| `PairSumTable::witnesses_fast_scan` | 6.60 × 10⁹ | 7.3% | 66.30 M |
| `Gf2::sqr` | 5.93 × 10⁹ | 6.5% | 0 |
| `clmul_u64` | 3.23 × 10⁹ | 3.6% | 0 |
| `PairSumTable::compact_contains` | 1.63 × 10⁹ | 1.8% | 32.99 M |
| `FrobeniusCanon::canon_in_place` | 1.39 × 10⁹ | 1.5% | 0.03 M |

**Reading it.**
- **The scan's canonical keys are the largest single cost.** `keys_of`
  is 46% of the set-up's instructions: about 565 a scanned summand, for
  74.4 M summands. At `icv1-f2m41-tm2308219-7f48b14a` it is about 390.
  It misses no cache.
  - This is the next lever after R02, and nothing in the plan's backlog
    named it. It is added as A11.
- **The scalar addition path is 38%.**
  - `Gf2::batch_inv`, `sqr`, `clmul_u64`, `inv` and `add_many_lazy`'s
    own body are the scalar path's.
  - At this size the AVX-512 kernel refuses the field, which is R02's
    hypothesis seen from the instruction side.
  - The matching split where the kernel runs needs a profile without the
    pricer's calibration. At `icv1-f2m41-tm2308219-7f48b14a` the first
    part is mostly calibration. R02's workflow profiles, under its
    amendment 2, give that split.
- **The cache misses sit in the filter and the compact table.**
  `witnesses_fast_scan` and `compact_contains` take 99 M last-level read
  misses against the simulated 2 MiB LL. The prefetches hide part of
  them. How much is a timing question callgrind cannot answer.

## v0's ledger row

| baseline | commit | class | `S` cold at the six top sizes (v0 unit) | online / rho (§23) | cold / rho (§23) | C suite | largest `n` at F0 | PR |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| v0 | `46ae2014` (`src/` tree `003badc2`) | accounting | 4.43, 7.65, 5.42, 7.53, 8.77, 14.49 | IC 8.8–17.5× faster | IC 4.9–31.6× slower | none yet | 61 | #1104 (declaration), this PR (results) |

The pin held, so §23's single-target figures are v0's. v0's binary is
`c776cc04…`; §23's was `6f58b075…`. They are built from `src/` trees that
differ only in two hyperelliptic files the pipeline does not call.

## Correction to the host manifest

`runs/host.json` (in `runs.tar.xz`) records the build commit as `4afd2990b6fd`. That was a
typing error at launch. The binary was built from
`4afd29903e4fcf65c1d7096083bbe9b9f5ec0a66`, whose `src/` tree is
`003badc259bcef6e56ad19ac8ae93df4221ffe6e`, as the same manifest
records. The manifest is left as it was written, and this note is the
correction.

## Files

| path | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration (#1104) |
| `run.py`, `analyse.py` | the steps and the figures |
| `analysis.json` | every figure above, derived from the run tree |
| `runs.tar.xz`, `runs.tar.xz.sha256` | the run tree, 1,282 files, as one archive: 2,083,252 bytes, SHA-256 `ae20e0bd…` |

The run tree is committed as one archive because its directories are
named by the suite's row ids. AGENTS.md §11, which landed after the run,
refuses those stems in new file names. It also says evidence keeps the
names it was written with. The archive does both: every name and byte
inside it is as the run wrote it. To audit the figures:

```
cd research/ic_tool_program/rounds/R01-baseline-v0
sha256sum -c runs.tar.xz.sha256
tar -xJf runs.tar.xz            # runs/, which .gitignore keeps out of git
python3 analyse.py | cmp - analysis.json
```

The archive is deterministic (`tar --sort=name --mtime='2026-10-01 00:00:00Z'
--owner=0 --group=0 --numeric-owner --mode='u+rwX,go+rX,go-w' --format=gnu`,
then `xz -6 -T1`), so packing the extracted tree again gives the same
bytes. Inside it:

| path in the archive | what it is |
|:--|:--|
| `runs/host.json` | the host manifest (see the correction) |
| `runs/pin/`, `runs/smoke/` | the pin against §23, and the smoke rows |
| `runs/profile/`, `runs/aa/`, `runs/thp/` | each process's report, isolation record and stderr; `-retry1` files are reruns of contended attempts, which are kept |
| `runs/calib/` | the chase's outputs, its source hash and its flags |
| `runs/callgrind/` | the per-phase splits (`*.phases.json`), the valgrind logs and the reports; the raw parts are in `callgrind-parts.tar.xz`, with its SHA-256 in `callgrind-parts.sha256` |
| `runs/run.log`, `runs/psi-waits.log`, `runs/*/refusals.log` | the run's log, the PSI waits and every refused start |
