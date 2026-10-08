# R07: the programme re-based on main's head

The `ic` tool programme's seventh round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md`, Track A). It asks
whether main's head computes the programme's answers, and what it costs,
size by size:
- the programme's baselines v1 and v2 were candidates on v0's lineage,
  and main has gained about 700 commits since;
- among them #1242, which reduces every `Gf2` product by two carry-less
  folds with the modulus's sparse tail;
- main's `ic` is what anyone building this repository runs.

It claims no lever of its own: any gain is main's. Its class is
**engineering**. The declaration is [`PROTOCOL.md`](PROTOCOL.md), merged
in #1395 before anything ran.

## Status: accepted, as baseline v3

The protocol accepts on the pin and the A/A bands, with no minimum gain.
Every declared condition holds:

| condition | declared | measured | met |
|:--|:--|:--|:--|
| the tests | `cargo test --release --lib -- koblitz_fast koblitz_index_calculus semaev_decomp` passes on the candidate's tree | 194 passed, none failed, 3 ignored by their own attributes ([`tests.log`](tests.log)) | yes |
| the pin | every pinned output of the candidate equals v0's, on all 90 rows, and the names agree | 90 of 90 rows, names agree | yes |
| no regression | no size's cold interval wholly below the round's own A/A interval | none: every size's cold interval lies above 1 | yes |
| callgrind | both arms recover the same logarithm at both profiled sizes | the same logarithm at `icv1-f2m53-tm56619371-dac20a85` and `icv1-f2m61-t158598901-ab42b6c5` | yes |

**What follows from the acceptance:**
- **v3 is the newest baseline:** main at `995ea207`, binary SHA-256
  `9c3320390d74ae48b74382d2ec3f8b4b5b1585e5e8aa2448c0d78369807a26b4`.
  `baselines.json` carries v3, with its `S` from this round's runs beside
  v2's from the same runs.
- **Its class is engineering.** Every pinned output is v0's, so the ratio
  to the floor is flat. The gain is main's, chiefly #1242's, not a lever
  of the programme's.
- **The rule's comparison runs at v3,** under rule v2's protocol with
  v3's binary. This pull request declares it as
  [`rule/v3`](../../rule/v3/PROTOCOL.md). Rule v2 never ran and stays on
  record, superseded.
- **R06 is declared on v3** in this pull request
  ([`R06-scan-key/PROTOCOL.md`](../R06-scan-key/PROTOCOL.md)), after this
  decision.
- **R02b's go/no-go exploration runs on v3**, as its amendment 3 fixed.

## What ran

R07 ran on one host class, through three container restarts.
- **Each timed process** ran isolated on CPU 2 with
  `RAYON_NUM_THREADS=1`, on the native tools only: `icprog run r07`
  (SHA-256 `70e87395…`) and `isolated_bench` (`83418020…`).
- **The arms:**
  - the base is v2's binary, `ic-r05-cand-on-30f6c153-edcb0bec`
    (SHA-256 `76a2a2fd…`), built from `edcb0bec`;
  - the candidate is main's `ic` at `995ea207`, built with v2's
    `Cargo.lock` (SHA-256 `9c332039…`).
- **The A/A bands are the round's own** (`aa-source.json`): this host's
  kernel build is `6.18.44-fc-v70`, and R01's was `-fc-v50`.

| step | when, UTC | processes | contended, each retried | failed |
|:--|:--|--:|--:|--:|
| the tests | 2026-10-05, before the manifest | — | — | 0 |
| the manifest and the pin, untimed | 19:52–19:56 | 90 | — | 0 |
| the A/A, the base against a byte-identical copy | 19:56–22:54 | 232 | 13 | 0 |
| the suite rows, rounds 1–5 | 22:54–02:55 | 424 | 25 | 0 |
| the holdouts, rounds 1–5 | 02:55–05:09 | 251 | 11 | 0 |
| the extension | 05:09 | 0 | — | — |
| the callgrind cross-check, untimed | 05:09–05:39 | 4 | — | 0 |

### Three container restarts

**The container restarted three times,** and each time the run tree
resumed on the same host class:
- about 20:42 UTC, during the A/A, resumed at 20:48:59;
- about 21:00, during the A/A, resumed at 21:08:31;
- about 03:59 on 2026-10-06, during the holdouts, resumed at 04:01:51.
  The runner had stopped at 03:07; the container was down for most of
  the hour between.

**Each resume checked the host first.** The CPU model, its flags, the
core count, the memory and the kernel build all matched `host.json`, as
the protocol requires. Each resume wrote a manifest (`host-resumed.json`,
`host-resumed-2.json`, `host-resumed-3.json`) and a line in
`resumes.log`, and the analysis lists them all.

**Three pairs straddle a restart.** Each is paired as R05's was: its first
arm before the restart, its second after. A check after the run, outside
the analysis:

| set | pair | its ratio | the set's geometric mean, with it | without it |
|:--|:--|--:|--:|--:|
| the A/A at `icv1-f2m19-tm797-9c54981b` | `M1-T02`, round 3 | 0.956 | 0.9880 | 0.9920 |
| the A/A at `icv1-f2m41-tm2308219-7f48b14a` | `M1-T01`, round 3 | 1.152 | 1.0146 | 1.0004 |
| the holdouts at `icv1-f2m53-tm56619371-dac20a85` | `M14-T119`, round 2 | 1.200 | 1.1704 | 1.1696 |

### Refused attempts: none was timed

**The isolation tool refused 182 attempts before they ran.** A refused
attempt runs nothing and writes no figure; the runner waits and tries
again. Every refusal and every wait for pressure is logged in the step's
`refusals.log` and `psi-waits.log`.

| step | refusals | runs refused at least once | the busiest neighbour |
|:--|--:|--:|:--|
| the A/A | 31 | 27 | this session's agent 28, a kernel worker 1, `git` 1, PSI alone 1 |
| the suite rows | 123 | 93 | the agent 108, kernel workers 9, `git` 4, `objdump` 1, PSI alone 1 |
| the holdouts | 28 | 21 | the agent 24, a kernel worker 2, `du` 2 |

### Contended runs: each retried

**49 timed runs were contended,** the busiest neighbour this session's
agent in all but two (`git` once in the A/A, `du` once in the holdouts),
at 0.10–1.03 CPU seconds from other processes. No contended run is in any
figure.
- **The runner retries a contended run at most twice,** and the first
  clean attempt is the row's figure; every attempt stays in the run tree.
- **Two pairs have no clean attempt** and are left out:
  - in the A/A at `icv1-f2m19-tm797-9c54981b`, the copy's run of
    `M1-T01`, round 1, so that band rests on nine pairs;
  - in the suite rows at `icv1-f2m53-tm56619371-dac20a85`, the
    candidate's run of `M4-T08`, round 3, so that set has 39 pairs.

## The comparison: v2 against main's head

**How to read the table.**
- **The ratio** is v2 over main's head, so above 1 means main is faster.
- **The rows.** The three target sizes have eight rows each (`M1`–`M4`,
  two targets each); the eight other sizes have `M1`'s two rows. Every
  row ran five rounds: no set's half-width exceeded 3%, so none was
  extended.
- **The intervals** are 95% `t` intervals of the geometric mean of the
  paired ratios.
- **Stage diagnostics.** Online, the collection stage and rho's online
  interval are stage diagnostics, not speedups.
- **`S`** is the median cold cost in v0's unit (R01's `unit_ns`), each
  arm from this round's own runs.
- **"Units a summand"** is the collection stage's median time per summand
  scanned, in the same unit.
- **The last column** is the round's own A/A: v2 against a byte-identical
  copy of itself.

| curve | log₂ r | rows | rounds | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate | collection, units a summand | the round's A/A, cold |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 2 | 5 | 1.157 [1.074, 1.246] | 1.430 [1.318, 1.551] | 1.256 [1.112, 1.419] | 1.110 [0.905, 1.361] | 133.55 → 112.53 | 30.44 → 23.13 | 0.988 [0.889, 1.097] |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 2 | 5 | 1.210 [1.148, 1.275] | 1.538 [1.464, 1.616] | 1.345 [1.263, 1.432] | 1.488 [1.397, 1.585] | 35.28 → 30.00 | 11.51 → 8.56 | 0.992 [0.956, 1.030] |
| `icv1-f2m45-tm6236725-40939294` | 24.8 | 2 | 5 | 1.209 [1.170, 1.249] | 1.679 [1.588, 1.774] | 1.474 [1.417, 1.532] | 1.818 [1.729, 1.912] | 21.59 → 17.91 | 7.82 → 5.45 | 0.991 [0.942, 1.043] |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 2 | 5 | 1.287 [1.232, 1.343] | 1.336 [1.270, 1.405] | 1.493 [1.417, 1.572] | 1.628 [1.565, 1.693] | 10.29 → 8.25 | 4.26 → 2.89 | 1.007 [0.967, 1.048] |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 2 | 5 | 1.251 [1.198, 1.306] | 1.267 [1.211, 1.325] | 1.314 [1.225, 1.409] | 1.465 [1.396, 1.539] | 6.00 → 4.79 | 2.27 → 1.66 | 1.003 [0.979, 1.027] |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 2 | 5 | 1.248 [1.070, 1.456] | 1.081 [0.887, 1.316] | 1.198 [1.018, 1.410] | 1.386 [1.280, 1.502] | 4.12 → 3.12 | 1.22 → 0.98 | 1.010 [0.972, 1.049] |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 2 | 5 | 1.231 [1.166, 1.299] | 1.359 [1.287, 1.435] | 1.195 [1.126, 1.269] | 1.446 [1.266, 1.652] | 4.37 → 3.67 | 1.17 → 0.99 | 1.022 [0.949, 1.101] |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 2 | 5 | 1.110 [1.022, 1.206] | 1.141 [0.901, 1.447] | 1.079 [0.983, 1.185] | 1.198 [1.044, 1.375] | 5.11 → 4.38 | 1.16 → 1.01 | 1.015 [0.920, 1.119] |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 8 | 5 | **1.135 [1.104, 1.167]** | 1.229 [1.181, 1.280] | 1.049 [1.017, 1.083] | 1.345 [1.301, 1.391] | 6.66 → 5.82 | 1.17 → 1.11 | 0.998 [0.960, 1.037] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 8 | 5 | **1.483 [1.452, 1.514]** | 1.713 [1.654, 1.774] | 1.627 [1.588, 1.666] | 1.342 [1.287, 1.399] | 8.11 → 5.48 | 1.85 → 1.15 | 1.032 [0.989, 1.078] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | 5 | **1.526 [1.507, 1.546]** | 1.583 [1.525, 1.643] | 1.577 [1.555, 1.599] | 1.312 [1.259, 1.367] | 12.26 → 8.09 | 1.80 → 1.14 | 1.043 [1.006, 1.082] |

**What the table says.**
- **Main's head is faster at every size,** and no cold interval reaches
  1: 1.110–1.287 at the eight smaller sizes, and 1.135, 1.483 and 1.526
  at `2^44.3`, `2^44.5` and `2^47.2`.
- **The wide-tail penalty is gone.** At `2^44.5` and `2^47.2`, the two
  fields whose polynomial has `n + deg t = 66`, a scanned summand falls
  from 1.85 and 1.80 units to 1.15 and 1.14, the 8-lane kernel's cost at
  `2^44.3` (1.11). The scan there runs the scalar subtraction, and #1242
  makes each of its products about 2.4 times cheaper. R02 and R02b
  targeted this gap with a kernel; main closed it underneath them.
- **`2^44.3` gains least,** 1.135 cold and 1.049 in collection: its scan
  runs the 8-lane kernel, which #1242 does not touch. The gain there comes
  mostly from the build and the field operations every phase shares.
- **`S` at `2^44.5` is now below `S` at `2^44.3`:** 5.48 against 5.82.
- **Rho's online interval is faster too,** 1.110–1.818: the reference
  walk shares the field arithmetic, and main rewrote it after v2 (#1334,
  #1360). The rule's comparison at v3 measures where that leaves the
  index calculus against rho.
- **The A/A bands** span 0.988–1.043 in their geometric means. One
  interval excludes 1: 1.043 [1.006, 1.082] at `2^47.2`.

## Against the prediction

The prediction came from the diagnostic before the declaration: one
process per arm, set-up time only.

| curve | log₂ r | set-up, diagnostic | cold, suite | cold, holdouts |
|:--|--:|--:|--:|--:|
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 1.32 | 1.248 [1.070, 1.456] | — |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | **1.05** | 1.135 [1.104, 1.167] | 1.170 [1.141, 1.201] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | **1.38** | 1.483 [1.452, 1.514] | 1.483 [1.445, 1.523] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.54** | 1.526 [1.507, 1.546] | 1.539 [1.521, 1.557] |

- **At `2^47.2` the diagnostic held** within 1%.
- **At `2^44.3` and `2^44.5` it was low,** by 8% and 7%.
- **The diagnostic's one surprise did not recur.** It read collection at
  `2^44.3` as slower on main (459 against 509 ms). The round reads it
  1.049 [1.017, 1.083] times faster, as the protocol expected a single
  process might mislead.
- **Elsewhere** the protocol predicted about 1.1–1.35; the eight smaller
  sizes read 1.110–1.287.

## The holdouts

**Fresh rows, five rounds each.**
- The frozen files in [`holdouts/`](holdouts/) (`SHA256SUMS`), at the
  three target sizes: recipe seeds 214–217, two targets each,
  `T119`–`T126`, drawn by `icprog holdouts`.
- None was extended: their half-widths after five rounds were 1.2–2.7%
  (the analysis's `extension_test`).

| curve | log₂ r | rows | rounds | pairs | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 8 | 5 | 40 | 1.170 [1.141, 1.201] | 1.265 [1.201, 1.332] | 1.081 [1.050, 1.112] | 1.358 [1.298, 1.421] | 7.10 → 5.96 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 8 | 5 | 40 | 1.483 [1.445, 1.523] | 1.715 [1.640, 1.794] | 1.635 [1.586, 1.686] | 1.367 [1.301, 1.435] | 8.12 → 5.40 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | 5 | 40 | 1.539 [1.521, 1.557] | 1.683 [1.623, 1.747] | 1.597 [1.572, 1.622] | 1.374 [1.327, 1.423] | 11.91 → 7.69 |

**The holdouts agree with the suite at every target size:** 1.170 against
1.135, 1.483 against 1.483, and 1.539 against 1.526.

## The callgrind cross-check

**What ran.** `ic workflow` under callgrind on `M1`'s first target at
`icv1-f2m61-t158598901-ab42b6c5` and `icv1-f2m53-tm56619371-dac20a85`,
both arms, untimed. Valgrind hides AVX-512 from both arms, so each runs
its portable paths; the carry-less multiply is still there, so #1242's
reduction is in the counts.

| curve | instructions, v2 over main's head | the same logarithm | functions that differ |
|:--|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.417 | yes | 37 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.365 | yes | 16 |

- **It decides one thing, and that held:** both arms recover the same
  logarithm at both sizes. The protocol makes it a cross-check (plan §5),
  not a control; the counts are reported, not tested.
- **Main runs 29% and 27% fewer instructions** on the portable paths. The
  functions that differ include the field's multiply, square, inversion
  and batched inversion, the curve's ladder and batched additions, and
  the key; the analysis lists every one with both counts.
- **The analysis's callgrind block now says so.** `icprog analyse`
  labelled every round's callgrind block with R02b's control rule. A
  cross-check's block now carries its own role and rule
  (`cross_check_held`), and R02b's control keeps its own; nothing else
  in the analysis changed.

## Files

- **`runs.tar.xz`** (`runs.tar.xz.sha256`): the run tree.
  - It holds every process's report, isolation record and stderr, every
    attempt of a retried run, the pin, the four host manifests, the
    callgrind profiles, and the refusal and PSI logs.
  - It reproduces `analysis.json` byte for byte; `tests/icprog.rs`
    checks this. From the repository root:

    ```
    tar -xJf research/ic_tool_program/rounds/R07-main-head/runs.tar.xz \
      -C research/ic_tool_program/rounds/R07-main-head
    cargo run --release --bin icprog -- analyse r07 \
      --runs research/ic_tool_program/rounds/R07-main-head/runs > analysis.json
    ```
- **`analysis.json`**: `icprog analyse r07`'s output from the archive.
- **The tables above** are `icprog table r07 analysis.json`'s output. The
  prediction table's measured columns are the same figures.
- **`tests.log`**: the candidate tree's tests, run before the timed steps.
- **`baseline.meta.json`**: v3's identity. `icprog baseline` appended v3's
  entry to `baselines.json` from it, from the analysis and from the run
  tree's host manifests.
