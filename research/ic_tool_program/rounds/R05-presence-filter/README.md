# R05: a sharper presence filter

The `ic` tool programme's fifth round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md`, Track A, A2). It
tests one lever on the collection scan's admitted stage:
- the presence filter set one bit a key under one hash, at four to eight
  bits a stored pair, and R04 found nearly every key it admitted at the
  top sizes to be a false positive;
- the candidate sets three bits a key in one 64-bit word, at eight to
  sixteen bits a stored pair, so a probe still reads one word.

Its class is **engineering**. The declaration is
[`PROTOCOL.md`](PROTOCOL.md), merged in #1164 before anything ran.

## Status: accepted, as baseline v2

Every declared condition holds:

| condition | declared | measured | met |
|:--|:--|:--|:--|
| the tests | `koblitz_index_calculus`'s library tests and `make -C gpu/ecc2k test` pass on the candidate's tree | 134 passed, none failed; every GPU host suite passed, ThreadSanitizer included ([`tests.log`](tests.log)) | yes |
| the pin | every pinned output of the candidate equals v0's, on all 90 rows | 90 of 90 rows, names agree | yes |
| `icv1-f2m53-tm56619371-dac20a85`, suite | cold interval above 1.10 | 1.293 [1.270, 1.317] | yes |
| `icv1-f2m53-tm56619371-dac20a85`, holdouts | cold interval above 1.10 | 1.276 [1.242, 1.311] | yes |
| `icv1-f2m59-tm943548413-98844ecc`, suite | cold interval above 1.10 | 1.238 [1.215, 1.262] | yes |
| `icv1-f2m59-tm943548413-98844ecc`, holdouts | cold interval above 1.10 | 1.240 [1.208, 1.273] | yes |
| `icv1-f2m61-t158598901-ab42b6c5`, suite | cold interval above 1.10 | 1.303 [1.281, 1.324] | yes |
| `icv1-f2m61-t158598901-ab42b6c5`, holdouts | cold interval above 1.10 | 1.319 [1.298, 1.340] | yes |
| no regression | no size's interval wholly below its A/A band | none | yes |

**What follows from the acceptance:**
- **v2 is the newest baseline:** `edcb0bec` on v1 `30f6c153`. This pull
  request applies `candidate.patch` to `src/`, `examples/` and
  `gpu/ecc2k/`, and `baselines.json` carries v2.
- **Its class is engineering.** Every pinned output is v0's. The filter
  decides only which keys reach a bucket probe, and it has no false
  negatives, so nothing the pipeline decides changes and the ratio to the
  floor is flat. `S` falls at every size from `2^36.6` up.
- **As declared, three things change besides time:**
  - the filter's memory doubles, to 4, 4 and 8 MiB at the target sizes;
  - the tables' byte estimates count the filter at eight bits a stored
    pair, not four;
  - `perfbench`'s pair-table fingerprint changes. Nothing pins it.
- **R02b runs on v2**, as its amendment 1 fixed.
- **The rule's comparison is owed** at the new baseline: §23's protocol
  at its six sizes, with 64 targets, on the candidate. It runs on the
  native port of §23's scripts (plan §10a, N3), which is not built yet.
  Until it runs, the ledger's online/rho and cold/rho columns read "not
  re-measured" for v2.

## What ran

R05 ran in two container sessions on one hardware class, because the
container was rebuilt during the holdouts (below).
- **Each timed process** ran isolated on CPU 2 with
  `RAYON_NUM_THREADS=1`, after PSI fell below 4.0.
- **The arms:**
  - the base is v1's binary, `ic-r03-cand-on-c1a2e5f8-30f6c153`
    (SHA-256 `23c83a47…`);
  - the candidate is built from `edcb0bec`, which is v1 plus
    `candidate.patch` and nothing else (SHA-256 `76a2a2fd…`).
- **The A/A bands are R01's.** The first session's manifest matches
  R01's (`aa-source.json`).

| step | runner | when, UTC | processes | contended, each rerun | failed |
|:--|:--|:--|--:|--:|--:|
| the tests | `cargo test`, `make -C gpu/ecc2k test` | 18:43–18:52 | — | — | 0 |
| the pin, untimed | declared | 18:52–18:58 | 90 | — | 0 |
| the suite rows, rounds 1–5 | declared | 18:58–21:56 | 412 | 12 | 0 |
| the holdouts, rounds 1–2 and most of round 3 | declared | 21:56–22:39 | 141 | 0 | 0 |
| the container rebuilt | — | about 22:47 | — | — | — |
| the holdouts, the rest of round 3 and rounds 4–5 | native | 23:00–23:26 | 101 | 2 | 0 |
| the extension: rounds 6–10 of the suite rows at two sizes | native | 23:26–23:53 | 160 | 0 | 0 |

- **The declared runner** is [`run.py`](run.py) on
  `harness/bench.py`, with `tools/isolated_bench.py`.
- **The native runner** is `icprog run r05`, with `isolated_bench`.
  - Both were built at 22:59 UTC from #1183's working tree, which was
    committed two minutes later as `0efdb104`. Its `src/bin/` is main's
    at `f7abedc0`.
  - Their hashes identify them, not a rebuild: a rebuild of `0efdb104`
    on this host gives other bytes. `icprog` is `f35d6c93…`
    ([`resume.log`](resume.log)), and `isolated_bench` is `219d1dfa…`
    (`host-resumed.json`).

Every row and round has a clean pair: the analysis found no missing pair
at any size.

### The rebuild and the native resume

**What happened.**
- At about 22:47 UTC the container was rebuilt, which killed the declared
  runner during the holdouts.
- By then rounds 1 and 2 were complete. Round 3 had 22 of its 24 pairs,
  and the base arm of a 23rd.
- AGENTS.md's rule against Python tooling had merged at 22:00 UTC. That
  was after the suite rows finished, while the holdouts ran. The rule
  asks for native tooling before a script's next use, so the rest ran on
  the native runner, not on a restarted `run.py`.

**How it resumed.**
- From 23:00 UTC, `icprog run r05 holdout` resumed the same run tree. It
  ran each missing attempt in the declared order and overwrote none.
- `icprog run r05 extend` then applied the extension rule.

**The rebuilt host.**
- It matches the first in CPU model and flags, cores, memory, THP,
  `rustc` and glibc.
- Its kernel build differs: `6.18.44-fc-v51` against `-v50`.
- Its manifest is `host-resumed.json`, beside `host.json`.
- Each process's isolation record names its tool. The analysis counts
  them (`isolation_by_tool`): in the suite rows 412 by the Python tool and
  160 natively, and in the holdouts 141 and 101.

**The split pair.** One pair straddles the rebuild:
`icv1-f2m61-t158598901-ab42b6c5`, row `M13-T117`, round 3.
- Its base arm ran at 22:38 UTC, on the first host and the declared
  runner. Its candidate arm ran at 23:00, on the rebuilt host and the
  native runner.
- It is paired as the declared runner pairs any round.
- **A check after the run, outside the analysis** (`jq` and `bc` on the
  archive):
  - its ratio, 1.309, lies inside its row's other four rounds,
    1.255–1.340;
  - without it, the set's geometric mean reads 1.3192, against 1.3189
    with it.

### Refused attempts: none was timed

**The isolation tools refused 82 attempts before they ran.** A refused
attempt runs no benchmark and writes no figure; the runner waited and
tried the same row again.
- **Every refusal was the machine busy:** other processes above 0.20 CPU
  s in 2 s, or PSI above 5.
- **Every attempt is logged** in its step's `refusals.log`, and every
  wait for pressure in `psi-waits.log`.

| step | refusals | runs refused at least once | the busiest neighbour | PSI alone |
|:--|--:|--:|:--|--:|
| the suite rows | 61 | 46 | this session's agent 58, a `python3` 3 | 9 |
| the holdouts | 21 | 10 | the agent 11, `git` 10 | 1 |

**Nine of the `git` refusals were this session's own.** They came at
23:02–23:04 UTC, from a `git worktree add` started in the background
while the resumed holdouts ran. The tenth came at 22:08.

### Contended runs: each rerun

**Fourteen timed runs were contended.** Each was rerun once, and every
rerun was clean. No contended run is in any figure.
- **In the suite rows, 12, all under the Python tool:**
  - in 11 the agent was the busiest neighbour, and other processes used
    0.16–0.72 CPU s in all;
  - in 1 a `python3` was, at 2.27 of 2.41 CPU s.
- **In the holdouts, 2, both under the native tool,** both at
  `icv1-f2m53-tm56619371-dac20a85`, round 4:
  - one by the same `git worktree add`, at 0.75 of 0.81 CPU s;
  - one by the agent, at 0.42 of 0.43 CPU s.

### One difference in the logs

**The two runners write refusals differently.** The native runner writes
each refusal's neighbour map as JSON, with double quotes. The Python
runner wrote Python's representation, with single quotes. No figure
reads these logs.

## The comparison: v1 against the candidate

**How to read the table.**
- **The ratio** is v1 over the candidate, so above 1 means the candidate
  is faster.
- **The rows.** The three target sizes have eight rows each (`M1`–`M4`,
  two targets each); the eight other sizes have `M1`'s two rows.
- **The rounds.** Every row ran five rounds, and ten at the two sizes the
  extension rule reached.
- **The intervals** are 95% `t` intervals of the geometric mean of the
  paired ratios.
- **Stage diagnostics.** Online, the collection stage and rho's online
  interval are stage diagnostics, not speedups.
- **`S`** is the median cold cost in v0's unit (R01's `unit_ns`), each
  arm from this round's own runs.
- **"Units a summand"** is the collection stage's median time per summand
  scanned, in the same unit.

| curve | log₂ r | rows | rounds | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate | collection, units a summand | R01's A/A, cold |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 2 | 5 | 0.939 [0.857, 1.028] | 0.933 [0.812, 1.073] | 1.003 [0.891, 1.130] | 0.936 [0.811, 1.081] | 130.53 → 136.54 | 29.89 → 30.86 | 1.020 [0.930, 1.118] |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 2 | 5 | 0.933 [0.861, 1.011] | 1.026 [0.962, 1.094] | 0.992 [0.936, 1.051] | 0.948 [0.879, 1.022] | 33.19 → 35.60 | 11.48 → 11.14 | 1.021 [0.984, 1.059] |
| `icv1-f2m45-tm6236725-40939294` | 24.8 | 2 | 5 | 1.012 [0.964, 1.063] | 1.076 [0.988, 1.171] | 1.111 [1.003, 1.232] | 1.038 [0.973, 1.108] | 22.67 → 22.05 | 8.63 → 7.98 | 0.985 [0.916, 1.060] |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 2 | 5 | 1.030 [0.946, 1.121] | 1.105 [1.043, 1.171] | 1.078 [0.995, 1.169] | 1.016 [0.937, 1.102] | 10.92 → 10.45 | 4.74 → 4.30 | 0.976 [0.903, 1.055] |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 2 | 5 | 1.023 [0.978, 1.070] | 1.052 [0.952, 1.163] | 1.095 [1.046, 1.146] | 0.946 [0.882, 1.015] | 6.63 → 6.45 | 2.67 → 2.35 | 1.028 [0.936, 1.129] |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 2 | 5 | 1.055 [0.923, 1.205] | 1.347 [1.178, 1.540] | 1.208 [1.050, 1.390] | 0.983 [0.878, 1.102] | 4.80 → 4.47 | 1.53 → 1.28 | 1.019 [0.922, 1.127] |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 2 | 5 | 1.139 [1.029, 1.261] | 1.210 [1.098, 1.332] | 1.293 [1.154, 1.448] | 0.952 [0.818, 1.108] | 4.96 → 4.45 | 1.53 → 1.22 | 1.012 [0.935, 1.095] |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 2 | 5 | 1.074 [0.973, 1.184] | 1.116 [1.028, 1.211] | 1.154 [1.055, 1.262] | 1.012 [0.904, 1.133] | 5.88 → 5.28 | 1.43 → 1.20 | 1.075 [1.005, 1.149] |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 8 | 10 | **1.293 [1.270, 1.317]** | 1.339 [1.294, 1.387] | 1.522 [1.486, 1.558] | 0.986 [0.957, 1.016] | 8.74 → 6.72 | 1.80 → 1.18 | 1.017 [0.951, 1.088] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 8 | 10 | **1.238 [1.215, 1.262]** | 1.250 [1.208, 1.293] | 1.371 [1.342, 1.402] | 0.985 [0.955, 1.016] | 10.12 → 8.13 | 2.55 → 1.87 | 0.990 [0.934, 1.049] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | 5 | **1.303 [1.281, 1.324]** | 1.294 [1.238, 1.353] | 1.356 [1.333, 1.380] | 1.005 [0.968, 1.044] | 16.74 → 12.77 | 2.56 → 1.88 | 0.986 [0.926, 1.049] |

**What the table says.**
- **At the three target sizes the filter is the gain.**
  - Collection runs 1.522, 1.371 and 1.356 times faster.
  - A scanned summand falls from 1.80, 2.55 and 2.56 units to 1.18,
    1.87 and 1.88.
  - Cold time is 1.293, 1.238 and 1.303 times faster.
  - `S` falls from 8.74, 10.12 and 16.74 to 6.72, 8.13 and 12.77.
- **The online interval gains too,** 1.250–1.339 at the target sizes.
  Its descent probes the same table and filter. Here it is a stage
  diagnostic; the rule's comparison will measure it against rho.
- **Rho's online interval does not move:** 0.985–1.005 at the target
  sizes. The patch does not change its code.
- **From `2^24.8` to `2^39.0` cold time reads 1.012–1.139.** Collection
  there runs 1.078–1.293 times faster. Every cold interval but one
  contains 1. The exception is `2^38.0`'s, 1.139 [1.029, 1.261].
- **The two smallest sizes read 0.939 [0.857, 1.028] and 0.933
  [0.861, 1.011].**
  - Both intervals reach into their A/A bands, [0.930, 1.118] and
    [0.984, 1.059], so neither is a regression by the declared test.
  - Collection is 6–11% of cold time there in v0's profile, and it did
    not move: 1.003 and 0.992.
  - Rho's online interval, whose code the patch does not change, reads
    0.936 and 0.948 in the same pairs. The shortfall is as large in code
    the patch did not change, so it does not point at the filter.
- **v1 reads slower in these runs than in R03's.** Its `S` at the six
  top sizes is 1–11% above its `S` in R03's runs: for example 8.74
  against 8.16 at `2^44.3`, and 16.74 against 15.61 at `2^47.2`. Runs
  hours apart on this host differ by several per cent, which is why the
  ledger compares rounds only by paired ratios.

## Against the prediction

The model is [`prediction.json`](prediction.json)'s. The declared
predictions at the target sizes also drew on the exploration.

| curve | log₂ r | collection, model | collection stage, measured | cold, predicted | cold, suite | cold, holdouts |
|:--|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 1.726 | 1.003 [0.891, 1.130] | 1.028 | 0.939 [0.857, 1.028] | — |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 1.422 | 0.992 [0.936, 1.051] | 1.033 | 0.933 [0.861, 1.011] | — |
| `icv1-f2m45-tm6236725-40939294` | 24.8 | 1.407 | 1.111 [1.003, 1.232] | 1.017 | 1.012 [0.964, 1.063] | — |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 1.111 | 1.078 [0.995, 1.169] | 1.029 | 1.030 [0.946, 1.121] | — |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 1.119 | 1.095 [1.046, 1.146] | 1.058 | 1.023 [0.978, 1.070] | — |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 1.207 | 1.208 [1.050, 1.390] | 1.075 | 1.055 [0.923, 1.205] | — |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 1.283 | 1.293 [1.154, 1.448] | 1.062 | 1.139 [1.029, 1.261] | — |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 1.208 | 1.154 [1.055, 1.262] | 1.133 | 1.074 [0.973, 1.184] | — |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.551 | 1.522 [1.486, 1.558] | **1.27** (model 1.307) | 1.293 [1.270, 1.317] | 1.276 [1.242, 1.311] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.385 | 1.371 [1.342, 1.402] | **1.25** (model 1.25) | 1.238 [1.215, 1.262] | 1.240 [1.208, 1.273] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.397 | 1.356 [1.333, 1.380] | **1.37** (model 1.336) | 1.303 [1.281, 1.324] | 1.319 [1.298, 1.340] |

**From `2^27.8` up the model held.**
- Its collection ratios are within 5% of the measured stage at every
  size there, and within 3% at the target sizes.
- Cold time came in 1–3% under the model at the target sizes. The model
  leaves out the larger filter's own cost, as the protocol said.
- The declared 1.37× at `2^47.2` leaned on the exploration's 1.374,
  which ran on Track B's stack, not on v1. The suite read 1.303 and the
  holdouts 1.319.

**At the three smallest sizes the model missed, and R04's counts suggest
why.**
- The model kept true hits, but it priced every admitted key at the
  admitted stage's mean cost.
- True hits are 0.10, 0.019 and 0.011 a summand at `2^18.0`, `2^22.0`
  and `2^24.8` (R04's `filter_check.json`), against under 0.0007 from
  `2^27.8` up.
- Where hits are that common, a hit's work dominates the stage, and the
  false positives the filter now rejects cost little.
- That is a reading of R04's counts, not a measurement.
- Collection is 6–11% of cold time at those sizes, so the miss moves
  cold time by a few per cent at most.

## The holdouts

**Fresh rows, run five rounds each.**
- The frozen files in [`holdouts/`](holdouts/) (`SHA256SUMS`), at the
  three target sizes.
- Recipe seeds 210–213, two targets each, `T111`–`T118`.
- None was extended: their half-widths after five rounds were 1.6–2.7%
  (the analysis's `extension_test`).

| curve | log₂ r | rows | rounds | pairs | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 8 | 5 | 40 | 1.276 [1.242, 1.311] | 1.333 [1.266, 1.403] | 1.499 [1.452, 1.547] | 1.010 [0.960, 1.063] | 8.10 → 6.63 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 8 | 5 | 40 | 1.240 [1.208, 1.273] | 1.231 [1.165, 1.300] | 1.365 [1.327, 1.405] | 0.997 [0.956, 1.039] | 10.03 → 8.07 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | 5 | 40 | 1.319 [1.298, 1.340] | 1.324 [1.281, 1.367] | 1.377 [1.354, 1.401] | 0.973 [0.937, 1.010] | 15.39 → 11.58 |

**The holdouts agree with the suite at every target size:** 1.276
against 1.293, 1.240 against 1.238, and 1.319 against 1.303.

## The extension

**After five rounds, two suite sets were wider than the declared 3%.**
- Their half-widths were 3.15% at `icv1-f2m53-tm56619371-dac20a85` and
  3.37% at `icv1-f2m59-tm943548413-98844ecc`.
- At `icv1-f2m61-t158598901-ab42b6c5` it was 1.68%.

**So those two got rounds 6–10,** and their figures pool all ten.
- `extended.json` in the run tree records the decision.
- The analysis re-tests the rule from the runs, and checks the record
  against the test.
- With five rounds they read 1.290 [1.251, 1.331] and 1.215
  [1.176, 1.256].
- With ten they read 1.293 [1.270, 1.317] and 1.238 [1.215, 1.262].

## Files

- **`runs.tar.xz`** (`runs.tar.xz.sha256`): the run tree.
  - It holds every process's report, isolation record and stderr, the
    pin, both host manifests, `extended.json`, and the refusal and PSI
    logs.
  - It reproduces `analysis.json` byte for byte; `tests/icprog.rs`
    checks this. From the repository root:

    ```
    tar -xJf research/ic_tool_program/rounds/R05-presence-filter/runs.tar.xz \
      -C research/ic_tool_program/rounds/R05-presence-filter
    cargo run --release --bin icprog -- analyse r05 > analysis.json
    ```
- **`analysis.json`**: `icprog analyse r05`'s output from the archive.
  `analyse.py`, written with the protocol, was not run: under AGENTS.md's
  rule the analysis is native (plan §10a).
- **The tables above** are `icprog table r05 analysis.json`'s output.
  The prediction table's measured columns are the same figures.
- **`tests.log`**: the candidate tree's tests, run before the timed
  steps.
- **`resume.sh` and `resume.log`**: the commands that resumed the run
  tree on the native runner after the rebuild, as run (the paths are the
  session's), and their console output, headed by both tools' hashes.
- **`baseline.meta.json`**: v2's identity. `icprog baseline` appended
  v2's entry to `baselines.json` from it, from the analysis and from the
  run tree's host manifests.
- **`candidate.patch`** (`candidate.patch.sha256`): the candidate, as
  frozen in the declaration. This pull request applies it.
