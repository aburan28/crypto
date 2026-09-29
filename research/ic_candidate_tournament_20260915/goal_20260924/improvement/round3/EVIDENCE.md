# Round three: provenance and interpretation

The third and final registered bounded attempt retained the qualified `pairinv`
incumbent. Development and selection chose `stop7_word`; it did not meet the
required 0.8 complete cold-instruction or cold-native-time ratio in either
confirmation or replay. Its online point estimate improved against the
separately qualified `ic_online` reference, but the familywise upper bounds
still included regression. This is a completed negative improvement round on
synthetic toy instances, not a production-curve or global-optimum result.

[The generated tables](README.md) list every variant in separate online,
cold-instruction and cold-native units. [RESULTS.json](RESULTS.json) retains
full precision, stage diagnostics, uncertainty, actual usable factor-base
sizes, folded columns, ranks and the frozen decision. [RUNS.csv](RUNS.csv)
retains all 3,480 canonical run keys and exclusive phase ledgers. All scheduled
pairs were verified; no failure, timeout or OOM row occurred in this attempt.

## Frozen execution

[PR 904](https://github.com/aburan28/crypto/pull/904) merged the
[third-round registration](../../improvement-v2/ROUND3-REGISTRATION.md),
[eleven-arm panel](../../improvement-v2/round3.json), portfolio correction,
runner and one-shot workflow at `f672329103dc2bffaecf95316b6354dbfde185c3`.
The panel's byte SHA-256 is
`6c87dafb98945bb49038f9389dc37881536700a9b9709d3862f9dbf79782763c`;
the exact candidate source-manifest SHA-256 is
`8582e4ab4b63e98696a0ff00ee296e2902923a2c39c325ab0f3e3950ffbb2b28`.
The main-branch [workflow run](https://github.com/aburan28/crypto/actions/runs/36463687634)
executed the registered panel once with `GITHUB_RUN_ATTEMPT=1` and seed
`2026092553`. Its fixed-vector controls passed 20/20 before the measurement;
the runner's frozen verifier completed successfully on Linux.

The contract used Linux x86-64, rustc 1.94.1, a musl static worker,
Valgrind 3.22.0, one pinned CPU and worker thread, an 8 GiB address-space cap,
and a 180-second child timeout. The retained host record identifies the
kernel and runner node but does not establish a specific CPU model.
All three rounds have distinct seeds and target exclusions. This attempt
restored and checked both prior decisions before preparing any new point,
then excluded their exposed targets and seven additional hash-pinned fixture
corpora. Its preflight verified 5,697 distinct points in the initial history
and supplemental corpora before adding the two prior rounds' exposures.

The primary workload supplied one previously unseen public point per case.
Reusable IC and rho preparation and fixture construction were outside the
online interval; target-dependent work through scalar replay was inside it.
Confirmation used 72 fresh points, 12 in each of six curve/subgroup cells,
with three process repetitions per point. Replay used the same 72 points in
fresh processes. Repetitions were not counted as additional targets.

| Stage | Native/profile pairs | Verified |
| --- | ---: | ---: |
| A/A | 30 | 30 |
| Smoke | 210 | 210 |
| Development | 630 | 630 |
| Selection | 450 | 450 |
| Confirmation | 1,080 | 1,080 |
| Replay | 1,080 | 1,080 |
| Total | 3,480 | 3,480 |

## Decision and comparison scope

Cold instruction and native-time candidate ratios use the qualified `pairinv`
incumbent. Online candidate ratios use the separately qualified `ic_online`
reference. The matched `rho_online` baseline solves the same single public
point; online rho / IC is computed directly from measured online costs, not
by dividing ratios with different IC denominators. Cold rho is a separate
qualified reference. Each final ratio is an equal-cell geometric aggregate of
per-point three-process medians. The nominal familywise one-sided rule spans
three attempts, both final stages and three metrics; its bootstrap coverage is
approximate, and replay repeats points in new processes.

| Selected `stop7_word` / qualified IC reference | Confirmation ratio | Familywise upper | Replay ratio | Familywise upper | Promotion gate |
| --- | ---: | ---: | ---: | ---: | --- |
| Single-target online time | 0.973335 | 1.033892 | 0.978750 | 1.040371 | Ratio ≤ 1; upper < 1 |
| Complete cold instructions | 0.990182 | 1.034953 | 0.990183 | 1.034954 | Ratio ≤ 0.8; upper < 1 |
| Complete cold native time | 0.985697 | 1.020563 | 0.980580 | 1.014370 | Ratio ≤ 0.8; upper < 1 |

The largest replay online cell ratio was 1.086595, below the 1.1 per-cell
limit, but all three familywise upper bounds exceeded one. Neither cold point
estimate approached 0.8. The frozen checker therefore retained the incumbent
with `promotion_eligible: false`; confirmation and replay were not used to
retune the panel. The incumbent's confirmation online cost was approximately
0.033926 ms against 0.252643 ms for the qualified online rho reference on the
matched toy panel, a descriptive rho / incumbent ratio of 7.446838. These
values do not establish a larger-curve crossover. On n31a0 the selected arm
actually admitted 434 subgroup-usable factor-base points and seven folded
columns, compared with 496 and eight for the incumbent; the labels alone are
not used as support-size evidence.

The reference's collision-table allocation and factor-base-only memory remain
unmeasured. Whole-process peak RSS is retained separately. F4, F5, SAT and
large sparse relation-LA backends were not executed in this specialized
pair-table panel. The broader generic adapters and correctness controls need
their own complete, source-bound comparative qualification before promotion.
The K-instruction floor in the generated report is specific to this full-rank
collector and is not a generic IC lower bound.

## Durable archive and independent replay

The [archive](../../../evidence/ic-improvement-round3-20260928.tar.zst) and
[manifest](../../../evidence/manifest.json) bind the raw sources, workers,
fixtures, admissions, native/profile outputs, mathematical certificates,
profiles, stage summaries, decision, command logs and all 20 candidate-policy
controls:

- SHA-256: `a2e88ae540ddfa3430b03873f5cd85df0d5b1ad12b8e59b7dea82a5b33df0a4f`.
- Compressed bytes: 37,223,911.
- Restored regular files: 361,985; original file bytes: 514,447,908.

Build caches, Python bytecode, operation locks and the restored prior/reference
input trees are excluded. The profiles reconstruct to the original gzip bytes
using recorded headers, trailers and SHA-256 hashes. A first local compression
attempt was interrupted by a transient full SSD; its incomplete file was
discarded before this archive and no manifest entry was made. No candidate
worker, fixture, timer, or measurement was rerun. [Fresh archive extraction](ARCHIVE_RESTORE.json)
passed the archive hash, every reconstructed-profile hash, file count and byte
count. The measured Linux frozen checker had already verified all 3,480 pairs.
Exact cross-platform floating-summary replay is Linux/Python 3.12-specific;
a local macOS replay rejected a stage-summary/provisional-selection value.
The evidence PR's Linux/Python 3.12 archive audit, byte-exact table export and
scoreboard check are the independent transport gate; a successful replay does
not convert the rejected challenger into a win.

From the repository root on Linux with Python 3.12 and zstd, choose an empty
output directory and run:

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-improvement-round3-20260928 --out /tmp/ic-round3-replay
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/audit_archive.py \
  --bundle /tmp/ic-round3-replay/ic-improvement-round3 \
  --out /tmp/ic-round3-audit.json
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/export_report.py \
  --round /tmp/ic-round3-replay/ic-improvement-round3/round/tournament \
  --out research/ic_candidate_tournament_20260915/goal_20260924/improvement/round3 --verify
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/render_report.py \
  --results research/ic_candidate_tournament_20260915/goal_20260924/improvement/round3/RESULTS.json \
  --markdown research/ic_candidate_tournament_20260915/goal_20260924/improvement/round3/README.md \
  --scoreboard docs/index-calculus-scoreboard.html --verify
```

These are read-only evidence replays. Do not redispatch any of the three
registered attempts. This result exhausts the three-round bounded goal without
a promoted challenger; future tests need a new predeclared protocol and fresh
targets.
