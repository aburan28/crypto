# Round two: provenance and interpretation

The registered version-two tournament retained the qualified `pairinv` incumbent.
The selected `stop5_word` challenger improved complete cold instructions only
modestly (confirmation ratio 0.922725) and did not meet the 0.8 cold gates; its
single-target online time against the separately qualified `ic_online` reference
regressed (confirmation ratio 1.093655). Replay agreed. This is a completed
negative improvement round on synthetic toy instances. It is not completion of
the three-round goal, a production-curve result, or a global-optimum claim.

[The generated tables](README.md) contain all variants, separately for online
native time, cold instructions and cold native time. [RESULTS.json](RESULTS.json)
retains full precision, phase diagnostics, per-cell counts and the exact
decision. [RUNS.csv](RUNS.csv) retains all 3,480 canonical run keys and their
exclusive ledgers. No scheduled pair failed or timed out in this round; the
exporter still preserves those statuses when present.

## Reporting correction

The original aggregate exporter divided rho/incumbent by challenger/ic_online.
Those denominators differ, so that quotient was not the matched rho/challenger
speedup. The corrected report divides the measured online costs directly and
names each row's IC denominator. It also includes the auxiliary IC reference
comparisons already retained by the frozen evaluator. Historical values remain
in [the complete before/after record](REPORTING-CORRECTION.json).

| Selected challenger stage | Previously exported online rho / challenger | Corrected matched ratio |
| --- | ---: | ---: |
| Confirmation | 5.058918983932264 | 5.0219292766882875 |
| Replay | 5.1760622769942755 | 5.130908564435435 |

This is an accounting correction. The raw archive, per-target paired rows,
measured costs, RUNS.csv, frozen evaluator, confirmation/replay gates and
negative promotion decision are unchanged. No worker was rerun. Linux evidence
CI must independently reproduce the correction from the retained archive.

## Frozen execution

[PR 892](https://github.com/aburan28/crypto/pull/892) merged the
[version-two protocol](../../improvement-v2/PROTOCOL.md),
[eleven-arm registry](../../improvement-v2/round2.json), candidate patch, runner
and workflow at `9a34071721bdfaea8c64d281f3412482a49b238d`. The measured head was
that same merge; the registered candidate source SHA-256 is
`8582e4ab4b63e98696a0ff00ee296e2902923a2c39c325ab0f3e3950ffbb2b28`. GitHub
`workflow_dispatch` with `run_round_two=true` was denied (HTTP 403) in this
environment, so the calibrated measurement ran once locally under the same
contract the workflow would have used: Linux amd64, rustc 1.94.1, the musl
static worker, Valgrind 3.22.0, seed `2026092552`, attempt 2, and
`run_improvement_v2.py`. No measured process, source, fixture or promotion rule
was retried or changed after the first complete prepare/run/verify.

Workers used one thread pinned to a single CPU, an 8 GiB address-space cap and a
180-second child timeout. The schedule bound was 3,480 native/profile pairs. The
recorded host/kernel is in `round/tournament/contract.json`; an exact CPU model
was not captured, so it is unknown.

Five development cells were `n17a1`, `n19a0`, `n23a0`, `n23a1`, `n31a0`; the
final panel added holdout `n29a1`. The `n` label is the binary-field degree, not
subgroup bits. Confirmation used 12 fresh points per cell, 72 total, and three
processes per point. Replay repeated those same 72 points in fresh processes.
Exclusions covered the restored round-one prior plus the readiness,
generic-reference-qualification and adapter-control fixture corpora. The round
preserves every generated fixture and admission, including candidates that were
not measured in the final stages.

| Stage | Native/profile pairs | Verified |
|---|---:|---:|
| A/A | 30 | 30 |
| Smoke | 210 | 210 |
| Development | 630 | 630 |
| Selection | 450 | 450 |
| Confirmation | 1,080 | 1,080 |
| Replay | 1,080 | 1,080 |
| Total | 3,480 | 3,480 |

The preceding fixed-vector candidate controls independently verified 20/20
pairs (ten challengers × two archived toy vectors). Those controls do not supply
natural-query yield or performance selection.

## Decision and interpretation

The primary interval starts with the supplied point after reusable preparation
and ends after scalar replay. Fixture generation and reusable preparation are
excluded. Complete cold costs remain additional acceptance gates; three process
repetitions never become three independent targets. Online ratios use the
qualified `ic_online` reference; cold instruction and cold native ratios use the
qualified incumbent.

`stop5_word` confirmation ratios were 1.093655 for online time, 0.922725 for
cold Ir and 0.987022 for cold native time. Familywise one-sided upper bounds
were 1.185021, 0.967882 and 1.026999 respectively. Replay ratios were 1.082607,
0.922726 and 1.004137, with uppers 1.172905, 0.967883 and 1.040030. Largest
per-cell online ratios exceeded 1.1 on `n23a1` in both stages. Neither cold
point estimate reached 0.8. The frozen checker therefore retained the incumbent
without promotion. Classification is accounting: the round completed the
registered measurement and recorded a negative promotion outcome; it does not
move the floor ratio.

The incumbent's measured confirmation online time was 0.0300712 ms against
0.166375 ms for the separately qualified online rho reference. Cold instruction
point estimates continue to beat the cold rho reference on this toy panel;
`beats_rho_strict` is true for the retained incumbent under the decision's
paired cell definition, and remains a descriptive toy-panel observation, not a
production-curve crossover.

Development and selection are stage diagnostics only. Selection's provisional
challenger was `stop5_word`. Confirmation and replay were not used to retune
candidates. Any next-round design must use development evidence only and exclude
all prior generated points, including every fixture from this round.

The K-instruction floor is a weak implementation-specific full-rank collector
floor, not a generic IC bound. Changing the factor-base policy is explicit;
actual B/K values are retained. Base-only memory was not instrumented and remains
unknown; peak process RSS is separate. F4/F5/SAT and large sparse LA were not
executed in this round. The tested mechanisms are pair-table PDP policies
including half/cover modes, Frobenius-orbit base construction and scalar-field
Gaussian row kernels.

## Durable evidence and independent replay

The [archive](../../../evidence/ic-improvement-round2-20260928.tar.zst) is
committed with [its manifest](../../../evidence/manifest.json):

- SHA-256: `020e0cef5f044e01613625da93d9cedaa48e875a59523f2859477c484f23b6f0`.
- Compressed bytes: 36,473,587.
- Restored regular files: 352,644; restored file bytes: 468,890,707.

It retains the original sources, worker binaries, fixtures, admissions,
native/profile output, mathematical certificates, all raw profiles, stage
summaries, decision, command logs, local completion provenance and the 20
candidate-policy controls. Build caches, Python bytecode, operation locks and
the restored prior-round / reference-evidence trees used only as inputs are
excluded. Expanded profiles in the zstd container reconstruct to their original
gzip bytes, which the restorer hashes.

[Fresh archive extraction](ARCHIVE_RESTORE.json) passed the archive hash,
reconstructed-profile hashes, file count and byte count. The measured Linux
checker verified all 3,480 records and 1,878 files across the source snapshots,
plus 20/20 policy controls. [Linux transport replay](TRANSPORT_LINUX.json)
records that independent audit. A fresh Linux/Python 3.12 restore, full audit,
byte-exact table export and scoreboard check are required by the evidence PR's
CI. Successful verification does not turn the rejected challenger into a winner.

From the repository root on Linux, with Python 3.12 and zstd available, choose a
new output directory and run:

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-improvement-round2-20260928 --out /tmp/ic-round2-replay
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/audit_archive.py \
  --bundle /tmp/ic-round2-replay/ic-improvement-round2 \
  --out /tmp/ic-round2-audit.json
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/export_report.py \
  --round /tmp/ic-round2-replay/ic-improvement-round2/round/tournament \
  --out research/ic_candidate_tournament_20260915/goal_20260924/improvement/round2 --verify
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/render_report.py \
  --results research/ic_candidate_tournament_20260915/goal_20260924/improvement/round2/RESULTS.json \
  --markdown research/ic_candidate_tournament_20260915/goal_20260924/improvement/round2/README.md \
  --scoreboard docs/index-calculus-scoreboard.html --verify
```

These are read-only evidence replays, not new measurements. Do not dispatch the
second round again. Do not redispatch round one. The repository skill and goal
status point to this result and keep the remaining one-round budget open.
