# PR1853: ordinary conflict integration

This is CI/review follow-up for the delivered research PR, with no campaign
rerun, target resampling, parameter search or new session. The previous
published head is 6351a6a3df065c5ffc7b81aba5d4a82157b022e6. Actual remote main
is f2c9202ca48ed16f086815fda9c0ab73f2aabcda. The PR API still reports the old
a80 base while correctly reporting a dirty merge state.

The fresh Copilot request6186854498 received no response or pushed resolution
after about 65 minutes. Per current AGENTS.md's reasonable-wait fallback,
[comment6187909093](https://github.com/aburan28/crypto/pull/1853#issuecomment-6187909093)
announces local integration and ends the pending local wait. Main is merged into
the branch normally, without rebasing or force-pushing. Conductor is unavailable
(exit127); conductor.log records the unavailable command.

## Resolution and preservation

The only textual conflict is benchmarks/README.md. Retain every row from both
parents verbatim. No frozen evidence, ledger row or claim conflicts. The
ecbench methods and module files merge automatically: incoming field support
coexists with ic.orbit_rank_v1. combined-dispatch.patch.gz records the result
against main; benchmark-resolution.patch.gz records the index resolution.
Incoming library implementations remain their exact main-parent objects;
incoming-groebner.patch.gz documents that imported engine change.

preserve.mjs compares repository object IDs, without executing research code.
Its passing receipt checks all 1104 prior campaign/study/benchmark artifacts
and 84188 incoming evidence files against their respective parent commits.
Every benchmark index row is retained and the original measured follow-up
binary still hashes to
50a5739fadcb08ab8f170d4fe03ba3427fb0fd387c1a5918843fcda4d68ff1fb.
The initial local Node process-permission failure is recorded separately;
the unchanged script passes in the network-enabled execution context.

The completed 6351a6a3d catalog/publication repair remains intact. Main changes
no published browser/leaderboard/catalog inputs, so the exact historical
metadata baseline remains a80. The existing strict native freshness check
passes; no derived view needs regeneration in this integration. No
assertion, verifier, workflow or timeout is disabled.

## Validation

validate.sh records exact commands, complete deterministic-gzip logs and exit
statuses. It runs only integration tests, compilation, lint and publication
freshness. The heavy validation uses the existing isolated_bench.py busy lock.
CLI tests execute their temporary smoke fixtures; no research session, campaign
or target resampling is performed. Full-library validation is
justified by the incoming library/dispatcher changes. Prior repair validation
is retained in ../pr-followup-20261015T0647Z/.

The single full release library invocation passes 4116 tests, with zero
failures and 117 existing ignored tests. Its complete output and exit0 are
retained in library.log.gz and library.exit.

| Check | Result | Receipt |
|---|---|---|
| Release library | 4116 pass, zero fail, 117 existing ignored | library.log.gz |
| CLI/integration and six original examples | 27 pass | cli-examples.log.gz |
| Incoming touched example | 7 pass | incoming-example.log.gz |
| All-target release compilation | pass | all-targets.log.gz |
| Warnings-as-errors clippy | pass | clippy.log.gz |
| Strict locked native metadata freshness | pass | native-freshness.log.gz |
| Site tests | 48 pass | site.log.gz |
| Workflow parser | 4 pass | workflow-parser.log.gz |
| Whole-repository Rust1.98 format | pass | format.log |

The ecbench binary, ecbench_leads_findings and ecbench_open_questions_run have
zero unit tests and are compiled targets. The other four original examples
contribute six passing tests; integration targets contribute 21. No zero-test
target is counted as a test case. validate.sh's introductory comment is clarified
after execution to distinguish smoke fixtures from scientific sessions; its
commands are unchanged.

Whole-repository Rust1.98 formatting now passes: main already formatted the
two examples that failed on a80. Their old failed logs remain unchanged.
The leaf-m10 archive and its actual workflow are byte-identical between a80
and this incoming main (archive-workflow-comparison.exit0). The first
archive-base-comparison established archive equality; the second comparison
names the actual outcome workflow. The previously
reproduced exact-file-set failure remains a merge gate; no archive, manifest,
verifier or workflow is repaired here. See the existing
[baseline-CI problem](../../../problems/ecbench-ten-question-baseline-ci.md).

All recorded integration validation passes and both-parent diffs are reviewed.
Every unsuccessful check remains recorded. Original 15 audits, 1203 exact
verified recoveries and 107 original-binary failure replays remain immutable;
this integration does not claim new audit executions or independent replay.
No graphable research result changed, so original graphs/PDFs and benchmark
copies stay unchanged.

## Remote gates

gates-before-merge.json records the single scheduled inspection of head635:
only three service checks exist, including neutral failed-to-run Bugbot; the
successful review/approval automation jobs do not supply review approval.
The latest Cursor submission at 07:24:50Z explicitly says Not approved and
requires human review. There are zero inline threads and no approval.
The new merge commit needs its own applicable checks; old greens do not
transfer. Leave the PR open until all required/applicable checks and review
permit merge, with no blocking conflict or feedback.

SHA256SUMS covers every file in this evidence directory except itself. Raw
commands and patches are deterministic gzip files; no whitespace in original
logs is normalized. The captured first merge returns1 because of the ordinary
index conflict, and the initial process-permission check returns1; both failed
statuses are retained separately from successful validation.
