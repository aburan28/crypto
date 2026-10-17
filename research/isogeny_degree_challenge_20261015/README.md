# Prime Degree Frontier: native challenge and admission evidence

## Current release pointer — 2026-10-16

The dedicated tool is now published as
[v0.1.2390](https://github.com/aburan28/crypto/releases/tag/v0.1.2390),
source `bd7eb117f094a336d3e3e8f31facb0c6b86f95bc`.
Use the [download and installation instructions](../../docs/challenges/isogeny-prime-degree/README.md#run-the-dedicated-executable)
and [verified publication receipts](../isogeny_world_frontier_20261016/replay-optimization-v1/official-release-v0.1.2390/README.md).

The packet below is the preserved source-`9bb7cbce2` checkpoint. Its season,
test totals and publication-pending row describe that historical packet.
The newer release has its own version/policy receipts; its source is not
substituted into those frozen runs. Newer campaign outcomes and representation
comparisons are linked in the [dated campaign checkpoint](../isogeny_world_frontier_20261016/record-comparison-20261016/README.md).

## Frozen implementation checkpoint

The dedicated `isogeny` executable now supplies a largest-prime-degree challenge
for **all 54 registered standard prime models**, ranked separately by exact ICV1.
The standard roster uses nonempty `standard_names`; registered research descendants
and benchmark models remain outside this track. The pinned catalogue has 47
construction/replay adapters in this roster. Seven models require checked coordinate
conversions and stay visible without a score. Three source models were exercised
in the final starting leaderboard; adapter availability is not a claim that all
47 have been run in this packet.

The native commands are `challenge policy`, `plan`, `run`, `judge` and `leaderboard`.
Each judgement snapshots bounded raw files and runs a fresh independent verifier.
A complete explicit normalized map and monic kernel are required. Equal prime
degrees share rank within a source curve. The judge verifies mathematical maps
rather than trusting a stored PASS or score. Every leaderboard build replays all
of its submitted raw artifacts again.

The final corrected executable is frozen at source
`9bb7cbce2b84151eff6846b02da842a8ef019343`. The earlier ff71f5772 packet
exposed an overstated COMPLETE label when a valid kernel was duplicated. The
corrected packet rejects that ninth adversarial case and retains the earlier
failure separately. The season is
`19f4503758a3dc2d9764f2283b4c92c7a35651e75ee7157f546b2563fcbf5c30`.
The frozen experimental executable is archived as
`store:evidence/isogeny-degree-challenge-20261015/isogeny-9bb7cbce2-macos-arm64`.
It is an evidence artifact. Official binary publication uses the existing ReleaseMe
workflow after the PR lands and its release checks pass.

## Results and requirements

| Requirement | Status and evidence |
|---|---|
| All registered standard prime models ranked separately | Verified: 54 exact source rows in the policy, plan and starting board; canonical roster cross-checked against the full embedded registry |
| Working native executable | Verified: installed-only directory, 74 package tests, retained command/exit/byte records; no Cargo or source checkout required at runtime |
| Independent map admission | Verified: P192/13 and P521/7 positives, all nine adversarial cases rejected, forged stored scores rejected |
| Resource lifecycle | Verified in final packet: 20/20 fresh short-map judgements accepted with positive RSS samples; one-second high-degree control censored with no score |
| Higher-degree candidate ladder | Verified as a screen: 157 curve-degree pairs from 10454 through 65537, order filter at most 16, descending per source; candidates earn no score |
| Historical starting scores | Freshly verified: P192/10453 and P224/1471 each have one eigenline, coverage PARTIAL; P521/7 has both eigenlines, coverage COMPLETE |
| Construction on every standard model | Partial: seven coordinate-conversion adapters remain unavailable; only three models exercised here |
| Root library release gate | Verified: 4121 passed, zero failed, 117 existing ignored; complete terminal receipt retained |
| PR and published release | PR publication follows the local gate; actual new ReleaseMe asset publication remains pending |

The [canonical starting board](../../docs/challenges/isogeny-prime-degree/leaderboard.json)
contains all 54 rows, three with verified starting entries and 51 unscored.
These are re-admissions of historical public maps. **New discoveries in this
packet: 0.** The earlier additional-map and CPU optimization campaigns remain
separate unfinished workstreams.

The highest screened candidates for P192 and P224 are 50341 and 28537 respectively.
They are construction candidates, not new scores. The source/target ICV1 identities,
map hashes, independent certificates and twenty public scalar transports per map
are linked in the [final native observations](reports/run-20261016T0029Z/observations.md).

## Retained runs

| Run | Observed outcome | Records and visual report |
|---|---|---|
| 20261015T2317Z, source 7e96187d3 | Installed controls passed; the separate full suite exposed a valid short replay with no RSS sample | [Observations](reports/run-20261015T2317Z/observations.md), [raw](reports/run-20261015T2317Z/raw.tar.gz), [PDF](reports/run-20261015T2317Z/report.pdf) |
| 20261015T2325Z, source 66b5b7ce0 | First installed positive censored after startup; later controls not reached in this run | [Raw](reports/run-20261015T2325Z/raw.tar.gz), [failed-control PDF](reports/run-20261015T2325Z/report.pdf) |
| 20261015T2333Z, source d72229b34 | Full controls passed; resource stress accepted 19/20 and censored one during OS exit cleanup | [Observations](reports/run-20261015T2333Z/observations.md), [raw](reports/run-20261015T2333Z/raw.tar.gz), [PDF](reports/run-20261015T2333Z/report.pdf) |
| 20261015T2343Z, source ff71f5772 | Full controls and 20/20 stress passed; short budget censored; all three starting entries freshly accepted | [Observations](reports/run-20261015T2343Z/observations.md), [raw](reports/run-20261015T2343Z/raw.tar.gz), [PDF](reports/run-20261015T2343Z/report.pdf) |
| 20261016T0022Z, baseline ff71f5772 | Additional duplicate-kernel regression failed: valid degree 13, unsupported COMPLETE label accepted | [Observation](reports/run-20261016T0022Z-duplicate/observations.md), [raw](reports/run-20261016T0022Z-duplicate/raw.tar.gz), [PDF](reports/run-20261016T0022Z-duplicate/report.pdf) |
| 20261016T0029Z, corrected source 9bb7cbce2 | All nine negative cases, 20/20 stress and three fresh starting admissions passed; short budget censored | [Observations](reports/run-20261016T0029Z/observations.md), [raw](reports/run-20261016T0029Z/raw.tar.gz), [PDF](reports/run-20261016T0029Z/report.pdf) |

Replay now waits for a positive RSS sample before parsing. It writes a completion
marker only after parsing, algebra and failure handling finish. RSS polling ends
at that marker; the process deadline continues through output and operating-system
cleanup. Missing telemetry during active work still fails closed. Earlier failures
and the 19/20 intermediate outcome stay in their own packets.

Limits are separate: constructor 1800 seconds and 8 GiB sampled RSS; replay at most
3600 seconds and 8 GiB; input map 256 MiB, summary 1 MiB, stderr 8 MiB. The bundled
constructor stops at degree 1000000; larger explicit u64-prime submissions may be
judged within the artifact/replay limits. External constructor costs are retained
but not independently attested, and are not used as a degree tiebreak. Shared-host
operational receipts establish no isolated timing comparison.

## Run and reproduce

Install the dedicated ReleaseMe `isogeny` executable once the new command has been
published. For a locally retained experimental binary, use the source/hash named
in the packet's `freeze.json` and `custody.json`.

```sh
isogeny challenge policy > policy.json
isogeny challenge plan --curve all --from 10454 --to 65537 > ladder.json
isogeny challenge run --curve p192 --ell 13 --out fresh-control
isogeny challenge run --curve p192 --ell 50341 --out fresh-next-rung
isogeny challenge judge --input submission --out fresh-judgement
isogeny challenge leaderboard --input submissions.json --out fresh-board
```

To replay an archived starting entry, extract its `raw.tar.gz` into a fresh
workspace and judge `anchor-board/entry-0/artifact` (P192), `entry-1/artifact`
(P224) or `entry-2/artifact` (P521). Original command paths stay in receipts;
replay resolves the extracted raw directory directly. Never submit a saved
`entry.json` or `leaderboard.json` as a proof. The one-second control and full-budget
judgements use separate directories.

Maintainer build commands and the embedded rules are documented in
[the tool README](../../tools/isogeny-cli/README.md) and
[the challenge contract](../../docs/challenges/isogeny-prime-degree/README.md).
The native [audit.rs](audit.rs) reproduces counts, source/binary binding and report
observations. [report.typ](report.typ) and [failed-report.typ](failed-report.typ)
render from the committed JSON and diagrams. `store_cli.rs` mounts the repository's
unchanged native artifact-store implementation for archival; it is not a companion
program required by the installed `isogeny` tool.

The IC/rho dashboard, progress timeline, IC leaderboard and canonical curve-edge
charts were checked and remain unchanged: this work adds challenge admission and
replays already registered maps, with no new IC ratio, registry model or curve edge.
