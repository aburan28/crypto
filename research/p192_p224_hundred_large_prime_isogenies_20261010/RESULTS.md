# One hundred large prime-degree isogeny maps

The campaign independently verified 100 distinct maps: 50 consecutive edges in each public P192/P224 isogeny class. There are two distinct prime degrees, 2341 and 1117. The first maps start on the published sources; later maps start on their independently verified codomains.

The [frozen protocol](PROTOCOL.md) defines map-level uniqueness and retains one-eigenline coverage PARTIAL at every step. Existing P192/10453 and P224/1471 maps and repeated replays are excluded. This meets a hundred-map target with repeated prime degrees allowed; a hundred different prime degrees remains a separate task.

| Public class | Prime degree | New maps | Constructor seconds | Independent replay seconds | Summed invocation seconds |
|---|---:|---:|---:|---:|---:|
| p192 | 2341 | 50 | 257.592 | 2466.415 | 2724.952 |
| p224 | 1117 | 50 | 124.170 | 667.004 | 792.191 |

These are complete operational intervals on a shared Apple M4 Pro host, with two concurrent walks. Constructors have an 1800-second limit and 8 GiB RSS cap. Independent replay is outside that constructor cap. These measurements do not establish a matched speedup. All 100 exact kernel, codomain and rational curve-equation checks passed, alongside subgroup-order checks and 20 public scalar-transport checks per map (2000 total). Replay uses separate field and polynomial implementations on the same host.

[SUMMARY.json](SUMMARY.json) lists all exact source/destination ICV1 and EC1 records, degrees and map hashes. [P192 walk](evidence/p192/walk.json) and [P224 walk](evidence/p224/walk.json) link every raw kernel and rational map, resource and stderr receipt, mapped public generator and independent replay. The embedded source revision was clean; both walks retain the binary digest.

![Verified routes and coverage](SEARCH.svg)

The screen covered every prime from 1010 through 1000000 on both sources. Orders at most 8 selected 2 P192 and 3 P224 candidates; orders at most 64 selected 40 and 38. These are structural candidates. The generalized constructor supports checked extension degrees 2 through 16. All enabled degrees passed field-law tests, and the original degree-3/5 selections remain unchanged.

## Validation and delivery

The CLI passed 62 tests plus 3 limb-field reference tests. Standalone algorithm suites passed 114 tests. The initial required root compilation failed with 643 errors. Three existing upstream type-repair commits restored compilation. The first complete runtime attempt exposed catalogue parser gaps, a serde-feature-dependent fixture, one sparse fixture omission and a one MiB digest buffer overflowing the default test stack. Their original diagnostics are retained, and the parser, fixture reference and buffer allocation were repaired without changing map acceptance checks. An additional strict CLI Clippy diagnostic retains 35 errors from inherited/generated modules; that check is not presented as passing. The required algorithm Clippy job is checked separately. PR #1606 remains open and the release remains unpublished.

Canonical updates append the 100 verified models and mapped generators, rebuild covers and links, and refresh aliases, leaderboard rosters and browser joins. The complete route graph is in [ROUTES.mmd](ROUTES.mmd). The IC boundary ledger, performance scoreboard and timeline retain their existing measurements because this round does not measure an IC improvement. Earlier frozen studies retain their historical findings.

Run the actual executable with: isogeny walk --curve p224 --ell 1117 --count 50 --out ./p224-walk. Replay it with: isogeny verify-walk --input ./p224-walk. Both arithmetic implementations and pinned public anchors are embedded. Runtime Cargo and a checkout are unnecessary. See [tool instructions](../../tools/isogeny-cli/README.md), [commands](COMMANDS.md) and [manifest](MANIFEST.json).

The final required library suite is running. Push/publication remains gated on its result.

The unsandboxed full library run has at least 10 observed failures outside the dedicated map verification path. Its retained [log prefix](evidence/validation/root-complete-prefix.log) is a checkpoint, not a completed passing suite. PR push is awaiting the owner's explicit choice to continue broader repairs or publish with merge/release still blocked.

The rebuilt 449-model executable passed 62 CLI tests and 3 independent limb-field reference tests. It also passed a two-step saved-chain replay and complete two-map P384/19 construction/replay from a directory without a checkout. The [Apple silicon executable archive](evidence/isogeny-aarch64-apple-darwin-3c39c5f23.tar.gz) includes BUILD.json and source commit 3c39c5f23e0bd1590d220fb32386ad9ed79bf399. Its source state was clean.
