# Fresh photographed-curve volcano campaign

**New run authorized after the historical artifacts became unavailable.** The
old experiment is not recovered or pooled. This replacement uses the same exact
photographed source curve with new securely sampled fixtures and a new frozen
identity. The original complete workload remains required.

## Original native corpus complete, 2026-10-17 18:21 UTC

All **4600 rho and 4600 BSGS cases** are verified: **100 problems on each of 46 curves**. An independent PARI replay verifies all **9200 recovered-log equations**. There are no failed final native cases.

[Complete native min/mean/max report](reports/native_complete_1823_v2_visuals/README.md) · [Full setup/query and retained-attempt ledger](reports/costs_complete_1821/README.md). The ledger retains **5062 attempts across 5060 logical jobs**.

For each paired key, BSGS online operation counts agree over every curve and generator; all 46 curve averages are **71837588.1 counted online group operations**. Table construction is separate. The separately frozen field-counter replay is still running. Full random-corpus IC recovery remains unresolved, so this completes the original generic-solver corpus, not the entire requested experiment. Earlier snapshots remain unchanged.

## Progress snapshot, 2026-10-17 18:00 UTC

All **4600 rho cases** and **4520/4600 BSGS cases** have verified answers.
The separately frozen native arithmetic replay is still running; its timings
do not replace the original measurements. All **46 curve cardinalities** and
**9200 point orders** match independent PARI checks. The full-path transport
phase verifies **23000 actual P/Q pair evaluations**. Full random-corpus IC
recovery remains unresolved; see [the precise gap](IC_RECOVERY_GAP.md).

New reports: [transport](reports/transport_1802/README.md),
[setup/query cost ledger](reports/costs_1801/README.md), and
[orders](reports/orders_1525_v2/README.md). Older dated snapshots below retain
their historical counts. The full experiment is not complete.

## Reconstructed and validated

- Complete selected rational 2-isogeny component: **46 curves, 90 directed maps,
  five levels with populations 1/3/6/12/24**. This is not the whole isogeny class
  over every prime degree. New C labels differ from historical labels.
- Exact source: p=72057980207847001, a=18456521102804593,
  b=19958176276690439; #E=72057979853368804, prime subgroup r=18014494963342201.
- Ten distinct generator points modulo sign and ten distinct nonzero uniform k
  values, shared across generators: 100 problems/curve and **4600 transported
  problem/curve cases**. Linux getrandom and unbiased rejection sampling;
  every raw word and accepted/rejected result is retained.
- PARI/GP independently counted points and group invariant factors on all 46
  curves and verified all 90 rational curve/dual identities symbolically.
  Native C++ independently checked primality, all 4600 scalar/order fixtures,
  9000 directed-map fixture pairs and 900 homomorphism/dual controls.
- Unmodified pinned native C backend tests passed: group 376415 checks, BSGS
  2063 checks and rho 134 checks. CMake and Cargo are unavailable; the library
  was built directly with GCC, and no whole-Rust-library test pass is claimed.

The first full-subgroup rho solve and first ten-target BSGS batch passed
independent scalar checks. Those are initial cases, not the full campaign.
Low-level field instrumentation, full-case solver completion and actual full
IC recovery remain open. No fastest-tool, speedup, hardness or complete-cost
claim is made. Local wall and process-CPU times are descriptive L0 records.

## Recoverable inputs and executable source

[artifacts/PUBLIC_MANIFEST.json](artifacts/PUBLIC_MANIFEST.json) identifies the
[public input/source archive](artifacts/public_inputs.tar.xz). It contains all
4600 answer-free P/Q fixtures, exact curve/map records, public rho walk seeds,
native source snapshots and protocol. Expected/solved scalars and raw getrandom
words are excluded. The complete local archive retains them; publication of
that complete archive requires explicit approval and has not occurred.

The fresh solver-fixture SHA-256 is:

`980d931ce1a804684b160bdab9b213b87688ab9100f0cdf1cd14c390e809acb8`

`sha256sum -c artifacts/public_inputs.tar.xz.sha256` checks the public archive.
Inspect its entries with `tar -tJf`, then extract into a fresh directory.
PUBLIC_MANIFEST.json inside lists exact file hashes and excluded private files.
The public solver input contains no expected scalar. Source is available for
review and reconstruction; the original worker's optional checker-data
preflight and raw-RNG audit cannot be restored from this public bundle alone.
No full-private-artifact restoration is claimed.

PARI runtime Debian URLs and hashes are retained locally for certificate replay.
For rebuilding the native library, use the preserved native snapshot as
upstream: `ln -s source/native_snapshot upstream`, then
`bash source/build_native.sh`. Preserve old source and binary bindings when
rebuilding. No timing or replay record is regenerated merely by extraction.

## Checkpointed execution

The [frozen protocol](PROTOCOL.md) preserves V1 and a pre-execution V2 memory
sizing amendment. BSGS uses m=2^27 and a measured 4-GiB table per generator,
reused across ten targets. Setup/query/verification/destruction costs are kept
separate; any derived or amortized sums retain their own labels. Rho uses the
pinned native default negation-map walk and a secure per-instance seed shared
across curves. Neither solver receives the checker scalar.

The native worker saves every target and all returned statistics/resources,
including failures, as atomic records. Complete cases are resumed only with
the same frozen identity. An interrupted BSGS table is rebuilt and that new
setup is charged as a retained new attempt. Field counters and per-step event
streams are not currently collected and are not imputed as zeros.

`source/launch.sh` performs the frozen checksum preflight and runs two BSGS
workers plus two rho workers. It first spreads coverage across curves, then
continues the full 100-problem-per-curve workload. One invocation per campaign
is required: do not start competing launchers against the same job paths.
Creation of STOP requests a checkpointed stop after the current job unit.
Check actual launcher receipts and process state; this file alone does not
claim the full campaign has been launched or completed.

Reports, per-curve and per-problem summaries, heatmaps and compact raw checkpoint
archives are appended as separate snapshots. No missing curve measurement is
invented. Full IC remains a separate implementation/recovery obligation;
relation-yield/regularity diagnostics cannot discharge it.

Historical recovery metadata remains in
[the earlier handoff](../photo_volcano_recovery_20261017/README.md). Its missing
archive problem remains open for that old evidence, but no longer prevents
this newly authorized replacement run.

## Audited diagnostic and visualization update

The [source-linked topology report](reports/topology_restart/REPORT.md),
[inspected two-page PDF](reports/topology_restart/REPORT.pdf) and
[vector volcano](reports/topology_restart/volcano.svg) retain all 46 vertices
and 90 directed maps. Its colors indicate level only.
The reusable [diagram skill](../../.agents/skills/isogeny-volcano-diagrams/SKILL.md)
supports arbitrary recorded depth, surface cycles, loops and parallel maps.
Its native Node tests passed at depths 0, 1, 4, 8 and 20; this verifies the
renderer, not mathematical discovery of those test graphs.

Actual Semaev/factor-base Gröbner probes ran on every curve: 9200 random
setup trials, 4600 target trials and 920 deliberately planted positive
controls. Singular independently verified equality of the returned ideals
for **all 14720 systems**. Random setup and target relation yields were both
zero under the frozen B=12, strict signed two-sum configuration; setup rank
remained zero and **no random-corpus full IC recovery occurred**. A separate
planted full-rank control passes GB → relations → linear algebra → verified
point equation and is explicitly excluded from those recoveries. See
[protocol](IC_PROTOCOL.md), [audit](validation/singular_audit.json), and
[diagnostic summary](reports/ic_audited_snapshot/summary.json).

The measured target top-forms Hilbert cutoff is 13; maximum temporary
accumulator degree is 23 and processed-pair sugar degree is 26. These are
separate quantities. None is presented as the CM regularity of the
homogenized ideal. The exact support audit has 288 nonzero points per
factor base out of r−1 possible nonzero subgroup points. Paired transports
across curves preserve relations, so repeated curve copies are not
independent yield samples. No bounded zero-yield run proves impossibility.

For answer-free independent replay, compile the added wrapper against the
extracted public native source snapshot and library:

```
g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror source/replay_public.cpp -Lnative -lcryptanalysis -Wl,-rpath,'$ORIGIN/native' -pthread -o replay_public
./replay_public "$PWD" 0 0 rho
```

It reads only the public P/Q inputs and public walk seed, checks the returned
point equation independently, and exports no recovered scalar. Its BSGS mode
is an independent cold setup+query solve, not a retained-table query timing.
For a private resumed campaign use `source/resume.sh` only after checking
that no campaign worker is active; its lock and process preflight reject a
competing launcher and its logs have a new per-launch namespace.

Full 4600-case rho/BSGS completion, full random-corpus IC recovery, rho/BSGS
field-level counters and complete per-step event streams remain unresolved.
Raw records, failed attempts and private randomness remain preserved locally.
The new sources do not modify the already frozen running workers or fixtures.

## 2026-10-17 14:50 UTC measured snapshot

**All 4600 rho cases are now independently verified, with zero failed cases.**
BSGS has 3297 verified cases at this snapshot; the interrupted execution
session was no longer connected, so remaining jobs were resumed with the
same frozen identity and fresh per-launch logs. All interrupted attempts are
retained and rebuilt table setup is recorded again. No completed job was
replaced. Check LIVE_STATUS.json locally for later actual progress.

The [quantitative report](reports/quantitative_1452/REPORT.md) includes an
[inspected 13-page run PDF](reports/quantitative_1452/RUN_REPORT.pdf),
[inspected eight-page results PDF](reports/quantitative_1452/RESULTS_REPORT.pdf),
seven full-graph heatmaps and per-curve/per-generator/per-problem min/mean/max
summaries. Every one of the 100 rho problems now has all 46 curve results.
Curve-average rho counts range from 107895789.39 to 134494184.67 group
operations, ratio 1.246519; individual cases range from 3759524 to 429007765.
These descriptive spreads do not establish an intrinsic hardness difference.
BSGS's paired online operation counts agree across every observed curve and
generator. All 4600 IC target GB stages count exactly 217410 field
multiplications, despite wall times from 11.519667 to 266.013969 ms under L0.

The [public IC metrics/replay archive](artifacts/public_ic_metrics_and_replay.tar.xz)
contains all 14720 polynomial-system/ideal replays, public factor-base point
coordinates and aggregate metrics. It excludes scalar witnesses, known random
setup coefficients, raw getrandom journals and raw native case records.
Its checksum is beside it. Source data and generated visuals remain separate
from the full private checkpoint, which has not been externally exported.

The [benchmark mirror](../../benchmarks/20261017-photo-volcano-restart/) copies
the same PDF/SVG/150 dpi PNG bytes. Missing field counters, successful full IC
recovery and full BSGS completion are still open; no fastest-tool or complete
cold-cost claim is made. `source/monitor.mjs` observes the real native workers,
appends local status, and generates a new local report on native completion.
It sends no data externally; that final report still requires PDF inspection
and a later PR update. It never restarts a solver or declares full IC complete.

## Orders and arithmetic replay, 2026-10-17 continuation

The [separately frozen order phase](ORDER_PROTOCOL.md) independently recomputed
all 46 curve cardinalities and checked **9200 native point orders against
9200 PARI point orders**, including both P and Q of every problem. All match.
Its [inspected report and heatmaps](reports/orders_1525_v2/README.md) include
full-group invariant factors, setup and point-check costs, field/group counters,
per-curve/per-problem statistics, timer resolution, and the actual GP default
thread count (5). Reference point-order reduction counts 12381 field
multiplications including squares per point; this is not an ECDLP-solver count.
Every order record and the compiled bundle remain preserved locally. The
readable source and full per-curve/per-problem order tables are published;
scalar witnesses and RNG journals are excluded.

A [new, separately frozen native arithmetic replay](FIELD_PROTOCOL.md) now
runs on the exact full corpus. Its instrumented backend passes the same
376415 group, 2063 BSGS and 134 rho checks, as well as counter classification
and reset controls. The first actual rho replay matches the original answer
and all nonclock solver statistics. Two rho counting workers are active;
the two BSGS counting workers wait for all original BSGS batches to complete.
The baseline binaries, data and timing records are unchanged. No replay
coverage is claimed before its case records exist.

The readable [instrumented backend](field_backend/INSTRUMENTATION.json) and
[source wrapper](source/field_replay_public.cpp) provide a tested
`field_replay_public` program that reads only public P/Q
fixtures and walk seeds, verifies the returned point equation, and exports no
scalar. It reproduces the first case's complete counter vector. Its BSGS mode
counts cold setup+query, while the checkpoint worker counts setup/query/free
separately. Rebuilding uses `-O3 -march=native -ftls-model=initial-exec` for
the instrumented shared library, loaded at process startup. Instrumented
wall time includes counter overhead and must not replace original timings.

Counters distinguish F_p from scalar-mod-r arithmetic and record actual
modular API invocations, Montgomery multiply/square/reduction, inverse and
Euclid divisions, and EC coordinate complements. They are not retired machine
instructions or a full state stream. Montgomery reductions overlap with the
multiplication/square categories; do not add them as independent costs.
Full coverage, full random-corpus IC recovery and unavailable hardware/state
metrics remain open. Check FIELD_LIVE_STATUS.json locally for actual progress.

The new compressed order/replay bundles are local: automatic upload review
rejected them because it could not verify an opaque archive against the
approved public subset. This change publishes readable source, aggregate
order tables and inspected visuals; it uploads no new compiled binary or raw
measurement journal. Existing approved public input archives supply the same
P/Q fixtures. Build an independent counter replay without touching the
running original files: `bash source/build_field_public.sh NEW_OUTPUT_DIR`,
then `NEW_OUTPUT_DIR/field_replay_public ROOT 0 0 rho`.

A direct own-process hardware-counter probe returned ENOENT for cycles,
instructions, cache misses and branch misses. These metrics remain null;
the primitive counters do not stand in for unavailable hardware events.
