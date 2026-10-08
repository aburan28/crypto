# n37 descendant-native six-sum batch against the same-Q strong rho

**Decision: fixed-cell batch diagnostic; no admitted speedup or no-go.** Two
native schedules each recovered and verified all 1,024 published b01 target
logs in ten IC and ten signed-Frobenius rho arms. All 20 IC outputs also
passed the separate source-curve replay, which rebuilt rank 42 and checked
every direct/shifted support decision and scalar. This closes the missing
same-Q **batch** comparator for the earlier [bounded residual result](../n37_native_residual_20261002/RESULT.md).
The medians of five interleaved blockwise IC/rho **complete-process CPU**
ratios were **10.298** in the first native schedule and **10.767** in its
formatting-only replication. Their separate descriptive log-t intervals
were **[9.728, 10.886]** and **[10.222, 10.952]**. Two first-schedule and
one repeat-schedule blocks exceeded the frozen 10% duplicate-drift limit; the
Apple M4 Pro host had no auditable exclusive CPU partition. These ratios are
exploratory. They are neither an isolated CPU speedup nor the primary
one-target online result, and say nothing by themselves about n41, n53,
degree 131, or the ECC2K-130 challenge.

The [protocol](PROTOCOL.md) froze the same public b01 Q, unchanged native IC
and rho producer sources and binaries, five alternating four-arm blocks,
180-second wall and 4 GiB RSS caps. The initial 20-arm pilot used a Python
process driver. Because this repository requires native research harnesses,
its files and source snapshots are retained as
[historical provenance](evidence/python-pilot/README.md), but **none of its
timings enters the decision below**. The protocol amendment and Rust driver
were committed as `c64a8ad2f` before the first native schedule. A
formatting-only amendment and Rust auditor were committed as `b52dc1c64`
before the second native schedule. Both are retained; neither schedule was
selected by its observed cost. The auditor reads only archived outputs and
published fixtures.
Every measured process started fresh, with empty in-process caches. The IC
arm rebuilt its K42 descendant-native base, signed 3+3 half table, relations,
rank and 16 frozen shifts; rho used 32 parallel signed-Frobenius walks,
normal-basis canonicalization, DP8, and a shared table **within the same
1,024-target batch**. Neither method read published fixture scalars during
the timed solve. The same public point file has SHA-256
`78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717`.

| Schedule / block | IC / rho child CPU | IC duplicate drift | Rho duplicate drift | Frozen noise gate |
|:--|--:|--:|--:|:--|
| First / 0 | 9.896 | 1.047 | 1.061 | pass |
| First / 1 | 10.399 | 1.115 | 1.175 | fail |
| First / 2 | 11.032 | 1.035 | 1.019 | pass |
| First / 3 | 9.870 | 1.195 | 1.092 | fail |
| First / 4 | 10.298 | 1.046 | 1.016 | pass |
| Formatted repeat / 0 | 10.810 | 1.016 | 1.015 | pass |
| Formatted repeat / 1 | 10.143 | 1.009 | 1.063 | pass |
| Formatted repeat / 2 | 10.767 | 1.025 | 1.039 | pass |
| Formatted repeat / 3 | 10.428 | 1.010 | 1.113 | fail |
| Formatted repeat / 4 | 10.773 | 1.004 | 1.012 | pass |

Each ratio divides the geometric mean of the block's two IC child CPU times
by the geometric mean of its two rho child CPU times. The table retains
all ten scheduled blocks, including noisy ones; no post-hoc trimming was
used. Median child CPU was 1,734.853/170.121 ms (IC/rho) in the first
schedule and 1,629.357/155.304 ms in the repeat. Median process wall was
1,952.864/216.918 ms and 1,820.165/198.907 ms respectively. All 40
processes exited zero, stayed below both resource caps, and passed their
native audits. Peak RSS across both schedules ranged 26.07–26.15 MB for IC
and 2.52–2.56 MB for rho. These process costs include startup, setup, all target work, output
and OS overhead; they cannot substitute for the one-target online interval.

| Fixed IC stage or control | First / formatted-repeat evidence |
|:--|--:|
| Descendant-native signed log columns | 42; independently replayed rank 42 |
| Signed half-table raw / distinct sums | 105,995 / 102,391 |
| Target-blind relation probes / misses | 49 / 7 |
| Direct target hits / shifted recoveries | 939 / 85 per arm |
| Target oracle queries / verified logs | 1,115 / 1,024 per arm |
| Setup and isogeny transport, median | 18.811 / 20.232 ms |
| Half-table construction, median | 54.271 / 49.448 ms |
| Relations, median | 75.115 / 83.805 ms |
| Target query and verification, median | 1,777.317 / 1,635.493 ms |
| Internal time excluding output, median | 1,933.659 / 1,800.933 ms |
| IC logical addition requests / explicit scalar multiplications | 33,647,429 / 4,289 per arm |
| Rho walk steps / group additions / explicit scalar multiplications | 596,818 / 632,240 / 3,104 per arm |

The IC addition requests and rho additions are produced by different
arithmetic paths; scalar internals, transport, field conversion and rho's
normal-basis work have no common calibrated group-addition price. Thus
`S = operations / sqrt(r)` and a counted IC/rho cost ratio are **unknown**.
The fixed policy's 10× descriptive process gap makes another local K42
batch-only micro-optimization a weak next bet; that prioritization is an
inference from this noisy engineering panel, not a degree-131 scaling law.

The compact [first](evidence/native/raw-panel-first.tar.gz) and
[formatted-repeat](evidence/native/raw-panel-repeat.tar.gz) native raw
archives are 1,707,441 and 1,707,388 bytes, SHA-256 respectively
`e4282a7efda282503be873561bdc3c6d9b3684d9bf32dcbf674bd6bf254a503d`
and
`e9adce3d3c78d46588bcd0df97361eeec400f68beaa40ba5d1ac9800542a08ed`.
Together they retain all 40 arm commands, output digests, exit and resource
records, full producer outputs, and 20 independent native replay receipts.
The [first audit receipt](evidence/native/VERIFICATION-first.json) and
[repeat audit receipt](NATIVE_VERIFICATION.json) were computed by the
[Rust auditor](../../../../examples/n37_native_batch_rho_audit.rs), which
checks the SHA-pinned executed sources, inputs and binaries; every row and
output digest; 20,480 recovered scalar outputs per schedule against public Q
and the post-run scalar fixture; replay rank and target support; phase sums;
and the frozen schedule and caps. Each receipt reports
`PASS_20_NATIVE_ARMS_20480_POINT_SOLVES`,
`all_duplicate_drift_pass=false`,
`timing_eligible_for_speed_claim=false`, and
`primary_one_target_online_speedup=null`. Extracting the committed archive
and rerunning that auditor reproduced each receipt byte-for-byte at its
own source snapshot. The first Rust driver and auditor source SHA-256 digests
are `281154a944fa85f0e6de14fa502c2e3a48865dfd3d94fd9e8a629f508315198c`
and `0b940fa2896fec272a6066f63a9989163822961c8dab37606d5ddf67eabda691`;
their exact source snapshots and protocol are preserved in
[`evidence/native/`](evidence/native/). The formatted driver and auditor
source SHA-256 digests are
`b1c56c19f63cbc6084c2a78c49f6820ba141fcdd64a787ef50d903e341316cbf`
and `ace47cc206554fd0659bc4b64b5d367a0841c6ab1b74b6c03363b0ccb6eabe3d`.

Reproduce from this PR head with the pinned
[`Cargo.lock`](../n37_native_m6_mitm_20261002/Cargo.lock) copied to the
checkout root (SHA-256 `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`):

```sh
cargo build --release --locked --offline --example n37_native_m6_residual --example n37_native_m6_residual_replay --example koblitz_rho_batch_ks_v3 --example n37_native_batch_rho_panel --example n37_native_batch_rho_audit
mkdir -p /tmp/n37-native-batch-rho-replay
tar -xzf research/notes/ecc2k130/n37_native_residual_rho_20261005/evidence/native/raw-panel-repeat.tar.gz -C /tmp/n37-native-batch-rho-replay
target/release/examples/n37_native_batch_rho_audit /tmp/n37-native-batch-rho-replay/native-raw-panel-repeat /tmp/n37-native-batch-rho-replay/VERIFICATION.json
cmp /tmp/n37-native-batch-rho-replay/VERIFICATION.json research/notes/ecc2k130/n37_native_residual_rho_20261005/NATIVE_VERIFICATION.json
```

**Next gate.** Use the established one-target same-Q contract to compare
the K8/K16 source, transported, descendant-native and pullback bases on
fresh n41/n53 points, with each base's *actual* usable points, folded
columns, setup, natural relation rank, target descent and verification
accounted for. Freeze the source-to-leaf route and endomorphism-order status
separately; a crater-depth label is not a proved conductor. Before any
performance promotion, acquire an isolated host receipt and a common
operation-count unit for IC and rho. Carry only a complete, verified route
with measured natural-query coverage toward degree 131.
