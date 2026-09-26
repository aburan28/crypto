# Cycle escape v3: correctness receipts

The original review was against `1000ae5`; this repair starts from `e07fd850`.
The newer baseline had already repaired the missing four-term tau checks and
measured v2 merge loss. V3 changes the exit decision, not that historical record.

| Check | Result | Scope |
|---|---|---|
| Pairwise tag patterns through 8 steps; both tau4 families | Pass | Existing hint predicate, not a complete cycle proof |
| Synthetic cyclic entry phases, cold/warm histories | Pass | Common anchor for detected cycles |
| A raw cycle containing a distinguished point | Pass | No escape bypasses a report |
| F131 `[1184]P` two-cycle, both entries and automorphisms | 6 common-orbit exits | Actual curve; constructed states |
| Packed vs reference on the actual cold path | 15 matching transitions per layout | Host execution of packed code, including hybrid layout |
| Sigma and table production tests | Pass | Verified collision recovery, I/O fault isolation, bounded restarts |
| Table checkpoint round trip / old version rejection | Pass | History preserved; old rule refused |
| Offline merge pair -> table client | Pass | Independent `[k]P = Q` verification |
| Worker/spool tests | 26 passed | Includes v3 framing |
| Ingest tests | 110 passed | Includes v3 header stripping |

The first table production test exposed a freshly written header still buffered
when the corpus was reopened for reload. The writer now flushes the header before
opening the second handle. The failing receipt is retained in
`production-table.log`. The first ingest run lacked boto3; after installing the
existing test dependency, all tests passed. Neither failure was hidden.

No GPU executed here. CUDA compilation and the Rust sigma certificate are CI
gates; their results are not represented by the local C++ receipts. No speedup
or campaign-scale zero-loss claim is made. New cold-path work invalidates reuse
of old v2 throughput figures for v3. Rare/unrecognized cycles and guard-induced
trail truncation remain bounded failure modes, described in CYCLE-ESCAPE-V3.md.

Commands are the Make targets `test-cycle-rule`, `test-cycle-escape`,
`test-table-walk-host TABLE_PROBE_POINTS=16`; `build/test-production`; the same
production translation unit compiled with `-DECC_WALK_TABLE=1 -DECC_BATCH=2`;
and unittest discovery for `test_worker_spool.py` and `test_dp_ingest.py`.
The workflow repeats the focused gates and adds CUDA compilation without a GPU.
