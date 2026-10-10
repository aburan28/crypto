# Same-target S3 worker screen on the N53 C1800 index

The frozen one-target screen found that **one target-query worker was faster
than four on seven of eight public targets**. Across T205–T212, the geometric
mean of one-worker/four-worker online wall time was **0.565×** (target-resampled
95% bootstrap interval **[0.411, 0.781]×**). Equivalently, four workers took
**1.77×** the one-worker online time on this local panel. All 16 fresh-process
solves reached prepared rank 1,800, passed native scalar and phase checks,
stayed within the 32 GiB observed-RSS envelope, and passed independent
checked-Sage scalar replay. The predeclared 1.25× advance gate failed.

Both arms use the same [Rust source](../../examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow.rs),
binary hash, factor-base set, 14 preparation workers, C1200-first/C1800-full
search schedule, eight-state inversion batch, root-slot prefetch, projective
scalar check, and full-rank row-elision path. The baseline has one
target-query worker; the candidate has four. The base contains **190,800**
distinct subgroup-usable points before orbit folding and **1,800** relation
columns. Each worker starts at a different index-state offset; the first
verified relation triggers cooperative cancellation, with all completed S3
calls and worker drain charged to that target. Candidate identities and exact
flags are in the [pair freeze](PAIR_FREEZE.json).

| Public target | One worker online ms | Four workers online ms | One/four ratio | Charged S3 calls, one/four |
| --- | ---: | ---: | ---: | ---: |
| T205 | 0.097 | 0.248 | 0.390× | 128 / 866 |
| T206 | 0.089 | 0.229 | 0.388× | 128 / 710 |
| T207 | 0.086 | 0.315 | 0.272× | 128 / 1,880 |
| T208 | 0.104 | 0.194 | 0.539× | 160 / 272 |
| T209 | 0.217 | 0.205 | 1.058× | 560 / 810 |
| T210 | 0.107 | 0.215 | 0.500× | 224 / 1,042 |
| T211 | 0.236 | 0.243 | 0.972× | 624 / 660 |
| T212 | 0.217 | 0.237 | 0.914× | 672 / 992 |

The online clock starts with the public target available and the reusable
index and base logs prepared. It includes target query, charged target PDP,
relation check, descent, and scalar recovery check; process launch, input
loading, and reusable setup are outside it. The paired geometric one/four
ratio for target query plus charged PDP was **0.517×** and for target descent
was **0.780×**. The largest observed process RSS was **25,140,281,344 bytes**.
The T025 controls independently verified two and four workers before this
screen; their archived native receipts and compressed reports are under
[`local/control`](local/control).

[`PAIR_FREEZE.json`](PAIR_FREEZE.json) binds the eight previously frozen
public points and one-target workload IDs, candidate manifests, source and
binary hashes, alternating arm order, one repetition per arm and point, the
180-second per-process cap, runner hashes, and the advance threshold. Raw
fresh-process receipts, compressed reports, stdout, stderr, ledger, host
record, and completion status are under [`local/paired_screen`](local/paired_screen).
[`analysis.json`](analysis.json), [`phase_summary.json`](phase_summary.json),
[`sage_runtime_info.json`](sage_runtime_info.json), and
[`sage_replay.json`](sage_replay.json) retain the calculations and independent
replay. The raw ledger SHA-256 is recorded in `analysis.json`.

This local host has no auditable exclusive CPU partition receipt; the
controlled-speedup field remains `null`. T205–T212 are development targets
already used in earlier kernel selection. The fresh T213–T246 confirmation
panel remains reserved for the C1200/C1800 comparison, whose execution plan
is in the [34-point freeze](../koblitz-s3-c1800-confirmation34-20261010/PANEL_FREEZE.json).
