# N83 intermediate root11 screen selects K2816 for the fresh-point holdout

The preregistered 128-point local screen selects **K2816/root11** as the
fourth arm of the sealed twenty-point Tier 1 comparison. Its geometric mean
target-online interval was **18.443 ms**, the smallest of the five tested
factor bases. Each arm solved all 128 public points; an independent checker
recomputed every fixture scalar, every enumerated base, and every four-point
relation. The exact K2816 candidate has **467,456 usable factor-base points**
and **2,816 folded relation columns**:
`IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b`.

| Root11 K | Usable base points | Folded columns | Geometric mean online, ms | Arithmetic mean online, ms | Complete shared-setup process, s | Allocated root state, GiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2304 | 382,464 | 2,304 | 32.102 | 56.733 | 323.363 | 2.257 |
| 2543 | 422,138 | 2,543 | 23.811 | 40.317 | 374.125 | 2.749 |
| 2688 | 446,208 | 2,688 | 25.106 | 40.405 | 461.860 | 3.072 |
| **2816** | **467,456** | **2,816** | **18.443** | **33.808** | **549.916** | **3.371** |
| 3072 | 509,952 | 3,072 | 21.962 | 44.997 | 514.623 | 4.012 |

K2816's geometric mean is 1.191 times shorter than K3072's on these same
128 points; it is shorter on 70 individual points. This is the frozen
fourth-arm selection statistic, not the primary one-target Tier 1 comparison.
The batch reused each arm's factor base, root index and log table across its
128 queries on an ordinary macOS ARM64 host with 48 GiB RAM. The complete
process times in the table include that shared preparation. The five arms ran
sequentially, so those process times are descriptive; no same-point rho
process ran in this secondary screen.

The [selection](results/batch_selection.json),
[full independent replay](results/batch_full_replay.json), and
[per-arm records](results/) retain all 640 target intervals, probes, base
hashes, source and raw-output digests. The public batch target SHA-256 is
`ac14c0326532b81041d357e921d99f86ef7ff4f6cf5e475b9c3a2bcc8997ac7c`;
the answer-sidecar SHA-256 is
`1fb66a7bf9bc1fd6f75214d1a633711d64c44e289aeffb86e5b9fb44cf7e46a2`.
The executed IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
The raw base dumps, target records, command and environment receipts, and
release binary remain under
`/Volumes/SSD990/crypto/artifacts/n83_k_intermediate_screen_20261017/`;
their individual SHA-256 digests appear in the selection and replay records.

The original prewritten checker stopped with `KeyError` because it requested
one combined PDP/relation-check timing field. The producer records
`target_pdp_ms` and `target_relation_check_ms` separately. The successor
checker summed those two recorded phases, retaining the same five-phase
online interval and every scalar, base and relation check. Both verifier
hashes and the original failure are retained in the
[correction receipt](results/verifier_correction.json); the original
[verifier plan](results/verify_plan.json) and [batch plan](results/batch_plan.json)
remain frozen.

The [twenty-point protocol](PROTOCOL.md) now runs K2304, K2543, K3072 and
K2816 on four distinct Tier 1 cloud allocations. Each point will have a
separately seeded same-point rho solve, alternating execution order,
before/after calibration and independent scalar replay. Its paired online
intervals and exact sign-flip tests will determine whether this local
selection transfers to fresh points and the cloud resource envelope.

The first four-host dispatch, [run 39037685265](https://github.com/aburan28/crypto/actions/runs/39037685265),
stopped in the preflight selector test before either producer was built or
measured. On the Linux runner, `import select` resolved Python's standard
module instead of this experiment's `select.py`; all four jobs reported the
same `AttributeError`. The [failed-dispatch copy](results/failed_dispatch_39037685265.json)
and [provider failure receipt](results/preflight_failure_39037685265.json)
retain the original head, four artifact digests, job-log hash, and failure
classification. The four provider ZIPs and full job receipts remain under
`/Volumes/SSD990/crypto/artifacts/n83_intermediate_holdout_20261017/failed_run_39037685265/`.
The selector test now loads its exact sibling file by path; a successor
dispatch will bind that correction while preserving the original failed run.
