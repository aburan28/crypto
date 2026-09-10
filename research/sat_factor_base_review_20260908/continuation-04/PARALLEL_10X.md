**Tenfold end-to-end SAT index-calculus speedup**

This continuation parallelizes the dominant relation-collection stage of the complete toy `K_1/F_(2^19)` index-calculus run. Each SAT target is independent, so a deterministic batch of public `(a,b)` targets can be solved concurrently and consumed in the original sampler order. The final policy uses ten Rayon threads, batches of ten targets, the selected 152-point factor base, signed Frobenius relation columns, verified rank stopping, and the previously selected 10,000-conflict per-target cutoff.

Every admitted time starts before curve construction and ends only after factor-base materialization, SAT relation collection, modular linear algebra, unknown-scalar recovery, and verification that the scalar reproduces the public target point. Offline factor-base discovery is common parameter tuning and is excluded from both arms; materialization of the fixed base is included in every run.

**Matched control.** The control records are the completing serial arms from continuation 03: identical curve, factor-base points, target scalar, sampler seed, target order, native-XOR `S_4` equations, exact domain trie, trace constraint, 10,000-conflict cutoff, model cap, and trial cap. The control retains eight ordinary Frobenius columns and a fixed two-relation surplus. The optimized arm combines `P` and `-P` into four signed columns, attempts a verified solve after five relations, and launches independent SAT targets in deterministic batches.

| Run seed | Unknown scalar | Parallel optimized | Serial Frobenius control | End-to-end speedup |
|---:|---:|---:|---:|---:|
| 1 | 4,242 | 24.014539 s | 291.720155 s | 12.147648x |
| 2 | 9,001 | 25.637970 s | 276.528004 s | 10.785878x |
| 3 | 17,001 | 16.539749 s | 194.972048 s | 11.788090x |
| 4 | 30,001 | 22.142284 s | 154.822686 s | 6.992173x |

All eight arms recovered and independently reverified their unknown scalar. No admitted run used a direct-relation shortcut, produced an invalid SAT model, or treated `UNKNOWN` as a relation or refutation.

| Aggregate | Parallel optimized | Serial control | Ratio |
|---|---:|---:|---:|
| End-to-end time | 88.334542 s | 918.042893 s | **10.392796x speedup** |
| Target trials | 110 | 184 | 1.672727x fewer |
| SAT conflicts | 1,018,353 | 1,669,784 | 1.639691x fewer |
| Accepted relations | 21 | 40 | 1.904762x fewer |

The median paired speedup is 11.286984x, the geometric mean is 10.194145x, and the minimum is 6.992173x. Seed 3 was used to select the batch size. Excluding it, the held-out seeds 1, 2, and 4 total 71.794793 seconds optimized versus 723.070844 seconds control, a **10.071355x holdout aggregate speedup**.

The result is a wall-time speedup on ten CPU threads. Parallel batching may perform speculative work beyond the first five usable relations; all launched attempts are charged in the reported trials, conflicts, CPU consumption, and end-to-end wall interval.

**Batch-size selection.** A bounded screen on seed 3 held the solver, conflict cutoff, target stream, and thread count fixed while changing the deterministic batch size:

| Batch size | Trials launched | Relations | End-to-end time |
|---:|---:|---:|---:|
| 8 | 24 | 5 | 25.154888 s |
| **10** | **20** | **5** | **16.539749 s** |
| 12 | 24 | 5 | 21.697368 s |
| 16 | 32 | 7 | 31.195652 s |

Batch 10 was fixed before running seeds 1, 2, and 4. Larger batches overproduced relations; smaller batches paid more barrier overhead. Seed-4 experiments with batches 14 and 30 were slower than batch 10 and are retained as non-admitted tuning evidence.

The lower 5,000-conflict cutoff was also tested on seed 4 with batch 10. It reached the 64-trial bound with only four relations and failed recovery. That run is preserved as an operational failure and is excluded from every speed aggregate.

**Implementation.** `KoblitzIcOptions::relation_batch_size` controls independent SAT targets per batch. Non-SAT strategies remain serial. Target pairs are generated from the seeded sampler before dispatch. Rayon evaluates each indexed target independently, and its indexed collection preserves sampler order. Solver statistics account for every launched target. Any invalid SAT model rejects the entire run. At a batch boundary, relation rows are appended in deterministic order and the modular solve is attempted only when enough equations exist; a candidate stops collection only after reproducing the public target.

`KoblitzIcReport` records the requested batch size and number of launched batches. The end-to-end example records the Rayon thread count, all SAT work, relation counts, stage times, and final verification status. A small-curve regression test exercises the parallel driver and checks recovery against the serial path.

`verify_parallel_10x.py` fails closed on mismatched seeds, targets, factor bases, conflict budgets, column modes, relation counts, batch sizes, thread counts, direct shortcuts, invalid models, or failed scalar verification. It checks internal time against external `/usr/bin/time`, requires every pair to improve, requires at least 10x full and holdout aggregate speedup, and separately verifies the batch-size screen and failed lower-cutoff arm.

Reproduction commands:

```sh
cp research/sat_factor_base_review_20260908/continuation-01/source_snapshots/Cargo.lock Cargo.lock
cargo build --release --example koblitz_selected_e2e --locked
RAYON_NUM_THREADS=10 target/release/examples/koblitz_selected_e2e optimized 1 4242 10000 10
RAYON_NUM_THREADS=10 target/release/examples/koblitz_selected_e2e optimized 2 9001 10000 10
RAYON_NUM_THREADS=10 target/release/examples/koblitz_selected_e2e optimized 3 17001 10000 10
RAYON_NUM_THREADS=10 target/release/examples/koblitz_selected_e2e optimized 4 30001 10000 10
python3 research/sat_factor_base_review_20260908/continuation-04/verify_parallel_10x.py
```

The serial control commands and artifacts are recorded in continuation 03.

**Evidence boundary.** This is a repeated, measured tenfold aggregate wall-time speedup for complete toy degree-19 SAT index-calculus runs on the recorded ten-thread host. It is a constant-factor parallel and relation-accounting improvement. It does not establish an exponent improvement, a single-thread 10x speedup, superiority over Pollard rho, or a threat to deployed binary curves.
