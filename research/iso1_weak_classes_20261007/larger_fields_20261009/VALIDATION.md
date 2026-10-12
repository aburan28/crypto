# Validation receipt, 2026-10-09

The expanded controls verify 17 ordinary source/target cardinality pairs
through log2(Q)=252, 49 geometric fixtures, 34 exact ICV1 identities,
and the prime-degree orbit theorem. Every larger-field cell completed
inside its registered cap. The complete-prime continuation is active
at p=61 under an independent finite local service; p=59 has completed.

| Check | Result and evidence |
| --- | --- |
| Norm, fourth root, explicit source-point map, two independent target halves | 49/49 pass; 17 geometry processes, exit 0 and empty stderr |
| Independent source and target counts, mod-16, Hasse, ordinary trace | 17/17 pairs pass; 34 counts; all count processes exit 0 and have empty stderr |
| Geometry/count fixture agreement | All field moduli and first lambda vectors agree between the two modes |
| Exact model identities and raw receipt hashes | Native checker passes 34 ICV1 identities and 68 raw-output SHA-256 checks |
| Short-model coordinates and j-invariants | Separate PARI ellchangecurve replay matches all 17 recorded model pairs exactly |
| Cleared orbit-sum polynomial | IDC1h90d58cc0e0c48fe3 accepts 32 recorded + 64 extra evaluations; constant mutation refused |
| Exact integer expansion | PARI/GP returns zero; empty stderr |
| Cyclic orbit populations | Eight complete audits, 6,928,970 nonidentity parameters; every size population matches |
| Stabilizer classification | 45 p/n cases match exact gcds and exceptional size-2 orbits |
| Focused report package | 8 tests pass, including four prior replay tests and four new orbit/model tests |
| Native queue package | 13 queue/SHA/storage tests pass after the relocated-manifest change |
| Degree-5 workflow smoke | Two curves over F_(13^10) pass the exact final workflow command |
| Formatting and dashboard | Rustfmt and staged whitespace checks; canonical progress JSON equals embedded data; every prior cost point retained |
| Full root library suite | Exit 101: sparse worktree excludes src/lib.rs; raw output retained |
| PDFs and figures | Four-page report and both figure PDFs generated; every final page and figure inspected |

Native commands from the isolated repository root:

```sh
rustc --edition=2021 -O research/iso1_weak_classes_20261007/larger_fields_20261009/run_controls.rs -o research/iso1_weak_classes_20261007/live-census/run_large_controls
research/iso1_weak_classes_20261007/live-census/run_large_controls research/iso1_weak_classes_20261007/larger_fields_20261009 research/iso1_weak_classes_20261007/larger_fields_20261009/evidence_run1
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bins
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/census_runtime/Cargo.toml --bin iso1_census_queue
cargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin iso1_prime_degree_check -- research/iso1_weak_classes_20261007/larger_fields_20261009/orbit_identity_certificate.json
cargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin iso1_large_records_check -- research/iso1_weak_classes_20261007/larger_fields_20261009/curve_records.json
```

PARI/GP 2.17.3 supplies finite-field arithmetic and independent cardinality
calls. Each process starts with a 256 MiB stack; the geometry cap is
30 seconds, and the cardinality cap is 120 seconds per field. The count
phase takes one source/target pair per field; the geometry phase takes
49 fixed-seed cases. Each seed is 20261009+i, and each exact irreducible
modulus is retained. Timings are stage resource accounting on the
ordinary host; no isolation record is assigned to this round.

The first native-driver build used an incorrect relative SHA source path;
the next used an unavailable SHA function name. Both compiler failures
were corrected before any larger-field control was launched. The final
native driver and GP source hashes are retained in
[evidence_run1/source_freeze.txt](evidence_run1/source_freeze.txt).
The p13 development checks precede the larger-field run; final CI uses
the retained degree5_smoke outputs and the registered final source.

The earlier temporary p59 queue was absent at this round's start.
The first launchd attempt could not read the external-volume script;
its OS error transcript is retained. The internal runtime bundle was
verified by source/executable hashes. A provisional legacy submit job
was stopped before completion because it restarts on exit. The final
LaunchAgent uses RunAtLoad=true and KeepAlive=false, runs independently
of the chat, and preserves failures or completion without relaunch.

At the recorded snapshot, queue PID 72083 is a child of launchd PID 1,
and its p59 census child is PID 72144. Live output is
`/Users/adamburan/Library/Application Support/crypto-iso1/census-20261009-run3`.
The source freeze and status snapshot are retained beside the launch
configuration. The frozen census executable remains byte-identical to
the p53 executable. Only the queue's explicit runtime-manifest location
changed; arithmetic kernels are unchanged.

The preceding published commit c12ceeb71546700aa329e5a72791f3f53e3d7e80
passed the independent Linux ISO-1 workflow in
[run 37892969993](https://github.com/aburan28/crypto/actions/runs/37892969993).
This round extends that workflow with prime-degree replay, model/hash
checks, and a degree-5 point-count smoke. Remote results for the new
commit remain separate from local validation.

The previous October 8 hash seal belongs to commit c12ceeb71. It remains
historical evidence. This round adds its own evidence seal for the new
sources, records, artifacts, and current dashboard context.

Readable transcript copies normalize trailing whitespace; adjacent .raw.gz
archives retain the exact original bytes and were checked by decompression.
Normalized copies: native_tests.stdout, census_restart_source_freeze.txt, orbit_tests.stdout.

The p59 archived result completed during this round. Its imported zstd
archive was independently decompressed and audited: all 205,380 rows,
24,241,684 weighted representatives, 100,408 weak ordinary rows,
0/100,950 depth-1 positives, and 540/100,948 higher-depth zeros agree.
All 5,000 controls are positive, with 4,565 distinct absolute traces.
The queue has advanced to p61; the earlier PID snapshot remains
historical and the final status snapshot records the completed transition.

Every class-panel value was replayed from the six census CSVs using
header-named fields, including frobenius_2_depth; all counts match.
The initial documentary check selected two_split by position and was
corrected before acceptance. The plotted values were unchanged.
