# C2800 shared-index S3 lookup diagnostic

Late root-slot prefetch reduced repeated N53 C2800 target lookup wall time by
**1.15×** on the local host (18 balanced triplets, bootstrap 95% interval
1.10–1.21×). Issuing the prefetch earlier increased lookup cost. All 54
queries recovered the frozen T025 scalar and the same first relation witness.
This result selects late placement for further candidates; the owner’s 2×
single-target advancement threshold was not reached by this lookup change.

## Method and evidence

The frozen source and binary hashes, three exact candidate manifests, T025
public point, one-target workload ID, and 54-mode schedule are in
`DIAGNOSTIC_FREEZE.json`. The program prepares one exact C2800 root index,
solves T025 once with the late mode, then measures the secondary lookup
diagnostic on the same public point. Each of the 18 balanced triplets contains
one query per mode, with a 64 MiB cache flush before each query. It records
wall and main-thread CPU time, witness, scalar, S3 calls, and relation check.
The previously known fixture scalar is used only for correctness checking.

`run_diagnostic.py` verifies hashes before launch and checks the primary solve,
all diagnostic rows, phase accounting, reference report, and 32 GiB memory
envelope. `analyze.py` checks the receipt, raw hashes, balanced schedule,
scalars, witnesses, and calls, then resamples the 18 paired triplets 20,000
times. `analysis.json` records the full numerical results. The raw primary
report, ledger, stdout, stderr, and receipt are under `local/verified_t025/`
as deterministic gzip files, bound by `archive_receipt.json`.

| Lookup mode | Median wall ms | Median thread CPU ms | Paired wall ratio against late |
| --- | ---: | ---: | ---: |
| No prefetch | 0.138 | 0.138 | 1.148× [1.099, 1.205] |
| Late prefetch | 0.121 | 0.122 | 1.000× |
| Early prefetch | 0.160 | 0.161 | 0.741× [0.663, 0.808] |

The primary T025 solve used the late candidate, reached prepared rank 2800,
verified the scalar and group relation, and had an online interval of
0.511 ms. Peak RSS was 25,823,084,544 bytes (24.05 GiB). The diagnostic
queries share one prepared index and repeatedly use a known point, so their
paired ratios are stage diagnostics, not the primary previously unseen
one-target online metric. No host-level exclusive CPU partition receipt was
available on this machine; `controlled_speedup` remains `null`.

To replay the analysis from the archive, decompress the five `*.gz` files to
one directory, verify each raw and gzip SHA-256 against `archive_receipt.json`,
then run `python3 analyze.py --run-dir <directory> --output <new-output.json>`.
The analyzer refuses to overwrite an existing output.
