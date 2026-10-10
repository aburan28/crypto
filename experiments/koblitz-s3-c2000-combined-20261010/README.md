# N53 C2000 combined S3 screen

The C2000 root index completed eight frozen N53 single-target solves, with
native scalar checks and independent checked-Sage replay. It made exactly the
same number of target S3 calls, inspected the same number of states, and found
the same first relation witness as C1800 on every point. The geometric
C1800/C2000 S3-call ratio is **1.00**, below the frozen **1.10** advancement
gate. C1800 remains the selected candidate for timed confirmation.

| Point | C1800 S3 calls | C2000 S3 calls | C2000 observed peak RSS, GiB |
| --- | ---: | ---: | ---: |
| T205 | 128 | 128 | 26.36 |
| T206 | 128 | 128 | 26.39 |
| T207 | 128 | 128 | 26.45 |
| T208 | 160 | 160 | 26.43 |
| T209 | 560 | 560 | 26.42 |
| T210 | 224 | 224 | 26.37 |
| T211 | 624 | 624 | 26.39 |
| T212 | 672 | 672 | 26.41 |

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) binds the source and binary,
candidate and workload identities, the C1800 reference analysis, public
points, resource envelope, and advancement threshold. The
[`analysis`](analysis.json), [`run ledger`](local/screen/ledger.jsonl),
[`completion`](local/screen/completion.json), and per-point compressed raw
reports retain the full screen. The independent
[`Sage replay`](sage_replay.json) used the saved
[`runtime receipt`](sage_runtime_info.json). Each run had 14 preparation
workers, one target worker, a 600-second cap, and a 32 GiB RSS envelope.
The C1800 and C2000 bases differ: C1800 has 190,800 actual usable points,
while C2000 has 212,000. Online wall times came from separate local runs
without a host isolation receipt and are diagnostics only.

## Initial control

The C2000 base and C1200-prefix target schedule reached complete rank 2,000
before the public T025 point, recovered its frozen scalar, and fit the
predeclared 32 GiB observed RSS envelope. The base contains **212,000**
distinct subgroup-usable points before sign/Frobenius folding, digest
`83907cda970c2ff8f087f884a3dc557ee32ad5684c3529f534d5f1d7944e1aa5`.
Its 212,000,000 retained S3 states produced 211,880,531 distinct root keys
in 536,870,912 packed 13-byte slots. Observed child peak RSS was
24,370,053,120 bytes. The checked scalar is 429728658711.

The T025 target required **512 S3 calls** and recovered the same first
relation indices `[55035, 67241, 57152, 189185]` as the C1800 combined
control. The C2000 target online interval was 0.267 ms; C1800's was
0.250 ms on separate, contended local executions. These single-point times
are diagnostics. The C2000 reusable setup took 78.972 s in this run and is
outside the single-target online interval.

[`CONTROL_FREEZE.json`](CONTROL_FREEZE.json) binds the source and release
binary hashes, checked public point, 14 preparation workers, one target
worker, rank threshold, 600-second cap, and 32 GiB memory envelope. The
[`receipt`](local/control/receipt.json) and
[`raw report`](local/control/report.json.gz) retain the full result. The
source uses the same late root-slot prefetch and projective scalar checks as
the C1800 combined candidate, while the root table includes the larger
C2000 base.
