# N53 C2000 combined S3 control

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
C2000 base. The next stage compares deterministic S3-call counts on T205–T212
using the same frozen points as the C1800 panel; any online wall ratio from
separate local days remains exploratory.
