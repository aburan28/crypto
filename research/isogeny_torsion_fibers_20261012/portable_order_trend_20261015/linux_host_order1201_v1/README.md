# Order-1,201 independent Linux artifact, V1

The frozen [`order1000` host protocol](../order1000/INDEPENDENT_HOST_PROTOCOL.md)
ran as [GitHub Actions 38838302462](https://github.com/aburan28/crypto/actions/runs/38838302462)
on commit `a5d81db3d47abd2ea9bdf3cfb471e357b1f4f12d` at
2026-10-16 02:08–02:13 UTC. The exact uploaded artifact is copied to
[`RAW_ARTIFACT.tar.gz`](RAW_ARTIFACT.tar.gz), 497,778 bytes, SHA-256
`b4ad7314d73b79b5e6573e5890a6a06066c80cafa10698737dc9115c812f773f`.
The download contains 515 files. Every file's SHA-256 is in
[`DOWNLOADED_ARTIFACT_FILES.sha256`](DOWNLOADED_ARTIFACT_FILES.sha256),
and all 515 passed the recorded
[`MANIFEST_VERIFICATION.txt`](MANIFEST_VERIFICATION.txt) check. The
archive and verification manifest are retained together so the run can
be reproduced independently of GitHub's artifact retention period.

The [host receipt](HOST.txt) records an `ubuntu-24.04` Azure VM with
four Intel Xeon 6973P-C vCPUs, two hardware threads per core, Linux
x86-64, and Rust 1.99.1. The virtual machine has no reserved physical
core. The [build hash receipt](BUILD_SHA256.txt) ties its binaries,
source, public input, and script to the run. The [independent checker
receipt](REPLAY_1237-1201.json) reports 84/84 successful children,
1,785/1,785 scalar replays, and 12,486 relation-witness checks. The
native [cross-host comparison](COMPARISON_V1.json) verifies that all 84
counted arithmetic vectors match the frozen Darwin records before
computing process-CPU ratios. The maximum child RSS is 189,353,984
bytes, below the 512 MiB ceiling.

Download and verification commands:

```sh
gh run download 38838302462 -n isogeny-fiber-order1201-linux-v1 -D /private/tmp/isogeny-fiber-order1201-linux-v1
(cd /private/tmp/isogeny-fiber-order1201-linux-v1 && shasum -a 256 -c /absolute/path/to/DOWNLOADED_ARTIFACT_FILES.sha256)
cargo run --release --locked --example different_j_order1201_host_performance -- /private/tmp/isogeny-fiber-order1201-linux-v1 research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000 NEW_COMPARISON.json NEW_COMPARISON.svg
```

The source of truth for descriptive timing ratios is the comparison
JSON, with its graph in [`COMPARISON_V1.svg`](COMPARISON_V1.svg).
The [run report](../reports/run-20261016T0550Z-linux-order1201/REPORT.md)
states the observed ratios and open pricing/isolation gates. This is an
independent-host correctness and CPU replay, not an isolated wall-time
session or a complete portable `S` price.
