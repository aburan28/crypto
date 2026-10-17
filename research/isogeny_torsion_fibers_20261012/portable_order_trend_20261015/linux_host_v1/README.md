# Independent Linux replay artifact, run 38831567671

The [GitHub Actions run](https://github.com/aburan28/crypto/actions/runs/38831567671)
started its `independent-replay` job at 2026-10-16 00:50 UTC and completed
at 00:54 UTC on commit `de340fa2161088a5188b841ec82278f5067f8aed`.
It used the frozen [four-order run spec](../BATCH_RUN_SPEC_V1.json),
public targets, native runner and independent checker. The workflow restored
the pinned runner/checker source bytes from that commit before building.
The [host receipt](HOST.txt) records Linux x86-64, four virtual AMD EPYC
9V74 CPUs on Microsoft Azure, and Rust 1.99.1. The [build hashes](BUILD_SHA256.txt)
identify the Linux binaries and sources.

The complete downloaded artifact is preserved in
[`RAW_ARTIFACT.tar.gz`](RAW_ARTIFACT.tar.gz): **1,330,324 bytes**, SHA-256
`b6c59b8e3c1efecd3b5695f472aee3b4d570ce3fcf7d4447eba7f60f9c9f1d88`.
It contains every child stdout, stderr, meter, exit status, host receipt,
and per-cell independent replay. The fresh
[`DOWNLOADED_ARTIFACT_FILES.sha256`](DOWNLOADED_ARTIFACT_FILES.sha256)
checks all **1,533 files** relative to the extracted artifact root; all
1,533 entries passed after archive extraction. Reproduce the check with:

```sh
mkdir -p /tmp/isogeny-four-check
artifact_manifest=$(realpath DOWNLOADED_ARTIFACT_FILES.sha256)
tar -xzf RAW_ARTIFACT.tar.gz -C /tmp/isogeny-four-check
(cd /tmp/isogeny-four-check/isogeny-four-linux-38831567671 && sha256sum -c "$artifact_manifest")
```

The artifact's original `sha256sums-before-manifest.txt` includes a hash
of itself because the first Linux script opened that file before running
`find`. Preserve it as a historical receipt, but use the fresh manifest
above for verification. The order-1,201 follow-on script excludes its own
manifest and uses relative paths.

The independent checker passed in each cell:
[`409/433`](REPLAY_409-433.json),
[`421/457`](REPLAY_421-457.json), and
[`457/421`](REPLAY_457-421.json). The native
[`different_j_host_performance` analyzer](../../../../examples/different_j_host_performance.rs)
then matched every Linux child arithmetic/scalar record against the
committed Darwin record before writing [the 252-child comparison](COMPARISON_V1.json).
Its comparison JSON is 17,962 bytes, SHA-256
`2aeada935e9e2715ed64ab5390ca070fa8eef5d98bca9699844169a619d7539e`.
The visual/report is [beside the run records](../reports/run-20261016T0050Z-linux-host/REPORT.md).

This hosted virtual machine has no exclusive physical-core reservation.
Its process CPU ratios are descriptive cross-host evidence; the separate
L2 isolated-wall and complete portable-price gates remain open.
