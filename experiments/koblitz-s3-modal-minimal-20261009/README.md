# Minimal-binary Modal S3 pilot

The frozen N41 C244/C800 and N53 C396/C1200 studies already bind two Rust
example sources, one shared library source-tree ID, candidate manifests,
factor-base digests, and public one-target fixtures. The first Modal pilot
will compare those same algorithms on one fresh public point per curve in
ABBA/BAAB order. It will test Linux correctness, capture each verified online
interval and its five exclusive phases, and collect guest CPU/NUMA/pressure
and process perf diagnostics where available. A later pilot can use twelve
independent targets per curve to estimate variance and a fixed confirmation
size.

Build only these two examples for `x86_64-unknown-linux-gnu` from the
checked-out source and `Cargo.lock`:

1. `s3_partial_prepared_parallel_targetspan` for the C244/C396 baselines.
2. `s3_partial_prepared_fullrank_targetspan` for the C800/C1200 candidates.

[`build_linux.py`](build_linux.py) checks both frozen source hashes and the
library tree ID, invokes Cargo offline through a Zig linker targeting glibc
2.28 (compatible with the Debian guest),
and writes its build log and binary hashes outside the repository. Run it
after the active local timed panel ends, so compilation does not add CPU
contention to that panel.
The [build receipts and logs](build/) preserve the successful GNU cross-build
and the earlier local toolchain and musl-link attempts. The successful
executables have SHA-256 values
`e999feba1bba11f8628493a6fe171f5d948d632aa87af05a78f06fd6fbce97d1`
and `18186665ee35c837ac8248ce818ee74d3513de50d8909eb625e553b9870afca5`.

[`prepare_bundle.py`](prepare_bundle.py) freezes the two public T181 points,
workload IDs, candidate IDs, and expected factor-base digests in one 2,411-byte
JSON file. The known-scalar fixtures remain local and
[`analyze.py`](analyze.py) compares them with the recovered scalars after the
guest receipt returns. The Modal image transfer inventory contains exactly
the two Linux executables and that JSON bundle. The
[`remote_pilot.py`](remote_pilot.py) text is passed as a Python command; it
captures per-run affinity, CPU and memory pressure, NUMA topology and sampled
page placement, cgroup throttling, steal ticks, interrupts, memory, and
process IPC where available. The analyzer assigns a transparent 16-check
guest diagnostic score, counting unavailable checks as zero, and records each
failed or unavailable check. IPC is retained as a whole-process mechanism
diagnostic and does not enter the fidelity score. The guest receives no
repository tree, Cargo cache, Sage installation, private credentials, or
existing result archive.

[`pilot_modal.py`](pilot_modal.py) verifies a frozen transfer manifest before
one Sandbox creation attempt. The requested and capped allocation is eight
physical CPU cores (16 vCPUs), with 16 logical execution CPUs selected in the
guest, 32 GiB requested and capped memory, and a 1,200-second Sandbox
timeout. The workload has a 1,080-second deadline and a 240-second cap per
process. The runner blocks guest networking and terminates the
Sandbox after retrieving its compressed raw receipt. An attempt marker is
written before creation so an ambiguous create is never retried automatically.
The maximum billed CPU-and-memory envelope at [Modal's published
rate](https://modal.com/pricing) is roughly $0.64 for the full 1,200 seconds,
plus image-build or storage charges if applicable.
[`PILOT_FREEZE.json`](PILOT_FREEZE.json) binds the exact 3,444,075-byte
transfer, the build and runner hashes, the public bundle, and the resource cap.
Its v2 resource fields correct an earlier prelaunch manifest that labeled
Modal's physical-core request as vCPUs; no Sandbox was created under v1.

The previous [source-free VM probe](../koblitz-s3-modal-host-probe-20261009/RESULT.md)
requested 16 physical cores and observed 32 logical CPUs. The new guest
runner requires 16 selected logical CPUs for its eight-core request. The
probe also verified guest process perf access and affinity but found no host-wide
exclusive CPU partition receipt. Modal's
[Sandbox resource guide](https://modal.com/docs/guide/sandbox-resources)
describes CPU requests as minimum reservations with possible bursting, and
its [VM guide](https://modal.com/docs/guide/vm-sandboxes) calls VM CPU
provisioning elastic. Guest pinning, PSI, and IPC will therefore be retained
as diagnostics; a controlled twofold CPU online claim still requires the
repository's auditable host-isolation gate and passing A/A noise checks.
