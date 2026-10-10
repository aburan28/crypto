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

Build only these two examples for `x86_64-unknown-linux-musl` from the
checked-out source and `Cargo.lock`:

1. `s3_partial_prepared_parallel_targetspan` for the C244/C396 baselines.
2. `s3_partial_prepared_fullrank_targetspan` for the C800/C1200 candidates.

[`build_linux.py`](build_linux.py) checks both frozen source hashes and the
library tree ID, invokes Cargo offline through a portable Zig musl linker,
and writes its build log and binary hashes outside the repository. Run it
after the active local timed panel ends, so compilation does not add CPU
contention to that panel.

The Modal image transfer inventory will contain those two Linux executables
and one canonical JSON bundle of the selected public points, workload IDs,
candidate IDs, and scalar replay fixtures. The runner code is a separate
small Python command. The remote VM receives no repository tree, Cargo
cache, Sage installation, private credentials, or existing result archive.
Freeze the exact binary and bundle SHA-256 values and the Modal resource and
timeout envelope before creating one Sandbox. Do not repeat an ambiguous
Sandbox creation automatically. Terminate the Sandbox after retrieving its
receipt.

The previous [source-free VM probe](../koblitz-s3-modal-host-probe-20261009/RESULT.md)
verified guest process perf access and affinity but found no host-wide
exclusive CPU partition receipt. Modal's
[Sandbox resource guide](https://modal.com/docs/guide/sandbox-resources)
describes CPU requests as minimum reservations with possible bursting, and
its [VM guide](https://modal.com/docs/guide/vm-sandboxes) calls VM CPU
provisioning elastic. Guest pinning, PSI, and IPC will therefore be retained
as diagnostics; a controlled twofold CPU online claim still requires the
repository's auditable host-isolation gate and passing A/A noise checks.
