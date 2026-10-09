# Source-free Modal VM CPU fidelity probe

This probe measures which CPU isolation and noise signals a Modal VM exposes before any S3 timing panel is considered there. It sends only the generic `remote_probe.py` program as a `python3 -c` argument. It does not mount or upload cryptanalysis source, curve fixtures, targets, candidate tables, or a repository archive.

## Frozen protocol

- Create **one** Modal VM Sandbox using the installed Modal SDK version in `FREEZE.json`, a Debian slim image, a request for 16 physical CPU cores and 32 GiB memory, a 120-second Sandbox timeout, and `block_network=True`. Do not retry an ambiguous creation.
- Wait three seconds idle, pin the process to the first allowed logical CPU, then execute a five-second deterministic integer loop. Record snapshots before and after both intervals.
- Record guest-visible CPU model and topology, affinity, NUMA nodes and process page policy, cgroup quota and throttling, CPU and memory PSI, selected-CPU steal ticks and interrupts, frequency, process scheduling counters, and process perf hardware cycles and instructions. Record perf permission errors when counters are unavailable.
- Terminate the Sandbox after collecting output. Preserve the raw JSON receipt, including failures. The local runner checks its own SHA-256, the remote code SHA-256, the analyzer SHA-256, and the expected resource/transfer fields against `FREEZE.json` before creation.

The seven-item **guest observability score** counts whether pinning held, topology was visible, NUMA policy was visible, cgroup CPU controls were visible, both system and cgroup PSI were visible, steal/interrupt counters were visible, and process perf IPC was available. It is a coverage count, not a performance score. The analyzer also reports idle/busy PSI, throttling, steal and interrupt deltas. The score and individual checks are fixed before the remote observation.

Host-wide CPU exclusivity, fixed frequency, IRQ routing, physical SMT separation, and contention outside the VM cannot be certified by guest observations alone. Therefore `host_exclusive_partition_verified` is always false and `controlled_speedup` is always null. A future timing panel requires an auditable host isolation receipt and the repository's paired-run noise gates before a CPU wall-time ratio can be promoted.

## Commands

From this directory, after committing the freeze:

```sh
python3 probe_modal_vm.py --preview
python3 probe_modal_vm.py --launch --output receipt.json
python3 analyze_receipt.py receipt.json --output analysis.json
```

The launch command attempts Sandbox creation once. If its result is ambiguous, inspect Modal's Sandbox state before any later manual attempt. The probe does not require the repository Sage runtime because it performs no curve computation.

Current published Modal pricing lists Sandbox CPU at $0.00003942 per core-second and memory at $0.00000667 per GiB-second. At the requested 16 cores and 32 GiB, 120 seconds of requested compute is approximately $0.1013, excluding any other charges. Modal documents CPU allocation as elastic and says the requested CPU count is a minimum; its resource description does not establish exclusive host cores. See [Sandbox resources](https://modal.com/docs/guide/sandbox-resources), [VM Sandboxes](https://frontend.modal.com/docs/guide/vm-sandboxes), and [pricing](https://modal.com/pricing).
