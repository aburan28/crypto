# Bounded Modal pilot for the N53 C1800 combined target path

This prepared pilot runs C1200 packed and C1800 prefix plus prefetch/projective
target code on the same public T205 point in ABBA order. It checks Linux
x86-64 execution, exact scalar recovery, factor-base and rank identities,
five exclusive online phases, and guest CPU/memory fidelity counters before
a larger confirmation panel. The point belongs to the committed eight-point
local panel; this four-run pilot is an environment and correctness check.

The [transfer manifest](TRANSFER_MANIFEST.json) inventories exactly two
source-derived Linux x86-64 executables and one 1,446-byte public-point JSON
bundle, totaling **3,571,686 bytes** in the image. The 19,419-byte guest
runner source is also sent as the Sandbox command text. The bundle contains
candidate/workload IDs and expected public outputs but no fixture scalar.
The repository, local run archive, caches, and credentials are outside the
transfer. [`build/receipt.json`](build/receipt.json) records the source,
compiler, lockfile, binary, and ELF hashes; the binaries are currently in
`/private/tmp/s3-modal-combined-build` and can be rebuilt from the frozen
source if that temporary directory is cleared.

[`PILOT_FREEZE.json`](PILOT_FREEZE.json) fixes one Sandbox creation attempt,
eight physical CPU cores requested and capped, 32 GiB memory requested and
capped, a 1,200-second Sandbox limit, a 1,080-second workload deadline,
240 seconds per run, no explicit region, VM runtime, and blocked guest
networking. The guest records affinity, physical-core/SMT topology, NUMA
placement, PSI, cgroup throttling, steal ticks, interrupts, memory, and
instructions/cycles when perf access is available. Each inherited perf event
also retains its enabled and running times. IPC is exposed only when both
events ran for at least 99% of their enabled intervals; a scaled diagnostic
remains in the raw receipt when coverage is lower. The analyzer records
multiplexing as a failed fidelity check and unavailable perf access as
unavailable. Three synthetic tests cover scaled counts, short reads, and
the fidelity gate. Seven analyzer-contract tests pass four archived T205
reports through a synthetic guest receipt, check complete, invalid-scalar,
and failed-producer outcomes, and exercise CPU selection and affinity
validation. The guest selects complete SMT sibling groups and prefers one
NUMA node when available. The analyzer checks the child's initial affinity
and sampled thread masks without assuming the selected CPU list is sorted.
Before this correction, an interleaved sibling selection failed the affinity
check, and a guest exposing extra cores could select 16 distinct physical
cores instead of eight whole cores. These are schema checks, not remote timing
evidence. [Modal's Sandbox resource guide](https://modal.com/docs/guide/sandbox-resources)
describes a CPU request as a guaranteed minimum allocation and allows
bursting when host capacity is available. That allocation does not establish
the host-wide exclusive CPU partition required by the cryptanalysis isolated
benchmark rule, so this pilot retains `controlled_speedup: null` unless an
independent host-level isolation receipt is supplied.

At [Modal's published Sandbox rates](https://modal.com/pricing), the fixed
request/limit has a CPU-and-memory base rate of
`8 × $0.00003942 + 32 × $0.00000667 = $0.00052880` per second, or
**$0.63456** at 1,200 seconds, with any separately billed image, storage,
or other usage outside that figure. [Modal's resource guide](https://modal.com/docs/guide/sandbox-resources)
defines CPU requests as physical cores and describes request/limit billing.

`python3 pilot_modal.py --preview` validates every frozen hash and transfer
byte count without creating a Sandbox. The preview passed. No Sandbox has
been created and no source-derived file or guest runner has been sent to
Modal. The owner authorization required by the prior automatic approval
review is pending.

After the CPU-selection and affinity correction, the `PILOT_FREEZE.json` SHA-256 is
`7dfc53e9764e4a5e88966e1f32389d0916d56dba4820299d320d167a265033ae`;
the exact `TRANSFER_MANIFEST.json` SHA-256 is
`9b9e1a679f6705d4edb0b10bd9ea29d1176a018ac0e48a92c8c66d2e50300286`.
The guest runner source SHA-256 is
`31ecf18ac7f6c6bc0c1a64da877133851a4d287f03b0851e8a4fdee10ce3a2af`.
The three image files, their hashes, the transfer byte count, and Sandbox resource cap above are unchanged. This
revised manifest is the reviewable transfer identity for the C1800 pilot.
