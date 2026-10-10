# Capped Modal pilot for the combined S3 target path

This prepared pilot compares the packed S3 reference with late root-slot
prefetch plus projective target checks on one frozen public T181 point each
for N41 C800 and N53 C1200. The two Linux x86-64 executables were built from
the accepted source tree with the frozen Cargo lockfile and a Zig glibc-2.28
cross-linker. `build/receipt.json` records the compiler, source hashes,
library-tree ID, build log hash, binary hashes, and exact ELF paths.

The only proposed Modal image transfer is the three-file inventory in
`TRANSFER_MANIFEST.json`: the two source-derived Linux executables and a
2,416-byte JSON bundle containing two public points, candidate IDs, workload
IDs, factor-base expectations, and run order. Its total is **3,572,376
bytes**. Fixture scalars, the repository, caches, credentials, and prior
result archives are outside that inventory. `pilot_modal.py --preview`
rechecks all file hashes, the bundle against the parent freeze, the source
receipts, and the transfer size without creating a Sandbox.

| Curve | Folded columns / usable points | Packed candidate | Combined candidate |
| --- | ---: | --- | --- |
| N41 | 800 / 65,600 | `IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0h161970671e93` | `IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0h700c174acd1e` |
| N53 | 1,200 / 127,200 | `IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hf8b850f65ece` | `IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h7d64b2ad00e8` |

`PILOT_FREEZE.json` fixes one Sandbox creation attempt, an eight-physical-core
request and limit (16 vCPUs), 32 GiB memory request and limit, a 1,200-second
Sandbox timeout, a 1,080-second workload deadline, a 240-second per-process
cap, no explicit region, and blocked guest networking. The eight jobs run
ABBA on N41 and BAAB on N53, with 14 preparation workers and one target
worker. The guest runner records exact output and all five online phases,
CPU affinity and topology, NUMA page placement, cgroup CPU/memory/io pressure,
throttling, steal ticks, interrupts, memory, and process instructions and
cycles when perf access is available. The analyzer replays local scalar
fixtures and assigns the frozen 16-check guest diagnostic score. It retains
unavailable checks as unavailable rather than treating them as passes.

At [Modal's published Sandbox rates](https://modal.com/pricing), the fixed
request/limit costs `8 × $0.00003942 + 32 × $0.00000667 = $0.00052880` per
second for CPU and memory, or **$0.63456** at the full 1,200-second timeout
before any separate image, storage, or other charges. [Modal's Sandbox
resource guide](https://modal.com/docs/guide/sandbox-resources) defines the
CPU request as physical cores and the upper limit as the burst cap.

The source-free [Modal host probe](../koblitz-s3-modal-host-probe-20261009/RESULT.md)
observed guest perf and affinity access but could not establish a host-wide
exclusive CPU partition. This pilot measures Linux correctness and guest
fidelity first. A controlled CPU speedup requires a qualifying host-level
isolation receipt and passing paired noise gates under the repository rule.
No Sandbox has been created for this pilot, and none of its three files has
been transferred to Modal.
