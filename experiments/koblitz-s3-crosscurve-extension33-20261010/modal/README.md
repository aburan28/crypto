# Bounded Modal pilot for the frozen N41/N53 rank-243 S3 comparison

This pilot is ready to compare the same source snapshot and four frozen public
points on Linux x86-64: T025 and T028 on N41 and N53. It schedules 16
single-target solves, two reference and two prepared runs per point, in the
frozen ABBA/BAAB order. T028 covers the observed N41 online fallback path.
The Modal repetitions use R3/R4, distinct from the local R1/R2 records.

The [transfer manifest](TRANSFER_MANIFEST.json) inventories **3,329,083
bytes** in exactly two source-derived ELF executables and one 6,259-byte
public-point bundle. The 15,550-byte [guest runner](remote_pilot.py) and its
375-byte argument JSON are sent as Sandbox command text, for **3,345,008
bytes of total input payload**. The public bundle contains candidate and
workload IDs, curve IDs, factor-base sizes, and public points; fixture scalars, the
repository source tree, local run archives, caches, and credentials stay local.
The [build receipt](build_receipt.json) and [build log](build.log) bind the
original Git snapshot `b562cdaed14254d4f6c5f2a4640b302592e8882b`,
Rust/Zig toolchain, source archive, and both Linux binary hashes. The
[pilot freeze](PILOT_FREEZE.json) binds those artifacts, runner, analyzer,
run order, and resource limits.

The one-attempt pilot requests and caps **16 physical cores and 32 GiB** for
at most **1,200 seconds**. It has a 1,080-second work deadline, a 120-second
per-run cap, blocked guest networking, no fixed region, and terminates the
Sandbox in a `finally` block. At [Modal's published Sandbox
rates](https://modal.com/pricing), the CPU and memory base rate is
`16 × $0.00003942 + 32 × $0.00000667 = $0.00084416` per second, or
**$1.012992** over the full 1,200 seconds. Image, storage, and other
separately billed usage are outside that base calculation. Modal's
[Sandbox resource guide](https://modal.com/docs/guide/sandbox-resources)
specifies physical cores and request/limit billing.

The guest samples CPU affinity and thread masks, physical core and SMT
topology, NUMA pages, cgroup quota/throttling, PSI, steal ticks, interrupts,
RSS, and process cycles/instructions when perf access is available. The
[analyzer](analyze.py) checks the recovered scalars against local fixtures,
every charged online phase, candidate/workload/run IDs, and the bundle and
binary identities. The online S3 operation count adds target-dependent rank
calls to final target calls in the reference arm, and fallback calls to final
target calls in the prepared arm. Its quality score retains failed and
unavailable counters separately. Passing the pilot requires 16 verified
solves, ≥2× exploratory
online wall gain on both curves, and the declared A/A noise gate before a
larger Modal panel is considered. The prior [guest-only Modal host
probe](../../koblitz-s3-modal-host-probe-20261009/RESULT.md) could not attest
host-wide exclusive CPU partitioning. This pilot records
`controlled_speedup: null`; promotion needs a separate auditable host-isolation
receipt and analysis.

Local preparation completed: the Linux build passed, `pilot_modal.py
--preview` verified all hashes and bytes, and a 16-row synthetic guest
receipt exercised the analyzer's successful and wrong-scalar paths. The
synthetic receipt is an analyzer contract check, not a Modal measurement.
The automatic approval review previously rejected transferring
private-repository-derived binaries without explicit owner approval of the
exact payload. **No Modal Sandbox has been created and no payload has been
uploaded.** The approved launch command, once authorized, is:

```sh
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/modal/pilot_modal.py \
  --launch --output /private/tmp/s3-crosscurve-modal-pilot-receipt.json
```

The launch script writes an attempt marker before creation and refuses a
second creation for the same output path. The analyzer can then produce a
local replay ledger for the checked Sage launcher:

```sh
python3 experiments/koblitz-s3-crosscurve-extension33-20261010/modal/analyze.py \
  /private/tmp/s3-crosscurve-modal-pilot-receipt.json \
  --output /private/tmp/s3-crosscurve-modal-pilot-analysis.json \
  --replay-ledger /private/tmp/s3-crosscurve-modal-pilot-replay.jsonl
```
