# Native Modal retry for the N53 S3 target path

The native guest records the VM's actual CPU affinity and topology before it
runs the frozen T205 public point in baseline, expanded, expanded, baseline
order. It selects at most eight observed physical cores and sixteen logical
CPUs. Complete SMT sibling groups on one NUMA node are preferred; a partial
mask is recorded and remains usable for an exploratory run. This addresses
the first Modal pilot's topology assertion, which stopped before any timed
row. The [first pilot receipt](../modal/pilot-T205-20261010.json) records that
single creation attempt and its remote error.

The [freeze](PILOT_FREEZE.json) fixes four transferred files totaling
**6,185,082 bytes**: two pinned Linux x86-64 Rust executables, the 1,446-byte
[public bundle](../modal/public_bundle.json), and one static Linux x86-64 Go
guest. The bundle carries the public point and candidate/workload identities;
the known scalar stays in the local independent analyzer. The [guest build
receipt](BUILD_RECEIPT.json) records the Go toolchain, target, source, module,
binary, and ELF hashes. The Rust executables reproduce the hashes in the
[original rebuild receipt](../modal/rebuild-20261010/receipt.json). The
repository, private fixture, local run archive, and credentials are not
transferred.

The freeze permits **one** new Sandbox creation attempt after separate owner
approval. It requests and limits eight physical CPU cores and 32 GiB of
memory, uses a 1,200-second VM timeout, limits the four-run workload to 1,080
seconds and each run to 240 seconds, blocks guest networking, and selects no
explicit region. The control program writes an attempt marker before calling
`Sandbox.Create`, preserves the Sandbox ID as soon as it receives one, copies
the guest's partial receipt after any guest exit, and terminates the Sandbox.
An ambiguous create is never retried automatically. At the published
[Modal Sandbox rates](https://modal.com/pricing), the CPU and memory base
compute figure for the full timeout is **$0.63456**; any separately billed
usage is outside that figure. [Modal's resource guide](https://modal.com/docs/guide/sandbox-resources)
defines CPU requests as minimum physical-core allocations and describes
request/limit behavior.

The guest records allowed and selected CPUs, core/SMT and NUMA topology,
thread affinity, NUMA pages, cgroup quota and throttling, CPU/memory/I/O PSI,
interrupts, steal ticks, frequency settings, peak resident memory, and perf
instructions/cycles with enabled and running times. Perf IPC is available
only when both events ran for at least 99% of their enabled time. Each native
report's original bytes are included in compressed form and checked against
their SHA-256 digest. The independent analyzer verifies the T205 fixture,
point-file hash, candidate and run identities, scalar, point replay, base and
rank, target schedule, and the five exclusive online phases. Failures and
missing metrics stay visible in the analysis; four verified rows are required
for an exploratory paired ratio. The guest score counts 17 prespecified
checks, with failed and unavailable checks listed separately.

The pilot is an environment and correctness check on a point already used in
the local panel. Its four runs cannot establish the requested at-least-2×
result across fresh N41 and N53 targets. [Modal's Sandbox allocation
contract](https://modal.com/docs/guide/sandbox-resources) also does not
provide the host-wide exclusive CPU partition receipt required by this
repository's CPU measurement rule, so the analyzer leaves
`controlled_speedup` null. A controlled comparison requires an auditable
host-level isolation receipt and a fresh, paired one-target panel with
independent scalar verification.

## Local verification and launch interface

From this directory, the no-cost preview checks every frozen source and
transfer hash before a remote create:

```sh
GOTMPDIR=/private/tmp GOCACHE=/private/tmp/s3-go-cache go test ./internal/topology ./cmd/analyze
GOTMPDIR=/private/tmp GOCACHE=/private/tmp/s3-go-cache go build -o /private/tmp/s3-modal-native-launch ./cmd/launch
/private/tmp/s3-modal-native-launch --preview
```

The guest build command is recorded verbatim in `BUILD_RECEIPT.json`. It
cross-compiles a static Linux x86-64 binary on this macOS ARM64 host; local
tests use archived T205 reports and synthetic topology. Remote Linux guest
execution has not yet occurred. After a new approval, a single paid attempt
uses `--launch --output ABSOLUTE_RECEIPT_PATH`, followed by the native
analyzer with `--output ABSOLUTE_ANALYSIS_PATH ABSOLUTE_RECEIPT_PATH`.

The current freeze SHA-256 is
`bb1d0635bca2360bcd578dd15222ccf09978c71c35642467d1bd91babfa0131b`.
The first pilot's one approved creation attempt has already been consumed;
this retry has not been launched.
