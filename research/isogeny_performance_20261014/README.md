# Native isogeny performance campaign

This campaign inventories the repository's isogeny algorithms, freezes native
CPU baselines, and measures successive optimization candidates. The dedicated
`isogeny-perf` binary runs without Cargo or a checkout. Start with
[COMMANDS.md](COMMANDS.md); the ReleaseMe archive integration is implemented,
while an actual new release asset remains unverified.

| Requested work | Current evidence and remaining gap |
|---|---|
| Comprehensive algorithm inventory | [INVENTORY.md](INVENTORY.md) separates constructor, finder, polynomial, path, action, model and higher-dimensional boundaries; many series remain open |
| Native baseline on several CPUs | Four 96-case output matrices pass; AMD has a complete admitted timing series; ARM's first series is censored; two macOS profiles have counts only |
| Iterative BMSS optimization | First three rungs have exact outputs and qualified AMD data; fourth candidate has 96/96 exact controls and unchanged counters, with new timing comparisons still required |
| Other constructors | Seven-method kernel adapter has two 84-case matched baseline runs on Apple M4 Pro; other CPU timing and optimization rounds remain open |
| Remaining families | Linux exploratory outputs retain all 24 legacy groups and verification flags; complete result fingerprints and matched controls are still required |
| Validation and delivery | Full root release gate: 4,104 passed, zero failed, 117 existing ignored; standalone suites/Clippy/fmt/workflow lint pass; ready PR1856 remains the publication vehicle |
| Named P192/P224 whole-process performance and additional maps | Separate requirements remain open; synthetic constructor fixtures are not counted as search discoveries |

The [CPU report](cpu-report-20261015/REPORT.md),
[reviewed PDF](cpu-report-20261015/REPORT.pdf), and
[all 64 timing cells](cpu-report-20261015/RESULTS.csv) show constructor-only
results. The AMD 224-bit/degree-31 linear-algebra paired median ratio is
16.9189, observed range 16.8594–16.9618. Slower small-degree Atkin cases are
retained. [Immutable CPU artifacts and native audits](cpu-evidence-38766895864/)
record the actual CPUs, binary/source identities and failed ARM admission.

The [fourth-candidate control](rung4-report-20261015/REPORT.md) preserves its
unverified timing hypothesis separately from those earlier AMD results.
[Kernel baseline protocol](KERNEL_BASELINE_PROTOCOL.md) distinguishes supplied
inputs and output materialization. [VALIDATION.md](VALIDATION.md),
[ATTEMPTS.md](ATTEMPTS.md) and [PROTOCOL.md](PROTOCOL.md) retain verification,
infrastructure failures and measurement boundaries. Canonical DLP/rho and
curve charts receive no value from these constructor-only measurements.

The [exploratory family census](family-report-20261015/REPORT.md) audits all
395 legacy records on AMD, ARM and local M4, retaining 394 true flags and one
false path-search flag on every profile. It does not establish a new timing
ratio. The local executable/source/raw/resource receipts are sealed in
`LEGACY_M4_RUN.tar.gz`; the first hosted wave's missing legacy executable
remains an explicit provenance gap. Future waves now retain that binary.
