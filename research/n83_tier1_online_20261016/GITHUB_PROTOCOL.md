# GitHub-hosted Tier 1 route after the rejected AWS allocation

The [AWS allocation result](RESULT.md) preserved a provider-level block
before instance creation. This route runs the already frozen four-arm
same-instance selector and five-point panel in one `ubuntu-24.04` GitHub
Actions job. Its exact source, five fresh public points, root-index flags,
same-point rho inputs, IC/rho caps, A/A thresholds, acceptance interval
and independent replay rules are those in [PROTOCOL.md](PROTOCOL.md).
The hosted route is separately frozen before its first execution.

The job first builds and tests the exact Rust producers and runs K3072 A/B,
K2543 A/B and same-point rho on the previously exposed point. The selector
pins four distinct physical cores on one visible NUMA node. It records the
host, CPU affinity, cgroup quota/throttle/OOM, steal, resource usage, base
and index identities, source/binary hashes and all raw outcomes. Independent
verification and the strict A/A and separation gate must pass before an
exact IC1 manifest and selection are written. A noisy or incomplete selector
stops before the five fresh public points are loaded.

The passing selector's raw rows, binaries, manifest and selection are
uploaded as a distinct `actions/upload-artifact@v4` object. The next step
requires a nonempty artifact ID and SHA-256 digest, fetches that exact ID
through GitHub's Actions artifact API, and checks its run ID, source commit,
name, digest, size and unexpired state. The resulting seal receipt binds the
selected manifest, analysis, source, host and four CPUs. The panel loader
requires that receipt and verifies the same run identity and CPU set before
reading the five-point fixture. The object is retained for 14 days; its raw
ZIP and local replay receipts must be archived outside GitHub before expiry.

The panel runs the five fresh same-point IC/rho pairs in alternating order,
with calibrations before and after each pair, under the same job and CPU
allocation. `verify.py` checks the hosted run identity, exact selection and
seal, input points, base/relation/scalar replay, telemetry, pair order and
paired confidence interval. The final artifact retains every selector and
panel row, including capped, failed and noisy outcomes. Local collection
must download both exact artifact IDs, verify their digests and contents,
and recreate the cloud analyses byte for byte on the launch commit.

The hardware tier is `tier1_github_hosted_runner` only when all exposed
topology, affinity, cgroup, throttling, steal and calibration gates pass.
Physical host isolation is not required for this tier. The accepted result
must be described with its observed runner envelope; if the runner lacks
four distinct cores or adequate memory, the job records a host-capacity
failure and the online speedup remains unknown. K3114's separate
exposed-point capacity-frontier screen does not change this already frozen
K3072/K2543 selector.
