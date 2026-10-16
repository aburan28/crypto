# N83 five-point Tier 1 single-target online panel

The [fixture manifest](fixture_manifest.json) freezes five distinct fresh
public points under domain `n83-tier1-online-five-fresh-20261016-v1`.
`targets.jsonl` has SHA-256
`4976679ac6c9e3858cd26312ec2e50b8eb30971041ea38286c9582dd9217a1c1`;
the independent answer sidecar `answers.jsonl` has SHA-256
`41004f94d7257c26a434a291767a06b87d41d590d952e5551edca1aa62c5545d`.
The generator replays the N83 subgroup and excludes all prior public points
in the pinned cold-screen prior manifest, the cold-screen point, and its four
reserve points. No producer receives a scalar or answer sidecar.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059` and cofactor
`1128547018`. The candidate source, exact factor-base manifest, binary,
flags, and rho reference must be frozen in a versioned selection record
before any of these five points is executed. The candidate is selected from
earlier exposed-point and reserve-point screens; an unknown or unverified
selection prevents the panel from starting.

The approved allocation is one `g5g.metal` Spot instance in `us-east-1b`,
for at most 165 minutes and $5 total. Before launch, check current Spot
price and capacity, use a one-time request with a $0.50/hour maximum bid,
cap the root volume at 50 GiB, enable instance-initiated termination, and
start an in-instance termination watchdog from the EC2 pending timestamp
before installing or building. Abort boot if the pending timestamp cannot
be read or the watchdog cannot be installed.
Use a dedicated minimum-permission SSM role and a security group with no
inbound rule. Retain the launch, pricing, boot, termination and cost
receipts. The panel's resource envelope pins the IC and rho processes to
the same four exposed CPUs, retains a 15.25 GiB process RSS guard, and
records topology, affinity, cgroup quota and memory limit, OS, compiler,
source and binary hashes. The bare-metal instance is reported as Tier 1
until the separate Tier 0 isolation and noise requirements are checked.
The bootstrap uploads its complete logs, binaries, producer rows and
verifier output to the dedicated private S3 bucket before shutting down.
`collect.py` archives that exact object locally, replays the analysis on the
same source commit, retains both analysis hashes and requests termination
if the instance remains active. Delete the temporary IAM, network and S3
resources after the local archive is independently checked.

Each point is one primary single-target workload. On points 1, 3 and 5,
run `IC,rho`; on points 2 and 4, run `rho,IC`. Use a fresh process, empty
factor base, index and factor-log system for each IC arm; retain its
target-independent setup time separately. Rho uses a fresh same-point
walk, no shared distinguished-point table and a point-specific pinned
seed. Cap each IC process at 900 seconds and rho at 120 seconds.

Before and after every pair, run the same fresh-process fixed-seed rho
calibration on the already exposed cold-screen point. It is outside both
online intervals. Require its two intervals to differ by at most 10% of
their mean. Preserve all ten calibration receipts. Record for every arm
UTC start/stop, monotonic outer wall, user/system CPU, RSS, major faults,
involuntary context switches, `/proc/stat` steal and interrupt deltas,
cgroup throttle/OOM deltas, CPU frequency where exposed, and source/point
hashes. A missing telemetry field stays explicit. Reject a pair if answer
replay fails, a cap is crossed, cgroup throttling or OOM occurs, steal time
exceeds 1% of the measured arm wall, or calibration fails. Do not replace
excluded points; a retry is a separately archived allocation.

The primary IC online interval starts at the first target-dependent query
after reusable preparation and stops after producer verification. Charge
all target-dependent attempts and report the five exclusive target phases.
The paired rho online interval starts at its first target-dependent walk
step and stops after its producer verification. Independently replay both
answers afterward and record the replay duration outside the producer
intervals. Declare one online speedup only if all five frozen pairs pass,
IC is faster on all five, and the two-sided paired 95% lower bound for the
mean log ratio `log(rho_online_ms / IC_online_ms)` is positive. Report all
five ratios, the median ratio, arithmetic mean intervals, uncertainty,
resource receipts and exact Tier 1 instance identity. Any missing pair or
phase keeps the aggregate result unknown.
