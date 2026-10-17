# Prospective N83 four-byte ID30 Tier 1 holdout

The 24 public points in `targets.jsonl` were generated before the new
four-byte ID30 producer was timed. Their SHA-256 is
`9382792536e697b8e1799e2b8a9fc7c00f3912fc2d871b409aa93ebbf35d7307`;
the answer sidecar SHA-256 is
`05890e113df6df51d1f76b0d7b741b61d0e7e455e74cf386820c6da4871ccd78`.
The deterministic domain is `n83-indexed4-id30-holdout24-20261017-v1`.
The points are distinct from prior N83 target fixtures checked at sealing.
The curve ID is `EC1N83Ckb1a1haacb69ee9f37`.

Dispatch this holdout only if the frozen 128-point K4096 local screen returns
complete independent base, rank, four-point-relation, scalar and phase replay
for both five-byte ID30 and four-byte ID30; the observed geometric mean
target-online interval for four-byte ID30 must be lower than five-byte ID30.
Otherwise retain the points unmeasured. Exact IC1 candidate IDs remain null
until their bases and executed configuration are materialized and checked.

For a dispatch, use four distinct Tier 1 cloud instances or Pods, each with
six point indices. Run both exact candidates on every point in fresh IC
processes with empty base, index and log caches. Pair each with an
independently solved strong rho instance on the identical public point.
Rotate candidate order and alternate IC/rho order; calibrate rho before and
after each pair. Record provider identity, source/binary hashes, topology,
affinity, cgroup limits, throttling, steal, memory, all target phases, raw
relations, and independent scalar replay. Bound each IC process by 2,400
seconds and rho/calibration by 120 seconds; size the producer/job caps for
all twelve IC processes per host before dispatch.

Admit a host only after every exact replay, source, resource, calibration,
provider and order gate passes. On 24 admitted paired points, the four-byte
online winner requires a geometric five-byte/four-byte online ratio above
one and exact two-sided paired sign-flip p below 0.05. Report all per-point
IC/rho ratios and compare setup-inclusive cold process times separately.
Preserve every timeout, noisy pair and failed host as a censored row.
