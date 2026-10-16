# Alternate host isolation preflight, 2026-10-09 UTC

This is a **read-only rejection probe**, separate from the Modal timing panels.
It checks whether an already-configured SSH host can supply the strict
physical-host receipt required by [`../PROTOCOL.md`](../PROTOCOL.md). It did
not build a solver, generate a target, start a timed run, or change host
settings. Its result is not an `isolab` strict receipt.

The existing `runpod` SSH alias returned exit code 255 to a batch-mode
`true` liveness check with a five-second connect timeout. Its quiet-mode
stderr was empty, so the specific connection or authentication failure is
unknown. The existing `runpod-isolated-blush-stork` alias returned exit 0
for the same check. Its [read-only JSON probe](alternate_host_20261009.json),
SHA-256 `654584490fac94416137568abd915b6fade40dd66c791240ae03db297c4fe451`,
was captured at `2026-10-09T03:52:30.245847+00:00`. The
[`host_isolation_preflight.py`](../host_isolation_preflight.py) source has
SHA-256 `fa8e319b4c3918586a2cc9e6a0dc0697026ee78997d5bd4c82be19a04b67cc1f`.

| Visible fact | Probe result | Strict-gate consequence |
| --- | --- | --- |
| CPU affinity and topology | 16 logical CPUs on one NUMA node, eight physical cores with both SMT siblings visible | Fewer than 15 separate cores if sibling-idle placement is required; the visible mask alone cannot certify the protocol's SMT policy |
| cgroups | Docker cgroup v1; cgroup v2 controllers absent | No cgroup v2 isolated partition or exclusive-cpuset readback |
| CPU quota | v1 `cpu.cfs_quota_us = -1` | No quota observed, but this does not establish host exclusivity |
| Hardware counters | `perf` unavailable; `perf_event_paranoid = 4` | Required PMU evidence is unavailable in the visible environment |
| Host-wide routing and contention | Not exposed by this container probe | IRQ routing and absence of foreign runnable work cannot be certified |

The strict preflight **fails**, and `strict_host_receipt` and
`controlled_speedup` remain `null` in the JSON. The DMI string reports a
physical server model, but the process is in a Docker cgroup; that string
does not turn the container's affinity mask into a host-wide isolation
receipt. The CPU/memory PSI values in the probe are point-in-time readbacks,
not a quiet-settle or timed-run interval.

The exact non-mutating commands used from the repository checkout were:

```sh
ssh -q -o BatchMode=yes -o ConnectTimeout=5 -o StrictHostKeyChecking=yes runpod true
ssh -q -o BatchMode=yes -o ConnectTimeout=5 -o StrictHostKeyChecking=yes runpod-isolated-blush-stork true
ssh -q -o BatchMode=yes -o ConnectTimeout=5 -o StrictHostKeyChecking=yes \
  runpod-isolated-blush-stork python3 - \
  < experiments/koblitz-s3-fidelity-20261007/host_isolation_preflight.py \
  > experiments/koblitz-s3-fidelity-20261007/host_preflight/alternate_host_20261009.json
```

A qualified next route is a physical Linux `isolab` worker with at least 15
execution CPUs on one NUMA node, separate housekeeping cores, a cgroup v2
isolated partition, SMT sibling and IRQ control, hardware PMU access, and a
passing quiet A/A pilot. This session has no configured `ISOLAB_URL` or
`isolab` command; no such worker was reached here. If that access is
provided, the read-only `isolab doctor` and strict calibration precede any
new timing run. The confirmatory target count remains unset.
