# N83 root11 finalist comparison

This experiment compares the two root11 factor bases whose target-online
timings were closest in the intermediate K sweep: K2816 and K3072. It fixes
32 new public points before cloud execution and tests their one-target online
intervals on four separate Tier 1 allocations. The exact candidate IDs are
resolved from their immutable manifests by `arms.py`; the IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.

The field is GF(2^83), the curve is the Koblitz `a=1,b=1` model identified by
`EC1N83Ckb1a1haacb69ee9f37`, and the subgroup order is
`8569786107849059`. K2816 has 467,456 actual usable factor-base points and
2,816 folded columns. K3072 has 509,952 usable points and 3,072 folded
columns. The two arms use the same four-summand indexed root decomposition,
single target thread, relation policy, final algebra and target handling;
only the factor-base size and its exact manifest differ.

`generate.py` fixes 32 distinct points by
`1 + little_endian(SHA256(domain || u32_le(index))) mod (r-1)` with domain
`n83-root11-finalists-32-20261017-v1`. The public targets have SHA-256
`4eb139e73c99f0921041a911d30698997913d1e71155a8e09035ffa647e02a58`;
the answer sidecar has SHA-256
`016da838827c59e79b0b0f2c5296004bc94b73c5f3f3d9a78b4222f826e84139`.
These points do not overlap any of the 55 other `targets.jsonl` files in the
frozen source checkout. The scalar sidecar is used only for independent
replay, outside the timed interval.

Each of four GitHub Actions `ubuntu-24.04` runners receives eight point
indices. Each point has one cold process for each IC arm, a separately seeded
same-point rho solve for each arm, and three fixed-point calibration rho solves
before and after each IC/rho pair. Candidate order rotates with point index;
IC/rho order alternates by point and candidate. The cloud verifier checks
fixture scalars, exact base dumps, decompositions, scalar recovery, command
flags, phase timings, resource receipts, event order, topology, cgroup limits,
steal time and calibration medians. Each before/after median pair must differ
by at most 10% of their mean. An independent replay must reproduce every
analysis from the provider ZIP. The four provider job IDs, runner IDs and
overlapping intervals must bind to the archived host receipts. Every failed,
timed-out, noisy or missing pair stays recorded and prevents a full-panel
selection.

The primary metric is the target-online interval from first target-dependent
work through scalar recovery and verification, with reusable factor-base and
index preparation outside that interval. The supplementary cold metric is
the full process wall time including preparation. `select.py` compares the
32 paired online values with a two-sided exact sign test at alpha 0.05. A
selection requires 23 or more strict wins by one arm, no ties, complete
independent replay and all Tier 1 host gates. It also reports the geometric
online ratio and a paired 95% interval. The same rule is applied separately
to cold process time. The selector never promotes a producer-only timing or
an incomplete panel. Each candidate is paired with a verified rho solve on
the identical public point and resource envelope; per-host rho/IC ratios are
reported separately from factor-base selection.

The dispatcher hashes every input, candidate manifest, source receipt,
verifier, selector, workflow and resource cap before the first cloud run.
Each IC process has a 900-second cap; each rho process has a 120-second cap;
the per-host producer cap is 270 minutes and the workflow cap is 330 minutes.
