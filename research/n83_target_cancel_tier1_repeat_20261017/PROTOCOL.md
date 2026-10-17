# N83 cooperative cancellation: fresh Tier 1 confirmation

The two-worker cancellation implementation returns after the first verified
four-point target relation and cooperatively stops the other worker. This
sealed confirmation compares it with two-worker join on 24 new public N83
points. Both candidates use the same 164,672 usable factor-base points,
992 folded columns, root11 index, rank construction, and source SHA-256
`8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829`.
Only `KIC_TARGET_CANCEL` differs. The two versioned candidate manifests
specify their complete identities and flags.

The public target file SHA-256 is
`aa02c031dedd97a5f4bb214838e550e4597f5af44b556c3f24f7e5f3eefd99f5`;
the independent checker's answer sidecar SHA-256 is
`c0134002e32ff898c8799f2024649ef286d0fbf881dc5a1a3eff5b589e4ce7d2`.
The deterministic target domain is
`n83-target-cancel-tier1-repeat24-20261017-v1`. The 24 points are distinct
from the first cancellation panel. The producer sees one public point at a
time; the answer sidecar is used by independent replay.

Four distinct, concurrent cloud instances or Pods each receive six points:
strata 0–3 cover indices 0–5, 6–11, 12–17, and 18–23. Rotate join and
cancellation by point index and alternate IC and independently solved
same-point rho order. For **each candidate and point**, run three
independent solves of a fixed calibration rho point before the IC/rho pair,
then three after. Compare the median of the three before times with the
median of the three after times. The relative difference divided by the
mean of those medians must be at most 10%. All six calibration solves must
independently replay and pass resource checks. This repetition rule and
threshold are frozen before cloud dispatch.

Each IC process starts with empty base, index, and log caches. Each IC run
has a 900-second cap; each rho or calibration run has a 120-second cap; each
host producer has a 300-minute cap. Freeze code, build, commands, flags,
public point, resource envelope, exposed topology, affinity, cgroup limits,
throttling, steal time, provider identity, and raw artifacts. Record every
timeout, failure, and noisy pair. The timing comparison uses the one-target
online interval, including target query, decomposition, relation check,
descent, and producer recovery check. Report cold process cost separately.

Independent replay checks all 48 IC relations and scalars, 48 matching rho
solves, 288 calibration solves, exact commands and event order, source and
binary digests, host/resource receipts, and provider allocation identity.
The 10% median calibration gate and all other Tier 1 gates must pass on all
four strata before the full selector can choose a winner. A stratum with
verified computations but noisy calibration remains `BOUNDED_UNKNOWN` and
`NOT_ADMITTED`; its evidence is preserved. The selector never converts a
noisy host into a failed mathematical replay.

The frozen online success rule requires all four admitted strata, all 24
defined paired ratios, geometric mean join/cancellation above one, and an
exact two-sided paired sign-flip p value below 0.05. Report the p value,
95% descriptive paired log-ratio interval, win count, probes, per-point
same-target rho/IC ratios, and all cold costs. The same exact selector is
applied to cold process cost separately. A fresh confirmation panel is
needed because the first 24-point panel had four calibration pairs over
10% on two hosts despite independently verified point computations.
