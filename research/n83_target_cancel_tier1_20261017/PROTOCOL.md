# N83 two-worker cancellation: prospective Tier 1 confirmation

The implementation contribution is cooperative cancellation after the first
verified four-point target relation. This protocol compares that policy with
the existing two-worker join policy on 24 fresh N83 public points. The exact
factor base has 164,672 usable points and 992 signed-Frobenius columns.
Both policies use the same source SHA-256
`8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829`,
root11 indexed search, rank construction, and two target workers. Only
`KIC_TARGET_CANCEL` changes between them. The exact candidate manifests are
`candidate_IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h33c9c3b2f933.json`
for join and
`candidate_IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h58ec337f71f7.json`
for cancellation.

The sealed public fixture has SHA-256
`b18623d86a131cef7b40ebd83c0465598f6a55554fa83bfa3d4f5634f34c74ba`.
Its answer sidecar has SHA-256
`270b48e9d60974c3c7c7bb39449ffb8bff17b28f33565432c1f2733c017ebf5a`.
The deterministic domain is `n83-target-cancel-tier1-holdout24-20261017-v1`.
These 24 points are distinct from all prior study `targets.jsonl` files in
the repository at sealing. The sidecar is available only to the independent
checker, and the producer receives only one public point at a time.

Four independent cloud instances or Pods each receive six points: strata
0–3 get point indices 0–5, 6–11, 12–17, and 18–23. Rotate the two policies
by point index. For each policy and point, run a calibration rho solve, one
IC/rho same-point pair in alternating order, and a second calibration.
Start every IC run in a fresh process with empty base, index, and log caches.
Freeze source, binary, commands, environment, point, resource envelope, host
topology and provider job/runner identity in raw receipts. Each IC process
has a 900-second cap; each rho or calibration has a 120-second cap; each host
producer has a 300-minute cap. Preserve every timeout, failed replay and
noisy pair.

The primary comparison is one-target online wall time for cancellation versus
join on each identical point, with the independently solved rho time for that
point alongside both arms. Record exclusive target phases, target probes,
setup-inclusive cold process time, memory, base/index/rank, and independent
group-law scalar replay. Before any timing comparison, require all 48 IC
relations and scalars, all 48 matching rho solves, the before/after
calibrations, exact launch flags and source, provider identity, resource caps,
event order and exposed topology to pass. Tier 1 permits cloud instances or
Pods; report observed affinity, cgroup, throttling, steal and calibration
noise. A failed stratum remains a censored result.

The frozen success rule for the cancellation speedup is: all four strata
admitted; all 24 paired online ratios defined; geometric mean join/cancel
above one; and the exact two-sided paired sign-flip p value below 0.05.
Report the exact p value, 95% descriptive paired log-ratio interval, win
count, probe reduction and all point rows. Apply the same calculation to
cold process cost separately. Rho/IC online ratios are reported by policy
and point, never substituted for the direct cancellation comparison.
