# Prospective K3072 cancellation Tier 1 holdout

The 24 public N83 points in `targets.jsonl` were generated and sealed before
the K3072 cancellation screen. Their SHA-256 is
`432d72bc03faba869d09124628a22b3e0ab5acc692c1e2f09a52a90efb22fb9a`;
the verifier-only answer SHA-256 is
`5fc1b3a259a4fe964c869e3bf0ba8559094933d07032cd9d36cd64c544ce8471`.
The deterministic domain is `n83-k3072-target-cancel-holdout24-20261017-v1`.
The points are distinct from the local 128-point batch and the two other
checked N83 cancellation/index holdouts.

Dispatch only if both frozen local policies complete independent replay and
the paired serial/cancellation geometric online ratio exceeds one. The two
exact IC1 manifests were materialized before the screen: serial
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0h49e5d0fe6271`
and cancellation
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hf5484487a579`.
Their manifest SHA-256 values are recorded in the local `candidate_ids.json`.

For a dispatch, use four distinct Tier 1 cloud allocations with six point
indices each. Run both K3072 policies in fresh IC processes on every point,
with the same public point separately solved by strong rho. Rotate policy
order, alternate IC and rho order, and calibrate rho before and after each
pair. Record provider identity, source and binary hashes, topology, exposed
affinity, cgroup limits, throttling, steal time, memory, calibration noise,
all exclusive target phases, raw relations and independent scalar replay.
Cap each IC process at 1,200 seconds and each rho/calibration at 120 seconds;
bound the whole host job for all twelve IC processes before dispatch.

Admit a host only after its exact replay, source, resource, calibration and
order gates pass. On 24 admitted paired points, require a geometric
serial/cancellation target-online ratio above one and an exact two-sided
paired sign-flip p below 0.05. Report every IC/rho online ratio and the cold
process costs separately. Preserve timeout, noisy and failed pairs as
censored records.
